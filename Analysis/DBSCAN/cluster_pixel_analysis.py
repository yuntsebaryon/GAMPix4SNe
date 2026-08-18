#!/usr/bin/env python3
"""
Pixel-first clustering and R-ratio analysis for waveform-based GAMPixPy samples.

Main change relative to the tile-first notebook:
  old logic: one tile hit defines one analysis object, then nearby pixels are matched;
  new logic: pixel hits are clustered first, then corresponding tile hits are matched.

The script keeps the full uncut cluster table as the primary output.  Optional
quality-cut tables are written as extra files, never as replacements.

Expected hit HDF5 layout:
  top-level datasets/groups: pixels, tiles, optionally meta
  pixel fields: event id/event_id, pixel x, pixel y, start t or trig t, waveform,
                label, attribution
  tile fields:  event id/event_id, tile x, tile y, start t or trig t, waveform,
                label, attribution

Expected truth HDF5 layout:
  dataset: segments
  fields: event_id/event id, segment_id, dE, and either x/y/z or *_start/*_end fields.
  Geometry-based drift length is computed from --detector-yaml drift_volumes.

Example:
  python3 cluster_pixel_first_analysis.py \
    --hit-file ../detsim_sample/gampixpy_fullgeoanatruth-vd-reduced_g4_00_2Mhz_segmentlabel_lowtrig_5mmpitch.h5 \
    --truth-file ../g4_cv1_sample/fullgeoanatruth-vd-reduced_g4_00.h5 \
    --output-dir ../detsim_sample/cluster_pixel_first \
    --max-events 100
"""

from __future__ import annotations

import argparse
import json
import math
import os
import warnings
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import h5py
import numpy as np
import pandas as pd

try:
    from sklearn.cluster import DBSCAN
except Exception as exc:  # pragma: no cover
    raise RuntimeError("This script needs scikit-learn: pip install scikit-learn") from exc

try:
    import matplotlib.pyplot as plt
except Exception:  # pragma: no cover
    plt = None

try:
    from scipy.optimize import curve_fit
except Exception:  # pragma: no cover
    curve_fit = None

try:
    from gampixpy import config as gampix_config
except Exception:  # pragma: no cover
    gampix_config = None

try:
    import torch
except Exception:  # pragma: no cover
    torch = None


INVALID_SEGMENT_ID = -9999


# -----------------------------------------------------------------------------
# Argument handling
# -----------------------------------------------------------------------------
@dataclass
class Config:
    hit_file: str
    truth_file: Optional[str]
    output_dir: str
    tag: str

    max_events: Optional[int]
    max_pixel_hits_per_event: Optional[int]

    pixel_pitch: float
    tile_size_pixels: float
    tick_size: float
    waveform_ticks: int
    pre_samples: int

    cluster_eps_space: float
    cluster_eps_time: float
    dbscan_eps: float
    dbscan_min_samples: int
    min_pixel_charge: float
    min_cluster_pixel_hits: int
    min_cluster_charge: float

    noise: float
    threshold_sigma: float
    anode_x: float
    detector_yaml: Optional[str]
    drift_mode: str
    truth_lookup: str
    truth_chunk_size: int
    make_plots: bool

    fit_waveforms: bool
    plot_fit_examples: bool
    fit_example_clusters: int
    fit_example_hits_per_cluster: int

    @property
    def threshold(self) -> float:
        return self.noise * self.threshold_sigma

    @property
    def tile_size(self) -> float:
        return self.tile_size_pixels * self.pixel_pitch

    @property
    def tile_half_size(self) -> float:
        return 0.5 * self.tile_size

    @property
    def pixel_half_size(self) -> float:
        return 0.5 * self.pixel_pitch

    @property
    def readout_window(self) -> float:
        return self.waveform_ticks * self.tick_size


def parse_args() -> Config:
    ap = argparse.ArgumentParser(
        description="Cluster pixel hits first, then match tile hits and truth segments."
    )
    ap.add_argument("--hit-file", required=True, help="GAMPixPy hit/readout HDF5 file")
    ap.add_argument("--truth-file", default=None, help="edep-sim truth HDF5 file with segments")
    ap.add_argument("--output-dir", required=True, help="Directory for CSV and plots")
    ap.add_argument("--tag", default="pixel_first", help="Prefix for output files")

    ap.add_argument("--max-events", type=int, default=None, help="Limit number of events for testing")
    ap.add_argument("--max-pixel-hits-per-event", type=int, default=None,
                    help="Optional safety cap; keeps largest-charge pixel hits per event")

    ap.add_argument("--pixel-pitch", type=float, default=5.0, help="Pixel pitch in spatial units")
    ap.add_argument("--tile-size-pixels", type=float, default=20.0,
                    help="Tile side length in units of pixel pitch; default 20 means 20x20 pixels")
    ap.add_argument("--tick-size", type=float, default=0.5,
                    help="Time spacing between waveform samples")
    ap.add_argument("--waveform-ticks", type=int, default=20,
                    help="Number of waveform ticks after activation")
    ap.add_argument("--pre-samples", type=int, default=0,
                    help="Number of waveform samples before the stored start/trig time, if applicable")

    ap.add_argument("--cluster-eps-space", type=float, default=1.5,
                    help="DBSCAN distance scale in pixel-pitch units")
    ap.add_argument("--cluster-eps-time", type=float, default=4.0,
                    help="DBSCAN distance scale in time units")
    ap.add_argument("--dbscan-eps", type=float, default=1.0,
                    help="DBSCAN eps after coordinate normalization")
    ap.add_argument("--dbscan-min-samples", type=int, default=2,
                    help="DBSCAN min_samples; use 1 to keep isolated pixel hits as clusters")
    ap.add_argument("--min-pixel-charge", type=float, default=0.0,
                    help="Only cluster pixel hits with summed waveform charge above this value")
    ap.add_argument("--min-cluster-pixel-hits", type=int, default=1)
    ap.add_argument("--min-cluster-charge", type=float, default=0.0)

    ap.add_argument("--noise", type=float, default=50.0, help="Noise value used for threshold definitions")
    ap.add_argument("--threshold-sigma", type=float, default=3.0,
                    help="Threshold is noise * threshold_sigma")
    ap.add_argument("--anode-x", type=float, default=325.0,
                    help="For diagnostic compatibility with notebook drift = anode_x - truth_x")
    ap.add_argument("--detector-yaml",
                    default="/home/yboxun/NeutrinoGAMPix/detsim_prediction/depth/far_detector_vd.yaml",
                    help="Detector YAML used to compute geometry-based drift length")
    ap.add_argument("--drift-mode", default="signed_volume0",
                    choices=["signed_volume0", "min_abs", "abs_volume0"],
                    help=("Geometry drift convention. signed_volume0 exactly follows the older "
                          "script formula dot(segment_midpoint - anode_center, -drift_axis) "
                          "using drift_volumes/volume_0. No abs() correction is applied."))
    ap.add_argument("--truth-lookup", default="segment_id",
                    choices=["segment_id", "event_segment", "event_segment_then_segment_id"],
                    help=("How detector labels are matched to truth segments. segment_id reproduces "
                          "the older script: segment_drift_map[segment_id]. event_segment uses "
                          "(event_id, segment_id)."))
    ap.add_argument("--truth-chunk-size", type=int, default=250000,
                    help="Number of truth segment rows to read at a time. Smaller uses less memory.")
    ap.add_argument("--make-plots", action="store_true", help="Write basic diagnostic plots")
    ap.add_argument("--fit-waveforms", action="store_true",
                    help=("Fit individual pixel waveforms and summed cluster waveforms with a Gaussian. "
                          "This enables Gaussian-fit charge columns, including pixel_fit_charge_3sd."))
    ap.add_argument("--plot-fit-examples", action="store_true",
                    help=("When --fit-waveforms is enabled, save example individual-pixel Gaussian-fit plots "
                          "for the first few clusters."))
    ap.add_argument("--fit-example-clusters", type=int, default=2,
                    help="Number of earliest clusters for which to save example fit plots")
    ap.add_argument("--fit-example-hits-per-cluster", type=int, default=10,
                    help="Maximum number of individual pixel-hit fit plots to save per example cluster")

    a = ap.parse_args()
    return Config(**vars(a))


# -----------------------------------------------------------------------------
# HDF5 and field-name helpers
# -----------------------------------------------------------------------------
def names(ds: h5py.Dataset) -> Tuple[str, ...]:
    if ds.dtype.names is None:
        raise ValueError(f"Dataset {ds.name} does not have named fields")
    return tuple(ds.dtype.names)


def pick_name(ds: h5py.Dataset, candidates: Sequence[str]) -> str:
    available = set(names(ds))
    for c in candidates:
        if c in available:
            return c
    raise KeyError(f"Could not find any of {candidates} in {ds.name}. Available: {sorted(available)}")


def read_field(ds: h5py.Dataset, candidates: Sequence[str], idx: Optional[np.ndarray] = None):
    name = pick_name(ds, candidates)
    arr = ds[name]
    return arr[:] if idx is None else arr[idx]


def event_slices(event_ids_sorted: np.ndarray) -> Iterable[Tuple[int, int, int]]:
    if len(event_ids_sorted) == 0:
        return
    starts = np.r_[0, np.nonzero(event_ids_sorted[1:] != event_ids_sorted[:-1])[0] + 1]
    stops = np.r_[starts[1:], len(event_ids_sorted)]
    for s, e in zip(starts, stops):
        yield int(event_ids_sorted[s]), int(s), int(e)


def safe_array(a, dtype=float):
    return np.asarray(a, dtype=dtype)


# -----------------------------------------------------------------------------
# Charge and waveform helpers
# -----------------------------------------------------------------------------
def waveform_sum(wf: np.ndarray) -> np.ndarray:
    wf = np.asarray(wf, dtype=float)
    if wf.ndim == 1:
        return np.asarray([float(np.nansum(wf))])
    return np.nansum(wf, axis=1)


def waveform_sum_above(wf: np.ndarray, threshold: float, subtract_threshold: bool = False) -> np.ndarray:
    wf = np.asarray(wf, dtype=float)
    y = np.maximum(wf - threshold, 0.0) if subtract_threshold else np.where(wf > threshold, wf, 0.0)
    if y.ndim == 1:
        return np.asarray([float(np.nansum(y))])
    return np.nansum(y, axis=1)


def waveform_active_mask(wf: np.ndarray) -> np.ndarray:
    wf = np.asarray(wf, dtype=float)
    if wf.ndim == 1:
        wf = wf[None, :]
    return np.isfinite(wf).any(axis=1) & (np.nansum(np.abs(wf), axis=1) > 0)


def sum_waveforms(wf: np.ndarray) -> np.ndarray:
    wf = np.asarray(wf, dtype=float)
    if wf.size == 0:
        return np.array([], dtype=float)
    if wf.ndim == 1:
        return wf.astype(float)
    return np.nansum(wf, axis=0)


def weighted_mean(x: np.ndarray, w: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    w = np.asarray(w, dtype=float)
    good = np.isfinite(x) & np.isfinite(w) & (w > 0)
    if not np.any(good):
        return np.nan
    return float(np.sum(x[good] * w[good]) / np.sum(w[good]))


def weighted_rms(x: np.ndarray, w: np.ndarray) -> float:
    m = weighted_mean(x, w)
    if not np.isfinite(m):
        return np.nan
    x = np.asarray(x, dtype=float)
    w = np.asarray(w, dtype=float)
    good = np.isfinite(x) & np.isfinite(w) & (w > 0)
    return float(np.sqrt(np.sum(w[good] * (x[good] - m) ** 2) / np.sum(w[good])))


def span(x: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    return float(np.max(x) - np.min(x)) if len(x) else np.nan


def gaussian_no_offset(t, amp, mu, sigma):
    sigma = np.maximum(sigma, 1e-12)
    return amp * np.exp(-0.5 * ((t - mu) / sigma) ** 2)


def gaussian_charge_total(amp: float, sigma: float, tick_size: float) -> float:
    """Integral of A exp(-0.5*((t-mu)/sigma)^2), divided by tick size."""
    if not (np.isfinite(amp) and np.isfinite(sigma)) or amp <= 0 or sigma <= 0:
        return np.nan
    return float(amp * abs(sigma) * math.sqrt(2.0 * math.pi) / tick_size)


def gaussian_charge_above_threshold(
    amp: float,
    sigma: float,
    threshold: float,
    tick_size: float,
    subtract_threshold: bool = False,
) -> float:
    """
    Gaussian charge in the region where the fitted signal is above threshold.

    If subtract_threshold=False, this is
        integral Gaussian(t) dt over Gaussian(t) > threshold
    divided by tick size.  This matches the raw thresholded waveform convention
    used in pixel_charge_3sd: sum waveform bins whose value is above threshold.

    If subtract_threshold=True, this is
        integral max(Gaussian(t) - threshold, 0) dt
    divided by tick size.
    """
    if not (np.isfinite(amp) and np.isfinite(sigma) and np.isfinite(threshold)):
        return np.nan
    amp = float(amp)
    sigma = abs(float(sigma))
    threshold = float(threshold)
    if amp <= 0 or sigma <= 0:
        return np.nan
    if threshold <= 0:
        return gaussian_charge_total(amp, sigma, tick_size)
    if amp <= threshold:
        return 0.0

    half_width = sigma * math.sqrt(2.0 * math.log(amp / threshold))
    area_inside = (
        amp * sigma * math.sqrt(2.0 * math.pi)
        * math.erf(half_width / (math.sqrt(2.0) * sigma))
    )
    if subtract_threshold:
        area_inside -= 2.0 * half_width * threshold
        area_inside = max(area_inside, 0.0)
    return float(area_inside / tick_size)


def fit_summed_waveform(wf_sum: np.ndarray, cfg: Config) -> Dict[str, float]:
    out = {
        "fit_amp": np.nan,
        "fit_mu": np.nan,
        "fit_sigma": np.nan,
        "fit_charge": np.nan,
        "fit_charge_3sd": np.nan,
        "fit_charge_above_threshold_subtracted": np.nan,
        "fit_rmse": np.nan,
        "fit_success": False,
    }
    y = np.asarray(wf_sum, dtype=float)
    if curve_fit is None or len(y) < 4 or not np.isfinite(y).any() or np.nanmax(y) <= 0:
        return out
    t = (np.arange(len(y), dtype=float) - cfg.pre_samples) * cfg.tick_size
    amp0 = float(np.nanmax(y))
    mu0 = weighted_mean(t, np.maximum(y, 0.0))
    sig0 = weighted_rms(t, np.maximum(y, 0.0))
    if not np.isfinite(mu0):
        mu0 = float(t[np.nanargmax(y)])
    if not np.isfinite(sig0) or sig0 <= 0:
        sig0 = max(cfg.tick_size, 1e-3)
    try:
        popt, _ = curve_fit(
            gaussian_no_offset, t, y,
            p0=[amp0, mu0, sig0],
            bounds=([0.0, float(np.min(t) - cfg.readout_window), 1e-9],
                    [np.inf, float(np.max(t) + cfg.readout_window), np.inf]),
            maxfev=10000,
        )
        amp, mu, sigma = [float(v) for v in popt]
        yhat = gaussian_no_offset(t, amp, mu, sigma)
        sigma_abs = abs(sigma)
        charge = gaussian_charge_total(amp, sigma_abs, cfg.tick_size)
        charge_3sd = gaussian_charge_above_threshold(
            amp, sigma_abs, cfg.threshold, cfg.tick_size, subtract_threshold=False
        )
        charge_sub = gaussian_charge_above_threshold(
            amp, sigma_abs, cfg.threshold, cfg.tick_size, subtract_threshold=True
        )
        out.update({
            "fit_amp": amp,
            "fit_mu": mu,
            "fit_sigma": sigma_abs,
            "fit_charge": float(charge),
            "fit_charge_3sd": float(charge_3sd),
            "fit_charge_above_threshold_subtracted": float(charge_sub),
            "fit_rmse": float(np.sqrt(np.mean((y - yhat) ** 2))),
            "fit_success": True,
        })
    except Exception:
        pass
    return out


def fit_pixel_waveforms_individually(waveforms: np.ndarray, cfg: Config) -> Dict[str, object]:
    """
    Fit each selected pixel-hit waveform individually and sum the fitted charges.

    The main new charge for this request is pixel_fit_charge_3sd:
      sum over pixel hits of the fitted Gaussian integral in the region where
      Gaussian(t) > 3*sigma_noise.

    This can be slower than simple waveform sums, so it is only run when
    --fit-waveforms is enabled.
    """
    wf = np.asarray(waveforms, dtype=float)
    if wf.size == 0:
        return {
            "individual_fit_charge": np.nan,
            "individual_fit_charge_3sd": np.nan,
            "individual_fit_charge_above_threshold_subtracted": np.nan,
            "individual_fit_success_count": 0,
            "individual_fit_fail_count": 0,
            "individual_fit_results": [],
        }
    if wf.ndim == 1:
        wf = wf[None, :]

    results: List[Dict[str, float]] = []
    total_charge = 0.0
    total_charge_3sd = 0.0
    total_charge_sub = 0.0
    success = 0

    for row in wf:
        fit = fit_summed_waveform(row, cfg)
        results.append(fit)
        if fit.get("fit_success", False):
            success += 1
            if np.isfinite(fit.get("fit_charge", np.nan)):
                total_charge += float(fit["fit_charge"])
            if np.isfinite(fit.get("fit_charge_3sd", np.nan)):
                total_charge_3sd += float(fit["fit_charge_3sd"])
            if np.isfinite(fit.get("fit_charge_above_threshold_subtracted", np.nan)):
                total_charge_sub += float(fit["fit_charge_above_threshold_subtracted"])

    return {
        "individual_fit_charge": float(total_charge) if success > 0 else np.nan,
        "individual_fit_charge_3sd": float(total_charge_3sd) if success > 0 else np.nan,
        "individual_fit_charge_above_threshold_subtracted": float(total_charge_sub) if success > 0 else np.nan,
        "individual_fit_success_count": int(success),
        "individual_fit_fail_count": int(len(wf) - success),
        "individual_fit_results": results,
    }


def save_pixel_fit_example_plots(
    cluster_id: int,
    event_id: int,
    pixel_waveforms: np.ndarray,
    individual_fit_results: List[Dict[str, float]],
    cfg: Config,
) -> None:
    """
    Save a small number of individual pixel-waveform Gaussian-fit plots.

    This intentionally only saves examples for the earliest clusters, because
    plotting every pixel hit would produce too many files.
    """
    if plt is None or not cfg.plot_fit_examples or not cfg.fit_waveforms:
        return
    if cluster_id >= cfg.fit_example_clusters:
        return

    wf = np.asarray(pixel_waveforms, dtype=float)
    if wf.ndim == 1:
        wf = wf[None, :]

    n_plot = min(len(wf), cfg.fit_example_hits_per_cluster)
    if n_plot <= 0:
        return

    plot_dir = os.path.join(cfg.output_dir, "plots", "gaussian_fit_examples")
    os.makedirs(plot_dir, exist_ok=True)

    t = (np.arange(wf.shape[1], dtype=float) - cfg.pre_samples) * cfg.tick_size

    for i in range(n_plot):
        y = wf[i]
        fit = individual_fit_results[i] if i < len(individual_fit_results) else {}

        plt.figure(figsize=(8, 5))
        plt.plot(t, y, marker="o", linestyle="-", label="pixel waveform")
        plt.axhline(cfg.threshold, linestyle="--", linewidth=1.2, label=f"3sigma threshold = {cfg.threshold:g}")

        if fit.get("fit_success", False):
            amp = float(fit["fit_amp"])
            mu = float(fit["fit_mu"])
            sigma = float(fit["fit_sigma"])
            tt = np.linspace(float(np.min(t)), float(np.max(t)), 300)
            yy = gaussian_no_offset(tt, amp, mu, sigma)
            plt.plot(tt, yy, linewidth=1.5, label="Gaussian fit")
            title = (
                f"cluster {cluster_id}, event {event_id}, pixel hit {i}\\n"
                f"Q_fit={fit.get('fit_charge', np.nan):.3g}, "
                f"Q_fit_3sd={fit.get('fit_charge_3sd', np.nan):.3g}, "
                f"sigma={sigma:.3g}"
            )
        else:
            title = f"cluster {cluster_id}, event {event_id}, pixel hit {i}\\nGaussian fit failed"

        plt.xlabel("time")
        plt.ylabel("charge")
        plt.title(title)
        plt.legend()
        plt.tight_layout()

        out = os.path.join(
            plot_dir,
            f"cluster_{cluster_id:04d}_event_{event_id}_pixelhit_{i:02d}_gaussfit.png",
        )
        plt.savefig(out, dpi=160)
        plt.close()


# -----------------------------------------------------------------------------
# Truth handling
# -----------------------------------------------------------------------------
def _to_numpy(x) -> np.ndarray:
    """Convert detector-config values, including torch tensors, to numpy arrays."""
    if torch is not None and torch.is_tensor(x):
        x = x.detach().cpu().numpy()
    return np.asarray(x, dtype=float)


def load_drift_volumes(detector_yaml: Optional[str]) -> List[Dict[str, object]]:
    """
    Load drift-volume geometry from the detector YAML.

    Each returned entry contains anode_center and drift_axis.  The drift distance
    is later computed by projecting the segment midpoint onto the axis normal to
    the anode, rather than using the notebook's hard-coded 325 - x expression.
    """
    if detector_yaml is None or detector_yaml == "":
        warnings.warn("No detector YAML given; geometry drift columns will be NaN")
        return []
    if not os.path.exists(detector_yaml):
        warnings.warn(f"Detector YAML not found: {detector_yaml}; geometry drift columns will be NaN")
        return []
    if gampix_config is None:
        warnings.warn("Could not import gampixpy.config; geometry drift columns will be NaN")
        return []

    detector_config = gampix_config.DetectorConfig(detector_yaml)
    try:
        drift_volumes = detector_config["drift_volumes"]
    except Exception as exc:
        warnings.warn(f"Could not read drift_volumes from {detector_yaml}: {exc}; geometry drift columns will be NaN")
        return []

    out: List[Dict[str, object]] = []
    for name, vol in drift_volumes.items():
        try:
            anode_center = _to_numpy(vol["anode_center"])
            drift_axis = _to_numpy(vol["drift_axis"])
            norm = float(np.linalg.norm(drift_axis))
            if not np.isfinite(norm) or norm <= 0:
                warnings.warn(f"Skipping drift volume {name!r}: invalid drift_axis={drift_axis}")
                continue
            drift_axis = drift_axis / norm
            out.append({
                "name": str(name),
                "anode_center": anode_center,
                "drift_axis": drift_axis,
            })
        except Exception as exc:
            warnings.warn(f"Skipping drift volume {name!r}: {exc}")
    return out


def geometry_drift_from_position(
    pos: np.ndarray,
    drift_volumes: Sequence[Dict[str, object]],
    mode: str,
) -> Tuple[float, float, str]:
    """
    Return geometry-based drift value for one truth-segment midpoint.

    The default signed_volume0 is intentionally the same formula used in the
    older script supplied by the user:

        drift = dot(segment_midpoint - anode_center, -drift_axis)

    where anode_center and drift_axis are taken from
    detector_config["drift_volumes"]["volume_0"].

    No abs() is applied in signed_volume0.  abs_volume0 and min_abs are kept only
    as diagnostic options, not as the default.
    """
    pos = np.asarray(pos, dtype=float)
    if len(drift_volumes) == 0 or pos.shape[0] != 3 or not np.all(np.isfinite(pos)):
        return np.nan, np.nan, ""

    signed_values: List[float] = []
    volume_names: List[str] = []
    for vol in drift_volumes:
        center = np.asarray(vol["anode_center"], dtype=float)
        axis = np.asarray(vol["drift_axis"], dtype=float)
        signed = float(np.dot(pos - center, -axis))
        signed_values.append(signed)
        volume_names.append(str(vol["name"]))

    signed_arr = np.asarray(signed_values, dtype=float)
    if mode == "signed_volume0":
        return signed_values[0], signed_values[0], volume_names[0]
    if mode == "abs_volume0":
        return abs(signed_values[0]), signed_values[0], volume_names[0]
    if mode == "min_abs":
        i = int(np.nanargmin(np.abs(signed_arr)))
        return abs(signed_values[i]), signed_values[i], volume_names[i]
    raise ValueError(f"Unknown drift_mode={mode!r}")


def load_truth_segments(
    path: Optional[str],
    cfg: Config,
    wanted_event_ids: Optional[set[int]] = None,
    wanted_segment_ids: Optional[set[int]] = None,
) -> Dict[Tuple[int, int], Dict[str, float]]:
    """
    Load only the truth segments that can be used by the selected detector events.

    Earlier versions read the entire truth file before processing detector events.
    That is wasteful for --max-events tests.  This version still scans the truth
    file in chunks, but only stores rows whose segment_id appears in the selected
    pixel labels.  For --truth-lookup event_segment it also requires the event_id
    to be one of the detector events being processed.

    With the default --truth-lookup segment_id, the stored (-1, segment_id)
    record reproduces the legacy script's dict(zip(segment_id, drift)) behavior:
    if a segment_id appears multiple times, the last matching truth row wins.
    """
    if path is None:
        return {}
    if not os.path.exists(path):
        warnings.warn(f"Truth file not found: {path}; truth columns will be NaN")
        return {}

    if wanted_segment_ids is not None:
        wanted_segment_ids = {int(s) for s in wanted_segment_ids if int(s) != INVALID_SEGMENT_ID}
        if len(wanted_segment_ids) == 0:
            print("No valid segment labels found in selected detector events; truth columns will be NaN")
            return {}

    if wanted_event_ids is not None:
        wanted_event_ids = {int(e) for e in wanted_event_ids}

    drift_volumes = load_drift_volumes(cfg.detector_yaml)
    if drift_volumes:
        print(f"Loaded {len(drift_volumes)} drift volume(s) from {cfg.detector_yaml}")
        print(f"Geometry drift mode: {cfg.drift_mode}")
    else:
        print("No drift-volume geometry loaded; geometry drift columns will be NaN")

    out: Dict[Tuple[int, int], Dict[str, float]] = {}
    n_scanned = 0
    n_kept_rows = 0

    with h5py.File(path, "r") as f:
        if "segments" not in f:
            warnings.warn(f"No 'segments' dataset in {path}; truth columns will be NaN")
            return {}

        S = f["segments"]
        nseg = len(S)
        event_name = pick_name(S, ["event_id", "event id"])
        seg_name = pick_name(S, ["segment_id", "segment id"])
        dE_name = pick_name(S, ["dE", "de", "energy", "E"])
        fields = names(S)

        chunk = max(1, int(cfg.truth_chunk_size))
        for i0 in range(0, nseg, chunk):
            i1 = min(i0 + chunk, nseg)
            sl = slice(i0, i1)

            ev = S[event_name][sl].astype(int)
            sid = S[seg_name][sl].astype(int)

            keep = np.ones(len(sid), dtype=bool)
            if wanted_segment_ids is not None:
                keep &= np.isin(sid, list(wanted_segment_ids))
            if cfg.truth_lookup == "event_segment" and wanted_event_ids is not None:
                keep &= np.isin(ev, list(wanted_event_ids))
            elif cfg.truth_lookup == "event_segment_then_segment_id" and wanted_event_ids is not None:
                # Keep all wanted segment IDs so the segment_id fallback remains possible.
                # Also store event-specific records for the selected detector events.
                pass

            n_scanned += len(sid)
            if not np.any(keep):
                continue

            ev_k = ev[keep]
            sid_k = sid[keep]
            dE_k = S[dE_name][sl].astype(float)[keep]

            if "x" in fields:
                x_k = S["x"][sl].astype(float)[keep]
            elif "x_start" in fields and "x_end" in fields:
                x_k = 0.5 * (S["x_start"][sl].astype(float)[keep] + S["x_end"][sl].astype(float)[keep])
            else:
                x_k = np.full(len(sid_k), np.nan)

            if "y" in fields:
                y_k = S["y"][sl].astype(float)[keep]
            elif "y_start" in fields and "y_end" in fields:
                y_k = 0.5 * (S["y_start"][sl].astype(float)[keep] + S["y_end"][sl].astype(float)[keep])
            else:
                y_k = np.full(len(sid_k), np.nan)

            if "z" in fields:
                z_k = S["z"][sl].astype(float)[keep]
            elif "z_start" in fields and "z_end" in fields:
                z_k = 0.5 * (S["z_start"][sl].astype(float)[keep] + S["z_end"][sl].astype(float)[keep])
            else:
                z_k = np.full(len(sid_k), np.nan)

            for e, s, de, xx, yy, zz in zip(ev_k, sid_k, dE_k, x_k, y_k, z_k):
                pos = np.asarray([xx, yy, zz], dtype=float)
                drift_geom, drift_geom_signed, drift_volume = geometry_drift_from_position(
                    pos, drift_volumes, cfg.drift_mode
                )
                rec = {
                    "event_id": int(e),
                    "segment_id": int(s),
                    "dE": float(de),
                    "x": float(xx),
                    "y": float(yy),
                    "z": float(zz),
                    "drift_325_minus_x": float(cfg.anode_x - xx) if np.isfinite(xx) else np.nan,
                    "drift_geometry": float(drift_geom) if np.isfinite(drift_geom) else np.nan,
                    "drift_geometry_signed": float(drift_geom_signed) if np.isfinite(drift_geom_signed) else np.nan,
                    "drift_volume": drift_volume,
                }

                # Event-specific lookup.
                if cfg.truth_lookup in ("event_segment", "event_segment_then_segment_id"):
                    if wanted_event_ids is None or int(e) in wanted_event_ids:
                        out[(int(e), int(s))] = rec

                # Legacy segment_id-only lookup.  Last matching row wins,
                # exactly like dict(zip(s_ids, segment_drift_values)).
                if cfg.truth_lookup in ("segment_id", "event_segment_then_segment_id"):
                    out[(-1, int(s))] = rec

                n_kept_rows += 1

            if wanted_segment_ids is not None:
                print(f"  truth scan {i1}/{nseg}: kept rows so far = {n_kept_rows}, stored records = {len(out)}", end="\r")

    if wanted_segment_ids is not None:
        print()
    print(f"Truth segments scanned: {n_scanned}; matching truth rows kept: {n_kept_rows}")
    return out


def get_truth_record(
    event_id: int,
    segment_id: int,
    truth: Dict[Tuple[int, int], Dict[str, float]],
    cfg: Config,
) -> Optional[Dict[str, float]]:
    """Look up one truth segment according to the requested lookup convention."""
    sid = int(segment_id)
    if cfg.truth_lookup == "segment_id":
        return truth.get((-1, sid))
    if cfg.truth_lookup == "event_segment":
        return truth.get((int(event_id), sid))
    if cfg.truth_lookup == "event_segment_then_segment_id":
        rec = truth.get((int(event_id), sid))
        return rec if rec is not None else truth.get((-1, sid))
    raise ValueError(f"Unknown truth_lookup={cfg.truth_lookup!r}")


def normalize_label_attr(label: np.ndarray, attribution: np.ndarray, n_ticks: int) -> Tuple[np.ndarray, np.ndarray]:
    """Return labels shape (n_hit, n_label), attr shape (n_hit, n_tick, n_label)."""
    lab = np.asarray(label)
    attr = np.asarray(attribution, dtype=float)

    if lab.ndim == 1:
        lab = lab[:, None]
    elif lab.ndim > 2:
        # Rare defensive case: flatten all non-hit axes except tick-like axis is not expected for labels.
        lab = lab.reshape((lab.shape[0], -1))

    if attr.ndim == 1:
        attr = attr[:, None, None]
    elif attr.ndim == 2:
        # Could be (n_hit, n_label) or (n_hit, n_tick).  Prefer n_label if it matches labels.
        if attr.shape[1] == lab.shape[1]:
            attr = np.repeat(attr[:, None, :], n_ticks, axis=1)
        else:
            attr = attr[:, :, None]
    elif attr.ndim >= 3:
        attr = attr.reshape((attr.shape[0], attr.shape[1], -1))

    # Harmonize the number of labels.
    n_label = min(lab.shape[1], attr.shape[2])
    lab = lab[:, :n_label]
    attr = attr[:, :, :n_label]
    return lab.astype(int), attr.astype(float)


def segment_weights_from_hits(
    labels: np.ndarray,
    attribution: np.ndarray,
    waveforms: np.ndarray,
    hit_indices: np.ndarray,
) -> Dict[int, float]:
    """
    Charge*attribution segment weights for selected channel-level hits.

    IMPORTANT:
      Only positive attribution entries are counted.  This means that if one hit has

          label       = [1,   2,   3]
          attribution = [0.5, 0.5, 0]       or [0.5, 0.5, -1] or [0.5, 0.5, -999]

      then only labels 1 and 2 contribute to the cluster truth weight.  The label
      paired with attribution <= 0 is treated as empty/sentinel information and
      is not allowed to enter the drift-length average.
    """
    seg_w: Dict[int, float] = defaultdict(float)
    if len(hit_indices) == 0:
        return seg_w

    wf = np.asarray(waveforms, dtype=float)[hit_indices]
    n_ticks = wf.shape[1] if wf.ndim == 2 else 1
    lab, attr = normalize_label_attr(
        np.asarray(labels)[hit_indices],
        np.asarray(attribution)[hit_indices],
        n_ticks,
    )

    if wf.ndim == 1:
        wf = wf[:, None]

    ticks = min(wf.shape[1], attr.shape[1])
    wf = wf[:, :ticks]
    attr = attr[:, :ticks, :]

    # Accumulate segment weights using only:
    #   valid segment id, finite positive waveform charge, finite positive attribution.
    # Negative or zero attribution values are excluded before the sum, rather than
    # being allowed to reduce/cancel another valid contribution.
    for k in range(lab.shape[1]):
        ids = lab[:, k]

        for i, sid in enumerate(ids):
            sid = int(sid)
            if sid == INVALID_SEGMENT_ID:
                continue

            valid_tick = (
                np.isfinite(wf[i, :])
                & np.isfinite(attr[i, :, k])
                & (wf[i, :] > 0.0)
                & (attr[i, :, k] > 0.0)
            )

            if not np.any(valid_tick):
                continue

            w = float(np.sum(wf[i, valid_tick] * attr[i, valid_tick, k]))
            if np.isfinite(w) and w > 0.0:
                seg_w[sid] += w

    return seg_w


def _weighted_avg_from_records(
    event_id: int,
    seg_w: Dict[int, float],
    truth: Dict[Tuple[int, int], Dict[str, float]],
    value_key: str,
    weight_kind: str,
    cfg: Config,
) -> Tuple[float, float]:
    """
    Weighted average over matched truth records.

    weight_kind='attr' uses the detector attribution weights accumulated from
    pixel charge.  weight_kind='dE' uses the truth segment dE for every segment
    that appears in the cluster labels.
    """
    num = 0.0
    den = 0.0
    for sid, attr_w in seg_w.items():
        rec = get_truth_record(event_id, int(sid), truth, cfg)
        if rec is None:
            continue
        val = rec.get(value_key, np.nan)
        if not np.isfinite(val):
            continue
        if weight_kind == "attr":
            w = float(attr_w)
        elif weight_kind == "dE":
            w = float(rec.get("dE", np.nan))
        else:
            raise ValueError(weight_kind)
        if np.isfinite(w) and w > 0:
            num += w * val
            den += w
    return (float(num / den), float(den)) if den > 0 else (np.nan, np.nan)


def truth_summary(event_id: int, seg_w: Dict[int, float], truth: Dict[Tuple[int, int], Dict[str, float]], cfg: Config) -> Dict[str, object]:
    empty = {
        "n_truth_segments": 0,
        "dominant_segment_id": INVALID_SEGMENT_ID,
        "dominant_segment_weight": np.nan,
        "truth_weight_sum": 0.0,
        "truth_dE_sum": np.nan,
        "truth_x_attr_avg": np.nan,
        "truth_y_attr_avg": np.nan,
        "truth_z_attr_avg": np.nan,
        "truth_x_dE_avg": np.nan,
        "truth_y_dE_avg": np.nan,
        "truth_z_dE_avg": np.nan,
        # Main drift columns are now geometry-based.
        "drift_length_attr_avg": np.nan,
        "drift_length_dE_avg": np.nan,
        "drift_geometry_attr_avg": np.nan,
        "drift_geometry_dE_avg": np.nan,
        # Diagnostic old notebook-style drift columns.
        "drift_325_minus_x_attr_avg": np.nan,
        "drift_325_minus_x_dE_avg": np.nan,
        "drift_geometry_signed_attr_avg": np.nan,
        "drift_geometry_signed_dE_avg": np.nan,
        "segment_ids": "",
        "segment_weights": "",
    }
    if not seg_w:
        return empty

    dominant_segment_id, dominant_weight = max(seg_w.items(), key=lambda kv: kv[1])
    truth_weight_sum = float(sum(seg_w.values()))

    x_attr, attr_den = _weighted_avg_from_records(event_id, seg_w, truth, "x", "attr", cfg)
    y_attr, _ = _weighted_avg_from_records(event_id, seg_w, truth, "y", "attr", cfg)
    z_attr, _ = _weighted_avg_from_records(event_id, seg_w, truth, "z", "attr", cfg)
    x_de, de_den = _weighted_avg_from_records(event_id, seg_w, truth, "x", "dE", cfg)
    y_de, _ = _weighted_avg_from_records(event_id, seg_w, truth, "y", "dE", cfg)
    z_de, _ = _weighted_avg_from_records(event_id, seg_w, truth, "z", "dE", cfg)

    drift_geom_attr, _ = _weighted_avg_from_records(event_id, seg_w, truth, "drift_geometry", "attr", cfg)
    drift_geom_de, _ = _weighted_avg_from_records(event_id, seg_w, truth, "drift_geometry", "dE", cfg)
    drift_325_attr, _ = _weighted_avg_from_records(event_id, seg_w, truth, "drift_325_minus_x", "attr", cfg)
    drift_325_de, _ = _weighted_avg_from_records(event_id, seg_w, truth, "drift_325_minus_x", "dE", cfg)
    drift_signed_attr, _ = _weighted_avg_from_records(event_id, seg_w, truth, "drift_geometry_signed", "attr", cfg)
    drift_signed_de, _ = _weighted_avg_from_records(event_id, seg_w, truth, "drift_geometry_signed", "dE", cfg)

    return {
        "n_truth_segments": int(len(seg_w)),
        "dominant_segment_id": int(dominant_segment_id),
        "dominant_segment_weight": float(dominant_weight),
        "truth_weight_sum": truth_weight_sum,
        "truth_dE_sum": float(de_den) if np.isfinite(de_den) else np.nan,
        "truth_x_attr_avg": x_attr,
        "truth_y_attr_avg": y_attr,
        "truth_z_attr_avg": z_attr,
        "truth_x_dE_avg": x_de,
        "truth_y_dE_avg": y_de,
        "truth_z_dE_avg": z_de,
        # Preserve the column names used by the old plotting code, but switch the
        # values to the geometry-based drift length.
        "drift_length_attr_avg": drift_geom_attr,
        "drift_length_dE_avg": drift_geom_de,
        "drift_geometry_attr_avg": drift_geom_attr,
        "drift_geometry_dE_avg": drift_geom_de,
        "drift_325_minus_x_attr_avg": drift_325_attr,
        "drift_325_minus_x_dE_avg": drift_325_de,
        "drift_geometry_signed_attr_avg": drift_signed_attr,
        "drift_geometry_signed_dE_avg": drift_signed_de,
        "segment_ids": ";".join(str(int(sid)) for sid in seg_w.keys()),
        "segment_weights": ";".join(f"{float(w):.8g}" for w in seg_w.values()),
    }


# -----------------------------------------------------------------------------
# Matching and row building
# -----------------------------------------------------------------------------
def match_tiles_to_pixel_cluster(
    px: np.ndarray,
    py: np.ndarray,
    p_start: np.ndarray,
    tile_x: np.ndarray,
    tile_y: np.ndarray,
    tile_start: np.ndarray,
    cfg: Config,
) -> np.ndarray:
    """
    Match tile hits whose 20x20-pixel active area overlaps the pixel cluster.

    Space: tile is treated as a square of side 20*pixel_pitch centered on tile x/y.
    Time: pixel and tile hits are treated as active for waveform_ticks after activation;
          a tile is matched if its active interval overlaps the cluster active interval.
    """
    if len(tile_x) == 0 or len(px) == 0:
        return np.zeros(len(tile_x), dtype=bool)

    x_min = float(np.nanmin(px) - cfg.pixel_half_size)
    x_max = float(np.nanmax(px) + cfg.pixel_half_size)
    y_min = float(np.nanmin(py) - cfg.pixel_half_size)
    y_max = float(np.nanmax(py) + cfg.pixel_half_size)
    t_min = float(np.nanmin(p_start) - cfg.pre_samples * cfg.tick_size)
    t_max = float(np.nanmax(p_start) + cfg.readout_window)

    tile_x_min = tile_x - cfg.tile_half_size
    tile_x_max = tile_x + cfg.tile_half_size
    tile_y_min = tile_y - cfg.tile_half_size
    tile_y_max = tile_y + cfg.tile_half_size
    tile_t_min = tile_start - cfg.pre_samples * cfg.tick_size
    tile_t_max = tile_start + cfg.readout_window

    spatial_overlap = (tile_x_max >= x_min) & (tile_x_min <= x_max) & (tile_y_max >= y_min) & (tile_y_min <= y_max)
    time_overlap = (tile_t_max >= t_min) & (tile_t_min <= t_max)
    return spatial_overlap & time_overlap


def build_cluster_row(
    cluster_id: int,
    event_id: int,
    cluster_local_indices: np.ndarray,
    p_global_idx: np.ndarray,
    px: np.ndarray,
    py: np.ndarray,
    pt: np.ndarray,
    pwf: np.ndarray,
    plabel: np.ndarray,
    pattr: np.ndarray,
    tile_global_idx: np.ndarray,
    tx: np.ndarray,
    ty: np.ndarray,
    tt: np.ndarray,
    twf: np.ndarray,
    truth: Dict[Tuple[int, int], Dict[str, float]],
    cfg: Config,
) -> Optional[Dict[str, object]]:
    psel = cluster_local_indices
    pcharge_each = waveform_sum(pwf[psel])
    pixel_charge = float(np.sum(pcharge_each))

    if len(psel) < cfg.min_cluster_pixel_hits or pixel_charge <= cfg.min_cluster_charge:
        return None

    tmask = match_tiles_to_pixel_cluster(px[psel], py[psel], pt[psel], tx, ty, tt, cfg)
    tsel = np.nonzero(tmask)[0]
    tcharge_each = waveform_sum(twf[tsel]) if len(tsel) else np.array([], dtype=float)
    tile_charge = float(np.sum(tcharge_each)) if len(tsel) else 0.0

    pixel_charge_3sd = float(np.sum(waveform_sum_above(pwf[psel], cfg.threshold, subtract_threshold=False)))
    pixel_charge_above_thr_sub = float(np.sum(waveform_sum_above(pwf[psel], cfg.threshold, subtract_threshold=True)))
    tile_charge_3sd = float(np.sum(waveform_sum_above(twf[tsel], cfg.threshold, subtract_threshold=False))) if len(tsel) else 0.0
    tile_charge_above_thr_sub = float(np.sum(waveform_sum_above(twf[tsel], cfg.threshold, subtract_threshold=True))) if len(tsel) else 0.0

    pixel_wf_sum = sum_waveforms(pwf[psel])
    tile_wf_sum = sum_waveforms(twf[tsel]) if len(tsel) else np.array([], dtype=float)

    pixel_fit = fit_summed_waveform(pixel_wf_sum, cfg) if cfg.fit_waveforms else {}
    tile_fit = fit_summed_waveform(tile_wf_sum, cfg) if cfg.fit_waveforms and len(tsel) else {}
    pixel_individual_fit = (
        fit_pixel_waveforms_individually(pwf[psel], cfg) if cfg.fit_waveforms else {}
    )
    if cfg.fit_waveforms and cfg.plot_fit_examples:
        save_pixel_fit_example_plots(
            cluster_id=cluster_id,
            event_id=event_id,
            pixel_waveforms=pwf[psel],
            individual_fit_results=pixel_individual_fit.get("individual_fit_results", []),
            cfg=cfg,
        )

    seg_w = segment_weights_from_hits(plabel, pattr, pwf, psel)
    ts = truth_summary(event_id, seg_w, truth, cfg)

    row: Dict[str, object] = {
        "cluster_id": int(cluster_id),
        "event_id": int(event_id),

        "n_pixel_hits": int(len(psel)),
        "n_tile_hits": int(len(tsel)),
        "pixel_global_rows": ";".join(map(str, p_global_idx[psel].tolist())),
        "matched_tile_global_rows": ";".join(map(str, tile_global_idx[tsel].tolist())) if len(tsel) else "",

        "pixel_charge": pixel_charge,
        "tile_charge": tile_charge,
        "R_pixel_over_tile": pixel_charge / tile_charge if tile_charge > 0 else np.nan,

        "pixel_charge_3sd": pixel_charge_3sd,
        "tile_charge_3sd": tile_charge_3sd,
        "R_3sd_pixel_over_tile_rawtile": pixel_charge_3sd / tile_charge if tile_charge > 0 else np.nan,
        "R_3sd_pixel_over_tile_3sdtile": pixel_charge_3sd / tile_charge_3sd if tile_charge_3sd > 0 else np.nan,

        "pixel_charge_above_threshold_subtracted": pixel_charge_above_thr_sub,
        "tile_charge_above_threshold_subtracted": tile_charge_above_thr_sub,
        "R_above_threshold_subtracted": (
            pixel_charge_above_thr_sub / tile_charge_above_thr_sub if tile_charge_above_thr_sub > 0 else np.nan
        ),
        "threshold": float(cfg.threshold),

        "pixel_centroid_x": weighted_mean(px[psel], pcharge_each),
        "pixel_centroid_y": weighted_mean(py[psel], pcharge_each),
        "pixel_centroid_t": weighted_mean(pt[psel], pcharge_each),
        "pixel_width_x": weighted_rms(px[psel], pcharge_each),
        "pixel_width_y": weighted_rms(py[psel], pcharge_each),
        "pixel_width_t": weighted_rms(pt[psel], pcharge_each),
        "pixel_span_x": span(px[psel]),
        "pixel_span_y": span(py[psel]),
        "pixel_span_t": span(pt[psel]),
        "pixel_time_min": float(np.nanmin(pt[psel])),
        "pixel_time_max": float(np.nanmax(pt[psel]) + cfg.readout_window),

        "tile_centroid_x": weighted_mean(tx[tsel], tcharge_each) if len(tsel) else np.nan,
        "tile_centroid_y": weighted_mean(ty[tsel], tcharge_each) if len(tsel) else np.nan,
        "tile_centroid_t": weighted_mean(tt[tsel], tcharge_each) if len(tsel) else np.nan,
        "tile_width_x": weighted_rms(tx[tsel], tcharge_each) if len(tsel) else np.nan,
        "tile_width_y": weighted_rms(ty[tsel], tcharge_each) if len(tsel) else np.nan,
        "tile_width_t": weighted_rms(tt[tsel], tcharge_each) if len(tsel) else np.nan,
        "tile_span_x": span(tx[tsel]) if len(tsel) else np.nan,
        "tile_span_y": span(ty[tsel]) if len(tsel) else np.nan,
        "tile_span_t": span(tt[tsel]) if len(tsel) else np.nan,
        "tile_time_min": float(np.nanmin(tt[tsel])) if len(tsel) else np.nan,
        "tile_time_max": float(np.nanmax(tt[tsel]) + cfg.readout_window) if len(tsel) else np.nan,
    }

    row.update(ts)

    if cfg.fit_waveforms:
        for prefix, d in [("pixel", pixel_fit), ("tile", tile_fit)]:
            for key, val in d.items():
                row[f"{prefix}_{key}"] = val

        row["pixel_individual_fit_charge"] = pixel_individual_fit.get("individual_fit_charge", np.nan)
        row["pixel_individual_fit_charge_3sd"] = pixel_individual_fit.get("individual_fit_charge_3sd", np.nan)
        row["pixel_individual_fit_charge_above_threshold_subtracted"] = pixel_individual_fit.get(
            "individual_fit_charge_above_threshold_subtracted", np.nan
        )
        row["pixel_individual_fit_success_count"] = pixel_individual_fit.get("individual_fit_success_count", 0)
        row["pixel_individual_fit_fail_count"] = pixel_individual_fit.get("individual_fit_fail_count", 0)

        # Ratios using the Gaussian-fitted pixel charge above 3sigma.
        # The raw-tile version is the closest analog of R_3sd_pixel_over_tile_rawtile,
        # but with pixel charge taken from fitted Gaussian waveforms.
        pf = row.get("pixel_fit_charge", np.nan)
        tf = row.get("tile_fit_charge", np.nan)
        row["R_fit_pixel_over_tile"] = pf / tf if np.isfinite(pf) and np.isfinite(tf) and tf > 0 else np.nan

        pif3 = row.get("pixel_individual_fit_charge_3sd", np.nan)
        row["R_pixel_individual_fit3sd_over_rawtile"] = (
            pif3 / tile_charge if np.isfinite(pif3) and tile_charge > 0 else np.nan
        )
        row["R_pixel_individual_fit3sd_over_3sdtile"] = (
            pif3 / tile_charge_3sd if np.isfinite(pif3) and tile_charge_3sd > 0 else np.nan
        )

        # Also provide the summed-cluster-waveform fit above threshold for comparison.
        psf3 = row.get("pixel_fit_charge_3sd", np.nan)
        row["R_pixel_summed_fit3sd_over_rawtile"] = (
            psf3 / tile_charge if np.isfinite(psf3) and tile_charge > 0 else np.nan
        )
        row["R_pixel_summed_fit3sd_over_3sdtile"] = (
            psf3 / tile_charge_3sd if np.isfinite(psf3) and tile_charge_3sd > 0 else np.nan
        )

    return row



def collect_selected_events_and_segment_ids(cfg: Config) -> Tuple[List[int], set[int]]:
    """
    Scan only the detector events that will be processed and collect their pixel
    segment labels.  This lets the truth loader store only relevant segments.
    """
    selected_event_ids: List[int] = []
    segment_ids: set[int] = set()

    with h5py.File(cfg.hit_file, "r") as f:
        if "pixels" not in f:
            raise KeyError("Hit file must contain top-level 'pixels' dataset")
        P = f["pixels"]

        p_event = read_field(P, ["event id", "event_id"]).astype(int)
        p_order = np.argsort(p_event, kind="stable")
        p_event_sorted = p_event[p_order]

        for iev, (event_id, ps0, ps1) in enumerate(event_slices(p_event_sorted), start=1):
            if cfg.max_events is not None and iev > cfg.max_events:
                break

            selected_event_ids.append(int(event_id))
            pidx = p_order[ps0:ps1]

            # Apply the same basic active-hit selection as the main processing,
            # so quick tests do not collect labels from irrelevant empty hits.
            pwf = read_field(P, ["waveform"], pidx).astype(float)
            pcharge = waveform_sum(pwf)
            active = waveform_active_mask(pwf) & (pcharge > cfg.min_pixel_charge)
            local_keep = np.nonzero(active)[0]

            if cfg.max_pixel_hits_per_event is not None and len(local_keep) > cfg.max_pixel_hits_per_event:
                qkeep = pcharge[local_keep]
                chosen = np.argpartition(qkeep, -cfg.max_pixel_hits_per_event)[-cfg.max_pixel_hits_per_event:]
                local_keep = local_keep[chosen]

            if len(local_keep) == 0:
                continue

            plabel = read_field(P, ["label", "labels"], pidx[local_keep]).astype(int)
            pattr = read_field(P, ["attribution", "attributions"], pidx[local_keep]).astype(float)
            wf_keep = pwf[local_keep]

            # Collect only labels with positive attribution.  This matches the
            # truth-weighting logic used later and avoids loading truth segments
            # for empty/sentinel label slots such as attribution 0, -1, or -999.
            n_ticks = wf_keep.shape[1] if wf_keep.ndim == 2 else 1
            lab, attr = normalize_label_attr(plabel, pattr, n_ticks)
            if wf_keep.ndim == 1:
                wf_keep = wf_keep[:, None]
            ticks = min(wf_keep.shape[1], attr.shape[1])
            wf_keep = wf_keep[:, :ticks]
            attr = attr[:, :ticks, :]

            for k in range(lab.shape[1]):
                positive = (
                    (lab[:, k] != INVALID_SEGMENT_ID)
                    & np.any(
                        np.isfinite(attr[:, :, k])
                        & np.isfinite(wf_keep)
                        & (attr[:, :, k] > 0.0)
                        & (wf_keep > 0.0),
                        axis=1,
                    )
                )
                if np.any(positive):
                    segment_ids.update(int(s) for s in np.unique(lab[positive, k]))

    return selected_event_ids, segment_ids


# -----------------------------------------------------------------------------
# Outputs
# -----------------------------------------------------------------------------
def write_outputs(df: pd.DataFrame, cfg: Config) -> None:
    os.makedirs(cfg.output_dir, exist_ok=True)
    all_path = os.path.join(cfg.output_dir, f"{cfg.tag}_clusters_all.csv")
    df.to_csv(all_path, index=False)
    print(f"Wrote all clusters: {all_path}  ({len(df)} rows)")

    # Keep cuts separate.  The all-cluster table is the authoritative output.
    cut_specs = {
        "cut_R_0_1": df["R_pixel_over_tile"].between(0, 1, inclusive="right") if "R_pixel_over_tile" in df else pd.Series(False, index=df.index),
        "cut_R3sd_rawtile_0_1": df["R_3sd_pixel_over_tile_rawtile"].between(0, 1, inclusive="right") if "R_3sd_pixel_over_tile_rawtile" in df else pd.Series(False, index=df.index),
        "matched_tiles_only": (df["n_tile_hits"] > 0) if "n_tile_hits" in df else pd.Series(False, index=df.index),
        "cut_Rfit3sd_rawtile_0_1": (
            df["R_pixel_individual_fit3sd_over_rawtile"].between(0, 1, inclusive="right")
            if "R_pixel_individual_fit3sd_over_rawtile" in df else pd.Series(False, index=df.index)
        ),
    }
    for name, mask in cut_specs.items():
        path = os.path.join(cfg.output_dir, f"{cfg.tag}_{name}.csv")
        df.loc[mask].to_csv(path, index=False)
        print(f"Wrote {name}: {path}  ({int(mask.sum())} rows)")

    summary = {
        "n_clusters_all": int(len(df)),
        "n_clusters_with_matched_tiles": int((df["n_tile_hits"] > 0).sum()) if "n_tile_hits" in df else 0,
        "n_clusters_with_truth": int((df["n_truth_segments"] > 0).sum()) if "n_truth_segments" in df else 0,
        "config": cfg.__dict__,
    }
    with open(os.path.join(cfg.output_dir, f"{cfg.tag}_summary.json"), "w") as fp:
        json.dump(summary, fp, indent=2)

    if cfg.make_plots:
        make_plots(df, cfg)


def make_plots(df: pd.DataFrame, cfg: Config) -> None:
    if plt is None:
        print("matplotlib not available; skipping plots")
        return
    plot_dir = os.path.join(cfg.output_dir, "plots")
    os.makedirs(plot_dir, exist_ok=True)

    def save_hist(col: str, bins=100, rng=None):
        if col not in df:
            return
        vals = pd.to_numeric(df[col], errors="coerce").replace([np.inf, -np.inf], np.nan).dropna()
        if len(vals) == 0:
            return
        plt.figure()
        plt.hist(vals, bins=bins, range=rng)
        plt.xlabel(col)
        plt.ylabel("clusters")
        plt.tight_layout()
        plt.savefig(os.path.join(plot_dir, f"hist_{col}.png"), dpi=160)
        plt.close()

    save_hist("n_tile_hits", bins=range(0, int(pd.to_numeric(df["n_tile_hits"], errors="coerce").max()) + 2) if "n_tile_hits" in df and len(df) else 50)
    save_hist("R_pixel_over_tile", rng=(0, 2))
    save_hist("R_3sd_pixel_over_tile_rawtile", rng=(0, 2))
    save_hist("R_pixel_individual_fit3sd_over_rawtile", rng=(0, 2))
    save_hist("R_pixel_individual_fit3sd_over_3sdtile", rng=(0, 2))
    save_hist("pixel_individual_fit_success_count")
    save_hist("pixel_individual_fit_fail_count")
    save_hist("drift_length_attr_avg")
    save_hist("drift_325_minus_x_attr_avg")
    save_hist("drift_geometry_signed_attr_avg")

    if "drift_length_attr_avg" in df and "R_pixel_over_tile" in df:
        x = pd.to_numeric(df["drift_length_attr_avg"], errors="coerce")
        y = pd.to_numeric(df["R_pixel_over_tile"], errors="coerce")
        good = np.isfinite(x) & np.isfinite(y)
        if good.sum() > 0:
            plt.figure()
            plt.scatter(x[good], y[good], s=4, alpha=0.4)
            plt.xlabel("drift_length_attr_avg")
            plt.ylabel("R_pixel_over_tile")
            plt.tight_layout()
            plt.savefig(os.path.join(plot_dir, "R_vs_drift_length_attr_avg.png"), dpi=160)
            plt.close()

    if "drift_length_attr_avg" in df and "R_pixel_individual_fit3sd_over_rawtile" in df:
        x = pd.to_numeric(df["drift_length_attr_avg"], errors="coerce")
        y = pd.to_numeric(df["R_pixel_individual_fit3sd_over_rawtile"], errors="coerce")
        good = np.isfinite(x) & np.isfinite(y)
        if good.sum() > 0:
            plt.figure()
            plt.scatter(x[good], y[good], s=4, alpha=0.4)
            plt.xlabel("drift_length_attr_avg")
            plt.ylabel("R_pixel_individual_fit3sd_over_rawtile")
            plt.tight_layout()
            plt.savefig(os.path.join(plot_dir, "R_fit3sd_rawtile_vs_drift_length_attr_avg.png"), dpi=160)
            plt.close()


# -----------------------------------------------------------------------------
# Main loop
# -----------------------------------------------------------------------------
def main() -> None:
    cfg = parse_args()
    os.makedirs(cfg.output_dir, exist_ok=True)

    print("Pixel-first clustering analysis")
    print(f"Hit file:   {cfg.hit_file}")
    print(f"Truth file: {cfg.truth_file}")
    print(f"Output dir: {cfg.output_dir}")
    print(f"Tile size:  {cfg.tile_size:g} = {cfg.tile_size_pixels:g} * pixel_pitch")
    print(f"Readout window: {cfg.readout_window:g} time units = {cfg.waveform_ticks} ticks")
    print(f"Detector YAML: {cfg.detector_yaml}")
    print(f"Drift mode: {cfg.drift_mode}")
    print(f"Truth lookup: {cfg.truth_lookup}")
    print("Truth weights: require attribution > 0 and waveform charge > 0 per tick")
    if cfg.fit_waveforms:
        print("Gaussian fits: enabled")
        if cfg.plot_fit_examples:
            print(
                f"Gaussian fit examples: first {cfg.fit_example_clusters} clusters, "
                f"{cfg.fit_example_hits_per_cluster} pixel hits each"
            )

    selected_event_ids, wanted_segment_ids = collect_selected_events_and_segment_ids(cfg)
    print(f"Detector events selected for processing: {len(selected_event_ids)}")
    print(f"Unique detector segment labels needed: {len(wanted_segment_ids)}")

    truth = load_truth_segments(
        cfg.truth_file,
        cfg,
        wanted_event_ids=set(selected_event_ids),
        wanted_segment_ids=wanted_segment_ids,
    )
    print(f"Loaded truth segment map entries: {len(truth)}")

    rows: List[Dict[str, object]] = []
    cluster_id = 0
    noise_hits_total = 0
    events_seen = 0

    with h5py.File(cfg.hit_file, "r") as f:
        if "pixels" not in f or "tiles" not in f:
            raise KeyError("Hit file must contain top-level 'pixels' and 'tiles' datasets")
        P = f["pixels"]
        T = f["tiles"]

        p_event = read_field(P, ["event id", "event_id"]).astype(int)
        t_event = read_field(T, ["event id", "event_id"]).astype(int)
        p_order = np.argsort(p_event, kind="stable")
        t_order = np.argsort(t_event, kind="stable")
        p_event_sorted = p_event[p_order]
        t_event_sorted = t_event[t_order]
        tile_slices = {eid: (s, e) for eid, s, e in event_slices(t_event_sorted)}

        for iev, (event_id, ps0, ps1) in enumerate(event_slices(p_event_sorted), start=1):
            if cfg.max_events is not None and iev > cfg.max_events:
                break
            events_seen += 1

            pidx = p_order[ps0:ps1]
            px = read_field(P, ["pixel x", "pixel_x", "x"], pidx).astype(float)
            py = read_field(P, ["pixel y", "pixel_y", "y"], pidx).astype(float)
            pt = read_field(P, ["start t", "trig t", "start_t", "trig_t", "hit t", "hit_t"], pidx).astype(float)
            pwf = read_field(P, ["waveform"], pidx).astype(float)
            plabel = read_field(P, ["label", "labels"], pidx).astype(int)
            pattr = read_field(P, ["attribution", "attributions"], pidx).astype(float)

            pcharge = waveform_sum(pwf)
            active = waveform_active_mask(pwf) & np.isfinite(px) & np.isfinite(py) & np.isfinite(pt) & (pcharge > cfg.min_pixel_charge)
            if not np.any(active):
                continue

            local_keep = np.nonzero(active)[0]
            if cfg.max_pixel_hits_per_event is not None and len(local_keep) > cfg.max_pixel_hits_per_event:
                # Keep the largest-charge pixel hits for quick tests only.
                qkeep = pcharge[local_keep]
                chosen = np.argpartition(qkeep, -cfg.max_pixel_hits_per_event)[-cfg.max_pixel_hits_per_event:]
                local_keep = local_keep[chosen]

            px_k = px[local_keep]
            py_k = py[local_keep]
            pt_k = pt[local_keep]
            pwf_k = pwf[local_keep]
            plabel_k = plabel[local_keep]
            pattr_k = pattr[local_keep]
            pidx_k = pidx[local_keep]
            pcharge_k = waveform_sum(pwf_k)

            if event_id in tile_slices:
                ts0, ts1 = tile_slices[event_id]
                tidx = t_order[ts0:ts1]
                tx = read_field(T, ["tile x", "tile_x", "x"], tidx).astype(float)
                ty = read_field(T, ["tile y", "tile_y", "y"], tidx).astype(float)
                tt = read_field(T, ["start t", "trig t", "start_t", "trig_t", "hit t", "hit_t"], tidx).astype(float)
                twf = read_field(T, ["waveform"], tidx).astype(float)
            else:
                tidx = np.array([], dtype=int)
                tx = ty = tt = np.array([], dtype=float)
                twf = np.empty((0, cfg.waveform_ticks), dtype=float)

            # Pixel-hit clustering: channel/hit level, not waveform-sample level.
            # Coordinates are normalized so eps=1 corresponds roughly to the two
            # requested scales: cluster_eps_space pixel pitches and cluster_eps_time time units.
            features = np.column_stack([
                px_k / (cfg.cluster_eps_space * cfg.pixel_pitch),
                py_k / (cfg.cluster_eps_space * cfg.pixel_pitch),
                pt_k / cfg.cluster_eps_time,
            ])
            db = DBSCAN(eps=cfg.dbscan_eps, min_samples=cfg.dbscan_min_samples).fit(features)
            labels = db.labels_
            noise_hits_total += int(np.sum(labels < 0))

            for cid in sorted(set(labels)):
                if cid < 0:
                    continue
                csel = np.nonzero(labels == cid)[0]
                row = build_cluster_row(
                    cluster_id=cluster_id,
                    event_id=event_id,
                    cluster_local_indices=csel,
                    p_global_idx=pidx_k,
                    px=px_k,
                    py=py_k,
                    pt=pt_k,
                    pwf=pwf_k,
                    plabel=plabel_k,
                    pattr=pattr_k,
                    tile_global_idx=tidx,
                    tx=tx,
                    ty=ty,
                    tt=tt,
                    twf=twf,
                    truth=truth,
                    cfg=cfg,
                )
                if row is None:
                    continue
                rows.append(row)
                cluster_id += 1

            if iev % 100 == 0:
                print(f"Processed {iev} events; clusters so far = {len(rows)}")

    df = pd.DataFrame(rows)
    print("Done processing")
    print(f"Events processed: {events_seen}")
    print(f"Clusters found:    {len(df)}")
    print(f"DBSCAN noise pixel hits: {noise_hits_total}")
    write_outputs(df, cfg)


if __name__ == "__main__":
    main()
