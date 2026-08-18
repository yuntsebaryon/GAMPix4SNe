First drift distance analysis: Work Log
=======================================

## First Iteration

### Detector simulation outputs

- GAMPixPy sample: `/Users/yuntse/data/lartpc_rd/gampix/detsim/radiologicals/gampixpy_fullgeoanatruth-vd-reduced_g4_00_2Mhz_segmentlabel_lowtrig_5mmpitch.h5`
    - 2026/4/7
    - Triggered when the sum over the true signals within the tile readout window above 3&times;50e<sup>-</sup>
- First look at the GAMPixPy output: `Validation/detsim/checkReadout_20260413.ipynb`
- Test tile-pixel matching algorithms: `Validation/detsim/testMatchingAlgorithm.ipynb`
- Test the Gaussian fit algorithms: `Validation/detsim/testGaussianFitAlgos.ipynb`
- Fit the waveforms with Gaussian: `Validation/detsim/fitGaussian.ipynb`
- Systematically study the output variables and their correlation: `Validation/detsim/studyReadout1.ipynb`

### Neural network analysis

- `v0`: Bahrudin's parameters of the NN hidden layer sizes: `Analysis/TileCluster/driftDistanceNN0.ipynb`
- `v0.1`: The default NN parameters: `Analysis/TileCluster/driftDistanceNN0.1.ipynb`
    - Slightly better results, but not significant by eyes
    - Without any tile readout to justify that we need/want the tiles: `Analysis/TileCluster/driftDistanceNN0.1_noTile.ipynb`
- `v1`: Add one more variable in the NN inputs: `Analysis/TileCluster/driftDistanceNN1.ipynb`
- `v1.1`: Add a different variable in the NN inputs: `Analysis/TileCluster/driftDistanceNN1.1.ipynb`

--------------------------------------------------------------------------------------------------------------

## Second Iteration

### Detector simulation outputs

- GAMPixPy sample: `/Users/yuntse/data/lartpc_rd/gampix/detsim/radiologicals/gampixpy_fullgeoanatruth-vd-reduced_g4_00_2Mhz_segmentlabel_lowtrig_5mmpitch_presample.h5`
    - 2026/5/14
    - Triggered when
        - the sum over the true signals within the tile readout window above 3&times;50e<sup>-</sup> AND
        - the true signal at the trigger time above 3&times;50e<sup>-</sup>
- Fit the waveforms with Gaussian: `Validation/detsim/fitGaussian1.ipynb`
- Systematically study the output variables and their correlation: `Validation/detsim/studyReadout1.1.ipynb`

### Neural network analysis

- `v0.1.1`: The same as `v0.1` but using the new sample: `Analysis/TileCluster/driftDistanceNN0.1.1.ipynb`
    - Without any tile readout to justify that we need/want the tiles: `Analysis/TileCluster/driftDistanceNN0.1.1_noTile.ipynb`

### DUNE Collaboration Meeting Presentation: 2026.5.20

- Use the `v0.1.1` results