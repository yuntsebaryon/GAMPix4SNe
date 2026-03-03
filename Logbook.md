Logbook
=======

## Energy Deposition issue: 2025.12.9

In August, John found that dE in the dumpTree stores the total energy (GetEnergyDeposit() in GEANT4), while we should use only the ionization energy in the following simulation.
He made `dumpTreeNew.py` correcting `dE` and adding `dE2`
> segment[iHit]["dE"] = (hitSegment.GetEnergyDeposit() - hitSegment.GetSecondaryDeposit()) \
> segment[iHit]["dE2"] = hitSegment.GetSecondaryDeposit()

2025.12.10

In `edep-sim`, there are two options on `/edep/phys/ionizationModel`. \
`/edep/phys/ionizationModel 0`: supposedly use NEST, but in SegmentHits (`src/EDepSimHitSegment.cc`), non-ionizing energy is always 0. \
`/edep/phys/ionizationModel 0`: supposedly use Birk model, but in SegmentHits, non-ionizing energy is always identical to the total energy.  Also, it doesn't seem to go through `src/EDepSimDokeBirks.cc`.

----------------------------------------------------------------------------------

## GEANT4 notes

### Neutral Paricles: Energy Deposition: 2025.1.7

Check how the neutral particles, neutrons and photons, deposit energy.  I only checked a handful of neutron events from my BRN samples used in the background study of the nue-Ar CC measurement in COHERENT:
- neutrons deposit energy in LAr when the end process is hadElastic, which is the hadron elastic scattering.
- photons deposit energy in LAr when the end process is phot, which I think is photo-electric effect
- not 100% of neutrons with the end process of hadElastic deposit energy in LAr, a small fraction of them do not.
- also a very small fraction of photons doesn't deposit energy with the end process phot
- All the energy deposition in each step in GEANT4 above refers to sub-MeV level or much smaller

__Update on 2025.1.15__

- Two variables: `G4Step::fTotalEnergyDeposit ` and `G4Step::fNonIonizingEnergyDeposit`.  The `dE` variable in my G4 code is `G4Step::fTotalEnergyDeposit`
> // Member data
> G4double fTotalEnergyDeposit;
>   // Accummulated total energy desposit in the current Step
>
> G4double fNonIonizingEnergyDeposit;
>  // Accummulated non-ionizing energy desposit in the current Step
- In the case of neutrons ending with `hadElastic`, the two variables `fTotalEnergyDeposit` and `fNonIonizingEnergyDeposit` have identical values, which means `dE` is non-ionizing energy deposit
- In the G4 hadron elastic process, https://apc.u-paris.fr/~franco/g4doxy/html/classG4HadronElasticProcess.html
>  edep = GetLocalEnergyDeposit + (if GetEnergyChange <= lowestEnergy) GetEnergyChange + (if p->GetKineticEnergy() <= tcut),
  where p is a secondary particle/recoil, and lowestEnergy default is 1 keV (or 1.e-6*eV for neutrons).  tcut depends on G4ProductionCutsTable, https://apc.u-paris.fr/~franco/g4doxy/html/G4ProductionCutsTable_8hh-source.html.
- In short, the non-zero dE from neutrons ending with hadElastic process represent the energy change below a threshold and the kinetic energy of the secondary particles below another threshold, as Hiro suggested.  It is not the ionizing energy.
- The dE from photons ending with photoelectric are ionizing energy.  They are typically very small (a few keV level or below) and can be safely ignored.

----------------------------------------------------------------------------------