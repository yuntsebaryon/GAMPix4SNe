import argparse
import os

if __name__ == "__main__":

    parser = argparse.ArgumentParser( description = 'Process GEANT4 simulation using edep-sim.' )

    parser.add_argument('--dir', type = str, required = True, help = 'directory of input and output files')

    args = parser.parse_args()

    sample = 'nurad'
    prefix = 'nueArCC_garching_nh_mxpypzDir-rad-vd-reduced'
    nFiles = 28
    nEventsPerFile = 1000

    gdml = '/Users/yuntse/work/lartpc_rd/GAMPix4SNe/MCProd/g4/gdml/InfiniteLAr.gdml'
    
    # executable = '/Users/yuntse/opt/edep-sim/edep-gcc-17.0.0-arm64-apple-darwin24.5.0/bin/edep-sim'
    executable = '/Users/yuntse/opt/edep-sim-origin/edep-gcc-17.0.0-arm64-apple-darwin24.5.0/bin/edep-sim'

    inDir = f'{args.dir}/gen/{sample}'
    if not os.path.isdir( inDir ):
        raise FileNotFoundError(f"Input directory '{inDir}' does not exist")

    outDir = f'{args.dir}/g4/{sample}'
    # if os.path.exists( outDir ):
    #     raise FileExistsError(f"Output directory '{outDir}' already exists.")
    # else:
    #     os.makedirs( f'{outDir}/mac')

    for iFile in range( 2, nFiles ):
        macfile = f'{outDir}/mac/{prefix}_g4_{iFile:02d}.mac'
        infile = f'{inDir}/{prefix}_{iFile:02d}.hepevt'
        outfile = f'{outDir}/{prefix}_g4_{iFile:02d}.root'

        with open( macfile, 'w') as f:
            f.write(
f'''
/edep/gdml/read {gdml}

/edep/db/set/neutronThreshold 0 MeV
/edep/db/set/lengthThreshold 0 mm
/edep/db/set/gammaThreshold 0 MeV
/edep/db/open {outfile}

/edep/hitSagitta LArVol 0.1 mm
/edep/hitLength LArVol 0.1 mm

/edep/hitSeparation LArVol -1 mm

/edep/update

/generator/kinematics/set hepevt

/generator/kinematics/hepevt/input {infile}
/generator/kinematics/hepevt/flavor marley

## Set the counter of the vertices in an event
/generator/count/set hepevt
/generator/count/hepevt/input {infile}

## Distribute the events based on the vertex in the rooTracker file.
/generator/position/set free

## Distribute the event times based onn the time in the rooTracker file.
/generator/time/set free

## Make sure EDEPSIM updates the kinematics generator.
/generator/add

/run/beamOn {nEventsPerFile}
'''
            )
        
        cmd = f'{executable} {macfile}'
        print( cmd )
        os.system( cmd )
