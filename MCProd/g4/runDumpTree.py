#!/usr/bin/env python
import argparse
import dumpTreeNew

if __name__ == "__main__":

    parser = argparse.ArgumentParser( description = 'Process GEANT4 simulation using edep-sim.' )

    parser.add_argument('--dir', type = str, required = True, help = 'directory of input and output files')

    args = parser.parse_args()

    sample = 'radiologicals'
    prefix = 'fullgeoanatruth-vd-reduced'
    nFiles = 28
    nEventsPerFile = 1000
    
    for iFile in range( 0, nFiles ):
        rootFile = f'{args.dir}/g4/{sample}/{prefix}_g4_{iFile:02d}.root'
        h5File = f'{args.dir}/g4/{sample}/{prefix}_g4_{iFile:02d}.h5'

        print(f'Dumping {rootFile}...')
        dumpTreeNew.dump(rootFile, h5File, True)
