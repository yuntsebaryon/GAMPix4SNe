#!/usr/bin/env python

import os

if __name__ == "__main__":

    dir = '/Users/yuntse/data/lartpc_rd/gampix/mixed/garching/nh'
    prefix = 'rad_nueArCC_garching_nh_mxDir'

    nSubruns = 2
    x = [ 315, 205, 105, 0, -105, -205, -305 ]

    for iSubrun in range( nSubruns ):
        for ix in x:

            infile = f'{dir}/temp/{prefix}_xpos{ix:03d}cm_g4_{iSubrun:04d}.csv'
            outfile = f'{dir}/g4/{prefix}_xpos{ix:03d}cm_g4_{iSubrun:04d}.csv'

            with open(outfile, 'w') as f:
                print('run,subrun,event,isSignal,pdgCode,trackID,motherID,startE,dE,startX,startY,startZ,startT,endX,endY,endZ,endT', file = f)

            cmd = f'cat {infile} >> {outfile}'
            os.system( cmd )