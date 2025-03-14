#!/usr/bin/env python

import os

if __name__ == "__main__":

    dir = '/Users/yuntse/data/lartpc_rd/gampix/sn/garching/nh/marley'
    name = 'nueArCC_garching_nh_mxDir'
    nFiles = 56

    for iFile in range( nFiles ):
        infile = f'{dir}/{name}_{iFile:02d}.ascii'
        outfile = f'{dir}/{name}_{iFile:02d}.root'

        cmd = f'root -q -b \'MARLEYToROOT.C("{infile}", "{outfile}")\''
        print( cmd )
        os.system( cmd )
