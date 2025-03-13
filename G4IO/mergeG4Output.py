import ROOT
import pandas as pd
import os

def writeCSV( infileName, outfileName, treeName, runNo, subrunNo, isSignal):

    print( f'Writing {infileName} to a CSV file...')
    infile = ROOT.TFile( infileName )
    tree = infile.Get( treeName )

    with open( outfileName, 'w') as f:
        if isSignal:
            print('run/I,subrun/I,event/I,isSignal/I,pdgCode/I,trackID/I,motherID/I,startE/D,dE/D,startX/D,startY/D,startZ/D,startT/D,endX/D,endY/D,endZ/D,endT/D', file = f)
        for Part, Evt in enumerate( tree ):
            print( f'{runNo},{subrunNo},{Evt.event},{isSignal},{Evt.pdg},{Evt.trackID},{Evt.motherID},{Evt.startE},{Evt.dE},{Evt.startX},{Evt.startY},{Evt.startZ},{Evt.startT},{Evt.endX},{Evt.endY},{Evt.endZ},{Evt.endT}', file = f)

    return
# def writeCSV()

def sortEvents( signalF, radF, outF ):

    print( f'Sorting and creating {outF}...')
    cmd = f'sort -t , -k 3,3 -s {signalF} {radF} > {outF}'
    os.system( cmd )
    return
# def sortEvents

if __name__ == "__main__":

    sigDir = '/Users/yuntse/data/lartpc_rd/gampix/sn/garching/nh/g4'
    radDir = '/Users/yuntse/data/lartpc_rd/gampix/radiologicals/g4'
    outDir = '/Users/yuntse/data/lartpc_rd/gampix/mixed/garching/nh/g4'
    runNo = 20000047
    subrunRange = 1

    xLoc = 0
    treeName = 'edep'
    sigName = f'nueArCC_garching_nh_mxDir_x{xLoc}_g4'
    radName = 'fullgeoanatruth-vd-reduced_g4'
    outName = f'rad_nueArCC_garching_nh_mxDir_x{xLoc}_g4'

    for iSubrun in range( subrunRange ):
    
        # SN neutrino signal file
        sigInfile = f'{sigDir}/{sigName}_{iSubrun:04d}.root'
        sigOutfile = f'{sigDir}/{sigName}_{iSubrun:04d}.csv'
        writeCSV( sigInfile, sigOutfile, treeName, runNo, iSubrun, 1)

        # Radiological background file
        radInfile = f'{radDir}/{radName}_{iSubrun:04d}.root'
        radOutfile = f'{radDir}/{radName}_{iSubrun:04d}.csv'
        writeCSV( radInfile, radOutfile, treeName, runNo, iSubrun, 0)

        outfile = f'{outDir}/{outName}_{iSubrun:04d}.csv'
        sortEvents( sigOutfile, radOutfile, outfile )