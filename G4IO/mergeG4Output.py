import ROOT
import pandas as pd
import os

def writeCSV( infileName, outfileName, treeName, runNo, subrunNo, isSignal):

    print( f'Writing {infileName} to a CSV file...')
    infile = ROOT.TFile( infileName )
    tree = infile.Get( treeName )

    with open( outfileName, 'w') as f:
        # if isSignal:
        #     print('run,subrun,event,isSignal,pdgCode,trackID,motherID,startE,dE,startX,startY,startZ,startT,endX,endY,endZ,endT', file = f)
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
    outDir = '/Users/yuntse/data/lartpc_rd/gampix/mixed/garching/nh/temp'
    runNo = 20000047
    nSubruns = 2

    x = [ 315, 205, 105, 0, -105, -205, -305 ]
    treeName = 'edep'
    sigName = f'nueArCC_garching_nh_mxDir'
    radName = 'fullgeoanatruth-vd-reduced_g4'
    outName = f'rad_nueArCC_garching_nh_mxDir'

    for iSubrun in range( nSubruns ):
    
        # Radiological background file
        radInfile = f'{radDir}/{radName}_{iSubrun:04d}.root'
        radOutfile = f'{radDir}/{radName}_{iSubrun:04d}.csv'
        writeCSV( radInfile, radOutfile, treeName, runNo, iSubrun, 0)
        
        for ix in x:
            # SN neutrino signal file
            sigInfile = f'{sigDir}/{sigName}_xpos{ix:03d}cm_g4_{iSubrun:04d}.root'
            sigOutfile = f'{sigDir}/{sigName}_xpos{ix:03d}cm_g4_{iSubrun:04d}.csv'
            writeCSV( sigInfile, sigOutfile, treeName, runNo, iSubrun, 1)

            outfile = f'{outDir}/{outName}_xpos{ix:03d}cm_g4_{iSubrun:04d}.csv'
            sortEvents( sigOutfile, radOutfile, outfile )