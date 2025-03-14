#!/usr/bin/env python

import ROOT
import pandas as pd

def divideAndwriteCVS( tRad, prefix, subrun, ymin, ymax, zmin, zmax, nEvents):

    csvName = f'{prefix}_{subrun:04d}.csv'
    yCenter = ( ymin + ymax )/2.
    zCenter = ( zmin + zmax )/2. - (-107+2198.88)/2

    with open(csvName, 'w') as f0:
        print('run/I,subrun/I,event/I,pdgCode/I,x/D,y/D,z/D,t/D,px/D,py/D,pz/D,E/D,mass/D', file = f0)
        for iPart, iEvt in enumerate( tRad ):
            if iEvt.event > nEvents: continue
            if ( iEvt.y > ymax or iEvt.y < ymin or iEvt.z > zmax or iEvt.z < zmin ): continue
            print(f'{iEvt.run},{subrun},{iEvt.event-1},{iEvt.pdg_code},{iEvt.x},{iEvt.y-yCenter},{iEvt.z-zCenter},{iEvt.t},{iEvt.px},{iEvt.py},{iEvt.pz},{iEvt.Etot},{iEvt.mass}', file = f0)

    csvName = f'{prefix}_{subrun+1:04d}.csv'

    with open(csvName, 'w') as f1:
        print('run/I,subrun/I,event/I,pdgCode/I,x/D,y/D,z/D,t/D,px/D,py/D,pz/D,E/D,mass/D', file = f1)
        for iPart, iEvt in enumerate( tRad ):
            if iEvt.event <= nEvents: continue
            if ( iEvt.y > ymax or iEvt.y < ymin or iEvt.z > zmax or iEvt.z < zmin ): continue
            EvtNo = iEvt.event - 1 - nEvents
            print(f'{iEvt.run},{subrun},{EvtNo},{iEvt.pdg_code},{iEvt.x},{iEvt.y-yCenter},{iEvt.z-zCenter},{iEvt.t},{iEvt.px},{iEvt.py},{iEvt.pz},{iEvt.Etot},{iEvt.mass}', file = f1)

    return
# def divideAndwriteCVS()

def writeROOT( rootName, csvName):

    fDest = ROOT.TFile( rootName, "RECREATE")
    newTree = ROOT.TTree("SNEvent", "SNEvent")
    newTree.ReadFile( csvName, "", ",")
    newTree.SetDirectory(fDest)
    newTree.Write()
    fDest.Close()

    return
# def writeROOT()

if __name__ == "__main__":

    Dir = '/Users/yuntse/data/lartpc_rd/gampix/radiologicals'
    infile = f'{Dir}/dune/fullgeoanatruth-vd1x8x14-1000.root'
    prefix = f'{Dir}/gen/fullgeoanatruth-vd-reduced'

    # Division
    nY = 4
    nZ = 7
    nEvents = 500
    nSub = 2

    # cm
    YLow = -600.
    ZLow = 0.
    YLength = 300.
    ZLength = 300.

    fRad = ROOT.TFile( infile )
    tRad = fRad.Get("fullgeoanatruth/FullGeoAnaTruth")

    for iy in range( nY ):
        for iz in range( nZ ):

            ymin = YLength*float(iy) + YLow
            ymax = YLength*float(iy+1) + YLow
            zmin = ZLength*float(iz) + ZLow
            zmax = ZLength*float(iz+1) + ZLow
            subrun = (iy*nZ + iz)*2

            print( f'Processing Subrun {subrun}, {ymin} < y < {ymax}cm, {zmin} < z < {zmax}cm...')
            divideAndwriteCVS( tRad, prefix, subrun, ymin, ymax, zmin, zmax, nEvents)

            for isub in range( nSub ):

                csvName = f'{prefix}_{subrun+isub:04d}.csv'
                rootName = f'{prefix}_{subrun+isub:04d}.root'

                writeROOT( rootName, csvName )