#!/usr/bin/env python

import ROOT
import pandas as pd

def writeCSV( tSN, csvName, run, subrun, x, y, z, t):

    iEvt = 0
    with open(csvName, 'w') as fOut:
        print('run/I,subrun/I,event/I,pdgCode/I,x/D,y/D,z/D,t/D,px/D,py/D,pz/D,E/D,mass/D', file = fOut)
        for Part, Evt in enumerate( tSN ):
            for pdg, E, px, py, pz, m in zip(Evt.pdg, Evt.E, Evt.px, Evt.py, Evt.pz, Evt.m):
                print(f'{run},{subrun},{iEvt},{pdg},{x},{y},{z},{t},{px/1000.},{py/1000.},{pz/1000.},{E/1000.},{m/1000.}', file = fOut)
            iEvt += 1
    return
# def writeCSV()

def writeROOT( rootName, csvName):

    fDest = ROOT.TFile(rootName, "RECREATE")
    newTree = ROOT.TTree("SNEvent", "SNEvent")
    newTree.ReadFile(csvName, "", ",")
    newTree.SetDirectory(fDest)
    newTree.Write()
    fDest.Close()

    return
# def writeROOT()

if __name__ == "__main__":

    Dir = '/Users/yuntse/data/lartpc_rd/gampix/sn/garching/nh'
    prefix = 'nueArCC_garching_nh_mxDir'

    run = 20000047
    nSubrun = 56
    x = [ 315, 205, 105, 0, -105, -205, -305 ]
    y = 0.
    z = (-107+2198.88)/2.
    t = 0.

    for subrun in range( nSubrun ):

        print( f'Processing subrun {subrun}...')

        infile = f'{Dir}/marley/{prefix}_{subrun:02d}.root'
        fSN = ROOT.TFile( infile )
        tSN = fSN.Get("kin")

        for ix in x:

            csvName = f'{Dir}/gen/{prefix}_xpos{ix:03d}cm_{subrun:04d}.csv'
            rootName = f'{Dir}/gen/{prefix}_xpos{ix:03d}cm_{subrun:04d}.root'

            writeCSV( tSN, csvName, run, subrun, ix, y, z, t)
            writeROOT( rootName, csvName )