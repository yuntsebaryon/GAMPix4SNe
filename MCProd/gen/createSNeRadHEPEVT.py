import ROOT
import numpy as np

def readHepevt(filename, EvtNo):
    events = []
    with open(filename, 'r') as f:
        while True:
            header = f.readline()
            if not header:
                break  # End of file
            parts = header.strip().split()
            if len(parts) < 2:
                continue
            n_particles = int(parts[1])
            event_id = EvtNo

            particles = []
            for _ in range(n_particles):
                line = f.readline()
                data = line.strip().split()
                if len(data) < 15:
                    continue
                particle = {
                    "ISTHEP": int(data[0]),
                    "IDHEP": int(data[1]),
                    "JMOHEP1": int(data[2]),
                    "JMOHEP2": int(data[3]),
                    "JDAHEP1": int(data[4]),
                    "JDAHEP2": int(data[5]),
                    "PHEP1": float(data[6]),
                    "PHEP2": float(data[7]),
                    "PHEP3": float(data[8]),
                    "PHEP4": float(data[9]),
                    "PHEP5": float(data[10]),
                    "VHEP1": float(data[11]),
                    "VHEP2": float(data[12]),
                    "VHEP3": float(data[13]),
                    "VHEP4": float(data[14]),
                }
                particles.append(particle)

            events.append((event_id, particles))
            EvtNo += 1
    return events
# def readHepevt(filename, EvtNo)

def assignXYZT(events, x, y ,z, t, firstEvtNo):
    lengths = [ len(events), len(x), len(y), len(z), len(t) ]
    if len(set(lengths)) != 1:
        raise ValueError( f'Array lengths are not identical!')
        
    for event_id, particles in events:
        ievt = event_id - firstEvtNo
        for p in particles:
            p['VHEP1'] = x[ievt]
            p['VHEP2'] = y[ievt]
            p['VHEP3'] = z[ievt]
            p['VHEP4'] = t[ievt]
    return events
# def assignXYZT(events, x, y ,z, t, firstEvtNo)

def xyzInBox(xmin, xmax, ymin, ymax, zmin, zmax, n):
    rng = np.random.default_rng()
    x = rng.uniform(xmin, xmax, n)
    y = rng.uniform(ymin, ymax, n)
    z = rng.uniform(zmin, zmax, n)
    return x, y, z
# def xyzInBox

def writeHEPEVT( tRad, snEvents, prefix, subrun, ymin, ymax, zmin, zmax):

    hepevtName = f'{prefix}_{subrun:02d}.hepevt'
    nParticlePerVtx = 1

    with open(hepevtName, 'w') as f:
        lastEvt = None
        iVertex = 0
        for iPart, iEvt in enumerate( tRad ):
            if ( iEvt.y > ymax or iEvt.y < ymin or iEvt.z > zmax or iEvt.z < zmin ): continue

            thisEvt = iEvt.event-1

            # First vertex is always the supernova neutrino interaction
            if thisEvt != lastEvt:    
                iVertex = 0
                snEvt = snEvents[thisEvt]
                nParticlesSN = len(snEvt[1])
                # print( snEvt[1])
                f.write( f'{thisEvt} {iVertex} {nParticlesSN}\n')
                for p in snEvt[1]:
                    f.write(f"{p['ISTHEP']} {p['IDHEP']} {p['JMOHEP1']} {p['JMOHEP2']} {p['JDAHEP1']} {p['JDAHEP2']} "
                            f"{p['PHEP1']:.6f} {p['PHEP2']:.6f} {p['PHEP3']:.6f} {p['PHEP4']:.6f} {p['PHEP5']:.6f} "
                            f"{p['VHEP1']:.6f} {p['VHEP2']:.6f} {p['VHEP3']:.6f} {p['VHEP4']}\n")

            # From the second vertex (iVertex = 1), start filling the radiological background
            else:
                f.write( f'{thisEvt} {iVertex} {nParticlePerVtx}\n')
            
                # ISTHEP IDHEP JMOHEP1 JMOHEP2 JDAHEP1 JDAHEP2 PHEP1 PHEP2 PHEP3 PHEP4 PHEP5 VHEP1 VHEP2 VHEP3 VHEP4
                # final-state particle
                ISTHEP = 1
                IDHEP = iEvt.pdg_code
                # The JMOHEP1, JMOHEP2, JDAHEP1, and JDAHEP2 entries record the indices (between 1 and NHEP, inclusive) 
                # of particles in the event record that correspond to the first mother, second mother, first daughter, 
                # and last daughter of the current particle, respectively. 
                JMOHEP1 = 0
                JMOHEP2 = 0
                JDAHEP1 = 0
                JDAHEP2 = 0

                PHEP1 = iEvt.px
                PHEP2 = iEvt.py
                PHEP3 = iEvt.pz
                PHEP4 = iEvt.Etot
                PHEP5 = iEvt.mass
                VHEP1 = iEvt.x
                VHEP2 = iEvt.y
                VHEP3 = iEvt.z
                VHEP4 = iEvt.t
                f.write( f'{ISTHEP} {IDHEP} {JMOHEP1} {JMOHEP2} {JDAHEP1} {JDAHEP2} {PHEP1} {PHEP2} {PHEP3} {PHEP4} {PHEP5} {VHEP1} {VHEP2} {VHEP3} {VHEP4}\n')
                
            iVertex += 1
            lastEvt = thisEvt
    return
# def writeHEPEVT

if __name__ == "__main__":

    # Radiological files
    radDir = '/Users/yuntse/data/lartpc_rd/gampix/gen'
    radFile = f'{radDir}/dune/fullgeoanatruth-vd1x8x14-1000.root'
    outPrefix = f'{radDir}/nurad/nueArCC_garching_nh_mxpypzDir-rad-vd-reduced'

    fRad = ROOT.TFile( radFile )
    tRad = fRad.Get("fullgeoanatruth/FullGeoAnaTruth")

    # Division
    nY = 4
    nZ = 7
    nEvents = 1000

    # cm
    YLow = -600.
    ZLow = 0.
    YLength = 300.
    ZLength = 300.
    xCathode = -325.
    xCRP = 325.

    # Supernova neutrinos
    snDir = '/Users/yuntse/data/lartpc_rd/gampix/gen/marley_raw'
    snInPrefix = 'nueArCC_garching_nh_mxpypzDir'

    for iy in range( nY ):
        for iz in range( nZ ):

            ymin = YLength*float(iy) + YLow
            ymax = YLength*float(iy+1) + YLow
            zmin = ZLength*float(iz) + ZLow
            zmax = ZLength*float(iz+1) + ZLow
            subrun = iy*nZ + iz

            snFile = f'{snDir}/{snInPrefix}_{subrun:02d}.hepevt'

            print( f'Processing Subrun {subrun}, radiologicals in {ymin} < y < {ymax}cm, {zmin} < z < {zmax}cm...')
            marleyEvts = readHepevt( snFile, 0)
            x, y, z = xyzInBox(xCathode, xCRP, ymin, ymax, zmin, zmax, nEvents)
            t = np.full(nEvents, 0)
            snEvents = assignXYZT(marleyEvts, x, y ,z, t, 0)
            writeHEPEVT( tRad, snEvents, outPrefix, subrun, ymin, ymax, zmin, zmax)