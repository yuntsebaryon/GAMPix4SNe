GAMPix for Supernova Neutrino Measurements in DUNE
==================================================

## Directories

- `MCProd`: produce the MC samples
    - `gen`: scripts producing the events at the level of primary particles
    - `g4`: scripts running `edep-sim` (a GEANT4 wrapper)

## Environment setups

### Edep-sim

```shell
   source /Users/yuntse/opt/root-v6.36.10/bin/thisroot.sh
   source /Users/yuntse/opt/geant4-v10.7.4/bin/geant4.sh
   cd source/edep-sim-origin
   source setup.sh
   export export CMAKE_PREFIX_PATH=${EDEP_ROOT}/${EDEP_TARGET}
   cd /Users/yuntse/work/lartpc_rd/GAMPix4SNe/MCProd/g4
   python runG4.py --dir /Users/yuntse/data/lartpc_rd/gampix
```

### Dump the edep-sim output to a h5 file

```shell
   source /Users/yuntse/opt/root-v6.36.10/bin/thisroot.sh
   export LD_LIBRARY_PATH="/Users/yuntse/opt/edep-sim-origin/edep-gcc-21.0.0-arm64-apple-darwin25.3.0/lib:/Users/yuntse/opt/root-v6.36.10/lib"
   python dumpTreeNew.py <input.root> <output.h5> True
```
alternatively, in the end
```shell
   python runDumpTree.py --dir /Users/yuntse/data/lartpc_rd/gampix
```

--------------------------------------------------------------------------------------------

## Install the packages

To install `edep-sim`, you will need `ROOT` and `GEANT4` installed.
In the instruction below, most of my installation is at `/Users/yuntse/opt`.

### Install ROOT

```shell
   cmake -DCMAKE_INSTALL_PREFIX=/Users/yuntse/opt/root-v6.36.10 -Dbuiltin_vdt=ON /Users/yuntse/source/root
   cmake --build . -- -j 16
   make install
```

### Compile GEANT4 10.7.4

- Note that you can only install `GEANT4 v10`.  `GEANT4 v11` is not compatible with `edep-sim`.
- I need to install `xerces-c` to parse the gdml files (detector geometry).  I also install `expat` and `zlib`.

```shell
    cmake -DCMAKE_INSTALL_PREFIX=/Users/yuntse/opt/geant4-v10.7.4 -DCMAKE_POLICY_DEFAULT_CMP0135=NEW -DGEANT4_INSTALL_DATA=ON -DGEANT4_USE_GDML=ON -DGEANT4_USE_QT=ON -DGEANT4_USE_RAYTRACER_X11=ON -DGEANT4_USE_OPENGL_X11=ON -DGEANT4_BUILD_MULTITHREADED=ON -DCMAKE_PREFIX_PATH=/Users/yuntse/opt/xerces-c -DEXPAT_INCLUDE_DIR=/Users/yuntse/opt/expat-2.6.2/include -DEXPAT_LIBRARY=/Users/yuntse/opt/expat-2.6.2/lib/libexpat.dylib -DGEANT4_USE_SYSTEM_ZLIB=ON /Users/yuntse/source/geant4-v10.7.4
    make -j8
    make install
```

### Compile edep-sim

- To keep the same version as mine, check out my fork of [edep-sim](https://github.com/yuntsebaryon/edep-sim-origin) and use the tag `cv2.0`.

```shell
   source /Users/yuntse/opt/root-v6.36.10/bin/thisroot.sh
   source /Users/yuntse/opt/geant4-v10.7.4/bin/geant4.sh
   cd source/edep-sim-origin
   source setup.sh
   export CMAKE_PREFIX_PATH=/Users/yuntse/opt/edep-sim-origin/${EDEP_TARGET}
   cd ../edep-sim-origin-build/
   cmake -DCMAKE_POLICY_VERSION_MINIMUM=3.5 -DCMAKE_INSTALL_PREFIX=/Users/yuntse/opt/edep-sim-origin/edep-gcc-21.0.0-arm64-apple-darwin25.3.0 -DXercesC_INCLUDE_DIR=/Users/yuntse/opt/xerces-c/include -DXercesC_LIBRARY_RELEASE=/Users/yuntse/opt/xerces-c/lib/libxerces-c.dylib /Users/yuntse/source/edep-sim-origin
   make
   make doc
   make install
```

