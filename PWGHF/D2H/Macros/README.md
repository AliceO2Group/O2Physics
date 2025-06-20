# Invariant mass fitter
Invariant mass fitter is implemented in the class `HFInvMassFitter` (`HFInvMassFitter.cxx/h` files), and its run is executed via `runMassFitter.C` macro.
Fitter is configured in `config_massfitter.json`.\
The fitter is **not** a part of O2Physics source code.

## Dependencies
1. ROOT

2. RapidJSON. Download the header-only (no compilation needed) RapidJSON library, see the link https://rapidjson.org.

## How to run
### As a ROOT macro
The `runMassFitter.C` can be compiled as ROOT macro.
```
source path-to-root-install/bin/thisroot.sh
export ROOT_INCLUDE_PATH=$ROOT_INCLUDE_PATH:path-to-json-include
root -l -x -b -q "HFInvMassFitter.cxx" "runMassFitter.C(\"config_massfitter.json\")"
```
If you have O2Physics compilation and enter into its environment there is no need to set environment variables (skip first two lines above).

### As a CMake project
It is also possible to compile the fitter as a CMake project or insert it into existing one if any.
Use the `CMakeLists_HFInvMassFitter.txt` (rename it into `CMakeLists.txt` before usage).\
Compile the fitter with the following steps:
```
mkdir build
cd build
source path-to-root-install/bin/thisroot.sh
cmake -DHFFITTER_RAPIDJSON_INCLUDE_DIRS=path-to-json-include ..
make
```
and run the fitter:
```
./runMassFitter ../config_massfitter.json
```
### Directly from the terminal
Compile the fitter with the following steps:
```
mkdir build
cd build
source path-to-root-install/bin/thisroot.sh

# Generate ROOT dictionary:
rootcling -f G__HFInvMassFitter.cxx -c ../HFInvMassFitter.h ../HFInvMassFitterLinkDef.h

# Compile source code:
g++ -fPIC -I$(root-config --incdir) -I path-to-json-include -c ../HFInvMassFitter.cxx ../runMassFitter.C G__HFInvMassFitter.cxx

# Link the executable:
g++ -o runMassFitter HFInvMassFitter.o runMassFitter.o G__HFInvMassFitter.o $(root-config --libs)  -lRooFit -lRooFitCore -lEG
```
and run the fitter:
```
./runMassFitter ../config_massfitter.json
```
