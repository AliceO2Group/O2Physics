# This is the base for the PWGLF tutorials for the O2AT

The tutorial (14-18 Oct 2024) can be used as a reference for the LF analysis.
It is built as a set of steps. Each step adds a level of complexity and is built in a separate executable.

The tutorial is divided into two directories depending on the collision system:
## pp
This repository contains the codes to analyse V0 and cascade particles in proton-proton collisions:
* `CMakeLists.txt` here are defined the source files to compile
* `strangeness_step0.cxx` Starting point: loop over all V0s and fill invariant mass histogram
* `strangeness_step1.cxx` Apply selections on topological variables of V0s
* `strangeness_step2.cxx` Apply PID selections on V0 daughter tracks
* `strangeness_step3.cxx` Check the MC information of the V0s and verify with the PID information of daughter tracks

The introduction and hands-on sessions can be found here: https://indico.cern.ch/event/1326201/

## PbPb
This repository contains the code to analyse cascade particles in Pb-Pb collisions.
The tutorial revolves around two aspects:
i) the derived data production in the directory DerivedDataProduction,
ii) the analysis of the derived data in the directory Analysis, containing:
* `CMakeLists.txt` here are defined the source files to compile
* `strangeness_pbpb_step0.cxx` Starting point: loop over all cascades and fill invariant mass histogram
* `strangeness_pbpb_step1.cxx` Apply selections on topological variables of cascades
* `strangeness_pbpb_step2.cxx` Apply TPC PID selections on cascade daughter tracks
* `strangeness_pbpb_step3.cxx` Apply TOF PID selections on cascade daughter tracks (if available)
* `strangeness_pbpb_step4.cxx` Check the MC information of the cascade

The introduction and hands-on sessions can be found here: https://indico.cern.ch/event/1425820/