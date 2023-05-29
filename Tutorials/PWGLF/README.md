# This is the base for the PWGLF tutorials for the O2AT

The tutorial (17-28 Apr 2023) can be still used and is a reference for the LF analyses.
It is built as a set of steps. Each step adds a level of complexity and is built in a separate executable.
The executables are defined in the `CMakeLists.txt`.

## Files
* `README.md` this readme
* `CMakeLists.txt` here are defined the source files to compile
* `strangeness_step0.cxx` Starting point: loop over all V0s and fill invariant mass histogram
* `strangeness_step1.cxx` Apply selections on topological variables of V0s
* `strangeness_step2.cxx` Apply PID selections on V0 daughter tracks
* `strangeness_step3.cxx` Check the MC information of the V0s and verify with the PID information of daughter tracks
