# This is the base for the PWGLF tutorials for the O2AT

The tutorial can be still used and is a reference for the LF Resonance analyses.
It is built as a set of steps. Each step adds a level of complexity and is built in a separate executable.
The executables are defined in the `CMakeLists.txt`.

## Files
* `README.md` this readme
* `CMakeLists.txt` here are defined the source files to compile
* `resonance_step0.cxx` Read the resonance table and run basic track loop
* `resonance_step1.cxx` Producing same event invariant mass distribution
* `resonance_step2.cxx` Producing Mixed event invariant mass distribution
* `resonance_step3.cxx` Starting point for MC: loop over all Generated MC particles
* `resonance_step4.cxx` Producing histograms (pt distributions etc.) for Generated MCs
* `resonance_step5.cxx` Loop over all MC Tracks and produce basic QA histograms
* `resonance_step6.cxx` Resonance reconstruction
