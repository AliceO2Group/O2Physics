# Multiplicity/centrality tools in O2/O2Physics

This directory serves to aggregate all files necessary for multiplicity
and centrality calibration going from pp to Pb-Pb. It offers functionality
to perform simple slicing in percentiles, such as what is done in
proton-proton collisions, as well as the Glauber + analytical NBD fitter
used to perform anchoring and hadronic event distribution estimates in
nucleus-nucleus collisions.

## Files
* `README.md` this readme
* `CMakeLists.txt` definition of source files that need compiling
* `multCalibrator.cxx/h` a class to do percentile slicing of a given histogram. Used for all systems.
* `multMCCalibrator.cxx/h` a class to perform data-to-mc matching of average Nch.
* `multGlauberNBDFitter.cxx/h` a class to do glauber fits.
* `multModule.h` a class to perform calculations of multiplicity and centrality tables for analysis. Meant to be used inside the main core service wagon 'multcenttable'.

