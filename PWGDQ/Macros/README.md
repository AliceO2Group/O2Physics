# DQFitter class
Collection of methods for the fit to the dilepton invariant mass distributions. The class is based on RooFit and allows to fit both binned and unbinned distributions.
- DQFitter.py: class collecting all the methods. It allows to fit single invariant mass spectra or perform a multi-trial fit
- fit_library: directory containing a collection of PDFs that are used in dilepton analyses. To add another PDF the user has to follow the same template as those already included
- tutorial.py: simple script which generates the tutorial sample and fits it
- configFit.json: configuration of the fit. Contains all the parameters and PDFs which should be used in the fit

## Tutorial
Tutorial showing how to use the DQFitter class:
- Generate the tutorial data sample (hisogram and tree):
  ```ruby
  python tutorial.py configFit.json --gen_tutorial
  ```
- Run the fit on the tutorial data sample:
  ```ruby
  python tutorial.py configFit.json --run_fit
  ```

The results of the fits are saved into `FitResults.root`