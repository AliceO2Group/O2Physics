# Cut-variation method
In order to compute the (non-)prompt fraction of charm hadrons with the cut-variation method, the `compute_fraction_cutvar.py` script can be used.
It requires as inputs `.root` files with transverse-momentum differential raw yields and efficiencies of prompt and non-prompt production stored in `ROOT::TH1F` objects. These have to be passed to the script via a `JSON` config file, such as `config_cutvar_example.json`. The script has to be executed via
```python
python3 compute_fraction_cutvar.py config_cutvar_example.json
```
It produces a `.root` file with output histograms and four `.pdf` files with fitted distributions, efficiencies, covariance matrices, and resulting fractions.
