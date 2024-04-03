#!/usr/bin/env python

import pandas as pd
import numpy as np
import uproot3
import matplotlib.pyplot as plt
import os
from matplotlib import colors

# input_dir = "data_20_03/"
# input_dir = "data_02_06/"
# input_dir = "data_18_08/"
input_dir = "/root/alice/O2Physics/Tools/PIDML/mlModelHists"
#input_files = os.listdir(input_dir)
#input_files = [file for file in input_files if "pidtracksmc" in file]
input_files = ["AnalysisResults.root"]
input_files
dataframes = []
for input_file in input_files:
    file = uproot3.open(input_dir+input_file)
    for dirname in file:
        dirname = dirname.decode("utf-8")
        pure_dirname = dirname.split(";")[0]
        #print("Pure dirname: ", pure_dirname)
        if pure_dirname.startswith("DF_"):
            tree_data = file["%s/O2pidtracksmc" % (dirname)].pandas.df()
            dataframes.append(tree_data)

data = pd.concat(dataframes,ignore_index=True)
print(data.head())
print(data.columns)

#p = np.sqrt(data.fPx ** 2 + data.fPy ** 2 + data.fPz ** 2)
#data["P"] = p
data["fBeta"].mask(np.isclose(data["fBeta"],-999),inplace=True)
data["fTOFSignal"].mask(np.isclose(data["fTOFSignal"],-999),inplace=True)
data["fTRDPattern"].mask(np.isclose(data["fTRDPattern"],0),inplace=True)
#data["fTRDSignal"].mask(np.isclose(data["fTRDSignal"],-999),inplace=True)
data = data[data["fTPCSignal"]>0]

data.to_csv("LHC18g4_023.csv")
