#!/usr/bin/env python3

# Copyright 2019-2020 CERN and copyright holders of ALICE O2.
# See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
# All rights not expressly granted are reserved.

# This software is distributed under the terms of the GNU General Public
# License v3 (GPL Version 3), copied verbatim in the file "COPYING".

# In applying this license CERN does not waive the privileges and immunities
# granted to it by virtue of its status as an Intergovernmental Organization
# or submit itself to any jurisdiction.


"""
file: prepare_samples.py
brief: script to prepare data samples for ML training
usage: python3 prepare_samples.py CONFIG
author: Alexandre Bigot <alexandre.bigot@cern.ch>, Strasbourg University
"""

import argparse
import os
import sys

import yaml  # pylint: disable=import-error

# pylint: disable=import-error,no-name-in-module
from ROOT import EnableImplicitMT, RDataFrame, TChain, TFile

ID_BKG = 0  # index of background class in labeling configuration


def get_path_input_trees(name_input_files, name_tree):
    """
    Helper method to get the path of the input TTrees

    Parameters
    -----------------
    - name_input_files: list with input file names
    - name_tree: the (common) name of the TTree

    Returns
    -----------------
    - path_input_trees: path of input trees
    """

    name_df = "DF"
    delimiter = ";"
    path_input_trees = []
    for name_file in name_input_files:
        myfile = TFile.Open(name_file, "READ")
        name_dirs = []
        keys = myfile.GetListOfKeys()
        for key in keys:
            # fill list only if DF in the dir name (no parent files dir)
            name = key.GetName()
            if name_df in name:
                name_dirs.append(name)
        myfile.Close()

        if len(name_dirs) > 2:
            print(
                "\033[91mERROR: too many DFs in <name of root file>!"
                "Run o2-aod-merger with --max-size <great number> to get only one DF.\033[0m"
            )
            sys.exit()
        if len(name_dirs) == 2:
            # check if the two DFs are are the same (o2-aod-merger artefact)
            name0 = name_dirs[0].split(delimiter)[0]
            name1 = name_dirs[1].split(delimiter)[0]
            if name0 != name1:
                print(
                    "\033[91mERROR: too many DFs in <name of root file>!"
                    "Run o2-aod-merger with --max-size <great number> to get only one DF.\033[0m"
                )
                sys.exit()

        path_input_trees.append(name_file + "/" + name_dirs[0] + "/" + name_tree)

    return path_input_trees


# pylint: disable=too-many-locals
def main(cfg):
    """
    Main function

    Parameters
    -----------------
    - config: dictionary with config read from a yaml file
    """

    channel = cfg["channel"]
    labels = cfg["labels"]
    col_tag = cfg["col_tag"]

    # input configurables
    name_in_files = cfg["prepare_samples"]["input"]["files"]
    name_in_tree = cfg["prepare_samples"]["input"]["tree_name"]

    # output configurables
    inv_mass_sidebands = cfg["filt_bkg_mass"]
    preselections = cfg["preselections"]
    cols_to_remove = cfg["cols_to_remove"]
    name_out_dirs = cfg["output"]["dirs"]
    name_out_tree = cfg["output"]["tree_name"]

    # configure access to the input TTrees
    path_input_trees = get_path_input_trees(name_in_files, name_in_tree)

    # define preselections, if enabled
    presels = ""
    if preselections is not None:
        for i, presel in enumerate(preselections):
            presels += presel
            if i < len(preselections) - 1:
                presels += " && "

    # extract dataframes from input
    counter_outdir = 0
    for path_input_tree in path_input_trees:
        print("\033[32mExtracting dataframes from input " f"{name_in_files[counter_outdir]}\033[0m")
        chain = TChain()
        chain.Add(path_input_tree)

        # define dataframe from the input TTrees
        EnableImplicitMT(8)  # tell ROOT to go parallel
        df = RDataFrame(chain)

        cols_to_keep = list(df.GetColumnNames())
        for col in cols_to_remove:
            cols_to_keep.remove(col)

        # apply preselections, if enabled
        if preselections is not None:
            df = df.Filter(presels)

        # divide dataframe into classes and save them in flagged .root files
        for idx, _ in enumerate(labels):
            name_out_file = f"{name_out_dirs[counter_outdir]}/{labels[idx]}_{channel}.root"
            if os.path.isfile(name_out_file):
                print(f"\033[93mWARNING: Output file {name_out_file} already exists, overwrite ongoing!\033[0m")
            if idx == ID_BKG:
                df.Filter(f"{col_tag} == {idx}").Filter(inv_mass_sidebands).Snapshot(
                    name_out_tree, name_out_file, cols_to_keep
                )
            else:
                df.Filter(f"{col_tag} == {idx}").Snapshot(name_out_tree, name_out_file, cols_to_keep)

    counter_outdir += 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument("config", metavar="text", default="config.yml", help="config file name for ml")
    args = parser.parse_args()

    print("Loading configuration: ...", end="\r")
    with open(args.config, "r", encoding="utf-8") as yml_cfg:
        configuration = yaml.load(yml_cfg, yaml.FullLoader)
    print("Loading configuration: Done!")

    main(configuration)
