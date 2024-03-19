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
brief: script for the preparation of the samples for the training of ML models to be used in PWGHF/D2H
note: inspired by EventFiltering/PWGHF/Macros/prepare_samples.py
usage: python3 prepare_samples.py CONFIG
author: Alexandre Bigot <alexandre.bigot@cern.ch>, Strasbourg University
"""

import argparse
import os

import numpy as np  # pylint: disable=import-error
import uproot  # pylint: disable=import-error
import yaml  # pylint: disable=import-error
from alive_progress import alive_bar  # pylint: disable=import-error

MAX_FILES = 1000


def enforce_list(x):
    """
    Helper method to enforce list type

    Parameters
    ----------
    - x: a string or a list of string

    Returns
    ----------
    - x_list if x was not a list (and not None), x itself otherwise
    """

    if not isinstance(x, list):
        # handle possible spaces in config file entry
        x = x.split(",")
        for i, element in enumerate(x):
            x[i] = element.strip()

    return x


# pylint: disable=too-many-locals
def combine_tpc_or_tof(df, mass_hypos, n_prongs):
    """
    Method to remove track if it has no TOF nor TPC information
    and to compute combined nSigma
    Parameters
    -----------------
    - df: pandas dataframe containing all candidates with fNSigTPC/TOF columns
    - mass_hypos: list of strings containing names of mass hypothesis
    - n_prongs: int number of prongs in the studied decay channel

    Outputs
    -----------------
    - df: pandas dataframe containing all candidates with fNSigCombined columns
    """

    nsigma_tolerance = 0.1
    nsigma_dummy = -999.0 + nsigma_tolerance
    logic = "or"

    nsigma_cols = []
    cols_to_remove = []
    nsigma_combined_cols = {}

    # filter rows with TPC or TOF information
    for mass_hypo in mass_hypos:
        for prong in range(n_prongs):
            suffix = mass_hypo + f"{prong}"
            nsigma_tpc_col = "fNSigTpc" + suffix
            nsigma_tof_col = "fNSigTof" + suffix
            nsigma_cols.append([nsigma_tpc_col, nsigma_tof_col])
            cols_to_remove.append(nsigma_tpc_col)
            cols_to_remove.append(nsigma_tof_col)
            nsigma_combined_cols[suffix] = []
            # apply TPC OR TOF logic
            df.query(f"{nsigma_tpc_col} > {nsigma_dummy} {logic} {nsigma_tof_col} > {nsigma_dummy}", inplace=True)

    for row in list(df.index):
        for mass_hypo in mass_hypos:
            for prong in range(n_prongs):
                suffix = mass_hypo + f"{prong}"
                nsigma_tpc_col = "fNSigTpc" + suffix
                nsigma_tof_col = "fNSigTof" + suffix

                nsigma_tpc = df.at[row, nsigma_tpc_col]
                nsigma_tof = df.at[row, nsigma_tof_col]

                if nsigma_tpc > nsigma_dummy:
                    if nsigma_tof > nsigma_dummy:  # TPC and TOF
                        nsigma_combined_cols[suffix].append(
                            np.sqrt(0.5 * (nsigma_tpc * nsigma_tpc + nsigma_tof * nsigma_tof))
                        )
                    else:  # only TPC
                        nsigma_combined_cols[suffix].append(np.abs(nsigma_tpc))
                else:  # only TOF
                    if nsigma_tof > nsigma_dummy:
                        nsigma_combined_cols[suffix].append(np.abs(nsigma_tof))

    # remove columns with nsigma information
    df.drop(columns=cols_to_remove, inplace=True)
    # insert column with combined nsigma
    for mass_hypo in mass_hypos:
        for prong in range(n_prongs):
            suffix = mass_hypo + f"{prong}"
            df.insert(0, "fNSigComb" + suffix, nsigma_combined_cols[suffix])

    return df


def divide_df_for_origin(df, cols_to_remove=None):
    """
    Method to divide a dataframe in three (prompt, non-prompt, bkg)

    Parameters
    -----------------
    - df: pandas dataframe containing all candidates with fOriginMcRec column
    - cols_to_remove: columns to be removed from output dataframes

    Outputs
    -----------------
    - df_bkg: pandas dataframe containing only background candidates
    - df_prompt: pandas dataframe containing only prompt signal
    - df_nonprompt: pandas dataframe containing only non-prompt signal
    """

    if cols_to_remove is None:
        cols_to_remove = ["fOriginMcRec"]

    df_bkg = df.query("fOriginMcRec == 0")
    df_prompt = df.query("fOriginMcRec == 1")
    df_nonprompt = df.query("fOriginMcRec == 2")

    cols_to_keep = list(df_prompt.columns)
    for col in cols_to_remove:
        cols_to_keep.remove(col)
    df_bkg = df_bkg[cols_to_keep]
    df_prompt = df_prompt[cols_to_keep]
    df_nonprompt = df_nonprompt[cols_to_keep]

    return df_bkg, df_prompt, df_nonprompt


def main(cfg):
    """
    Main function

    Parameters
    -----------------
    - config: dictionary with configs
    """

    # import configurables
    channel = cfg["channel"]
    labels = enforce_list(cfg["labels"])
    input_files = enforce_list(cfg["prepare_samples"]["input"]["files"])
    input_tree_name = cfg["prepare_samples"]["input"]["tree_name"]
    combine_pid_vars = cfg["prepare_samples"]["pid"]["combine_vars"]
    downscale_bkg = cfg["prepare_samples"]["downscale_bkg"]
    inv_mass_sidebands = cfg["prepare_samples"]["filt_bkg_mass"]
    seed_split = cfg["seed_split"]
    force = cfg["prepare_samples"]["output"]["force"]

    mass_hypos, n_prongs = None, None
    if combine_pid_vars:
        mass_hypos = enforce_list(cfg["prepare_samples"]["pid"]["mass_hypos"])
        n_prongs = cfg["prepare_samples"]["pid"]["n_prongs"]

    df_tot = None

    with alive_bar(len(input_files[:MAX_FILES])) as bar_alive:
        for file in input_files[:MAX_FILES]:
            print(f"\033[32mExtracting dataframes from input " f"{file}\033[0m")

            file_root = uproot.open(file)
            indir = os.path.split(file)[0]

            is_filtered = False

            for exfile in os.listdir(indir):
                if f"{channel}.parquet.gzip" in exfile:
                    is_filtered = True
                    break
            if not is_filtered or force:
                list_of_df = []
                for tree_name in file_root.keys():
                    if input_tree_name in tree_name:
                        list_of_df.append(f"{file}:{tree_name}")
                df_tot = uproot.concatenate(list_of_df, library="pd")

                if combine_pid_vars:
                    df_tot = combine_tpc_or_tof(df_tot, mass_hypos, n_prongs)

                df_bkg, df_prompt, df_nonprompt = divide_df_for_origin(df_tot)

                df_bkg = df_bkg.sample(frac=downscale_bkg, random_state=seed_split)
                df_bkg.query(inv_mass_sidebands, inplace=True)  # select bkg candidates
                df_bkg.to_parquet(os.path.join(indir, f"{labels[0]}_{channel}.parquet.gzip"), compression="gzip")
                df_prompt.to_parquet(os.path.join(indir, f"{labels[1]}_{channel}.parquet.gzip"), compression="gzip")
                if not df_nonprompt.empty:
                    df_nonprompt.to_parquet(
                        os.path.join(indir, f"{labels[2]}_{channel}.parquet.gzip"), compression="gzip"
                    )

                df_tot = None

            bar_alive()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument("config", metavar="text", default="config.yml", help="config file name for ml")
    args = parser.parse_args()

    print("Loading analysis configuration: ...", end="\r")
    with open(args.config, "r", encoding="utf-8") as yml_cfg:
        config = yaml.load(yml_cfg, yaml.FullLoader)
    print("Loading analysis configuration: Done!")

    main(config)
