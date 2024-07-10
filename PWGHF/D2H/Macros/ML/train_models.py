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
file: train_models.py
brief: script for the training of ML models to be used in D2H
note: inspired by EventFiltering/PWGHF/Macros/train_hf_triggers.py and Run2 macros
usage: python3 train_models.py CONFIG
author: Alexandre Bigot <alexandre.bigot@cern.ch>, Strasbourg University
author: Mingyu Zhang <mingyu.zang@cern.ch>, Central China Normal University
"""

import argparse
import os
import pickle
import sys

# pylint: disable=import-error
try:
    from hipe4ml import plot_utils
    from hipe4ml.model_handler import ModelHandler
    from hipe4ml.tree_handler import TreeHandler
except ModuleNotFoundError:
    print("Module 'hipe4ml' is not installed. Please install it to run this macro")

import matplotlib.pyplot as plt  # pylint: disable=import-error
import numpy as np  # pylint: disable=import-error
import pandas as pd  # pylint: disable=import-error
import xgboost as xgb  # pylint: disable=import-error
import yaml  # pylint: disable=import-error
from sklearn.model_selection import train_test_split  # pylint: disable=import-error

try:
    from hipe4ml_converter.h4ml_converter import H4MLConverter
except ModuleNotFoundError:
    print("Module 'hipe4ml_converter' is not installed. Please install it to run this macro")

LABEL_BKG = 0
LABEL_PROMPT = 1
LABEL_NONPROMPT = 2

MAX_BKG_FRAC = 0.4  # max of bkg fraction to keep for training


def enforce_list(x):
    """
    Helper method to enforce list type

    Parameters
    -----------------
    - x: a string or a list of string

    Returns
    -----------------
    - x_list if x was not a list (and not None), x itself otherwise
    """

    if not isinstance(x, list):
        # handle possible spaces in config file entry
        x = x.split(",")
        for i, element in enumerate(x):
            x[i] = element.strip()

    return x


# pylint: disable= too-many-instance-attributes, too-few-public-methods
class MlTrainer:
    """
    Class for ML training and testing
    """

    def __init__(self, config):
        """
        Init method

        Parameters
        -----------------
        - config: dictionary with config read from a yaml file
        """

        # input
        self.channel = config["channel"]
        self.seed_split = config["seed_split"]
        self.labels = enforce_list(config["labels"])
        pt_bins_limits = enforce_list(config["data_prep"]["pt_bins_limits"])
        self.pt_bins = [[a, b] for a, b in zip(pt_bins_limits[:-1], pt_bins_limits[1:])]

        self.indirs = config["data_prep"]["indirs"]
        self.tree_name = config["data_prep"]["tree_name"]
        self.infile_extension = "parquet.gzip" if self.tree_name is None else "root"
        self.binary = self.indirs["Nonprompt"] is None
        self.file_lists = {}

        # (hyper)parameters
        self.share = config["data_prep"]["class_balance"]["share"]
        self.bkg_factor = config["data_prep"]["class_balance"]["bkg_factor"]
        self.downsample_bkg_factor = config["data_prep"]["downsample_bkg_factor"]
        self.training_vars = enforce_list(config["ml"]["training_vars"])
        self.test_frac = config["data_prep"]["test_fraction"]
        self.vars_to_draw = config["plots"]["extra_columns"] + self.training_vars
        self.name_pt_var = config["data_prep"]["name_pt_var"]

        self.raw_output = config["ml"]["raw_output"]
        self.roc_auc_approach = config["ml"]["roc_auc_approach"]
        self.roc_auc_average = config["ml"]["roc_auc_average"]
        self.score_metric = "roc_auc"
        self.hyper_pars = config["ml"]["hyper_pars"]
        self.hyper_pars_opt = config["ml"]["hyper_pars_opt"]

        # output
        self.outdir = config["output"]["dir"]
        self.extension = config["plots"]["extension"]
        self.column_to_save_list = config["output"]["column_to_save_list"]
        self.log_file = config["output"]["log_file"]

    def __check_input_consistency(self):
        """
        Helper method to check self consistency of inputs
        """

        # class balance
        if self.share not in ("equal", "all_signal"):
            print(f"\033[91mERROR: class_balance option {self.share} not implemented\033[0m")
            sys.exit()
        if self.share == "all_signal" and len(self.bkg_factor) != len(self.pt_bins):
            print("\033[91mERROR: bkg_factor must be defined for each pT bin!\033[0m")
            sys.exit()
        # training
        if self.training_vars is None:
            print("\033[91mERROR: training columns must be defined!\033[0m")
            sys.exit()
        # hyper-parameters options
        if not isinstance(self.hyper_pars, list):
            print("\033[91mERROR: hyper-parameters must be defined or be a list containing an empty dict!\033[0m")
            sys.exit()
        if not isinstance(self.hyper_pars[0], dict):
            print("\033[91mERROR: hyper-parameters must be a list of dict!\033[0m")
            sys.exit()
        if len(self.hyper_pars) != len(self.pt_bins):
            print("\033[91mERROR: hyper-parameters definition does not match pT binning!\033[0m")
            sys.exit()
        if not isinstance(self.hyper_pars_opt["hyper_par_ranges"], dict):
            print("\033[91mERROR: hyper_pars_opt_config must be defined!\033[0m")
            sys.exit()
        if not self.binary:
            if not (self.roc_auc_average in ["macro", "weighted"] and self.roc_auc_approach in ["ovo", "ovr"]):
                print("\033[91mERROR: selected ROC configuration is not valid!\033[0m")
                sys.exit()

            if self.roc_auc_average == "weighted":
                self.score_metric += f"_{self.roc_auc_approach}_{self.roc_auc_average}"
            else:
                self.score_metric += f"_{self.roc_auc_approach}"

    def __fill_list_input_files(self):
        """
        Helper method to fill a dictionary with lists of input files for each class
        """

        file_lists = {}
        for cand_type in self.indirs:
            file_lists[cand_type] = []
            if self.indirs[cand_type] is None:
                continue
            for indir in self.indirs[cand_type]:
                file = os.path.join(indir, f"{cand_type}_{self.channel}.{self.infile_extension}")
                if os.path.isfile(file):
                    file_lists[cand_type].append(file)
                else:
                    print(
                        "\033[91mERROR: missing file, did you prepare the samples?\n"
                        "If not, do python3 prepare_samples.py CONFIG\033[0m"
                    )
                    sys.exit()

        self.file_lists = file_lists

    def __get_sliced_dfs(self):
        """
        Helper method to get pT-sliced dataframes for each class

        Returns
        -----------------
        - hdl_bkg: pandas dataframe containing only background candidates
        - hdl_prompt: pandas dataframe containing only prompt signal
        - hdl_nonprompt: pandas dataframe containing only non-prompt signal
        """

        print("Loading and preparing data files: ...", end="\r")

        self.__fill_list_input_files()

        hdl_bkg = TreeHandler(file_name=self.file_lists[self.labels[0]], tree_name=self.tree_name)
        hdl_prompt = TreeHandler(file_name=self.file_lists[self.labels[1]], tree_name=self.tree_name)
        hdl_nonprompt = (
            None if self.binary else TreeHandler(file_name=self.file_lists[self.labels[2]], tree_name=self.tree_name)
        )

        hdl_prompt.slice_data_frame(self.name_pt_var, self.pt_bins, True)
        if hdl_nonprompt is not None:
            hdl_nonprompt.slice_data_frame(self.name_pt_var, self.pt_bins, True)
        hdl_bkg.slice_data_frame(self.name_pt_var, self.pt_bins, True)

        print("Loading and preparing data files: Done!")
        return hdl_bkg, hdl_prompt, hdl_nonprompt

    # pylint: disable=too-many-statements, too-many-branches, too-many-arguments, too-many-locals
    def __data_prep(self, df_bkg, df_prompt, df_nonprompt, pt_bin, out_dir, bkg_factor):
        """
        Helper method for pt-dependent data preparation

        Parameters
        -----------------
        - df_bkg: pandas dataframe containing only background candidates
        - df_prompt: pandas dataframe containing only prompt signal
        - df_nonprompt: pandas dataframe containing only non-prompt signal
        - pt_bin: pT bin
        - out_dir: output directory
        - bkg_factor: multiplier for (n_prompt + n_nonprompt) used to determine n_cand_bkg in the 'all_signal' option

        Returns
        -----------------
        - train_test_data: list containing train/test sets and the associated model predictions
        """

        n_prompt = len(df_prompt)
        n_nonprompt = len(df_nonprompt)
        n_bkg = len(df_bkg)
        if self.binary:
            log_available_cands = (
                f"\nNumber of available candidates "
                f"in {pt_bin[0]} < pT < {pt_bin[1]} GeV/c: \n   "
                f"Signal: {n_prompt}\n   Bkg: {n_bkg}"
            )
        else:
            log_available_cands = (
                f"\nNumber of available candidates "
                f"in {pt_bin[0]} < pT < {pt_bin[1]} GeV/c: \n   "
                f"Prompt: {n_prompt}\n   Nonprompt: {n_nonprompt}\n   Bkg: {n_bkg}"
            )
        print(log_available_cands)

        if self.share == "equal":
            if self.binary:
                n_cand_min = min([n_prompt, n_bkg])
                bkg_fraction = n_cand_min / n_bkg
                n_bkg = n_prompt = n_cand_min
            else:
                n_cand_min = min([n_prompt, n_nonprompt, n_bkg])
                bkg_fraction = n_cand_min / n_bkg
                n_bkg = n_prompt = n_nonprompt = n_cand_min
            log_share = (
                "\nKeep the same number of candidates for each class, "
                "chosen as the minimal number of candidates among all classes."
            )

        elif self.share == "all_signal":
            n_cand_bkg = int(min([n_bkg, (n_prompt + n_nonprompt) * bkg_factor]))
            if self.binary:
                log_share = (
                    f"\nKeep all signal and use {n_cand_bkg} bkg candidates "
                    f"for training and testing ({1 - self.test_frac}-{self.test_frac})"
                )
            else:
                log_share = (
                    f"\nKeep all prompt and nonprompt and use {n_cand_bkg} "
                    "bkg candidates for training and testing "
                    f"({1 - self.test_frac}-{self.test_frac})"
                )
            bkg_fraction = n_cand_bkg / n_bkg
            n_bkg = n_cand_bkg

        else:
            print(f"\033[91mERROR: class_balance option {self.share} not implemented\033[0m")
            sys.exit()

        print(log_share)

        log_bkg_fraction = (
            "\nFraction of original (i.e. from original dataset) bkg candidates used for ML: "
            f"{100*bkg_fraction*self.downsample_bkg_factor:.2f}%"
        )
        if (1 - self.test_frac) * bkg_fraction * self.downsample_bkg_factor > MAX_BKG_FRAC:
            log_bkg_fraction += (
                f"\n\033[93m\nWARNING: using more than {100*MAX_BKG_FRAC:.0f}% "
                "of original (i.e. from original dataset) bkg available for training!\033[0m"
            )
        print(log_bkg_fraction)

        if self.binary:
            log_training_cands = (
                "\nNumber of candidates used for training and testing: \n   " f"Signal: {n_prompt}\n   Bkg: {n_bkg}\n"
            )
        else:
            log_training_cands = (
                "\nNumber of candidates used for training and testing: \n   "
                f"Prompt: {n_prompt}\n   Nonprompt: {n_nonprompt}\n   Bkg: {n_bkg}\n"
            )

        print(log_training_cands)

        # write logs in log file
        with open(os.path.join(out_dir, self.log_file), "w", encoding="utf-8") as file:
            file.write(log_available_cands)
            file.write(log_share)
            file.write(log_bkg_fraction)
            file.write(log_training_cands)

        df_tot = pd.concat([df_bkg[:n_bkg], df_prompt[:n_prompt], df_nonprompt[:n_nonprompt]], sort=True)

        labels_array = np.array([LABEL_BKG] * n_bkg + [LABEL_PROMPT] * n_prompt + [LABEL_NONPROMPT] * n_nonprompt)
        if 0 < self.test_frac < 1:
            train_set, test_set, y_train, y_test = train_test_split(
                df_tot, labels_array, test_size=self.test_frac, random_state=self.seed_split
            )
        else:
            print("ERROR: test_fraction must belong to ]0,1[")
            sys.exit(0)

        train_test_data = [train_set, y_train, test_set, y_test]
        del df_tot  # release memory

        # safety
        if len(np.unique(train_test_data[3])) != len(self.labels):
            print(
                "\033[91mERROR: The number of labels defined does not match"
                "the number of classes! \nCheck the CONFIG file\033[0m"
            )
            sys.exit()

        # plots
        df_list = [df_bkg, df_prompt] if self.binary else [df_bkg, df_prompt, df_nonprompt]

        # _____________________________________________
        plot_utils.plot_distr(
            df_list, self.vars_to_draw, 100, self.labels, figsize=(12, 7), alpha=0.3, log=True, grid=False, density=True
        )
        plt.subplots_adjust(left=0.06, bottom=0.06, right=0.99, top=0.96, hspace=0.55, wspace=0.55)
        for ext in self.extension:
            plt.savefig(f"{out_dir}/DistributionsAll_pT_{pt_bin[0]}_{pt_bin[1]}.{ext}")
        plt.close("all")
        # _____________________________________________
        corr_matrix_fig = plot_utils.plot_corr(df_list, self.vars_to_draw, self.labels)
        for fig, lab in zip(corr_matrix_fig, self.labels):
            plt.figure(fig.number)
            plt.subplots_adjust(left=0.2, bottom=0.25, right=0.95, top=0.9)
            for ext in self.extension:
                fig.savefig(f"{out_dir}/CorrMatrix{lab}_pT_{pt_bin[0]}_{pt_bin[1]}.{ext}")

        return train_test_data

    # pylint: disable=too-many-statements, too-many-branches
    def __train_test(self, train_test_data, hyper_pars, pt_bin, out_dir):
        """
        Helper method for model training and testing

        Parameters
        -----------------
        - train_test_data: list containing train/test sets and the associated model predictions
        - hyper_pars: default hyper-parameters (can be modified if Optuna enabled)
        - pt_bin: pT bin
        - out_dir: output directory
        """

        n_classes = len(np.unique(train_test_data[3]))
        model_clf = xgb.XGBClassifier(use_label_encoder=False)
        model_hdl = ModelHandler(model_clf, self.training_vars, hyper_pars)

        # hyperparams optimization
        if self.hyper_pars_opt["activate"]:
            print("Performing optuna hyper-parameters optimisation: ...", end="\r")

            with open(os.path.join(out_dir, self.log_file), "a", encoding="utf-8") as file:
                file.write("\nOptuna hyper-parameters optimisation:")
                sys.stdout = file
                model_hdl.optimize_params_optuna(
                    train_test_data,
                    self.hyper_pars_opt["hyper_par_ranges"],
                    cross_val_scoring=self.score_metric,
                    timeout=self.hyper_pars_opt["timeout"],
                    n_jobs=self.hyper_pars_opt["njobs"],
                    n_trials=self.hyper_pars_opt["ntrials"],
                    direction="maximize",
                )
            sys.stdout = sys.__stdout__
            print("Performing optuna hyper-parameters optimisation: Done!")
            print(f"Optuna hyper-parameters:\n{model_hdl.get_model_params()}")
        else:
            model_hdl.set_model_params(hyper_pars)

        # store final hyperparameters in info file
        with open(os.path.join(out_dir, self.log_file), "a", encoding="utf-8") as file:
            file.write(f"\nModel hyperparameters:\n {model_hdl.get_model_params()}")

        # train and test the model with the updated hyper-parameters
        y_pred_test = model_hdl.train_test_model(
            train_test_data,
            True,
            output_margin=self.raw_output,
            average=self.roc_auc_average,
            multi_class_opt=self.roc_auc_approach,
        )

        y_pred_train = model_hdl.predict(train_test_data[0], self.raw_output)

        # Save applied model to test set
        test_set_df = train_test_data[2]
        test_set_df = test_set_df.loc[:, self.column_to_save_list]
        test_set_df["Labels"] = train_test_data[3]

        for pred, lab in enumerate(self.labels):
            if self.binary:
                test_set_df["ML_output"] = y_pred_test
            else:
                test_set_df[f"ML_output_{lab}"] = y_pred_test[:, pred]

        test_set_df.to_parquet(f"{out_dir}/{self.channel}_ModelApplied" f"_pT_{pt_bin[0]}_{pt_bin[1]}.parquet.gzip")

        # save model
        if os.path.isfile(f"{out_dir}/ModelHandler_{self.channel}.pickle"):
            os.remove(f"{out_dir}/ModelHandler_{self.channel}.pickle")
        if os.path.isfile(f"{out_dir}/ModelHandler_onnx_{self.channel}.onnx"):
            os.remove(f"{out_dir}/ModelHandler_onnx_{self.channel}.onnx")

        model_hdl.dump_model_handler(f"{out_dir}/ModelHandler_{self.channel}" f"_pT_{pt_bin[0]}_{pt_bin[1]}.pickle")
        model_conv = H4MLConverter(model_hdl)
        model_conv.convert_model_onnx(1)
        model_conv.dump_model_onnx(f"{out_dir}/ModelHandler_onnx_{self.channel}" f"_pT_{pt_bin[0]}_{pt_bin[1]}.onnx")

        # plots
        # _____________________________________________
        plt.rcParams["figure.figsize"] = (10, 7)
        fig_ml_output = plot_utils.plot_output_train_test(
            model_hdl, train_test_data, 80, self.raw_output, self.labels, True, density=True
        )
        if n_classes > 2:
            for fig, lab in zip(fig_ml_output, self.labels):
                for ext in self.extension:
                    fig.savefig(f"{out_dir}/MLOutputDistr{lab}_pT_{pt_bin[0]}_{pt_bin[1]}.{ext}")
        else:
            for ext in self.extension:
                fig_ml_output.savefig(f"{out_dir}/MLOutputDistr_pT_{pt_bin[0]}_{pt_bin[1]}.{ext}")
        # _____________________________________________
        plt.rcParams["figure.figsize"] = (10, 9)
        fig_roc_curve = plot_utils.plot_roc(
            train_test_data[3], y_pred_test, None, self.labels, self.roc_auc_average, self.roc_auc_approach
        )
        for ext in self.extension:
            fig_roc_curve.savefig(f"{out_dir}/ROCCurveAll_pT_{pt_bin[0]}_{pt_bin[1]}.{ext}")
        pickle.dump(fig_roc_curve, open(f"{out_dir}/ROCCurveAll_pT_{pt_bin[0]}_{pt_bin[1]}.pkl", "wb"))
        # _____________________________________________
        plt.rcParams["figure.figsize"] = (10, 9)
        fig_roc_curve_tt = plot_utils.plot_roc_train_test(
            train_test_data[3],
            y_pred_test,
            train_test_data[1],
            y_pred_train,
            None,
            self.labels,
            self.roc_auc_average,
            self.roc_auc_approach,
        )

        fig_roc_curve_tt.savefig(f"{out_dir}/ROCCurveTrainTest_pT_{pt_bin[0]}_{pt_bin[1]}.pdf")
        # _____________________________________________
        precision_recall_fig = plot_utils.plot_precision_recall(train_test_data[3], y_pred_test, self.labels)
        precision_recall_fig.savefig(f"{out_dir}/PrecisionRecallAll_pT_{pt_bin[0]}_{pt_bin[1]}.pdf")
        # _____________________________________________
        plt.rcParams["figure.figsize"] = (12, 7)
        fig_feat_importance = plot_utils.plot_feature_imp(
            train_test_data[2][train_test_data[0].columns], train_test_data[3], model_hdl, self.labels
        )
        n_plot = n_classes if n_classes > 2 else 1
        for i_fig, fig in enumerate(fig_feat_importance):
            if i_fig < n_plot:
                lab = self.labels[i_fig] if n_classes > 2 else ""
                for ext in self.extension:
                    fig.savefig(f"{out_dir}/FeatureImportance_{lab}_{self.channel}.{ext}")
            else:
                for ext in self.extension:
                    fig.savefig(f"{out_dir}/FeatureImportanceAll_{self.channel}.{ext}")

    def process(self):
        """
        Process function of the class, performing data preparation,
        training, testing, saving the model and important plots
        """

        self.__check_input_consistency()
        df_bkg, df_prompt, df_nonprompt = self.__get_sliced_dfs()

        for i_pt, pt_bin in enumerate(self.pt_bins):
            print(f"\n\033[94mStarting ML analysis --- {pt_bin[0]} < pT < {pt_bin[1]} GeV/c\033[0m")

            out_dir_pt = os.path.join(os.path.expanduser(self.outdir), f"pt{pt_bin[0]}_{pt_bin[1]}")
            if os.path.isdir(out_dir_pt):
                print(
                    (
                        f"\033[93mWARNING: Output directory '{out_dir_pt}' already exists,"
                        " overwrites possibly ongoing!\033[0m"
                    )
                )
            else:
                os.makedirs(out_dir_pt)

            df_pt_nonprompt = pd.DataFrame() if df_nonprompt is None else df_nonprompt.get_slice(i_pt)
            if self.share == "all_signal":
                bkg_factor = self.bkg_factor[i_pt]
            else:
                bkg_factor = None

            train_test_data = self.__data_prep(
                df_bkg.get_slice(i_pt), df_prompt.get_slice(i_pt), df_pt_nonprompt, pt_bin, out_dir_pt, bkg_factor
            )
            self.__train_test(train_test_data, self.hyper_pars[i_pt], pt_bin, out_dir_pt)


def main(cfg):
    """
    Main function

    Parameters
    -----------------
    - config: dictionary with config read from a yaml file
    """
    MlTrainer(cfg).process()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument("config", metavar="text", default="config.yml", help="config file name for ml")
    args = parser.parse_args()

    print("Loading analysis configuration: ...", end="\r")
    with open(args.config, "r", encoding="utf-8") as yml_cfg:
        configuration = yaml.load(yml_cfg, yaml.FullLoader)
    print("Loading analysis configuration: Done!")

    main(configuration)
