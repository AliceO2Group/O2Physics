"""
Script for the training of ML models to be used in HF triggers
It requires hipe4ml and hipe4ml_converter to be installed:
  pip install hipe4ml
  pip install hipe4ml_converter

\author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN
\author Alexandre Bigot <alexandre.bigot@cern.ch>, Strasbourg University
\author Biao Zhang <biao.zhang@cern.ch>, CCNU
"""

import os
import sys
import argparse
import pickle
import numpy as np
import pandas as pd
import xgboost as xgb
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt
import yaml

from hipe4ml import plot_utils
from hipe4ml.model_handler import ModelHandler
from hipe4ml.tree_handler import TreeHandler
from hipe4ml_converter.h4ml_converter import H4MLConverter


def is_selected_mass_kk(cand, mass_cut):
    """
    Function that returns True if the mass of the KK is in the selected mass range
    """
    if cand["fDeltaMassKKFirst"] < mass_cut or cand["fDeltaMassKKSecond"] < mass_cut:
        return 1
    return 0


def is_selected_proton_pid(cand, nsigma_tpc, nsigma_tof):
    """
    Function that returns True if the proton PID is in the selected range
    """
    is_proton_first = cand["fNsigmaPrTPC1"] < nsigma_tpc and \
        (cand["fNsigmaPrTOF1"] < nsigma_tof or cand["fNsigmaPrTOF1"] == np.nan)
    is_proton_third = cand["fNsigmaPrTPC3"] < nsigma_tpc and \
        (cand["fNsigmaPrTOF3"] < nsigma_tof or cand["fNsigmaPrTOF3"] == np.nan)
    if is_proton_first or is_proton_third:
        return 1

    return 0


def get_list_input_files(indirs, channel):
    """
    function that returns the list of files

    Parameters
    -----------------
    - indirs: dictionary with lists of input directories for prompt, nonprompt, and bkg
    - channel: decay channel, options:
        D0ToKPi, DPlusToPiKPi, DsToKKPi, LcToPKPi, XicToPKPi

    Outputs
    -----------------
    - file_lists: dictionary with lists of input files for prompt, nonprompt, and bkg
    """

    if channel not in ["D0ToKPi", "DplusToPiKPi", "DsToKKPi", "LcToPKPi", "XicToPKPi"]:
        print(f"ERROR: channel {channel} not implemented, return None")
        return {"Prompt": None, "Nonprompt": None, "Bkg": None}

    file_lists = {}
    for cand_type in indirs:  # pylint: disable=too-many-nested-blocks
        file_lists[cand_type] = []
        for indir in indirs[cand_type]:
            subdirs = os.listdir(indir)
            for subdir in subdirs:
                subdir = os.path.join(indir, subdir)
                if os.path.isdir(subdir):
                    for subsubdir in os.listdir(subdir):
                        subsubdir = os.path.join(subdir, subsubdir)
                        if os.path.isdir(subsubdir):
                            file = os.path.join(
                                subsubdir, f"{cand_type}_{channel}.parquet.gzip")
                            if os.path.isfile(file):
                                file_lists[cand_type].append(file)
                            else:
                                for subsubsubdir in os.listdir(subsubdir):
                                    subsubsubdir = os.path.join(subsubdir, subsubsubdir)
                                    if os.path.isdir(subsubsubdir):
                                        file = os.path.join(
                                            subsubsubdir, f"{cand_type}_{channel}.parquet.gzip")
                                        if os.path.isfile(file):
                                            file_lists[cand_type].append(file)

    return file_lists


# pylint: disable=too-many-statements, too-many-branches, too-many-locals
def data_prep(config):
    """
    function for data preparation

    Parameters
    -----------------
    - config: dictionary with config read from a yaml file
    """

    input_dirs = config["data_prep"]["dirs"]
    channel = config["data_prep"]["channel"]
    test_f = config["data_prep"]["test_fraction"]
    seed_split = config["data_prep"]["seed_split"]
    out_dir = config["output"]["directory"]
    training_vars = config["ml"]["training_vars"]
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    file_lists = get_list_input_files(input_dirs, channel)
    for file_list in file_lists:
        if file_list is None:
            sys.exit()

    hdl_prompt = TreeHandler(file_lists["Prompt"])
    hdl_nonprompt = TreeHandler(file_lists["Nonprompt"])
    hdl_bkg = TreeHandler(file_lists["Bkg"])

    df_prompt = hdl_prompt.get_data_frame()
    df_nonprompt = hdl_nonprompt.get_data_frame()
    df_bkg = hdl_bkg.get_data_frame()

    if config["data_prep"]["preselection"]["enable"]:
        if channel in ["D0ToKPi", "DPlusToPiKPi"]:
            print("\nWARNING: No preselection for D0 and D+ to be applied, skip it\n")
        else:
            if channel == "DsToKKPi":
                mass_cut = config["data_prep"]["preselection"]["delta_mass_kk"]
                df_prompt["fIsSel"] = df_prompt.apply(
                    lambda row: is_selected_mass_kk(row, mass_cut), axis=1)
                df_nonprompt["fIsSel"] = df_nonprompt.apply(
                    lambda row: is_selected_mass_kk(row, mass_cut), axis=1)
                df_bkg["fIsSel"] = df_bkg.apply(
                    lambda row: is_selected_mass_kk(row, mass_cut), axis=1)
            elif channel in ["LcToPKPi", "XicToPKPi"]:
                nsigma_tpc = config["data_prep"]["preselection"]["nsigma_tpc_proton"]
                nsigma_tof = config["data_prep"]["preselection"]["nsigma_tof_proton"]
                df_prompt["fIsSel"] = df_prompt.apply(
                    lambda row: is_selected_proton_pid(row, nsigma_tpc, nsigma_tof), axis=1)
                df_nonprompt["fIsSel"] = df_nonprompt.apply(
                    lambda row: is_selected_proton_pid(row, nsigma_tpc, nsigma_tof), axis=1)
                df_bkg["fIsSel"] = df_bkg.apply(
                    lambda row: is_selected_proton_pid(row, nsigma_tpc, nsigma_tof), axis=1)
            df_prompt = df_prompt.query("fIsSel > 0.5")
            df_nonprompt = df_nonprompt.query("fIsSel > 0.5")
            df_bkg = df_bkg.query("fIsSel > 0.5")

    use_pid = False
    for var in training_vars:
        if "NSigma" in var:
            use_pid = True
    if use_pid:
        df_prompt.replace([np.inf, -np.inf], np.nan, inplace=True)
        df_nonprompt.replace([np.inf, -np.inf], np.nan, inplace=True)
        df_bkg.replace([np.inf, -np.inf], np.nan, inplace=True)
        df_prompt.dropna(inplace=True)
        df_nonprompt.dropna(inplace=True)
        df_bkg.dropna(inplace=True)

    n_prompt = len(df_prompt)
    n_nonprompt = len(df_nonprompt)
    n_bkg = len(df_bkg)
    print("\nNumber of available candidates: \n     "
          f"prompt: {n_prompt}\n     nonprompt: {n_nonprompt}\n     bkg: {n_bkg}\n")

    n_cand_min = min([n_prompt, n_nonprompt, n_bkg])
    share = config["data_prep"]["class_balance"]["share"]
    if share == "equal":
        n_bkg = n_prompt = n_nonprompt = n_cand_min
    elif share == "all_signal":
        n_bkg = min(
            [n_bkg, (n_prompt + n_nonprompt) * config["data_prep"]["class_balance"]["bkg_factor"]])
    else:
        print(f"ERROR: class_balance option {share} not implemented")
        sys.exit()

    print("\nNumber of candidates used for training and test: \n     "
          f"prompt: {n_prompt}\n     nonprompt: {n_nonprompt}\n     bkg: {n_bkg}\n")

    df_tot = pd.concat(
        [df_bkg[:n_bkg],
         df_prompt[:n_prompt],
         df_nonprompt[:n_nonprompt]],
        sort=True
    )

    labels_array = np.array([0]*n_bkg + [1]*n_prompt + [2]*n_nonprompt)
    if 0 < test_f < 1:
        train_set, test_set, y_train, y_test = train_test_split(
            df_tot, labels_array, test_size=test_f, random_state=seed_split
        )
    else:
        print("ERROR: test_fraction must belong to ]0,1[")
        sys.exit(0)

    train_test_data = [train_set, y_train, test_set, y_test]
    del df_tot

    df_list = [df_bkg, df_prompt, df_nonprompt]
    leg_labels = ["bkg", "prompt", "nonprompt"]

    # _____________________________________________
    plot_utils.plot_distr(df_list, training_vars, 100, leg_labels,
                          figsize=(12, 7), alpha=0.3, log=True, grid=False, density=True)
    plt.subplots_adjust(left=0.06, bottom=0.06, right=0.99,
                        top=0.96, hspace=0.55, wspace=0.55)
    plt.savefig(f"{out_dir}/DistributionsAll_{channel}.pdf")
    plt.savefig(f"{out_dir}/DistributionsAll_{channel}.svg")
    plt.close("all")

    # _____________________________________________
    corr_matrix_fig = plot_utils.plot_corr(df_list, training_vars, leg_labels)
    for fig, lab in zip(corr_matrix_fig, leg_labels):
        plt.figure(fig.number)
        plt.subplots_adjust(left=0.2, bottom=0.25, right=0.95, top=0.9)
        fig.savefig(f"{out_dir}/CorrMatrix_{channel}_{lab}.pdf")
        fig.savefig(f"{out_dir}/CorrMatrix_{channel}_{lab}.svg")

    return train_test_data


def train(config, train_test_data):  # pylint: disable=too-many-locals
    """
    Function for the training

    Parameters
    -----------------
    - config: dictionary with config read from a yaml file
    - train_test_data: list with training and test data
    """

    out_dir = config["output"]["directory"]
    channel = config["data_prep"]["channel"]
    n_classes = len(np.unique(train_test_data[3]))
    model_clf = xgb.XGBClassifier(use_label_encoder=False)
    training_vars = config["ml"]["training_vars"]
    hyper_pars = config["ml"]["hyper_pars"]
    model_hdl = ModelHandler(model_clf, training_vars, hyper_pars)

    # hyperparameters optimization
    if config["ml"]["hyper_pars_opt"]["activate"]:
        model_hdl.optimize_params_optuna(
            train_test_data,
            config["ml"]["hyper_pars_opt"]["hyper_par_ranges"],
            cross_val_scoring="roc_auc_ovo",
            timeout=config['ml']['hyper_pars_opt']['timeout'],
            n_jobs=config['ml']['hyper_pars_opt']['njobs'],
            n_trials=config['ml']['hyper_pars_opt']['ntrials'],
            direction='maximize'
        )
    else:
        model_hdl.set_model_params(hyper_pars)

    # train and test the model with the updated hyper-parameters
    y_pred_test = model_hdl.train_test_model(
        train_test_data,
        True,
        output_margin=config["ml"]["raw_output"],
        average=config["ml"]["roc_auc_average"],
        multi_class_opt=config["ml"]["roc_auc_approach"]
    )

    output_labels = [config["output"]["out_labels"]["Bkg"], config["output"]
                     ["out_labels"]["Prompt"], config["output"]["out_labels"]["Nonprompt"]]
    df_column_to_save_list = config["output"]["column_to_save_list"]

    test_set_df = train_test_data[2]
    test_set_df = test_set_df.loc[:, df_column_to_save_list]
    test_set_df[f'Labels'] = train_test_data[3]

    for pred, lab in enumerate(output_labels):
        test_set_df[f'ML_output_{lab}'] = y_pred_test[:, pred]
    test_set_df.to_parquet(f"{out_dir}/{channel}_ModelApplied.parquet.gzip")

    # save model
    if os.path.isfile(f"{out_dir}/ModelHandler_{channel}.pickle"):
        os.remove(f"{out_dir}/ModelHandler_{channel}.pickle")
    if os.path.isfile(f"{out_dir}/ModelHandler_onnx_{channel}.onnx"):
        os.remove(f"{out_dir}/ModelHandler_onnx_{channel}.onnx")
    if os.path.isfile(f"{out_dir}/ModelHandler_onnx_hummingbird_{channel}"):
        os.remove(f"{out_dir}/ModelHandler_onnx_hummingbird_{channel}")

    model_hdl.dump_model_handler(f"{out_dir}/ModelHandler_{channel}.pickle")
    model_conv = H4MLConverter(model_hdl)
    model_conv.convert_model_onnx(1)
    model_conv.dump_model_onnx(f"{out_dir}/ModelHandler_onnx_{channel}.onnx")
    model_conv.convert_model_hummingbird("onnx", 1)
    model_conv.dump_model_hummingbird(
        f"{out_dir}/ModelHandler_onnx_hummingbird_{channel}")

    # plots
    leg_labels = ["bkg", "prompt", "nonprompt"]

    # _____________________________________________
    plt.rcParams["figure.figsize"] = (10, 7)
    fig_ml_output = plot_utils.plot_output_train_test(
        model_hdl,
        train_test_data,
        80,
        config['ml']['raw_output'],
        leg_labels,
        True,
        density=True
    )

    if n_classes > 2:
        for fig, lab in zip(fig_ml_output, leg_labels):
            fig.savefig(f'{out_dir}/MLOutputDistr_{lab}_{channel}.pdf')
            fig.savefig(f'{out_dir}/MLOutputDistr_{lab}_{channel}.svg')
    else:
        fig_ml_output.savefig(f'{out_dir}/MLOutputDistr_{channel}.pdf')
        fig_ml_output.savefig(f'{out_dir}/MLOutputDistr_{channel}.svg')

    # _____________________________________________
    plt.rcParams["figure.figsize"] = (10, 9)
    fig_roc_curve = plot_utils.plot_roc(
        train_test_data[3],
        y_pred_test,
        None,
        leg_labels,
        config['ml']['roc_auc_average'],
        config['ml']['roc_auc_approach']
    )
    fig_roc_curve.savefig(f'{out_dir}/ROCCurveAll_{channel}.pdf')
    fig_roc_curve.savefig(f'{out_dir}/ROCCurveAll_{channel}.svg')
    pickle.dump(fig_roc_curve, open(
        f'{out_dir}/ROCCurveAll_{channel}.pkl', 'wb'))

    # _____________________________________________
    plt.rcParams["figure.figsize"] = (12, 7)
    fig_feat_importance = plot_utils.plot_feature_imp(
        train_test_data[2][train_test_data[0].columns],
        train_test_data[3],
        model_hdl,
        leg_labels
    )
    n_plot = n_classes if n_classes > 2 else 1
    for i_fig, fig in enumerate(fig_feat_importance):
        if i_fig < n_plot:
            lab = leg_labels[i_fig] if n_classes > 2 else ''
            fig.savefig(f'{out_dir}/FeatureImportance_{lab}_{channel}.pdf')
            fig.savefig(f'{out_dir}/FeatureImportance_{lab}_{channel}.svg')
        else:
            fig.savefig(f'{out_dir}/FeatureImportanceAll_{channel}.pdf')
            fig.savefig(f'{out_dir}/FeatureImportanceAll_{channel}.svg')


def main(config):
    """
    Main function

    Parameters
    -----------------
    - config: dictionary with config read from a yaml file
    """
    train_test_data = data_prep(config)
    train(config, train_test_data)

    os._exit(0)  # pylint: disable=protected-access


if __name__ == "__main__":
    PARSER = argparse.ArgumentParser(description="Arguments")
    PARSER.add_argument("config", metavar="text", default="config_training.yml",
                        help="config file for training")
    ARGS = PARSER.parse_args()

    with open(ARGS.config, "r") as yml_cfg:  # pylint: disable=bad-option-value
        CFG = yaml.load(yml_cfg, yaml.FullLoader)

    main(CFG)
