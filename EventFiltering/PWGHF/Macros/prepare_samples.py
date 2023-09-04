"""
Script for the preparation of the samples for the
training of ML models to be used in the HF triggers

\author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN
\author Alexandre Bigot <alexandre.bigot@cern.ch>, Strasbourg University
\author Biao Zhang <biao.zhang@cern.ch>, CCNU
"""

import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
import uproot
from alive_progress import alive_bar
from ROOT import TFile, gRandom  # pylint: disable=no-name-in-module

BITS3P = {"DplusToPiKPi": 0,
          "LcToPKPi": 1,
          "DsToKKPi": 2,
          "XicToPKPi": 3}

CHANNELS3P = {"DplusToPiKPi": 1,
              "DsToKKPi": 2,
              "LcToPKPi": 3,
              "XicToPKPi": 4}


# pylint: disable=too-many-locals
def do_dca_smearing(input_df, n_prongs=None):
    """
    Method to do the DCA smearing for 2 prongs and 3 prongs

    Parameters
    -----------------
    - input_df: pandas dataframe containing all candidates with fFlagOrigin column
    - n_prongs: option to 2 prongs or 3 prongs

    Outputs
    -----------------
    - input_df: New dataframe with the smeared DCA columns
    """

    print(f"Start to do the smearing for {n_prongs} prong")
    # Open the input files
    file_names_data = ["New_DCA_smear/results_dcaXY_LHC22m_pass4_run523308.root",
                       "New_DCA_smear/results_dcaZ_LHC22m_pass4_run523308.root"]
    file_names_mc = ["New_DCA_smear/output_DCAxy_all.root",
                     "New_DCA_smear/output_DCAz_all.root"]
    input_files_data = [TFile.Open(name) for name in file_names_data]
    input_files_mc = [TFile.Open(name) for name in file_names_mc]
    input_files_mean = [TFile.Open(name) for name in file_names_data]

    # Extract DCA resolution histograms
    dca_reso_data, dca_reso_mc, dca_reso_data_meanshift = {}, {}, {}
    for i, par in enumerate(["XY", "Z"]):
        dca_reso_data[par] = input_files_data[i].Get("tge_DCA_res_withoutPVrefit_all")
        dca_reso_mc[par] = input_files_mc[i].Get("tge_DCA_res_withoutPVrefit_all")
        dca_reso_data_meanshift[par] = input_files_mean[i].Get("tge_DCA_res_withoutPVrefit_all")

    # Add smeared DCA columns to the dataframe
    smear_cols = ["XY1", "XY2", "Z1", "Z2"]
    if n_prongs == 3:
        smear_cols.extend(["XY3", "Z3"])
    for col in smear_cols:
        dca_col = f"fDCAPrim{col}"
        pt_col = f"fPT{col[-1]}"

        smear_values = []
        for dca, pt in zip(input_df[dca_col], input_df[pt_col]):
            if pt < 7:
                dca_reso_data_mean_shift = dca_reso_data_meanshift[col[:-1]].Eval(pt) * 1e-4
                dca_reso_data_val = dca_reso_data[col[:-1]].Eval(pt)
                dca_reso_mc_val = dca_reso_mc[col[:-1]].Eval(pt)
            else:
                dca_reso_data_mean_shift = dca_reso_data_meanshift[col[:-1]].Eval(7) * 1e-4
                dca_reso_data_val = dca_reso_data[col[:-1]].Eval(7)
                dca_reso_mc_val = dca_reso_mc[col[:-1]].Eval(7)

            smear_value = gRandom.Gaus(dca + dca_reso_data_mean_shift,
                                       np.sqrt(dca_reso_data_val**2 - dca_reso_mc_val**2) * 1e-4)
            smear_values.append(smear_value)

        input_df[f"{dca_col}_SMEAR"] = smear_values

    # Close the input files
    for file_data, file_mc, file_mean in zip(input_files_data, input_files_mc, input_files_mean):
        file_data.Close()
        file_mc.Close()
        file_mean.Close()

    # Make a figure comparing the DCA variables before and after smearing
    num_cols = len(smear_cols)
    num_rows = num_cols // 2 + num_cols % 2
    fig, axs = plt.subplots(num_rows, 2, figsize=(10, 8))
    axs = axs.flatten()

    for i, col in enumerate(smear_cols):
        dca_col = f"fDCAPrim{col}"
        smear_col = f"{dca_col}_SMEAR"
        axs[i].hist(input_df[dca_col], bins=500, alpha=0.5, label="Before Smearing")
        axs[i].hist(input_df[smear_col], bins=500, alpha=0.3, label="After Smearing")
        axs[i].set_xlabel(col)
        axs[i].set_xlim(-0.05, 0.05)
        axs[i].set_ylim(10e2, None)
        axs[i].set_yscale('log')
        axs[i].legend()

    plt.tight_layout()
    fig.savefig(f"dca_comparison_{n_prongs}prong.png")
    plt.close(fig)

    return input_df


def divide_df_for_origin(input_df, cols_to_remove=None, channel=None):
    """
    Method to divide a dataframe in three (prompt, non-prompt, bkg)

    Parameters
    -----------------
    - input_df: pandas dataframe containing all candidates with fFlagOrigin column
    - cols_to_remove: columns to be removed from output dataframes
    - channel: integer corresponding a specific fChannel for signal

    Outputs
    -----------------
    - df_prompt: pandas dataframe containing only prompt signal
    - df_nonprompt: pandas dataframe containing only non-prompt signal
    - df_bkg: pandas dataframe containing only background candidates
    """

    if cols_to_remove is None:
        cols_to_remove = ['fFlagOrigin']

    df_prompt = input_df.query("fFlagOrigin == 1")
    df_nonprompt = input_df.query("fFlagOrigin == 2")
    if channel is not None and "fChannel" in df_prompt.columns:
        df_prompt = df_prompt.query(f"fChannel == {channel}")
        df_nonprompt = df_nonprompt.query(f"fChannel == {channel}")
    df_bkg = input_df.query("fFlagOrigin == 0")
    cols_to_keep = list(df_prompt.columns)
    for col in cols_to_remove:
        cols_to_keep.remove(col)
    df_prompt = df_prompt[cols_to_keep]
    df_nonprompt = df_nonprompt[cols_to_keep]
    df_bkg = df_bkg[cols_to_keep]

    return df_prompt, df_nonprompt, df_bkg


# pylint: disable=too-many-locals,too-many-branches
# pylint: disable=too-many-nested-blocks,too-many-statements
def main(input_dir, max_files=1000, downscale_bkg=1., force=False, do_smearing=False):
    """
    Main function

    Parameters
    -----------------
    - config: dictionary with configs
    """

    input_files = []
    for subdir in os.listdir(input_dir):
        if os.path.isdir(os.path.join(input_dir, subdir)):
            for subsubdir in os.listdir(os.path.join(input_dir, subdir)):
                if os.path.isdir(os.path.join(input_dir, subdir, subsubdir)):
                    for file in os.listdir(os.path.join(input_dir, subdir, subsubdir)):
                        if "AO2D.root" in file:
                            input_files.append(os.path.join(
                                input_dir, subdir, subsubdir, file))
                        elif os.path.isdir(os.path.join(input_dir, subdir, subsubdir, file)):
                            for file2 in os.listdir(os.path.join(
                                    input_dir, subdir, subsubdir, file)):
                                if "AO2D.root" in file2:
                                    input_files.append(os.path.join(
                                        input_dir, subdir, subsubdir, file, file2))

    df_2p, df_3p = None, None
    with alive_bar(len(input_files[:max_files])) as bar_alive:
        for file in input_files[:max_files]:
            print(f"\033[32mExtracting dataframes from input "
                  f"{file}\033[0m")

            file_root = uproot.open(file)
            indir = os.path.split(file)[0]

            # 2-prongs --> only D0
            is_d0_filtered = False
            for exfile in os.listdir(indir):
                if "D0ToKPi.parquet.gzip" in exfile:
                    is_d0_filtered = True
                    break
            if not is_d0_filtered or force:
                list_of_2p_df = []
                for tree_name in file_root.keys():
                    if "O2hftrigtrain2p" in tree_name:
                        list_of_2p_df.append(f"{file}:{tree_name}")
                df_2p = uproot.concatenate(list_of_2p_df, library="pd")
                if do_smearing:
                    df_2p = do_dca_smearing(df_2p, 2)

                df_2p_prompt, df_2p_nonprompt, df_2p_bkg = divide_df_for_origin(df_2p)
                df_2p_bkg = df_2p_bkg.sample(frac=downscale_bkg, random_state=42)
                df_2p_prompt.to_parquet(
                    os.path.join(indir, "Prompt_D0ToKPi.parquet.gzip"),
                    compression="gzip"
                )
                df_2p_nonprompt.to_parquet(
                    os.path.join(indir, "Nonprompt_D0ToKPi.parquet.gzip"),
                    compression="gzip"
                )
                df_2p_bkg.to_parquet(
                    os.path.join(indir, "Bkg_D0ToKPi.parquet.gzip"),
                    compression="gzip"
                )
                df_2p = None

            # 3-prongs --> D+, Ds+, Lc+, Xic+
            is_3p_filtered = False
            for channel_3p in BITS3P:
                for exfile in os.listdir(indir):
                    if f"{channel_3p}.parquet.gzip" in exfile:
                        is_3p_filtered = True
                        break

            list_of_3p_df = []
            if not is_3p_filtered or force:
                for tree_name in file_root.keys():
                    if "O2hftrigtrain3p" in tree_name:
                        list_of_3p_df.append(f"{file}:{tree_name}")
                df_3p = uproot.concatenate(list_of_3p_df, library="pd")
                if do_smearing:
                    df_3p = do_dca_smearing(df_3p, 3)

                for channel_3p in BITS3P:
                    flags = df_3p["fHFSelBit"].astype(
                        int) & 2**BITS3P[channel_3p]
                    df_channel_3p = df_3p[flags.astype("bool").to_numpy()]
                    df_3p_prompt, df_3p_nonprompt, df_3p_bkg = divide_df_for_origin(
                        df_channel_3p,
                        ["fFlagOrigin", "fChannel", "fHFSelBit"],
                        channel=CHANNELS3P[channel_3p]
                    )
                    df_3p_bkg = df_3p_bkg.sample(
                        frac=downscale_bkg, random_state=42)
                    df_3p_prompt.to_parquet(
                        os.path.join(
                            indir, f"Prompt_{channel_3p}.parquet.gzip"),
                        compression="gzip"
                    )
                    df_3p_nonprompt.to_parquet(
                        os.path.join(
                            indir, f"Nonprompt_{channel_3p}.parquet.gzip"),
                        compression="gzip"
                    )
                    df_3p_bkg.to_parquet(
                        os.path.join(indir, f"Bkg_{channel_3p}.parquet.gzip"),
                        compression="gzip"
                    )

                df_3p = None

            bar_alive()


if __name__ == "__main__":
    PARSER = argparse.ArgumentParser(description="Arguments")
    PARSER.add_argument("input_dir", metavar="text", default="AO2D/trains",
                        help="input directory with AO2D input files")
    PARSER.add_argument("--max_files", type=int, default=1000,
                        help="max input files to be processed")
    PARSER.add_argument("--downscale_bkg", type=float, default=1.,
                        help="fraction of bkg to be kept")
    PARSER.add_argument("--force", action="store_true", default=False,
                        help="force re-creation of output files")
    PARSER.add_argument("--dosmearing", action="store_true", default=False,
                        help="do smearing on the dca of daughter tracks ")
    ARGS = PARSER.parse_args()

    main(ARGS.input_dir, ARGS.max_files, ARGS.downscale_bkg, ARGS.force, ARGS.dosmearing)
