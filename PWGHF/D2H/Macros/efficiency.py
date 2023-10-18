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
file: efficiency.py
brief: script for (acceptance-times-)efficiency computation
usage: python3 efficiency.py CONFIG
author: Alexandre Bigot <alexandre.bigot@cern.ch>, Strasbourg University
"""

import argparse
import os
import sys
from typing import Dict, List, Optional, Union

import numpy as np  # pylint: disable=import-error
import yaml  # pylint: disable=import-error

# pylint: disable=import-error, no-name-in-module
from ROOT import (
    TH1,
    TH1F,
    TH2,
    TCanvas,
    TEfficiency,
    TFile,
    TLatex,
    TLegend,
    kAzure,
    kBlack,
    kFullCircle,
    kFullSquare,
    kOpenSquare,
    kRed,
)
from style_formatter import set_global_style, set_object_style

# options for object style
COLOR = {"All": kBlack, "Prompt": kRed + 1, "Nonprompt": kAzure + 4}
MARKERSTYLE = {"All": kFullSquare, "Prompt": kFullCircle, "Nonprompt": kOpenSquare}
LINEWIDTH = {"All": 1, "Prompt": 1, "Nonprompt": 2}

LINESTYLE = {"All": 1, "Prompt": 2, "Nonprompt": 7}


def enforce_list(x: Union[str, List[str]]) -> List[str]:
    """
    Helper method to enforce list type

    Parameters
    ----------
    - x: a string or a list of string

    Returns
    ----------
    - x_list if x was not a list, x itself otherwise
    """

    if not isinstance(x, list):
        # handle possible whitespaces in config file entry
        x_list = x.split(",")
        for i, element in enumerate(x_list):
            x_list[i] = element.strip()  # remove possible whitespaces
        return x_list

    return x


def enforce_trailing_slash(path: str) -> str:
    """
    Helper method to enforce '/' at the and of directory name

    Parameters
    ----------
    - path: some path

    Returns
    ----------
    - path with a trailing slash at the end if it was not there yet
    """

    if path is not None and path[-1] != "/":
        path += "/"

    return path


# pylint: disable=too-many-arguments
def configure_canvas(
    name_canvas: str,
    pt_min: Union[int, float],
    y_axis_min: Union[int, float],
    pt_max: Union[int, float],
    y_axis_max: Union[int, float],
    title: str,
    log_y_axis: bool,
) -> TCanvas:
    """
    Helper method to configure canvas

    Parameters
    ----------
    - name_canvas: name of the canvas
    - pt_min: lower limit of x axis
    - y_axis_min: lower limit of y axis
    - pt_max: upper limit of x axis
    - y_axis_max: upper limit of y axis
    - title: title of the canvas
    - log_y_axis: switch for log scale along y axis

    Returns
    ----------
    - c_eff: TCanvas instance
    """

    c_eff = TCanvas(name_canvas, "", 800, 800)
    c_eff.DrawFrame(pt_min, y_axis_min, pt_max, y_axis_max, title)
    if log_y_axis:
        c_eff.SetLogy()

    return c_eff


def save_canvas(canvas: TCanvas, out_dir: str, name_file: str, extension: List[str]) -> None:
    """
    Save canvas in formats chosen by extension

    Parameters
    ----------
    - canvas: a TCanvas instance
    - out_dir: output directory where canvas will be saved
    - name_file: name of the output file
    - extension: file format
    """

    for ext in extension:
        canvas.SaveAs(out_dir + name_file + "." + ext)


def __set_object_style(obj: Union[TEfficiency, TH1], key: str) -> None:
    """
    Helper method to set style of TEfficiency or TH1 object with local options

    Parameters
    ----------
    - obj: TEfficiency or TH1 object
    - key: key of options dictionary
    """

    set_object_style(
        obj, color=COLOR[key], markerstyle=MARKERSTYLE[key], linewidth=LINEWIDTH[key], linestyle=LINESTYLE[key]
    )


def compute_efficiency(
    h_rec: Union[TH1, TH2],
    h_gen: Union[TH1, TH2],
    axis_rapidity: str = "Y",
    rapidity_cut: Optional[float] = None,
    pt_bins_limits: Optional[np.ndarray] = None,
) -> TEfficiency:
    """
    Helper method to compute the efficiency as function of the feature in axis.

    Parameters
    ----------
    - h_rec: histogram with the reconstructed information
    - h_gen: histogram with the generated information
    - rapidity_cut: applies the selection |y| < y_cut to the efficiency
    - axis_rapidity: histogram axis containting rapidity

    Returns
    ----------
    - efficiency: a TEfficiency instance
    """

    epsilon = 0.0001

    is_reco_2d, is_gen_2d = False, False
    if h_rec.GetDimension() == 2:
        is_reco_2d = True
    if h_gen.GetDimension() == 2:
        is_gen_2d = True

    if is_reco_2d:
        if axis_rapidity == "Y":
            h_rec = h_rec.ProjectionX()
        else:
            h_rec = h_rec.ProjectionY()

    if is_gen_2d:
        if axis_rapidity == "Y":
            if rapidity_cut is not None:
                h_gen.GetYaxis().SetRangeUser(-1.0 * rapidity_cut + epsilon, rapidity_cut - epsilon)
            h_gen = h_gen.ProjectionX()
        else:
            if rapidity_cut is not None:
                h_gen.GetXaxis().SetRangeUser(-1.0 * rapidity_cut + epsilon, rapidity_cut - epsilon)
            h_gen = h_gen.ProjectionY()

    # rebin histograms, if enabled
    if pt_bins_limits is not None:
        n_bins = pt_bins_limits.shape[0] - 1
        h_rec = h_rec.Rebin(n_bins, "hRec", pt_bins_limits)
        h_gen = h_gen.Rebin(n_bins, "hGen", pt_bins_limits)

    efficiency = TEfficiency(h_rec, h_gen)

    return efficiency


def get_th1_from_tefficiency(teff: TEfficiency, h_eff: TH1) -> TH1:
    """
    Helper method to convert TEfficiency object to TH1

    Parameters
    ----------
    - teff: TEfficiency instance
    - h_eff: clone of histogram with reconstructed information

    Returns
    ----------
    - h_eff: histogram containing teff information
    """

    n_bins = h_eff.GetXaxis().GetNbins()
    for i_bin in range(1, n_bins + 1):
        h_eff.SetBinContent(i_bin, teff.GetEfficiency(i_bin))
        if (
            abs(teff.GetEfficiencyErrorLow(i_bin) - teff.GetEfficiencyErrorUp(i_bin))
            / teff.GetEfficiencyErrorLow(i_bin)
            > 1e-02
        ):
            print(
                f"\033[93mWARNING: efficiency error is asymmetric in bin {i_bin},"
                " setting the maximum error in the histogram!\033[0m"
            )
            # TH1 can't handle asymmetric errors so we take the max
            err = max(teff.GetEfficiencyErrorLow(i_bin), teff.GetEfficiencyErrorUp(i_bin))
        else:
            err = teff.GetEfficiencyErrorLow(i_bin)

        h_eff.SetBinError(i_bin, err)

    return h_eff


# pylint: disable=too-many-locals, too-many-branches, too-many-statements
def main(cfg: Dict, batch: bool) -> None:
    """
    Main function

    Parameters
    -----------------
    - config: dictionary with config read from a yaml file
    - batch: bool to suppress video output
    """

    set_global_style(padleftmargin=0.14, padbottommargin=0.12, titlesize=0.045, labelsize=0.04)

    # input configuration
    name_file = cfg["input"]["filename"]
    name_tree = cfg["input"]["treename"]
    file = TFile.Open(name_file)
    name_tree = enforce_trailing_slash(name_tree) if name_tree is not None else ""

    # fill dictionaries with histograms
    dic_rec: Dict[str, Optional[Union[TH1, TH2]]] = {}
    dic_gen: Dict[str, Optional[Union[TH1, TH2]]] = {}
    for key, histoname in cfg["input"]["histoname"]["reconstructed"].items():
        if histoname is None:
            dic_rec[key] = None
        else:
            dic_rec[key] = file.Get(name_tree + histoname)
    for key, histoname in cfg["input"]["histoname"]["generated"].items():
        if histoname is None:
            dic_gen[key] = None
        else:
            dic_gen[key] = file.Get(name_tree + histoname)

    # configure pt binning
    pt_bins_limits: Optional[np.ndarray] = None
    pt_min, pt_max = 0.0, 36.0  # default values
    is_retrieved_pt_interval = False
    if cfg["pt_bins_limits"] is not None:
        pt_bins_limits = np.array(enforce_list(cfg["pt_bins_limits"]), "d")

    if pt_bins_limits is not None:
        pt_min = pt_bins_limits[0]
        pt_max = pt_bins_limits[-1]
        is_retrieved_pt_interval = True

    #  configure possible cut on rapidity
    rapidity_cut = cfg["rapidity"]["cut"]
    axis_rapidity = cfg["rapidity"]["axis"]

    # output configuration
    out_dir = enforce_trailing_slash(cfg["output"]["dir"])
    out_label = cfg["output"]["plots"]["label"]
    name_axis = cfg["output"]["plots"]["y_axis"]["name"]
    y_axis_min = cfg["output"]["plots"]["y_axis"]["min"]
    y_axis_max = cfg["output"]["plots"]["y_axis"]["max"]
    log_y_axis = cfg["output"]["plots"]["y_axis"]["log_scale"]
    if name_axis is None:
        name_axis = "#varepsilon"  # default value
    title = ";#it{p}_{T} (GeV/#it{c});" + name_axis + ";"

    overlap: Union[str, List[str]] = str()
    if cfg["output"]["plots"]["overlap"] is not None:
        overlap = enforce_list(cfg["output"]["plots"]["overlap"])

    # output save options
    save_tefficiency = cfg["output"]["save"]["TEfficiency"]
    save_th1 = cfg["output"]["save"]["TH1"]
    save_tcanvas_individual = cfg["output"]["save"]["TCanvas"]["individual"]
    save_tcanvas_overlap = cfg["output"]["save"]["TCanvas"]["overlap"]
    extension: List[str] = ["pdf"]  # default value
    if cfg["output"]["save"]["TCanvas"]["extension"] is None:
        if save_tcanvas_individual or save_tcanvas_overlap:
            print(
                "\033[93mWARNING: No extension provided for saving canvas in extra file,"
                " '.pdf' set as default.\033[0m"
            )
    else:
        extension = enforce_list(cfg["output"]["save"]["TCanvas"]["extension"])

    if os.path.isdir(out_dir):
        print((f"\033[93mWARNING: Output directory '{out_dir}' already exists," " overwrites possibly ongoing!\033[0m"))
    else:
        os.makedirs(out_dir)

    # plots watermark configuration
    latex = TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.035)
    latex.SetTextAlign(13)  # align at top
    latex.SetTextFont(42)
    watermark = cfg["output"]["plots"]["watermark"]

    # open output file
    out_file = TFile(out_dir + out_label + ".root", "recreate")

    # compute efficiencies
    efficiencies = {}
    for (key, h_rec), (_, h_gen) in zip(dic_rec.items(), dic_gen.items()):
        if h_rec is not None and h_gen is not None:
            # compute efficiency, store it in TEfficiency instance and configure it
            efficiencies[key] = compute_efficiency(h_rec, h_gen, axis_rapidity, rapidity_cut, pt_bins_limits)
            efficiencies[key].SetName("TEfficiency_" + out_label + key)
            efficiencies[key].SetTitle(title)
            __set_object_style(efficiencies[key], key)

            # save TEfficiency instance in output file, if enabled
            if save_tefficiency:
                efficiencies[key].Write()

            # retrieve pt_min and pt_max from imported histogram
            if not is_retrieved_pt_interval:
                pt_min, pt_max = h_rec.GetXaxis().GetXmin(), h_rec.GetXaxis().GetXmax()
                is_retrieved_pt_interval = True

            # plot efficiency on canvas
            c_eff = configure_canvas(f"c{out_label + key}", pt_min, y_axis_min, pt_max, y_axis_max, title, log_y_axis)
            efficiencies[key].Draw("same")  # to draw TEfficiency in TCanvas
            latex.DrawLatex(0.18, 0.92, watermark)
            c_eff.Write()
            # save canvas in separate file, if enabled
            if save_tcanvas_individual:
                save_canvas(c_eff, out_dir, out_label + key, extension)

            # convert TEfficiency instance to TH1 histogram in output file, if save enabled
            if save_th1:
                if pt_bins_limits is not None:
                    h_eff = TH1F("h" + out_label + key, title, len(pt_bins_limits) - 1, pt_bins_limits)
                    h_eff = get_th1_from_tefficiency(efficiencies[key], h_eff)
                else:
                    h_eff = get_th1_from_tefficiency(efficiencies[key], h_rec.Clone("h" + out_label + key))
                    h_eff.SetTitle(title)
                __set_object_style(h_eff, key)
                h_eff.Write()

        elif h_rec != h_gen:  # one histogram name is None, the other is not
            print(
                f"\033[91mERROR: efficiency for {key} could not be computed,"
                " one of the histogram names is None!\033[0m"
            )
            sys.exit()
        else:  # both histogram names are None
            efficiencies[key] = None
            print(f"\033[94mEfficiency for {key} not computed.\033[0m")

    # overlap plots, if enabled
    if overlap:
        c_overlap = configure_canvas(
            f"cOverlap{''.join(overlap)}", pt_min, y_axis_min, pt_max, y_axis_max, title, log_y_axis
        )

        leg = TLegend(0.6, 0.2, 0.8, 0.4)
        leg.SetTextSize(0.045)
        leg.SetFillStyle(0)

        for key, teff in efficiencies.items():
            if key in overlap:
                if teff is None:
                    print(f"\033[91mERROR: efficiency for {key} could not be computed," " cannot overlap it!\033[0m")
                    sys.exit()

                leg.AddEntry(teff, key, "p")
                teff.Draw("same")

        leg.Draw()
        latex.DrawLatex(0.18, 0.92, watermark)

        if save_tcanvas_overlap:
            save_canvas(c_overlap, out_dir, out_label + "Overlap", extension)
        c_overlap.Write()

    # close output file
    out_file.Close()

    if not batch:
        input("Press enter to exit")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument("config", metavar="text", default="config.yml", help="config file name for ml")
    parser.add_argument("--batch", help="suppress video output", action="store_true")
    args = parser.parse_args()

    print("Loading analysis configuration: ...", end="\r")
    with open(args.config, "r", encoding="utf-8") as yml_cfg:
        configuration = yaml.load(yml_cfg, yaml.FullLoader)
    print("Loading analysis configuration: Done!")

    main(configuration, args.batch)
