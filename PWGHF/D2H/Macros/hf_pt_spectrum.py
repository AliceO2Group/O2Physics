#!/usr/bin/env python3

"""
file: hf_pt_spectrum.py
brief: script for computation of pT-differential yields (cross sections)
usage: python3 HfPtSpectrum.py CONFIG
author: Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN
"""

import argparse
import os
import sys

import numpy as np  # pylint: disable=import-error
import yaml  # pylint: disable=import-error
from hf_analysis_utils import (
    compute_crosssection,
    compute_fraction_fc,
    compute_fraction_nb,
    get_hist_binlimits,
)
from ROOT import (  # pylint: disable=import-error,no-name-in-module
    TH1,
    TH1F,
    TCanvas,
    TFile,
    TGraphAsymmErrors,
    TStyle,
    gPad,
    gROOT,
    kAzure,
    kFullCircle,
)


# pylint: disable=too-many-locals,too-many-statements
def load_inputs(input_cfg):
    """
    Helper method to load inputs

    Parameters
    ----------
    - input_cfg: dictionary from yaml config file

    Returns
    ----------
    - histos: dictionary with input histos, keys (rawyields, acceffp, acceffnp)
    - norm: dictionary with normalisations, keys (BR, events, sigmaMB)
    """

    observable = input_cfg["observable"]
    if observable not in ["dsigmadpt", "dNdpt"]:
        print(f"\033[91mERROR: observable {observable} not supported. Exit\033[0m")
        sys.exit(1)

    channel = input_cfg["channel"]
    if channel not in [
        "D0toKpi",
        "DplustoKpipi",
        "DstoKpipi",
        "DstartoD0pi",
        "LctopKpi",
        "LctopK0S",
    ]:
        print(f"\033[91mERROR: channel {channel} not supported. Exit\033[0m")
        sys.exit(2)

    system = input_cfg["system"]
    if system not in ["pp", "pPb", "PbPb"]:
        print(f"\033[91mERROR: channel {channel} not supported. Exit\033[0m")
        sys.exit(3)
    if system in ["pPb", "PbPb"] and observable == "dsigmadpt":
        print("\033[93mWARNING: switching from dsigmadpt to dNdpt\033[0m")

    energy = input_cfg["energy"]
    if energy not in ["5TeV", "13TeV"]:
        print(f"\033[91mERROR: energy {energy} not supported. Exit\033[0m")
        sys.exit(4)

    frac_method = input_cfg["fraction"]
    if frac_method not in ["Nb", "fc"]:
        print(f"\033[91mERROR: method to subtract nonprompt" f" {frac_method} not supported. Exit\033[0m")
        sys.exit(5)

    rawy_file_name = input_cfg["rawyield"]["filename"]
    rawy_hist_name = input_cfg["rawyield"]["rawyieldhist"]
    norm_hist_name = input_cfg["rawyield"]["normhist"]

    eff_file_name = input_cfg["acceff"]["filename"]
    effprompt_hist_name = input_cfg["acceff"]["prompthist"]
    effnonprompt_hist_name = input_cfg["acceff"]["nonprompthist"]

    pred_file_name = input_cfg["FONLL"]

    # load histos from root files
    histos = {}
    infile_rawy = TFile.Open(rawy_file_name)
    histos["rawyields"] = infile_rawy.Get(rawy_hist_name)
    if not histos["rawyields"]:
        print(f"\033[91mERROR: raw-yield histo {rawy_hist_name}" f" not found in {rawy_file_name}. Exit\033[0m")
        sys.exit(6)
    histos["rawyields"].SetDirectory(0)
    h_events = infile_rawy.Get(norm_hist_name)
    if not h_events:
        print(f"\033[91mERROR: normalisation histo {norm_hist_name}" f" not found in {rawy_file_name}. Exit\033[0m")
        sys.exit(7)
    h_events.SetDirectory(0)
    infile_rawy.Close()

    infile_eff = TFile.Open(eff_file_name)
    histos["acceffp"] = infile_eff.Get(effprompt_hist_name)
    if not histos["acceffp"]:
        print(
            f"\033[91mERROR: prompt (acc x eff) histo {effprompt_hist_name}"
            f" not found in {eff_file_name}. Exit\033[0m"
        )
        sys.exit(8)
    histos["acceffp"].SetDirectory(0)
    histos["acceffnp"] = infile_eff.Get(effnonprompt_hist_name)
    if not histos["acceffnp"]:
        print(
            f"\033[91mERROR: nonprompt (acc x eff) histo {effprompt_hist_name}"
            f"not found in {eff_file_name}. Exit\033[0m"
        )
        sys.exit(9)
    histos["acceffnp"].SetDirectory(0)
    infile_eff.Close()

    fonll_hist_name = {
        "D0toKpi": "hD0Kpi",
        "DplustoKpipi": "hDpluskpipi",
        "DstoKKpi": "hDsPhipitoKkpi",
        "DstartoD0pi": "hDstarD0pi",
        "LctopKpi": "hLcpkpi",
        "LctopK0S": "hLcK0sp",
    }
    histos["FONLL"] = {"prompt": {}, "nonprompt": {}}
    infile_fonll = TFile.Open(pred_file_name)
    for pred in ("central", "min", "max"):
        histos["FONLL"]["nonprompt"][pred] = infile_fonll.Get(f"{fonll_hist_name[channel]}fromBpred_{pred}_corr")
        histos["FONLL"]["nonprompt"][pred].SetDirectory(0)
        if frac_method == "fc":
            histos["FONLL"]["prompt"][pred] = infile_fonll.Get(f"{fonll_hist_name[channel]}pred_{pred}")
            histos["FONLL"]["prompt"][pred].SetDirectory(0)
    infile_fonll.Close()

    # load normalisation info from common database
    norm = {}
    with open("config/norm_database.yml", "r") as yml_norm_db:
        norm_db = yaml.safe_load(yml_norm_db)
    norm["BR"] = norm_db["BR"][channel]["value"]
    norm["events"] = h_events.GetBinContent(1)
    norm["sigmaMB"] = norm_db["sigma"]["Run2"][system][energy] if observable == "dsigmadpt" else 1.0

    return histos, norm


def main():
    """
    Main function
    """

    parser = argparse.ArgumentParser(description="Arguments to pass")
    parser.add_argument(
        "configfile_name",
        metavar="text",
        default="config_Dplus_pp5TeV.yml",
        help="input yaml config file",
    )
    parser.add_argument("--batch", action="store_true", default=False, help="suppress video output")
    args = parser.parse_args()

    if args.batch:
        gROOT.SetBatch(True)

    # final plots style settings
    style_hist = TStyle("style_hist", "Histo graphics style")
    style_hist.SetOptStat("n")
    style_hist.SetMarkerColor(kAzure + 4)
    style_hist.SetMarkerStyle(kFullCircle)
    style_hist.SetMarkerSize(1)
    style_hist.SetHistLineColor(kAzure + 4)
    style_hist.SetHistLineWidth(2)
    style_hist.SetLabelSize(0.030)
    style_hist.SetLabelOffset(0.010)
    style_hist.SetTitleXOffset(1.3)
    style_hist.SetTitleYOffset(1.3)
    style_hist.SetDrawOption("AP")
    gROOT.SetStyle("style_hist")
    gROOT.ForceStyle()

    # load info from config file
    with open(args.configfile_name, "r") as yml_configfile:
        cfg = yaml.safe_load(yml_configfile)
    frac_method = cfg["fraction"]

    histos, norm = load_inputs(cfg)

    # consistency check of bins
    ptlims = {}
    for histo in ["rawyields", "acceffp", "acceffnp"]:
        ptlims[histo] = get_hist_binlimits(histos[histo])
        if histo != "rawyields" and not np.equal(ptlims[histo], ptlims["rawyields"]).all():
            print("\033[91mERROR: histo binning not consistent. Exit\033[0m")
            sys.exit(10)

    # compute cross section
    axistit_cross = "d#sigma/d#it{p}_{T} (pb GeV^{-1} #it{c})"
    axistit_cross_times_br = "d#sigma/d#it{p}_{T} #times BR (pb GeV^{-1} #it{c})"
    axistit_pt = "#it{p}_{T} (GeV/#it{c})"
    axistit_fprompt = "#if{f}_{prompt}"
    gfraction = TGraphAsymmErrors()
    gfraction.SetNameTitle("gfraction", f";{axistit_pt};{axistit_fprompt}")

    hptspectrum = TH1F(
        "hptspectrum",
        f";{axistit_pt};{axistit_cross}",
        len(ptlims["rawyields"]) - 1,
        ptlims["rawyields"],
    )
    hptspectrum_wo_br = TH1F(
        "hptspectrum_wo_br",
        f";{axistit_pt};{axistit_cross_times_br}",
        len(ptlims["rawyields"]) - 1,
        ptlims["rawyields"],
    )

    for i_pt, (ptmin, ptmax) in enumerate(zip(ptlims["rawyields"][:-1], ptlims["rawyields"][1:])):
        pt_cent = (ptmax + ptmin) / 2
        pt_delta = ptmax - ptmin
        rawy = histos["rawyields"].GetBinContent(i_pt + 1)
        rawy_unc = histos["rawyields"].GetBinError(i_pt + 1)
        eff_times_acc_prompt = histos["acceffp"].GetBinContent(i_pt + 1)
        eff_times_acc_nonprompt = histos["acceffnp"].GetBinContent(i_pt + 1)
        ptmin_fonll = histos["FONLL"]["nonprompt"]["central"].GetXaxis().FindBin(ptmin * 1.0001)
        ptmax_fonll = histos["FONLL"]["nonprompt"]["central"].GetXaxis().FindBin(ptmax * 0.9999)
        crosssec_nonprompt_fonll = [
            histos["FONLL"]["nonprompt"][pred].Integral(ptmin_fonll, ptmax_fonll, "width") / (ptmax - ptmin)
            for pred in histos["FONLL"]["nonprompt"]
        ]

        # compute prompt fraction
        frac = [0.0, 0.0, 0.0]
        if frac_method == "Nb":
            frac = compute_fraction_nb(  # BR already included in FONLL prediction
                rawy,
                eff_times_acc_prompt,
                eff_times_acc_nonprompt,
                crosssec_nonprompt_fonll,
                pt_delta,
                1.0,
                1.0,
                norm["events"],
                norm["sigmaMB"],
            )
        elif frac_method == "fc":
            crosssec_prompt_fonll = [
                histos["FONLL"]["prompt"][pred].Integral(ptmin_fonll, ptmax_fonll, "width") / (ptmax - ptmin)
                for pred in histos["FONLL"]["prompt"]
            ]
            frac, _ = compute_fraction_fc(
                eff_times_acc_prompt,
                eff_times_acc_nonprompt,
                crosssec_prompt_fonll,
                crosssec_nonprompt_fonll,
            )

        # compute cross section times BR
        crosssec, crosssec_unc = compute_crosssection(
            rawy,
            rawy_unc,
            frac[0],
            eff_times_acc_prompt,
            ptmax - ptmin,
            1.0,
            norm["sigmaMB"],
            norm["events"],
            1.0,
            frac_method,
        )

        hptspectrum.SetBinContent(i_pt + 1, crosssec / norm["BR"])
        hptspectrum.SetBinError(i_pt + 1, crosssec_unc / norm["BR"])
        hptspectrum_wo_br.SetBinContent(i_pt + 1, crosssec)
        hptspectrum_wo_br.SetBinError(i_pt + 1, crosssec_unc)
        gfraction.SetPoint(i_pt, pt_cent, frac[0])
        gfraction.SetPointError(i_pt, pt_delta / 2, pt_delta / 2, frac[0] - frac[1], frac[2] - frac[0])

    c = TCanvas("c", "c", 600, 800)
    c.Divide(1, 2)
    c.cd(1)
    gPad.SetLogy(True)
    hptspectrum.Draw()
    c.cd(2)
    gPad.SetLogy(True)
    hptspectrum_wo_br.Draw()

    output_dir = cfg["output"]["directory"]
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    c.Print(os.path.join(output_dir, f'{cfg["output"]["filename"]}.pdf'))

    # save output file
    output_file = TFile(os.path.join(output_dir, f'{cfg["output"]["filename"]}.root'), "recreate")
    hptspectrum.Write()
    hptspectrum_wo_br.Write()
    gfraction.Write()
    for hist in histos:
        if isinstance(histos[hist], TH1):
            histos[hist].Write()
        else:
            for flav in histos[hist]:
                for pred in histos[hist][flav]:
                    histos[hist][flav][pred].Write()
    output_file.Close()


# call main function
main()
