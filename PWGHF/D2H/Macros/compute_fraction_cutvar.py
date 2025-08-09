"""
Script for the (non-)prompt fraction calculation with the cut-variation method

\author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN
\author Fabio Catalano <fabio.catalano@cern.ch>, Politecnico and INFN Torino
\author Stefano Politanò <stefano.politano@cern.ch>, Politecnico and INFN Torino
\author Daniel Battistini <daniel.battistini@cern.ch>, TUM
"""

import argparse
import json
import os
import sys

import numpy as np  # pylint: disable=import-error
import ROOT  # pylint: disable=import-error
sys.path.insert(0, '..')
from cut_variation import CutVarMinimiser
from style_formatter import set_object_style

# pylint: disable=no-member,too-many-locals,too-many-statements


def main(config):
    """
    Main function
    """

    ROOT.gROOT.SetBatch(True)
    ROOT.TH1.AddDirectory(False)

    with open(config, encoding="utf8") as fil:
        cfg = json.load(fil)

    hist_rawy, hist_effp, hist_effnp = ([] for _ in range(3))
    for filename_rawy, filename_eff in zip(cfg["rawyields"]["inputfiles"], cfg["efficiencies"]["inputfiles"]):
        infile_rawy = ROOT.TFile.Open(os.path.join(cfg["rawyields"]["inputdir"], filename_rawy))
        hist_rawy_name = cfg["rawyields"]["histoname"]
        hist_rawy.append(infile_rawy.Get(hist_rawy_name))
        if(hist_rawy[-1] is None):
            sys.exit(f"Fatal error: Histogram with raw yield \"{hist_rawy_name}\" is absent. Exit.")
        hist_rawy[-1].SetDirectory(0)
        infile_rawy.Close()

        infile_eff = ROOT.TFile.Open(os.path.join(cfg["efficiencies"]["inputdir"], filename_eff))
        hist_effp_name = cfg["efficiencies"]["histonames"]["prompt"]
        hist_effnp_name = cfg["efficiencies"]["histonames"]["nonprompt"]
        hist_effp.append(infile_eff.Get(hist_effp_name))
        hist_effnp.append(infile_eff.Get(hist_effnp_name))
        if(hist_effp[-1] is None):
            sys.exit(f"Fatal error: Histogram with efficiency for prompt \"{hist_effp_name}\" is absent. Exit.")
        if(hist_effnp[-1] is None):
            sys.exit(f"Fatal error: Histogram with efficiency for nonprompt \"{hist_effnp}\" is absent. Exit.")
        hist_effp[-1].SetDirectory(0)
        hist_effnp[-1].SetDirectory(0)
        infile_eff.Close()

    pt_bin_to_process = cfg.get("pt_bin_to_process", -1)
    if not isinstance(pt_bin_to_process, int):
        sys.exit("Fatal error: pt_bin_to_process must be an integer value. Exit.")
    if (pt_bin_to_process != -1 and pt_bin_to_process < 1) or pt_bin_to_process > hist_rawy[0].GetNbinsX():
        sys.exit("Fatal error: pt_bin_to_process must be a positive value up to number of bins in raw yield histogram. Exit.")

    if cfg["central_efficiency"]["computerawfrac"]:
        infile_name = os.path.join(cfg["central_efficiency"]["inputdir"], cfg["central_efficiency"]["inputfile"])
        infile_central_eff = ROOT.TFile.Open(infile_name)
        hist_central_effp_name = cfg["central_efficiency"]["histonames"]["prompt"]
        hist_central_effp = infile_central_eff.Get(hist_central_effp_name)
        if(hist_central_effp is None):
            sys.exit(f"Fatal error: Histogram with central efficiency for prompt \"{hist_central_effp_name}\" is absent. Exit.")
        hist_central_effnp_name = cfg["central_efficiency"]["histonames"]["nonprompt"]
        hist_central_effnp = infile_central_eff.Get(hist_central_effnp_name)
        if(hist_central_effnp is None):
            sys.exit(f"Fatal error: Histogram with central efficiency for nonprompt \"{hist_central_effnp_name}\" is absent. Exit.")
        hist_central_effp.SetDirectory(0)
        hist_central_effnp.SetDirectory(0)
        infile_central_eff.Close()

        # Check if the histograms have the same binning, assuming the two sets of inputs are consistent
        for ibin in range(1, hist_rawy[0].GetNbinsX() + 2):  # +2 to include the upper edge of last bin
            if hist_rawy[0].GetBinLowEdge(ibin) != hist_central_effp.GetBinLowEdge(ibin):
                raise ValueError("Histograms have different binning, check the input files")

    hist_corry_prompt = hist_rawy[0].Clone("hCorrYieldsPrompt")
    hist_corry_nonprompt = hist_rawy[0].Clone("hCorrYieldsNonPrompt")
    hist_covariance_pnp = hist_rawy[0].Clone("hCovPromptNonPrompt")
    hist_covariance_pp = hist_rawy[0].Clone("hCovPromptPrompt")
    hist_covariance_npnp = hist_rawy[0].Clone("hCovNonPromptNonPrompt")
    hist_corrfrac_prompt = hist_rawy[0].Clone("hCorrFracPrompt")
    hist_corrfrac_nonprompt = hist_rawy[0].Clone("hCorrFracNonPrompt")
    for histo in hist_corry_prompt, hist_corry_nonprompt, hist_covariance_pnp, hist_covariance_pp, hist_covariance_npnp, hist_corrfrac_prompt, hist_corrfrac_nonprompt:
        histo.Reset()
    hist_corry_prompt.GetYaxis().SetTitle("corrected yields prompt")
    hist_corry_nonprompt.GetYaxis().SetTitle("corrected yields non-prompt")
    hist_covariance_pnp.GetYaxis().SetTitle("#sigma(prompt, non-prompt)")
    hist_covariance_pp.GetYaxis().SetTitle("#sigma(prompt, prompt)")
    hist_covariance_npnp.GetYaxis().SetTitle("#sigma(non-prompt, non-prompt)")
    hist_corrfrac_prompt.GetYaxis().SetTitle("corrected fraction prompt")
    hist_corrfrac_nonprompt.GetYaxis().SetTitle("corrected fraction non-prompt")
    set_object_style(
        hist_corry_prompt,
        color=ROOT.kRed + 1,
        fillstyle=0,
        markerstyle=ROOT.kFullCircle,
    )
    set_object_style(
        hist_corry_nonprompt,
        color=ROOT.kAzure + 4,
        fillstyle=0,
        markerstyle=ROOT.kFullSquare,
    )
    set_object_style(hist_covariance_pnp)
    set_object_style(hist_covariance_pp)
    set_object_style(hist_covariance_npnp)
    set_object_style(
        hist_corrfrac_prompt,
        color=ROOT.kRed + 1,
        fillstyle=0,
        markerstyle=ROOT.kFullCircle,
    )
    set_object_style(
        hist_corrfrac_nonprompt,
        color=ROOT.kAzure + 4,
        fillstyle=0,
        markerstyle=ROOT.kFullSquare,
    )
    if cfg["central_efficiency"]["computerawfrac"]:
        hist_frac_raw_prompt = hist_rawy[0].Clone("hRawFracPrompt")
        hist_frac_raw_nonprompt = hist_rawy[0].Clone("hRawFracNonPrompt")
        for histo in hist_frac_raw_prompt, hist_frac_raw_nonprompt:
            histo.Reset()
        hist_frac_raw_prompt.GetYaxis().SetTitle("raw fraction prompt")
        hist_frac_raw_nonprompt.GetYaxis().SetTitle("raw fraction non-prompt")
        set_object_style(
            hist_frac_raw_prompt,
            color=ROOT.kRed + 1,
            fillstyle=0,
            markerstyle=ROOT.kFullSquare,
        )

        set_object_style(
            hist_frac_raw_nonprompt,
            color=ROOT.kAzure + 4,
            fillstyle=0,
            markerstyle=ROOT.kFullSquare,
        )

    pt_bin_to_process_name_suffix = ""
    if pt_bin_to_process != -1:
        pt_bin_to_process_name_suffix = "_bin_" + str(pt_bin_to_process)

    output_name_template = cfg['output']['file'].replace(".root", "") + pt_bin_to_process_name_suffix + ".root"
    output = ROOT.TFile(os.path.join(cfg["output"]["directory"], output_name_template), "recreate")
    n_sets = len(hist_rawy)
    pt_axis_title = hist_rawy[0].GetXaxis().GetTitle()
    for ipt in range(hist_rawy[0].GetNbinsX()):
        if pt_bin_to_process !=-1 and ipt+1 != pt_bin_to_process:
            continue
        pt_min = hist_rawy[0].GetXaxis().GetBinLowEdge(ipt + 1)
        pt_max = hist_rawy[0].GetXaxis().GetBinUpEdge(ipt + 1)
        print(f"\n\nINFO: processing pt range {ipt+1} from {pt_min} to {pt_max} {pt_axis_title}")

        rawy, effp, effnp, unc_rawy, unc_effp, unc_effnp = (np.zeros(n_sets) for _ in range(6))
        for iset, (hrawy, heffp, heffnp) in enumerate(zip(hist_rawy, hist_effp, hist_effnp)):
            rawy[iset] = hrawy.GetBinContent(ipt + 1)
            effp[iset] = heffp.GetBinContent(ipt + 1)
            effnp[iset] = heffnp.GetBinContent(ipt + 1)
            unc_rawy[iset] = hrawy.GetBinError(ipt + 1)
            unc_effp[iset] = heffp.GetBinError(ipt + 1)
            unc_effnp[iset] = heffnp.GetBinError(ipt + 1)

        if cfg["minimisation"]["correlated"]:
            if not (np.all(rawy[1:] > rawy[:-1]) or np.all(rawy[1:] < rawy[:-1])):
                print("WARNING! main(): the raw yield vector is not monotonous. Check the input for stability.")
                print(f"raw yield vector elements = {rawy}\n")
            if not (np.all(unc_rawy[1:] > unc_rawy[:-1]) or np.all(unc_rawy[1:] < unc_rawy[:-1])):
                print("WARNING! main(): the raw yield uncertainties vector is not monotonous. Check the input for stability.")
                print(f"raw yield uncertainties vector elements = {unc_rawy}\n")

        minimiser = CutVarMinimiser(rawy, effp, effnp, unc_rawy, unc_effp, unc_effnp)
        status = minimiser.minimise_system(cfg["minimisation"]["correlated"])

        if status:
            hist_corry_prompt.SetBinContent(ipt + 1, minimiser.get_prompt_yield_and_error()[0])
            hist_corry_prompt.SetBinError(ipt + 1, minimiser.get_prompt_yield_and_error()[1])
            hist_corry_nonprompt.SetBinContent(ipt + 1, minimiser.get_nonprompt_yield_and_error()[0])
            hist_corry_nonprompt.SetBinError(ipt + 1, minimiser.get_nonprompt_yield_and_error()[1])
            hist_covariance_pnp.SetBinContent(ipt + 1, minimiser.get_prompt_nonprompt_cov())
            hist_covariance_pnp.SetBinError(ipt + 1, 0)
            hist_covariance_pp.SetBinContent(ipt + 1, minimiser.get_prompt_prompt_cov())
            hist_covariance_pp.SetBinError(ipt + 1, 0)
            hist_covariance_npnp.SetBinContent(ipt + 1, minimiser.get_nonprompt_nonprompt_cov())
            hist_covariance_npnp.SetBinError(ipt + 1, 0)
            corr_frac_prompt = minimiser.get_corr_prompt_fraction()
            corr_frac_nonprompt = minimiser.get_corr_nonprompt_fraction()
            hist_corrfrac_prompt.SetBinContent(ipt + 1, corr_frac_prompt[0])
            hist_corrfrac_prompt.SetBinError(ipt + 1, corr_frac_prompt[1])
            hist_corrfrac_nonprompt.SetBinContent(ipt + 1, corr_frac_nonprompt[0])
            hist_corrfrac_nonprompt.SetBinError(ipt + 1, corr_frac_nonprompt[1])
            if cfg["central_efficiency"]["computerawfrac"]:
                raw_frac_prompt = minimiser.get_raw_prompt_fraction(
                    hist_central_effp.GetBinContent(ipt + 1), hist_central_effnp.GetBinContent(ipt + 1)
                )
                raw_frac_nonprompt = minimiser.get_raw_nonprompt_fraction(
                    hist_central_effp.GetBinContent(ipt + 1), hist_central_effnp.GetBinContent(ipt + 1)
                )
                hist_frac_raw_prompt.SetBinContent(ipt + 1, raw_frac_prompt[0])
                hist_frac_raw_prompt.SetBinError(ipt + 1, raw_frac_prompt[1])
                hist_frac_raw_nonprompt.SetBinContent(ipt + 1, raw_frac_nonprompt[0])
                hist_frac_raw_nonprompt.SetBinError(ipt + 1, raw_frac_nonprompt[1])

            hist_bin_title = f"bin # {ipt+1}; {pt_axis_title}#in ({pt_min}; {pt_max})"

            canv_rawy, histos_rawy, leg_r = minimiser.plot_result(f"_pt{pt_min}_{pt_max}", hist_bin_title)
            output.cd()
            canv_rawy.Write()
            for _, hist in histos_rawy.items():
                hist.Write()

            canv_unc, histos_unc, leg_unc = minimiser.plot_uncertainties(f"_pt{pt_min}_{pt_max}", hist_bin_title)
            output.cd()
            canv_unc.Write()
            for _, hist in histos_unc.items():
                hist.Write()

            canv_eff, histos_eff, leg_e = minimiser.plot_efficiencies(f"_pt{pt_min}_{pt_max}", hist_bin_title)
            output.cd()
            canv_eff.Write()
            for _, hist in histos_eff.items():
                hist.Write()

            canv_frac, histos_frac, leg_f = minimiser.plot_fractions(f"_pt{pt_min}_{pt_max}", hist_bin_title)
            output.cd()
            canv_frac.Write()
            for _, hist in histos_frac.items():
                hist.Write()

            canv_cov, histo_cov = minimiser.plot_cov_matrix(True, f"_pt{pt_min}_{pt_max}", hist_bin_title)
            output.cd()
            canv_cov.Write()
            histo_cov.Write()
        else:
            print(f"Minimization for pT {pt_min}, {pt_max} not successful")
            canv_rawy = ROOT.TCanvas("c_rawy_minimization_error", "Minimization error", 500, 500)
            canv_eff = ROOT.TCanvas("c_eff_minimization_error", "Minimization error", 500, 500)
            canv_frac = ROOT.TCanvas("c_frac_minimization_error", "Minimization error", 500, 500)
            canv_cov = ROOT.TCanvas("c_conv_minimization_error", "Minimization error", 500, 500)

        canv_combined = ROOT.TCanvas(f"canv_combined_{ipt}", "", 1000, 1000)
        canv_combined.Divide(2, 2)
        canv_combined.cd(1)
        canv_rawy.DrawClonePad()
        canv_combined.cd(2)
        canv_eff.DrawClonePad()
        canv_combined.cd(3)
        canv_frac.DrawClonePad()
        canv_combined.cd(4)
        canv_cov.DrawClonePad()

        output_name_template = output_name_template.replace('.root', '.pdf')

        output_name_rawy_pdf = f"Distr_{output_name_template}"
        output_name_eff_pdf = f"Eff_{output_name_template}"
        output_name_frac_pdf = f"Frac_{output_name_template}"
        output_name_covmat_pdf = f"CovMatrix_{output_name_template}"
        output_name_unc_pdf = f"Unc_{output_name_template}"
        output_name_pdf = f"{output_name_template}"

        if hist_rawy[0].GetNbinsX() == 1 or pt_bin_to_process != -1:
            print_bracket = ""
        elif ipt == 0:
            print_bracket = "("
        elif ipt == hist_rawy[0].GetNbinsX() - 1:
            print_bracket = ")"
        else:
            print_bracket = ""
        canv_rawy.Print(f"{os.path.join(cfg['output']['directory'], output_name_rawy_pdf)}{print_bracket}")
        canv_eff.Print(f"{os.path.join(cfg['output']['directory'], output_name_eff_pdf)}{print_bracket}")
        canv_frac.Print(f"{os.path.join(cfg['output']['directory'], output_name_frac_pdf)}{print_bracket}")
        canv_cov.Print(f"{os.path.join(cfg['output']['directory'], output_name_covmat_pdf)}{print_bracket}")
        canv_combined.Print(f"{os.path.join(cfg['output']['directory'], output_name_pdf)}{print_bracket}")
        canv_unc.Print(f"{os.path.join(cfg['output']['directory'], output_name_unc_pdf)}{print_bracket}")

    output.cd()
    hist_corry_prompt.Write()
    hist_corry_nonprompt.Write()
    hist_covariance_pnp.Write()
    hist_covariance_pp.Write()
    hist_covariance_npnp.Write()
    hist_corrfrac_prompt.Write()
    hist_corrfrac_nonprompt.Write()
    if cfg["central_efficiency"]["computerawfrac"]:
        hist_frac_raw_prompt.Write()
        hist_frac_raw_nonprompt.Write()
    output.Close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument(
        "config",
        metavar="text",
        default="config_cutvar_example.json",
        help="JSON config file",
    )
    args = parser.parse_args()

    main(args.config)
