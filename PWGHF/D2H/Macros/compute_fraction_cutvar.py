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

import numpy as np
import ROOT
from cut_variation import CutVarMinimiser
from style_formatter import set_object_style

# pylint: disable=no-member,too-many-locals,too-many-statements


def main(config):
    """
    Main function
    """

    ROOT.gROOT.SetBatch(True)

    with open(config, encoding="utf8") as fil:
        cfg = json.load(fil)

    hist_rawy, hist_effp, hist_effnp = ([] for _ in range(3))
    for filename_rawy, filename_eff in zip(
        cfg["rawyields"]["inputfiles"], cfg["efficiencies"]["inputfiles"]
    ):
        infile_rawy = ROOT.TFile.Open(
            os.path.join(cfg["rawyields"]["inputdir"], filename_rawy)
        )
        hist_rawy.append(infile_rawy.Get(cfg["rawyields"]["histoname"]))
        hist_rawy[-1].SetDirectory(0)
        infile_rawy.Close()

        infile_eff = ROOT.TFile.Open(
            os.path.join(cfg["efficiencies"]["inputdir"], filename_eff)
        )
        hist_effp.append(infile_eff.Get(cfg["efficiencies"]["histonames"]["prompt"]))
        hist_effnp.append(
            infile_eff.Get(cfg["efficiencies"]["histonames"]["nonprompt"])
        )
        hist_effp[-1].SetDirectory(0)
        hist_effnp[-1].SetDirectory(0)
        infile_eff.Close()

    hist_corry_prompt = hist_rawy[0].Clone("hCorrYieldsPrompt")
    hist_corry_nonprompt = hist_rawy[0].Clone("hCorrYieldsNonPrompt")
    hist_covariance = hist_rawy[0].Clone("hCovPromptNonPrompt")
    hist_corry_prompt.GetYaxis().SetTitle("corrected yields prompt")
    hist_corry_nonprompt.GetYaxis().SetTitle("corrected yields non-prompt")
    hist_covariance.GetYaxis().SetTitle("#sigma(prompt, non-prompt)")
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
    set_object_style(hist_covariance)

    output = ROOT.TFile(
        os.path.join(cfg["output"]["directory"], cfg["output"]["file"]), "recreate"
    )
    n_sets = len(hist_rawy)
    for ipt in range(hist_rawy[0].GetNbinsX()):
        pt_min = hist_rawy[0].GetXaxis().GetBinLowEdge(ipt + 1)
        pt_max = hist_rawy[0].GetXaxis().GetBinUpEdge(ipt + 1)

        rawy, effp, effnp, unc_rawy, unc_effp, unc_effnp = (
            np.zeros(n_sets) for _ in range(6)
        )
        for iset, (hrawy, heffp, heffnp) in enumerate(
            zip(hist_rawy, hist_effp, hist_effnp)
        ):
            rawy.itemset(iset, hrawy.GetBinContent(ipt + 1))
            effp.itemset(iset, heffp.GetBinContent(ipt + 1))
            effnp.itemset(iset, heffnp.GetBinContent(ipt + 1))
            unc_rawy.itemset(iset, hrawy.GetBinError(ipt + 1))
            unc_effp.itemset(iset, heffp.GetBinError(ipt + 1))
            unc_effnp.itemset(iset, heffnp.GetBinError(ipt + 1))

        minimiser = CutVarMinimiser(rawy, effp, effnp, unc_rawy, unc_effp, unc_effnp)
        minimiser.minimise_system(cfg["minimisation"]["correlated"])

        hist_corry_prompt.SetBinContent(
            ipt + 1, minimiser.get_prompt_yield_and_error()[0]
        )
        hist_corry_prompt.SetBinError(
            ipt + 1, minimiser.get_prompt_yield_and_error()[1]
        )
        hist_corry_nonprompt.SetBinContent(
            ipt + 1, minimiser.get_nonprompt_yield_and_error()[0]
        )
        hist_corry_nonprompt.SetBinError(
            ipt + 1, minimiser.get_nonprompt_yield_and_error()[1]
        )
        hist_covariance.SetBinContent(ipt + 1, minimiser.get_prompt_nonprompt_cov())
        hist_covariance.SetBinError(ipt + 1, 0)

        canv_rawy, histos_rawy, leg_r = minimiser.plot_result(
            f"_pt{pt_min:.0f}_{pt_max:.0f}"
        )
        output.cd()
        canv_rawy.Write()
        for _, hist in histos_rawy.items():
            hist.Write()

        canv_eff, histos_eff, leg_e = minimiser.plot_efficiencies(
            f"_pt{pt_min:.0f}_{pt_max:.0f}"
        )
        output.cd()
        canv_eff.Write()
        for _, hist in histos_eff.items():
            hist.Write()

        canv_frac, histos_frac, leg_f = minimiser.plot_fractions(
            f"_pt{pt_min:.0f}_{pt_max:.0f}"
        )
        output.cd()
        canv_frac.Write()
        for _, hist in histos_frac.items():
            hist.Write()

        canv_cov, histo_cov = minimiser.plot_cov_matrix(f"_pt{pt_min:.0f}_{pt_max:.0f}")
        output.cd()
        canv_cov.Write()
        histo_cov.Write()

        output_name_rawy_pdf = f"Distr_{cfg['output']['file'].replace('.root', '.pdf')}"
        output_name_eff_pdf = f"Eff_{cfg['output']['file'].replace('.root', '.pdf')}"
        output_name_frac_pdf = f"Frac_{cfg['output']['file'].replace('.root', '.pdf')}"
        output_name_covmat_pdf = (
            f"CovMatrix_{cfg['output']['file'].replace('.root', '.pdf')}"
        )
        if ipt == 0:
            canv_rawy.SaveAs(
                f"{os.path.join(cfg['output']['directory'], output_name_rawy_pdf)}["
            )
            canv_eff.SaveAs(
                f"{os.path.join(cfg['output']['directory'], output_name_eff_pdf)}["
            )
            canv_frac.SaveAs(
                f"{os.path.join(cfg['output']['directory'], output_name_frac_pdf)}["
            )
            canv_cov.SaveAs(
                f"{os.path.join(cfg['output']['directory'], output_name_covmat_pdf)}["
            )
        canv_rawy.SaveAs(
            f"{os.path.join(cfg['output']['directory'], output_name_rawy_pdf)}"
        )
        canv_eff.SaveAs(
            f"{os.path.join(cfg['output']['directory'], output_name_eff_pdf)}"
        )
        canv_frac.SaveAs(
            f"{os.path.join(cfg['output']['directory'], output_name_frac_pdf)}"
        )
        canv_cov.SaveAs(
            f"{os.path.join(cfg['output']['directory'], output_name_covmat_pdf)}"
        )
        if ipt == hist_rawy[0].GetNbinsX() - 1:
            canv_rawy.SaveAs(
                f"{os.path.join(cfg['output']['directory'], output_name_rawy_pdf)}]"
            )
            canv_eff.SaveAs(
                f"{os.path.join(cfg['output']['directory'], output_name_eff_pdf)}]"
            )
            canv_frac.SaveAs(
                f"{os.path.join(cfg['output']['directory'], output_name_frac_pdf)}]"
            )
            canv_cov.SaveAs(
                f"{os.path.join(cfg['output']['directory'], output_name_covmat_pdf)}]"
            )

    output.cd()
    hist_corry_prompt.Write()
    hist_corry_nonprompt.Write()
    hist_covariance.Write()

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
