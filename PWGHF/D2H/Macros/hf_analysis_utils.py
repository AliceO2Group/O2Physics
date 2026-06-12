#!/usr/bin/env python3

"""
file: hf_analysis_utils.py
brief: script with miscellanea utils methods for the HF analyses
author: Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN
"""

import numpy as np  # pylint: disable=import-error


def make_list(object) -> list:
    """
    Returns the object as a list if it is not a list already.
    """
    return object if isinstance(object, list) else [object]


# pylint: disable=too-many-arguments
def compute_crosssection(
    rawy,
    rawy_unc,
    frac,
    eff_times_acc,
    delta_pt,
    delta_y,
    sigma_mb,
    n_events,
    br,
    method_frac="Nb",
):
    """
    Method to compute cross section and its statistical uncertainty
    Only the statistical uncertainty on the raw yield and prompt (non-prompt)
    fraction are considered (the others are systematics)

    Parameters
    ----------
    - rawy: raw yield
    - rawy_unc: raw-yield statistical uncertainty
    - frac: either prompt or non-prompt fraction
    - eff_times_acc: efficiency times acceptance for prompt or non-prompt
    - delta_pt: pT interval
    - delta_y: Y interval
    - sigma_mb: hadronic cross section for MB
    - n_events: number of events
    - br: branching ratio of the decay channel
    - method_frac: method used to compute frac needed to propoer compute uncertainty

    Returns
    ----------
    - crosssection: cross section
    - crosssec_unc: cross-section statistical uncertainty
    """

    crosssection = rawy * frac * sigma_mb / (2 * delta_pt * delta_y * eff_times_acc * n_events * br)
    if method_frac == "Nb":
        crosssec_unc = rawy_unc / (rawy * frac) * crosssection
    else:
        crosssec_unc = rawy_unc / rawy * crosssection

    return crosssection, crosssec_unc


# pylint: disable=too-many-branches,too-many-arguments,too-many-locals
def compute_fraction_fc(
    acc_eff_prompt,
    acc_eff_fd,
    cross_sec_prompt,
    cross_sec_fd,
    raa_prompt=1.0,
    raa_fd=1.0,
) -> "tuple[list[float], list[float]]":
    """
    Method to get fraction of prompt / FD fraction with fc method

    Parameters
    ----------
    - acc_eff_prompt: efficiency times acceptance of prompt D
    - acc_eff_fd: efficiency times acceptance of non-prompt D
    - cross_sec_prompt: list of production cross sections (cent, min, max)
                        of prompt D in pp collisions from theory
    - cross_sec_fd: list of production cross sections (cent, min, max)
                    of non-prompt D in pp collisions from theory
    - raa_prompt: list of nuclear modification factors (cent, min, max) of prompt D from theory
    - raa_fd: list of nuclear modification factors of (cent, min, max) non-prompt D from theory

    Returns
    ----------
    - frac_prompt: list of fraction of prompt D (central, min, max)
    - frac_fd: list of fraction of non-prompt D (central, min, max)
    """

    cross_sec_prompt = make_list(cross_sec_prompt)
    cross_sec_fd = make_list(cross_sec_fd)
    raa_prompt = make_list(raa_prompt)
    raa_fd = make_list(raa_fd)

    frac_prompt: list[float] = []
    frac_fd: list[float] = []
    if acc_eff_prompt == 0:
        frac_fd_cent = 1.0
        frac_prompt_cent = 0.0
        frac_prompt = [frac_prompt_cent, frac_prompt_cent, frac_prompt_cent]
        frac_fd = [frac_fd_cent, frac_fd_cent, frac_fd_cent]
        return frac_prompt, frac_fd
    if acc_eff_fd == 0:
        frac_fd_cent = 0.0
        frac_prompt_cent = 1.0
        frac_prompt = [frac_prompt_cent, frac_prompt_cent, frac_prompt_cent]
        frac_fd = [frac_fd_cent, frac_fd_cent, frac_fd_cent]
        return frac_prompt, frac_fd

    for i_sigma, (sigma_p, sigma_f) in enumerate(zip(cross_sec_prompt, cross_sec_fd)):
        for i_raa, (raa_p, raa_f) in enumerate(zip(raa_prompt, raa_fd)):
            if i_sigma == 0 and i_raa == 0:
                frac_prompt_cent = 1.0 / (1 + acc_eff_fd / acc_eff_prompt * sigma_f / sigma_p * raa_f / raa_p)
                frac_fd_cent = 1.0 / (1 + acc_eff_prompt / acc_eff_fd * sigma_p / sigma_f * raa_p / raa_f)
            else:
                frac_prompt.append(1.0 / (1 + acc_eff_fd / acc_eff_prompt * sigma_f / sigma_p * raa_f / raa_p))
                frac_fd.append(1.0 / (1 + acc_eff_prompt / acc_eff_fd * sigma_p / sigma_f * raa_p / raa_f))

    if frac_prompt and frac_fd:
        frac_prompt.sort()
        frac_fd.sort()
        frac_prompt = [frac_prompt_cent, frac_prompt[0], frac_prompt[-1]]
        frac_fd = [frac_fd_cent, frac_fd[0], frac_fd[-1]]
    else:
        frac_prompt = [frac_prompt_cent, frac_prompt_cent, frac_prompt_cent]
        frac_fd = [frac_fd_cent, frac_fd_cent, frac_fd_cent]

    return frac_prompt, frac_fd


# pylint: disable=too-many-branches,too-many-arguments,too-many-locals,invalid-name
def compute_fraction_nb(
    rawy,
    acc_eff_same,
    acc_eff_other,
    crosssection,
    delta_pt,
    delta_y,
    br,
    n_events,
    sigma_mb,
    raa_ratio=1.0,
    taa=1.0,
) -> "list[float]":
    """
    Method to get fraction of prompt / FD fraction with Nb method

    Parameters
    ----------
    - acc_eff_same: efficiency times acceptance of prompt (non-prompt) D
    - acc_eff_other: efficiency times acceptance of non-prompt (prompt) D
    - crosssection: list of production cross sections (cent, min, max) of non-prompt (prompt)
                D in pp collisions from theory
    - delta_pt: width of pT interval
    - delta_y: width of Y interval
    - br: branching ratio for the chosen decay channel
    - n_events: number of events corresponding to the raw yields
    - sigma_mb: MB cross section
    - raa_ratio: list of D nuclear modification factor ratios
                non-prompt / prompt (prompt / non-prompt) (cent, min, max) (=1 in case of pp)
    - taa: average nuclear overlap function (=1 in case of pp)

    Returns
    ----------
    - frac: list of fraction of prompt (non-prompt) D (central, min, max)
    """

    crosssection = make_list(crosssection)
    raa_ratio = make_list(raa_ratio)

    frac: list[float] = []
    for i_sigma, sigma in enumerate(crosssection):
        for i_raa_ratio, raa_rat in enumerate(raa_ratio):
            raa_other = 1.0
            if i_sigma == 0 and i_raa_ratio == 0:
                if raa_rat == 1.0 and taa == 1.0:  # pp
                    frac_cent = 1 - sigma * delta_pt * delta_y * acc_eff_other * br * n_events * 2 / rawy / sigma_mb
                else:  # p-Pb or Pb-Pb: iterative evaluation of Raa needed
                    delta_raa = 1.0
                    while delta_raa > 1.0e-3:
                        raw_fd = (
                            taa * raa_rat * raa_other * sigma * delta_pt * delta_y * acc_eff_other * br * n_events * 2
                        )
                        frac_cent = 1 - raw_fd / rawy
                        raa_other_old = raa_other
                        raa_other = frac_cent * rawy * sigma_mb / 2 / acc_eff_same / delta_pt / delta_y / br / n_events
                        delta_raa = abs((raa_other - raa_other_old) / raa_other)
            else:
                if raa_rat == 1.0 and taa == 1.0:  # pp
                    frac.append(1 - sigma * delta_pt * delta_y * acc_eff_other * br * n_events * 2 / rawy / sigma_mb)
                else:  # p-Pb or Pb-Pb: iterative evaluation of Raa needed
                    delta_raa = 1.0
                    frac_tmp = 1.0
                    while delta_raa > 1.0e-3:
                        raw_fd = (
                            taa * raa_rat * raa_other * sigma * delta_pt * delta_y * acc_eff_other * br * n_events * 2
                        )
                        frac_tmp = 1 - raw_fd / rawy
                        raa_other_old = raa_other
                        raa_other = frac_tmp * rawy * sigma_mb / 2 / acc_eff_same / delta_pt / delta_y / br / n_events
                        delta_raa = abs((raa_other - raa_other_old) / raa_other)
                    frac.append(frac_tmp)

    if frac:
        frac.sort()
        frac = [frac_cent, frac[0], frac[-1]]
    else:
        frac = [frac_cent, frac_cent, frac_cent]

    return frac


def get_hist_binlimits(histo):
    """
    Method to retrieve bin limits of ROOT.TH1

    Parameters
    ----------
    - histo: ROOT.TH1

    Returns
    ----------
    - bin_limits: numpy array of bin limits
    """

    if np.array(histo.GetXaxis().GetXbins(), "d").any():  # variable binning
        bin_limits = np.array(histo.GetXaxis().GetXbins(), "d")
    else:  # constant binning
        n_limits = histo.GetNbinsX() + 1
        low_edge = histo.GetBinLowEdge(1)
        bin_width = histo.GetBinWidth(1)
        bin_limits = np.array([low_edge + i_bin * bin_width for i_bin in range(n_limits)], "d")

    return bin_limits
