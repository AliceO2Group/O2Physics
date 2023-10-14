#!/usr/bin/env python3

"""
file: efficiency.py
brief: script for (acceptance-times-)efficiency computation
usage: python3 efficiency.py CONFIG
author: Alexandre Bigot <alexandre.bigot@cern.ch>, Strasbourg University
"""

import argparse
import os
import sys
from typing import Union

import yaml
#pylint: disable=no-name-in-module
from ROOT import (TH1, TH2, TCanvas, TEfficiency, TFile, TLatex, TLegend,
                  kAzure, kBlack, kFullCircle, kFullSquare, kOpenSquare, kRed)
from style_formatter import set_global_style, set_object_style

# options for object style
COLOR = {'All': kBlack,
         'Prompt': kRed+1,
         'Nonprompt': kAzure+4}
MARKERSTYLE = {'All': kFullSquare,
               'Prompt': kFullCircle,
               'Nonprompt': kOpenSquare}
LINEWIDTH = {'All': 1,
         'Prompt': 1,
         'Nonprompt': 2}

LINESTYLE = {'All': 1,
         'Prompt': 2,
         'Nonprompt': 7}

def enforce_list(x: Union[str, list, None]) -> Union[list, None]:
    """
    Helper method to enforce list type

    Parameters
    ----------
    - x: a string or a list of string

    Returns
    ----------
    - x_list if x was not a list (and not None), x itself otherwise
    """

    if not isinstance(x, list) and x is not None:
        # handle possible whitespaces in config file entry
        x_list = x.split(',')
        for i, element in enumerate(x_list):
            x_list[i] = element.strip() # remove possible whitespaces
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

    if path[-1] != '/':
        path += '/'

    return path

def save_canvas(
    canvas: TCanvas,
    out_dir: str,
    name_file: str,
    extension: list) -> None:
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
        canvas.SaveAs(out_dir + name_file + '.' + ext)

def compute_efficiency(
    h_rec: Union[TH1, TH2],
    h_gen: Union[TH1, TH2],
    rapidity_cut: float = None,
    axis_rapidity: str = 'Y') -> TEfficiency:
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
        if axis_rapidity == 'Y':
            h_rec = h_rec.ProjectionX()
        else:
            h_rec = h_rec.ProjectionY()

    if is_gen_2d:
        if rapidity_cut is not None:
            h_gen.GetYaxis().SetRangeUser(-1.0 * rapidity_cut + epsilon, rapidity_cut - epsilon)
        if axis_rapidity == 'Y':
            h_gen = h_gen.ProjectionX()
        else:
            h_gen = h_gen.ProjectionY()

    efficiency = TEfficiency(h_rec, h_gen)

    return efficiency

def get_th1_from_tefficiency(
    teff: TEfficiency,
    h_eff: TH1) -> TH1:
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
    for i_bin in range(1, n_bins+1):
        h_eff.SetBinContent(i_bin, teff.GetEfficiency(i_bin))
        # TH1 can't handle asymmetric errors so we take the max
        err_max = max(teff.GetEfficiencyErrorLow(i_bin), teff.GetEfficiencyErrorUp(i_bin))
        h_eff.SetBinError(i_bin, err_max)

    return h_eff

#pylint: disable=too-many-arguments
def configure_canvas(
    name_canvas: str,
    pt_min: Union[int, float],
    min_y_axis: Union[int, float],
    pt_max: Union[int, float],
    max_y_axis: Union[int, float],
    title: str,
    log_y_axis: bool) -> TCanvas:
    """
    Helper method to configure canvas

    Parameters
    ----------
    - name_canvas: name of the canvas
    - pt_min: lower limit of x axis
    - min_y_axis: lower limit of y axis
    - pt_max: upper limit of x axis
    - max_y_axis: upper limit of y axis
    - title: title of the canvas
    - log_y_axis: switch for log scale along y axis

    Returns
    ----------
    - c_eff: TCanvas instance
    """

    c_eff = TCanvas(name_canvas, '', 800, 800)
    c_eff.DrawFrame(pt_min, min_y_axis, pt_max, max_y_axis, title)
    if log_y_axis:
        c_eff.SetLogy()

    return c_eff

def __set_object_style(obj: Union[TEfficiency, TH1], key: str) -> None:
    """
    Helper method to set style of TEfficiency or TH1 object with local options

    Parameters
    ----------
    - obj: TEfficiency or TH1 object
    - key: key of options dictionary
    """

    set_object_style(
                obj,
                color=COLOR[key],
                markerstyle=MARKERSTYLE[key],
                linewidth=LINEWIDTH[key],
                linestyle=LINESTYLE[key])

#pylint: disable=too-many-locals, too-many-branches, too-many-statements
def main(
    cfg: dict,
    batch: bool) -> None:
    """
    Main function

    Parameters
    -----------------
    - config: dictionary with config read from a yaml file
    - batch: bool to suppress video output
    """

    set_global_style(padleftmargin=0.14, padbottommargin=0.12, titlesize=0.045, labelsize=0.04)

    # input configuration
    name_file = cfg['input']['filename']
    name_tree = cfg['input']['treename']
    file = TFile.Open(name_file)
    name_tree = enforce_trailing_slash(name_tree)

    # fill dictionaries with histograms
    dic_rec, dic_gen = {}, {}
    for key, histoname in cfg['input']['histoname']['reconstructed'].items():
        if histoname is None:
            dic_rec[key] = None
        else:
            dic_rec[key] = file.Get(name_tree + histoname)
    for key, histoname in cfg['input']['histoname']['generated'].items():
        if histoname is None:
            dic_gen[key] = None
        else:
            dic_gen[key] = file.Get(name_tree + histoname)

    #  configure possible cut on rapidity
    rapidity_cut = cfg['rapidity']['cut']
    axis_rapidity = cfg['rapidity']['axis']

    # output configuration
    out_dir = enforce_trailing_slash(cfg['output']['dir'])
    out_label = cfg['output']['plots']['label']
    name_axis = cfg['output']['plots']['y_axis']['name']
    min_y_axis = cfg['output']['plots']['y_axis']['min']
    max_y_axis = cfg['output']['plots']['y_axis']['max']
    log_y_axis = cfg['output']['plots']['y_axis']['log_scale']
    if name_axis is None:
        name_axis = '#varepsilon' # default value
    pt_min, pt_max = 0., 35. # default values
    title = ';#it{p}_{T} (GeV/#it{c});' + name_axis + ';'
    overlap = enforce_list(cfg['output']['plots']['overlap'])

    # output save options
    save_tefficiency = cfg['output']['save']['TEfficiency']
    save_th1 = cfg['output']['save']['TH1']
    save_tcanvas_individual = cfg['output']['save']['TCanvas']['individual']
    save_tcanvas_overlap = cfg['output']['save']['TCanvas']['overlap']
    extension = enforce_list(cfg['output']['save']['TCanvas']['extension'])

    if os.path.isdir(out_dir):
        print((f'\033[93mWARNING: Output directory \'{out_dir}\' already exists,'
                    ' overwrites possibly ongoing!\033[0m'))
    else:
        os.makedirs(out_dir)

    # plots watermark configuration
    latex = TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.035)
    latex.SetTextAlign(13)  # align at top
    latex.SetTextFont(42)
    watermark = cfg['output']['plots']['watermark']

    # open output file
    out_file = TFile(out_dir + out_label + '.root', 'recreate')

    # compute efficiencies
    efficiencies = {}
    is_retrieved_pt_interval = False
    for (key, h_rec), (_, h_gen) in zip(dic_rec.items(), dic_gen.items()):
        if h_rec is not None and h_gen is not None:
            # retrieve pt_min and pt_max from imported histogram
            if not is_retrieved_pt_interval:
                pt_min, pt_max = h_rec.GetXaxis().GetXmin(), h_rec.GetXaxis().GetXmax()
                is_retrieved_pt_interval = True

            # compute efficiency, store it in TEfficiency instance and configure it
            efficiencies[key] = compute_efficiency(
                h_rec, h_gen, rapidity_cut, axis_rapidity)
            efficiencies[key].SetName('TEfficiency_' + out_label + key)
            efficiencies[key].SetTitle(title)
            __set_object_style(efficiencies[key], key)

            # save TEfficiency instance in output file, if enabled
            if save_tefficiency:
                efficiencies[key].Write()

            # plot efficiency on canvas
            c_eff = configure_canvas(f'c{out_label + key}',
                                    pt_min, min_y_axis, pt_max, max_y_axis, title, log_y_axis)
            efficiencies[key].Draw('same') # to draw TEfficiency in TCanvas
            latex.DrawLatex(0.18, 0.92, watermark)
            c_eff.Write()
            # save canvas in separate file, if enabled
            if save_tcanvas_individual:
                save_canvas(c_eff, out_dir, out_label + key, extension)

            # convert TEfficiency instance to TH1 histogram in output file, if save enabled
            if save_th1:
                h_eff = get_th1_from_tefficiency(
                    efficiencies[key], h_rec.Clone('h' + out_label + key))
                h_eff.SetTitle(title)
                __set_object_style(h_eff, key)
                h_eff.Write()

        elif h_rec != h_gen: # one histogram name is None, the other is not
            print(f'\033[91mERROR: efficiency for {key} could not be computed,' \
                ' one of the histogram names is None!\033[0m')
            sys.exit()
        else: # both histogram names are None
            efficiencies[key] = None
            print(f'\033[94mEfficiency for {key} not computed.\033[0m')

    # overlap plots, if enabled
    if overlap is not None:
        c_overlap = configure_canvas(
            f"cOverlap{''.join(overlap)}",
            pt_min,
            min_y_axis,
            pt_max,
            max_y_axis,
            title,
            log_y_axis)

        leg = TLegend(0.6, 0.2, 0.8, 0.4)
        leg.SetTextSize(0.045)
        leg.SetFillStyle(0)

        for key, teff in efficiencies.items():
            if key in overlap:
                if teff is None:
                    print(f'\033[91mERROR: efficiency for {key} could not be computed,' \
                    ' cannot overlap it!\033[0m')
                    sys.exit()

                leg.AddEntry(teff, key, 'p')
                teff.Draw('same')

        leg.Draw()
        latex.DrawLatex(0.18, 0.92, watermark)

        if save_tcanvas_overlap:
            save_canvas(c_overlap, out_dir, out_label + 'Overlap', extension)
        c_overlap.Write()

    # close output file
    out_file.Close()

    if not batch:
        input('Press enter to exit')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Arguments')
    parser.add_argument('config', metavar='text', default='config.yml',
                        help='config file name for ml')
    parser.add_argument('--batch', help='suppress video output', action='store_true')
    args = parser.parse_args()

    print('Loading analysis configuration: ...', end='\r')
    with open(args.config, "r", encoding="utf-8") as yml_cfg:
        configuration = yaml.load(yml_cfg, yaml.FullLoader)
    print('Loading analysis configuration: Done!')

    main(configuration, args.batch)
