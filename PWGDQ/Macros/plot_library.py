"""
\author Luca Micheletti <luca.micheletti@cern.ch>, CERN
"""
import os
from os import path

import ROOT
from ROOT import (
    TCanvas,
    gStyle
)

def SetLatex(latex):
    latex.SetTextSize(0.035)
    latex.SetNDC()
    latex.SetTextFont(42)

def SetLegend(legend):
    legend.SetBorderSize(0)
    legend.SetFillColor(10)
    legend.SetFillStyle(1)
    legend.SetLineStyle(0)
    legend.SetLineColor(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.03)

def LoadStyle():
    gStyle.SetPadLeftMargin(0.15)
    gStyle.SetPadBottomMargin(0.15)
    gStyle.SetPadTopMargin(0.05)
    gStyle.SetPadRightMargin(0.05)
    gStyle.SetEndErrorSize(0.0)
    gStyle.SetTitleSize(0.05,"X")
    gStyle.SetTitleSize(0.045,"Y")
    gStyle.SetLabelSize(0.045,"X")
    gStyle.SetLabelSize(0.045,"Y")
    gStyle.SetTitleOffset(1.2,"X")
    gStyle.SetTitleOffset(1.35,"Y")

def DoResidualPlot(rooPlot, rooVar, trialName):
    rooHistResidual = rooPlot.residHist()
    rooPlotResidual = rooVar.frame(ROOT.RooFit.Title("Residual Distribution"))
    rooPlotResidual.addPlotable(rooHistResidual,"P")
    canvasResidual = TCanvas("residual_plot_{}".format(trialName), "residual_plot_{}".format(trialName), 600, 600)
    canvasResidual.SetLeftMargin(0.15)
    rooPlotResidual.GetYaxis().SetTitleOffset(1.4)
    rooPlotResidual.Draw()
    return canvasResidual

