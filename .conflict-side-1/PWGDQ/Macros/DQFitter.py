"""
\author Luca Micheletti <luca.micheletti@cern.ch>, CERN
\author Victor Valencia <valencia@subatech.in2p3.fr>, subatech
"""
from os.path import exists
from plot_library import DoResidualPlot, DoPullPlot, DoCorrMatPlot, DoPropagandaPlot

import ROOT
from ROOT import (
    TH1F,
    RooArgSet,
    RooDataHist,
    RooDataSet,
    RooRealVar,
    RooWorkspace,
    TCanvas,
    TFile,
    gPad,
    gROOT,
)


class DQFitter:
    def __init__(self, fInName, fInputName, fOutPath, fMinDataRange, fMaxDataRange):
        if not exists(fInName):
            print("The input file does not exist, exit...")
            exit()
        self.fPdfDict = {}
        self.fOutPath = fOutPath
        self.fFileOut = TFile("{}{}.root".format(fOutPath, fInputName), "RECREATE")
        self.fInputName = fInputName
        self.fFileIn = TFile.Open(fInName)
        self.fInput = self.fFileIn.Get(fInputName)
        self.fRooWorkspace = RooWorkspace("w", "workspace")
        self.fParNames = []
        self.fFitRangeMin = []
        self.fFitRangeMax = []
        self.fTrialName = ""
        self.fRooMass = RooRealVar("m", "#it{M} (GeV/#it{c}^{2})", fMinDataRange, fMaxDataRange)

    def SetFitConfig(self, pdfDict):
        """
        Method to set the fit PDFs
        """
        self.fPdfDict = pdfDict
        self.fFitRangeMin = pdfDict["fitRangeMin"]
        self.fFitRangeMax = pdfDict["fitRangeMax"]

        pdfList = []
        for pdf in self.fPdfDict["pdf"][:-1]:
            self.fTrialName = self.fTrialName + pdf + "_"

        for i in range(0, len(self.fPdfDict["pdf"])):
            if not self.fPdfDict["pdf"][i] == "SUM":
                gROOT.ProcessLineSync(
                    ".x fit_library/{}Pdf.cxx+".format(self.fPdfDict["pdf"][i])
                )

        for i in range(0, len(self.fPdfDict["pdf"])):
            parVal = self.fPdfDict["parVal"][i]
            parLimMin = self.fPdfDict["parLimMin"][i]
            parLimMax = self.fPdfDict["parLimMax"][i]
            parName = self.fPdfDict["parName"][i]

            if not self.fPdfDict["pdf"][i] == "SUM":
                for j in range(0, len(parVal)):
                    if ("sum" in parName[j]) or ("prod" in parName[j]):
                        self.fRooWorkspace.factory("{}".format(parName[j]))
                        r1 = parName[j].find("::") + 2
                        r2 = parName[j].find("(", r1)
                        parName[j] = parName[j][r1:r2]
                        if (parLimMin == parLimMax):
                            self.fRooWorkspace.factory("{}[{}]".format(parName[j], parVal[j]))
                    else:
                        if (parLimMin == parLimMax):
                            self.fRooWorkspace.factory("{}[{}]".format(parName[j], parVal[j]))
                        else:
                            self.fRooWorkspace.factory("{}[{},{},{}]".format(parName[j], parVal[j], parLimMin[j], parLimMax[j]))
                        self.fParNames.append(parName[j])
                nameFunc = self.fPdfDict["pdf"][i]
                nameFunc += "Pdf::{}Pdf(m[{},{}]".format(self.fPdfDict["pdfName"][i],self.fPdfDict["fitRangeMin"][0],self.fPdfDict["fitRangeMax"][0])
                pdfList.append(self.fPdfDict["pdfName"][i])
                for j in range(0, len(parVal)):
                    nameFunc += ",{}".format(parName[j])
                nameFunc += ")"
                self.fRooWorkspace.factory(nameFunc)
            else:
                nameFunc = self.fPdfDict["pdf"][i]
                nameFunc += "::sum("
                for j in range(0, len(pdfList)):
                    nameFunc += "{}[{},{},{}]*{}Pdf".format(parName[j], parVal[j], parLimMin[j], parLimMax[j], pdfList[j])
                    self.fParNames.append(parName[j])
                    if not j == len(pdfList) - 1:
                        nameFunc += ","
                nameFunc += ")"
                self.fRooWorkspace.factory(nameFunc)



    def FitInvMassSpectrum(self, fitRangeMin, fitRangeMax):
        """
        Method to perform binned / unbinned fit to a ROOT histogram / tree
        """
        trialName = self.fTrialName + "_" + str(fitRangeMin) + "_" + str(fitRangeMax)
        self.fRooWorkspace.Print()
        pdf = self.fRooWorkspace.pdf("sum")
        self.fRooMass.setRange("range", fitRangeMin, fitRangeMax)
        if "TTree" in self.fInput.ClassName():
            print("Perform unbinned fit")
            rooDs = RooDataSet(
                "data",
                "data",
                RooArgSet(self.fRooMass),
                ROOT.RooFit.Import(self.fInput),
            )
        else:
            print("Perform binned fit")
            rooDs = RooDataHist(
                "data",
                "data",
                RooArgSet(self.fRooMass),
                ROOT.RooFit.Import(self.fInput),
            )
        rooFitRes = ROOT.RooFitResult(pdf.fitTo(rooDs, ROOT.RooFit.Save(True)))
        fRooPlot = self.fRooMass.frame(ROOT.RooFit.Title(trialName))
        fRooPlotCopy = self.fRooMass.frame(ROOT.RooFit.Title(trialName))
        rooDs.plotOn(fRooPlot, ROOT.RooFit.MarkerStyle(20), ROOT.RooFit.MarkerSize(0.6), ROOT.RooFit.Range(fitRangeMin, fitRangeMax))
        pdf.plotOn(fRooPlot, ROOT.RooFit.LineColor(ROOT.kRed+1), ROOT.RooFit.LineWidth(2), ROOT.RooFit.Range(fitRangeMin, fitRangeMax))

        index = 1
        histResults = TH1F(
            "fit_results_{}".format(trialName),
            "fit_results_{}".format(trialName),
            len(self.fParNames),
            0.0,
            len(self.fParNames),
        )
        for parName in self.fParNames:
            histResults.GetXaxis().SetBinLabel(index, parName)
            histResults.SetBinContent(index, self.fRooWorkspace.var(parName).getVal())
            histResults.SetBinContent(index, self.fRooWorkspace.var(parName).getError())
            index += 1

        for i in range(0, len(self.fPdfDict["pdf"])):
            if not self.fPdfDict["pdfName"][i] == "SUM":
                pdf.plotOn(fRooPlot, ROOT.RooFit.Components("{}Pdf".format(self.fPdfDict["pdfName"][i])), ROOT.RooFit.LineColor(self.fPdfDict["pdfColor"][i]), ROOT.RooFit.LineStyle(self.fPdfDict["pdfStyle"][i]), ROOT.RooFit.LineWidth(2), ROOT.RooFit.Range(fitRangeMin, fitRangeMax))

        extraText = []
        paveText = ROOT.TPaveText(0.60, 0.45, 0.99, 0.94, "brNDC")
        paveText.SetTextFont(42)
        paveText.SetTextSize(0.025)
        paveText.SetFillColor(ROOT.kWhite)
        for parName in self.fParNames:
            paveText.AddText("{} = {:.4f} #pm {:.4f}".format(parName, self.fRooWorkspace.var(parName).getVal(), self.fRooWorkspace.var(parName).getError()))
            if self.fPdfDict["parForPropagandaPlot"].count(parName) > 0:
                text = self.fPdfDict["parNameForPropagandaPlot"][self.fPdfDict["parForPropagandaPlot"].index(parName)]
                if "sig" in parName:
                    extraText.append("{} = {:.0f} #pm {:.0f}".format(text, self.fRooWorkspace.var(parName).getVal(), self.fRooWorkspace.var(parName).getError()))
                else:
                    extraText.append("{} = {:.3f} #pm {:.3f}".format(text, self.fRooWorkspace.var(parName).getVal(), self.fRooWorkspace.var(parName).getError()))
            for i in range(0, len(self.fPdfDict["pdfName"])):
                if self.fPdfDict["pdfName"][i] in parName:
                    (paveText.GetListOfLines().Last()).SetTextColor(self.fPdfDict["pdfColor"][i])

        nPars = rooFitRes.floatParsFinal().getSize()
        if "TTree" in self.fInput.ClassName():
            # Convert RooDataSet into RooDataHist to extract the Chi2 value
            rooDh = RooDataHist("rooDh", "binned version of rooDs", RooArgSet(self.fRooMass), rooDs)
            chi2 = ROOT.RooChi2Var("chi2", "chi2", pdf, rooDh)
            nbinsperGev = rooDh.numEntries() / (self.fPdfDict["fitRangeMax"][0] - self.fPdfDict["fitRangeMin"][0])
            nBins = (fitRangeMax - fitRangeMin) * nbinsperGev
            ndof = nBins - nPars
            reduced_chi2 = chi2.getVal() / ndof
            paveText.AddText("#bf{#chi^{2}/dof = %3.2f}" % (reduced_chi2))
            extraText.append("#chi^{2}/dof = %3.2f" % reduced_chi2)
        else:
            # To Do : Find a way to get the number of bins differently. The following is a temparary solution.
            # WARNING : The largest fit range has to come first in the config file otherwise it does not work
            chi2 = ROOT.RooChi2Var("chi2", "chi2", pdf, rooDs)
            nbinsperGev = rooDs.numEntries() / (self.fPdfDict["fitRangeMax"][0] - self.fPdfDict["fitRangeMin"][0])
            nBins = (fitRangeMax - fitRangeMin) * nbinsperGev
            ndof = nBins - nPars
            reduced_chi2 = chi2.getVal() / ndof
            paveText.AddText("n Par = %3.2f" % (nPars))
            paveText.AddText("n Bins = %3.2f" % (nBins))
            paveText.AddText("#bf{#chi^{2}/dof = %3.2f}" % reduced_chi2)
            fRooPlot.addObject(paveText)
            extraText.append("#chi^{2}/dof = %3.2f" % reduced_chi2)

        # Fit plot
        canvasFit = TCanvas(
            "fit_plot_{}".format(trialName), "fit_plot_{}".format(trialName), 600, 600
        )
        canvasFit.SetLeftMargin(0.15)
        gPad.SetLeftMargin(0.15)
        fRooPlot.GetYaxis().SetTitleOffset(1.4)
        fRooPlot.Draw()

        # Residual plot
        canvasResidual = DoResidualPlot(fRooPlot, self.fRooMass, trialName)

        # Pull plot
        canvasPull = DoPullPlot(fRooPlot, self.fRooMass, trialName)

        # Correlation matrix plot
        canvasCorrMat = DoCorrMatPlot(rooFitRes, trialName)

        # Propaganda plot
        if self.fPdfDict["doPropagandaPlot"]:
            DoPropagandaPlot(rooDs, pdf, fRooPlotCopy, self.fPdfDict, self.fInputName, trialName, self.fOutPath, extraText)

        # Save results
        self.fFileOut.cd()
        canvasFit.Write()
        canvasResidual.Write()
        canvasPull.Write()
        canvasCorrMat.Write()
        histResults.Write()

    def MultiTrial(self):
        """
        Method to run multiple fits of the same invariant mass distribution
        """
        for iRange in range(0, len(self.fFitRangeMin)):
            self.FitInvMassSpectrum(
                self.fFitRangeMin[iRange], self.fFitRangeMax[iRange]
            )
        self.fFileOut.Close()
