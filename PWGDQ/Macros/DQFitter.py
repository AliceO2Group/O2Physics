"""
\author Luca Micheletti <luca.micheletti@cern.ch>, CERN
"""
from os.path import exists
from traceback import print_tb
import ROOT
from ROOT import TCanvas, TFile, TH1F, RooRealVar, RooDataSet, RooWorkspace, RooDataHist, RooArgSet
from ROOT import gPad, gROOT

class DQFitter:
    def __init__(self, fInName, fInputName):
        if not exists(fInName):
            print("The input file does not exist, exit...")
            exit()
        self.fPdfDict          = {}
        self.fFileOut          = TFile("FitResults.root", "RECREATE")
        self.fFileIn           = TFile.Open(fInName)
        self.fInput            = self.fFileIn.Get(fInputName)
        self.fRooWorkspace     = RooWorkspace('w','workspace')
        self.fParNames         = []
        self.fFitRangeMin      = []
        self.fFitRangeMax      = []
        self.fTrialName        = ""
        self.fRooMass          = RooRealVar("m", "#it{M} (GeV/#it{c}^{2})", 2, 5)

    def SetFitConfig(self, pdfDict):
        '''
        Method to set the fit PDFs
        '''
        self.fPdfDict = pdfDict
        self.fFitRangeMin = pdfDict["fitRangeMin"]
        self.fFitRangeMax = pdfDict["fitRangeMax"]

        pdfList = []
        for pdf in self.fPdfDict["pdf"][:-1]:
            self.fTrialName = self.fTrialName + pdf + "_"

        for i in range(0, len(self.fPdfDict["pdf"])):
            parVal = self.fPdfDict["parVal"][i]
            parLimMin = self.fPdfDict["parLimMin"][i]
            parLimMax = self.fPdfDict["parLimMax"][i]
            parName = self.fPdfDict["parName"][i]

            if not self.fPdfDict["pdf"][i] == "SUM":
                gROOT.ProcessLineSync(".x fit_library/{}Pdf.cxx+".format(self.fPdfDict["pdf"][i]))
                nameFunc = self.fPdfDict["pdf"][i]
                nameFunc += "Pdf::{}Pdf(m[2,5]".format(self.fPdfDict["pdf"][i])
                pdfList.append(self.fPdfDict["pdf"][i])

                for j in range(0, len(parVal)):
                    nameFunc += ",{}[{},{},{}]".format(parName[j], parVal[j], parLimMin[j], parLimMax[j])
                    self.fParNames.append(parName[j])
                nameFunc += ")"
                print(nameFunc)
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
                print(nameFunc)
                self.fRooWorkspace.factory(nameFunc)

    def FitInvMassSpectrum(self, fitRangeMin, fitRangeMax):
        '''
        Method to perform binned / unbinned fit to a ROOT histogram / tree
        '''
        trialName = self.fTrialName + "_" + str(fitRangeMin) + "_" + str(fitRangeMax)
        self.fRooWorkspace.Print()
        pdf = self.fRooWorkspace.pdf("sum")
        self.fRooMass.setRange("range", fitRangeMin, fitRangeMax)
        fRooPlot = self.fRooMass.frame(ROOT.RooFit.Title(trialName))
        if "TTree" in self.fInput.ClassName():
            print("Perform unbinned fit")
            rooDs = RooDataSet("data", "data", RooArgSet(self.fRooMass), ROOT.RooFit.Import(self.fInput))
        else:
            print("Perform binned fit")
            rooDs = RooDataHist("data", "data", RooArgSet(self.fRooMass), ROOT.RooFit.Import(self.fInput))
        pdf.fitTo(rooDs)

        index = 1
        histResults = TH1F("fit_results_{}".format(trialName), "fit_results_{}".format(trialName), len(self.fParNames), 0., len(self.fParNames))
        for parName in self.fParNames:
            histResults.GetXaxis().SetBinLabel(index, parName)
            histResults.SetBinContent(index, self.fRooWorkspace.var(parName).getVal())
            histResults.SetBinContent(index, self.fRooWorkspace.var(parName).getError())
            index += 1

        rooDs.plotOn(fRooPlot)
        pdf.plotOn(fRooPlot)
        pdf.paramOn(fRooPlot, ROOT.RooFit.Layout(0.55))

        # Fit plot
        canvasFit = TCanvas("fit_plot_{}".format(trialName), "fit_plot_{}".format(trialName), 600, 600)
        canvasFit.SetLeftMargin(0.15)
        gPad.SetLeftMargin(0.15)
        fRooPlot.GetYaxis().SetTitleOffset(1.4)
        fRooPlot.Draw()

        # Ratio plot
        rooHistRatio = fRooPlot.residHist()
        rooPlotRatio = self.fRooMass.frame(ROOT.RooFit.Title("Residual Distribution"))
        rooPlotRatio.addPlotable(rooHistRatio,"P")
        canvasRatio = TCanvas("ratio_plot_{}".format(trialName), "ratio_plot_{}".format(trialName), 600, 600)
        canvasRatio.SetLeftMargin(0.15)
        rooPlotRatio.GetYaxis().SetTitleOffset(1.4)
        rooPlotRatio.Draw()

        # Pull plot
        rooHistPull = fRooPlot.pullHist()
        rooPlotPull = self.fRooMass.frame(ROOT.RooFit.Title("Pull Distribution"))
        rooPlotPull.addPlotable(rooHistPull,"P")
        canvasPull = TCanvas("pull_plot_{}".format(trialName), "pull_plot_{}".format(trialName), 600, 600)
        canvasPull.SetLeftMargin(0.15)
        rooPlotPull.GetYaxis().SetTitleOffset(1.4)
        rooPlotPull.Draw()

        # Save results
        self.fFileOut.cd()
        canvasFit.Write()
        canvasRatio.Write()
        canvasPull.Write()
        histResults.Write()


    def MultiTrial(self):
        '''
        Method to run multiple fits of the same invariant mass distribution
        '''
        for iRange in range(0, len(self.fFitRangeMin)):
            self.FitInvMassSpectrum(self.fFitRangeMin[iRange], self.fFitRangeMax[iRange])
        self.fFileOut.Close()