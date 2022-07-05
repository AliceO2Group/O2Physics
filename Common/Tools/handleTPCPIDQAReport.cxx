// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

///
/// \file handleTPCPIDQAReport.cxx
/// \author Tiantian Cheng
/// \brief exec for root file for the TPC PID and qa task

#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TLine.h"
#include "TList.h"
#include "TMath.h"
#include "TPaveText.h"
#include "TObjArray.h"
#include "TString.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TF1.h"
#include "TKey.h"
#include "TPDF.h"
#include "TColor.h"
#include <array>
#include <fstream>
#include <sstream>
#include <string>
#include "TFile.h"
#include "Algorithm/RangeTokenizer.h"
#include "handleParamBase.h"

void SetupStyle();
TH2* Get2DHistogramfromList(TDirectoryFile* pidqalist, const char* subdirname, TObject* histoname);
void AddFit(TH2* h2d);
void PublishCanvas(TDirectoryFile* qaList, const char* subdirname);
void SetupPadStyle();

bool initOptionsAndParse(bpo::options_description& options, int argc, char* argv[])
{
  options.add_options()("inputFileName", bpo::value<std::string>()->default_value("AnalysisResults.root"), "The name of input file")("outputFileName", bpo::value<std::string>()->default_value("TPCPIDQA.pdf"), "The name of output file")("help,h", "Produce help message");
  try {
    bpo::store(parse_command_line(argc, argv, options), arguments);
    // help
    if (arguments.count("help")) {
      LOG(info) << options;
      return false;
    }
    bpo::notify(arguments);
  } catch (const bpo::error& e) {
    LOG(error) << e.what() << "\n";
    LOG(error) << "Error parsing command line arguments; Available options:";
    LOG(error) << options;
    return false;
  }
  return true;
} // initOptionsAndParse

TCanvas* fCanvas = 0x0;

int main(int argc, char* argv[])
{

  bpo::options_description options("Allowed options");
  if (!initOptionsAndParse(options, argc, argv)) {
    return 1;
  }

  const std::string inputFile = arguments["inputFileName"].as<std::string>();
  const std::string outputFile = arguments["outputFileName"].as<std::string>();

  SetupStyle();

  TFile f(inputFile.data(), "READ");
  if (!f.IsOpen()) {
    LOG(error) << "Input file " << inputFile << "could not be read";
    return 1;
  }

  fCanvas = new TCanvas;
  TPDF p(outputFile.c_str());

  std::string dirname;
  TDirectoryFile* qaList = nullptr;

  //-- tpc-pid-full-qa
  dirname = "tpc-pid-full-qa";
  qaList = (TDirectoryFile*)f.Get(dirname.c_str());
  if (!qaList) {
    printf("Could not find directory '%s' in file '%s' \n", dirname.c_str(), f.GetName());
    return false;
  }
  PublishCanvas(qaList, "nsigma");
  delete qaList;

  //--- tpc-pid-qa
  dirname = "tpc-pid-qa";
  qaList = (TDirectoryFile*)f.Get(dirname.c_str());
  if (!qaList) {
    printf("Could not find directory '%s' in file '%s' \n", dirname.c_str(), f.GetName());
    return false;
  }
  PublishCanvas(qaList, "nsigma");
  delete qaList;

  //--- qa-tpc-tof
  dirname = "qa-tpc-tof";
  qaList = (TDirectoryFile*)f.Get(dirname.c_str());
  if (!qaList) {
    printf("Could not find directory '%s' in file '%s' \n", dirname.c_str(), f.GetName());
    return false;
  }
  PublishCanvas(qaList, "nsigmaTPCTOF");
  delete qaList;

  //--- qa-tpc-v0
  dirname = "qa-tpc-v0";
  qaList = (TDirectoryFile*)f.Get(dirname.c_str());
  if (!qaList) {
    printf("Could not find directory '%s' in file '%s' \n", dirname.c_str(), f.GetName());
    return false;
  }
  PublishCanvas(qaList, "nsigmaTPCV0");
  PublishCanvas(qaList, "nsigmaTPCV0VsEta");
  delete qaList;

  p.Close();
  delete fCanvas;

} // main()

void SetupStyle()
{
  const Int_t NCont = 255;

  TStyle* st = new TStyle("mystyle", "mystyle");
  gROOT->GetStyle("Plain")->Copy((*st));
  st->SetTitleX(0.1);
  st->SetTitleW(0.8);
  st->SetTitleH(0.08);
  st->SetStatX(.9);
  st->SetStatY(.9);
  st->SetNumberContours(NCont);
  st->SetPalette(1, 0);
  st->SetOptStat("erm");
  st->SetOptFit(0);
  st->SetGridColor(kGray + 1);
  st->SetPadGridX(kTRUE);
  st->SetPadGridY(kTRUE);
  st->SetPadTickX(kTRUE);
  st->SetPadTickY(kTRUE);
  st->cd();

  const Int_t NRGBs = 5;
  Double_t stops[NRGBs] = {0.00, 0.34, 0.61, 0.84, 1.00};
  Double_t red[NRGBs] = {0.00, 0.00, 0.87, 1.00, 0.51};
  Double_t green[NRGBs] = {0.00, 0.81, 1.00, 0.20, 0.00};
  Double_t blue[NRGBs] = {0.51, 1.00, 0.12, 0.00, 0.00};

  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);

} // SetupStyle

TH2* Get2DHistogramfromList(TDirectoryFile* pidqalist, const char* subdirname, TObject* histoname)
{

  TDirectoryFile* histolist = (TDirectoryFile*)pidqalist->Get(subdirname);
  if (!histolist) {
    printf(" directory not found \n");
    return 0x0;
  }

  TH2* histo = (TH2*)histolist->FindObject(histoname);
  if (!histo) {
    printf(" histogram not found \n");
    return 0x0;
  }

  return histo;
}

void AddFit(TH2* h2d)
{
  //
  // Fit in slices and draw mean and sigma
  //

  TF1* f1 = new TF1("f1", "gaus");
  f1->SetRange(-1.5, 1.5);
  TObjArray aSlices;

  h2d->FitSlicesY(f1, 0, -1, 0, "QNR", &aSlices);
  aSlices.SetOwner(1);

  TH1* hMean = (TH1*)aSlices.At(1);
  TH1* hSigma = (TH1*)aSlices.At(2);
  TH1* hChi2 = (TH1*)aSlices.At(3);

  hChi2->Scale(1. / 10.);
  aSlices.AddAt(0x0, 1);
  aSlices.AddAt(0x0, 2);
  aSlices.AddAt(0x0, 3);

  hMean->SetMarkerStyle(20);
  hMean->SetMarkerSize(0.3);
  hMean->SetOption("same");
  h2d->GetListOfFunctions()->Add(hMean);

  hSigma->SetMarkerStyle(20);
  hSigma->SetMarkerSize(0.3);
  hSigma->SetOption("same");
  hSigma->SetMarkerColor(kMagenta);
  h2d->GetListOfFunctions()->Add(hSigma);

  hChi2->SetOption("same");
  hChi2->SetMarkerColor(kMagenta + 2);
  hChi2->SetLineColor(kMagenta + 2);
  h2d->GetListOfFunctions()->Add(hChi2);

  TLine* l = 0x0;
  l = new TLine(h2d->GetXaxis()->GetXmin(), 0, h2d->GetXaxis()->GetXmax(), 0);
  l->SetLineStyle(2);
  h2d->GetListOfFunctions()->Add(l);
  l = new TLine(h2d->GetXaxis()->GetXmin(), 1, h2d->GetXaxis()->GetXmax(), 1);
  l->SetLineStyle(2);
  h2d->GetListOfFunctions()->Add(l);

} // AddFit
void PublishCanvas(TDirectoryFile* qaList, const char* subdirname)
{
  TObjArray arrHistos;

  TPaveText pt(.1, .1, .9, .9, "NDC");
  pt.SetBorderSize(1);
  pt.SetFillColor(0);
  pt.SetTextSizePixels(16);

  pt.AddText(Form("%s %s", qaList->GetName(), subdirname));
  arrHistos.Add(&pt);

  TDirectoryFile* QA = (TDirectoryFile*)qaList->Get(subdirname);
  TList* listry = QA->GetListOfKeys();

  TIter next(listry);
  TKey* key;
  TObject* obj;

  while ((key = (TKey*)next())) {
    obj = key->ReadObj();
    TH2* h = Get2DHistogramfromList(qaList, subdirname, obj);
    h->SetOption("colz");
    AddFit(h);
    arrHistos.Add(h);
  }

  Int_t nPads = arrHistos.GetEntriesFast();
  Int_t nCols = (Int_t)TMath::Ceil(TMath::Sqrt(nPads));
  Int_t nRows = (Int_t)TMath::Ceil((Double_t)nPads / (Double_t)nCols);
  fCanvas->Divide(nCols, nRows);

  for (Int_t i = 0; i < nPads; ++i) {
    fCanvas->cd(i + 1);
    SetupPadStyle();
    if (strcmp(subdirname, "nsigmaTPCV0VsEta") == 0)
      gPad->SetLogx(kFALSE);
    arrHistos.At(i)->Draw();
  }

  fCanvas->Update();
  fCanvas->Clear();

} // PublishCanvas
void SetupPadStyle()
{
  gPad->SetLogx();
  gPad->SetLogz();
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetTickx();
  gPad->SetTicky();

} // SetupPadStyle
