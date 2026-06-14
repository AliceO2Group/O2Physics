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
//
// This code calculates a centrality calibration based on output
// from the multiplicityQa task.
//
// Comments, suggestions, questions? Please write to:
// - victor.gonzalez@cern.ch
// - david.dobrigkeit.chinellato@cern.ch
//
#include "multCalibrator.h"

#include <TDirectory.h>
#include <TFile.h>
#include <TH1.h>
#include <TList.h>
#include <TMathBase.h>
#include <TNamed.h>
#include <TString.h>

#include <RtypesCore.h>

#include <cmath>
#include <iostream> // FIXME

using namespace std;

const TString multCalibrator::fCentEstimName[kNCentEstim] = {
  "RawV0M", "RawT0M", "RawFDD", "RawNTracks",
  "ZeqV0M", "ZeqT0M", "ZeqFDD", "ZeqNTracks"};

multCalibrator::multCalibrator() : TNamed(),
                                   lDesiredBoundaries(0),
                                   lNDesiredBoundaries(0),
                                   fkPrecisionWarningThreshold(1000.0),
                                   fInputFileName("AnalysisResults.root"),
                                   fOutputFileName("CCDB-objects.root"),
                                   fAnchorPointValue(-1),
                                   fAnchorPointPercentage(100),
                                   fCalibHists(0x0),
                                   fPrecisionHistogram(0x0)
{
  // Constructor
  // Make sure the TList owns its objects
  fCalibHists = new TList();
  fCalibHists->SetOwner(kTRUE);
}

multCalibrator::multCalibrator(const char* name, const char* title) : TNamed(name, title),
                                                                      lDesiredBoundaries(0),
                                                                      lNDesiredBoundaries(0),
                                                                      fkPrecisionWarningThreshold(1000.0),
                                                                      fInputFileName("AnalysisResults.root"),
                                                                      fOutputFileName("CCDB-objects.root"),
                                                                      fAnchorPointValue(-1),
                                                                      fAnchorPointPercentage(90),
                                                                      fCalibHists(0x0),
                                                                      fPrecisionHistogram(0x0)
{
  // Named Constructor
  // Make sure the TList owns its objects
  fCalibHists = new TList();
  fCalibHists->SetOwner(kTRUE);
}
//________________________________________________________________
multCalibrator::~multCalibrator()
{
  // Destructor
  if (fCalibHists) {
    delete fCalibHists;
    fCalibHists = 0x0;
  }
  if (fPrecisionHistogram) {
    delete fPrecisionHistogram;
    fPrecisionHistogram = 0x0;
  }
}

//________________________________________________________________
Bool_t multCalibrator::Calibrate()
{
  // Function meant to generate calibration OADB
  //
  // --- input : fInputFileName, containing a TTree object
  // --- output: fOutputFileName, containing OABD object
  //

  cout << "=== STARTING CALIBRATION PROCEDURE ===" << endl;
  cout << " * Input File.....: " << fInputFileName.Data() << endl;
  cout << " * Output File....: " << fOutputFileName.Data() << endl;
  cout << endl;

  TFile* fileInput = new TFile(fInputFileName.Data(), "READ");
  if (!fileInput) {
    cout << "Input file not found!" << endl;
    return kFALSE;
  }

  // Step 1: verify if input file contains desired histograms
  TH1D* hRaw[kNCentEstim];
  for (Int_t iv = 0; iv < kNCentEstim; iv++) {
    hRaw[iv] = reinterpret_cast<TH1D*>(fileInput->Get(Form("multiplicity-qa/multiplicityQa/h%s", fCentEstimName[iv].Data())));
    if (!hRaw[iv]) {
      cout << Form("File does not contain histogram h%s, which is necessary for calibration!", fCentEstimName[iv].Data()) << endl;
      return kFALSE;
    }
  }

  cout << "Histograms loaded! Will now calibrate..." << endl;

  // Create output file
  TFile* fOut = new TFile(fOutputFileName.Data(), "RECREATE");
  TH1F* hCalib[kNCentEstim];
  for (Int_t iv = 0; iv < kNCentEstim; iv++) {
    cout << Form("Calibrating estimator: %s", fCentEstimName[iv].Data()) << endl;
    hCalib[iv] = GetCalibrationHistogram(hRaw[iv], Form("hCalib%s", fCentEstimName[iv].Data()));
    hCalib[iv]->Write();
  }

  cout << "Saving calibration file..." << endl;
  fOut->Write();
  cout << "Done! Enjoy!" << endl;
  return kTRUE;
}

Double_t multCalibrator::GetRawMax(TH1* histo)
{
  // This function gets the max X value (right edge) which is filled.
  for (Int_t ii = histo->GetNbinsX(); ii > 0; ii--) {
    if (histo->GetBinContent(ii) < 1e-10)
      return histo->GetBinLowEdge(ii + 1);
  }
  return 1e+6;
}

Double_t multCalibrator::GetBoundaryForPercentile(TH1* histo, Double_t lPercentileRequested, Double_t& lPrecisionEstimate)
{
  // This function returns the boundary for a specific percentile.
  // It uses a linear interpolation in an attempt to get more precision
  // than the binning of the histogram used for quantiling.
  //
  // It also estimates a certain level of precision of the procedure
  // by explicitly comparing the bin content of the bins around the boundary
  // with the entire cross section, effectively reporting back a percentage
  // that corresponds to those bins. If this percentage is O(percentile bin
  // width requested), then the user should worry and we print out a warning.

  const Double_t lPrecisionConstant = 2.0;

  Double_t lRawMax = GetRawMax(histo);

  if (lPercentileRequested < 1e-7)
    return lRawMax; // safeguard
  if (lPercentileRequested > 100 - 1e-7)
    return 0.0; // safeguard

  Double_t lReturnValue = 0.0;
  Double_t lPercentile = 100.0 - lPercentileRequested;
  Double_t lPercentileAnchor = 100.0 - fAnchorPointPercentage;
  lPrecisionEstimate = -1;
  if (lPercentile < lPercentileAnchor + 1e-7)
    return fAnchorPointValue;
  Double_t lCount = 0;

  // Anchor point changes: if anchored, start at the first bin that includes that
  Long_t lFirstBin = 1;
  Double_t lHadronicTotal = histo->Integral(1, histo->GetNbinsX()); // histo->GetEntries();
  if (fAnchorPointValue > 0) {
    lFirstBin = histo->FindBin(fAnchorPointValue + 1e-6);
    Double_t lAbove = histo->Integral(lFirstBin, histo->GetNbinsX());
    lHadronicTotal = lAbove * 100.0 / (fAnchorPointPercentage);
    lCount = lHadronicTotal - lAbove; // the relevant anchored-out part
  }

  const Long_t lNBins = histo->GetNbinsX();
  Double_t lCountDesired = lPercentile * lHadronicTotal / 100;

  for (Long_t ibin = lFirstBin; ibin < lNBins; ibin++) {
    lCount += histo->GetBinContent(ibin);
    if (lCount >= lCountDesired) {
      // Found bin I am looking for!
      Double_t lWidth = histo->GetBinWidth(ibin);
      Double_t lLeftPercentile = 100. * (lCount - histo->GetBinContent(ibin)) / lHadronicTotal;
      Double_t lRightPercentile = 100. * lCount / lHadronicTotal;
      lPrecisionEstimate = (lRightPercentile - lLeftPercentile) / lPrecisionConstant;

      Double_t lProportion = (lPercentile - lLeftPercentile) / (lRightPercentile - lLeftPercentile);

      lReturnValue = histo->GetBinLowEdge(ibin) + lProportion * lWidth;
      break;
    }
  }
  return lReturnValue;
}

//________________________________________________________________
void multCalibrator::SetStandardAdaptiveBoundaries()
{
  // Function to set standard adaptive boundaries
  // Typically used in pp, goes to 0.001% binning for highest multiplicity
  lNDesiredBoundaries = 0;
  lDesiredBoundaries = new Double_t[1100];
  lDesiredBoundaries[0] = 100;
  // From Low To High Multiplicity
  for (Int_t ib = 1; ib < 91; ib++) {
    lNDesiredBoundaries++;
    lDesiredBoundaries[lNDesiredBoundaries] = lDesiredBoundaries[lNDesiredBoundaries - 1] - 1.0;
  }
  for (Int_t ib = 1; ib < 91; ib++) {
    lNDesiredBoundaries++;
    lDesiredBoundaries[lNDesiredBoundaries] = lDesiredBoundaries[lNDesiredBoundaries - 1] - 0.1;
  }
  for (Int_t ib = 1; ib < 91; ib++) {
    lNDesiredBoundaries++;
    lDesiredBoundaries[lNDesiredBoundaries] = lDesiredBoundaries[lNDesiredBoundaries - 1] - 0.01;
  }
  for (Int_t ib = 1; ib < 101; ib++) {
    lNDesiredBoundaries++;
    lDesiredBoundaries[lNDesiredBoundaries] = lDesiredBoundaries[lNDesiredBoundaries - 1] - 0.001;
  }
  lNDesiredBoundaries++;
  cout << "Set standard adaptive percentile boundaries! Nboundaries: " << lNDesiredBoundaries << endl;
}

//________________________________________________________________
void multCalibrator::SetRun3AdaptiveBoundaries()
{
  // Function to set standard adaptive boundaries
  // Run 3 exclusive: goes to 0.00001% binning for highest multiplicity
  // Warning: this setting requires well filled input histograms
  //
  // -> at 0.00001%, one needs 10 million events to get a single one in the highest
  // category. In Run 3, at 500 kHz, this corresponds to 20 seconds of dataking;
  // any run lasting over 5 minutes should have enough to calibrate to the highest
  // boundary.
  //
  // QA when using hyperfine binning is still advised:
  // Selecting very fine percentages may select on spurious situations
  // such as beam-background, pileup, etc.

  lNDesiredBoundaries = 0;
  lDesiredBoundaries = new Double_t[1100];
  lDesiredBoundaries[0] = 100;
  // From Low To High Multiplicity
  for (Int_t ib = 1; ib < 91; ib++) {
    lNDesiredBoundaries++;
    lDesiredBoundaries[lNDesiredBoundaries] = lDesiredBoundaries[lNDesiredBoundaries - 1] - 1.0;
  }
  for (Int_t ib = 1; ib < 91; ib++) {
    lNDesiredBoundaries++;
    lDesiredBoundaries[lNDesiredBoundaries] = lDesiredBoundaries[lNDesiredBoundaries - 1] - 0.1;
  }
  for (Int_t ib = 1; ib < 91; ib++) {
    lNDesiredBoundaries++;
    lDesiredBoundaries[lNDesiredBoundaries] = lDesiredBoundaries[lNDesiredBoundaries - 1] - 0.01;
  }
  for (Int_t ib = 1; ib < 91; ib++) {
    lNDesiredBoundaries++;
    lDesiredBoundaries[lNDesiredBoundaries] = lDesiredBoundaries[lNDesiredBoundaries - 1] - 0.001;
  }
  for (Int_t ib = 1; ib < 91; ib++) {
    lNDesiredBoundaries++;
    lDesiredBoundaries[lNDesiredBoundaries] = lDesiredBoundaries[lNDesiredBoundaries - 1] - 0.0001;
  }
  for (Int_t ib = 1; ib < 101; ib++) {
    lNDesiredBoundaries++;
    lDesiredBoundaries[lNDesiredBoundaries] = lDesiredBoundaries[lNDesiredBoundaries - 1] - 0.00001;
  }
  lNDesiredBoundaries++;
  cout << "Set run 3 adaptive percentile boundaries! Nboundaries: " << lNDesiredBoundaries << endl;
}

//________________________________________________________________
void multCalibrator::SetStandardOnePercentBoundaries()
{
  // Function to set standard adaptive boundaries
  // Typically used in pp, goes to 0.001% binning for highest multiplicity
  lNDesiredBoundaries = 101;
  lDesiredBoundaries = new Double_t[101];
  lDesiredBoundaries[0] = 100;
  // From Low To High Multiplicity
  for (Int_t ib = 1; ib < 101; ib++)
    lDesiredBoundaries[ib] = lDesiredBoundaries[ib - 1] - 1.0;
  cout << "Set standard 1%-wide percentile boundaries! Nboundaries: " << lNDesiredBoundaries << endl;
}

//________________________________________________________________
bool multCalibrator::IsBinningSane(TH1* histogram)
{
  // check if binning that will be attempted is too fine for this histogram
  double minimumBinWidth = 1000.0f;
  for (int ib = 0; ib < lNDesiredBoundaries - 1; ib++) {
    if (std::abs(lDesiredBoundaries[ib + 1] - lDesiredBoundaries[ib]) < minimumBinWidth) {
      minimumBinWidth = std::abs(lDesiredBoundaries[ib + 1] - lDesiredBoundaries[ib]);
    }
  }
  if (minimumBinWidth < 1e-9) {
    cout << "Excessively fine binning requested: minimum bin width registers as " << minimumBinWidth << endl;
    return false; // not reasonable
  }
  if (histogram->GetEntries() < 1000.0 / minimumBinWidth) {
    cout << "Histogram " << histogram->GetName() << " does not have enough entries (" << histogram->GetEntries() << ") to calibrate!" << endl;
    return false;
  }
  if (std::abs(histogram->GetMean()) < 1e-9) {
    cout << "Histogram " << histogram->GetName() << " has suspicious mean: " << histogram->GetMean() << endl;
    return false;
  }
  cout << "Sanity check: " << histogram->GetName() << ", entries = " << histogram->GetEntries() << ", 1/(min bin width) = " << 100.0 / minimumBinWidth << " -> OK ! " << endl;
  return true;
}

//________________________________________________________________
TH1F* multCalibrator::GetCalibrationHistogram(TH1* histoRaw, TString lHistoName)
{
  // This function returns a calibration histogram
  //(pp or p-Pb like, no anchor point considered)

  // Reset + recreate precision histogram
  ResetPrecisionHistogram();

  // Consistency check
  if (fAnchorPointValue > 0 && TMath::Abs(fAnchorPointPercentage - lDesiredBoundaries[0]) > 1e-6) {
    cout << "PROBLEM WITH ANCHOR POINT SETTINGS! " << endl;
    cout << "Anchor point percentage requested: " << fAnchorPointPercentage << endl;
    cout << "Last boundary: " << lDesiredBoundaries[0] << endl;
  }

  // binning check
  if (!IsBinningSane(histoRaw)) {
    cout << "Requested binning is not viable for input histogram named " << histoRaw->GetName() << "! Will return 100.5 for all requests." << endl;
    Double_t lDummyBounds[2];
    lDummyBounds[0] = 0;
    lDummyBounds[1] = 1.0;
    TH1F* hCalib = new TH1F(lHistoName.Data(), "", 1, lDummyBounds);
    hCalib->SetBinContent(0, 100.5);
    hCalib->SetBinContent(1, 100.5);
    hCalib->SetBinContent(2, 100.5);
    return hCalib;
  }

  // Aux vars
  Double_t lMiddleOfBins[1000];
  for (Long_t lB = 1; lB < lNDesiredBoundaries; lB++) {
    // place squarely at the middle to ensure it's all fine
    lMiddleOfBins[lB - 1] = 0.5 * (lDesiredBoundaries[lB] + lDesiredBoundaries[lB - 1]);
  }
  Double_t lBounds[lNDesiredBoundaries + 1];
  Double_t lPrecision[lNDesiredBoundaries + 1];

  if (fAnchorPointValue > 0) {
    lBounds[0] = 0;
  }

  for (Int_t ii = 0; ii < lNDesiredBoundaries; ii++) {
    Int_t lDisplacedii = ii;
    if (fAnchorPointValue > 0)
      lDisplacedii++;
    lBounds[lDisplacedii] = GetBoundaryForPercentile(histoRaw, lDesiredBoundaries[ii], lPrecision[ii]);
    bool warnUser = false;
    TString lPrecisionString = "(Precision OK)";
    if (ii != 0 && ii != lNDesiredBoundaries - 1) {
      // check precision, please
      if ((lPrecision[ii] / TMath::Abs(lDesiredBoundaries[ii + 1] - lDesiredBoundaries[ii])) > fkPrecisionWarningThreshold) {
        lPrecisionString = "(WARNING: BINNING MAY LEAD TO IMPRECISION!)";
        warnUser = true;
      }
      if ((lPrecision[ii] / TMath::Abs(lDesiredBoundaries[ii - 1] - lDesiredBoundaries[ii])) > fkPrecisionWarningThreshold) {
        lPrecisionString = "(WARNING: BINNING MAY LEAD TO IMPRECISION!)";
        warnUser = true;
      }
    }
    if (warnUser) {
      cout << histoRaw->GetName() << " boundaries, percentile: " << lDesiredBoundaries[ii] << "%\t Signal value = " << lBounds[lDisplacedii] << "\tprecision = " << lPrecision[ii] << "% " << lPrecisionString.Data() << endl;
    }
  }
  TH1F* hCalib = new TH1F(lHistoName.Data(), "", fAnchorPointValue < 0 ? lNDesiredBoundaries - 1 : lNDesiredBoundaries, lBounds);
  hCalib->SetDirectory(0);
  hCalib->SetBinContent(0, 100.5);
  hCalib->SetBinContent(1, 100.5);
  for (Long_t ibin = fAnchorPointValue < 0 ? 1 : 2; ibin < lNDesiredBoundaries; ibin++) {
    hCalib->SetBinContent(ibin, lMiddleOfBins[fAnchorPointValue < 0 ? ibin - 1 : ibin - 2]);
    fPrecisionHistogram->SetBinContent(lNDesiredBoundaries - ibin + 1, std::hypot(lPrecision[ibin - 1], lPrecision[ibin]));
  }
  return hCalib;
}

//________________________________________________________________
void multCalibrator::ResetPrecisionHistogram()
{
  if (fPrecisionHistogram) {
    delete fPrecisionHistogram;
    fPrecisionHistogram = 0x0;
  }
  if (lNDesiredBoundaries > 0) { // only if initialized
    // invert boundaries, please
    Double_t lInverseDesiredBoundaries[1100];
    for (Int_t ii = 0; ii < lNDesiredBoundaries; ii++) {
      lInverseDesiredBoundaries[ii] = lDesiredBoundaries[lNDesiredBoundaries - (ii + 1)];
    }
    fPrecisionHistogram = new TH1D("hPrecisionHistogram", "", lNDesiredBoundaries - 1, lInverseDesiredBoundaries);
  }
}
