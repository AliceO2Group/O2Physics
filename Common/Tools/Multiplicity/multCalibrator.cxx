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
#include "TList.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TStopwatch.h"
#include "TArrayL64.h"
#include "TArrayF.h"
#include "multCalibrator.h"

multCalibrator::multCalibrator() : TNamed(),
                                   lDesiredBoundaries(0),
                                   lNDesiredBoundaries(0),
                                   fkPrecisionWarningThreshold(1.0),
                                   fInputFileName("AnalysisResults.root"),
                                   fOutputFileName("CCDB-objects.root"),
                                   fCalibHists(0x0)
{
  // Constructor
  // Make sure the TList owns its objects
  fCalibHists = new TList();
  fCalibHists->SetOwner(kTRUE);
}

multCalibrator::multCalibrator(const char* name, const char* title) : TNamed(name, title),
                                                                      lDesiredBoundaries(0),
                                                                      lNDesiredBoundaries(0),
                                                                      fkPrecisionWarningThreshold(1.0),
                                                                      fInputFileName("AnalysisResults.root"),
                                                                      fOutputFileName("CCDB-objects.root"),
                                                                      fCalibHists(0x0)
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

  //Step 1: verify if input file contains desired histograms
  TH1D* hRawV0 = (TH1D*)fileInput->Get("multiplicity-qa/multiplicityQa/hRawV0");
  TH1D* hRawT0 = (TH1D*)fileInput->Get("multiplicity-qa/multiplicityQa/hRawT0");
  TH1D* hRawFDD = (TH1D*)fileInput->Get("multiplicity-qa/multiplicityQa/hRawFDD");
  TH1D* hRawNTracks = (TH1D*)fileInput->Get("multiplicity-qa/multiplicityQa/hRawNTracks");
  TH1D* hZeqV0 = (TH1D*)fileInput->Get("multiplicity-qa/multiplicityQa/hZeqV0");
  TH1D* hZeqT0 = (TH1D*)fileInput->Get("multiplicity-qa/multiplicityQa/hZeqT0");
  TH1D* hZeqFDD = (TH1D*)fileInput->Get("multiplicity-qa/multiplicityQa/hZeqFDD");
  TH1D* hZeqNTracks = (TH1D*)fileInput->Get("multiplicity-qa/multiplicityQa/hZeqNTracks");

  if (!hRawV0 || !hRawT0 || !hRawFDD || !hRawNTracks || !hZeqV0 || !hZeqT0 || !hZeqFDD || !hZeqNTracks) {
    cout << "File does not contain histograms necessary for calibration!" << endl;
    return kFALSE;
  }
  cout << "Histograms loaded! Will now calibrate..." << endl;
  Double_t lRawV0Bounds[lNDesiredBoundaries];
  Double_t lRawT0Bounds[lNDesiredBoundaries];
  Double_t lRawFDDBounds[lNDesiredBoundaries];
  Double_t lRawNTracksBounds[lNDesiredBoundaries];
  Double_t lZeqV0Bounds[lNDesiredBoundaries];
  Double_t lZeqT0Bounds[lNDesiredBoundaries];
  Double_t lZeqFDDBounds[lNDesiredBoundaries];
  Double_t lZeqNTracksBounds[lNDesiredBoundaries];

  //Aux vars
  Double_t lMiddleOfBins[1000];
  for (Long_t lB = 1; lB < lNDesiredBoundaries; lB++) {
    //place squarely at the middle to ensure it's all fine
    lMiddleOfBins[lB - 1] = 0.5 * (lDesiredBoundaries[lB] + lDesiredBoundaries[lB - 1]);
  }

  //Create output file
  TFile* fOut = new TFile(fOutputFileName.Data(), "RECREATE");

  //_________________________________________________________________________
  cout << "Raw V0 calibration" << endl;
  for (Int_t ii = 0; ii < lNDesiredBoundaries; ii++) {
    Double_t lPrecision = 0;
    lRawV0Bounds[ii] = GetBoundaryForPercentile(hRawV0, lDesiredBoundaries[ii], lPrecision);

    TString lPrecisionString = "(Precision OK)";
    if (ii != 0 && ii != lNDesiredBoundaries - 1) {
      //check precision, please
      if (lPrecision / TMath::Abs(lDesiredBoundaries[ii + 1] - lDesiredBoundaries[ii]) > fkPrecisionWarningThreshold)
        lPrecisionString = "(WARNING: BINNING MAY LEAD TO IMPRECISION!)";
      if (lPrecision / TMath::Abs(lDesiredBoundaries[ii - 1] - lDesiredBoundaries[ii]) > fkPrecisionWarningThreshold)
        lPrecisionString = "(WARNING: BINNING MAY LEAD TO IMPRECISION!)";
    }
    cout << "Raw V0 bound percentile: " << lDesiredBoundaries[ii] << "%\t Signal value = " << lRawV0Bounds[ii] << "\tprecision = " << lPrecision << "% " << lPrecisionString.Data() << endl;
  }
  TH1F* hCalibRawV0 = new TH1F("hCalibRawV0", "", lNDesiredBoundaries - 1, lRawV0Bounds);
  hCalibRawV0->SetDirectory(0);
  hCalibRawV0->SetBinContent(0, 100.5);
  for (Long_t ibin = 1; ibin < lNDesiredBoundaries; ibin++) {
    hCalibRawV0->SetBinContent(ibin, lMiddleOfBins[ibin - 1]);
  }
  hCalibRawV0->Write();
  //_________________________________________________________________________
  cout << "Raw T0 calibration" << endl;
  for (Int_t ii = 0; ii < lNDesiredBoundaries; ii++) {
    Double_t lPrecision = 0;
    lRawT0Bounds[ii] = GetBoundaryForPercentile(hRawT0, lDesiredBoundaries[ii], lPrecision);
    TString lPrecisionString = "(Precision OK)";
    if (ii != 0 && ii != lNDesiredBoundaries - 1) {
      //check precision, please
      if (lPrecision / TMath::Abs(lDesiredBoundaries[ii + 1] - lDesiredBoundaries[ii]) > fkPrecisionWarningThreshold)
        lPrecisionString = "(WARNING: BINNING MAY LEAD TO IMPRECISION!)";
      if (lPrecision / TMath::Abs(lDesiredBoundaries[ii - 1] - lDesiredBoundaries[ii]) > fkPrecisionWarningThreshold)
        lPrecisionString = "(WARNING: BINNING MAY LEAD TO IMPRECISION!)";
    }
    cout << "Raw T0 bound percentile: " << lDesiredBoundaries[ii] << "%\t Signal value = " << lRawT0Bounds[ii] << "\tprecision = " << lPrecision << "% " << lPrecisionString.Data() << endl;
  }
  TH1F* hCalibRawT0 = new TH1F("hCalibRawT0", "", lNDesiredBoundaries - 1, lRawT0Bounds);
  hCalibRawT0->SetDirectory(0);
  hCalibRawT0->SetBinContent(0, 100.5);
  for (Long_t ibin = 1; ibin < lNDesiredBoundaries; ibin++) {
    hCalibRawT0->SetBinContent(ibin, lMiddleOfBins[ibin - 1]);
  }
  hCalibRawT0->Write();
  //_________________________________________________________________________
  cout << "Raw FDD calibration" << endl;
  for (Int_t ii = 0; ii < lNDesiredBoundaries; ii++) {
    Double_t lPrecision = 0;
    lRawFDDBounds[ii] = GetBoundaryForPercentile(hRawFDD, lDesiredBoundaries[ii], lPrecision);
    TString lPrecisionString = "(Precision OK)";
    if (ii != 0 && ii != lNDesiredBoundaries - 1) {
      //check precision, please
      if (lPrecision / TMath::Abs(lDesiredBoundaries[ii + 1] - lDesiredBoundaries[ii]) > fkPrecisionWarningThreshold)
        lPrecisionString = "(WARNING: BINNING MAY LEAD TO IMPRECISION!)";
      if (lPrecision / TMath::Abs(lDesiredBoundaries[ii - 1] - lDesiredBoundaries[ii]) > fkPrecisionWarningThreshold)
        lPrecisionString = "(WARNING: BINNING MAY LEAD TO IMPRECISION!)";
    }
    cout << "Raw FDD bound percentile: " << lDesiredBoundaries[ii] << "%\t Signal value = " << lRawFDDBounds[ii] << "\tprecision = " << lPrecision << "% " << lPrecisionString.Data() << endl;
  }
  TH1F* hCalibRawFDD = new TH1F("hCalibRawFDD", "", lNDesiredBoundaries - 1, lRawFDDBounds);
  hCalibRawFDD->SetDirectory(0);
  hCalibRawFDD->SetBinContent(0, 100.5);
  for (Long_t ibin = 1; ibin < lNDesiredBoundaries; ibin++) {
    hCalibRawFDD->SetBinContent(ibin, lMiddleOfBins[ibin - 1]);
  }
  hCalibRawFDD->Write();
  //_________________________________________________________________________
  cout << "Raw NTracks calibration" << endl;
  for (Int_t ii = 0; ii < lNDesiredBoundaries; ii++) {
    Double_t lPrecision = 0;
    lRawNTracksBounds[ii] = GetBoundaryForPercentile(hRawNTracks, lDesiredBoundaries[ii], lPrecision);
    TString lPrecisionString = "(Precision OK)";
    if (ii != 0 && ii != lNDesiredBoundaries - 1) {
      //check precision, please
      if (lPrecision / TMath::Abs(lDesiredBoundaries[ii + 1] - lDesiredBoundaries[ii]) > fkPrecisionWarningThreshold)
        lPrecisionString = "(WARNING: BINNING MAY LEAD TO IMPRECISION!)";
      if (lPrecision / TMath::Abs(lDesiredBoundaries[ii - 1] - lDesiredBoundaries[ii]) > fkPrecisionWarningThreshold)
        lPrecisionString = "(WARNING: BINNING MAY LEAD TO IMPRECISION!)";
    }
    cout << "Raw NTracks bound percentile: " << lDesiredBoundaries[ii] << "%\t Signal value = " << lRawNTracksBounds[ii] << "\tprecision = " << lPrecision << "% " << lPrecisionString.Data() << endl;
  }
  TH1F* hCalibRawNTracks = new TH1F("hCalibRawNTracks", "", lNDesiredBoundaries - 1, lRawNTracksBounds);
  hCalibRawNTracks->SetDirectory(0);
  hCalibRawNTracks->SetBinContent(0, 100.5); //Just in case correction functions screw up the values
  for (Long_t ibin = 1; ibin < lNDesiredBoundaries; ibin++) {
    hCalibRawNTracks->SetBinContent(ibin, lMiddleOfBins[ibin - 1]);
  }
  hCalibRawNTracks->Write();

  //_________________________________________________________________________
  cout << "Vertex-Z equalized V0 calibration" << endl;
  for (Int_t ii = 0; ii < lNDesiredBoundaries; ii++) {
    Double_t lPrecision = 0;
    lZeqV0Bounds[ii] = GetBoundaryForPercentile(hZeqV0, lDesiredBoundaries[ii], lPrecision);

    TString lPrecisionString = "(Precision OK)";
    if (ii != 0 && ii != lNDesiredBoundaries - 1) {
      //check precision, please
      if (lPrecision / TMath::Abs(lDesiredBoundaries[ii + 1] - lDesiredBoundaries[ii]) > fkPrecisionWarningThreshold)
        lPrecisionString = "(WARNING: BINNING MAY LEAD TO IMPRECISION!)";
      if (lPrecision / TMath::Abs(lDesiredBoundaries[ii - 1] - lDesiredBoundaries[ii]) > fkPrecisionWarningThreshold)
        lPrecisionString = "(WARNING: BINNING MAY LEAD TO IMPRECISION!)";
    }
    cout << "Vertex-Z equalized V0 bound percentile: " << lDesiredBoundaries[ii] << "%\t Signal value = " << lZeqV0Bounds[ii] << "\tprecision = " << lPrecision << "% " << lPrecisionString.Data() << endl;
  }
  TH1F* hCalibZeqV0 = new TH1F("hCalibZeqV0", "", lNDesiredBoundaries - 1, lZeqV0Bounds);
  hCalibZeqV0->SetDirectory(0);
  hCalibZeqV0->SetBinContent(0, 100.5);
  for (Long_t ibin = 1; ibin < lNDesiredBoundaries; ibin++) {
    hCalibZeqV0->SetBinContent(ibin, lMiddleOfBins[ibin - 1]);
  }
  hCalibZeqV0->Write();
  //_________________________________________________________________________
  cout << "Vertex-Z equalized T0 calibration" << endl;
  for (Int_t ii = 0; ii < lNDesiredBoundaries; ii++) {
    Double_t lPrecision = 0;
    lZeqT0Bounds[ii] = GetBoundaryForPercentile(hZeqT0, lDesiredBoundaries[ii], lPrecision);
    TString lPrecisionString = "(Precision OK)";
    if (ii != 0 && ii != lNDesiredBoundaries - 1) {
      //check precision, please
      if (lPrecision / TMath::Abs(lDesiredBoundaries[ii + 1] - lDesiredBoundaries[ii]) > fkPrecisionWarningThreshold)
        lPrecisionString = "(WARNING: BINNING MAY LEAD TO IMPRECISION!)";
      if (lPrecision / TMath::Abs(lDesiredBoundaries[ii - 1] - lDesiredBoundaries[ii]) > fkPrecisionWarningThreshold)
        lPrecisionString = "(WARNING: BINNING MAY LEAD TO IMPRECISION!)";
    }
    cout << "Vertex-Z equalized T0 bound percentile: " << lDesiredBoundaries[ii] << "%\t Signal value = " << lZeqT0Bounds[ii] << "\tprecision = " << lPrecision << "% " << lPrecisionString.Data() << endl;
  }
  TH1F* hCalibZeqT0 = new TH1F("hCalibZeqT0", "", lNDesiredBoundaries - 1, lZeqT0Bounds);
  hCalibZeqT0->SetDirectory(0);
  hCalibZeqT0->SetBinContent(0, 100.5);
  for (Long_t ibin = 1; ibin < lNDesiredBoundaries; ibin++) {
    hCalibZeqT0->SetBinContent(ibin, lMiddleOfBins[ibin - 1]);
  }
  hCalibZeqT0->Write();
  //_________________________________________________________________________
  cout << "Vertex-Z equalized FDD calibration" << endl;
  for (Int_t ii = 0; ii < lNDesiredBoundaries; ii++) {
    Double_t lPrecision = 0;
    lZeqFDDBounds[ii] = GetBoundaryForPercentile(hZeqFDD, lDesiredBoundaries[ii], lPrecision);
    TString lPrecisionString = "(Precision OK)";
    if (ii != 0 && ii != lNDesiredBoundaries - 1) {
      //check precision, please
      if (lPrecision / TMath::Abs(lDesiredBoundaries[ii + 1] - lDesiredBoundaries[ii]) > fkPrecisionWarningThreshold)
        lPrecisionString = "(WARNING: BINNING MAY LEAD TO IMPRECISION!)";
      if (lPrecision / TMath::Abs(lDesiredBoundaries[ii - 1] - lDesiredBoundaries[ii]) > fkPrecisionWarningThreshold)
        lPrecisionString = "(WARNING: BINNING MAY LEAD TO IMPRECISION!)";
    }
    cout << "Vertex-Z equalized FDD bound percentile: " << lDesiredBoundaries[ii] << "%\t Signal value = " << lZeqFDDBounds[ii] << "\tprecision = " << lPrecision << "% " << lPrecisionString.Data() << endl;
  }
  TH1F* hCalibZeqFDD = new TH1F("hCalibZeqFDD", "", lNDesiredBoundaries - 1, lZeqFDDBounds);
  hCalibZeqFDD->SetDirectory(0);
  hCalibZeqFDD->SetBinContent(0, 100.5);
  for (Long_t ibin = 1; ibin < lNDesiredBoundaries; ibin++) {
    hCalibZeqFDD->SetBinContent(ibin, lMiddleOfBins[ibin - 1]);
  }
  hCalibZeqFDD->Write();
  //_________________________________________________________________________
  cout << "Vertex-Z equalized NTracks calibration" << endl;
  for (Int_t ii = 0; ii < lNDesiredBoundaries; ii++) {
    Double_t lPrecision = 0;
    lZeqNTracksBounds[ii] = GetBoundaryForPercentile(hZeqNTracks, lDesiredBoundaries[ii], lPrecision);
    TString lPrecisionString = "(Precision OK)";
    if (ii != 0 && ii != lNDesiredBoundaries - 1) {
      //check precision, please
      if (lPrecision / TMath::Abs(lDesiredBoundaries[ii + 1] - lDesiredBoundaries[ii]) > fkPrecisionWarningThreshold)
        lPrecisionString = "(WARNING: BINNING MAY LEAD TO IMPRECISION!)";
      if (lPrecision / TMath::Abs(lDesiredBoundaries[ii - 1] - lDesiredBoundaries[ii]) > fkPrecisionWarningThreshold)
        lPrecisionString = "(WARNING: BINNING MAY LEAD TO IMPRECISION!)";
    }
    cout << "Vertex-Z equalized NTracks bound percentile: " << lDesiredBoundaries[ii] << "%\t Signal value = " << lZeqNTracksBounds[ii] << "\tprecision = " << lPrecision << "% " << lPrecisionString.Data() << endl;
  }
  TH1F* hCalibZeqNTracks = new TH1F("hCalibZeqNTracks", "", lNDesiredBoundaries - 1, lZeqNTracksBounds);
  hCalibZeqNTracks->SetDirectory(0);
  hCalibZeqNTracks->SetBinContent(0, 100.5); //Just in case correction functions screw up the values
  for (Long_t ibin = 1; ibin < lNDesiredBoundaries; ibin++) {
    hCalibZeqNTracks->SetBinContent(ibin, lMiddleOfBins[ibin - 1]);
  }
  hCalibZeqNTracks->Write();

  cout << "Saving calibration file..." << endl;
  fOut->Write();
  cout << "Done! Enjoy!" << endl;
  return kTRUE;
}

Double_t multCalibrator::GetBoundaryForPercentile(TH1D* histo, Double_t lPercentileRequested, Double_t& lPrecisionEstimate)
{
  //This function returns the boundary for a specific percentile.
  //It uses a linear interpolation in an attempt to get more precision
  //than the binning of the histogram used for quantiling.
  //
  //It also estimates a certain level of precision of the procedure
  //by explicitly comparing the bin content of the bins around the boundary
  //with the entire cross section, effectively reporting back a percentage
  //that corresponds to those bins. If this percentage is O(percentile bin
  //width requested), then the user should worry and we print out a warning.

  //if( lPercentileRequested < 1e-7 ) return 1e+6; //safeguard
  if (lPercentileRequested > 100 - 1e-7)
    return 0.0; //safeguard

  Double_t lReturnValue = 0.0;
  Double_t lPercentile = 100.0 - lPercentileRequested;
  lPrecisionEstimate = -1;

  const Long_t lNBins = histo->GetNbinsX();
  Double_t lCountDesired = lPercentile * histo->GetEntries() / 100;
  Long_t lCount = 0;
  for (Long_t ibin = 1; ibin < lNBins; ibin++) {
    lCount += histo->GetBinContent(ibin);
    if (lCount >= lCountDesired) {
      //Found bin I am looking for!
      Double_t lWidth = histo->GetBinWidth(ibin);
      Double_t lLeftPercentile = 100. * (lCount - histo->GetBinContent(ibin)) / histo->GetEntries();
      Double_t lRightPercentile = 100. * lCount / histo->GetEntries();
      lPrecisionEstimate = (lRightPercentile - lLeftPercentile) / 2;

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
  //Function to set standard adaptive boundaries
  //Typically used in pp, goes to 0.001% binning for highest multiplicity
  lNDesiredBoundaries = 0;
  lDesiredBoundaries = new Double_t[1100];
  lDesiredBoundaries[0] = 100;
  //From Low To High Multiplicity
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
  lDesiredBoundaries[lNDesiredBoundaries - 1] = 0;
}
