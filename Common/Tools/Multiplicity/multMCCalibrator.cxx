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
#include "multMCCalibrator.h"

#include "multCalibrator.h"

#include <TDirectory.h>
#include <TF1.h>
#include <TFile.h>
#include <TList.h>
#include <TNamed.h>
#include <TProfile.h>
#include <TString.h>

#include <RtypesCore.h>

#include <iostream> // FIXME

using namespace std;

multMCCalibrator::multMCCalibrator() : TNamed(),
                                       fDataInputFileName("AnalysisResults.root"),
                                       fSimInputFileName("AnalysisResultsMC.root"),
                                       fOutputFileName("CCDB-objects.root"),
                                       fCalibHists(0x0)
{
  // Constructor
  // Make sure the TList owns its objects
  fCalibHists = new TList();
  fCalibHists->SetOwner(kTRUE);
}

multMCCalibrator::multMCCalibrator(const char* name, const char* title) : TNamed(name, title),
                                                                          fDataInputFileName("AnalysisResults.root"),
                                                                          fSimInputFileName("AnalysisResultsMC.root"),
                                                                          fOutputFileName("CCDB-objects.root"),
                                                                          fCalibHists(0x0)
{
  // Named Constructor
  // Make sure the TList owns its objects
  fCalibHists = new TList();
  fCalibHists->SetOwner(kTRUE);
}
//________________________________________________________________
multMCCalibrator::~multMCCalibrator()
{
  // Destructor
  if (fCalibHists) {
    delete fCalibHists;
    fCalibHists = 0x0;
  }
}

//________________________________________________________________
Bool_t multMCCalibrator::Calibrate()
{
  // Function meant to generate calibration OADB
  //
  // --- input : fInputFileName, containing a TTree object
  // --- output: fOutputFileName, containing OABD object
  //

  cout << "=== STARTING CALIBRATION PROCEDURE ===" << endl;
  cout << " * Data Input File........: " << fDataInputFileName.Data() << endl;
  cout << " * Simulation Input File..: " << fSimInputFileName.Data() << endl;
  cout << " * Output File............: " << fOutputFileName.Data() << endl;
  cout << endl;

  // Opening data and simulation file...
  TFile* fileData = new TFile(fDataInputFileName.Data(), "READ");
  TFile* fileSim = new TFile(fSimInputFileName.Data(), "READ");

  // Step 1: verify if input file contains desired histograms
  TProfile* hProfData[multCalibrator::kNCentEstim];
  TProfile* hProfSim[multCalibrator::kNCentEstim];
  cout << " * acquiring input profiles..." << endl;
  for (Int_t iv = 0; iv < multCalibrator::kNCentEstim; iv++) {
    hProfData[iv] = reinterpret_cast<TProfile*>(fileData->Get(Form("multiplicity-qa/multiplicityQa/hProf%s", multCalibrator::fCentEstimName[iv].Data())));
    if (!hProfData[iv]) {
      cout << Form("Data file does not contain histogram h%s, which is necessary for calibration!", multCalibrator::fCentEstimName[iv].Data()) << endl;
      return kFALSE;
    }
    hProfData[iv]->SetName(Form("hProfData_%s", multCalibrator::fCentEstimName[iv].Data()));
    hProfSim[iv] = reinterpret_cast<TProfile*>(fileSim->Get(Form("multiplicity-qa/multiplicityQa/hProf%s", multCalibrator::fCentEstimName[iv].Data())));
    if (!hProfSim[iv]) {
      cout << Form("Sim file does not contain histogram h%s, which is necessary for calibration!", multCalibrator::fCentEstimName[iv].Data()) << endl;
      return kFALSE;
    }
    hProfSim[iv]->SetName(Form("hProfSim_%s", multCalibrator::fCentEstimName[iv].Data()));
  }
  cout << " * fitting profiles..." << endl;

  TF1* hFitData[multCalibrator::kNCentEstim];
  TF1* hFitSim[multCalibrator::kNCentEstim];
  TF1* hMapping[multCalibrator::kNCentEstim];
  for (Int_t iv = 0; iv < multCalibrator::kNCentEstim; iv++) {
    hFitData[iv] = GetFit(hProfData[iv]);
    fCalibHists->Add(hFitData[iv]);
  }
  for (Int_t iv = 0; iv < multCalibrator::kNCentEstim; iv++) {
    hFitSim[iv] = GetFit(hProfSim[iv]);
    fCalibHists->Add(hFitSim[iv]);
  }

  cout << " * creating maps..." << endl;
  for (Int_t iv = 0; iv < multCalibrator::kNCentEstim; iv++) {
    TString lTempDef = Form("TMath::Power(( (%.10f + %.10f * TMath::Power(%s,%.10f)) - %.10f )/%.10f, 1./%.10f)",
                            hFitSim[iv]->GetParameter(0),
                            hFitSim[iv]->GetParameter(1),
                            "x",
                            hFitSim[iv]->GetParameter(2),
                            hFitData[iv]->GetParameter(0),
                            hFitData[iv]->GetParameter(1),
                            hFitData[iv]->GetParameter(2));
    hMapping[iv] = new TF1(Form("hMapping_%s", multCalibrator::fCentEstimName[iv].Data()), lTempDef.Data(), hProfData[iv]->GetBinLowEdge(1), hProfData[iv]->GetBinLowEdge(hProfData[iv]->GetNbinsX()));
    fCalibHists->Add(hMapping[iv]);
  }
  cout << " * done! Number of objects in calibration list: " << fCalibHists->GetEntries() << endl;
  return kTRUE;
}

//________________________________________________________________
TF1* multMCCalibrator::GetFit(TProfile* fProf, Bool_t lQuadratic)
{
  TString fFormula = "[0]*x"; // old/deprecated (avoid if possible, please)
  if (lQuadratic)
    fFormula = "[0]+[1]*TMath::Power(x,[2])";

  // Function to return fit function to profile for posterior inversion
  TF1* fit = new TF1(Form("%s_fit", fProf->GetName()), fFormula.Data(), fProf->GetBinLowEdge(1), fProf->GetBinLowEdge(fProf->GetNbinsX()));

  // Guesstimate inclination from data points in profile
  Double_t lMeanInclination = 0;
  Long_t lInclinationCount = 0;
  for (Int_t ii = 2; ii < fProf->GetNbinsX(); ii++) {
    if (fProf->GetBinContent(ii) < 1e-10)
      continue;
    if (fProf->GetBinError(ii) / fProf->GetBinContent(ii) > 0.1)
      continue;
    lMeanInclination = fProf->GetBinContent(ii) / fProf->GetBinCenter(ii);
    lInclinationCount++;
  }
  if (lInclinationCount < 5)
    lMeanInclination = 1;
  if (lInclinationCount >= 5)
    lMeanInclination /= lInclinationCount;

  // Give it a little nudge, cause life's hard
  fit->SetParameter(0, 0.0);
  fit->SetParameter(1, lMeanInclination);
  fit->SetParameter(2, 1.0);

  fProf->Fit(fit->GetName(), "IREM0");
  return fit;
}
