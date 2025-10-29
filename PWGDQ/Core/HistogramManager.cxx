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

#include "PWGDQ/Core/HistogramManager.h"

#include <iostream>
#include <memory>
#include <fstream>
#include <list>
#include <vector>
#include <algorithm>
#include "Framework/Logger.h"
using namespace std;

#include <TObject.h>
#include <TObjArray.h>
#include <THashList.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TProfile3D.h>
#include <THn.h>
#include <THnSparse.h>
#include <TIterator.h>
#include <TClass.h>

ClassImp(HistogramManager);

//_______________________________________________________________________________
HistogramManager::HistogramManager() : TNamed("", ""),
                                       fMainList(nullptr),
                                       fNVars(0),
                                       fUsedVars(nullptr),
                                       fVariablesMap(),
                                       fUseDefaultVariableNames(false),
                                       fBinsAllocated(0),
                                       fVariableNames(nullptr),
                                       fVariableUnits(nullptr)
{
  //
  // Constructor
  //
}

//_______________________________________________________________________________
HistogramManager::HistogramManager(const char* name, const char* title, const int maxNVars) : TNamed(name, title),
                                                                                              fMainList(),
                                                                                              fNVars(maxNVars),
                                                                                              fUsedVars(),
                                                                                              fVariablesMap(),
                                                                                              fUseDefaultVariableNames(kFALSE),
                                                                                              fBinsAllocated(0),
                                                                                              fVariableNames(),
                                                                                              fVariableUnits()
{
  //
  // Constructor
  //
  fMainList = new THashList;
  fMainList->SetOwner(kTRUE);
  fMainList->SetName(name);
  fUsedVars = new bool[maxNVars];
  for (int i = 0; i < maxNVars; ++i) {
    fUsedVars[i] = false;
  }
  fVariableNames = new TString[maxNVars];
  fVariableUnits = new TString[maxNVars];
}

//_______________________________________________________________________________
HistogramManager::~HistogramManager()
{
  //
  // De-constructor
  //
  delete fMainList;
  delete[] fUsedVars;
}

//_______________________________________________________________________________
void HistogramManager::SetDefaultVarNames(TString* vars, TString* units)
{
  //
  // Set default variable names
  //
  for (int i = 0; i < fNVars; ++i) {
    fVariableNames[i] = vars[i];
    fVariableUnits[i] = units[i];
  }
};

//__________________________________________________________________
void HistogramManager::AddHistClass(const char* histClass)
{
  //
  // Add a new histogram list
  //
  if (fMainList->FindObject(histClass)) {
    LOG(warn) << "HistogramManager::AddHistClass(): Cannot add histogram class " << histClass
              << " because it already exists.";
    return;
  }
  auto* hList = new TList;
  hList->SetOwner(kTRUE);
  hList->SetName(histClass);
  fMainList->Add(hList);
  std::list<std::vector<int>> varList;
  fVariablesMap[histClass] = varList;
}

//_________________________________________________________________
void HistogramManager::AddHistogram(const char* histClass, const char* hname, const char* title, bool isProfile,
                                    int nXbins, double xmin, double xmax, int varX,
                                    int nYbins, double ymin, double ymax, int varY,
                                    int nZbins, double zmin, double zmax, int varZ,
                                    const char* xLabels, const char* yLabels, const char* zLabels,
                                    int varT, int varW, bool isdouble, bool isFillLabelx)
{
  //
  // add a histogram  (this function can define TH1F,TH2F,TH3F,TProfile,TProfile2D, and TProfile3D)
  //
  // TODO: replace the cout warning messages with LOG (same for all the other functions)

  // get the list to which the histogram should be added
  auto* hList = reinterpret_cast<TList*>(fMainList->FindObject(histClass));
  if (!hList) {
    LOG(warn) << "HistogramManager::AddHistogram(): Histogram list " << histClass << " not found!";
    LOG(warn) << "         Histogram not created";
    return;
  }
  // check whether this histogram name was used before
  if (hList->FindObject(hname)) {
    LOG(warn) << "HistogramManager::AddHistogram(): Histogram " << hname << " already exists in class " << histClass;
    return;
  }

  // deduce the dimension of the histogram from parameters
  // NOTE: in case of profile histograms, one extra variable is needed
  int dimension = 1;
  if (varY > kNothing) {
    dimension = 2;
  }
  if (varZ > kNothing) {
    dimension = 3;
  }

  // tokenize the title string; the user may include in it axis titles which will overwrite the defaults
  TString titleStr(title);
  std::unique_ptr<TObjArray> arr(titleStr.Tokenize(";"));
  // mark required variables as being used
  if (varX > kNothing) {
    fUsedVars[varX] = kTRUE;
  }
  if (varY > kNothing) {
    fUsedVars[varY] = kTRUE;
  }
  if (varZ > kNothing) {
    fUsedVars[varZ] = kTRUE;
  }
  if (varT > kNothing) {
    fUsedVars[varT] = kTRUE;
  }
  if (varW > kNothing) {
    fUsedVars[varW] = kTRUE;
  }

  // encode needed variable identifiers in a vector and push it to the std::list corresponding to the current histogram list
  std::vector<int> varVector;
  varVector.push_back(isProfile ? 1 : 0); // whether the histogram is a profile
  varVector.push_back(0);                 // whether it is a THn
  varVector.push_back(varW);              // variable used for weighting
  varVector.push_back(varX);              // variables on each axis
  varVector.push_back(varY);
  varVector.push_back(varZ);
  varVector.push_back(varT);                 // variable used for profiling in case of TProfile3D
  varVector.push_back(isFillLabelx ? 1 : 0); // whether to fill with the x-axis labels
  std::list varList = fVariablesMap[histClass];
  varList.push_back(varVector);
  fVariablesMap[histClass] = varList;

  // create and configure histograms according to required options
  TH1* h = nullptr;
  switch (dimension) {
    case 1: // TH1F
      if (!isdouble) {
        h = new TH1F(hname, (arr->At(0) ? arr->At(0)->GetName() : ""), nXbins, xmin, xmax);
      } else {
        h = new TH1D(hname, (arr->At(0) ? arr->At(0)->GetName() : ""), nXbins, xmin, xmax);
      }
      fBinsAllocated += nXbins + 2;
      // TODO: possibly make the call of Sumw2() optional for all histograms
      h->Sumw2();
      if (fVariableNames[varX][0]) {
        h->GetXaxis()->SetTitle(Form("%s %s", fVariableNames[varX].Data(),
                                     (fVariableUnits[varX][0] ? Form("(%s)", fVariableUnits[varX].Data()) : "")));
      }
      if (arr->At(1)) {
        h->GetXaxis()->SetTitle(arr->At(1)->GetName());
      }
      if (xLabels[0] != '\0') {
        MakeAxisLabels(h->GetXaxis(), xLabels);
      }
      hList->Add(h);
      h->SetDirectory(nullptr);
      break;

    case 2: // either TH2F or TProfile
      if (isProfile) {
        h = new TProfile(hname, (arr->At(0) ? arr->At(0)->GetName() : ""), nXbins, xmin, xmax);
        fBinsAllocated += nXbins + 2;
        h->Sumw2();
        // if requested, build the profile using the profile widths instead of stat errors
        // TODO: make this option more transparent to the user ?
        if (titleStr.Contains("--s--")) {
          (reinterpret_cast<TProfile*>(h))->BuildOptions(0., 0., "s");
        }
      } else {
        if (!isdouble) {
          h = new TH2F(hname, (arr->At(0) ? arr->At(0)->GetName() : ""), nXbins, xmin, xmax, nYbins, ymin, ymax);
        } else {
          h = new TH2D(hname, (arr->At(0) ? arr->At(0)->GetName() : ""), nXbins, xmin, xmax, nYbins, ymin, ymax);
        }
        fBinsAllocated += (nXbins + 2) * (nYbins + 2);
        h->Sumw2();
      }
      if (fVariableNames[varX][0]) {
        h->GetXaxis()->SetTitle(Form("%s %s", fVariableNames[varX].Data(),
                                     (fVariableUnits[varX][0] ? Form("(%s)", fVariableUnits[varX].Data()) : "")));
      }
      if (arr->At(1)) {
        h->GetXaxis()->SetTitle(arr->At(1)->GetName());
      }
      if (xLabels[0] != '\0') {
        MakeAxisLabels(h->GetXaxis(), xLabels);
      }

      if (fVariableNames[varY][0]) {
        h->GetYaxis()->SetTitle(Form("%s %s", fVariableNames[varY].Data(),
                                     (fVariableUnits[varY][0] ? Form("(%s)", fVariableUnits[varY].Data()) : "")));
      }
      if (fVariableNames[varY][0] && isProfile) {
        h->GetYaxis()->SetTitle(Form("<%s> %s", fVariableNames[varY].Data(),
                                     (fVariableUnits[varY][0] ? Form("(%s)", fVariableUnits[varY].Data()) : "")));
      }
      if (arr->At(2)) {
        h->GetYaxis()->SetTitle(arr->At(2)->GetName());
      }
      if (yLabels[0] != '\0') {
        MakeAxisLabels(h->GetYaxis(), yLabels);
      }
      hList->Add(h);
      h->SetDirectory(nullptr);
      break;

    case 3: // TH3F, TProfile2D or TProfile3D
      if (isProfile) {
        if (varT > kNothing) { // TProfile3D
          h = new TProfile3D(hname, (arr->At(0) ? arr->At(0)->GetName() : ""), nXbins, xmin, xmax, nYbins, ymin, ymax, nZbins, zmin, zmax);
          fBinsAllocated += (nXbins + 2) * (nYbins + 2) * (nZbins + 2);
          h->Sumw2();
          if (titleStr.Contains("--s--")) {
            (reinterpret_cast<TProfile3D*>(h))->BuildOptions(0., 0., "s");
          }
        } else { // TProfile2D
          h = new TProfile2D(hname, (arr->At(0) ? arr->At(0)->GetName() : ""), nXbins, xmin, xmax, nYbins, ymin, ymax);
          fBinsAllocated += (nXbins + 2) * (nYbins + 2);
          h->Sumw2();
          if (titleStr.Contains("--s--")) {
            (reinterpret_cast<TProfile2D*>(h))->BuildOptions(0., 0., "s");
          }
        }
      } else { // TH3F
        if (!isdouble) {
          h = new TH3F(hname, (arr->At(0) ? arr->At(0)->GetName() : ""), nXbins, xmin, xmax, nYbins, ymin, ymax, nZbins, zmin, zmax);
        } else {
          h = new TH3D(hname, (arr->At(0) ? arr->At(0)->GetName() : ""), nXbins, xmin, xmax, nYbins, ymin, ymax, nZbins, zmin, zmax);
        }
        fBinsAllocated += (nXbins + 2) * (nYbins + 2) * (nZbins + 2);
        h->Sumw2();
      }
      if (fVariableNames[varX][0]) {
        h->GetXaxis()->SetTitle(Form("%s %s", fVariableNames[varX].Data(),
                                     (fVariableUnits[varX][0] ? Form("(%s)", fVariableUnits[varX].Data()) : "")));
      }
      if (arr->At(1)) {
        h->GetXaxis()->SetTitle(arr->At(1)->GetName());
      }
      if (xLabels[0] != '\0') {
        MakeAxisLabels(h->GetXaxis(), xLabels);
      }
      if (fVariableNames[varY][0]) {
        h->GetYaxis()->SetTitle(Form("%s %s", fVariableNames[varY].Data(),
                                     (fVariableUnits[varY][0] ? Form("(%s)", fVariableUnits[varY].Data()) : "")));
      }
      if (arr->At(2)) {
        h->GetYaxis()->SetTitle(arr->At(2)->GetName());
      }
      if (yLabels[0] != '\0') {
        MakeAxisLabels(h->GetYaxis(), yLabels);
      }
      if (fVariableNames[varZ][0]) {
        h->GetZaxis()->SetTitle(Form("%s %s", fVariableNames[varZ].Data(),
                                     (fVariableUnits[varZ][0] ? Form("(%s)", fVariableUnits[varZ].Data()) : "")));
      }
      if (fVariableNames[varZ][0] && isProfile && varT < 0) { // for TProfile2D
        h->GetZaxis()->SetTitle(Form("<%s> %s", fVariableNames[varZ].Data(),
                                     (fVariableUnits[varZ][0] ? Form("(%s)", fVariableUnits[varZ].Data()) : "")));
      }
      if (arr->At(3)) {
        h->GetZaxis()->SetTitle(arr->At(3)->GetName());
      }
      if (zLabels[0] != '\0') {
        MakeAxisLabels(h->GetZaxis(), zLabels);
      }
      h->SetDirectory(nullptr);
      hList->Add(h);
      break;
  } // end switch
}

//_________________________________________________________________
void HistogramManager::AddHistogram(const char* histClass, const char* hname, const char* title, bool isProfile,
                                    int nXbins, double* xbins, int varX,
                                    int nYbins, double* ybins, int varY,
                                    int nZbins, double* zbins, int varZ,
                                    const char* xLabels, const char* yLabels, const char* zLabels,
                                    int varT, int varW, bool isdouble, bool isFillLabelx)
{
  //
  // add a histogram
  //

  // get the list to which the histogram should be added
  auto* hList = reinterpret_cast<TList*>(fMainList->FindObject(histClass));
  if (!hList) {
    LOG(warn) << "HistogramManager::AddHistogram(): Histogram list " << histClass << " not found!";
    LOG(warn) << "         Histogram not created";
    return;
  }
  // check whether this histogram name was used before
  if (hList->FindObject(hname)) {
    LOG(warn) << "HistogramManager::AddHistogram(): Histogram " << hname << " already exists";
    return;
  }

  // deduce the dimension of the histogram from parameters
  // NOTE: in case of profile histograms, one extra variable is needed
  int dimension = 1;
  if (varY > kNothing) {
    dimension = 2;
  }
  if (varZ > kNothing) {
    dimension = 3;
  }

  // mark required variables as being used
  if (varX > kNothing) {
    fUsedVars[varX] = kTRUE;
  }
  if (varY > kNothing) {
    fUsedVars[varY] = kTRUE;
  }
  if (varZ > kNothing) {
    fUsedVars[varZ] = kTRUE;
  }
  if (varT > kNothing) {
    fUsedVars[varT] = kTRUE;
  }
  if (varW > kNothing) {
    fUsedVars[varW] = kTRUE;
  }

  // tokenize the title string; the user may include in it axis titles which will overwrite the defaults
  TString titleStr(title);
  std::unique_ptr<TObjArray> arr(titleStr.Tokenize(";"));

  // encode needed variable identifiers in a vector and push it to the std::list corresponding to the current histogram list
  std::vector<int> varVector;
  varVector.push_back(isProfile ? 1 : 0); // whether the histogram is a profile
  varVector.push_back(0);                 // whether it is a THn
  varVector.push_back(varW);              // variable used for weighting
  varVector.push_back(varX);              // variables on each axis
  varVector.push_back(varY);
  varVector.push_back(varZ);
  varVector.push_back(varT);                 // variable used for profiling in case of TProfile3D
  varVector.push_back(isFillLabelx ? 1 : 0); // whether to fill with the x-axis labels
  std::list varList = fVariablesMap[histClass];
  varList.push_back(varVector);
  fVariablesMap[histClass] = varList;

  TH1* h = nullptr;
  switch (dimension) {
    case 1:
      if (!isdouble) {
        h = new TH1F(hname, (arr->At(0) ? arr->At(0)->GetName() : ""), nXbins, xbins);
      } else {
        h = new TH1D(hname, (arr->At(0) ? arr->At(0)->GetName() : ""), nXbins, xbins);
      }
      fBinsAllocated += nXbins + 2;
      h->Sumw2();
      if (fVariableNames[varX][0]) {
        h->GetXaxis()->SetTitle(Form("%s %s", fVariableNames[varX].Data(),
                                     (fVariableUnits[varX][0] ? Form("(%s)", fVariableUnits[varX].Data()) : "")));
      }
      if (arr->At(1)) {
        h->GetXaxis()->SetTitle(arr->At(1)->GetName());
      }
      if (xLabels[0] != '\0') {
        MakeAxisLabels(h->GetXaxis(), xLabels);
      }
      h->SetDirectory(nullptr);
      hList->Add(h);
      break;

    case 2:
      if (isProfile) {
        h = new TProfile(hname, (arr->At(0) ? arr->At(0)->GetName() : ""), nXbins, xbins);
        fBinsAllocated += nXbins + 2;
        h->Sumw2();
        if (titleStr.Contains("--s--")) {
          (reinterpret_cast<TProfile*>(h))->BuildOptions(0., 0., "s");
        }
      } else {
        if (!isdouble) {
          h = new TH2F(hname, (arr->At(0) ? arr->At(0)->GetName() : ""), nXbins, xbins, nYbins, ybins);
        } else {
          h = new TH2D(hname, (arr->At(0) ? arr->At(0)->GetName() : ""), nXbins, xbins, nYbins, ybins);
        }
        fBinsAllocated += (nXbins + 2) * (nYbins + 2);
        h->Sumw2();
      }
      if (fVariableNames[varX][0]) {
        h->GetXaxis()->SetTitle(Form("%s %s", fVariableNames[varX].Data(),
                                     (fVariableUnits[varX][0] ? Form("(%s)", fVariableUnits[varX].Data()) : "")));
      }
      if (arr->At(1)) {
        h->GetXaxis()->SetTitle(arr->At(1)->GetName());
      }
      if (xLabels[0] != '\0') {
        MakeAxisLabels(h->GetXaxis(), xLabels);
      }
      if (fVariableNames[varY][0]) {
        h->GetYaxis()->SetTitle(Form("%s %s", fVariableNames[varY].Data(),
                                     (fVariableUnits[varY][0] ? Form("(%s)", fVariableUnits[varY].Data()) : "")));
      }
      if (fVariableNames[varY][0] && isProfile) {
        h->GetYaxis()->SetTitle(Form("<%s> %s", fVariableNames[varY].Data(),
                                     (fVariableUnits[varY][0] ? Form("(%s)", fVariableUnits[varY].Data()) : "")));
      }

      if (arr->At(2)) {
        h->GetYaxis()->SetTitle(arr->At(2)->GetName());
      }
      if (yLabels[0] != '\0') {
        MakeAxisLabels(h->GetYaxis(), yLabels);
      }
      h->SetDirectory(nullptr);
      hList->Add(h);
      break;

    case 3:
      if (isProfile) {
        if (varT > kNothing) {
          h = new TProfile3D(hname, (arr->At(0) ? arr->At(0)->GetName() : ""), nXbins, xbins, nYbins, ybins, nZbins, zbins);
          fBinsAllocated += (nXbins + 2) * (nYbins + 2) * (nZbins + 2);
          h->Sumw2();
          if (titleStr.Contains("--s--")) {
            (reinterpret_cast<TProfile3D*>(h))->BuildOptions(0., 0., "s");
          }
        } else {
          h = new TProfile2D(hname, (arr->At(0) ? arr->At(0)->GetName() : ""), nXbins, xbins, nYbins, ybins);
          fBinsAllocated += (nXbins + 2) * (nYbins + 2);
          h->Sumw2();
          if (titleStr.Contains("--s--")) {
            (reinterpret_cast<TProfile2D*>(h))->BuildOptions(0., 0., "s");
          }
        }
      } else {
        if (!isdouble) {
          h = new TH3F(hname, (arr->At(0) ? arr->At(0)->GetName() : ""), nXbins, xbins, nYbins, ybins, nZbins, zbins);
        } else {
          h = new TH3D(hname, (arr->At(0) ? arr->At(0)->GetName() : ""), nXbins, xbins, nYbins, ybins, nZbins, zbins);
        }
        fBinsAllocated += (nXbins + 2) * (nYbins + 2) * (nZbins + 2);
        h->Sumw2();
      }
      if (fVariableNames[varX][0]) {
        h->GetXaxis()->SetTitle(Form("%s %s", fVariableNames[varX].Data(),
                                     (fVariableUnits[varX][0] ? Form("(%s)", fVariableUnits[varX].Data()) : "")));
      }
      if (arr->At(1)) {
        h->GetXaxis()->SetTitle(arr->At(1)->GetName());
      }
      if (xLabels[0] != '\0') {
        MakeAxisLabels(h->GetXaxis(), xLabels);
      }
      if (fVariableNames[varY][0]) {
        h->GetYaxis()->SetTitle(Form("%s %s", fVariableNames[varY].Data(),
                                     (fVariableUnits[varY][0] ? Form("(%s)", fVariableUnits[varY].Data()) : "")));
      }
      if (arr->At(2)) {
        h->GetYaxis()->SetTitle(arr->At(2)->GetName());
      }
      if (yLabels[0] != '\0') {
        MakeAxisLabels(h->GetYaxis(), yLabels);
      }
      if (fVariableNames[varZ][0]) {
        h->GetZaxis()->SetTitle(Form("%s %s", fVariableNames[varZ].Data(),
                                     (fVariableUnits[varZ][0] ? Form("(%s)", fVariableUnits[varZ].Data()) : "")));
      }
      if (fVariableNames[varZ][0] && isProfile && varT < 0) { // TProfile2D
        h->GetZaxis()->SetTitle(Form("<%s> %s", fVariableNames[varZ].Data(),
                                     (fVariableUnits[varZ][0] ? Form("(%s)", fVariableUnits[varZ].Data()) : "")));
      }

      if (arr->At(3)) {
        h->GetZaxis()->SetTitle(arr->At(3)->GetName());
      }
      if (zLabels[0] != '\0') {
        MakeAxisLabels(h->GetZaxis(), zLabels);
      }
      hList->Add(h);
      break;
  } // end switch(dimension)
}

//_________________________________________________________________
void HistogramManager::AddHistogram(const char* histClass, const char* hname, const char* title,
                                    int nDimensions, int* vars, int* nBins, double* xmin, double* xmax,
                                    TString* axLabels, int varW, bool useSparse, bool isdouble)
{
  //
  // add a multi-dimensional histogram THnF or THnFSparseF
  //

  // get the list to which the histogram should be added
  auto* hList = reinterpret_cast<TList*>(fMainList->FindObject(histClass));
  if (!hList) {
    LOG(warn) << "HistogramManager::AddHistogram(): Histogram list " << histClass << " not found!";
    LOG(warn) << "         Histogram not created";
    return;
  }
  // check whether this histogram name was used before
  if (hList->FindObject(hname)) {
    LOG(warn) << "HistogramManager::AddHistogram(): Histogram " << hname << " already exists";
    return;
  }

  // tokenize the title string; the user may include in it axis titles which will overwrite the defaults
  TString titleStr(title);
  std::unique_ptr<TObjArray> arr(titleStr.Tokenize(";"));

  if (varW > kNothing) {
    fUsedVars[varW] = kTRUE;
  }

  for (int i = 0; i < nDimensions; i++) {
    if (xmax[i] <= xmin[i]) {
      LOG(warn) << "HistogramManager::AddHistogram(): Histogram " << hname << " has wrong axes ranges for dimension " << i
                << ", (xmin/xmax): " << xmin[i] << " / " << xmax[i];
    }
  }

  // encode needed variable identifiers in a vector and push it to the std::list corresponding to the current histogram list
  std::vector<int> varVector;
  varVector.push_back(0);           // whether the histogram is a profile
  varVector.push_back(nDimensions); // number of dimensions
  varVector.push_back(varW);        // variable used for weighting
  for (int idim = 0; idim < nDimensions; ++idim) {
    varVector.push_back(vars[idim]); // axes variables
  }
  std::list varList = fVariablesMap[histClass];
  varList.push_back(varVector);
  fVariablesMap[histClass] = varList;

  uint32_t nbins = 1;
  THnBase* h = nullptr;
  if (!isdouble) {
    if (useSparse) {
      h = new THnSparseF(hname, (arr->At(0) ? arr->At(0)->GetName() : ""), nDimensions, nBins, xmin, xmax);
    } else {
      h = new THnF(hname, (arr->At(0) ? arr->At(0)->GetName() : ""), nDimensions, nBins, xmin, xmax);
    }
  } else {
    if (useSparse) {
      h = new THnSparseD(hname, (arr->At(0) ? arr->At(0)->GetName() : ""), nDimensions, nBins, xmin, xmax);
    } else {
      h = new THnD(hname, (arr->At(0) ? arr->At(0)->GetName() : ""), nDimensions, nBins, xmin, xmax);
    }
  }
  h->Sumw2();

  // configure the THn histogram and count the allocated bins
  for (int idim = 0; idim < nDimensions; ++idim) {
    nbins *= (nBins[idim] + 2);
    TAxis* axis = h->GetAxis(idim);
    if (fVariableNames[vars[idim]][0]) {
      axis->SetTitle(Form("%s %s", fVariableNames[vars[idim]].Data(),
                          (fVariableUnits[vars[idim]][0] ? Form("(%s)", fVariableUnits[vars[idim]].Data()) : "")));
    }
    if (arr->At(1 + idim)) {
      axis->SetTitle(arr->At(1 + idim)->GetName());
    }
    if (axLabels && !axLabels[idim].IsNull()) {
      MakeAxisLabels(axis, axLabels[idim].Data());
    }

    fUsedVars[vars[idim]] = kTRUE;
  }
  if (!isdouble) {
    if (useSparse) {
      hList->Add(reinterpret_cast<THnSparseF*>(h));
    } else {
      hList->Add(reinterpret_cast<THnF*>(h));
    }
  } else {
    if (useSparse) {
      hList->Add(reinterpret_cast<THnSparseD*>(h));
    } else {
      hList->Add(reinterpret_cast<THnD*>(h));
    }
  }

  fBinsAllocated += nbins;
}

//_________________________________________________________________
void HistogramManager::AddHistogram(const char* histClass, const char* hname, const char* title,
                                    int nDimensions, int* vars, TArrayD* binLimits,
                                    TString* axLabels, int varW, bool useSparse, bool isdouble)
{
  //
  // add a multi-dimensional histogram THnF or THnSparseF with equal or variable bin widths
  //

  // get the list to which the histogram should be added
  auto* hList = reinterpret_cast<TList*>(fMainList->FindObject(histClass));
  if (!hList) {
    LOG(warn) << "HistogramManager::AddHistogram(): Histogram list " << histClass << " not found!";
    LOG(warn) << "         Histogram not created";
    return;
  }
  // check whether this histogram name was used before
  if (hList->FindObject(hname)) {
    LOG(warn) << "HistogramManager::AddHistogram(): Histogram " << hname << " already exists";
    return;
  }

  // tokenize the title string; the user may include in it axis titles which will overwrite the defaults
  TString titleStr(title);
  std::unique_ptr<TObjArray> arr(titleStr.Tokenize(";"));

  if (varW > kNothing) {
    fUsedVars[varW] = kTRUE;
  }

  // encode needed variable identifiers in a vector and push it to the std::list corresponding to the current histogram list
  std::vector<int> varVector;
  varVector.push_back(0);           // whether the histogram is a profile
  varVector.push_back(nDimensions); // number of dimensions
  varVector.push_back(varW);        // variable used for weighting
  for (int idim = 0; idim < nDimensions; ++idim) {
    varVector.push_back(vars[idim]); // axes variables
  }
  std::list varList = fVariablesMap[histClass];
  varList.push_back(varVector);
  fVariablesMap[histClass] = varList;

  // get the min and max for each axis
  auto* xmin = new double[nDimensions];
  auto* xmax = new double[nDimensions];
  int* nBins = new int[nDimensions];
  for (int idim = 0; idim < nDimensions; ++idim) {
    nBins[idim] = binLimits[idim].GetSize() - 1;
    xmin[idim] = binLimits[idim][0];
    xmax[idim] = binLimits[idim][nBins[idim]];
  }

  // initialize the THn with equal spaced bins
  THnBase* h = nullptr;
  if (!isdouble) {
    if (useSparse) {
      h = new THnSparseF(hname, (arr->At(0) ? arr->At(0)->GetName() : ""), nDimensions, nBins, xmin, xmax);
    } else {
      h = new THnF(hname, (arr->At(0) ? arr->At(0)->GetName() : ""), nDimensions, nBins, xmin, xmax);
    }
  } else {
    if (useSparse) {
      h = new THnSparseD(hname, (arr->At(0) ? arr->At(0)->GetName() : ""), nDimensions, nBins, xmin, xmax);
    } else {
      h = new THnD(hname, (arr->At(0) ? arr->At(0)->GetName() : ""), nDimensions, nBins, xmin, xmax);
    }
  }
  // rebin the axes according to the user requested binning
  for (int idim = 0; idim < nDimensions; ++idim) {
    TAxis* axis = h->GetAxis(idim);
    axis->Set(nBins[idim], binLimits[idim].GetArray());
  }
  h->Sumw2();

  uint32_t bins = 1;
  for (int idim = 0; idim < nDimensions; ++idim) {
    bins *= (nBins[idim] + 2);
    TAxis* axis = h->GetAxis(idim);
    if (fVariableNames[vars[idim]][0]) {
      axis->SetTitle(Form("%s %s", fVariableNames[vars[idim]].Data(),
                          (fVariableUnits[vars[idim]][0] ? Form("(%s)", fVariableUnits[vars[idim]].Data()) : "")));
    }
    if (arr->At(1 + idim)) {
      axis->SetTitle(arr->At(1 + idim)->GetName());
    }
    if (axLabels && !axLabels[idim].IsNull()) {
      MakeAxisLabels(axis, axLabels[idim].Data());
    }
    fUsedVars[vars[idim]] = kTRUE;
  }
  if (!isdouble) {
    if (useSparse) {
      hList->Add(reinterpret_cast<THnSparseF*>(h));
    } else {
      hList->Add(reinterpret_cast<THnF*>(h));
    }
  } else {
    if (useSparse) {
      hList->Add(reinterpret_cast<THnSparseD*>(h));
    } else {
      hList->Add(reinterpret_cast<THnD*>(h));
    }
  }
  fBinsAllocated += bins;
}

//__________________________________________________________________
void HistogramManager::FillHistClass(const char* className, Float_t* values)
{
  //
  //  fill a class of histograms
  //

  // get the needed histogram list
  auto* hList = reinterpret_cast<TList*>(fMainList->FindObject(className));
  if (!hList) {
    // TODO: add some meaningfull error message
    /*LOG(warn) << "HistogramManager::FillHistClass(): Histogram list " << className << " not found!";
    LOG(warn) << "         Histogram list not filled" << endl; */
    return;
  }

  // get the corresponding std::list containng identifiers to the needed variables to be filled
  auto const& varList = fVariablesMap[className];

  TIter next(hList);

  TObject* h = nullptr;
  bool isProfile;
  bool isTHn;
  int dimension = 0;
  bool isSparse = kFALSE;
  bool isFillLabelx = kFALSE;
  // TODO: At the moment, maximum 20 dimensions are foreseen for the THn histograms. We should make this more dynamic
  //       But maybe its better to have it like to avoid dynamically allocating this array in the histogram loop
  double fillValues[20] = {0.0};
  int varX = -1, varY = -1, varZ = -1, varT = -1, varW = -1;

  // loop over the histogram and std::list
  // NOTE: these two should contain the same number of elements and be synchronized, otherwise its a mess
  for (auto varIter = varList.begin(); varIter != varList.end(); varIter++) {
    h = next(); // get the histogram
    // decode information from the vector of indices
    isProfile = ((*varIter)[0] == 1 ? true : false);
    isTHn = ((*varIter)[1] > 0 ? true : false);
    if (isTHn) {
      dimension = (*varIter)[1];
    } else {
      dimension = (reinterpret_cast<TH1*>(h))->GetDimension();
    }

    // get the various variable indices
    varW = (*varIter)[2];
    if (isTHn) {
      for (int i = 0; i < dimension; i++) {
        fillValues[i] = values[(*varIter)[3 + i]];
      }
    } else {
      varX = (*varIter)[3];
      varY = (*varIter)[4];
      varZ = (*varIter)[5];
      varT = (*varIter)[6];
      isFillLabelx = ((*varIter)[7] == 1 ? true : false);
    }

    if (!isTHn) {
      switch (dimension) {
        case 1:
          if (isProfile) {
            if (varW > kNothing) {
              if (isFillLabelx) {
                (reinterpret_cast<TProfile*>(h))->Fill(Form("%d", static_cast<int>(values[varX])), values[varY], values[varW]);
              } else {
                (reinterpret_cast<TProfile*>(h))->Fill(values[varX], values[varY], values[varW]);
              }
            } else {
              if (isFillLabelx) {
                (reinterpret_cast<TProfile*>(h))->Fill(Form("%d", static_cast<int>(values[varX])), values[varY]);
              } else {
                (reinterpret_cast<TProfile*>(h))->Fill(values[varX], values[varY]);
              }
            }
          } else {
            if (varW > kNothing) {
              if (isFillLabelx) {
                (reinterpret_cast<TH1*>(h))->Fill(Form("%d", static_cast<int>(values[varX])), values[varW]);
              } else {
                (reinterpret_cast<TH1*>(h))->Fill(values[varX], values[varW]);
              }
            } else {
              if (isFillLabelx) {
                (reinterpret_cast<TH1*>(h))->Fill(Form("%d", static_cast<int>(values[varX])), 1.);
              } else {
                (reinterpret_cast<TH1*>(h))->Fill(values[varX]);
              }
            }
          }
          break;
        case 2:
          if (isProfile) {
            if (varW > kNothing) {
              (reinterpret_cast<TProfile2D*>(h))->Fill(values[varX], values[varY], values[varZ], values[varW]);
            } else {
              (reinterpret_cast<TProfile2D*>(h))->Fill(values[varX], values[varY], values[varZ]);
            }
          } else {
            if (varW > kNothing) {
              if (isFillLabelx) {
                (reinterpret_cast<TH2*>(h))->Fill(Form("%d", static_cast<int>(values[varX])), values[varY], values[varW]);
              } else {
                (reinterpret_cast<TH2*>(h))->Fill(values[varX], values[varY], values[varW]);
              }
            } else {
              if (isFillLabelx) {
                (reinterpret_cast<TH2*>(h))->Fill(Form("%d", static_cast<int>(values[varX])), values[varY], 1.);
              } else {
                (reinterpret_cast<TH2*>(h))->Fill(values[varX], values[varY]);
              }
            }
          }
          break;
        case 3:
          if (isProfile) {
            if (varW > kNothing) {
              (reinterpret_cast<TProfile3D*>(h))->Fill(values[varX], values[varY], values[varZ], values[varT], values[varW]);
            } else {
              (reinterpret_cast<TProfile3D*>(h))->Fill(values[varX], values[varY], values[varZ], values[varT]);
            }
          } else {
            if (varW > kNothing) {
              (reinterpret_cast<TH3*>(h))->Fill(values[varX], values[varY], values[varZ], values[varW]);
            } else {
              (reinterpret_cast<TH3*>(h))->Fill(values[varX], values[varY], values[varZ]);
            }
          }
          break;

        default:
          break;
      } // end switch
      // end if(!isTHn)
    } else {
      if (varW > kNothing) {
        if (isSparse) {
          (reinterpret_cast<THnSparse*>(h))->Fill(fillValues, values[varW]);
        } else {
          (reinterpret_cast<THn*>(h))->Fill(fillValues, values[varW]);
        }
      } else {
        if (isSparse) {
          (reinterpret_cast<THnSparse*>(h))->Fill(fillValues);
        } else {
          (reinterpret_cast<THn*>(h))->Fill(fillValues);
        }
      }
    } // end else
  } // end loop over histograms
}

//____________________________________________________________________________________
void HistogramManager::MakeAxisLabels(TAxis* ax, const char* labels)
{
  //
  // add bin labels to an axis
  //
  TString labelsStr(labels);
  std::unique_ptr<TObjArray> arr(labelsStr.Tokenize(";"));
  for (int ib = 1; ib <= ax->GetNbins(); ++ib) {
    if (ib >= arr->GetEntries() + 1) {
      break;
    }
    ax->SetBinLabel(ib, arr->At(ib - 1)->GetName());
  }
}

//____________________________________________________________________________________
void HistogramManager::Print(Option_t*) const
{
  //
  // Print the defined histograms
  //
  cout << "###################################################################" << endl;
  cout << "HistogramManager:: " << fMainList->GetName() << endl;
  for (int i = 0; i < fMainList->GetEntries(); ++i) {
    auto* list = reinterpret_cast<TList*>(fMainList->At(i));
    cout << "************** List " << list->GetName() << endl;
    for (int j = 0; j < list->GetEntries(); ++j) {
      TObject* obj = list->At(j);
      cout << obj->GetName() << ": " << obj->IsA()->GetName() << endl;
    }
  }
}
