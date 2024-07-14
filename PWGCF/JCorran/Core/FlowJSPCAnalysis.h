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

// \brief   Calculation class for the SPC-related analyses.
// \author  Maxim Virta (maxim.virta@cern.fi), Cindy Mordasini (cindy.mordasini@cern.ch)

#ifndef PWGCF_JCORRAN_CORE_FLOWJSPCANALYSIS_H_
#define PWGCF_JCORRAN_CORE_FLOWJSPCANALYSIS_H_

/* Header files. */
#include <iostream>
#include <array>
#include <vector>
#include <TComplex.h>
#include <TProfile.h>

// O2 headers. //
#include "Framework/HistogramRegistry.h"
#include "PWGCF/JCorran/Core/JQVectors.h"

using namespace o2;
using namespace o2::framework;
using namespace std;

class FlowJSPCAnalysis
{
 public:
  FlowJSPCAnalysis() = default;

  void SetHistRegistry(HistogramRegistry* histReg) { mHistRegistry = histReg; }
  Int_t GetCentBin(float cValue);

  using JQVectorsT = JQVectors<TComplex, 113, 15, false>;
  inline void SetQvectors(const JQVectorsT* _qvecs) { qvecs = _qvecs; }
  void Correlation(Int_t c_nPart, Int_t c_nHarmo, Int_t* harmo, Double_t* correlData);
  void CalculateCorrelators(const Int_t fCentBin);
  void FillHistograms(const Int_t fCentBin, Int_t ind, Double_t cNum, Double_t cDenom, Double_t wNum, Double_t wDenom);
  void FillQAHistograms(const Int_t fCentBin, Double_t phi, Double_t phiWeight);
  TComplex Recursion(int n, int* harmonic, int mult, int skip);
  TComplex Q(const Int_t harmN, const Int_t p);

  void CreateHistos()
  {
    if (!mHistRegistry) {
      LOGF(error, "QA histogram registry missing. Quitting...");
      return;
    }
    mHistRegistry->add("FullCentrality", "FullCentrality", HistType::kTH1D, {{100, 0., 100.}}, true);
    mHistRegistry->add("Centrality_0/fResults", "Numerators and denominators", {HistType::kTProfile, {{24, 0., 24.}}}, true);
    mHistRegistry->add("Centrality_0/fCovResults", "Covariance N*D", {HistType::kTProfile, {{48, 0., 48.}}}, true);
    mHistRegistry->add("Centrality_0/phiBefore", "Phi before", {HistType::kTH1D, {{100, 0., TMath::TwoPi()}}}, true);
    mHistRegistry->add("Centrality_0/phiAfter", "Phi after", {HistType::kTH1D, {{100, 0., TMath::TwoPi()}}}, true);

    for (UInt_t i = 1; i < 8; i++) {
      mHistRegistry->addClone("Centrality_0/", Form("Centrality_%u/", i));
    }
  }

  void SetCorrSet(Int_t obsInd, Int_t harmo[8])
  {
    for (int i = 0; i < 8; i++) {
      fHarmosArray[obsInd][i] = harmo[i];
    }
  }
  void SetFullCorrSet(Int_t harmo[12][8])
  {
    memcpy(fHarmosArray, harmo, sizeof(Int_t) * 12 * 8);
  }

  static constexpr std::string_view mCentClasses[] = {
    "Centrality_0/",
    "Centrality_1/",
    "Centrality_2/",
    "Centrality_3/",
    "Centrality_4/",
    "Centrality_5/",
    "Centrality_6/",
    "Centrality_7/",
    "Centrality_8/",
    "Centrality_9/"};

 private:
  const Int_t mNqHarmos = 113; ///< Highest harmo for Q(n,p): (v8*14part)+1.
  const Int_t mNqPowers = 15;  ///< Max power for Q(n,p): 14part+1.
  const JQVectorsT* qvecs;

  HistogramRegistry* mHistRegistry = nullptr;

  Int_t fHarmosArray[12][8];

  Double_t fCorrelDenoms[14];

  ClassDefNV(FlowJSPCAnalysis, 1);
};
#endif // PWGCF_JCORRAN_CORE_FLOWJSPCANALYSIS_H_
