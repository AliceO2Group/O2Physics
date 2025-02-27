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
#include <array>
#include <vector>
#include <TComplex.h>
#include <TProfile.h>

// O2 headers. //
#include "Framework/HistogramRegistry.h"
#include "PWGCF/JCorran/Core/JQVectors.h"
#include "CommonConstants/MathConstants.h"

class FlowJSPCAnalysis
{
 public:
  FlowJSPCAnalysis() = default;

  void setHistRegistry(o2::framework::HistogramRegistry* histReg) { mHistRegistry = histReg; }
  int getCentBin(float cValue);

  using JQVectorsT = JQVectors<TComplex, 113, 15, false>;
  inline void setQvectors(const JQVectorsT* _qvecs) { qvecs = _qvecs; }
  void correlation(int c_nPart, int c_nHarmo, int* harmo, double* correlData);
  void calculateCorrelators(const int fCentBin);
  void fillHistograms(const int fCentBin, int ind, double cNum, double cDenom, double wNum, double wDenom);
  void fillQAHistograms(const int fCentBin, double phi, double phiWeight);
  TComplex recursion(int n, int* harmonic, int mult, int skip);
  TComplex q(const int harmN, const int p);

  void createHistos()
  {
    if (!mHistRegistry) {
      LOGF(error, "QA histogram registry missing. Quitting...");
      return;
    }
    mHistRegistry->add("FullCentrality", "FullCentrality", o2::framework::HistType::kTH1D, {{100, 0., 100.}}, true);
    mHistRegistry->add("Centrality_0/fResults", "Numerators and denominators", {o2::framework::HistType::kTProfile, {{24, 0., 24.}}}, true);
    mHistRegistry->add("Centrality_0/fCovResults", "Covariance N*D", {o2::framework::HistType::kTProfile, {{48, 0., 48.}}}, true);
    mHistRegistry->add("Centrality_0/phiBefore", "Phi before", {o2::framework::HistType::kTH1D, {{100, 0., o2::constants::math::TwoPI}}}, true);
    mHistRegistry->add("Centrality_0/phiAfter", "Phi after", {o2::framework::HistType::kTH1D, {{100, 0., o2::constants::math::TwoPI}}}, true);

    for (uint i = 1; i < 9; i++) {
      mHistRegistry->addClone("Centrality_0/", Form("Centrality_%u/", i));
    }
  }

  void setCorrSet(int obsInd, int harmo[8])
  {
    for (int i = 0; i < 8; i++) {
      fHarmosArray[obsInd][i] = harmo[i];
    }
  }
  void setFullCorrSet(int harmo[12][8])
  {
    memcpy(fHarmosArray, harmo, sizeof(int) * 12 * 8);
  }

  static constexpr std::string_view MCentClasses[] = {
    "Centrality_0/",
    "Centrality_1/",
    "Centrality_2/",
    "Centrality_3/",
    "Centrality_4/",
    "Centrality_5/",
    "Centrality_6/",
    "Centrality_7/",
    "Centrality_8/"};

 private:
  const int mNqHarmos = 113; ///< Highest harmo for Q(n,p): (v8*14part)+1.
  const int mNqPowers = 15;  ///< Max power for Q(n,p): 14part+1.
  const JQVectorsT* qvecs;

  o2::framework::HistogramRegistry* mHistRegistry = nullptr;

  int fHarmosArray[12][8];

  double fCorrelDenoms[14];

  ClassDefNV(FlowJSPCAnalysis, 1);
};
#endif // PWGCF_JCORRAN_CORE_FLOWJSPCANALYSIS_H_
