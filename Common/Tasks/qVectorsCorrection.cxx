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
/// \file   qVectorsCorrection.cxx
/// \author Cindy Mordasini <cindy.mordasini@cern.ch>
/// \author Anna Önnerstad <anna.onnerstad@cern.ch>
///
/// \brief  ...
///

// C++/ROOT includes.
#include <chrono>
#include <string>
#include <vector>
#include <TComplex.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TMath.h>

// o2Physics includes.
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/StaticFor.h"

#include "Common/DataModel/Qvectors.h"
#include "Common/Core/EventPlaneHelper.h"

// o2 includes.

using namespace o2;
using namespace o2::framework;

namespace qV
{
static constexpr std::string_view centClasses[8] = {
  "Centrality_0-5/", "Centrality_5-10/", "Centrality_10-20/", "Centrality_20-30/",
  "Centrality_30-40/", "Centrality_40-50/", "Centrality_50-60/", "Centrality_60-80/"};
} // namespace qV

struct qVectorsCorrection {
  // No correction = recenter, recentered Qvectors = twist, twisted Qvectors = rescale.
  // NOTE: As of no, the twist gets both twist and rescale correction constants.

  // Histogram registry for the output QA figures and list of centrality classes for it.
  // Objects are NOT saved in alphabetical orders, and registry names are NOT saved
  // as TDirectoryFile.
  HistogramRegistry histosQA{"histosQA", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  Configurable<int> cfgMinTPCTracks{"cfgMinTPCTracks", 20, "minimum TPC tracks participating in Q-vector reconstruction"};
  Configurable<int> cfgnMod{"cfgnMod", 2, "Modulation of interest"};

  // Helper variables.
  EventPlaneHelper helperEP;

  void init(InitContext const&)
  {
    // Fill the registry with the needed objects.
    const AxisSpec axisCent{110, 0., 110.};
    const AxisSpec axisQvec{1000, -5, 5};
    const AxisSpec axisConst{12, 0., 12.}; // 4 constants x 3 detectors.
    const AxisSpec axisEvtPl{360, -constants::math::PI, constants::math::PI};

    histosQA.add("histCentFull", "Centrality distribution for valid events",
                 HistType::kTH1F, {axisCent});

    histosQA.add("Centrality_0-5/histCent", "Centrality distribution",
                 HistType::kTH1F, {axisCent});

    histosQA.add("Centrality_0-5/histQvecFT0CUncor", "", {HistType::kTH2F, {axisQvec, axisQvec}});
    histosQA.add("Centrality_0-5/histQvecFT0CRectr", "", {HistType::kTH2F, {axisQvec, axisQvec}});
    histosQA.add("Centrality_0-5/histQvecFT0CTwist", "", {HistType::kTH2F, {axisQvec, axisQvec}});
    histosQA.add("Centrality_0-5/histQvecFT0CFinal", "", {HistType::kTH2F, {axisQvec, axisQvec}});

    histosQA.add("Centrality_0-5/histQvecBPosUncor", "", {HistType::kTH2F, {axisQvec, axisQvec}});
    histosQA.add("Centrality_0-5/histQvecBPosRectr", "", {HistType::kTH2F, {axisQvec, axisQvec}});
    histosQA.add("Centrality_0-5/histQvecBPosTwist", "", {HistType::kTH2F, {axisQvec, axisQvec}});
    histosQA.add("Centrality_0-5/histQvecBPosFinal", "", {HistType::kTH2F, {axisQvec, axisQvec}});

    histosQA.add("Centrality_0-5/histQvecBNegUncor", "", {HistType::kTH2F, {axisQvec, axisQvec}});
    histosQA.add("Centrality_0-5/histQvecBNegRectr", "", {HistType::kTH2F, {axisQvec, axisQvec}});
    histosQA.add("Centrality_0-5/histQvecBNegTwist", "", {HistType::kTH2F, {axisQvec, axisQvec}});
    histosQA.add("Centrality_0-5/histQvecBNegFinal", "", {HistType::kTH2F, {axisQvec, axisQvec}});

    histosQA.add("Centrality_0-5/histEvtPlFT0CUncor", "", {HistType::kTH1F, {axisEvtPl}});
    histosQA.add("Centrality_0-5/histEvtPlFT0CRectr", "", {HistType::kTH1F, {axisEvtPl}});
    histosQA.add("Centrality_0-5/histEvtPlFT0CTwist", "", {HistType::kTH1F, {axisEvtPl}});
    histosQA.add("Centrality_0-5/histEvtPlFT0CFinal", "", {HistType::kTH1F, {axisEvtPl}});

    histosQA.add("Centrality_0-5/histEvtPlBPosUncor", "", {HistType::kTH1F, {axisEvtPl}});
    histosQA.add("Centrality_0-5/histEvtPlBPosRectr", "", {HistType::kTH1F, {axisEvtPl}});
    histosQA.add("Centrality_0-5/histEvtPlBPosTwist", "", {HistType::kTH1F, {axisEvtPl}});
    histosQA.add("Centrality_0-5/histEvtPlBPosFinal", "", {HistType::kTH1F, {axisEvtPl}});

    histosQA.add("Centrality_0-5/histEvtPlBNegUncor", "", {HistType::kTH1F, {axisEvtPl}});
    histosQA.add("Centrality_0-5/histEvtPlBNegRectr", "", {HistType::kTH1F, {axisEvtPl}});
    histosQA.add("Centrality_0-5/histEvtPlBNegTwist", "", {HistType::kTH1F, {axisEvtPl}});
    histosQA.add("Centrality_0-5/histEvtPlBNegFinal", "", {HistType::kTH1F, {axisEvtPl}});

    histosQA.add("Centrality_0-5/histEvtPlResolution", "", {HistType::kTH1F, {axisEvtPl}});

    for (int iBin = 1; iBin < 8; iBin++) {
      histosQA.addClone("Centrality_0-5/", qV::centClasses[iBin].data());
    }

  } // End void init(InitContext const&)

  // Definition of all the needed template functions.
  template <int cBin, typename T>
  void fillHistosQvec(const T& vec)
  {
    histosQA.fill(HIST(qV::centClasses[cBin]) + HIST("histCent"), vec.cent());

    if (vec.sumAmplFT0C() > 1e-8) {
      histosQA.fill(HIST(qV::centClasses[cBin]) + HIST("histQvecFT0CUncor"), vec.qvecFT0CUncorRe(), vec.qvecFT0CUncorIm());
      histosQA.fill(HIST(qV::centClasses[cBin]) + HIST("histQvecFT0CRectr"), vec.qvecFT0CRectrRe(), vec.qvecFT0CRectrIm());
      histosQA.fill(HIST(qV::centClasses[cBin]) + HIST("histQvecFT0CTwist"), vec.qvecFT0CTwistRe(), vec.qvecFT0CTwistIm());
      histosQA.fill(HIST(qV::centClasses[cBin]) + HIST("histQvecFT0CFinal"), vec.qvecFT0CRe(), vec.qvecFT0CIm());

      histosQA.fill(HIST(qV::centClasses[cBin]) + HIST("histEvtPlFT0CUncor"), helperEP.GetEventPlane(vec.qvecFT0CUncorRe(), vec.qvecFT0CUncorIm()));
      histosQA.fill(HIST(qV::centClasses[cBin]) + HIST("histEvtPlFT0CRectr"), helperEP.GetEventPlane(vec.qvecFT0CRectrRe(), vec.qvecFT0CRectrIm()));
      histosQA.fill(HIST(qV::centClasses[cBin]) + HIST("histEvtPlFT0CTwist"), helperEP.GetEventPlane(vec.qvecFT0CTwistRe(), vec.qvecFT0CTwistIm()));
      histosQA.fill(HIST(qV::centClasses[cBin]) + HIST("histEvtPlFT0CFinal"), helperEP.GetEventPlane(vec.qvecFT0CRe(), vec.qvecFT0CIm()));
    }

    if (vec.nTrkBPos() > cfgMinTPCTracks) {
      histosQA.fill(HIST(qV::centClasses[cBin]) + HIST("histQvecBPosUncor"), vec.qvecBPosUncorRe(), vec.qvecBPosUncorIm());
      histosQA.fill(HIST(qV::centClasses[cBin]) + HIST("histQvecBPosRectr"), vec.qvecBPosRectrRe(), vec.qvecBPosRectrIm());
      histosQA.fill(HIST(qV::centClasses[cBin]) + HIST("histQvecBPosTwist"), vec.qvecBPosTwistRe(), vec.qvecBPosTwistIm());
      histosQA.fill(HIST(qV::centClasses[cBin]) + HIST("histQvecBPosFinal"), vec.qvecBPosRe(), vec.qvecBPosIm());

      histosQA.fill(HIST(qV::centClasses[cBin]) + HIST("histEvtPlBPosUncor"), helperEP.GetEventPlane(vec.qvecBPosUncorRe(), vec.qvecBPosUncorIm()));
      histosQA.fill(HIST(qV::centClasses[cBin]) + HIST("histEvtPlBPosRectr"), helperEP.GetEventPlane(vec.qvecBPosRectrRe(), vec.qvecBPosRectrIm()));
      histosQA.fill(HIST(qV::centClasses[cBin]) + HIST("histEvtPlBPosTwist"), helperEP.GetEventPlane(vec.qvecBPosTwistRe(), vec.qvecBPosTwistIm()));
      histosQA.fill(HIST(qV::centClasses[cBin]) + HIST("histEvtPlBPosFinal"), helperEP.GetEventPlane(vec.qvecBPosRe(), vec.qvecBPosIm()));
    }

    if (vec.nTrkBNeg() > cfgMinTPCTracks) {
      histosQA.fill(HIST(qV::centClasses[cBin]) + HIST("histQvecBNegUncor"), vec.qvecBNegUncorRe(), vec.qvecBNegUncorIm());
      histosQA.fill(HIST(qV::centClasses[cBin]) + HIST("histQvecBNegRectr"), vec.qvecBNegRectrRe(), vec.qvecBNegRectrIm());
      histosQA.fill(HIST(qV::centClasses[cBin]) + HIST("histQvecBNegTwist"), vec.qvecBNegTwistRe(), vec.qvecBNegTwistIm());
      histosQA.fill(HIST(qV::centClasses[cBin]) + HIST("histQvecBNegFinal"), vec.qvecBNegRe(), vec.qvecBNegIm());

      histosQA.fill(HIST(qV::centClasses[cBin]) + HIST("histEvtPlBNegUncor"), helperEP.GetEventPlane(vec.qvecBNegUncorRe(), vec.qvecBNegUncorIm()));
      histosQA.fill(HIST(qV::centClasses[cBin]) + HIST("histEvtPlBNegRectr"), helperEP.GetEventPlane(vec.qvecBNegRectrRe(), vec.qvecBNegRectrIm()));
      histosQA.fill(HIST(qV::centClasses[cBin]) + HIST("histEvtPlBNegTwist"), helperEP.GetEventPlane(vec.qvecBNegTwistRe(), vec.qvecBNegTwistIm()));
      histosQA.fill(HIST(qV::centClasses[cBin]) + HIST("histEvtPlBNegFinal"), helperEP.GetEventPlane(vec.qvecBNegRe(), vec.qvecBNegIm()));
    }

    if (vec.nTrkBPos() > cfgMinTPCTracks && vec.nTrkBNeg() > cfgMinTPCTracks && vec.sumAmplFT0C() > 1e-8) {
      histosQA.fill(HIST(qV::centClasses[cBin]) + HIST("histEvtPlResolution"), helperEP.GetResolution(
                                                                                 helperEP.GetEventPlane(vec.qvecBPosRe(), vec.qvecBPosIm()),
                                                                                 helperEP.GetEventPlane(vec.qvecBNegRe(), vec.qvecBNegIm()),
                                                                                 helperEP.GetEventPlane(vec.qvecBPosRe(), vec.qvecBPosIm()), 2));
    }
  }

  void process(aod::Qvector const& qVec)
  {
    // Get the centrality bin, and fill the centrality QA histograms.
    int centBin = helperEP.GetCentBin(qVec.cent());
    histosQA.fill(HIST("histCentFull"), qVec.cent());
    switch (centBin) { // LOKI: See if we can replace that with a const something like below.
      case 0:
        fillHistosQvec<0>(qVec);
        break;
      case 1:
        fillHistosQvec<1>(qVec);
        break;
      case 2:
        fillHistosQvec<2>(qVec);
        break;
      case 3:
        fillHistosQvec<3>(qVec);
        break;
      case 4:
        fillHistosQvec<4>(qVec);
        break;
      case 5:
        fillHistosQvec<5>(qVec);
        break;
      case 6:
        fillHistosQvec<6>(qVec);
        break;
      case 7:
        fillHistosQvec<7>(qVec);
        break;
    }   // End switch(centBin)
  }     // End void process(...)
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<qVectorsCorrection>(cfgc)};
}
