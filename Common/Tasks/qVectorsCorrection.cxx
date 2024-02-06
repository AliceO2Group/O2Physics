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

  Configurable<int> cfgnMod{"cfgnMod", 2, "Modulation of interest"};

  Configurable<std::string> cfgDetName{"cfgDetName", "FT0C", "The name of detector to be analyzed"};
  Configurable<std::string> cfgRefAName{"cfgRefAName", "BPos", "The name of detector for reference A"};
  Configurable<std::string> cfgRefBName{"cfgRefBName", "BNeg", "The name of detector for reference B"};

  // Helper variables.
  EventPlaneHelper helperEP;

  int DetId;
  int RefAId;
  int RefBId;

  template <typename T>
  int GetDetId(const T& name)
  {
    if (name.value == "FT0C") {
      return 0;
    } else if (name.value == "FT0A") {
      return 1;
    } else if (name.value == "FT0M") {
      return 2;
    } else if (name.value == "FV0A") {
      return 3;
    } else if (name.value == "BPos") {
      return 4;
    } else if (name.value == "BNeg") {
      return 5;
    } else {
      return 0;
    }
  }

  void init(InitContext const&)
  {
    DetId = GetDetId(cfgDetName);
    RefAId = GetDetId(cfgRefAName);
    RefBId = GetDetId(cfgRefBName);

    if (DetId == RefAId || DetId == RefBId || RefAId == RefBId) {
      LOGF(info, "Wrong detector configuration \n The FT0C will be used to get Q-Vector \n The BPos and BNeg will be used as reference systems");
      DetId = 0;
      RefAId = 4;
      RefBId = 5;
    }

    // Fill the registry with the needed objects.
    const AxisSpec axisCent{110, 0., 110.};
    const AxisSpec axisQvec{1000, -3, 3};
    const AxisSpec axisQvecF{1000, -0.3, 0.3};
    const AxisSpec axisConst{12, 0., 12.}; // 4 constants x 3 detectors.
    const AxisSpec axisEvtPl{360, -constants::math::PI, constants::math::PI};

    histosQA.add("histCentFull", "Centrality distribution for valid events",
                 HistType::kTH1F, {axisCent});

    histosQA.add("Centrality_0-5/histCent", "Centrality distribution",
                 HistType::kTH1F, {axisCent});

    histosQA.add("Centrality_0-5/histQvecUncor", "", {HistType::kTH2F, {axisQvecF, axisQvecF}});
    histosQA.add("Centrality_0-5/histQvecRectr", "", {HistType::kTH2F, {axisQvecF, axisQvecF}});
    histosQA.add("Centrality_0-5/histQvecTwist", "", {HistType::kTH2F, {axisQvecF, axisQvecF}});
    histosQA.add("Centrality_0-5/histQvecFinal", "", {HistType::kTH2F, {axisQvec, axisQvec}});

    histosQA.add("Centrality_0-5/histQvecRefAUncor", "", {HistType::kTH2F, {axisQvecF, axisQvecF}});
    histosQA.add("Centrality_0-5/histQvecRefARectr", "", {HistType::kTH2F, {axisQvecF, axisQvecF}});
    histosQA.add("Centrality_0-5/histQvecRefATwist", "", {HistType::kTH2F, {axisQvecF, axisQvecF}});
    histosQA.add("Centrality_0-5/histQvecRefAFinal", "", {HistType::kTH2F, {axisQvec, axisQvec}});

    histosQA.add("Centrality_0-5/histQvecRefBUncor", "", {HistType::kTH2F, {axisQvecF, axisQvecF}});
    histosQA.add("Centrality_0-5/histQvecRefBRectr", "", {HistType::kTH2F, {axisQvecF, axisQvecF}});
    histosQA.add("Centrality_0-5/histQvecRefBTwist", "", {HistType::kTH2F, {axisQvecF, axisQvecF}});
    histosQA.add("Centrality_0-5/histQvecRefBFinal", "", {HistType::kTH2F, {axisQvec, axisQvec}});

    histosQA.add("Centrality_0-5/histEvtPlUncor", "", {HistType::kTH1F, {axisEvtPl}});
    histosQA.add("Centrality_0-5/histEvtPlRectr", "", {HistType::kTH1F, {axisEvtPl}});
    histosQA.add("Centrality_0-5/histEvtPlTwist", "", {HistType::kTH1F, {axisEvtPl}});
    histosQA.add("Centrality_0-5/histEvtPlFinal", "", {HistType::kTH1F, {axisEvtPl}});

    histosQA.add("Centrality_0-5/histEvtPlRefAUncor", "", {HistType::kTH1F, {axisEvtPl}});
    histosQA.add("Centrality_0-5/histEvtPlRefARectr", "", {HistType::kTH1F, {axisEvtPl}});
    histosQA.add("Centrality_0-5/histEvtPlRefATwist", "", {HistType::kTH1F, {axisEvtPl}});
    histosQA.add("Centrality_0-5/histEvtPlRefAFinal", "", {HistType::kTH1F, {axisEvtPl}});

    histosQA.add("Centrality_0-5/histEvtPlRefBUncor", "", {HistType::kTH1F, {axisEvtPl}});
    histosQA.add("Centrality_0-5/histEvtPlRefBRectr", "", {HistType::kTH1F, {axisEvtPl}});
    histosQA.add("Centrality_0-5/histEvtPlRefBTwist", "", {HistType::kTH1F, {axisEvtPl}});
    histosQA.add("Centrality_0-5/histEvtPlRefBFinal", "", {HistType::kTH1F, {axisEvtPl}});

    histosQA.add("Centrality_0-5/histEvtPlRes_SigRefA", "", {HistType::kTH1F, {axisEvtPl}});
    histosQA.add("Centrality_0-5/histEvtPlRes_SigRefB", "", {HistType::kTH1F, {axisEvtPl}});
    histosQA.add("Centrality_0-5/histEvtPlRes_RefARefB", "", {HistType::kTH1F, {axisEvtPl}});

    for (int iBin = 1; iBin < 8; iBin++) {
      histosQA.addClone("Centrality_0-5/", qV::centClasses[iBin].data());
    }

  } // End void init(InitContext const&)

  // Definition of all the needed template functions.
  template <int cBin, typename T>
  void fillHistosQvec(const T& vec)
  {
    histosQA.fill(HIST(qV::centClasses[cBin]) + HIST("histCent"), vec.cent());

    if (vec.qvecAmp()[DetId] > 1e-8) {
      histosQA.fill(HIST(qV::centClasses[cBin]) + HIST("histQvecUncor"), vec.qvecRe()[DetId * 4], vec.qvecIm()[DetId * 4]);
      histosQA.fill(HIST(qV::centClasses[cBin]) + HIST("histQvecRectr"), vec.qvecRe()[DetId * 4 + 1], vec.qvecIm()[DetId * 4 + 1]);
      histosQA.fill(HIST(qV::centClasses[cBin]) + HIST("histQvecTwist"), vec.qvecRe()[DetId * 4 + 2], vec.qvecIm()[DetId * 4 + 2]);
      histosQA.fill(HIST(qV::centClasses[cBin]) + HIST("histQvecFinal"), vec.qvecRe()[DetId * 4 + 3], vec.qvecIm()[DetId * 4 + 3]);

      histosQA.fill(HIST(qV::centClasses[cBin]) + HIST("histEvtPlUncor"), helperEP.GetEventPlane(vec.qvecRe()[DetId * 4], vec.qvecIm()[DetId * 4], cfgnMod));
      histosQA.fill(HIST(qV::centClasses[cBin]) + HIST("histEvtPlRectr"), helperEP.GetEventPlane(vec.qvecRe()[DetId * 4 + 1], vec.qvecIm()[DetId * 4 + 1], cfgnMod));
      histosQA.fill(HIST(qV::centClasses[cBin]) + HIST("histEvtPlTwist"), helperEP.GetEventPlane(vec.qvecRe()[DetId * 4 + 2], vec.qvecIm()[DetId * 4 + 2], cfgnMod));
      histosQA.fill(HIST(qV::centClasses[cBin]) + HIST("histEvtPlFinal"), helperEP.GetEventPlane(vec.qvecRe()[DetId * 4 + 3], vec.qvecIm()[DetId * 4 + 3], cfgnMod));
    }

    if (vec.qvecAmp()[RefAId] > 1e-8) {
      histosQA.fill(HIST(qV::centClasses[cBin]) + HIST("histQvecRefAUncor"), vec.qvecRe()[RefAId * 4], vec.qvecIm()[RefAId * 4]);
      histosQA.fill(HIST(qV::centClasses[cBin]) + HIST("histQvecRefARectr"), vec.qvecRe()[RefAId * 4 + 1], vec.qvecIm()[RefAId * 4 + 1]);
      histosQA.fill(HIST(qV::centClasses[cBin]) + HIST("histQvecRefATwist"), vec.qvecRe()[RefAId * 4 + 2], vec.qvecIm()[RefAId * 4 + 2]);
      histosQA.fill(HIST(qV::centClasses[cBin]) + HIST("histQvecRefAFinal"), vec.qvecRe()[RefAId * 4 + 3], vec.qvecIm()[RefAId * 4 + 3]);

      histosQA.fill(HIST(qV::centClasses[cBin]) + HIST("histEvtPlRefAUncor"), helperEP.GetEventPlane(vec.qvecRe()[RefAId * 4], vec.qvecIm()[RefAId * 4], cfgnMod));
      histosQA.fill(HIST(qV::centClasses[cBin]) + HIST("histEvtPlRefARectr"), helperEP.GetEventPlane(vec.qvecRe()[RefAId * 4 + 1], vec.qvecIm()[RefAId * 4 + 1], cfgnMod));
      histosQA.fill(HIST(qV::centClasses[cBin]) + HIST("histEvtPlRefATwist"), helperEP.GetEventPlane(vec.qvecRe()[RefAId * 4 + 2], vec.qvecIm()[RefAId * 4 + 2], cfgnMod));
      histosQA.fill(HIST(qV::centClasses[cBin]) + HIST("histEvtPlRefAFinal"), helperEP.GetEventPlane(vec.qvecRe()[RefAId * 4 + 3], vec.qvecIm()[RefAId * 4 + 3], cfgnMod));
    }

    if (vec.qvecAmp()[RefBId] > 1e-8) {
      histosQA.fill(HIST(qV::centClasses[cBin]) + HIST("histQvecRefBUncor"), vec.qvecRe()[RefBId * 4], vec.qvecIm()[RefBId * 4]);
      histosQA.fill(HIST(qV::centClasses[cBin]) + HIST("histQvecRefBRectr"), vec.qvecRe()[RefBId * 4 + 1], vec.qvecIm()[RefBId * 4 + 1]);
      histosQA.fill(HIST(qV::centClasses[cBin]) + HIST("histQvecRefBTwist"), vec.qvecRe()[RefBId * 4 + 2], vec.qvecIm()[RefBId * 4 + 2]);
      histosQA.fill(HIST(qV::centClasses[cBin]) + HIST("histQvecRefBFinal"), vec.qvecRe()[RefBId * 4 + 3], vec.qvecIm()[RefBId * 4 + 3]);

      histosQA.fill(HIST(qV::centClasses[cBin]) + HIST("histEvtPlRefBUncor"), helperEP.GetEventPlane(vec.qvecRe()[RefBId * 4], vec.qvecIm()[RefBId * 4], cfgnMod));
      histosQA.fill(HIST(qV::centClasses[cBin]) + HIST("histEvtPlRefBRectr"), helperEP.GetEventPlane(vec.qvecRe()[RefBId * 4 + 1], vec.qvecIm()[RefBId * 4 + 1], cfgnMod));
      histosQA.fill(HIST(qV::centClasses[cBin]) + HIST("histEvtPlRefBTwist"), helperEP.GetEventPlane(vec.qvecRe()[RefBId * 4 + 2], vec.qvecIm()[RefBId * 4 + 2], cfgnMod));
      histosQA.fill(HIST(qV::centClasses[cBin]) + HIST("histEvtPlRefBFinal"), helperEP.GetEventPlane(vec.qvecRe()[RefBId * 4 + 3], vec.qvecIm()[RefBId * 4 + 3], cfgnMod));
    }

    if (vec.qvecAmp()[DetId] > 1e-8 && vec.qvecAmp()[RefAId] > 1e-8 && vec.qvecAmp()[RefBId] > 1e-8) {
      histosQA.fill(HIST(qV::centClasses[cBin]) + HIST("histEvtPlRes_SigRefA"), helperEP.GetResolution(helperEP.GetEventPlane(vec.qvecRe()[DetId * 4 + 3], vec.qvecIm()[DetId * 4 + 3], cfgnMod), helperEP.GetEventPlane(vec.qvecRe()[RefAId * 4 + 3], vec.qvecIm()[RefAId * 4 + 3], cfgnMod), cfgnMod));
      histosQA.fill(HIST(qV::centClasses[cBin]) + HIST("histEvtPlRes_SigRefB"), helperEP.GetEventPlane(helperEP.GetEventPlane(vec.qvecRe()[DetId * 4 + 3], vec.qvecIm()[DetId * 4 + 3], cfgnMod), helperEP.GetEventPlane(vec.qvecRe()[RefBId * 4 + 3], vec.qvecIm()[RefBId * 4 + 3], cfgnMod), cfgnMod));
      histosQA.fill(HIST(qV::centClasses[cBin]) + HIST("histEvtPlRes_RefARefB"), helperEP.GetEventPlane(helperEP.GetEventPlane(vec.qvecRe()[RefAId * 4 + 3], vec.qvecIm()[RefAId * 4 + 3], cfgnMod), helperEP.GetEventPlane(vec.qvecRe()[RefBId * 4 + 3], vec.qvecIm()[RefBId * 4 + 3], cfgnMod), cfgnMod));
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
    } // End switch(centBin)
  }   // End void process(...)
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<qVectorsCorrection>(cfgc)};
}
