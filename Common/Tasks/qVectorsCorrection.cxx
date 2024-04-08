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

  ConfigurableAxis cfgaxisQvecF{"cfgaxisQvecF", {100, -1, 1}, ""};
  ConfigurableAxis cfgaxisQvec{"cfgaxisQvec", {100, -3, 3}, ""};
  ConfigurableAxis cfgaxisCent{"cfgaxisCent", {90, 0, 90}, ""};

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
    AxisSpec axisCent{cfgaxisCent, "centrality"};
    AxisSpec axisQvec{cfgaxisQvec, "Q"};
    AxisSpec axisQvecF{cfgaxisQvecF, "Q"};
    AxisSpec axisEvtPl{360, -constants::math::PI, constants::math::PI};

    histosQA.add("histCentFull", "Centrality distribution for valid events",
                 HistType::kTH1F, {axisCent});

    histosQA.add("histQvecUncor", "", {HistType::kTH3F, {axisQvecF, axisQvecF, axisCent}});
    histosQA.add("histQvecRectr", "", {HistType::kTH3F, {axisQvecF, axisQvecF, axisCent}});
    histosQA.add("histQvecTwist", "", {HistType::kTH3F, {axisQvecF, axisQvecF, axisCent}});
    histosQA.add("histQvecFinal", "", {HistType::kTH3F, {axisQvec, axisQvec, axisCent}});

    histosQA.add("histQvecRefAUncor", "", {HistType::kTH3F, {axisQvecF, axisQvecF, axisCent}});
    histosQA.add("histQvecRefARectr", "", {HistType::kTH3F, {axisQvecF, axisQvecF, axisCent}});
    histosQA.add("histQvecRefATwist", "", {HistType::kTH3F, {axisQvecF, axisQvecF, axisCent}});
    histosQA.add("histQvecRefAFinal", "", {HistType::kTH3F, {axisQvec, axisQvec, axisCent}});

    histosQA.add("histQvecRefBUncor", "", {HistType::kTH3F, {axisQvecF, axisQvecF, axisCent}});
    histosQA.add("histQvecRefBRectr", "", {HistType::kTH3F, {axisQvecF, axisQvecF, axisCent}});
    histosQA.add("histQvecRefBTwist", "", {HistType::kTH3F, {axisQvecF, axisQvecF, axisCent}});
    histosQA.add("histQvecRefBFinal", "", {HistType::kTH3F, {axisQvec, axisQvec, axisCent}});

    histosQA.add("histEvtPlUncor", "", {HistType::kTH2F, {axisEvtPl, axisCent}});
    histosQA.add("histEvtPlRectr", "", {HistType::kTH2F, {axisEvtPl, axisCent}});
    histosQA.add("histEvtPlTwist", "", {HistType::kTH2F, {axisEvtPl, axisCent}});
    histosQA.add("histEvtPlFinal", "", {HistType::kTH2F, {axisEvtPl, axisCent}});

    histosQA.add("histEvtPlRefAUncor", "", {HistType::kTH2F, {axisEvtPl, axisCent}});
    histosQA.add("histEvtPlRefARectr", "", {HistType::kTH2F, {axisEvtPl, axisCent}});
    histosQA.add("histEvtPlRefATwist", "", {HistType::kTH2F, {axisEvtPl, axisCent}});
    histosQA.add("histEvtPlRefAFinal", "", {HistType::kTH2F, {axisEvtPl, axisCent}});

    histosQA.add("histEvtPlRefBUncor", "", {HistType::kTH2F, {axisEvtPl, axisCent}});
    histosQA.add("histEvtPlRefBRectr", "", {HistType::kTH2F, {axisEvtPl, axisCent}});
    histosQA.add("histEvtPlRefBTwist", "", {HistType::kTH2F, {axisEvtPl, axisCent}});
    histosQA.add("histEvtPlRefBFinal", "", {HistType::kTH2F, {axisEvtPl, axisCent}});

    histosQA.add("histEvtPlRes_SigRefA", "", {HistType::kTH2F, {axisEvtPl, axisCent}});
    histosQA.add("histEvtPlRes_SigRefB", "", {HistType::kTH2F, {axisEvtPl, axisCent}});
    histosQA.add("histEvtPlRes_RefARefB", "", {HistType::kTH2F, {axisEvtPl, axisCent}});

  } // End void init(InitContext const&)

  // Definition of all the needed template functions.
  template <typename T>
  void fillHistosQvec(const T& vec)
  {
    if (vec.qvecAmp()[DetId] > 1e-8) {
      histosQA.fill(HIST("histQvecUncor"), vec.qvecRe()[DetId * 4], vec.qvecIm()[DetId * 4], vec.cent());
      histosQA.fill(HIST("histQvecRectr"), vec.qvecRe()[DetId * 4 + 1], vec.qvecIm()[DetId * 4 + 1], vec.cent());
      histosQA.fill(HIST("histQvecTwist"), vec.qvecRe()[DetId * 4 + 2], vec.qvecIm()[DetId * 4 + 2], vec.cent());
      histosQA.fill(HIST("histQvecFinal"), vec.qvecRe()[DetId * 4 + 3], vec.qvecIm()[DetId * 4 + 3], vec.cent());

      histosQA.fill(HIST("histEvtPlUncor"), helperEP.GetEventPlane(vec.qvecRe()[DetId * 4], vec.qvecIm()[DetId * 4], cfgnMod), vec.cent());
      histosQA.fill(HIST("histEvtPlRectr"), helperEP.GetEventPlane(vec.qvecRe()[DetId * 4 + 1], vec.qvecIm()[DetId * 4 + 1], cfgnMod), vec.cent());
      histosQA.fill(HIST("histEvtPlTwist"), helperEP.GetEventPlane(vec.qvecRe()[DetId * 4 + 2], vec.qvecIm()[DetId * 4 + 2], cfgnMod), vec.cent());
      histosQA.fill(HIST("histEvtPlFinal"), helperEP.GetEventPlane(vec.qvecRe()[DetId * 4 + 3], vec.qvecIm()[DetId * 4 + 3], cfgnMod), vec.cent());
    }

    if (vec.qvecAmp()[RefAId] > 1e-8) {
      histosQA.fill(HIST("histQvecRefAUncor"), vec.qvecRe()[RefAId * 4], vec.qvecIm()[RefAId * 4], vec.cent());
      histosQA.fill(HIST("histQvecRefARectr"), vec.qvecRe()[RefAId * 4 + 1], vec.qvecIm()[RefAId * 4 + 1], vec.cent());
      histosQA.fill(HIST("histQvecRefATwist"), vec.qvecRe()[RefAId * 4 + 2], vec.qvecIm()[RefAId * 4 + 2], vec.cent());
      histosQA.fill(HIST("histQvecRefAFinal"), vec.qvecRe()[RefAId * 4 + 3], vec.qvecIm()[RefAId * 4 + 3], vec.cent());

      histosQA.fill(HIST("histEvtPlRefAUncor"), helperEP.GetEventPlane(vec.qvecRe()[RefAId * 4], vec.qvecIm()[RefAId * 4], cfgnMod), vec.cent());
      histosQA.fill(HIST("histEvtPlRefARectr"), helperEP.GetEventPlane(vec.qvecRe()[RefAId * 4 + 1], vec.qvecIm()[RefAId * 4 + 1], cfgnMod), vec.cent());
      histosQA.fill(HIST("histEvtPlRefATwist"), helperEP.GetEventPlane(vec.qvecRe()[RefAId * 4 + 2], vec.qvecIm()[RefAId * 4 + 2], cfgnMod), vec.cent());
      histosQA.fill(HIST("histEvtPlRefAFinal"), helperEP.GetEventPlane(vec.qvecRe()[RefAId * 4 + 3], vec.qvecIm()[RefAId * 4 + 3], cfgnMod), vec.cent());
    }

    if (vec.qvecAmp()[RefBId] > 1e-8) {
      histosQA.fill(HIST("histQvecRefBUncor"), vec.qvecRe()[RefBId * 4], vec.qvecIm()[RefBId * 4], vec.cent());
      histosQA.fill(HIST("histQvecRefBRectr"), vec.qvecRe()[RefBId * 4 + 1], vec.qvecIm()[RefBId * 4 + 1], vec.cent());
      histosQA.fill(HIST("histQvecRefBTwist"), vec.qvecRe()[RefBId * 4 + 2], vec.qvecIm()[RefBId * 4 + 2], vec.cent());
      histosQA.fill(HIST("histQvecRefBFinal"), vec.qvecRe()[RefBId * 4 + 3], vec.qvecIm()[RefBId * 4 + 3], vec.cent());

      histosQA.fill(HIST("histEvtPlRefBUncor"), helperEP.GetEventPlane(vec.qvecRe()[RefBId * 4], vec.qvecIm()[RefBId * 4], cfgnMod), vec.cent());
      histosQA.fill(HIST("histEvtPlRefBRectr"), helperEP.GetEventPlane(vec.qvecRe()[RefBId * 4 + 1], vec.qvecIm()[RefBId * 4 + 1], cfgnMod), vec.cent());
      histosQA.fill(HIST("histEvtPlRefBTwist"), helperEP.GetEventPlane(vec.qvecRe()[RefBId * 4 + 2], vec.qvecIm()[RefBId * 4 + 2], cfgnMod), vec.cent());
      histosQA.fill(HIST("histEvtPlRefBFinal"), helperEP.GetEventPlane(vec.qvecRe()[RefBId * 4 + 3], vec.qvecIm()[RefBId * 4 + 3], cfgnMod), vec.cent());
    }

    if (vec.qvecAmp()[DetId] > 1e-8 && vec.qvecAmp()[RefAId] > 1e-8 && vec.qvecAmp()[RefBId] > 1e-8) {
      histosQA.fill(HIST("histEvtPlRes_SigRefA"), helperEP.GetResolution(helperEP.GetEventPlane(vec.qvecRe()[DetId * 4 + 3], vec.qvecIm()[DetId * 4 + 3], cfgnMod), helperEP.GetEventPlane(vec.qvecRe()[RefAId * 4 + 3], vec.qvecIm()[RefAId * 4 + 3], cfgnMod), cfgnMod), vec.cent());
      histosQA.fill(HIST("histEvtPlRes_SigRefB"), helperEP.GetResolution(helperEP.GetEventPlane(vec.qvecRe()[DetId * 4 + 3], vec.qvecIm()[DetId * 4 + 3], cfgnMod), helperEP.GetEventPlane(vec.qvecRe()[RefBId * 4 + 3], vec.qvecIm()[RefBId * 4 + 3], cfgnMod), cfgnMod), vec.cent());
      histosQA.fill(HIST("histEvtPlRes_RefARefB"), helperEP.GetResolution(helperEP.GetEventPlane(vec.qvecRe()[RefAId * 4 + 3], vec.qvecIm()[RefAId * 4 + 3], cfgnMod), helperEP.GetEventPlane(vec.qvecRe()[RefBId * 4 + 3], vec.qvecIm()[RefBId * 4 + 3], cfgnMod), cfgnMod), vec.cent());
    }
  }

  void process(aod::Qvector const& qVec)
  {
    histosQA.fill(HIST("histCentFull"), qVec.cent());
    fillHistosQvec(qVec);
  } // End void process(...)
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<qVectorsCorrection>(cfgc)};
}
