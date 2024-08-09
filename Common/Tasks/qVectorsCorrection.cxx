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
#include "Common/DataModel/EventSelection.h"
#include "Common/Core/EventPlaneHelper.h"

// o2 includes.

using namespace o2;
using namespace o2::framework;

using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Qvectors>;

struct qVectorsCorrection {
  // No correction = recenter, recentered Qvectors = twist, twisted Qvectors = rescale.
  // NOTE: As of no, the twist gets both twist and rescale correction constants.

  // Histogram registry for the output QA figures and list of centrality classes for it.
  // Objects are NOT saved in alphabetical orders, and registry names are NOT saved
  // as TDirectoryFile.
  HistogramRegistry histosQA{"histosQA", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  Configurable<std::vector<int>> cfgnMods{"cfgnMods", {2, 3}, "Modulation of interest"};

  Configurable<std::string> cfgDetName{"cfgDetName", "FT0C", "The name of detector to be analyzed"};
  Configurable<std::string> cfgRefAName{"cfgRefAName", "TPCpos", "The name of detector for reference A"};
  Configurable<std::string> cfgRefBName{"cfgRefBName", "TPCneg", "The name of detector for reference B"};
  Configurable<bool> cfgAddEvtSel{"cfgAddEvtSel", true, "event selection"};

  Configurable<int> cfgnTotalSystem{"cfgnTotalSystem", 7, "total qvector number"};

  ConfigurableAxis cfgaxisQvecF{"cfgaxisQvecF", {300, -1, 1}, ""};
  ConfigurableAxis cfgaxisQvec{"cfgaxisQvec", {100, -3, 3}, ""};
  ConfigurableAxis cfgaxisCent{"cfgaxisCent", {90, 0, 90}, ""};

  Configurable<bool> cfgQAAll{"cfgQAAll", false, "draw all q-vector steps"};
  Configurable<bool> cfgQAFinal{"cfgQAFinal", false, "draw final q-vector steps"};

  // Helper variables.
  EventPlaneHelper helperEP;

  int DetId;
  int RefAId;
  int RefBId;

  template <typename T>
  int GetDetId(const T& name)
  {
    if (name.value == "BPos" || name.value == "BNeg" || name.value == "BTot") {
      LOGF(warning, "Using deprecated label: %s. Please use TPCpos, TPCneg, TPCall instead.", name.value);
    }
    if (name.value == "FT0C") {
      return 0;
    } else if (name.value == "FT0A") {
      return 1;
    } else if (name.value == "FT0M") {
      return 2;
    } else if (name.value == "FV0A") {
      return 3;
    } else if (name.value == "TPCpos" || name.value == "BPos") {
      return 4;
    } else if (name.value == "TPCneg" || name.value == "BNeg") {
      return 5;
    } else if (name.value == "TPCall" || name.value == "BTot") {
      return 6;
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
      LOGF(info, "Wrong detector configuration \n The FT0C will be used to get Q-Vector \n The TPCpos and TPCneg will be used as reference systems");
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

    for (uint i = 0; i < cfgnMods->size(); i++) {
      histosQA.add(Form("histQvecUncorV%d", cfgnMods->at(i)), "", {HistType::kTH3F, {axisQvecF, axisQvecF, axisCent}});
      histosQA.add(Form("histQvecRefAUncorV%d", cfgnMods->at(i)), "", {HistType::kTH3F, {axisQvecF, axisQvecF, axisCent}});
      histosQA.add(Form("histQvecRefBUncorV%d", cfgnMods->at(i)), "", {HistType::kTH3F, {axisQvecF, axisQvecF, axisCent}});

      histosQA.add(Form("histEvtPlUncorV%d", cfgnMods->at(i)), "", {HistType::kTH2F, {axisEvtPl, axisCent}});
      histosQA.add(Form("histEvtPlRefAUncorV%d", cfgnMods->at(i)), "", {HistType::kTH2F, {axisEvtPl, axisCent}});
      histosQA.add(Form("histEvtPlRefBUncorV%d", cfgnMods->at(i)), "", {HistType::kTH2F, {axisEvtPl, axisCent}});

      if (cfgQAFinal) {
        histosQA.add(Form("histQvecFinalV%d", cfgnMods->at(i)), "", {HistType::kTH3F, {axisQvec, axisQvec, axisCent}});
        histosQA.add(Form("histQvecRefAFinalV%d", cfgnMods->at(i)), "", {HistType::kTH3F, {axisQvec, axisQvec, axisCent}});
        histosQA.add(Form("histQvecRefBFinalV%d", cfgnMods->at(i)), "", {HistType::kTH3F, {axisQvec, axisQvec, axisCent}});

        histosQA.add(Form("histEvtPlFinalV%d", cfgnMods->at(i)), "", {HistType::kTH2F, {axisEvtPl, axisCent}});
        histosQA.add(Form("histEvtPlRefAFinalV%d", cfgnMods->at(i)), "", {HistType::kTH2F, {axisEvtPl, axisCent}});
        histosQA.add(Form("histEvtPlRefBFinalV%d", cfgnMods->at(i)), "", {HistType::kTH2F, {axisEvtPl, axisCent}});

        histosQA.add(Form("histEvtPlRes_SigRefAV%d", cfgnMods->at(i)), "", {HistType::kTH2F, {axisEvtPl, axisCent}});
        histosQA.add(Form("histEvtPlRes_SigRefBV%d", cfgnMods->at(i)), "", {HistType::kTH2F, {axisEvtPl, axisCent}});
        histosQA.add(Form("histEvtPlRes_RefARefBV%d", cfgnMods->at(i)), "", {HistType::kTH2F, {axisEvtPl, axisCent}});

        if (cfgQAAll) {
          histosQA.add(Form("histQvecRectrV%d", cfgnMods->at(i)), "", {HistType::kTH3F, {axisQvecF, axisQvecF, axisCent}});
          histosQA.add(Form("histQvecTwistV%d", cfgnMods->at(i)), "", {HistType::kTH3F, {axisQvecF, axisQvecF, axisCent}});

          histosQA.add(Form("histQvecRefARectrV%d", cfgnMods->at(i)), "", {HistType::kTH3F, {axisQvecF, axisQvecF, axisCent}});
          histosQA.add(Form("histQvecRefATwistV%d", cfgnMods->at(i)), "", {HistType::kTH3F, {axisQvecF, axisQvecF, axisCent}});

          histosQA.add(Form("histQvecRefBRectrV%d", cfgnMods->at(i)), "", {HistType::kTH3F, {axisQvecF, axisQvecF, axisCent}});
          histosQA.add(Form("histQvecRefBTwistV%d", cfgnMods->at(i)), "", {HistType::kTH3F, {axisQvecF, axisQvecF, axisCent}});

          histosQA.add(Form("histEvtPlRectrV%d", cfgnMods->at(i)), "", {HistType::kTH2F, {axisEvtPl, axisCent}});
          histosQA.add(Form("histEvtPlTwistV%d", cfgnMods->at(i)), "", {HistType::kTH2F, {axisEvtPl, axisCent}});

          histosQA.add(Form("histEvtPlRefARectrV%d", cfgnMods->at(i)), "", {HistType::kTH2F, {axisEvtPl, axisCent}});
          histosQA.add(Form("histEvtPlRefATwistV%d", cfgnMods->at(i)), "", {HistType::kTH2F, {axisEvtPl, axisCent}});

          histosQA.add(Form("histEvtPlRefBRectrV%d", cfgnMods->at(i)), "", {HistType::kTH2F, {axisEvtPl, axisCent}});
          histosQA.add(Form("histEvtPlRefBTwistV%d", cfgnMods->at(i)), "", {HistType::kTH2F, {axisEvtPl, axisCent}});
        }
      }
    }
  } // End void init(InitContext const&)

  // Definition of all the needed template functions.
  template <typename T>
  void fillHistosQvec(const T& vec, int nmode)
  {
    int DetInd = DetId * 4 + cfgnTotalSystem * 4 * (nmode - 2);
    int RefAInd = RefAId * 4 + cfgnTotalSystem * 4 * (nmode - 2);
    int RefBInd = RefBId * 4 + cfgnTotalSystem * 4 * (nmode - 2);
    if (nmode == 2) {
      if (vec.qvecAmp()[DetId] > 1e-8) {
        histosQA.fill(HIST("histQvecUncorV2"), vec.qvecRe()[DetInd], vec.qvecIm()[DetInd], vec.cent());
        histosQA.fill(HIST("histEvtPlUncorV2"), helperEP.GetEventPlane(vec.qvecRe()[DetInd], vec.qvecIm()[DetInd], nmode), vec.cent());
        if (cfgQAFinal) {
          histosQA.fill(HIST("histQvecFinalV2"), vec.qvecRe()[DetInd + 3], vec.qvecIm()[DetInd + 3], vec.cent());
          histosQA.fill(HIST("histEvtPlFinalV2"), helperEP.GetEventPlane(vec.qvecRe()[DetInd + 3], vec.qvecIm()[DetInd + 3], nmode), vec.cent());
          if (cfgQAAll) {
            histosQA.fill(HIST("histQvecRectrV2"), vec.qvecRe()[DetInd + 1], vec.qvecIm()[DetInd + 1], vec.cent());
            histosQA.fill(HIST("histQvecTwistV2"), vec.qvecRe()[DetInd + 2], vec.qvecIm()[DetInd + 2], vec.cent());

            histosQA.fill(HIST("histEvtPlRectrV2"), helperEP.GetEventPlane(vec.qvecRe()[DetInd + 1], vec.qvecIm()[DetInd + 1], nmode), vec.cent());
            histosQA.fill(HIST("histEvtPlTwistV2"), helperEP.GetEventPlane(vec.qvecRe()[DetInd + 2], vec.qvecIm()[DetInd + 2], nmode), vec.cent());
          }
        }
      }
      if (vec.qvecAmp()[RefAId] > 1e-8) {
        histosQA.fill(HIST("histQvecRefAUncorV2"), vec.qvecRe()[RefAInd], vec.qvecIm()[RefAInd], vec.cent());
        histosQA.fill(HIST("histEvtPlRefAUncorV2"), helperEP.GetEventPlane(vec.qvecRe()[RefAInd], vec.qvecIm()[RefAInd], nmode), vec.cent());
        if (cfgQAFinal) {
          histosQA.fill(HIST("histQvecRefAFinalV2"), vec.qvecRe()[RefAInd + 3], vec.qvecIm()[RefAInd + 3], vec.cent());
          histosQA.fill(HIST("histEvtPlRefAFinalV2"), helperEP.GetEventPlane(vec.qvecRe()[RefAInd + 3], vec.qvecIm()[RefAInd + 3], nmode), vec.cent());
          if (cfgQAAll) {
            histosQA.fill(HIST("histQvecRefARectrV2"), vec.qvecRe()[RefAInd + 1], vec.qvecIm()[RefAInd + 1], vec.cent());
            histosQA.fill(HIST("histQvecRefATwistV2"), vec.qvecRe()[RefAInd + 2], vec.qvecIm()[RefAInd + 2], vec.cent());

            histosQA.fill(HIST("histEvtPlRefARectrV2"), helperEP.GetEventPlane(vec.qvecRe()[RefAInd + 1], vec.qvecIm()[RefAInd + 1], nmode), vec.cent());
            histosQA.fill(HIST("histEvtPlRefATwistV2"), helperEP.GetEventPlane(vec.qvecRe()[RefAInd + 2], vec.qvecIm()[RefAInd + 2], nmode), vec.cent());
          }
        }
      }
      if (vec.qvecAmp()[RefBId] > 1e-8) {
        histosQA.fill(HIST("histQvecRefBUncorV2"), vec.qvecRe()[RefBInd], vec.qvecIm()[RefBInd], vec.cent());
        histosQA.fill(HIST("histEvtPlRefBUncorV2"), helperEP.GetEventPlane(vec.qvecRe()[RefBInd], vec.qvecIm()[RefBInd], nmode), vec.cent());
        if (cfgQAFinal) {
          histosQA.fill(HIST("histQvecRefBFinalV2"), vec.qvecRe()[RefBInd + 3], vec.qvecIm()[RefBInd + 3], vec.cent());
          histosQA.fill(HIST("histEvtPlRefBFinalV2"), helperEP.GetEventPlane(vec.qvecRe()[RefBInd + 3], vec.qvecIm()[RefBInd + 3], nmode), vec.cent());
          if (cfgQAAll) {
            histosQA.fill(HIST("histQvecRefBRectrV2"), vec.qvecRe()[RefBInd + 1], vec.qvecIm()[RefBInd + 1], vec.cent());
            histosQA.fill(HIST("histQvecRefBTwistV2"), vec.qvecRe()[RefBInd + 2], vec.qvecIm()[RefBInd + 2], vec.cent());

            histosQA.fill(HIST("histEvtPlRefBRectrV2"), helperEP.GetEventPlane(vec.qvecRe()[RefBInd + 1], vec.qvecIm()[RefBInd + 1], nmode), vec.cent());
            histosQA.fill(HIST("histEvtPlRefBTwistV2"), helperEP.GetEventPlane(vec.qvecRe()[RefBInd + 2], vec.qvecIm()[RefBInd + 2], nmode), vec.cent());
          }
        }
      }
      if (vec.qvecAmp()[DetId] > 1e-8 && vec.qvecAmp()[RefAId] > 1e-8 && vec.qvecAmp()[RefBId] > 1e-8 && cfgQAFinal) {
        histosQA.fill(HIST("histEvtPlRes_SigRefAV2"), helperEP.GetResolution(helperEP.GetEventPlane(vec.qvecRe()[DetInd + 3], vec.qvecIm()[DetInd + 3], nmode), helperEP.GetEventPlane(vec.qvecRe()[RefAInd + 3], vec.qvecIm()[RefAInd + 3], nmode), nmode), vec.cent());
        histosQA.fill(HIST("histEvtPlRes_SigRefBV2"), helperEP.GetResolution(helperEP.GetEventPlane(vec.qvecRe()[DetInd + 3], vec.qvecIm()[DetInd + 3], nmode), helperEP.GetEventPlane(vec.qvecRe()[RefBInd + 3], vec.qvecIm()[RefBInd + 3], nmode), nmode), vec.cent());
        histosQA.fill(HIST("histEvtPlRes_RefARefBV2"), helperEP.GetResolution(helperEP.GetEventPlane(vec.qvecRe()[RefAInd + 3], vec.qvecIm()[RefAInd + 3], nmode), helperEP.GetEventPlane(vec.qvecRe()[RefBInd + 3], vec.qvecIm()[RefBInd + 3], nmode), nmode), vec.cent());
      }
    } else if (nmode == 3) {
      if (vec.qvecAmp()[DetId] > 1e-8) {
        histosQA.fill(HIST("histQvecUncorV3"), vec.qvecRe()[DetInd], vec.qvecIm()[DetInd], vec.cent());
        histosQA.fill(HIST("histEvtPlUncorV3"), helperEP.GetEventPlane(vec.qvecRe()[DetInd], vec.qvecIm()[DetInd], nmode), vec.cent());
        if (cfgQAFinal) {
          histosQA.fill(HIST("histQvecFinalV3"), vec.qvecRe()[DetInd + 3], vec.qvecIm()[DetInd + 3], vec.cent());
          histosQA.fill(HIST("histEvtPlFinalV3"), helperEP.GetEventPlane(vec.qvecRe()[DetInd + 3], vec.qvecIm()[DetInd + 3], nmode), vec.cent());
          if (cfgQAAll) {
            histosQA.fill(HIST("histQvecRectrV3"), vec.qvecRe()[DetInd + 1], vec.qvecIm()[DetInd + 1], vec.cent());
            histosQA.fill(HIST("histQvecTwistV3"), vec.qvecRe()[DetInd + 2], vec.qvecIm()[DetInd + 2], vec.cent());

            histosQA.fill(HIST("histEvtPlRectrV3"), helperEP.GetEventPlane(vec.qvecRe()[DetInd + 1], vec.qvecIm()[DetInd + 1], nmode), vec.cent());
            histosQA.fill(HIST("histEvtPlTwistV3"), helperEP.GetEventPlane(vec.qvecRe()[DetInd + 2], vec.qvecIm()[DetInd + 2], nmode), vec.cent());
          }
        }
      }
      if (vec.qvecAmp()[RefAId] > 1e-8) {
        histosQA.fill(HIST("histQvecRefAUncorV3"), vec.qvecRe()[RefAInd], vec.qvecIm()[RefAInd], vec.cent());
        histosQA.fill(HIST("histEvtPlRefAUncorV3"), helperEP.GetEventPlane(vec.qvecRe()[RefAInd], vec.qvecIm()[RefAInd], nmode), vec.cent());
        if (cfgQAFinal) {
          histosQA.fill(HIST("histQvecRefAFinalV3"), vec.qvecRe()[RefAInd + 3], vec.qvecIm()[RefAInd + 3], vec.cent());
          histosQA.fill(HIST("histEvtPlRefAFinalV3"), helperEP.GetEventPlane(vec.qvecRe()[RefAInd + 3], vec.qvecIm()[RefAInd + 3], nmode), vec.cent());
          if (cfgQAAll) {
            histosQA.fill(HIST("histQvecRefARectrV3"), vec.qvecRe()[RefAInd + 1], vec.qvecIm()[RefAInd + 1], vec.cent());
            histosQA.fill(HIST("histQvecRefATwistV3"), vec.qvecRe()[RefAInd + 2], vec.qvecIm()[RefAInd + 2], vec.cent());

            histosQA.fill(HIST("histEvtPlRefARectrV3"), helperEP.GetEventPlane(vec.qvecRe()[RefAInd + 1], vec.qvecIm()[RefAInd + 1], nmode), vec.cent());
            histosQA.fill(HIST("histEvtPlRefATwistV3"), helperEP.GetEventPlane(vec.qvecRe()[RefAInd + 2], vec.qvecIm()[RefAInd + 2], nmode), vec.cent());
          }
        }
      }
      if (vec.qvecAmp()[RefBId] > 1e-8) {
        histosQA.fill(HIST("histQvecRefBUncorV3"), vec.qvecRe()[RefBInd], vec.qvecIm()[RefBInd], vec.cent());
        histosQA.fill(HIST("histEvtPlRefBUncorV3"), helperEP.GetEventPlane(vec.qvecRe()[RefBInd], vec.qvecIm()[RefBInd], nmode), vec.cent());
        if (cfgQAFinal) {
          histosQA.fill(HIST("histQvecRefBFinalV3"), vec.qvecRe()[RefBInd + 3], vec.qvecIm()[RefBInd + 3], vec.cent());
          histosQA.fill(HIST("histEvtPlRefBFinalV3"), helperEP.GetEventPlane(vec.qvecRe()[RefBInd + 3], vec.qvecIm()[RefBInd + 3], nmode), vec.cent());
          if (cfgQAAll) {
            histosQA.fill(HIST("histQvecRefBRectrV3"), vec.qvecRe()[RefBInd + 1], vec.qvecIm()[RefBInd + 1], vec.cent());
            histosQA.fill(HIST("histQvecRefBTwistV3"), vec.qvecRe()[RefBInd + 2], vec.qvecIm()[RefBInd + 2], vec.cent());

            histosQA.fill(HIST("histEvtPlRefBRectrV3"), helperEP.GetEventPlane(vec.qvecRe()[RefBInd + 1], vec.qvecIm()[RefBInd + 1], nmode), vec.cent());
            histosQA.fill(HIST("histEvtPlRefBTwistV3"), helperEP.GetEventPlane(vec.qvecRe()[RefBInd + 2], vec.qvecIm()[RefBInd + 2], nmode), vec.cent());
          }
        }
      }
      if (vec.qvecAmp()[DetId] > 1e-8 && vec.qvecAmp()[RefAId] > 1e-8 && vec.qvecAmp()[RefBId] > 1e-8 && cfgQAFinal) {
        histosQA.fill(HIST("histEvtPlRes_SigRefAV3"), helperEP.GetResolution(helperEP.GetEventPlane(vec.qvecRe()[DetInd + 3], vec.qvecIm()[DetInd + 3], nmode), helperEP.GetEventPlane(vec.qvecRe()[RefAInd + 3], vec.qvecIm()[RefAInd + 3], nmode), nmode), vec.cent());
        histosQA.fill(HIST("histEvtPlRes_SigRefBV3"), helperEP.GetResolution(helperEP.GetEventPlane(vec.qvecRe()[DetInd + 3], vec.qvecIm()[DetInd + 3], nmode), helperEP.GetEventPlane(vec.qvecRe()[RefBInd + 3], vec.qvecIm()[RefBInd + 3], nmode), nmode), vec.cent());
        histosQA.fill(HIST("histEvtPlRes_RefARefBV3"), helperEP.GetResolution(helperEP.GetEventPlane(vec.qvecRe()[RefAInd + 3], vec.qvecIm()[RefAInd + 3], nmode), helperEP.GetEventPlane(vec.qvecRe()[RefBInd + 3], vec.qvecIm()[RefBInd + 3], nmode), nmode), vec.cent());
      }
    } else if (nmode == 4) {
      if (vec.qvecAmp()[DetId] > 1e-8) {
        histosQA.fill(HIST("histQvecUncorV4"), vec.qvecRe()[DetInd], vec.qvecIm()[DetInd], vec.cent());
        histosQA.fill(HIST("histEvtPlUncorV4"), helperEP.GetEventPlane(vec.qvecRe()[DetInd], vec.qvecIm()[DetInd], nmode), vec.cent());
        if (cfgQAFinal) {
          histosQA.fill(HIST("histQvecFinalV4"), vec.qvecRe()[DetInd + 3], vec.qvecIm()[DetInd + 3], vec.cent());
          histosQA.fill(HIST("histEvtPlFinalV4"), helperEP.GetEventPlane(vec.qvecRe()[DetInd + 3], vec.qvecIm()[DetInd + 3], nmode), vec.cent());
          if (cfgQAAll) {
            histosQA.fill(HIST("histQvecRectrV4"), vec.qvecRe()[DetInd + 1], vec.qvecIm()[DetInd + 1], vec.cent());
            histosQA.fill(HIST("histQvecTwistV4"), vec.qvecRe()[DetInd + 2], vec.qvecIm()[DetInd + 2], vec.cent());

            histosQA.fill(HIST("histEvtPlRectrV4"), helperEP.GetEventPlane(vec.qvecRe()[DetInd + 1], vec.qvecIm()[DetInd + 1], nmode), vec.cent());
            histosQA.fill(HIST("histEvtPlTwistV4"), helperEP.GetEventPlane(vec.qvecRe()[DetInd + 2], vec.qvecIm()[DetInd + 2], nmode), vec.cent());
          }
        }
      }
      if (vec.qvecAmp()[RefAId] > 1e-8) {
        histosQA.fill(HIST("histQvecRefAUncorV4"), vec.qvecRe()[RefAInd], vec.qvecIm()[RefAInd], vec.cent());
        histosQA.fill(HIST("histEvtPlRefAUncorV4"), helperEP.GetEventPlane(vec.qvecRe()[RefAInd], vec.qvecIm()[RefAInd], nmode), vec.cent());
        if (cfgQAFinal) {
          histosQA.fill(HIST("histQvecRefAFinalV4"), vec.qvecRe()[RefAInd + 3], vec.qvecIm()[RefAInd + 3], vec.cent());
          histosQA.fill(HIST("histEvtPlRefAFinalV4"), helperEP.GetEventPlane(vec.qvecRe()[RefAInd + 3], vec.qvecIm()[RefAInd + 3], nmode), vec.cent());
          if (cfgQAAll) {
            histosQA.fill(HIST("histQvecRefARectrV4"), vec.qvecRe()[RefAInd + 1], vec.qvecIm()[RefAInd + 1], vec.cent());
            histosQA.fill(HIST("histQvecRefATwistV4"), vec.qvecRe()[RefAInd + 2], vec.qvecIm()[RefAInd + 2], vec.cent());

            histosQA.fill(HIST("histEvtPlRefARectrV4"), helperEP.GetEventPlane(vec.qvecRe()[RefAInd + 1], vec.qvecIm()[RefAInd + 1], nmode), vec.cent());
            histosQA.fill(HIST("histEvtPlRefATwistV4"), helperEP.GetEventPlane(vec.qvecRe()[RefAInd + 2], vec.qvecIm()[RefAInd + 2], nmode), vec.cent());
          }
        }
      }
      if (vec.qvecAmp()[RefBId] > 1e-8) {
        histosQA.fill(HIST("histQvecRefBUncorV4"), vec.qvecRe()[RefBInd], vec.qvecIm()[RefBInd], vec.cent());
        histosQA.fill(HIST("histEvtPlRefBUncorV4"), helperEP.GetEventPlane(vec.qvecRe()[RefBInd], vec.qvecIm()[RefBInd], nmode), vec.cent());
        if (cfgQAFinal) {
          histosQA.fill(HIST("histQvecRefBFinalV4"), vec.qvecRe()[RefBInd + 3], vec.qvecIm()[RefBInd + 3], vec.cent());
          histosQA.fill(HIST("histEvtPlRefBFinalV4"), helperEP.GetEventPlane(vec.qvecRe()[RefBInd + 3], vec.qvecIm()[RefBInd + 3], nmode), vec.cent());
          if (cfgQAAll) {
            histosQA.fill(HIST("histQvecRefBRectrV4"), vec.qvecRe()[RefBInd + 1], vec.qvecIm()[RefBInd + 1], vec.cent());
            histosQA.fill(HIST("histQvecRefBTwistV4"), vec.qvecRe()[RefBInd + 2], vec.qvecIm()[RefBInd + 2], vec.cent());

            histosQA.fill(HIST("histEvtPlRefBRectrV4"), helperEP.GetEventPlane(vec.qvecRe()[RefBInd + 1], vec.qvecIm()[RefBInd + 1], nmode), vec.cent());
            histosQA.fill(HIST("histEvtPlRefBTwistV4"), helperEP.GetEventPlane(vec.qvecRe()[RefBInd + 2], vec.qvecIm()[RefBInd + 2], nmode), vec.cent());
          }
        }
      }
      if (vec.qvecAmp()[DetId] > 1e-8 && vec.qvecAmp()[RefAId] > 1e-8 && vec.qvecAmp()[RefBId] > 1e-8 && cfgQAFinal) {
        histosQA.fill(HIST("histEvtPlRes_SigRefAV4"), helperEP.GetResolution(helperEP.GetEventPlane(vec.qvecRe()[DetInd + 3], vec.qvecIm()[DetInd + 3], nmode), helperEP.GetEventPlane(vec.qvecRe()[RefAInd + 3], vec.qvecIm()[RefAInd + 3], nmode), nmode), vec.cent());
        histosQA.fill(HIST("histEvtPlRes_SigRefBV4"), helperEP.GetResolution(helperEP.GetEventPlane(vec.qvecRe()[DetInd + 3], vec.qvecIm()[DetInd + 3], nmode), helperEP.GetEventPlane(vec.qvecRe()[RefBInd + 3], vec.qvecIm()[RefBInd + 3], nmode), nmode), vec.cent());
        histosQA.fill(HIST("histEvtPlRes_RefARefBV4"), helperEP.GetResolution(helperEP.GetEventPlane(vec.qvecRe()[RefAInd + 3], vec.qvecIm()[RefAInd + 3], nmode), helperEP.GetEventPlane(vec.qvecRe()[RefBInd + 3], vec.qvecIm()[RefBInd + 3], nmode), nmode), vec.cent());
      }
    }
  }
  void process(MyCollisions::iterator const& qVec)
  {
    histosQA.fill(HIST("histCentFull"), qVec.cent());
    if (cfgAddEvtSel && (!qVec.sel8() ||
                         !qVec.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV) ||
                         !qVec.selection_bit(aod::evsel::kNoSameBunchPileup))) {
      return;
    }
    for (uint i = 0; i < cfgnMods->size(); i++) {
      fillHistosQvec(qVec, cfgnMods->at(i));
    }
  } // End void process(...)
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<qVectorsCorrection>(cfgc)};
}
