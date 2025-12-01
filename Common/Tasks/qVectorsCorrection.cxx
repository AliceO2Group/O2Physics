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
#include <TString.h>
#include <TVector2.h>

#include <sys/types.h>

#include <string>
#include <vector>

// o2Physics includes.
#include "Common/CCDB/EventSelectionParams.h"
#include "Common/Core/EventPlaneHelper.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Qvectors.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CommonConstants/MathConstants.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

// o2 includes.

using namespace o2;
using namespace o2::framework;

using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Qvectors>;
using MyCollisionsWithSC = soa::Join<aod::Collisions, aod::EvSels, aod::QvectorsShifteds>;
using MyTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TrackSelectionExtension>;

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

  Configurable<int> cfgEvtSel{"cfgEvtSel", 0, "Event selection flags\n0: Sel8\n1: Sel8+kIsGoodZvtxFT0vsPV+kNoSameBunchPileup\n2: Sel8+kIsGoodZvtxFT0vsPV+kNoSameBunchPileup+kNoCollInTimeRangeStandard\n3: Sel8+kNoSameBunchPileup"};

  Configurable<int> cfgnTotalSystem{"cfgnTotalSystem", 7, "total qvector number"};
  Configurable<int> cfgNbinsEP{"cfgNbinsEP", 360, "nbins for EP histograms"};

  Configurable<bool> cfgQAAll{"cfgQAAll", false, "draw all q-vector steps"};
  Configurable<bool> cfgQAFinal{"cfgQAFinal", false, "draw final q-vector steps"};
  Configurable<bool> cfgQAFlowStudy{"cfgQAFlowStudy", false, "configurable for flow study"};
  Configurable<bool> cfgQAOccupancyStudy{"cfgQAOccupancyStudy", false, "configurable for occupancy study"};
  Configurable<bool> cfgAddEvtSelPileup{"cfgAddEvtSelPileup", false, "configurable for pileup selection"};
  Configurable<bool> cfgShiftCorPrep{"cfgShiftCorPrep", false, "configurable for shift correction"};

  Configurable<float> cfgMinPt{"cfgMinPt", 0.15, "Minimum transverse momentum for charged track"};
  Configurable<float> cfgMaxEta{"cfgMaxEta", 0.8, "Maximum pseudorapidiy for charged track"};
  Configurable<float> cfgMaxDCArToPVcut{"cfgMaxDCArToPVcut", 0.1, "Maximum transverse DCA"};
  Configurable<float> cfgMaxDCAzToPVcut{"cfgMaxDCAzToPVcut", 1.0, "Maximum longitudinal DCA"};

  Configurable<int> cfgMaxOccupancy{"cfgMaxOccupancy", 999999, "maximum occupancy of tracks in neighbouring collisions in a given time range"};
  Configurable<int> cfgMinOccupancy{"cfgMinOccupancy", 0, "maximum occupancy of tracks in neighbouring collisions in a given time range"};

  ConfigurableAxis cfgaxisQvecF{"cfgaxisQvecF", {300, -1, 1}, ""};
  ConfigurableAxis cfgaxisQvec{"cfgaxisQvec", {100, -3, 3}, ""};
  ConfigurableAxis cfgaxisCent{"cfgaxisCent", {100, 0, 100}, ""};

  ConfigurableAxis cfgaxiscos{"cfgaxiscos", {102, -1.02, 1.02}, ""};
  ConfigurableAxis cfgaxispt{"cfgaxispt", {100, 0, 10}, ""};
  ConfigurableAxis cfgaxisCentMerged{"cfgaxisCentMerged", {20, 0, 100}, ""};
  ConfigurableAxis cfgaxisAzimuth{"cfgaxisAzimuth", {72, 0, 2.0 * constants::math::PI}, ""};
  ConfigurableAxis cfgaxisOccupancy{"cfgaxisOccupancy", {VARIABLE_WIDTH, -1, 0, 100, 500, 1000, 2000, 3000, 4000, 5000, 10000, 99999}, ""};

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

  template <typename TrackType>
  bool SelTrack(const TrackType track)
  {
    if (track.pt() < 0.15)
      return false;
    if (!track.passedITSNCls())
      return false;
    if (!track.passedITSChi2NDF())
      return false;
    if (!track.passedITSHits())
      return false;
    if (!track.passedTPCCrossedRowsOverNCls())
      return false;
    if (!track.passedTPCChi2NDF())
      return false;
    if (!track.passedDCAxy())
      return false;
    if (!track.passedDCAz())
      return false;

    return true;
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
    AxisSpec axisEvtPl{cfgNbinsEP, -constants::math::PI, constants::math::PI};

    AxisSpec axisCos{cfgaxiscos, "angle function"};
    AxisSpec axisPt{cfgaxispt, "trasverse momentum"};
    AxisSpec axisCentMerged{cfgaxisCentMerged, "merged centrality"};
    AxisSpec axisAzimuth{cfgaxisAzimuth, "relative azimuthal angle"};
    AxisSpec axisOccupancy{cfgaxisOccupancy, "Occupancy"};

    AxisSpec axisShift = {10, 0, 10, "shift"};
    AxisSpec axisBasis = {20, 0, 20, "basis"};
    AxisSpec axisVertex = {220, -11, 11, "vertex"};

    histosQA.add("histCentFull", "Centrality distribution for valid events", HistType::kTH1F, {axisCent});
    histosQA.add("histCentSelected", "Centrality distribution for valid events", HistType::kTH1F, {axisCent});
    histosQA.add("histVtxSelected", "Centrality distribution for valid events", HistType::kTH1F, {axisVertex});

    for (uint i = 0; i < cfgnMods->size(); i++) {
      histosQA.add(Form("histQvecUncorV%d", cfgnMods->at(i)), "", {HistType::kTH3F, {axisQvecF, axisQvecF, axisCent}});
      histosQA.add(Form("histQvecRefAUncorV%d", cfgnMods->at(i)), "", {HistType::kTH3F, {axisQvecF, axisQvecF, axisCent}});
      histosQA.add(Form("histQvecRefBUncorV%d", cfgnMods->at(i)), "", {HistType::kTH3F, {axisQvecF, axisQvecF, axisCent}});

      histosQA.add(Form("histEvtPlUncorV%d", cfgnMods->at(i)), "", {HistType::kTH2F, {axisEvtPl, axisCent}});
      histosQA.add(Form("histEvtPlRefAUncorV%d", cfgnMods->at(i)), "", {HistType::kTH2F, {axisEvtPl, axisCent}});
      histosQA.add(Form("histEvtPlRefBUncorV%d", cfgnMods->at(i)), "", {HistType::kTH2F, {axisEvtPl, axisCent}});

      if (cfgQAOccupancyStudy) {
        histosQA.add(Form("histQvecOccUncorV%d", cfgnMods->at(i)), "", {HistType::kTHnSparseF, {axisQvecF, axisQvecF, axisCent, axisOccupancy}});
        histosQA.add(Form("histQvecRefAOccUncorV%d", cfgnMods->at(i)), "", {HistType::kTHnSparseF, {axisQvecF, axisQvecF, axisCent, axisOccupancy}});
        histosQA.add(Form("histQvecRefBOccUncorV%d", cfgnMods->at(i)), "", {HistType::kTHnSparseF, {axisQvecF, axisQvecF, axisCent, axisOccupancy}});
      }

      if (cfgQAFinal) {
        histosQA.add(Form("histQvecFinalV%d", cfgnMods->at(i)), "", {HistType::kTH3F, {axisQvec, axisQvec, axisCent}});
        histosQA.add(Form("histQvecRefAFinalV%d", cfgnMods->at(i)), "", {HistType::kTH3F, {axisQvec, axisQvec, axisCent}});
        histosQA.add(Form("histQvecRefBFinalV%d", cfgnMods->at(i)), "", {HistType::kTH3F, {axisQvec, axisQvec, axisCent}});

        if (cfgShiftCorPrep) {
          histosQA.add(Form("histShiftV%d", cfgnMods->at(i)), "", {HistType::kTProfile3D, {axisCent, axisBasis, axisShift}});
        }

        if (cfgQAOccupancyStudy) {
          histosQA.add(Form("histQvecOccFinalV%d", cfgnMods->at(i)), "", {HistType::kTHnSparseF, {axisQvecF, axisQvecF, axisCent, axisOccupancy}});
          histosQA.add(Form("histQvecRefAOccFinalV%d", cfgnMods->at(i)), "", {HistType::kTHnSparseF, {axisQvecF, axisQvecF, axisCent, axisOccupancy}});
          histosQA.add(Form("histQvecRefBOccFinalV%d", cfgnMods->at(i)), "", {HistType::kTHnSparseF, {axisQvecF, axisQvecF, axisCent, axisOccupancy}});
        }

        histosQA.add(Form("histQvecRes_SigRefAV%d", cfgnMods->at(i)), "", {HistType::kTH2F, {axisQvecF, axisCent}});
        histosQA.add(Form("histQvecRes_SigRefBV%d", cfgnMods->at(i)), "", {HistType::kTH2F, {axisQvecF, axisCent}});
        histosQA.add(Form("histQvecRes_RefARefBV%d", cfgnMods->at(i)), "", {HistType::kTH2F, {axisQvecF, axisCent}});

        histosQA.add(Form("histEvtPlFinalV%d", cfgnMods->at(i)), "", {HistType::kTH2F, {axisEvtPl, axisCent}});
        histosQA.add(Form("histEvtPlRefAFinalV%d", cfgnMods->at(i)), "", {HistType::kTH2F, {axisEvtPl, axisCent}});
        histosQA.add(Form("histEvtPlRefBFinalV%d", cfgnMods->at(i)), "", {HistType::kTH2F, {axisEvtPl, axisCent}});

        histosQA.add(Form("histEvtPlRes_SigRefAV%d", cfgnMods->at(i)), "", {HistType::kTH2F, {axisEvtPl, axisCent}});
        histosQA.add(Form("histEvtPlRes_SigRefBV%d", cfgnMods->at(i)), "", {HistType::kTH2F, {axisEvtPl, axisCent}});
        histosQA.add(Form("histEvtPlRes_RefARefBV%d", cfgnMods->at(i)), "", {HistType::kTH2F, {axisEvtPl, axisCent}});

        histosQA.add(Form("hist_EP_cos_Det_v%d", cfgnMods->at(i)), "", {HistType::kTH3F, {axisCentMerged, axisPt, axisCos}});
        histosQA.add(Form("hist_EP_sin_Det_v%d", cfgnMods->at(i)), "", {HistType::kTH3F, {axisCentMerged, axisPt, axisCos}});
        histosQA.add(Form("hist_EP_azimuth_Det_v%d", cfgnMods->at(i)), "", {HistType::kTH3F, {axisCentMerged, axisPt, axisAzimuth}});

        histosQA.add(Form("hist_EP_cos_RefA_v%d", cfgnMods->at(i)), "", {HistType::kTH3F, {axisCentMerged, axisPt, axisCos}});
        histosQA.add(Form("hist_EP_sin_RefA_v%d", cfgnMods->at(i)), "", {HistType::kTH3F, {axisCentMerged, axisPt, axisCos}});
        histosQA.add(Form("hist_EP_azimuth_RefA_v%d", cfgnMods->at(i)), "", {HistType::kTH3F, {axisCentMerged, axisPt, axisAzimuth}});

        histosQA.add(Form("hist_EP_cos_RefB_v%d", cfgnMods->at(i)), "", {HistType::kTH3F, {axisCentMerged, axisPt, axisCos}});
        histosQA.add(Form("hist_EP_sin_RefB_v%d", cfgnMods->at(i)), "", {HistType::kTH3F, {axisCentMerged, axisPt, axisCos}});
        histosQA.add(Form("hist_EP_azimuth_RefB_v%d", cfgnMods->at(i)), "", {HistType::kTH3F, {axisCentMerged, axisPt, axisAzimuth}});

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

  template <typename T>
  void fillHistosShiftCor(const T& vec, int nmode)
  {
    int DetInd = DetId * 4 + cfgnTotalSystem * 4 * (nmode - 2) + 3;
    int RefAInd = RefAId * 4 + cfgnTotalSystem * 4 * (nmode - 2) + 3;
    int RefBInd = RefBId * 4 + cfgnTotalSystem * 4 * (nmode - 2) + 3;

    if (vec.qvecAmp()[DetId] < 1e-8 || vec.qvecAmp()[RefAId] < 1e-8 || vec.qvecAmp()[RefBId] < 1e-8) {
      return;
    }

    if (nmode == 2) {
      for (int ishift = 1; ishift <= 10; ishift++) {
        histosQA.fill(HIST("histShiftV2"), vec.cent(), 2.0 * DetId + 0.5, ishift - 0.5, TMath::Sin(ishift * static_cast<float>(nmode) * TMath::ATan2(vec.qvecIm()[DetInd], vec.qvecRe()[DetInd]) / static_cast<float>(nmode)));
        histosQA.fill(HIST("histShiftV2"), vec.cent(), 2.0 * DetId + 1.5, ishift - 0.5, TMath::Cos(ishift * static_cast<float>(nmode) * TMath::ATan2(vec.qvecIm()[DetInd], vec.qvecRe()[DetInd]) / static_cast<float>(nmode)));

        histosQA.fill(HIST("histShiftV2"), vec.cent(), 2.0 * RefAId + 0.5, ishift - 0.5, TMath::Sin(ishift * static_cast<float>(nmode) * TMath::ATan2(vec.qvecIm()[RefAInd], vec.qvecRe()[RefAInd]) / static_cast<float>(nmode)));
        histosQA.fill(HIST("histShiftV2"), vec.cent(), 2.0 * RefAId + 1.5, ishift - 0.5, TMath::Cos(ishift * static_cast<float>(nmode) * TMath::ATan2(vec.qvecIm()[RefAInd], vec.qvecRe()[RefAInd]) / static_cast<float>(nmode)));

        histosQA.fill(HIST("histShiftV2"), vec.cent(), 2.0 * RefBId + 0.5, ishift - 0.5, TMath::Sin(ishift * static_cast<float>(nmode) * TMath::ATan2(vec.qvecIm()[RefBInd], vec.qvecRe()[RefBInd]) / static_cast<float>(nmode)));
        histosQA.fill(HIST("histShiftV2"), vec.cent(), 2.0 * RefBId + 1.5, ishift - 0.5, TMath::Cos(ishift * static_cast<float>(nmode) * TMath::ATan2(vec.qvecIm()[RefBInd], vec.qvecRe()[RefBInd]) / static_cast<float>(nmode)));
      }
    } else if (nmode == 3) {
      for (int ishift = 1; ishift <= 10; ishift++) {
        histosQA.fill(HIST("histShiftV3"), vec.cent(), 2.0 * DetId + 0.5, ishift - 0.5, TMath::Sin(ishift * static_cast<float>(nmode) * TMath::ATan2(vec.qvecIm()[DetInd], vec.qvecRe()[DetInd]) / static_cast<float>(nmode)));
        histosQA.fill(HIST("histShiftV3"), vec.cent(), 2.0 * DetId + 1.5, ishift - 0.5, TMath::Cos(ishift * static_cast<float>(nmode) * TMath::ATan2(vec.qvecIm()[DetInd], vec.qvecRe()[DetInd]) / static_cast<float>(nmode)));

        histosQA.fill(HIST("histShiftV3"), vec.cent(), 2.0 * RefAId + 0.5, ishift - 0.5, TMath::Sin(ishift * static_cast<float>(nmode) * TMath::ATan2(vec.qvecIm()[RefAInd], vec.qvecRe()[RefAInd]) / static_cast<float>(nmode)));
        histosQA.fill(HIST("histShiftV3"), vec.cent(), 2.0 * RefAId + 1.5, ishift - 0.5, TMath::Cos(ishift * static_cast<float>(nmode) * TMath::ATan2(vec.qvecIm()[RefAInd], vec.qvecRe()[RefAInd]) / static_cast<float>(nmode)));

        histosQA.fill(HIST("histShiftV3"), vec.cent(), 2.0 * RefBId + 0.5, ishift - 0.5, TMath::Sin(ishift * static_cast<float>(nmode) * TMath::ATan2(vec.qvecIm()[RefBInd], vec.qvecRe()[RefBInd]) / static_cast<float>(nmode)));
        histosQA.fill(HIST("histShiftV3"), vec.cent(), 2.0 * RefBId + 1.5, ishift - 0.5, TMath::Cos(ishift * static_cast<float>(nmode) * TMath::ATan2(vec.qvecIm()[RefBInd], vec.qvecRe()[RefBInd]) / static_cast<float>(nmode)));
      }
    } else if (nmode == 4) {
      for (int ishift = 1; ishift <= 10; ishift++) {
        histosQA.fill(HIST("histShiftV4"), vec.cent(), 2.0 * DetId + 0.5, ishift - 0.5, TMath::Sin(ishift * static_cast<float>(nmode) * TMath::ATan2(vec.qvecIm()[DetInd], vec.qvecRe()[DetInd]) / static_cast<float>(nmode)));
        histosQA.fill(HIST("histShiftV4"), vec.cent(), 2.0 * DetId + 1.5, ishift - 0.5, TMath::Cos(ishift * static_cast<float>(nmode) * TMath::ATan2(vec.qvecIm()[DetInd], vec.qvecRe()[DetInd]) / static_cast<float>(nmode)));

        histosQA.fill(HIST("histShiftV4"), vec.cent(), 2.0 * RefAId + 0.5, ishift - 0.5, TMath::Sin(ishift * static_cast<float>(nmode) * TMath::ATan2(vec.qvecIm()[RefAInd], vec.qvecRe()[RefAInd]) / static_cast<float>(nmode)));
        histosQA.fill(HIST("histShiftV4"), vec.cent(), 2.0 * RefAId + 1.5, ishift - 0.5, TMath::Cos(ishift * static_cast<float>(nmode) * TMath::ATan2(vec.qvecIm()[RefAInd], vec.qvecRe()[RefAInd]) / static_cast<float>(nmode)));

        histosQA.fill(HIST("histShiftV4"), vec.cent(), 2.0 * RefBId + 0.5, ishift - 0.5, TMath::Sin(ishift * static_cast<float>(nmode) * TMath::ATan2(vec.qvecIm()[RefBInd], vec.qvecRe()[RefBInd]) / static_cast<float>(nmode)));
        histosQA.fill(HIST("histShiftV4"), vec.cent(), 2.0 * RefBId + 1.5, ishift - 0.5, TMath::Cos(ishift * static_cast<float>(nmode) * TMath::ATan2(vec.qvecIm()[RefBInd], vec.qvecRe()[RefBInd]) / static_cast<float>(nmode)));
      }
    }
  }

  template <typename CollType, typename TrackType>
  void fillHistosFlowWithSC(const CollType& coll, const TrackType& track, int nmode)
  {
    int DetInd = DetId + cfgnTotalSystem * (nmode - 2);
    int RefAInd = RefAId + cfgnTotalSystem * (nmode - 2);
    int RefBInd = RefBId + cfgnTotalSystem * (nmode - 2);

    if (coll.qvecAmp()[DetId] < 1e-8 || coll.qvecAmp()[RefAId] < 1e-8 || coll.qvecAmp()[RefBId] < 1e-8) {
      return;
    }

    for (auto& trk : track) {
      if (!SelTrack(trk)) {
        continue;
      }

      if (std::abs(trk.eta()) > 0.8) {
        continue;
      }
      if (nmode == 2) {
        histosQA.fill(HIST("hist_EP_cos_Det_v2"), coll.cent(), trk.pt(), std::cos(2.0 * (trk.phi() - helperEP.GetEventPlane(coll.qvecShiftedRe()[DetInd], coll.qvecShiftedIm()[DetInd], nmode))));
        histosQA.fill(HIST("hist_EP_sin_Det_v2"), coll.cent(), trk.pt(), std::sin(2.0 * (trk.phi() - helperEP.GetEventPlane(coll.qvecShiftedRe()[DetInd], coll.qvecShiftedIm()[DetInd], nmode))));
        histosQA.fill(HIST("hist_EP_azimuth_Det_v2"), coll.cent(), trk.pt(), TVector2::Phi_0_2pi(2.0 * (trk.phi() - helperEP.GetEventPlane(coll.qvecShiftedRe()[DetInd], coll.qvecShiftedIm()[DetInd], nmode))));

        histosQA.fill(HIST("hist_EP_cos_RefA_v2"), coll.cent(), trk.pt(), std::cos(2.0 * (trk.phi() - helperEP.GetEventPlane(coll.qvecShiftedRe()[RefAInd], coll.qvecShiftedIm()[RefAInd], nmode))));
        histosQA.fill(HIST("hist_EP_sin_RefA_v2"), coll.cent(), trk.pt(), std::sin(2.0 * (trk.phi() - helperEP.GetEventPlane(coll.qvecShiftedRe()[RefAInd], coll.qvecShiftedIm()[RefAInd], nmode))));
        histosQA.fill(HIST("hist_EP_azimuth_RefA_v2"), coll.cent(), trk.pt(), TVector2::Phi_0_2pi(2.0 * (trk.phi() - helperEP.GetEventPlane(coll.qvecShiftedRe()[RefAInd], coll.qvecShiftedIm()[RefAInd], nmode))));

        histosQA.fill(HIST("hist_EP_cos_RefB_v2"), coll.cent(), trk.pt(), std::cos(2.0 * (trk.phi() - helperEP.GetEventPlane(coll.qvecShiftedRe()[RefBInd], coll.qvecShiftedIm()[RefBInd], nmode))));
        histosQA.fill(HIST("hist_EP_sin_RefB_v2"), coll.cent(), trk.pt(), std::sin(2.0 * (trk.phi() - helperEP.GetEventPlane(coll.qvecShiftedRe()[RefBInd], coll.qvecShiftedIm()[RefBInd], nmode))));
        histosQA.fill(HIST("hist_EP_azimuth_RefB_v2"), coll.cent(), trk.pt(), TVector2::Phi_0_2pi(2.0 * (trk.phi() - helperEP.GetEventPlane(coll.qvecShiftedRe()[RefBInd], coll.qvecShiftedIm()[RefBInd], nmode))));
      } else if (nmode == 3) {
        histosQA.fill(HIST("hist_EP_cos_Det_v3"), coll.cent(), trk.pt(), std::cos(3.0 * (trk.phi() - helperEP.GetEventPlane(coll.qvecShiftedRe()[DetInd], coll.qvecShiftedIm()[DetInd], nmode))));
        histosQA.fill(HIST("hist_EP_sin_Det_v3"), coll.cent(), trk.pt(), std::sin(3.0 * (trk.phi() - helperEP.GetEventPlane(coll.qvecShiftedRe()[DetInd], coll.qvecShiftedIm()[DetInd], nmode))));
        histosQA.fill(HIST("hist_EP_azimuth_Det_v3"), coll.cent(), trk.pt(), TVector2::Phi_0_2pi(3.0 * (trk.phi() - helperEP.GetEventPlane(coll.qvecShiftedRe()[DetInd], coll.qvecShiftedIm()[DetInd], nmode))));

        histosQA.fill(HIST("hist_EP_cos_RefA_v3"), coll.cent(), trk.pt(), std::cos(3.0 * (trk.phi() - helperEP.GetEventPlane(coll.qvecShiftedRe()[RefAInd], coll.qvecShiftedIm()[RefAInd], nmode))));
        histosQA.fill(HIST("hist_EP_sin_RefA_v3"), coll.cent(), trk.pt(), std::sin(3.0 * (trk.phi() - helperEP.GetEventPlane(coll.qvecShiftedRe()[RefAInd], coll.qvecShiftedIm()[RefAInd], nmode))));
        histosQA.fill(HIST("hist_EP_azimuth_RefA_v3"), coll.cent(), trk.pt(), TVector2::Phi_0_2pi(3.0 * (trk.phi() - helperEP.GetEventPlane(coll.qvecShiftedRe()[RefAInd], coll.qvecShiftedIm()[RefAInd], nmode))));

        histosQA.fill(HIST("hist_EP_cos_RefB_v3"), coll.cent(), trk.pt(), std::cos(3.0 * (trk.phi() - helperEP.GetEventPlane(coll.qvecShiftedRe()[RefBInd], coll.qvecShiftedIm()[RefBInd], nmode))));
        histosQA.fill(HIST("hist_EP_sin_RefB_v3"), coll.cent(), trk.pt(), std::sin(3.0 * (trk.phi() - helperEP.GetEventPlane(coll.qvecShiftedRe()[RefBInd], coll.qvecShiftedIm()[RefBInd], nmode))));
        histosQA.fill(HIST("hist_EP_azimuth_RefB_v3"), coll.cent(), trk.pt(), TVector2::Phi_0_2pi(3.0 * (trk.phi() - helperEP.GetEventPlane(coll.qvecShiftedRe()[RefBInd], coll.qvecShiftedIm()[RefBInd], nmode))));
      } else if (nmode == 4) {
        histosQA.fill(HIST("hist_EP_cos_Det_v4"), coll.cent(), trk.pt(), std::cos(4.0 * (trk.phi() - helperEP.GetEventPlane(coll.qvecShiftedRe()[DetInd], coll.qvecShiftedIm()[DetInd], nmode))));
        histosQA.fill(HIST("hist_EP_sin_Det_v4"), coll.cent(), trk.pt(), std::sin(4.0 * (trk.phi() - helperEP.GetEventPlane(coll.qvecShiftedRe()[DetInd], coll.qvecShiftedIm()[DetInd], nmode))));
        histosQA.fill(HIST("hist_EP_azimuth_Det_v4"), coll.cent(), trk.pt(), TVector2::Phi_0_2pi(4.0 * (trk.phi() - helperEP.GetEventPlane(coll.qvecShiftedRe()[DetInd], coll.qvecShiftedIm()[DetInd], nmode))));

        histosQA.fill(HIST("hist_EP_cos_RefA_v4"), coll.cent(), trk.pt(), std::cos(4.0 * (trk.phi() - helperEP.GetEventPlane(coll.qvecShiftedRe()[RefAInd], coll.qvecShiftedIm()[RefAInd], nmode))));
        histosQA.fill(HIST("hist_EP_sin_RefA_v4"), coll.cent(), trk.pt(), std::sin(4.0 * (trk.phi() - helperEP.GetEventPlane(coll.qvecShiftedRe()[RefAInd], coll.qvecShiftedIm()[RefAInd], nmode))));
        histosQA.fill(HIST("hist_EP_azimuth_RefA_v4"), coll.cent(), trk.pt(), TVector2::Phi_0_2pi(4.0 * (trk.phi() - helperEP.GetEventPlane(coll.qvecShiftedRe()[RefAInd], coll.qvecShiftedIm()[RefAInd], nmode))));

        histosQA.fill(HIST("hist_EP_cos_RefB_v4"), coll.cent(), trk.pt(), std::cos(4.0 * (trk.phi() - helperEP.GetEventPlane(coll.qvecShiftedRe()[RefBInd], coll.qvecShiftedIm()[RefBInd], nmode))));
        histosQA.fill(HIST("hist_EP_sin_RefB_v4"), coll.cent(), trk.pt(), std::sin(4.0 * (trk.phi() - helperEP.GetEventPlane(coll.qvecShiftedRe()[RefBInd], coll.qvecShiftedIm()[RefBInd], nmode))));
        histosQA.fill(HIST("hist_EP_azimuth_RefB_v4"), coll.cent(), trk.pt(), TVector2::Phi_0_2pi(4.0 * (trk.phi() - helperEP.GetEventPlane(coll.qvecShiftedRe()[RefBInd], coll.qvecShiftedIm()[RefBInd], nmode))));
      }
    }
  }

  template <typename CollType, typename TrackType>
  void fillHistosFlow(const CollType& coll, const TrackType& track, int nmode)
  {
    if (coll.qvecAmp()[DetId] < 1e-8 || coll.qvecAmp()[RefAId] < 1e-8 || coll.qvecAmp()[RefBId] < 1e-8) {
      return;
    }
    int DetInd = DetId * 4 + cfgnTotalSystem * 4 * (nmode - 2);
    int RefAInd = RefAId * 4 + cfgnTotalSystem * 4 * (nmode - 2);
    int RefBInd = RefBId * 4 + cfgnTotalSystem * 4 * (nmode - 2);

    for (auto& trk : track) {
      if (!SelTrack(trk)) {
        continue;
      }

      if (std::abs(trk.eta()) > 0.8) {
        continue;
      }

      if (nmode == 2) {
        histosQA.fill(HIST("hist_EP_cos_Det_v2"), coll.cent(), trk.pt(), std::cos(2.0 * (trk.phi() - helperEP.GetEventPlane(coll.qvecRe()[DetInd + 3], coll.qvecIm()[DetInd + 3], nmode))));
        histosQA.fill(HIST("hist_EP_sin_Det_v2"), coll.cent(), trk.pt(), std::sin(2.0 * (trk.phi() - helperEP.GetEventPlane(coll.qvecRe()[DetInd + 3], coll.qvecIm()[DetInd + 3], nmode))));
        histosQA.fill(HIST("hist_EP_azimuth_Det_v2"), coll.cent(), trk.pt(), TVector2::Phi_0_2pi(2.0 * (trk.phi() - helperEP.GetEventPlane(coll.qvecRe()[DetInd + 3], coll.qvecIm()[DetInd + 3], nmode))));

        histosQA.fill(HIST("hist_EP_cos_RefA_v2"), coll.cent(), trk.pt(), std::cos(2.0 * (trk.phi() - helperEP.GetEventPlane(coll.qvecRe()[RefAInd + 3], coll.qvecIm()[RefAInd + 3], nmode))));
        histosQA.fill(HIST("hist_EP_sin_RefA_v2"), coll.cent(), trk.pt(), std::sin(2.0 * (trk.phi() - helperEP.GetEventPlane(coll.qvecRe()[RefAInd + 3], coll.qvecIm()[RefAInd + 3], nmode))));
        histosQA.fill(HIST("hist_EP_azimuth_RefA_v2"), coll.cent(), trk.pt(), TVector2::Phi_0_2pi(2.0 * (trk.phi() - helperEP.GetEventPlane(coll.qvecRe()[RefAInd + 3], coll.qvecIm()[RefAInd + 3], nmode))));

        histosQA.fill(HIST("hist_EP_cos_RefB_v2"), coll.cent(), trk.pt(), std::cos(2.0 * (trk.phi() - helperEP.GetEventPlane(coll.qvecRe()[RefBInd + 3], coll.qvecIm()[RefBInd + 3], nmode))));
        histosQA.fill(HIST("hist_EP_sin_RefB_v2"), coll.cent(), trk.pt(), std::sin(2.0 * (trk.phi() - helperEP.GetEventPlane(coll.qvecRe()[RefBInd + 3], coll.qvecIm()[RefBInd + 3], nmode))));
        histosQA.fill(HIST("hist_EP_azimuth_RefB_v2"), coll.cent(), trk.pt(), TVector2::Phi_0_2pi(2.0 * (trk.phi() - helperEP.GetEventPlane(coll.qvecRe()[RefBInd + 3], coll.qvecIm()[RefBInd + 3], nmode))));
      } else if (nmode == 3) {
        histosQA.fill(HIST("hist_EP_cos_Det_v3"), coll.cent(), trk.pt(), std::cos(3.0 * (trk.phi() - helperEP.GetEventPlane(coll.qvecRe()[DetInd + 3], coll.qvecIm()[DetInd + 3], nmode))));
        histosQA.fill(HIST("hist_EP_sin_Det_v3"), coll.cent(), trk.pt(), std::sin(3.0 * (trk.phi() - helperEP.GetEventPlane(coll.qvecRe()[DetInd + 3], coll.qvecIm()[DetInd + 3], nmode))));
        histosQA.fill(HIST("hist_EP_azimuth_Det_v3"), coll.cent(), trk.pt(), TVector2::Phi_0_2pi(3.0 * (trk.phi() - helperEP.GetEventPlane(coll.qvecRe()[DetInd + 3], coll.qvecIm()[DetInd + 3], nmode))));

        histosQA.fill(HIST("hist_EP_cos_RefA_v3"), coll.cent(), trk.pt(), std::cos(3.0 * (trk.phi() - helperEP.GetEventPlane(coll.qvecRe()[RefAInd + 3], coll.qvecIm()[RefAInd + 3], nmode))));
        histosQA.fill(HIST("hist_EP_sin_RefA_v3"), coll.cent(), trk.pt(), std::sin(3.0 * (trk.phi() - helperEP.GetEventPlane(coll.qvecRe()[RefAInd + 3], coll.qvecIm()[RefAInd + 3], nmode))));
        histosQA.fill(HIST("hist_EP_azimuth_RefA_v3"), coll.cent(), trk.pt(), TVector2::Phi_0_2pi(3.0 * (trk.phi() - helperEP.GetEventPlane(coll.qvecRe()[RefAInd + 3], coll.qvecIm()[RefAInd + 3], nmode))));

        histosQA.fill(HIST("hist_EP_cos_RefB_v3"), coll.cent(), trk.pt(), std::cos(3.0 * (trk.phi() - helperEP.GetEventPlane(coll.qvecRe()[RefBInd + 3], coll.qvecIm()[RefBInd + 3], nmode))));
        histosQA.fill(HIST("hist_EP_sin_RefB_v3"), coll.cent(), trk.pt(), std::sin(3.0 * (trk.phi() - helperEP.GetEventPlane(coll.qvecRe()[RefBInd + 3], coll.qvecIm()[RefBInd + 3], nmode))));
        histosQA.fill(HIST("hist_EP_azimuth_RefB_v3"), coll.cent(), trk.pt(), TVector2::Phi_0_2pi(3.0 * (trk.phi() - helperEP.GetEventPlane(coll.qvecRe()[RefBInd + 3], coll.qvecIm()[RefBInd + 3], nmode))));
      } else if (nmode == 4) {
        histosQA.fill(HIST("hist_EP_cos_Det_v4"), coll.cent(), trk.pt(), std::cos(4.0 * (trk.phi() - helperEP.GetEventPlane(coll.qvecRe()[DetInd + 3], coll.qvecIm()[DetInd + 3], nmode))));
        histosQA.fill(HIST("hist_EP_sin_Det_v4"), coll.cent(), trk.pt(), std::sin(4.0 * (trk.phi() - helperEP.GetEventPlane(coll.qvecRe()[DetInd + 3], coll.qvecIm()[DetInd + 3], nmode))));
        histosQA.fill(HIST("hist_EP_azimuth_Det_v4"), coll.cent(), trk.pt(), TVector2::Phi_0_2pi(4.0 * (trk.phi() - helperEP.GetEventPlane(coll.qvecRe()[DetInd + 3], coll.qvecIm()[DetInd + 3], nmode))));

        histosQA.fill(HIST("hist_EP_cos_RefA_v4"), coll.cent(), trk.pt(), std::cos(4.0 * (trk.phi() - helperEP.GetEventPlane(coll.qvecRe()[RefAInd + 3], coll.qvecIm()[RefAInd + 3], nmode))));
        histosQA.fill(HIST("hist_EP_sin_RefA_v4"), coll.cent(), trk.pt(), std::sin(4.0 * (trk.phi() - helperEP.GetEventPlane(coll.qvecRe()[RefAInd + 3], coll.qvecIm()[RefAInd + 3], nmode))));
        histosQA.fill(HIST("hist_EP_azimuth_RefA_v4"), coll.cent(), trk.pt(), TVector2::Phi_0_2pi(4.0 * (trk.phi() - helperEP.GetEventPlane(coll.qvecRe()[RefAInd + 3], coll.qvecIm()[RefAInd + 3], nmode))));

        histosQA.fill(HIST("hist_EP_cos_RefB_v4"), coll.cent(), trk.pt(), std::cos(4.0 * (trk.phi() - helperEP.GetEventPlane(coll.qvecRe()[RefBInd + 3], coll.qvecIm()[RefBInd + 3], nmode))));
        histosQA.fill(HIST("hist_EP_sin_RefB_v4"), coll.cent(), trk.pt(), std::sin(4.0 * (trk.phi() - helperEP.GetEventPlane(coll.qvecRe()[RefBInd + 3], coll.qvecIm()[RefBInd + 3], nmode))));
        histosQA.fill(HIST("hist_EP_azimuth_RefB_v4"), coll.cent(), trk.pt(), TVector2::Phi_0_2pi(4.0 * (trk.phi() - helperEP.GetEventPlane(coll.qvecRe()[RefBInd + 3], coll.qvecIm()[RefBInd + 3], nmode))));
      }
    }
  }

  template <typename T>
  void fillHistosQvecWithSC(const T& vec, int nmode)
  {
    int DetInd = DetId + cfgnTotalSystem * (nmode - 2);
    int RefAInd = RefAId + cfgnTotalSystem * (nmode - 2);
    int RefBInd = RefBId + cfgnTotalSystem * (nmode - 2);

    if (vec.qvecAmp()[DetId] < 1e-8 || vec.qvecAmp()[RefAId] < 1e-8 || vec.qvecAmp()[RefBId] < 1e-8) {
      return;
    }

    if (nmode == 2) {
      histosQA.fill(HIST("histQvecFinalV2"), vec.qvecShiftedRe()[DetInd], vec.qvecShiftedIm()[DetInd], vec.cent());
      histosQA.fill(HIST("histEvtPlFinalV2"), helperEP.GetEventPlane(vec.qvecShiftedRe()[DetInd], vec.qvecShiftedIm()[DetInd], nmode), vec.cent());

      histosQA.fill(HIST("histQvecRefAFinalV2"), vec.qvecShiftedRe()[RefAInd], vec.qvecShiftedIm()[RefAInd], vec.cent());
      histosQA.fill(HIST("histEvtPlRefAFinalV2"), helperEP.GetEventPlane(vec.qvecShiftedRe()[RefAInd], vec.qvecShiftedIm()[RefAInd], nmode), vec.cent());

      histosQA.fill(HIST("histQvecRefBFinalV2"), vec.qvecShiftedRe()[RefBInd], vec.qvecShiftedIm()[RefBInd], vec.cent());
      histosQA.fill(HIST("histEvtPlRefBFinalV2"), helperEP.GetEventPlane(vec.qvecShiftedRe()[RefBInd], vec.qvecShiftedIm()[RefBInd], nmode), vec.cent());

      histosQA.fill(HIST("histEvtPlRes_SigRefAV2"), helperEP.GetResolution(helperEP.GetEventPlane(vec.qvecShiftedRe()[DetInd], vec.qvecShiftedIm()[DetInd], nmode), helperEP.GetEventPlane(vec.qvecShiftedRe()[RefAInd], vec.qvecShiftedIm()[RefAInd], nmode), nmode), vec.cent());
      histosQA.fill(HIST("histEvtPlRes_SigRefBV2"), helperEP.GetResolution(helperEP.GetEventPlane(vec.qvecShiftedRe()[DetInd], vec.qvecShiftedIm()[DetInd], nmode), helperEP.GetEventPlane(vec.qvecShiftedRe()[RefBInd], vec.qvecShiftedIm()[RefBInd], nmode), nmode), vec.cent());
      histosQA.fill(HIST("histEvtPlRes_RefARefBV2"), helperEP.GetResolution(helperEP.GetEventPlane(vec.qvecShiftedRe()[RefAInd], vec.qvecShiftedIm()[RefAInd], nmode), helperEP.GetEventPlane(vec.qvecShiftedRe()[RefBInd], vec.qvecShiftedIm()[RefBInd], nmode), nmode), vec.cent());
    } else if (nmode == 3) {
      histosQA.fill(HIST("histQvecFinalV3"), vec.qvecShiftedRe()[DetInd], vec.qvecShiftedIm()[DetInd], vec.cent());
      histosQA.fill(HIST("histEvtPlFinalV3"), helperEP.GetEventPlane(vec.qvecShiftedRe()[DetInd], vec.qvecShiftedIm()[DetInd], nmode), vec.cent());

      histosQA.fill(HIST("histQvecRefAFinalV3"), vec.qvecShiftedRe()[RefAInd], vec.qvecShiftedIm()[RefAInd], vec.cent());
      histosQA.fill(HIST("histEvtPlRefAFinalV3"), helperEP.GetEventPlane(vec.qvecShiftedRe()[RefAInd], vec.qvecShiftedIm()[RefAInd], nmode), vec.cent());

      histosQA.fill(HIST("histQvecRefBFinalV3"), vec.qvecShiftedRe()[RefBInd], vec.qvecShiftedIm()[RefBInd], vec.cent());
      histosQA.fill(HIST("histEvtPlRefBFinalV3"), helperEP.GetEventPlane(vec.qvecShiftedRe()[RefBInd], vec.qvecShiftedIm()[RefBInd], nmode), vec.cent());

      histosQA.fill(HIST("histEvtPlRes_SigRefAV3"), helperEP.GetResolution(helperEP.GetEventPlane(vec.qvecShiftedRe()[DetInd], vec.qvecShiftedIm()[DetInd], nmode), helperEP.GetEventPlane(vec.qvecShiftedRe()[RefAInd], vec.qvecShiftedIm()[RefAInd], nmode), nmode), vec.cent());
      histosQA.fill(HIST("histEvtPlRes_SigRefBV3"), helperEP.GetResolution(helperEP.GetEventPlane(vec.qvecShiftedRe()[DetInd], vec.qvecShiftedIm()[DetInd], nmode), helperEP.GetEventPlane(vec.qvecShiftedRe()[RefBInd], vec.qvecShiftedIm()[RefBInd], nmode), nmode), vec.cent());
      histosQA.fill(HIST("histEvtPlRes_RefARefBV3"), helperEP.GetResolution(helperEP.GetEventPlane(vec.qvecShiftedRe()[RefAInd], vec.qvecShiftedIm()[RefAInd], nmode), helperEP.GetEventPlane(vec.qvecShiftedRe()[RefBInd], vec.qvecShiftedIm()[RefBInd], nmode), nmode), vec.cent());
    } else if (nmode == 4) {
      histosQA.fill(HIST("histQvecFinalV4"), vec.qvecShiftedRe()[DetInd], vec.qvecShiftedIm()[DetInd], vec.cent());
      histosQA.fill(HIST("histEvtPlFinalV4"), helperEP.GetEventPlane(vec.qvecShiftedRe()[DetInd], vec.qvecShiftedIm()[DetInd], nmode), vec.cent());

      histosQA.fill(HIST("histQvecRefAFinalV4"), vec.qvecShiftedRe()[RefAInd], vec.qvecShiftedIm()[RefAInd], vec.cent());
      histosQA.fill(HIST("histEvtPlRefAFinalV4"), helperEP.GetEventPlane(vec.qvecShiftedRe()[RefAInd], vec.qvecShiftedIm()[RefAInd], nmode), vec.cent());

      histosQA.fill(HIST("histQvecRefBFinalV4"), vec.qvecShiftedRe()[RefBInd], vec.qvecShiftedIm()[RefBInd], vec.cent());
      histosQA.fill(HIST("histEvtPlRefBFinalV4"), helperEP.GetEventPlane(vec.qvecShiftedRe()[RefBInd], vec.qvecShiftedIm()[RefBInd], nmode), vec.cent());

      histosQA.fill(HIST("histEvtPlRes_SigRefAV4"), helperEP.GetResolution(helperEP.GetEventPlane(vec.qvecShiftedRe()[DetInd], vec.qvecShiftedIm()[DetInd], nmode), helperEP.GetEventPlane(vec.qvecShiftedRe()[RefAInd], vec.qvecShiftedIm()[RefAInd], nmode), nmode), vec.cent());
      histosQA.fill(HIST("histEvtPlRes_SigRefBV4"), helperEP.GetResolution(helperEP.GetEventPlane(vec.qvecShiftedRe()[DetInd], vec.qvecShiftedIm()[DetInd], nmode), helperEP.GetEventPlane(vec.qvecShiftedRe()[RefBInd], vec.qvecShiftedIm()[RefBInd], nmode), nmode), vec.cent());
      histosQA.fill(HIST("histEvtPlRes_RefARefBV4"), helperEP.GetResolution(helperEP.GetEventPlane(vec.qvecShiftedRe()[RefAInd], vec.qvecShiftedIm()[RefAInd], nmode), helperEP.GetEventPlane(vec.qvecShiftedRe()[RefBInd], vec.qvecShiftedIm()[RefBInd], nmode), nmode), vec.cent());
    }
  }

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
        if (cfgQAOccupancyStudy) {
          histosQA.fill(HIST("histQvecOccUncorV2"), vec.qvecRe()[DetInd], vec.qvecIm()[DetInd], vec.cent(), vec.trackOccupancyInTimeRange());
        }
        if (cfgQAFinal) {
          histosQA.fill(HIST("histQvecFinalV2"), vec.qvecRe()[DetInd + 3], vec.qvecIm()[DetInd + 3], vec.cent());
          histosQA.fill(HIST("histEvtPlFinalV2"), helperEP.GetEventPlane(vec.qvecRe()[DetInd + 3], vec.qvecIm()[DetInd + 3], nmode), vec.cent());
          if (cfgQAOccupancyStudy) {
            histosQA.fill(HIST("histQvecOccFinalV2"), vec.qvecRe()[DetInd + 3], vec.qvecIm()[DetInd + 3], vec.cent(), vec.trackOccupancyInTimeRange());
          }
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
        if (cfgQAOccupancyStudy) {
          histosQA.fill(HIST("histQvecRefAOccUncorV2"), vec.qvecRe()[RefAInd], vec.qvecIm()[RefAInd], vec.cent(), vec.trackOccupancyInTimeRange());
        }
        if (cfgQAFinal) {
          histosQA.fill(HIST("histQvecRefAFinalV2"), vec.qvecRe()[RefAInd + 3], vec.qvecIm()[RefAInd + 3], vec.cent());
          histosQA.fill(HIST("histEvtPlRefAFinalV2"), helperEP.GetEventPlane(vec.qvecRe()[RefAInd + 3], vec.qvecIm()[RefAInd + 3], nmode), vec.cent());
          if (cfgQAOccupancyStudy) {
            histosQA.fill(HIST("histQvecRefAOccFinalV2"), vec.qvecRe()[RefAInd + 3], vec.qvecIm()[RefAInd + 3], vec.cent(), vec.trackOccupancyInTimeRange());
          }
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
        if (cfgQAOccupancyStudy) {
          histosQA.fill(HIST("histQvecRefBOccUncorV2"), vec.qvecRe()[RefBInd], vec.qvecIm()[RefBInd], vec.cent(), vec.trackOccupancyInTimeRange());
        }
        if (cfgQAFinal) {
          histosQA.fill(HIST("histQvecRefBFinalV2"), vec.qvecRe()[RefBInd + 3], vec.qvecIm()[RefBInd + 3], vec.cent());
          histosQA.fill(HIST("histEvtPlRefBFinalV2"), helperEP.GetEventPlane(vec.qvecRe()[RefBInd + 3], vec.qvecIm()[RefBInd + 3], nmode), vec.cent());
          if (cfgQAOccupancyStudy) {
            histosQA.fill(HIST("histQvecRefBOccFinalV2"), vec.qvecRe()[RefBInd + 3], vec.qvecIm()[RefBInd + 3], vec.cent(), vec.trackOccupancyInTimeRange());
          }
          if (cfgQAAll) {
            histosQA.fill(HIST("histQvecRefBRectrV2"), vec.qvecRe()[RefBInd + 1], vec.qvecIm()[RefBInd + 1], vec.cent());
            histosQA.fill(HIST("histQvecRefBTwistV2"), vec.qvecRe()[RefBInd + 2], vec.qvecIm()[RefBInd + 2], vec.cent());

            histosQA.fill(HIST("histEvtPlRefBRectrV2"), helperEP.GetEventPlane(vec.qvecRe()[RefBInd + 1], vec.qvecIm()[RefBInd + 1], nmode), vec.cent());
            histosQA.fill(HIST("histEvtPlRefBTwistV2"), helperEP.GetEventPlane(vec.qvecRe()[RefBInd + 2], vec.qvecIm()[RefBInd + 2], nmode), vec.cent());
          }
        }
      }
      if (vec.qvecAmp()[DetId] > 1e-8 && vec.qvecAmp()[RefAId] > 1e-8 && vec.qvecAmp()[RefBId] > 1e-8 && cfgQAFinal) {
        histosQA.fill(HIST("histQvecRes_SigRefAV2"), vec.qvecRe()[DetInd + 3] * vec.qvecRe()[RefAInd + 3] + vec.qvecIm()[DetInd + 3] * vec.qvecIm()[RefAInd + 3], vec.cent());
        histosQA.fill(HIST("histQvecRes_SigRefBV2"), vec.qvecRe()[DetInd + 3] * vec.qvecRe()[RefBInd + 3] + vec.qvecIm()[DetInd + 3] * vec.qvecIm()[RefBInd + 3], vec.cent());
        histosQA.fill(HIST("histQvecRes_RefARefBV2"), vec.qvecRe()[RefAInd + 3] * vec.qvecRe()[RefBInd + 3] + vec.qvecIm()[RefAInd + 3] * vec.qvecIm()[RefBInd + 3], vec.cent());

        histosQA.fill(HIST("histEvtPlRes_SigRefAV2"), helperEP.GetResolution(helperEP.GetEventPlane(vec.qvecRe()[DetInd + 3], vec.qvecIm()[DetInd + 3], nmode), helperEP.GetEventPlane(vec.qvecRe()[RefAInd + 3], vec.qvecIm()[RefAInd + 3], nmode), nmode), vec.cent());
        histosQA.fill(HIST("histEvtPlRes_SigRefBV2"), helperEP.GetResolution(helperEP.GetEventPlane(vec.qvecRe()[DetInd + 3], vec.qvecIm()[DetInd + 3], nmode), helperEP.GetEventPlane(vec.qvecRe()[RefBInd + 3], vec.qvecIm()[RefBInd + 3], nmode), nmode), vec.cent());
        histosQA.fill(HIST("histEvtPlRes_RefARefBV2"), helperEP.GetResolution(helperEP.GetEventPlane(vec.qvecRe()[RefAInd + 3], vec.qvecIm()[RefAInd + 3], nmode), helperEP.GetEventPlane(vec.qvecRe()[RefBInd + 3], vec.qvecIm()[RefBInd + 3], nmode), nmode), vec.cent());
      }
    } else if (nmode == 3) {
      if (vec.qvecAmp()[DetId] > 1e-8) {
        histosQA.fill(HIST("histQvecUncorV3"), vec.qvecRe()[DetInd], vec.qvecIm()[DetInd], vec.cent());
        histosQA.fill(HIST("histEvtPlUncorV3"), helperEP.GetEventPlane(vec.qvecRe()[DetInd], vec.qvecIm()[DetInd], nmode), vec.cent());
        if (cfgQAOccupancyStudy) {
          histosQA.fill(HIST("histQvecOccUncorV3"), vec.qvecRe()[DetInd], vec.qvecIm()[DetInd], vec.cent(), vec.trackOccupancyInTimeRange());
        }
        if (cfgQAFinal) {
          if (cfgQAOccupancyStudy) {
            histosQA.fill(HIST("histQvecOccFinalV3"), vec.qvecRe()[DetInd + 3], vec.qvecIm()[DetInd + 3], vec.cent(), vec.trackOccupancyInTimeRange());
          }
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
        if (cfgQAOccupancyStudy) {
          histosQA.fill(HIST("histQvecRefAOccUncorV3"), vec.qvecRe()[RefAInd], vec.qvecIm()[RefAInd], vec.cent(), vec.trackOccupancyInTimeRange());
        }
        if (cfgQAFinal) {
          histosQA.fill(HIST("histQvecRefAFinalV3"), vec.qvecRe()[RefAInd + 3], vec.qvecIm()[RefAInd + 3], vec.cent());
          histosQA.fill(HIST("histEvtPlRefAFinalV3"), helperEP.GetEventPlane(vec.qvecRe()[RefAInd + 3], vec.qvecIm()[RefAInd + 3], nmode), vec.cent());
          if (cfgQAOccupancyStudy) {
            histosQA.fill(HIST("histQvecRefAOccFinalV3"), vec.qvecRe()[RefAInd + 3], vec.qvecIm()[RefAInd + 3], vec.cent(), vec.trackOccupancyInTimeRange());
          }
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
        if (cfgQAOccupancyStudy) {
          histosQA.fill(HIST("histQvecRefBOccUncorV3"), vec.qvecRe()[RefBInd], vec.qvecIm()[RefBInd], vec.cent(), vec.trackOccupancyInTimeRange());
        }
        if (cfgQAFinal) {
          histosQA.fill(HIST("histQvecRefBFinalV3"), vec.qvecRe()[RefBInd + 3], vec.qvecIm()[RefBInd + 3], vec.cent());
          histosQA.fill(HIST("histEvtPlRefBFinalV3"), helperEP.GetEventPlane(vec.qvecRe()[RefBInd + 3], vec.qvecIm()[RefBInd + 3], nmode), vec.cent());
          if (cfgQAOccupancyStudy) {
            histosQA.fill(HIST("histQvecRefBOccFinalV3"), vec.qvecRe()[RefBInd + 3], vec.qvecIm()[RefBInd + 3], vec.cent(), vec.trackOccupancyInTimeRange());
          }
          if (cfgQAAll) {
            histosQA.fill(HIST("histQvecRefBRectrV3"), vec.qvecRe()[RefBInd + 1], vec.qvecIm()[RefBInd + 1], vec.cent());
            histosQA.fill(HIST("histQvecRefBTwistV3"), vec.qvecRe()[RefBInd + 2], vec.qvecIm()[RefBInd + 2], vec.cent());

            histosQA.fill(HIST("histEvtPlRefBRectrV3"), helperEP.GetEventPlane(vec.qvecRe()[RefBInd + 1], vec.qvecIm()[RefBInd + 1], nmode), vec.cent());
            histosQA.fill(HIST("histEvtPlRefBTwistV3"), helperEP.GetEventPlane(vec.qvecRe()[RefBInd + 2], vec.qvecIm()[RefBInd + 2], nmode), vec.cent());
          }
        }
      }
      if (vec.qvecAmp()[DetId] > 1e-8 && vec.qvecAmp()[RefAId] > 1e-8 && vec.qvecAmp()[RefBId] > 1e-8 && cfgQAFinal) {
        histosQA.fill(HIST("histQvecRes_SigRefAV3"), vec.qvecRe()[DetInd + 3] * vec.qvecRe()[RefAInd + 3] + vec.qvecIm()[DetInd + 3] * vec.qvecIm()[RefAInd + 3], vec.cent());
        histosQA.fill(HIST("histQvecRes_SigRefBV3"), vec.qvecRe()[DetInd + 3] * vec.qvecRe()[RefBInd + 3] + vec.qvecIm()[DetInd + 3] * vec.qvecIm()[RefBInd + 3], vec.cent());
        histosQA.fill(HIST("histQvecRes_RefARefBV3"), vec.qvecRe()[RefAInd + 3] * vec.qvecRe()[RefBInd + 3] + vec.qvecIm()[RefAInd + 3] * vec.qvecIm()[RefBInd + 3], vec.cent());

        histosQA.fill(HIST("histEvtPlRes_SigRefAV3"), helperEP.GetResolution(helperEP.GetEventPlane(vec.qvecRe()[DetInd + 3], vec.qvecIm()[DetInd + 3], nmode), helperEP.GetEventPlane(vec.qvecRe()[RefAInd + 3], vec.qvecIm()[RefAInd + 3], nmode), nmode), vec.cent());
        histosQA.fill(HIST("histEvtPlRes_SigRefBV3"), helperEP.GetResolution(helperEP.GetEventPlane(vec.qvecRe()[DetInd + 3], vec.qvecIm()[DetInd + 3], nmode), helperEP.GetEventPlane(vec.qvecRe()[RefBInd + 3], vec.qvecIm()[RefBInd + 3], nmode), nmode), vec.cent());
        histosQA.fill(HIST("histEvtPlRes_RefARefBV3"), helperEP.GetResolution(helperEP.GetEventPlane(vec.qvecRe()[RefAInd + 3], vec.qvecIm()[RefAInd + 3], nmode), helperEP.GetEventPlane(vec.qvecRe()[RefBInd + 3], vec.qvecIm()[RefBInd + 3], nmode), nmode), vec.cent());
      }
    } else if (nmode == 4) {
      if (vec.qvecAmp()[DetId] > 1e-8) {
        histosQA.fill(HIST("histQvecUncorV4"), vec.qvecRe()[DetInd], vec.qvecIm()[DetInd], vec.cent());
        histosQA.fill(HIST("histEvtPlUncorV4"), helperEP.GetEventPlane(vec.qvecRe()[DetInd], vec.qvecIm()[DetInd], nmode), vec.cent());
        if (cfgQAOccupancyStudy) {
          histosQA.fill(HIST("histQvecOccUncorV4"), vec.qvecRe()[DetInd], vec.qvecIm()[DetInd], vec.cent(), vec.trackOccupancyInTimeRange());
        }
        if (cfgQAFinal) {
          histosQA.fill(HIST("histQvecFinalV4"), vec.qvecRe()[DetInd + 3], vec.qvecIm()[DetInd + 3], vec.cent());
          histosQA.fill(HIST("histEvtPlFinalV4"), helperEP.GetEventPlane(vec.qvecRe()[DetInd + 3], vec.qvecIm()[DetInd + 3], nmode), vec.cent());
          if (cfgQAOccupancyStudy) {
            histosQA.fill(HIST("histQvecOccFinalV4"), vec.qvecRe()[DetInd + 3], vec.qvecIm()[DetInd + 3], vec.cent(), vec.trackOccupancyInTimeRange());
          }
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
        if (cfgQAOccupancyStudy) {
          histosQA.fill(HIST("histQvecRefAOccUncorV4"), vec.qvecRe()[RefAInd], vec.qvecIm()[RefAInd], vec.cent(), vec.trackOccupancyInTimeRange());
        }
        if (cfgQAFinal) {
          histosQA.fill(HIST("histQvecRefAFinalV4"), vec.qvecRe()[RefAInd + 3], vec.qvecIm()[RefAInd + 3], vec.cent());
          histosQA.fill(HIST("histEvtPlRefAFinalV4"), helperEP.GetEventPlane(vec.qvecRe()[RefAInd + 3], vec.qvecIm()[RefAInd + 3], nmode), vec.cent());
          if (cfgQAOccupancyStudy) {
            histosQA.fill(HIST("histQvecRefAOccFinalV4"), vec.qvecRe()[RefAInd + 3], vec.qvecIm()[RefAInd + 3], vec.cent(), vec.trackOccupancyInTimeRange());
          }
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
        if (cfgQAOccupancyStudy) {
          histosQA.fill(HIST("histQvecRefBOccUncorV4"), vec.qvecRe()[RefBInd], vec.qvecIm()[RefBInd], vec.cent(), vec.trackOccupancyInTimeRange());
        }
        if (cfgQAFinal) {
          histosQA.fill(HIST("histQvecRefBFinalV4"), vec.qvecRe()[RefBInd + 3], vec.qvecIm()[RefBInd + 3], vec.cent());
          histosQA.fill(HIST("histEvtPlRefBFinalV4"), helperEP.GetEventPlane(vec.qvecRe()[RefBInd + 3], vec.qvecIm()[RefBInd + 3], nmode), vec.cent());
          if (cfgQAOccupancyStudy) {
            histosQA.fill(HIST("histQvecRefBOccFinalV4"), vec.qvecRe()[RefBInd + 3], vec.qvecIm()[RefBInd + 3], vec.cent(), vec.trackOccupancyInTimeRange());
          }
          if (cfgQAAll) {
            histosQA.fill(HIST("histQvecRefBRectrV4"), vec.qvecRe()[RefBInd + 1], vec.qvecIm()[RefBInd + 1], vec.cent());
            histosQA.fill(HIST("histQvecRefBTwistV4"), vec.qvecRe()[RefBInd + 2], vec.qvecIm()[RefBInd + 2], vec.cent());

            histosQA.fill(HIST("histEvtPlRefBRectrV4"), helperEP.GetEventPlane(vec.qvecRe()[RefBInd + 1], vec.qvecIm()[RefBInd + 1], nmode), vec.cent());
            histosQA.fill(HIST("histEvtPlRefBTwistV4"), helperEP.GetEventPlane(vec.qvecRe()[RefBInd + 2], vec.qvecIm()[RefBInd + 2], nmode), vec.cent());
          }
        }
      }
      if (vec.qvecAmp()[DetId] > 1e-8 && vec.qvecAmp()[RefAId] > 1e-8 && vec.qvecAmp()[RefBId] > 1e-8 && cfgQAFinal) {
        histosQA.fill(HIST("histQvecRes_SigRefAV4"), vec.qvecRe()[DetInd + 3] * vec.qvecRe()[RefAInd + 3] + vec.qvecIm()[DetInd + 3] * vec.qvecIm()[RefAInd + 3], vec.cent());
        histosQA.fill(HIST("histQvecRes_SigRefBV4"), vec.qvecRe()[DetInd + 3] * vec.qvecRe()[RefBInd + 3] + vec.qvecIm()[DetInd + 3] * vec.qvecIm()[RefBInd + 3], vec.cent());
        histosQA.fill(HIST("histQvecRes_RefARefBV4"), vec.qvecRe()[RefAInd + 3] * vec.qvecRe()[RefBInd + 3] + vec.qvecIm()[RefAInd + 3] * vec.qvecIm()[RefBInd + 3], vec.cent());

        histosQA.fill(HIST("histEvtPlRes_SigRefAV4"), helperEP.GetResolution(helperEP.GetEventPlane(vec.qvecRe()[DetInd + 3], vec.qvecIm()[DetInd + 3], nmode), helperEP.GetEventPlane(vec.qvecRe()[RefAInd + 3], vec.qvecIm()[RefAInd + 3], nmode), nmode), vec.cent());
        histosQA.fill(HIST("histEvtPlRes_SigRefBV4"), helperEP.GetResolution(helperEP.GetEventPlane(vec.qvecRe()[DetInd + 3], vec.qvecIm()[DetInd + 3], nmode), helperEP.GetEventPlane(vec.qvecRe()[RefBInd + 3], vec.qvecIm()[RefBInd + 3], nmode), nmode), vec.cent());
        histosQA.fill(HIST("histEvtPlRes_RefARefBV4"), helperEP.GetResolution(helperEP.GetEventPlane(vec.qvecRe()[RefAInd + 3], vec.qvecIm()[RefAInd + 3], nmode), helperEP.GetEventPlane(vec.qvecRe()[RefBInd + 3], vec.qvecIm()[RefBInd + 3], nmode), nmode), vec.cent());
      }
    }
  }
  void processDefault(MyCollisions::iterator const& qVec, MyTracks const& tracks)
  {
    histosQA.fill(HIST("histCentFull"), qVec.cent());
    if (cfgAddEvtSel) {
      if (std::abs(qVec.posZ()) > 10.)
        return;
      switch (cfgEvtSel) {
        case 0: // Sel8
          if (!qVec.sel8())
            return;
          break;
        case 1: // PbPb standard
          if (!qVec.sel8() || !qVec.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV) || !qVec.selection_bit(aod::evsel::kNoSameBunchPileup))
            return;
          break;
        case 2: // PbPb with pileup
          if (!qVec.sel8() || !qVec.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard) ||
              !qVec.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV) || !qVec.selection_bit(aod::evsel::kNoSameBunchPileup))
            return;
          break;
        case 3: // Small systems (OO, NeNe, pp)
          if (!qVec.sel8() || !qVec.selection_bit(aod::evsel::kNoSameBunchPileup))
            return;
          break;
        default:
          LOGF(warning, "Event selection flag was not found, continuing without basic event selections!\n");
      }
      // Check occupancy
      if (qVec.trackOccupancyInTimeRange() > cfgMaxOccupancy || qVec.trackOccupancyInTimeRange() < cfgMinOccupancy)
        return;
    }
    histosQA.fill(HIST("histCentSelected"), qVec.cent());
    histosQA.fill(HIST("histVtxSelected"), qVec.posZ());

    if (cfgShiftCorPrep) {
      for (uint i = 0; i < cfgnMods->size(); i++) {
        fillHistosShiftCor(qVec, cfgnMods->at(i));
      }
    }

    for (uint i = 0; i < cfgnMods->size(); i++) {
      fillHistosQvec(qVec, cfgnMods->at(i));
      if (cfgQAFinal && cfgQAFlowStudy) {
        fillHistosFlow(qVec, tracks, cfgnMods->at(i));
      }
    }
  }
  PROCESS_SWITCH(qVectorsCorrection, processDefault, "default process", true);

  void processWithSC(MyCollisionsWithSC::iterator const& qVec, MyTracks const& tracks)
  {
    histosQA.fill(HIST("histCentFull"), qVec.cent());
    if (cfgAddEvtSel) {
      if (std::abs(qVec.posZ()) > 10.)
        return;
      switch (cfgEvtSel) {
        case 0: // Sel8
          if (!qVec.sel8())
            return;
          break;
        case 1: // PbPb standard
          if (!qVec.sel8() || !qVec.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV) || !qVec.selection_bit(aod::evsel::kNoSameBunchPileup))
            return;
          break;
        case 2: // PbPb with pileup
          if (!qVec.sel8() || !qVec.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard) ||
              !qVec.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV) || !qVec.selection_bit(aod::evsel::kNoSameBunchPileup))
            return;
          break;
        case 3: // Small systems (OO, NeNe, pp)
          if (!qVec.sel8() || !qVec.selection_bit(aod::evsel::kNoSameBunchPileup))
            return;
          break;
        default:
          LOGF(warning, "Event selection flag was not found, continuing without basic event selections!\n");
      }
      // Check occupancy
      if (qVec.trackOccupancyInTimeRange() > cfgMaxOccupancy || qVec.trackOccupancyInTimeRange() < cfgMinOccupancy)
        return;
    }
    histosQA.fill(HIST("histCentSelected"), qVec.cent());
    histosQA.fill(HIST("histVtxSelected"), qVec.posZ());

    for (uint i = 0; i < cfgnMods->size(); i++) {
      fillHistosQvecWithSC(qVec, cfgnMods->at(i));
      if (cfgQAFinal && cfgQAFlowStudy) {
        fillHistosFlowWithSC(qVec, tracks, cfgnMods->at(i));
      }
    }
  }
  PROCESS_SWITCH(qVectorsCorrection, processWithSC, "process with shift correction", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<qVectorsCorrection>(cfgc)};
}
