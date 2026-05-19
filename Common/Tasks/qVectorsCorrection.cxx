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

#include <TMath.h>
#include <TString.h>
#include <TVector2.h>

#include <sys/types.h>

#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;

using MySpCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Qvectors>;
using MyEseCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::EseQvectors>;
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

  template <typename TColl>
  float getQVecRe(TColl const& collision, int detId)
  {
    if constexpr (std::derived_from<TColl, MyEseCollisions::iterator>) {
      return collision.eseQvecRe()[detId];
    } else {
      return collision.qvecRe()[detId];
    }
  }

  template <typename TColl>
  float getQVecIm(TColl const& collision, int detId)
  {
    if constexpr (std::derived_from<TColl, MyEseCollisions::iterator>) {
      return collision.eseQvecIm()[detId];
    } else {
      return collision.qvecIm()[detId];
    }
  }

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

    histosQA.add("histCentFull", "Centrality distribution for all events", HistType::kTH1F, {axisCent});
    histosQA.add("histCentSelected", "Centrality distribution for selected events", HistType::kTH1F, {axisCent});
    histosQA.add("histVtxSelected", "Z vtx coordinate distribution for valid events", HistType::kTH1F, {axisVertex});

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

    auto aTan2DetInd = TMath::ATan2(getQVecIm(vec, DetInd), getQVecRe(vec, DetInd));
    auto aTan2RefAInd = TMath::ATan2(getQVecIm(vec, RefAInd), getQVecRe(vec, RefAInd));
    auto aTan2RefBInd = TMath::ATan2(getQVecIm(vec, RefBInd), getQVecRe(vec, RefBInd));

    if (nmode == 2) {
      for (int ishift = 1; ishift <= 10; ishift++) {
        histosQA.fill(HIST("histShiftV2"), vec.cent(), 2.0 * DetId + 0.5, ishift - 0.5, TMath::Sin(ishift * aTan2DetInd));
        histosQA.fill(HIST("histShiftV2"), vec.cent(), 2.0 * DetId + 1.5, ishift - 0.5, TMath::Cos(ishift * aTan2DetInd));

        histosQA.fill(HIST("histShiftV2"), vec.cent(), 2.0 * RefAId + 0.5, ishift - 0.5, TMath::Sin(ishift * aTan2RefAInd));
        histosQA.fill(HIST("histShiftV2"), vec.cent(), 2.0 * RefAId + 1.5, ishift - 0.5, TMath::Cos(ishift * aTan2RefAInd));

        histosQA.fill(HIST("histShiftV2"), vec.cent(), 2.0 * RefBId + 0.5, ishift - 0.5, TMath::Sin(ishift * aTan2RefBInd));
        histosQA.fill(HIST("histShiftV2"), vec.cent(), 2.0 * RefBId + 1.5, ishift - 0.5, TMath::Cos(ishift * aTan2RefBInd));
      }
    } else if (nmode == 3) {
      for (int ishift = 1; ishift <= 10; ishift++) {
        histosQA.fill(HIST("histShiftV3"), vec.cent(), 2.0 * DetId + 0.5, ishift - 0.5, TMath::Sin(ishift * aTan2DetInd));
        histosQA.fill(HIST("histShiftV3"), vec.cent(), 2.0 * DetId + 1.5, ishift - 0.5, TMath::Cos(ishift * aTan2DetInd));

        histosQA.fill(HIST("histShiftV3"), vec.cent(), 2.0 * RefAId + 0.5, ishift - 0.5, TMath::Sin(ishift * aTan2RefAInd));
        histosQA.fill(HIST("histShiftV3"), vec.cent(), 2.0 * RefAId + 1.5, ishift - 0.5, TMath::Cos(ishift * aTan2RefAInd));

        histosQA.fill(HIST("histShiftV3"), vec.cent(), 2.0 * RefBId + 0.5, ishift - 0.5, TMath::Sin(ishift * aTan2RefBInd));
        histosQA.fill(HIST("histShiftV3"), vec.cent(), 2.0 * RefBId + 1.5, ishift - 0.5, TMath::Cos(ishift * aTan2RefBInd));
      }
    } else if (nmode == 4) {
      for (int ishift = 1; ishift <= 10; ishift++) {
        histosQA.fill(HIST("histShiftV4"), vec.cent(), 2.0 * DetId + 0.5, ishift - 0.5, TMath::Sin(ishift * aTan2DetInd));
        histosQA.fill(HIST("histShiftV4"), vec.cent(), 2.0 * DetId + 1.5, ishift - 0.5, TMath::Cos(ishift * aTan2DetInd));

        histosQA.fill(HIST("histShiftV4"), vec.cent(), 2.0 * RefAId + 0.5, ishift - 0.5, TMath::Sin(ishift * aTan2RefAInd));
        histosQA.fill(HIST("histShiftV4"), vec.cent(), 2.0 * RefAId + 1.5, ishift - 0.5, TMath::Cos(ishift * aTan2RefAInd));

        histosQA.fill(HIST("histShiftV4"), vec.cent(), 2.0 * RefBId + 0.5, ishift - 0.5, TMath::Sin(ishift * aTan2RefBInd));
        histosQA.fill(HIST("histShiftV4"), vec.cent(), 2.0 * RefBId + 1.5, ishift - 0.5, TMath::Cos(ishift * aTan2RefBInd));
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

      auto trkPhi = trk.phi();
      auto trkPt = trk.pt();
      auto cent = coll.cent();
      auto epDetInd = helperEP.GetEventPlane(getQVecRe(coll, DetInd + 3), getQVecIm(coll, DetInd + 3), nmode);
      auto epRefAInd = helperEP.GetEventPlane(getQVecRe(coll, RefAInd + 3), getQVecIm(coll, RefAInd + 3), nmode);
      auto epRefBInd = helperEP.GetEventPlane(getQVecRe(coll, RefBInd + 3), getQVecIm(coll, RefBInd + 3), nmode);

      if (nmode == 2) {
        histosQA.fill(HIST("hist_EP_cos_Det_v2"), cent, trkPt, std::cos(2.0 * (trkPhi - epDetInd)));
        histosQA.fill(HIST("hist_EP_sin_Det_v2"), cent, trkPt, std::sin(2.0 * (trkPhi - epDetInd)));
        histosQA.fill(HIST("hist_EP_azimuth_Det_v2"), cent, trkPt, TVector2::Phi_0_2pi(2.0 * (trkPhi - epDetInd)));

        histosQA.fill(HIST("hist_EP_cos_RefA_v2"), cent, trkPt, std::cos(2.0 * (trkPhi - epRefAInd)));
        histosQA.fill(HIST("hist_EP_sin_RefA_v2"), cent, trkPt, std::sin(2.0 * (trkPhi - epRefAInd)));
        histosQA.fill(HIST("hist_EP_azimuth_RefA_v2"), cent, trkPt, TVector2::Phi_0_2pi(2.0 * (trkPhi - epRefAInd)));

        histosQA.fill(HIST("hist_EP_cos_RefB_v2"), cent, trkPt, std::cos(2.0 * (trkPhi - epRefBInd)));
        histosQA.fill(HIST("hist_EP_sin_RefB_v2"), cent, trkPt, std::sin(2.0 * (trkPhi - epRefBInd)));
        histosQA.fill(HIST("hist_EP_azimuth_RefB_v2"), cent, trkPt, TVector2::Phi_0_2pi(2.0 * (trkPhi - epRefBInd)));
      } else if (nmode == 3) {
        histosQA.fill(HIST("hist_EP_cos_Det_v3"), cent, trkPt, std::cos(3.0 * (trkPhi - epDetInd)));
        histosQA.fill(HIST("hist_EP_sin_Det_v3"), cent, trkPt, std::sin(3.0 * (trkPhi - epDetInd)));
        histosQA.fill(HIST("hist_EP_azimuth_Det_v3"), cent, trkPt, TVector2::Phi_0_2pi(3.0 * (trkPhi - epDetInd)));

        histosQA.fill(HIST("hist_EP_cos_RefA_v3"), cent, trkPt, std::cos(3.0 * (trkPhi - epRefAInd)));
        histosQA.fill(HIST("hist_EP_sin_RefA_v3"), cent, trkPt, std::sin(3.0 * (trkPhi - epRefAInd)));
        histosQA.fill(HIST("hist_EP_azimuth_RefA_v3"), cent, trkPt, TVector2::Phi_0_2pi(3.0 * (trkPhi - epRefAInd)));

        histosQA.fill(HIST("hist_EP_cos_RefB_v3"), cent, trkPt, std::cos(3.0 * (trkPhi - epRefBInd)));
        histosQA.fill(HIST("hist_EP_sin_RefB_v3"), cent, trkPt, std::sin(3.0 * (trkPhi - epRefBInd)));
        histosQA.fill(HIST("hist_EP_azimuth_RefB_v3"), cent, trkPt, TVector2::Phi_0_2pi(3.0 * (trkPhi - epRefBInd)));
      } else if (nmode == 4) {
        histosQA.fill(HIST("hist_EP_cos_Det_v4"), cent, trkPt, std::cos(4.0 * (trkPhi - epDetInd)));
        histosQA.fill(HIST("hist_EP_sin_Det_v4"), cent, trkPt, std::sin(4.0 * (trkPhi - epDetInd)));
        histosQA.fill(HIST("hist_EP_azimuth_Det_v4"), cent, trkPt, TVector2::Phi_0_2pi(4.0 * (trkPhi - epDetInd)));

        histosQA.fill(HIST("hist_EP_cos_RefA_v4"), cent, trkPt, std::cos(4.0 * (trkPhi - epRefAInd)));
        histosQA.fill(HIST("hist_EP_sin_RefA_v4"), cent, trkPt, std::sin(4.0 * (trkPhi - epRefAInd)));
        histosQA.fill(HIST("hist_EP_azimuth_RefA_v4"), cent, trkPt, TVector2::Phi_0_2pi(4.0 * (trkPhi - epRefAInd)));

        histosQA.fill(HIST("hist_EP_cos_RefB_v4"), cent, trkPt, std::cos(4.0 * (trkPhi - epRefBInd)));
        histosQA.fill(HIST("hist_EP_sin_RefB_v4"), cent, trkPt, std::sin(4.0 * (trkPhi - epRefBInd)));
        histosQA.fill(HIST("hist_EP_azimuth_RefB_v4"), cent, trkPt, TVector2::Phi_0_2pi(4.0 * (trkPhi - epRefBInd)));
      }
    }
  }

  // Definition of all the needed template functions.
  template <typename T>
  void fillHistosQvec(const T& vec, int nmode)
  {
    int DetInd = DetId * 4 + cfgnTotalSystem * 4 * (nmode - 2);
    int RefAInd = RefAId * 4 + cfgnTotalSystem * 4 * (nmode - 2);
    int RefBInd = RefBId * 4 + cfgnTotalSystem * 4 * (nmode - 2);

    auto epDetInd = helperEP.GetEventPlane(getQVecRe(vec, DetInd), getQVecIm(vec, DetInd), nmode);
    auto epDetIndRectr = helperEP.GetEventPlane(getQVecRe(vec, DetInd + 1), getQVecIm(vec, DetInd + 1), nmode);
    auto epDetIndTwist = helperEP.GetEventPlane(getQVecRe(vec, DetInd + 2), getQVecIm(vec, DetInd + 2), nmode);
    auto epDetIndFinal = helperEP.GetEventPlane(getQVecRe(vec, DetInd + 3), getQVecIm(vec, DetInd + 3), nmode);
    auto epRefAInd = helperEP.GetEventPlane(getQVecRe(vec, RefAInd), getQVecIm(vec, RefAInd), nmode);
    auto epRefARectrInd = helperEP.GetEventPlane(getQVecRe(vec, RefAInd + 1), getQVecIm(vec, RefAInd + 1), nmode);
    auto epRefATwistInd = helperEP.GetEventPlane(getQVecRe(vec, RefAInd + 2), getQVecIm(vec, RefAInd + 2), nmode);
    auto epRefAFinalInd = helperEP.GetEventPlane(getQVecRe(vec, RefAInd + 3), getQVecIm(vec, RefAInd + 3), nmode);
    auto epRefBInd = helperEP.GetEventPlane(getQVecRe(vec, RefBInd), getQVecIm(vec, RefBInd), nmode);
    auto epRefBRectrInd = helperEP.GetEventPlane(getQVecRe(vec, RefBInd + 1), getQVecIm(vec, RefBInd + 1), nmode);
    auto epRefBTwistInd = helperEP.GetEventPlane(getQVecRe(vec, RefBInd + 2), getQVecIm(vec, RefBInd + 2), nmode);
    auto epRefBFinalInd = helperEP.GetEventPlane(getQVecRe(vec, RefBInd + 3), getQVecIm(vec, RefBInd + 3), nmode);

    auto occ = vec.trackOccupancyInTimeRange();
    auto cent = vec.cent();

    if (nmode == 2) {
      if (vec.qvecAmp()[DetId] > 1e-8) {
        histosQA.fill(HIST("histQvecUncorV2"), getQVecRe(vec, DetInd), getQVecIm(vec, DetInd), cent);
        histosQA.fill(HIST("histEvtPlUncorV2"), epDetInd, cent);
        if (cfgQAOccupancyStudy) {
          histosQA.fill(HIST("histQvecOccUncorV2"), getQVecRe(vec, DetInd), getQVecIm(vec, DetInd), cent, occ);
        }
        if (cfgQAFinal) {
          histosQA.fill(HIST("histQvecFinalV2"), getQVecRe(vec, DetInd + 3), getQVecIm(vec, DetInd + 3), cent);
          histosQA.fill(HIST("histEvtPlFinalV2"), epDetIndFinal, cent);
          if (cfgQAOccupancyStudy) {
            histosQA.fill(HIST("histQvecOccFinalV2"), getQVecRe(vec, DetInd + 3), getQVecIm(vec, DetInd + 3), cent, occ);
          }
          if (cfgQAAll) {
            histosQA.fill(HIST("histQvecRectrV2"), getQVecRe(vec, DetInd + 1), getQVecIm(vec, DetInd + 1), cent);
            histosQA.fill(HIST("histQvecTwistV2"), getQVecRe(vec, DetInd + 2), getQVecIm(vec, DetInd + 2), cent);

            histosQA.fill(HIST("histEvtPlRectrV2"), epDetIndRectr, cent);
            histosQA.fill(HIST("histEvtPlTwistV2"), epDetIndTwist, cent);
          }
        }
      }
      if (vec.qvecAmp()[RefAId] > 1e-8) {
        histosQA.fill(HIST("histQvecRefAUncorV2"), getQVecRe(vec, RefAInd), getQVecIm(vec, RefAInd), cent);
        histosQA.fill(HIST("histEvtPlRefAUncorV2"), epRefAInd, cent);
        if (cfgQAOccupancyStudy) {
          histosQA.fill(HIST("histQvecRefAOccUncorV2"), getQVecRe(vec, RefAInd), getQVecIm(vec, RefAInd), cent, occ);
        }
        if (cfgQAFinal) {
          histosQA.fill(HIST("histQvecRefAFinalV2"), getQVecRe(vec, RefAInd + 3), getQVecIm(vec, RefAInd + 3), cent);
          histosQA.fill(HIST("histEvtPlRefAFinalV2"), epRefAFinalInd, cent);
          if (cfgQAOccupancyStudy) {
            histosQA.fill(HIST("histQvecRefAOccFinalV2"), getQVecRe(vec, RefAInd + 3), getQVecIm(vec, RefAInd + 3), cent, occ);
          }
          if (cfgQAAll) {
            histosQA.fill(HIST("histQvecRefARectrV2"), getQVecRe(vec, RefAInd + 1), getQVecIm(vec, RefAInd + 1), cent);
            histosQA.fill(HIST("histQvecRefATwistV2"), getQVecRe(vec, RefAInd + 2), getQVecIm(vec, RefAInd + 2), cent);

            histosQA.fill(HIST("histEvtPlRefARectrV2"), epRefARectrInd, cent);
            histosQA.fill(HIST("histEvtPlRefATwistV2"), epRefATwistInd, cent);
          }
        }
      }
      if (vec.qvecAmp()[RefBId] > 1e-8) {
        histosQA.fill(HIST("histQvecRefBUncorV2"), getQVecRe(vec, RefBInd), getQVecIm(vec, RefBInd), cent);
        histosQA.fill(HIST("histEvtPlRefBUncorV2"), epRefBInd, cent);
        if (cfgQAOccupancyStudy) {
          histosQA.fill(HIST("histQvecRefBOccUncorV2"), getQVecRe(vec, RefBInd), getQVecIm(vec, RefBInd), cent, occ);
        }
        if (cfgQAFinal) {
          histosQA.fill(HIST("histQvecRefBFinalV2"), getQVecRe(vec, RefBInd + 3), getQVecIm(vec, RefBInd + 3), cent);
          histosQA.fill(HIST("histEvtPlRefBFinalV2"), epRefBFinalInd, cent);
          if (cfgQAOccupancyStudy) {
            histosQA.fill(HIST("histQvecRefBOccFinalV2"), getQVecRe(vec, RefBInd + 3), getQVecIm(vec, RefBInd + 3), cent, occ);
          }
          if (cfgQAAll) {
            histosQA.fill(HIST("histQvecRefBRectrV2"), getQVecRe(vec, RefBInd + 1), getQVecIm(vec, RefBInd + 1), cent);
            histosQA.fill(HIST("histQvecRefBTwistV2"), getQVecRe(vec, RefBInd + 2), getQVecIm(vec, RefBInd + 2), cent);

            histosQA.fill(HIST("histEvtPlRefBRectrV2"), epRefBRectrInd, cent);
            histosQA.fill(HIST("histEvtPlRefBTwistV2"), epRefBTwistInd, cent);
          }
        }
      }
      if (vec.qvecAmp()[DetId] > 1e-8 && vec.qvecAmp()[RefAId] > 1e-8 && vec.qvecAmp()[RefBId] > 1e-8 && cfgQAFinal) {
        histosQA.fill(HIST("histQvecRes_SigRefAV2"), getQVecRe(vec, DetInd + 3) * getQVecRe(vec, RefAInd + 3) + getQVecIm(vec, DetInd + 3) * getQVecIm(vec, RefAInd + 3), cent);
        histosQA.fill(HIST("histQvecRes_SigRefBV2"), getQVecRe(vec, DetInd + 3) * getQVecRe(vec, RefBInd + 3) + getQVecIm(vec, DetInd + 3) * getQVecIm(vec, RefBInd + 3), cent);
        histosQA.fill(HIST("histQvecRes_RefARefBV2"), getQVecRe(vec, RefAInd + 3) * getQVecRe(vec, RefBInd + 3) + getQVecIm(vec, RefAInd + 3) * getQVecIm(vec, RefBInd + 3), cent);

        histosQA.fill(HIST("histEvtPlRes_SigRefAV2"), helperEP.GetResolution(epDetIndFinal, epRefAFinalInd, nmode), cent);
        histosQA.fill(HIST("histEvtPlRes_SigRefBV2"), helperEP.GetResolution(epDetIndFinal, epRefBFinalInd, nmode), cent);
        histosQA.fill(HIST("histEvtPlRes_RefARefBV2"), helperEP.GetResolution(epRefAFinalInd, epRefBFinalInd, nmode), cent);
      }
    } else if (nmode == 3) {
      if (vec.qvecAmp()[DetId] > 1e-8) {
        histosQA.fill(HIST("histQvecUncorV3"), getQVecRe(vec, DetInd), getQVecIm(vec, DetInd), cent);
        histosQA.fill(HIST("histEvtPlUncorV3"), epDetInd, cent);
        if (cfgQAOccupancyStudy) {
          histosQA.fill(HIST("histQvecOccUncorV3"), getQVecRe(vec, DetInd), getQVecIm(vec, DetInd), cent, occ);
        }
        if (cfgQAFinal) {
          if (cfgQAOccupancyStudy) {
            histosQA.fill(HIST("histQvecOccFinalV3"), getQVecRe(vec, DetInd + 3), getQVecIm(vec, DetInd + 3), cent, occ);
          }
          histosQA.fill(HIST("histQvecFinalV3"), getQVecRe(vec, DetInd + 3), getQVecIm(vec, DetInd + 3), cent);
          histosQA.fill(HIST("histEvtPlFinalV3"), epDetIndFinal, cent);
          if (cfgQAAll) {
            histosQA.fill(HIST("histQvecRectrV3"), getQVecRe(vec, DetInd + 1), getQVecIm(vec, DetInd + 1), cent);
            histosQA.fill(HIST("histQvecTwistV3"), getQVecRe(vec, DetInd + 2), getQVecIm(vec, DetInd + 2), cent);

            histosQA.fill(HIST("histEvtPlRectrV3"), epDetIndRectr, cent);
            histosQA.fill(HIST("histEvtPlTwistV3"), epDetIndTwist, cent);
          }
        }
      }
      if (vec.qvecAmp()[RefAId] > 1e-8) {
        histosQA.fill(HIST("histQvecRefAUncorV3"), getQVecRe(vec, RefAInd), getQVecIm(vec, RefAInd), cent);
        histosQA.fill(HIST("histEvtPlRefAUncorV3"), epRefAInd, cent);
        if (cfgQAOccupancyStudy) {
          histosQA.fill(HIST("histQvecRefAOccUncorV3"), getQVecRe(vec, RefAInd), getQVecIm(vec, RefAInd), cent, occ);
        }
        if (cfgQAFinal) {
          histosQA.fill(HIST("histQvecRefAFinalV3"), getQVecRe(vec, RefAInd + 3), getQVecIm(vec, RefAInd + 3), cent);
          histosQA.fill(HIST("histEvtPlRefAFinalV3"), epRefAFinalInd, cent);
          if (cfgQAOccupancyStudy) {
            histosQA.fill(HIST("histQvecRefAOccFinalV3"), getQVecRe(vec, RefAInd + 3), getQVecIm(vec, RefAInd + 3), cent, occ);
          }
          if (cfgQAAll) {
            histosQA.fill(HIST("histQvecRefARectrV3"), getQVecRe(vec, RefAInd + 1), getQVecIm(vec, RefAInd + 1), cent);
            histosQA.fill(HIST("histQvecRefATwistV3"), getQVecRe(vec, RefAInd + 2), getQVecIm(vec, RefAInd + 2), cent);

            histosQA.fill(HIST("histEvtPlRefARectrV3"), epRefARectrInd, cent);
            histosQA.fill(HIST("histEvtPlRefATwistV3"), epRefATwistInd, cent);
          }
        }
      }
      if (vec.qvecAmp()[RefBId] > 1e-8) {
        histosQA.fill(HIST("histQvecRefBUncorV3"), getQVecRe(vec, RefBInd), getQVecIm(vec, RefBInd), cent);
        histosQA.fill(HIST("histEvtPlRefBUncorV3"), epRefBInd, cent);
        if (cfgQAOccupancyStudy) {
          histosQA.fill(HIST("histQvecRefBOccUncorV3"), getQVecRe(vec, RefBInd), getQVecIm(vec, RefBInd), cent, occ);
        }
        if (cfgQAFinal) {
          histosQA.fill(HIST("histQvecRefBFinalV3"), getQVecRe(vec, RefBInd + 3), getQVecIm(vec, RefBInd + 3), cent);
          histosQA.fill(HIST("histEvtPlRefBFinalV3"), epRefBFinalInd, cent);
          if (cfgQAOccupancyStudy) {
            histosQA.fill(HIST("histQvecRefBOccFinalV3"), getQVecRe(vec, RefBInd + 3), getQVecIm(vec, RefBInd + 3), cent, occ);
          }
          if (cfgQAAll) {
            histosQA.fill(HIST("histQvecRefBRectrV3"), getQVecRe(vec, RefBInd + 1), getQVecIm(vec, RefBInd + 1), cent);
            histosQA.fill(HIST("histQvecRefBTwistV3"), getQVecRe(vec, RefBInd + 2), getQVecIm(vec, RefBInd + 2), cent);

            histosQA.fill(HIST("histEvtPlRefBRectrV3"), epRefBRectrInd, cent);
            histosQA.fill(HIST("histEvtPlRefBTwistV3"), epRefBTwistInd, cent);
          }
        }
      }
      if (vec.qvecAmp()[DetId] > 1e-8 && vec.qvecAmp()[RefAId] > 1e-8 && vec.qvecAmp()[RefBId] > 1e-8 && cfgQAFinal) {
        histosQA.fill(HIST("histQvecRes_SigRefAV3"), getQVecRe(vec, DetInd + 3) * getQVecRe(vec, RefAInd + 3) + getQVecIm(vec, DetInd + 3) * getQVecIm(vec, RefAInd + 3), cent);
        histosQA.fill(HIST("histQvecRes_SigRefBV3"), getQVecRe(vec, DetInd + 3) * getQVecRe(vec, RefBInd + 3) + getQVecIm(vec, DetInd + 3) * getQVecIm(vec, RefBInd + 3), cent);
        histosQA.fill(HIST("histQvecRes_RefARefBV3"), getQVecRe(vec, RefAInd + 3) * getQVecRe(vec, RefBInd + 3) + getQVecIm(vec, RefAInd + 3) * getQVecIm(vec, RefBInd + 3), cent);

        histosQA.fill(HIST("histEvtPlRes_SigRefAV3"), helperEP.GetResolution(epDetIndFinal, epRefAFinalInd, nmode), cent);
        histosQA.fill(HIST("histEvtPlRes_SigRefBV3"), helperEP.GetResolution(epDetIndFinal, epRefBFinalInd, nmode), cent);
        histosQA.fill(HIST("histEvtPlRes_RefARefBV3"), helperEP.GetResolution(epRefAFinalInd, epRefBFinalInd, nmode), cent);
      }
    } else if (nmode == 4) {
      if (vec.qvecAmp()[DetId] > 1e-8) {
        histosQA.fill(HIST("histQvecUncorV4"), getQVecRe(vec, DetInd), getQVecIm(vec, DetInd), cent);
        histosQA.fill(HIST("histEvtPlUncorV4"), epDetInd, cent);
        if (cfgQAOccupancyStudy) {
          histosQA.fill(HIST("histQvecOccUncorV4"), getQVecRe(vec, DetInd), getQVecIm(vec, DetInd), cent, occ);
        }
        if (cfgQAFinal) {
          histosQA.fill(HIST("histQvecFinalV4"), getQVecRe(vec, DetInd + 3), getQVecIm(vec, DetInd + 3), cent);
          histosQA.fill(HIST("histEvtPlFinalV4"), epDetIndFinal, cent);
          if (cfgQAOccupancyStudy) {
            histosQA.fill(HIST("histQvecOccFinalV4"), getQVecRe(vec, DetInd + 3), getQVecIm(vec, DetInd + 3), cent, occ);
          }
          if (cfgQAAll) {
            histosQA.fill(HIST("histQvecRectrV4"), getQVecRe(vec, DetInd + 1), getQVecIm(vec, DetInd + 1), cent);
            histosQA.fill(HIST("histQvecTwistV4"), getQVecRe(vec, DetInd + 2), getQVecIm(vec, DetInd + 2), cent);

            histosQA.fill(HIST("histEvtPlRectrV4"), epDetIndRectr, cent);
            histosQA.fill(HIST("histEvtPlTwistV4"), epDetIndTwist, cent);
          }
        }
      }
      if (vec.qvecAmp()[RefAId] > 1e-8) {
        histosQA.fill(HIST("histQvecRefAUncorV4"), getQVecRe(vec, RefAInd), getQVecIm(vec, RefAInd), cent);
        histosQA.fill(HIST("histEvtPlRefAUncorV4"), epRefAInd, cent);
        if (cfgQAOccupancyStudy) {
          histosQA.fill(HIST("histQvecRefAOccUncorV4"), getQVecRe(vec, RefAInd), getQVecIm(vec, RefAInd), cent, occ);
        }
        if (cfgQAFinal) {
          histosQA.fill(HIST("histQvecRefAFinalV4"), getQVecRe(vec, RefAInd + 3), getQVecIm(vec, RefAInd + 3), cent);
          histosQA.fill(HIST("histEvtPlRefAFinalV4"), epRefAFinalInd, cent);
          if (cfgQAOccupancyStudy) {
            histosQA.fill(HIST("histQvecRefAOccFinalV4"), getQVecRe(vec, RefAInd + 3), getQVecIm(vec, RefAInd + 3), cent, occ);
          }
          if (cfgQAAll) {
            histosQA.fill(HIST("histQvecRefARectrV4"), getQVecRe(vec, RefAInd + 1), getQVecIm(vec, RefAInd + 1), cent);
            histosQA.fill(HIST("histQvecRefATwistV4"), getQVecRe(vec, RefAInd + 2), getQVecIm(vec, RefAInd + 2), cent);

            histosQA.fill(HIST("histEvtPlRefARectrV4"), epRefARectrInd, cent);
            histosQA.fill(HIST("histEvtPlRefATwistV4"), epRefATwistInd, cent);
          }
        }
      }
      if (vec.qvecAmp()[RefBId] > 1e-8) {
        histosQA.fill(HIST("histQvecRefBUncorV4"), getQVecRe(vec, RefBInd), getQVecIm(vec, RefBInd), cent);
        histosQA.fill(HIST("histEvtPlRefBUncorV4"), epRefBInd, cent);
        if (cfgQAOccupancyStudy) {
          histosQA.fill(HIST("histQvecRefBOccUncorV4"), getQVecRe(vec, RefBInd), getQVecIm(vec, RefBInd), cent, occ);
        }
        if (cfgQAFinal) {
          histosQA.fill(HIST("histQvecRefBFinalV4"), getQVecRe(vec, RefBInd + 3), getQVecIm(vec, RefBInd + 3), cent);
          histosQA.fill(HIST("histEvtPlRefBFinalV4"), epRefBFinalInd, cent);
          if (cfgQAOccupancyStudy) {
            histosQA.fill(HIST("histQvecRefBOccFinalV4"), getQVecRe(vec, RefBInd + 3), getQVecIm(vec, RefBInd + 3), cent, occ);
          }
          if (cfgQAAll) {
            histosQA.fill(HIST("histQvecRefBRectrV4"), getQVecRe(vec, RefBInd + 1), getQVecIm(vec, RefBInd + 1), cent);
            histosQA.fill(HIST("histQvecRefBTwistV4"), getQVecRe(vec, RefBInd + 2), getQVecIm(vec, RefBInd + 2), cent);

            histosQA.fill(HIST("histEvtPlRefBRectrV4"), epRefBRectrInd, cent);
            histosQA.fill(HIST("histEvtPlRefBTwistV4"), epRefBTwistInd, cent);
          }
        }
      }
      if (vec.qvecAmp()[DetId] > 1e-8 && vec.qvecAmp()[RefAId] > 1e-8 && vec.qvecAmp()[RefBId] > 1e-8 && cfgQAFinal) {
        histosQA.fill(HIST("histQvecRes_SigRefAV4"), getQVecRe(vec, DetInd + 3) * getQVecRe(vec, RefAInd + 3) + getQVecIm(vec, DetInd + 3) * getQVecIm(vec, RefAInd + 3), cent);
        histosQA.fill(HIST("histQvecRes_SigRefBV4"), getQVecRe(vec, DetInd + 3) * getQVecRe(vec, RefBInd + 3) + getQVecIm(vec, DetInd + 3) * getQVecIm(vec, RefBInd + 3), cent);
        histosQA.fill(HIST("histQvecRes_RefARefBV4"), getQVecRe(vec, RefAInd + 3) * getQVecRe(vec, RefBInd + 3) + getQVecIm(vec, RefAInd + 3) * getQVecIm(vec, RefBInd + 3), cent);

        histosQA.fill(HIST("histEvtPlRes_SigRefAV4"), helperEP.GetResolution(epDetIndFinal, epRefAFinalInd, nmode), cent);
        histosQA.fill(HIST("histEvtPlRes_SigRefBV4"), helperEP.GetResolution(epDetIndFinal, epRefBFinalInd, nmode), cent);
        histosQA.fill(HIST("histEvtPlRes_RefARefBV4"), helperEP.GetResolution(epRefAFinalInd, epRefBFinalInd, nmode), cent);
      }
    }
  }

  template <typename TCollision>
  void processColls(TCollision const& flowQVec, MyTracks const& tracks)
  {
    histosQA.fill(HIST("histCentFull"), flowQVec.cent());
    if (cfgAddEvtSel) {
      if (std::abs(flowQVec.posZ()) > 10.)
        return;
      switch (cfgEvtSel) {
        case 0: // Sel8
          if (!flowQVec.sel8())
            return;
          break;
        case 1: // PbPb standard
          if (!flowQVec.sel8() || !flowQVec.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV) || !flowQVec.selection_bit(aod::evsel::kNoSameBunchPileup))
            return;
          break;
        case 2: // PbPb with pileup
          if (!flowQVec.sel8() || !flowQVec.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard) ||
              !flowQVec.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV) || !flowQVec.selection_bit(aod::evsel::kNoSameBunchPileup))
            return;
          break;
        case 3: // Small systems (OO, NeNe, pp)
          if (!flowQVec.sel8() || !flowQVec.selection_bit(aod::evsel::kNoSameBunchPileup))
            return;
          break;
        default:
          LOGF(warning, "Event selection flag was not found, continuing without basic event selections!\n");
      }
      // Check occupancy
      if (flowQVec.trackOccupancyInTimeRange() > cfgMaxOccupancy || flowQVec.trackOccupancyInTimeRange() < cfgMinOccupancy)
        return;
    }
    histosQA.fill(HIST("histCentSelected"), flowQVec.cent());
    histosQA.fill(HIST("histVtxSelected"), flowQVec.posZ());

    if (cfgShiftCorPrep) {
      for (uint i = 0; i < cfgnMods->size(); i++) {
        fillHistosShiftCor(flowQVec, cfgnMods->at(i));
      }
    }
    for (uint i = 0; i < cfgnMods->size(); i++) {
      fillHistosQvec(flowQVec, cfgnMods->at(i));
      if (cfgQAFinal && cfgQAFlowStudy) {
        fillHistosFlow(flowQVec, tracks, cfgnMods->at(i));
      }
    }
  }

  void processSp(MySpCollisions::iterator const& spQVec, MyTracks const& tracks)
  {
    processColls(spQVec, tracks);
  }
  PROCESS_SWITCH(qVectorsCorrection, processSp, "process SP", true);

  void processEse(MyEseCollisions::iterator const& eseQVec, MyTracks const& tracks)
  {
    processColls(eseQVec, tracks);
  }

  PROCESS_SWITCH(qVectorsCorrection, processEse, "process Ese", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<qVectorsCorrection>(cfgc)};
}
