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

/// \file flowEseTask.cxx
/// \brief Task for flow and event shape engineering correlation with other observation.
/// \author Alice Collaboration
/// \since 2023-05-15
/// \version 1.0
///
/// This task calculates flow and event shape engineering
/// using Q-vector and event plane methods.

// C++/ROOT includes.
#include <TComplex.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TMath.h>
#include <TVector2.h>

#include <chrono>
#include <string>
#include <vector>

// o2Physics includes.
#include "Common/Core/EventPlaneHelper.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Qvectors.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/StaticFor.h"
#include "Framework/runDataProcessing.h"

// o2 includes.

using namespace o2;
using namespace o2::framework;

using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentFV0As, aod::QvectorFT0CVecs, aod::QvectorTPCposVecs, aod::QvectorTPCnegVecs>;
using MyTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TrackSelectionExtension>;
using BCsWithRun3Matchings = soa::Join<aod::BCs, aod::Timestamps, aod::Run3MatchedToBCSparse>;

struct FlowEseTask {
  HistogramRegistry histosQA{"histosQA", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  Configurable<std::vector<int>> cfgNmods{"cfgNmods", {2}, "Modulation of interest"};
  Configurable<std::string> cfgDetName{"cfgDetName", "FT0C", "The name of detector to be analyzed"};
  Configurable<std::string> cfgRefAName{"cfgRefAName", "TPCpos", "The name of detector for reference A"};
  Configurable<std::string> cfgRefBName{"cfgRefBName", "TPCneg", "The name of detector for reference B"};

  Configurable<float> cfgMinPt{"cfgMinPt", 0.15f, "Minimum transverse momentum for charged track"};
  Configurable<float> cfgMaxEta{"cfgMaxEta", 0.8f, "Maximum pseudorapidiy for charged track"};
  Configurable<float> cfgMaxDCArToPVcut{"cfgMaxDCArToPVcut", 0.1f, "Maximum transverse DCA"};
  Configurable<float> cfgMaxDCAzToPVcut{"cfgMaxDCAzToPVcut", 1.0f, "Maximum longitudinal DCA"};

  ConfigurableAxis cfgAxisQvecF{"cfgAxisQvecF", {300, -1, 1}, ""};
  ConfigurableAxis cfgAxisQvec{"cfgAxisQvec", {100, -3, 3}, ""};
  ConfigurableAxis cfgAxisCent{"cfgAxisCent", {100, 0, 100}, ""};

  ConfigurableAxis cfgAxisCos{"cfgAxisCos", {102, -1.02, 1.02}, ""};
  ConfigurableAxis cfgAxisPt{"cfgAxisPt", {100, 0, 10}, ""};
  ConfigurableAxis cfgAxisCentMerged{"cfgAxisCentMerged", {20, 0, 100}, ""};
  ConfigurableAxis cfgAxisMultNum{"cfgAxisMultNum", {300, 0, 2700}, ""};
  ConfigurableAxis cfgaxisQ{"cfgaxisQ", {1000, 0, 1000}, ""};

  static constexpr float kMinAmplitudeThreshold = 1e-4f;
  static constexpr int kDefaultModulation = 2;
  static constexpr std::array<double, 11> kCent = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
  static constexpr std::array<double, 20> kSeparator = {125, 204, 110, 172, 92, 143, 74, 116, 57, 92, 43, 70, 31, 50, 21, 34, 14, 22, 5, 9};

  EventPlaneHelper helperEP;

  void init(InitContext const&)
  {
    AxisSpec axisCent{cfgAxisCent, "centrality"};
    AxisSpec axisQvec{cfgAxisQvec, "Q"};
    AxisSpec axisQvecF{cfgAxisQvecF, "Q"};
    AxisSpec axisEvtPl = {100, -1.0 * constants::math::PI, constants::math::PI};

    AxisSpec axisCos{cfgAxisCos, "angle function"};
    AxisSpec axisPt{cfgAxisPt, "trasverse momentum"};
    AxisSpec axisCentMerged{cfgAxisCentMerged, "merged centrality"};
    AxisSpec axisMultNum{cfgAxisMultNum, "statistic of mult"};
    AxisSpec axisQ{cfgaxisQ, "result of q2"};

    histosQA.add(Form("histQvecV2"), "", {HistType::kTH3F, {axisQvecF, axisQvecF, axisCent}});
    histosQA.add(Form("histQvecCent"), "", {HistType::kTH2F, {axisQ, axisCent}});
    histosQA.add(Form("histEvtPlV2"), "", {HistType::kTH2F, {axisEvtPl, axisCent}});
    histosQA.add(Form("histQvecRes_SigRefAV2"), "", {HistType::kTH2F, {axisQvecF, axisCent}});
    histosQA.add(Form("histCosDetV2"), "", {HistType::kTH3F, {axisCentMerged, axisPt, axisCos}});
    histosQA.add(Form("histCosDetV2_0010Left"), "", {HistType::kTH3F, {axisQvecF, axisPt, axisCos}});
    histosQA.add(Form("histCosDetV2_0010Mid"), "", {HistType::kTH3F, {axisQvecF, axisPt, axisCos}});
    histosQA.add(Form("histCosDetV2_0010Right"), "", {HistType::kTH3F, {axisQvecF, axisPt, axisCos}});
    histosQA.add(Form("histCosDetV2_1020Left"), "", {HistType::kTH3F, {axisQvecF, axisPt, axisCos}});
    histosQA.add(Form("histCosDetV2_1020Mid"), "", {HistType::kTH3F, {axisQvecF, axisPt, axisCos}});
    histosQA.add(Form("histCosDetV2_1020Right"), "", {HistType::kTH3F, {axisQvecF, axisPt, axisCos}});
    histosQA.add(Form("histCosDetV2_2030Left"), "", {HistType::kTH3F, {axisQvecF, axisPt, axisCos}});
    histosQA.add(Form("histCosDetV2_2030Mid"), "", {HistType::kTH3F, {axisQvecF, axisPt, axisCos}});
    histosQA.add(Form("histCosDetV2_2030Right"), "", {HistType::kTH3F, {axisQvecF, axisPt, axisCos}});
    histosQA.add(Form("histCosDetV2_3040Left"), "", {HistType::kTH3F, {axisQvecF, axisPt, axisCos}});
    histosQA.add(Form("histCosDetV2_3040Mid"), "", {HistType::kTH3F, {axisQvecF, axisPt, axisCos}});
    histosQA.add(Form("histCosDetV2_3040Right"), "", {HistType::kTH3F, {axisQvecF, axisPt, axisCos}});
    histosQA.add(Form("histCosDetV2_4050Left"), "", {HistType::kTH3F, {axisQvecF, axisPt, axisCos}});
    histosQA.add(Form("histCosDetV2_4050Mid"), "", {HistType::kTH3F, {axisQvecF, axisPt, axisCos}});
    histosQA.add(Form("histCosDetV2_4050Right"), "", {HistType::kTH3F, {axisQvecF, axisPt, axisCos}});
    histosQA.add(Form("histCosDetV2_5060Left"), "", {HistType::kTH3F, {axisQvecF, axisPt, axisCos}});
    histosQA.add(Form("histCosDetV2_5060Mid"), "", {HistType::kTH3F, {axisQvecF, axisPt, axisCos}});
    histosQA.add(Form("histCosDetV2_5060Right"), "", {HistType::kTH3F, {axisQvecF, axisPt, axisCos}});
    histosQA.add(Form("histCosDetV2_6070Left"), "", {HistType::kTH3F, {axisQvecF, axisPt, axisCos}});
    histosQA.add(Form("histCosDetV2_6070Mid"), "", {HistType::kTH3F, {axisQvecF, axisPt, axisCos}});
    histosQA.add(Form("histCosDetV2_6070Right"), "", {HistType::kTH3F, {axisQvecF, axisPt, axisCos}});
    histosQA.add(Form("histCosDetV2_7080Left"), "", {HistType::kTH3F, {axisQvecF, axisPt, axisCos}});
    histosQA.add(Form("histCosDetV2_7080Mid"), "", {HistType::kTH3F, {axisQvecF, axisPt, axisCos}});
    histosQA.add(Form("histCosDetV2_7080Right"), "", {HistType::kTH3F, {axisQvecF, axisPt, axisCos}});
    histosQA.add(Form("histCosDetV2_8090Left"), "", {HistType::kTH3F, {axisQvecF, axisPt, axisCos}});
    histosQA.add(Form("histCosDetV2_8090Mid"), "", {HistType::kTH3F, {axisQvecF, axisPt, axisCos}});
    histosQA.add(Form("histCosDetV2_8090Right"), "", {HistType::kTH3F, {axisQvecF, axisPt, axisCos}});
    histosQA.add(Form("histCosDetV2_9010Left"), "", {HistType::kTH3F, {axisQvecF, axisPt, axisCos}});
    histosQA.add(Form("histCosDetV2_9010Mid"), "", {HistType::kTH3F, {axisQvecF, axisPt, axisCos}});
    histosQA.add(Form("histCosDetV2_9010Right"), "", {HistType::kTH3F, {axisQvecF, axisPt, axisCos}});
    histosQA.add(Form("histMult_Cent"), "", {HistType::kTH2F, {axisMultNum, axisCent}});
  }

  template <typename CollType>
  bool selectEvent(CollType const& collision)
  {
    if (!collision.sel8()) {
      return false;
    }
    if (!collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    if (!collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    if (!collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      return false;
    }

    return true;
  }

  template <typename TrackType>
  bool selectTrack(TrackType const& track)
  {
    if (track.pt() < cfgMinPt) {
      return false;
    }
    if (std::abs(track.eta()) > cfgMaxEta) {
      return false;
    }
    if (!track.passedITSNCls()) {
      return false;
    }
    if (!track.passedITSChi2NDF()) {
      return false;
    }
    if (!track.passedITSHits()) {
      return false;
    }
    if (!track.passedTPCCrossedRowsOverNCls()) {
      return false;
    }
    if (!track.passedTPCChi2NDF()) {
      return false;
    }
    if (!track.passedDCAxy()) {
      return false;
    }
    if (!track.passedDCAz()) {
      return false;
    }

    return true;
  }

  template <typename CollType>
  void fillHistosQvec(CollType const& collision, int nmode)
  {
    if (nmode == kDefaultModulation) {
      histosQA.fill(HIST("histQvecV2"), collision.qvecFT0CReVec()[0], collision.qvecFT0CImVec()[0],
                    collision.centFT0C());
      histosQA.fill(HIST("histQvecCent"), std::sqrt(collision.qvecFT0CReVec()[0] * collision.qvecFT0CReVec()[0] + collision.qvecFT0CImVec()[0] * collision.qvecFT0CImVec()[0]) * std::sqrt(collision.sumAmplFT0C()), collision.centFT0C());
      histosQA.fill(HIST("histEvtPlV2"),
                    helperEP.GetEventPlane(collision.qvecFT0CReVec()[0], collision.qvecFT0CImVec()[0], nmode),
                    collision.centFT0C());
      histosQA.fill(HIST("histQvecRes_SigRefAV2"), helperEP.GetResolution(helperEP.GetEventPlane(collision.qvecFT0CReVec()[0], collision.qvecFT0CImVec()[0], nmode), helperEP.GetEventPlane(collision.qvecTPCposReVec()[0], collision.qvecTPCposImVec()[0], nmode), nmode),
                    collision.centFT0C());
      histosQA.fill(HIST("histMult_Cent"), collision.sumAmplFT0C(), collision.centFT0C());
    }
  }

  template <typename CollType, typename TrackType>
  void fillHistosFlow(CollType const& collision, TrackType const& tracks, int nmode)
  {
    double q2 = std::sqrt(collision.qvecFT0CReVec()[0] * collision.qvecFT0CReVec()[0] + collision.qvecFT0CImVec()[0] * collision.qvecFT0CImVec()[0]) * std::sqrt(collision.sumAmplFT0C());
    if (collision.sumAmplFT0C() < kMinAmplitudeThreshold) {
      return;
    }
    for (auto const& trk : tracks) {
      if (!selectTrack(trk)) {
        continue;
      }
      if (nmode == kDefaultModulation) {
        histosQA.fill(HIST("histCosDetV2"), collision.centFT0C(), trk.pt(),
                      std::cos(static_cast<float>(nmode) * (trk.phi() -
                                                            helperEP.GetEventPlane(collision.qvecFT0CReVec()[0],
                                                                                   collision.qvecFT0CImVec()[0],
                                                                                   nmode))));
        if (collision.centFT0C() > kCent[0] && collision.centFT0C() <= kCent[1] && q2 < kSeparator[0]) {
          histosQA.fill(HIST("histCosDetV2_0010Left"), helperEP.GetResolution(helperEP.GetEventPlane(collision.qvecFT0CReVec()[0], collision.qvecFT0CImVec()[0], nmode), helperEP.GetEventPlane(collision.qvecTPCposReVec()[0], collision.qvecTPCposImVec()[0], nmode), nmode), trk.pt(),
                        std::cos(static_cast<float>(nmode) * (trk.phi() -
                                                              helperEP.GetEventPlane(collision.qvecFT0CReVec()[0],
                                                                                     collision.qvecFT0CImVec()[0],
                                                                                     nmode))));
        }
        if (collision.centFT0C() > kCent[0] && collision.centFT0C() <= kCent[1] && q2 > kSeparator[0] && q2 <= kSeparator[1]) {
          histosQA.fill(HIST("histCosDetV2_0010Mid"), helperEP.GetResolution(helperEP.GetEventPlane(collision.qvecFT0CReVec()[0], collision.qvecFT0CImVec()[0], nmode), helperEP.GetEventPlane(collision.qvecTPCposReVec()[0], collision.qvecTPCposImVec()[0], nmode), nmode), trk.pt(),
                        std::cos(static_cast<float>(nmode) * (trk.phi() -
                                                              helperEP.GetEventPlane(collision.qvecFT0CReVec()[0],
                                                                                     collision.qvecFT0CImVec()[0],
                                                                                     nmode))));
        }
        if (collision.centFT0C() > kCent[0] && collision.centFT0C() <= kCent[1] && q2 > kSeparator[1]) {
          histosQA.fill(HIST("histCosDetV2_0010Right"), helperEP.GetResolution(helperEP.GetEventPlane(collision.qvecFT0CReVec()[0], collision.qvecFT0CImVec()[0], nmode), helperEP.GetEventPlane(collision.qvecTPCposReVec()[0], collision.qvecTPCposImVec()[0], nmode), nmode), trk.pt(),
                        std::cos(static_cast<float>(nmode) * (trk.phi() -
                                                              helperEP.GetEventPlane(collision.qvecFT0CReVec()[0],
                                                                                     collision.qvecFT0CImVec()[0],
                                                                                     nmode))));
        }
        if (collision.centFT0C() > kCent[1] && collision.centFT0C() <= kCent[2] && q2 < kSeparator[2]) {
          histosQA.fill(HIST("histCosDetV2_1020Left"), helperEP.GetResolution(helperEP.GetEventPlane(collision.qvecFT0CReVec()[0], collision.qvecFT0CImVec()[0], nmode), helperEP.GetEventPlane(collision.qvecTPCposReVec()[0], collision.qvecTPCposImVec()[0], nmode), nmode), trk.pt(),
                        std::cos(static_cast<float>(nmode) * (trk.phi() -
                                                              helperEP.GetEventPlane(collision.qvecFT0CReVec()[0],
                                                                                     collision.qvecFT0CImVec()[0],
                                                                                     nmode))));
        }
        if (collision.centFT0C() > kCent[1] && collision.centFT0C() <= kCent[2] && q2 > kSeparator[2] && q2 <= kSeparator[3]) {
          histosQA.fill(HIST("histCosDetV2_1020Mid"), helperEP.GetResolution(helperEP.GetEventPlane(collision.qvecFT0CReVec()[0], collision.qvecFT0CImVec()[0], nmode), helperEP.GetEventPlane(collision.qvecTPCposReVec()[0], collision.qvecTPCposImVec()[0], nmode), nmode), trk.pt(),
                        std::cos(static_cast<float>(nmode) * (trk.phi() -
                                                              helperEP.GetEventPlane(collision.qvecFT0CReVec()[0],
                                                                                     collision.qvecFT0CImVec()[0],
                                                                                     nmode))));
        }
        if (collision.centFT0C() > kCent[1] && collision.centFT0C() <= kCent[2] && q2 > kSeparator[3]) {
          histosQA.fill(HIST("histCosDetV2_1020Right"), helperEP.GetResolution(helperEP.GetEventPlane(collision.qvecFT0CReVec()[0], collision.qvecFT0CImVec()[0], nmode), helperEP.GetEventPlane(collision.qvecTPCposReVec()[0], collision.qvecTPCposImVec()[0], nmode), nmode), trk.pt(),
                        std::cos(static_cast<float>(nmode) * (trk.phi() -
                                                              helperEP.GetEventPlane(collision.qvecFT0CReVec()[0],
                                                                                     collision.qvecFT0CImVec()[0],
                                                                                     nmode))));
        }
        if (collision.centFT0C() > kCent[2] && collision.centFT0C() <= kCent[3] && q2 < kSeparator[4]) {
          histosQA.fill(HIST("histCosDetV2_2030Left"), helperEP.GetResolution(helperEP.GetEventPlane(collision.qvecFT0CReVec()[0], collision.qvecFT0CImVec()[0], nmode), helperEP.GetEventPlane(collision.qvecTPCposReVec()[0], collision.qvecTPCposImVec()[0], nmode), nmode), trk.pt(),
                        std::cos(static_cast<float>(nmode) * (trk.phi() -
                                                              helperEP.GetEventPlane(collision.qvecFT0CReVec()[0],
                                                                                     collision.qvecFT0CImVec()[0],
                                                                                     nmode))));
        }
        if (collision.centFT0C() > kCent[2] && collision.centFT0C() <= kCent[3] && q2 > kSeparator[4] && q2 <= kSeparator[5]) {
          histosQA.fill(HIST("histCosDetV2_2030Mid"), helperEP.GetResolution(helperEP.GetEventPlane(collision.qvecFT0CReVec()[0], collision.qvecFT0CImVec()[0], nmode), helperEP.GetEventPlane(collision.qvecTPCposReVec()[0], collision.qvecTPCposImVec()[0], nmode), nmode), trk.pt(),
                        std::cos(static_cast<float>(nmode) * (trk.phi() -
                                                              helperEP.GetEventPlane(collision.qvecFT0CReVec()[0],
                                                                                     collision.qvecFT0CImVec()[0],
                                                                                     nmode))));
        }
        if (collision.centFT0C() > kCent[2] && collision.centFT0C() <= kCent[3] && q2 > kSeparator[5]) {
          histosQA.fill(HIST("histCosDetV2_2030Right"), helperEP.GetResolution(helperEP.GetEventPlane(collision.qvecFT0CReVec()[0], collision.qvecFT0CImVec()[0], nmode), helperEP.GetEventPlane(collision.qvecTPCposReVec()[0], collision.qvecTPCposImVec()[0], nmode), nmode), trk.pt(),
                        std::cos(static_cast<float>(nmode) * (trk.phi() -
                                                              helperEP.GetEventPlane(collision.qvecFT0CReVec()[0],
                                                                                     collision.qvecFT0CImVec()[0],
                                                                                     nmode))));
        }
        if (collision.centFT0C() > kCent[3] && collision.centFT0C() <= kCent[4] && q2 < kSeparator[6]) {
          histosQA.fill(HIST("histCosDetV2_3040Left"), helperEP.GetResolution(helperEP.GetEventPlane(collision.qvecFT0CReVec()[0], collision.qvecFT0CImVec()[0], nmode), helperEP.GetEventPlane(collision.qvecTPCposReVec()[0], collision.qvecTPCposImVec()[0], nmode), nmode), trk.pt(),
                        std::cos(static_cast<float>(nmode) * (trk.phi() -
                                                              helperEP.GetEventPlane(collision.qvecFT0CReVec()[0],
                                                                                     collision.qvecFT0CImVec()[0],
                                                                                     nmode))));
        }
        if (collision.centFT0C() > kCent[3] && collision.centFT0C() <= kCent[4] && q2 > kSeparator[6] && q2 <= kSeparator[7]) {
          histosQA.fill(HIST("histCosDetV2_3040Mid"), helperEP.GetResolution(helperEP.GetEventPlane(collision.qvecFT0CReVec()[0], collision.qvecFT0CImVec()[0], nmode), helperEP.GetEventPlane(collision.qvecTPCposReVec()[0], collision.qvecTPCposImVec()[0], nmode), nmode), trk.pt(),
                        std::cos(static_cast<float>(nmode) * (trk.phi() -
                                                              helperEP.GetEventPlane(collision.qvecFT0CReVec()[0],
                                                                                     collision.qvecFT0CImVec()[0],
                                                                                     nmode))));
        }
        if (collision.centFT0C() > kCent[3] && collision.centFT0C() <= kCent[4] && q2 > kSeparator[7]) {
          histosQA.fill(HIST("histCosDetV2_3040Right"), helperEP.GetResolution(helperEP.GetEventPlane(collision.qvecFT0CReVec()[0], collision.qvecFT0CImVec()[0], nmode), helperEP.GetEventPlane(collision.qvecTPCposReVec()[0], collision.qvecTPCposImVec()[0], nmode), nmode), trk.pt(),
                        std::cos(static_cast<float>(nmode) * (trk.phi() -
                                                              helperEP.GetEventPlane(collision.qvecFT0CReVec()[0],
                                                                                     collision.qvecFT0CImVec()[0],
                                                                                     nmode))));
        }
        if (collision.centFT0C() > kCent[4] && collision.centFT0C() <= kCent[5] && q2 < kSeparator[8]) {
          histosQA.fill(HIST("histCosDetV2_4050Left"), helperEP.GetResolution(helperEP.GetEventPlane(collision.qvecFT0CReVec()[0], collision.qvecFT0CImVec()[0], nmode), helperEP.GetEventPlane(collision.qvecTPCposReVec()[0], collision.qvecTPCposImVec()[0], nmode), nmode), trk.pt(),
                        std::cos(static_cast<float>(nmode) * (trk.phi() -
                                                              helperEP.GetEventPlane(collision.qvecFT0CReVec()[0],
                                                                                     collision.qvecFT0CImVec()[0],
                                                                                     nmode))));
        }
        if (collision.centFT0C() > kCent[4] && collision.centFT0C() <= kCent[5] && q2 > kSeparator[8] && q2 <= kSeparator[9]) {
          histosQA.fill(HIST("histCosDetV2_4050Mid"), helperEP.GetResolution(helperEP.GetEventPlane(collision.qvecFT0CReVec()[0], collision.qvecFT0CImVec()[0], nmode), helperEP.GetEventPlane(collision.qvecTPCposReVec()[0], collision.qvecTPCposImVec()[0], nmode), nmode), trk.pt(),
                        std::cos(static_cast<float>(nmode) * (trk.phi() -
                                                              helperEP.GetEventPlane(collision.qvecFT0CReVec()[0],
                                                                                     collision.qvecFT0CImVec()[0],
                                                                                     nmode))));
        }
        if (collision.centFT0C() > kCent[4] && collision.centFT0C() <= kCent[5] && q2 > kSeparator[9]) {
          histosQA.fill(HIST("histCosDetV2_4050Right"), helperEP.GetResolution(helperEP.GetEventPlane(collision.qvecFT0CReVec()[0], collision.qvecFT0CImVec()[0], nmode), helperEP.GetEventPlane(collision.qvecTPCposReVec()[0], collision.qvecTPCposImVec()[0], nmode), nmode), trk.pt(),
                        std::cos(static_cast<float>(nmode) * (trk.phi() -
                                                              helperEP.GetEventPlane(collision.qvecFT0CReVec()[0],
                                                                                     collision.qvecFT0CImVec()[0],
                                                                                     nmode))));
        }
        if (collision.centFT0C() > kCent[5] && collision.centFT0C() <= kCent[6] && q2 < kSeparator[10]) {
          histosQA.fill(HIST("histCosDetV2_5060Left"), helperEP.GetResolution(helperEP.GetEventPlane(collision.qvecFT0CReVec()[0], collision.qvecFT0CImVec()[0], nmode), helperEP.GetEventPlane(collision.qvecTPCposReVec()[0], collision.qvecTPCposImVec()[0], nmode), nmode), trk.pt(),
                        std::cos(static_cast<float>(nmode) * (trk.phi() -
                                                              helperEP.GetEventPlane(collision.qvecFT0CReVec()[0],
                                                                                     collision.qvecFT0CImVec()[0],
                                                                                     nmode))));
        }
        if (collision.centFT0C() > kCent[5] && collision.centFT0C() <= kCent[6] && q2 > kSeparator[10] && q2 <= kSeparator[11]) {
          histosQA.fill(HIST("histCosDetV2_5060Mid"), helperEP.GetResolution(helperEP.GetEventPlane(collision.qvecFT0CReVec()[0], collision.qvecFT0CImVec()[0], nmode), helperEP.GetEventPlane(collision.qvecTPCposReVec()[0], collision.qvecTPCposImVec()[0], nmode), nmode), trk.pt(),
                        std::cos(static_cast<float>(nmode) * (trk.phi() -
                                                              helperEP.GetEventPlane(collision.qvecFT0CReVec()[0],
                                                                                     collision.qvecFT0CImVec()[0],
                                                                                     nmode))));
        }
        if (collision.centFT0C() > kCent[5] && collision.centFT0C() <= kCent[6] && q2 > kSeparator[11]) {
          histosQA.fill(HIST("histCosDetV2_5060Right"), helperEP.GetResolution(helperEP.GetEventPlane(collision.qvecFT0CReVec()[0], collision.qvecFT0CImVec()[0], nmode), helperEP.GetEventPlane(collision.qvecTPCposReVec()[0], collision.qvecTPCposImVec()[0], nmode), nmode), trk.pt(),
                        std::cos(static_cast<float>(nmode) * (trk.phi() -
                                                              helperEP.GetEventPlane(collision.qvecFT0CReVec()[0],
                                                                                     collision.qvecFT0CImVec()[0],
                                                                                     nmode))));
        }
        if (collision.centFT0C() > kCent[6] && collision.centFT0C() <= kCent[7] && q2 < kSeparator[12]) {
          histosQA.fill(HIST("histCosDetV2_6070Left"), helperEP.GetResolution(helperEP.GetEventPlane(collision.qvecFT0CReVec()[0], collision.qvecFT0CImVec()[0], nmode), helperEP.GetEventPlane(collision.qvecTPCposReVec()[0], collision.qvecTPCposImVec()[0], nmode), nmode), trk.pt(),
                        std::cos(static_cast<float>(nmode) * (trk.phi() -
                                                              helperEP.GetEventPlane(collision.qvecFT0CReVec()[0],
                                                                                     collision.qvecFT0CImVec()[0],
                                                                                     nmode))));
        }
        if (collision.centFT0C() > kCent[6] && collision.centFT0C() <= kCent[7] && q2 > kSeparator[12] && q2 <= kSeparator[13]) {
          histosQA.fill(HIST("histCosDetV2_6070Mid"), helperEP.GetResolution(helperEP.GetEventPlane(collision.qvecFT0CReVec()[0], collision.qvecFT0CImVec()[0], nmode), helperEP.GetEventPlane(collision.qvecTPCposReVec()[0], collision.qvecTPCposImVec()[0], nmode), nmode), trk.pt(),
                        std::cos(static_cast<float>(nmode) * (trk.phi() -
                                                              helperEP.GetEventPlane(collision.qvecFT0CReVec()[0],
                                                                                     collision.qvecFT0CImVec()[0],
                                                                                     nmode))));
        }
        if (collision.centFT0C() > kCent[6] && collision.centFT0C() <= kCent[7] && q2 > kSeparator[13]) {
          histosQA.fill(HIST("histCosDetV2_6070Right"), helperEP.GetResolution(helperEP.GetEventPlane(collision.qvecFT0CReVec()[0], collision.qvecFT0CImVec()[0], nmode), helperEP.GetEventPlane(collision.qvecTPCposReVec()[0], collision.qvecTPCposImVec()[0], nmode), nmode), trk.pt(),
                        std::cos(static_cast<float>(nmode) * (trk.phi() -
                                                              helperEP.GetEventPlane(collision.qvecFT0CReVec()[0],
                                                                                     collision.qvecFT0CImVec()[0],
                                                                                     nmode))));
        }
        if (collision.centFT0C() > kCent[7] && collision.centFT0C() <= kCent[8] && q2 < kSeparator[14]) {
          histosQA.fill(HIST("histCosDetV2_7080Left"), helperEP.GetResolution(helperEP.GetEventPlane(collision.qvecFT0CReVec()[0], collision.qvecFT0CImVec()[0], nmode), helperEP.GetEventPlane(collision.qvecTPCposReVec()[0], collision.qvecTPCposImVec()[0], nmode), nmode), trk.pt(),
                        std::cos(static_cast<float>(nmode) * (trk.phi() -
                                                              helperEP.GetEventPlane(collision.qvecFT0CReVec()[0],
                                                                                     collision.qvecFT0CImVec()[0],
                                                                                     nmode))));
        }
        if (collision.centFT0C() > kCent[7] && collision.centFT0C() <= kCent[8] && q2 > kSeparator[14] && q2 <= kSeparator[15]) {
          histosQA.fill(HIST("histCosDetV2_7080Mid"), helperEP.GetResolution(helperEP.GetEventPlane(collision.qvecFT0CReVec()[0], collision.qvecFT0CImVec()[0], nmode), helperEP.GetEventPlane(collision.qvecTPCposReVec()[0], collision.qvecTPCposImVec()[0], nmode), nmode), trk.pt(),
                        std::cos(static_cast<float>(nmode) * (trk.phi() -
                                                              helperEP.GetEventPlane(collision.qvecFT0CReVec()[0],
                                                                                     collision.qvecFT0CImVec()[0],
                                                                                     nmode))));
        }
        if (collision.centFT0C() > kCent[7] && collision.centFT0C() <= kCent[8] && q2 > kSeparator[15]) {
          histosQA.fill(HIST("histCosDetV2_7080Right"), helperEP.GetResolution(helperEP.GetEventPlane(collision.qvecFT0CReVec()[0], collision.qvecFT0CImVec()[0], nmode), helperEP.GetEventPlane(collision.qvecTPCposReVec()[0], collision.qvecTPCposImVec()[0], nmode), nmode), trk.pt(),
                        std::cos(static_cast<float>(nmode) * (trk.phi() -
                                                              helperEP.GetEventPlane(collision.qvecFT0CReVec()[0],
                                                                                     collision.qvecFT0CImVec()[0],
                                                                                     nmode))));
        }
        if (collision.centFT0C() > kCent[8] && collision.centFT0C() <= kCent[9] && q2 < kSeparator[16]) {
          histosQA.fill(HIST("histCosDetV2_8090Left"), helperEP.GetResolution(helperEP.GetEventPlane(collision.qvecFT0CReVec()[0], collision.qvecFT0CImVec()[0], nmode), helperEP.GetEventPlane(collision.qvecTPCposReVec()[0], collision.qvecTPCposImVec()[0], nmode), nmode), trk.pt(),
                        std::cos(static_cast<float>(nmode) * (trk.phi() -
                                                              helperEP.GetEventPlane(collision.qvecFT0CReVec()[0],
                                                                                     collision.qvecFT0CImVec()[0],
                                                                                     nmode))));
        }
        if (collision.centFT0C() > kCent[8] && collision.centFT0C() <= kCent[9] && q2 > kSeparator[16] && q2 <= kSeparator[17]) {
          histosQA.fill(HIST("histCosDetV2_8090Mid"), helperEP.GetResolution(helperEP.GetEventPlane(collision.qvecFT0CReVec()[0], collision.qvecFT0CImVec()[0], nmode), helperEP.GetEventPlane(collision.qvecTPCposReVec()[0], collision.qvecTPCposImVec()[0], nmode), nmode), trk.pt(),
                        std::cos(static_cast<float>(nmode) * (trk.phi() -
                                                              helperEP.GetEventPlane(collision.qvecFT0CReVec()[0],
                                                                                     collision.qvecFT0CImVec()[0],
                                                                                     nmode))));
        }
        if (collision.centFT0C() > kCent[8] && collision.centFT0C() <= kCent[9] && q2 > kSeparator[17]) {
          histosQA.fill(HIST("histCosDetV2_8090Right"), helperEP.GetResolution(helperEP.GetEventPlane(collision.qvecFT0CReVec()[0], collision.qvecFT0CImVec()[0], nmode), helperEP.GetEventPlane(collision.qvecTPCposReVec()[0], collision.qvecTPCposImVec()[0], nmode), nmode), trk.pt(),
                        std::cos(static_cast<float>(nmode) * (trk.phi() -
                                                              helperEP.GetEventPlane(collision.qvecFT0CReVec()[0],
                                                                                     collision.qvecFT0CImVec()[0],
                                                                                     nmode))));
        }
        if (collision.centFT0C() > kCent[9] && collision.centFT0C() <= kCent[10] && q2 < kSeparator[18]) {
          histosQA.fill(HIST("histCosDetV2_9010Left"), helperEP.GetResolution(helperEP.GetEventPlane(collision.qvecFT0CReVec()[0], collision.qvecFT0CImVec()[0], nmode), helperEP.GetEventPlane(collision.qvecTPCposReVec()[0], collision.qvecTPCposImVec()[0], nmode), nmode), trk.pt(),
                        std::cos(static_cast<float>(nmode) * (trk.phi() -
                                                              helperEP.GetEventPlane(collision.qvecFT0CReVec()[0],
                                                                                     collision.qvecFT0CImVec()[0],
                                                                                     nmode))));
        }
        if (collision.centFT0C() > kCent[9] && collision.centFT0C() <= kCent[10] && q2 > kSeparator[18] && q2 <= kSeparator[19]) {
          histosQA.fill(HIST("histCosDetV2_9010Mid"), helperEP.GetResolution(helperEP.GetEventPlane(collision.qvecFT0CReVec()[0], collision.qvecFT0CImVec()[0], nmode), helperEP.GetEventPlane(collision.qvecTPCposReVec()[0], collision.qvecTPCposImVec()[0], nmode), nmode), trk.pt(),
                        std::cos(static_cast<float>(nmode) * (trk.phi() -
                                                              helperEP.GetEventPlane(collision.qvecFT0CReVec()[0],
                                                                                     collision.qvecFT0CImVec()[0],
                                                                                     nmode))));
        }
        if (collision.centFT0C() > kCent[9] && collision.centFT0C() <= kCent[10] && q2 > kSeparator[19]) {
          histosQA.fill(HIST("histCosDetV2_9010Right"), helperEP.GetResolution(helperEP.GetEventPlane(collision.qvecFT0CReVec()[0], collision.qvecFT0CImVec()[0], nmode), helperEP.GetEventPlane(collision.qvecTPCposReVec()[0], collision.qvecTPCposImVec()[0], nmode), nmode), trk.pt(),
                        std::cos(static_cast<float>(nmode) * (trk.phi() -
                                                              helperEP.GetEventPlane(collision.qvecFT0CReVec()[0],
                                                                                     collision.qvecFT0CImVec()[0],
                                                                                     nmode))));
        }
      }
    }
  }

  void process(MyCollisions::iterator const& collision, MyTracks const& tracks)
  {
    if (!selectEvent(collision)) {
      return;
    }
    for (std::size_t i = 0; i < cfgNmods->size(); i++) {
      fillHistosQvec(collision, cfgNmods->at(i));
      fillHistosFlow(collision, tracks, cfgNmods->at(i));
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FlowEseTask>(cfgc)};
}
