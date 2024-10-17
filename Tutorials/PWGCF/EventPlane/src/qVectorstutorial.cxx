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

// C++/ROOT includes.
#include <chrono>
#include <string>
#include <vector>
#include <TComplex.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TMath.h>
#include <TVector2.h>

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
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Centrality.h"
#include "Common/Core/EventPlaneHelper.h"
#include "Common/Core/TrackSelection.h"

#include "CommonConstants/PhysicsConstants.h"

// o2 includes.

using namespace o2;
using namespace o2::framework;

using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::QvectorFT0CVecs, aod::QvectorTPCposVecs, aod::QvectorTPCnegVecs>;
using MyTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TrackSelectionExtension>;

struct qVectorstutorial {
  HistogramRegistry histosQA{"histosQA", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  Configurable<std::vector<int>> cfgnMods{"cfgnMods", {2}, "Modulation of interest"};
  Configurable<std::string> cfgDetName{"cfgDetName", "FT0C", "The name of detector to be analyzed"};
  Configurable<std::string> cfgRefAName{"cfgRefAName", "TPCpos", "The name of detector for reference A"};
  Configurable<std::string> cfgRefBName{"cfgRefBName", "TPCneg", "The name of detector for reference B"};

  Configurable<float> cfgMinPt{"cfgMinPt", 0.15, "Minimum transverse momentum for charged track"};
  Configurable<float> cfgMaxEta{"cfgMaxEta", 0.8, "Maximum pseudorapidiy for charged track"};
  Configurable<float> cfgMaxDCArToPVcut{"cfgMaxDCArToPVcut", 0.1, "Maximum transverse DCA"};
  Configurable<float> cfgMaxDCAzToPVcut{"cfgMaxDCAzToPVcut", 1.0, "Maximum longitudinal DCA"};

  ConfigurableAxis cfgaxisQvecF{"cfgaxisQvecF", {300, -1, 1}, ""};
  ConfigurableAxis cfgaxisQvec{"cfgaxisQvec", {100, -3, 3}, ""};
  ConfigurableAxis cfgaxisCent{"cfgaxisCent", {100, 0, 100}, ""};

  ConfigurableAxis cfgaxiscos{"cfgaxiscos", {102, -1.02, 1.02}, ""};
  ConfigurableAxis cfgaxispt{"cfgaxispt", {100, 0, 10}, ""};
  ConfigurableAxis cfgaxisCentMerged{"cfgaxisCentMerged", {20, 0, 100}, ""};

  EventPlaneHelper helperEP;

  void init(InitContext const&)
  {
    AxisSpec axisCent{cfgaxisCent, "centrality"};
    AxisSpec axisQvec{cfgaxisQvec, "Q"};
    AxisSpec axisQvecF{cfgaxisQvecF, "Q"};
    AxisSpec axisEvtPl = {100, -1.0 * constants::math::PI, constants::math::PI};

    AxisSpec axisCos{cfgaxiscos, "angle function"};
    AxisSpec axisPt{cfgaxispt, "trasverse momentum"};
    AxisSpec axisCentMerged{cfgaxisCentMerged, "merged centrality"};

    histosQA.add(Form("histQvecV2"), "", {HistType::kTH3F, {axisQvecF, axisQvecF, axisCent}});
    histosQA.add(Form("histEvtPlV2"), "", {HistType::kTH2F, {axisEvtPl, axisCent}});
    histosQA.add(Form("histQvecRes_SigRefAV2"), "", {HistType::kTH2F, {axisQvecF, axisCent}});
    histosQA.add(Form("histCosDetV2"), "", {HistType::kTH3F, {axisCentMerged, axisPt, axisCos}});
  }

  template <typename CollType>
  bool SelEvent(const CollType& collision)
  {
    if (!collision.sel8()) {
      return 0;
    }
    if (!collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return 0;
    }
    if (!collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      return 0;
    }
    if (!collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      return 0;
    }

    return 1;
  }

  template <typename TrackType>
  bool SelTrack(const TrackType track)
  {
    if (track.pt() < cfgMinPt)
      return false;
    if (std::abs(track.eta()) > cfgMaxEta)
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

  template <typename CollType>
  void fillHistosQvec(const CollType& collision, int nmode)
  {
    if (nmode == 2) {
      histosQA.fill(HIST("histQvecV2"), collision.qvecFT0CReVec()[0], collision.qvecFT0CImVec()[0], collision.centFT0C());
      histosQA.fill(HIST("histEvtPlV2"), helperEP.GetEventPlane(collision.qvecFT0CReVec()[0], collision.qvecFT0CImVec()[0], nmode), collision.centFT0C());
      histosQA.fill(HIST("histQvecRes_SigRefAV2"), helperEP.GetResolution(helperEP.GetEventPlane(collision.qvecFT0CReVec()[0], collision.qvecFT0CImVec()[0], nmode), helperEP.GetEventPlane(collision.qvecTPCposReVec()[0], collision.qvecTPCposImVec()[0], nmode), nmode), collision.centFT0C());
    }
  }

  template <typename CollType, typename TrackType>
  void fillHistosFlow(const CollType& collision, const TrackType& track, int nmode)
  {
    if (collision.sumAmplFT0C() < 1e-4) {
      return;
    }
    for (auto& trk : track) {
      if (!SelTrack(trk)) {
        continue;
      }
      if (nmode == 2) {
        histosQA.fill(HIST("histCosDetV2"), collision.centFT0C(), trk.pt(),
                      std::cos(static_cast<float>(nmode) * (trk.phi() - helperEP.GetEventPlane(collision.qvecFT0CReVec()[0], collision.qvecFT0CImVec()[0], nmode))));
      }
    }
  }

  void process(MyCollisions::iterator const& collision, MyTracks const& tracks)
  {
    if (!SelEvent(collision)) {
      return;
    }
    for (auto i = 0; i < cfgnMods->size(); i++) {
      fillHistosQvec(collision, cfgnMods->at(i));
      fillHistosFlow(collision, tracks, cfgnMods->at(i));
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<qVectorstutorial>(cfgc)};
}
