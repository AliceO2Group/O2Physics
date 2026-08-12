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

/// \file kaDaughterProducer.cxx
/// \brief Table producer for shared kaon/pion daughter derived data (K*(892), phi)
///
/// \author sarjeeta.gami@cern.ch

#include "PWGLF/DataModel/EPCalibrationTables.h"
#include "PWGLF/DataModel/LFKaonDaughterTables.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/RCTSelectionFlags.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CommonConstants/MathConstants.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::rctsel;

struct KaDaughterProducer {
  Produces<aod::KaDaughterEvents> kaEvent;
  Produces<aod::KaDaughterTracks> kaTrack;

  struct RCTCut : ConfigurableGroup {
    Configurable<bool> requireRCTFlagChecker{"requireRCTFlagChecker", true, "Check event quality in run condition table"};
    Configurable<std::string> cfgEvtRCTFlagCheckerLabel{"cfgEvtRCTFlagCheckerLabel", "CBT_hadronPID", "Evt sel: RCT flag checker label"};
    Configurable<bool> cfgEvtRCTFlagCheckerZDCCheck{"cfgEvtRCTFlagCheckerZDCCheck", false, "Evt sel: RCT flag checker ZDC check"};
    Configurable<bool> cfgEvtRCTFlagCheckerLimitAcceptAsBad{"cfgEvtRCTFlagCheckerLimitAcceptAsBad", true, "Evt sel: RCT flag checker treat Limited Acceptance As Bad"};
    RCTFlagsChecker rctChecker;
  };
  RCTCut rctCut;

  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutCentrality{"cfgCutCentrality", 80.0f, "Accepted maximum Centrality"};
  Configurable<float> cfgCutPT{"cfgCutPT", 0.2, "PT cut on daughter track"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8, "Eta cut on daughter track"};
  Configurable<float> cfgCutDCAxy{"cfgCutDCAxy", 2.0f, "DCAxy range for tracks"};
  Configurable<float> cfgCutDCAz{"cfgCutDCAz", 2.0f, "DCAz range for tracks"};
  Configurable<bool> useGlobalTrack{"useGlobalTrack", true, "use Global track"};
  Configurable<int> cfgITScluster{"cfgITScluster", 0, "Number of ITS cluster"};
  Configurable<int> cfgTPCcluster{"cfgTPCcluster", 70, "Number of TPC cluster"};
  Configurable<bool> additionalEvSel1{"additionalEvSel1", true, "Additional evsel1"};
  Configurable<bool> additionalEvSel2{"additionalEvSel2", true, "Additional evsel2"};
  Configurable<bool> additionalEvSel3{"additionalEvSel3", true, "Additional evsel3"};
  Configurable<bool> additionalEvSel4{"additionalEvSel4", true, "Additional evsel4"};
  Configurable<bool> fillOccupancy{"fillOccupancy", false, "fill Occupancy"};
  Configurable<int> cfgOccupancyCut{"cfgOccupancyCut", 500, "Occupancy cut"};
  Configurable<bool> additionalQAplots1{"additionalQAplots1", true, "Additional QA plots (event-plane resolution etc.)"};

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter centralityFilter = nabs(aod::cent::centFT0C) < cfgCutCentrality;
  Filter acceptanceFilter = (nabs(aod::track::eta) < cfgCutEta && nabs(aod::track::pt) > cfgCutPT);
  Filter dcacutFilter = (nabs(aod::track::dcaXY) < cfgCutDCAxy) && (nabs(aod::track::dcaZ) < cfgCutDCAz);

  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs,
                                                  aod::EPCalibrationTables, aod::Mults>>;
  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA,
                                                  aod::TrackSelection, aod::pidTOFbeta,
                                                  aod::pidTPCFullKa, aod::pidTOFFullKa,
                                                  aod::pidTPCFullPi, aod::pidTOFFullPi>>;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
    rctCut.rctChecker.init(rctCut.cfgEvtRCTFlagCheckerLabel, rctCut.cfgEvtRCTFlagCheckerZDCCheck, rctCut.cfgEvtRCTFlagCheckerLimitAcceptAsBad);

    histos.add("hEvtSelInfo", "hEvtSelInfo", kTH1F, {{10, 0, 10.0}});
    if (additionalQAplots1) {
      std::vector<double> occupancyBinning = {0.0, 500.0, 1000.0, 1500.0, 2000.0, 3000.0, 4000.0, 5000.0, 50000.0};
      AxisSpec phiAxis = {500, -6.28, 6.28, "phi"};
      AxisSpec resAxis = {6000, -30, 30, "Res"};
      AxisSpec centAxis = {8, 0, 80, "V0M (%)"};
      AxisSpec occupancyAxis = {occupancyBinning, "Occupancy"};
      histos.add("hFTOCvsTPCSelected", "Mult correlation FT0C vs. TPC after selection", kTH2F, {{80, 0.0f, 80.0f}, {100, -0.5f, 5999.5f}});
      histos.add("hCentrality", "Centrality distribution", kTH1F, {{200, 0.0, 200.0}});
      histos.add("hOccupancy", "Occupancy distribution", kTH1F, {occupancyAxis});
      histos.add("hVtxZ", "Vertex distribution in Z;Z (cm)", kTH1F, {{400, -20.0, 20.0}});
      histos.add("hPsiFT0C", "PsiFT0C", kTH2F, {centAxis, phiAxis});
      histos.add("hPsiFT0A", "PsiFT0A", kTH2F, {centAxis, phiAxis});
      histos.add("hPsiTPC", "PsiTPC", kTH2F, {centAxis, phiAxis});
      histos.add("ResFT0CTPC", "ResFT0CTPC", kTH2F, {centAxis, resAxis});
      histos.add("ResFT0CFT0A", "ResFT0CFT0A", kTH2F, {centAxis, resAxis});
      histos.add("ResFT0ATPC", "ResFT0ATPC", kTH2F, {centAxis, resAxis});
      histos.add("ResFT0CTPCSP", "ResFT0CTPCSP", kTH2F, {centAxis, resAxis});
      histos.add("ResFT0CFT0ASP", "ResFT0CFT0ASP", kTH2F, {centAxis, resAxis});
      histos.add("ResFT0ATPCSP", "ResFT0ATPCSP", kTH2F, {centAxis, resAxis});
      histos.add("ResTrackSPFT0CTPC", "ResTrackSPFT0CTPC", kTH3F, {centAxis, occupancyAxis, resAxis});
      histos.add("ResTrackSPFT0CFT0A", "ResTrackSPFT0CFT0A", kTH3F, {centAxis, occupancyAxis, resAxis});
      histos.add("ResTrackSPFT0ATPC", "ResTrackSPFT0ATPC", kTH3F, {centAxis, occupancyAxis, resAxis});
    }
  }

  template <typename T>
  bool selectionTrack(const T& candidate)
  {
    if (useGlobalTrack && !(candidate.isGlobalTrack() && candidate.isPVContributor() && candidate.itsNCls() > cfgITScluster && candidate.tpcNClsFound() > cfgTPCcluster)) {
      return false;
    }
    if (!useGlobalTrack && !(candidate.isPVContributor() && candidate.itsNCls() > cfgITScluster)) {
      return false;
    }
    return true;
  }

  void processData(EventCandidates::iterator const& collision, TrackCandidates const& tracks, aod::BCsWithTimestamps const&)
  {
    histos.fill(HIST("hEvtSelInfo"), 0.5);
    if (rctCut.requireRCTFlagChecker && !rctCut.rctChecker(collision)) return;
    histos.fill(HIST("hEvtSelInfo"), 1.5);
    if (!collision.sel8()) return;
    if (!collision.triggereventep()) return;
    if (additionalEvSel1 && !collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) return;
    if (additionalEvSel2 && !collision.selection_bit(aod::evsel::kNoITSROFrameBorder)) return;
    if (additionalEvSel3 && !collision.selection_bit(aod::evsel::kNoSameBunchPileup)) return;
    if (additionalEvSel4 && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) return;
    histos.fill(HIST("hEvtSelInfo"), 2.5);

    auto centrality = collision.centFT0C();
    int occupancy = collision.trackOccupancyInTimeRange();
    auto psiFT0C = collision.psiFT0C();
    auto qFT0C = collision.qFT0C();

    if (fillOccupancy && occupancy > cfgOccupancyCut) return;
    histos.fill(HIST("hEvtSelInfo"), 3.5);

    if (additionalQAplots1) {
      auto multTPC = collision.multNTracksPV();
      auto psiFT0A = collision.psiFT0A();
      auto psiTPC = collision.psiTPC();
      auto qFT0A = collision.qFT0A();
      auto qTPC = collision.qTPC();
      histos.fill(HIST("hFTOCvsTPCSelected"), centrality, multTPC);
      histos.fill(HIST("hPsiFT0C"), centrality, psiFT0C);
      histos.fill(HIST("hPsiFT0A"), centrality, psiFT0A);
      histos.fill(HIST("hPsiTPC"), centrality, psiTPC);
      histos.fill(HIST("ResFT0CTPC"), centrality, std::cos(2.0 * (psiFT0C - psiTPC)));
      histos.fill(HIST("ResFT0CFT0A"), centrality, std::cos(2.0 * (psiFT0C - psiFT0A)));
      histos.fill(HIST("ResFT0ATPC"), centrality, std::cos(2.0 * (psiTPC - psiFT0A)));
      histos.fill(HIST("ResFT0CTPCSP"), centrality, qFT0C * qTPC * std::cos(2.0 * (psiFT0C - psiTPC)));
      histos.fill(HIST("ResFT0CFT0ASP"), centrality, qFT0C * qFT0A * std::cos(2.0 * (psiFT0C - psiFT0A)));
      histos.fill(HIST("ResFT0ATPCSP"), centrality, qTPC * qFT0A * std::cos(2.0 * (psiTPC - psiFT0A)));
      histos.fill(HIST("hCentrality"), centrality);
      histos.fill(HIST("hOccupancy"), occupancy);
      histos.fill(HIST("hVtxZ"), collision.posZ());
      // filled once per event now, instead of once per track-pair as in the original processSE
      histos.fill(HIST("ResTrackSPFT0CTPC"), centrality, occupancy, qFT0C * qTPC * std::cos(2.0 * (psiFT0C - psiTPC)));
      histos.fill(HIST("ResTrackSPFT0CFT0A"), centrality, occupancy, qFT0C * qFT0A * std::cos(2.0 * (psiFT0C - psiFT0A)));
      histos.fill(HIST("ResTrackSPFT0ATPC"), centrality, occupancy, qTPC * qFT0A * std::cos(2.0 * (psiTPC - psiFT0A)));
    }

    struct Stored {
      float px, py, pz, dcaXY, dcaZ, tpcKa, tofKa, tpcPi, tofPi, beta;
      int8_t sign;
      int64_t id;
      bool hasTOF;
    };
    std::vector<Stored> sel;
    for (const auto& t : tracks) {
      if (!selectionTrack(t)) continue;
      sel.push_back({t.px(), t.py(), t.pz(), t.dcaXY(), t.dcaZ(),
                     t.tpcNSigmaKa(), t.tofNSigmaKa(), t.tpcNSigmaPi(), t.tofNSigmaPi(),
                     t.hasTOF() ? t.beta() : -999.f, static_cast<int8_t>(t.sign()), t.globalIndex(), t.hasTOF()});
    }
    if (sel.empty()) return;

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    kaEvent(centrality, collision.posZ(), occupancy, psiFT0C, qFT0C, /*psiZDCC=*/0.f,
           collision.bcId(), bc.runNumber(), bc.timestamp());
    const int64_t idx = kaEvent.lastIndex();
    for (const auto& s : sel) {
      kaTrack(idx, s.px, s.py, s.pz, s.sign, s.id, s.dcaXY, s.dcaZ,
             s.tpcKa, s.tofKa, s.tpcPi, s.tofPi, s.hasTOF, s.beta);
    }
  }
  PROCESS_SWITCH(KaDaughterProducer, processData, "Produce shared kaon/pion daughter candidates", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<KaDaughterProducer>(cfgc)};
}
