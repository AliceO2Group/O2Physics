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
/// \file ptmultCorr.cxx
///
/// \brief task for analysis of charged-particle RAA at midrapidity in light-ion collisions
/// \author Abhi Modak (abhi.modak@cern.ch)
/// \since October 01, 2025

#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGMM/Mult/DataModel/Index.h"
#include "PWGMM/Mult/DataModel/bestCollisionTable.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/MathConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/Configurable.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/GlobalTrackID.h"
#include "ReconstructionDataFormats/Track.h"

#include <TPDGCode.h>

#include <cmath>
#include <cstdlib>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::track;
using namespace o2::aod::evsel;

using ColDataTablePbPb = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::PVMults, aod::CentFT0Cs>;
using ColDataTablepp = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::PVMults>;
using ColMCRecTablePbPb = soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions, aod::EvSels, aod::Mults, aod::PVMults, aod::CentFT0Cs>>;
using ColMCRecTablepp = soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions, aod::EvSels, aod::Mults, aod::PVMults>>;
using ColMCTrueTable = soa::Join<aod::McCollisions, aod::MultMCExtras, aod::McCollsExtra>;
using TrackDataTable = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>;
using FilTrackDataTable = soa::Filtered<TrackDataTable>;
using TrackMCRecTable = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::McTrackLabels, aod::TrackSelection>;
using FilTrackMCRecTable = soa::Filtered<TrackMCRecTable>;
using TrackMCTrueTable = aod::McParticles;

enum {
  kTrackTypebegin = 0,
  kGlobalplusITS = 1,
  kGlobalonly,
  kITSonly,
  kTrackTypeend
};

enum {
  kGenpTbegin = 0,
  kNoGenpTVar = 1,
  kGenpTup,
  kGenpTdown,
  kGenpTend
};

enum {
  kGenTrkTypebegin = 0,
  kGenAll = 1,
  kGenPion,
  kGenKaon,
  kGenProton,
  kGenOther,
  kGenTrkTypeend
};

enum {
  kRecTrkTypebegin = 0,
  kRecoAll = 1,
  kRecoHasmc,
  kRecoPrimary,
  kRecoPion,
  kRecoKaon,
  kRecoProton,
  kRecoOther,
  kRecoSecondary,
  kRecoWeakDecay,
  kRecoFake,
  kRecoBkg,
  kRecTrkTypeend
};

static constexpr TrackSelectionFlags::flagtype TrackSelectionIts =
  TrackSelectionFlags::kITSNCls | TrackSelectionFlags::kITSChi2NDF |
  TrackSelectionFlags::kITSHits;
static constexpr TrackSelectionFlags::flagtype TrackSelectionTpc =
  TrackSelectionFlags::kTPCNCls |
  TrackSelectionFlags::kTPCCrossedRowsOverNCls |
  TrackSelectionFlags::kTPCChi2NDF;
static constexpr TrackSelectionFlags::flagtype TrackSelectionDca =
  TrackSelectionFlags::kDCAz | TrackSelectionFlags::kDCAxy;
static constexpr TrackSelectionFlags::flagtype TrackSelectionDcaxyOnly =
  TrackSelectionFlags::kDCAxy;

AxisSpec axisEvent{15, 0.5, 15.5, "#Event", "EventAxis"};
AxisSpec axisVtxZ{40, -20, 20, "Vertex Z", "VzAxis"};
AxisSpec axisEta{40, -2, 2, "#eta", "EtaAxis"};
AxisSpec axisPhi{{0, o2::constants::math::PIQuarter, o2::constants::math::PIHalf, o2::constants::math::PIQuarter * 3., o2::constants::math::PI, o2::constants::math::PIQuarter * 5., o2::constants::math::PIHalf * 3., o2::constants::math::PIQuarter * 7., o2::constants::math::TwoPI}, "#phi", "PhiAxis"};
AxisSpec axisPhi2{629, 0, o2::constants::math::TwoPI, "#phi"};
AxisSpec axisCent{100, 0, 100, "#Cent"};
AxisSpec axisDeltaPt{50, -1.0, +1.0, "#Delta(pT)"};
AxisSpec axisGenPtVary = {kGenpTend - 1, +kGenpTbegin + 0.5, +kGenpTend - 0.5, "", "GenpTVaryAxis"};
AxisSpec axisGenTrkType = {kGenTrkTypeend - 1, +kGenTrkTypebegin + 0.5, +kGenTrkTypeend - 0.5, "", "GenTrackTypeAxis"};
AxisSpec axisRecTrkType = {kRecTrkTypeend - 1, +kRecTrkTypebegin + 0.5, +kRecTrkTypeend - 0.5, "", "RecTrackTypeAxis"};
AxisSpec axisTrackType = {kTrackTypeend - 1, +kTrackTypebegin + 0.5, +kTrackTypeend - 0.5, "", "TrackTypeAxis"};
auto static constexpr KminCharge = 3.f;
auto static constexpr KminPtCut = 0.1f;

struct PtmultCorr {

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Service<o2::framework::O2DatabasePDG> pdg;
  Preslice<TrackMCRecTable> perCollision = aod::track::collisionId;
  Configurable<float> etaRange{"etaRange", 1.0f, "Eta range to consider"};
  Configurable<float> vtxRange{"vtxRange", 10.0f, "Vertex Z range to consider"};
  Configurable<float> occuRange{"occuRange", 500.0f, "Occupancy range to consider"};
  Configurable<float> dcaZ{"dcaZ", 0.2f, "Custom DCA Z cut (ignored if negative)"};
  Configurable<float> cfgPtCutMin{"cfgPtCutMin", 0.15f, "minimum accepted track pT"};
  Configurable<float> extraphicut1{"extraphicut1", 3.07666f, "Extra Phi cut 1"};
  Configurable<float> extraphicut2{"extraphicut2", 3.12661f, "Extra Phi cut 2"};
  Configurable<float> extraphicut3{"extraphicut3", 0.03f, "Extra Phi cut 3"};
  Configurable<float> extraphicut4{"extraphicut4", 6.253f, "Extra Phi cut 4"};
  ConfigurableAxis ptHistBin{"ptHistBin", {200, 0., 20.}, ""};
  ConfigurableAxis centralityBinning{"centralityBinning", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100}, ""};
  ConfigurableAxis binsImpactPar{"binsImpactPar", {VARIABLE_WIDTH, 0.0, 3.00065, 4.28798, 6.14552, 7.6196, 8.90942, 10.0897, 11.2002, 12.2709, 13.3167, 14.4173, 23.2518}, "Binning of the impact parameter axis"};

  Configurable<bool> isApplySameBunchPileup{"isApplySameBunchPileup", false, "Enable SameBunchPileup cut"};
  Configurable<bool> isApplyGoodZvtxFT0vsPV{"isApplyGoodZvtxFT0vsPV", false, "Enable GoodZvtxFT0vsPV cut"};
  Configurable<bool> isApplyExtraPhiCut{"isApplyExtraPhiCut", false, "Enable extra phi cut"};
  Configurable<bool> isApplyNoCollInTimeRangeStandard{"isApplyNoCollInTimeRangeStandard", true, "Enable NoCollInTimeRangeStandard cut"};
  Configurable<bool> isApplyNoCollInRofStandard{"isApplyNoCollInRofStandard", false, "Enable NoCollInRofStandard cut"};
  Configurable<bool> isApplyNoHighMultCollInPrevRof{"isApplyNoHighMultCollInPrevRof", false, "Enable NoHighMultCollInPrevRof cut"};
  Configurable<bool> isApplyFT0CbasedOccupancy{"isApplyFT0CbasedOccupancy", false, "Enable FT0CbasedOccupancy cut"};
  Configurable<bool> isApplyInelgt0{"isApplyInelgt0", false, "Enable INEL > 0 condition"};
  Configurable<bool> isApplyOccuCut{"isApplyOccuCut", false, "Enable occupancy selection"};

  void init(InitContext const&)
  {
    AxisSpec centAxis = {centralityBinning, "Centrality", "CentralityAxis"};
    AxisSpec axisPt = {ptHistBin, "pT", "pTAxis"};
    AxisSpec impactParAxis = {binsImpactPar, "Impact Parameter"};

    histos.add("EventHist", "EventHist", kTH1D, {axisEvent}, false);
    histos.add("VtxZHist", "VtxZHist", kTH1D, {axisVtxZ}, false);
    histos.add("CentPercentileHist", "CentPercentileHist", kTH1D, {axisCent}, false);
    histos.add("CentPercentileMCRecHist", "CentPercentileMCRecHist", kTH1D, {axisCent}, false);
    histos.add("PhiVsEtaHistNoCut", "PhiVsEtaHistNoCut", kTH2D, {axisPhi2, axisEta}, false);
    histos.add("PhiVsEtaHistWithCut", "PhiVsEtaHistWithCut", kTH2D, {axisPhi2, axisEta}, false);

    if (doprocessDataPbPb) {
      histos.add("hdatazvtxcent", "hdatazvtxcent", kTH2D, {axisVtxZ, centAxis}, false);
      histos.add("hdatahistPbPb", "hdatahistPbPb", kTHnSparseD, {axisVtxZ, centAxis, axisPt, axisPhi, axisTrackType}, false);
    }

    if (doprocessDatapp) {
      histos.add("hdatahistpp", "hdatahistpp", kTHnSparseD, {axisVtxZ, axisPt, axisPhi, axisTrackType}, false);
    }

    if (doprocessMCeffPbPb) {
      histos.add("hPbPbGenMCvtxz", "hPbPbGenMCvtxz", kTH1D, {axisVtxZ}, false);
      histos.add("hPbPbGenMCvtxzcent", "hPbPbGenMCvtxzcent", kTH2D, {axisVtxZ, centAxis}, false);
      histos.add("hPbPbGenMCAssoRecvtxz", "hPbPbGenMCAssoRecvtxz", kTH1D, {axisVtxZ}, false);
      histos.add("hPbPbGenMCAssoRecvtxzcent", "hPbPbGenMCAssoRecvtxzcent", kTH2D, {axisVtxZ, centAxis}, false);
      histos.add("hPbPbGenMCdndpt", "hPbPbGenMCdndpt", kTHnSparseD, {axisVtxZ, centAxis, axisPt, axisPhi}, false);
      histos.add("hPbPbGenMCAssoRecdndpt", "hPbPbGenMCAssoRecdndpt", kTHnSparseD, {axisVtxZ, centAxis, axisPt, axisPhi, axisGenTrkType, axisGenPtVary}, false);

      histos.add("hPbPbRecMCvtxz", "hPbPbRecMCvtxz", kTH1D, {axisVtxZ}, false);
      histos.add("hPbPbRecMCvtxzcent", "hPbPbRecMCvtxzcent", kTH2D, {axisVtxZ, centAxis}, false);
      histos.add("hPbPbRecMCcent", "hPbPbRecMCcent", kTH1D, {axisCent}, false);
      histos.add("hPbPbRecMCdndpt", "hPbPbRecMCdndpt", kTHnSparseD, {axisVtxZ, centAxis, axisPt, axisPhi, axisRecTrkType}, false);
      histos.add("hPbPbEtaReso", "hPbPbEtaReso", kTH2D, {axisPt, axisDeltaPt});
    }

    if (doprocessMCeffpp) {
      histos.add("hppGenMCvtxz", "hppGenMCvtxz", kTH1D, {axisVtxZ}, false);
      histos.add("hppGenMCAssoRecvtxz", "hppGenMCAssoRecvtxz", kTH1D, {axisVtxZ}, false);
      histos.add("hppGenMCdndpt", "hppGenMCdndpt", kTHnSparseD, {axisVtxZ, axisPt, axisPhi}, false);
      histos.add("hppGenMCAssoRecdndpt", "hppGenMCAssoRecdndpt", kTHnSparseD, {axisVtxZ, axisPt, axisPhi, axisGenTrkType, axisGenPtVary}, false);

      histos.add("hppRecMCvtxz", "hppRecMCvtxz", kTH1D, {axisVtxZ}, false);
      histos.add("hppRecMCdndpt", "hppRecMCdndpt", kTHnSparseD, {axisVtxZ, axisPt, axisPhi, axisRecTrkType}, false);
      histos.add("hppEtaReso", "hppEtaReso", kTH2D, {axisPt, axisDeltaPt});
    }

    if (doprocessEvtLossSigLossMCpp || doprocessEvtLossSigLossMCPbPb) {
      histos.add("MCEventHist", "MCEventHist", kTH1F, {axisEvent}, false);
      auto hstat = histos.get<TH1>(HIST("MCEventHist"));
      auto* x = hstat->GetXaxis();
      x->SetBinLabel(1, "All MC events");
      x->SetBinLabel(2, "MC events with atleast one reco event");
      histos.add("hgenptBeforeEvtSel", "hgenptBeforeEvtSel", kTH1F, {axisPt}, false);
      histos.add("hgenptAfterEvtSel", "hgenptAfterEvtSel", kTH1F, {axisPt}, false);
      histos.add("hgenptBeforeEvtSelPbPb", "hgenptBeforeEvtSelPbPb", kTH2F, {axisPt, impactParAxis}, false);
      histos.add("hgenptAfterEvtSelPbPb", "hgenptAfterEvtSelPbPb", kTH2F, {axisPt, impactParAxis}, false);
      histos.add("hImpactParameterGen", "Impact parameter of generated MC events", kTH1F, {impactParAxis});
      histos.add("hImpactParameterRec", "Impact parameter of selected MC events", kTH1F, {impactParAxis});
      histos.add("hImpactParvsCentrRec", "Impact parameter of selected MC events vs centrality", kTH2F, {axisCent, impactParAxis});
    }

    auto hstat = histos.get<TH1>(HIST("EventHist"));
    auto* x = hstat->GetXaxis();
    x->SetBinLabel(1, "All events");
    x->SetBinLabel(2, "sel8");
    x->SetBinLabel(3, "kNoSameBunchPileup"); // reject collisions in case of pileup with another collision in the same foundBC
    x->SetBinLabel(4, "kIsGoodZvtxFT0vsPV"); // small difference between z-vertex from PV and from FT0
    x->SetBinLabel(5, "ApplyNoCollInTimeRangeStandard");
    x->SetBinLabel(6, "ApplyNoCollInRofStandard");
    x->SetBinLabel(7, "ApplyNoHighMultCollInPrevRof");
    x->SetBinLabel(8, "INEL > 0");
    x->SetBinLabel(9, "|vz|<10");
    x->SetBinLabel(10, "Occupancy<500");
  }

  template <typename CheckCol>
  bool isEventSelected(CheckCol const& col)
  {
    histos.fill(HIST("EventHist"), 1);

    if (!col.sel8()) {
      return false;
    }
    histos.fill(HIST("EventHist"), 2);

    if (isApplySameBunchPileup && !col.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    histos.fill(HIST("EventHist"), 3);

    if (isApplyGoodZvtxFT0vsPV && !col.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    histos.fill(HIST("EventHist"), 4);

    if (isApplyNoCollInTimeRangeStandard && !col.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      return false;
    }
    histos.fill(HIST("EventHist"), 5);

    if (isApplyNoCollInRofStandard && !col.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
      return false;
    }
    histos.fill(HIST("EventHist"), 6);

    if (isApplyNoHighMultCollInPrevRof && !col.selection_bit(o2::aod::evsel::kNoHighMultCollInPrevRof)) {
      return false;
    }
    histos.fill(HIST("EventHist"), 7);

    if (isApplyInelgt0 && !col.isInelGt0()) {
      return false;
    }
    histos.fill(HIST("EventHist"), 8);

    if (std::abs(col.posZ()) >= vtxRange) {
      return false;
    }
    histos.fill(HIST("EventHist"), 9);

    auto occu = isApplyFT0CbasedOccupancy ? col.ft0cOccupancyInTimeRange() : col.trackOccupancyInTimeRange();
    if (isApplyOccuCut && occu > occuRange) {
      return false;
    }
    histos.fill(HIST("EventHist"), 10);
    return true;
  }

  template <typename CheckTrack>
  bool isTrackSelected(CheckTrack const& track)
  {
    if (std::abs(track.eta()) >= etaRange) {
      return false;
    }
    histos.fill(HIST("PhiVsEtaHistNoCut"), track.phi(), track.eta());
    if (isApplyExtraPhiCut && ((track.phi() > extraphicut1 && track.phi() < extraphicut2) || track.phi() <= extraphicut3 || track.phi() >= extraphicut4)) {
      return false;
    }
    histos.fill(HIST("PhiVsEtaHistWithCut"), track.phi(), track.eta());
    return true;
  }

  template <typename CheckGenTrack>
  bool isGenTrackSelected(CheckGenTrack const& track)
  {
    if (!track.isPhysicalPrimary()) {
      return false;
    }
    if (!track.producedByGenerator()) {
      return false;
    }
    auto pdgTrack = pdg->GetParticle(track.pdgCode());
    if (pdgTrack == nullptr) {
      return false;
    }
    if (std::abs(pdgTrack->Charge()) < KminCharge) {
      return false;
    }
    if (std::abs(track.eta()) >= etaRange) {
      return false;
    }
    if (isApplyExtraPhiCut && ((track.phi() > extraphicut1 && track.phi() < extraphicut2) || track.phi() <= extraphicut3 || track.phi() >= extraphicut4)) {
      return false;
    }
    return true;
  }

  Filter fTrackSelectionITS = ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::ITS) &&
                              ncheckbit(aod::track::trackCutFlag, TrackSelectionIts);
  Filter fTrackSelectionTPC = ifnode(ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::TPC),
                                     ncheckbit(aod::track::trackCutFlag, TrackSelectionTpc), true);
  Filter fTrackSelectionDCA = ifnode(dcaZ.node() > 0.f, nabs(aod::track::dcaZ) <= dcaZ && ncheckbit(aod::track::trackCutFlag, TrackSelectionDcaxyOnly),
                                     ncheckbit(aod::track::trackCutFlag, TrackSelectionDca));
  Filter fTracksPt = aod::track::pt > cfgPtCutMin;

  void processDataPbPb(ColDataTablePbPb::iterator const& cols, FilTrackDataTable const& tracks)
  {
    if (!isEventSelected(cols)) {
      return;
    }
    histos.fill(HIST("VtxZHist"), cols.posZ());
    histos.fill(HIST("CentPercentileHist"), cols.centFT0C());
    histos.fill(HIST("hdatazvtxcent"), cols.posZ(), cols.centFT0C());

    for (const auto& track : tracks) {
      if (!isTrackSelected(track)) {
        continue;
      }
      histos.fill(HIST("hdatahistPbPb"), cols.posZ(), cols.centFT0C(), track.pt(), track.phi(), kGlobalplusITS);
      if (track.hasTPC()) {
        histos.fill(HIST("hdatahistPbPb"), cols.posZ(), cols.centFT0C(), track.pt(), track.phi(), kGlobalonly);
      } else {
        histos.fill(HIST("hdatahistPbPb"), cols.posZ(), cols.centFT0C(), track.pt(), track.phi(), kITSonly);
      }
    }
  }

  void processDatapp(ColDataTablepp::iterator const& cols, FilTrackDataTable const& tracks)
  {
    if (!isEventSelected(cols)) {
      return;
    }
    histos.fill(HIST("VtxZHist"), cols.posZ());

    for (const auto& track : tracks) {
      if (!isTrackSelected(track)) {
        continue;
      }
      histos.fill(HIST("hdatahistpp"), cols.posZ(), track.pt(), track.phi(), kGlobalplusITS);
      if (track.hasTPC()) {
        histos.fill(HIST("hdatahistpp"), cols.posZ(), track.pt(), track.phi(), kGlobalonly);
      } else {
        histos.fill(HIST("hdatahistpp"), cols.posZ(), track.pt(), track.phi(), kITSonly);
      }
    }
  }

  void processMCeffPbPb(ColMCTrueTable::iterator const& mcCollision, ColMCRecTablePbPb const& RecCols, TrackMCTrueTable const& GenParticles, FilTrackMCRecTable const& RecTracks)
  {
    auto gencent = -999;
    bool atLeastOne = false;
    for (const auto& RecCol : RecCols) {
      if (!isEventSelected(RecCol)) {
        continue;
      }
      if (RecCol.globalIndex() != mcCollision.bestCollisionIndex()) {
        continue;
      }
      atLeastOne = true;
      gencent = RecCol.centFT0C();
    }
    histos.fill(HIST("hPbPbGenMCvtxz"), mcCollision.posZ());
    histos.fill(HIST("hPbPbGenMCvtxzcent"), mcCollision.posZ(), gencent);
    if (atLeastOne) {
      histos.fill(HIST("hPbPbGenMCAssoRecvtxz"), mcCollision.posZ());
      histos.fill(HIST("hPbPbGenMCAssoRecvtxzcent"), mcCollision.posZ(), gencent);
    }
    for (const auto& particle : GenParticles) {
      if (!isGenTrackSelected(particle)) {
        continue;
      }
      histos.fill(HIST("hPbPbGenMCdndpt"), mcCollision.posZ(), gencent, particle.pt(), particle.phi());
      if (atLeastOne) {
        histos.fill(HIST("hPbPbGenMCAssoRecdndpt"), mcCollision.posZ(), gencent, particle.pt(), particle.phi(), static_cast<double>(kGenAll), kNoGenpTVar);
        if (particle.pt() < KminPtCut) {
          histos.fill(HIST("hPbPbGenMCAssoRecdndpt"), mcCollision.posZ(), gencent, particle.pt(), particle.phi(), static_cast<double>(kGenAll), kGenpTup, -10.0 * particle.pt() + 2);
          histos.fill(HIST("hPbPbGenMCAssoRecdndpt"), mcCollision.posZ(), gencent, particle.pt(), particle.phi(), static_cast<double>(kGenAll), kGenpTdown, 5.0 * particle.pt() + 0.5);
        } else {
          histos.fill(HIST("hPbPbGenMCAssoRecdndpt"), mcCollision.posZ(), gencent, particle.pt(), particle.phi(), static_cast<double>(kGenAll), kGenpTup);
          histos.fill(HIST("hPbPbGenMCAssoRecdndpt"), mcCollision.posZ(), gencent, particle.pt(), particle.phi(), static_cast<double>(kGenAll), kGenpTdown);
        }
        int pid = 0;
        switch (std::abs(particle.pdgCode())) {
          case PDG_t::kPiPlus:
            pid = kGenPion;
            break;
          case PDG_t::kKPlus:
            pid = kGenKaon;
            break;
          case PDG_t::kProton:
            pid = kGenProton;
            break;
          default:
            pid = kGenOther;
            break;
        }
        histos.fill(HIST("hPbPbGenMCAssoRecdndpt"), mcCollision.posZ(), gencent, particle.pt(), particle.phi(), static_cast<double>(pid), kNoGenpTVar);
      } // Associated with reco col
    } // track (mcgen) loop

    for (const auto& RecCol : RecCols) {
      if (!isEventSelected(RecCol)) {
        continue;
      }
      if (RecCol.globalIndex() != mcCollision.bestCollisionIndex()) {
        continue;
      }
      histos.fill(HIST("hPbPbRecMCvtxz"), RecCol.posZ());
      histos.fill(HIST("hPbPbRecMCcent"), RecCol.centFT0C());
      histos.fill(HIST("hPbPbRecMCvtxzcent"), RecCol.posZ(), RecCol.centFT0C());
      auto recTracksPart = RecTracks.sliceBy(perCollision, RecCol.globalIndex());
      std::vector<int> mclabels;
      for (const auto& Rectrack : recTracksPart) {
        if (!isTrackSelected(Rectrack)) {
          continue;
        }
        histos.fill(HIST("hPbPbRecMCdndpt"), RecCol.posZ(), RecCol.centFT0C(), Rectrack.pt(), Rectrack.phi(), static_cast<double>(kRecoAll));
        if (Rectrack.has_mcParticle()) {
          int pid = 0;
          auto mcpart = Rectrack.mcParticle();
          histos.fill(HIST("hPbPbEtaReso"), Rectrack.pt(), Rectrack.pt() - mcpart.pt());
          histos.fill(HIST("hPbPbRecMCdndpt"), RecCol.posZ(), RecCol.centFT0C(), mcpart.pt(), mcpart.phi(), static_cast<double>(kRecoHasmc));
          if (mcpart.isPhysicalPrimary()) {
            histos.fill(HIST("hPbPbRecMCdndpt"), RecCol.posZ(), RecCol.centFT0C(), mcpart.pt(), mcpart.phi(), static_cast<double>(kRecoPrimary));
            switch (std::abs(mcpart.pdgCode())) {
              case PDG_t::kPiPlus:
                pid = kRecoPion;
                break;
              case PDG_t::kKPlus:
                pid = kRecoKaon;
                break;
              case PDG_t::kProton:
                pid = kRecoProton;
                break;
              default:
                pid = kRecoOther;
                break;
            }
          } else {
            pid = kRecoSecondary;
          }
          if (mcpart.has_mothers()) {
            auto mcpartMother = mcpart.template mothers_as<aod::McParticles>().front();
            if (mcpartMother.pdgCode() == PDG_t::kK0Short || std::abs(mcpartMother.pdgCode()) == PDG_t::kLambda0) {
              pid = kRecoWeakDecay;
            }
          }
          if (find(mclabels.begin(), mclabels.end(), Rectrack.mcParticleId()) != mclabels.end()) {
            pid = kRecoFake;
          }
          mclabels.push_back(Rectrack.mcParticleId());
          histos.fill(HIST("hPbPbRecMCdndpt"), RecCol.posZ(), RecCol.centFT0C(), mcpart.pt(), mcpart.phi(), static_cast<double>(pid));
        } else {
          histos.fill(HIST("hPbPbRecMCdndpt"), RecCol.posZ(), RecCol.centFT0C(), Rectrack.pt(), Rectrack.phi(), static_cast<double>(kRecoBkg));
        }
      } // track (mcrec) loop
    } // collision loop
  }

  void processMCeffpp(ColMCTrueTable::iterator const& mcCollision, ColMCRecTablepp const& RecCols, TrackMCTrueTable const& GenParticles, FilTrackMCRecTable const& RecTracks)
  {
    bool atLeastOne = false;
    for (const auto& RecCol : RecCols) {
      if (!isEventSelected(RecCol)) {
        continue;
      }
      if (RecCol.globalIndex() != mcCollision.bestCollisionIndex()) {
        continue;
      }
      atLeastOne = true;
    }
    histos.fill(HIST("hppGenMCvtxz"), mcCollision.posZ());
    if (atLeastOne) {
      histos.fill(HIST("hppGenMCAssoRecvtxz"), mcCollision.posZ());
    }
    for (const auto& particle : GenParticles) {
      if (!isGenTrackSelected(particle)) {
        continue;
      }
      histos.fill(HIST("hppGenMCdndpt"), mcCollision.posZ(), particle.pt(), particle.phi());
      if (atLeastOne) {
        histos.fill(HIST("hppGenMCAssoRecdndpt"), mcCollision.posZ(), particle.pt(), particle.phi(), static_cast<double>(kGenAll), kNoGenpTVar);
        if (particle.pt() < KminPtCut) {
          histos.fill(HIST("hppGenMCAssoRecdndpt"), mcCollision.posZ(), particle.pt(), particle.phi(), static_cast<double>(kGenAll), kGenpTup, -10.0 * particle.pt() + 2);
          histos.fill(HIST("hppGenMCAssoRecdndpt"), mcCollision.posZ(), particle.pt(), particle.phi(), static_cast<double>(kGenAll), kGenpTdown, 5.0 * particle.pt() + 0.5);
        } else {
          histos.fill(HIST("hppGenMCAssoRecdndpt"), mcCollision.posZ(), particle.pt(), particle.phi(), static_cast<double>(kGenAll), kGenpTup);
          histos.fill(HIST("hppGenMCAssoRecdndpt"), mcCollision.posZ(), particle.pt(), particle.phi(), static_cast<double>(kGenAll), kGenpTdown);
        }
        int pid = 0;
        switch (std::abs(particle.pdgCode())) {
          case PDG_t::kPiPlus:
            pid = kGenPion;
            break;
          case PDG_t::kKPlus:
            pid = kGenKaon;
            break;
          case PDG_t::kProton:
            pid = kGenProton;
            break;
          default:
            pid = kGenOther;
            break;
        }
        histos.fill(HIST("hppGenMCAssoRecdndpt"), mcCollision.posZ(), particle.pt(), particle.phi(), static_cast<double>(pid), kNoGenpTVar);
      } // Associated with reco col
    } // track (mcgen) loop

    for (const auto& RecCol : RecCols) {
      if (!isEventSelected(RecCol)) {
        continue;
      }
      if (RecCol.globalIndex() != mcCollision.bestCollisionIndex()) {
        continue;
      }
      histos.fill(HIST("hppRecMCvtxz"), RecCol.posZ());
      auto recTracksPart = RecTracks.sliceBy(perCollision, RecCol.globalIndex());
      std::vector<int> mclabels;
      for (const auto& Rectrack : recTracksPart) {
        if (!isTrackSelected(Rectrack)) {
          continue;
        }
        histos.fill(HIST("hppRecMCdndpt"), RecCol.posZ(), Rectrack.pt(), Rectrack.phi(), static_cast<double>(kRecoAll));
        if (Rectrack.has_mcParticle()) {
          int pid = 0;
          auto mcpart = Rectrack.mcParticle();
          histos.fill(HIST("hppEtaReso"), Rectrack.pt(), Rectrack.pt() - mcpart.pt());
          histos.fill(HIST("hppRecMCdndpt"), RecCol.posZ(), mcpart.pt(), mcpart.phi(), static_cast<double>(kRecoHasmc));
          if (mcpart.isPhysicalPrimary()) {
            histos.fill(HIST("hppRecMCdndpt"), RecCol.posZ(), mcpart.pt(), mcpart.phi(), static_cast<double>(kRecoPrimary));
            switch (std::abs(mcpart.pdgCode())) {
              case PDG_t::kPiPlus:
                pid = kRecoPion;
                break;
              case PDG_t::kKPlus:
                pid = kRecoKaon;
                break;
              case PDG_t::kProton:
                pid = kRecoProton;
                break;
              default:
                pid = kRecoOther;
                break;
            }
          } else {
            pid = kRecoSecondary;
          }
          if (mcpart.has_mothers()) {
            auto mcpartMother = mcpart.template mothers_as<aod::McParticles>().front();
            if (mcpartMother.pdgCode() == PDG_t::kK0Short || std::abs(mcpartMother.pdgCode()) == PDG_t::kLambda0) {
              pid = kRecoWeakDecay;
            }
          }
          if (find(mclabels.begin(), mclabels.end(), Rectrack.mcParticleId()) != mclabels.end()) {
            pid = kRecoFake;
          }
          mclabels.push_back(Rectrack.mcParticleId());
          histos.fill(HIST("hppRecMCdndpt"), RecCol.posZ(), mcpart.pt(), mcpart.phi(), static_cast<double>(pid));
        } else {
          histos.fill(HIST("hppRecMCdndpt"), RecCol.posZ(), Rectrack.pt(), Rectrack.phi(), static_cast<double>(kRecoBkg));
        }
      } // track (mcrec) loop
    } // collision loop
  }

  void processEvtLossSigLossMCpp(ColMCTrueTable::iterator const& mcCollision, ColMCRecTablepp const& RecCols, TrackMCTrueTable const& GenParticles)
  {
    if (isApplyInelgt0 && !mcCollision.isInelGt0()) {
      return;
    }
    if (std::abs(mcCollision.posZ()) >= vtxRange) {
      return;
    }
    bool atLeastOne = false;
    for (const auto& RecCol : RecCols) {
      if (!isEventSelected(RecCol)) {
        continue;
      }
      if (RecCol.globalIndex() != mcCollision.bestCollisionIndex()) {
        continue;
      }
      atLeastOne = true;
    }
    // All generated events
    histos.fill(HIST("MCEventHist"), 1);
    if (atLeastOne) {
      histos.fill(HIST("MCEventHist"), 2);
    }
    for (const auto& particle : GenParticles) {
      if (!isGenTrackSelected(particle)) {
        continue;
      }
      // All generated particles
      histos.fill(HIST("hgenptBeforeEvtSel"), particle.pt());
      if (atLeastOne) {
        // All generated particles with at least one reconstructed collision (signal loss estimation)
        histos.fill(HIST("hgenptAfterEvtSel"), particle.pt());
      }
    }
  }

  void processEvtLossSigLossMCPbPb(ColMCTrueTable::iterator const& mcCollision, ColMCRecTablePbPb const& RecCols, TrackMCTrueTable const& GenParticles)
  {
    if (isApplyInelgt0 && !mcCollision.isInelGt0()) {
      return;
    }
    if (std::abs(mcCollision.posZ()) >= vtxRange) {
      return;
    }
    bool atLeastOne = false;
    auto centrality = -999.;
    for (const auto& RecCol : RecCols) {
      if (!isEventSelected(RecCol)) {
        continue;
      }
      if (RecCol.globalIndex() != mcCollision.bestCollisionIndex()) {
        continue;
      }
      centrality = RecCol.centFT0C();
      atLeastOne = true;
    }
    // All generated events
    histos.fill(HIST("MCEventHist"), 1);
    histos.fill(HIST("hImpactParameterGen"), mcCollision.impactParameter());
    if (atLeastOne) {
      histos.fill(HIST("MCEventHist"), 2);
      histos.fill(HIST("hImpactParameterRec"), mcCollision.impactParameter());
      histos.fill(HIST("hImpactParvsCentrRec"), centrality, mcCollision.impactParameter());
    }
    for (const auto& particle : GenParticles) {
      if (!isGenTrackSelected(particle)) {
        continue;
      }
      // All generated particles
      histos.fill(HIST("hgenptBeforeEvtSelPbPb"), particle.pt(), mcCollision.impactParameter());
      if (atLeastOne) {
        // All generated particles with at least one reconstructed collision (signal loss estimation)
        histos.fill(HIST("hgenptAfterEvtSelPbPb"), particle.pt(), mcCollision.impactParameter());
      }
    }
  }

  PROCESS_SWITCH(PtmultCorr, processDataPbPb, "process data heavy-ion", false);
  PROCESS_SWITCH(PtmultCorr, processDatapp, "process data pp", false);
  PROCESS_SWITCH(PtmultCorr, processMCeffPbPb, "process MC heavy-ion", false);
  PROCESS_SWITCH(PtmultCorr, processMCeffpp, "process MC pp", false);
  PROCESS_SWITCH(PtmultCorr, processEvtLossSigLossMCpp, "process Signal Loss, Event Loss", false);
  PROCESS_SWITCH(PtmultCorr, processEvtLossSigLossMCPbPb, "process Signal Loss, Event Loss", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<PtmultCorr>(cfgc)};
}
