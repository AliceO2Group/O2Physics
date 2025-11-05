// Copyright 2019-2025 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
//
/// \file cascqaanalysis.cxx
/// \brief Analysis of cascades in pp collisions
/// \author Roman Nepeivoda (roman.nepeivoda@cern.ch)

#include "PWGLF/DataModel/cascqaanalysis.h"

#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/mcCentrality.h"
#include "PWGLF/Utils/inelGt.h"

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"

#include "TRandom2.h"
#include <TPDGCode.h>

#include <algorithm>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::rctsel;

// using DauTracks = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA, aod::pidTPCPi, aod::pidTPCPr, aod::pidTPCKa, aod::pidTOFPi, aod::pidTOFPr>;
using TrkPidInfo = soa::Join<aod::pidTPCFullPi, aod::pidTPCFullPr, aod::pidTPCFullKa, aod::pidTOFPi, aod::pidTOFPr, aod::pidTOFKa>;
using DauTracks = soa::Join<aod::TracksIU, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA, TrkPidInfo>;
using LabeledCascades = soa::Join<aod::CascDataExt, aod::McCascLabels>;

struct Cascqaanalysis {

  // Tables to produce
  Produces<aod::MyCascades> mycascades;
  Produces<aod::MyMCCascades> myMCcascades;

  HistogramRegistry registry{"registry"};

  // Axes
  ConfigurableAxis ptAxis{"ptAxis", {200, 0.0f, 10.0f}, "#it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis rapidityAxis{"rapidityAxis", {200, -2.0f, 2.0f}, "y"};
  ConfigurableAxis centFT0MAxis{"centFT0MAxis",
                                {VARIABLE_WIDTH, 0., 0.01, 0.05, 0.1, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 105.5},
                                "FT0M (%)"};
  ConfigurableAxis centFV0AAxis{"centFV0AAxis",
                                {VARIABLE_WIDTH, 0., 0.01, 0.05, 0.1, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 105.5},
                                "FV0A (%)"};
  ConfigurableAxis eventTypeAxis{"eventTypeAxis", {3, -0.5f, 2.5f}, "Event Type"};

  ConfigurableAxis nAssocCollAxis{"nAssocCollAxis", {5, -0.5f, 4.5f}, "N_{assoc.}"};
  ConfigurableAxis nChargedFT0MGenAxis{"nChargedFT0MGenAxis", {300, 0, 300}, "N_{FT0M, gen.}"};
  ConfigurableAxis nChargedFV0AGenAxis{"nChargedFV0AGenAxis", {300, 0, 300}, "N_{FV0A, gen.}"};
  ConfigurableAxis multNTracksAxis{"multNTracksAxis", {500, 0, 500}, "N_{tracks}"};
  ConfigurableAxis signalFT0MAxis{"signalFT0MAxis", {10000, 0, 40000}, "FT0M amplitude"};
  ConfigurableAxis signalFV0AAxis{"signalFV0AAxis", {10000, 0, 40000}, "FV0A amplitude"};
  ConfigurableAxis nCandidates{"nCandidates", {30, -0.5, 29.5}, "N_{cand.}"};

  // Event selection criteria
  Configurable<float> cutzvertex{"cutzvertex", 10.0f, "Accepted z-vertex range (cm)"};
  Configurable<bool> isVertexITSTPC{"isVertexITSTPC", 0, "Select collisions with at least one ITS-TPC track"};
  Configurable<bool> isNoSameBunchPileup{"isNoSameBunchPileup", 0, "Same found-by-T0 bunch crossing rejection"};
  Configurable<bool> isGoodZvtxFT0vsPV{"isGoodZvtxFT0vsPV", 0, "z of PV by tracks and z of PV from FT0 A-C time difference cut"};
  Configurable<bool> isVertexTOFmatched{"isVertexTOFmatched", 0, "Is Vertex TOF matched"};

  Configurable<bool> isTriggerTVX{"isTriggerTVX", 1, "TVX trigger"};
  Configurable<bool> isNoTimeFrameBorder{"isNoTimeFrameBorder", 1, "TF border cut"};
  Configurable<bool> isNoITSROFrameBorder{"isNoITSROFrameBorder", 1, "ITS ROF border cut"};
  Configurable<bool> isNoCollInTimeRangeNarrow{"isNoCollInTimeRangeNarrow", 1, "No collisions in +-2us window"};

  Configurable<bool> requireRCTFlagChecker{"requireRCTFlagChecker", true, "Check event quality in run condition table"};
  Configurable<std::string> cfgEvtRCTFlagCheckerLabel{"cfgEvtRCTFlagCheckerLabel", "CBT_hadronPID", "Evt sel: RCT flag checker label"};
  Configurable<bool> cfgEvtRCTFlagCheckerZDCCheck{"cfgEvtRCTFlagCheckerZDCCheck", false, "Evt sel: RCT flag checker ZDC check"};
  Configurable<bool> cfgEvtRCTFlagCheckerFV0Check{"cfgEvtRCTFlagCheckerFV0Check", false, "Evt sel: RCT flag checker FV0 check"};
  Configurable<bool> cfgEvtRCTFlagCheckerLimitAcceptAsBad{"cfgEvtRCTFlagCheckerLimitAcceptAsBad", false, "Evt sel: RCT flag checker treat Limited Acceptance As Bad"};

  RCTFlagsChecker rctChecker;

  // Cascade selection criteria
  Configurable<float> scalefactor{"scalefactor", 1.0, "Scaling factor"};
  Configurable<double> casccospa{"casccospa", 0.97, "Casc CosPA"};
  Configurable<double> v0cospa{"v0cospa", 0.97, "V0 CosPA"};
  Configurable<float> dcacascdau{"dcacascdau", 2.0, "DCA Casc Daughters"};
  Configurable<float> dcav0dau{"dcav0dau", 2.0, "DCA V0 Daughters"};
  Configurable<float> dcanegtopv{"dcanegtopv", 0.0, "DCA Neg To PV"};
  Configurable<float> dcapostopv{"dcapostopv", 0.0, "DCA Pos To PV"};
  Configurable<float> dcabachtopv{"dcabachtopv", 0.0, "DCA Bach To PV"};
  Configurable<float> v0radius{"v0radius", 0.0, "V0 Radius"};
  Configurable<float> cascradius{"cascradius", 0.0, "Casc Radius"};
  Configurable<float> etadau{"etadau", 0.8, "Eta Daughters"};

  // Switch between Data/MC-dedicated histograms
  Configurable<bool> isMC{"isMC", 0, "0 - data, 1 - MC"};

  // QA histograms for the multiplicity estimation
  Configurable<bool> multQA{"multQA", 0, "0 - not to do QA, 1 - do the QA"};

  // QA histograms for cascade rec.
  Configurable<bool> candidateQA{"candidateQA", 1, "0 - not to do QA, 1 - do the QA"};

  // Necessary for particle charges
  Service<o2::framework::O2DatabasePDG> pdgDB;

  SliceCache cache;

  // Random number generator for event scaling
  TRandom2* fRand = new TRandom2();

  // Struct to select on event type
  typedef struct CollisionIndexAndType {
    int64_t index;
    uint8_t typeFlag;
  } CollisionIndexAndType;

  void init(InitContext const&)
  {
    TString hCandidateCounterLabels[4] = {"All candidates", "passed topo cuts", "has associated MC particle", "associated with Xi(Omega)"};
    TString hNEventsMCLabels[6] = {"All", "z vrtx", "INEL", "INEL>0", "INEL>1", "Associated with rec. collision"};
    TString hNEventsLabels[14] = {"All", "kIsTriggerTVX", "kNoTimeFrameBorder", "kNoITSROFrameBorder", "kIsVertexITSTPC", "kNoSameBunchPileup", "kIsGoodZvtxFT0vsPV", "isVertexTOFmatched", "kNoCollInTimeRangeNarrow", "z vrtx", "RCTFlagsChecker", "INEL", "INEL>0", "INEL>1"};

    registry.add("hNEvents", "hNEvents", {HistType::kTH1D, {{14, 0.f, 14.f}}});

    for (int n = 1; n <= registry.get<TH1>(HIST("hNEvents"))->GetNbinsX(); n++) {
      registry.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(n, hNEventsLabels[n - 1]);
    }
    registry.add("hZCollision", "hZCollision", {HistType::kTH1D, {{200, -20.f, 20.f}}});

    registry.add("hCandidateCounter", "hCandidateCounter", {HistType::kTH1D, {{4, 0.0f, 4.0f}}});
    for (int n = 1; n <= registry.get<TH1>(HIST("hCandidateCounter"))->GetNbinsX(); n++) {
      registry.get<TH1>(HIST("hCandidateCounter"))->GetXaxis()->SetBinLabel(n, hCandidateCounterLabels[n - 1]);
    }
    if (isMC) {
      // Rec. lvl
      registry.add("fakeEvents", "fakeEvents", {HistType::kTH1F, {{1, -0.5f, 0.5f}}});
      // Gen. lvl
      registry.add("hNEventsMC", "hNEventsMC", {HistType::kTH1D, {{6, 0.0f, 6.0f}}});
      for (int n = 1; n <= registry.get<TH1>(HIST("hNEventsMC"))->GetNbinsX(); n++) {
        registry.get<TH1>(HIST("hNEventsMC"))->GetXaxis()->SetBinLabel(n, hNEventsMCLabels[n - 1]);
      }
      registry.add("hZCollisionGen", "hZCollisionGen", {HistType::kTH1D, {{200, -20.f, 20.f}}});
      registry.add("hNchFT0MNAssocMCCollisions", "hNchFT0MNAssocMCCollisions", {HistType::kTH3D, {nChargedFT0MGenAxis, nAssocCollAxis, eventTypeAxis}});
      registry.add("hNchFT0MNAssocMCCollisionsSameType", "hNchFT0MNAssocMCCollisionsSameType", {HistType::kTH3D, {nChargedFT0MGenAxis, nAssocCollAxis, eventTypeAxis}});
      registry.add("hNContributorsCorrelation", "hNContributorsCorrelation", {HistType::kTH2F, {{250, -0.5f, 249.5f, "Secondary Contributor"}, {250, -0.5f, 249.5f, "Main Contributor"}}});
      registry.add("hNchFT0MGenEvType", "hNchFT0MGenEvType", {HistType::kTH2D, {nChargedFT0MGenAxis, eventTypeAxis}});
      registry.add("hNchFV0AGenEvType", "hNchFV0AGenEvType", {HistType::kTH2D, {nChargedFV0AGenAxis, eventTypeAxis}});
      registry.add("hCentFT0M_genMC", "hCentFT0M_genMC", {HistType::kTH2D, {centFT0MAxis, eventTypeAxis}});
    }

    registry.add("hCentFT0M_rec", "hCentFT0M_rec", {HistType::kTH2D, {centFT0MAxis, eventTypeAxis}});

    if (candidateQA) {
      registry.add("hNcandidates", "hNcandidates", {HistType::kTH3D, {nCandidates, centFT0MAxis, {2, -0.5f, 1.5f}}});
    }

    if (multQA) {
      if (isMC) {
        // Rec. lvl
        registry.add("hNchFT0Mglobal", "hNchFT0Mglobal", {HistType::kTH3D, {nChargedFT0MGenAxis, multNTracksAxis, eventTypeAxis}});
        registry.add("hNchFT0MPVContr", "hNchFT0MPVContr", {HistType::kTH3D, {nChargedFT0MGenAxis, multNTracksAxis, eventTypeAxis}});
        registry.add("hNchFV0APVContr", "hNchFV0APVContr", {HistType::kTH3D, {nChargedFV0AGenAxis, multNTracksAxis, eventTypeAxis}});
      }
      registry.add("hFT0MpvContr", "hFT0MpvContr", {HistType::kTH3D, {centFT0MAxis, multNTracksAxis, eventTypeAxis}});
      registry.add("hFV0ApvContr", "hFV0ApvContr", {HistType::kTH3D, {centFV0AAxis, multNTracksAxis, eventTypeAxis}});
      registry.add("hFT0Mglobal", "hFT0Mglobal", {HistType::kTH3D, {centFT0MAxis, multNTracksAxis, eventTypeAxis}});
      registry.add("hFV0AFT0M", "hFV0AFT0M", {HistType::kTH3D, {centFV0AAxis, centFT0MAxis, eventTypeAxis}});
      registry.add("hFT0MFV0Asignal", "hFT0MFV0Asignal", {HistType::kTH2D, {signalFT0MAxis, signalFV0AAxis}});
      registry.add("hFT0MsignalPVContr", "hFT0MsignalPVContr", {HistType::kTH3D, {signalFT0MAxis, multNTracksAxis, eventTypeAxis}});
    }

    rctChecker.init(cfgEvtRCTFlagCheckerLabel, cfgEvtRCTFlagCheckerZDCCheck, cfgEvtRCTFlagCheckerLimitAcceptAsBad);
    if (cfgEvtRCTFlagCheckerFV0Check) {
      rctChecker.set(o2::aod::rctsel::kFV0Bad);
    }
  }

  Filter preFilter =
    (nabs(aod::cascdata::dcapostopv) > dcapostopv &&
     nabs(aod::cascdata::dcanegtopv) > dcanegtopv &&
     nabs(aod::cascdata::dcabachtopv) > dcabachtopv &&
     aod::cascdata::dcaV0daughters < dcav0dau &&
     aod::cascdata::dcacascdaughters < dcacascdau);

  Partition<DauTracks> pvContribTracksIUEta1 = (nabs(aod::track::eta) < 1.0f) && ((aod::track::flags & static_cast<uint32_t>(o2::aod::track::PVContributor)) == static_cast<uint32_t>(o2::aod::track::PVContributor));
  Partition<DauTracks> globalTracksIUEta05 = (nabs(aod::track::eta) < 0.5f) && (requireGlobalTrackInFilter());

  template <class TCascTracksTo, typename TCascade>
  bool acceptCascCandidate(TCascade const& cascCand, float const& pvx, float const& pvy, float const& pvz)
  {
    // Access daughter tracks
    auto posdau = cascCand.template posTrack_as<TCascTracksTo>();
    auto negdau = cascCand.template negTrack_as<TCascTracksTo>();
    auto bachelor = cascCand.template bachelor_as<TCascTracksTo>();

    // Basic set of selections
    if (cascCand.cascradius() > cascradius &&
        cascCand.v0radius() > v0radius &&
        cascCand.casccosPA(pvx, pvy, pvz) > casccospa &&
        cascCand.v0cosPA(pvx, pvy, pvz) > v0cospa &&
        std::fabs(posdau.eta()) < etadau &&
        std::fabs(negdau.eta()) < etadau &&
        std::fabs(bachelor.eta()) < etadau) {
      return true;
    } else {
      return false;
    }
  }

  template <typename TMcParticles>
  uint16_t getGenNchInFT0Mregion(TMcParticles particles)
  {
    // Particle counting in FITFT0: -3.3<η<-2.1; 3.5<η<4.9
    uint16_t nchFT0 = 0;
    for (const auto& mcParticle : particles) {
      if (!mcParticle.isPhysicalPrimary()) {
        continue;
      }
      const auto& pdgInfo = pdgDB->GetParticle(mcParticle.pdgCode());
      if (!pdgInfo) {
        continue;
      }
      if (pdgInfo->Charge() == 0) {
        continue;
      }
      if (mcParticle.eta() < -3.3 || mcParticle.eta() > 4.9 || (mcParticle.eta() > -2.1 && mcParticle.eta() < 3.5)) {
        continue; // select on T0M Nch region
      }
      nchFT0++; // increment
    }
    return nchFT0;
  }

  template <typename TMcParticles>
  uint16_t getGenNchInFV0Aregion(TMcParticles particles)
  {
    // Particle counting in FV0A: 2.2<η<5.1
    uint16_t nchFV0A = 0;
    for (const auto& mcParticle : particles) {
      if (!mcParticle.isPhysicalPrimary()) {
        continue;
      }
      const auto& pdgInfo = pdgDB->GetParticle(mcParticle.pdgCode());
      if (!pdgInfo) {
        continue;
      }
      if (pdgInfo->Charge() == 0) {
        continue;
      }
      if (mcParticle.eta() < 2.2 || mcParticle.eta() > 5.1) {
        continue; // select on V0A Nch region
      }
      nchFV0A++; // increment
    }
    return nchFV0A;
  }

  template <typename TCollision>
  int getEventTypeFlag(TCollision const& collision)
  {
    // 0 - INEL, 1 - INEL>0, 2 - INEL>1
    int evFlag = 0;
    registry.fill(HIST("hNEvents"), 11.5); // INEL
    if (collision.isInelGt0()) {
      evFlag += 1;
      registry.fill(HIST("hNEvents"), 12.5); // INEL>0
    }
    if (collision.isInelGt1()) {
      evFlag += 1;
      registry.fill(HIST("hNEvents"), 13.5); // INEL>1
    }
    return evFlag;
  }

  template <typename TCollision>
  bool acceptEvent(TCollision const& collision, bool isFillEventSelectionQA)
  {
    if (isFillEventSelectionQA) {
      registry.fill(HIST("hNEvents"), 0.5);
    }

    // kIsTriggerTVX selection
    if (isTriggerTVX && !collision.selection_bit(aod::evsel::kIsTriggerTVX)) {
      return false;
    }

    if (isFillEventSelectionQA) {
      registry.fill(HIST("hNEvents"), 1.5);
    }

    // kNoTimeFrameBorder selection
    if (isNoTimeFrameBorder && !collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
      return false;
    }

    if (isFillEventSelectionQA) {
      registry.fill(HIST("hNEvents"), 2.5);
    }

    // kNoITSROFrameBorder selection
    if (isNoITSROFrameBorder && !collision.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
      return false;
    }

    if (isFillEventSelectionQA) {
      registry.fill(HIST("hNEvents"), 3.5);
    }

    // kIsVertexITSTPC selection
    if (isVertexITSTPC && !collision.selection_bit(aod::evsel::kIsVertexITSTPC)) {
      return false;
    }
    if (isFillEventSelectionQA) {
      registry.fill(HIST("hNEvents"), 4.5);
    }
    // kNoSameBunchPileup selection
    if (isNoSameBunchPileup && !collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    if (isFillEventSelectionQA) {
      registry.fill(HIST("hNEvents"), 5.5);
    }
    // kIsGoodZvtxFT0vsPV selection
    if (isGoodZvtxFT0vsPV && !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    if (isFillEventSelectionQA) {
      registry.fill(HIST("hNEvents"), 6.5);
    }
    // isVertexTOFmatched selection
    if (isVertexTOFmatched && !collision.selection_bit(aod::evsel::kIsVertexTOFmatched)) {
      return false;
    }
    if (isFillEventSelectionQA) {
      registry.fill(HIST("hNEvents"), 7.5);
    }
    // kNoCollInTimeRangeNarrow selection
    if (isNoCollInTimeRangeNarrow && !collision.selection_bit(aod::evsel::kNoCollInTimeRangeNarrow)) {
      return false;
    }
    if (isFillEventSelectionQA) {
      registry.fill(HIST("hNEvents"), 8.5);
    }

    // Z vertex selection
    if (std::fabs(collision.posZ()) > cutzvertex) {
      return false;
    }
    if (isFillEventSelectionQA) {
      registry.fill(HIST("hNEvents"), 9.5);
      registry.fill(HIST("hZCollision"), collision.posZ());
    }

    // RCTFlagChecker selection
    if (requireRCTFlagChecker && !rctChecker(collision)) {
      return false;
    }
    if (isFillEventSelectionQA) {
      registry.fill(HIST("hNEvents"), 10.5);
    }

    return true;
  }

  void processData(soa::Join<aod::Collisions, aod::EvSels,
                             aod::PVMults, aod::FT0Mults, aod::FV0Mults,
                             aod::CentFT0Ms, aod::CentFV0As>::iterator const& collision,
                   soa::Filtered<aod::CascDataExt> const& Cascades,
                   aod::V0Datas const&,
                   DauTracks const&)
  {
    if (!acceptEvent(collision, 1)) {
      return;
    }

    int evType = getEventTypeFlag(collision);

    auto tracksGroupedPVcontr = pvContribTracksIUEta1->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    int nTracksPVcontr = tracksGroupedPVcontr.size();

    auto tracksGroupedGlobal = globalTracksIUEta05->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    int nTracksGlobal = tracksGroupedGlobal.size();

    registry.fill(HIST("hCentFT0M_rec"), collision.centFT0M(), evType);

    if (multQA) {
      registry.fill(HIST("hFT0MpvContr"), collision.centFT0M(), nTracksPVcontr, evType);
      registry.fill(HIST("hFV0ApvContr"), collision.centFV0A(), nTracksPVcontr, evType);
      registry.fill(HIST("hFT0Mglobal"), collision.centFT0M(), nTracksGlobal, evType);
      registry.fill(HIST("hFV0AFT0M"), collision.centFV0A(), collision.centFT0M(), evType);
      registry.fill(HIST("hFT0MFV0Asignal"), collision.multFT0A() + collision.multFT0C(), collision.multFV0A());
      registry.fill(HIST("hFT0MsignalPVContr"), collision.multFT0A() + collision.multFT0C(), nTracksPVcontr, evType);
    }

    float lEventScale = scalefactor;
    int nCandSel = 0;
    int nCandAll = 0;

    for (const auto& casc : Cascades) {              // loop over Cascades
      registry.fill(HIST("hCandidateCounter"), 0.5); // all candidates
      nCandAll++;
      if (acceptCascCandidate<DauTracks>(casc, collision.posX(), collision.posY(), collision.posZ())) {
        registry.fill(HIST("hCandidateCounter"), 1.5); // passed topo cuts
        nCandSel++;
        // Fill table
        if (fRand->Rndm() < lEventScale) {
          auto posdau = casc.posTrack_as<DauTracks>();
          auto negdau = casc.negTrack_as<DauTracks>();
          auto bachelor = casc.bachelor_as<DauTracks>();

          // ITS N hits
          int posITSNhits = 0, negITSNhits = 0, bachITSNhits = 0;
          for (unsigned int i = 0; i < 7; i++) {
            if (posdau.itsClusterMap() & (1 << i)) {
              posITSNhits++;
            }
            if (negdau.itsClusterMap() & (1 << i)) {
              negITSNhits++;
            }
            if (bachelor.itsClusterMap() & (1 << i)) {
              bachITSNhits++;
            }
          }

          uint8_t evFlag = 0;
          evFlag |= o2::aod::mycascades::EvFlags::EvINEL;
          if (collision.multNTracksPVeta1() > 0) {
            evFlag |= o2::aod::mycascades::EvFlags::EvINELgt0;
          }
          if (collision.multNTracksPVeta1() > 1) {
            evFlag |= o2::aod::mycascades::EvFlags::EvINELgt1;
          }

          // c x tau
          float cascpos = std::hypot(casc.x() - collision.posX(), casc.y() - collision.posY(), casc.z() - collision.posZ());
          float cascptotmom = std::hypot(casc.px(), casc.py(), casc.pz());
          float ctauXi = o2::constants::physics::MassXiMinus * cascpos / (cascptotmom + 1e-13);
          float ctauOmega = o2::constants::physics::MassOmegaMinus * cascpos / (cascptotmom + 1e-13);

          mycascades(collision.posZ(),
                     collision.centFT0M(), collision.centFV0A(),
                     collision.multFT0M(), collision.multFV0A(),
                     casc.sign(), casc.pt(), casc.yXi(), casc.yOmega(), casc.eta(),
                     casc.mXi(), casc.mOmega(), casc.mLambda(), casc.cascradius(), casc.v0radius(),
                     casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()), casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()),
                     casc.dcapostopv(), casc.dcanegtopv(), casc.dcabachtopv(), casc.dcacascdaughters(), casc.dcaV0daughters(), casc.dcav0topv(collision.posX(), collision.posY(), collision.posZ()),
                     posdau.eta(), negdau.eta(), bachelor.eta(), posITSNhits, negITSNhits, bachITSNhits,
                     ctauXi, ctauOmega, negdau.tpcNSigmaPr(), posdau.tpcNSigmaPr(), negdau.tpcNSigmaPi(), posdau.tpcNSigmaPi(), bachelor.tpcNSigmaPi(), bachelor.tpcNSigmaKa(),
                     negdau.tofNSigmaPr(), posdau.tofNSigmaPr(), negdau.tofNSigmaPi(), posdau.tofNSigmaPi(), bachelor.tofNSigmaPi(), bachelor.tofNSigmaKa(),
                     posdau.tpcNClsFound(), negdau.tpcNClsFound(), bachelor.tpcNClsFound(),
                     posdau.tpcNClsCrossedRows(), negdau.tpcNClsCrossedRows(), bachelor.tpcNClsCrossedRows(),
                     posdau.hasTOF(), negdau.hasTOF(), bachelor.hasTOF(),
                     posdau.pt(), negdau.pt(), bachelor.pt(), -1, -1, casc.bachBaryonCosPA(), casc.bachBaryonDCAxyToPV(), evFlag, 1e3, 1e3);
        }
      }
    }

    if (candidateQA) {
      registry.fill(HIST("hNcandidates"), nCandAll, collision.centFT0M(), 0);
      registry.fill(HIST("hNcandidates"), nCandSel, collision.centFT0M(), 1);
    }
  }

  Preslice<aod::McParticles> perMcCollision = aod::mcparticle::mcCollisionId;

  void processMCrec(soa::Join<aod::Collisions, o2::aod::McCollisionLabels, aod::EvSels,
                              aod::PVMults, aod::FT0Mults, aod::FV0Mults,
                              aod::CentFT0Ms, aod::CentFV0As>::iterator const& collision,
                    aod::V0Datas const&,
                    soa::Filtered<LabeledCascades> const& Cascades,
                    DauTracks const&,
                    soa::Join<aod::McCollisions, aod::McCentFT0Ms> const&, // aod::McCentFV0As to be added
                    aod::McParticles const& mcParticles)
  {
    if (!acceptEvent(collision, 1)) {
      return;
    }

    if (!collision.has_mcCollision()) {
      registry.fill(HIST("fakeEvents"), 0); // no assoc. MC collisions
      return;
    }

    const auto& mcCollision = collision.mcCollision_as<soa::Join<aod::McCollisions, aod::McCentFT0Ms>>(); // aod::McCentFV0As to be added

    auto tracksGroupedPVcontr = pvContribTracksIUEta1->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    int nTracksPVcontr = tracksGroupedPVcontr.size();

    auto tracksGroupedGlobal = globalTracksIUEta05->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    int nTracksGlobal = tracksGroupedGlobal.size();

    // N charged in FT0M region in corresponding gen. MC collision
    auto mcPartSlice = mcParticles.sliceBy(perMcCollision, collision.mcCollision_as<soa::Join<aod::McCollisions, aod::McCentFT0Ms>>().globalIndex()); // mcCollision.centFV0A() to be added
    uint16_t nchFT0 = getGenNchInFT0Mregion(mcPartSlice);
    uint16_t nchFV0 = getGenNchInFV0Aregion(mcPartSlice);

    int evType = 0;
    registry.fill(HIST("hNEvents"), 11.5); // INEL
    // Rec. collision associated with INEL>0 gen. one
    if (pwglf::isINELgtNmc(mcPartSlice, 0, pdgDB)) {
      registry.fill(HIST("hNEvents"), 12.5); // INEL
      evType++;
    }
    // Rec. collision associated with INEL>1 gen. one
    if (pwglf::isINELgtNmc(mcPartSlice, 1, pdgDB)) {
      registry.fill(HIST("hNEvents"), 13.5); // INEL
      evType++;
    }

    registry.fill(HIST("hCentFT0M_rec"), mcCollision.centFT0M(), evType);

    if (multQA) {
      registry.fill(HIST("hNchFT0MPVContr"), nchFT0, nTracksPVcontr, evType);
      registry.fill(HIST("hNchFV0APVContr"), nchFV0, nTracksPVcontr, evType);
      registry.fill(HIST("hFT0MpvContr"), mcCollision.centFT0M(), nTracksPVcontr, evType);
      registry.fill(HIST("hFV0ApvContr"), 0, nTracksPVcontr, evType); // mcCollision.centFV0A() to be added
      registry.fill(HIST("hFT0Mglobal"), mcCollision.centFT0M(), nTracksGlobal, evType);
      registry.fill(HIST("hFV0AFT0M"), 0, mcCollision.centFT0M(), evType); // mcCollision.centFV0A() to be added
      registry.fill(HIST("hNchFT0Mglobal"), nchFT0, nTracksGlobal, evType);
      registry.fill(HIST("hFT0MFV0Asignal"), collision.multFT0A() + collision.multFT0C(), collision.multFV0A());
      registry.fill(HIST("hFT0MsignalPVContr"), collision.multFT0A() + collision.multFT0C(), nTracksPVcontr, evType);
    }

    float lEventScale = scalefactor;
    int nCandSel = 0;
    int nCandAll = 0;

    for (const auto& casc : Cascades) {              // loop over Cascades
      registry.fill(HIST("hCandidateCounter"), 0.5); // all candidates
      nCandAll++;
      if (acceptCascCandidate<DauTracks>(casc, collision.posX(), collision.posY(), collision.posZ())) {
        registry.fill(HIST("hCandidateCounter"), 1.5); // passed topo cuts
        nCandSel++;
        // Check mc association
        float lPDG = 1e3;
        float genPt = 1e3;
        float genY = 1e3;
        float isPrimary = -1;
        if (casc.has_mcParticle()) {
          registry.fill(HIST("hCandidateCounter"), 2.5); // has associated MC particle
          auto cascmc = casc.mcParticle();
          if (std::abs(cascmc.pdgCode()) == PDG_t::kXiMinus || std::abs(cascmc.pdgCode()) == PDG_t::kOmegaMinus) {
            registry.fill(HIST("hCandidateCounter"), 3.5); // associated with Xi or Omega
            lPDG = cascmc.pdgCode();
            isPrimary = cascmc.isPhysicalPrimary() ? 1 : 0;
            genPt = cascmc.pt();
            genY = cascmc.y();
          }
        }
        if (fRand->Rndm() < lEventScale) {
          // Fill table
          auto posdau = casc.posTrack_as<DauTracks>();
          auto negdau = casc.negTrack_as<DauTracks>();
          auto bachelor = casc.bachelor_as<DauTracks>();

          // ITS N hits
          int posITSNhits = 0, negITSNhits = 0, bachITSNhits = 0;
          for (unsigned int i = 0; i < 7; i++) {
            if (posdau.itsClusterMap() & (1 << i)) {
              posITSNhits++;
            }
            if (negdau.itsClusterMap() & (1 << i)) {
              negITSNhits++;
            }
            if (bachelor.itsClusterMap() & (1 << i)) {
              bachITSNhits++;
            }
          }

          // Event type flag
          uint8_t evFlag = 0;
          evFlag |= o2::aod::mycascades::EvFlags::EvINEL;
          if (collision.multNTracksPVeta1() > 0) {
            evFlag |= o2::aod::mycascades::EvFlags::EvINELgt0;
          }
          if (collision.multNTracksPVeta1() > 1) {
            evFlag |= o2::aod::mycascades::EvFlags::EvINELgt1;
          }

          // c x tau
          float cascpos = std::hypot(casc.x() - collision.posX(), casc.y() - collision.posY(), casc.z() - collision.posZ());
          float cascptotmom = std::hypot(casc.px(), casc.py(), casc.pz());
          float ctauXi = o2::constants::physics::MassXiMinus * cascpos / (cascptotmom + 1e-13);
          float ctauOmega = o2::constants::physics::MassOmegaMinus * cascpos / (cascptotmom + 1e-13);

          mycascades(collision.posZ(),
                     mcCollision.centFT0M(), 0, // mcCollision.centFV0A() to be added
                     collision.multFT0M(), collision.multFV0A(),
                     casc.sign(), casc.pt(), casc.yXi(), casc.yOmega(), casc.eta(),
                     casc.mXi(), casc.mOmega(), casc.mLambda(), casc.cascradius(), casc.v0radius(),
                     casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()), casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()),
                     casc.dcapostopv(), casc.dcanegtopv(), casc.dcabachtopv(), casc.dcacascdaughters(), casc.dcaV0daughters(), casc.dcav0topv(collision.posX(), collision.posY(), collision.posZ()),
                     posdau.eta(), negdau.eta(), bachelor.eta(), posITSNhits, negITSNhits, bachITSNhits,
                     ctauXi, ctauOmega, negdau.tpcNSigmaPr(), posdau.tpcNSigmaPr(), negdau.tpcNSigmaPi(), posdau.tpcNSigmaPi(), bachelor.tpcNSigmaPi(), bachelor.tpcNSigmaKa(),
                     negdau.tofNSigmaPr(), posdau.tofNSigmaPr(), negdau.tofNSigmaPi(), posdau.tofNSigmaPi(), bachelor.tofNSigmaPi(), bachelor.tofNSigmaKa(),
                     posdau.tpcNClsFound(), negdau.tpcNClsFound(), bachelor.tpcNClsFound(),
                     posdau.tpcNClsCrossedRows(), negdau.tpcNClsCrossedRows(), bachelor.tpcNClsCrossedRows(),
                     posdau.hasTOF(), negdau.hasTOF(), bachelor.hasTOF(),
                     posdau.pt(), negdau.pt(), bachelor.pt(), lPDG, isPrimary, casc.bachBaryonCosPA(), casc.bachBaryonDCAxyToPV(), evFlag, genPt, genY);
        }
      }
    }

    if (candidateQA) {
      registry.fill(HIST("hNcandidates"), nCandAll, mcCollision.centFT0M(), 0);
      registry.fill(HIST("hNcandidates"), nCandSel, mcCollision.centFT0M(), 1);
    }
  }

  void processMCgen(soa::Join<aod::McCollisions, aod::McCentFT0Ms>::iterator const& mcCollision, // mcCollision.centFV0A() to be added
                    aod::McParticles const& mcParticles,
                    const soa::SmallGroups<o2::soa::Join<o2::aod::Collisions, o2::aod::McCollisionLabels, o2::aod::EvSels, aod::PVMults, aod::FT0Mults, aod::CentFT0Ms, aod::CentFV0As>>& collisions,
                    DauTracks const&)
  {
    // All generated collisions
    registry.fill(HIST("hNEventsMC"), 0.5);

    // Generated with accepted z vertex
    if (std::fabs(mcCollision.posZ()) > cutzvertex) {
      return;
    }
    registry.fill(HIST("hZCollisionGen"), mcCollision.posZ());
    registry.fill(HIST("hNEventsMC"), 1.5);

    // Define the type of generated MC collision
    int evType = 0;
    uint8_t flagsGen = 0;
    flagsGen |= o2::aod::myMCcascades::EvFlags::EvINEL;
    registry.fill(HIST("hNEventsMC"), 2.5);
    // Generated collision is INEL>0
    if (pwglf::isINELgtNmc(mcParticles, 0, pdgDB)) {
      flagsGen |= o2::aod::myMCcascades::EvFlags::EvINELgt0;
      evType++;
      registry.fill(HIST("hNEventsMC"), 3.5);
    }
    // Generated collision is INEL>1
    if (pwglf::isINELgtNmc(mcParticles, 1, pdgDB)) {
      flagsGen |= o2::aod::myMCcascades::EvFlags::EvINELgt1;
      evType++;
      registry.fill(HIST("hNEventsMC"), 4.5);
    }

    registry.fill(HIST("hCentFT0M_genMC"), mcCollision.centFT0M(), evType);

    uint16_t nchFT0 = getGenNchInFT0Mregion(mcParticles);
    uint16_t nchFV0 = getGenNchInFV0Aregion(mcParticles);
    registry.fill(HIST("hNchFT0MGenEvType"), nchFT0, evType);
    registry.fill(HIST("hNchFV0AGenEvType"), nchFV0, evType);

    std::vector<CollisionIndexAndType> selectedEvents(collisions.size());
    std::vector<int64_t> numberOfContributors;
    int nevts = 0;
    int nAssocColl = 0;
    for (const auto& collision : collisions) {
      CollisionIndexAndType collWithType = {0, 0x0};
      if (!acceptEvent(collision, 0)) {
        continue;
      }
      collWithType.index = collision.mcCollision_as<soa::Join<aod::McCollisions, aod::McCentFT0Ms>>().globalIndex(); // mcCollision.centFV0A() to be added
      collWithType.typeFlag |= o2::aod::myMCcascades::EvFlags::EvINEL;

      if (collision.isInelGt0()) {
        collWithType.typeFlag |= o2::aod::myMCcascades::EvFlags::EvINELgt0;
      }
      if (collision.isInelGt1()) {
        collWithType.typeFlag |= o2::aod::myMCcascades::EvFlags::EvINELgt1;
      }

      selectedEvents[nevts++] = collWithType;
      if (collision.mcCollision_as<soa::Join<aod::McCollisions, aod::McCentFT0Ms>>().globalIndex() == mcCollision.globalIndex()) { // mcCollision.centFV0A() to be added
        nAssocColl++;
        numberOfContributors.push_back(collision.numContrib());
      }
    }
    selectedEvents.resize(nevts);

    registry.fill(HIST("hNchFT0MNAssocMCCollisions"), nchFT0, nAssocColl, evType);

    if (numberOfContributors.size() == 2) {
      std::sort(numberOfContributors.begin(), numberOfContributors.end());
      registry.fill(HIST("hNContributorsCorrelation"), numberOfContributors[0], numberOfContributors[1]);
    }

    auto isAssocToINEL = [&mcCollision](CollisionIndexAndType i) { return (i.index == mcCollision.globalIndex()) && ((i.typeFlag & o2::aod::myMCcascades::EvFlags::EvINEL) == o2::aod::myMCcascades::EvFlags::EvINEL); };
    auto isAssocToINELgt0 = [&mcCollision](CollisionIndexAndType i) { return (i.index == mcCollision.globalIndex()) && ((i.typeFlag & o2::aod::myMCcascades::EvFlags::EvINELgt0) == o2::aod::myMCcascades::EvFlags::EvINELgt0); };
    auto isAssocToINELgt1 = [&mcCollision](CollisionIndexAndType i) { return (i.index == mcCollision.globalIndex()) && ((i.typeFlag & o2::aod::myMCcascades::EvFlags::EvINELgt1) == o2::aod::myMCcascades::EvFlags::EvINELgt1); };
    // number of reconstructed INEL events that have the same global index as mcCollision
    const auto evtReconstructedAndINEL = std::count_if(selectedEvents.begin(), selectedEvents.end(), isAssocToINEL);
    // number of reconstructed INEL > 0 events that have the same global index as mcCollision
    const auto evtReconstructedAndINELgt0 = std::count_if(selectedEvents.begin(), selectedEvents.end(), isAssocToINELgt0);
    // number of reconstructed INEL > 1 events that have the same global index as mcCollision
    const auto evtReconstructedAndINELgt1 = std::count_if(selectedEvents.begin(), selectedEvents.end(), isAssocToINELgt1);

    switch (evType) {
      case 0: {
        registry.fill(HIST("hNchFT0MNAssocMCCollisionsSameType"), nchFT0, evtReconstructedAndINEL, evType);
        break;
      }
      case 1: {
        registry.fill(HIST("hNchFT0MNAssocMCCollisionsSameType"), nchFT0, evtReconstructedAndINELgt0, evType);
        break;
      }
      case 2: {
        registry.fill(HIST("hNchFT0MNAssocMCCollisionsSameType"), nchFT0, evtReconstructedAndINELgt1, evType);
        break;
      }
      default:
        LOGF(fatal, "incorrect evType in cascqaanalysis task");
        break;
    }

    uint8_t flagsAssoc = 0;
    if (evtReconstructedAndINEL) {
      flagsAssoc |= o2::aod::myMCcascades::EvFlags::EvINEL;
      registry.fill(HIST("hNEventsMC"), 5.5);
    }
    if (evtReconstructedAndINELgt0) {
      flagsAssoc |= o2::aod::myMCcascades::EvFlags::EvINELgt0;
    }
    if (evtReconstructedAndINELgt1) {
      flagsAssoc |= o2::aod::myMCcascades::EvFlags::EvINELgt1;
    }

    for (const auto& mcParticle : mcParticles) {
      float sign = 0;
      if (mcParticle.pdgCode() == PDG_t::kXiPlusBar || mcParticle.pdgCode() == PDG_t::kOmegaPlusBar) {
        sign = 1;
      } else if (mcParticle.pdgCode() == PDG_t::kXiMinus || mcParticle.pdgCode() == PDG_t::kOmegaMinus) {
        sign = -1;
      } else {
        continue;
      }
      myMCcascades(mcCollision.posZ(), sign, mcParticle.pdgCode(),
                   mcParticle.y(), mcParticle.eta(), mcParticle.phi(), mcParticle.pt(),
                   mcParticle.isPhysicalPrimary(), nAssocColl,
                   nchFT0, nchFV0,
                   mcCollision.centFT0M(), 0, // mcCollision.centFV0A() to be added
                   flagsAssoc,
                   flagsGen);
    }
  }

  PROCESS_SWITCH(Cascqaanalysis, processData, "Process Run 3 data", true);
  PROCESS_SWITCH(Cascqaanalysis, processMCrec, "Process Run 3 mc, reconstructed", false);
  PROCESS_SWITCH(Cascqaanalysis, processMCgen, "Process Run 3 mc, genereated", false);
};

struct MyCascades {

  HistogramRegistry registry{"registry"};

  // QA
  Configurable<bool> doQA{"doQA", 0, "Fill QA histograms"};

  void init(InitContext const&)
  {
    TString pdglabels[3] = {"Unknown", "3312", "3334"};
    registry.add("hPt", "hPt", {HistType::kTH1F, {{100, 0.0f, 10.0f}}});
    registry.add("hMassXi", "hMassXi", {HistType::kTH1F, {{1000, 1.0f, 2.0f}}});
    registry.add("hMassOmega", "hMassOmega", {HistType::kTH1F, {{1000, 1.0f, 2.0f}}});
    if (doQA) {
      registry.add("hCascRadius", "hCascRadius", {HistType::kTH1D, {{100, 0.0f, 40.0f}}});
      registry.add("hV0Radius", "hV0Radius", {HistType::kTH1D, {{100, 0.0f, 40.0f}}});
      registry.add("hCascCosPA", "hCascCosPA", {HistType::kTH1F, {{100, 0.9f, 1.0f}}});
      registry.add("hV0CosPA", "hV0CosPA", {HistType::kTH1F, {{100, 0.9f, 1.0f}}});
      registry.add("hDCANegToPV", "hDCANegToPV", {HistType::kTH1F, {{100, -1.0f, 1.0f}}});
      registry.add("hDCAPosToPV", "hDCAPosToPV", {HistType::kTH1F, {{100, -1.0f, 1.0f}}});
      registry.add("hDCABachToPV", "hDCABachToPV", {HistType::kTH1F, {{100, -1.0f, 1.0f}}});
      registry.add("hDCACascDaughters", "hDCACascDaughters", {HistType::kTH1F, {{55, 0.0f, 2.20f}}});
      registry.add("hDCAV0Daughters", "hDCAV0Daughters", {HistType::kTH1F, {{55, 0.0f, 2.20f}}});
      registry.add("hCtauXi", "hCtauXi", {HistType::kTH1F, {{100, 0.0f, 40.0f}}});
      registry.add("hCtauOmega", "hCtauOmega", {HistType::kTH1F, {{100, 0.0f, 40.0f}}});
      registry.add("hTPCNSigmaPosPi", "hTPCNSigmaPosPi", {HistType::kTH1F, {{100, -10.0f, 10.0f}}});
      registry.add("hTPCNSigmaNegPi", "hTPCNSigmaNegPi", {HistType::kTH1F, {{100, -10.0f, 10.0f}}});
      registry.add("hTPCNSigmaPosPr", "hTPCNSigmaPosPr", {HistType::kTH1F, {{100, -10.0f, 10.0f}}});
      registry.add("hTPCNSigmaNegPr", "hTPCNSigmaNegPr", {HistType::kTH1F, {{100, -10.0f, 10.0f}}});
      registry.add("hTPCNSigmaBachPi", "hTPCNSigmaBachPi", {HistType::kTH1F, {{100, -10.0f, 10.0f}}});
      registry.add("hTOFNSigmaPosPi", "hTOFNSigmaPosPi", {HistType::kTH1F, {{100, -10.0f, 10.0f}}});
      registry.add("hTOFNSigmaNegPi", "hTOFNSigmaNegPi", {HistType::kTH1F, {{100, -10.0f, 10.0f}}});
      registry.add("hTOFNSigmaPosPr", "hTOFNSigmaPosPr", {HistType::kTH1F, {{100, -10.0f, 10.0f}}});
      registry.add("hTOFNSigmaNegPr", "hTOFNSigmaNegPr", {HistType::kTH1F, {{100, -10.0f, 10.0f}}});
      registry.add("hTOFNSigmaBachPi", "hTOFNSigmaBachPi", {HistType::kTH1F, {{100, -10.0f, 10.0f}}});
      registry.add("hPosITSHits", "hPosITSHits", {HistType::kTH1F, {{8, -0.5f, 7.5f}}});
      registry.add("hNegITSHits", "hNegITSHits", {HistType::kTH1F, {{8, -0.5f, 7.5f}}});
      registry.add("hBachITSHits", "hBachITSHits", {HistType::kTH1F, {{8, -0.5f, 7.5f}}});
      registry.add("hIsPrimary", "hIsPrimary", {HistType::kTH1F, {{3, -1.5f, 1.5f}}});
      registry.add("hPDGcode", "hPDGcode", {HistType::kTH1F, {{3, -1.5f, 1.5f}}});
      for (int n = 1; n <= registry.get<TH1>(HIST("hPDGcode"))->GetNbinsX(); n++) {
        registry.get<TH1>(HIST("hPDGcode"))->GetXaxis()->SetBinLabel(n, pdglabels[n - 1]);
      }
      registry.add("hBachBaryonCosPA", "hBachBaryonCosPA", {HistType::kTH1F, {{100, 0.0f, 1.0f}}});
      registry.add("hBachBaryonDCAxyToPV", "hBachBaryonDCAxyToPV", {HistType::kTH1F, {{300, -3.0f, 3.0f}}});
    }
  }

  void process(aod::MyCascades const& mycascades)
  {
    for (const auto& candidate : mycascades) {

      registry.fill(HIST("hMassXi"), candidate.massxi());
      registry.fill(HIST("hMassOmega"), candidate.massomega());
      registry.fill(HIST("hPt"), candidate.pt());
      if (doQA) {
        registry.fill(HIST("hCascRadius"), candidate.cascradius());
        registry.fill(HIST("hV0Radius"), candidate.v0radius());
        registry.fill(HIST("hCascCosPA"), candidate.casccospa());
        registry.fill(HIST("hV0CosPA"), candidate.v0cospa());
        registry.fill(HIST("hDCANegToPV"), candidate.dcanegtopv());
        registry.fill(HIST("hDCAPosToPV"), candidate.dcapostopv());
        registry.fill(HIST("hDCABachToPV"), candidate.dcabachtopv());
        registry.fill(HIST("hDCACascDaughters"), candidate.dcacascdaughters());
        registry.fill(HIST("hDCAV0Daughters"), candidate.dcav0daughters());
        registry.fill(HIST("hCtauXi"), candidate.ctauxi());
        registry.fill(HIST("hCtauOmega"), candidate.ctauomega());
        registry.fill(HIST("hTPCNSigmaPosPi"), candidate.ntpcsigmapospi());
        registry.fill(HIST("hTPCNSigmaNegPi"), candidate.ntpcsigmanegpi());
        registry.fill(HIST("hTPCNSigmaPosPr"), candidate.ntpcsigmapospr());
        registry.fill(HIST("hTPCNSigmaNegPr"), candidate.ntpcsigmanegpr());
        registry.fill(HIST("hTPCNSigmaBachPi"), candidate.ntpcsigmabachpi());
        registry.fill(HIST("hTOFNSigmaPosPi"), candidate.ntofsigmapospi());
        registry.fill(HIST("hTOFNSigmaNegPi"), candidate.ntofsigmanegpi());
        registry.fill(HIST("hTOFNSigmaPosPr"), candidate.ntofsigmapospr());
        registry.fill(HIST("hTOFNSigmaNegPr"), candidate.ntofsigmanegpr());
        registry.fill(HIST("hTOFNSigmaBachPi"), candidate.ntofsigmabachpi());
        registry.fill(HIST("hPosITSHits"), candidate.positshits());
        registry.fill(HIST("hNegITSHits"), candidate.negitshits());
        registry.fill(HIST("hBachITSHits"), candidate.bachitshits());
        registry.fill(HIST("hIsPrimary"), candidate.isPrimary());
        registry.fill(HIST("hBachBaryonCosPA"), candidate.bachBaryonCosPA());
        registry.fill(HIST("hBachBaryonDCAxyToPV"), candidate.bachBaryonDCAxyToPV());

        if (std::abs(candidate.mcPdgCode()) == PDG_t::kXiMinus || std::abs(candidate.mcPdgCode()) == PDG_t::kOmegaMinus) {
          registry.fill(HIST("hPDGcode"), std::abs(candidate.mcPdgCode()) == PDG_t::kXiMinus ? 0 : 1); // 0 if Xi, 1 if Omega
        } else {
          registry.fill(HIST("hPDGcode"), -1); // -1 if unknown
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<Cascqaanalysis>(cfgc),
    adaptAnalysisTask<MyCascades>(cfgc)};
}
