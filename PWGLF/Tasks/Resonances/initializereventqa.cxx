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
/// \file initializereventqa.cxx
/// \brief QA for the event-loss and signal-loss correction at the generator level for the ResonanceInitializer in pp collisions (referred to TableProducer/Strangeness/cascqaanalysis.cxx)
///
///  Following the discussions at the two PAG meetings (https://indico.cern.ch/event/1518979, https://indico.cern.ch/event/1575984)
///  we have introduced an auxiliary task that, when the resonanceInitializer.cxx is used,
///  computes the event-loss and signal-loss correction factors at the generator level.
///  With minor configuration tuning for a truth-tagging,
///  we expect it to be applicable to most analyses that rely on the initializer.
///
/// \author Minjae Kim (minjae.kim@cern.ch)

#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/cascqaanalysis.h"
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
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using TrkPidInfo = soa::Join<aod::pidTPCFullPi, aod::pidTPCFullPr, aod::pidTPCFullKa, aod::pidTOFPi, aod::pidTOFPr, aod::pidTOFKa>;
using DauTracks = soa::Join<aod::TracksIU, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA, TrkPidInfo>;

struct Initializereventqa {

  // Templates used, new hedder wil be added
  Produces<aod::MyCascades> mycascades;
  Produces<aod::MyMCCascades> myMCcascades;

  HistogramRegistry registry{"registry"};

  // Axes
  ConfigurableAxis ptAxis{"ptAxis", {400, 0.0f, 20.0f}, "#it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis rapidityAxis{"rapidityAxis", {200, -2.0f, 2.0f}, "y"};
  ConfigurableAxis centFT0MAxis{"centFT0MAxis",
                                {VARIABLE_WIDTH, 0.0, 1.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0, 110.0},
                                "FT0M (%)"};
  ConfigurableAxis eventTypeAxis{"eventTypeAxis", {2, -0.5f, 1.5f}, "Event Type"};

  ConfigurableAxis nAssocCollAxis{"nAssocCollAxis", {5, -0.5f, 4.5f}, "N_{assoc.}"};
  ConfigurableAxis nChargedFT0MGenAxis{"nChargedFT0MGenAxis", {300, 0, 300}, "N_{FT0M, gen.}"};
  ConfigurableAxis multNTracksAxis{"multNTracksAxis", {500, 0, 500}, "N_{tracks}"};
  ConfigurableAxis signalFT0MAxis{"signalFT0MAxis", {4000, 0, 40000}, "FT0M amplitude"};
  ConfigurableAxis signalFV0AAxis{"signalFV0AAxis", {4000, 0, 40000}, "FV0A amplitude"};
  ConfigurableAxis nCandidates{"nCandidates", {30, -0.5, 29.5}, "N_{cand.}"};

  // Event selection criteria
  Configurable<float> cutzvertex{"cutzvertex", 10.0f, "Accepted z-vertex range (cm)"};
  Configurable<bool> isZvtxcut{"isZvtxcut", 1, "Select collisions with Accepted z-vertex"};
  Configurable<bool> isVertexITSTPC{"isVertexITSTPC", 0, "Select collisions with at least one ITS-TPC track"};
  Configurable<bool> isNoSameBunchPileup{"isNoSameBunchPileup", 0, "Same found-by-T0 bunch crossing rejection"};
  Configurable<bool> isGoodZvtxFT0vsPV{"isGoodZvtxFT0vsPV", 0, "z of PV by tracks and z of PV from FT0 A-C time difference cut"};
  Configurable<bool> isVertexTOFmatched{"isVertexTOFmatched", 0, "Is Vertex TOF matched"};

  Configurable<bool> isTriggerTVX{"isTriggerTVX", 1, "TVX trigger"};
  Configurable<bool> isNoTimeFrameBorder{"isNoTimeFrameBorder", 1, "TF border cut"};
  Configurable<bool> isNoITSROFrameBorder{"isNoITSROFrameBorder", 1, "ITS ROF border cut"};
  Configurable<bool> isNoCollInTimeRangeNarrow{"isNoCollInTimeRangeNarrow", 0, "No collisions in +-2us window"};

  // QA histograms for the multiplicity estimation
  Configurable<bool> multQA{"multQA", 1, "0 - not to do QA, 1 - do the QA"};

  // Selection for signal-loss corrections
  Configurable<bool> isDaughterCheck{"isDaughterCheck", 1, "Check if the candidate has the correct daughters when it is considered"};

  Configurable<float> cfgRapidityCut{"cfgRapidityCut", 0.5, "Rapidity cut for the truth particle"};

  Configurable<int> pdgTruthMother{"pdgTruthMother", 3324, "pdgcode for the truth mother particle, e.g. Xi(1530) (3324)"};
  Configurable<int> pdgTruthDaughter1{"pdgTruthDaughter1", 3312, "pdgcode for the first daughter particle, e.g. Xi-3312"};
  Configurable<int> pdgTruthDaughter2{"pdgTruthDaughter2", 211, "pdgcode for the second daughter particle, e.g. Xi-3312"};

  // Necessary for particle charges
  Service<o2::framework::O2DatabasePDG> pdgDB;

  SliceCache cache;

  // Struct to select on event type
  typedef struct CollisionIndexAndType {
    int64_t index;
    uint8_t typeFlag;
  } CollisionIndexAndType;

  void init(InitContext const&)
  {
    TString hNEventsMCLabels[5] = {"All", "z vrtx", "INEL", "INEL>0", "Associated with rec. collision"};
    TString hNEventsLabels[12] = {"All", "kIsTriggerTVX", "kNoTimeFrameBorder", "kNoITSROFrameBorder", "kIsVertexITSTPC", "kNoSameBunchPileup", "kIsGoodZvtxFT0vsPV", "isVertexTOFmatched", "kNoCollInTimeRangeNarrow", "z vrtx", "INEL", "INEL>0"};

    registry.add("hNEvents", "hNEvents", {HistType::kTH1D, {{12, 0.f, 12.f}}});

    for (int n = 1; n <= registry.get<TH1>(HIST("hNEvents"))->GetNbinsX(); n++) {
      registry.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(n, hNEventsLabels[n - 1]);
    }
    registry.add("hZCollision", "hZCollision", {HistType::kTH1D, {{200, -20.f, 20.f}}});

    registry.add("fakeEvents", "fakeEvents", {HistType::kTH1F, {{1, -0.5f, 0.5f}}});

    registry.add("hNEventsMC", "hNEventsMC", {HistType::kTH1D, {{5, 0.0f, 5.0f}}});
    for (int n = 1; n <= registry.get<TH1>(HIST("hNEventsMC"))->GetNbinsX(); n++) {
      registry.get<TH1>(HIST("hNEventsMC"))->GetXaxis()->SetBinLabel(n, hNEventsMCLabels[n - 1]);
    }
    registry.add("hZCollisionGen", "hZCollisionGen", {HistType::kTH1D, {{200, -20.f, 20.f}}});
    registry.add("hCentFT0MNAssocMCCollisions", "hCentFT0MNAssocMCCollisions", {HistType::kTH3D, {centFT0MAxis, nAssocCollAxis, eventTypeAxis}});
    registry.add("hCentFT0MNAssocMCCollisionsSameType", "hCentFT0MNAssocMCCollisionsSameType", {HistType::kTH3D, {centFT0MAxis, nAssocCollAxis, eventTypeAxis}});
    registry.add("hNchFT0MNAssocMCCollisions", "hNchFT0MNAssocMCCollisions", {HistType::kTH3D, {nChargedFT0MGenAxis, nAssocCollAxis, eventTypeAxis}});
    registry.add("hNchFT0MNAssocMCCollisionsSameType", "hNchFT0MNAssocMCCollisionsSameType", {HistType::kTH3D, {nChargedFT0MGenAxis, nAssocCollAxis, eventTypeAxis}});
    registry.add("hNContributorsCorrelation", "hNContributorsCorrelation", {HistType::kTH2F, {{250, -0.5f, 249.5f, "Secondary Contributor"}, {250, -0.5f, 249.5f, "Main Contributor"}}});
    registry.add("hNchFT0MGenEvType", "hNchFT0MGenEvType", {HistType::kTH2D, {nChargedFT0MGenAxis, eventTypeAxis}});
    registry.add("hCentFT0M_genMC", "hCentFT0M_genMC", {HistType::kTH2D, {centFT0MAxis, eventTypeAxis}});

    registry.add("hCentFT0M_rec", "hCentFT0M_rec", {HistType::kTH2D, {centFT0MAxis, eventTypeAxis}});
    registry.add("hCentFT0M_corr", "hCentFT0M_Corr", {HistType::kTH2D, {centFT0MAxis, centFT0MAxis}});

    if (multQA) {
      registry.add("hNchFT0Mglobal", "hNchFT0Mglobal", {HistType::kTH3D, {nChargedFT0MGenAxis, multNTracksAxis, eventTypeAxis}});
      registry.add("hNchFT0MPVContr", "hNchFT0MPVContr", {HistType::kTH3D, {nChargedFT0MGenAxis, multNTracksAxis, eventTypeAxis}});
      registry.add("hFT0MpvContr", "hFT0MpvContr", {HistType::kTH3D, {centFT0MAxis, multNTracksAxis, eventTypeAxis}});
      registry.add("hFT0Mglobal", "hFT0Mglobal", {HistType::kTH3D, {centFT0MAxis, multNTracksAxis, eventTypeAxis}});
      registry.add("hFT0MsignalPVContr", "hFT0MsignalPVContr", {HistType::kTH3D, {signalFT0MAxis, multNTracksAxis, eventTypeAxis}});
    }

    registry.add("h3ResonanceTruth", "pT distribution of True Resonance", kTHnSparseF, {eventTypeAxis, ptAxis, centFT0MAxis});
    registry.add("h3ResonanceTruthAnti", "pT distribution of True Resonance Anti", kTHnSparseF, {eventTypeAxis, ptAxis, centFT0MAxis});
  }
  float pvEta1 = 1.0f;
  float globalEta05 = 0.5f;

  Partition<DauTracks> pvContribTracksIUEta1 = (nabs(aod::track::eta) < pvEta1) && ((aod::track::flags & static_cast<uint32_t>(o2::aod::track::PVContributor)) == static_cast<uint32_t>(o2::aod::track::PVContributor));
  Partition<DauTracks> globalTracksIUEta05 = (nabs(aod::track::eta) < globalEta05) && (requireGlobalTrackInFilter());

  template <typename TMcParticles>
  uint16_t getGenNchInFT0Mregion(TMcParticles particles)
  {
    float region1FT0 = -3.3f;
    float region2FT0 = -2.1f;
    float region3FT0 = 3.5f;
    float region4FT0 = 4.9f;
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
      if (mcParticle.eta() < region1FT0 || mcParticle.eta() > region4FT0 || (mcParticle.eta() > region2FT0 && mcParticle.eta() < region3FT0)) {
        continue; // select on T0M Nch region
      }
      nchFT0++; // increment
    }
    return nchFT0;
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
    if (isZvtxcut && std::fabs(collision.posZ()) > cutzvertex) {
      return false;
    }
    if (isFillEventSelectionQA) {
      registry.fill(HIST("hNEvents"), 9.5);
      registry.fill(HIST("hZCollision"), collision.posZ());
    }

    return true;
  }

  template <typename TotalMCParts, typename MultMCGen, typename evtType>
  void fillMCParticles(TotalMCParts const& mcParticles, MultMCGen const& multiplicity, evtType const& eventType)
  {
    for (auto const& mcPart : mcParticles) {

      if (std::abs(mcPart.pdgCode()) != pdgTruthMother || std::abs(mcPart.y()) >= cfgRapidityCut)
        continue;
      std::vector<int> daughterPDGs;
      if (mcPart.has_daughters()) {
        auto daughter01 = mcParticles.rawIteratorAt(mcPart.daughtersIds()[0] - mcParticles.offset());
        auto daughter02 = mcParticles.rawIteratorAt(mcPart.daughtersIds()[1] - mcParticles.offset());
        daughterPDGs = {daughter01.pdgCode(), daughter02.pdgCode()};
      } else {
        daughterPDGs = {-1, -1};
      }

      if (isDaughterCheck) {
        bool pass1 = std::abs(daughterPDGs[0]) == pdgTruthDaughter1 || std::abs(daughterPDGs[1]) == pdgTruthDaughter1;
        bool pass2 = std::abs(daughterPDGs[0]) == pdgTruthDaughter2 || std::abs(daughterPDGs[1]) == pdgTruthDaughter2;
        if (!pass1 || !pass2)
          continue;
      }
      if (mcPart.pdgCode() > 0) // Consider INELt0 or INEL
        registry.fill(HIST("h3ResonanceTruth"), eventType, mcPart.pt(), multiplicity);
      else
        registry.fill(HIST("h3ResonanceTruthAnti"), eventType, mcPart.pt(), multiplicity);

      daughterPDGs.clear();
    }
  }
  void processData(soa::Join<aod::Collisions, aod::EvSels,
                             aod::PVMults, aod::FT0Mults, aod::FV0Mults,
                             aod::CentFT0Ms, aod::CentFV0As>::iterator const& collision,
                   DauTracks const&)
  {
    if (!acceptEvent(collision, 1)) {
      return;
    }

    int evType = 0;
    registry.fill(HIST("hNEvents"), 10.5); // INEL
    if (collision.isInelGt0()) {
      evType += 1;
      registry.fill(HIST("hNEvents"), 11.5); // INEL>0
    }

    auto tracksGroupedPVcontr = pvContribTracksIUEta1->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    int nTracksPVcontr = tracksGroupedPVcontr.size();

    auto tracksGroupedGlobal = globalTracksIUEta05->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    int nTracksGlobal = tracksGroupedGlobal.size();

    registry.fill(HIST("hCentFT0M_rec"), collision.centFT0M(), evType);

    if (multQA) {
      registry.fill(HIST("hFT0MpvContr"), collision.centFT0M(), nTracksPVcontr, evType);
      registry.fill(HIST("hFT0Mglobal"), collision.centFT0M(), nTracksGlobal, evType);
      registry.fill(HIST("hFT0MsignalPVContr"), collision.multFT0A() + collision.multFT0C(), nTracksPVcontr, evType);
    }
  }

  Preslice<aod::McParticles> perMcCollision = aod::mcparticle::mcCollisionId;
  void processMCrec(soa::Join<aod::Collisions, o2::aod::McCollisionLabels, aod::EvSels,
                              aod::PVMults, aod::FT0Mults, aod::FV0Mults,
                              aod::CentFT0Ms, aod::CentFV0As>::iterator const& collision,
                    soa::Join<aod::McCollisions, aod::McCentFT0Ms> const&,
                    DauTracks const&,
                    aod::McParticles const& mcParticles)
  {
    if (!acceptEvent(collision, 1)) {
      return;
    }

    if (!collision.has_mcCollision()) {
      registry.fill(HIST("fakeEvents"), 0); // no assoc. MC collisions
      return;
    }

    const auto& mcCollision = collision.mcCollision_as<soa::Join<aod::McCollisions, aod::McCentFT0Ms>>();

    auto tracksGroupedPVcontr = pvContribTracksIUEta1->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    int nTracksPVcontr = tracksGroupedPVcontr.size();

    auto tracksGroupedGlobal = globalTracksIUEta05->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    int nTracksGlobal = tracksGroupedGlobal.size();

    // N charged in FT0M region in corresponding gen. MC collision
    auto mcPartSlice = mcParticles.sliceBy(perMcCollision, collision.mcCollision_as<soa::Join<aod::McCollisions, aod::McCentFT0Ms>>().globalIndex());
    uint16_t nchFT0 = getGenNchInFT0Mregion(mcPartSlice);

    int evType = 0;
    registry.fill(HIST("hNEvents"), 10.5); // reco INEL
    if (collision.isInelGt0()) {
      evType += 1;
      registry.fill(HIST("hNEvents"), 11.5); // reco INEL>0
    }

    registry.fill(HIST("hCentFT0M_rec"), mcCollision.centFT0M(), evType); // correction only reco level in this stage
    registry.fill(HIST("hCentFT0M_corr"), mcCollision.centFT0M(), mcCollision.centFT0M(), evType);

    if (multQA) {
      registry.fill(HIST("hNchFT0MPVContr"), nchFT0, nTracksPVcontr, evType);
      registry.fill(HIST("hFT0MpvContr"), mcCollision.centFT0M(), nTracksPVcontr, evType);
      registry.fill(HIST("hFT0Mglobal"), mcCollision.centFT0M(), nTracksGlobal, evType);
      registry.fill(HIST("hNchFT0Mglobal"), nchFT0, nTracksGlobal, evType);
      registry.fill(HIST("hFT0MsignalPVContr"), collision.multFT0A() + collision.multFT0C(), nTracksPVcontr, evType);
    }
  }

  void processMCgen(soa::Join<aod::McCollisions, aod::McCentFT0Ms>::iterator const& mcCollision,
                    aod::McParticles const& mcParticles,
                    const soa::SmallGroups<o2::soa::Join<o2::aod::Collisions, o2::aod::McCollisionLabels,
                                                         o2::aod::EvSels, aod::PVMults, aod::FV0Mults, aod::FT0Mults, aod::CentFT0Ms, aod::CentFV0As>>& collisions)
  {
    auto cent = mcCollision.centFT0M();

    registry.fill(HIST("hNEventsMC"), 0.5);

    if (isZvtxcut && std::fabs(mcCollision.posZ()) > cutzvertex) {
      return;
    }
    registry.fill(HIST("hZCollisionGen"), mcCollision.posZ());
    registry.fill(HIST("hNEventsMC"), 1.5);

    int evType = 0;
    registry.fill(HIST("hNEventsMC"), 2.5);
    if (pwglf::isINELgtNmc(mcParticles, 0, pdgDB)) { // Truth INEL>0
      evType++;
      registry.fill(HIST("hNEventsMC"), 3.5);
    }

    fillMCParticles(mcParticles, cent, evType);

    registry.fill(HIST("hCentFT0M_genMC"), cent, evType);

    uint16_t nchFT0 = getGenNchInFT0Mregion(mcParticles);
    registry.fill(HIST("hNchFT0MGenEvType"), nchFT0, evType);

    std::vector<CollisionIndexAndType> selectedEvents(collisions.size());
    std::vector<int64_t> numberOfContributors;
    int nevts = 0;
    int nAssocColl = 0;
    const int nContSize = 2;
    for (const auto& collision : collisions) {
      CollisionIndexAndType collWithType = {0, 0x0};
      if (!acceptEvent(collision, 0)) {
        continue;
      }
      collWithType.index = collision.mcCollision_as<soa::Join<aod::McCollisions, aod::McCentFT0Ms>>().globalIndex();
      collWithType.typeFlag |= o2::aod::myMCcascades::EvFlags::EvINEL;

      if (collision.isInelGt0()) { // reco INEL>0
        collWithType.typeFlag |= o2::aod::myMCcascades::EvFlags::EvINELgt0;
      }
      selectedEvents[nevts++] = collWithType;
      if (collision.mcCollision_as<soa::Join<aod::McCollisions, aod::McCentFT0Ms>>().globalIndex() == mcCollision.globalIndex()) {
        nAssocColl++;
        numberOfContributors.push_back(collision.numContrib());
      }
    }
    selectedEvents.resize(nevts);

    registry.fill(HIST("hCentFT0MNAssocMCCollisions"), cent, nAssocColl, evType);
    registry.fill(HIST("hNchFT0MNAssocMCCollisions"), nchFT0, nAssocColl, evType);

    if (numberOfContributors.size() == nContSize) {
      std::sort(numberOfContributors.begin(), numberOfContributors.end());
      registry.fill(HIST("hNContributorsCorrelation"), numberOfContributors[0], numberOfContributors[1]);
    }

    auto isAssocToINEL = [&mcCollision](CollisionIndexAndType i) { return (i.index == mcCollision.globalIndex()) && ((i.typeFlag & o2::aod::myMCcascades::EvFlags::EvINEL) == o2::aod::myMCcascades::EvFlags::EvINEL); };
    auto isAssocToINELgt0 = [&mcCollision](CollisionIndexAndType i) { return (i.index == mcCollision.globalIndex()) && ((i.typeFlag & o2::aod::myMCcascades::EvFlags::EvINELgt0) == o2::aod::myMCcascades::EvFlags::EvINELgt0); };
    // number of reconstructed INEL events that have the same global index as mcCollision
    const auto evtReconstructedAndINEL = std::count_if(selectedEvents.begin(), selectedEvents.end(), isAssocToINEL);
    // number of reconstructed INEL > 0 events that have the same global index as mcCollision
    const auto evtReconstructedAndINELgt0 = std::count_if(selectedEvents.begin(), selectedEvents.end(), isAssocToINELgt0);
    switch (evType) {
      case 0: {
        registry.fill(HIST("hCentFT0MNAssocMCCollisionsSameType"), cent, evtReconstructedAndINEL, evType);
        registry.fill(HIST("hNchFT0MNAssocMCCollisionsSameType"), nchFT0, evtReconstructedAndINEL, evType);
        break;
      }
      case 1: {
        registry.fill(HIST("hCentFT0MNAssocMCCollisionsSameType"), cent, evtReconstructedAndINELgt0, evType);
        registry.fill(HIST("hNchFT0MNAssocMCCollisionsSameType"), nchFT0, evtReconstructedAndINELgt0, evType);
        break;
      }
      default:
        LOGF(fatal, "incorrect evType in event task");
        break;
    }

    if (evtReconstructedAndINELgt0) { // N INEL>0 reconstructed events associated with the MC collision
      registry.fill(HIST("hNEventsMC"), 4.5);
    }
  }
  PROCESS_SWITCH(Initializereventqa, processData, "Process Run 3 data", false);
  PROCESS_SWITCH(Initializereventqa, processMCrec, "Process Run 3 mc, Reconstructed", true);
  PROCESS_SWITCH(Initializereventqa, processMCgen, "Process Run 3 mc, genereated", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<Initializereventqa>(cfgc),
  };
}
