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
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include <CCDB/BasicCCDBManager.h>
#include "Framework/GroupSlicer.h"
#include "Framework/StepTHn.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "CommonConstants/MathConstants.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Centrality.h"
#include "PWGCF/DataModel/CorrelationsDerived.h"
#include "PWGCF/Core/CorrelationContainer.h"
#include "PWGCF/Core/PairCuts.h"
#include "DataFormatsParameters/GRPObject.h"

#include <TH1F.h>
#include <cmath>
#include <TDirectory.h>
#include <THn.h>

namespace o2::aod
{
namespace hash
{
DECLARE_SOA_COLUMN(Bin, bin, int);
} // namespace hash
DECLARE_SOA_TABLE(Hashes, "AOD", "HASH", hash::Bin);

using Hash = Hashes::iterator;
} // namespace o2::aod

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace constants::math;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

static constexpr float cfgPairCutDefaults[1][5] = {{-1, -1, -1, -1, -1}};

struct CorrelationTask {

  // Configuration
  O2_DEFINE_CONFIGURABLE(cfgCutVertex, float, 7.0f, "Accepted z-vertex range")
  O2_DEFINE_CONFIGURABLE(cfgCutPt, float, 0.5f, "Minimal pT for tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutEta, float, 0.8f, "Eta range for tracks")

  O2_DEFINE_CONFIGURABLE(cfgPtOrder, int, 1, "Only consider pairs for which pT,1 < pT,2 (0 = OFF, 1 = ON)");
  O2_DEFINE_CONFIGURABLE(cfgTriggerCharge, int, 0, "Select on charge of trigger particle: 0 = all; 1 = positive; -1 = negative");
  O2_DEFINE_CONFIGURABLE(cfgAssociatedCharge, int, 0, "Select on charge of associated particle: 0 = all; 1 = positive; -1 = negative");
  O2_DEFINE_CONFIGURABLE(cfgPairCharge, int, 0, "Select on charge of particle pair: 0 = all; 1 = like sign; -1 = unlike sign");

  O2_DEFINE_CONFIGURABLE(cfgTwoTrackCut, float, -1, "Two track cut: -1 = off; >0 otherwise distance value (suggested: 0.02)");
  O2_DEFINE_CONFIGURABLE(cfgTwoTrackCutMinRadius, float, 0.8f, "Two track cut: radius in m from which two track cuts are applied");

  // Suggested values: Photon: 0.004; K0 and Lambda: 0.005
  Configurable<LabeledArray<float>> cfgPairCut{"cfgPairCut", {cfgPairCutDefaults[0], 5, {"Photon", "K0", "Lambda", "Phi", "Rho"}}, "Pair cuts on various particles"};

  O2_DEFINE_CONFIGURABLE(cfgEfficiencyTrigger, std::string, "", "CCDB path to efficiency object for trigger particles")
  O2_DEFINE_CONFIGURABLE(cfgEfficiencyAssociated, std::string, "", "CCDB path to efficiency object for associated particles")

  O2_DEFINE_CONFIGURABLE(cfgNoMixedEvents, int, 5, "Number of mixed events per event")

  ConfigurableAxis axisVertex{"axisVertex", {7, -7, 7}, "vertex axis for histograms"};
  ConfigurableAxis axisDeltaPhi{"axisDeltaPhi", {72, -PIHalf, PIHalf * 3}, "delta phi axis for histograms"};
  ConfigurableAxis axisDeltaEta{"axisDeltaEta", {40, -2, 2}, "delta eta axis for histograms"};
  ConfigurableAxis axisPtTrigger{"axisPtTrigger", {VARIABLE_WIDTH, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 10.0}, "pt trigger axis for histograms"};
  ConfigurableAxis axisPtAssoc{"axisPtAssoc", {VARIABLE_WIDTH, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0}, "pt associated axis for histograms"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 100.1}, "multiplicity / centrality axis for histograms"};

  ConfigurableAxis axisVertexEfficiency{"axisVertexEfficiency", {10, -10, 10}, "vertex axis for efficiency histograms"};
  ConfigurableAxis axisEtaEfficiency{"axisEtaEfficiency", {20, -1.0, 1.0}, "eta axis for efficiency histograms"};
  ConfigurableAxis axisPtEfficiency{"axisPtEfficiency", {VARIABLE_WIDTH, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0}, "pt axis for efficiency histograms"};

  // This filter is applied to AOD and derived data (column names are identical)
  Filter collisionZVtxFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  // This filter is only applied to AOD
  Filter collisionVertexTypeFilter = (aod::collision::flags & (uint16_t)aod::collision::CollisionFlagsRun2::Run2VertexerTracks) == (uint16_t)aod::collision::CollisionFlagsRun2::Run2VertexerTracks;

  // Track filters
  Filter trackFilter = (nabs(aod::track::eta) < cfgCutEta) && (aod::track::pt > cfgCutPt) && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true));
  Filter cfTrackFilter = (nabs(aod::cftrack::eta) < cfgCutEta) && (aod::cftrack::pt > cfgCutPt);

  using aodTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection>>;
  using derivedTracks = soa::Filtered<aod::CFTracks>;

  // Output definitions
  OutputObj<CorrelationContainer> same{"sameEvent"};
  OutputObj<CorrelationContainer> mixed{"mixedEvent"};

  struct Config {
    bool mPairCuts = false;
    THn* mEfficiencyTrigger = nullptr;
    THn* mEfficiencyAssociated = nullptr;
  } cfg;

  HistogramRegistry registry{"registry"};
  PairCuts mPairCuts;

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  using aodCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms>>;
  using derivedCollisions = soa::Filtered<aod::CFCollisions>;
  using BinningType = BinningPolicy<aod::collision::PosZ, aod::cent::CentRun2V0M>;

  void init(o2::framework::InitContext&)
  {
    registry.add("yields", "centrality vs pT vs eta", {HistType::kTH3F, {{100, 0, 100, "centrality"}, {40, 0, 20, "p_{T}"}, {100, -2, 2, "#eta"}}});
    registry.add("etaphi", "centrality vs eta vs phi", {HistType::kTH3F, {{100, 0, 100, "centrality"}, {100, -2, 2, "#eta"}, {200, 0, 2 * M_PI, "#varphi"}}});

    const int maxMixBin = axisMultiplicity->size() * axisVertex->size();
    registry.add("eventcount", "bin", {HistType::kTH1F, {{maxMixBin + 2, -2.5, -0.5 + maxMixBin, "bin"}}});

    mPairCuts.SetHistogramRegistry(&registry);

    if (cfgPairCut->get("Photon") > 0 || cfgPairCut->get("K0") > 0 || cfgPairCut->get("Lambda") > 0 || cfgPairCut->get("Phi") > 0 || cfgPairCut->get("Rho") > 0) {
      mPairCuts.SetPairCut(PairCuts::Photon, cfgPairCut->get("Photon"));
      mPairCuts.SetPairCut(PairCuts::K0, cfgPairCut->get("K0"));
      mPairCuts.SetPairCut(PairCuts::Lambda, cfgPairCut->get("Lambda"));
      mPairCuts.SetPairCut(PairCuts::Phi, cfgPairCut->get("Phi"));
      mPairCuts.SetPairCut(PairCuts::Rho, cfgPairCut->get("Rho"));
      cfg.mPairCuts = true;
    }

    if (cfgTwoTrackCut > 0) {
      mPairCuts.SetTwoTrackCuts(cfgTwoTrackCut, cfgTwoTrackCutMinRadius);
    }

    // --- OBJECT INIT ---

    std::vector<AxisSpec> axisList = {{axisDeltaEta, "#Delta#eta"},
                                      {axisPtAssoc, "p_{T} (GeV/c)"},
                                      {axisPtTrigger, "p_{T} (GeV/c)"},
                                      {axisMultiplicity, "multiplicity / centrality"},
                                      {axisDeltaPhi, "#Delta#varphi (rad)"},
                                      {axisVertex, "z-vtx (cm)"},
                                      {axisEtaEfficiency, "#eta"},
                                      {axisPtEfficiency, "p_{T} (GeV/c)"},
                                      {axisVertexEfficiency, "z-vtx (cm)"}};
    same.setObject(new CorrelationContainer("sameEvent", "sameEvent", axisList));
    mixed.setObject(new CorrelationContainer("mixedEvent", "mixedEvent", axisList));

    same->setTrackEtaCut(cfgCutEta);
    mixed->setTrackEtaCut(cfgCutEta);

    // o2-ccdb-upload -p Users/jgrosseo/correlations/LHC15o -f /tmp/correction_2011_global.root -k correction

    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    long now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now); // TODO must become global parameter from the train creation time
  }

  int getMagneticField(uint64_t timestamp)
  {
    // TODO done only once (and not per run). Will be replaced by CCDBConfigurable
    static o2::parameters::GRPObject* grpo = nullptr;
    if (grpo == nullptr) {
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>("GLO/GRP/GRP", timestamp);
      if (grpo == nullptr) {
        LOGF(fatal, "GRP object not found for timestamp %llu", timestamp);
        return 0;
      }
      LOGF(info, "Retrieved GRP for timestamp %llu with magnetic field of %d kG", timestamp, grpo->getNominalL3Field());
    }
    return grpo->getNominalL3Field();
  }

  template <typename TCollision, typename TTracks>
  void fillQA(TCollision collision, float centrality, TTracks tracks)
  {
    for (auto& track1 : tracks) {
      registry.fill(HIST("yields"), centrality, track1.pt(), track1.eta());
      registry.fill(HIST("etaphi"), centrality, track1.eta(), track1.phi());
    }
  }

  template <typename TTarget, typename TCollision>
  bool fillCollisionAOD(TTarget target, TCollision collision, float centrality)
  {
    target->fillEvent(centrality, CorrelationContainer::kCFStepAll);

    if (!collision.alias()[kINT7] || !collision.sel7()) {
      return false;
    }

    target->fillEvent(centrality, CorrelationContainer::kCFStepReconstructed);

    return true;
  }

  template <typename TTarget, typename TTracks>
  void fillCorrelations(TTarget target, TTracks tracks1, TTracks tracks2, float centrality, float posZ, int magField)
  {
    // Cache efficiency for particles (too many FindBin lookups)
    float* efficiencyAssociated = nullptr;
    if (cfg.mEfficiencyAssociated) {
      efficiencyAssociated = new float[tracks2.size()];
      int i = 0;
      for (auto& track : tracks2) {
        efficiencyAssociated[i++] = getEfficiency(cfg.mEfficiencyAssociated, track.eta(), track.pt(), centrality, posZ);
      }
    }

    for (auto& track1 : tracks1) {
      // LOGF(info, "Track %f | %f | %f  %d %d", track1.eta(), track1.phi(), track1.pt(), track1.isGlobalTrack(), track1.isGlobalTrackSDD());

      if (cfgTriggerCharge != 0 && cfgTriggerCharge * track1.sign() < 0) {
        continue;
      }

      float triggerWeight = 1.0;
      if (cfg.mEfficiencyTrigger) {
        triggerWeight = getEfficiency(cfg.mEfficiencyTrigger, track1.eta(), track1.pt(), centrality, posZ);
      }

      target->getTriggerHist()->Fill(CorrelationContainer::kCFStepReconstructed, track1.pt(), centrality, posZ, triggerWeight);

      int i = -1;
      for (auto& track2 : tracks2) {
        i++; // HACK
        if (track1 == track2) {
          continue;
        }

        if (cfgPtOrder != 0 && track2.pt() >= track1.pt()) {
          continue;
        }

        if (cfgAssociatedCharge != 0 && cfgAssociatedCharge * track2.sign() < 0) {
          continue;
        }
        if (cfgPairCharge != 0 && cfgPairCharge * track1.sign() * track2.sign() < 0) {
          continue;
        }

        if (cfg.mPairCuts && mPairCuts.conversionCuts(track1, track2)) {
          continue;
        }

        if (cfgTwoTrackCut > 0 && mPairCuts.twoTrackCut(track1, track2, magField)) {
          continue;
        }

        float associatedWeight = 1.0;
        if (cfg.mEfficiencyAssociated) {
          associatedWeight = efficiencyAssociated[i];
        }

        float deltaPhi = track1.phi() - track2.phi();
        if (deltaPhi > 1.5f * PI) {
          deltaPhi -= TwoPI;
        }
        if (deltaPhi < -PIHalf) {
          deltaPhi += TwoPI;
        }

        target->getPairHist()->Fill(CorrelationContainer::kCFStepReconstructed,
                                    track1.eta() - track2.eta(), track2.pt(), track1.pt(), centrality, deltaPhi, posZ,
                                    triggerWeight * associatedWeight);
      }
    }

    delete[] efficiencyAssociated;
  }

  // Version with explicit nested loop
  void processSameAOD(aodCollisions::iterator const& collision, aod::BCsWithTimestamps const&, aodTracks const& tracks)
  {
    // TODO will go to CCDBConfigurable
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();

    if (cfgEfficiencyTrigger.value.empty() == false) {
      cfg.mEfficiencyTrigger = ccdb->getForTimeStamp<THnT<float>>(cfgEfficiencyTrigger, bc.timestamp());
      if (cfg.mEfficiencyTrigger == nullptr) {
        LOGF(fatal, "Could not load efficiency histogram for trigger particles from %s", cfgEfficiencyTrigger.value.c_str());
      }
      LOGF(info, "Loaded efficiency histogram for trigger particles from %s (%p)", cfgEfficiencyTrigger.value.c_str(), (void*)cfg.mEfficiencyTrigger);
    }
    if (cfgEfficiencyAssociated.value.empty() == false) {
      cfg.mEfficiencyAssociated = ccdb->getForTimeStamp<THnT<float>>(cfgEfficiencyAssociated, bc.timestamp());
      if (cfg.mEfficiencyAssociated == nullptr) {
        LOGF(fatal, "Could not load efficiency histogram for associated particles from %s", cfgEfficiencyAssociated.value.c_str());
      }
      LOGF(info, "Loaded efficiency histogram for associated particles from %s (%p)", cfgEfficiencyAssociated.value.c_str(), (void*)cfg.mEfficiencyAssociated);
    }

    LOGF(info, "processSameAOD: Tracks for collision: %d | Vertex: %.1f | INT7: %d | V0M: %.1f", tracks.size(), collision.posZ(), collision.sel7(), collision.centRun2V0M());

    const auto centrality = collision.centRun2V0M();

    if (fillCollisionAOD(same, collision, centrality) == false) {
      return;
    }
    registry.fill(HIST("eventcount"), -2);
    fillQA(collision, centrality, tracks);
    fillCorrelations(same, tracks, tracks, centrality, collision.posZ(), getMagneticField(bc.timestamp()));
  }
  PROCESS_SWITCH(CorrelationTask, processSameAOD, "Process same event on AOD", true);

  void processSameDerived(derivedCollisions::iterator const& collision, soa::Filtered<aod::CFTracks> const& tracks)
  {
    LOGF(info, "processSameDerived: Tracks for collision: %d | Vertex: %.1f | V0M: %.1f", tracks.size(), collision.posZ(), collision.centRun2V0M());

    const auto centrality = collision.centRun2V0M();

    same->fillEvent(centrality, CorrelationContainer::kCFStepReconstructed);
    registry.fill(HIST("eventcount"), -2);
    fillQA(collision, centrality, tracks);
    fillCorrelations(same, tracks, tracks, centrality, collision.posZ(), getMagneticField(collision.timestamp()));
  }
  PROCESS_SWITCH(CorrelationTask, processSameDerived, "Process same event on derived data", false);

  NoBinningPolicy<aod::hash::Bin> hashBin;
  void processMixedAOD(soa::Filtered<soa::Join<aod::Collisions, aod::Hashes, aod::EvSels, aod::CentRun2V0Ms>>& collisions, aodTracks const& tracks, aod::BCsWithTimestamps const&)
  {
    // TODO loading of efficiency histogram missing here, because it will happen somehow in the CCDBConfigurable

    collisions.bindExternalIndices(&tracks);
    auto tracksTuple = std::make_tuple(tracks);
    GroupSlicer slicer(collisions, tracksTuple);

    // Strictly upper categorised collisions, for cfgNoMixedEvents combinations per bin, skipping those in entry -1
    for (auto& [collision1, collision2] : selfCombinations(hashBin, cfgNoMixedEvents, -1, collisions, collisions)) {

      LOGF(info, "processMixedAOD: Mixed collisions bin: %d pair: %d (%.3f, %.3f), %d (%.3f, %.3f)", collision1.bin(), collision1.globalIndex(), collision1.posZ(), collision1.centRun2V0M(), collision2.globalIndex(), collision2.posZ(), collision2.centRun2V0M());

      // TODO in principle these should be already checked on hash level, because in this way we don't check collision 2
      // TODO not correct because event-level histograms on collision1 are filled for each pair (important :))
      if (fillCollisionAOD(mixed, collision1, collision1.centRun2V0M()) == false) {
        continue;
      }
      registry.fill(HIST("eventcount"), collision1.bin());

      auto it1 = slicer.begin();
      auto it2 = slicer.begin();
      for (auto& slice : slicer) {
        if (slice.groupingElement().index() == collision1.index()) {
          it1 = slice;
          break;
        }
      }
      for (auto& slice : slicer) {
        if (slice.groupingElement().index() == collision2.index()) {
          it2 = slice;
          break;
        }
      }

      auto tracks1 = std::get<aodTracks>(it1.associatedTables());
      tracks1.bindExternalIndices(&collisions);
      auto tracks2 = std::get<aodTracks>(it2.associatedTables());
      tracks2.bindExternalIndices(&collisions);

      auto bc = collision1.bc_as<aod::BCsWithTimestamps>();

      // LOGF(info, "Tracks: %d and %d entries", tracks1.size(), tracks2.size());

      // TODO mixed event weight missing
      fillCorrelations(mixed, tracks1, tracks2, collision1.centRun2V0M(), collision1.posZ(), getMagneticField(bc.timestamp()));
    }
  }
  PROCESS_SWITCH(CorrelationTask, processMixedAOD, "Process mixed events on AOD", true);

  void processMixedDerived(soa::Filtered<soa::Join<aod::CFCollisions, aod::Hashes>>& collisions, derivedTracks const& tracks)
  {
    // TODO loading of efficiency histogram missing here, because it will happen somehow in the CCDBConfigurable

    collisions.bindExternalIndices(&tracks);
    auto tracksTuple = std::make_tuple(tracks);
    GroupSlicer slicer(collisions, tracksTuple);

    // Strictly upper categorised collisions, for cfgNoMixedEvents combinations per bin, skipping those in entry -1
    for (auto& [collision1, collision2] : selfCombinations(hashBin, cfgNoMixedEvents, -1, collisions, collisions)) {

      LOGF(info, "processMixedDerived: Mixed collisions bin: %d pair: %d (%.3f, %.3f), %d (%.3f, %.3f)", collision1.bin(), collision1.globalIndex(), collision1.posZ(), collision1.centRun2V0M(), collision2.globalIndex(), collision2.posZ(), collision2.centRun2V0M());

      registry.fill(HIST("eventcount"), collision1.bin());
      mixed->fillEvent(collision1.centRun2V0M(), CorrelationContainer::kCFStepReconstructed);

      auto it1 = slicer.begin();
      auto it2 = slicer.begin();
      for (auto& slice : slicer) {
        if (slice.groupingElement().index() == collision1.index()) {
          it1 = slice;
          break;
        }
      }
      for (auto& slice : slicer) {
        if (slice.groupingElement().index() == collision2.index()) {
          it2 = slice;
          break;
        }
      }

      auto tracks1 = std::get<derivedTracks>(it1.associatedTables());
      tracks1.bindExternalIndices(&collisions);
      auto tracks2 = std::get<derivedTracks>(it2.associatedTables());
      tracks2.bindExternalIndices(&collisions);

      // LOGF(info, "Tracks: %d and %d entries", tracks1.size(), tracks2.size());

      fillCorrelations(mixed, tracks1, tracks2, collision1.centRun2V0M(), collision1.posZ(), getMagneticField(collision1.timestamp()));
    }
  }
  PROCESS_SWITCH(CorrelationTask, processMixedDerived, "Process mixed events on derived data", false);

  BinningType pairBinning{{axisVertex, axisMultiplicity}, true};                               // true is for 'ignore overflows' (true by default). Underflows and overflows will have bin -1.
  SameKindPair<aodCollisions, aodTracks, BinningType> pair{pairBinning, cfgNoMixedEvents, -1}; // -1 is the number of the bin to skip
  void processMixedAODNoHash(aodCollisions& collisions, aodTracks const& tracks, aod::BCsWithTimestamps const&)
  {
    // TODO loading of efficiency histogram missing here, because it will happen somehow in the CCDBConfigurable

    // Strictly upper categorised collisions, for cfgNoMixedEvents combinations per bin, skipping those in entry -1
    for (auto& [collision1, tracks1, collision2, tracks2] : pair) {
      int bin = pairBinning.getBin({collision1.posZ(), collision1.centRun2V0M()});
      LOGF(info, "processMixedAODNoHash: Mixed collisions bin: %d pair: %d (%.3f, %.3f), %d (%.3f, %.3f)", bin, collision1.globalIndex(), collision1.posZ(), collision1.centRun2V0M(), collision2.globalIndex(), collision2.posZ(), collision2.centRun2V0M());

      // TODO in principle these should be already checked on hash level, because in this way we don't check collision 2
      // TODO not correct because event-level histograms on collision1 are filled for each pair (important :))
      if (fillCollisionAOD(mixed, collision1, collision1.centRun2V0M()) == false) {
        continue;
      }
      registry.fill(HIST("eventcount"), bin);

      auto bc = collision1.bc_as<aod::BCsWithTimestamps>();

      // LOGF(info, "Tracks: %d and %d entries", tracks1.size(), tracks2.size());

      // TODO mixed event weight missing
      fillCorrelations(mixed, tracks1, tracks2, collision1.centRun2V0M(), collision1.posZ(), getMagneticField(bc.timestamp()));
    }
  }
  PROCESS_SWITCH(CorrelationTask, processMixedAODNoHash, "Process mixed events on AOD", false);

  void processMixedAODNoHashInProcess(aodCollisions& collisions, aodTracks const& tracks, aod::BCsWithTimestamps const&)
  {
    // TODO loading of efficiency histogram missing here, because it will happen somehow in the CCDBConfigurable

    // Strictly upper categorised collisions, for cfgNoMixedEvents combinations per bin, skipping those in entry -1
    BinningType configurableBinning{{axisVertex, axisMultiplicity}, true}; // true is for 'ignore overflows' (true by default). Underflows and overflows will have bin -1.
    auto tracksTuple = std::make_tuple(tracks);
    SameKindPair<aodCollisions, aodTracks, BinningType> pairInProcess{configurableBinning, cfgNoMixedEvents, -1, collisions, tracksTuple}; // -1 is the number of the bin to skip

    for (auto& [collision1, tracks1, collision2, tracks2] : pairInProcess) {
      int bin = configurableBinning.getBin({collision1.posZ(), collision1.centRun2V0M()});
      LOGF(info, "processMixedAODNoHashInProcess: Mixed collisions bin: %d pair: %d (%.3f, %.3f), %d (%.3f, %.3f)", bin, collision1.globalIndex(), collision1.posZ(), collision1.centRun2V0M(), collision2.globalIndex(), collision2.posZ(), collision2.centRun2V0M());

      // TODO in principle these should be already checked on hash level, because in this way we don't check collision 2
      // TODO not correct because event-level histograms on collision1 are filled for each pair (important :))
      if (fillCollisionAOD(mixed, collision1, collision1.centRun2V0M()) == false) {
        continue;
      }
      registry.fill(HIST("eventcount"), bin);

      auto bc = collision1.bc_as<aod::BCsWithTimestamps>();

      // LOGF(info, "Tracks: %d and %d entries", tracks1.size(), tracks2.size());

      // TODO mixed event weight missing
      fillCorrelations(mixed, tracks1, tracks2, collision1.centRun2V0M(), collision1.posZ(), getMagneticField(bc.timestamp()));
    }
  }
  PROCESS_SWITCH(CorrelationTask, processMixedAODNoHashInProcess, "Process mixed events on AOD", false);

  // Version with combinations
  void processWithCombinations(soa::Join<aod::Collisions, aod::CentRun2V0Ms>::iterator const& collision, aod::BCsWithTimestamps const&, soa::Filtered<aod::Tracks> const& tracks)
  {
    LOGF(info, "Tracks for collision (Combination run): %d", tracks.size());

    // TODO will go to CCDBConfigurable
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();

    const auto centrality = collision.centRun2V0M();

    for (auto track1 = tracks.begin(); track1 != tracks.end(); ++track1) {

      if (cfgTriggerCharge != 0 && cfgTriggerCharge * track1.sign() < 0) {
        continue;
      }

      //       LOGF(info, "TRACK %f %f | %f %f | %f %f", track1.eta(), track1.eta(), track1.phi(), track1.phi2(), track1.pt(), track1.pt());

      same->getTriggerHist()->Fill(CorrelationContainer::kCFStepReconstructed, track1.pt(), centrality, collision.posZ());
      // mixed->getTriggerHist()->Fill(eventValues, CorrelationContainer::kCFStepReconstructed);
    }

    for (auto& [track1, track2] : combinations(tracks, tracks)) {
      // LOGF(info, "Combination %d %d", track1.index(), track2.index());

      if (cfgTriggerCharge != 0 && cfgTriggerCharge * track1.sign() < 0) {
        continue;
      }
      if (cfgAssociatedCharge != 0 && cfgAssociatedCharge * track2.sign() < 0) {
        continue;
      }
      if (cfgPairCharge != 0 && cfgPairCharge * track1.sign() * track2.sign() < 0) {
        continue;
      }

      if (cfg.mPairCuts && mPairCuts.conversionCuts(track1, track2)) {
        continue;
      }

      if (cfgTwoTrackCut > 0 && mPairCuts.twoTrackCut(track1, track2, getMagneticField(bc.timestamp()))) {
        continue;
      }

      float deltaPhi = track1.phi() - track2.phi();
      if (deltaPhi > 1.5f * PI) {
        deltaPhi -= TwoPI;
      }
      if (deltaPhi < -PIHalf) {
        deltaPhi += TwoPI;
      }

      same->getPairHist()->Fill(CorrelationContainer::kCFStepReconstructed,
                                track1.eta() - track2.eta(), track2.pt(), track1.pt(), centrality, deltaPhi, collision.posZ());
      // mixed->getPairHist()->Fill(values, CorrelationContainer::kCFStepReconstructed);
    }
  }

  PROCESS_SWITCH(CorrelationTask, processWithCombinations, "Process same event on AOD with combinations", false);

  double getEfficiency(THn* eff, float eta, float pt, float centrality, float posZ)
  {
    int effVars[4];
    effVars[0] = eff->GetAxis(0)->FindBin(eta);
    effVars[1] = eff->GetAxis(1)->FindBin(pt);
    effVars[2] = eff->GetAxis(2)->FindBin(centrality);
    effVars[3] = eff->GetAxis(3)->FindBin(posZ);
    return eff->GetBinContent(effVars);
  }
};

struct CorrelationHashTask {
  std::vector<float> vtxBins;
  std::vector<float> multBins;

  Produces<aod::Hashes> hashes;

  void fillArray(int length, double* source, std::vector<float>& target)
  {
    // Expand binning from Configurable. Can we let some code in AxisSpec do this?

    target.clear();
    if (source[0] == VARIABLE_WIDTH) {
      for (int i = 1; i < length; i++) {
        target.push_back(source[i]);
      }
    } else {
      for (int i = 0; i <= source[0]; i++) {
        target.push_back(source[1] + (source[2] - source[1]) / source[0] * i);
      }
    }
  }

  void init(o2::framework::InitContext& initContext)
  {
    // get own suffix. Is there a better way?
    auto& deviceSpec = initContext.services().get<DeviceSpec const>();
    std::string suffix(deviceSpec.name);
    suffix.replace(0, strlen("correlation-hash-task"), "");

    // get axis config from CorrelationTask
    auto& workflows = initContext.services().get<RunningWorkflowInfo const>();
    for (DeviceSpec device : workflows.devices) {
      if (device.name == "correlation-task" + suffix) {
        for (auto option : device.options) {
          if (option.name == "axisVertex") {
            fillArray(option.defaultValue.size(), option.defaultValue.get<double*>(), vtxBins);
            LOGF(info, "Initialized vertex binning for mixing from configurable %s", option.name);
          }
          if (option.name == "axisMultiplicity") {
            fillArray(option.defaultValue.size(), option.defaultValue.get<double*>(), multBins);
            LOGF(info, "Initialized multiplicity binning for mixing from configurable %s", option.name);
          }
        }
      }
    }

    if (vtxBins.size() == 0) {
      LOGF(fatal, "vtxBins not configured. Check configuration.");
    }
    if (multBins.size() == 0) {
      LOGF(fatal, "multBins not configured. Check configuration.");
    }
  }

  // Calculate hash for an element based on 2 properties and their bins.
  int getHash(float vtx, float mult)
  {
    // underflow
    if (vtx < vtxBins[0]) {
      return -1;
    }
    if (mult < multBins[0]) {
      return -1;
    }

    for (unsigned int i = 1; i < vtxBins.size(); i++) {
      if (vtx < vtxBins[i]) {
        for (unsigned int j = 1; j < multBins.size(); j++) {
          if (mult < multBins[j]) {
            return i + j * (vtxBins.size() + 1);
          }
        }
      }
    }
    // overflow
    return -1;
  }

  void processAOD(soa::Join<aod::Collisions, aod::CentRun2V0Ms> const& collisions)
  {
    for (auto& collision : collisions) {
      int hash = getHash(collision.posZ(), collision.centRun2V0M());
      // LOGF(info, "Collision: %d (%f, %f) hash: %d", collision.index(), collision.posZ(), collision.centRun2V0M(), hash);
      hashes(hash);
    }
  }
  PROCESS_SWITCH(CorrelationHashTask, processAOD, "Create hashes for mixing on AOD", true);

  void processDerived(aod::CFCollisions const& collisions)
  {
    for (auto& collision : collisions) {
      int hash = getHash(collision.posZ(), collision.centRun2V0M());
      // LOGF(info, "Collision: %d (%f, %f) hash: %d", collision.index(), collision.posZ(), collision.centRun2V0M(), hash);
      hashes(hash);
    }
  }
  PROCESS_SWITCH(CorrelationHashTask, processDerived, "Create hashes for mixing on AOD", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<CorrelationHashTask>(cfgc),
    adaptAnalysisTask<CorrelationTask>(cfgc)};
}
