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
#include "DataFormatsParameters/GRPMagField.h"

#include <TH1F.h>
#include <cmath>
#include <TDirectory.h>
#include <THn.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace constants::math;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

// NOTE This is a nice idea but will again make it impossible to use subwagon configurations...
// namespace o2::aod
// {
// namespace cfcorr
// {
// DECLARE_SOA_COLUMN(Correction, correction, float);          //! Correction factor for this track
// } // namespace cfcorr
// DECLARE_SOA_TABLE(CFCorrections, "AOD", "CFCORRECTIONS", //! Table which attaches efficiency correction factor to tracks
//                   cfcorreff::Correction);
// } // namespace o2::aod

static constexpr float cfgPairCutDefaults[1][5] = {{-1, -1, -1, -1, -1}};

struct CorrelationTask {
  SliceCache cache;

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

  O2_DEFINE_CONFIGURABLE(cfgVerbosity, int, 1, "Verbosity level (0 = major, 1 = per collision)")

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

  // MC filters
  Filter cfMCCollisionFilter = nabs(aod::mccollision::posZ) < cfgCutVertex;
  Filter cfMCParticleFilter = (nabs(aod::cfmcparticle::eta) < cfgCutEta) && (aod::cfmcparticle::pt > cfgCutPt) && (aod::cfmcparticle::sign != 0);

  // Output definitions
  OutputObj<CorrelationContainer> same{"sameEvent"};
  OutputObj<CorrelationContainer> mixed{"mixedEvent"};

  struct Config {
    bool mPairCuts = false;
    THn* mEfficiencyTrigger = nullptr;
    THn* mEfficiencyAssociated = nullptr;
    bool efficiencyLoaded = false;
  } cfg;

  HistogramRegistry registry{"registry"};
  PairCuts mPairCuts;

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  using aodCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms>>;
  using aodTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection>>;

  using derivedTracks = soa::Filtered<aod::CFTracks>;
  using derivedCollisions = soa::Filtered<aod::CFCollisions>;

  void init(o2::framework::InitContext&)
  {
    registry.add("yields", "multiplicity/centrality vs pT vs eta", {HistType::kTH3F, {{100, 0, 100, "/multiplicity/centrality"}, {40, 0, 20, "p_{T}"}, {100, -2, 2, "#eta"}}});
    registry.add("etaphi", "multiplicity/centrality vs eta vs phi", {HistType::kTH3F, {{100, 0, 100, "multiplicity/centrality"}, {100, -2, 2, "#eta"}, {200, 0, 2 * M_PI, "#varphi"}}});

    const int maxMixBin = AxisSpec(axisMultiplicity).getNbins() * AxisSpec(axisVertex).getNbins();
    registry.add("eventcount_same", "bin", {HistType::kTH1F, {{maxMixBin + 2, -2.5, -0.5 + maxMixBin, "bin"}}});
    registry.add("eventcount_mixed", "bin", {HistType::kTH1F, {{maxMixBin + 2, -2.5, -0.5 + maxMixBin, "bin"}}});

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

    std::vector<AxisSpec> corrAxis = {{axisDeltaEta, "#Delta#eta"},
                                      {axisPtAssoc, "p_{T} (GeV/c)"},
                                      {axisPtTrigger, "p_{T} (GeV/c)"},
                                      {axisMultiplicity, "multiplicity / centrality"},
                                      {axisDeltaPhi, "#Delta#varphi (rad)"},
                                      {axisVertex, "z-vtx (cm)"}};
    std::vector<AxisSpec> effAxis = {{axisEtaEfficiency, "#eta"},
                                     {axisPtEfficiency, "p_{T} (GeV/c)"},
                                     {axisVertexEfficiency, "z-vtx (cm)"}};
    same.setObject(new CorrelationContainer("sameEvent", "sameEvent", corrAxis, effAxis, {}));
    mixed.setObject(new CorrelationContainer("mixedEvent", "mixedEvent", corrAxis, effAxis, {}));

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
    // static o2::parameters::GRPMagField* grpo = nullptr;
    if (grpo == nullptr) {
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>("GLO/GRP/GRP", timestamp);
      // grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", timestamp);
      if (grpo == nullptr) {
        LOGF(fatal, "GRP object not found for timestamp %llu", timestamp);
        return 0;
      }
      LOGF(info, "Retrieved GRP for timestamp %llu with magnetic field of %d kG", timestamp, grpo->getNominalL3Field());
    }
    return grpo->getNominalL3Field();
  }

  template <typename TCollision, typename TTracks>
  void fillQA(TCollision collision, float multiplicity, TTracks tracks)
  {
    for (auto& track1 : tracks) {
      registry.fill(HIST("yields"), multiplicity, track1.pt(), track1.eta());
      registry.fill(HIST("etaphi"), multiplicity, track1.eta(), track1.phi());
    }
  }

  template <typename TTarget, typename TCollision>
  bool fillCollisionAOD(TTarget target, TCollision collision, float multiplicity)
  {
    target->fillEvent(multiplicity, CorrelationContainer::kCFStepAll);

    if (!collision.alias_bit(kINT7) || !collision.sel7()) {
      return false;
    }

    target->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);

    return true;
  }

  template <CorrelationContainer::CFStep step, typename TTrack>
  bool checkObject(TTrack& track)
  {
    if constexpr (step <= CorrelationContainer::kCFStepAnaTopology) {
      return track.isPhysicalPrimary();
    } else if constexpr (step == CorrelationContainer::kCFStepTrackedOnlyPrim) {
      return track.isPhysicalPrimary() && (track.flags() & aod::cfmcparticle::kReconstructed);
    } else if constexpr (step == CorrelationContainer::kCFStepTracked) {
      return (track.flags() & aod::cfmcparticle::kReconstructed);
    }

    return true;
  }

  template <CorrelationContainer::CFStep step, typename TTarget, typename TTracks>
  void fillCorrelations(TTarget target, TTracks& tracks1, TTracks& tracks2, float multiplicity, float posZ, int magField, float eventWeight)
  {
    // Cache efficiency for particles (too many FindBin lookups)
    float* efficiencyAssociated = nullptr;
    if constexpr (step == CorrelationContainer::kCFStepCorrected) {
      if (cfg.mEfficiencyAssociated) {
        efficiencyAssociated = new float[tracks2.size()];
        int i = 0;
        for (auto& track : tracks2) {
          efficiencyAssociated[i++] = getEfficiencyCorrection(cfg.mEfficiencyAssociated, track.eta(), track.pt(), multiplicity, posZ);
        }
      }
    }

    for (auto& track1 : tracks1) {
      // LOGF(info, "Track %f | %f | %f  %d %d", track1.eta(), track1.phi(), track1.pt(), track1.isGlobalTrack(), track1.isGlobalTrackSDD());

      if constexpr (step <= CorrelationContainer::kCFStepTracked) {
        if (!checkObject<step>(track1)) {
          continue;
        }
      }

      if (cfgTriggerCharge != 0 && cfgTriggerCharge * track1.sign() < 0) {
        continue;
      }

      float triggerWeight = eventWeight;
      if constexpr (step == CorrelationContainer::kCFStepCorrected) {
        if (cfg.mEfficiencyTrigger) {
          triggerWeight *= getEfficiencyCorrection(cfg.mEfficiencyTrigger, track1.eta(), track1.pt(), multiplicity, posZ);
        }
      }

      target->getTriggerHist()->Fill(step, track1.pt(), multiplicity, posZ, triggerWeight);

      for (auto& track2 : tracks2) {
        if (track1.globalIndex() == track2.globalIndex()) {
          // LOGF(info, "Track identical: %f | %f | %f || %f | %f | %f", track1.eta(), track1.phi(), track1.pt(),  track2.eta(), track2.phi(), track2.pt());
          continue;
        }

        if constexpr (step <= CorrelationContainer::kCFStepTracked) {
          if (!checkObject<step>(track2)) {
            continue;
          }
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

        if constexpr (step >= CorrelationContainer::kCFStepReconstructed) {
          if (cfg.mPairCuts && mPairCuts.conversionCuts(track1, track2)) {
            continue;
          }

          if (cfgTwoTrackCut > 0 && mPairCuts.twoTrackCut(track1, track2, magField)) {
            continue;
          }
        }

        float associatedWeight = triggerWeight;
        if constexpr (step == CorrelationContainer::kCFStepCorrected) {
          if (cfg.mEfficiencyAssociated) {
            associatedWeight *= efficiencyAssociated[track2.filteredIndex()];
          }
        }

        float deltaPhi = track1.phi() - track2.phi();
        if (deltaPhi > 1.5f * PI) {
          deltaPhi -= TwoPI;
        }
        if (deltaPhi < -PIHalf) {
          deltaPhi += TwoPI;
        }

        target->getPairHist()->Fill(step,
                                    track1.eta() - track2.eta(), track2.pt(), track1.pt(), multiplicity, deltaPhi, posZ, associatedWeight);
      }
    }

    delete[] efficiencyAssociated;
  }

  void loadEfficiency(uint64_t timestamp)
  {
    if (cfg.efficiencyLoaded) {
      return;
    }
    if (cfgEfficiencyTrigger.value.empty() == false) {
      cfg.mEfficiencyTrigger = ccdb->getForTimeStamp<THnT<float>>(cfgEfficiencyTrigger, timestamp);
      if (cfg.mEfficiencyTrigger == nullptr) {
        LOGF(fatal, "Could not load efficiency histogram for trigger particles from %s", cfgEfficiencyTrigger.value.c_str());
      }
      LOGF(info, "Loaded efficiency histogram for trigger particles from %s (%p)", cfgEfficiencyTrigger.value.c_str(), (void*)cfg.mEfficiencyTrigger);
    }
    if (cfgEfficiencyAssociated.value.empty() == false) {
      cfg.mEfficiencyAssociated = ccdb->getForTimeStamp<THnT<float>>(cfgEfficiencyAssociated, timestamp);
      if (cfg.mEfficiencyAssociated == nullptr) {
        LOGF(fatal, "Could not load efficiency histogram for associated particles from %s", cfgEfficiencyAssociated.value.c_str());
      }
      LOGF(info, "Loaded efficiency histogram for associated particles from %s (%p)", cfgEfficiencyAssociated.value.c_str(), (void*)cfg.mEfficiencyAssociated);
    }
    cfg.efficiencyLoaded = true;
  }

  double getEfficiencyCorrection(THn* eff, float eta, float pt, float multiplicity, float posZ)
  {
    int effVars[4];
    effVars[0] = eff->GetAxis(0)->FindBin(eta);
    effVars[1] = eff->GetAxis(1)->FindBin(pt);
    effVars[2] = eff->GetAxis(2)->FindBin(multiplicity);
    effVars[3] = eff->GetAxis(3)->FindBin(posZ);
    return eff->GetBinContent(effVars);
  }

  // Version with explicit nested loop
  void processSameAOD(aodCollisions::iterator const& collision, aod::BCsWithTimestamps const&, aodTracks const& tracks)
  {
    // NOTE legacy function for O2 integration tests. Full version needs derived data

    if (cfgVerbosity > 0) {
      LOGF(info, "processSameAOD: Tracks for collision: %d | Vertex: %.1f | INT7: %d | V0M: %.1f", tracks.size(), collision.posZ(), collision.sel7(), collision.centRun2V0M());
    }

    // TODO will go to CCDBConfigurable
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    loadEfficiency(bc.timestamp());

    const auto multiplicity = collision.centRun2V0M();

    if (fillCollisionAOD(same, collision, multiplicity) == false) {
      return;
    }
    registry.fill(HIST("eventcount_same"), -2);
    fillQA(collision, multiplicity, tracks);
    fillCorrelations<CorrelationContainer::kCFStepReconstructed>(same, tracks, tracks, multiplicity, collision.posZ(), getMagneticField(bc.timestamp()), 1.0f);
  }
  PROCESS_SWITCH(CorrelationTask, processSameAOD, "Process same event on AOD", true);

  void processSameDerived(derivedCollisions::iterator const& collision, soa::Filtered<aod::CFTracks> const& tracks)
  {
    if (cfgVerbosity > 0) {
      LOGF(info, "processSameDerived: Tracks for collision: %d | Vertex: %.1f | Multiplicity/Centrality: %.1f", tracks.size(), collision.posZ(), collision.multiplicity());
    }
    loadEfficiency(collision.timestamp());

    const auto multiplicity = collision.multiplicity();
    int field = 0;
    if (cfgTwoTrackCut > 0) {
      field = getMagneticField(collision.timestamp());
    }

    int bin = configurableBinningDerived.getBin({collision.posZ(), collision.multiplicity()});
    registry.fill(HIST("eventcount_same"), bin);
    fillQA(collision, multiplicity, tracks);

    same->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
    fillCorrelations<CorrelationContainer::kCFStepReconstructed>(same, tracks, tracks, multiplicity, collision.posZ(), field, 1.0f);

    if (cfg.mEfficiencyAssociated || cfg.mEfficiencyTrigger) {
      same->fillEvent(multiplicity, CorrelationContainer::kCFStepCorrected);
      fillCorrelations<CorrelationContainer::kCFStepCorrected>(same, tracks, tracks, multiplicity, collision.posZ(), field, 1.0f);
    }
  }
  PROCESS_SWITCH(CorrelationTask, processSameDerived, "Process same event on derived data", false);

  using BinningTypeAOD = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentRun2V0M>;
  void processMixedAOD(aodCollisions& collisions, aodTracks const& tracks, aod::BCsWithTimestamps const&)
  {
    // NOTE legacy function for O2 integration tests. Full version needs derived data

    // Strictly upper categorised collisions, for cfgNoMixedEvents combinations per bin, skipping those in entry -1
    BinningTypeAOD configurableBinning{{axisVertex, axisMultiplicity}, true}; // true is for 'ignore overflows' (true by default). Underflows and overflows will have bin -1.
    auto tracksTuple = std::make_tuple(tracks);
    SameKindPair<aodCollisions, aodTracks, BinningTypeAOD> pairs{configurableBinning, cfgNoMixedEvents, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip

    int skipID = -1;
    for (auto it = pairs.begin(); it != pairs.end(); it++) {
      auto& [collision1, tracks1, collision2, tracks2] = *it;
      int bin = configurableBinning.getBin({collision1.posZ(), collision1.centRun2V0M()});
      if (cfgVerbosity > 0) {
        LOGF(info, "processMixedAOD: Mixed collisions bin: %d pair: %d (%.3f, %.3f), %d (%.3f, %.3f)", bin, collision1.globalIndex(), collision1.posZ(), collision1.centRun2V0M(), collision2.globalIndex(), collision2.posZ(), collision2.centRun2V0M());
      }
      if (collision1.globalIndex() == skipID) {
        continue;
      }

      if (it.isNewWindow()) {
        skipID = -1;
        if (fillCollisionAOD(mixed, collision1, collision1.centRun2V0M()) == false) {
          skipID = collision1.globalIndex();
          continue;
        }
      }
      if (!collision2.alias_bit(kINT7) || !collision2.sel7()) {
        continue;
      }

      registry.fill(HIST("eventcount_mixed"), bin);

      auto bc = collision1.bc_as<aod::BCsWithTimestamps>();

      // LOGF(info, "Tracks: %d and %d entries", tracks1.size(), tracks2.size());

      fillCorrelations<CorrelationContainer::kCFStepReconstructed>(mixed, tracks1, tracks2, collision1.centRun2V0M(), collision1.posZ(), getMagneticField(bc.timestamp()), 1.0f / it.currentWindowNeighbours());
    }
  }
  PROCESS_SWITCH(CorrelationTask, processMixedAOD, "Process mixed events on AOD", false);

  using BinningTypeDerived = ColumnBinningPolicy<aod::collision::PosZ, aod::cfcollision::Multiplicity>;
  BinningTypeDerived configurableBinningDerived{{axisVertex, axisMultiplicity}, true}; // true is for 'ignore overflows' (true by default). Underflows and overflows will have bin -1.
  void processMixedDerived(derivedCollisions& collisions, derivedTracks const& tracks)
  {
    // Strictly upper categorised collisions, for cfgNoMixedEvents combinations per bin, skipping those in entry -1
    auto tracksTuple = std::make_tuple(tracks);
    SameKindPair<derivedCollisions, derivedTracks, BinningTypeDerived> pairs{configurableBinningDerived, cfgNoMixedEvents, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip

    for (auto it = pairs.begin(); it != pairs.end(); it++) {
      auto& [collision1, tracks1, collision2, tracks2] = *it;
      int bin = configurableBinningDerived.getBin({collision1.posZ(), collision1.multiplicity()});
      float eventWeight = 1.0f / it.currentWindowNeighbours();
      int field = 0;
      if (cfgTwoTrackCut > 0) {
        field = getMagneticField(collision1.timestamp());
      }

      if (cfgVerbosity > 0) {
        LOGF(info, "processMixedDerived: Mixed collisions bin: %d pair: [%d, %d] %d (%.3f, %.3f), %d (%.3f, %.3f)", bin, it.isNewWindow(), it.currentWindowNeighbours(), collision1.globalIndex(), collision1.posZ(), collision1.multiplicity(), collision2.globalIndex(), collision2.posZ(), collision2.multiplicity());
      }

      if (it.isNewWindow()) {
        loadEfficiency(collision1.timestamp());

        mixed->fillEvent(collision1.multiplicity(), CorrelationContainer::kCFStepReconstructed);
      }

      // LOGF(info, "Tracks: %d and %d entries", tracks1.size(), tracks2.size());

      registry.fill(HIST("eventcount_mixed"), bin);
      fillCorrelations<CorrelationContainer::kCFStepReconstructed>(mixed, tracks1, tracks2, collision1.multiplicity(), collision1.posZ(), field, eventWeight);

      if (cfg.mEfficiencyAssociated || cfg.mEfficiencyTrigger) {
        if (it.isNewWindow()) {
          mixed->fillEvent(collision1.multiplicity(), CorrelationContainer::kCFStepCorrected);
        }
        fillCorrelations<CorrelationContainer::kCFStepCorrected>(mixed, tracks1, tracks2, collision1.multiplicity(), collision1.posZ(), field, eventWeight);
      }
    }
  }
  PROCESS_SWITCH(CorrelationTask, processMixedDerived, "Process mixed events on derived data", false);

  // Version with combinations
  /*void processWithCombinations(soa::Join<aod::Collisions, aod::CentRun2V0Ms>::iterator const& collision, aod::BCsWithTimestamps const&, soa::Filtered<aod::Tracks> const& tracks)
  {
    LOGF(info, "Tracks for collision (Combination run): %d", tracks.size());

    // TODO will go to CCDBConfigurable
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();

    const auto multiplicity = collision.centRun2V0M();

    for (auto track1 = tracks.begin(); track1 != tracks.end(); ++track1) {

      if (cfgTriggerCharge != 0 && cfgTriggerCharge * track1.sign() < 0) {
        continue;
      }

      //       LOGF(info, "TRACK %f %f | %f %f | %f %f", track1.eta(), track1.eta(), track1.phi(), track1.phi2(), track1.pt(), track1.pt());

      same->getTriggerHist()->Fill(CorrelationContainer::kCFStepReconstructed, track1.pt(), multiplicity, collision.posZ());
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
                                track1.eta() - track2.eta(), track2.pt(), track1.pt(), multiplicity, deltaPhi, collision.posZ());
      // mixed->getPairHist()->Fill(values, CorrelationContainer::kCFStepReconstructed);
    }
  }
  PROCESS_SWITCH(CorrelationTask, processWithCombinations, "Process same event on AOD with combinations", false);*/

  int GetSpecies(int pdgCode)
  {
    switch (pdgCode) {
      case 211: // pion
      case -211:
        return 0;
      case 321: // Kaon
      case -321:
        return 1;
      case 2212: // proton
      case -2212:
        return 2;
      default:
        return 3;
    }
  }

  // NOTE SmallGroups includes soa::Filtered always
  Preslice<aod::CFTracksWithLabel> perCollision = aod::cftrack::cfCollisionId;
  void processMCEfficiency(soa::Filtered<aod::CFMcCollisions>::iterator const& mcCollision, aod::CFMcParticles const& mcParticles, soa::SmallGroups<aod::CFCollisionsWithLabel> const& collisions, aod::CFTracksWithLabel const& tracks)
  {
    if (cfgVerbosity > 0) {
      LOGF(info, "MC collision at vtx-z = %f with %d mc particles and %d reconstructed collisions", mcCollision.posZ(), mcParticles.size(), collisions.size());
    }

    // Primaries
    auto multiplicity = mcCollision.multiplicity();
    for (auto& mcParticle : mcParticles) {
      if (mcParticle.isPhysicalPrimary() && mcParticle.sign() != 0) {
        same->getTrackHistEfficiency()->Fill(CorrelationContainer::MC, mcParticle.eta(), mcParticle.pt(), GetSpecies(mcParticle.pdgCode()), multiplicity, mcCollision.posZ());
      }
    }

    for (auto& collision : collisions) {
      auto groupedTracks = tracks.sliceBy(perCollision, collision.globalIndex());
      if (cfgVerbosity > 0) {
        LOGF(info, "  Reconstructed collision at vtx-z = %f", collision.posZ());
        LOGF(info, "  which has %d tracks", groupedTracks.size());
      }

      for (auto& track : groupedTracks) {
        if (track.has_cfMCParticle()) {
          const auto& mcParticle = track.cfMCParticle();
          if (mcParticle.isPhysicalPrimary()) {
            same->getTrackHistEfficiency()->Fill(CorrelationContainer::RecoPrimaries, mcParticle.eta(), mcParticle.pt(), GetSpecies(mcParticle.pdgCode()), multiplicity, mcCollision.posZ());
          }
          same->getTrackHistEfficiency()->Fill(CorrelationContainer::RecoAll, mcParticle.eta(), mcParticle.pt(), GetSpecies(mcParticle.pdgCode()), multiplicity, mcCollision.posZ());
          // LOGF(info, "Filled track %d", track.globalIndex());
        } else {
          // fake track
          same->getTrackHistEfficiency()->Fill(CorrelationContainer::Fake, track.eta(), track.pt(), 0, multiplicity, mcCollision.posZ());
        }
      }
    }
  }
  PROCESS_SWITCH(CorrelationTask, processMCEfficiency, "MC: Extract efficiencies", false);

  // NOTE SmallGroups includes soa::Filtered always
  void processMCSameDerived(soa::Filtered<aod::CFMcCollisions>::iterator const& mcCollision, soa::Filtered<aod::CFMcParticles> const& mcParticles, soa::SmallGroups<aod::CFCollisionsWithLabel> const& collisions)
  {
    if (cfgVerbosity > 0) {
      LOGF(info, "processMCSameDerived. MC collision: %d, particles: %d, collisions: %d", mcCollision.globalIndex(), mcParticles.size(), collisions.size());
    }

    same->fillEvent(mcCollision.multiplicity(), CorrelationContainer::kCFStepAll);
    fillCorrelations<CorrelationContainer::kCFStepAll>(same, mcParticles, mcParticles, mcCollision.multiplicity(), mcCollision.posZ(), 0, 1.0f);

    if (collisions.size() == 0) {
      return;
    }

    same->fillEvent(mcCollision.multiplicity(), CorrelationContainer::kCFStepVertex);
    fillCorrelations<CorrelationContainer::kCFStepVertex>(same, mcParticles, mcParticles, mcCollision.multiplicity(), mcCollision.posZ(), 0, 1.0f);

    same->fillEvent(mcCollision.multiplicity(), CorrelationContainer::kCFStepTrackedOnlyPrim);
    fillCorrelations<CorrelationContainer::kCFStepTrackedOnlyPrim>(same, mcParticles, mcParticles, mcCollision.multiplicity(), mcCollision.posZ(), 0, 1.0f);

    same->fillEvent(mcCollision.multiplicity(), CorrelationContainer::kCFStepTracked);
    fillCorrelations<CorrelationContainer::kCFStepTracked>(same, mcParticles, mcParticles, mcCollision.multiplicity(), mcCollision.posZ(), 0, 1.0f);

    // NOTE kCFStepReconstructed and kCFStepCorrected are filled in processSameDerived
    //      This also means that if a MC collision had several reconstructed vertices (collisions), all of them are filled
  }
  PROCESS_SWITCH(CorrelationTask, processMCSameDerived, "Process MC same event on derived data", false);

  using BinningTypeMCDerived = ColumnBinningPolicy<aod::mccollision::PosZ, aod::cfmccollision::Multiplicity>;
  PresliceUnsorted<aod::CFCollisionsWithLabel> collisionPerMCCollision = aod::cfcollision::cfMcCollisionId;
  void processMCMixedDerived(soa::Filtered<aod::CFMcCollisions>& mcCollisions, soa::Filtered<aod::CFMcParticles> const& mcParticles, soa::Filtered<aod::CFCollisionsWithLabel> const& collisions)
  {
    // Strictly upper categorised collisions, for cfgNoMixedEvents combinations per bin, skipping those in entry -1
    BinningTypeMCDerived configurableBinning{{axisVertex, axisMultiplicity}, true}; // true is for 'ignore overflows' (true by default). Underflows and overflows will have bin -1.
    auto tuple = std::make_tuple(mcParticles);
    SameKindPair<soa::Filtered<aod::CFMcCollisions>, soa::Filtered<aod::CFMcParticles>, BinningTypeMCDerived> pairs{configurableBinning, cfgNoMixedEvents, -1, mcCollisions, tuple, &cache}; // -1 is the number of the bin to skip

    for (auto it = pairs.begin(); it != pairs.end(); it++) {
      auto& [collision1, tracks1, collision2, tracks2] = *it;
      float eventWeight = 1.0f / it.currentWindowNeighbours();

      if (cfgVerbosity > 0) {
        int bin = configurableBinning.getBin({collision1.posZ(), collision1.multiplicity()});
        LOGF(info, "processMCMixedDerived: Mixed collisions bin: %d pair: [%d, %d] %d (%.3f, %.3f), %d (%.3f, %.3f)", bin, it.isNewWindow(), it.currentWindowNeighbours(), collision1.globalIndex(), collision1.posZ(), collision1.multiplicity(), collision2.globalIndex(), collision2.posZ(), collision2.multiplicity());
      }

      // STEP 0
      if (it.isNewWindow()) {
        mixed->fillEvent(collision1.multiplicity(), CorrelationContainer::kCFStepAll);
      }
      fillCorrelations<CorrelationContainer::kCFStepAll>(mixed, tracks1, tracks2, collision1.multiplicity(), collision1.posZ(), 0, eventWeight);

      // check if collision1 has at least one reconstructed collision
      auto groupedCollisions = collisions.sliceBy(collisionPerMCCollision, collision1.globalIndex());
      if (cfgVerbosity > 0) {
        LOGF(info, "Found %d related collisions", groupedCollisions.size());
      }
      if (groupedCollisions.size() == 0) {
        continue;
      }

      // STEP 2, 4, 5
      if (it.isNewWindow()) {
        mixed->fillEvent(collision1.multiplicity(), CorrelationContainer::kCFStepVertex);
        mixed->fillEvent(collision1.multiplicity(), CorrelationContainer::kCFStepTrackedOnlyPrim);
        mixed->fillEvent(collision1.multiplicity(), CorrelationContainer::kCFStepTracked);
      }
      fillCorrelations<CorrelationContainer::kCFStepVertex>(mixed, tracks1, tracks2, collision1.multiplicity(), collision1.posZ(), 0, eventWeight);
      fillCorrelations<CorrelationContainer::kCFStepTrackedOnlyPrim>(mixed, tracks1, tracks2, collision1.multiplicity(), collision1.posZ(), 0, eventWeight);
      fillCorrelations<CorrelationContainer::kCFStepTracked>(mixed, tracks1, tracks2, collision1.multiplicity(), collision1.posZ(), 0, eventWeight);

      // NOTE kCFStepReconstructed and kCFStepCorrected are filled in processMixedDerived
      //      This also means that if a MC collision had several reconstructed vertices (collisions), all of them are filled
    }
  }
  PROCESS_SWITCH(CorrelationTask, processMCMixedDerived, "Process MC mixed events on derived data", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<CorrelationTask>(cfgc)};
}
