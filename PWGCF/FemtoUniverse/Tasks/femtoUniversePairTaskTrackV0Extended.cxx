// Copyright 2019-2022 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \brief Tasks that build pairs of track particles and v0s
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de
/// \author Zuzanna Chochulska, WUT Warsaw, zuzanna.chochulska.stud@pw.edu.pl
/// \author Shirajum Monira, WUT Warsaw, shirajum.monira.dokt@pw.edu.pl

#include <vector>
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/StepTHn.h"
#include "PWGCF/FemtoUniverse/DataModel/FemtoDerived.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseParticleHisto.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseEventHisto.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniversePairCleaner.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseContainer.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseDetaDphiStar.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUtils.h"

using namespace o2;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::analysis::femtoUniverse;
using namespace o2::aod::pidutils;

struct femtoUniversePairTaskTrackV0Extended {

  SliceCache cache;
  using FemtoFullParticles = soa::Join<aod::FDParticles, aod::FDExtParticles>;
  Preslice<FemtoFullParticles> perCol = aod::femtouniverseparticle::fdCollisionId;

  /// Particle 1 (track)
  Configurable<int> ConfTrkPDGCodePartOne{"ConfTrkPDGCodePartOne", 211, "Particle 1 (Track) - PDG code"};
  Configurable<int> ConfTrackChoicePartOne{"ConfTrackChoicePartOne", 1, "0:Proton, 1:Pion, 2:Kaon"};
  ConfigurableAxis ConfTrkTempFitVarBins{"ConfTrkDTempFitVarBins", {300, -0.15, 0.15}, "binning of the TempFitVar in the pT vs. TempFitVar plot"};
  ConfigurableAxis ConfTrkTempFitVarpTBins{"ConfTrkTempFitVarpTBins", {20, 0.5, 4.05}, "pT binning of the pT vs. TempFitVar plot"};
  Configurable<int> ConfChargePart1{"ConfChargePart1", 0, "sign of particle 1"};
  Configurable<float> ConfHPtPart1{"ConfHPtPart1", 4.0f, "higher limit for pt of particle 1"};
  Configurable<float> ConfLPtPart1{"ConfLPtPart1", 0.3f, "lower limit for pt of particle 1"};
  Configurable<float> Confmom{"Confmom", 0.5, "momentum threshold for particle identification using TOF"};
  Configurable<float> ConfNsigmaTPCParticle{"ConfNsigmaTPCParticle", 3.0, "TPC Sigma for particle momentum < Confmom"};
  Configurable<float> ConfNsigmaCombinedParticle{"ConfNsigmaCombinedParticle", 3.0, "TPC and TOF Sigma (combined) for particle momentum > Confmom"};

  /// Partition for particle 1
  Partition<FemtoFullParticles> partsOne = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kTrack)) && aod::femtouniverseparticle::sign == ConfChargePart1 && aod::femtouniverseparticle::pt < ConfHPtPart1 && aod::femtouniverseparticle::pt > ConfLPtPart1;

  /// Histogramming for particle 1
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kTrack, 3> trackHistoPartOnePos;
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kTrack, 4> trackHistoPartOneNeg;

  /// Particle 2 (V0)
  Configurable<int> ConfV0PDGCodePartTwo{"ConfV0PDGCodePartTwo", 3122, "Particle 2 (V0) - PDG code"};
  ConfigurableAxis ConfV0TempFitVarBins{"ConfV0TempFitVarBins", {300, 0.95, 1.}, "V0: binning of the TempFitVar in the pT vs. TempFitVar plot"};
  ConfigurableAxis ConfV0TempFitVarpTBins{"ConfV0TempFitVarpTBins", {20, 0.5, 4.05}, "V0: pT binning of the pT vs. TempFitVar plot"};
  Configurable<int> ConfV0PosChildType{"ConfV0PosChildType", -1, "identifier for positive v0 daughter"};
  ConfigurableAxis ConfChildTempFitVarBins{"ConfChildTempFitVarBins", {300, -0.15, 0.15}, "V0 child: binning of the TempFitVar in the pT vs. TempFitVar plot"};
  ConfigurableAxis ConfChildTempFitVarpTBins{"ConfChildTempFitVarpTBins", {20, 0.5, 4.05}, "V0 child: pT binning of the pT vs. TempFitVar plot"};
  Configurable<float> ConfHPtPart2{"ConfHPtPart2", 4.0f, "higher limit for pt of particle 2"};
  Configurable<float> ConfLPtPart2{"ConfLPtPart2", 0.3f, "lower limit for pt of particle 2"};

  /// Partition for particle 2
  Partition<FemtoFullParticles> partsTwo = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kV0)) && aod::femtouniverseparticle::pt < ConfHPtPart2 && aod::femtouniverseparticle::pt > ConfLPtPart2;

  /// Histogramming for particle 2
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kV0, 2> trackHistoPartTwo;
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kV0Child, 3> posChildHistos;
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kV0Child, 4> negChildHistos;

  /// Histogramming for Event
  FemtoUniverseEventHisto eventHisto;

  /// Correlation part
  // Configurable<int> ConfTrackChoicePartTwo{"ConfTrackChoicePartTwo", 1, "0:Proton, 1:Pion, 2:Kaon"}; //not used
  Configurable<bool> ConfIsMC{"ConfIsMC", false, "Enable additional Histogramms in the case of a MonteCarlo Run"};
  Configurable<bool> ConfUse3D{"ConfUse3D", false, "Enable three dimensional histogramms (to be used only for analysis with high statistics): k* vs mT vs multiplicity"};
  ConfigurableAxis ConfMultBins{"ConfMultBins", {VARIABLE_WIDTH, 0.0f, 20.0f, 40.0f, 60.0f, 80.0f, 100.0f, 200.0f, 99999.f}, "Mixing bins - multiplicity"};
  ConfigurableAxis ConfVtxBins{"ConfVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis ConfkstarBins{"ConfkstarBins", {1500, 0., 6.}, "binning kstar"};
  ConfigurableAxis ConfkTBins{"ConfkTBins", {150, 0., 9.}, "binning kT"};
  ConfigurableAxis ConfmTBins{"ConfmTBins", {225, 0., 7.5}, "binning mT"};
  // Configurable<int> ConfNEventsMix{"ConfNEventsMix", 5, "Number of events for mixing"}; // not used
  Configurable<bool> ConfIsCPR{"ConfIsCPR", true, "Close Pair Rejection"};
  Configurable<bool> ConfCPRPlotPerRadii{"ConfCPRPlotPerRadii", false, "Plot CPR per radii"};
  Configurable<float> ConfCPRdeltaPhiMax{"ConfCPRdeltaPhiMax", 0.01, "Max. Delta Phi for Close Pair Rejection"};
  Configurable<float> ConfCPRdeltaEtaMax{"ConfCPRdeltaEtaMax", 0.01, "Max. Delta Eta for Close Pair Rejection"};
  Configurable<float> ConfCPRdeltaPhiCut{"ConfCPRdeltaPhiCut", 0.0, "Delta Phi cut for Close Pair Rejection"};
  Configurable<float> ConfCPRdeltaEtaCut{"ConfCPRdeltaEtaCut", 0.0, "Delta Eta cut for Close Pair Rejection"};
  Configurable<int> ConfPhiBins{"ConfPhiBins", 29, "Number of phi bins in deta dphi"};
  Configurable<int> ConfEtaBins{"ConfEtaBins", 29, "Number of eta bins in deta dphi"};
  ConfigurableAxis ConfmTBins3D{"ConfmTBins3D", {VARIABLE_WIDTH, 1.02f, 1.14f, 1.20f, 1.26f, 1.38f, 1.56f, 1.86f, 4.50f}, "mT Binning for the 3Dimensional plot: k* vs multiplicity vs mT (set <<ConfUse3D>> to true in order to use)"};
  ConfigurableAxis ConfmultBins3D{"ConfMultBins3D", {VARIABLE_WIDTH, 0.0f, 20.0f, 30.0f, 40.0f, 99999.0f}, "multiplicity Binning for the 3Dimensional plot: k* vs multiplicity vs mT (set <<ConfUse3D>> to true in order to use)"};

  FemtoUniverseContainer<femtoUniverseContainer::EventType::same, femtoUniverseContainer::Observable::kstar> sameEventCont;
  FemtoUniverseContainer<femtoUniverseContainer::EventType::mixed, femtoUniverseContainer::Observable::kstar> mixedEventCont;
  FemtoUniversePairCleaner<aod::femtouniverseparticle::ParticleType::kTrack, aod::femtouniverseparticle::ParticleType::kV0> pairCleaner;
  FemtoUniverseDetaDphiStar<aod::femtouniverseparticle::ParticleType::kTrack, aod::femtouniverseparticle::ParticleType::kV0> pairCloseRejection;
  /// Histogram output
  HistogramRegistry qaRegistry{"TrackQA", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry resultRegistry{"Correlations", {}, OutputObjHandlingPolicy::AnalysisObject};

  bool IsNSigmaCombined(float mom, float nsigmaTPCParticle, float nsigmaTOFParticle)
  {
    if (mom <= Confmom) {
      if (TMath::Abs(nsigmaTPCParticle) < ConfNsigmaTPCParticle) {
        return true;
      } else {
        return false;
      }
    } else {
      if (TMath::Hypot(nsigmaTOFParticle, nsigmaTPCParticle) < ConfNsigmaCombinedParticle) {
        return true;
      } else {
        return false;
      }
    }
  }

  bool IsNSigmaTPC(float nsigmaTPCParticle)
  {
    if (TMath::Abs(nsigmaTPCParticle) < ConfNsigmaTPCParticle) {
      return true;
    } else {
      return false;
    }
  }

  template <typename T>
  bool IsParticleCombined(const T& part, int id)
  {
    const float tpcNSigmas[3] = {unPackInTable<aod::pidtpc_tiny::binning>(part.tpcNSigmaStorePr()), unPackInTable<aod::pidtpc_tiny::binning>(part.tpcNSigmaStorePi()), unPackInTable<aod::pidtpc_tiny::binning>(part.tpcNSigmaStoreKa())};
    // const float tofNSigmas[3] = {part.tofNSigmaPr(), part.tofNSigmaPi(), part.tofNSigmaKa()};
    const float tofNSigmas[3] = {unPackInTable<aod::pidtof_tiny::binning>(part.tofNSigmaStorePr()), unPackInTable<aod::pidtof_tiny::binning>(part.tofNSigmaStorePi()), unPackInTable<aod::pidtof_tiny::binning>(part.tofNSigmaStoreKa())};

    return IsNSigmaCombined(part.p(), tpcNSigmas[id], tofNSigmas[id]);
  }

  template <typename T>
  bool IsParticleTPC(const T& part, int id)
  {
    const float tpcNSigmas[3] = {unPackInTable<aod::pidtpc_tiny::binning>(part.tpcNSigmaStorePr()), unPackInTable<aod::pidtpc_tiny::binning>(part.tpcNSigmaStorePi()), unPackInTable<aod::pidtpc_tiny::binning>(part.tpcNSigmaStoreKa())};

    return IsNSigmaTPC(tpcNSigmas[id]);
  }

  void init(InitContext&)
  {
    eventHisto.init(&qaRegistry);
    qaRegistry.add("Tracks_pos/nSigmaTPC", "; #it{p} (GeV/#it{c}); n#sigma_{TPC}", kTH2F, {{100, 0, 10}, {100, -5, 5}});
    qaRegistry.add("Tracks_pos/nSigmaTOF", "; #it{p} (GeV/#it{c}); n#sigma_{TOF}", kTH2F, {{100, 0, 10}, {100, -5, 5}});
    qaRegistry.add("Tracks_neg/nSigmaTPC", "; #it{p} (GeV/#it{c}); n#sigma_{TPC}", kTH2F, {{100, 0, 10}, {100, -5, 5}});
    qaRegistry.add("Tracks_neg/nSigmaTOF", "; #it{p} (GeV/#it{c}); n#sigma_{TOF}", kTH2F, {{100, 0, 10}, {100, -5, 5}});
    trackHistoPartOnePos.init(&qaRegistry, ConfTrkTempFitVarpTBins, ConfTrkTempFitVarBins, ConfIsMC, ConfTrkPDGCodePartOne);
    trackHistoPartOneNeg.init(&qaRegistry, ConfTrkTempFitVarpTBins, ConfTrkTempFitVarBins, ConfIsMC, ConfTrkPDGCodePartOne);
    trackHistoPartTwo.init(&qaRegistry, ConfV0TempFitVarpTBins, ConfV0TempFitVarBins, ConfIsMC, ConfV0PDGCodePartTwo);
    posChildHistos.init(&qaRegistry, ConfChildTempFitVarpTBins, ConfChildTempFitVarBins, false, false);
    negChildHistos.init(&qaRegistry, ConfChildTempFitVarpTBins, ConfChildTempFitVarBins, false, false);

    sameEventCont.init(&resultRegistry, ConfkstarBins, ConfMultBins, ConfkTBins, ConfmTBins, ConfmultBins3D, ConfmTBins3D, ConfEtaBins, ConfPhiBins, ConfIsMC, ConfUse3D);
    sameEventCont.setPDGCodes(ConfTrkPDGCodePartOne, ConfV0PDGCodePartTwo);
    mixedEventCont.init(&resultRegistry, ConfkstarBins, ConfMultBins, ConfkTBins, ConfmTBins, ConfmultBins3D, ConfmTBins3D, ConfEtaBins, ConfPhiBins, ConfIsMC, ConfUse3D);
    mixedEventCont.setPDGCodes(ConfTrkPDGCodePartOne, ConfV0PDGCodePartTwo);

    pairCleaner.init(&qaRegistry);
    if (ConfIsCPR.value) {
      pairCloseRejection.init(&resultRegistry, &qaRegistry, ConfCPRdeltaPhiMax.value, ConfCPRdeltaEtaMax.value, ConfCPRdeltaPhiCut.value, ConfCPRdeltaEtaCut.value, ConfCPRPlotPerRadii.value);
    }
  }
  /// This function processes the same event and takes care of all the histogramming
  void processSameEvent(o2::aod::FDCollision& col, FemtoFullParticles& parts)
  {
    const auto& magFieldTesla = col.magField();

    auto groupPartsOne = partsOne->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto groupPartsTwo = partsTwo->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    const int multCol = col.multNtr();

    eventHisto.fillQA(col);

    /// Histogramming same event
    for (auto& part : groupPartsTwo) {
      const auto& posChild = parts.iteratorAt(part.index() - 2);
      const auto& negChild = parts.iteratorAt(part.index() - 1);
      // printf("-- V0 d- %d d+ %d\n",negChild.globalIndex(),posChild.globalIndex());
      if (ConfV0PosChildType >= 0 && !IsParticleTPC(posChild, ConfV0PosChildType))
        continue;

      trackHistoPartTwo.fillQA<false, false>(part);
      posChildHistos.fillQA<false, false>(posChild);
      negChildHistos.fillQA<false, false>(negChild);
    }

    for (auto& part : groupPartsOne) {
      /// PID plot for particle 1
      const float tpcNSigmas[3] = {unPackInTable<aod::pidtpc_tiny::binning>(part.tpcNSigmaStorePr()), unPackInTable<aod::pidtpc_tiny::binning>(part.tpcNSigmaStorePi()), unPackInTable<aod::pidtpc_tiny::binning>(part.tpcNSigmaStoreKa())};
      const float tofNSigmas[3] = {unPackInTable<aod::pidtof_tiny::binning>(part.tofNSigmaStorePr()), unPackInTable<aod::pidtof_tiny::binning>(part.tofNSigmaStorePi()), unPackInTable<aod::pidtof_tiny::binning>(part.tofNSigmaStoreKa())};

      if (!IsNSigmaCombined(part.p(), tpcNSigmas[ConfTrackChoicePartOne], tofNSigmas[ConfTrackChoicePartOne]))
        continue;
      if (part.sign() > 0) {
        qaRegistry.fill(HIST("Tracks_pos/nSigmaTPC"), part.p(), tpcNSigmas[ConfTrackChoicePartOne]);
        qaRegistry.fill(HIST("Tracks_pos/nSigmaTOF"), part.p(), tofNSigmas[ConfTrackChoicePartOne]);
        trackHistoPartOnePos.fillQA<false, false>(part);
      } else if (part.sign() < 0) {
        qaRegistry.fill(HIST("Tracks_neg/nSigmaTPC"), part.p(), tpcNSigmas[ConfTrackChoicePartOne]);
        qaRegistry.fill(HIST("Tracks_neg/nSigmaTOF"), part.p(), tofNSigmas[ConfTrackChoicePartOne]);
        trackHistoPartOneNeg.fillQA<false, false>(part);
      }
    }

    /// Now build the combinations
    for (auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(groupPartsOne, groupPartsTwo))) {
      if (ConfIsCPR.value) {
        if (pairCloseRejection.isClosePair(p1, p2, parts, magFieldTesla)) {
          continue;
        }
      }
      // track cleaning
      if (!pairCleaner.isCleanPair(p1, p2, parts)) {
        continue;
      }
      /// PID using stored binned nsigma
      if (!IsParticleCombined(p1, ConfTrackChoicePartOne))
        continue;
      const auto& posChild = parts.iteratorAt(p2.index() - 2);
      if (ConfV0PosChildType >= 0 && !IsParticleTPC(posChild, ConfV0PosChildType))
        continue;
      sameEventCont.setPair<false>(p1, p2, multCol, ConfUse3D);
    }
  }

  PROCESS_SWITCH(femtoUniversePairTaskTrackV0Extended, processSameEvent, "Enable processing same event", true);

  /// This function processes the mixed event
  void processMixedEvent(o2::aod::FDCollisions& cols, FemtoFullParticles& parts)
  {
    ColumnBinningPolicy<aod::collision::PosZ, aod::femtouniversecollision::MultNtr> colBinning{{ConfVtxBins, ConfMultBins}, true};

    for (auto& [collision1, collision2] : soa::selfCombinations(colBinning, 5, -1, cols, cols)) {

      const int multCol = collision1.multNtr();

      auto groupPartsOne = partsOne->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);
      auto groupPartsTwo = partsTwo->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision2.globalIndex(), cache);

      const auto& magFieldTesla1 = collision1.magField();
      const auto& magFieldTesla2 = collision2.magField();

      if (magFieldTesla1 != magFieldTesla2) {
        continue;
      }

      for (auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(groupPartsOne, groupPartsTwo))) {
        if (ConfIsCPR.value) {
          if (pairCloseRejection.isClosePair(p1, p2, parts, magFieldTesla1)) {
            continue;
          }
        }
        // track cleaning
        if (!pairCleaner.isCleanPair(p1, p2, parts)) {
          continue;
        }
        /// PID using stored binned nsigma
        if (!IsParticleCombined(p1, ConfTrackChoicePartOne))
          continue;
        const auto& posChild = parts.iteratorAt(p2.index() - 2);
        if (ConfV0PosChildType >= 0 && !IsParticleTPC(posChild, ConfV0PosChildType))
          continue;
        mixedEventCont.setPair<false>(p1, p2, multCol, ConfUse3D);
      }
    }
  }

  PROCESS_SWITCH(femtoUniversePairTaskTrackV0Extended, processMixedEvent, "Enable processing mixed events", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<femtoUniversePairTaskTrackV0Extended>(cfgc),
  };
  return workflow;
}
