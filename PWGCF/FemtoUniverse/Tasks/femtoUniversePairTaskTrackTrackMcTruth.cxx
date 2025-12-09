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

/// \file femtoUniversePairTaskTrackTrackMcTruth.cxx
/// \brief Tasks that reads the track tables used for the pairing and builds pairs of two tracks
/// \author Malgorzata Janik, WUT, majanik@cern.ch
/// \author Zuzanna Chochulska, WUT Warsaw & CTU Prague, zchochul@cern.ch

#include "PWGCF/FemtoUniverse/Core/FemtoUniverseContainer.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseEventHisto.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniversePairCleaner.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseParticleHisto.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseTrackSelection.h"
#include "PWGCF/FemtoUniverse/Core/femtoUtils.h"
#include "PWGCF/FemtoUniverse/DataModel/FemtoDerived.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/StepTHn.h"
#include "Framework/runDataProcessing.h"

#include <random>
#include <vector>

using namespace o2;
using namespace o2::analysis::femto_universe;
using namespace o2::constants::physics;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

struct FemtoUniversePairTaskTrackTrackMcTruth {
  SliceCache cache;
  Preslice<aod::FDParticles> perCol = aod::femtouniverseparticle::fdCollisionId;

  /// Particle selection part

  /// Configurables for both particles
  Configurable<bool> confUse3D{"confUse3D", false, "Enable three dimensional histogramms (to be used only for analysis with high statistics): k* vs mT vs multiplicity"};
  Configurable<float> confEtaMax{"confEtaMax", 0.8f, "Higher limit for |Eta| (the same for both particles)"};

  /// Particle 1
  Configurable<int> confPDGCodePartOne{"confPDGCodePartOne", 2212, "Particle 1 - PDG code"};
  Configurable<bool> confNoPDGPartOne{"confNoPDGPartOne", false, "0: selecting part by PDG, 1: no PID selection"};
  Configurable<float> confPtLowPart1{"confPtLowPart1", 0.2, "Lower limit for Pt for the first particle"};
  Configurable<float> confPtHighPart1{"confPtHighPart1", 2.5, "Higher limit for Pt for the first particle"};

  /// Partition for particle 1
  Partition<aod::FDParticles> partsOne = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kMCTruthTrack)) &&
                                         (aod::femtouniverseparticle::pt < confPtHighPart1) && (aod::femtouniverseparticle::pt > confPtLowPart1) && (nabs(aod::femtouniverseparticle::eta) < confEtaMax);

  /// Histogramming for particle 1
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kMCTruthTrack, 1> trackHistoPartOne;

  /// Particle 2
  Configurable<bool> confIsSame{"confIsSame", false, "Pairs of the same particle"};
  Configurable<int> confPDGCodePartTwo{"confPDGCodePartTwo", 333, "Particle 2 - PDG code"};
  Configurable<bool> confNoPDGPartTwo{"confNoPDGPartTwo", false, "0: selecting part by PDG, 1: no PID selection"};
  Configurable<float> confPtLowPart2{"confPtLowPart2", 0.2, "Lower limit for Pt for the second particle"};
  Configurable<float> confPtHighPart2{"confPtHighPart2", 2.5, "Higher limit for Pt for the second particle"};

  /// Partition for particle 2
  Partition<aod::FDParticles> partsTwo = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kMCTruthTrack)) &&
                                         (aod::femtouniverseparticle::pt < confPtHighPart2) && (aod::femtouniverseparticle::pt > confPtLowPart2) && (nabs(aod::femtouniverseparticle::eta) < confEtaMax);

  /// Histogramming for particle 2
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kMCTruthTrack, 2> trackHistoPartTwo;

  /// Histogramming for Event
  FemtoUniverseEventHisto eventHisto;

  /// The configurables need to be passed to an std::vector
  int vPIDPartOne, vPIDPartTwo;
  std::vector<float> kNsigma;

  // D0/D0bar options
  Configurable<bool> confActiveD0OriginCheck{"confActiveD0OriginCheck", false, "If true - calculate correlation for D0/D0bar mesons with a given origin"};
  Configurable<int8_t> confD0OriginFlag{"confD0OriginFlag", 1, "D0/D0bar origin: 0 - none, 1 - prompt, 2 - non-prompt"};

  /// particle part
  ConfigurableAxis confTempFitVarpTBins{"confTempFitVarpTBins", {20, 0.5, 4.05}, "pT binning of the pT vs. TempFitVar plot"};
  ConfigurableAxis confTempFitVarPDGBins{"confTempFitVarPDGBins", {6000, -2300, 2300}, "binning of the TempFitVar in the pT vs. TempFitVar plot"};

  /// Correlation part
  ConfigurableAxis confMultBins{"confMultBins", {VARIABLE_WIDTH, 0.0f, 4.0f, 8.0f, 12.0f, 16.0f, 20.0f, 24.0f, 28.0f, 32.0f, 36.0f, 40.0f, 44.0f, 48.0f, 52.0f, 56.0f, 60.0f, 64.0f, 68.0f, 72.0f, 76.0f, 80.0f, 84.0f, 88.0f, 92.0f, 96.0f, 100.0f, 200.0f, 99999.f}, "Mixing bins - multiplicity"}; // \todo to be obtained from the hash task
  // ConfigurableAxis confMultBins{"CfgMultBins", {VARIABLE_WIDTH, 0.0f, 20.0f, 40.0f, 60.0f, 80.0f, 100.0f, 200.0f, 99999.f}, "Mixing bins - multiplicity"};
  ConfigurableAxis confVtxBins{"confVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};

  ConfigurableAxis confmTBins3D{"confmTBins3D", {VARIABLE_WIDTH, 1.02f, 1.14f, 1.20f, 1.26f, 1.38f, 1.56f, 1.86f, 4.50f}, "mT Binning for the 3Dimensional plot: k* vs multiplicity vs mT (set <<confUse3D>> to true in order to use)"};
  ConfigurableAxis confMultBins3D{"confMultBins3D", {VARIABLE_WIDTH, 0.0f, 20.0f, 30.0f, 40.0f, 99999.0f}, "multiplicity Binning for the 3Dimensional plot: k* vs multiplicity vs mT (set <<confUse3D>> to true in order to use)"};

  ColumnBinningPolicy<aod::collision::PosZ, aod::femtouniversecollision::MultNtr> colBinning{{confVtxBins, confMultBins}, true};

  ConfigurableAxis confkstarBins{"confkstarBins", {1500, 0., 6.}, "binning kstar"};
  ConfigurableAxis confkTBins{"confkTBins", {150, 0., 9.}, "binning kT"};
  ConfigurableAxis confmTBins{"confmTBins", {225, 0., 7.5}, "binning mT"};
  Configurable<int> confNEventsMix{"confNEventsMix", 5, "Number of events for mixing"};
  Configurable<bool> confIsCPR{"confIsCPR", true, "Close Pair Rejection"};
  Configurable<bool> confCPRPlotPerRadii{"confCPRPlotPerRadii", false, "Plot CPR per radii"};
  Configurable<int> confPhiBins{"confPhiBins", 29, "Number of phi bins in deta dphi"};
  Configurable<int> confEtaBins{"confEtaBins", 29, "Number of eta bins in deta dphi"};

  FemtoUniverseContainer<femto_universe_container::EventType::same, femto_universe_container::Observable::kstar> sameEventCont;
  FemtoUniverseContainer<femto_universe_container::EventType::mixed, femto_universe_container::Observable::kstar> mixedEventCont;
  FemtoUniversePairCleaner<aod::femtouniverseparticle::ParticleType::kMCTruthTrack, aod::femtouniverseparticle::ParticleType::kMCTruthTrack> pairCleaner;
  /// Histogram output
  HistogramRegistry qaRegistry{"TrackQA", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry resultRegistry{"Correlations", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry mixQaRegistry{"mixQaRegistry", {}, OutputObjHandlingPolicy::AnalysisObject};

  /// @brief Counter for particle swapping
  int fNeventsProcessed = 0;

  void init(InitContext&)
  {

    eventHisto.init(&qaRegistry);
    trackHistoPartOne.init(&qaRegistry, confTempFitVarpTBins, confTempFitVarPDGBins, 0, confPDGCodePartOne, false);
    if (!confIsSame) {
      trackHistoPartTwo.init(&qaRegistry, confTempFitVarpTBins, confTempFitVarPDGBins, 0, confPDGCodePartTwo, false);
    }

    mixQaRegistry.add("MixingQA/hSECollisionBins", ";bin;Entries", kTH1F, {{120, -0.5, 119.5}});
    mixQaRegistry.add("MixingQA/hMECollisionBins", ";bin;Entries", kTH1F, {{120, -0.5, 119.5}});

    sameEventCont.init(&resultRegistry, confkstarBins, confMultBins, confkTBins, confmTBins, confMultBins3D, confmTBins3D, confEtaBins, confPhiBins, 0, confUse3D);
    mixedEventCont.init(&resultRegistry, confkstarBins, confMultBins, confkTBins, confmTBins, confMultBins3D, confmTBins3D, confEtaBins, confPhiBins, 0, confUse3D);
    sameEventCont.setPDGCodes(confPDGCodePartOne, confPDGCodePartTwo);
    mixedEventCont.setPDGCodes(confPDGCodePartOne, confPDGCodePartTwo);
    pairCleaner.init(&qaRegistry);
  }

  template <typename CollisionType>
  void fillCollision(CollisionType col)
  {
    mixQaRegistry.fill(HIST("MixingQA/hSECollisionBins"), colBinning.getBin({col.posZ(), col.multNtr()}));
    eventHisto.fillQA(col);
  }

  /// This function processes the same event and takes care of all the histogramming
  /// \todo the trivial loops over the tracks should be factored out since they will be common to all combinations of T-T, T-V0, V0-V0, ...
  /// @tparam PartitionType
  /// @tparam PartType
  /// @tparam isMC: enables Monte Carlo truth specific histograms
  /// @param groupPartsOne partition for the first particle passed by the process function
  /// @param groupPartsTwo partition for the second particle passed by the process function
  /// @param parts femtoUniverseParticles table (in case of Monte Carlo joined with FemtoUniverseMCLabels)
  /// @param magFieldTesla magnetic field of the collision
  /// @param multCol multiplicity of the collision
  template <bool isMC, typename PartitionType, typename PartType>
  void doSameEvent(PartitionType groupPartsOne, PartitionType groupPartsTwo, PartType parts, float /*magFieldTesla*/, int multCol)
  {
    bool swpart = fNeventsProcessed % 2;
    fNeventsProcessed++;

    /// Histogramming same event
    for (auto const& part : groupPartsOne) {
      if (!confNoPDGPartOne && part.tempFitVar() != confPDGCodePartOne) {
        continue;
      }
      if (static_cast<int>(part.tempFitVar()) == static_cast<int>(Pdg::kD0) && confActiveD0OriginCheck && part.mLambda() != confD0OriginFlag) {
        continue;
      }
      trackHistoPartOne.fillQA<isMC, false>(part);
    }

    if (!confIsSame) {
      for (auto const& part : groupPartsTwo) {
        if (!confNoPDGPartTwo && part.tempFitVar() != confPDGCodePartTwo) {
          continue;
        }
        if (static_cast<int>(part.tempFitVar()) == static_cast<int>(Pdg::kD0) && confActiveD0OriginCheck && part.mLambda() != confD0OriginFlag) {
          continue;
        }
        trackHistoPartTwo.fillQA<isMC, false>(part);
      }
    }
    /// Now build the combinations
    if (!confIsSame) {
      // Build the combinations for pairs of non-identical particles
      for (auto const& [p1, p2] : combinations(CombinationsFullIndexPolicy(groupPartsOne, groupPartsTwo))) {
        // track cleaning
        if (!pairCleaner.isCleanPair(p1, p2, parts)) {
          continue;
        }
        if ((!confNoPDGPartOne && static_cast<int>(p1.tempFitVar()) != confPDGCodePartOne) || (!confNoPDGPartTwo && static_cast<int>(p2.tempFitVar()) != confPDGCodePartTwo)) {
          continue;
        }

        if (swpart) {
          sameEventCont.setPair<isMC>(p1, p2, multCol, confUse3D);
        } else {
          sameEventCont.setPair<isMC>(p2, p1, multCol, confUse3D);
        }

        swpart = !swpart;
      }
    } else {
      // Build the combinations for pairs of identical pairs
      for (auto const& [p1, p2] : combinations(CombinationsStrictlyUpperIndexPolicy(groupPartsOne, groupPartsTwo))) {
        // track cleaning
        if (!pairCleaner.isCleanPair(p1, p2, parts)) {
          continue;
        }
        if ((!confNoPDGPartOne && p2.tempFitVar() != confPDGCodePartOne) || (!confNoPDGPartTwo && p1.tempFitVar() != confPDGCodePartTwo)) {
          continue;
        }
        if (swpart)
          sameEventCont.setPair<isMC>(p1, p2, multCol, confUse3D);
        else
          sameEventCont.setPair<isMC>(p2, p1, multCol, confUse3D);

        swpart = !swpart;
      }
    }
  }
  /// process function for to call doSameEvent with Data
  /// \param col subscribe to the collision table (Data)
  /// \param parts subscribe to the femtoUniverseParticleTable
  void processSameEvent(o2::aod::FdCollision const& col,
                        o2::aod::FDParticles const& parts)
  {
    fillCollision(col);

    auto thegroupPartsOne = partsOne->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto thegroupPartsTwo = partsTwo->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);

    doSameEvent<false>(thegroupPartsOne, thegroupPartsTwo, parts, col.magField(), col.multNtr());
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackTrackMcTruth, processSameEvent, "Enable processing same event", true);

  /// This function processes the mixed event
  /// \todo the trivial loops over the collisions and tracks should be factored out since they will be common to all combinations of T-T, T-V0, V0-V0, ...
  /// \tparam PartitionType
  /// \tparam PartType
  /// \tparam isMC: enables Monte Carlo truth specific histograms
  /// \param groupPartsOne partition for the first particle passed by the process function
  /// \param groupPartsTwo partition for the second particle passed by the process function
  /// \param parts femtoUniverseParticles table (in case of Monte Carlo joined with FemtoUniverseMCLabels)
  /// \param magFieldTesla magnetic field of the collision
  /// \param multCol multiplicity of the collision
  template <bool isMC, typename PartitionType, typename PartType>
  void doMixedEvent(PartitionType groupPartsOne, PartitionType groupPartsTwo, PartType /*parts*/, float /*magFieldTesla*/, int multCol)
  {
    bool swpart = fNeventsProcessed % 2;
    fNeventsProcessed++;

    for (auto const& [p1, p2] : combinations(CombinationsFullIndexPolicy(groupPartsOne, groupPartsTwo))) {
      if ((!confNoPDGPartOne && static_cast<int>(p1.tempFitVar()) != confPDGCodePartOne) || (!confNoPDGPartTwo && static_cast<int>(p2.tempFitVar()) != confPDGCodePartTwo)) {
        continue;
      }
      if (swpart)
        mixedEventCont.setPair<isMC>(p1, p2, multCol, confUse3D);
      else
        mixedEventCont.setPair<isMC>(p2, p1, multCol, confUse3D);

      swpart = !swpart;
    }
  }

  /// process function for to call doMixedEvent with Data
  /// @param cols subscribe to the collisions table (Data)
  /// @param parts subscribe to the femtoUniverseParticleTable
  void processMixedEvent(o2::aod::FdCollisions const& cols,
                         o2::aod::FDParticles const& parts)
  {
    for (auto const& [collision1, collision2] : soa::selfCombinations(colBinning, 5, -1, cols, cols)) {

      const int multiplicityCol = collision1.multNtr();
      mixQaRegistry.fill(HIST("MixingQA/hMECollisionBins"), colBinning.getBin({collision1.posZ(), multiplicityCol}));

      auto groupPartsOne = partsOne->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);
      auto groupPartsTwo = partsTwo->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision2.globalIndex(), cache);

      const auto& magFieldTesla1 = collision1.magField();
      const auto& magFieldTesla2 = collision2.magField();

      if (magFieldTesla1 != magFieldTesla2) {
        continue;
      }
      /// \todo before mixing we should check whether both collisions contain a pair of particles!
      // if (partsOne.size() == 0 || nPart2Evt1 == 0 || nPart1Evt2 == 0 || partsTwo.size() == 0 ) continue;

      doMixedEvent<false>(groupPartsOne, groupPartsTwo, parts, magFieldTesla1, multiplicityCol);
    }
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackTrackMcTruth, processMixedEvent, "Enable processing mixed events", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<FemtoUniversePairTaskTrackTrackMcTruth>(cfgc),
  };
  return workflow;
}
