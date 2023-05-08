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

/// \file femtoDreamPairTaskTrackTrack.cxx
/// \brief Tasks that reads the track tables used for the pairing and builds pairs of two tracks
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de
/// \author Georgios Mantzaridis, TU München, georgios.mantzaridis@tum.de
/// \author Anton Riedel, TU München, anton.riedel@tum.de

#include <vector>
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/StepTHn.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "TDatabasePDG.h"

#include "PWGCF/DataModel/FemtoDerived.h"
#include "FemtoDreamParticleHisto.h"
#include "FemtoDreamEventHisto.h"
#include "FemtoDreamPairCleaner.h"
#include "FemtoDreamContainer.h"
#include "FemtoDreamDetaDphiStar.h"
#include "FemtoUtils.h"

using namespace o2;
using namespace o2::analysis::femtoDream;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

namespace
{
static constexpr int nPart = 2;
static constexpr int nCuts = 5;
static const std::vector<std::string> partNames{"PartOne", "PartTwo"};
static const std::vector<std::string> cutNames{"MaxPt", "PIDthr", "nSigmaTPC", "nSigmaTPCTOF", "MaxP"};
static const float cutsTable[nPart][nCuts]{
  {4.05f, 1.f, 3.f, 3.f, 100.f},
  {4.05f, 1.f, 3.f, 3.f, 100.f}};
} // namespace

struct femtoDreamPairTaskTrackTrack {
  Service<O2DatabasePDG> pdg;
  SliceCache cache;
  Preslice<aod::FemtoDreamDebugParticles> perCol = aod::femtodreamparticle::femtoDreamCollisionId;

  /// Particle selection part

  /// Table for both particles
  Configurable<LabeledArray<float>> cfgCutTable{"cfgCutTable", {cutsTable[0], nPart, nCuts, partNames, cutNames}, "Particle selections"};
  Configurable<int> cfgNspecies{"ccfgNspecies", 4, "Number of particle spieces with PID info"};
  Configurable<std::vector<float>> ConfPIDnSigmaMax{"ConfPIDnSigmaMax", std::vector<float>{3.5f, 3.f, 2.5f}, "This configurable needs to be the same as the one used in the producer task"};
  Configurable<bool> ConfIsMC{"ConfIsMC", false, "Enable additional Histogramms in the case of a MonteCarlo Run"};

  /// Particle 1
  Configurable<int> ConfPDGCodePartOne{"ConfPDGCodePartOne", 2212, "Particle 1 - PDG code"};
  Configurable<uint32_t> ConfCutPartOne{"ConfCutPartOne", 5542474, "Particle 1 - Selection bit from cutCulator"};
  Configurable<std::vector<int>> ConfPIDPartOne{"ConfPIDPartOne", std::vector<int>{2}, "Particle 1 - Read from cutCulator"}; // we also need the possibility to specify whether the bit is true/false ->std>>vector<std::pair<int, int>>int>>

  /// Partition for particle 1
  Partition<aod::FemtoDreamParticles> partsOne = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack)) && ((aod::femtodreamparticle::cut & ConfCutPartOne) == ConfCutPartOne);
  Partition<soa::Join<aod::FemtoDreamParticles, aod::FemtoDreamMCLabels>> partsOneMC = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack)) && ((aod::femtodreamparticle::cut & ConfCutPartOne) == ConfCutPartOne);

  /// Histogramming for particle 1
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kTrack, 1> trackHistoPartOne;

  /// Particle 2
  Configurable<bool> ConfIsSame{"ConfIsSame", false, "Pairs of the same particle"};
  Configurable<int> ConfPDGCodePartTwo{"ConfPDGCodePartTwo", 2212, "Particle 2 - PDG code"};
  Configurable<uint32_t> ConfCutPartTwo{"ConfCutPartTwo", 5542474, "Particle 2 - Selection bit"};
  Configurable<std::vector<int>> ConfPIDPartTwo{"ConfPIDPartTwo", std::vector<int>{2}, "Particle 2 - Read from cutCulator"}; // we also need the possibility to specify whether the bit is true/false ->std>>vector<std::pair<int, int>>

  /// Partition for particle 2
  Partition<aod::FemtoDreamParticles> partsTwo = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack)) &&
                                                 //  (aod::femtodreamparticle::pt < cfgCutTable->get("PartTwo", "MaxPt")) &&
                                                 ((aod::femtodreamparticle::cut & ConfCutPartTwo) == ConfCutPartTwo);
  Partition<soa::Join<aod::FemtoDreamParticles, aod::FemtoDreamMCLabels>> partsTwoMC = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack)) &&
                                                                                       ((aod::femtodreamparticle::cut & ConfCutPartTwo) == ConfCutPartTwo);

  /// Histogramming for particle 2
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kTrack, 2> trackHistoPartTwo;

  /// Histogramming for Event
  FemtoDreamEventHisto eventHisto;

  /// The configurables need to be passed to an std::vector
  std::vector<int> vPIDPartOne, vPIDPartTwo;
  std::vector<float> kNsigma;

  /// particle part
  ConfigurableAxis CfgTempFitVarBins{"CfgDTempFitVarBins", {300, -0.15, 0.15}, "binning of the TempFitVar in the pT vs. TempFitVar plot"};
  ConfigurableAxis CfgTempFitVarpTBins{"CfgTempFitVarpTBins", {20, 0.5, 4.05}, "pT binning of the pT vs. TempFitVar plot"};

  /// Correlation part
  ConfigurableAxis CfgMultBins{"CfgMultBins", {VARIABLE_WIDTH, 0.0f, 4.0f, 8.0f, 12.0f, 16.0f, 20.0f, 24.0f, 28.0f, 32.0f, 36.0f, 40.0f, 44.0f, 48.0f, 52.0f, 56.0f, 60.0f, 64.0f, 68.0f, 72.0f, 76.0f, 80.0f, 84.0f, 88.0f, 92.0f, 96.0f, 100.0f, 200.0f, 99999.f}, "Mixing bins - multiplicity"}; // \todo to be obtained from the hash task
  // ConfigurableAxis CfgMultBins{"CfgMultBins", {VARIABLE_WIDTH, 0.0f, 20.0f, 40.0f, 60.0f, 80.0f, 100.0f, 200.0f, 99999.f}, "Mixing bins - multiplicity"};
  ConfigurableAxis CfgVtxBins{"CfgVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};

  ColumnBinningPolicy<aod::collision::PosZ, aod::femtodreamcollision::MultNtr> colBinning{{CfgVtxBins, CfgMultBins}, true};

  ConfigurableAxis CfgkstarBins{"CfgkstarBins", {1500, 0., 6.}, "binning kstar"};
  ConfigurableAxis CfgkTBins{"CfgkTBins", {150, 0., 9.}, "binning kT"};
  ConfigurableAxis CfgmTBins{"CfgmTBins", {225, 0., 7.5}, "binning mT"};
  Configurable<int> ConfNEventsMix{"ConfNEventsMix", 5, "Number of events for mixing"};
  Configurable<bool> ConfIsCPR{"ConfIsCPR", true, "Close Pair Rejection"};
  Configurable<bool> ConfCPRPlotPerRadii{"ConfCPRPlotPerRadii", false, "Plot CPR per radii"};

  FemtoDreamContainer<femtoDreamContainer::EventType::same, femtoDreamContainer::Observable::kstar> sameEventCont;
  FemtoDreamContainer<femtoDreamContainer::EventType::mixed, femtoDreamContainer::Observable::kstar> mixedEventCont;
  FemtoDreamPairCleaner<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kTrack> pairCleaner;
  FemtoDreamDetaDphiStar<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kTrack> pairCloseRejection;
  /// Histogram output
  HistogramRegistry qaRegistry{"TrackQA", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry resultRegistry{"Correlations", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry MixQaRegistry{"MixQaRegistry", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
    eventHisto.init(&qaRegistry);
    trackHistoPartOne.init(&qaRegistry, CfgTempFitVarpTBins, CfgTempFitVarBins, ConfIsMC);
    if (!ConfIsSame) {
      trackHistoPartTwo.init(&qaRegistry, CfgTempFitVarpTBins, CfgTempFitVarBins, ConfIsMC);
    }

    MixQaRegistry.add("MixingQA/hSECollisionBins", ";bin;Entries", kTH1F, {{120, -0.5, 119.5}});
    MixQaRegistry.add("MixingQA/hMECollisionBins", ";bin;Entries", kTH1F, {{120, -0.5, 119.5}});

    sameEventCont.init(&resultRegistry, CfgkstarBins, CfgMultBins, CfgkTBins, CfgmTBins, ConfIsMC);
    mixedEventCont.init(&resultRegistry, CfgkstarBins, CfgMultBins, CfgkTBins, CfgmTBins, ConfIsMC);
    sameEventCont.setPDGCodes(ConfPDGCodePartOne, ConfPDGCodePartTwo);
    mixedEventCont.setPDGCodes(ConfPDGCodePartOne, ConfPDGCodePartTwo);
    pairCleaner.init(&qaRegistry);
    if (ConfIsCPR) {
      pairCloseRejection.init(&resultRegistry, &qaRegistry, 0.01, 0.01, ConfCPRPlotPerRadii); /// \todo add config for Δη and ΔΦ cut values
    }

    vPIDPartOne = ConfPIDPartOne;
    vPIDPartTwo = ConfPIDPartTwo;
    kNsigma = ConfPIDnSigmaMax;
  }

  template <typename CollisionType>
  void fillCollision(CollisionType col)
  {
    MixQaRegistry.fill(HIST("MixingQA/hSECollisionBins"), colBinning.getBin({col.posZ(), col.multNtr()}));
    eventHisto.fillQA(col);
  }

  /// This function processes the same event and takes care of all the histogramming
  /// \todo the trivial loops over the tracks should be factored out since they will be common to all combinations of T-T, T-V0, V0-V0, ...
  /// @tparam PartitionType
  /// @tparam PartType
  /// @tparam isMC: enables Monte Carlo truth specific histograms
  /// @param groupPartsOne partition for the first particle passed by the process function
  /// @param groupPartsTwo partition for the second particle passed by the process function
  /// @param parts femtoDreamParticles table (in case of Monte Carlo joined with FemtoDreamMCLabels)
  /// @param magFieldTesla magnetic field of the collision
  /// @param multCol multiplicity of the collision
  template <bool isMC, typename PartitionType, typename PartType>
  void doSameEvent(PartitionType groupPartsOne, PartitionType groupPartsTwo, PartType parts, float magFieldTesla, int multCol)
  {

    /// Histogramming same event
    for (auto& part : groupPartsOne) {
      if (part.p() > cfgCutTable->get("PartOne", "MaxP") || part.pt() > cfgCutTable->get("PartOne", "MaxPt")) {
        continue;
      }
      if (!isFullPIDSelected(part.pidcut(),
                             part.p(),
                             cfgCutTable->get("PartOne", "PIDthr"),
                             vPIDPartOne,
                             cfgNspecies,
                             kNsigma,
                             cfgCutTable->get("PartOne", "nSigmaTPC"),
                             cfgCutTable->get("PartOne", "nSigmaTPCTOF"))) {
        continue;
      }

      trackHistoPartOne.fillQA<isMC>(part);
    }

    if (!ConfIsSame) {
      for (auto& part : groupPartsTwo) {
        if (part.p() > cfgCutTable->get("PartTwo", "MaxP") || part.pt() > cfgCutTable->get("PartTwo", "MaxPt")) {
          continue;
        }
        if (!isFullPIDSelected(part.pidcut(),
                               part.p(),
                               cfgCutTable->get("PartTwo", "PIDthr"),
                               vPIDPartTwo,
                               cfgNspecies,
                               kNsigma,
                               cfgCutTable->get("PartTwo", "nSigmaTPC"),
                               cfgCutTable->get("PartTwo", "nSigmaTPCTOF"))) {
          continue;
        }
        trackHistoPartTwo.fillQA<isMC>(part);
      }
    }
    /// Now build the combinations
    for (auto& [p1, p2] : combinations(CombinationsStrictlyUpperIndexPolicy(groupPartsOne, groupPartsTwo))) {
      if (p1.p() > cfgCutTable->get("PartOne", "MaxP") || p1.pt() > cfgCutTable->get("PartOne", "MaxPt") || p2.p() > cfgCutTable->get("PartTwo", "MaxP") || p2.pt() > cfgCutTable->get("PartTwo", "MaxPt")) {
        continue;
      }
      if (!isFullPIDSelected(p1.pidcut(),
                             p1.p(),
                             cfgCutTable->get("PartOne", "PIDthr"),
                             vPIDPartOne,
                             cfgNspecies,
                             kNsigma,
                             cfgCutTable->get("PartOne", "nSigmaTPC"),
                             cfgCutTable->get("PartOne", "nSigmaTPCTOF")) ||
          !isFullPIDSelected(p2.pidcut(),
                             p2.p(),
                             cfgCutTable->get("PartTwo", "PIDthr"),
                             vPIDPartTwo,
                             cfgNspecies,
                             kNsigma,
                             cfgCutTable->get("PartTwo", "nSigmaTPC"),
                             cfgCutTable->get("PartTwo", "nSigmaTPCTOF"))) {
        continue;
      }

      if (ConfIsCPR) {
        if (pairCloseRejection.isClosePair(p1, p2, parts, magFieldTesla)) {
          continue;
        }
      }

      // track cleaning
      if (!pairCleaner.isCleanPair(p1, p2, parts)) {
        continue;
      }
      sameEventCont.setPair<isMC>(p1, p2, multCol);
    }
  }

  /// process function for to call doSameEvent with Data
  /// \param col subscribe to the collision table (Data)
  /// \param parts subscribe to the femtoDreamParticleTable
  void processSameEvent(o2::aod::FemtoDreamCollision& col,
                        o2::aod::FemtoDreamParticles& parts)
  {
    fillCollision(col);

    auto thegroupPartsOne = partsOne->sliceByCached(aod::femtodreamparticle::femtoDreamCollisionId, col.globalIndex(), cache);
    auto thegroupPartsTwo = partsTwo->sliceByCached(aod::femtodreamparticle::femtoDreamCollisionId, col.globalIndex(), cache);

    doSameEvent<false>(thegroupPartsOne, thegroupPartsTwo, parts, col.magField(), col.multNtr());
  }
  PROCESS_SWITCH(femtoDreamPairTaskTrackTrack, processSameEvent, "Enable processing same event", true);

  /// process function for to call doSameEvent with Monte Carlo
  /// \param col subscribe to the collision table (Monte Carlo Reconstructed reconstructed)
  /// \param parts subscribe to joined table FemtoDreamParticles and FemtoDreamMCLables to access Monte Carlo truth
  /// \param FemtoDreamMCParticles subscribe to the Monte Carlo truth table
  void processSameEventMC(o2::aod::FemtoDreamCollision& col,
                          soa::Join<o2::aod::FemtoDreamParticles, o2::aod::FemtoDreamMCLabels>& parts,
                          o2::aod::FemtoDreamMCParticles&)
  {
    fillCollision(col);

    auto thegroupPartsOne = partsOneMC->sliceByCached(aod::femtodreamparticle::femtoDreamCollisionId, col.globalIndex(), cache);
    auto thegroupPartsTwo = partsTwoMC->sliceByCached(aod::femtodreamparticle::femtoDreamCollisionId, col.globalIndex(), cache);

    doSameEvent<true>(thegroupPartsOne, thegroupPartsTwo, parts, col.magField(), col.multNtr());
  }
  PROCESS_SWITCH(femtoDreamPairTaskTrackTrack, processSameEventMC, "Enable processing same event for Monte Carlo", false);

  /// This function processes the mixed event
  /// \todo the trivial loops over the collisions and tracks should be factored out since they will be common to all combinations of T-T, T-V0, V0-V0, ...
  /// \tparam PartitionType
  /// \tparam PartType
  /// \tparam isMC: enables Monte Carlo truth specific histograms
  /// \param groupPartsOne partition for the first particle passed by the process function
  /// \param groupPartsTwo partition for the second particle passed by the process function
  /// \param parts femtoDreamParticles table (in case of Monte Carlo joined with FemtoDreamMCLabels)
  /// \param magFieldTesla magnetic field of the collision
  /// \param multCol multiplicity of the collision
  template <bool isMC, typename PartitionType, typename PartType>
  void doMixedEvent(PartitionType groupPartsOne, PartitionType groupPartsTwo, PartType parts, float magFieldTesla, int multCol)
  {

    for (auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(groupPartsOne, groupPartsTwo))) {
      if (p1.p() > cfgCutTable->get("PartOne", "MaxP") || p1.pt() > cfgCutTable->get("PartOne", "MaxPt") || p2.p() > cfgCutTable->get("PartTwo", "MaxP") || p2.pt() > cfgCutTable->get("PartTwo", "MaxPt")) {
        continue;
      }
      if (!isFullPIDSelected(p1.pidcut(),
                             p1.p(),
                             cfgCutTable->get("PartOne", "PIDthr"),
                             vPIDPartOne,
                             cfgNspecies,
                             kNsigma,
                             cfgCutTable->get("PartOne", "nSigmaTPC"),
                             cfgCutTable->get("PartOne", "nSigmaTPCTOF")) ||
          !isFullPIDSelected(p2.pidcut(),
                             p2.p(),
                             cfgCutTable->get("PartTwo", "PIDthr"),
                             vPIDPartTwo,
                             cfgNspecies,
                             kNsigma,
                             cfgCutTable->get("PartTwo", "nSigmaTPC"),
                             cfgCutTable->get("PartTwo", "nSigmaTPCTOF"))) {
        continue;
      }

      if (ConfIsCPR) {
        if (pairCloseRejection.isClosePair(p1, p2, parts, magFieldTesla)) {
          continue;
        }
      }

      mixedEventCont.setPair<isMC>(p1, p2, multCol);
    }
  }

  /// process function for to call doMixedEvent with Data
  /// @param cols subscribe to the collisions table (Data)
  /// @param parts subscribe to the femtoDreamParticleTable
  void processMixedEvent(o2::aod::FemtoDreamCollisions& cols,
                         o2::aod::FemtoDreamParticles& parts)
  {
    for (auto& [collision1, collision2] : soa::selfCombinations(colBinning, 5, -1, cols, cols)) {

      const int multiplicityCol = collision1.multNtr();
      MixQaRegistry.fill(HIST("MixingQA/hMECollisionBins"), colBinning.getBin({collision1.posZ(), multiplicityCol}));

      auto groupPartsOne = partsOne->sliceByCached(aod::femtodreamparticle::femtoDreamCollisionId, collision1.globalIndex(), cache);
      auto groupPartsTwo = partsTwo->sliceByCached(aod::femtodreamparticle::femtoDreamCollisionId, collision2.globalIndex(), cache);

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
  PROCESS_SWITCH(femtoDreamPairTaskTrackTrack, processMixedEvent, "Enable processing mixed events", true);

  /// brief process function for to call doMixedEvent with Monte Carlo
  /// @param cols subscribe to the collisions table (Monte Carlo Reconstructed reconstructed)
  /// @param parts subscribe to joined table FemtoDreamParticles and FemtoDreamMCLables to access Monte Carlo truth
  /// @param FemtoDreamMCParticles subscribe to the Monte Carlo truth table
  void processMixedEventMC(o2::aod::FemtoDreamCollisions& cols,
                           soa::Join<o2::aod::FemtoDreamParticles, o2::aod::FemtoDreamMCLabels>& parts,
                           o2::aod::FemtoDreamMCParticles&)
  {
    for (auto& [collision1, collision2] : soa::selfCombinations(colBinning, 5, -1, cols, cols)) {

      const int multiplicityCol = collision1.multNtr();
      MixQaRegistry.fill(HIST("MixingQA/hMECollisionBins"), colBinning.getBin({collision1.posZ(), multiplicityCol}));

      auto groupPartsOne = partsOneMC->sliceByCached(aod::femtodreamparticle::femtoDreamCollisionId, collision1.globalIndex(), cache);
      auto groupPartsTwo = partsTwoMC->sliceByCached(aod::femtodreamparticle::femtoDreamCollisionId, collision2.globalIndex(), cache);

      const auto& magFieldTesla1 = collision1.magField();
      const auto& magFieldTesla2 = collision2.magField();

      if (magFieldTesla1 != magFieldTesla2) {
        continue;
      }
      /// \todo before mixing we should check whether both collisions contain a pair of particles!
      // if (partsOne.size() == 0 || nPart2Evt1 == 0 || nPart1Evt2 == 0 || partsTwo.size() == 0 ) continue;

      doMixedEvent<true>(groupPartsOne, groupPartsTwo, parts, magFieldTesla1, multiplicityCol);
    }
  }
  PROCESS_SWITCH(femtoDreamPairTaskTrackTrack, processMixedEventMC, "Enable processing mixed events MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<femtoDreamPairTaskTrackTrack>(cfgc),
  };
  return workflow;
}
