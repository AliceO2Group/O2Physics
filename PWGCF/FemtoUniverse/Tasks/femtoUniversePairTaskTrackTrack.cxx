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

/// \file femtoUniversePairTaskTrackTrack.cxx
/// \brief Tasks that reads the track tables used for the pairing and builds pairs of two tracks
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de
/// \author Georgios Mantzaridis, TU München, georgios.mantzaridis@tum.de
/// \author Anton Riedel, TU München, anton.riedel@tum.de
/// \author Zuzanna Chochulska, WUT Warsaw & CTU Prague, zchochul@cern.ch

#include "PWGCF/FemtoUniverse/Core/FemtoUniverseAngularContainer.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseDetaDphiStar.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseEventHisto.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseFemtoContainer.h"
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

#include <vector>

using namespace o2;
using namespace o2::analysis::femto_universe;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

namespace
{
static constexpr int NPart = 2;
static constexpr int NCuts = 5;
static const std::vector<std::string> partNames{"PartOne", "PartTwo"};
static const std::vector<std::string> cutNames{"MaxPt", "PIDthr", "nSigmaTPC", "nSigmaTPCTOF", "MaxP"};
static const float cutsTable[NPart][NCuts]{
  {4.05f, 1.f, 3.f, 3.f, 100.f},
  {4.05f, 1.f, 3.f, 3.f, 100.f}};
} // namespace

struct FemtoUniversePairTaskTrackTrack {
  SliceCache cache;
  Preslice<aod::FDParticles> perCol = aod::femtouniverseparticle::fdCollisionId;

  /// Particle selection part

  /// Table for both particles
  Configurable<LabeledArray<float>> confCutTable{"confCutTable", {cutsTable[0], NPart, NCuts, partNames, cutNames}, "Particle selections"};
  Configurable<int> confNspecies{"confNspecies", 2, "Number of particle spieces with PID info"};
  Configurable<bool> confIsMC{"confIsMC", false, "Enable additional Histogramms in the case of a MonteCarlo Run"};
  Configurable<std::vector<float>> confTrkPIDnSigmaMax{"confTrkPIDnSigmaMax", std::vector<float>{4.f, 3.f, 2.f}, "This configurable needs to be the same as the one used in the producer task"};
  Configurable<bool> confUse3D{"confUse3D", false, "Enable three dimensional histogramms (to be used only for analysis with high statistics): k* vs mT vs multiplicity"};

  /// Particle 1
  Configurable<int> confPDGCodePartOne{"confPDGCodePartOne", 2212, "Particle 1 - PDG code"};
  Configurable<uint32_t> confCutPartOne{"confCutPartOne", 5542474, "Particle 1 - Selection bit from cutCulator"};
  Configurable<int> confPIDPartOne{"confPIDPartOne", 2, "Particle 1 - Read from cutCulator"}; // we also need the possibility to specify whether the bit is true/false ->std>>vector<std::pair<int, int>>int>>

  /// Partition for particle 1
  Partition<aod::FDParticles> partsOne = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kTrack)) && ((aod::femtouniverseparticle::cut & confCutPartOne) == confCutPartOne);
  Partition<soa::Join<aod::FDParticles, aod::FDMCLabels>> partsOneMC = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kTrack)) && ((aod::femtouniverseparticle::cut & confCutPartOne) == confCutPartOne);

  /// Histogramming for particle 1
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kTrack, 1> trackHistoPartOne;

  /// Particle 2
  Configurable<bool> confIsSame{"confIsSame", false, "Pairs of the same particle"};
  Configurable<int> confPDGCodePartTwo{"confPDGCodePartTwo", 2212, "Particle 2 - PDG code"};
  Configurable<uint32_t> confCutPartTwo{"confCutPartTwo", 5542474, "Particle 2 - Selection bit"};
  Configurable<int> confPIDPartTwo{"confPIDPartTwo", 2, "Particle 2 - Read from cutCulator"}; // we also need the possibility to specify whether the bit is true/false ->std>>vector<std::pair<int, int>>

  /// Partition for particle 2
  Partition<aod::FDParticles> partsTwo = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kTrack)) &&
                                         //  (aod::femtouniverseparticle::pt < cfgCutTable->get("PartTwo", "MaxPt")) &&
                                         ((aod::femtouniverseparticle::cut & confCutPartTwo) == confCutPartTwo);
  Partition<soa::Join<aod::FDParticles, aod::FDMCLabels>> partsTwoMC = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kTrack)) &&
                                                                       ((aod::femtouniverseparticle::cut & confCutPartTwo) == confCutPartTwo);

  /// Histogramming for particle 2
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kTrack, 2> trackHistoPartTwo;

  /// Histogramming for Event
  FemtoUniverseEventHisto eventHisto;

  /// The configurables need to be passed to an std::vector
  int vPIDPartOne, vPIDPartTwo;
  std::vector<float> kNsigma;

  /// particle part
  ConfigurableAxis confTempFitVarBins{"confTempFitVarBins", {300, -0.15, 0.15}, "binning of the TempFitVar in the pT vs. TempFitVar plot"};
  ConfigurableAxis confTempFitVarpTBins{"confTempFitVarpTBins", {20, 0.5, 4.05}, "pT binning of the pT vs. TempFitVar plot"};

  /// Correlation part
  ConfigurableAxis confMultBins{"confMultBins", {VARIABLE_WIDTH, 0.0f, 4.0f, 8.0f, 12.0f, 16.0f, 20.0f, 24.0f, 28.0f, 32.0f, 36.0f, 40.0f, 44.0f, 48.0f, 52.0f, 56.0f, 60.0f, 64.0f, 68.0f, 72.0f, 76.0f, 80.0f, 84.0f, 88.0f, 92.0f, 96.0f, 100.0f, 200.0f, 99999.f}, "Mixing bins - multiplicity"}; // \todo to be obtained from the hash task
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
  Configurable<float> confCPRdeltaPhiCutMax{"confCPRdeltaPhiCutMax", 0.0, "Delta Phi max cut for Close Pair Rejection"};
  Configurable<float> confCPRdeltaPhiCutMin{"confCPRdeltaPhiCutMin", 0.0, "Delta Phi min cut for Close Pair Rejection"};
  Configurable<float> confCPRdeltaEtaCutMax{"confCPRdeltaEtaCutMax", 0.0, "Delta Eta max cut for Close Pair Rejection"};
  Configurable<float> confCPRdeltaEtaCutMin{"confCPRdeltaEtaCutMin", 0.0, "Delta Eta min cut for Close Pair Rejection"};
  Configurable<float> confCPRChosenRadii{"confCPRChosenRadii", 0.80, "Delta Eta cut for Close Pair Rejection"};
  Configurable<int> confPhiBins{"confPhiBins", 29, "Number of phi bins in deta dphi"};
  Configurable<int> confEtaBins{"confEtaBins", 29, "Number of eta bins in deta dphi"};

  FemtoUniverseFemtoContainer<femto_universe_femto_container::EventType::same, femto_universe_femto_container::Observable::kstar> sameEventFemtoCont;
  FemtoUniverseFemtoContainer<femto_universe_femto_container::EventType::mixed, femto_universe_femto_container::Observable::kstar> mixedEventFemtoCont;
  FemtoUniverseAngularContainer<femto_universe_angular_container::EventType::same, femto_universe_angular_container::Observable::kstar> sameEventAngularCont;
  FemtoUniverseAngularContainer<femto_universe_angular_container::EventType::mixed, femto_universe_angular_container::Observable::kstar> mixedEventAngularCont;
  FemtoUniversePairCleaner<aod::femtouniverseparticle::ParticleType::kTrack, aod::femtouniverseparticle::ParticleType::kTrack> pairCleaner;
  FemtoUniverseDetaDphiStar<aod::femtouniverseparticle::ParticleType::kTrack, aod::femtouniverseparticle::ParticleType::kTrack> pairCloseRejection;

  /// Histogram output
  HistogramRegistry qaRegistry{"TrackQA", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry resultRegistry{"Correlations", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry mixQaRegistry{"mixQaRegistry", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {

    eventHisto.init(&qaRegistry);
    trackHistoPartOne.init(&qaRegistry, confTempFitVarpTBins, confTempFitVarBins, confIsMC, confPDGCodePartOne, false);
    if (!confIsSame) {
      trackHistoPartTwo.init(&qaRegistry, confTempFitVarpTBins, confTempFitVarBins, confIsMC, confPDGCodePartTwo, false);
    }

    mixQaRegistry.add("MixingQA/hSECollisionBins", ";bin;Entries", kTH1F, {{120, -0.5, 119.5}});
    mixQaRegistry.add("MixingQA/hMECollisionBins", ";bin;Entries", kTH1F, {{120, -0.5, 119.5}});

    sameEventFemtoCont.init(&resultRegistry, confkstarBins, confMultBins, confkTBins, confmTBins, confMultBins3D, confmTBins3D, confIsMC, confUse3D);
    mixedEventFemtoCont.init(&resultRegistry, confkstarBins, confMultBins, confkTBins, confmTBins, confMultBins3D, confmTBins3D, confIsMC, confUse3D);
    sameEventAngularCont.init(&resultRegistry, confkstarBins, confMultBins, confkTBins, confmTBins, confMultBins3D, confmTBins3D, confEtaBins, confPhiBins, confIsMC, confUse3D);
    mixedEventAngularCont.init(&resultRegistry, confkstarBins, confMultBins, confkTBins, confmTBins, confMultBins3D, confmTBins3D, confEtaBins, confPhiBins, confIsMC, confUse3D);

    sameEventFemtoCont.setPDGCodes(confPDGCodePartOne, confPDGCodePartTwo);
    mixedEventFemtoCont.setPDGCodes(confPDGCodePartOne, confPDGCodePartTwo);
    sameEventAngularCont.setPDGCodes(confPDGCodePartOne, confPDGCodePartTwo);
    mixedEventAngularCont.setPDGCodes(confPDGCodePartOne, confPDGCodePartTwo);

    pairCleaner.init(&qaRegistry);
    if (confIsCPR.value) {
      pairCloseRejection.init(&resultRegistry, &qaRegistry, confCPRdeltaPhiCutMin.value, confCPRdeltaPhiCutMax.value, confCPRdeltaEtaCutMin.value, confCPRdeltaEtaCutMax.value, confCPRChosenRadii.value, confCPRPlotPerRadii.value);
    }

    vPIDPartOne = confPIDPartOne.value;
    vPIDPartTwo = confPIDPartTwo.value;
    kNsigma = confTrkPIDnSigmaMax.value;
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
  void doSameEvent(PartitionType groupPartsOne, PartitionType groupPartsTwo, PartType parts, float magFieldTesla, int multCol)
  {

    /// Histogramming same event
    for (const auto& part : groupPartsOne) {
      if (part.p() > confCutTable->get("PartOne", "MaxP") || part.pt() > confCutTable->get("PartOne", "MaxPt")) {
        continue;
      }
      if (!isFullPIDSelected(part.pidCut(),
                             part.p(),
                             confCutTable->get("PartOne", "PIDthr"),
                             vPIDPartOne,
                             confNspecies,
                             kNsigma,
                             confCutTable->get("PartOne", "nSigmaTPC"),
                             confCutTable->get("PartOne", "nSigmaTPCTOF"))) {
        continue;
      }

      trackHistoPartOne.fillQA<isMC, false>(part);
    }

    if (!confIsSame) {
      for (const auto& part : groupPartsTwo) {
        if (part.p() > confCutTable->get("PartTwo", "MaxP") || part.pt() > confCutTable->get("PartTwo", "MaxPt")) {
          continue;
        }
        if (!isFullPIDSelected(part.pidCut(),
                               part.p(),
                               confCutTable->get("PartTwo", "PIDthr"),
                               vPIDPartTwo,
                               confNspecies,
                               kNsigma,
                               confCutTable->get("PartTwo", "nSigmaTPC"),
                               confCutTable->get("PartTwo", "nSigmaTPCTOF"))) {
          continue;
        }
        trackHistoPartTwo.fillQA<isMC, false>(part);
      }
    }
    /// Now build the combinations
    for (const auto& [p1, p2] : combinations(CombinationsStrictlyUpperIndexPolicy(groupPartsOne, groupPartsTwo))) {
      if (p1.p() > confCutTable->get("PartOne", "MaxP") || p1.pt() > confCutTable->get("PartOne", "MaxPt") || p2.p() > confCutTable->get("PartTwo", "MaxP") || p2.pt() > confCutTable->get("PartTwo", "MaxPt")) {
        continue;
      }
      if (!isFullPIDSelected(p1.pidCut(),
                             p1.p(),
                             confCutTable->get("PartOne", "PIDthr"),
                             vPIDPartOne,
                             confNspecies,
                             kNsigma,
                             confCutTable->get("PartOne", "nSigmaTPC"),
                             confCutTable->get("PartOne", "nSigmaTPCTOF")) ||
          !isFullPIDSelected(p2.pidCut(),
                             p2.p(),
                             confCutTable->get("PartTwo", "PIDthr"),
                             vPIDPartTwo,
                             confNspecies,
                             kNsigma,
                             confCutTable->get("PartTwo", "nSigmaTPC"),
                             confCutTable->get("PartTwo", "nSigmaTPCTOF"))) {
        continue;
      }

      if (confIsCPR.value) {
        if (pairCloseRejection.isClosePair(p1, p2, parts, magFieldTesla, femto_universe_container::EventType::same)) {
          continue;
        }
      }

      // track cleaning
      if (!pairCleaner.isCleanPair(p1, p2, parts)) {
        continue;
      }

      sameEventFemtoCont.setPair<isMC>(p1, p2, multCol, confUse3D);
      sameEventAngularCont.setPair<isMC>(p1, p2, multCol, confUse3D);
    }
  }

  /// process function for to call doSameEvent with Data
  /// \param col subscribe to the collision table (Data)
  /// \param parts subscribe to the femtoUniverseParticleTable
  void processSameEvent(const o2::aod::FdCollision& col,
                        const o2::aod::FDParticles& parts)
  {
    fillCollision(col);

    auto thegroupPartsOne = partsOne->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto thegroupPartsTwo = partsTwo->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);

    doSameEvent<false>(thegroupPartsOne, thegroupPartsTwo, parts, col.magField(), col.multNtr());
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackTrack, processSameEvent, "Enable processing same event", true);

  /// process function for to call doSameEvent with Monte Carlo
  /// \param col subscribe to the collision table (Monte Carlo Reconstructed reconstructed)
  /// \param parts subscribe to joined table FemtoUniverseParticles and FemtoUniverseMCLables to access Monte Carlo truth
  /// \param FemtoUniverseMCParticles subscribe to the Monte Carlo truth table
  void processSameEventMC(const o2::aod::FdCollision& col,
                          const soa::Join<o2::aod::FDParticles, o2::aod::FDMCLabels>& parts,
                          const o2::aod::FdMCParticles&)
  {
    fillCollision(col);

    auto thegroupPartsOne = partsOneMC->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto thegroupPartsTwo = partsTwoMC->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);

    doSameEvent<true>(thegroupPartsOne, thegroupPartsTwo, parts, col.magField(), col.multNtr());
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackTrack, processSameEventMC, "Enable processing same event for Monte Carlo", false);

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
  void doMixedEvent(PartitionType groupPartsOne, PartitionType groupPartsTwo, PartType parts, float magFieldTesla, int multCol)
  {

    for (const auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(groupPartsOne, groupPartsTwo))) {
      if (p1.p() > confCutTable->get("PartOne", "MaxP") || p1.pt() > confCutTable->get("PartOne", "MaxPt") || p2.p() > confCutTable->get("PartTwo", "MaxP") || p2.pt() > confCutTable->get("PartTwo", "MaxPt")) {
        continue;
      }
      if (!isFullPIDSelected(p1.pidCut(),
                             p1.p(),
                             confCutTable->get("PartOne", "PIDthr"),
                             vPIDPartOne,
                             confNspecies,
                             kNsigma,
                             confCutTable->get("PartOne", "nSigmaTPC"),
                             confCutTable->get("PartOne", "nSigmaTPCTOF")) ||
          !isFullPIDSelected(p2.pidCut(),
                             p2.p(),
                             confCutTable->get("PartTwo", "PIDthr"),
                             vPIDPartTwo,
                             confNspecies,
                             kNsigma,
                             confCutTable->get("PartTwo", "nSigmaTPC"),
                             confCutTable->get("PartTwo", "nSigmaTPCTOF"))) {
        continue;
      }

      if (confIsCPR.value) {
        if (pairCloseRejection.isClosePair(p1, p2, parts, magFieldTesla, femto_universe_container::EventType::mixed)) {
          continue;
        }
      }

      mixedEventFemtoCont.setPair<isMC>(p1, p2, multCol, confUse3D);
      mixedEventAngularCont.setPair<isMC>(p1, p2, multCol, confUse3D);
    }
  }

  /// process function for to call doMixedEvent with Data
  /// @param cols subscribe to the collisions table (Data)
  /// @param parts subscribe to the femtoUniverseParticleTable
  void processMixedEvent(const o2::aod::FdCollisions& cols,
                         const o2::aod::FDParticles& parts)
  {
    for (const auto& [collision1, collision2] : soa::selfCombinations(colBinning, 5, -1, cols, cols)) {

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
  PROCESS_SWITCH(FemtoUniversePairTaskTrackTrack, processMixedEvent, "Enable processing mixed events", true);

  /// brief process function for to call doMixedEvent with Monte Carlo
  /// @param cols subscribe to the collisions table (Monte Carlo Reconstructed reconstructed)
  /// @param parts subscribe to joined table FemtoUniverseParticles and FemtoUniverseMCLables to access Monte Carlo truth
  /// @param FemtoUniverseMCParticles subscribe to the Monte Carlo truth table
  void processMixedEventMC(const o2::aod::FdCollisions& cols,
                           const soa::Join<o2::aod::FDParticles, o2::aod::FDMCLabels>& parts,
                           const o2::aod::FdMCParticles&)
  {
    for (const auto& [collision1, collision2] : soa::selfCombinations(colBinning, 5, -1, cols, cols)) {

      const int multiplicityCol = collision1.multNtr();
      mixQaRegistry.fill(HIST("MixingQA/hMECollisionBins"), colBinning.getBin({collision1.posZ(), multiplicityCol}));

      auto groupPartsOne = partsOneMC->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);
      auto groupPartsTwo = partsTwoMC->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision2.globalIndex(), cache);

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
  PROCESS_SWITCH(FemtoUniversePairTaskTrackTrack, processMixedEventMC, "Enable processing mixed events MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<FemtoUniversePairTaskTrackTrack>(cfgc),
  };
  return workflow;
}
