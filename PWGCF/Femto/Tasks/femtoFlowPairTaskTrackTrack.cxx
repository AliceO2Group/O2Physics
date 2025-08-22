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
/// \author Wenya Wu, TUM, wenya.wu@cern.ch
/// \NOTE:  The femtoflow framework borrows and copies the framework of femtouniverse and femtodream

#include <vector>
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/StepTHn.h"
#include "Framework/O2DatabasePDGPlugin.h"

#include "PWGCF/Femto/DataModel/FemtoDerived.h"
#include "PWGCF/Femto/Core/FemtoFlowParticleHisto.h"
#include "PWGCF/Femto/Core/FemtoFlowEventHisto.h"
#include "PWGCF/Femto/Core/FemtoFlowPairCleaner.h"
#include "PWGCF/Femto/Core/FemtoFlowFemtoContainer.h"
#include "PWGCF/Femto/Core/FemtoFlowAngularContainer.h"
#include "PWGCF/Femto/Core/FemtoFlowDetaDphiStar.h"
#include "PWGCF/Femto/Core/femtoUtils.h"
#include "PWGCF/Femto/Core/FemtoFlowTrackSelection.h"

using namespace o2;
using namespace o2::analysis::femto_flow;
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

struct FemtoFlowPairTaskTrackTrack {
  SliceCache cache;
  Preslice<aod::FDParticles> perCol = aod::femtoflowparticle::fdCollisionId;

  using FemtoFullParticles = soa::Join<aod::FDParticles, aod::FDExtParticles>;
  Preslice<FemtoFullParticles> perColReco = aod::femtoflowparticle::fdCollisionId;

  /// Particle selection part

  /// Table for both particles
  Configurable<LabeledArray<float>> confCutTable{"confCutTable", {cutsTable[0], NPart, NCuts, partNames, cutNames}, "Particle selections"};
  Configurable<int> confNspecies{"confNspecies", 1, "Number of particle spieces with PID info"};
  Configurable<bool> confIsMC{"confIsMC", false, "Enable additional Histogramms in the case of a MonteCarlo Run"};
  Configurable<std::vector<float>> confTrkPIDnSigmaMax{"confTrkPIDnSigmaMax", std::vector<float>{3.5f, 3.f, 2.5f}, "This configurable needs to be the same as the one used in the producer task"};
  Configurable<bool> confUse3D{"confUse3D", false, "Enable three dimensional histogramms (to be used only for analysis with high statistics): k* vs mT vs multiplicity"};
  Configurable<bool> confUseqnDivide{"confUseqnDivide", true, "Enable histogramms of k* vs mT in separated qn bins"};
  Configurable<bool> confIsDebug{"confIsDebug", true, "Enable Debug tables"};

  /// Particle 1
  Configurable<int> confPDGCodePartOne{"confPDGCodePartOne", 2212, "Particle 1 - PDG code"};
  Configurable<uint32_t> confCutPartOne{"confCutPartOne", 170220070, "Particle 1 - Selection bit from cutCulator"};
  Configurable<int> confPIDPartOne{"confPIDPartOne", 0, "Particle 1 - Read from cutCulator"}; // Should be bit-mask from cutCulator
  // we also need the possibility to specify whether the bit is true/false ->std>>vector<std::pair<int, int>>int>>

  /// Partition for particle 1
  Partition<aod::FDParticles> partsOne = (aod::femtoflowparticle::partType == uint8_t(aod::femtoflowparticle::ParticleType::kTrack)) && ((aod::femtoflowparticle::cut & confCutPartOne) == confCutPartOne);
  Partition<soa::Join<aod::FDParticles, aod::FDMCLabels>> partsOneMC = (aod::femtoflowparticle::partType == uint8_t(aod::femtoflowparticle::ParticleType::kTrack)) && ((aod::femtoflowparticle::cut & confCutPartOne) == confCutPartOne);

  /// Histogramming for particle 1
  FemtoFlowParticleHisto<aod::femtoflowparticle::ParticleType::kTrack, 1> trackHistoPartOne;

  /// Particle 2
  Configurable<bool> confIsSame{"confIsSame", false, "Pairs of the same particle"};
  Configurable<int> confPDGCodePartTwo{"confPDGCodePartTwo", 2212, "Particle 2 - PDG code"};
  Configurable<uint32_t> confCutPartTwo{"confCutPartTwo", 170220070, "Particle 2 - Selection bit"};
  Configurable<int> confPIDPartTwo{"confPIDPartTwo", 0, "Particle 2 - Read from cutCulator"}; // we also need the possibility to specify whether the bit is true/false ->std>>vector<std::pair<int, int>>

  /// Partition for particle 2
  Partition<aod::FDParticles> partsTwo = (aod::femtoflowparticle::partType == uint8_t(aod::femtoflowparticle::ParticleType::kTrack)) &&
                                         //  (aod::femtoflowparticle::pt < cfgCutTable->get("PartTwo", "MaxPt")) &&
                                         ((aod::femtoflowparticle::cut & confCutPartTwo) == confCutPartTwo);
  Partition<soa::Join<aod::FDParticles, aod::FDMCLabels>> partsTwoMC = (aod::femtoflowparticle::partType == uint8_t(aod::femtoflowparticle::ParticleType::kTrack)) &&
                                                                       ((aod::femtoflowparticle::cut & confCutPartTwo) == confCutPartTwo);

  /// Histogramming for particle 2
  FemtoFlowParticleHisto<aod::femtoflowparticle::ParticleType::kTrack, 2> trackHistoPartTwo;

  // Particle for debug
  Partition<FemtoFullParticles> partsDebug = (aod::femtoflowparticle::partType == uint8_t(aod::femtoflowparticle::ParticleType::kTrack)) && 
                                             (aod::femtoflowparticle::sign == int8_t(1)) &&
                                             ((aod::femtoflowparticle::cut & confCutPartOne) == confCutPartOne); // Bit_mask for partOne or partTwo

  /// Histogramming for debug particle 
  FemtoFlowParticleHisto<aod::femtoflowparticle::ParticleType::kTrack> trackHistoPartDebug; //FolderSuffixType name set as 0 = "_debug"

  /// Histogramming for Event
  FemtoFlowEventHisto eventHisto;

  /// The configurables need to be passed to an std::vector
  int vPIDPartOne, vPIDPartTwo;
  std::vector<float> kNsigma;

  /// particle part
  ConfigurableAxis confTempFitVarBins{"confTempFitVarBins", {300, -0.15, 0.15}, "binning of the TempFitVar in the pT vs. TempFitVar plot"};
  ConfigurableAxis confTempFitVarpTBins{"confTempFitVarpTBins", {20, 0.5, 4.05}, "pT binning of the pT vs. TempFitVar plot"};

  /// Correlation part
  ConfigurableAxis confMultBins{"confMultBins", {VARIABLE_WIDTH, 0.0f, 50.0f, 100.0f, 200.0f, 300.0f, 500.0f, 800.0f, 1000.0f, 1500.0f, 2000.0f, 2500.0f, 3000.0f, 4000.f}, "Mixing bins - multiplicity"}; // \todo to be obtained from the hash task
  ConfigurableAxis confVtxBins{"confVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};

  ConfigurableAxis confmTBins3D{"confmTBins3D", {VARIABLE_WIDTH, 1.02f, 1.14f, 1.20f, 1.26f, 1.38f, 1.56f, 1.86f, 4.50f}, "mT Binning for the 3Dimensional plot: k* vs multiplicity vs mT (set <<confUse3D>> to true in order to use)"};
  ConfigurableAxis confMultBins3D{"confMultBins3D", {VARIABLE_WIDTH, 0.0f, 20.0f, 30.0f, 40.0f, 99999.0f}, "multiplicity Binning for the 3Dimensional plot: k* vs multiplicity vs mT (set <<confUse3D>> to true in order to use)"};

  ColumnBinningPolicy<aod::collision::PosZ, aod::femtoflowcollision::MultNtr> colBinning{{confVtxBins, confMultBins}, true};

  ConfigurableAxis confkstarBins{"confkstarBins", {1500, 0., 6.}, "binning kstar"};
  ConfigurableAxis confkTBins{"confkTBins", {150, 0., 9.}, "binning kT"};
  ConfigurableAxis confmTBins{"confmTBins", {225, 0., 7.5}, "binning mT"};
  Configurable<int> confNEventsMix{"confNEventsMix", 5, "Number of events for mixing"};
  Configurable<bool> confIsCPR{"confIsCPR", false, "Close Pair Rejection"};
  Configurable<bool> confCPRPlotPerRadii{"confCPRPlotPerRadii", false, "Plot CPR per radii"};
  Configurable<float> confCPRdeltaPhiCutMax{"confCPRdeltaPhiCutMax", 0.0, "Delta Phi max cut for Close Pair Rejection"};
  Configurable<float> confCPRdeltaPhiCutMin{"confCPRdeltaPhiCutMin", 0.0, "Delta Phi min cut for Close Pair Rejection"};
  Configurable<float> confCPRdeltaEtaCutMax{"confCPRdeltaEtaCutMax", 0.0, "Delta Eta max cut for Close Pair Rejection"};
  Configurable<float> confCPRdeltaEtaCutMin{"confCPRdeltaEtaCutMin", 0.0, "Delta Eta min cut for Close Pair Rejection"};
  Configurable<float> confCPRChosenRadii{"confCPRChosenRadii", 0.01, "Delta Eta cut for Close Pair Rejection"};
  Configurable<int> confPhiBins{"confPhiBins", 29, "Number of phi bins in deta dphi"};
  Configurable<int> confEtaBins{"confEtaBins", 29, "Number of eta bins in deta dphi"};

  FemtoFlowFemtoContainer<femto_flow_femto_container::EventType::same, femto_flow_femto_container::Observable::kstar> sameEventFemtoCont;
  FemtoFlowFemtoContainer<femto_flow_femto_container::EventType::mixed, femto_flow_femto_container::Observable::kstar> mixedEventFemtoCont;
  FemtoFlowAngularContainer<femto_flow_angular_container::EventType::same, femto_flow_angular_container::Observable::kstar> sameEventAngularCont;
  FemtoFlowAngularContainer<femto_flow_angular_container::EventType::mixed, femto_flow_angular_container::Observable::kstar> mixedEventAngularCont;
  FemtoFlowPairCleaner<aod::femtoflowparticle::ParticleType::kTrack, aod::femtoflowparticle::ParticleType::kTrack> pairCleaner;
  FemtoFlowDetaDphiStar<aod::femtoflowparticle::ParticleType::kTrack, aod::femtoflowparticle::ParticleType::kTrack> pairCloseRejection;

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

    if (confIsDebug) trackHistoPartDebug.init(&qaRegistry, confTempFitVarpTBins, confTempFitVarBins, confIsMC, confPDGCodePartOne, true);

    mixQaRegistry.add("MixingQA/hSECollisionBins", ";bin;Entries", kTH1F, {{120, -0.5, 119.5}});
    mixQaRegistry.add("MixingQA/hMECollisionBins", ";bin;Entries", kTH1F, {{120, -0.5, 119.5}});

    sameEventFemtoCont.init(&resultRegistry, confkstarBins, confMultBins, confkTBins, confmTBins, confMultBins3D, confmTBins3D, confIsMC, confUse3D, confUseqnDivide);
    mixedEventFemtoCont.init(&resultRegistry, confkstarBins, confMultBins, confkTBins, confmTBins, confMultBins3D, confmTBins3D, confIsMC, confUse3D, false);
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
  /// @param parts femtoUniverseParticles table (in case of Monte Carlo joined with FemtoFlowMCLabels)
  /// @param magFieldTesla magnetic field of the collision
  /// @param multCol multiplicity of the collision
  template <bool isMC, typename PartitionType, typename PartType>
  void doSameEvent(PartitionType groupPartsOne, PartitionType groupPartsTwo, PartType parts, float magFieldTesla, int multCol, int myqnBin)
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
        if (pairCloseRejection.isClosePair(p1, p2, parts, magFieldTesla, femto_flow_femto_container::EventType::same)) {
          continue;
        }
      }

      // track cleaning
      if (!pairCleaner.isCleanPair(p1, p2, parts)) {
        continue;
      }
      sameEventFemtoCont.setPair<isMC>(p1, p2, multCol, confUse3D, confUseqnDivide, myqnBin);
      sameEventAngularCont.setPair<isMC>(p1, p2, multCol, confUse3D);
    }
  }

  /// process function for to call doSameEvent with Data
  /// \param col subscribe to the collision table (Data)
  /// \param parts subscribe to the femtoUniverseParticleTable
  void processSameEvent(const o2::aod::FDCollision& col,
                        const o2::aod::FDParticles& parts,
                        FemtoFullParticles&)
  {
    fillCollision(col);

    auto thegroupPartsOne = partsOne->sliceByCached(aod::femtoflowparticle::fdCollisionId, col.globalIndex(), cache);
    auto thegroupPartsTwo = partsTwo->sliceByCached(aod::femtoflowparticle::fdCollisionId, col.globalIndex(), cache);

    doSameEvent<false>(thegroupPartsOne, thegroupPartsTwo, parts, col.magField(), col.multNtr(), col.qnbin());
  
    if(confIsDebug){   
      auto thegroupPartsDebug = partsDebug->sliceByCached(aod::femtoflowparticle::fdCollisionId, col.globalIndex(), cache);
      for (const auto& partDebug : thegroupPartsDebug) {
        if (partDebug.p() > confCutTable->get("PartTwo", "MaxP") || partDebug.pt() > confCutTable->get("PartTwo", "MaxPt")) {
          continue;
        }
        if (!isFullPIDSelected(partDebug.pidCut(),
                               partDebug.p(),
                               confCutTable->get("PartOne", "PIDthr"),
                               vPIDPartOne,
                               confNspecies,
                               kNsigma,
                               confCutTable->get("PartOne", "nSigmaTPC"),
                               confCutTable->get("PartOne", "nSigmaTPCTOF"))) {
          continue;
        }
        trackHistoPartDebug.fillQA<false, true>(partDebug);
      }
    }
  }
  PROCESS_SWITCH(FemtoFlowPairTaskTrackTrack, processSameEvent, "Enable processing same event", true);

  /// process function for to call doSameEvent with Monte Carlo
  /// \param col subscribe to the collision table (Monte Carlo Reconstructed reconstructed)
  /// \param parts subscribe to joined table FemtoFlowParticles and FemtoFlowMCLables to access Monte Carlo truth
  /// \param FemtoFlowMCParticles subscribe to the Monte Carlo truth table
  void processSameEventMC(const o2::aod::FDCollision& col,
                          const soa::Join<o2::aod::FDParticles, o2::aod::FDMCLabels>& parts,
                          const o2::aod::FDMCParticles&)
  {
    fillCollision(col);

    auto thegroupPartsOne = partsOneMC->sliceByCached(aod::femtoflowparticle::fdCollisionId, col.globalIndex(), cache);
    auto thegroupPartsTwo = partsTwoMC->sliceByCached(aod::femtoflowparticle::fdCollisionId, col.globalIndex(), cache);

    doSameEvent<true>(thegroupPartsOne, thegroupPartsTwo, parts, col.magField(), col.multNtr(), col.qnbin());
  }
  PROCESS_SWITCH(FemtoFlowPairTaskTrackTrack, processSameEventMC, "Enable processing same event for Monte Carlo", false);

  /// This function processes the mixed event
  /// \todo the trivial loops over the collisions and tracks should be factored out since they will be common to all combinations of T-T, T-V0, V0-V0, ...
  /// \tparam PartitionType
  /// \tparam PartType
  /// \tparam isMC: enables Monte Carlo truth specific histograms
  /// \param groupPartsOne partition for the first particle passed by the process function
  /// \param groupPartsTwo partition for the second particle passed by the process function
  /// \param parts femtoUniverseParticles table (in case of Monte Carlo joined with FemtoFlowMCLabels)
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
        if (pairCloseRejection.isClosePair(p1, p2, parts, magFieldTesla, femto_flow_femto_container::EventType::mixed)) {
          continue;
        }
      }

      mixedEventFemtoCont.setPair<isMC>(p1, p2, multCol, confUse3D, false, -999);
      mixedEventAngularCont.setPair<isMC>(p1, p2, multCol, confUse3D);
    }
  }

  /// process function for to call doMixedEvent with Data
  /// @param cols subscribe to the collisions table (Data)
  /// @param parts subscribe to the femtoFlowParticleTable
  void processMixedEvent(const o2::aod::FDCollisions& cols,
                         const o2::aod::FDParticles& parts)
  {
    for (const auto& [collision1, collision2] : soa::selfCombinations(colBinning, 5, -1, cols, cols)) {

      const int multiplicityCol = collision1.multNtr();
      mixQaRegistry.fill(HIST("MixingQA/hMECollisionBins"), colBinning.getBin({collision1.posZ(), multiplicityCol}));

      auto groupPartsOne = partsOne->sliceByCached(aod::femtoflowparticle::fdCollisionId, collision1.globalIndex(), cache);
      auto groupPartsTwo = partsTwo->sliceByCached(aod::femtoflowparticle::fdCollisionId, collision2.globalIndex(), cache);

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
  PROCESS_SWITCH(FemtoFlowPairTaskTrackTrack, processMixedEvent, "Enable processing mixed events", true);

  /// brief process function for to call doMixedEvent with Monte Carlo
  /// @param cols subscribe to the collisions table (Monte Carlo Reconstructed reconstructed)
  /// @param parts subscribe to joined table FemtoFlowParticles and FemtoFlowMCLables to access Monte Carlo truth
  /// @param FemtoFlowMCParticles subscribe to the Monte Carlo truth table
  void processMixedEventMC(const o2::aod::FDCollisions& cols,
                           const soa::Join<o2::aod::FDParticles, o2::aod::FDMCLabels>& parts,
                           const o2::aod::FDMCParticles&)
  {
    for (const auto& [collision1, collision2] : soa::selfCombinations(colBinning, 5, -1, cols, cols)) {

      const int multiplicityCol = collision1.multNtr();
      mixQaRegistry.fill(HIST("MixingQA/hMECollisionBins"), colBinning.getBin({collision1.posZ(), multiplicityCol}));

      auto groupPartsOne = partsOneMC->sliceByCached(aod::femtoflowparticle::fdCollisionId, collision1.globalIndex(), cache);
      auto groupPartsTwo = partsTwoMC->sliceByCached(aod::femtoflowparticle::fdCollisionId, collision2.globalIndex(), cache);

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
  PROCESS_SWITCH(FemtoFlowPairTaskTrackTrack, processMixedEventMC, "Enable processing mixed events MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<FemtoFlowPairTaskTrackTrack>(cfgc),
  };
  return workflow;
}
