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

/// \file femtoDreamPairTaskTrackTrackTrack.cxx
/// \brief Tasks that reads the track tables and creates track triplets; should be currently used only for identical tracks
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de
/// \author Georgios Mantzaridis, TU München, georgios.mantzaridis@tum.de
/// \author Anton Riedel, TU München, anton.riedel@tum.de
/// \author Laura Serksnyte, TU München, laura.serksnyte@tum.de

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
#include "FemtoDreamContainerThreeBody.h"
#include "FemtoDreamDetaDphiStar.h"
#include "FemtoUtils.h"

using namespace o2;
using namespace o2::analysis::femtoDream;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

namespace
{
static constexpr int nPart = 3;
static constexpr int nCuts = 5;
static const std::vector<std::string> partNames{"PartOne", "PartTwo", "PartThree"};
static const std::vector<std::string> cutNames{"MaxPt", "PIDthr", "nSigmaTPC", "nSigmaTPCTOF", "MaxP"};
static const float cutsTable[nPart][nCuts]{
  {4.05f, 1.f, 3.f, 3.f, 100.f},
  {4.05f, 1.f, 3.f, 3.f, 100.f},
  {4.05f, 1.f, 3.f, 3.f, 100.f}};
} // namespace

struct femtoDreamPairTaskTrackTrackTrack {
  SliceCache cache;
  Preslice<aod::FDParticles> perCol = aod::femtodreamparticle::fdCollisionId;

  /// Particle selection part

  /// Table for both particles
  Configurable<LabeledArray<float>> ConfCutTable{"ConfCutTable", {cutsTable[0], nPart, nCuts, partNames, cutNames}, "Particle selections"};
  Configurable<int> ConfNspecies{"ConfNspecies", 2, "Number of particle spieces with PID info"};
  Configurable<bool> ConfIsMC{"ConfIsMC", false, "Enable additional Histogramms in the case of a MonteCarlo Run"};
  Configurable<std::vector<float>> ConfTrkPIDnSigmaMax{"ConfTrkPIDnSigmaMax", std::vector<float>{4.f, 3.f, 2.f}, "This configurable needs to be the same as the one used in the producer task"};
  Configurable<bool> ConfUse3D{"ConfUse3D", false, "Enable three dimensional histogramms (to be used only for analysis with high statistics): k* vs mT vs multiplicity"};

  // Which particles to analyse; currently support only for same species and cuts triplets
  Configurable<int> ConfPDGCodePart{"ConfPDGCodePart", 2212, "Particle PDG code"};
  Configurable<int> ConfPIDPart{"ConfPIDPart", 2, "Particles - Read from cutCulator"};
  Configurable<uint32_t> ConfCutPart{"ConfCutPart", 5542474, "Particles - Selection bit from cutCulator"};

  // Keep this format as next code update will include also different particle species and cuts:

  /// Particle 1

  /// Partition for particle 1
  Partition<aod::FDParticles> partsOne = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack)) && ((aod::femtodreamparticle::cut & ConfCutPart) == ConfCutPart);
  Partition<soa::Join<aod::FDParticles, aod::FDMCLabels>> partsOneMC = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack)) && ((aod::femtodreamparticle::cut & ConfCutPart) == ConfCutPart);

  /// Histogramming for particle 1
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kTrack, 1> trackHistoPartOne;

  /// Particle 2

  /// Partition for particle 2
  Partition<aod::FDParticles> partsTwo = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack)) &&
                                         //  (aod::femtodreamparticle::pt < cfgCutTable->get("PartTwo", "MaxPt")) &&
                                         ((aod::femtodreamparticle::cut & ConfCutPart) == ConfCutPart);
  Partition<soa::Join<aod::FDParticles, aod::FDMCLabels>> partsTwoMC = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack)) &&
                                                                       ((aod::femtodreamparticle::cut & ConfCutPart) == ConfCutPart);

  /// Particle 3

  /// Partition for particle 3
  Partition<aod::FDParticles> partsThree = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack)) &&
                                           //  (aod::femtodreamparticle::pt < cfgCutTable->get("PartTwo", "MaxPt")) &&
                                           ((aod::femtodreamparticle::cut & ConfCutPart) == ConfCutPart);
  Partition<soa::Join<aod::FDParticles, aod::FDMCLabels>> partsThreeMC = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack)) &&
                                                                         ((aod::femtodreamparticle::cut & ConfCutPart) == ConfCutPart);

  /// Histogramming for Event
  FemtoDreamEventHisto eventHisto;

  /// The configurables need to be passed to an std::vector
  int vPIDPart;
  std::vector<float> kNsigma;

  /// particle part
  ConfigurableAxis ConfTempFitVarBins{"ConfDTempFitVarBins", {300, -0.15, 0.15}, "binning of the TempFitVar in the pT vs. TempFitVar plot"};
  ConfigurableAxis ConfTempFitVarpTBins{"ConfTempFitVarpTBins", {20, 0.5, 4.05}, "pT binning of the pT vs. TempFitVar plot"};

  /// Correlation part
  // ConfigurableAxis ConfMultBins{"ConfMultBins", {VARIABLE_WIDTH, 0.0f, 4.0f, 8.0f, 12.0f, 16.0f, 20.0f, 24.0f, 28.0f, 32.0f, 36.0f, 40.0f, 44.0f, 48.0f, 52.0f, 56.0f, 60.0f, 64.0f, 68.0f, 72.0f, 76.0f, 80.0f, 84.0f, 88.0f, 92.0f, 96.0f, 100.0f, 200.0f, 99999.f}, "Mixing bins - multiplicity"}; // \todo to be obtained from the hash task
  ConfigurableAxis ConfMultBins{"CfgMultBins", {VARIABLE_WIDTH, 0.0f, 20.0f, 40.0f, 60.0f, 80.0f, 100.0f, 200.0f, 99999.f}, "Mixing bins - multiplicity"};
  ConfigurableAxis ConfVtxBins{"ConfVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};

  ColumnBinningPolicy<aod::collision::PosZ, aod::femtodreamcollision::MultNtr> colBinning{{ConfVtxBins, ConfMultBins}, true};

  ConfigurableAxis ConfQ3Bins{"ConfQ3Bins", {2000, 0., 8.}, "binning Q3"};
  Configurable<int> ConfNEventsMix{"ConfNEventsMix", 5, "Number of events for mixing"};
  Configurable<bool> ConfIsCPR{"ConfIsCPR", true, "Close Pair Rejection"};
  Configurable<bool> ConfCPRPlotPerRadii{"ConfCPRPlotPerRadii", false, "Plot CPR per radii"};
  Configurable<float> ConfCPRdeltaPhiMax{"ConfCPRdeltaPhiMax", 0.01, "Max. Delta Phi for Close Pair Rejection"};
  Configurable<float> ConfCPRdeltaEtaMax{"ConfCPRdeltaEtaMax", 0.01, "Max. Delta Eta for Close Pair Rejection"};
  ConfigurableAxis ConfDummy{"ConfDummy", {1, 0, 1}, "Dummy axis"};

  FemtoDreamContainerThreeBody<femtoDreamContainerThreeBody::EventType::same, femtoDreamContainerThreeBody::Observable::Q3> sameEventCont;
  FemtoDreamContainerThreeBody<femtoDreamContainerThreeBody::EventType::mixed, femtoDreamContainerThreeBody::Observable::Q3> mixedEventCont;
  FemtoDreamPairCleaner<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kTrack> pairCleaner;
  FemtoDreamDetaDphiStar<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kTrack> pairCloseRejection;
  /// Histogram output
  HistogramRegistry qaRegistry{"TrackQA", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry resultRegistry{"Correlations", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry MixQaRegistry{"MixQaRegistry", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {

    eventHisto.init(&qaRegistry);
    trackHistoPartOne.init(&qaRegistry, ConfTempFitVarpTBins, ConfTempFitVarBins, ConfDummy, ConfDummy, ConfDummy, ConfDummy, ConfIsMC, ConfPDGCodePart);

    MixQaRegistry.add("MixingQA/hSECollisionBins", ";bin;Entries", kTH1F, {{120, -0.5, 119.5}});
    MixQaRegistry.add("MixingQA/hMECollisionBins", ";bin;Entries", kTH1F, {{120, -0.5, 119.5}});

    sameEventCont.init(&resultRegistry, ConfQ3Bins, ConfMultBins, ConfIsMC);
    mixedEventCont.init(&resultRegistry, ConfQ3Bins, ConfMultBins, ConfIsMC);
    sameEventCont.setPDGCodes(ConfPDGCodePart, ConfPDGCodePart, ConfPDGCodePart);
    mixedEventCont.setPDGCodes(ConfPDGCodePart, ConfPDGCodePart, ConfPDGCodePart);
    pairCleaner.init(&qaRegistry); // SERKSNYTE : later check if init should be updated to have 3 separate histos
    if (ConfIsCPR.value) {
      pairCloseRejection.init(&resultRegistry, &qaRegistry, ConfCPRdeltaPhiMax.value, ConfCPRdeltaEtaMax.value, ConfCPRPlotPerRadii.value);
    }

    // CURRENTLY do only for one species
    vPIDPart = ConfPIDPart.value;
    kNsigma = ConfTrkPIDnSigmaMax.value;
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
  void doSameEvent(PartitionType groupPartsOne, PartitionType groupPartsTwo, PartitionType groupPartsThree, PartType parts, float magFieldTesla, int multCol)
  {
    /// Histogramming same event
    int NGoodTracks = 0;
    for (auto& part : groupPartsOne) {
      if (part.p() > ConfCutTable->get("PartOne", "MaxP") || part.pt() > ConfCutTable->get("PartOne", "MaxPt")) {
        continue;
      }
      if (!isFullPIDSelected(part.pidcut(),
                             part.p(),
                             ConfCutTable->get("PartOne", "PIDthr"),
                             vPIDPart,
                             ConfNspecies,
                             kNsigma,
                             ConfCutTable->get("PartOne", "nSigmaTPC"),
                             ConfCutTable->get("PartOne", "nSigmaTPCTOF"))) {
        continue;
      }
      NGoodTracks++;

      trackHistoPartOne.fillQA<isMC, false>(part, aod::femtodreamparticle::kPt);
    }

    /// Now build the combinations
    for (auto& [p1, p2, p3] : combinations(CombinationsStrictlyUpperIndexPolicy(groupPartsOne, groupPartsTwo, groupPartsThree))) {
      if (p1.p() > ConfCutTable->get("PartOne", "MaxP") || p1.pt() > ConfCutTable->get("PartOne", "MaxPt") || p2.p() > ConfCutTable->get("PartTwo", "MaxP") || p2.pt() > ConfCutTable->get("PartTwo", "MaxPt") || p3.p() > ConfCutTable->get("PartThree", "MaxP") || p3.pt() > ConfCutTable->get("PartThree", "MaxPt")) {
        continue;
      }
      if (!isFullPIDSelected(p1.pidcut(),
                             p1.p(),
                             ConfCutTable->get("PartOne", "PIDthr"),
                             vPIDPart,
                             ConfNspecies,
                             kNsigma,
                             ConfCutTable->get("PartOne", "nSigmaTPC"),
                             ConfCutTable->get("PartOne", "nSigmaTPCTOF")) ||
          !isFullPIDSelected(p2.pidcut(),
                             p2.p(),
                             ConfCutTable->get("PartTwo", "PIDthr"),
                             vPIDPart,
                             ConfNspecies,
                             kNsigma,
                             ConfCutTable->get("PartTwo", "nSigmaTPC"),
                             ConfCutTable->get("PartTwo", "nSigmaTPCTOF")) ||
          !isFullPIDSelected(p3.pidcut(),
                             p3.p(),
                             ConfCutTable->get("PartThree", "PIDthr"),
                             vPIDPart,
                             ConfNspecies,
                             kNsigma,
                             ConfCutTable->get("PartThree", "nSigmaTPC"),
                             ConfCutTable->get("PartThree", "nSigmaTPCTOF"))) {
        continue;
      }

      if (ConfIsCPR.value) {
        if (pairCloseRejection.isClosePair(p1, p2, parts, magFieldTesla)) {
          continue;
        }
        if (pairCloseRejection.isClosePair(p2, p3, parts, magFieldTesla)) {
          continue;
        }
        if (pairCloseRejection.isClosePair(p1, p3, parts, magFieldTesla)) {
          continue;
        }
      }

      // track cleaning
      if (!pairCleaner.isCleanPair(p1, p2, parts)) {
        continue;
      }
      if (!pairCleaner.isCleanPair(p2, p3, parts)) {
        continue;
      }
      if (!pairCleaner.isCleanPair(p1, p3, parts)) {
        continue;
      }
      sameEventCont.setTriplet<isMC>(p1, p2, p3, multCol);
    }
  }

  /// process function to call doSameEvent with Data
  /// \param col subscribe to the collision table (Data)
  /// \param parts subscribe to the femtoDreamParticleTable
  void processSameEvent(o2::aod::FDCollision& col,
                        o2::aod::FDParticles& parts)
  {
    fillCollision(col);

    auto thegroupPartsOne = partsOne->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    auto thegroupPartsTwo = partsTwo->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    auto thegroupPartsThree = partsThree->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);

    doSameEvent<false>(thegroupPartsOne, thegroupPartsTwo, thegroupPartsTwo, parts, col.magField(), col.multNtr());
  }
  PROCESS_SWITCH(femtoDreamPairTaskTrackTrackTrack, processSameEvent, "Enable processing same event", true);

  /// process function for to call doSameEvent with Monte Carlo
  /// \param col subscribe to the collision table (Monte Carlo Reconstructed reconstructed)
  /// \param parts subscribe to joined table FemtoDreamParticles and FemtoDreamMCLables to access Monte Carlo truth
  /// \param FemtoDreamMCParticles subscribe to the Monte Carlo truth table
  void processSameEventMC(o2::aod::FDCollision& col,
                          soa::Join<o2::aod::FDParticles, o2::aod::FDMCLabels>& parts,
                          o2::aod::FDMCParticles&)
  {
    fillCollision(col);

    auto thegroupPartsOne = partsOneMC->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    auto thegroupPartsTwo = partsTwoMC->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    auto thegroupPartsThree = partsThreeMC->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);

    doSameEvent<true>(thegroupPartsOne, thegroupPartsTwo, thegroupPartsThree, parts, col.magField(), col.multNtr());
  }
  PROCESS_SWITCH(femtoDreamPairTaskTrackTrackTrack, processSameEventMC, "Enable processing same event for Monte Carlo", false);

  /// This function processes the mixed event
  /// \todo the trivial loops over the collisions and tracks should be factored out since they will be common to all combinations of T-T, T-V0, V0-V0, ...
  /// \tparam PartitionType
  /// \tparam PartType
  /// \tparam isMC: enables Monte Carlo truth specific histograms
  /// \param groupPartsOne partition for the first particle passed by the process function
  /// \param groupPartsTwo partition for the second particle passed by the process function
  /// \param groupPartsThree partition for the third particle passed by the process function
  /// \param parts femtoDreamParticles table (in case of Monte Carlo joined with FemtoDreamMCLabels)
  /// \param magFieldTesla magnetic field of the collision
  /// \param multCol multiplicity of the collision
  template <bool isMC, typename PartitionType, typename PartType>
  void doMixedEvent(PartitionType groupPartsOne, PartitionType groupPartsTwo, PartitionType groupPartsThree, PartType parts, float magFieldTesla, int multCol)
  {
    for (auto& [p1, p2, p3] : combinations(CombinationsFullIndexPolicy(groupPartsOne, groupPartsTwo, groupPartsThree))) {
      if (p1.p() > ConfCutTable->get("PartOne", "MaxP") || p1.pt() > ConfCutTable->get("PartOne", "MaxPt") || p2.p() > ConfCutTable->get("PartTwo", "MaxP") || p2.pt() > ConfCutTable->get("PartTwo", "MaxPt") || p3.p() > ConfCutTable->get("PartThree", "MaxP") || p3.pt() > ConfCutTable->get("PartThree", "MaxPt")) {
        continue;
      }

      if (!isFullPIDSelected(p1.pidcut(),
                             p1.p(),
                             ConfCutTable->get("PartOne", "PIDthr"),
                             vPIDPart,
                             ConfNspecies,
                             kNsigma,
                             ConfCutTable->get("PartOne", "nSigmaTPC"),
                             ConfCutTable->get("PartOne", "nSigmaTPCTOF")) ||
          !isFullPIDSelected(p2.pidcut(),
                             p2.p(),
                             ConfCutTable->get("PartTwo", "PIDthr"),
                             vPIDPart,
                             ConfNspecies,
                             kNsigma,
                             ConfCutTable->get("PartTwo", "nSigmaTPC"),
                             ConfCutTable->get("PartTwo", "nSigmaTPCTOF")) ||
          !isFullPIDSelected(p3.pidcut(),
                             p3.p(),
                             ConfCutTable->get("PartThree", "PIDthr"),
                             vPIDPart,
                             ConfNspecies,
                             kNsigma,
                             ConfCutTable->get("PartThree", "nSigmaTPC"),
                             ConfCutTable->get("PartThree", "nSigmaTPCTOF"))) {
        continue;
      }

      if (ConfIsCPR.value) {
        if (pairCloseRejection.isClosePair(p1, p2, parts, magFieldTesla)) {
          continue;
        }
        if (pairCloseRejection.isClosePair(p2, p3, parts, magFieldTesla)) {
          continue;
        }

        if (pairCloseRejection.isClosePair(p1, p3, parts, magFieldTesla)) {
          continue;
        }
      }
      mixedEventCont.setTriplet<isMC>(p1, p2, p3, multCol);
    }
  }

  /// process function for to call doMixedEvent with Data
  /// @param cols subscribe to the collisions table (Data)
  /// @param parts subscribe to the femtoDreamParticleTable
  void processMixedEvent(o2::aod::FDCollisions& cols,
                         o2::aod::FDParticles& parts)
  {
    for (auto& [collision1, collision2, collision3] : soa::selfCombinations(colBinning, 5, -1, cols, cols, cols)) {
      const int multiplicityCol = collision1.multNtr();
      MixQaRegistry.fill(HIST("MixingQA/hMECollisionBins"), colBinning.getBin({collision1.posZ(), multiplicityCol}));

      auto groupPartsOne = partsOne->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision1.globalIndex(), cache);
      auto groupPartsTwo = partsTwo->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision2.globalIndex(), cache);
      auto groupPartsThree = partsThree->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision3.globalIndex(), cache);

      const auto& magFieldTesla1 = collision1.magField();
      const auto& magFieldTesla2 = collision2.magField();
      const auto& magFieldTesla3 = collision3.magField();

      if ((magFieldTesla1 != magFieldTesla2) || (magFieldTesla2 != magFieldTesla3) || (magFieldTesla1 != magFieldTesla3)) {
        continue;
      }
      // CONSIDER testing different strategies to which events to use

      doMixedEvent<false>(groupPartsOne, groupPartsTwo, groupPartsThree, parts, magFieldTesla1, multiplicityCol);
    }
  }
  PROCESS_SWITCH(femtoDreamPairTaskTrackTrackTrack, processMixedEvent, "Enable processing mixed events", true);

  /// brief process function for to call doMixedEvent with Monte Carlo
  /// @param cols subscribe to the collisions table (Monte Carlo Reconstructed reconstructed)
  /// @param parts subscribe to joined table FemtoDreamParticles and FemtoDreamMCLables to access Monte Carlo truth
  /// @param FemtoDreamMCParticles subscribe to the Monte Carlo truth table
  void processMixedEventMC(o2::aod::FDCollisions& cols,
                           soa::Join<o2::aod::FDParticles, o2::aod::FDMCLabels>& parts,
                           o2::aod::FDMCParticles&)
  {
    for (auto& [collision1, collision2, collision3] : soa::selfCombinations(colBinning, 5, -1, cols, cols, cols)) {

      const int multiplicityCol = collision1.multNtr();
      MixQaRegistry.fill(HIST("MixingQA/hMECollisionBins"), colBinning.getBin({collision1.posZ(), multiplicityCol}));

      auto groupPartsOne = partsOneMC->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision1.globalIndex(), cache);
      auto groupPartsTwo = partsTwoMC->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision2.globalIndex(), cache);
      auto groupPartsThree = partsThreeMC->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision3.globalIndex(), cache);

      const auto& magFieldTesla1 = collision1.magField();
      const auto& magFieldTesla2 = collision2.magField();
      const auto& magFieldTesla3 = collision3.magField();

      if ((magFieldTesla1 != magFieldTesla2) || (magFieldTesla2 != magFieldTesla3) || (magFieldTesla1 != magFieldTesla3)) {
        continue;
      }
      // CONSIDER testing different strategies to which events to use

      doMixedEvent<true>(groupPartsOne, groupPartsTwo, groupPartsThree, parts, magFieldTesla1, multiplicityCol);
    }
  }
  PROCESS_SWITCH(femtoDreamPairTaskTrackTrackTrack, processMixedEventMC, "Enable processing mixed events MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<femtoDreamPairTaskTrackTrackTrack>(cfgc),
  };
  return workflow;
}
