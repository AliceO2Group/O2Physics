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

/// \file femtoDreamPairTaskTrackTrack.cxx
/// \brief Tasks that reads the track tables used for the pairing and builds pairs of two tracks
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de
/// \author Georgios Mantzaridis, TU München, georgios.mantzaridis@tum.de
/// \author Anton Riedel, TU München, anton.riedel@tum.de

#include <Framework/Configurable.h>
#include <Framework/Expressions.h>
#include <cstdint>
#include <vector>
#include <bitset>
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

using namespace o2::aod;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::analysis::femtoDream;

struct femtoDreamPairTaskTrackTrack {
  SliceCache cache;
  Preslice<aod::FDParticles> perCol = aod::femtodreamparticle::fdCollisionId;

  using MaskedCollisions = soa::Join<FDCollisions, FDColMasks>;
  using MaskedCollision = MaskedCollisions::iterator;
  femtodreamcollision::BitMaskType MaskBit = -1;

  /// Table for both particles
  Configurable<bool> ConfIsMC{"ConfIsMC", false, "Enable additional Histogramms in the case of a MonteCarlo Run"};
  Configurable<bool> ConfUse3D{"ConfUse3D", false, "Enable three dimensional histogramms (to be used only for analysis with high statistics): k* vs mT vs multiplicity"};
  Configurable<bool> ConfExtendedPlots{"ConfExtendedPlots", false, "Enable additional three dimensional histogramms. High memory consumption. Use for debugging"};
  Configurable<float> ConfHighkstarCut{"ConfHighkstarCut", -1., "Set a cut for high k*, above which the pairs are rejected. Set it to -1 to deactivate it"};

  /// Particle selection part

  /// Particle 1
  Configurable<int> ConfTrk1_PDGCode{"ConfTrk1_PDGCode", 2212, "PDG code of particle 1 (Track)"};
  Configurable<femtodreamparticle::cutContainerType> ConfTrk1_CutBit{"ConfTrk1_CutBit", 3191978, "Selection bit from cutCulator for particle 1 (Track)"};
  Configurable<femtodreamparticle::cutContainerType> ConfTrk1_TPCBit{"ConfTrk1_TPCBit", 4, "PID TPC bit from cutCulator for particle 1 (Track)"};
  Configurable<femtodreamparticle::cutContainerType> ConfTrk1_TPCTOFBit{"ConfTrk1_TPCTOFBit", 2, "PID TPCTOF bit from cutCulator for particle 1 (Track)"};
  Configurable<float> ConfTrk1_PIDThres{"ConfTrk1_PIDThres", 0.75, "Momentum threshold for PID selection for particle 1 (Track)"};
  Configurable<float> ConfTrk1_minPt{"ConfTrk1_minPt", 0., "Minimum pT of partricle 1 (Track)"};
  Configurable<float> ConfTrk1_maxPt{"ConfTrk1_maxPt", 999., "Maximum pT of partricle 1 (Track)"};
  Configurable<float> ConfTrk1_minEta{"ConfTrk1_minEta", -10., "Minimum eta of partricle 1 (Track)"};
  Configurable<float> ConfTrk1_maxEta{"ConfTrk1_maxEta", 10., "Maximum eta of partricle 1 (Track)"};

  /// Partition for particle 1
  Partition<aod::FDParticles> PartitionTrk1 = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack)) &&
                                              (ncheckbit(aod::femtodreamparticle::cut, ConfTrk1_CutBit)) &&
                                              ifnode(aod::femtodreamparticle::pt * (nexp(aod::femtodreamparticle::eta) + nexp(-1.f * aod::femtodreamparticle::eta)) / 2.f <= ConfTrk1_PIDThres, ncheckbit(aod::femtodreamparticle::pidcut, ConfTrk1_TPCBit), ncheckbit(aod::femtodreamparticle::pidcut, ConfTrk1_TPCTOFBit)) &&
                                              (aod::femtodreamparticle::pt > ConfTrk1_minPt) &&
                                              (aod::femtodreamparticle::pt < ConfTrk1_maxPt) &&
                                              (aod::femtodreamparticle::eta > ConfTrk1_minEta) &&
                                              (aod::femtodreamparticle::eta < ConfTrk1_maxEta);

  Partition<soa::Join<aod::FDParticles, aod::FDMCLabels>> PartitionMCTrk1 = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack)) &&
                                                                            (ncheckbit(aod::femtodreamparticle::cut, ConfTrk1_CutBit)) &&
                                                                            ifnode(aod::femtodreamparticle::pt * (nexp(aod::femtodreamparticle::eta) + nexp(-1.f * aod::femtodreamparticle::eta)) / 2.f <= ConfTrk1_PIDThres, ncheckbit(aod::femtodreamparticle::pidcut, ConfTrk1_TPCBit), ncheckbit(aod::femtodreamparticle::pidcut, ConfTrk1_TPCTOFBit)) &&
                                                                            (aod::femtodreamparticle::pt > ConfTrk1_minPt) &&
                                                                            (aod::femtodreamparticle::pt < ConfTrk1_maxPt) &&
                                                                            (aod::femtodreamparticle::eta > ConfTrk1_minEta) &&
                                                                            (aod::femtodreamparticle::eta < ConfTrk1_maxEta);

  /// Histogramming for particle 1
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kTrack, 1> trackHistoPartOne;

  /// Particle 2
  Configurable<bool> ConfIsSame{"ConfIsSame", false, "Set to true if particle 1 and particle 2 are the same species"};
  Configurable<int> ConfTrk2_PDGCode{"ConfTrk2_PDGCode", 2212, "PDG code of particle 2 (Track)"};
  Configurable<femtodreamparticle::cutContainerType> ConfTrk2_CutBit{"ConfTrk2_CutBit", 3191978, "Selection bit from cutCulator for particle 2 (Track)"};
  Configurable<femtodreamparticle::cutContainerType> ConfTrk2_TPCBit{"ConfTrk2_TPCBit", 4, "PID TPC bit from cutCulator for particle 2 (Track)"};
  Configurable<femtodreamparticle::cutContainerType> ConfTrk2_TPCTOFBit{"ConfTrk2_TPCTOFBit", 2, "PID TPCTOF bit from cutCulator for particle 2 (Track)"};
  Configurable<float> ConfTrk2_PIDThres{"ConfTrk2_PIDThres", 0.75, "Momentum threshold for PID selection for particle 2 (Track)"};
  Configurable<float> ConfTrk2_minPt{"ConfTrk2_minPt", 0., "Minimum pT of partricle 2 (Track)"};
  Configurable<float> ConfTrk2_maxPt{"ConfTrk2_maxPt", 999., "Maximum pT of partricle 2 (Track)"};
  Configurable<float> ConfTrk2_minEta{"ConfTrk2_minEta", -10., "Minimum eta of partricle 2 (Track)"};
  Configurable<float> ConfTrk2_maxEta{"ConfTrk2_maxEta", 10., "Maximum eta of partricle 2 (Track)"};

  /// Partition for particle 2
  Partition<aod::FDParticles> PartitionTrk2 = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack)) &&
                                              (ncheckbit(aod::femtodreamparticle::cut, ConfTrk2_CutBit)) &&
                                              ifnode(aod::femtodreamparticle::pt * (nexp(aod::femtodreamparticle::eta) + nexp(-1.f * aod::femtodreamparticle::eta)) / 2.f <= ConfTrk2_PIDThres, ncheckbit(aod::femtodreamparticle::pidcut, ConfTrk2_TPCBit), ncheckbit(aod::femtodreamparticle::pidcut, ConfTrk2_TPCTOFBit)) &&
                                              (aod::femtodreamparticle::pt > ConfTrk2_minPt) &&
                                              (aod::femtodreamparticle::pt < ConfTrk2_maxPt) &&
                                              (aod::femtodreamparticle::eta > ConfTrk2_minEta) &&
                                              (aod::femtodreamparticle::eta < ConfTrk2_maxEta);

  Partition<soa::Join<aod::FDParticles, aod::FDMCLabels>> PartitionMCTrk2 = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack)) &&
                                                                            (ncheckbit(aod::femtodreamparticle::cut, ConfTrk2_CutBit)) &&
                                                                            ifnode(aod::femtodreamparticle::pt * (nexp(aod::femtodreamparticle::eta) + nexp(-1.f * aod::femtodreamparticle::eta)) / 2.f <= ConfTrk2_PIDThres, ncheckbit(aod::femtodreamparticle::pidcut, ConfTrk2_TPCBit), ncheckbit(aod::femtodreamparticle::pidcut, ConfTrk2_TPCTOFBit)) &&
                                                                            (aod::femtodreamparticle::pt > ConfTrk2_minPt) &&
                                                                            (aod::femtodreamparticle::pt < ConfTrk2_maxPt) &&
                                                                            (aod::femtodreamparticle::eta > ConfTrk2_minEta) &&
                                                                            (aod::femtodreamparticle::eta < ConfTrk2_maxEta);

  /// Histogramming for particle 2
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kTrack, 2> trackHistoPartTwo;

  /// Histogramming for Event
  FemtoDreamEventHisto eventHisto;

  /// particle part
  ConfigurableAxis ConfTempFitVarBins{"ConfTempFitVarBins", {300, -0.15, 0.15}, "binning of the TempFitVar in the pT vs. TempFitVar plot"};
  ConfigurableAxis ConfTempFitVarpTBins{"ConfTempFitVarpTBins", {20, 0.5, 4.05}, "pT binning of the pT vs. TempFitVar plot"};

  /// Correlation part
  ConfigurableAxis ConfMultBins{"ConfMultBins", {VARIABLE_WIDTH, 0.0f, 4.0f, 8.0f, 12.0f, 16.0f, 20.0f, 24.0f, 28.0f, 32.0f, 36.0f, 40.0f, 44.0f, 48.0f, 52.0f, 56.0f, 60.0f, 64.0f, 68.0f, 72.0f, 76.0f, 80.0f, 84.0f, 88.0f, 92.0f, 96.0f, 100.0f, 200.0f, 99999.f}, "Mixing bins - multiplicity"}; // \todo to be obtained from the hash task
  ConfigurableAxis ConfVtxBins{"ConfVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};

  ConfigurableAxis ConfmTBins3D{"ConfmTBins3D", {VARIABLE_WIDTH, 1.02f, 1.14f, 1.20f, 1.26f, 1.38f, 1.56f, 1.86f, 4.50f}, "mT Binning for the 3Dimensional plot: k* vs multiplicity vs mT (set <<ConfUse3D>> to true in order to use)"};
  ConfigurableAxis ConfmultBins3D{"ConfmultBins3D", {VARIABLE_WIDTH, 0.0f, 20.0f, 30.0f, 40.0f, 99999.0f}, "multiplicity Binning for the 3Dimensional plot: k* vs multiplicity vs mT (set <<ConfUse3D>> to true in order to use)"};
  ColumnBinningPolicy<aod::collision::PosZ, aod::femtodreamcollision::MultNtr> colBinning{{ConfVtxBins, ConfMultBins}, true};

  ConfigurableAxis ConfkstarBins{"ConfkstarBins", {1500, 0., 6.}, "binning kstar"};
  ConfigurableAxis ConfkTBins{"ConfkTBins", {150, 0., 9.}, "binning kT"};
  ConfigurableAxis ConfmTBins{"ConfmTBins", {225, 0., 7.5}, "binning mT"};
  Configurable<int> ConfNEventsMix{"ConfNEventsMix", 5, "Number of events for mixing"};
  Configurable<bool> ConfIsCPR{"ConfIsCPR", true, "Close Pair Rejection"};
  Configurable<bool> ConfCPRPlotPerRadii{"ConfCPRPlotPerRadii", false, "Plot CPR per radii"};
  Configurable<float> ConfCPRdeltaPhiMax{"ConfCPRdeltaPhiMax", 0.01, "Max. Delta Phi for Close Pair Rejection"};
  Configurable<float> ConfCPRdeltaEtaMax{"ConfCPRdeltaEtaMax", 0.01, "Max. Delta Eta for Close Pair Rejection"};
  ConfigurableAxis ConfDummy{"ConfDummy", {1, 0, 1}, "Dummy axis"};

  FemtoDreamContainer<femtoDreamContainer::EventType::same, femtoDreamContainer::Observable::kstar> sameEventCont;
  FemtoDreamContainer<femtoDreamContainer::EventType::mixed, femtoDreamContainer::Observable::kstar> mixedEventCont;
  FemtoDreamPairCleaner<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kTrack> pairCleaner;
  FemtoDreamDetaDphiStar<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kTrack> pairCloseRejection;
  /// Histogram output
  HistogramRegistry qaRegistry{"TrackQA", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry resultRegistry{"Correlations", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry MixQaRegistry{"MixQaRegistry", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext& context)
  {
    eventHisto.init(&qaRegistry);
    trackHistoPartOne.init(&qaRegistry, ConfTempFitVarpTBins, ConfTempFitVarBins, ConfDummy, ConfDummy, ConfDummy, ConfDummy, ConfIsMC, ConfTrk1_PDGCode);
    if (!ConfIsSame) {
      trackHistoPartTwo.init(&qaRegistry, ConfTempFitVarpTBins, ConfTempFitVarBins, ConfDummy, ConfDummy, ConfDummy, ConfDummy, ConfIsMC, ConfTrk2_PDGCode);
    }

    MixQaRegistry.add("MixingQA/hSECollisionBins", ";bin;Entries", kTH1F, {{120, -0.5, 119.5}});
    MixQaRegistry.add("MixingQA/hMECollisionBins", ";bin;Entries", kTH1F, {{120, -0.5, 119.5}});

    sameEventCont.init(&resultRegistry, ConfkstarBins, ConfMultBins, ConfkTBins, ConfmTBins, ConfmultBins3D, ConfmTBins3D, ConfIsMC, ConfUse3D, ConfExtendedPlots, ConfHighkstarCut);
    mixedEventCont.init(&resultRegistry, ConfkstarBins, ConfMultBins, ConfkTBins, ConfmTBins, ConfmultBins3D, ConfmTBins3D, ConfIsMC, ConfUse3D, ConfExtendedPlots, ConfHighkstarCut);
    sameEventCont.setPDGCodes(ConfTrk1_PDGCode, ConfTrk2_PDGCode);
    mixedEventCont.setPDGCodes(ConfTrk1_PDGCode, ConfTrk2_PDGCode);
    pairCleaner.init(&qaRegistry);
    if (ConfIsCPR.value) {
      pairCloseRejection.init(&resultRegistry, &qaRegistry, ConfCPRdeltaPhiMax.value, ConfCPRdeltaEtaMax.value, ConfCPRPlotPerRadii.value);
    }

    // get bit for the collision mask
    std::bitset<8 * sizeof(femtodreamcollision::BitMaskType)> mask;
    int index = 0;
    auto& workflows = context.services().get<RunningWorkflowInfo const>();
    for (DeviceSpec const& device : workflows.devices) {
      if (device.name.find("femto-dream-pair-task-track-track") != std::string::npos) {
        if (containsNameValuePair(device.options, "ConfTrk1_CutBit", ConfTrk1_CutBit.value) &&
            containsNameValuePair(device.options, "ConfTrk1_TPCBit", ConfTrk1_TPCBit.value) &&
            containsNameValuePair(device.options, "ConfTrk1_TPCTOFBit", ConfTrk1_TPCTOFBit.value) &&
            containsNameValuePair(device.options, "ConfTrk1_PIDThres", ConfTrk1_PIDThres.value) &&
            containsNameValuePair(device.options, "ConfTrk1_minPt", ConfTrk1_minPt.value) &&
            containsNameValuePair(device.options, "ConfTrk1_maxPt", ConfTrk1_maxPt.value) &&
            containsNameValuePair(device.options, "ConfTrk1_minEta", ConfTrk1_minEta.value) &&
            containsNameValuePair(device.options, "ConfTrk1_maxEta", ConfTrk1_maxEta.value) &&
            containsNameValuePair(device.options, "ConfTrk2_CutBit", ConfTrk2_CutBit.value) &&
            containsNameValuePair(device.options, "ConfTrk2_TPCBit", ConfTrk2_TPCBit.value) &&
            containsNameValuePair(device.options, "ConfTrk2_TPCTOFBit", ConfTrk2_TPCTOFBit.value) &&
            containsNameValuePair(device.options, "ConfTrk2_PIDThres", ConfTrk2_PIDThres.value) &&
            containsNameValuePair(device.options, "ConfTrk2_minPt", ConfTrk2_minPt.value) &&
            containsNameValuePair(device.options, "ConfTrk2_maxPt", ConfTrk2_maxPt.value) &&
            containsNameValuePair(device.options, "ConfTrk2_minEta", ConfTrk2_minEta.value) &&
            containsNameValuePair(device.options, "ConfTrk2_maxEta", ConfTrk2_maxEta.value)) {
          mask.set(index);
          MaskBit = static_cast<femtodreamcollision::BitMaskType>(mask.to_ulong());
          LOG(info) << "Device name matched: " << device.name;
          LOG(info) << "Bitmask for collisions: " << mask.to_string();
          break;
        } else {
          index++;
        }
      }
    }
    if ((doprocessSameEvent && doprocessSameEventMasked) ||
        (doprocessMixedEvent && doprocessMixedEventMasked) ||
        (doprocessSameEventMC && doprocessSameEventMCMasked) ||
        (doprocessMixedEventMC && doprocessMixedEventMCMasked)) {
      LOG(fatal) << "Normal and masked processing cannot be activated simultaneously!";
    }
  };

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
  void doSameEvent(PartitionType SliceTrk1, PartitionType SliceTrk2, PartType parts, float magFieldTesla, int multCol)
  {
    /// Histogramming same event
    for (auto& part : SliceTrk1) {
      trackHistoPartOne.fillQA<isMC, false>(part, aod::femtodreamparticle::kPt);
    }

    if (!ConfIsSame.value) {
      for (auto& part : SliceTrk2) {
        trackHistoPartTwo.fillQA<isMC, false>(part, aod::femtodreamparticle::kPt);
      }
    }

    /// Now build the combinations
    if (ConfIsSame.value) {
      for (auto& [p1, p2] : combinations(CombinationsStrictlyUpperIndexPolicy(SliceTrk1, SliceTrk2))) {
        if (ConfIsCPR.value) {
          if (pairCloseRejection.isClosePair(p1, p2, parts, magFieldTesla)) {
            continue;
          }
        }
        // track cleaning
        if (!pairCleaner.isCleanPair(p1, p2, parts)) {
          continue;
        }
        sameEventCont.setPair<isMC>(p1, p2, multCol, ConfUse3D, ConfExtendedPlots);
      }
    } else {
      for (auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(SliceTrk1, SliceTrk2))) {
        if (ConfIsCPR.value) {
          if (pairCloseRejection.isClosePair(p1, p2, parts, magFieldTesla)) {
            continue;
          }
        }
        // track cleaning
        if (!pairCleaner.isCleanPair(p1, p2, parts)) {
          continue;
        }
        sameEventCont.setPair<isMC>(p1, p2, multCol, ConfUse3D, ConfExtendedPlots);
      }
    }
  }

  /// process function for to call doSameEvent with Data
  /// \param col subscribe to the collision table (Data)
  /// \param parts subscribe to the femtoDreamParticleTable
  void processSameEvent(o2::aod::FDCollision& col, o2::aod::FDParticles& parts)
  {
    fillCollision(col);
    auto SliceTrk1 = PartitionTrk1->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    auto SliceTrk2 = PartitionTrk2->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    if (SliceTrk1.size() == 0 && SliceTrk2.size() == 0) {
      return;
    }
    doSameEvent<false>(SliceTrk1, SliceTrk2, parts, col.magField(), col.multNtr());
  }
  PROCESS_SWITCH(femtoDreamPairTaskTrackTrack, processSameEvent, "Enable processing same event", true);

  void processSameEventMasked(MaskedCollision& col, o2::aod::FDParticles& parts)
  {
    if (ConfIsSame.value) {
      if ((col.bitmaskTrackOne() & MaskBit) != MaskBit) {
        return;
      }
    } else {
      if ((col.bitmaskTrackOne() & MaskBit) != MaskBit && (col.bitmaskTrackTwo() & MaskBit) != MaskBit) {
        return;
      }
    }
    fillCollision(col);
    auto SliceTrk1 = PartitionTrk1->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    auto SliceTrk2 = PartitionTrk2->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    doSameEvent<false>(SliceTrk1, SliceTrk2, parts, col.magField(), col.multNtr());
  }
  PROCESS_SWITCH(femtoDreamPairTaskTrackTrack, processSameEventMasked, "Enable processing same event with masks", false);

  /// process function for to call doSameEvent with Monte Carlo
  /// \param col subscribe to the collision table (Monte Carlo Reconstructed reconstructed)
  /// \param parts subscribe to joined table FemtoDreamParticles and FemtoDreamMCLables to access Monte Carlo truth
  /// \param FemtoDreamMCParticles subscribe to the Monte Carlo truth table
  void processSameEventMC(o2::aod::FDCollision& col, soa::Join<o2::aod::FDParticles, o2::aod::FDMCLabels>& parts,
                          o2::aod::FDMCParticles&)
  {
    fillCollision(col);
    auto SliceMCTrk1 = PartitionMCTrk1->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    auto SliceMCTrk2 = PartitionMCTrk2->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    if (SliceMCTrk1.size() == 0 && SliceMCTrk2.size() == 0) {
      return;
    }
    doSameEvent<true>(SliceMCTrk1, SliceMCTrk2, parts, col.magField(), col.multNtr());
  }
  PROCESS_SWITCH(femtoDreamPairTaskTrackTrack, processSameEventMC, "Enable processing same event for Monte Carlo", false);

  void processSameEventMCMasked(MaskedCollision& col, soa::Join<o2::aod::FDParticles, o2::aod::FDMCLabels>& parts,
                                o2::aod::FDMCParticles&)
  {
    if (ConfIsSame.value) {
      if ((col.bitmaskTrackOne() & MaskBit) != MaskBit) {
        return;
      }
    } else {
      if ((col.bitmaskTrackOne() & MaskBit) != MaskBit && (col.bitmaskTrackTwo() & MaskBit) != MaskBit) {
        return;
      }
    }
    fillCollision(col);
    auto SliceMCTrk1 = PartitionMCTrk1->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    auto SliceMCTrk2 = PartitionMCTrk2->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    doSameEvent<true>(SliceMCTrk1, SliceMCTrk2, parts, col.magField(), col.multNtr());
  }
  PROCESS_SWITCH(femtoDreamPairTaskTrackTrack, processSameEventMCMasked, "Enable processing same event for Monte Carlo with masked collisions", false);

  /// This function processes the mixed event
  /// \tparam PartitionType
  /// \tparam PartType
  /// \tparam isMC: enables Monte Carlo truth specific histograms
  /// \param groupPartsOne partition for the first particle passed by the process function
  /// \param groupPartsTwo partition for the second particle passed by the process function
  /// \param parts femtoDreamParticles table (in case of Monte Carlo joined with FemtoDreamMCLabels)
  /// \param magFieldTesla magnetic field of the collision
  /// \param multCol multiplicity of the collision
  template <bool isMC, typename PartitionType, typename PartType>
  void doMixedEvent(PartitionType SliceTrk1, PartitionType SliceTrk2, PartType parts, float magFieldTesla, int multCol)
  {
    for (auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(SliceTrk1, SliceTrk2))) {
      if (ConfIsCPR.value) {
        if (pairCloseRejection.isClosePair(p1, p2, parts, magFieldTesla)) {
          continue;
        }
      }
      mixedEventCont.setPair<isMC>(p1, p2, multCol, ConfUse3D, ConfExtendedPlots);
    }
  }

  /// process function for to call doMixedEvent with Data
  /// @param cols subscribe to the collisions table (Data)
  /// @param parts subscribe to the femtoDreamParticleTable
  void processMixedEvent(o2::aod::FDCollisions& cols, o2::aod::FDParticles& parts)
  {
    for (auto const& [collision1, collision2] : soa::selfCombinations(colBinning, ConfNEventsMix, -1, cols, cols)) {
      const int multiplicityCol = collision1.multNtr();
      MixQaRegistry.fill(HIST("MixingQA/hMECollisionBins"), colBinning.getBin({collision1.posZ(), multiplicityCol}));
      auto SliceTrk1 = PartitionTrk1->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision1.globalIndex(), cache);
      auto SliceTrk2 = PartitionTrk2->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision2.globalIndex(), cache);
      if (SliceTrk1.size() == 0 || SliceTrk2.size() == 0) {
        continue;
      }
      const auto magFieldTesla1 = collision1.magField();
      doMixedEvent<false>(SliceTrk1, SliceTrk2, parts, magFieldTesla1, multiplicityCol);
    }
  }
  PROCESS_SWITCH(femtoDreamPairTaskTrackTrack, processMixedEvent, "Enable processing mixed events", true);

  void processMixedEventMasked(MaskedCollisions& cols, o2::aod::FDParticles& parts)
  {
    Partition<MaskedCollisions> PartitionMaskedCol1 = (aod::femtodreamcollision::bitmaskTrackOne & MaskBit) == MaskBit;
    Partition<MaskedCollisions> PartitionMaskedCol2 = (aod::femtodreamcollision::bitmaskTrackTwo & MaskBit) == MaskBit;

    PartitionMaskedCol1.bindTable(cols);
    PartitionMaskedCol2.bindTable(cols);

    for (auto const& [collision1, collision2] : combinations(soa::CombinationsBlockUpperIndexPolicy(colBinning, ConfNEventsMix.value, -1, PartitionMaskedCol1, PartitionMaskedCol2))) {
      const auto multiplicityCol = collision1.multNtr();
      const auto magFieldTesla1 = collision1.magField();
      MixQaRegistry.fill(HIST("MixingQA/hMECollisionBins"), colBinning.getBin({collision1.posZ(), multiplicityCol}));
      auto SliceTrk1 = PartitionTrk1->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision1.globalIndex(), cache);
      auto SliceTrk2 = PartitionTrk2->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision2.globalIndex(), cache);
      doMixedEvent<false>(SliceTrk1, SliceTrk2, parts, magFieldTesla1, multiplicityCol);
    }
  }
  PROCESS_SWITCH(femtoDreamPairTaskTrackTrack, processMixedEventMasked, "Enable processing mixed events", false);

  /// brief process function for to call doMixedEvent with Monte Carlo
  /// @param cols subscribe to the collisions table (Monte Carlo Reconstructed reconstructed)
  /// @param parts subscribe to joined table FemtoDreamParticles and FemtoDreamMCLables to access Monte Carlo truth
  /// @param FemtoDreamMCParticles subscribe to the Monte Carlo truth table
  void processMixedEventMC(o2::aod::FDCollisions& cols, soa::Join<o2::aod::FDParticles, o2::aod::FDMCLabels>& parts, o2::aod::FDMCParticles&)
  {
    for (auto& [collision1, collision2] : soa::selfCombinations(colBinning, ConfNEventsMix, -1, cols, cols)) {
      const int multiplicityCol = collision1.multNtr();
      MixQaRegistry.fill(HIST("MixingQA/hMECollisionBins"), colBinning.getBin({collision1.posZ(), multiplicityCol}));
      auto SliceMCTrk1 = PartitionMCTrk1->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision1.globalIndex(), cache);
      auto SliceMCTrk2 = PartitionMCTrk2->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision2.globalIndex(), cache);
      if (SliceMCTrk1.size() == 0 || SliceMCTrk2.size() == 0) {
        continue;
      }
      const auto magFieldTesla1 = collision1.magField();
      doMixedEvent<true>(SliceMCTrk1, SliceMCTrk2, parts, magFieldTesla1, multiplicityCol);
    }
  }
  PROCESS_SWITCH(femtoDreamPairTaskTrackTrack, processMixedEventMC, "Enable processing mixed events MC", false);

  void processMixedEventMCMasked(MaskedCollisions& cols, soa::Join<o2::aod::FDParticles, o2::aod::FDMCLabels>& parts, o2::aod::FDMCParticles&)
  {
    Partition<MaskedCollisions> PartitionMaskedCol1 = (aod::femtodreamcollision::bitmaskTrackOne & MaskBit) == MaskBit;
    Partition<MaskedCollisions> PartitionMaskedCol2 = (aod::femtodreamcollision::bitmaskTrackTwo & MaskBit) == MaskBit;

    PartitionMaskedCol1.bindTable(cols);
    PartitionMaskedCol2.bindTable(cols);

    for (auto const& [collision1, collision2] : combinations(soa::CombinationsBlockUpperIndexPolicy(colBinning, ConfNEventsMix.value, -1, PartitionMaskedCol1, PartitionMaskedCol2))) {
      const auto multiplicityCol = collision1.multNtr();
      const auto& magFieldTesla1 = collision1.magField();
      MixQaRegistry.fill(HIST("MixingQA/hMECollisionBins"), colBinning.getBin({collision1.posZ(), multiplicityCol}));
      auto SliceMCTrk1 = PartitionMCTrk1->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision1.globalIndex(), cache);
      auto SliceMCTrk2 = PartitionMCTrk2->sliceByCached(aod::femtodreamparticle::fdCollisionId, collision2.globalIndex(), cache);
      doMixedEvent<true>(SliceMCTrk1, SliceMCTrk2, parts, magFieldTesla1, multiplicityCol);
    }
  }
  PROCESS_SWITCH(femtoDreamPairTaskTrackTrack, processMixedEventMCMasked, "Enable processing mixed events MC with masked collisions", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<femtoDreamPairTaskTrackTrack>(cfgc),
  };
  return workflow;
}
