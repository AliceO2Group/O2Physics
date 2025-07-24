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

/// \file femtoWorldPairTaskPionPion.cxx
/// \brief Tasks that reads the track tables used for the pairing and builds pairs of two tracks - pions
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de
/// \author Zuzanna Chochulska, WUT Warsaw, zchochul@cern.ch
/// \author Alicja Plachta, WUT Warsaw

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/StepTHn.h"

#include "PWGCF/FemtoWorld/DataModel/FemtoWorldDerived.h"
#include "PWGCF/FemtoWorld/Core/FemtoWorldParticleHisto.h"
#include "PWGCF/FemtoWorld/Core/FemtoWorldEventHisto.h"
#include "PWGCF/FemtoWorld/Core/FemtoWorldPairCleaner.h"
#include "PWGCF/FemtoWorld/Core/FemtoWorldPionContainer.h"
#include "PWGCF/FemtoWorld/Core/FemtoWorldDetaDphiStar.h"
#include "PWGCF/FemtoWorld/Core/FemtoWorldUtils.h"

using namespace o2;
using namespace o2::analysis::femtoWorld;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

namespace o2::aod
{
using FemtoWorldParticlesMerged = aod::FemtoWorldParticles;
using FemtoWorldParticleMerged = FemtoWorldParticlesMerged::iterator;
} // namespace o2::aod

struct femtoWorldPairTaskPionPion {

  SliceCache cache;
  Preslice<aod::FemtoWorldParticlesMerged> perCol = aod::femtoworldparticle::femtoWorldCollisionId;

  /// Particle selection part

  Configurable<float> ConfNsigmaTPCPion{"ConfNsigmaTPCPion", 3.0, "TPC Pion Sigma for momentum < 0.5"};
  Configurable<float> ConfNsigmaCombinedPion{"ConfNsigmaCombinedPion", 3.0, "TPC and TOF Pion Sigma (combined) for momentum > 0.5"};

  // Configurables for cuts
  // First particle
  Configurable<float> cfgPtLowPart1{"cfgPtLowPart1", 0.14, "Lower limit for Pt for the first particle"};
  Configurable<float> cfgPtHighPart1{"cfgPtHighPart1", 1.5, "Higher limit for Pt for the first particle"};
  Configurable<float> cfgEtaLowPart1{"cfgEtaLowPart1", -0.8, "Lower limit for Eta for the first particle"};
  Configurable<float> cfgEtaHighPart1{"cfgEtaHighPart1", 0.8, "Higher limit for Eta for the first particle"};
  Configurable<float> cfgDcaXYPart1{"cfgDcaXYPart1", 2.4, "Value for DCA_XY for the first particle"};
  Configurable<float> cfgDcaZPart1{"cfgDcaZPart1", 3.2, "Value for DCA_Z for the first particle"};
  Configurable<int> cfgTpcClPart1{"cfgTpcClPart1", 88, "Number of tpc clasters for the first particle"};             // min number of found TPC clusters
  Configurable<int> cfgTpcCrosRoPart1{"cfgTpcCrosRoPart1", 70, "Number of tpc crossed rows for the first particle"}; // min number of crossed rows
  Configurable<float> cfgChi2TpcPart1{"cfgChi2TpcPart1", 4.0, "Chi2 / cluster for the TPC track segment for the first particle"};
  Configurable<float> cfgChi2ItsPart1{"cfgChi2ItsPart1", 36.0, "Chi2 / cluster for the ITS track segment for the first particle"};
  // Second particle
  Configurable<float> cfgPtLowPart2{"cfgPtLowPart2", 0.14, "Lower limit for Pt for the second particle"};
  Configurable<float> cfgPtHighPart2{"cfgPtHighPart2", 1.5, "Higher limit for Pt for the second particle"};
  Configurable<float> cfgEtaLowPart2{"cfgEtaLowPart2", -0.8, "Lower limit for Eta for the second particle"};
  Configurable<float> cfgEtaHighPart2{"cfgEtaHighPart2", 0.8, "Higher limit for Eta for the second particle"};
  Configurable<float> cfgDcaXYPart2{"cfgDcaXYPart2", 2.4, "Value for DCA_XY for the second particle"};
  Configurable<float> cfgDcaZPart2{"cfgDcaZPart2", 3.2, "Value for DCA_Z for the second particle"};
  Configurable<int> cfgTpcClPart2{"cfgTpcClPart2", 88, "Number of tpc clasters for the second particle"}; // min number of found TPC clusters

  Configurable<int> cfgTpcCrosRoPart2{"cfgTpcCrosRoPart2", 70, "Number of tpc crossed rows for the second particle"}; // min number of crossed rows
  Configurable<float> cfgChi2TpcPart2{"cfgChi2TpcPart2", 4.0, "Chi2 / cluster for the TPC track segment for the second particle"};
  Configurable<float> cfgChi2ItsPart2{"cfgChi2ItsPart2", 36.0, "Chi2 / cluster for the ITS track segment for the second particle"};

  /// Particle 1 Pion+
  Configurable<int> ConfPDGCodePartOne{"ConfPDGCodePartOne", 211, "Particle 1 - PDG code"};
  Configurable<std::vector<int>> ConfPIDPartOne{"ConfPIDPartOne", std::vector<int>{2}, "Particle 1 - Read from json file"}; // we also need the possibility to specify whether the bit is true/false ->std>>vector<std::pair<int, int>>int>>

  /// Partition for particle 1
  Partition<aod::FemtoWorldParticlesMerged> partsOne = (aod::femtoworldparticle::partType == uint8_t(aod::femtoworldparticle::ParticleType::kTrack))                         // particle type cut
                                                       && (aod::femtoworldparticle::pt < cfgPtHighPart1) && (aod::femtoworldparticle::pt > cfgPtLowPart1)                    // simple pT cuts
                                                       && (aod::femtoworldparticle::eta < cfgEtaHighPart1) && (aod::femtoworldparticle::eta > cfgEtaLowPart1)                // Eta cuts
                                                       && (o2::aod::track::dcaXY < cfgDcaXYPart1) && (o2::aod::track::dcaZ < cfgDcaZPart1)                                   // DCA cuts for XY and Z
                                                       && (aod::femtoworldparticle::tpcNClsFound > (uint8_t)cfgTpcClPart1)                                                   // Number of found TPC clusters
                                                       && (aod::femtoworldparticle::tpcNClsCrossedRows > (uint8_t)cfgTpcCrosRoPart1)                                         // Crossed rows TPC
                                                       && (aod::femtoworldparticle::itsChi2NCl < cfgChi2ItsPart1) && (aod::femtoworldparticle::tpcChi2NCl < cfgChi2TpcPart1) //&& // chi2 cuts
                                                       && (aod::femtoworldparticle::sign > int8_t(0))                                                                        // Sign (Pion+)
    // todo: globaltrack cut
    ;

  Partition<aod::FemtoWorldParticlesMerged> partsOneFailed =                                                              //~(aod::femtoworldparticle::partType == uint8_t(aod::femtoworldparticle::ParticleType::kTrack)) ||                   // particle type cut
    ((aod::femtoworldparticle::pt > cfgPtHighPart1) || (aod::femtoworldparticle::pt < cfgPtLowPart1)) ||                  // pT cuts
    ((aod::femtoworldparticle::eta > cfgEtaHighPart1) || (aod::femtoworldparticle::eta < cfgEtaLowPart1)) ||              // Eta cuts
    (o2::aod::track::dcaXY > cfgDcaXYPart1) || (o2::aod::track::dcaZ > cfgDcaZPart1) ||                                   // DCA cuts for XY and Z
    (aod::femtoworldparticle::tpcNClsFound < (uint8_t)cfgTpcClPart1) ||                                                   // Number of found TPC clusters
    (aod::femtoworldparticle::tpcNClsCrossedRows < (uint8_t)cfgTpcCrosRoPart1) ||                                         // Crossed rows TPC
    (aod::femtoworldparticle::itsChi2NCl > cfgChi2ItsPart1) || (aod::femtoworldparticle::tpcChi2NCl > cfgChi2TpcPart1) || //&& // chi2 cuts
    (aod::femtoworldparticle::sign < int8_t(0));
  /// Histogramming for particle 1
  FemtoWorldParticleHisto<aod::femtoworldparticle::ParticleType::kTrack, 1> trackHistoPartOne;
  FemtoWorldParticleHisto<aod::femtoworldparticle::ParticleType::kTrack, 2> trackHistoPartOneFailed;

  /// Particle 2 Pion-
  Configurable<bool> ConfIsSame{"ConfIsSame", false, "Pairs of the same particle"};
  Configurable<int> ConfPDGCodePartTwo{"ConfPDGCodePartTwo", 211, "Particle 2 - PDG code"};
  Configurable<std::vector<int>> ConfPIDPartTwo{"ConfPIDPartTwo", std::vector<int>{2}, "Particle 2 - Read from json file"}; // we also need the possibility to specify whether the bit is true/false ->std>>vector<std::pair<int, int>>

  /// Partition for particle 2
  Partition<aod::FemtoWorldParticlesMerged> partsTwo = (aod::femtoworldparticle::partType == uint8_t(aod::femtoworldparticle::ParticleType::kTrack))                         // particle type cut
                                                       && (aod::femtoworldparticle::pt < cfgPtHighPart2) && (aod::femtoworldparticle::pt > cfgPtLowPart2)                    // pT cuts
                                                       && (aod::femtoworldparticle::eta < cfgEtaHighPart2) && (aod::femtoworldparticle::eta > cfgEtaLowPart2)                // Eta cuts
                                                       && (o2::aod::track::dcaXY < cfgDcaXYPart2) && (o2::aod::track::dcaZ < cfgDcaZPart2)                                   // DCA cuts for XY and Z
                                                       && (aod::femtoworldparticle::tpcNClsFound > (uint8_t)cfgTpcClPart2)                                                   // Number of found TPC clusters
                                                       && (aod::femtoworldparticle::tpcNClsCrossedRows > (uint8_t)cfgTpcCrosRoPart2)                                         // Crossed rows TPC
                                                       && (aod::femtoworldparticle::itsChi2NCl < cfgChi2ItsPart2) && (aod::femtoworldparticle::tpcChi2NCl < cfgChi2TpcPart2) //  // chi2 cuts
                                                       && (aod::femtoworldparticle::sign < int8_t(0))                                                                        // Sign (K-)
                                                                                                                                                                             //  todo: globaltrack cut
    ;
  Partition<aod::FemtoWorldParticlesMerged> partsTwoFailed =                                                           //(aod::femtoworldparticle::partType == uint8_t(aod::femtoworldparticle::ParticleType::kTrack)) ||                   // particle type cut
    ((aod::femtoworldparticle::pt > cfgPtHighPart2) || (aod::femtoworldparticle::pt < cfgPtLowPart2)) ||               // pT cuts
    ((aod::femtoworldparticle::eta > cfgEtaHighPart2) || (aod::femtoworldparticle::eta < cfgEtaLowPart2)) ||           // Eta cuts
    (o2::aod::track::dcaXY > cfgDcaXYPart2) || (o2::aod::track::dcaZ > cfgDcaZPart2) ||                                // DCA cuts for XY and Z
    (aod::femtoworldparticle::tpcNClsFound < (uint8_t)cfgTpcClPart2) ||                                                // Number of found TPC clusters
    (aod::femtoworldparticle::tpcNClsCrossedRows < (uint8_t)cfgTpcCrosRoPart2) ||                                      // Crossed rows TPC
    (aod::femtoworldparticle::itsChi2NCl > cfgChi2ItsPart2) || (aod::femtoworldparticle::tpcChi2NCl > cfgChi2TpcPart2) //&& // chi2 cuts
    ;
  /// Histogramming for particle 2
  FemtoWorldParticleHisto<aod::femtoworldparticle::ParticleType::kTrack, 3> trackHistoPartTwo;
  FemtoWorldParticleHisto<aod::femtoworldparticle::ParticleType::kTrack, 4> trackHistoPartTwoFailed;

  /// Histogramming for Event
  FemtoWorldEventHisto eventHisto;

  /// The configurables need to be passed to an std::vector
  std::vector<int> vPIDPartOne, vPIDPartTwo;

  /// Correlation part
  ConfigurableAxis CfgMultBins{"CfgMultBins", {VARIABLE_WIDTH, 0.0f, 20.0f, 40.0f, 60.0f, 80.0f, 100.0f, 200.0f, 99999.f}, "Mixing bins - multiplicity"};
  ConfigurableAxis CfgVtxBins{"CfgVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};

  ColumnBinningPolicy<aod::collision::PosZ, aod::femtoworldcollision::MultV0M> colBinning{{CfgVtxBins, CfgMultBins}, true};

  ConfigurableAxis CfgkstarBins{"CfgkstarBins", {1500, 0., 6.}, "binning kstar"};
  ConfigurableAxis CfgkTBins{"CfgkTBins", {150, 0., 9.}, "binning kT"};
  ConfigurableAxis CfgmTBins{"CfgmTBins", {225, 0., 7.5}, "binning mT"};
  Configurable<int> ConfNEventsMix{"ConfNEventsMix", 5, "Number of events for mixing"};
  Configurable<bool> ConfIsCPR{"ConfIsCPR", true, "Close Pair Rejection"};
  Configurable<bool> ConfCPRPlotPerRadii{"ConfCPRPlotPerRadii", false, "Plot CPR per radii"};
  Configurable<int> ConfPhiBins{"ConfPhiBins", 15, "Number of phi bins in deta dphi"};
  Configurable<int> ConfEtaBins{"ConfEtaBins", 15, "Number of eta bins in deta dphi"};
  Configurable<int> ConfMInvBins{"ConfMInvBins", 1000, "Number of bins in mInv distribution"};

  Configurable<bool> ConfPlusMinus{"ConfPlusMinus", true, "Process pi+pi- pions pairs"};
  Configurable<bool> ConfPlusPlus{"ConfPlusPlus", true, "Process pi+pi+ pion pairs"};
  Configurable<bool> ConfMinusMinus{"ConfMinusMinus", true, "Process pi-pi- pion pairs"};
  FemtoWorldPionContainer<femtoWorldPionContainer::EventType::samePM, femtoWorldPionContainer::Observable::kstar> sameEventPlusMinCont;
  FemtoWorldPionContainer<femtoWorldPionContainer::EventType::samePP, femtoWorldPionContainer::Observable::kstar> sameEventPlusPlusCont;
  FemtoWorldPionContainer<femtoWorldPionContainer::EventType::sameMM, femtoWorldPionContainer::Observable::kstar> sameEventMinMinCont;
  FemtoWorldPionContainer<femtoWorldPionContainer::EventType::mixedPM, femtoWorldPionContainer::Observable::kstar> mixedEventPlusMinCont;
  FemtoWorldPionContainer<femtoWorldPionContainer::EventType::mixedPP, femtoWorldPionContainer::Observable::kstar> mixedEventPlusPlusCont;
  FemtoWorldPionContainer<femtoWorldPionContainer::EventType::mixedMM, femtoWorldPionContainer::Observable::kstar> mixedEventMinMinCont;

  FemtoWorldPairCleaner<aod::femtoworldparticle::ParticleType::kTrack, aod::femtoworldparticle::ParticleType::kTrack> pairCleaner;
  FemtoWorldDetaDphiStar<aod::femtoworldparticle::ParticleType::kTrack, aod::femtoworldparticle::ParticleType::kTrack> pairCloseRejection;
  /// Histogram output
  HistogramRegistry qaRegistry{"TrackQA", {}, OutputObjHandlingPolicy::AnalysisObject};
  // HistogramRegistry qaRegistryFail{"TrackQAFailed", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry resultRegistryPlusMin{"CorrelationsPlusMin", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry resultRegistryPlusPlus{"CorrelationsPlusPlus", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry resultRegistryMinMin{"CorrelationsMinMin", {}, OutputObjHandlingPolicy::AnalysisObject};

  HistogramRegistry MixQaRegistry{"MixQaRegistry", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
    eventHisto.init(&qaRegistry);
    trackHistoPartOne.init(&qaRegistry);
    trackHistoPartOneFailed.init(&qaRegistry);
    if (!ConfIsSame) {
      trackHistoPartTwo.init(&qaRegistry);
      trackHistoPartTwoFailed.init(&qaRegistry);
    }

    MixQaRegistry.add("MixingQA/hSECollisionBins", ";bin;Entries", kTH1F, {{120, -0.5, 119.5}});
    MixQaRegistry.add("MixingQA/hMECollisionBins", ";bin;Entries", kTH1F, {{120, -0.5, 119.5}});

    if (ConfPlusMinus)
      sameEventPlusMinCont.init(&resultRegistryPlusMin, CfgkstarBins, CfgMultBins, CfgkTBins, CfgmTBins, ConfPhiBins, ConfEtaBins, ConfMInvBins);
    if (ConfPlusMinus)
      sameEventPlusMinCont.setPDGCodes(ConfPDGCodePartOne, ConfPDGCodePartTwo);
    if (ConfPlusPlus)
      sameEventPlusPlusCont.init(&resultRegistryPlusPlus, CfgkstarBins, CfgMultBins, CfgkTBins, CfgmTBins, ConfPhiBins, ConfEtaBins, ConfMInvBins);
    if (ConfPlusPlus)
      sameEventPlusPlusCont.setPDGCodes(ConfPDGCodePartOne, ConfPDGCodePartTwo);
    sameEventMinMinCont.init(&resultRegistryMinMin, CfgkstarBins, CfgMultBins, CfgkTBins, CfgmTBins, ConfPhiBins, ConfEtaBins, ConfMInvBins);
    sameEventMinMinCont.setPDGCodes(ConfPDGCodePartOne, ConfPDGCodePartTwo);

    if (ConfPlusMinus)
      mixedEventPlusMinCont.init(&resultRegistryPlusMin, CfgkstarBins, CfgMultBins, CfgkTBins, CfgmTBins, ConfPhiBins, ConfEtaBins, ConfMInvBins);
    if (ConfPlusMinus)
      mixedEventPlusMinCont.setPDGCodes(ConfPDGCodePartOne, ConfPDGCodePartTwo);
    if (ConfPlusPlus)
      mixedEventPlusPlusCont.init(&resultRegistryPlusPlus, CfgkstarBins, CfgMultBins, CfgkTBins, CfgmTBins, ConfPhiBins, ConfEtaBins, ConfMInvBins);
    if (ConfPlusPlus)
      mixedEventPlusPlusCont.setPDGCodes(ConfPDGCodePartOne, ConfPDGCodePartTwo);
    if (ConfMinusMinus)
      mixedEventMinMinCont.init(&resultRegistryMinMin, CfgkstarBins, CfgMultBins, CfgkTBins, CfgmTBins, ConfPhiBins, ConfEtaBins, ConfMInvBins);
    if (ConfMinusMinus)
      mixedEventMinMinCont.setPDGCodes(ConfPDGCodePartOne, ConfPDGCodePartTwo);

    pairCleaner.init(&qaRegistry);

    if (ConfIsCPR) {
      pairCloseRejection.init(&resultRegistryPlusMin, &qaRegistry, 0.01, 0.01, ConfCPRPlotPerRadii); /// \todo add config for Δη and ΔΦ cut values
      pairCloseRejection.init(&resultRegistryPlusPlus, &qaRegistry, 0.01, 0.01, ConfCPRPlotPerRadii);
      pairCloseRejection.init(&resultRegistryMinMin, &qaRegistry, 0.01, 0.01, ConfCPRPlotPerRadii);
    }

    vPIDPartOne = ConfPIDPartOne;
    vPIDPartTwo = ConfPIDPartTwo;
  }

  // PID
  // femtoWorldPairTaskTrackPhi.cxx
  bool IsPionNSigma(float mom, float nsigmaTPCPi, float nsigmaTOFPi)
  {
    //|nsigma_TPC| < 3 for p < 0.5 GeV/c
    //|nsigma_combined| < 3 for p > 0.5

    // using configurables:
    // ConfNsigmaTPCPion -> TPC Kaon Sigma for momentum < 0.5
    // ConfNsigmaCombinedPion -> TPC and TOF Pion Sigma (combined) for momentum > 0.5
    if (true) {
      if (mom < 0.5) {
        if (TMath::Abs(nsigmaTPCPi) < ConfNsigmaTPCPion) {
          return true;
        } else {
          return false;
        }
      } else if (mom > 0.5) {
        if (TMath::Hypot(nsigmaTOFPi, nsigmaTPCPi) < ConfNsigmaCombinedPion) {
          return true;
        } else {
          return false;
        }
      }
    }
    return false;
  }

  // Function to build combinations
  template <typename T1, typename T2, typename T3, typename T4>
  void CombineParticles(T1 groupPartsOne, T1 groupPartsTwo, T2 cont, T3 parts, T4 magFieldTesla, int multCol, int sameOrMixed)
  {
    if (sameOrMixed == 1) {
      for (auto& [p1, p2] : combinations(groupPartsOne, groupPartsTwo)) {
        if (!(IsPionNSigma(p1.p(), p1.tpcNSigmaPi(), p1.tofNSigmaPi()))) {
          continue;
        }
        if (!(IsPionNSigma(p2.p(), p2.tpcNSigmaPi(), p2.tofNSigmaPi()))) {
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
        cont.setPair(p1, p2, multCol);
      }
    } else if (sameOrMixed == 2) {
      for (auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(groupPartsOne, groupPartsTwo))) {
        if (!(IsPionNSigma(p1.p(), p1.tpcNSigmaPi(), p1.tofNSigmaPi()))) {
          continue;
        }
        if (!(IsPionNSigma(p2.p(), p2.tpcNSigmaPi(), p2.tofNSigmaPi()))) {
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
        cont.setPair(p1, p2, multCol);
      }
    }
  }

  /// This function processes the same event and takes care of all the histogramming
  void processSameEvent(o2::aod::FemtoWorldCollision& col,
                        o2::aod::FemtoWorldParticlesMerged& parts)
  {
    auto groupPartsOne = partsOne->sliceByCached(aod::femtoworldparticle::femtoWorldCollisionId, col.globalIndex(), cache);
    auto groupPartsTwo = partsTwo->sliceByCached(aod::femtoworldparticle::femtoWorldCollisionId, col.globalIndex(), cache);
    auto groupPartsOneFailed = partsOneFailed->sliceByCached(aod::femtoworldparticle::femtoWorldCollisionId, col.globalIndex(), cache);
    auto groupPartsTwoFailed = partsTwoFailed->sliceByCached(aod::femtoworldparticle::femtoWorldCollisionId, col.globalIndex(), cache);

    const auto& magFieldTesla = col.magField();
    const int multCol = col.multV0M();
    eventHisto.fillQA(col);
    MixQaRegistry.fill(HIST("MixingQA/hSECollisionBins"), colBinning.getBin({col.posZ(), col.multV0M()}));
    /// Histogramming same event
    for (auto& part : groupPartsOne) {
      if (!(IsPionNSigma(part.p(), part.tpcNSigmaPi(), part.tofNSigmaPi()))) {
        continue;
      }
      trackHistoPartOne.fillQA(part);
    }
    for (auto& part : groupPartsOneFailed) {

      trackHistoPartOneFailed.fillQA(part);
    }
    if (!ConfIsSame) {
      for (auto& part : groupPartsTwo) {
        if (!(IsPionNSigma(part.p(), part.tpcNSigmaPi(), part.tofNSigmaPi()))) {
          continue;
        }
        trackHistoPartTwo.fillQA(part);
      }
      for (auto& part : groupPartsTwoFailed) {

        trackHistoPartTwoFailed.fillQA(part);
      }
    }
    /// Build the combinations
    if (ConfPlusMinus)
      CombineParticles(groupPartsOne, groupPartsTwo, sameEventPlusMinCont, parts, magFieldTesla, multCol, 1);
    if (ConfPlusPlus)
      CombineParticles(groupPartsOne, groupPartsOne, sameEventPlusPlusCont, parts, magFieldTesla, multCol, 1);
    if (ConfMinusMinus)
      CombineParticles(groupPartsTwo, groupPartsTwo, sameEventMinMinCont, parts, magFieldTesla, multCol, 1);
  }

  PROCESS_SWITCH(femtoWorldPairTaskPionPion, processSameEvent, "Enable processing same event", true);

  /// This function processes the mixed event
  void processMixedEvent(o2::aod::FemtoWorldCollisions& cols,
                         o2::aod::FemtoWorldParticlesMerged& parts)
  {

    for (auto& [collision1, collision2] : soa::selfCombinations(colBinning, 5, -1, cols, cols)) {

      MixQaRegistry.fill(HIST("MixingQA/hMECollisionBins"), colBinning.getBin({collision1.posZ(), collision1.multV0M()}));

      const auto& magFieldTesla1 = collision1.magField();
      const auto& magFieldTesla2 = collision2.magField();

      if (magFieldTesla1 != magFieldTesla2) {
        continue;
      }

      const int multCol1 = collision1.multV0M();

      /// \todo before mixing we should check whether both collisions contain a pair of particles!
      // if (partsOne.size() == 0 || nPart2Evt1 == 0 || nPart1Evt2 == 0 || partsTwo.size() == 0 ) continue;

      if (ConfPlusMinus) {
        auto groupPartsOne = partsOne->sliceByCached(aod::femtoworldparticle::femtoWorldCollisionId, collision1.globalIndex(), cache);
        auto groupPartsTwo = partsTwo->sliceByCached(aod::femtoworldparticle::femtoWorldCollisionId, collision2.globalIndex(), cache);
        CombineParticles(groupPartsOne, groupPartsTwo, mixedEventPlusMinCont, parts, magFieldTesla1, multCol1, 2);
      }
      if (ConfPlusPlus) {
        auto groupPartsOne = partsOne->sliceByCached(aod::femtoworldparticle::femtoWorldCollisionId, collision1.globalIndex(), cache);
        auto groupPartsTwo = partsOne->sliceByCached(aod::femtoworldparticle::femtoWorldCollisionId, collision2.globalIndex(), cache);
        CombineParticles(groupPartsOne, groupPartsTwo, mixedEventPlusPlusCont, parts, magFieldTesla1, multCol1, 2);
      }
      if (ConfMinusMinus) {
        auto groupPartsOne = partsTwo->sliceByCached(aod::femtoworldparticle::femtoWorldCollisionId, collision1.globalIndex(), cache);
        auto groupPartsTwo = partsTwo->sliceByCached(aod::femtoworldparticle::femtoWorldCollisionId, collision2.globalIndex(), cache);
        CombineParticles(groupPartsOne, groupPartsTwo, mixedEventMinMinCont, parts, magFieldTesla1, multCol1, 2);
      }
    }
  }

  PROCESS_SWITCH(femtoWorldPairTaskPionPion, processMixedEvent, "Enable processing mixed events", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<femtoWorldPairTaskPionPion>(cfgc),
  };
  return workflow;
}
