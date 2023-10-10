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

/// \file FemtoWorldPionAllPair.cxx
/// \brief Tasks that reads the track tables used for the pairing and builds pairs of two tracks
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de
/// \author Zuzanna Chochulska, WUT Warsaw, zchochul@cern.ch
/// \author Deependra Sharma, IITB, deependra.sharma@cern.ch

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/StepTHn.h"
#include "TDatabasePDG.h"

#include "PWGCF/FemtoWorld/DataModel/FemtoWorldDerived.h"
#include "PWGCF/FemtoWorld/Core/FemtoWorldParticleHisto.h"
#include "PWGCF/FemtoWorld/Core/FemtoWorldEventHisto.h"
#include "PWGCF/FemtoWorld/Core/FemtoWorldPairCleaner.h"
#include "PWGCF/FemtoWorld/Core/FemtoWorldContainer.h"
#include "PWGCF/FemtoWorld/Core/FemtoWorldDetaDphiStar.h"
#include "PWGCF/FemtoWorld/Core/FemtoWorldUtils.h"
#include "PWGCF/FemtoWorld/Core/FemtoWorldMath.h"
#include "PWGCF/FemtoWorld/Core/FemtoWorldPairWithCentrality.h"

using namespace o2;
using namespace o2::analysis::femtoWorld;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

namespace o2::aod
{
// using FemtoWorldParticlesMerged = soa::Join<aod::FemtoWorldParticles,aod::FemtoWorldDebugParticles>;
using FemtoWorldParticlesMerged = aod::FemtoWorldParticles;
using FemtoWorldParticleMerged = FemtoWorldParticlesMerged::iterator;
} // namespace o2::aod

struct FemtoWorldIdenticalPionPair {
  SliceCache cache;
  Preslice<aod::FemtoWorldParticlesMerged> perCol = aod::femtoworldparticle::femtoWorldCollisionId;

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

  /// Particle 1 Pi+
  Configurable<int> ConfPDGCodePartOne{"ConfPDGCodePartOne", 211, "Particle 1 - PDG code"};

  /// Particle 2 Pi-
  Configurable<bool> ConfIsSame{"ConfIsSame", false, "Pairs of the same particle"};
  Configurable<int> ConfPDGCodePartTwo{"ConfPDGCodePartTwo", 211, "Particle 2 - PDG code"};
  Configurable<int> cfgSignPartOne{"cfgSignPartOne", 1, "Sign of Firts particle"};
  Configurable<int> cfgSignPartTwo{"cfgSignPartTwo", 1, "Sign of Second particle"};

  /// Partition for particle 1
  Partition<aod::FemtoWorldParticlesMerged> partsOne = (aod::femtoworldparticle::partType == uint8_t(aod::femtoworldparticle::ParticleType::kTrack))                          // particle type cut
                                                       && (aod::femtoworldparticle::pt < cfgPtHighPart1) && (aod::femtoworldparticle::pt > cfgPtLowPart1)                     // simple pT cuts
                                                       && (aod::femtoworldparticle::eta < cfgEtaHighPart1) && (aod::femtoworldparticle::eta > cfgEtaLowPart1)                 // Eta cuts
                                                       && (nabs(o2::aod::track::dcaXY) < cfgDcaXYPart1) && (nabs(o2::aod::track::dcaZ) < cfgDcaZPart1)                        // DCA cuts for XY and Z
                                                       && (aod::femtoworldparticle::tpcNClsFound > (uint8_t)cfgTpcClPart1)                                                    // Number of found TPC clusters
                                                       && (aod::femtoworldparticle::tpcNClsCrossedRows > (uint8_t)cfgTpcCrosRoPart1)                                          // Crossed rows TPC
                                                       && (aod::femtoworldparticle::itsChi2NCl < cfgChi2ItsPart1) && (aod::femtoworldparticle::tpcChi2NCl < cfgChi2TpcPart1); //&& // chi2 cuts

  /// Histogramming for particle 1
  FemtoWorldParticleHisto<aod::femtoworldparticle::ParticleType::kTrack, 1> trackHistoPartOne;

  /// Histogramming for particle 2
  FemtoWorldParticleHisto<aod::femtoworldparticle::ParticleType::kTrack, 2> trackHistoPartTwo;

  /// Histogramming for Event
  FemtoWorldEventHisto eventHisto;

  /// Correlation part
  // ConfigurableAxis CfgMultBins{"CfgMultBins", {VARIABLE_WIDTH, 0.0f, 4.0f, 8.0f, 12.0f, 16.0f, 20.0f, 24.0f, 28.0f, 32.0f, 36.0f, 40.0f, 44.0f, 48.0f, 52.0f, 56.0f, 60.0f, 64.0f, 68.0f, 72.0f, 76.0f, 80.0f, 84.0f, 88.0f, 92.0f, 96.0f, 100.0f, 200.0f, 99999.f}, "Mixing bins - multiplicity"}; // \todo to be obtained from the hash task
  ConfigurableAxis CfgMultBins{"CfgMultBins", {VARIABLE_WIDTH, 0.0f, 20.0f, 40.0f, 60.0f, 80.0f, 100.0f, 200.0f, 99999.f}, "Mixing bins - multiplicity"};
  ConfigurableAxis CfgVtxBins{"CfgVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis CfgCentBins{"CfgCentBins", {VARIABLE_WIDTH, 0.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f}, "Centrality Bins"};
  ConfigurableAxis CfgCentBinsMixing{"CfgCentBinsMixing", {VARIABLE_WIDTH, 0.0f, 4.0f, 8.0f, 12.0f, 16.0f, 20.0f, 24.0f, 28.0f, 32.0f, 36.0f, 40.0f, 44.0f, 48.0f, 52.0f, 56.0f, 60.0f, 64.0f, 68.0f, 72.0f, 76.0f, 80.0f, 84.0f, 88.0f, 92.0f, 96.0f, 100.0f}, "Mixing-Bins Centrality Bins"};

  // Configurable<std::vector<float>> CfgMultBins{"CfgMultBins",{0.0f, 20.0f, 40.0f, 60.0f, 80.0f, 100.0f, 200.0f, 99999.f},"Mixing bins - multiplicity"};
  // Configurable<std::vector<float>> CfgVtxBins{"CfgVtxBins",{-10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f},"Mixing bins - z-vertex"};
  // Configurable<std::vector<float>> CfgCentBins{"CfgCentBins",{5.0f,10.0f,20.0f,30.0f,40.0f,50.0f},"Centrality Bins"};

  ColumnBinningPolicy<aod::collision::PosZ, aod::femtoworldcollision::MultV0M> colBinning{{CfgVtxBins, CfgMultBins}, true};
  ColumnBinningPolicy<aod::collision::PosZ, aod::femtoworldcollision::RunCentrality> colBinning2{{CfgVtxBins, CfgCentBinsMixing}, true};

  ConfigurableAxis CfgkstarBins{"CfgkstarBins", {1500, 0., 6.}, "binning kstar"};
  ConfigurableAxis CfgkTBins{"CfgkTBins", {150, 0., 9.}, "binning kT"};
  ConfigurableAxis CfgmTBins{"CfgmTBins", {225, 0., 7.5}, "binning mT"};
  Configurable<int> ConfNEventsMix{"ConfNEventsMix", 5, "Number of events for mixing"};
  Configurable<bool> ConfIsCPR{"ConfIsCPR", true, "Close Pair Rejection"};
  Configurable<bool> ConfCPRPlotPerRadii{"ConfCPRPlotPerRadii", false, "Plot CPR per radii"};
  Configurable<int> ConfPhiBins{"ConfPhiBins", 15, "Number of phi bins in deta dphi"};
  Configurable<int> ConfEtaBins{"ConfEtaBins", 15, "Number of eta bins in deta dphi"};
  Configurable<int> ConfMInvBins{"ConfMInvBins", 1000, "Number of bins in mInv distribution"};

  FemtoWorldContainer<femtoWorldContainer::EventType::same, femtoWorldContainer::Observable::kstar> sameEventCont;
  FemtoWorldContainer<femtoWorldContainer::EventType::mixed, femtoWorldContainer::Observable::kstar> mixedEventCont;
  FemtoWorldPairCleaner<aod::femtoworldparticle::ParticleType::kTrack, aod::femtoworldparticle::ParticleType::kTrack> pairCleaner;
  FemtoWorldDetaDphiStar<aod::femtoworldparticle::ParticleType::kTrack, aod::femtoworldparticle::ParticleType::kTrack> pairCloseRejection;
  PairWithCentrality SameEvCorrWithCent;
  PairWithCentrality MixedEvCorrWithCent;

  /// Histogram output
  HistogramRegistry qaRegistry{"TrackQA", {}, OutputObjHandlingPolicy::AnalysisObject};
  // HistogramRegistry qaRegistryFail{"TrackQAFailed", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry resultRegistry{"Correlations", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry MixQaRegistry{"MixQaRegistry", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry SameEvPairWithCentReg{"SameEvPairWithCentReg", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry MixedEvPairWithCentReg{"MixedEvPairWithCentReg", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  void init(InitContext&)
  {
    qaRegistry.add("ChargeInfo/TrackOne", "Sign of track one", kTH1F, {{4, -2, 2}});
    qaRegistry.add("ChargeInfo/TrackTwo", "Sign of track two", kTH1F, {{4, -2, 2}});

    mMassOne = TDatabasePDG::Instance()->GetParticle(ConfPDGCodePartOne)->Mass();
    mMassTwo = TDatabasePDG::Instance()->GetParticle(ConfPDGCodePartTwo)->Mass();

    SameEvCorrWithCent.init(&SameEvPairWithCentReg, CfgkstarBins, CfgCentBins);
    MixedEvCorrWithCent.init(&MixedEvPairWithCentReg, CfgkstarBins, CfgCentBins);

    eventHisto.init(&qaRegistry);
    trackHistoPartOne.init(&qaRegistry);

    MixQaRegistry.add("MixingQA/hSECollisionBins", ";bin;Entries", kTH1F, {{120, -0.5, 119.5}});
    MixQaRegistry.add("MixingQA/hMECollisionBins", ";bin;Entries", kTH1F, {{120, -0.5, 119.5}});

    sameEventCont.init(&resultRegistry, CfgkstarBins, CfgMultBins, CfgkTBins, CfgmTBins, ConfPhiBins, ConfEtaBins, ConfMInvBins);
    sameEventCont.setPDGCodes(ConfPDGCodePartOne, ConfPDGCodePartTwo);
    mixedEventCont.init(&resultRegistry, CfgkstarBins, CfgMultBins, CfgkTBins, CfgmTBins, ConfPhiBins, ConfEtaBins, ConfMInvBins);
    mixedEventCont.setPDGCodes(ConfPDGCodePartOne, ConfPDGCodePartTwo);
    pairCleaner.init(&qaRegistry);
    if (ConfIsCPR) {
      pairCloseRejection.init(&resultRegistry, &qaRegistry, 0.01, 0.01, ConfCPRPlotPerRadii); /// \todo add config for Δη and ΔΦ cut values
    }
  }

  bool IsPionNSigma(float nsigmap)
  {
    if (TMath::Abs(nsigmap) < 3.0)
      return true;
    return false;
  }

  void processSameEvent(o2::aod::FemtoWorldCollision& col, o2::aod::FemtoWorldParticlesMerged& parts)
  {

    auto groupPartsOne = partsOne->sliceByCached(aod::femtoworldparticle::femtoWorldCollisionId, col.globalIndex(), cache);

    // auto groupPartsOne=groupParts.sliceByCached(aod::femtoworldparticle::sign,cfgSignPartOne);
    // auto groupPartsTwo=groupParts.sliceByCached(aod::femtoworldparticle::sign,cfgSignPartTwo);

    const auto& magFieldTesla = col.magField();
    const int multCol = col.multV0M();
    eventHisto.fillQA(col);
    MixQaRegistry.fill(HIST("MixingQA/hSECollisionBins"), colBinning2.getBin({col.posZ(), col.runCent()}));

    /// Histogramming same event
    for (auto& part : groupPartsOne) {
      if ((part.p() > 0.50)) {

        float NSigmaPion = TMath::Sqrt(2.0) * TMath::Sqrt(TMath::Sq(part.tpcNSigmaPi()) + TMath::Sq(part.tofNSigmaPi()));
        if (!IsPionNSigma(NSigmaPion)) {
          continue;
        }

      } else if ((part.p() <= 0.50)) {
        float NSigmaPion = part.tpcNSigmaPi();
        if (!IsPionNSigma(NSigmaPion)) {
          continue;
        }
      }

      if (part.sign() == cfgSignPartOne) {
        qaRegistry.fill(HIST("ChargeInfo/TrackOne"), part.sign());
        trackHistoPartOne.fillQA(part);
      }
      if (part.sign() == cfgSignPartTwo) {
        qaRegistry.fill(HIST("ChargeInfo/TrackTwo"), part.sign());
        trackHistoPartTwo.fillQA(part);
      }
    }

    // Pair Correlation
    for (auto& [p1, p2] : combinations(CombinationsStrictlyUpperIndexPolicy(groupPartsOne, groupPartsOne))) {
      if ((p1.sign() != cfgSignPartOne) || (p2.sign() != cfgSignPartTwo))
        continue;
      if ((p1.p() > 0.50)) {
        float NSigmaPion1 = TMath::Sqrt(2.0) * TMath::Sqrt(TMath::Sq(p1.tpcNSigmaPi()) + TMath::Sq(p1.tofNSigmaPi()));
        if (!IsPionNSigma(NSigmaPion1)) {
          continue;
        }

      } else if ((p1.p() <= 0.50)) {
        float NSigmaPion1 = p1.tpcNSigmaPi();
        if (!IsPionNSigma(NSigmaPion1)) {
          continue;
        }
      }
      if ((p2.p() > 0.50)) {
        float NSigmaPion2 = TMath::Sqrt(2.0) * TMath::Sqrt(TMath::Sq(p2.tpcNSigmaPi()) + TMath::Sq(p2.tofNSigmaPi()));
        if (!IsPionNSigma(NSigmaPion2)) {
          continue;
        }

      } else if ((p2.p() <= 0.50)) {
        float NSigmaPion2 = p2.tpcNSigmaPi();
        if (!IsPionNSigma(NSigmaPion2)) {
          continue;
        }
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

      sameEventCont.setPair(p1, p2, multCol);

      float kstar = FemtoWorldMath::getkstar(p1, mMassOne, p2, mMassTwo);
      // float kT=FemtoWorldMath::getkT(p1, mMassOne, p2, mMassTwo)
      float v0mCent = col.runCent();
      SameEvCorrWithCent.fill<float>(kstar, v0mCent);
    }
  }

  PROCESS_SWITCH(FemtoWorldIdenticalPionPair, processSameEvent, "Enable processing same event", true);

  /// This function processes the mixed event
  void processMixedEvent(o2::aod::FemtoWorldCollisions& cols, o2::aod::FemtoWorldParticlesMerged& parts)
  {
    for (auto& [collision1, collision2] : soa::selfCombinations(colBinning2, 5, -1, cols, cols)) {

      MixQaRegistry.fill(HIST("MixingQA/hMECollisionBins"), colBinning2.getBin({collision1.posZ(), collision1.runCent()}));

      auto groupPartsOne = partsOne->sliceByCached(aod::femtoworldparticle::femtoWorldCollisionId, collision1.globalIndex(), cache);
      auto groupPartsTwo = partsOne->sliceByCached(aod::femtoworldparticle::femtoWorldCollisionId, collision2.globalIndex(), cache);

      const auto& magFieldTesla1 = collision1.magField();
      const auto& magFieldTesla2 = collision2.magField();

      if (magFieldTesla1 != magFieldTesla2) {
        continue;
      }

      /// \todo before mixing we should check whether both collisions contain a pair of particles!
      // if (partsOne.size() == 0 || nPart2Evt1 == 0 || nPart1Evt2 == 0 || partsTwo.size() == 0 ) continue;

      for (auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(groupPartsOne, groupPartsTwo))) {
        if ((p1.sign() != cfgSignPartOne) || (p2.sign() != cfgSignPartTwo))
          continue;

        if ((p1.p() > 0.50)) {
          float NSigmaPion1 = TMath::Sqrt(2.0) * TMath::Sqrt(TMath::Sq(p1.tpcNSigmaPi()) + TMath::Sq(p1.tofNSigmaPi()));
          if (!IsPionNSigma(NSigmaPion1)) {
            continue;
          }

        } else if ((p1.p() <= 0.50)) {
          float NSigmaPion1 = p1.tpcNSigmaPi();
          if (!IsPionNSigma(NSigmaPion1)) {
            continue;
          }
        }
        if ((p2.p() > 0.50)) {
          float NSigmaPion2 = TMath::Sqrt(2.0) * TMath::Sqrt(TMath::Sq(p2.tpcNSigmaPi()) + TMath::Sq(p2.tofNSigmaPi()));
          if (!IsPionNSigma(NSigmaPion2)) {
            continue;
          }

        } else if ((p2.p() <= 0.50)) {
          float NSigmaPion2 = p2.tpcNSigmaPi();
          if (!IsPionNSigma(NSigmaPion2)) {
            continue;
          }
        }

        if (ConfIsCPR) {
          if (pairCloseRejection.isClosePair(p1, p2, parts, magFieldTesla1)) {
            continue;
          }
        }
        mixedEventCont.setPair(p1, p2, collision1.multV0M());

        float kstar = FemtoWorldMath::getkstar(p1, mMassOne, p2, mMassTwo);
        float v0mCent = collision1.runCent();
        // float kT=FemtoWorldMath::getkT(p1, mMassOne, p2, mMassTwo)
        MixedEvCorrWithCent.fill<float>(kstar, v0mCent);
      }
    }
  }
  PROCESS_SWITCH(FemtoWorldIdenticalPionPair, processMixedEvent, "Enable processing mixed events", true);

  float mMassOne = -999., mMassTwo = -999.;
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<FemtoWorldIdenticalPionPair>(cfgc),
  };
  return workflow;
}
