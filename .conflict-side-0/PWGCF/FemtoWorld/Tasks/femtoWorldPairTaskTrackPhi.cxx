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

/// \file femtoWorldPairTaskTrackTrack.cxx
/// \brief Tasks that reads the track tables used for the pairing and builds pairs of two tracks
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de
/// \author Zuzanna Chochulska, WUT Warsaw, zchochul@cern.ch

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
#include "PWGCF/FemtoWorld/Core/FemtoWorldContainer.h"
#include "PWGCF/FemtoWorld/Core/FemtoWorldDetaDphiStar.h"
#include "PWGCF/FemtoWorld/Core/FemtoWorldUtils.h"

using namespace o2;
using namespace o2::analysis::femtoWorld;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

namespace
{
static constexpr int nPart = 2;                                                                         // number of particle types (for us it will be proton and phi)
static constexpr int nCuts = 5;                                                                         // number of cuts
static const std::vector<std::string> partNames{"PartOne", "PartTwo"};                                  // names of the rows
static const std::vector<std::string> cutNames{"MaxPt", "PIDthr", "nSigmaTPC", "nSigmaTPCTOF", "MaxP"}; // names of the columns
static const float cutsTable[nPart][nCuts]{                                                             // cuts table [rows = particles][columns = types of cuts]
                                           {1.5f, 1.f, 3.f, 3.f, 100.f},
                                           {1.5f, 1.f, 5.f, 5.f, 100.f}};

static const std::vector<float> kNsigma = {3.5f, 3.f, 2.5f};

} // namespace

struct femtoWorldPairTaskTrackPhi {
  SliceCache cache;
  Preslice<aod::FemtoWorldParticles> perCol = aod::femtoworldparticle::femtoWorldCollisionId;

  /// Particle selection part
  Configurable<int> ConfTrackChoice{"ConfTrackChoice", 0, "Type of particle (track1): {0:Proton, 1:Pion, 2:Kaon}"};
  Configurable<float> ConfNsigmaCombinedKaon{"ConfNsigmaCombinedKaon", 3.0, "TPC and TOF Kaon Sigma (combined) for momentum > 0.4"};
  Configurable<float> ConfNsigmaTPCKaon{"ConfNsigmaTPCKaon", 3.0, "TPC Kaon Sigma for momentum < 0.4"};
  Configurable<float> ConfNsigmaCombinedProton{"ConfNsigmaCombinedProton", 3.0, "TPC and TOF Proton Sigma (combined) for momentum > 0.5"};
  Configurable<float> ConfNsigmaTPCProton{"ConfNsigmaTPCProton", 3.0, "TPC Proton Sigma for momentum < 0.5"};
  Configurable<float> ConfNsigmaCombinedPion{"ConfNsigmaCombinedPion", 3.0, "TPC and TOF Pion Sigma (combined) for momentum > 0.5"};
  Configurable<float> ConfNsigmaTPCPion{"ConfNsigmaTPCPion", 3.0, "TPC Pion Sigma for momentum < 0.5"};
  /// Table for both particles
  Configurable<LabeledArray<float>> cfgCutTable{"cfgCutTable", {cutsTable[0], nPart, nCuts, partNames, cutNames}, "Particle selections"};
  Configurable<int> cfgNspecies{"ccfgNspecies", 4, "Number of particle spieces with PID info"};

  /// Particle 1 (track)
  Configurable<int> ConfPDGCodePartOne{"ConfPDGCodePartOne", 2212, "Particle 1 - PDG code"};                                 // proton (2212)
  Configurable<std::vector<int>> ConfPIDPartOne{"ConfPIDPartOne", std::vector<int>{2}, "Particle 1 - Read from cutCulator"}; // we also need the possibility to specify whether the bit is true/false ->std>>vector<std::pair<int, int>>int>>
  Configurable<float> cfgPtLowPart1{"cfgPtLowPart1", 0.5, "Lower limit for Pt for the first particle"};                      // change according to wrzesa cuts
  Configurable<float> cfgPtHighPart1{"cfgPtHighPart1", 1.5, "Higher limit for Pt for the first particle"};
  Configurable<float> cfgEtaLowPart1{"cfgEtaLowPart1", -0.8, "Lower limit for Eta for the first particle"};
  Configurable<float> cfgEtaHighPart1{"cfgEtaHighPart1", 0.8, "Higher limit for Eta for the first particle"};
  Configurable<float> cfgDcaXYPart1{"cfgDcaXYPart1", 2.4, "Value for DCA_XY for the first particle"};
  Configurable<float> cfgDcaZPart1{"cfgDcaZPart1", 3.2, "Value for DCA_Z for the first particle"};
  Configurable<int> cfgTpcClPart1{"cfgTpcClPart1", 88, "Number of tpc clasters for the first particle"};             // min number of found TPC clusters
  Configurable<int> cfgTpcCrosRoPart1{"cfgTpcCrosRoPart1", 70, "Number of tpc crossed rows for the first particle"}; // min number of crossed rows
  Configurable<float> cfgChi2TpcPart1{"cfgChi2TpcPart1", 4.0, "Chi2 / cluster for the TPC track segment for the first particle"};
  Configurable<float> cfgChi2ItsPart1{"cfgChi2ItsPart1", 36.0, "Chi2 / cluster for the ITS track segment for the first particle"};

  /// Partition for particle 1 (proton)
  Partition<aod::FemtoWorldParticles> partsOne = (aod::femtoworldparticle::partType == uint8_t(aod::femtoworldparticle::ParticleType::kTrack)) && (aod::femtoworldparticle::pt < cfgPtHighPart1) && (aod::femtoworldparticle::pt > cfgPtLowPart1) // simple pT cuts
                                                 && (aod::femtoworldparticle::eta < cfgEtaHighPart1) && (aod::femtoworldparticle::eta > cfgEtaLowPart1)                                                                                           // Eta cuts
                                                 && (o2::aod::track::dcaXY < cfgDcaXYPart1) && (o2::aod::track::dcaZ < cfgDcaZPart1)                                                                                                              // DCA cuts for XY and Z
                                                 && (aod::femtoworldparticle::tpcNClsFound > (uint8_t)cfgTpcClPart1)                                                                                                                              // Number of found TPC clusters
                                                 && (aod::femtoworldparticle::tpcNClsCrossedRows > (uint8_t)cfgTpcCrosRoPart1)                                                                                                                    // Crossed rows TPC
                                                 && (aod::femtoworldparticle::itsChi2NCl < cfgChi2ItsPart1) && (aod::femtoworldparticle::tpcChi2NCl < cfgChi2TpcPart1)                                                                            //&& // chi2 cuts
                                                 && (aod::femtoworldparticle::sign > int8_t(0));

  /// Histogramming for particle 1
  FemtoWorldParticleHisto<aod::femtoworldparticle::ParticleType::kTrack, 0> trackHistoPartOne;

  /// Particle 2 (Phi)
  Configurable<int> ConfPDGCodePartTwo{"ConfPDGCodePartTwo", 333, "Particle 1 - PDG code"}; // phi meson (333)
  Configurable<uint32_t> ConfCutPartTwo{"ConfCutPartTwo", 338, "Particle 2 - Selection bit"};
  Configurable<float> cfgPtLowPart2{"cfgPtLowPart2", 0.14, "Lower limit for Pt for the second particle"};
  Configurable<float> cfgPtHighPart2{"cfgPtHighPart2", 5.0, "Higher limit for Pt for the second particle"};
  Configurable<float> cfgEtaLowPart2{"cfgEtaLowPart2", -0.8, "Lower limit for Eta for the second particle"};
  Configurable<float> cfgEtaHighPart2{"cfgEtaHighPart2", 0.8, "Higher limit for Eta for the second particle"};

  /// Partition for particle 2
  Partition<aod::FemtoWorldParticles> partsTwo = (aod::femtoworldparticle::partType == uint8_t(aod::femtoworldparticle::ParticleType::kPhi)) && (aod::femtoworldparticle::pt < cfgPtHighPart2) && (aod::femtoworldparticle::pt > cfgPtLowPart2) // pT cuts
                                                 && (aod::femtoworldparticle::eta < cfgEtaHighPart2) && (aod::femtoworldparticle::eta > cfgEtaLowPart2);                                                                                        // Eta cuts

  /// Histogramming for particle 2
  FemtoWorldParticleHisto<aod::femtoworldparticle::ParticleType::kPhi, 0> trackHistoPartTwo;

  // Partition for particle 3 (Phi daughters (K+ K-))
  Partition<aod::FemtoWorldParticles> partsThree = (aod::femtoworldparticle::partType == uint8_t(aod::femtoworldparticle::ParticleType::kPhiChild));

  /// Histogramming for particle 3
  FemtoWorldParticleHisto<aod::femtoworldparticle::ParticleType::kPhiChild, 0> trackHistoPartThree;

  /// Histogramming for Event
  FemtoWorldEventHisto eventHisto;

  /// The configurables need to be passed to an std::vector
  std::vector<int> vPIDPartOne;

  /// Correlation part
  ConfigurableAxis CfgMultBins{"CfgMultBins", {VARIABLE_WIDTH, 0.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.0f, 150.0f, 200.0f, 99999.f}, "Mixing bins - multiplicity"};
  ConfigurableAxis CfgVtxBins{"CfgVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis CfgkstarBins{"CfgkstarBins", {1500, 0., 6.}, "binning kstar"};
  ConfigurableAxis CfgkTBins{"CfgkTBins", {150, 0., 9.}, "binning kT"};
  ConfigurableAxis CfgmTBins{"CfgmTBins", {225, 0., 7.5}, "binning mT"};
  Configurable<int> ConfPhiBins{"ConfPhiBins", 29, "Number of phi bins in deta dphi"};
  Configurable<int> ConfEtaBins{"ConfEtaBins", 29, "Number of eta bins in deta dphi"};
  Configurable<int> ConfMInvBins{"ConfMInvBins", 1000, "Number of bins in mInv distribution"};
  Configurable<int> ConfNEventsMix{"ConfNEventsMix", 5, "Number of events for mixing"};
  Configurable<bool> ConfIsCPR{"ConfIsCPR", true, "Close Pair Rejection"};
  Configurable<bool> ConfCPRPlotPerRadii{"ConfCPRPlotPerRadii", false, "Plot CPR per radii"};

  FemtoWorldContainer<femtoWorldContainer::EventType::same, femtoWorldContainer::Observable::kstar> sameEventCont;
  FemtoWorldContainer<femtoWorldContainer::EventType::mixed, femtoWorldContainer::Observable::kstar> mixedEventCont;
  FemtoWorldPairCleaner<aod::femtoworldparticle::ParticleType::kTrack, aod::femtoworldparticle::ParticleType::kPhi> pairCleaner;
  FemtoWorldDetaDphiStar<aod::femtoworldparticle::ParticleType::kTrack, aod::femtoworldparticle::ParticleType::kPhi> pairCloseRejection;
  /// Histogram output
  HistogramRegistry qaRegistry{"TrackQA", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry resultRegistry{"Correlations", {}, OutputObjHandlingPolicy::AnalysisObject};

  // PID for protons
  bool IsProtonNSigma(float mom, float nsigmaTPCPr, float nsigmaTOFPr) // previous version from: https://github.com/alisw/AliPhysics/blob/master/PWGCF/FEMTOSCOPY/AliFemtoUser/AliFemtoMJTrackCut.cxx
  {
    //|nsigma_TPC| < 3 for p < 0.5 GeV/c
    //|nsigma_combined| < 3 for p > 0.5

    // using configurables:
    // ConfNsigmaTPCProton -> TPC Kaon Sigma for momentum < 0.5
    // ConfNsigmaCombinedProton -> TPC and TOF Kaon Sigma (combined) for momentum > 0.5

    if (mom < 0.5) {
      if (TMath::Abs(nsigmaTPCPr) < ConfNsigmaTPCProton) {
        return true;
      } else {
        return false;
      }
    } else if (mom > 0.4) {
      if (TMath::Hypot(nsigmaTOFPr, nsigmaTPCPr) < ConfNsigmaCombinedProton) {
        return true;
      } else {
        return false;
      }
    }
    return false;
  }

  bool IsKaonNSigma(float mom, float nsigmaTPCK, float nsigmaTOFK)
  {
    //|nsigma_TPC| < 3 for p < 0.5 GeV/c
    //|nsigma_combined| < 3 for p > 0.5

    // using configurables:
    // ConfNsigmaTPCTOFKaon -> are we doing TPC TOF PID for Kaons? (boolean)
    // ConfNsigmaTPCKaon -> TPC Kaon Sigma for momentum < 0.5
    // ConfNsigmaCombinedKaon -> TPC and TOF Kaon Sigma (combined) for momentum > 0.5
    if (true) {
      if (mom < 0.5) {
        if (TMath::Abs(nsigmaTPCK) < ConfNsigmaTPCKaon) {
          return true;
        } else {
          return false;
        }
      } else if (mom > 0.5) {
        if (TMath::Hypot(nsigmaTOFK, nsigmaTPCK) < ConfNsigmaCombinedKaon) {
          return true;
        } else {
          return false;
        }
      }
    }
    return false;
  }

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

  bool IsParticleNSigma(float mom, float nsigmaTPCPr, float nsigmaTOFPr, float nsigmaTPCPi, float nsigmaTOFPi, float nsigmaTPCK, float nsigmaTOFK)
  {
    switch (ConfTrackChoice) {
      case 0: // Proton
        return IsProtonNSigma(mom, nsigmaTPCPr, nsigmaTOFPr);
        break;
      case 1: // Pion
        return IsPionNSigma(mom, nsigmaTPCPi, nsigmaTOFPi);
        break;
      case 2: // Kaon
        return IsKaonNSigma(mom, nsigmaTPCK, nsigmaTOFK);
        break;
      default:
        return false;
    }
    return false;
  }

  void init(InitContext&)
  {
    eventHisto.init(&qaRegistry);
    trackHistoPartOne.init(&qaRegistry);
    trackHistoPartTwo.init(&qaRegistry);
    trackHistoPartThree.init(&qaRegistry);

    sameEventCont.init(&resultRegistry, CfgkstarBins, CfgMultBins, CfgkTBins, CfgmTBins, ConfPhiBins, ConfEtaBins, ConfMInvBins);
    sameEventCont.setPDGCodes(ConfPDGCodePartOne, ConfPDGCodePartTwo);
    mixedEventCont.init(&resultRegistry, CfgkstarBins, CfgMultBins, CfgkTBins, CfgmTBins, ConfPhiBins, ConfEtaBins, ConfMInvBins);
    mixedEventCont.setPDGCodes(ConfPDGCodePartOne, ConfPDGCodePartTwo);
    pairCleaner.init(&qaRegistry);
    if (ConfIsCPR) {
      pairCloseRejection.init(&resultRegistry, &qaRegistry, 0.01, 0.01, ConfCPRPlotPerRadii); /// \todo add config for Δη and ΔΦ cut values
    }

    vPIDPartOne = ConfPIDPartOne;
  }

  /// This function processes the same event and takes care of all the histogramming
  /// \todo the trivial loops over the tracks should be factored out since they will be common to all combinations of T-T, T-Phi, Phi-Phi, ...
  void processSameEvent(o2::aod::FemtoWorldCollision& col,
                        o2::aod::FemtoWorldParticles&)
  {
    // const auto& magFieldTesla = col.magField();
    auto groupPartsOne = partsOne->sliceByCached(aod::femtoworldparticle::femtoWorldCollisionId, col.globalIndex(), cache);
    auto groupPartsTwo = partsTwo->sliceByCached(aod::femtoworldparticle::femtoWorldCollisionId, col.globalIndex(), cache);
    auto groupPartsThree = partsThree->sliceByCached(aod::femtoworldparticle::femtoWorldCollisionId, col.globalIndex(), cache);

    const int multCol = col.multV0M();
    eventHisto.fillQA(col);

    /// Histogramming same event
    for (auto& part : groupPartsOne) {
      if (!(IsParticleNSigma(part.p(), part.tpcNSigmaPr(), part.tofNSigmaPr(), part.tpcNSigmaPi(), part.tofNSigmaPi(), part.tpcNSigmaKa(), part.tofNSigmaKa()))) {
        continue;
      }
      trackHistoPartOne.fillQA(part);
    }
    for (auto& part : groupPartsTwo) {
      trackHistoPartTwo.fillQAMult(part, multCol);
    }
    for (auto& part : groupPartsThree) {
      trackHistoPartThree.fillQA(part);
    }

    /// Now build the combinations
    for (auto& [p1, p2] : combinations(groupPartsOne, groupPartsTwo)) {
      if (!(IsParticleNSigma(p1.p(), p1.tpcNSigmaPr(), p1.tofNSigmaPr(), p1.tpcNSigmaPi(), p1.tofNSigmaPi(), p1.tpcNSigmaKa(), p1.tofNSigmaKa()))) {
        continue;
      }
      // TODO: Include pairCloseRejection and pairCleaner
      /*
      if (ConfIsCPR) {
        if (pairCloseRejection.isClosePair(p1, p2, parts, magFieldTesla)) {
          continue;
        }
      }
      */
      // track cleaning
      /*if (!pairCleaner.isCleanPair(p1, p2, parts)) {
        continue;
      }*/

      sameEventCont.setPair(p1, p2, multCol);
    }
  }

  PROCESS_SWITCH(femtoWorldPairTaskTrackPhi, processSameEvent, "Enable processing same event", true);

  /// This function processes the mixed event
  /// \todo the trivial loops over the collisions and tracks should be factored out since they will be common to all combinations of T-T, T-Phi, Phi-Phi, ...
  void processMixedEvent(o2::aod::FemtoWorldCollisions& cols,
                         o2::aod::FemtoWorldParticles&)
  {
    ColumnBinningPolicy<aod::collision::PosZ, aod::femtoworldcollision::MultV0M> colBinning{{CfgVtxBins, CfgMultBins}, true};

    for (auto& [collision1, collision2] : soa::selfCombinations(colBinning, ConfNEventsMix, -1, cols, cols)) {

      auto groupPartsOne = partsOne->sliceByCached(aod::femtoworldparticle::femtoWorldCollisionId, collision1.globalIndex(), cache);
      auto groupPartsTwo = partsTwo->sliceByCached(aod::femtoworldparticle::femtoWorldCollisionId, collision2.globalIndex(), cache);
      auto groupPartsThree = partsThree->sliceByCached(aod::femtoworldparticle::femtoWorldCollisionId, collision2.globalIndex(), cache);

      const auto& magFieldTesla1 = collision1.magField();
      const auto& magFieldTesla2 = collision2.magField();

      if (magFieldTesla1 != magFieldTesla2) {
        continue;
      }

      for (auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(groupPartsOne, groupPartsTwo))) {
        if (!(IsParticleNSigma(p1.p(), p1.tpcNSigmaPr(), p1.tofNSigmaPr(), p1.tpcNSigmaPi(), p1.tofNSigmaPi(), p1.tpcNSigmaKa(), p1.tofNSigmaKa()))) {
          continue;
        }
        // TODO: Include pairCloseRejection and pairCleaner
        /*
        if (ConfIsCPR) {
          if (pairCloseRejection.isClosePair(p1, p2, parts, magFieldTesla1)) {
            continue;
          }
        }*/
        mixedEventCont.setPair(p1, p2, collision1.multV0M());
      }
    }
  }

  PROCESS_SWITCH(femtoWorldPairTaskTrackPhi, processMixedEvent, "Enable processing mixed events", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<femtoWorldPairTaskTrackPhi>(cfgc),
  };
  return workflow;
}
