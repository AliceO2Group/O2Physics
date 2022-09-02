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

namespace o2::aod
{
// using FemtoWorldParticlesMerged = soa::Join<aod::FemtoWorldParticles,aod::FemtoWorldDebugParticles>;
using FemtoWorldParticlesMerged = aod::FemtoWorldParticles;
using FemtoWorldParticleMerged = FemtoWorldParticlesMerged::iterator;
} // namespace o2::aod

struct femtoWorldPairTaskTrackTrack {
  // for filling phi candidates table
  // Produces<aod::FemtoWorldCollisions> outputCollision;
  // Produces<aod::FemtoWorldPhiCandidates> outputPhiCan;
  // std::vector<int> tmpIDtrack; // this vector keeps track of the matching of the primary track table row <-> aod::track table global index
  /// Particle selection part

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
  Configurable<int> cfgTpcClPart2{"cfgTpcClPart2", 88, "Number of tpc clasters for the second particle"};             // min number of found TPC clusters
  Configurable<int> cfgTpcCrosRoPart2{"cfgTpcCrosRoPart2", 70, "Number of tpc crossed rows for the second particle"}; // min number of crossed rows
  Configurable<float> cfgChi2TpcPart2{"cfgChi2TpcPart2", 4.0, "Chi2 / cluster for the TPC track segment for the second particle"};
  Configurable<float> cfgChi2ItsPart2{"cfgChi2ItsPart2", 36.0, "Chi2 / cluster for the ITS track segment for the second particle"};

  /// Particle 1 K+
  Configurable<int> ConfPDGCodePartOne{"ConfPDGCodePartOne", 321, "Particle 1 - PDG code"};
  Configurable<std::vector<int>> ConfPIDPartOne{"ConfPIDPartOne", std::vector<int>{2}, "Particle 1 - Read from json file"}; // we also need the possibility to specify whether the bit is true/false ->std>>vector<std::pair<int, int>>int>>

  /// Partition for particle 1
  Partition<aod::FemtoWorldParticlesMerged> partsOne = (aod::femtoworldparticle::partType == uint8_t(aod::femtoworldparticle::ParticleType::kTrack))                         // particle type cut
                                                       && (aod::femtoworldparticle::pt < cfgPtHighPart1) && (aod::femtoworldparticle::pt > cfgPtLowPart1)                    // simple pT cuts
                                                       && (aod::femtoworldparticle::eta < cfgEtaHighPart1) && (aod::femtoworldparticle::eta > cfgEtaLowPart1)                // Eta cuts
                                                       && (o2::aod::track::dcaXY < cfgDcaXYPart1) && (o2::aod::track::dcaZ < cfgDcaZPart1)                                   // DCA cuts for XY and Z
                                                       && (aod::femtoworldparticle::tpcNClsFound > (uint8_t)cfgTpcClPart1)                                                   // Number of found TPC clusters
                                                       && (aod::femtoworldparticle::tpcNClsCrossedRows > (uint8_t)cfgTpcCrosRoPart1)                                         // Crossed rows TPC
                                                       && (aod::femtoworldparticle::itsChi2NCl < cfgChi2ItsPart1) && (aod::femtoworldparticle::tpcChi2NCl < cfgChi2TpcPart1) //&& // chi2 cuts
                                                       && (aod::femtoworldparticle::sign > int8_t(0))                                                                        // Sign (K+)
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

  /// Particle 2 K-
  Configurable<bool> ConfIsSame{"ConfIsSame", false, "Pairs of the same particle"};
  Configurable<int> ConfPDGCodePartTwo{"ConfPDGCodePartTwo", 321, "Particle 2 - PDG code"};
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
  // ConfigurableAxis CfgMultBins{"CfgMultBins", {VARIABLE_WIDTH, 0.0f, 4.0f, 8.0f, 12.0f, 16.0f, 20.0f, 24.0f, 28.0f, 32.0f, 36.0f, 40.0f, 44.0f, 48.0f, 52.0f, 56.0f, 60.0f, 64.0f, 68.0f, 72.0f, 76.0f, 80.0f, 84.0f, 88.0f, 92.0f, 96.0f, 100.0f, 200.0f, 99999.f}, "Mixing bins - multiplicity"}; // \todo to be obtained from the hash task
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

  FemtoWorldContainer<femtoWorldContainer::EventType::same, femtoWorldContainer::Observable::kstar> sameEventCont;
  FemtoWorldContainer<femtoWorldContainer::EventType::mixed, femtoWorldContainer::Observable::kstar> mixedEventCont;
  FemtoWorldPairCleaner<aod::femtoworldparticle::ParticleType::kTrack, aod::femtoworldparticle::ParticleType::kTrack> pairCleaner;
  FemtoWorldDetaDphiStar<aod::femtoworldparticle::ParticleType::kTrack, aod::femtoworldparticle::ParticleType::kTrack> pairCloseRejection;
  /// Histogram output
  HistogramRegistry qaRegistry{"TrackQA", {}, OutputObjHandlingPolicy::AnalysisObject};
  // HistogramRegistry qaRegistryFail{"TrackQAFailed", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry resultRegistry{"Correlations", {}, OutputObjHandlingPolicy::AnalysisObject};
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

    sameEventCont.init(&resultRegistry, CfgkstarBins, CfgMultBins, CfgkTBins, CfgmTBins, ConfPhiBins, ConfEtaBins, ConfMInvBins);
    sameEventCont.setPDGCodes(ConfPDGCodePartOne, ConfPDGCodePartTwo);
    mixedEventCont.init(&resultRegistry, CfgkstarBins, CfgMultBins, CfgkTBins, CfgmTBins, ConfPhiBins, ConfEtaBins, ConfMInvBins);
    mixedEventCont.setPDGCodes(ConfPDGCodePartOne, ConfPDGCodePartTwo);
    pairCleaner.init(&qaRegistry);
    if (ConfIsCPR) {
      pairCloseRejection.init(&resultRegistry, &qaRegistry, 0.01, 0.01, ConfCPRPlotPerRadii); /// \todo add config for Δη and ΔΦ cut values
    }

    vPIDPartOne = ConfPIDPartOne;
    vPIDPartTwo = ConfPIDPartTwo;
  }

  // PID
  bool IsKaonTPCdEdxNSigma(float mom, float nsigmaK) // true if accepted, false if rejected
  {
    if (mom < 0.4 && TMath::Abs(nsigmaK) < 2.0)
      return true;
    if (mom >= 0.4 && mom < 0.5 && TMath::Abs(nsigmaK) < 1.0)
      return true;
    if (mom > 0.5 && TMath::Abs(nsigmaK) < 3.0)
      return true;

    return false;
  }

  bool IsKaonTOFNSigma(float mom, float nsigmaK) // true if accepted, false if rejected
  {
    if (mom >= 0.45 && mom < 0.8 && TMath::Abs(nsigmaK) < 2.0)
      return true;
    if (mom >= 0.8 && mom < 1.0 && TMath::Abs(nsigmaK) < 1.5)
      return true;
    if (mom > 1.0 && TMath::Abs(nsigmaK) < 1.0)
      return true;
    return false;
  }

  /*
    void processSameEventPhiCand(o2::aod::FemtoWorldCollision& col,
                          o2::aod::FemtoWorldParticlesMerged& parts){

      auto groupPartsOne = partsOne->sliceByCached(aod::femtoworldparticle::femtoWorldCollisionId, col.globalIndex());
      auto groupPartsTwo = partsTwo->sliceByCached(aod::femtoworldparticle::femtoWorldCollisionId, col.globalIndex());
      auto groupPartsOneFailed = partsOneFailed->sliceByCached(aod::femtoworldparticle::femtoWorldCollisionId, col.globalIndex());
      auto groupPartsTwoFailed = partsTwoFailed->sliceByCached(aod::femtoworldparticle::femtoWorldCollisionId, col.globalIndex());

      const auto& magFieldTesla = col.magField();
        /// Now build the combinations
      for (auto& [p1, p2] : combinations(groupPartsOne, groupPartsTwo)) {
        if ((p1.p() > (float)0.45)) {
          if (!((IsKaonTPCdEdxNSigma(p1.p(), p1.tpcNSigmaKaon())) && (IsKaonTOFNSigma(p1.p(), p1.tofNSigmaKaon())))) {
            continue;
          }

        } else if ((p1.p() <= (float)0.45)) {
          if (!(IsKaonTPCdEdxNSigma(p1.p(), p1.tpcNSigmaKaon()))) {
            continue;
          }
        }
        if ((p2.p() > (float)0.45)) {
          if (!((IsKaonTPCdEdxNSigma(p2.p(), p2.tpcNSigmaKaon())) && (IsKaonTOFNSigma(p2.p(), p2.tofNSigmaKaon())))) {
            continue;
          }

        } else if ((p2.p() <= (float)0.45)) {
          if (!(IsKaonTPCdEdxNSigma(p2.p(), p2.tpcNSigmaKaon()))) {
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
        float phiPx = p1.px() + p2.px();
        float phiPy = p1.py() + p2.py();
        float phiPz = p1.pz() + p2.pz();

        TLorentzVector part1Vec;
        TLorentzVector part2Vec;
        float mMassOne = TDatabasePDG::Instance()->GetParticle(ConfPDGCodePartOne)->Mass();
        float mMassTwo = TDatabasePDG::Instance()->GetParticle(ConfPDGCodePartTwo)->Mass();
        TLorentzVector sumVec(part1Vec);
        sumVec += part2Vec;

        part1Vec.SetPtEtaPhiM(p1.pt(), p1.eta(), p1.phi(), mMassOne);
        part2Vec.SetPtEtaPhiM(p2.pt(), p2.eta(), p2.phi(), mMassTwo);

        //float phiEta = p1.eta() + p2.eta();
        //float phiPhi = p1.phi() + p2.phi();
        float phiEta = sumVec.Eta();
        float phiPhi = sumVec.Phi();

         outputPhiCan(//outputCollision.lastIndex(),
                    phiPx,
                    phiPy,
                    phiPz,
                    phiEta,
                    phiPhi);
        tmpIDtrack.push_back(p1.globalIndex());
        //Partition<aod::FemtoWorldPhiCandidates> Phiparts = (aod::femtoworldparticle::eta < cfgEtaHighPart2) && (aod::femtoworldparticle::eta > cfgEtaLowPart2);
        //sameEventCont.setPair(p1, phi multCol);
      }

    }
    PROCESS_SWITCH(femtoWorldPairTaskTrackTrack, processSameEventPhiCand, "Enable processing same event and filling phi candidates", true);*/

  /// This function processes the same event and takes care of all the histogramming
  /// \todo the trivial loops over the tracks should be factored out since they will be common to all combinations of T-T, T-V0, V0-V0, ...
  void processSameEvent(o2::aod::FemtoWorldCollision& col,
                        o2::aod::FemtoWorldParticlesMerged& parts)
  {
    // Partition<aod::FemtoWorldPhiCandidates> Phiparts = (aod::femtoworldparticle::eta < cfgEtaHighPart2) && (aod::femtoworldparticle::eta > cfgEtaLowPart2)                // Eta cuts
    //;
    // auto groupPhiParts = Phiparts->sliceByCached(aod::femtoworldparticle::femtoWorldCollisionId, col.globalIndex());
    auto groupPartsOne = partsOne->sliceByCached(aod::femtoworldparticle::femtoWorldCollisionId, col.globalIndex());
    auto groupPartsTwo = partsTwo->sliceByCached(aod::femtoworldparticle::femtoWorldCollisionId, col.globalIndex());
    auto groupPartsOneFailed = partsOneFailed->sliceByCached(aod::femtoworldparticle::femtoWorldCollisionId, col.globalIndex());
    auto groupPartsTwoFailed = partsTwoFailed->sliceByCached(aod::femtoworldparticle::femtoWorldCollisionId, col.globalIndex());

    const auto& magFieldTesla = col.magField();
    const int multCol = col.multV0M();
    eventHisto.fillQA(col);
    MixQaRegistry.fill(HIST("MixingQA/hSECollisionBins"), colBinning.getBin({col.posZ(), col.multV0M()}));
    /// Histogramming same event
    for (auto& part : groupPartsOne) {
      if ((part.p() > (float)0.45)) {
        if (!((IsKaonTPCdEdxNSigma(part.p(), part.tpcNSigmaKaon())) && (IsKaonTOFNSigma(part.p(), part.tofNSigmaKaon())))) {
          continue;
        }

      } else if ((part.p() <= (float)0.45)) {
        if (!(IsKaonTPCdEdxNSigma(part.p(), part.tpcNSigmaKaon()))) {
          continue;
        }
      }
      trackHistoPartOne.fillQA(part);
    }
    for (auto& part : groupPartsOneFailed) {

      trackHistoPartOneFailed.fillQA(part);
    }
    if (!ConfIsSame) {
      for (auto& part : groupPartsTwo) {
        if ((part.p() > (float)0.45)) {
          if (!((IsKaonTPCdEdxNSigma(part.p(), part.tpcNSigmaKaon())) && (IsKaonTOFNSigma(part.p(), part.tofNSigmaKaon())))) {
            continue;
          }

        } else if ((part.p() <= (float)0.45)) {
          if (!(IsKaonTPCdEdxNSigma(part.p(), part.tpcNSigmaKaon()))) {
            continue;
          }
        }
        trackHistoPartTwo.fillQA(part);
      }
      for (auto& part : groupPartsTwoFailed) {

        trackHistoPartTwoFailed.fillQA(part);
      }
    }
    /// Now build the combinations
    for (auto& [p1, p2] : combinations(groupPartsOne, groupPartsTwo)) {
      if ((p1.p() > (float)0.45)) {
        if (!((IsKaonTPCdEdxNSigma(p1.p(), p1.tpcNSigmaKaon())) && (IsKaonTOFNSigma(p1.p(), p1.tofNSigmaKaon())))) {
          continue;
        }

      } else if ((p1.p() <= (float)0.45)) {
        if (!(IsKaonTPCdEdxNSigma(p1.p(), p1.tpcNSigmaKaon()))) {
          continue;
        }
      }
      if ((p2.p() > (float)0.45)) {
        if (!((IsKaonTPCdEdxNSigma(p2.p(), p2.tpcNSigmaKaon())) && (IsKaonTOFNSigma(p2.p(), p2.tofNSigmaKaon())))) {
          continue;
        }

      } else if ((p2.p() <= (float)0.45)) {
        if (!(IsKaonTPCdEdxNSigma(p2.p(), p2.tpcNSigmaKaon()))) {
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
    }
  }

  PROCESS_SWITCH(femtoWorldPairTaskTrackTrack, processSameEvent, "Enable processing same event", true);

  /// This function processes the mixed event
  /// \todo the trivial loops over the collisions and tracks should be factored out since they will be common to all combinations of T-T, T-V0, V0-V0, ...
  void processMixedEvent(o2::aod::FemtoWorldCollisions& cols,
                         o2::aod::FemtoWorldParticlesMerged& parts)
  {

    for (auto& [collision1, collision2] : soa::selfCombinations(colBinning, 5, -1, cols, cols)) {

      MixQaRegistry.fill(HIST("MixingQA/hMECollisionBins"), colBinning.getBin({collision1.posZ(), collision1.multV0M()}));

      auto groupPartsOne = partsOne->sliceByCached(aod::femtoworldparticle::femtoWorldCollisionId, collision1.globalIndex());
      auto groupPartsTwo = partsTwo->sliceByCached(aod::femtoworldparticle::femtoWorldCollisionId, collision2.globalIndex());

      const auto& magFieldTesla1 = collision1.magField();
      const auto& magFieldTesla2 = collision2.magField();

      if (magFieldTesla1 != magFieldTesla2) {
        continue;
      }

      /// \todo before mixing we should check whether both collisions contain a pair of particles!
      // if (partsOne.size() == 0 || nPart2Evt1 == 0 || nPart1Evt2 == 0 || partsTwo.size() == 0 ) continue;

      for (auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(groupPartsOne, groupPartsTwo))) {
        if ((p1.p() > (float)0.45)) {
          if (!((IsKaonTPCdEdxNSigma(p1.p(), p1.tpcNSigmaKaon())) && (IsKaonTOFNSigma(p1.p(), p1.tofNSigmaKaon())))) {
            continue;
          }
        } else if ((p1.p() <= (float)0.45)) {
          if (!(IsKaonTPCdEdxNSigma(p1.p(), p1.tpcNSigmaKaon()))) {
            continue;
          }
        }
        if ((p2.p() > (float)0.45)) {
          if (!((IsKaonTPCdEdxNSigma(p2.p(), p2.tpcNSigmaKaon())) && (IsKaonTOFNSigma(p2.p(), p2.tofNSigmaKaon())))) {
            continue;
          }
        } else if ((p2.p() <= (float)0.45)) {
          if (!(IsKaonTPCdEdxNSigma(p2.p(), p2.tpcNSigmaKaon()))) {
            continue;
          }
        }

        if (ConfIsCPR) {
          if (pairCloseRejection.isClosePair(p1, p2, parts, magFieldTesla1)) {
            continue;
          }
        }
        mixedEventCont.setPair(p1, p2, collision1.multV0M());
      }
    }
  }

  PROCESS_SWITCH(femtoWorldPairTaskTrackTrack, processMixedEvent, "Enable processing mixed events", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<femtoWorldPairTaskTrackTrack>(cfgc),
  };
  return workflow;
}
