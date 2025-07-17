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

/// \file femtoUniversePairTaskV0CascadeExtended.cxx
/// \brief Task for v0-cascade correlations and QA
/// \author Shirajum Monira, WUT Warsaw, shirajum.monira@cern.ch

#include "PWGCF/FemtoUniverse/Core/FemtoUniverseContainer.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseDetaDphiStar.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseEventHisto.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniversePairCleaner.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseParticleHisto.h"
#include "PWGCF/FemtoUniverse/Core/femtoUtils.h"
#include "PWGCF/FemtoUniverse/DataModel/FemtoDerived.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/StepTHn.h"
#include "Framework/runDataProcessing.h"

#include <set>
#include <vector>

using namespace o2;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::analysis::femto_universe;
using namespace o2::aod::pidutils;

struct FemtoUniversePairTaskV0CascadeExtended {

  SliceCache cache;
  using FemtoParticles = soa::Join<aod::FDParticles, aod::FDExtParticles>;
  Preslice<FemtoParticles> perCol = aod::femtouniverseparticle::fdCollisionId;

  /// applying narrow cut
  Configurable<float> confZVertexCut{"confZVertexCut", 10.f, "Event sel: Maximum z-Vertex (cm)"};

  Filter collisionFilter = (nabs(aod::collision::posZ) < confZVertexCut);
  using FilteredFDCollisions = soa::Filtered<o2::aod::FdCollisions>;
  using FilteredFDCollision = FilteredFDCollisions::iterator;

  /// particle 1 (v0)
  Configurable<float> confHPtPart1{"confHPtPart1", 4.0f, "higher limit for pt of particle 1"};
  Configurable<float> confLPtPart1{"confLPtPart1", 0.5f, "lower limit for pt of particle 1"};
  Configurable<int> confV0Type{"confV0Type", 0, "select one of the V0s (lambda = 0, anti-lambda = 1, k0 = 2)"};
  Configurable<float> confV0InvMassLowLimit{"confV0InvMassLowLimit", 1.10, "Lower limit of the V0 invariant mass"};
  Configurable<float> confV0InvMassUpLimit{"confV0InvMassUpLimit", 1.13, "Upper limit of the V0 invariant mass"};
  Configurable<int> confV0PDGCode{"confV0PDGCode", 3122, "Particle 1 (V0) - PDG code"};

  /// particle 2 (cascade)
  Configurable<float> confHPtPart2{"confHPtPart2", 4.0f, "higher limit for pt of particle 2"};
  Configurable<float> confLPtPart2{"confLPtPart2", 0.3f, "lower limit for pt of particle 2"};
  Configurable<int> confCascType{"confCascType", 0, "select one of the cascades (Omega = 0, Xi = 1, anti-Omega = 2, anti-Xi = 3)"};
  Configurable<float> confCascInvMassLowLimit{"confCascInvMassLowLimit", 1.315, "Lower limit of the cascade invariant mass"};
  Configurable<float> confCascInvMassUpLimit{"confCascInvMassUpLimit", 1.325, "Upper limit of the cascade invariant mass"};

  /// nSigma cuts
  Configurable<float> confmom{"confmom", 0.75, "momentum threshold for particle identification using TOF"};
  Configurable<float> confNsigmaTPCParticleChild{"confNsigmaTPCParticleChild", 3.0, "TPC Sigma for cascade (daugh & bach) momentum < Confmom"};
  Configurable<float> confNsigmaTOFParticleChild{"confNsigmaTOFParticleChild", 3.0, "TOF Sigma for cascade (daugh & bach) momentum > Confmom"};

  /// for correlation part
  Configurable<bool> confIsMC{"confIsMC", false, "Enable additional Histograms in the case of a MonteCarlo Run"};
  ConfigurableAxis confkstarBins{"confkstarBins", {1500, 0., 6.}, "binning kstar"};
  ConfigurableAxis confMultBins{"confMultBins", {VARIABLE_WIDTH, 0.0f, 20.0f, 40.0f, 60.0f, 80.0f, 100.0f, 200.0f, 99999.f}, "Mixing bins - multiplicity"};
  ConfigurableAxis confkTBins{"confkTBins", {150, 0., 9.}, "binning kT"};
  ConfigurableAxis confmTBins{"confmTBins", {225, 0., 7.5}, "binning mT"};
  ConfigurableAxis confMultBins3D{"confMultBins3D", {VARIABLE_WIDTH, 0.0f, 20.0f, 30.0f, 40.0f, 99999.0f}, "multiplicity Binning for the 3Dimensional plot: k* vs multiplicity vs mT (set <<confUse3D>> to true in order to use)"};
  ConfigurableAxis confmTBins3D{"confmTBins3D", {VARIABLE_WIDTH, 1.02f, 1.14f, 1.20f, 1.26f, 1.38f, 1.56f, 1.86f, 4.50f}, "mT Binning for the 3Dimensional plot: k* vs multiplicity vs mT (set <<confUse3D>> to true in order to use)"};
  Configurable<int> confEtaBins{"confEtaBins", 29, "Number of eta bins in deta dphi"};
  Configurable<int> confPhiBins{"confPhiBins", 29, "Number of phi bins in deta dphi"};
  Configurable<bool> confUse3D{"confUse3D", false, "Enable three dimensional histogramms (to be used only for analysis with high statistics): k* vs mT vs multiplicity"};
  ConfigurableAxis confVtxBins{"confVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis confTrkTempFitVarpTBins{"confTrkTempFitVarpTBins", {20, 0.5, 4.05}, "pT binning of the pT vs. TempFitVar plot"};
  ConfigurableAxis confTrkTempFitVarBins{"confTrkTempFitVarBins", {300, -0.15, 0.15}, "binning of the TempFitVar in the pT vs. TempFitVar plot"};
  ConfigurableAxis confChildTempFitVarBins{"confChildTempFitVarBins", {300, -0.15, 0.15}, "V0 child: binning of the TempFitVar in the pT vs. TempFitVar plot"};
  ConfigurableAxis confChildTempFitVarpTBins{"confChildTempFitVarpTBins", {20, 0.5, 4.05}, "V0 child: pT binning of the pT vs. TempFitVar plot"};
  ConfigurableAxis confV0TempFitVarpTBins{"confV0TempFitVarpTBins", {20, 0.5, 4.05}, "V0: pT binning of the pT vs. TempFitVar plot"};
  ConfigurableAxis confV0TempFitVarBins{"confV0TempFitVarBins", {300, 0.95, 1.}, "V0: binning of the TempFitVar in the pT vs. TempFitVar plot"};

  /// Partition for particle 1
  Partition<FemtoParticles> partsOne = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kV0)) && (aod::femtouniverseparticle::pt < confHPtPart1) && (aod::femtouniverseparticle::pt > confLPtPart1);

  /// Partition for particle 2
  Partition<FemtoParticles> partsTwo = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kCascade)) && (aod::femtouniverseparticle::pt < confHPtPart2) && (aod::femtouniverseparticle::pt > confLPtPart2);

  /// Histogramming for v0
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kV0, 2> trackHistoV0;
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kV0Child, 3> posChildV0;
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kV0Child, 4> negChildV0;

  /// Histogramming for cascade
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kV0Child, 3> posChildHistosCasc;
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kV0Child, 4> negChildHistosCasc;
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kCascadeBachelor, 0> bachHistosCasc;
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kCascade, 0> cascQAHistos;

  /// Histogramming for Event
  FemtoUniverseEventHisto eventHisto;

  FemtoUniverseContainer<femto_universe_container::EventType::same, femto_universe_container::Observable::kstar> sameEventCont;
  FemtoUniverseContainer<femto_universe_container::EventType::mixed, femto_universe_container::Observable::kstar> mixedEventCont;
  FemtoUniversePairCleaner<aod::femtouniverseparticle::ParticleType::kV0, aod::femtouniverseparticle::ParticleType::kCascade> pairCleaner;

  HistogramRegistry qaRegistry{"TrackQA", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry resultRegistry{"Correlations", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Table to select v0 daughters
  static constexpr unsigned int V0ChildTable[][2] = {{0, 1}, {1, 0}, {1, 1}};

  bool invMLambda(float invMassLambda, float invMassAntiLambda)
  {
    if ((invMassLambda < confV0InvMassLowLimit || invMassLambda > confV0InvMassUpLimit) && (invMassAntiLambda < confV0InvMassLowLimit || invMassAntiLambda > confV0InvMassUpLimit)) {
      return false;
    }
    return true;
  }

  // Table to select cascade daughters
  // Charges: = +--, +--, +-+, +-+
  static constexpr unsigned int CascChildTable[][3] = {{0, 1, 2}, {0, 1, 1}, {1, 0, 2}, {1, 0, 1}};

  bool invMCascade(float invMassXi, float invMassOmega, int cascType)
  {
    return (((cascType == 1 || cascType == 3) && (invMassXi > confCascInvMassLowLimit && invMassXi < confCascInvMassUpLimit)) || ((cascType == 0 || cascType == 2) && (invMassOmega > confCascInvMassLowLimit && invMassOmega < confCascInvMassUpLimit)));
  }

  bool isNSigmaTPC(float nsigmaTPCParticle)
  {
    if (std::abs(nsigmaTPCParticle) < confNsigmaTPCParticleChild) {
      return true;
    } else {
      return false;
    }
  }

  bool isNSigmaTOF(float mom, float nsigmaTOFParticle, float hasTOF)
  {
    // Cut only on daughter and bachelor tracks, that have TOF signal
    if (mom > confmom && hasTOF == 1) {
      if (std::abs(nsigmaTOFParticle) < confNsigmaTOFParticleChild) {
        return true;
      } else {
        return false;
      }
    } else {
      return true;
    }
  }

  template <typename T>
  bool isParticleTPC(const T& part, int id)
  {
    const float tpcNSigmas[3] = {unPackInTable<aod::pidtpc_tiny::binning>(part.tpcNSigmaStorePr()), unPackInTable<aod::pidtpc_tiny::binning>(part.tpcNSigmaStorePi()), unPackInTable<aod::pidtpc_tiny::binning>(part.tpcNSigmaStoreKa())};

    return isNSigmaTPC(tpcNSigmas[id]);
  }

  template <typename T>
  bool isParticleTOF(const T& part, int id)
  {
    const float tofNSigmas[3] = {unPackInTable<aod::pidtof_tiny::binning>(part.tofNSigmaStorePr()), unPackInTable<aod::pidtof_tiny::binning>(part.tofNSigmaStorePi()), unPackInTable<aod::pidtof_tiny::binning>(part.tofNSigmaStoreKa())};

    return isNSigmaTOF(part.p(), tofNSigmas[id], part.tempFitVar());
  }

  void init(InitContext const&)
  {
    trackHistoV0.init(&qaRegistry, confV0TempFitVarpTBins, confV0TempFitVarBins, confIsMC, confV0PDGCode, true, "trackHistoV0");
    posChildV0.init(&qaRegistry, confChildTempFitVarpTBins, confChildTempFitVarBins, false, 0, true, "posChildV0");
    negChildV0.init(&qaRegistry, confChildTempFitVarpTBins, confChildTempFitVarBins, false, 0, true, "negChildV0");

    posChildHistosCasc.init(&qaRegistry, confChildTempFitVarpTBins, confChildTempFitVarBins, false, 0, true, "posChildCasc");
    negChildHistosCasc.init(&qaRegistry, confChildTempFitVarpTBins, confChildTempFitVarBins, false, 0, true, "negChildCasc");
    bachHistosCasc.init(&qaRegistry, confChildTempFitVarpTBins, confChildTempFitVarBins, false, 0, true, "hBachelor");
    cascQAHistos.init(&qaRegistry, confChildTempFitVarpTBins, confChildTempFitVarBins, false, 0, true);

    sameEventCont.init(&resultRegistry, confkstarBins, confMultBins, confkTBins, confmTBins, confMultBins3D, confmTBins3D, confEtaBins, confPhiBins, confIsMC, confUse3D);
    mixedEventCont.init(&resultRegistry, confkstarBins, confMultBins, confkTBins, confmTBins, confMultBins3D, confmTBins3D, confEtaBins, confPhiBins, confIsMC, confUse3D);
    pairCleaner.init(&qaRegistry);
  }

  /// v0-cascade correlations same event
  void processSameEvent(const FilteredFDCollision& col, const FemtoParticles& parts)
  {
    auto groupPartsOne = partsOne->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto groupPartsTwo = partsTwo->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);

    eventHisto.fillQA(col);

    for (const auto& part : groupPartsTwo) {
      /// inv Mass check for cascade
      if (!invMCascade(part.mLambda(), part.mAntiLambda(), confCascType))
        continue;

      cascQAHistos.fillQA<false, true>(part);

      const auto& posChild1 = parts.iteratorAt(part.index() - 3);
      const auto& negChild1 = parts.iteratorAt(part.index() - 2);
      const auto& bachelor1 = parts.iteratorAt(part.index() - 1);
      /// Children of cascade must pass this condition to be selected
      if (!isParticleTPC(posChild1, CascChildTable[confCascType][0]) || !isParticleTPC(negChild1, CascChildTable[confCascType][1]) || !isParticleTPC(bachelor1, CascChildTable[confCascType][2]))
        continue;

      if (!isParticleTOF(posChild1, CascChildTable[confCascType][0]) || !isParticleTOF(negChild1, CascChildTable[confCascType][1]) || !isParticleTOF(bachelor1, CascChildTable[confCascType][2]))
        continue;

      posChildHistosCasc.fillQABase<false, true>(posChild1, HIST("posChildCasc"));
      negChildHistosCasc.fillQABase<false, true>(negChild1, HIST("negChildCasc"));
      bachHistosCasc.fillQABase<false, true>(bachelor1, HIST("hBachelor"));
    }
    for (const auto& part : groupPartsOne) {
      /// inv Mass check for V0s
      if (!invMLambda(part.mLambda(), part.mAntiLambda()))
        continue;
      const auto& posChild2 = parts.iteratorAt(part.index() - 2);
      const auto& negChild2 = parts.iteratorAt(part.index() - 1);
      /// Daughters of v0 must pass this condition to be selected
      if (!isParticleTPC(posChild2, V0ChildTable[confV0Type][0]) || !isParticleTPC(negChild2, V0ChildTable[confV0Type][1]))
        continue;

      trackHistoV0.fillQABase<false, true>(part, HIST("trackHistoV0"));
      posChildV0.fillQABase<false, true>(posChild2, HIST("posChildV0"));
      negChildV0.fillQABase<false, true>(negChild2, HIST("negChildV0"));
    }
    for (const auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(groupPartsOne, groupPartsTwo))) {
      if (!invMLambda(p1.mLambda(), p1.mAntiLambda()))
        continue;
      // Cascase inv Mass cut (mLambda stores Xi mass, mAntiLambda stored Omega mass)
      if (!invMCascade(p2.mLambda(), p2.mAntiLambda(), confCascType))
        continue;
      // track cleaning
      if (!pairCleaner.isCleanPair(p1, p2, parts)) {
        continue;
      }
      // V0
      const auto& posChild2 = parts.iteratorAt(p1.index() - 2);
      const auto& negChild2 = parts.iteratorAt(p1.index() - 1);
      /// Daughters of v0 must pass this condition to be selected
      if (!isParticleTPC(posChild2, V0ChildTable[confV0Type][0]) || !isParticleTPC(negChild2, V0ChildTable[confV0Type][1]))
        continue;
      // cascade
      const auto& posChild1 = parts.iteratorAt(p2.index() - 3);
      const auto& negChild1 = parts.iteratorAt(p2.index() - 2);
      const auto& bachelor1 = parts.iteratorAt(p2.index() - 1);
      /// Daughters of cascade must pass this condition to be selected
      if (!isParticleTPC(posChild1, CascChildTable[confCascType][0]) || !isParticleTPC(negChild1, CascChildTable[confCascType][1]) || !isParticleTPC(bachelor1, CascChildTable[confCascType][2]))
        continue;

      if (!isParticleTOF(posChild1, CascChildTable[confCascType][0]) || !isParticleTOF(negChild1, CascChildTable[confCascType][1]) || !isParticleTOF(bachelor1, CascChildTable[confCascType][2]))
        continue;

      sameEventCont.setPair<false>(p1, p2, col.multNtr(), confUse3D, 1.0f);
    }
  }
  PROCESS_SWITCH(FemtoUniversePairTaskV0CascadeExtended, processSameEvent, "Enable processing same event for v0 - cascade", false);

  /// v0-cascade correlations mixed event
  void processMixedEvent(const FilteredFDCollisions& cols, const FemtoParticles& parts)
  {
    ColumnBinningPolicy<aod::collision::PosZ, aod::femtouniversecollision::MultNtr> colBinning{{confVtxBins, confMultBins}, true};

    for (const auto& [collision1, collision2] : soa::selfCombinations(colBinning, 5, -1, cols, cols)) {
      const int multCol = collision1.multNtr();

      auto groupPartsOne = partsOne->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);
      auto groupPartsTwo = partsTwo->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision2.globalIndex(), cache);

      const auto& magFieldTesla1 = collision1.magField();
      const auto& magFieldTesla2 = collision2.magField();

      if (magFieldTesla1 != magFieldTesla2) {
        continue;
      }
      for (const auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(groupPartsOne, groupPartsTwo))) {
        if (!invMLambda(p1.mLambda(), p1.mAntiLambda()))
          continue;
        // Cascase inv Mass cut (mLambda stores Xi mass, mAntiLambda stored Omega mass)
        if (!invMCascade(p2.mLambda(), p2.mAntiLambda(), confCascType))
          continue;
        // V0
        const auto& posChild2 = parts.iteratorAt(p1.index() - 2);
        const auto& negChild2 = parts.iteratorAt(p1.index() - 1);
        /// Daughters of v0 must pass this condition to be selected
        if (!isParticleTPC(posChild2, V0ChildTable[confV0Type][0]) || !isParticleTPC(negChild2, V0ChildTable[confV0Type][1]))
          continue;
        // cascade
        const auto& posChild1 = parts.iteratorAt(p2.index() - 3);
        const auto& negChild1 = parts.iteratorAt(p2.index() - 2);
        const auto& bachelor1 = parts.iteratorAt(p2.index() - 1);
        /// Daughters of cascade must pass this condition to be selected
        if (!isParticleTPC(posChild1, CascChildTable[confCascType][0]) || !isParticleTPC(negChild1, CascChildTable[confCascType][1]) || !isParticleTPC(bachelor1, CascChildTable[confCascType][2]))
          continue;

        if (!isParticleTOF(posChild1, CascChildTable[confCascType][0]) || !isParticleTOF(negChild1, CascChildTable[confCascType][1]) || !isParticleTOF(bachelor1, CascChildTable[confCascType][2]))
          continue;

        // track cleaning
        if (!pairCleaner.isCleanPair(p1, p2, parts)) {
          continue;
        }

        mixedEventCont.setPair<false>(p1, p2, multCol, confUse3D, 1.0f);
      }
    }
  }
  PROCESS_SWITCH(FemtoUniversePairTaskV0CascadeExtended, processMixedEvent, "Enable processing mixed event for v0 - cascade", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<FemtoUniversePairTaskV0CascadeExtended>(cfgc),
  };
  return workflow;
}
