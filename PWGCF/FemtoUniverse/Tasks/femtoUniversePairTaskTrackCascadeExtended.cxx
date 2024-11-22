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

/// \brief Task for cascade QA; in the future: for cascade correlations
/// \author Barbara Chytla, WUT Warsaw, barbara.chytla@cern.ch
/// \author Shirajum Monira, WUT Warsaw, shirajum.monira@cern.ch

#include <vector>
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/StepTHn.h"
#include "PWGCF/FemtoUniverse/DataModel/FemtoDerived.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseParticleHisto.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseEventHisto.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniversePairCleaner.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseContainer.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseDetaDphiStar.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUtils.h"

using namespace o2;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::analysis::femtoUniverse;
using namespace o2::aod::pidutils;

struct femtoUniversePairTaskTrackCascadeExtended {

  SliceCache cache;
  using FemtoFullParticles = soa::Join<aod::FDCascParticles, aod::FDExtParticles>;
  Preslice<FemtoFullParticles> perCol = aod::femtouniverseparticle::fdCollisionId;

  Configurable<float> ConfZVertexCut{"ConfZVertexCut", 10.f, "Event sel: Maximum z-Vertex (cm)"};

  Filter collisionFilter = (nabs(aod::collision::posZ) < ConfZVertexCut);
  using FilteredFDCollisions = soa::Filtered<o2::aod::FDCollisions>;
  using FilteredFDCollision = FilteredFDCollisions::iterator;

  ConfigurableAxis ConfChildTempFitVarpTBins{"ConfChildTempFitVarpTBins", {20, 0.5, 4.05}, "V0 child: pT binning of the pT vs. TempFitVar plot"};
  ConfigurableAxis ConfChildTempFitVarBins{"ConfChildTempFitVarBins", {300, -0.15, 0.15}, "V0 child: binning of the TempFitVar in the pT vs. TempFitVar plot"};
  Configurable<float> ConfCascInvMassLowLimit{"ConfCascInvMassLowLimit", 1.315, "Lower limit of the Casc invariant mass"};
  Configurable<float> ConfCascInvMassUpLimit{"ConfCascInvMassUpLimit", 1.325, "Upper limit of the Casc invariant mass"};
  Configurable<float> ConfCascTranRad{"ConfCascTranRad", 0.5, "Cascade transverse radius"};

  Configurable<float> NSigmaTPCPion{"NSigmaTPCPion", 4, "NSigmaTPCPion"};
  Configurable<float> NSigmaTPCProton{"NSigmaTPCProton", 4, "NSigmaTPCProton"};

  // configs for correlation part
  Configurable<int> ConfTrackChoicePartOne{"ConfTrackChoicePartOne", 0, "0:Proton, 1:Pion, 2:Kaon"};
  Configurable<int> ConfTrkPDGCodePartOne{"ConfTrkPDGCodePartOne", 2212, "Particle 1 (Track) - PDG code"};
  Configurable<int> ConfCascType1{"ConfCascType1", 0, "select one of the V0s (Omega = 0, Xi = 1, anti-Omega = 2, anti-Xi = 3) for track-cascade combination"};
  Configurable<int> ConfChargePart1{"ConfChargePart1", 0, "sign of particle 1"};
  Configurable<float> ConfHPtPart1{"ConfHPtPart1", 4.0f, "higher limit for pt of particle 1"};
  Configurable<float> ConfLPtPart1{"ConfLPtPart1", 0.5f, "lower limit for pt of particle 1"};
  Configurable<float> ConfHPtPart2{"ConfHPtPart2", 4.0f, "higher limit for pt of particle 2"};
  Configurable<float> ConfLPtPart2{"ConfLPtPart2", 0.3f, "lower limit for pt of particle 2"};
  Configurable<float> Confmom{"Confmom", 0.75, "momentum threshold for particle identification using TOF"};
  Configurable<float> ConfNsigmaTPCParticle{"ConfNsigmaTPCParticle", 3.0, "TPC Sigma for particle momentum < Confmom"};
  Configurable<float> ConfNsigmaCombinedParticle{"ConfNsigmaCombinedParticle", 3.0, "TPC and TOF Sigma (combined) for particle momentum > Confmom"};

  ConfigurableAxis ConfkstarBins{"ConfkstarBins", {1500, 0., 6.}, "binning kstar"};
  ConfigurableAxis ConfMultBins{"ConfMultBins", {VARIABLE_WIDTH, 0.0f, 20.0f, 40.0f, 60.0f, 80.0f, 100.0f, 200.0f, 99999.f}, "Mixing bins - multiplicity"};
  ConfigurableAxis ConfkTBins{"ConfkTBins", {150, 0., 9.}, "binning kT"};
  ConfigurableAxis ConfmTBins{"ConfmTBins", {225, 0., 7.5}, "binning mT"};
  ConfigurableAxis ConfmultBins3D{"ConfMultBins3D", {VARIABLE_WIDTH, 0.0f, 20.0f, 30.0f, 40.0f, 99999.0f}, "multiplicity Binning for the 3Dimensional plot: k* vs multiplicity vs mT (set <<ConfUse3D>> to true in order to use)"};
  ConfigurableAxis ConfmTBins3D{"ConfmTBins3D", {VARIABLE_WIDTH, 1.02f, 1.14f, 1.20f, 1.26f, 1.38f, 1.56f, 1.86f, 4.50f}, "mT Binning for the 3Dimensional plot: k* vs multiplicity vs mT (set <<ConfUse3D>> to true in order to use)"};
  Configurable<int> ConfEtaBins{"ConfEtaBins", 29, "Number of eta bins in deta dphi"};
  Configurable<int> ConfPhiBins{"ConfPhiBins", 29, "Number of phi bins in deta dphi"};
  Configurable<bool> ConfIsMC{"ConfIsMC", false, "Enable additional Histograms in the case of a MonteCarlo Run"};
  Configurable<bool> ConfUse3D{"ConfUse3D", false, "Enable three dimensional histogramms (to be used only for analysis with high statistics): k* vs mT vs multiplicity"};
  ConfigurableAxis ConfVtxBins{"ConfVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis ConfTrkTempFitVarpTBins{"ConfTrkTempFitVarpTBins", {20, 0.5, 4.05}, "pT binning of the pT vs. TempFitVar plot"};
  ConfigurableAxis ConfTrkTempFitVarBins{"ConfTrkDTempFitVarBins", {300, -0.15, 0.15}, "binning of the TempFitVar in the pT vs. TempFitVar plot"};

  /// Partition for particle 1 (track)
  Partition<FemtoFullParticles> partsOne = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kTrack)) && (aod::femtouniverseparticle::sign == ConfChargePart1) && (aod::femtouniverseparticle::pt < ConfHPtPart1) && (aod::femtouniverseparticle::pt > ConfLPtPart1);

  /// Partition for particle 2 (cascade)
  Partition<FemtoFullParticles> partsTwo = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kCascade)) && (aod::femtouniverseparticle::pt < ConfHPtPart2) && (aod::femtouniverseparticle::pt > ConfLPtPart2);

  /// Partition for cascades
  Partition<FemtoFullParticles> cascs = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kCascade));

  /// Histogramming for track particle
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kTrack, 3> trackHistoPartOnePos;
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kTrack, 4> trackHistoPartOneNeg;

  /// Histogramming for cascade
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kV0Child, 3> posChildHistos;
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kV0Child, 4> negChildHistos;
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kCascadeBachelor, 2> bachHistos;

  FemtoUniverseContainer<femtoUniverseContainer::EventType::same, femtoUniverseContainer::Observable::kstar> sameEventCont;
  FemtoUniverseContainer<femtoUniverseContainer::EventType::mixed, femtoUniverseContainer::Observable::kstar> mixedEventCont;

  HistogramRegistry rXiQA{"xi", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry qaRegistry{"TrackQA", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry resultRegistry{"Correlations", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Table to select cascade daughters
  // Charges: = +--, +--, +-+, +-+
  static constexpr UInt_t CascChildTable[][3] = {{0, 1, 2}, {0, 1, 1}, {1, 0, 2}, {1, 0, 1}};

  bool invMCascade(float invMassCascade, float invMassAntiCascade)
  {
    if ((invMassCascade < ConfCascInvMassLowLimit || invMassCascade > ConfCascInvMassUpLimit) && (invMassAntiCascade < ConfCascInvMassLowLimit || invMassAntiCascade > ConfCascInvMassUpLimit)) {
      return false;
    }
    return true;
  }

  bool IsNSigmaTPC(float nsigmaTPCParticle)
  {
    if (TMath::Abs(nsigmaTPCParticle) < ConfNsigmaTPCParticle) {
      return true;
    } else {
      return false;
    }
  }

  bool IsNSigmaCombined(float mom, float nsigmaTPCParticle, float nsigmaTOFParticle)
  {
    if (mom <= Confmom) {
      return (TMath::Abs(nsigmaTPCParticle) < ConfNsigmaTPCParticle);
    } else {
      return (TMath::Hypot(nsigmaTOFParticle, nsigmaTPCParticle) < ConfNsigmaCombinedParticle);
    }
  }

  template <typename T>
  bool IsParticleTPC(const T& part, int id)
  {
    const float tpcNSigmas[3] = {unPackInTable<aod::pidtpc_tiny::binning>(part.tpcNSigmaStorePr()), unPackInTable<aod::pidtpc_tiny::binning>(part.tpcNSigmaStorePi()), unPackInTable<aod::pidtpc_tiny::binning>(part.tpcNSigmaStoreKa())};

    return IsNSigmaTPC(tpcNSigmas[id]);
  }

  template <typename T>
  bool IsParticleCombined(const T& part, int id)
  {
    const float tpcNSigmas[3] = {unPackInTable<aod::pidtpc_tiny::binning>(part.tpcNSigmaStorePr()), unPackInTable<aod::pidtpc_tiny::binning>(part.tpcNSigmaStorePi()), unPackInTable<aod::pidtpc_tiny::binning>(part.tpcNSigmaStoreKa())};
    const float tofNSigmas[3] = {unPackInTable<aod::pidtof_tiny::binning>(part.tofNSigmaStorePr()), unPackInTable<aod::pidtof_tiny::binning>(part.tofNSigmaStorePi()), unPackInTable<aod::pidtof_tiny::binning>(part.tofNSigmaStoreKa())};

    return IsNSigmaCombined(part.p(), tpcNSigmas[id], tofNSigmas[id]);
  }

  void init(InitContext const&)
  {
    // Axes
    AxisSpec XiMassAxis = {200, 1.28f, 1.36f, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec ptAxis = {100, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec etaAxis = {100, -2.0f, 2.0f, "#it{#eta}"};
    AxisSpec phiAxis = {100, 0.0f, 6.0f, "#it{#phi}"};
    AxisSpec DCADaughAxis = {1000, 0.0f, 2.0f, "DCA (cm)"};
    AxisSpec CPAAxis = {1000, 0.95f, 1.0f, "#it{cos #theta_{p}}"};
    AxisSpec tranRadAxis = {1000, 0.0f, 100.0f, "#it{r}_{xy} (cm)"};
    AxisSpec DCAToPVAxis = {1000, -10.0f, 10.0f, "DCA to PV (cm)"};

    // Histograms
    rXiQA.add("hMassXi", "hMassXi", {HistType::kTH1F, {XiMassAxis}});
    rXiQA.add("hMassXiSelected", "hMassXiSelected", {HistType::kTH1F, {XiMassAxis}});
    rXiQA.add("hPtXi", "hPtXi", {HistType::kTH1F, {{ptAxis}}});
    rXiQA.add("hEtaXi", "hEtaXi", {HistType::kTH1F, {{etaAxis}}});
    rXiQA.add("hPhiXi", "hPhiXi", {HistType::kTH1F, {{phiAxis}}});
    rXiQA.add("hDCAV0Daughters", "hDCAV0Daughters", {HistType::kTH1F, {DCADaughAxis}});
    rXiQA.add("hV0CosPA", "hV0CosPA", {HistType::kTH1F, {CPAAxis}});
    rXiQA.add("hV0TranRad", "hV0TranRad", {HistType::kTH1F, {tranRadAxis}});
    rXiQA.add("hDCACascDaughters", "hDCACascDaughters", {HistType::kTH1F, {DCADaughAxis}});
    rXiQA.add("hCascCosPA", "hCascCosPA", {HistType::kTH1F, {CPAAxis}});
    rXiQA.add("hCascTranRad", "hCascTranRad", {HistType::kTH1F, {tranRadAxis}});
    rXiQA.add("hDcaPostoPV", "hDcaPostoPV", {HistType::kTH1F, {DCAToPVAxis}});
    rXiQA.add("hDcaNegtoPV", "hDcaNegtoPV", {HistType::kTH1F, {DCAToPVAxis}});
    rXiQA.add("hDcaBachtoPV", "hDcaBachtoPV", {HistType::kTH1F, {DCAToPVAxis}});
    rXiQA.add("hDcaV0toPV", "hDcaV0toPV", {HistType::kTH1F, {DCAToPVAxis}});

    qaRegistry.add("Tracks_pos/nSigmaTPC", "; #it{p} (GeV/#it{c}); n#sigma_{TPC}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("Tracks_pos/nSigmaTOF", "; #it{p} (GeV/#it{c}); n#sigma_{TOF}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("Tracks_neg/nSigmaTPC", "; #it{p} (GeV/#it{c}); n#sigma_{TPC}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("Tracks_neg/nSigmaTOF", "; #it{p} (GeV/#it{c}); n#sigma_{TOF}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});

    trackHistoPartOnePos.init(&qaRegistry, ConfTrkTempFitVarpTBins, ConfTrkTempFitVarBins, ConfIsMC, ConfTrkPDGCodePartOne);
    trackHistoPartOneNeg.init(&qaRegistry, ConfTrkTempFitVarpTBins, ConfTrkTempFitVarBins, ConfIsMC, ConfTrkPDGCodePartOne);
    posChildHistos.init(&qaRegistry, ConfChildTempFitVarpTBins, ConfChildTempFitVarBins, false, 0, true);
    negChildHistos.init(&qaRegistry, ConfChildTempFitVarpTBins, ConfChildTempFitVarBins, false, 0, true);
    bachHistos.init(&qaRegistry, ConfChildTempFitVarpTBins, ConfChildTempFitVarBins, false, 0, true, "hBachelor");

    sameEventCont.init(&resultRegistry, ConfkstarBins, ConfMultBins, ConfkTBins, ConfmTBins, ConfmultBins3D, ConfmTBins3D, ConfEtaBins, ConfPhiBins, ConfIsMC, ConfUse3D);
    mixedEventCont.init(&resultRegistry, ConfkstarBins, ConfMultBins, ConfkTBins, ConfmTBins, ConfmultBins3D, ConfmTBins3D, ConfEtaBins, ConfPhiBins, ConfIsMC, ConfUse3D);
  }

  void processCascades(FilteredFDCollision& col, FemtoFullParticles& parts)
  {
    auto groupCascs = cascs->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    // const int multCol = col.multNtr();

    for (auto& casc : groupCascs) {
      rXiQA.fill(HIST("hMassXi"), casc.mLambda());

      // if (!invMCascade(casc.mLambda(), casc.mAntiLambda()))
      //   continue;

      const auto& posChild = parts.iteratorAt(casc.index() - 3);
      const auto& negChild = parts.iteratorAt(casc.index() - 2);
      const auto& bachelor = parts.iteratorAt(casc.index() - 1);

      // if (casc.transRadius() < ConfCascTranRad)
      //   continue;
      // std::cout<<std::endl;
      // std::cout<<"TYPE:"<<std::endl;
      // std::cout<<casc.partType()<<std::endl;
      //  nSigma selection for daughter and bachelor tracks
      if (casc.sign() < 0) {
        if (TMath::Abs(posChild.tpcNSigmaPr()) > NSigmaTPCProton) {
          continue;
        }
        if (TMath::Abs(negChild.tpcNSigmaPi()) > NSigmaTPCPion) {
          continue;
        }
      } else {
        if (TMath::Abs(negChild.tpcNSigmaPr()) > NSigmaTPCProton) {
          continue;
        }
        if (TMath::Abs(posChild.tpcNSigmaPi()) > NSigmaTPCPion) {
          continue;
        }
      }
      if (TMath::Abs(bachelor.tpcNSigmaPi()) > NSigmaTPCPion) {
        continue;
      }

      rXiQA.fill(HIST("hPtXi"), casc.pt());
      rXiQA.fill(HIST("hEtaXi"), casc.eta());
      rXiQA.fill(HIST("hPhiXi"), casc.phi());
      rXiQA.fill(HIST("hMassXiSelected"), casc.mLambda());
      rXiQA.fill(HIST("hDCAV0Daughters"), casc.dcaV0daughters());
      rXiQA.fill(HIST("hV0CosPA"), casc.cpav0());
      rXiQA.fill(HIST("hV0TranRad"), casc.v0radius());
      rXiQA.fill(HIST("hCascCosPA"), casc.cpaCasc());
      rXiQA.fill(HIST("hDCACascDaughters"), casc.dcacascdaughters());
      rXiQA.fill(HIST("hCascTranRad"), casc.cascradius());
      rXiQA.fill(HIST("hDcaPostoPV"), casc.dcapostopv());
      rXiQA.fill(HIST("hDcaNegtoPV"), casc.dcanegtopv());
      rXiQA.fill(HIST("hDcaBachtoPV"), casc.dcabachtopv());
      rXiQA.fill(HIST("hDcaV0toPV"), casc.dcav0topv());

      posChildHistos.fillQA<false, true>(posChild);
      negChildHistos.fillQA<false, true>(negChild);
      bachHistos.fillQABase<false, true>(bachelor, HIST("hBachelor"));
    }
  }
  PROCESS_SWITCH(femtoUniversePairTaskTrackCascadeExtended, processCascades, "Enable processing cascades", true);

  void processSameEvent(FilteredFDCollision& col, FemtoFullParticles& parts)
  {
    // const auto& magFieldTesla = col.magField();

    auto groupPartsOne = partsOne->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto groupPartsTwo = partsTwo->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);

    for (auto& part : groupPartsTwo) {
      if (!invMCascade(part.mLambda(), part.mAntiLambda()))
        continue;
      const auto& posChild = parts.iteratorAt(part.index() - 3);
      const auto& negChild = parts.iteratorAt(part.index() - 2);
      const auto& bachelor = parts.iteratorAt(part.index() - 1);
      /// Child particles must pass this condition to be selected
      /*if (!IsParticleTPC(posChild, CascChildTable[ConfCascType1][0]) || !IsParticleTPC(negChild, CascChildTable[ConfCascType1][1]) || !IsParticleTPC(bachelor, CascChildTable[ConfCascType1][2]))
        continue;*/

      posChildHistos.fillQA<false, true>(posChild);
      negChildHistos.fillQA<false, true>(negChild);
      bachHistos.fillQABase<false, true>(bachelor, HIST("hBachelor"));
    }

    for (auto& part : groupPartsOne) {
      /// PID plot for track particle
      const float tpcNSigmas[3] = {unPackInTable<aod::pidtpc_tiny::binning>(part.tpcNSigmaStorePr()), unPackInTable<aod::pidtpc_tiny::binning>(part.tpcNSigmaStorePi()), unPackInTable<aod::pidtpc_tiny::binning>(part.tpcNSigmaStoreKa())};
      const float tofNSigmas[3] = {unPackInTable<aod::pidtof_tiny::binning>(part.tofNSigmaStorePr()), unPackInTable<aod::pidtof_tiny::binning>(part.tofNSigmaStorePi()), unPackInTable<aod::pidtof_tiny::binning>(part.tofNSigmaStoreKa())};

      if (!IsNSigmaCombined(part.p(), tpcNSigmas[ConfTrackChoicePartOne], tofNSigmas[ConfTrackChoicePartOne]))
        continue;

      if (part.sign() > 0) {
        qaRegistry.fill(HIST("Tracks_pos/nSigmaTPC"), part.p(), tpcNSigmas[ConfTrackChoicePartOne]);
        qaRegistry.fill(HIST("Tracks_pos/nSigmaTOF"), part.p(), tofNSigmas[ConfTrackChoicePartOne]);
        trackHistoPartOnePos.fillQA<false, false>(part);
      } else if (part.sign() < 0) {
        qaRegistry.fill(HIST("Tracks_neg/nSigmaTPC"), part.p(), tpcNSigmas[ConfTrackChoicePartOne]);
        qaRegistry.fill(HIST("Tracks_neg/nSigmaTOF"), part.p(), tofNSigmas[ConfTrackChoicePartOne]);
        trackHistoPartOneNeg.fillQA<false, false>(part);
      }
    }

    for (auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(groupPartsOne, groupPartsTwo))) {
      if (!invMCascade(p2.mLambda(), p2.mAntiLambda()))
        continue;
      if (!IsParticleCombined(p1, ConfTrackChoicePartOne))
        continue;
      // const auto& posChild = parts.iteratorAt(p2.index() - 3);
      // const auto& negChild = parts.iteratorAt(p2.index() - 2);
      // const auto& bachelor = parts.iteratorAt(p2.index() - 1);

      /// Child particles must pass this condition to be selected
      /*if (!IsParticleTPC(posChild, CascChildTable[ConfCascType1][0]) || !IsParticleTPC(negChild, CascChildTable[ConfCascType1][1]) || !IsParticleTPC(bachelor, CascChildTable[ConfCascType1][2]))
        continue;*/

      sameEventCont.setPair<false>(p1, p2, col.multNtr(), ConfUse3D, 1.0f);
    }
  }
  PROCESS_SWITCH(femtoUniversePairTaskTrackCascadeExtended, processSameEvent, "Enable processing same event for track - cascade", true);

  void processMixedEvent(FilteredFDCollisions& cols, FemtoFullParticles& /*parts*/)
  {
    ColumnBinningPolicy<aod::collision::PosZ, aod::femtouniversecollision::MultNtr> colBinning{{ConfVtxBins, ConfMultBins}, true};

    for (auto& [collision1, collision2] : soa::selfCombinations(colBinning, 5, -1, cols, cols)) {

      auto groupPartsOne = partsOne->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);
      auto groupPartsTwo = partsTwo->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision2.globalIndex(), cache);

      const auto& magFieldTesla1 = collision1.magField();
      const auto& magFieldTesla2 = collision2.magField();

      if (magFieldTesla1 != magFieldTesla2) {
        continue;
      }
      for (auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(groupPartsOne, groupPartsTwo))) {
        if (!invMCascade(p2.mLambda(), p2.mAntiLambda()))
          continue;
        if (!IsParticleCombined(p1, ConfTrackChoicePartOne))
          continue;

        // const auto& posChild = parts.iteratorAt(p2.index() - 3);
        // const auto& negChild = parts.iteratorAt(p2.index() - 2);
        // const auto& bachelor = parts.iteratorAt(p2.index() - 1);

        /// Child particles must pass this condition to be selected
        /*if (!IsParticleTPC(posChild, CascChildTable[ConfCascType1][0]) || !IsParticleTPC(negChild, CascChildTable[ConfCascType1][1]) || !IsParticleTPC(bachelor, CascChildTable[ConfCascType1][2]))
          continue;*/

        mixedEventCont.setPair<false>(p1, p2, collision1.multNtr(), ConfUse3D, 1.0f);
      }
    }
  }
  PROCESS_SWITCH(femtoUniversePairTaskTrackCascadeExtended, processMixedEvent, "Enable processing mixed event for track - cascade", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<femtoUniversePairTaskTrackCascadeExtended>(cfgc),
  };
  return workflow;
}
