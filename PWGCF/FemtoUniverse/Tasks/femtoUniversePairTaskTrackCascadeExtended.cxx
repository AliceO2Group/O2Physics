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

/// \file femtoUniversePairTaskTrackCascadeExtended.cxx
/// \brief Task for cascade correlations and QA
/// \author Barbara Chytla, WUT Warsaw, barbara.chytla@cern.ch
/// \author Shirajum Monira, WUT Warsaw, shirajum.monira@cern.ch
// o2-linter: disable=name/workflow-file

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
#include "PWGCF/FemtoUniverse/Core/femtoUtils.h"
#include "Framework/O2DatabasePDGPlugin.h"

using namespace o2;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::analysis::femto_universe;
using namespace o2::aod::pidutils;

struct femtoUniversePairTaskTrackCascadeExtended { // o2-linter: disable=name/struct

  Service<o2::framework::O2DatabasePDG> pdgMC;
  SliceCache cache;
  using FemtoFullParticles = soa::Join<aod::FDCascParticles, aod::FDExtParticles>;
  Preslice<FemtoFullParticles> perCol = aod::femtouniverseparticle::fdCollisionId;

  using FemtoRecoParticles = soa::Join<aod::FDCascParticles, aod::FDExtParticles, aod::FDMCLabels>;
  Preslice<FemtoRecoParticles> perColReco = aod::femtouniverseparticle::fdCollisionId;

  ConfigurableAxis confChildTempFitVarpTBins{"ConfChildTempFitVarpTBins", {20, 0.5, 4.05}, "V0 child: pT binning of the pT vs. TempFitVar plot"};
  ConfigurableAxis confChildTempFitVarBins{"ConfChildTempFitVarBins", {300, -0.15, 0.15}, "V0 child: binning of the TempFitVar in the pT vs. TempFitVar plot"};
  Configurable<float> confCascInvMassLowLimit{"ConfCascInvMassLowLimit", 1.315, "Lower limit of the Casc invariant mass"};
  Configurable<float> confCascInvMassUpLimit{"ConfCascInvMassUpLimit", 1.325, "Upper limit of the Casc invariant mass"};
  Configurable<float> confCascTranRad{"ConfCascTranRad", 0.5, "Cascade transverse radius"};

  Configurable<float> confNSigmaTPCPion{"NSigmaTPCPion", 4, "NSigmaTPCPion"};
  Configurable<float> confNSigmaTPCProton{"NSigmaTPCProton", 4, "NSigmaTPCProton"};

  /// applying narrow cut
  Configurable<float> confZVertexCut{"ConfZVertexCut", 10.f, "Event sel: Maximum z-Vertex (cm)"};
  Configurable<float> confEta{"ConfEta", 0.8, "Eta cut for the global track"};

  // configurations for correlation part
  Configurable<int> confTrackChoicePartOne{"ConfTrackChoicePartOne", 0, "0:Proton, 1:Pion, 2:Kaon"};
  Configurable<int> confTrkPDGCodePartOne{"ConfTrkPDGCodePartOne", 2212, "Particle 1 (Track) - PDG code"};
  Configurable<int> confCascType1{"ConfCascType1", 0, "select one of the V0s (Omega = 0, Xi = 1, anti-Omega = 2, anti-Xi = 3) for track-cascade combination"};
  Configurable<int> confCascType2{"ConfCascType2", 0, "select one of the V0s (Omega = 0, Xi = 1, anti-Omega = 2, anti-Xi = 3) for cascade-cascade combination"};
  Configurable<bool> confIsCPR{"ConfIsCPR", false, "Close Pair Rejection"};
  Configurable<float> confCPRdeltaPhiCutMax{"ConfCPRdeltaPhiCutMax", 0.0, "Delta Phi max cut for Close Pair Rejection"};
  Configurable<float> confCPRdeltaPhiCutMin{"ConfCPRdeltaPhiCutMin", 0.0, "Delta Phi min cut for Close Pair Rejection"};
  Configurable<float> confCPRdeltaEtaCutMax{"ConfCPRdeltaEtaCutMax", 0.0, "Delta Eta max cut for Close Pair Rejection"};
  Configurable<float> confCPRdeltaEtaCutMin{"ConfCPRdeltaEtaCutMin", 0.0, "Delta Eta min cut for Close Pair Rejection"};
  Configurable<bool> confCPRPlotPerRadii{"ConfCPRPlotPerRadii", false, "Plot CPR per radii"};
  Configurable<float> confCPRChosenRadii{"ConfCPRChosenRadii", 0.0, "Delta Eta cut for Close Pair Rejection"};
  Configurable<int> confChargePart1{"ConfChargePart1", 1, "sign of particle 1"};
  Configurable<float> confHPtPart1{"ConfHPtPart1", 4.0f, "higher limit for pt of particle 1"};
  Configurable<float> confLPtPart1{"ConfLPtPart1", 0.5f, "lower limit for pt of particle 1"};
  Configurable<float> confHPtPart2{"ConfHPtPart2", 4.0f, "higher limit for pt of particle 2"};
  Configurable<float> confLPtPart2{"ConfLPtPart2", 0.3f, "lower limit for pt of particle 2"};
  Configurable<float> confmom{"Confmom", 0.75, "momentum threshold for particle identification using TOF"};
  Configurable<float> confNsigmaTPCParticle{"ConfNsigmaTPCParticle", 3.0, "TPC Sigma for particle momentum < Confmom"};
  Configurable<float> confNsigmaCombinedParticle{"ConfNsigmaCombinedParticle", 3.0, "TPC and TOF Sigma (combined) for particle momentum > Confmom"};

  ConfigurableAxis confkstarBins{"ConfkstarBins", {1500, 0., 6.}, "binning kstar"};
  ConfigurableAxis confMultBins{"ConfMultBins", {VARIABLE_WIDTH, 0.0f, 20.0f, 40.0f, 60.0f, 80.0f, 100.0f, 200.0f, 99999.f}, "Mixing bins - multiplicity"};
  ConfigurableAxis confkTBins{"ConfkTBins", {150, 0., 9.}, "binning kT"};
  ConfigurableAxis confmTBins{"ConfmTBins", {225, 0., 7.5}, "binning mT"};
  ConfigurableAxis confmultBins3D{"ConfMultBins3D", {VARIABLE_WIDTH, 0.0f, 20.0f, 30.0f, 40.0f, 99999.0f}, "multiplicity Binning for the 3Dimensional plot: k* vs multiplicity vs mT (set <<confUse3D>> to true in order to use)"};
  ConfigurableAxis confmTBins3D{"ConfmTBins3D", {VARIABLE_WIDTH, 1.02f, 1.14f, 1.20f, 1.26f, 1.38f, 1.56f, 1.86f, 4.50f}, "mT Binning for the 3Dimensional plot: k* vs multiplicity vs mT (set <<confUse3D>> to true in order to use)"};
  Configurable<int> confEtaBins{"ConfEtaBins", 29, "Number of eta bins in deta dphi"};
  Configurable<int> confPhiBins{"ConfPhiBins", 29, "Number of phi bins in deta dphi"};
  Configurable<bool> confIsMC{"ConfIsMC", false, "Enable additional Histograms in the case of a MonteCarlo Run"};
  Configurable<bool> confUse3D{"ConfUse3D", false, "Enable three dimensional histogramms (to be used only for analysis with high statistics): k* vs mT vs multiplicity"};
  Configurable<bool> confUseCent{"confUseCent", false, "Use centrality in place of multiplicity"};
  ConfigurableAxis confVtxBins{"ConfVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis confTrkTempFitVarpTBins{"ConfTrkTempFitVarpTBins", {20, 0.5, 4.05}, "pT binning of the pT vs. TempFitVar plot"};
  ConfigurableAxis confTrkTempFitVarBins{"ConfTrkDTempFitVarBins", {300, -0.15, 0.15}, "binning of the TempFitVar in the pT vs. TempFitVar plot"};

  Filter collisionFilter = (nabs(aod::collision::posZ) < confZVertexCut);
  using FilteredFDCollisions = soa::Filtered<o2::aod::FdCollisions>;
  using FilteredFDCollision = FilteredFDCollisions::iterator;

  /// Partition for particle 1 (track)
  Partition<FemtoFullParticles> partsOne = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kTrack)) && (aod::femtouniverseparticle::sign == confChargePart1) && (nabs(aod::femtouniverseparticle::eta) < confEta) && (aod::femtouniverseparticle::pt < confHPtPart1) && (aod::femtouniverseparticle::pt > confLPtPart1);
  Partition<FemtoFullParticles> partsOneMCgen = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kMCTruthTrack)) && (nabs(aod::femtouniverseparticle::eta) < confEta) && (aod::femtouniverseparticle::pt < confHPtPart1) && (aod::femtouniverseparticle::pt > confLPtPart1);
  Partition<FemtoRecoParticles> partsOneMCreco = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kTrack)) && (aod::femtouniverseparticle::sign == confChargePart1) && (nabs(aod::femtouniverseparticle::eta) < confEta) && (aod::femtouniverseparticle::pt < confHPtPart1) && (aod::femtouniverseparticle::pt > confLPtPart1);

  /// Partition for particle 2 (cascade)
  Partition<FemtoFullParticles> partsTwo = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kCascade)) && (aod::femtouniverseparticle::pt < confHPtPart2) && (aod::femtouniverseparticle::pt > confLPtPart2);
  Partition<FemtoFullParticles> partsTwoMCgen = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kMCTruthTrack)) && (aod::femtouniverseparticle::pt < confHPtPart2) && (aod::femtouniverseparticle::pt > confLPtPart2);
  Partition<FemtoRecoParticles> partsTwoMCreco = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kCascade)) && (aod::femtouniverseparticle::pt < confHPtPart2) && (aod::femtouniverseparticle::pt > confLPtPart2);

  /// Partition for cascades
  Partition<FemtoFullParticles> cascs = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kCascade));

  /// Histogramming for track particle
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kTrack, 3> trackHistoPartOnePos;
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kTrack, 4> trackHistoPartOneNeg;

  /// Histogramming for cascade
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kV0Child, 3> posChildHistos;
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kV0Child, 4> negChildHistos;
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kCascadeBachelor, 0> bachHistos;
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kCascade, 0> cascQAHistos;

  /// Histogramming for Event
  FemtoUniverseEventHisto eventHisto;

  FemtoUniverseContainer<femto_universe_container::EventType::same, femto_universe_container::Observable::kstar> sameEventCont;
  FemtoUniverseContainer<femto_universe_container::EventType::mixed, femto_universe_container::Observable::kstar> mixedEventCont;
  FemtoUniversePairCleaner<aod::femtouniverseparticle::ParticleType::kTrack, aod::femtouniverseparticle::ParticleType::kCascade> pairCleaner;
  FemtoUniversePairCleaner<aod::femtouniverseparticle::ParticleType::kCascade, aod::femtouniverseparticle::ParticleType::kCascade> pairCleanerCasc;
  FemtoUniverseDetaDphiStar<aod::femtouniverseparticle::ParticleType::kCascade, aod::femtouniverseparticle::ParticleType::kCascade> pairCloseRejection;

  HistogramRegistry rXiQA{"xi", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry qaRegistry{"TrackQA", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry resultRegistry{"Correlations", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registryMCgen{"MCgenHistos", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry registryMCreco{"MCrecoHistos", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};

  // Table to select cascade daughters
  // Charges: = +--, +--, +-+, +-+
  static constexpr unsigned int CascChildTable[][3] = {{0, 1, 2}, {0, 1, 1}, {1, 0, 2}, {1, 0, 1}};

  bool invMCascade(float invMassCascade, float invMassAntiCascade)
  {
    if ((invMassCascade < confCascInvMassLowLimit || invMassCascade > confCascInvMassUpLimit) && (invMassAntiCascade < confCascInvMassLowLimit || invMassAntiCascade > confCascInvMassUpLimit)) {
      return false;
    }
    return true;
  }

  bool isNSigmaTPC(float nsigmaTPCParticle)
  {
    if (std::abs(nsigmaTPCParticle) < confNsigmaTPCParticle) {
      return true;
    } else {
      return false;
    }
  }

  bool isNSigmaCombined(float mom, float nsigmaTPCParticle, float nsigmaTOFParticle)
  {
    if (mom <= confmom) {
      return (std::abs(nsigmaTPCParticle) < confNsigmaTPCParticle);
    } else {
      return (TMath::Hypot(nsigmaTOFParticle, nsigmaTPCParticle) < confNsigmaCombinedParticle); // o2-linter: disable=root-entity
    }
  }

  template <typename T>
  bool isParticleTPC(const T& part, int id)
  {
    const float tpcNSigmas[3] = {unPackInTable<aod::pidtpc_tiny::binning>(part.tpcNSigmaStorePr()), unPackInTable<aod::pidtpc_tiny::binning>(part.tpcNSigmaStorePi()), unPackInTable<aod::pidtpc_tiny::binning>(part.tpcNSigmaStoreKa())};

    return isNSigmaTPC(tpcNSigmas[id]);
  }

  template <typename T>
  bool isParticleCombined(const T& part, int id)
  {
    const float tpcNSigmas[3] = {unPackInTable<aod::pidtpc_tiny::binning>(part.tpcNSigmaStorePr()), unPackInTable<aod::pidtpc_tiny::binning>(part.tpcNSigmaStorePi()), unPackInTable<aod::pidtpc_tiny::binning>(part.tpcNSigmaStoreKa())};
    const float tofNSigmas[3] = {unPackInTable<aod::pidtof_tiny::binning>(part.tofNSigmaStorePr()), unPackInTable<aod::pidtof_tiny::binning>(part.tofNSigmaStorePi()), unPackInTable<aod::pidtof_tiny::binning>(part.tofNSigmaStoreKa())};

    return isNSigmaCombined(part.p(), tpcNSigmas[id], tofNSigmas[id]);
  }

  void init(InitContext const&)
  {
    // Axes
    AxisSpec aXiMassAxis = {200, 1.28f, 1.36f, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec ptAxis = {100, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec etaAxis = {100, -2.0f, 2.0f, "#it{#eta}"};
    AxisSpec phiAxis = {100, 0.0f, 6.0f, "#it{#phi}"};
    AxisSpec aDCADaughAxis = {1000, 0.0f, 2.0f, "DCA (cm)"};
    AxisSpec aCPAAxis = {1000, 0.95f, 1.0f, "#it{cos #theta_{p}}"};
    AxisSpec tranRadAxis = {1000, 0.0f, 100.0f, "#it{r}_{xy} (cm)"};
    AxisSpec aDCAToPVAxis = {1000, -10.0f, 10.0f, "DCA to PV (cm)"};

    // Histograms
    rXiQA.add("hMassXi", "hMassXi", {HistType::kTH1F, {aXiMassAxis}});
    rXiQA.add("hMassXiSelected", "hMassXiSelected", {HistType::kTH1F, {aXiMassAxis}});
    rXiQA.add("hPtXi", "hPtXi", {HistType::kTH1F, {{ptAxis}}});
    rXiQA.add("hEtaXi", "hEtaXi", {HistType::kTH1F, {{etaAxis}}});
    rXiQA.add("hPhiXi", "hPhiXi", {HistType::kTH1F, {{phiAxis}}});
    rXiQA.add("hDCAV0Daughters", "hDCAV0Daughters", {HistType::kTH1F, {aDCADaughAxis}});
    rXiQA.add("hV0CosPA", "hV0CosPA", {HistType::kTH1F, {aCPAAxis}});
    rXiQA.add("hV0TranRad", "hV0TranRad", {HistType::kTH1F, {tranRadAxis}});
    rXiQA.add("hDCACascDaughters", "hDCACascDaughters", {HistType::kTH1F, {aDCADaughAxis}});
    rXiQA.add("hCascCosPA", "hCascCosPA", {HistType::kTH1F, {aCPAAxis}});
    rXiQA.add("hCascTranRad", "hCascTranRad", {HistType::kTH1F, {tranRadAxis}});
    rXiQA.add("hDcaPostoPV", "hDcaPostoPV", {HistType::kTH1F, {aDCAToPVAxis}});
    rXiQA.add("hDcaNegtoPV", "hDcaNegtoPV", {HistType::kTH1F, {aDCAToPVAxis}});
    rXiQA.add("hDcaBachtoPV", "hDcaBachtoPV", {HistType::kTH1F, {aDCAToPVAxis}});
    rXiQA.add("hDcaV0toPV", "hDcaV0toPV", {HistType::kTH1F, {aDCAToPVAxis}});

    eventHisto.init(&qaRegistry);
    qaRegistry.add("Tracks_pos/nSigmaTPC", "; #it{p} (GeV/#it{c}); n#sigma_{TPC}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("Tracks_pos/nSigmaTOF", "; #it{p} (GeV/#it{c}); n#sigma_{TOF}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("Tracks_neg/nSigmaTPC", "; #it{p} (GeV/#it{c}); n#sigma_{TPC}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("Tracks_neg/nSigmaTOF", "; #it{p} (GeV/#it{c}); n#sigma_{TOF}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});

    // MC gen
    registryMCgen.add("plus/MCgenCasc", "MC gen cascades;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 5}, {400, -1.0, 1.0}}});
    registryMCgen.add("minus/MCgenCasc", "MC gen cascades;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 5}, {400, -1.0, 1.0}}});

    registryMCgen.add("plus/MCgenAllPt", "MC gen all;#it{p}_{T} (GeV/c); #eta", {HistType::kTH1F, {{500, 0, 5}}});
    registryMCgen.add("minus/MCgenAllPt", "MC gen all;#it{p}_{T} (GeV/c); #eta", {HistType::kTH1F, {{500, 0, 5}}});

    registryMCgen.add("plus/MCgenPr", "MC gen protons;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 5}, {400, -1.0, 1.0}}});
    registryMCgen.add("minus/MCgenPr", "MC gen protons;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 5}, {400, -1.0, 1.0}}});

    registryMCgen.add("plus/MCgenPrPt", "MC gen protons;#it{p}_{T} (GeV/c)", {HistType::kTH1F, {{500, 0, 5}}});
    registryMCgen.add("minus/MCgenPrPt", "MC gen protons;#it{p}_{T} (GeV/c)", {HistType::kTH1F, {{500, 0, 5}}});

    // MC reco
    registryMCreco.add("plus/MCrecoCascade", "MC reco Cascades;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 5}, {400, -1.0, 1.0}}});
    registryMCreco.add("minus/MCrecoCascade", "MC reco Cascades;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 5}, {400, -1.0, 1.0}}});

    registryMCreco.add("plus/MCrecoAllPt", "MC reco all;#it{p}_{T} (GeV/c); #eta", {HistType::kTH1F, {{500, 0, 5}}});
    registryMCreco.add("minus/MCrecoAllPt", "MC reco all;#it{p}_{T} (GeV/c); #eta", {HistType::kTH1F, {{500, 0, 5}}});

    registryMCreco.add("plus/MCrecoPr", "MC reco protons;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 5}, {400, -1.0, 1.0}}});
    registryMCreco.add("minus/MCrecoPr", "MC reco protons;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 5}, {400, -1.0, 1.0}}});

    registryMCreco.add("plus/MCrecoPrPt", "MC reco protons;#it{p}_{T} (GeV/c)", {HistType::kTH1F, {{500, 0, 5}}});
    registryMCreco.add("minus/MCrecoPrPt", "MC reco protons;#it{p}_{T} (GeV/c)", {HistType::kTH1F, {{500, 0, 5}}});

    trackHistoPartOnePos.init(&qaRegistry, confTrkTempFitVarpTBins, confTrkTempFitVarBins, confIsMC, confTrkPDGCodePartOne);
    trackHistoPartOneNeg.init(&qaRegistry, confTrkTempFitVarpTBins, confTrkTempFitVarBins, confIsMC, confTrkPDGCodePartOne);
    posChildHistos.init(&qaRegistry, confChildTempFitVarpTBins, confChildTempFitVarBins, false, 0, true);
    negChildHistos.init(&qaRegistry, confChildTempFitVarpTBins, confChildTempFitVarBins, false, 0, true);
    bachHistos.init(&qaRegistry, confChildTempFitVarpTBins, confChildTempFitVarBins, false, 0, true, "hBachelor");
    cascQAHistos.init(&qaRegistry, confChildTempFitVarpTBins, confChildTempFitVarBins, false, 0, true);

    sameEventCont.init(&resultRegistry, confkstarBins, confMultBins, confkTBins, confmTBins, confmultBins3D, confmTBins3D, confEtaBins, confPhiBins, confIsMC, confUse3D);
    mixedEventCont.init(&resultRegistry, confkstarBins, confMultBins, confkTBins, confmTBins, confmultBins3D, confmTBins3D, confEtaBins, confPhiBins, confIsMC, confUse3D);
    pairCleaner.init(&qaRegistry);
    if (confIsCPR.value) {
      pairCloseRejection.init(&resultRegistry, &qaRegistry, confCPRdeltaPhiCutMin.value, confCPRdeltaPhiCutMax.value, confCPRdeltaEtaCutMin.value, confCPRdeltaEtaCutMax.value, confCPRChosenRadii.value, confCPRPlotPerRadii.value);
    }
  }

  void processCascades(const FilteredFDCollision& col, const FemtoFullParticles& parts)
  {
    auto groupCascs = cascs->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    // const int multCol = col.multNtr();

    for (const auto& casc : groupCascs) {
      rXiQA.fill(HIST("hMassXi"), casc.mLambda());

      // if (!invMCascade(casc.mLambda(), casc.mAntiLambda()))
      //   continue;

      const auto& posChild = parts.iteratorAt(casc.index() - 3);
      const auto& negChild = parts.iteratorAt(casc.index() - 2);
      const auto& bachelor = parts.iteratorAt(casc.index() - 1);

      // if (casc.transRadius() < confCascTranRad)
      //   continue;
      // std::cout<<std::endl;
      // std::cout<<"TYPE:"<<std::endl;
      // std::cout<<casc.partType()<<std::endl;
      //  nSigma selection for daughter and bachelor tracks
      if (casc.sign() < 0) {
        if (std::abs(posChild.tpcNSigmaPr()) > confNSigmaTPCProton) {
          continue;
        }
        if (std::abs(negChild.tpcNSigmaPi()) > confNSigmaTPCPion) {
          continue;
        }
      } else {
        if (std::abs(negChild.tpcNSigmaPr()) > confNSigmaTPCProton) {
          continue;
        }
        if (std::abs(posChild.tpcNSigmaPi()) > confNSigmaTPCPion) {
          continue;
        }
      }
      if (std::abs(bachelor.tpcNSigmaPi()) > confNSigmaTPCPion) {
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
  PROCESS_SWITCH(femtoUniversePairTaskTrackCascadeExtended, processCascades, "Enable processing cascades", false);
  /// track - cascade
  void processSameEvent(const FilteredFDCollision& col, const FemtoFullParticles& parts)
  {
    auto groupPartsOne = partsOne->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto groupPartsTwo = partsTwo->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);

    eventHisto.fillQA(col);

    const int multCol = confUseCent ? col.multV0M() : col.multNtr();

    for (const auto& part : groupPartsTwo) {
      if (!invMCascade(part.mLambda(), part.mAntiLambda()))
        continue;

      cascQAHistos.fillQA<false, true>(part);

      const auto& posChild = parts.iteratorAt(part.index() - 3);
      const auto& negChild = parts.iteratorAt(part.index() - 2);
      const auto& bachelor = parts.iteratorAt(part.index() - 1);
      /// Child particles must pass this condition to be selected
      if (!isParticleTPC(posChild, CascChildTable[confCascType1][0]) || !isParticleTPC(negChild, CascChildTable[confCascType1][1]) || !isParticleTPC(bachelor, CascChildTable[confCascType1][2]))
        continue;

      posChildHistos.fillQA<false, true>(posChild);
      negChildHistos.fillQA<false, true>(negChild);
      bachHistos.fillQABase<false, true>(bachelor, HIST("hBachelor"));
    }

    for (const auto& part : groupPartsOne) {
      /// PID plot for track particle
      const float tpcNSigmas[3] = {unPackInTable<aod::pidtpc_tiny::binning>(part.tpcNSigmaStorePr()), unPackInTable<aod::pidtpc_tiny::binning>(part.tpcNSigmaStorePi()), unPackInTable<aod::pidtpc_tiny::binning>(part.tpcNSigmaStoreKa())};
      const float tofNSigmas[3] = {unPackInTable<aod::pidtof_tiny::binning>(part.tofNSigmaStorePr()), unPackInTable<aod::pidtof_tiny::binning>(part.tofNSigmaStorePi()), unPackInTable<aod::pidtof_tiny::binning>(part.tofNSigmaStoreKa())};

      if (!isNSigmaCombined(part.p(), tpcNSigmas[confTrackChoicePartOne], tofNSigmas[confTrackChoicePartOne]))
        continue;

      if (part.sign() > 0) {
        qaRegistry.fill(HIST("Tracks_pos/nSigmaTPC"), part.p(), tpcNSigmas[confTrackChoicePartOne]);
        qaRegistry.fill(HIST("Tracks_pos/nSigmaTOF"), part.p(), tofNSigmas[confTrackChoicePartOne]);
        trackHistoPartOnePos.fillQA<false, false>(part);
      } else if (part.sign() < 0) {
        qaRegistry.fill(HIST("Tracks_neg/nSigmaTPC"), part.p(), tpcNSigmas[confTrackChoicePartOne]);
        qaRegistry.fill(HIST("Tracks_neg/nSigmaTOF"), part.p(), tofNSigmas[confTrackChoicePartOne]);
        trackHistoPartOneNeg.fillQA<false, false>(part);
      }
    }

    for (const auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(groupPartsOne, groupPartsTwo))) {
      // Cascade invariant mass cut
      if (!invMCascade(p2.mLambda(), p2.mAntiLambda()))
        continue;
      // PID
      if (!isParticleCombined(p1, confTrackChoicePartOne))
        continue;
      // track cleaning
      if (!pairCleaner.isCleanPair(p1, p2, parts)) {
        continue;
      }
      const auto& posChild = parts.iteratorAt(p2.index() - 3);
      const auto& negChild = parts.iteratorAt(p2.index() - 2);
      const auto& bachelor = parts.iteratorAt(p2.index() - 1);
      /// Child particles must pass this condition to be selected
      if (!isParticleTPC(posChild, CascChildTable[confCascType1][0]) || !isParticleTPC(negChild, CascChildTable[confCascType1][1]) || !isParticleTPC(bachelor, CascChildTable[confCascType1][2]))
        continue;

      sameEventCont.setPair<false>(p1, p2, multCol, confUse3D, 1.0f);
    }
  }
  PROCESS_SWITCH(femtoUniversePairTaskTrackCascadeExtended, processSameEvent, "Enable processing same event for track - cascade", false);
  /// cascade - cascade
  void processSameEventCasc(const FilteredFDCollision& col, const FemtoFullParticles& parts)
  {
    const auto& magFieldTesla = col.magField();

    auto groupPartsTwo = partsTwo->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);

    eventHisto.fillQA(col);

    const int multCol = confUseCent ? col.multV0M() : col.multNtr();

    for (const auto& part : groupPartsTwo) {
      if (!invMCascade(part.mLambda(), part.mAntiLambda()))
        continue;

      cascQAHistos.fillQA<false, true>(part);

      const auto& posChild = parts.iteratorAt(part.index() - 3);
      const auto& negChild = parts.iteratorAt(part.index() - 2);
      const auto& bachelor = parts.iteratorAt(part.index() - 1);
      /// Check daughters of first cascade
      if (isParticleTPC(posChild, CascChildTable[confCascType1][0]) && isParticleTPC(negChild, CascChildTable[confCascType1][1]) && isParticleTPC(bachelor, CascChildTable[confCascType1][2])) {

        posChildHistos.fillQA<false, true>(posChild);
        negChildHistos.fillQA<false, true>(negChild);
        bachHistos.fillQABase<false, true>(bachelor, HIST("hBachelor"));
      }
      /// Check daughters of second cascade
      /*if (isParticleTPC(posChild, CascChildTable[confCascType2][0]) && isParticleTPC(negChild, CascChildTable[confCascType2][1]) && isParticleTPC(bachelor, CascChildTable[confCascType2][2])) {
      }*/
    }

    auto pairProcessFunc = [&](auto& p1, auto& p2) -> void {
      // Cascade invariant mass cut for p1
      if (!invMCascade(p1.mLambda(), p1.mAntiLambda()))
        return;
      // Cascade invariant mass cut for p2
      if (!invMCascade(p2.mLambda(), p2.mAntiLambda()))
        return;
      // track cleaning
      if (!pairCleanerCasc.isCleanPair(p1, p2, parts)) {
        return;
      }
      if (confIsCPR.value) {
        if (pairCloseRejection.isClosePair(p1, p2, parts, magFieldTesla, femto_universe_container::EventType::same)) {
          return;
        }
      }
      const auto& posChild1 = parts.iteratorAt(p1.index() - 3);
      const auto& negChild1 = parts.iteratorAt(p1.index() - 2);
      const auto& bachelor1 = parts.iteratorAt(p1.index() - 1);
      /// Child particles must pass this condition to be selected
      if (!isParticleTPC(posChild1, CascChildTable[confCascType1][0]) || !isParticleTPC(negChild1, CascChildTable[confCascType1][1]) || !isParticleTPC(bachelor1, CascChildTable[confCascType1][2]))
        return;
      const auto& posChild2 = parts.iteratorAt(p2.index() - 3);
      const auto& negChild2 = parts.iteratorAt(p2.index() - 2);
      const auto& bachelor2 = parts.iteratorAt(p2.index() - 1);
      /// Child particles must pass this condition to be selected
      if (!isParticleTPC(posChild2, CascChildTable[confCascType2][0]) || !isParticleTPC(negChild2, CascChildTable[confCascType2][1]) || !isParticleTPC(bachelor2, CascChildTable[confCascType2][2]))
        return;

      sameEventCont.setPair<false>(p1, p2, multCol, confUse3D, 1.0f);
    };
    if (confCascType1 == confCascType2) {
      /// Now build the combinations for identical cascades
      for (const auto& [p1, p2] : combinations(CombinationsStrictlyUpperIndexPolicy(groupPartsTwo, groupPartsTwo))) {
        pairProcessFunc(p1, p2);
      }
    } else {
      /// Now build the combinations for non-identical cascades
      for (const auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(groupPartsTwo, groupPartsTwo))) {
        pairProcessFunc(p1, p2);
      }
    }
  }
  PROCESS_SWITCH(femtoUniversePairTaskTrackCascadeExtended, processSameEventCasc, "Enable processing same event for cascade - cascade", false);
  /// track - cascade
  void processMixedEvent(const FilteredFDCollisions& cols, const FemtoFullParticles& parts)
  {
    ColumnBinningPolicy<aod::collision::PosZ, aod::femtouniversecollision::MultNtr> colBinning{{confVtxBins, confMultBins}, true};

    for (const auto& [collision1, collision2] : soa::selfCombinations(colBinning, 5, -1, cols, cols)) {
      const int multCol = confUseCent ? collision1.multV0M() : collision1.multNtr();

      auto groupPartsOne = partsOne->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);
      auto groupPartsTwo = partsTwo->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision2.globalIndex(), cache);

      const auto& magFieldTesla1 = collision1.magField();
      const auto& magFieldTesla2 = collision2.magField();

      if (magFieldTesla1 != magFieldTesla2) {
        continue;
      }
      for (const auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(groupPartsOne, groupPartsTwo))) {
        // Cascade invariant mass cut
        if (!invMCascade(p2.mLambda(), p2.mAntiLambda()))
          continue;
        // PID
        if (!isParticleCombined(p1, confTrackChoicePartOne))
          continue;

        const auto& posChild = parts.iteratorAt(p2.index() - 3);
        const auto& negChild = parts.iteratorAt(p2.index() - 2);
        const auto& bachelor = parts.iteratorAt(p2.index() - 1);
        /// Child particles must pass this condition to be selected
        if (!isParticleTPC(posChild, CascChildTable[confCascType1][0]) || !isParticleTPC(negChild, CascChildTable[confCascType1][1]) || !isParticleTPC(bachelor, CascChildTable[confCascType1][2]))
          continue;
        // track cleaning
        if (!pairCleaner.isCleanPair(p1, p2, parts)) {
          continue;
        }

        mixedEventCont.setPair<false>(p1, p2, multCol, confUse3D, 1.0f);
      }
    }
  }
  PROCESS_SWITCH(femtoUniversePairTaskTrackCascadeExtended, processMixedEvent, "Enable processing mixed event for track - cascade", false);
  /// cascade - cascade
  void processMixedEventCasc(const FilteredFDCollisions& cols, const FemtoFullParticles& parts)
  {
    ColumnBinningPolicy<aod::collision::PosZ, aod::femtouniversecollision::MultNtr> colBinning{{confVtxBins, confMultBins}, true};

    for (const auto& [collision1, collision2] : soa::selfCombinations(colBinning, 5, -1, cols, cols)) {
      const int multCol = confUseCent ? collision1.multV0M() : collision1.multNtr();

      auto groupPartsOne = partsTwo->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);
      auto groupPartsTwo = partsTwo->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision2.globalIndex(), cache);

      const auto& magFieldTesla1 = collision1.magField();
      const auto& magFieldTesla2 = collision2.magField();

      if (magFieldTesla1 != magFieldTesla2) {
        continue;
      }
      for (const auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(groupPartsOne, groupPartsTwo))) {
        // Cascade invariant mass cut for p1
        if (!invMCascade(p1.mLambda(), p1.mAntiLambda()))
          continue;
        // Cascade invariant mass cut for p2
        if (!invMCascade(p2.mLambda(), p2.mAntiLambda()))
          continue;

        const auto& posChild1 = parts.iteratorAt(p1.index() - 3);
        const auto& negChild1 = parts.iteratorAt(p1.index() - 2);
        const auto& bachelor1 = parts.iteratorAt(p1.index() - 1);
        /// Child particles must pass this condition to be selected
        if (!isParticleTPC(posChild1, CascChildTable[confCascType1][0]) || !isParticleTPC(negChild1, CascChildTable[confCascType1][1]) || !isParticleTPC(bachelor1, CascChildTable[confCascType1][2]))
          return;
        const auto& posChild2 = parts.iteratorAt(p2.index() - 3);
        const auto& negChild2 = parts.iteratorAt(p2.index() - 2);
        const auto& bachelor2 = parts.iteratorAt(p2.index() - 1);
        /// Child particles must pass this condition to be selected
        if (!isParticleTPC(posChild2, CascChildTable[confCascType2][0]) || !isParticleTPC(negChild2, CascChildTable[confCascType2][1]) || !isParticleTPC(bachelor2, CascChildTable[confCascType2][2]))
          return;
        // track cleaning
        if (!pairCleanerCasc.isCleanPair(p1, p2, parts)) {
          continue;
        }
        if (confIsCPR.value) {
          if (pairCloseRejection.isClosePair(p1, p2, parts, magFieldTesla1, femto_universe_container::EventType::mixed)) {
            continue;
          }
        }

        mixedEventCont.setPair<false>(p1, p2, multCol, confUse3D, 1.0f);
      }
    }
  }
  PROCESS_SWITCH(femtoUniversePairTaskTrackCascadeExtended, processMixedEventCasc, "Enable processing mixed event for cascade - cascade", false);
  // MC truth
  void processSameEventMCgen(const FilteredFDCollision& col, [[maybe_unused]] const FemtoFullParticles& parts)
  {
    auto groupPartsOne = partsOneMCgen->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto groupPartsTwo = partsTwoMCgen->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);

    eventHisto.fillQA(col);

    for (const auto& part : groupPartsTwo) {
      int pdgCode = static_cast<int>(part.pidCut());
      if ((confCascType1 == 0 && pdgCode != 3334) || (confCascType1 == 2 && pdgCode != -3334) || (confCascType1 == 1 && pdgCode != 3312) || (confCascType1 == 3 && pdgCode != -3312))
        continue;

      cascQAHistos.fillQA<false, true>(part);

      for (const auto& part : groupPartsOne) {
        int pdgCode = static_cast<int>(part.pidCut());
        if (pdgCode != confTrkPDGCodePartOne)
          continue;
        const auto& pdgTrackParticle = pdgMC->GetParticle(pdgCode);
        if (!pdgTrackParticle) {
          continue;
        }

        if (pdgTrackParticle->Charge() > 0) {
          trackHistoPartOnePos.fillQA<false, false>(part);
        } else if (pdgTrackParticle->Charge() < 0) {
          trackHistoPartOneNeg.fillQA<false, false>(part);
        }
      }

      for (const auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(groupPartsOne, groupPartsTwo))) {
        if (static_cast<int>(p1.pidCut()) != confTrkPDGCodePartOne)
          continue;
        int pdgCodeCasc = static_cast<int>(p2.pidCut());
        if ((confCascType1 == 0 && pdgCodeCasc != 3334) || (confCascType1 == 2 && pdgCodeCasc != -3334) || (confCascType1 == 1 && pdgCodeCasc != 3312) || (confCascType1 == 3 && pdgCodeCasc != -3312))
          continue;
        sameEventCont.setPair<false>(p1, p2, col.multNtr(), confUse3D, 1.0f);
      }
    }
  }
  PROCESS_SWITCH(femtoUniversePairTaskTrackCascadeExtended, processSameEventMCgen, "Enable processing same event MC truth for track - cascade", false);

  void processMixedEventMCgen(const FilteredFDCollisions& cols, [[maybe_unused]] const FemtoFullParticles& parts)
  {
    ColumnBinningPolicy<aod::collision::PosZ, aod::femtouniversecollision::MultNtr> colBinning{{confVtxBins, confMultBins}, true};

    for (const auto& [collision1, collision2] : soa::selfCombinations(colBinning, 5, -1, cols, cols)) {

      auto groupPartsOne = partsOneMCgen->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);
      auto groupPartsTwo = partsTwoMCgen->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision2.globalIndex(), cache);

      const auto& magFieldTesla1 = collision1.magField();
      const auto& magFieldTesla2 = collision2.magField();

      if (magFieldTesla1 != magFieldTesla2) {
        continue;
      }
      for (const auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(groupPartsOne, groupPartsTwo))) {
        if (static_cast<int>(p1.pidCut()) != confTrkPDGCodePartOne)
          continue;
        int pdgCodeCasc = static_cast<int>(p2.pidCut());
        if ((confCascType1 == 0 && pdgCodeCasc != 3334) || (confCascType1 == 2 && pdgCodeCasc != -3334) || (confCascType1 == 1 && pdgCodeCasc != 3312) || (confCascType1 == 3 && pdgCodeCasc != -3312))
          continue;
        mixedEventCont.setPair<false>(p1, p2, collision1.multNtr(), confUse3D, 1.0f);
      }
    }
  }
  PROCESS_SWITCH(femtoUniversePairTaskTrackCascadeExtended, processMixedEventMCgen, "Enable processing mixed event MC truth for track - cascade", false);

  /// This function fills MC truth particles from derived MC table
  void processMCgen(aod::FDCascParticles const& parts)
  {
    for (const auto& part : parts) {
      if (part.partType() != uint8_t(aod::femtouniverseparticle::ParticleType::kMCTruthTrack))
        continue;

      int pdgCode = static_cast<int>(part.pidCut());
      const auto& pdgParticle = pdgMC->GetParticle(pdgCode);
      if (!pdgParticle) {
        continue;
      }

      if ((confCascType1 == 0 && pdgCode == 3334) || (confCascType1 == 1 && pdgCode == 3312)) {
        registryMCgen.fill(HIST("plus/MCgenCasc"), part.pt(), part.eta());
        continue;
      } else if ((confCascType1 == 0 && pdgCode == -3334) || (confCascType1 == 1 && pdgCode == -3312)) {
        registryMCgen.fill(HIST("minus/MCgenCasc"), part.pt(), part.eta());
        continue;
      }

      if (pdgParticle->Charge() > 0.0) {
        registryMCgen.fill(HIST("plus/MCgenAllPt"), part.pt());
      }
      if (pdgCode == 2212) {
        registryMCgen.fill(HIST("plus/MCgenPr"), part.pt(), part.eta());
        registryMCgen.fill(HIST("plus/MCgenPrPt"), part.pt());
      }

      if (pdgParticle->Charge() < 0.0) {
        registryMCgen.fill(HIST("minus/MCgenAllPt"), part.pt());
      }
      if (pdgCode == -2212) {
        registryMCgen.fill(HIST("minus/MCgenPr"), part.pt(), part.eta());
        registryMCgen.fill(HIST("minus/MCgenPrPt"), part.pt());
      }
    }
  }

  PROCESS_SWITCH(femtoUniversePairTaskTrackCascadeExtended, processMCgen, "Process MC truth data for cascades", false);

  void processMCReco(FemtoRecoParticles const& parts, aod::FdMCParticles const& mcparts)
  {
    for (const auto& part : parts) {
      auto mcPartId = part.fdMCParticleId();
      if (mcPartId == -1)
        continue; // no MC particle
      const auto& mcpart = mcparts.iteratorAt(mcPartId);
      //
      if (part.partType() == aod::femtouniverseparticle::ParticleType::kCascade) {
        if ((confCascType1 == 0 && mcpart.pdgMCTruth() == 3334) || (confCascType1 == 1 && mcpart.pdgMCTruth() == 3312)) {
          const auto& posChild = parts.iteratorAt(part.index() - 3);
          const auto& negChild = parts.iteratorAt(part.index() - 2);
          const auto& bachelor = parts.iteratorAt(part.index() - 1);
          /// Daughters that do not pass this condition are not selected
          if (isParticleTPC(posChild, CascChildTable[confCascType1][0]) && isParticleTPC(negChild, CascChildTable[confCascType1][1]) && isParticleTPC(bachelor, CascChildTable[confCascType1][2])) {
            registryMCreco.fill(HIST("plus/MCrecoCascade"), mcpart.pt(), mcpart.eta());
          }
        } else if ((confCascType1 == 0 && mcpart.pdgMCTruth() == -3334) || (confCascType1 == 1 && mcpart.pdgMCTruth() == -3312)) {
          /// Daughters that do not pass this condition are not selected
          const auto& posChild = parts.iteratorAt(part.index() - 3);
          const auto& negChild = parts.iteratorAt(part.index() - 2);
          const auto& bachelor = parts.iteratorAt(part.index() - 1);
          if (isParticleTPC(posChild, CascChildTable[confCascType1 + 2][0]) && isParticleTPC(negChild, CascChildTable[confCascType1 + 2][1]) && isParticleTPC(bachelor, CascChildTable[confCascType1 + 2][2])) {
            registryMCreco.fill(HIST("minus/MCrecoCascade"), mcpart.pt(), mcpart.eta());
          }
        }
      } else if (part.partType() == aod::femtouniverseparticle::ParticleType::kTrack) {
        if (part.sign() > 0) {
          registryMCreco.fill(HIST("plus/MCrecoAllPt"), mcpart.pt());
          if (mcpart.pdgMCTruth() == 2212 && isNSigmaCombined(part.p(), unPackInTable<aod::pidtpc_tiny::binning>(part.tpcNSigmaStorePr()), unPackInTable<aod::pidtof_tiny::binning>(part.tofNSigmaStorePr()))) {
            registryMCreco.fill(HIST("plus/MCrecoPr"), mcpart.pt(), mcpart.eta());
            registryMCreco.fill(HIST("plus/MCrecoPrPt"), mcpart.pt());
          }
        }

        if (part.sign() < 0) {
          registryMCreco.fill(HIST("minus/MCrecoAllPt"), mcpart.pt());
          if (mcpart.pdgMCTruth() == -2212 && isNSigmaCombined(part.p(), unPackInTable<aod::pidtpc_tiny::binning>(part.tpcNSigmaStorePr()), unPackInTable<aod::pidtof_tiny::binning>(part.tofNSigmaStorePr()))) {
            registryMCreco.fill(HIST("minus/MCrecoPr"), mcpart.pt(), mcpart.eta());
            registryMCreco.fill(HIST("minus/MCrecoPrPt"), mcpart.pt());
          }
        }
      }
    }
  }

  PROCESS_SWITCH(femtoUniversePairTaskTrackCascadeExtended, processMCReco, "Process MC reco data for cascades", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<femtoUniversePairTaskTrackCascadeExtended>(cfgc),
  };
  return workflow;
}
