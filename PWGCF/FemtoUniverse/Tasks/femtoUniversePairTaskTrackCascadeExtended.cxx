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

/// \file femtoUniversePairTaskTrackCascadeExtended.cxx
/// \brief Task for cascade correlations and QA
/// \author Barbara Chytla, WUT Warsaw, barbara.chytla@cern.ch
/// \author Shirajum Monira, WUT Warsaw, shirajum.monira@cern.ch

#include <vector>
#include <set>
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

struct femtoUniversePairTaskTrackCascadeExtended {

  Service<o2::framework::O2DatabasePDG> pdgMC;
  SliceCache cache;
  using FemtoFullParticles = soa::Join<aod::FDParticles, aod::FDExtParticles>;
  Preslice<FemtoFullParticles> perCol = aod::femtouniverseparticle::fdCollisionId;

  using FemtoRecoParticles = soa::Join<aod::FDParticles, aod::FDExtParticles, aod::FDMCLabels>;
  Preslice<FemtoRecoParticles> perColReco = aod::femtouniverseparticle::fdCollisionId;

  ConfigurableAxis confChildTempFitVarpTBins{"confChildTempFitVarpTBins", {20, 0.5, 4.05}, "V0 child: pT binning of the pT vs. TempFitVar plot"};
  ConfigurableAxis confChildTempFitVarBins{"confChildTempFitVarBins", {300, -0.15, 0.15}, "V0 child: binning of the TempFitVar in the pT vs. TempFitVar plot"};
  Configurable<float> confCascInvMassLowLimit{"confCascInvMassLowLimit", 1.315, "Lower limit of the Casc invariant mass"};
  Configurable<float> confCascInvMassUpLimit{"confCascInvMassUpLimit", 1.325, "Upper limit of the Casc invariant mass"};
  Configurable<float> confCascTranRad{"confCascTranRad", 0.5, "Cascade transverse radius"};

  Configurable<float> confNSigmaTPCPion{"confNSigmaTPCPion", 4, "NSigmaTPCPion"};
  Configurable<float> confNSigmaTPCProton{"confNSigmaTPCProton", 4, "NSigmaTPCProton"};

  /// applying narrow cut
  Configurable<float> confZVertexCut{"confZVertexCut", 10.f, "Event sel: Maximum z-Vertex (cm)"};
  Configurable<float> confEta{"confEta", 0.8, "Eta cut for the global track"};

  // configurations for correlation part
  Configurable<int> confTrackChoicePartOne{"confTrackChoicePartOne", 0, "0:Proton, 1:Pion, 2:Kaon"};
  Configurable<int> confTrkPDGCodePartOne{"confTrkPDGCodePartOne", 2212, "Particle 1 (Track) - PDG code"};
  Configurable<int> confCascType1{"confCascType1", 0, "select one of the Cascades (Omega = 0, Xi = 1, anti-Omega = 2, anti-Xi = 3) for track-cascade combination"};
  Configurable<int> confCascType2{"confCascType2", 0, "select one of the Cascades (Omega = 0, Xi = 1, anti-Omega = 2, anti-Xi = 3) for cascade-cascade combination"};
  Configurable<bool> confIsCPR{"confIsCPR", false, "Close Pair Rejection"};
  Configurable<float> confCPRdeltaPhiCutMax{"confCPRdeltaPhiCutMax", 0.0, "Delta Phi max cut for Close Pair Rejection"};
  Configurable<float> confCPRdeltaPhiCutMin{"confCPRdeltaPhiCutMin", 0.0, "Delta Phi min cut for Close Pair Rejection"};
  Configurable<float> confCPRdeltaEtaCutMax{"confCPRdeltaEtaCutMax", 0.0, "Delta Eta max cut for Close Pair Rejection"};
  Configurable<float> confCPRdeltaEtaCutMin{"confCPRdeltaEtaCutMin", 0.0, "Delta Eta min cut for Close Pair Rejection"};
  Configurable<bool> confCPRPlotPerRadii{"confCPRPlotPerRadii", false, "Plot CPR per radii"};
  Configurable<float> confCPRChosenRadii{"confCPRChosenRadii", 0.0, "Delta Eta cut for Close Pair Rejection"};
  Configurable<bool> confIsSameSignCPR{"confIsSameSignCPR", false, "Close Pair Rejection for same sign children of cascades"};
  Configurable<int> confChargePart1{"confChargePart1", 1, "sign of particle 1"};
  Configurable<float> confHPtPart1{"confHPtPart1", 4.0f, "higher limit for pt of particle 1"};
  Configurable<float> confLPtPart1{"confLPtPart1", 0.5f, "lower limit for pt of particle 1"};
  Configurable<float> confHPtPart2{"confHPtPart2", 4.0f, "higher limit for pt of particle 2"};
  Configurable<float> confLPtPart2{"confLPtPart2", 0.3f, "lower limit for pt of particle 2"};
  Configurable<float> confmom{"confmom", 0.75, "momentum threshold for particle identification using TOF"};
  Configurable<float> confNsigmaTPCParticle{"confNsigmaTPCParticle", 3.0, "TPC Sigma for particle (track) momentum < Confmom"};
  Configurable<float> confNsigmaCombinedParticle{"confNsigmaCombinedParticle", 3.0, "TPC and TOF Sigma (combined) for particle (track) momentum > Confmom"};
  Configurable<float> confNsigmaTPCParticleChild{"confNsigmaTPCParticleChild", 3.0, "TPC Sigma for particle (daugh & bach) momentum < Confmom"};
  Configurable<float> confNsigmaTOFParticleChild{"confNsigmaTOFParticleChild", 3.0, "TOF Sigma for particle (daugh & bach) momentum > Confmom"};

  ConfigurableAxis confkstarBins{"confkstarBins", {1500, 0., 6.}, "binning kstar"};
  ConfigurableAxis confMultBins{"confMultBins", {VARIABLE_WIDTH, 0.0f, 20.0f, 40.0f, 60.0f, 80.0f, 100.0f, 200.0f, 99999.f}, "Mixing bins - multiplicity"};
  ConfigurableAxis confkTBins{"confkTBins", {150, 0., 9.}, "binning kT"};
  ConfigurableAxis confmTBins{"confmTBins", {225, 0., 7.5}, "binning mT"};
  ConfigurableAxis confMultBins3D{"confMultBins3D", {VARIABLE_WIDTH, 0.0f, 20.0f, 30.0f, 40.0f, 99999.0f}, "multiplicity Binning for the 3Dimensional plot: k* vs multiplicity vs mT (set <<confUse3D>> to true in order to use)"};
  ConfigurableAxis confmTBins3D{"confmTBins3D", {VARIABLE_WIDTH, 1.02f, 1.14f, 1.20f, 1.26f, 1.38f, 1.56f, 1.86f, 4.50f}, "mT Binning for the 3Dimensional plot: k* vs multiplicity vs mT (set <<confUse3D>> to true in order to use)"};
  Configurable<int> confEtaBins{"confEtaBins", 29, "Number of eta bins in deta dphi"};
  Configurable<int> confPhiBins{"confPhiBins", 29, "Number of phi bins in deta dphi"};
  Configurable<bool> confIsMC{"confIsMC", false, "Enable additional Histograms in the case of a MonteCarlo Run"};
  Configurable<bool> confUse3D{"confUse3D", false, "Enable three dimensional histogramms (to be used only for analysis with high statistics): k* vs mT vs multiplicity"};
  Configurable<bool> confUseCent{"confUseCent", false, "Use centrality in place of multiplicity"};
  ConfigurableAxis confVtxBins{"confVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis confTrkTempFitVarpTBins{"confTrkTempFitVarpTBins", {20, 0.5, 4.05}, "pT binning of the pT vs. TempFitVar plot"};
  ConfigurableAxis confTrkTempFitVarBins{"confTrkTempFitVarBins", {300, -0.15, 0.15}, "binning of the TempFitVar in the pT vs. TempFitVar plot"};

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

  std::set<int> cascDuplicates;

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

  bool isNSigmaCombined(float mom, float nsigmaTPCParticle, float nsigmaTOFParticle)
  {
    if (mom <= confmom) {
      return (std::abs(nsigmaTPCParticle) < confNsigmaTPCParticle);
    } else {
      return (TMath::Hypot(nsigmaTOFParticle, nsigmaTPCParticle) < confNsigmaCombinedParticle);
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

  template <typename T>
  bool isParticleCombined(const T& part, int id)
  {
    const float tpcNSigmas[3] = {unPackInTable<aod::pidtpc_tiny::binning>(part.tpcNSigmaStorePr()), unPackInTable<aod::pidtpc_tiny::binning>(part.tpcNSigmaStorePi()), unPackInTable<aod::pidtpc_tiny::binning>(part.tpcNSigmaStoreKa())};
    const float tofNSigmas[3] = {unPackInTable<aod::pidtof_tiny::binning>(part.tofNSigmaStorePr()), unPackInTable<aod::pidtof_tiny::binning>(part.tofNSigmaStorePi()), unPackInTable<aod::pidtof_tiny::binning>(part.tofNSigmaStoreKa())};

    return isNSigmaCombined(part.p(), tpcNSigmas[id], tofNSigmas[id]);
  }

  void init(InitContext const&)
  {
    std::vector<double> multBinning = {0.0, 5.0, 10.0, 20.0, 30.0f, 40.0, 50.0, 60.0f, 70.0, 80.0, 100.0, 200.0, 99999.0};
    // Axes
    AxisSpec aXiMassAxis = {200, 1.28f, 1.36f, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec ptAxis = {100, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec etaAxis = {100, -2.0f, 2.0f, "#it{#eta}"};
    AxisSpec phiAxis = {100, 0.0f, 6.0f, "#it{#phi}"};
    AxisSpec aDCADaughAxis = {1000, 0.0f, 2.0f, "DCA (cm)"};
    AxisSpec aCPAAxis = {1000, 0.95f, 1.0f, "#it{cos #theta_{p}}"};
    AxisSpec tranRadAxis = {1000, 0.0f, 100.0f, "#it{r}_{xy} (cm)"};
    AxisSpec aDCAToPVAxis = {1000, -10.0f, 10.0f, "DCA to PV (cm)"};
    AxisSpec multAxis = {multBinning, "Multiplicity"};

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
    rXiQA.add("hInvMpT", "hInvMpT", kTH2F, {{ptAxis}, {aXiMassAxis}});
    rXiQA.add("hInvMpTmult", "hInvMpTmult", kTH3F, {{ptAxis}, {aXiMassAxis}, {multAxis}});

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

    sameEventCont.init(&resultRegistry, confkstarBins, confMultBins, confkTBins, confmTBins, confMultBins3D, confmTBins3D, confEtaBins, confPhiBins, confIsMC, confUse3D);
    mixedEventCont.init(&resultRegistry, confkstarBins, confMultBins, confkTBins, confmTBins, confMultBins3D, confmTBins3D, confEtaBins, confPhiBins, confIsMC, confUse3D);
    pairCleaner.init(&qaRegistry);
    pairCleanerCasc.init(&qaRegistry);
    if (confIsCPR.value) {
      pairCloseRejection.init(&resultRegistry, &qaRegistry, confCPRdeltaPhiCutMin.value, confCPRdeltaPhiCutMax.value, confCPRdeltaEtaCutMin.value, confCPRdeltaEtaCutMax.value, confCPRChosenRadii.value, confCPRPlotPerRadii.value, 0, 0, confIsSameSignCPR.value);
    }
  }

  void processCascades([[maybe_unused]] const FilteredFDCollision& col, const FemtoFullParticles& parts, const aod::FDCascParticles& fdcascs)
  {
    // auto groupCascs = cascs->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    //  const int multCol = col.multNtr();

    for (const auto& casc : fdcascs) {
      const auto& part = casc.fdParticle_as<FemtoFullParticles>();
      rXiQA.fill(HIST("hMassXi"), part.mLambda());

      // if (!invMCascade(casc.mLambda(), casc.mAntiLambda()))
      //   continue;

      const auto& posChild = parts.iteratorAt(part.globalIndex() - 3 - parts.begin().globalIndex());
      const auto& negChild = parts.iteratorAt(part.globalIndex() - 2 - parts.begin().globalIndex());
      const auto& bachelor = parts.iteratorAt(part.globalIndex() - 1 - parts.begin().globalIndex());

      // if (casc.transRadius() < confCascTranRad)
      //   continue;
      // std::cout<<std::endl;
      // std::cout<<"TYPE:"<<std::endl;
      // std::cout<<casc.partType()<<std::endl;
      //  nSigma selection for daughter and bachelor tracks

      if (part.sign() < 0) {
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

      rXiQA.fill(HIST("hPtXi"), part.pt());
      rXiQA.fill(HIST("hEtaXi"), part.eta());
      rXiQA.fill(HIST("hPhiXi"), part.phi());
      rXiQA.fill(HIST("hMassXiSelected"), part.mLambda());
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
      rXiQA.fill(HIST("hInvMpT"), part.pt(), part.mLambda());

      posChildHistos.fillQA<false, true>(posChild);
      negChildHistos.fillQA<false, true>(negChild);
      bachHistos.fillQABase<false, true>(bachelor, HIST("hBachelor"));
    }
  }
  PROCESS_SWITCH(femtoUniversePairTaskTrackCascadeExtended, processCascades, "Enable processing cascades", false);
  /// track - cascade correlations
  void processSameEvent(const FilteredFDCollision& col, const FemtoFullParticles& parts)
  {
    auto groupPartsOne = partsOne->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto groupPartsTwo = partsTwo->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);

    eventHisto.fillQA(col);

    const int multCol = confUseCent ? col.multV0M() : col.multNtr();

    for (const auto& part : groupPartsTwo) {
      if (!invMCascade(part.mLambda(), part.mAntiLambda(), confCascType1)) /// mLambda stores Xi mass, mAntiLambda stores Omega mass
        continue;

      cascQAHistos.fillQA<false, true>(part);

      const auto& posChild = parts.iteratorAt(part.index() - 3);
      const auto& negChild = parts.iteratorAt(part.index() - 2);
      const auto& bachelor = parts.iteratorAt(part.index() - 1);
      /// Child particles must pass this condition to be selected
      if (!isParticleTPC(posChild, CascChildTable[confCascType1][0]) || !isParticleTPC(negChild, CascChildTable[confCascType1][1]) || !isParticleTPC(bachelor, CascChildTable[confCascType1][2]))
        continue;

      if (!isParticleTOF(posChild, CascChildTable[confCascType1][0]) || !isParticleTOF(negChild, CascChildTable[confCascType1][1]) || !isParticleTOF(bachelor, CascChildTable[confCascType1][2]))
        continue;

      posChildHistos.fillQA<false, true>(posChild);
      negChildHistos.fillQA<false, true>(negChild);
      bachHistos.fillQABase<false, true>(bachelor, HIST("hBachelor"));

      rXiQA.fill(HIST("hInvMpTmult"), part.pt(), part.mLambda(), multCol);
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
      // Cascade inv mass cut (mLambda stores Xi mass, mAntiLambda stores Omega mass)
      if (!invMCascade(p2.mLambda(), p2.mAntiLambda(), confCascType1))
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
      if (!isParticleTOF(posChild, CascChildTable[confCascType1][0]) || !isParticleTOF(negChild, CascChildTable[confCascType1][1]) || !isParticleTOF(bachelor, CascChildTable[confCascType1][2]))
        continue;

      sameEventCont.setPair<false>(p1, p2, multCol, confUse3D, 1.0f);
    }
  }
  PROCESS_SWITCH(femtoUniversePairTaskTrackCascadeExtended, processSameEvent, "Enable processing same event for track - cascade", false);
  /// cascade - cascade correlations
  void processSameEventCasc(const FilteredFDCollision& col, const FemtoFullParticles& parts)
  {
    const auto& magFieldTesla = col.magField();

    auto groupPartsTwo = partsTwo->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);

    eventHisto.fillQA(col);

    const int multCol = confUseCent ? col.multV0M() : col.multNtr();

    for (const auto& part : groupPartsTwo) {
      if (!invMCascade(part.mLambda(), part.mAntiLambda(), confCascType1)) /// mLambda stores Xi mass, mAntiLambda stores Omega mass
        continue;

      cascQAHistos.fillQA<false, true>(part);

      const auto& posChild = parts.iteratorAt(part.index() - 3);
      const auto& negChild = parts.iteratorAt(part.index() - 2);
      const auto& bachelor = parts.iteratorAt(part.index() - 1);
      /// Check daughters of first cascade
      if (!isParticleTPC(posChild, CascChildTable[confCascType1][0]) || !isParticleTPC(negChild, CascChildTable[confCascType1][1]) || !isParticleTPC(bachelor, CascChildTable[confCascType1][2]))
        continue;
      if (!isParticleTOF(posChild, CascChildTable[confCascType1][0]) || !isParticleTOF(negChild, CascChildTable[confCascType1][1]) || !isParticleTOF(bachelor, CascChildTable[confCascType1][2]))
        continue;

      posChildHistos.fillQA<false, true>(posChild);
      negChildHistos.fillQA<false, true>(negChild);
      bachHistos.fillQABase<false, true>(bachelor, HIST("hBachelor"));
      /// Check daughters of second cascade
      /*if (isParticleTPC(posChild, CascChildTable[confCascType2][0]) && isParticleTPC(negChild, CascChildTable[confCascType2][1]) && isParticleTPC(bachelor, CascChildTable[confCascType2][2])) {
      }*/
    }

    auto pairDuplicateCheckFunc = [&](auto& p1, auto& p2) -> void {
      // Cascade inv mass cut for p1 (mLambda stores Xi mass, mAntiLambda stores Omega mass)
      if (!invMCascade(p1.mLambda(), p1.mAntiLambda(), confCascType1))
        return;
      // Cascade inv mass cut for p2
      if (!invMCascade(p2.mLambda(), p2.mAntiLambda(), confCascType2))
        return;
      // track cleaning & checking for duplicate pairs
      if (!pairCleanerCasc.isCleanPair(p1, p2, parts)) {
        // mark for rejection the cascades that share a daughter with other cascades
        cascDuplicates.insert(p1.globalIndex());
        cascDuplicates.insert(p2.globalIndex());
      }
    };

    auto pairProcessFunc = [&](auto& p1, auto& p2) -> void {
      if (cascDuplicates.contains(p1.globalIndex()) || cascDuplicates.contains(p2.globalIndex()))
        return;
      if (!invMCascade(p1.mLambda(), p1.mAntiLambda(), confCascType1))
        return;
      if (!invMCascade(p2.mLambda(), p2.mAntiLambda(), confCascType2))
        return;
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
      if (!isParticleTOF(posChild1, CascChildTable[confCascType1][0]) || !isParticleTOF(negChild1, CascChildTable[confCascType1][1]) || !isParticleTOF(bachelor1, CascChildTable[confCascType1][2]))
        return;
      const auto& posChild2 = parts.iteratorAt(p2.index() - 3);
      const auto& negChild2 = parts.iteratorAt(p2.index() - 2);
      const auto& bachelor2 = parts.iteratorAt(p2.index() - 1);
      /// Child particles must pass this condition to be selected
      if (!isParticleTPC(posChild2, CascChildTable[confCascType2][0]) || !isParticleTPC(negChild2, CascChildTable[confCascType2][1]) || !isParticleTPC(bachelor2, CascChildTable[confCascType2][2]))
        return;
      if (!isParticleTOF(posChild2, CascChildTable[confCascType2][0]) || !isParticleTOF(negChild2, CascChildTable[confCascType2][1]) || !isParticleTOF(bachelor2, CascChildTable[confCascType2][2]))
        return;

      sameEventCont.setPair<false>(p1, p2, multCol, confUse3D, 1.0f);
    };
    cascDuplicates.clear();
    if (confCascType1 == confCascType2) {
      for (const auto& [p1, p2] : combinations(CombinationsStrictlyUpperIndexPolicy(groupPartsTwo, groupPartsTwo))) {
        pairDuplicateCheckFunc(p1, p2);
      }
      /// Now build the combinations for identical cascades
      for (const auto& [p1, p2] : combinations(CombinationsStrictlyUpperIndexPolicy(groupPartsTwo, groupPartsTwo))) {
        pairProcessFunc(p1, p2);
      }
    } else {
      for (const auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(groupPartsTwo, groupPartsTwo))) {
        pairDuplicateCheckFunc(p1, p2);
      }
      /// Now build the combinations for non-identical cascades
      for (const auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(groupPartsTwo, groupPartsTwo))) {
        pairProcessFunc(p1, p2);
      }
    }
  }
  PROCESS_SWITCH(femtoUniversePairTaskTrackCascadeExtended, processSameEventCasc, "Enable processing same event for cascade - cascade", false);
  /// track - cascade correlations
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
        // Cascade inv mass cut (mLambda stores Xi mass, mAntiLambda stores Omega mass)
        if (!invMCascade(p2.mLambda(), p2.mAntiLambda(), confCascType1))
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
        if (!isParticleTOF(posChild, CascChildTable[confCascType1][0]) || !isParticleTOF(negChild, CascChildTable[confCascType1][1]) || !isParticleTOF(bachelor, CascChildTable[confCascType1][2]))
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
  /// cascade - cascade correlations
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
        // Cascade inv mass cut for p1 (mLambda stores Xi mass, mAntiLambda stores Omega mass)
        if (!invMCascade(p1.mLambda(), p1.mAntiLambda(), confCascType1))
          continue;
        // Cascade inv mass cut for p2
        if (!invMCascade(p2.mLambda(), p2.mAntiLambda(), confCascType2))
          continue;

        const auto& posChild1 = parts.iteratorAt(p1.index() - 3);
        const auto& negChild1 = parts.iteratorAt(p1.index() - 2);
        const auto& bachelor1 = parts.iteratorAt(p1.index() - 1);
        /// Child particles must pass this condition to be selected
        if (!isParticleTPC(posChild1, CascChildTable[confCascType1][0]) || !isParticleTPC(negChild1, CascChildTable[confCascType1][1]) || !isParticleTPC(bachelor1, CascChildTable[confCascType1][2]))
          return;
        if (!isParticleTOF(posChild1, CascChildTable[confCascType1][0]) || !isParticleTOF(negChild1, CascChildTable[confCascType1][1]) || !isParticleTOF(bachelor1, CascChildTable[confCascType1][2]))
          return;
        const auto& posChild2 = parts.iteratorAt(p2.index() - 3);
        const auto& negChild2 = parts.iteratorAt(p2.index() - 2);
        const auto& bachelor2 = parts.iteratorAt(p2.index() - 1);
        /// Child particles must pass this condition to be selected
        if (!isParticleTPC(posChild2, CascChildTable[confCascType2][0]) || !isParticleTPC(negChild2, CascChildTable[confCascType2][1]) || !isParticleTPC(bachelor2, CascChildTable[confCascType2][2]))
          return;
        if (!isParticleTOF(posChild2, CascChildTable[confCascType2][0]) || !isParticleTOF(negChild2, CascChildTable[confCascType2][1]) || !isParticleTOF(bachelor2, CascChildTable[confCascType2][2]))
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
    const int multCol = confUseCent ? col.multV0M() : col.multNtr();

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
        sameEventCont.setPair<false>(p1, p2, multCol, confUse3D, 1.0f);
      }
    }
  }
  PROCESS_SWITCH(femtoUniversePairTaskTrackCascadeExtended, processSameEventMCgen, "Enable processing same event MC truth for track - cascade", false);

  void processMixedEventMCgen(const FilteredFDCollisions& cols, [[maybe_unused]] const FemtoFullParticles& parts)
  {
    ColumnBinningPolicy<aod::collision::PosZ, aod::femtouniversecollision::MultNtr> colBinning{{confVtxBins, confMultBins}, true};

    for (const auto& [collision1, collision2] : soa::selfCombinations(colBinning, 5, -1, cols, cols)) {
      const int multCol = confUseCent ? collision1.multV0M() : collision1.multNtr();

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
        mixedEventCont.setPair<false>(p1, p2, multCol, confUse3D, 1.0f);
      }
    }
  }
  PROCESS_SWITCH(femtoUniversePairTaskTrackCascadeExtended, processMixedEventMCgen, "Enable processing mixed event MC truth for track - cascade", false);

  /// This function fills MC truth particles from derived MC table
  void processMCgen(aod::FDParticles const& parts)
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
