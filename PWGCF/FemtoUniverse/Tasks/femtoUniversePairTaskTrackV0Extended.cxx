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

/// \file femtoUniversePairTaskTrackV0Extended.cxx
/// \brief Tasks that build pairs of track particles and v0s
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de
/// \author Zuzanna Chochulska, WUT Warsaw & CTU Prague, zchochul@cern.ch
/// \author Shirajum Monira, WUT Warsaw, shirajum.monira.dokt@pw.edu.pl

#include <vector>
#include <string>
#include <memory>
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
#include "Framework/O2DatabasePDGPlugin.h"
#include <TFile.h>
#include <TH1.h>

using namespace o2;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::analysis::femtoUniverse;
using namespace o2::aod::pidutils;
using namespace o2::track;

struct FemtoUniversePairTaskTrackV0Extended {

  Service<o2::framework::O2DatabasePDG> pdgMC;

  SliceCache cache;
  using FemtoFullParticles = soa::Join<aod::FDParticles, aod::FDExtParticles>;
  Preslice<FemtoFullParticles> perCol = aod::femtouniverseparticle::fdCollisionId;

  using FemtoRecoParticles = soa::Join<aod::FDParticles, aod::FDExtParticles, aod::FDMCLabels>;
  Preslice<FemtoRecoParticles> perColMC = aod::femtouniverseparticle::fdCollisionId;

  /// To apply narrow cut
  Configurable<float> confZVertexCut{"confZVertexCut", 10.f, "Event sel: Maximum z-Vertex (cm)"};
  Configurable<float> confEta{"confEta", 0.8, "Eta cut for the global track"};

  /// Particle 1 (track)
  Configurable<int> confTrkPDGCodePartOne{"confTrkPDGCodePartOne", 211, "Particle 1 (Track) - PDG code"};
  Configurable<int> confTrackChoicePartOne{"confTrackChoicePartOne", 1, "0:Proton, 1:Pion, 2:Kaon"};
  ConfigurableAxis confTrkTempFitVarBins{"confTrkTempFitVarBins", {300, -0.15, 0.15}, "binning of the TempFitVar in the pT vs. TempFitVar plot"};
  ConfigurableAxis confTrkTempFitVarpTBins{"confTrkTempFitVarpTBins", {20, 0.5, 4.05}, "pT binning of the pT vs. TempFitVar plot"};
  Configurable<int> confChargePart1{"confChargePart1", 0, "sign of particle 1"};
  Configurable<float> confHPtPart1{"confHPtPart1", 4.0f, "higher limit for pt of particle 1"};
  Configurable<float> confLPtPart1{"confLPtPart1", 0.3f, "lower limit for pt of particle 1"};
  Configurable<float> confmom{"confmom", 0.5, "momentum threshold for particle identification using TOF"};
  Configurable<float> confNsigmaTPCParticle{"confNsigmaTPCParticle", 3.0, "TPC Sigma for particle momentum < confmom"};
  Configurable<float> confNsigmaCombinedParticle{"confNsigmaCombinedParticle", 3.0, "TPC and TOF Sigma (combined) for particle momentum > confmom"};

  Filter collisionFilter = (nabs(aod::collision::posZ) < confZVertexCut);
  using FilteredFDCollisions = soa::Filtered<o2::aod::FDCollisions>;
  using FilteredFDCollision = FilteredFDCollisions::iterator;

  /// Partition for particle 1
  Partition<FemtoFullParticles> partsOne = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kTrack)) && (aod::femtouniverseparticle::sign == confChargePart1) && (nabs(aod::femtouniverseparticle::eta) < confEta) && (aod::femtouniverseparticle::pt < confHPtPart1) && (aod::femtouniverseparticle::pt > confLPtPart1);
  Partition<FemtoFullParticles> partsOneMC = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kMCTruthTrack)) && (nabs(aod::femtouniverseparticle::eta) < confEta) && (aod::femtouniverseparticle::pt < confHPtPart1) && (aod::femtouniverseparticle::pt > confLPtPart1);
  Partition<FemtoRecoParticles> partsOneMCReco = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kTrack)) && (aod::femtouniverseparticle::sign == confChargePart1) && (nabs(aod::femtouniverseparticle::eta) < confEta) && (aod::femtouniverseparticle::pt < confHPtPart1) && (aod::femtouniverseparticle::pt > confLPtPart1);

  /// Histogramming for particle 1
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kTrack, 3> trackHistoPartOnePos;
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kTrack, 4> trackHistoPartOneNeg;

  /// Particle 2 (V0)
  Configurable<int> confV0PDGCodePartTwo{"confV0PDGCodePartTwo", 3122, "Particle 2 (V0) - PDG code"};
  ConfigurableAxis confV0TempFitVarBins{"confV0TempFitVarBins", {300, 0.95, 1.}, "V0: binning of the TempFitVar in the pT vs. TempFitVar plot"};
  ConfigurableAxis confV0TempFitVarpTBins{"confV0TempFitVarpTBins", {20, 0.5, 4.05}, "V0: pT binning of the pT vs. TempFitVar plot"};
  Configurable<int> confV0Type1{"confV0Type1", 0, "select one of the V0s (lambda = 0, anti-lambda = 1, k0 = 2) for v0-v0 and Track-v0 combination"};
  Configurable<int> confV0Type2{"confV0Type2", 0, "select one of the V0s (lambda = 0, anti-lambda = 1, k0 = 2) for v0-v0 combination"};
  Configurable<float> confV0InvMassLowLimit{"confV0InvMassLowLimit", 1.10, "Lower limit of the V0 invariant mass"};
  Configurable<float> confV0InvMassUpLimit{"confV0InvMassUpLimit", 1.13, "Upper limit of the V0 invariant mass"};
  ConfigurableAxis confChildTempFitVarBins{"confChildTempFitVarBins", {300, -0.15, 0.15}, "V0 child: binning of the TempFitVar in the pT vs. TempFitVar plot"};
  ConfigurableAxis confChildTempFitVarpTBins{"confChildTempFitVarpTBins", {20, 0.5, 4.05}, "V0 child: pT binning of the pT vs. TempFitVar plot"};
  Configurable<float> confHPtPart2{"confHPtPart2", 4.0f, "higher limit for pt of particle 2"};
  Configurable<float> confLPtPart2{"confLPtPart2", 0.3f, "lower limit for pt of particle 2"};

  /// Partition for particle 2
  Partition<FemtoFullParticles> partsTwo = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kV0)) && (aod::femtouniverseparticle::pt < confHPtPart2) && (aod::femtouniverseparticle::pt > confLPtPart2);
  Partition<FemtoFullParticles> partsTwoMC = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kMCTruthTrack)) && (aod::femtouniverseparticle::pt < confHPtPart2) && (aod::femtouniverseparticle::pt > confLPtPart2);
  Partition<FemtoRecoParticles> partsTwoMCReco = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kV0)) && (aod::femtouniverseparticle::pt < confHPtPart2) && (aod::femtouniverseparticle::pt > confLPtPart2);

  /// Histogramming for particle 2
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kV0, 2> trackHistoPartTwo;
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kV0Child, 3> posChildHistos;
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kV0Child, 4> negChildHistos;

  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kV0, 2> trackHistoV0Type1;
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kV0Child, 3> posChildV0Type1;
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kV0Child, 4> negChildV0Type1;
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kV0, 2> trackHistoV0Type2;
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kV0Child, 3> posChildV0Type2;
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kV0Child, 4> negChildV0Type2;

  /// Histogramming for Event
  FemtoUniverseEventHisto eventHisto;

  /// Correlation part
  // Configurable<int> confTrackChoicePartTwo{"confTrackChoicePartTwo", 1, "0:Proton, 1:Pion, 2:Kaon"}; //not used
  Configurable<bool> confIsMC{"confIsMC", false, "Enable additional Histograms in the case of a MonteCarlo Run"};
  Configurable<bool> confUse3D{"confUse3D", false, "Enable three dimensional histogramms (to be used only for analysis with high statistics): k* vs mT vs multiplicity"};
  Configurable<bool> confUseCent{"confUseCent", false, "Use centrality in place of multiplicity"};
  ConfigurableAxis confMultBins{"confMultBins", {VARIABLE_WIDTH, 0.0f, 20.0f, 40.0f, 60.0f, 80.0f, 100.0f, 200.0f, 99999.f}, "Mixing bins - multiplicity"};
  ConfigurableAxis confVtxBins{"confVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  Configurable<int> confNEventsMix{"confNEventsMix", 5, "Number of events for mixing"};
  ConfigurableAxis confkstarBins{"confkstarBins", {1500, 0., 6.}, "binning kstar"};
  ConfigurableAxis confkTBins{"confkTBins", {150, 0., 9.}, "binning kT"};
  ConfigurableAxis confmTBins{"confmTBins", {225, 0., 7.5}, "binning mT"};
  Configurable<bool> confIsCPR{"confIsCPR", true, "Close Pair Rejection"};
  Configurable<bool> confCPRPlotPerRadii{"confCPRPlotPerRadii", false, "Plot CPR per radii"};
  Configurable<float> confCPRdeltaPhiCutMax{"confCPRdeltaPhiCutMax", 0.0, "Delta Phi max cut for Close Pair Rejection"};
  Configurable<float> confCPRdeltaPhiCutMin{"confCPRdeltaPhiCutMin", 0.0, "Delta Phi min cut for Close Pair Rejection"};
  Configurable<float> confCPRdeltaEtaCutMax{"confCPRdeltaEtaCutMax", 0.0, "Delta Eta max cut for Close Pair Rejection"};
  Configurable<float> confCPRdeltaEtaCutMin{"confCPRdeltaEtaCutMin", 0.0, "Delta Eta min cut for Close Pair Rejection"};
  Configurable<float> confCPRChosenRadii{"confCPRChosenRadii", 0.80, "Delta Eta cut for Close Pair Rejection"};
  Configurable<int> confPhiBins{"confPhiBins", 29, "Number of phi bins in deta dphi"};
  Configurable<int> confEtaBins{"confEtaBins", 29, "Number of eta bins in deta dphi"};
  ConfigurableAxis confmTBins3D{"confmTBins3D", {VARIABLE_WIDTH, 1.02f, 1.14f, 1.20f, 1.26f, 1.38f, 1.56f, 1.86f, 4.50f}, "mT Binning for the 3Dimensional plot: k* vs multiplicity vs mT (set <<confUse3D>> to true in order to use)"};
  ConfigurableAxis confMultBins3D{"confMultBins3D", {VARIABLE_WIDTH, 0.0f, 20.0f, 30.0f, 40.0f, 99999.0f}, "multiplicity Binning for the 3Dimensional plot: k* vs multiplicity vs mT (set <<confUse3D>> to true in order to use)"};

  // Efficiency
  Configurable<std::string> confLocalEfficiency{"confLocalEfficiency", "", "Local path to efficiency .root file"};

  static constexpr unsigned int V0ChildTable[][2] = {{0, 1}, {1, 0}, {1, 1}}; // Table to select the V0 children

  FemtoUniverseContainer<femtoUniverseContainer::EventType::same, femtoUniverseContainer::Observable::kstar> sameEventCont;
  FemtoUniverseContainer<femtoUniverseContainer::EventType::mixed, femtoUniverseContainer::Observable::kstar> mixedEventCont;
  FemtoUniversePairCleaner<aod::femtouniverseparticle::ParticleType::kTrack, aod::femtouniverseparticle::ParticleType::kV0> pairCleaner;
  FemtoUniversePairCleaner<aod::femtouniverseparticle::ParticleType::kV0, aod::femtouniverseparticle::ParticleType::kV0> pairCleanerV0;
  FemtoUniverseDetaDphiStar<aod::femtouniverseparticle::ParticleType::kTrack, aod::femtouniverseparticle::ParticleType::kV0> pairCloseRejection;
  FemtoUniverseDetaDphiStar<aod::femtouniverseparticle::ParticleType::kV0, aod::femtouniverseparticle::ParticleType::kV0> pairCloseRejectionV0;
  /// Histogram output
  HistogramRegistry qaRegistry{"TrackQA", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry resultRegistry{"Correlations", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registryMCtruth{"MCtruthHistos", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry registryMCreco{"MCrecoHistos", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};

  HistogramRegistry mixQaRegistry{"mixQaRegistry", {}, OutputObjHandlingPolicy::AnalysisObject};

  std::unique_ptr<TFile> plocalEffFile;
  std::unique_ptr<TH1> plocalEffp1;
  std::unique_ptr<TH1> plocalEffp2;

  bool isNSigmaCombined(float mom, float nsigmaTPCParticle, float nsigmaTOFParticle)
  {
    if (mom <= confmom) {
      return (std::abs(nsigmaTPCParticle) < confNsigmaTPCParticle);
    } else {
      return (std::hypot(nsigmaTOFParticle, nsigmaTPCParticle) < confNsigmaCombinedParticle);
    }
  }

  bool invMLambda(float invMassLambda, float invMassAntiLambda)
  {
    if ((invMassLambda < confV0InvMassLowLimit || invMassLambda > confV0InvMassUpLimit) && (invMassAntiLambda < confV0InvMassLowLimit || invMassAntiLambda > confV0InvMassUpLimit)) {
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

  template <typename T>
  bool isParticleCombined(const T& part, int id)
  {
    const float tpcNSigmas[3] = {unPackInTable<aod::pidtpc_tiny::binning>(part.tpcNSigmaStorePr()), unPackInTable<aod::pidtpc_tiny::binning>(part.tpcNSigmaStorePi()), unPackInTable<aod::pidtpc_tiny::binning>(part.tpcNSigmaStoreKa())};
    // const float tofNSigmas[3] = {part.tofNSigmaPr(), part.tofNSigmaPi(), part.tofNSigmaKa()};
    const float tofNSigmas[3] = {unPackInTable<aod::pidtof_tiny::binning>(part.tofNSigmaStorePr()), unPackInTable<aod::pidtof_tiny::binning>(part.tofNSigmaStorePi()), unPackInTable<aod::pidtof_tiny::binning>(part.tofNSigmaStoreKa())};

    return isNSigmaCombined(part.p(), tpcNSigmas[id], tofNSigmas[id]);
  }

  template <typename T>
  bool isParticleTPC(const T& part, int id)
  {
    const float tpcNSigmas[3] = {unPackInTable<aod::pidtpc_tiny::binning>(part.tpcNSigmaStorePr()), unPackInTable<aod::pidtpc_tiny::binning>(part.tpcNSigmaStorePi()), unPackInTable<aod::pidtpc_tiny::binning>(part.tpcNSigmaStoreKa())};

    return isNSigmaTPC(tpcNSigmas[id]);
  }

  void init(InitContext&)
  {
    eventHisto.init(&qaRegistry);
    qaRegistry.add("Tracks_pos/nSigmaTPC", "; #it{p} (GeV/#it{c}); n#sigma_{TPC}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("Tracks_pos/nSigmaTOF", "; #it{p} (GeV/#it{c}); n#sigma_{TOF}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("Tracks_neg/nSigmaTPC", "; #it{p} (GeV/#it{c}); n#sigma_{TPC}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("Tracks_neg/nSigmaTOF", "; #it{p} (GeV/#it{c}); n#sigma_{TOF}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    trackHistoPartOnePos.init(&qaRegistry, confTrkTempFitVarpTBins, confTrkTempFitVarBins, confIsMC, confTrkPDGCodePartOne);
    trackHistoPartOneNeg.init(&qaRegistry, confTrkTempFitVarpTBins, confTrkTempFitVarBins, confIsMC, confTrkPDGCodePartOne);
    trackHistoPartTwo.init(&qaRegistry, confV0TempFitVarpTBins, confV0TempFitVarBins, confIsMC, confV0PDGCodePartTwo, true);
    posChildHistos.init(&qaRegistry, confChildTempFitVarpTBins, confChildTempFitVarBins, false, 0, true);
    negChildHistos.init(&qaRegistry, confChildTempFitVarpTBins, confChildTempFitVarBins, false, 0, true);

    trackHistoV0Type1.init(&qaRegistry, confV0TempFitVarpTBins, confV0TempFitVarBins, confIsMC, confV0PDGCodePartTwo, true, "V0Type1");
    posChildV0Type1.init(&qaRegistry, confChildTempFitVarpTBins, confChildTempFitVarBins, false, 0, true, "posChildV0Type1");
    negChildV0Type1.init(&qaRegistry, confChildTempFitVarpTBins, confChildTempFitVarBins, false, 0, true, "negChildV0Type1");
    trackHistoV0Type2.init(&qaRegistry, confV0TempFitVarpTBins, confV0TempFitVarBins, confIsMC, confV0PDGCodePartTwo, true, "V0Type2");
    posChildV0Type2.init(&qaRegistry, confChildTempFitVarpTBins, confChildTempFitVarBins, false, 0, true, "posChildV0Type2");
    negChildV0Type2.init(&qaRegistry, confChildTempFitVarpTBins, confChildTempFitVarBins, false, 0, true, "negChildV0Type2");

    mixQaRegistry.add("MixingQA/hMECollisionBins", ";bin;Entries", kTH1F, {{120, -0.5, 119.5}});

    // MC truth
    registryMCtruth.add("plus/MCtruthLambda", "MC truth Lambdas;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 5}, {400, -1.0, 1.0}}});
    registryMCtruth.add("minus/MCtruthLambda", "MC truth Lambdas;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 5}, {400, -1.0, 1.0}}});

    registryMCtruth.add("plus/MCtruthAllPt", "MC truth all;#it{p}_{T} (GeV/c); #eta", {HistType::kTH1F, {{500, 0, 5}}});
    registryMCtruth.add("minus/MCtruthAllPt", "MC truth all;#it{p}_{T} (GeV/c); #eta", {HistType::kTH1F, {{500, 0, 5}}});

    registryMCtruth.add("plus/MCtruthPi", "MC truth pions;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 5}, {400, -1.0, 1.0}}});
    registryMCtruth.add("plus/MCtruthPr", "MC truth protons;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 5}, {400, -1.0, 1.0}}});

    registryMCtruth.add("minus/MCtruthPi", "MC truth pions;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 5}, {400, -1.0, 1.0}}});
    registryMCtruth.add("minus/MCtruthPr", "MC truth protons;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 5}, {400, -1.0, 1.0}}});

    registryMCtruth.add("plus/MCtruthPiPt", "MC truth pions;#it{p}_{T} (GeV/c)", {HistType::kTH1F, {{500, 0, 5}}});
    registryMCtruth.add("plus/MCtruthPrPt", "MC truth protons;#it{p}_{T} (GeV/c)", {HistType::kTH1F, {{500, 0, 5}}});
    registryMCtruth.add("minus/MCtruthPiPt", "MC truth pions;#it{p}_{T} (GeV/c)", {HistType::kTH1F, {{500, 0, 5}}});
    registryMCtruth.add("minus/MCtruthPrPt", "MC truth protons;#it{p}_{T} (GeV/c)", {HistType::kTH1F, {{500, 0, 5}}});

    // MC reco
    registryMCreco.add("plus/MCrecoLambda", "MC reco Lambdas;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 5}, {400, -1.0, 1.0}}});
    registryMCreco.add("plus/MCrecoLambdaChildPr", "MC reco Lambdas;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 5}, {400, -1.0, 1.0}}});
    registryMCreco.add("plus/MCrecoLambdaChildPi", "MC reco Lambdas;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 5}, {400, -1.0, 1.0}}});
    registryMCreco.add("minus/MCrecoLambda", "MC reco Lambdas;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 5}, {400, -1.0, 1.0}}});
    registryMCreco.add("minus/MCrecoLambdaChildPr", "MC reco Lambdas;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 5}, {400, -1.0, 1.0}}});
    registryMCreco.add("minus/MCrecoLambdaChildPi", "MC reco Lambdas;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 5}, {400, -1.0, 1.0}}});

    registryMCreco.add("plus/MCrecoAllPt", "MC reco all;#it{p}_{T} (GeV/c); #eta", {HistType::kTH1F, {{500, 0, 5}}});
    registryMCreco.add("minus/MCrecoAllPt", "MC reco all;#it{p}_{T} (GeV/c); #eta", {HistType::kTH1F, {{500, 0, 5}}});

    registryMCreco.add("plus/MCrecoPi", "MC reco pions;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 5}, {400, -1.0, 1.0}}});
    registryMCreco.add("plus/MCrecoPr", "MC reco protons;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 5}, {400, -1.0, 1.0}}});

    registryMCreco.add("minus/MCrecoPi", "MC reco pions;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 5}, {400, -1.0, 1.0}}});
    registryMCreco.add("minus/MCrecoPr", "MC reco protons;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 5}, {400, -1.0, 1.0}}});

    registryMCreco.add("plus/MCrecoPiPt", "MC reco pions;#it{p}_{T} (GeV/c)", {HistType::kTH1F, {{500, 0, 5}}});
    registryMCreco.add("plus/MCrecoPrPt", "MC reco protons;#it{p}_{T} (GeV/c)", {HistType::kTH1F, {{500, 0, 5}}});
    registryMCreco.add("minus/MCrecoPiPt", "MC reco pions;#it{p}_{T} (GeV/c)", {HistType::kTH1F, {{500, 0, 5}}});
    registryMCreco.add("minus/MCrecoPrPt", "MC reco protons;#it{p}_{T} (GeV/c)", {HistType::kTH1F, {{500, 0, 5}}});

    sameEventCont.init(&resultRegistry, confkstarBins, confMultBins, confkTBins, confmTBins, confMultBins3D, confmTBins3D, confEtaBins, confPhiBins, confIsMC, confUse3D);
    sameEventCont.setPDGCodes(confTrkPDGCodePartOne, confV0PDGCodePartTwo);
    mixedEventCont.init(&resultRegistry, confkstarBins, confMultBins, confkTBins, confmTBins, confMultBins3D, confmTBins3D, confEtaBins, confPhiBins, confIsMC, confUse3D);
    mixedEventCont.setPDGCodes(confTrkPDGCodePartOne, confV0PDGCodePartTwo);

    pairCleaner.init(&qaRegistry);
    pairCleanerV0.init(&qaRegistry);
    if (confIsCPR.value) {
      pairCloseRejection.init(&resultRegistry, &qaRegistry, confCPRdeltaPhiCutMin.value, confCPRdeltaPhiCutMax.value, confCPRdeltaEtaCutMin.value, confCPRdeltaEtaCutMax.value, confCPRChosenRadii.value, confCPRPlotPerRadii.value);
      pairCloseRejectionV0.init(&resultRegistry, &qaRegistry, confCPRdeltaPhiCutMin.value, confCPRdeltaPhiCutMax.value, confCPRdeltaEtaCutMin.value, confCPRdeltaEtaCutMax.value, confCPRChosenRadii.value, confCPRPlotPerRadii.value);
    }

    if (!confLocalEfficiency.value.empty()) {
      plocalEffFile = std::unique_ptr<TFile>(TFile::Open(confLocalEfficiency.value.c_str(), "read"));
      if (!plocalEffFile || plocalEffFile.get()->IsZombie())
        LOGF(fatal, "Could not load efficiency histogram from %s", confLocalEfficiency.value.c_str());
      if (doprocessSameEvent || doprocessMixedEvent) {
        plocalEffp1 = (confChargePart1 > 0) ? std::unique_ptr<TH1>(plocalEffFile.get()->Get<TH1>("PrPlus")) : std::unique_ptr<TH1>(plocalEffFile.get()->Get<TH1>("PrMinus")); // note: works only for protons for now
        plocalEffp2 = (confV0Type1 == 0) ? std::unique_ptr<TH1>(plocalEffFile.get()->Get<TH1>("Lambda")) : std::unique_ptr<TH1>(plocalEffFile.get()->Get<TH1>("AntiLambda"));
        LOGF(info, "Loaded efficiency histograms for track-V0.");
      } else if (doprocessSameEventV0 || doprocessMixedEventV0) {
        plocalEffp1 = (confV0Type1 == 0) ? std::unique_ptr<TH1>(plocalEffFile.get()->Get<TH1>("Lambda")) : std::unique_ptr<TH1>(plocalEffFile.get()->Get<TH1>("AntiLambda"));
        plocalEffp2 = (confV0Type2 == 0) ? std::unique_ptr<TH1>(plocalEffFile.get()->Get<TH1>("Lambda")) : std::unique_ptr<TH1>(plocalEffFile.get()->Get<TH1>("AntiLambda"));
        LOGF(info, "Loaded efficiency histograms for V0-V0.");
      }
    }
  }
  /// This function processes the same event for track - V0
  template <typename PartType, typename PartitionType, typename MCParticles = std::nullptr_t>
  void doSameEvent(FilteredFDCollision const& col, PartType const& parts, PartitionType& groupPartsOne, PartitionType& groupPartsTwo, [[maybe_unused]] MCParticles mcParts = nullptr)
  {
    const auto& magFieldTesla = col.magField();

    const int multCol = confUseCent ? col.multV0M() : col.multNtr();

    eventHisto.fillQA(col);

    /// Histogramming same event
    for (const auto& part : groupPartsTwo) {
      if (!invMLambda(part.mLambda(), part.mAntiLambda()))
        continue;
      const auto& posChild = parts.iteratorAt(part.index() - 2);
      const auto& negChild = parts.iteratorAt(part.index() - 1);
      /// Daughters that do not pass this condition are not selected
      if (!isParticleTPC(posChild, V0ChildTable[confV0Type1][0]) || !isParticleTPC(negChild, V0ChildTable[confV0Type1][1]))
        continue;

      trackHistoPartTwo.fillQA<false, true>(part);
      posChildHistos.fillQA<false, true>(posChild);
      negChildHistos.fillQA<false, true>(negChild);
    }

    for (const auto& part : groupPartsOne) {
      /// PID plot for particle 1
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

    /// Now build the combinations
    for (const auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(groupPartsOne, groupPartsTwo))) {
      // Lambda invariant mass cut
      if (!invMLambda(p2.mLambda(), p2.mAntiLambda()))
        continue;
      /// PID using stored binned nsigma
      if (!isParticleCombined(p1, confTrackChoicePartOne))
        continue;
      // track cleaning
      if (!pairCleaner.isCleanPair(p1, p2, parts)) {
        continue;
      }
      if (confIsCPR.value) {
        if (pairCloseRejection.isClosePair(p1, p2, parts, magFieldTesla, femtoUniverseContainer::EventType::same)) {
          continue;
        }
      }
      const auto& posChild = parts.iteratorAt(p2.index() - 2);
      const auto& negChild = parts.iteratorAt(p2.index() - 1);

      /// Daughters that do not pass this condition are not selected
      if (!isParticleTPC(posChild, V0ChildTable[confV0Type1][0]) || !isParticleTPC(negChild, V0ChildTable[confV0Type1][1]))
        continue;

      float weight = 1.0f;
      if (plocalEffp1)
        weight = plocalEffp1.get()->GetBinContent(plocalEffp1->FindBin(p1.pt(), p1.eta())) * plocalEffp2.get()->GetBinContent(plocalEffp2->FindBin(p2.pt(), p2.eta()));
      if constexpr (std::is_same<PartType, FemtoRecoParticles>::value)
        sameEventCont.setPair<true>(p1, p2, multCol, confUse3D, weight);
      else
        sameEventCont.setPair<false>(p1, p2, multCol, confUse3D, weight);
    }
  }

  void processSameEvent(FilteredFDCollision const& col, FemtoFullParticles const& parts)
  {
    auto groupPartsOne = partsOne->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto groupPartsTwo = partsTwo->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    doSameEvent(col, parts, groupPartsOne, groupPartsTwo);
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackV0Extended, processSameEvent, "Enable processing same event for track - V0", false);

  void processSameEventMCReco(FilteredFDCollision const& col, FemtoRecoParticles const& parts, aod::FDMCParticles const& mcparts)
  {
    auto groupPartsOne = partsOneMCReco->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto groupPartsTwo = partsTwoMCReco->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    doSameEvent(col, parts, groupPartsOne, groupPartsTwo, mcparts);
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackV0Extended, processSameEventMCReco, "Enable processing same event for track - V0 MC Reco", false);

  /// This function processes the same event for V0 - V0
  void processSameEventV0(FilteredFDCollision const& col, FemtoFullParticles const& parts)
  {
    const auto& magFieldTesla = col.magField();

    auto groupPartsTwo = partsTwo->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    const int multCol = confUseCent ? col.multV0M() : col.multNtr();

    eventHisto.fillQA(col);

    /// Histogramming same event
    for (const auto& part : groupPartsTwo) {
      if (!invMLambda(part.mLambda(), part.mAntiLambda()))
        continue;
      const auto& posChild = parts.iteratorAt(part.index() - 2);
      const auto& negChild = parts.iteratorAt(part.index() - 1);

      /// Check daughters of first V0 particle
      if (isParticleTPC(posChild, V0ChildTable[confV0Type1][0]) && isParticleTPC(negChild, V0ChildTable[confV0Type1][1])) {
        trackHistoV0Type1.fillQABase<false, true>(part, HIST("V0Type1"));
        posChildV0Type1.fillQABase<false, true>(posChild, HIST("posChildV0Type1"));
        negChildV0Type1.fillQABase<false, true>(negChild, HIST("negChildV0Type1"));
      }
      /// Check daughters of second V0 particle
      if (isParticleTPC(posChild, V0ChildTable[confV0Type2][0]) && isParticleTPC(negChild, V0ChildTable[confV0Type2][1])) {
        trackHistoV0Type2.fillQABase<false, true>(part, HIST("V0Type2"));
        posChildV0Type2.fillQABase<false, true>(posChild, HIST("posChildV0Type2"));
        negChildV0Type2.fillQABase<false, true>(negChild, HIST("negChildV0Type2"));
      }
    }

    auto pairProcessFunc = [&](auto& p1, auto& p2) -> void {
      // Lambda invariant mass cut for p1
      if (!invMLambda(p1.mLambda(), p1.mAntiLambda()))
        return;
      // Lambda invariant mass cut for p2
      if (!invMLambda(p2.mLambda(), p2.mAntiLambda()))
        return;
      // track cleaning
      if (!pairCleanerV0.isCleanPair(p1, p2, parts)) {
        return;
      }
      if (confIsCPR.value) {
        if (pairCloseRejectionV0.isClosePair(p1, p2, parts, magFieldTesla, femtoUniverseContainer::EventType::same)) {
          return;
        }
      }
      const auto& posChild1 = parts.iteratorAt(p1.index() - 2);
      const auto& negChild1 = parts.iteratorAt(p1.index() - 1);
      /// Daughters that do not pass this condition are not selected
      if (!isParticleTPC(posChild1, V0ChildTable[confV0Type1][0]) || !isParticleTPC(negChild1, V0ChildTable[confV0Type1][1]))
        return;

      const auto& posChild2 = parts.iteratorAt(p2.index() - 2);
      const auto& negChild2 = parts.iteratorAt(p2.index() - 1);
      /// Daughters that do not pass this condition are not selected
      if (!isParticleTPC(posChild2, V0ChildTable[confV0Type2][0]) || !isParticleTPC(negChild2, V0ChildTable[confV0Type2][1]))
        return;

      sameEventCont.setPair<false>(p1, p2, multCol, confUse3D);
    };
    if (confV0Type1 == confV0Type2) {
      /// Now build the combinations for identical V0s
      for (const auto& [p1, p2] : combinations(CombinationsStrictlyUpperIndexPolicy(groupPartsTwo, groupPartsTwo))) {
        pairProcessFunc(p1, p2);
      }
    } else {
      /// Now build the combinations for not identical identical V0s
      for (const auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(groupPartsTwo, groupPartsTwo))) {
        pairProcessFunc(p1, p2);
      }
    }
  }

  PROCESS_SWITCH(FemtoUniversePairTaskTrackV0Extended, processSameEventV0, "Enable processing same event for V0 - V0", false);

  /// This function processes MC same events for Track - V0
  void processMCSameEvent(FilteredFDCollision const& col, FemtoFullParticles const& parts)
  {
    const auto& magFieldTesla = col.magField();

    auto groupPartsOne = partsOneMC->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto groupPartsTwo = partsTwoMC->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    const int multCol = confUseCent ? col.multV0M() : col.multNtr();

    eventHisto.fillQA(col);

    /// Histogramming same event
    for (const auto& part : groupPartsTwo) {
      int pdgCode = static_cast<int>(part.pidcut());
      if ((confV0Type1 == 0 && pdgCode != 3122) || (confV0Type1 == 1 && pdgCode != -3122))
        continue;
      trackHistoPartTwo.fillQA<false, true>(part);
    }

    for (const auto& part : groupPartsOne) {
      int pdgCode = static_cast<int>(part.pidcut());
      if (pdgCode != confTrkPDGCodePartOne)
        continue;
      const auto& pdgParticle = pdgMC->GetParticle(pdgCode);
      if (!pdgParticle) {
        continue;
      }
      /// PID plot for particle 1
      if (pdgParticle->Charge() > 0.0) {
        trackHistoPartOnePos.fillQA<false, false>(part);
      } else if (pdgParticle->Charge() < 0.0) {
        trackHistoPartOneNeg.fillQA<false, false>(part);
      }
    }

    /// Now build the combinations
    for (const auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(groupPartsOne, groupPartsTwo))) {
      if (static_cast<int>(p1.pidcut()) != confTrkPDGCodePartOne)
        continue;
      int pdgCode2 = static_cast<int>(p2.pidcut());
      if ((confV0Type1 == 0 && pdgCode2 != 3122) || (confV0Type1 == 1 && pdgCode2 != -3122))
        continue;
      // track cleaning
      if (confIsCPR.value) {
        if (pairCloseRejection.isClosePair(p1, p2, parts, magFieldTesla, femtoUniverseContainer::EventType::same)) {
          continue;
        }
      }
      sameEventCont.setPair<false>(p1, p2, multCol, confUse3D);
    }
  }

  PROCESS_SWITCH(FemtoUniversePairTaskTrackV0Extended, processMCSameEvent, "Enable processing same event for MC truth track - V0", false);

  /// This function processes MC same events for V0 - V0
  void processMCSameEventV0(FilteredFDCollision const& col, FemtoFullParticles const& /*parts*/)
  {
    auto groupPartsTwo = partsTwoMC->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    const int multCol = confUseCent ? col.multV0M() : col.multNtr();

    eventHisto.fillQA(col);

    /// Histogramming same event
    for (const auto& part : groupPartsTwo) {
      int pdgCode = static_cast<int>(part.pidcut());
      if ((confV0Type1 == 0 && pdgCode != 3122) || (confV0Type1 == 1 && pdgCode != -3122))
        continue;
      trackHistoPartTwo.fillQA<false, true>(part);
    }

    auto pairProcessFunc = [&](auto& p1, auto& p2) -> void {
      int pdgCode1 = static_cast<int>(p1.pidcut());
      if ((confV0Type1 == 0 && pdgCode1 != 3122) || (confV0Type1 == 1 && pdgCode1 != -3122))
        return;
      int pdgCode2 = static_cast<int>(p2.pidcut());
      if ((confV0Type2 == 0 && pdgCode2 != 3122) || (confV0Type2 == 1 && pdgCode2 != -3122))
        return;
      sameEventCont.setPair<false>(p1, p2, multCol, confUse3D);
    };
    /// Now build the combinations
    if (confV0Type1 == confV0Type2) {
      /// Now build the combinations for identical V0s
      for (const auto& [p1, p2] : combinations(CombinationsStrictlyUpperIndexPolicy(groupPartsTwo, groupPartsTwo))) {
        pairProcessFunc(p1, p2);
      }
    } else {
      /// Now build the combinations for not identical identical V0s
      for (const auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(groupPartsTwo, groupPartsTwo))) {
        pairProcessFunc(p1, p2);
      }
    }
  }

  PROCESS_SWITCH(FemtoUniversePairTaskTrackV0Extended, processMCSameEventV0, "Enable processing same event for MC truth V0 - V0", false);

  /// This function processes the mixed event for track - V0
  template <typename PartType, typename PartitionType, typename MCParticles = std::nullptr_t>
  void doMixedEvent(FilteredFDCollisions const& cols, PartType const& parts, PartitionType& partitionOne, PartitionType& partitionTwo, [[maybe_unused]] MCParticles mcParts = nullptr)
  {
    ColumnBinningPolicy<aod::collision::PosZ, aod::femtouniversecollision::MultNtr> colBinningMult{{confVtxBins, confMultBins}, true};
    ColumnBinningPolicy<aod::collision::PosZ, aod::femtouniversecollision::MultV0M> colBinningCent{{confVtxBins, confMultBins}, true};

    auto mixedCollProcessFunc = [&](auto& collision1, auto& collision2) -> void {
      const int multCol = confUseCent ? collision1.multV0M() : collision1.multNtr();

      auto groupPartsOne = partitionOne->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);
      auto groupPartsTwo = partitionTwo->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision2.globalIndex(), cache);

      const auto& magFieldTesla1 = collision1.magField();
      const auto& magFieldTesla2 = collision2.magField();

      if (magFieldTesla1 != magFieldTesla2) {
        return;
      }

      for (const auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(groupPartsOne, groupPartsTwo))) {
        // Lambda invariant mass cut
        if (!invMLambda(p2.mLambda(), p2.mAntiLambda()))
          continue;
        /// PID using stored binned nsigma
        if (!isParticleCombined(p1, confTrackChoicePartOne))
          continue;

        const auto& posChild = parts.iteratorAt(p2.globalIndex() - 2);
        const auto& negChild = parts.iteratorAt(p2.globalIndex() - 1);
        /// Daughters that do not pass this condition are not selected
        if (!isParticleTPC(posChild, V0ChildTable[confV0Type1][0]) || !isParticleTPC(negChild, V0ChildTable[confV0Type1][1]))
          continue;

        // track cleaning
        if (!pairCleaner.isCleanPair(p1, p2, parts)) {
          continue;
        }
        if (confIsCPR.value) {
          if (pairCloseRejection.isClosePair(p1, p2, parts, magFieldTesla1, femtoUniverseContainer::EventType::mixed)) {
            continue;
          }
        }
        float weight = 1.0f;
        if (plocalEffp1)
          weight = plocalEffp1.get()->GetBinContent(plocalEffp1->FindBin(p1.pt(), p1.eta())) * plocalEffp2.get()->GetBinContent(plocalEffp2->FindBin(p2.pt(), p2.eta()));

        if constexpr (std::is_same<PartType, FemtoRecoParticles>::value)
          mixedEventCont.setPair<true>(p1, p2, multCol, confUse3D, weight);
        else
          mixedEventCont.setPair<false>(p1, p2, multCol, confUse3D, weight);
      }
    };

    if (confUseCent) {
      for (const auto& [collision1, collision2] : soa::selfCombinations(colBinningCent, confNEventsMix, -1, cols, cols)) {
        mixedCollProcessFunc(collision1, collision2);
        mixQaRegistry.fill(HIST("MixingQA/hMECollisionBins"), colBinningCent.getBin({collision1.posZ(), collision1.multV0M()}));
      }
    } else {
      for (const auto& [collision1, collision2] : soa::selfCombinations(colBinningMult, confNEventsMix, -1, cols, cols)) {
        mixedCollProcessFunc(collision1, collision2);
        mixQaRegistry.fill(HIST("MixingQA/hMECollisionBins"), colBinningMult.getBin({collision1.posZ(), collision1.multNtr()}));
      }
    }
  }

  void processMixedEvent(FilteredFDCollisions const& cols, FemtoFullParticles const& parts)
  {
    doMixedEvent(cols, parts, partsOne, partsTwo);
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackV0Extended, processMixedEvent, "Enable processing mixed event for track - V0", false);

  void processMixedEventMCReco(FilteredFDCollisions const& cols, FemtoRecoParticles const& parts, aod::FDMCParticles const& mcparts)
  {
    doMixedEvent(cols, parts, partsOneMCReco, partsTwoMCReco, mcparts);
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackV0Extended, processMixedEventMCReco, "Enable processing mixed event for track - V0 for MC Reco", false);

  /// This function processes the mixed event for V0 - V0
  void processMixedEventV0(FilteredFDCollisions const& cols, FemtoFullParticles const& parts)
  {
    ColumnBinningPolicy<aod::collision::PosZ, aod::femtouniversecollision::MultNtr> colBinningMult{{confVtxBins, confMultBins}, true};
    ColumnBinningPolicy<aod::collision::PosZ, aod::femtouniversecollision::MultV0M> colBinningCent{{confVtxBins, confMultBins}, true};

    auto mixedCollProcessFunc = [&](auto& collision1, auto& collision2) -> void {
      const int multCol = confUseCent ? collision1.multV0M() : collision1.multNtr();

      auto groupPartsOne = partsTwo->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);
      auto groupPartsTwo = partsTwo->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision2.globalIndex(), cache);

      const auto& magFieldTesla1 = collision1.magField();
      const auto& magFieldTesla2 = collision2.magField();

      if (magFieldTesla1 != magFieldTesla2) {
        return;
      }

      for (const auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(groupPartsOne, groupPartsTwo))) {
        // Lambda invariant mass cut for p1
        if (!invMLambda(p1.mLambda(), p1.mAntiLambda())) {
          continue;
        }
        // Lambda invariant mass cut for p2
        if (!invMLambda(p2.mLambda(), p2.mAntiLambda())) {
          continue;
        }

        const auto& posChild1 = parts.iteratorAt(p1.globalIndex() - 2);
        const auto& negChild1 = parts.iteratorAt(p1.globalIndex() - 1);
        /// Daughters that do not pass this condition are not selected
        if (!isParticleTPC(posChild1, V0ChildTable[confV0Type1][0]) || !isParticleTPC(negChild1, V0ChildTable[confV0Type1][1]))
          continue;

        const auto& posChild2 = parts.iteratorAt(p2.globalIndex() - 2);
        const auto& negChild2 = parts.iteratorAt(p2.globalIndex() - 1);
        /// Daughters that do not pass this condition are not selected
        if (!isParticleTPC(posChild2, V0ChildTable[confV0Type2][0]) || !isParticleTPC(negChild2, V0ChildTable[confV0Type2][1]))
          continue;

        // track cleaning
        if (!pairCleanerV0.isCleanPair(p1, p2, parts)) {
          continue;
        }
        if (confIsCPR.value) {
          if (pairCloseRejectionV0.isClosePair(p1, p2, parts, magFieldTesla1, femtoUniverseContainer::EventType::mixed)) {
            continue;
          }
        }
        mixedEventCont.setPair<false>(p1, p2, multCol, confUse3D);
      }
    };

    if (confUseCent) {
      for (const auto& [collision1, collision2] : soa::selfCombinations(colBinningCent, confNEventsMix, -1, cols, cols)) {
        mixedCollProcessFunc(collision1, collision2);
        mixQaRegistry.fill(HIST("MixingQA/hMECollisionBins"), colBinningCent.getBin({collision1.posZ(), collision1.multV0M()}));
      }
    } else {
      for (const auto& [collision1, collision2] : soa::selfCombinations(colBinningMult, confNEventsMix, -1, cols, cols)) {
        mixedCollProcessFunc(collision1, collision2);
        mixQaRegistry.fill(HIST("MixingQA/hMECollisionBins"), colBinningMult.getBin({collision1.posZ(), collision1.multNtr()}));
      }
    }
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackV0Extended, processMixedEventV0, "Enable processing mixed events for V0 - V0", false);

  /// This function processes MC mixed events for Track - V0
  void processMCMixedEvent(FilteredFDCollisions const& cols, FemtoFullParticles const& parts)
  {
    ColumnBinningPolicy<aod::collision::PosZ, aod::femtouniversecollision::MultNtr> colBinningMult{{confVtxBins, confMultBins}, true};
    ColumnBinningPolicy<aod::collision::PosZ, aod::femtouniversecollision::MultV0M> colBinningCent{{confVtxBins, confMultBins}, true};

    auto mixedCollProcessFunc = [&](auto& collision1, auto& collision2) -> void {
      const int multCol = confUseCent ? collision1.multV0M() : collision1.multNtr();

      auto groupPartsOne = partsOneMC->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);
      auto groupPartsTwo = partsTwoMC->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision2.globalIndex(), cache);

      const auto& magFieldTesla1 = collision1.magField();
      const auto& magFieldTesla2 = collision2.magField();

      if (magFieldTesla1 != magFieldTesla2) {
        return;
      }
      for (const auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(groupPartsOne, groupPartsTwo))) {
        if (static_cast<int>(p1.pidcut()) != confTrkPDGCodePartOne)
          continue;
        int pdgCode2 = static_cast<int>(p2.pidcut());
        if ((confV0Type1 == 0 && pdgCode2 != 3122) || (confV0Type1 == 1 && pdgCode2 != -3122))
          continue;
        if (confIsCPR.value) {
          if (pairCloseRejection.isClosePair(p1, p2, parts, magFieldTesla1, femtoUniverseContainer::EventType::mixed)) {
            continue;
          }
        }
        mixedEventCont.setPair<false>(p1, p2, multCol, confUse3D);
      }
    };

    if (confUseCent) {
      for (const auto& [collision1, collision2] : soa::selfCombinations(colBinningCent, confNEventsMix, -1, cols, cols)) {
        mixedCollProcessFunc(collision1, collision2);
        mixQaRegistry.fill(HIST("MixingQA/hMECollisionBins"), colBinningCent.getBin({collision1.posZ(), collision1.multV0M()}));
      }
    } else {
      for (const auto& [collision1, collision2] : soa::selfCombinations(colBinningMult, confNEventsMix, -1, cols, cols)) {
        mixedCollProcessFunc(collision1, collision2);
        mixQaRegistry.fill(HIST("MixingQA/hMECollisionBins"), colBinningMult.getBin({collision1.posZ(), collision1.multNtr()}));
      }
    }
  }

  PROCESS_SWITCH(FemtoUniversePairTaskTrackV0Extended, processMCMixedEvent, "Enable processing mixed events for MC truth track - V0", false);

  /// This function processes MC mixed events for V0 - V0
  void processMCMixedEventV0(FilteredFDCollisions const& cols, FemtoFullParticles const& /*parts*/)
  {
    ColumnBinningPolicy<aod::collision::PosZ, aod::femtouniversecollision::MultNtr> colBinningMult{{confVtxBins, confMultBins}, true};
    ColumnBinningPolicy<aod::collision::PosZ, aod::femtouniversecollision::MultV0M> colBinningCent{{confVtxBins, confMultBins}, true};

    auto mixedCollProcessFunc = [&](auto& collision1, auto& collision2) -> void {
      const int multCol = confUseCent ? collision1.multV0M() : collision1.multNtr();

      auto groupPartsOne = partsTwoMC->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);
      auto groupPartsTwo = partsTwoMC->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision2.globalIndex(), cache);

      for (const auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(groupPartsOne, groupPartsTwo))) {
        int pdgCode1 = static_cast<int>(p1.pidcut());
        if ((confV0Type1 == 0 && pdgCode1 != 3122) || (confV0Type1 == 1 && pdgCode1 != -3122))
          continue;
        int pdgCode2 = static_cast<int>(p2.pidcut());
        if ((confV0Type2 == 0 && pdgCode2 != 3122) || (confV0Type2 == 1 && pdgCode2 != -3122))
          continue;
        mixedEventCont.setPair<false>(p1, p2, multCol, confUse3D);
      }
    };

    if (confUseCent) {
      for (const auto& [collision1, collision2] : soa::selfCombinations(colBinningCent, confNEventsMix, -1, cols, cols)) {
        mixedCollProcessFunc(collision1, collision2);
        mixQaRegistry.fill(HIST("MixingQA/hMECollisionBins"), colBinningCent.getBin({collision1.posZ(), collision1.multV0M()}));
      }
    } else {
      for (const auto& [collision1, collision2] : soa::selfCombinations(colBinningMult, confNEventsMix, -1, cols, cols)) {
        mixedCollProcessFunc(collision1, collision2);
        mixQaRegistry.fill(HIST("MixingQA/hMECollisionBins"), colBinningMult.getBin({collision1.posZ(), collision1.multNtr()}));
      }
    }
  }

  PROCESS_SWITCH(FemtoUniversePairTaskTrackV0Extended, processMCMixedEventV0, "Enable processing mixed events for MC truth V0 - V0", false);
  ///--------------------------------------------MC-------------------------------------------------///

  /// This function fills MC truth particles from derived MC table
  void processMCTruth(aod::FDParticles const& parts)
  {
    for (const auto& part : parts) {
      if (part.partType() != uint8_t(aod::femtouniverseparticle::ParticleType::kMCTruthTrack))
        continue;

      int pdgCode = static_cast<int>(part.pidcut());
      const auto& pdgParticle = pdgMC->GetParticle(pdgCode);
      if (!pdgParticle) {
        continue;
      }

      if (pdgCode == 3122) {
        registryMCtruth.fill(HIST("plus/MCtruthLambda"), part.pt(), part.eta());
        continue;
      } else if (pdgCode == -3122) {
        registryMCtruth.fill(HIST("minus/MCtruthLambda"), part.pt(), part.eta());
        continue;
      }

      if (pdgParticle->Charge() > 0.0) {
        registryMCtruth.fill(HIST("plus/MCtruthAllPt"), part.pt());
      }
      if (pdgCode == 211) {
        registryMCtruth.fill(HIST("plus/MCtruthPi"), part.pt(), part.eta());
        registryMCtruth.fill(HIST("plus/MCtruthPiPt"), part.pt());
      }
      if (pdgCode == 2212) {
        registryMCtruth.fill(HIST("plus/MCtruthPr"), part.pt(), part.eta());
        registryMCtruth.fill(HIST("plus/MCtruthPrPt"), part.pt());
      }

      if (pdgParticle->Charge() < 0.0) {
        registryMCtruth.fill(HIST("minus/MCtruthAllPt"), part.pt());
      }
      if (pdgCode == -211) {
        registryMCtruth.fill(HIST("minus/MCtruthPi"), part.pt(), part.eta());
        registryMCtruth.fill(HIST("minus/MCtruthPiPt"), part.pt());
      }
      if (pdgCode == -2212) {
        registryMCtruth.fill(HIST("minus/MCtruthPr"), part.pt(), part.eta());
        registryMCtruth.fill(HIST("minus/MCtruthPrPt"), part.pt());
      }
    }
  }

  PROCESS_SWITCH(FemtoUniversePairTaskTrackV0Extended, processMCTruth, "Process MC truth data", false);

  void processMCReco(FemtoRecoParticles const& parts, aod::FDMCParticles const& mcparts)
  {
    for (const auto& part : parts) {
      auto mcPartId = part.fdMCParticleId();
      if (mcPartId == -1)
        continue; // no MC particle
      const auto& mcpart = mcparts.iteratorAt(mcPartId);
      //
      if (part.partType() == aod::femtouniverseparticle::ParticleType::kV0) {
        if (mcpart.pdgMCTruth() == 3122) {
          const auto& posChild = parts.iteratorAt(part.globalIndex() - 2);
          const auto& negChild = parts.iteratorAt(part.globalIndex() - 1);
          /// Daughters that do not pass this condition are not selected
          if (isParticleTPC(posChild, 0) && isParticleTPC(negChild, 1)) {
            registryMCreco.fill(HIST("plus/MCrecoLambda"), mcpart.pt(), mcpart.eta()); // lambda
            if (auto mcpartIdChild = posChild.fdMCParticleId(); mcpartIdChild != -1) {
              const auto& mcpartChild = mcparts.iteratorAt(mcpartIdChild);
              registryMCreco.fill(HIST("plus/MCrecoLambdaChildPr"), mcpartChild.pt(), mcpartChild.eta()); // lambda proton child
            }
            if (auto mcpartIdChild = negChild.fdMCParticleId(); mcpartIdChild != -1) {
              const auto& mcpartChild = mcparts.iteratorAt(mcpartIdChild);
              registryMCreco.fill(HIST("plus/MCrecoLambdaChildPi"), mcpartChild.pt(), mcpartChild.eta()); // lambda pion child
            }
          }
        } else if (mcpart.pdgMCTruth() == -3122) {
          const auto& posChild = parts.iteratorAt(part.globalIndex() - 2);
          const auto& negChild = parts.iteratorAt(part.globalIndex() - 1);
          /// Daughters that do not pass this condition are not selected
          if (isParticleTPC(posChild, 1) && isParticleTPC(negChild, 0)) {
            registryMCreco.fill(HIST("minus/MCrecoLambda"), mcpart.pt(), mcpart.eta()); // anti-lambda
            if (auto mcpartIdChild = posChild.fdMCParticleId(); mcpartIdChild != -1) {
              const auto& mcpartChild = mcparts.iteratorAt(mcpartIdChild);
              registryMCreco.fill(HIST("minus/MCrecoLambdaChildPi"), mcpartChild.pt(), mcpartChild.eta()); // anti-lambda pion child
            }
            if (auto mcpartIdChild = negChild.fdMCParticleId(); mcpartIdChild != -1) {
              const auto& mcpartChild = mcparts.iteratorAt(mcpartIdChild);
              registryMCreco.fill(HIST("minus/MCrecoLambdaChildPr"), mcpartChild.pt(), mcpartChild.eta()); // anti-lambda proton child
            }
          }
        }
      } else if (part.partType() == aod::femtouniverseparticle::ParticleType::kTrack) {
        if (part.sign() > 0) {
          registryMCreco.fill(HIST("plus/MCrecoAllPt"), mcpart.pt());
          if (mcpart.pdgMCTruth() == 211 && isNSigmaCombined(part.p(), unPackInTable<aod::pidtpc_tiny::binning>(part.tpcNSigmaStorePi()), unPackInTable<aod::pidtof_tiny::binning>(part.tofNSigmaStorePi()))) {
            registryMCreco.fill(HIST("plus/MCrecoPi"), mcpart.pt(), mcpart.eta());
            registryMCreco.fill(HIST("plus/MCrecoPiPt"), mcpart.pt());
          } else if (mcpart.pdgMCTruth() == 2212 && isNSigmaCombined(part.p(), unPackInTable<aod::pidtpc_tiny::binning>(part.tpcNSigmaStorePr()), unPackInTable<aod::pidtof_tiny::binning>(part.tofNSigmaStorePr()))) {
            registryMCreco.fill(HIST("plus/MCrecoPr"), mcpart.pt(), mcpart.eta());
            registryMCreco.fill(HIST("plus/MCrecoPrPt"), mcpart.pt());
          }
        }

        if (part.sign() < 0) {
          registryMCreco.fill(HIST("minus/MCrecoAllPt"), mcpart.pt());
          if (mcpart.pdgMCTruth() == -211 && isNSigmaCombined(part.p(), unPackInTable<aod::pidtpc_tiny::binning>(part.tpcNSigmaStorePi()), unPackInTable<aod::pidtof_tiny::binning>(part.tofNSigmaStorePi()))) {
            registryMCreco.fill(HIST("minus/MCrecoPi"), mcpart.pt(), mcpart.eta());
            registryMCreco.fill(HIST("minus/MCrecoPiPt"), mcpart.pt());
          } else if (mcpart.pdgMCTruth() == -2212 && isNSigmaCombined(part.p(), unPackInTable<aod::pidtpc_tiny::binning>(part.tpcNSigmaStorePr()), unPackInTable<aod::pidtof_tiny::binning>(part.tofNSigmaStorePr()))) {
            registryMCreco.fill(HIST("minus/MCrecoPr"), mcpart.pt(), mcpart.eta());
            registryMCreco.fill(HIST("minus/MCrecoPrPt"), mcpart.pt());
          }
        }
      } // partType
    }
  }

  PROCESS_SWITCH(FemtoUniversePairTaskTrackV0Extended, processMCReco, "Process MC reco data", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<FemtoUniversePairTaskTrackV0Extended>(cfgc),
  };
  return workflow;
}
