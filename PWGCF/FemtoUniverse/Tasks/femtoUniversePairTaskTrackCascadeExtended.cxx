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

#include "PWGCF/FemtoUniverse/Core/FemtoUniverseContainer.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseDetaDphiStar.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseEventHisto.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniversePairCleaner.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseParticleHisto.h"
#include "PWGCF/FemtoUniverse/Core/femtoUtils.h"
#include "PWGCF/FemtoUniverse/DataModel/FemtoDerived.h"

#include "CCDB/BasicCCDBManager.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/StepTHn.h"
#include "Framework/runDataProcessing.h"

#include <TPDGCode.h>

#include <memory>
#include <set>
#include <string>
#include <vector>

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
  using FemtoRecoFullParticles = soa::Join<aod::FDParticles, aod::FDExtParticles, aod::FDMCLabels>;
  using FemtoRecoBasicParticles = soa::Join<aod::FDParticles, aod::FDMCLabels>;

  ConfigurableAxis confChildTempFitVarpTBins{"confChildTempFitVarpTBins", {20, 0.5, 4.05}, "V0 child: pT binning of the pT vs. TempFitVar plot"};
  ConfigurableAxis confChildTempFitVarBins{"confChildTempFitVarBins", {300, -0.15, 0.15}, "V0 child: binning of the TempFitVar in the pT vs. TempFitVar plot"};
  Configurable<float> confCascInvMassLowLimit{"confCascInvMassLowLimit", 1.315, "Lower limit of the Casc invariant mass"};
  Configurable<float> confCascInvMassUpLimit{"confCascInvMassUpLimit", 1.325, "Upper limit of the Casc invariant mass"};

  // TODO: Add seperate selection for daughter particles
  // Configurable<float> confNSigmaTPCPion{"confNSigmaTPCPion", 4, "NSigmaTPCPion"};
  // Configurable<float> confNSigmaTPCProton{"confNSigmaTPCProton", 4, "NSigmaTPCProton"};

  /// applying narrow cut
  Configurable<float> confZVertexCut{"confZVertexCut", 10.f, "Event sel: Maximum z-Vertex (cm)"};
  Configurable<float> confEta{"confEta", 0.8, "Eta cut for the global track"};

  // configurations for correlation part
  Configurable<int> confTrackChoicePartOne{"confTrackChoicePartOne", 0, "0:Proton, 1:Pion, 2:Kaon"};
  Configurable<int> confTrkPDGCodePartOne{"confTrkPDGCodePartOne", 2212, "Particle 1 (Track) - PDG code"};
  Configurable<int> confCascPDGCodePartTwo{"confCascPDGCodePartTwo", 3312, "Particle 2 (Cascade) - PDG code"};
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
  Configurable<bool> confUseStrangenessTOF{"confUseStrangenessTOF", true, "Use strangeness TOF for cascade PID"};

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
  Configurable<int> confNEventsMix{"confNEventsMix", 5, "Number of events for mixing"};

  // Efficiency
  Configurable<std::string> confLocalEfficiency{"confLocalEfficiency", "", "Local path to efficiency .root file"};
  Configurable<std::string> confCCDBEfficiency{"confCCDBEfficiency", "", "CCDB path to efficiency object"};

  Filter collisionFilter = (nabs(aod::collision::posZ) < confZVertexCut);
  using FilteredFDCollisions = soa::Filtered<o2::aod::FdCollisions>;
  using FilteredFDCollision = FilteredFDCollisions::iterator;

  /// Partition for particle 1 using extended table (track)
  Partition<FemtoFullParticles> partsOneFull = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kTrack)) && (aod::femtouniverseparticle::mAntiLambda == confChargePart1) && (nabs(aod::femtouniverseparticle::eta) < confEta) && (aod::femtouniverseparticle::pt < confHPtPart1) && (aod::femtouniverseparticle::pt > confLPtPart1);

  /// Partition for particle 1 without extended table (track)
  Partition<aod::FDParticles> partsOneBasic = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kTrack)) && (aod::femtouniverseparticle::mAntiLambda == confChargePart1) && (nabs(aod::femtouniverseparticle::eta) < confEta) && (aod::femtouniverseparticle::pt < confHPtPart1) && (aod::femtouniverseparticle::pt > confLPtPart1);
  Partition<aod::FDParticles> partsOneMCgenBasic = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kMCTruthTrack)) && (nabs(aod::femtouniverseparticle::eta) < confEta) && (aod::femtouniverseparticle::pt < confHPtPart1) && (aod::femtouniverseparticle::pt > confLPtPart1);
  Partition<FemtoRecoBasicParticles> partsOneMCrecoBasic = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kTrack)) && (aod::femtouniverseparticle::mAntiLambda == confChargePart1) && (nabs(aod::femtouniverseparticle::eta) < confEta) && (aod::femtouniverseparticle::pt < confHPtPart1) && (aod::femtouniverseparticle::pt > confLPtPart1);

  /// Partition for particle 2 using extended table (cascade)
  Partition<FemtoFullParticles> partsTwoFull = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kCascade)) && (aod::femtouniverseparticle::pt < confHPtPart2) && (aod::femtouniverseparticle::pt > confLPtPart2);

  /// Partition for particle 2 without extended table (cascade)
  Partition<aod::FDParticles> partsTwoBasic = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kCascade)) && (aod::femtouniverseparticle::pt < confHPtPart2) && (aod::femtouniverseparticle::pt > confLPtPart2);
  Partition<aod::FDParticles> partsTwoMCgenBasic = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kMCTruthTrack)) && (aod::femtouniverseparticle::pt < confHPtPart2) && (aod::femtouniverseparticle::pt > confLPtPart2);
  Partition<FemtoRecoBasicParticles> partsTwoMCrecoBasic = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kCascade)) && (aod::femtouniverseparticle::pt < confHPtPart2) && (aod::femtouniverseparticle::pt > confLPtPart2);

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
  FemtoUniverseDetaDphiStar<aod::femtouniverseparticle::ParticleType::kTrack, aod::femtouniverseparticle::ParticleType::kCascade> pairCloseRejection;
  FemtoUniverseDetaDphiStar<aod::femtouniverseparticle::ParticleType::kCascade, aod::femtouniverseparticle::ParticleType::kCascade> pairCloseRejectionCasc;

  HistogramRegistry rXiQA{"xi", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry qaRegistry{"TrackQA", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry resultRegistry{"Correlations", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registryMCgen{"MCgenHistos", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry registryMCreco{"MCrecoHistos", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};

  std::set<int> cascDuplicates;

  std::unique_ptr<TFile> plocalEffFile;
  std::unique_ptr<TH1> pEffHistp1;
  std::unique_ptr<TH1> pEffHistp2;

  Service<o2::ccdb::BasicCCDBManager> ccdb;

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
  bool isParticleTPC(const T& part, int id, float* partSigma = 0)
  {
    const float tpcNSigmas[3] = {aod::pidtpc_tiny::binning::unPackInTable(part.tpcNSigmaStorePr()), aod::pidtpc_tiny::binning::unPackInTable(part.tpcNSigmaStorePi()), aod::pidtpc_tiny::binning::unPackInTable(part.tpcNSigmaStoreKa())};
    if (partSigma)
      *partSigma = tpcNSigmas[id];
    return isNSigmaTPC(tpcNSigmas[id]);
  }

  template <typename T>
  bool isParticleTOF(const T& part, int id, float* partSigma = 0)
  {
    const float tofNSigmas[3] = {aod::pidtof_tiny::binning::unPackInTable(part.tofNSigmaStorePr()), aod::pidtof_tiny::binning::unPackInTable(part.tofNSigmaStorePi()), aod::pidtof_tiny::binning::unPackInTable(part.tofNSigmaStoreKa())};
    if (partSigma)
      *partSigma = tofNSigmas[id];
    return isNSigmaTOF(part.p(), tofNSigmas[id], part.tempFitVar());
  }

  template <typename T>
  bool isParticleCombined(const T& part, int id)
  {
    const float tpcNSigmas[3] = {aod::pidtpc_tiny::binning::unPackInTable(part.tpcNSigmaStorePr()), aod::pidtpc_tiny::binning::unPackInTable(part.tpcNSigmaStorePi()), aod::pidtpc_tiny::binning::unPackInTable(part.tpcNSigmaStoreKa())};
    const float tofNSigmas[3] = {aod::pidtof_tiny::binning::unPackInTable(part.tofNSigmaStorePr()), aod::pidtof_tiny::binning::unPackInTable(part.tofNSigmaStorePi()), aod::pidtof_tiny::binning::unPackInTable(part.tofNSigmaStoreKa())};

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
    /// nSigma debug histograms for the selected particle species only i.e. not sigmas of all particles mixed together
    qaRegistry.add("Tracks_pos/nSigmaTPC", "; #it{p} (GeV/#it{c}); n#sigma_{TPC}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("Tracks_pos/nSigmaTOF", "; #it{p} (GeV/#it{c}); n#sigma_{TOF}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("Tracks_neg/nSigmaTPC", "; #it{p} (GeV/#it{c}); n#sigma_{TPC}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("Tracks_neg/nSigmaTOF", "; #it{p} (GeV/#it{c}); n#sigma_{TOF}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("MixingQA/hMECollisionBins", ";bin;Entries", kTH1F, {{120, -0.5, 119.5}});

    qaRegistry.add("V0Child_pos/nSigmaTPC", "; #it{p} (GeV/#it{c}); n#sigma_{TPC}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("V0Child_neg/nSigmaTPC", "; #it{p} (GeV/#it{c}); n#sigma_{TPC}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("hBachelor/nSigmaTPC", "; #it{p} (GeV/#it{c}); n#sigma_{TPC}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("V0Child_pos/nSigmaTOF", "; #it{p} (GeV/#it{c}); n#sigma_{TOF}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("V0Child_neg/nSigmaTOF", "; #it{p} (GeV/#it{c}); n#sigma_{TOF}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("hBachelor/nSigmaTOF", "; #it{p} (GeV/#it{c}); n#sigma_{TOF}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});

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
    sameEventCont.setPDGCodes(confTrkPDGCodePartOne, confCascPDGCodePartTwo);
    mixedEventCont.setPDGCodes(confTrkPDGCodePartOne, confCascPDGCodePartTwo);

    pairCleaner.init(&qaRegistry);
    pairCleanerCasc.init(&qaRegistry);
    if (confIsCPR.value) {
      if (doprocessSameEvent || doprocessSameEventBitmask || doprocessMixedEvent || doprocessMixedEventBitmask)
        pairCloseRejection.init(&resultRegistry, &qaRegistry, confCPRdeltaPhiCutMin.value, confCPRdeltaPhiCutMax.value, confCPRdeltaEtaCutMin.value, confCPRdeltaEtaCutMax.value, confCPRChosenRadii.value, confCPRPlotPerRadii.value, 0, 0, confIsSameSignCPR.value);
      if (doprocessSameEventCasc || doprocessSameEventCascBitmask || doprocessMixedEventCasc || doprocessMixedEventCascBitmask)
        pairCloseRejectionCasc.init(&resultRegistry, &qaRegistry, confCPRdeltaPhiCutMin.value, confCPRdeltaPhiCutMax.value, confCPRdeltaEtaCutMin.value, confCPRdeltaEtaCutMax.value, confCPRChosenRadii.value, confCPRPlotPerRadii.value, 0, 0, confIsSameSignCPR.value);
    }

    if (!confLocalEfficiency.value.empty()) {
      plocalEffFile = std::unique_ptr<TFile>(TFile::Open(confLocalEfficiency.value.c_str(), "read"));
      if (!plocalEffFile || plocalEffFile.get()->IsZombie())
        LOGF(fatal, "Could not load efficiency histogram from %s", confLocalEfficiency.value.c_str());
      if (doprocessSameEvent || doprocessSameEventBitmask || doprocessMixedEvent || doprocessMixedEventBitmask) {
        pEffHistp1 = (confChargePart1 > 0) ? std::unique_ptr<TH1>(plocalEffFile.get()->Get<TH1>("PrPlus")) : std::unique_ptr<TH1>(plocalEffFile.get()->Get<TH1>("PrMinus")); // note: works only for protons for now
        pEffHistp2 = (confCascType1 == 0 || confCascType1 == 1) ? std::unique_ptr<TH1>(plocalEffFile.get()->Get<TH1>("Cascade")) : std::unique_ptr<TH1>(plocalEffFile.get()->Get<TH1>("AntiCascade"));
        LOGF(info, "Loaded efficiency histograms for track-Cascade.");
      } else if (doprocessSameEventCasc || doprocessSameEventCascBitmask || doprocessMixedEventCasc || doprocessMixedEventCascBitmask) {
        pEffHistp1 = (confCascType1 == 0 || confCascType1 == 1) ? std::unique_ptr<TH1>(plocalEffFile.get()->Get<TH1>("Cascade")) : std::unique_ptr<TH1>(plocalEffFile.get()->Get<TH1>("AntiCascade"));
        pEffHistp2 = (confCascType2 == 0 || confCascType2 == 1) ? std::unique_ptr<TH1>(plocalEffFile.get()->Get<TH1>("Cascade")) : std::unique_ptr<TH1>(plocalEffFile.get()->Get<TH1>("AntiCascade"));
        LOGF(info, "Loaded efficiency histograms for Cascade-Cascade.");
      }
    } else if (!confCCDBEfficiency.value.empty()) {
      ccdb->setURL("http://alice-ccdb.cern.ch");
      ccdb->setCaching(true);
      ccdb->setLocalObjectValidityChecking();

      auto now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
      ccdb->setCreatedNotAfter(now);
      if (doprocessSameEvent || doprocessSameEventBitmask || doprocessMixedEvent || doprocessMixedEventBitmask) {
        pEffHistp1 = (confChargePart1 > 0) ? std::unique_ptr<TH1>(ccdb->getForTimeStamp<TH1>(confCCDBEfficiency.value + "/PrPlus", now)) : std::unique_ptr<TH1>(ccdb->getForTimeStamp<TH1>(confCCDBEfficiency.value + "/PrMinus", now)); // note: works only for protons for now
        pEffHistp2 = (confCascType1 == 0 || confCascType1 == 1) ? std::unique_ptr<TH1>(ccdb->getForTimeStamp<TH1>(confCCDBEfficiency.value + "/Cascade", now)) : std::unique_ptr<TH1>(ccdb->getForTimeStamp<TH1>(confCCDBEfficiency.value + "/AntiCascade", now));
        LOGF(info, "Loaded efficiency histograms for track-Cascade from CCDB");
      } else if (doprocessSameEventCasc || doprocessSameEventCascBitmask || doprocessMixedEventCasc || doprocessMixedEventCascBitmask) {
        pEffHistp1 = (confCascType1 == 0 || confCascType1 == 1) ? std::unique_ptr<TH1>(ccdb->getForTimeStamp<TH1>(confCCDBEfficiency.value + "/Cascade", now)) : std::unique_ptr<TH1>(ccdb->getForTimeStamp<TH1>(confCCDBEfficiency.value + "/AntiCascade", now));
        pEffHistp2 = (confCascType2 == 0 || confCascType2 == 1) ? std::unique_ptr<TH1>(ccdb->getForTimeStamp<TH1>(confCCDBEfficiency.value + "/Cascade", now)) : std::unique_ptr<TH1>(ccdb->getForTimeStamp<TH1>(confCCDBEfficiency.value + "/AntiCascade", now));
        LOGF(info, "Loaded efficiency histograms for Cascade-Cascade from CCDB.");
      }
    }
  }

  void processCascadeQA([[maybe_unused]] const FilteredFDCollision& col, const FemtoFullParticles& parts, const aod::FDCascParticles& fdcascs)
  {
    for (const auto& casc : fdcascs) {
      const auto& part = casc.fdParticle_as<FemtoFullParticles>();
      rXiQA.fill(HIST("hMassXi"), part.mLambda());

      const auto& posChild = parts.iteratorAt(part.globalIndex() - 3 - parts.begin().globalIndex());
      const auto& negChild = parts.iteratorAt(part.globalIndex() - 2 - parts.begin().globalIndex());
      const auto& bachelor = parts.iteratorAt(part.globalIndex() - 1 - parts.begin().globalIndex());

      float posChildTPC, negChildTPC, bachelorTPC, posChildTOF, negChildTOF, bachelorTOF;
      if (!isParticleTPC(posChild, CascChildTable[confCascType1][0], &posChildTPC) || !isParticleTPC(negChild, CascChildTable[confCascType1][1], &negChildTPC) || !isParticleTPC(bachelor, CascChildTable[confCascType1][2], &bachelorTPC))
        continue;

      if (!isParticleTOF(posChild, CascChildTable[confCascType1][0], &posChildTOF) || !isParticleTOF(negChild, CascChildTable[confCascType1][1], &negChildTOF) || !isParticleTOF(bachelor, CascChildTable[confCascType1][2], &bachelorTOF))
        continue;

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
    }
  }
  PROCESS_SWITCH(femtoUniversePairTaskTrackCascadeExtended, processCascadeQA, "Enable processing cascades", false);

  template <class T>
  using hasSigma = decltype(std::declval<T&>().tpcNSigmaStorePr());

  /// track - cascade correlations
  template <class TableType, typename PartitionType>
  void doSameEvent(const FilteredFDCollision& col, const TableType& parts, PartitionType& partsOne, PartitionType& partsTwo)
  {
    const auto& magFieldTesla = col.magField();

    auto groupPartsOne = partsOne->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto groupPartsTwo = partsTwo->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);

    eventHisto.fillQA(col);

    const int multCol = confUseCent ? col.multV0M() : col.multNtr();

    for (const auto& part : groupPartsTwo) {
      if (!invMCascade(part.mLambda(), part.mAntiLambda(), confCascType1)) /// mLambda stores Xi mass, mAntiLambda stores Omega mass
        continue;

      const auto& posChild = parts.iteratorAt(part.globalIndex() - 3 - parts.begin().globalIndex());
      const auto& negChild = parts.iteratorAt(part.globalIndex() - 2 - parts.begin().globalIndex());
      const auto& bachelor = parts.iteratorAt(part.globalIndex() - 1 - parts.begin().globalIndex());
      /// Child particles must pass this condition to be selected
      if constexpr (std::experimental::is_detected<hasSigma, typename TableType::iterator>::value) {
        float posChildTPC, negChildTPC, bachelorTPC, posChildTOF, negChildTOF, bachelorTOF;
        if (!isParticleTPC(posChild, CascChildTable[confCascType1][0], &posChildTPC) || !isParticleTPC(negChild, CascChildTable[confCascType1][1], &negChildTPC) || !isParticleTPC(bachelor, CascChildTable[confCascType1][2], &bachelorTPC))
          continue;

        if (!isParticleTOF(posChild, CascChildTable[confCascType1][0], &posChildTOF) || !isParticleTOF(negChild, CascChildTable[confCascType1][1], &negChildTOF) || !isParticleTOF(bachelor, CascChildTable[confCascType1][2], &bachelorTOF))
          continue;

        posChildHistos.fillQA<false, true>(posChild);
        negChildHistos.fillQA<false, true>(negChild);
        bachHistos.fillQABase<false, true>(bachelor, HIST("hBachelor"));

        qaRegistry.fill(HIST("V0Child_pos/nSigmaTPC"), posChild.p(), posChildTPC);
        qaRegistry.fill(HIST("V0Child_neg/nSigmaTPC"), negChild.p(), negChildTPC);
        qaRegistry.fill(HIST("hBachelor/nSigmaTPC"), bachelor.p(), bachelorTPC);
        qaRegistry.fill(HIST("V0Child_pos/nSigmaTOF"), posChild.p(), posChildTOF);
        qaRegistry.fill(HIST("V0Child_neg/nSigmaTOF"), negChild.p(), negChildTOF);
        qaRegistry.fill(HIST("hBachelor/nSigmaTOF"), bachelor.p(), bachelorTOF);

        cascQAHistos.fillQA<false, true>(part);

      } else {
        if ((posChild.pidCut() & (1u << CascChildTable[confCascType1][0])) == 0 || (negChild.pidCut() & (1u << CascChildTable[confCascType1][1])) == 0 || (bachelor.pidCut() & (1u << CascChildTable[confCascType1][2])) == 0)
          continue;

        if (confUseStrangenessTOF) {
          if (((confCascType1 == 1 || confCascType1 == 3) && (part.pidCut() & 7) != 7) || ((confCascType1 == 0 || confCascType1 == 2) && (part.pidCut() & 56) != 56))
            continue;
        } else {
          if ((posChild.pidCut() & (8u << CascChildTable[confCascType1][0])) == 0 || (negChild.pidCut() & (8u << CascChildTable[confCascType1][1])) == 0 || (bachelor.pidCut() & (8u << CascChildTable[confCascType1][2])) == 0)
            continue;
        }

        posChildHistos.fillQA<false, false>(posChild);
        negChildHistos.fillQA<false, false>(negChild);
        bachHistos.fillQABase<false, false>(bachelor, HIST("hBachelor"));
        cascQAHistos.fillQA<false, false>(part);
      }
      rXiQA.fill(HIST("hInvMpTmult"), part.pt(), part.mLambda(), multCol);
    }

    if constexpr (std::experimental::is_detected<hasSigma, typename TableType::iterator>::value) {
      for (const auto& part : groupPartsOne) {
        /// PID plot for track particle
        const float tpcNSigmas[3] = {aod::pidtpc_tiny::binning::unPackInTable(part.tpcNSigmaStorePr()), aod::pidtpc_tiny::binning::unPackInTable(part.tpcNSigmaStorePi()), aod::pidtpc_tiny::binning::unPackInTable(part.tpcNSigmaStoreKa())};
        const float tofNSigmas[3] = {aod::pidtof_tiny::binning::unPackInTable(part.tofNSigmaStorePr()), aod::pidtof_tiny::binning::unPackInTable(part.tofNSigmaStorePi()), aod::pidtof_tiny::binning::unPackInTable(part.tofNSigmaStoreKa())};

        if (!isNSigmaCombined(part.p(), tpcNSigmas[confTrackChoicePartOne], tofNSigmas[confTrackChoicePartOne]))
          continue;

        if (part.mAntiLambda() > 0) {
          qaRegistry.fill(HIST("Tracks_pos/nSigmaTPC"), part.p(), tpcNSigmas[confTrackChoicePartOne]);
          qaRegistry.fill(HIST("Tracks_pos/nSigmaTOF"), part.p(), tofNSigmas[confTrackChoicePartOne]);
          trackHistoPartOnePos.fillQA<false, false>(part);
        } else if (part.mAntiLambda() < 0) {
          qaRegistry.fill(HIST("Tracks_neg/nSigmaTPC"), part.p(), tpcNSigmas[confTrackChoicePartOne]);
          qaRegistry.fill(HIST("Tracks_neg/nSigmaTOF"), part.p(), tofNSigmas[confTrackChoicePartOne]);
          trackHistoPartOneNeg.fillQA<false, false>(part);
        }
      }
    }
    for (const auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(groupPartsOne, groupPartsTwo))) {
      // Cascade inv mass cut (mLambda stores Xi mass, mAntiLambda stores Omega mass)
      if (!invMCascade(p2.mLambda(), p2.mAntiLambda(), confCascType1))
        continue;
      // PID
      if constexpr (std::experimental::is_detected<hasSigma, typename TableType::iterator>::value) {
        if (!isParticleCombined(p1, confTrackChoicePartOne))
          continue;
      } else {
        if ((p1.pidCut() & (64u << confTrackChoicePartOne)) == 0)
          continue;
      }
      // track cleaning
      if (!pairCleaner.isCleanPair(p1, p2, parts)) {
        continue;
      }
      if (confIsCPR.value) {
        if (pairCloseRejection.isClosePair(p1, p2, parts, magFieldTesla, femto_universe_container::EventType::same)) {
          continue;
        }
      }

      const auto& posChild = parts.iteratorAt(p2.globalIndex() - 3 - parts.begin().globalIndex());
      const auto& negChild = parts.iteratorAt(p2.globalIndex() - 2 - parts.begin().globalIndex());
      const auto& bachelor = parts.iteratorAt(p2.globalIndex() - 1 - parts.begin().globalIndex());
      /// Child particles must pass this condition to be selected
      if constexpr (std::experimental::is_detected<hasSigma, typename TableType::iterator>::value) {
        if (!isParticleTPC(posChild, CascChildTable[confCascType1][0]) || !isParticleTPC(negChild, CascChildTable[confCascType1][1]) || !isParticleTPC(bachelor, CascChildTable[confCascType1][2]))
          continue;
        if (!isParticleTOF(posChild, CascChildTable[confCascType1][0]) || !isParticleTOF(negChild, CascChildTable[confCascType1][1]) || !isParticleTOF(bachelor, CascChildTable[confCascType1][2]))
          continue;
      } else {
        if ((posChild.pidCut() & (1u << CascChildTable[confCascType1][0])) == 0 || (negChild.pidCut() & (1u << CascChildTable[confCascType1][1])) == 0 || (bachelor.pidCut() & (1u << CascChildTable[confCascType1][2])) == 0)
          continue;
        if (confUseStrangenessTOF) {
          if (((confCascType1 == 1 || confCascType1 == 3) && (p2.pidCut() & 7) != 7) || ((confCascType1 == 0 || confCascType1 == 2) && (p2.pidCut() & 56) != 56))
            continue;
        } else {
          if ((posChild.pidCut() & (8u << CascChildTable[confCascType1][0])) == 0 || (negChild.pidCut() & (8u << CascChildTable[confCascType1][1])) == 0 || (bachelor.pidCut() & (8u << CascChildTable[confCascType1][2])) == 0)
            continue;
        }
      }
      float weight = 1.0f;
      if (pEffHistp1)
        weight = pEffHistp1.get()->GetBinContent(pEffHistp1->FindBin(p1.pt(), p1.eta())) * pEffHistp2.get()->GetBinContent(pEffHistp2->FindBin(p2.pt(), p2.eta()));
      sameEventCont.setPair<false>(p1, p2, multCol, confUse3D, weight);
    }
  }

  void processSameEvent(const FilteredFDCollision& col, const FemtoFullParticles& parts)
  {
    doSameEvent(col, parts, partsOneFull, partsTwoFull);
  }
  PROCESS_SWITCH(femtoUniversePairTaskTrackCascadeExtended, processSameEvent, "Enable processing same event for track - cascade", false);

  void processSameEventBitmask(const FilteredFDCollision& col, const aod::FDParticles& parts)
  {
    doSameEvent(col, parts, partsOneBasic, partsTwoBasic);
  }
  PROCESS_SWITCH(femtoUniversePairTaskTrackCascadeExtended, processSameEventBitmask, "Enable processing same event for track - cascade using bitmask for PID", false);

  /// cascade - cascade correlations
  Preslice<aod::FDCascParticles> perFDPartsCasc = aod::femtouniversecascparticle::fdParticleId;

  template <class TableType, typename PartitionType>
  void doSameEventCasc(const FilteredFDCollision& col, const TableType& parts, PartitionType& partsTwo, const aod::FDCascParticles& cascs)
  {
    const auto& magFieldTesla = col.magField();

    auto groupPartsTwo = partsTwo->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);

    eventHisto.fillQA(col);

    const int multCol = confUseCent ? col.multV0M() : col.multNtr();

    for (const auto& part : groupPartsTwo) {
      if (!invMCascade(part.mLambda(), part.mAntiLambda(), confCascType1)) /// mLambda stores Xi mass, mAntiLambda stores Omega mass
        continue;

      if constexpr (std::experimental::is_detected<hasSigma, typename TableType::iterator>::value)
        cascQAHistos.fillQA<false, true>(part);
      else
        cascQAHistos.fillQA<false, false>(part);

      const auto& posChild = parts.iteratorAt(part.globalIndex() - 3 - parts.begin().globalIndex());
      const auto& negChild = parts.iteratorAt(part.globalIndex() - 2 - parts.begin().globalIndex());
      const auto& bachelor = parts.iteratorAt(part.globalIndex() - 1 - parts.begin().globalIndex());
      /// Check daughters of first cascade
      if constexpr (std::experimental::is_detected<hasSigma, typename TableType::iterator>::value) {
        if (!isParticleTPC(posChild, CascChildTable[confCascType1][0]) || !isParticleTPC(negChild, CascChildTable[confCascType1][1]) || !isParticleTPC(bachelor, CascChildTable[confCascType1][2]))
          continue;
        if (!isParticleTOF(posChild, CascChildTable[confCascType1][0]) || !isParticleTOF(negChild, CascChildTable[confCascType1][1]) || !isParticleTOF(bachelor, CascChildTable[confCascType1][2]))
          continue;

        posChildHistos.fillQA<false, true>(posChild);
        negChildHistos.fillQA<false, true>(negChild);
        bachHistos.fillQABase<false, true>(bachelor, HIST("hBachelor"));
      } else {
        if ((posChild.pidCut() & (1u << CascChildTable[confCascType1][0])) == 0 || (negChild.pidCut() & (1u << CascChildTable[confCascType1][1])) == 0 || (bachelor.pidCut() & (1u << CascChildTable[confCascType1][2])) == 0)
          continue;
        if (confUseStrangenessTOF) {
          if (((confCascType1 == 1 || confCascType1 == 3) && (part.pidCut() & 7) != 7) || ((confCascType1 == 0 || confCascType1 == 2) && (part.pidCut() & 56) != 56))
            continue;
        } else {
          if ((posChild.pidCut() & (8u << CascChildTable[confCascType1][0])) == 0 || (negChild.pidCut() & (8u << CascChildTable[confCascType1][1])) == 0 || (bachelor.pidCut() & (8u << CascChildTable[confCascType1][2])) == 0)
            continue;
        }

        posChildHistos.fillQA<false, false>(posChild);
        negChildHistos.fillQA<false, false>(negChild);
        bachHistos.fillQABase<false, false>(bachelor, HIST("hBachelor"));
      }
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
        // mark for rejection the cascade that share a daughter with the other cascade and has a better cosPA value
        auto groupedCasc1 = cascs.sliceBy(perFDPartsCasc, p1.globalIndex());
        auto groupedCasc2 = cascs.sliceBy(perFDPartsCasc, p2.globalIndex());
        if (groupedCasc1.size() <= 0 || groupedCasc2.size() <= 0) {
          LOGF(warning, "Either cascade1 (%u) or cascade2 (%u) list is empty", groupedCasc1.size(), groupedCasc2.size()); // this should never happen but just for a sanity check
          return;
        }
        if (std::abs(groupedCasc1.begin().cpaCasc() - 1) < std::abs(groupedCasc2.begin().cpaCasc() - 1)) {
          cascDuplicates.insert(p1.globalIndex());
        } else {
          cascDuplicates.insert(p2.globalIndex());
        }
      }
    };

    auto pairProcessFunc = [&](auto& p1, auto& p2) -> bool {
      if (cascDuplicates.contains(p1.globalIndex()) || cascDuplicates.contains(p2.globalIndex()))
        return false;
      if (!invMCascade(p1.mLambda(), p1.mAntiLambda(), confCascType1))
        return false;
      if (!invMCascade(p2.mLambda(), p2.mAntiLambda(), confCascType2))
        return false;

      const auto& posChild1 = parts.iteratorAt(p1.globalIndex() - 3 - parts.begin().globalIndex());
      const auto& negChild1 = parts.iteratorAt(p1.globalIndex() - 2 - parts.begin().globalIndex());
      const auto& bachelor1 = parts.iteratorAt(p1.globalIndex() - 1 - parts.begin().globalIndex());
      /// Child particles must pass this condition to be selected
      if constexpr (std::experimental::is_detected<hasSigma, typename TableType::iterator>::value) {
        if (!isParticleTPC(posChild1, CascChildTable[confCascType1][0]) || !isParticleTPC(negChild1, CascChildTable[confCascType1][1]) || !isParticleTPC(bachelor1, CascChildTable[confCascType1][2]))
          return false;
        if (!isParticleTOF(posChild1, CascChildTable[confCascType1][0]) || !isParticleTOF(negChild1, CascChildTable[confCascType1][1]) || !isParticleTOF(bachelor1, CascChildTable[confCascType1][2]))
          return false;
      } else {
        if ((posChild1.pidCut() & (1u << CascChildTable[confCascType1][0])) == 0 || (negChild1.pidCut() & (1u << CascChildTable[confCascType1][1])) == 0 || (bachelor1.pidCut() & (1u << CascChildTable[confCascType1][2])) == 0)
          return false;
        if (confUseStrangenessTOF) {
          if (((confCascType1 == 1 || confCascType1 == 3) && (p1.pidCut() & 7) != 7) || ((confCascType1 == 0 || confCascType1 == 2) && (p1.pidCut() & 56) != 56))
            return false;
        } else {
          if ((posChild1.pidCut() & (8u << CascChildTable[confCascType1][0])) == 0 || (negChild1.pidCut() & (8u << CascChildTable[confCascType1][1])) == 0 || (bachelor1.pidCut() & (8u << CascChildTable[confCascType1][2])) == 0)
            return false;
        }
      }

      const auto& posChild2 = parts.iteratorAt(p2.globalIndex() - 3 - parts.begin().globalIndex());
      const auto& negChild2 = parts.iteratorAt(p2.globalIndex() - 2 - parts.begin().globalIndex());
      const auto& bachelor2 = parts.iteratorAt(p2.globalIndex() - 1 - parts.begin().globalIndex());
      /// Child particles must pass this condition to be selected
      if constexpr (std::experimental::is_detected<hasSigma, typename TableType::iterator>::value) {
        if (!isParticleTPC(posChild2, CascChildTable[confCascType2][0]) || !isParticleTPC(negChild2, CascChildTable[confCascType2][1]) || !isParticleTPC(bachelor2, CascChildTable[confCascType2][2]))
          return false;
        if (!isParticleTOF(posChild2, CascChildTable[confCascType2][0]) || !isParticleTOF(negChild2, CascChildTable[confCascType2][1]) || !isParticleTOF(bachelor2, CascChildTable[confCascType2][2]))
          return false;
      } else {
        if ((posChild2.pidCut() & (1u << CascChildTable[confCascType2][0])) == 0 || (negChild2.pidCut() & (1u << CascChildTable[confCascType2][1])) == 0 || (bachelor2.pidCut() & (1u << CascChildTable[confCascType2][2])) == 0)
          return false;
        if (confUseStrangenessTOF) {
          if (((confCascType2 == 1 || confCascType2 == 3) && (p2.pidCut() & 7) != 7) || ((confCascType2 == 0 || confCascType2 == 2) && (p2.pidCut() & 56) != 56))
            return false;
        } else {
          if ((posChild2.pidCut() & (8u << CascChildTable[confCascType2][0])) == 0 || (negChild2.pidCut() & (8u << CascChildTable[confCascType2][1])) == 0 || (bachelor2.pidCut() & (8u << CascChildTable[confCascType2][2])) == 0)
            return false;
        }
      }

      if (confIsCPR.value) {
        if (pairCloseRejectionCasc.isClosePair(p1, p2, parts, magFieldTesla, femto_universe_container::EventType::same)) {
          return false;
        }
      }

      float weight = 1.0f;
      if (pEffHistp1)
        weight = pEffHistp1.get()->GetBinContent(pEffHistp1->FindBin(p1.pt(), p1.eta())) * pEffHistp2.get()->GetBinContent(pEffHistp2->FindBin(p2.pt(), p2.eta()));
      sameEventCont.setPair<false>(p1, p2, multCol, confUse3D, weight);
      return true;
    };

    cascDuplicates.clear();
    for (const auto& [p1, p2] : combinations(CombinationsStrictlyUpperIndexPolicy(groupPartsTwo, groupPartsTwo))) {
      pairDuplicateCheckFunc(p1, p2);
    }
    /// Now build the combinations for cascades
    for (const auto& [p1, p2] : combinations(CombinationsStrictlyUpperIndexPolicy(groupPartsTwo, groupPartsTwo))) {
      if (!pairProcessFunc(p1, p2))
        pairProcessFunc(p2, p1);
    }
  }

  void processSameEventCasc(const FilteredFDCollision& col, const FemtoFullParticles& parts, const aod::FDCascParticles& cascs)
  {
    doSameEventCasc(col, parts, partsTwoFull, cascs);
  }
  PROCESS_SWITCH(femtoUniversePairTaskTrackCascadeExtended, processSameEventCasc, "Enable processing same event for cascade - cascade", false);

  void processSameEventCascBitmask(const FilteredFDCollision& col, const aod::FDParticles& parts, const aod::FDCascParticles& cascs)
  {
    doSameEventCasc(col, parts, partsTwoBasic, cascs);
  }
  PROCESS_SWITCH(femtoUniversePairTaskTrackCascadeExtended, processSameEventCascBitmask, "Enable processing same event for cascade - cascade using bitmask for PID", false);

  /// track - cascade correlations
  template <class TableType, typename PartitionType>
  void doMixedEvent(const FilteredFDCollisions& cols, const TableType& parts, PartitionType& partsOne, PartitionType& partsTwo)
  {
    ColumnBinningPolicy<aod::collision::PosZ, aod::femtouniversecollision::MultNtr> colBinningMult{{confVtxBins, confMultBins}, true};
    ColumnBinningPolicy<aod::collision::PosZ, aod::femtouniversecollision::MultV0M> colBinningCent{{confVtxBins, confMultBins}, true};

    auto mixedCollProcessFunc = [&](auto& collision1, auto& collision2) -> void {
      const int multCol = confUseCent ? collision1.multV0M() : collision1.multNtr();

      auto groupPartsOne = partsOne->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);
      auto groupPartsTwo = partsTwo->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision2.globalIndex(), cache);

      const auto& magFieldTesla1 = collision1.magField();
      const auto& magFieldTesla2 = collision2.magField();

      if (magFieldTesla1 != magFieldTesla2) {
        return;
      }
      for (const auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(groupPartsOne, groupPartsTwo))) {
        // Cascade inv mass cut (mLambda stores Xi mass, mAntiLambda stores Omega mass)
        if (!invMCascade(p2.mLambda(), p2.mAntiLambda(), confCascType1))
          continue;
        // PID
        if constexpr (std::experimental::is_detected<hasSigma, typename TableType::iterator>::value) {
          if (!isParticleCombined(p1, confTrackChoicePartOne))
            continue;
        } else {
          if ((p1.pidCut() & (64u << confTrackChoicePartOne)) == 0)
            continue;
        }

        const auto& posChild = parts.iteratorAt(p2.globalIndex() - 3 - parts.begin().globalIndex());
        const auto& negChild = parts.iteratorAt(p2.globalIndex() - 2 - parts.begin().globalIndex());
        const auto& bachelor = parts.iteratorAt(p2.globalIndex() - 1 - parts.begin().globalIndex());
        /// Child particles must pass this condition to be selected
        if constexpr (std::experimental::is_detected<hasSigma, typename TableType::iterator>::value) {
          if (!isParticleTPC(posChild, CascChildTable[confCascType1][0]) || !isParticleTPC(negChild, CascChildTable[confCascType1][1]) || !isParticleTPC(bachelor, CascChildTable[confCascType1][2]))
            continue;
          if (!isParticleTOF(posChild, CascChildTable[confCascType1][0]) || !isParticleTOF(negChild, CascChildTable[confCascType1][1]) || !isParticleTOF(bachelor, CascChildTable[confCascType1][2]))
            continue;
        } else {
          if ((posChild.pidCut() & (1u << CascChildTable[confCascType1][0])) == 0 || (negChild.pidCut() & (1u << CascChildTable[confCascType1][1])) == 0 || (bachelor.pidCut() & (1u << CascChildTable[confCascType1][2])) == 0)
            continue;
          if (confUseStrangenessTOF) {
            if (((confCascType1 == 1 || confCascType1 == 3) && (p2.pidCut() & 7) != 7) || ((confCascType1 == 0 || confCascType1 == 2) && (p2.pidCut() & 56) != 56))
              continue;
          } else {
            if ((posChild.pidCut() & (8u << CascChildTable[confCascType1][0])) == 0 || (negChild.pidCut() & (8u << CascChildTable[confCascType1][1])) == 0 || (bachelor.pidCut() & (8u << CascChildTable[confCascType1][2])) == 0)
              continue;
          }
        }

        // track cleaning
        if (!pairCleaner.isCleanPair(p1, p2, parts)) {
          continue;
        }
        if (confIsCPR.value) {
          if (pairCloseRejection.isClosePair(p1, p2, parts, magFieldTesla1, femto_universe_container::EventType::mixed)) {
            continue;
          }
        }

        float weight = 1.0f;
        if (pEffHistp1)
          weight = pEffHistp1.get()->GetBinContent(pEffHistp1->FindBin(p1.pt(), p1.eta())) * pEffHistp2.get()->GetBinContent(pEffHistp2->FindBin(p2.pt(), p2.eta()));
        mixedEventCont.setPair<false>(p1, p2, multCol, confUse3D, weight);
      }
    };

    if (confUseCent) {
      for (const auto& [collision1, collision2] : soa::selfCombinations(colBinningCent, confNEventsMix, -1, cols, cols)) {
        mixedCollProcessFunc(collision1, collision2);
        qaRegistry.fill(HIST("MixingQA/hMECollisionBins"), colBinningCent.getBin({collision1.posZ(), collision1.multV0M()}));
      }
    } else {
      for (const auto& [collision1, collision2] : soa::selfCombinations(colBinningMult, confNEventsMix, -1, cols, cols)) {
        mixedCollProcessFunc(collision1, collision2);
        qaRegistry.fill(HIST("MixingQA/hMECollisionBins"), colBinningMult.getBin({collision1.posZ(), collision1.multNtr()}));
      }
    }
  }

  void processMixedEvent(const FilteredFDCollisions& cols, const FemtoFullParticles& parts)
  {
    doMixedEvent(cols, parts, partsOneFull, partsTwoFull);
  }
  PROCESS_SWITCH(femtoUniversePairTaskTrackCascadeExtended, processMixedEvent, "Enable processing mixed event for track - cascade", false);

  void processMixedEventBitmask(const FilteredFDCollisions& cols, const aod::FDParticles& parts)
  {
    doMixedEvent(cols, parts, partsOneBasic, partsTwoBasic);
  }
  PROCESS_SWITCH(femtoUniversePairTaskTrackCascadeExtended, processMixedEventBitmask, "Enable processing mixed event for track - cascade using bitmask for PID", false);

  /// cascade - cascade correlations
  template <class TableType, typename PartitionType>
  void doMixedEventCasc(const FilteredFDCollisions& cols, const TableType& parts, PartitionType& partsTwo)
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

        const auto& posChild1 = parts.iteratorAt(p1.globalIndex() - 3 - parts.begin().globalIndex());
        const auto& negChild1 = parts.iteratorAt(p1.globalIndex() - 2 - parts.begin().globalIndex());
        const auto& bachelor1 = parts.iteratorAt(p1.globalIndex() - 1 - parts.begin().globalIndex());
        /// Child particles must pass this condition to be selected
        if constexpr (std::experimental::is_detected<hasSigma, typename TableType::iterator>::value) {
          if (!isParticleTPC(posChild1, CascChildTable[confCascType1][0]) || !isParticleTPC(negChild1, CascChildTable[confCascType1][1]) || !isParticleTPC(bachelor1, CascChildTable[confCascType1][2]))
            continue;
          if (!isParticleTOF(posChild1, CascChildTable[confCascType1][0]) || !isParticleTOF(negChild1, CascChildTable[confCascType1][1]) || !isParticleTOF(bachelor1, CascChildTable[confCascType1][2]))
            continue;
        } else {
          if ((posChild1.pidCut() & (1u << CascChildTable[confCascType1][0])) == 0 || (negChild1.pidCut() & (1u << CascChildTable[confCascType1][1])) == 0 || (bachelor1.pidCut() & (1u << CascChildTable[confCascType1][2])) == 0)
            continue;
          if (confUseStrangenessTOF) {
            if (((confCascType1 == 1 || confCascType1 == 3) && (p1.pidCut() & 7) != 7) || ((confCascType1 == 0 || confCascType1 == 2) && (p1.pidCut() & 56) != 56))
              continue;
          } else {
            if ((posChild1.pidCut() & (8u << CascChildTable[confCascType1][0])) == 0 || (negChild1.pidCut() & (8u << CascChildTable[confCascType1][1])) == 0 || (bachelor1.pidCut() & (8u << CascChildTable[confCascType1][2])) == 0)
              continue;
          }
        }

        const auto& posChild2 = parts.iteratorAt(p2.globalIndex() - 3 - parts.begin().globalIndex());
        const auto& negChild2 = parts.iteratorAt(p2.globalIndex() - 2 - parts.begin().globalIndex());
        const auto& bachelor2 = parts.iteratorAt(p2.globalIndex() - 1 - parts.begin().globalIndex());
        /// Child particles must pass this condition to be selected
        if constexpr (std::experimental::is_detected<hasSigma, typename TableType::iterator>::value) {
          if (!isParticleTPC(posChild2, CascChildTable[confCascType2][0]) || !isParticleTPC(negChild2, CascChildTable[confCascType2][1]) || !isParticleTPC(bachelor2, CascChildTable[confCascType2][2]))
            continue;
          if (!isParticleTOF(posChild2, CascChildTable[confCascType2][0]) || !isParticleTOF(negChild2, CascChildTable[confCascType2][1]) || !isParticleTOF(bachelor2, CascChildTable[confCascType2][2]))
            continue;
        } else {
          if ((posChild2.pidCut() & (1u << CascChildTable[confCascType2][0])) == 0 || (negChild2.pidCut() & (1u << CascChildTable[confCascType2][1])) == 0 || (bachelor2.pidCut() & (1u << CascChildTable[confCascType2][2])) == 0)
            continue;
          if (confUseStrangenessTOF) {
            if (((confCascType2 == 1 || confCascType2 == 3) && (p2.pidCut() & 7) != 7) || ((confCascType2 == 0 || confCascType2 == 2) && (p2.pidCut() & 56) != 56))
              continue;
          } else {
            if ((posChild2.pidCut() & (8u << CascChildTable[confCascType2][0])) == 0 || (negChild2.pidCut() & (8u << CascChildTable[confCascType2][1])) == 0 || (bachelor2.pidCut() & (8u << CascChildTable[confCascType2][2])) == 0)
              continue;
          }
        }
        // track cleaning
        if (!pairCleanerCasc.isCleanPair(p1, p2, parts)) {
          continue;
        }
        if (confIsCPR.value) {
          if (pairCloseRejectionCasc.isClosePair(p1, p2, parts, magFieldTesla1, femto_universe_container::EventType::mixed)) {
            continue;
          }
        }

        float weight = 1.0f;
        if (pEffHistp1)
          weight = pEffHistp1.get()->GetBinContent(pEffHistp1->FindBin(p1.pt(), p1.eta())) * pEffHistp2.get()->GetBinContent(pEffHistp2->FindBin(p2.pt(), p2.eta()));
        mixedEventCont.setPair<false>(p1, p2, multCol, confUse3D, weight);
      }
    }
  }

  void processMixedEventCasc(const FilteredFDCollisions& cols, const FemtoFullParticles& parts)
  {
    doMixedEventCasc(cols, parts, partsTwoFull);
  }
  PROCESS_SWITCH(femtoUniversePairTaskTrackCascadeExtended, processMixedEventCasc, "Enable processing mixed event for cascade - cascade", false);

  void processMixedEventCascBitmask(const FilteredFDCollisions& cols, const aod::FDParticles& parts)
  {
    doMixedEventCasc(cols, parts, partsTwoBasic);
  }
  PROCESS_SWITCH(femtoUniversePairTaskTrackCascadeExtended, processMixedEventCascBitmask, "Enable processing mixed event for cascade - cascade using bitmask for PID", false);

  // MC truth for track - cascade
  void processSameEventMCgen(const FilteredFDCollision& col, [[maybe_unused]] const aod::FDParticles& parts)
  {
    const int multCol = confUseCent ? col.multV0M() : col.multNtr();

    auto groupPartsOne = partsOneMCgenBasic->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto groupPartsTwo = partsTwoMCgenBasic->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);

    eventHisto.fillQA(col);

    for (const auto& part : groupPartsTwo) {
      int pdgCode = static_cast<int>(part.pidCut());
      if ((confCascType1 == 0 && pdgCode != kOmegaMinus) || (confCascType1 == 2 && pdgCode != kOmegaPlusBar) || (confCascType1 == 1 && pdgCode != kXiMinus) || (confCascType1 == 3 && pdgCode != kXiPlusBar))
        continue;

      cascQAHistos.fillQA<false, false>(part);

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
        if ((confCascType1 == 0 && pdgCodeCasc != kOmegaMinus) || (confCascType1 == 2 && pdgCodeCasc != kOmegaPlusBar) || (confCascType1 == 1 && pdgCodeCasc != kXiMinus) || (confCascType1 == 3 && pdgCodeCasc != kXiPlusBar))
          continue;
        sameEventCont.setPair<false>(p1, p2, multCol, confUse3D, 1.0f);
      }
    }
  }
  PROCESS_SWITCH(femtoUniversePairTaskTrackCascadeExtended, processSameEventMCgen, "Enable processing same event MC truth for track - cascade", false);

  // MC truth for cascade - cascade
  void processSameEventCascMCgen(const FilteredFDCollision& col, [[maybe_unused]] const aod::FDParticles& parts)
  {
    const int multCol = confUseCent ? col.multV0M() : col.multNtr();

    auto groupPartsTwo = partsTwoMCgenBasic->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);

    eventHisto.fillQA(col);

    for (const auto& part : groupPartsTwo) {
      int pdgCode = static_cast<int>(part.pidCut());
      if ((confCascType1 == 0 && pdgCode != kOmegaMinus) || (confCascType1 == 2 && pdgCode != kOmegaPlusBar) || (confCascType1 == 1 && pdgCode != kXiMinus) || (confCascType1 == 3 && pdgCode != kXiPlusBar))
        continue;

      cascQAHistos.fillQA<false, false>(part);

      auto pairProcessFunc = [&](auto& p1, auto& p2) -> void {
        int pdgCodeCasc1 = static_cast<int>(p1.pidCut());
        if ((confCascType1 == 0 && pdgCodeCasc1 != kOmegaMinus) || (confCascType1 == 2 && pdgCodeCasc1 != kOmegaPlusBar) || (confCascType1 == 1 && pdgCodeCasc1 != kXiMinus) || (confCascType1 == 3 && pdgCodeCasc1 != kXiPlusBar))
          return;
        int pdgCodeCasc2 = static_cast<int>(p2.pidCut());
        if ((confCascType2 == 0 && pdgCodeCasc2 != kOmegaMinus) || (confCascType2 == 2 && pdgCodeCasc2 != kOmegaPlusBar) || (confCascType2 == 1 && pdgCodeCasc2 != kXiMinus) || (confCascType2 == 3 && pdgCodeCasc2 != kXiPlusBar))
          return;
        sameEventCont.setPair<false>(p1, p2, multCol, confUse3D, 1.0f);
      };

      if (confCascType1 == confCascType2) {
        for (const auto& [p1, p2] : combinations(CombinationsStrictlyUpperIndexPolicy(groupPartsTwo, groupPartsTwo)))
          pairProcessFunc(p1, p2);
      } else {
        for (const auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(groupPartsTwo, groupPartsTwo)))
          pairProcessFunc(p1, p2);
      }
    }
  }
  PROCESS_SWITCH(femtoUniversePairTaskTrackCascadeExtended, processSameEventCascMCgen, "Enable processing same event MC truth for cascade - cascade", false);

  // MC truth for track - cascade
  void processMixedEventMCgen(const FilteredFDCollisions& cols, [[maybe_unused]] const aod::FDParticles& parts)
  {
    ColumnBinningPolicy<aod::collision::PosZ, aod::femtouniversecollision::MultNtr> colBinning{{confVtxBins, confMultBins}, true};

    for (const auto& [collision1, collision2] : soa::selfCombinations(colBinning, 5, -1, cols, cols)) {
      const int multCol = confUseCent ? collision1.multV0M() : collision1.multNtr();

      auto groupPartsOne = partsOneMCgenBasic->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);
      auto groupPartsTwo = partsTwoMCgenBasic->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision2.globalIndex(), cache);

      const auto& magFieldTesla1 = collision1.magField();
      const auto& magFieldTesla2 = collision2.magField();

      if (magFieldTesla1 != magFieldTesla2) {
        continue;
      }
      for (const auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(groupPartsOne, groupPartsTwo))) {
        if (static_cast<int>(p1.pidCut()) != confTrkPDGCodePartOne)
          continue;
        int pdgCodeCasc = static_cast<int>(p2.pidCut());
        if ((confCascType1 == 0 && pdgCodeCasc != kOmegaMinus) || (confCascType1 == 2 && pdgCodeCasc != kOmegaPlusBar) || (confCascType1 == 1 && pdgCodeCasc != kXiMinus) || (confCascType1 == 3 && pdgCodeCasc != kXiPlusBar))
          continue;
        mixedEventCont.setPair<false>(p1, p2, multCol, confUse3D, 1.0f);
      }
    }
  }
  PROCESS_SWITCH(femtoUniversePairTaskTrackCascadeExtended, processMixedEventMCgen, "Enable processing mixed event MC truth for track - cascade", false);

  // MC truth for cascade - cascade
  void processMixedEventCascMCgen(const FilteredFDCollisions& cols, [[maybe_unused]] const aod::FDParticles& parts)
  {
    ColumnBinningPolicy<aod::collision::PosZ, aod::femtouniversecollision::MultNtr> colBinning{{confVtxBins, confMultBins}, true};

    for (const auto& [collision1, collision2] : soa::selfCombinations(colBinning, 5, -1, cols, cols)) {
      const int multCol = confUseCent ? collision1.multV0M() : collision1.multNtr();

      auto groupPartsOne = partsTwoMCgenBasic->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);
      auto groupPartsTwo = partsTwoMCgenBasic->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision2.globalIndex(), cache);

      const auto& magFieldTesla1 = collision1.magField();
      const auto& magFieldTesla2 = collision2.magField();

      if (magFieldTesla1 != magFieldTesla2) {
        continue;
      }
      for (const auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(groupPartsOne, groupPartsTwo))) {
        int pdgCodeCasc1 = static_cast<int>(p1.pidCut());
        if ((confCascType1 == 0 && pdgCodeCasc1 != kOmegaMinus) || (confCascType1 == 2 && pdgCodeCasc1 != kOmegaPlusBar) || (confCascType1 == 1 && pdgCodeCasc1 != kXiMinus) || (confCascType1 == 3 && pdgCodeCasc1 != kXiPlusBar))
          continue;
        int pdgCodeCasc2 = static_cast<int>(p2.pidCut());
        if ((confCascType2 == 0 && pdgCodeCasc2 != kOmegaMinus) || (confCascType2 == 2 && pdgCodeCasc2 != kOmegaPlusBar) || (confCascType2 == 1 && pdgCodeCasc2 != kXiMinus) || (confCascType2 == 3 && pdgCodeCasc2 != kXiPlusBar))
          continue;
        mixedEventCont.setPair<false>(p1, p2, multCol, confUse3D, 1.0f);
      }
    }
  }
  PROCESS_SWITCH(femtoUniversePairTaskTrackCascadeExtended, processMixedEventCascMCgen, "Enable processing mixed event MC truth for cascade - cascade", false);

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

      if ((confCascType1 == 0 && pdgCode == kOmegaMinus) || (confCascType1 == 1 && pdgCode == kXiMinus)) {
        registryMCgen.fill(HIST("plus/MCgenCasc"), part.pt(), part.eta());
        continue;
      } else if ((confCascType1 == 0 && pdgCode == kOmegaPlusBar) || (confCascType1 == 1 && pdgCode == kXiPlusBar)) {
        registryMCgen.fill(HIST("minus/MCgenCasc"), part.pt(), part.eta());
        continue;
      }

      if (pdgParticle->Charge() > 0.0) {
        registryMCgen.fill(HIST("plus/MCgenAllPt"), part.pt());
      }
      if (pdgCode == kProton) {
        registryMCgen.fill(HIST("plus/MCgenPr"), part.pt(), part.eta());
        registryMCgen.fill(HIST("plus/MCgenPrPt"), part.pt());
      }

      if (pdgParticle->Charge() < 0.0) {
        registryMCgen.fill(HIST("minus/MCgenAllPt"), part.pt());
      }
      if (pdgCode == kProtonBar) {
        registryMCgen.fill(HIST("minus/MCgenPr"), part.pt(), part.eta());
        registryMCgen.fill(HIST("minus/MCgenPrPt"), part.pt());
      }
    }
  }

  PROCESS_SWITCH(femtoUniversePairTaskTrackCascadeExtended, processMCgen, "Process MC truth data for cascades", false);

  template <class TableType>
  void doMCReco(TableType const& parts, aod::FdMCParticles const& mcparts)
  {
    for (const auto& part : parts) {
      auto mcPartId = part.fdMCParticleId();
      if (mcPartId == -1)
        continue; // no MC particle
      const auto& mcpart = mcparts.iteratorAt(mcPartId);
      //
      if (part.partType() == aod::femtouniverseparticle::ParticleType::kCascade) {
        if (!invMCascade(part.mLambda(), part.mAntiLambda(), confCascType1)) /// mLambda stores Xi mass, mAntiLambda stores Omega mass
          continue;
        if ((confCascType1 == 0 && mcpart.pdgMCTruth() == kOmegaMinus) || (confCascType1 == 1 && mcpart.pdgMCTruth() == kXiMinus)) {
          const auto& posChild = parts.iteratorAt(part.globalIndex() - 3 - parts.begin().globalIndex());
          const auto& negChild = parts.iteratorAt(part.globalIndex() - 2 - parts.begin().globalIndex());
          const auto& bachelor = parts.iteratorAt(part.globalIndex() - 1 - parts.begin().globalIndex());
          /// Daughters that do not pass this condition are not selected
          if constexpr (std::experimental::is_detected<hasSigma, typename TableType::iterator>::value) {
            if (!isParticleTPC(posChild, CascChildTable[confCascType1][0]) || !isParticleTPC(negChild, CascChildTable[confCascType1][1]) || !isParticleTPC(bachelor, CascChildTable[confCascType1][2]))
              continue;
            if (!isParticleTOF(posChild, CascChildTable[confCascType1][0]) || !isParticleTOF(negChild, CascChildTable[confCascType1][1]) || !isParticleTOF(bachelor, CascChildTable[confCascType1][2]))
              continue;
          } else {
            if ((posChild.pidCut() & (1u << CascChildTable[confCascType1][0])) == 0 || (negChild.pidCut() & (1u << CascChildTable[confCascType1][1])) == 0 || (bachelor.pidCut() & (1u << CascChildTable[confCascType1][2])) == 0)
              continue;
            if (confUseStrangenessTOF) {
              if (((confCascType1 == 1) && (part.pidCut() & 7) != 7) || ((confCascType1 == 0) && (part.pidCut() & 56) != 56))
                continue;
            } else {
              if ((posChild.pidCut() & (8u << CascChildTable[confCascType1][0])) == 0 || (negChild.pidCut() & (8u << CascChildTable[confCascType1][1])) == 0 || (bachelor.pidCut() & (8u << CascChildTable[confCascType1][2])) == 0)
                continue;
            }
          }
          registryMCreco.fill(HIST("plus/MCrecoCascade"), mcpart.pt(), mcpart.eta());

        } else if ((confCascType1 == 0 && mcpart.pdgMCTruth() == kOmegaPlusBar) || (confCascType1 == 1 && mcpart.pdgMCTruth() == kXiPlusBar)) {
          const auto& posChild = parts.iteratorAt(part.globalIndex() - 3 - parts.begin().globalIndex());
          const auto& negChild = parts.iteratorAt(part.globalIndex() - 2 - parts.begin().globalIndex());
          const auto& bachelor = parts.iteratorAt(part.globalIndex() - 1 - parts.begin().globalIndex());
          if constexpr (std::experimental::is_detected<hasSigma, typename TableType::iterator>::value) {
            if (!isParticleTPC(posChild, CascChildTable[confCascType1 + 2][0]) && !isParticleTPC(negChild, CascChildTable[confCascType1 + 2][1]) && !isParticleTPC(bachelor, CascChildTable[confCascType1 + 2][2]))
              continue;
            if (!isParticleTOF(posChild, CascChildTable[confCascType1 + 2][0]) || !isParticleTOF(negChild, CascChildTable[confCascType1 + 2][1]) || !isParticleTOF(bachelor, CascChildTable[confCascType1 + 2][2]))
              continue;
          } else {
            if ((posChild.pidCut() & (1u << CascChildTable[confCascType1 + 2][0])) == 0 || (negChild.pidCut() & (1u << CascChildTable[confCascType1 + 2][1])) == 0 || (bachelor.pidCut() & (1u << CascChildTable[confCascType1 + 2][2])) == 0)
              continue;
            if (confUseStrangenessTOF) {
              if (((confCascType1 == 1) && (part.pidCut() & 7) != 7) || ((confCascType1 == 0) && (part.pidCut() & 56) != 56))
                continue;
            } else {
              if ((posChild.pidCut() & (8u << CascChildTable[confCascType1 + 2][0])) == 0 || (negChild.pidCut() & (8u << CascChildTable[confCascType1 + 2][1])) == 0 || (bachelor.pidCut() & (8u << CascChildTable[confCascType1 + 2][2])) == 0)
                continue;
            }
          }
          registryMCreco.fill(HIST("minus/MCrecoCascade"), mcpart.pt(), mcpart.eta());
        }

      } else if (part.partType() == aod::femtouniverseparticle::ParticleType::kTrack) {
        if (part.mAntiLambda() > 0) {
          registryMCreco.fill(HIST("plus/MCrecoAllPt"), mcpart.pt());
          if (mcpart.pdgMCTruth() != kProton)
            continue;
          if constexpr (std::experimental::is_detected<hasSigma, typename TableType::iterator>::value) {
            if (!isNSigmaCombined(part.p(), aod::pidtpc_tiny::binning::unPackInTable(part.tpcNSigmaStorePr()), aod::pidtof_tiny::binning::unPackInTable(part.tofNSigmaStorePr())))
              continue;
          } else {
            if ((part.pidCut() & 64u) == 0)
              continue;
          }
          registryMCreco.fill(HIST("plus/MCrecoPr"), mcpart.pt(), mcpart.eta());
          registryMCreco.fill(HIST("plus/MCrecoPrPt"), mcpart.pt());
        } else if (part.mAntiLambda() < 0) {
          registryMCreco.fill(HIST("minus/MCrecoAllPt"), mcpart.pt());
          if (mcpart.pdgMCTruth() != kProtonBar)
            continue;
          if constexpr (std::experimental::is_detected<hasSigma, typename TableType::iterator>::value) {
            if (!isNSigmaCombined(part.p(), aod::pidtpc_tiny::binning::unPackInTable(part.tpcNSigmaStorePr()), aod::pidtof_tiny::binning::unPackInTable(part.tofNSigmaStorePr())))
              continue;
          } else {
            if ((part.pidCut() & 64u) == 0)
              continue;
          }
          registryMCreco.fill(HIST("minus/MCrecoPr"), mcpart.pt(), mcpart.eta());
          registryMCreco.fill(HIST("minus/MCrecoPrPt"), mcpart.pt());
        }
      }
    }
  }

  void processMCReco(FemtoRecoFullParticles const& parts, aod::FdMCParticles const& mcparts)
  {
    doMCReco(parts, mcparts);
  }
  PROCESS_SWITCH(femtoUniversePairTaskTrackCascadeExtended, processMCReco, "Process MC reco data for cascades using nSigma for PID", false);

  void processMCRecoBitmask(FemtoRecoBasicParticles const& parts, aod::FdMCParticles const& mcparts)
  {
    doMCReco(parts, mcparts);
  }
  PROCESS_SWITCH(femtoUniversePairTaskTrackCascadeExtended, processMCRecoBitmask, "Process MC reco data for cascades using Bitmask for PID", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<femtoUniversePairTaskTrackCascadeExtended>(cfgc),
  };
  return workflow;
}
