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

/// \file femtoUniversePairTaskTrackV0Extended.cxx
/// \brief Tasks that build pairs of track particles and v0s
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de
/// \author Zuzanna Chochulska, WUT Warsaw & CTU Prague, zchochul@cern.ch
/// \author Shirajum Monira, WUT Warsaw, shirajum.monira.dokt@pw.edu.pl

#include "PWGCF/FemtoUniverse/Core/FemtoUniverseContainer.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseDetaDphiStar.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseEfficiencyCorrection.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseEventHisto.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniversePairCleaner.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseParticleHisto.h"

#include "CCDB/BasicCCDBManager.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"

#include <TFile.h>
#include <TH1.h>
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
using namespace o2::track;
using namespace o2::analysis::femto_universe::efficiency_correction;

struct FemtoUniversePairTaskTrackV0Extended {

  Service<o2::framework::O2DatabasePDG> pdgMC;

  SliceCache cache;
  using FemtoFullParticles = soa::Join<aod::FDParticles, aod::FDExtParticles>;
  Preslice<FemtoFullParticles> perCol = aod::femtouniverseparticle::fdCollisionId;

  using FemtoRecoParticles = soa::Join<aod::FDParticles, aod::FDExtParticles, aod::FDMCLabels, aod::FDExtMCParticles>;
  Preslice<FemtoRecoParticles> perColMC = aod::femtouniverseparticle::fdCollisionId;

  /// To apply narrow cut
  Configurable<float> confZVertexCut{"confZVertexCut", 10.f, "Event sel: Maximum z-Vertex (cm)"};
  Configurable<float> confEta{"confEta", 0.8, "Eta cut for the global track"};

  /// Particle 1 (track)
  struct : o2::framework::ConfigurableGroup {
    Configurable<int> confTrkPDGCodePartOne{"confTrkPDGCodePartOne", 211, "Particle 1 (Track) - PDG code"};
    Configurable<int> confTrackChoicePartOne{"confTrackChoicePartOne", 1, "0:Proton, 1:Pion, 2:Kaon"};
    ConfigurableAxis confTrkTempFitVarBins{"confTrkTempFitVarBins", {300, -0.15, 0.15}, "binning of the TempFitVar in the pT vs. TempFitVar plot"};
    ConfigurableAxis confTrkTempFitVarpTBins{"confTrkTempFitVarpTBins", {20, 0.5, 4.05}, "pT binning of the pT vs. TempFitVar plot"};
    Configurable<int> confChargePart1{"confChargePart1", 0, "sign of particle 1"};
    Configurable<float> confHPtPart1{"confHPtPart1", 4.0f, "higher limit for pt of particle 1"};
    Configurable<float> confLPtPart1{"confLPtPart1", 0.3f, "lower limit for pt of particle 1"};
  } ConfTrkSelection;

  Configurable<float> confmom{"confmom", 0.5, "momentum threshold for particle identification using TOF"};
  Configurable<float> confNsigmaTPCParticle{"confNsigmaTPCParticle", 3.0, "TPC Sigma for particle momentum < confmom"};
  Configurable<float> confNsigmaTPCDaughter{"confNsigmaTPCDaughter", 3.0, "TPC Sigma for daughter"};

  Configurable<float> confNsigmaTOFParticle{"confNsigmaTOFParticle", 3.0, "TOF Sigma for particle (daugh & bach) momentum > Confmom"};
  Configurable<float> confNsigmaCombinedParticle{"confNsigmaCombinedParticle", 3.0, "TPC and TOF Sigma (combined) for particle momentum > confmom"};

  Filter collisionFilter = (nabs(aod::collision::posZ) < confZVertexCut);
  using FilteredFDCollisions = soa::Filtered<o2::aod::FdCollisions>;
  using FilteredFDCollision = FilteredFDCollisions::iterator;

  /// Partition for particle 1 using extended table (track)
  Partition<FemtoFullParticles> partsOneFull = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kTrack)) && (aod::femtouniverseparticle::sign == as<int8_t>(ConfTrkSelection.confChargePart1)) && (nabs(aod::femtouniverseparticle::eta) < confEta) && (aod::femtouniverseparticle::pt < ConfTrkSelection.confHPtPart1) && (aod::femtouniverseparticle::pt > ConfTrkSelection.confLPtPart1);
  Partition<FemtoFullParticles> partsOneMCFull = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kMCTruthTrack)) && (nabs(aod::femtouniverseparticle::eta) < confEta) && (aod::femtouniverseparticle::pt < ConfTrkSelection.confHPtPart1) && (aod::femtouniverseparticle::pt > ConfTrkSelection.confLPtPart1);
  Partition<FemtoRecoParticles> partsOneMCRecoFull = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kTrack)) && (aod::femtouniverseparticle::sign == as<int8_t>(ConfTrkSelection.confChargePart1)) && (nabs(aod::femtouniverseparticle::eta) < confEta) && (aod::femtouniverseparticle::pt < ConfTrkSelection.confHPtPart1) && (aod::femtouniverseparticle::pt > ConfTrkSelection.confLPtPart1);

  /// Partition for particle 1 without extended table (track)
  Partition<aod::FDParticles> partsOneBasic = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kTrack)) && (aod::femtouniverseparticle::mAntiLambda == as<int8_t>(ConfTrkSelection.confChargePart1)) && (nabs(aod::femtouniverseparticle::eta) < confEta) && (aod::femtouniverseparticle::pt < ConfTrkSelection.confHPtPart1) && (aod::femtouniverseparticle::pt > ConfTrkSelection.confLPtPart1);
  Partition<aod::FDParticles> partsOneMCBasic = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kMCTruthTrack)) && (nabs(aod::femtouniverseparticle::eta) < confEta) && (aod::femtouniverseparticle::pt < ConfTrkSelection.confHPtPart1) && (aod::femtouniverseparticle::pt > ConfTrkSelection.confLPtPart1);
  Partition<aod::FDParticles> partsOneMCRecoBasic = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kTrack)) && (aod::femtouniverseparticle::mAntiLambda == as<int8_t>(ConfTrkSelection.confChargePart1)) && (nabs(aod::femtouniverseparticle::eta) < confEta) && (aod::femtouniverseparticle::pt < ConfTrkSelection.confHPtPart1) && (aod::femtouniverseparticle::pt > ConfTrkSelection.confLPtPart1);

  /// Histogramming for particle 1
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kTrack, 3> trackHistoPartOnePos;
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kTrack, 4> trackHistoPartOneNeg;

  /// Particle 2 (V0)
  struct : o2::framework::ConfigurableGroup {
    Configurable<int> confV0PDGCodePartTwo{"confV0PDGCodePartTwo", 3122, "Particle 2 (V0) - PDG code"};
    Configurable<int> confV0Type1{"confV0Type1", 0, "select one of the V0s (lambda = 0, anti-lambda = 1, k0 = 2) for v0-v0 and Track-v0 combination"};
    Configurable<int> confV0Type2{"confV0Type2", 0, "select one of the V0s (lambda = 0, anti-lambda = 1, k0 = 2) for v0-v0 combination"};
    ConfigurableAxis confV0TempFitVarBins{"confV0TempFitVarBins", {300, 0.95, 1.}, "V0: binning of the TempFitVar in the pT vs. TempFitVar plot"};
    ConfigurableAxis confV0TempFitVarpTBins{"confV0TempFitVarpTBins", {20, 0.5, 4.05}, "V0: pT binning of the pT vs. TempFitVar plot"};
    Configurable<float> confV0InvMassLowLimit{"confV0InvMassLowLimit", 1.10, "Lower limit of the V0 invariant mass"};
    Configurable<float> confV0InvMassUpLimit{"confV0InvMassUpLimit", 1.13, "Upper limit of the V0 invariant mass"};
    ConfigurableAxis confChildTempFitVarBins{"confChildTempFitVarBins", {300, -0.15, 0.15}, "V0 child: binning of the TempFitVar in the pT vs. TempFitVar plot"};
    ConfigurableAxis confChildTempFitVarpTBins{"confChildTempFitVarpTBins", {20, 0.5, 4.05}, "V0 child: pT binning of the pT vs. TempFitVar plot"};
    Configurable<float> confHPtPart2{"confHPtPart2", 4.0f, "higher limit for pt of particle 2"};
    Configurable<float> confLPtPart2{"confLPtPart2", 0.3f, "lower limit for pt of particle 2"};
    Configurable<bool> confUseStrangenessTOF{"confUseStrangenessTOF", true, "Use strangeness TOF for cascade PID"};
    Configurable<float> confHPtChildProton{"confHPtChildProton", 10.0f, "Higher limit for pt of children protons/antiprotons"};
    Configurable<float> confLPtChildProton{"confLPtChildProton", 0.f, "Lower limit for pt of children protons/antiprotons"};
    Configurable<float> confHPtChildPion{"confHPtChildPion", 10.0f, "Higher limit for pt of children pions"};
    Configurable<float> confLPtChildPion{"confLPtChildPion", 0.f, "Lower limit for pt of children pions"};
    Configurable<bool> confV0DuplCosPA{"confV0DuplCosPA", false, "Use cosPA instead of inv. mass as a deciding factor in rejecting a V0 in V0V0 pairs"};
    Configurable<bool> confSeparateInvMassCheck{"confSeparateInvMassCheck", false, "Apply additional cut separate for mLambda and mAntiLambda"};
  } ConfV0Selection;

  /// Partition for particle 2 using extended table
  Partition<FemtoFullParticles> partsTwoFull = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kV0)) && (aod::femtouniverseparticle::pt < ConfV0Selection.confHPtPart2) && (aod::femtouniverseparticle::pt > ConfV0Selection.confLPtPart2);
  Partition<FemtoFullParticles> partsTwoMCFull = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kMCTruthTrack)) && (aod::femtouniverseparticle::pt < ConfV0Selection.confHPtPart2) && (aod::femtouniverseparticle::pt > ConfV0Selection.confLPtPart2);
  Partition<FemtoRecoParticles> partsTwoMCRecoFull = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kV0)) && (aod::femtouniverseparticle::pt < ConfV0Selection.confHPtPart2) && (aod::femtouniverseparticle::pt > ConfV0Selection.confLPtPart2);

  /// Partition for particle 2 without extended table
  Partition<aod::FDParticles> partsTwoBasic = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kV0)) && (aod::femtouniverseparticle::pt < ConfV0Selection.confHPtPart2) && (aod::femtouniverseparticle::pt > ConfV0Selection.confLPtPart2);
  Partition<aod::FDParticles> partsTwoMCBasic = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kMCTruthTrack)) && (aod::femtouniverseparticle::pt < ConfV0Selection.confHPtPart2) && (aod::femtouniverseparticle::pt > ConfV0Selection.confLPtPart2);
  Partition<aod::FDParticles> partsTwoMCRecoBasic = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kV0)) && (aod::femtouniverseparticle::pt < ConfV0Selection.confHPtPart2) && (aod::femtouniverseparticle::pt > ConfV0Selection.confLPtPart2);

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
  Configurable<bool> confIsDebug{"confIsDebug", false, "Enable additional histograms (e.g. three-momentum)"};
  Configurable<bool> confUse3D{"confUse3D", false, "Enable three dimensional histogramms (to be used only for analysis with high statistics): k* vs mT vs multiplicity"};
  Configurable<bool> confUseCent{"confUseCent", false, "Use centrality in place of multiplicity"};
  ConfigurableAxis confMultBins{"confMultBins", {VARIABLE_WIDTH, 0.0f, 20.0f, 40.0f, 60.0f, 80.0f, 100.0f, 200.0f, 99999.f}, "Mixing bins - multiplicity"};
  ConfigurableAxis confVtxBins{"confVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  Configurable<int> confNEventsMix{"confNEventsMix", 5, "Number of events for mixing"};
  ConfigurableAxis confkstarBins{"confkstarBins", {1500, 0., 6.}, "binning kstar"};
  ConfigurableAxis confkTBins{"confkTBins", {150, 0., 9.}, "binning kT"};
  ConfigurableAxis confmTBins{"confmTBins", {225, 0., 7.5}, "binning mT"};
  Configurable<int> confPhiBins{"confPhiBins", 29, "Number of phi bins in deta dphi"};
  Configurable<int> confEtaBins{"confEtaBins", 29, "Number of eta bins in deta dphi"};
  ConfigurableAxis confmTBins3D{"confmTBins3D", {VARIABLE_WIDTH, 1.02f, 1.14f, 1.20f, 1.26f, 1.38f, 1.56f, 1.86f, 4.50f}, "mT Binning for the 3Dimensional plot: k* vs multiplicity vs mT (set <<confUse3D>> to true in order to use)"};
  ConfigurableAxis confMultBins3D{"confMultBins3D", {VARIABLE_WIDTH, 0.0f, 20.0f, 30.0f, 40.0f, 99999.0f}, "multiplicity Binning for the 3Dimensional plot: k* vs multiplicity vs mT (set <<confUse3D>> to true in order to use)"};

  struct : o2::framework::ConfigurableGroup {
    Configurable<bool> confIsCPR{"confIsCPR", true, "Close Pair Rejection"};
    Configurable<bool> confRectV0V0CPR{"confRectV0V0CPR", true, "Enable rectangular CPR cut for V0-V0 pairs"};
    Configurable<bool> confCPRPlotPerRadii{"confCPRPlotPerRadii", false, "Plot CPR per radii"};
    Configurable<float> confCPRdeltaPhiCutMax{"confCPRdeltaPhiCutMax", 0.0, "Delta Phi max cut for Close Pair Rejection"};
    Configurable<float> confCPRdeltaPhiCutMin{"confCPRdeltaPhiCutMin", 0.0, "Delta Phi min cut for Close Pair Rejection"};
    Configurable<float> confCPRdeltaEtaCutMax{"confCPRdeltaEtaCutMax", 0.0, "Delta Eta max cut for Close Pair Rejection"};
    Configurable<float> confCPRdeltaEtaCutMin{"confCPRdeltaEtaCutMin", 0.0, "Delta Eta min cut for Close Pair Rejection"};
    Configurable<float> confCPRChosenRadii{"confCPRChosenRadii", 0.80, "Delta Eta cut for Close Pair Rejection"};
  } ConfCPR;

  // Efficiency
  Configurable<std::string> confLocalEfficiency{"confLocalEfficiency", "", "Local path to efficiency .root file"};

  EffCorConfigurableGroup effCorConfGroup;
  EfficiencyCorrection effCorrection{&effCorConfGroup};

  static constexpr unsigned int V0ChildTable[][2] = {{0, 1}, {1, 0}, {1, 1}}; // Table to select the V0 children
  static constexpr double v0InvMass[] = {1.115, 1.115, 0.497};                // Table to select invariant mass of V0s

  FemtoUniverseContainer<femto_universe_container::EventType::same, femto_universe_container::Observable::kstar> sameEventCont;
  FemtoUniverseContainer<femto_universe_container::EventType::mixed, femto_universe_container::Observable::kstar> mixedEventCont;
  FemtoUniversePairCleaner<aod::femtouniverseparticle::ParticleType::kTrack, aod::femtouniverseparticle::ParticleType::kV0> pairCleaner;
  FemtoUniversePairCleaner<aod::femtouniverseparticle::ParticleType::kV0, aod::femtouniverseparticle::ParticleType::kV0> pairCleanerV0;
  FemtoUniverseDetaDphiStar<aod::femtouniverseparticle::ParticleType::kTrack, aod::femtouniverseparticle::ParticleType::kV0> pairCloseRejection;
  FemtoUniverseDetaDphiStar<aod::femtouniverseparticle::ParticleType::kV0, aod::femtouniverseparticle::ParticleType::kV0> pairCloseRejectionV0;
  /// Histogram output
  HistogramRegistry qaRegistry{"TrackQA", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry resultRegistry{"Correlations", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registryMCtruth{"MCtruthHistos", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry registryMCreco{"MCrecoHistos", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};

  std::set<int> v0Duplicates;

  std::unique_ptr<TFile> plocalEffFile;
  std::unique_ptr<TH1> pEffHistp1;
  std::unique_ptr<TH1> pEffHistp2;
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  bool isNSigmaCombined(float mom, float nsigmaTPCParticle, float nsigmaTOFParticle)
  {
    if (mom <= confmom) {
      return (std::abs(nsigmaTPCParticle) < confNsigmaTPCParticle);
    } else {
      return (std::hypot(nsigmaTOFParticle, nsigmaTPCParticle) < confNsigmaCombinedParticle);
    }
  }

  bool invMLambda(float invMassLambda, float invMassAntiLambda, int V0Type)
  {
    if (ConfV0Selection.confSeparateInvMassCheck) {
      const float pMass = V0Type ? invMassAntiLambda : invMassLambda;
      if (pMass < ConfV0Selection.confV0InvMassLowLimit || pMass > ConfV0Selection.confV0InvMassUpLimit) {
        return false;
      }
    } else {
      if ((invMassLambda < ConfV0Selection.confV0InvMassLowLimit || invMassLambda > ConfV0Selection.confV0InvMassUpLimit) && (invMassAntiLambda < ConfV0Selection.confV0InvMassLowLimit || invMassAntiLambda > ConfV0Selection.confV0InvMassUpLimit)) {
        return false;
      }
    }
    return true;
  }

  bool isNSigmaTPC(float nsigmaTPCParticle)
  {
    if (std::abs(nsigmaTPCParticle) < confNsigmaTPCDaughter) {
      return true;
    } else {
      return false;
    }
  }

  bool isNSigmaTOF(float mom, float nsigmaTOFParticle, float hasTOF)
  {
    // Cut only on daughter tracks, that have TOF signal
    if (mom > confmom && hasTOF == 1) {
      if (std::abs(nsigmaTOFParticle) < confNsigmaTOFParticle) {
        return true;
      } else {
        return false;
      }
    } else {
      return true;
    }
  }

  template <typename T>
  bool isParticleCombined(const T& part, int id)
  {
    const float tpcNSigmas[3] = {aod::pidtpc_tiny::binning::unPackInTable(part.tpcNSigmaStorePr()), aod::pidtpc_tiny::binning::unPackInTable(part.tpcNSigmaStorePi()), aod::pidtpc_tiny::binning::unPackInTable(part.tpcNSigmaStoreKa())};
    const float tofNSigmas[3] = {aod::pidtof_tiny::binning::unPackInTable(part.tofNSigmaStorePr()), aod::pidtof_tiny::binning::unPackInTable(part.tofNSigmaStorePi()), aod::pidtof_tiny::binning::unPackInTable(part.tofNSigmaStoreKa())};

    return isNSigmaCombined(part.p(), tpcNSigmas[id], tofNSigmas[id]);
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
    return isNSigmaTOF(part.p(), tofNSigmas[id], (part.pidCut() & 512u) != 0);
  }

  void init(InitContext&)
  {
    eventHisto.init(&qaRegistry);
    qaRegistry.add("Tracks_pos/nSigmaTPC", "; #it{p} (GeV/#it{c}); n#sigma_{TPC}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("Tracks_pos/nSigmaTOF", "; #it{p} (GeV/#it{c}); n#sigma_{TOF}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("Tracks_neg/nSigmaTPC", "; #it{p} (GeV/#it{c}); n#sigma_{TPC}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("Tracks_neg/nSigmaTOF", "; #it{p} (GeV/#it{c}); n#sigma_{TOF}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    trackHistoPartOnePos.init(&qaRegistry, ConfTrkSelection.confTrkTempFitVarpTBins, ConfTrkSelection.confTrkTempFitVarBins, confIsMC, ConfTrkSelection.confTrkPDGCodePartOne);
    trackHistoPartOneNeg.init(&qaRegistry, ConfTrkSelection.confTrkTempFitVarpTBins, ConfTrkSelection.confTrkTempFitVarBins, confIsMC, ConfTrkSelection.confTrkPDGCodePartOne);
    trackHistoPartTwo.init(&qaRegistry, ConfV0Selection.confV0TempFitVarpTBins, ConfV0Selection.confV0TempFitVarBins, confIsMC, ConfV0Selection.confV0PDGCodePartTwo, true);
    posChildHistos.init(&qaRegistry, ConfV0Selection.confChildTempFitVarpTBins, ConfV0Selection.confChildTempFitVarBins, false, 0, true);
    negChildHistos.init(&qaRegistry, ConfV0Selection.confChildTempFitVarpTBins, ConfV0Selection.confChildTempFitVarBins, false, 0, true);

    qaRegistry.add("V0Type1/hInvMassLambdaVsCent", "; Centrality; M_{#Lambda}; Entries", kTH2F, {confMultBins, {2000, 1.f, 3.f}});
    qaRegistry.add("V0Type2/hInvMassLambdaVsCent", "; Centrality; M_{#Lambda}; Entries", kTH2F, {confMultBins, {2000, 1.f, 3.f}});
    qaRegistry.add("V0Type1/hInvMassAntiLambdaVsCent", "; Centrality; M_{#Lambda}; Entries", kTH2F, {confMultBins, {2000, 1.f, 3.f}});
    qaRegistry.add("V0Type2/hInvMassAntiLambdaVsCent", "; Centrality; M_{#Lambda}; Entries", kTH2F, {confMultBins, {2000, 1.f, 3.f}});

    if (confIsDebug) {
      qaRegistry.add("SameEvent/hPtPosDaugh", ";  #it{p}_{T}^{1} (GeV/c);  #it{p}_{T}^{2} (GeV/c)", kTH2F, {{500, 0, 5}, {500, 0, 5}});
      qaRegistry.add("SameEvent/hPtNegDaugh", ";  #it{p}_{T}^{1} (GeV/c);  #it{p}_{T}^{2} (GeV/c)", kTH2F, {{500, 0, 5}, {500, 0, 5}});
      qaRegistry.add("SameEvent/hDaughMomPart1", "; #it{p}_{T}^{+} (GeV/c);  #it{p}_{T}^{-} (GeV/c)", kTH2F, {{500, 0, 5}, {500, 0, 5}});
      qaRegistry.add("SameEvent/hDaughMomPart2", "; #it{p}_{T}^{+} (GeV/c);  #it{p}_{T}^{-} (GeV/c)", kTH2F, {{500, 0, 5}, {500, 0, 5}});
    }

    trackHistoV0Type1.init(&qaRegistry, ConfV0Selection.confV0TempFitVarpTBins, ConfV0Selection.confV0TempFitVarBins, confIsMC, ConfV0Selection.confV0PDGCodePartTwo, true, "V0Type1");
    posChildV0Type1.init(&qaRegistry, ConfV0Selection.confChildTempFitVarpTBins, ConfV0Selection.confChildTempFitVarBins, false, 0, true, "posChildV0Type1");
    negChildV0Type1.init(&qaRegistry, ConfV0Selection.confChildTempFitVarpTBins, ConfV0Selection.confChildTempFitVarBins, false, 0, true, "negChildV0Type1");
    trackHistoV0Type2.init(&qaRegistry, ConfV0Selection.confV0TempFitVarpTBins, ConfV0Selection.confV0TempFitVarBins, confIsMC, ConfV0Selection.confV0PDGCodePartTwo, true, "V0Type2");
    posChildV0Type2.init(&qaRegistry, ConfV0Selection.confChildTempFitVarpTBins, ConfV0Selection.confChildTempFitVarBins, false, 0, true, "posChildV0Type2");
    negChildV0Type2.init(&qaRegistry, ConfV0Selection.confChildTempFitVarpTBins, ConfV0Selection.confChildTempFitVarBins, false, 0, true, "negChildV0Type2");

    qaRegistry.add("MixingQA/hMECollisionBins", ";bin;Entries", kTH1F, {{120, -0.5, 119.5}});

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

    if (doprocessPairFractionsMCTruthV0 || doprocessPairFractionsMCTruth) {
      registryMCtruth.add("mothersTruth/motherParticle", "pair fractions;part1 mother PDG;part2 mother PDG", {HistType::kTH2F, {{8001, -4000, 4000}, {8001, -4000, 4000}}});
    }

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

    if (doprocessPairFractions || doprocessPairFractionsV0) {
      registryMCreco.add("mothersReco/motherParticle", "pair fractions;part1 mother PDG;part2 mother PDG", {HistType::kTH2F, {{8001, -4000, 4000}, {8001, -4000, 4000}}});
      registryMCreco.add("mothersReco/motherParticlePDGCheck", "pair fractions;part1 mother PDG;part2 mother PDG", {HistType::kTH2F, {{8001, -4000, 4000}, {8001, -4000, 4000}}});
    }
    sameEventCont.init(&resultRegistry, confkstarBins, confMultBins, confkTBins, confmTBins, confMultBins3D, confmTBins3D, confEtaBins, confPhiBins, confIsMC, confUse3D);
    sameEventCont.setPDGCodes(ConfTrkSelection.confTrkPDGCodePartOne, ConfV0Selection.confV0PDGCodePartTwo);
    mixedEventCont.init(&resultRegistry, confkstarBins, confMultBins, confkTBins, confmTBins, confMultBins3D, confmTBins3D, confEtaBins, confPhiBins, confIsMC, confUse3D);
    mixedEventCont.setPDGCodes(ConfTrkSelection.confTrkPDGCodePartOne, ConfV0Selection.confV0PDGCodePartTwo);

    pairCleaner.init(&qaRegistry);
    pairCleanerV0.init(&qaRegistry);
    if (ConfCPR.confIsCPR.value) {
      pairCloseRejection.init(&resultRegistry, &qaRegistry, ConfCPR.confCPRdeltaPhiCutMin.value, ConfCPR.confCPRdeltaPhiCutMax.value, ConfCPR.confCPRdeltaEtaCutMin.value, ConfCPR.confCPRdeltaEtaCutMax.value, ConfCPR.confCPRChosenRadii.value, ConfCPR.confCPRPlotPerRadii.value);
      pairCloseRejectionV0.init(&resultRegistry, &qaRegistry, ConfCPR.confCPRdeltaPhiCutMin.value, ConfCPR.confCPRdeltaPhiCutMax.value, ConfCPR.confCPRdeltaEtaCutMin.value, ConfCPR.confCPRdeltaEtaCutMax.value, ConfCPR.confCPRChosenRadii.value, ConfCPR.confCPRPlotPerRadii.value);
    }

    if (!confLocalEfficiency.value.empty()) {
      plocalEffFile = std::unique_ptr<TFile>(TFile::Open(confLocalEfficiency.value.c_str(), "read"));
      if (!plocalEffFile || plocalEffFile.get()->IsZombie())
        LOGF(fatal, "Could not load efficiency histogram from %s", confLocalEfficiency.value.c_str());
      if (doprocessSameEvent || doprocessSameEventBitmask || doprocessMixedEvent || doprocessMixedEventBitmask) {
        pEffHistp1 = (ConfTrkSelection.confChargePart1 > 0) ? std::unique_ptr<TH1>(plocalEffFile.get()->Get<TH1>("PrPlus")) : std::unique_ptr<TH1>(plocalEffFile.get()->Get<TH1>("PrMinus")); // note: works only for protons for now
        pEffHistp2 = (ConfV0Selection.confV0Type1 == 0) ? std::unique_ptr<TH1>(plocalEffFile.get()->Get<TH1>("Lambda")) : std::unique_ptr<TH1>(plocalEffFile.get()->Get<TH1>("AntiLambda"));
        LOGF(info, "Loaded efficiency histograms for track-V0.");
      } else if (doprocessSameEventV0 || doprocessSameEventV0Bitmask || doprocessMixedEventV0 || doprocessMixedEventV0Bitmask) {
        pEffHistp1 = (ConfV0Selection.confV0Type1 == 0) ? std::unique_ptr<TH1>(plocalEffFile.get()->Get<TH1>("Lambda")) : std::unique_ptr<TH1>(plocalEffFile.get()->Get<TH1>("AntiLambda"));
        pEffHistp2 = (ConfV0Selection.confV0Type2 == 0) ? std::unique_ptr<TH1>(plocalEffFile.get()->Get<TH1>("Lambda")) : std::unique_ptr<TH1>(plocalEffFile.get()->Get<TH1>("AntiLambda"));
        LOGF(info, "Loaded efficiency histograms for V0-V0.");
      }
    } else if (!effCorConfGroup.confEffCorCCDBPath.value.empty()) {
      ccdb->setURL("http://alice-ccdb.cern.ch");
      ccdb->setCaching(true);
      ccdb->setLocalObjectValidityChecking();

      auto now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
      ccdb->setCreatedNotAfter(now);
      if (doprocessSameEvent || doprocessSameEventBitmask || doprocessMixedEvent || doprocessMixedEventBitmask) {
        pEffHistp1 = (ConfTrkSelection.confChargePart1 > 0) ? std::unique_ptr<TH1>(ccdb->getForTimeStamp<TH1>(effCorConfGroup.confEffCorCCDBPath.value + "/PrPlus", now)) : std::unique_ptr<TH1>(ccdb->getForTimeStamp<TH1>(effCorConfGroup.confEffCorCCDBPath.value + "/PrMinus", now)); // note: works only for protons for now
        pEffHistp2 = (ConfV0Selection.confV0Type1 == 0) ? std::unique_ptr<TH1>(ccdb->getForTimeStamp<TH1>(effCorConfGroup.confEffCorCCDBPath.value + "/Lambda", now)) : std::unique_ptr<TH1>(ccdb->getForTimeStamp<TH1>(effCorConfGroup.confEffCorCCDBPath.value + "/AntiLambda", now));
        LOGF(info, "Loaded efficiency histograms for track-V0 from CCDB.");
      } else if (doprocessSameEventV0 || doprocessSameEventV0Bitmask || doprocessMixedEventV0 || doprocessMixedEventV0Bitmask) {
        pEffHistp1 = (ConfV0Selection.confV0Type1 == 0) ? std::unique_ptr<TH1>(ccdb->getForTimeStamp<TH1>(effCorConfGroup.confEffCorCCDBPath.value + "/Lambda", now)) : std::unique_ptr<TH1>(ccdb->getForTimeStamp<TH1>(effCorConfGroup.confEffCorCCDBPath.value + "/AntiLambda", now));
        pEffHistp2 = (ConfV0Selection.confV0Type2 == 0) ? std::unique_ptr<TH1>(ccdb->getForTimeStamp<TH1>(effCorConfGroup.confEffCorCCDBPath.value + "/Lambda", now)) : std::unique_ptr<TH1>(ccdb->getForTimeStamp<TH1>(effCorConfGroup.confEffCorCCDBPath.value + "/AntiLambda", now));
        LOGF(info, "Loaded efficiency histograms for V0-V0 from CCDB.");
      }
    }

    effCorrection.init(&qaRegistry, {static_cast<framework::AxisSpec>(ConfV0Selection.confV0TempFitVarpTBins), {confEtaBins, -2, 2}, confMultBins});
  }

  template <class T>
  using hasSigma = decltype(std::declval<T&>().tpcNSigmaStorePr());

  /// This function processes the same event for track - V0
  template <typename PartType, typename PartitionType, typename MCParticles = std::nullptr_t>
  void doSameEvent(FilteredFDCollision const& col, PartType const& parts, PartitionType& partsOne, PartitionType& partsTwo, [[maybe_unused]] MCParticles mcParts = nullptr)
  {
    const auto& magFieldTesla = col.magField();

    auto groupPartsOne = partsOne->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto groupPartsTwo = partsTwo->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);

    const int multCol = confUseCent ? col.multV0M() : col.multNtr();

    eventHisto.fillQA(col);

    float v0DaughPtLowTable[3][2] = {{ConfV0Selection.confLPtChildProton, ConfV0Selection.confLPtChildPion}, {ConfV0Selection.confLPtChildPion, ConfV0Selection.confLPtChildProton}, {ConfV0Selection.confLPtChildPion, ConfV0Selection.confLPtChildPion}};
    float v0DaughPtHighTable[3][2] = {{ConfV0Selection.confHPtChildProton, ConfV0Selection.confHPtChildPion}, {ConfV0Selection.confHPtChildPion, ConfV0Selection.confHPtChildProton}, {ConfV0Selection.confHPtChildPion, ConfV0Selection.confHPtChildPion}};

    /// Histogramming same event
    for (const auto& part : groupPartsTwo) {
      if (!invMLambda(part.mLambda(), part.mAntiLambda(), ConfV0Selection.confV0Type1))
        continue;
      const auto& posChild = parts.iteratorAt(part.globalIndex() - 2 - parts.begin().globalIndex());
      const auto& negChild = parts.iteratorAt(part.globalIndex() - 1 - parts.begin().globalIndex());

      if (posChild.pt() < v0DaughPtLowTable[ConfV0Selection.confV0Type1][0] || negChild.pt() < v0DaughPtLowTable[ConfV0Selection.confV0Type1][1] || posChild.pt() > v0DaughPtHighTable[ConfV0Selection.confV0Type1][0] || negChild.pt() > v0DaughPtHighTable[ConfV0Selection.confV0Type1][1]) {
        continue;
      }

      /// Daughters that do not pass this condition are not selected
      if constexpr (std::experimental::is_detected<hasSigma, typename PartType::iterator>::value) {
        if (!isParticleTPC(posChild, V0ChildTable[ConfV0Selection.confV0Type1][0]) || !isParticleTPC(negChild, V0ChildTable[ConfV0Selection.confV0Type1][1]))
          continue;

        if (!isParticleTOF(posChild, V0ChildTable[ConfV0Selection.confV0Type1][0]) || !isParticleTOF(negChild, V0ChildTable[ConfV0Selection.confV0Type1][1]))
          continue;

        trackHistoPartTwo.fillQA<false, true>(part);
        posChildHistos.fillQA<false, true>(posChild);
        negChildHistos.fillQA<false, true>(negChild);

      } else {
        if ((posChild.pidCut() & (1u << V0ChildTable[ConfV0Selection.confV0Type1][0])) == 0 || (negChild.pidCut() & (1u << V0ChildTable[ConfV0Selection.confV0Type1][1])) == 0)
          continue;

        if (ConfV0Selection.confUseStrangenessTOF) {
          if (((ConfV0Selection.confV0Type1 == 0) && (part.pidCut() & 3) != 3) || ((ConfV0Selection.confV0Type1 == 1) && (part.pidCut() & 12) != 12) || ((ConfV0Selection.confV0Type1 == 2) && (part.pidCut() & 48) != 48))
            continue;
        } else {
          if ((posChild.pidCut() & (8u << V0ChildTable[ConfV0Selection.confV0Type1][0])) == 0 || (negChild.pidCut() & (8u << V0ChildTable[ConfV0Selection.confV0Type1][1])) == 0)
            continue;
        }
        trackHistoPartTwo.fillQA<false, false>(part);
        posChildHistos.fillQA<false, false>(posChild);
        negChildHistos.fillQA<false, false>(negChild);
      }
    }

    for (const auto& part : groupPartsOne) {
      if constexpr (std::experimental::is_detected<hasSigma, typename PartType::iterator>::value) {
        /// PID plot for particle 1
        const float tpcNSigmas[3] = {aod::pidtpc_tiny::binning::unPackInTable(part.tpcNSigmaStorePr()), aod::pidtpc_tiny::binning::unPackInTable(part.tpcNSigmaStorePi()), aod::pidtpc_tiny::binning::unPackInTable(part.tpcNSigmaStoreKa())};
        const float tofNSigmas[3] = {aod::pidtof_tiny::binning::unPackInTable(part.tofNSigmaStorePr()), aod::pidtof_tiny::binning::unPackInTable(part.tofNSigmaStorePi()), aod::pidtof_tiny::binning::unPackInTable(part.tofNSigmaStoreKa())};

        if (!isNSigmaCombined(part.p(), tpcNSigmas[ConfTrkSelection.confTrackChoicePartOne], tofNSigmas[ConfTrkSelection.confTrackChoicePartOne]))
          continue;
        if (part.sign() > 0) {
          qaRegistry.fill(HIST("Tracks_pos/nSigmaTPC"), part.p(), tpcNSigmas[ConfTrkSelection.confTrackChoicePartOne]);
          qaRegistry.fill(HIST("Tracks_pos/nSigmaTOF"), part.p(), tofNSigmas[ConfTrkSelection.confTrackChoicePartOne]);
          trackHistoPartOnePos.fillQA<false, false>(part);
        } else if (part.sign() < 0) {
          qaRegistry.fill(HIST("Tracks_neg/nSigmaTPC"), part.p(), tpcNSigmas[ConfTrkSelection.confTrackChoicePartOne]);
          qaRegistry.fill(HIST("Tracks_neg/nSigmaTOF"), part.p(), tofNSigmas[ConfTrkSelection.confTrackChoicePartOne]);
          trackHistoPartOneNeg.fillQA<false, false>(part);
        }
      } else {
        if ((part.pidCut() & (64u << ConfTrkSelection.confTrackChoicePartOne)) == 0)
          continue;
        if (ConfTrkSelection.confChargePart1 > 0)
          trackHistoPartOnePos.fillQA<false, false>(part);
        if (ConfTrkSelection.confChargePart1 < 0)
          trackHistoPartOnePos.fillQA<false, false>(part);
      }
    }

    /// Now build the combinations
    for (const auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(groupPartsOne, groupPartsTwo))) {
      // Lambda invariant mass cut
      if (!invMLambda(p2.mLambda(), p2.mAntiLambda(), ConfV0Selection.confV0Type1))
        continue;
      /// PID using stored binned nsigma
      if constexpr (std::experimental::is_detected<hasSigma, typename PartType::iterator>::value) {
        if (!isParticleCombined(p1, ConfTrkSelection.confTrackChoicePartOne))
          continue;
      } else {
        if ((p1.pidCut() & (64u << ConfTrkSelection.confTrackChoicePartOne)) == 0)
          continue;
      }
      // track cleaning
      if (!pairCleaner.isCleanPair(p1, p2, parts)) {
        continue;
      }

      const auto& posChild = parts.iteratorAt(p2.globalIndex() - 2 - parts.begin().globalIndex());
      const auto& negChild = parts.iteratorAt(p2.globalIndex() - 1 - parts.begin().globalIndex());
      if (posChild.pt() < v0DaughPtLowTable[ConfV0Selection.confV0Type1][0] || negChild.pt() < v0DaughPtLowTable[ConfV0Selection.confV0Type1][1] || posChild.pt() > v0DaughPtHighTable[ConfV0Selection.confV0Type1][0] || negChild.pt() > v0DaughPtHighTable[ConfV0Selection.confV0Type1][1]) {
        continue;
      }

      /// Daughters that do not pass this condition are not selected
      if constexpr (std::experimental::is_detected<hasSigma, typename PartType::iterator>::value) {
        if (!isParticleTPC(posChild, V0ChildTable[ConfV0Selection.confV0Type1][0]) || !isParticleTPC(negChild, V0ChildTable[ConfV0Selection.confV0Type1][1]))
          continue;

        if (!isParticleTOF(posChild, V0ChildTable[ConfV0Selection.confV0Type1][0]) || !isParticleTOF(negChild, V0ChildTable[ConfV0Selection.confV0Type1][1]))
          continue;

      } else {
        if ((posChild.pidCut() & (1u << V0ChildTable[ConfV0Selection.confV0Type1][0])) == 0 || (negChild.pidCut() & (1u << V0ChildTable[ConfV0Selection.confV0Type1][1])) == 0)
          continue;

        if (ConfV0Selection.confUseStrangenessTOF) {
          if (((ConfV0Selection.confV0Type1 == 0) && (p2.pidCut() & 3) != 3) || ((ConfV0Selection.confV0Type1 == 1) && (p2.pidCut() & 12) != 12) || ((ConfV0Selection.confV0Type1 == 2) && (p2.pidCut() & 48) != 48))
            continue;
        } else {
          if ((posChild.pidCut() & (8u << V0ChildTable[ConfV0Selection.confV0Type1][0])) == 0 || (negChild.pidCut() & (8u << V0ChildTable[ConfV0Selection.confV0Type1][1])) == 0)
            continue;
        }
      }

      if (ConfCPR.confIsCPR.value) {
        if (pairCloseRejection.isClosePair(p1, p2, parts, magFieldTesla, femto_universe_container::EventType::same)) {
          continue;
        }
      }

      float weight = 1.0f;
      if (pEffHistp1)
        weight = pEffHistp1.get()->GetBinContent(pEffHistp1->FindBin(p1.pt(), p1.eta())) * pEffHistp2.get()->GetBinContent(pEffHistp2->FindBin(p2.pt(), p2.eta()));
      if constexpr (std::is_same<PartType, FemtoRecoParticles>::value)
        sameEventCont.setPair<true>(p1, p2, multCol, confUse3D, weight);
      else
        sameEventCont.setPair<false>(p1, p2, multCol, confUse3D, weight);
    }
  }
  /// This function processes the same event for V0 - V0
  template <bool isMC, typename PartType, typename PartitionType, typename MCParticles = std::nullptr_t>
  void doSameEventV0(FilteredFDCollision const& col, PartType const& parts, PartitionType& groupPartsTwo, [[maybe_unused]] MCParticles mcParts = nullptr)
  {
    const auto& magFieldTesla = col.magField();

    const int multCol = confUseCent ? col.multV0M() : col.multNtr();

    eventHisto.fillQA(col);

    float v0DaughPtLowTable[3][2] = {{ConfV0Selection.confLPtChildProton, ConfV0Selection.confLPtChildPion}, {ConfV0Selection.confLPtChildPion, ConfV0Selection.confLPtChildProton}, {ConfV0Selection.confLPtChildPion, ConfV0Selection.confLPtChildPion}};
    float v0DaughPtHighTable[3][2] = {{ConfV0Selection.confHPtChildProton, ConfV0Selection.confHPtChildPion}, {ConfV0Selection.confHPtChildPion, ConfV0Selection.confHPtChildProton}, {ConfV0Selection.confHPtChildPion, ConfV0Selection.confHPtChildPion}};

    /// Histogramming same event for first V0 particle
    for (const auto& part : groupPartsTwo) {
      if (!invMLambda(part.mLambda(), part.mAntiLambda(), ConfV0Selection.confV0Type1))
        continue;
      const auto& posChild = parts.iteratorAt(part.globalIndex() - 2 - parts.begin().globalIndex());
      const auto& negChild = parts.iteratorAt(part.globalIndex() - 1 - parts.begin().globalIndex());
      if (posChild.pt() < v0DaughPtLowTable[ConfV0Selection.confV0Type1][0] || negChild.pt() < v0DaughPtLowTable[ConfV0Selection.confV0Type1][1] || posChild.pt() > v0DaughPtHighTable[ConfV0Selection.confV0Type1][0] || negChild.pt() > v0DaughPtHighTable[ConfV0Selection.confV0Type1][1]) {
        continue;
      }

      /// Check daughters of first V0 particle
      if constexpr (std::experimental::is_detected<hasSigma, typename PartType::iterator>::value) {
        if (!isParticleTPC(posChild, V0ChildTable[ConfV0Selection.confV0Type1][0]) || !isParticleTPC(negChild, V0ChildTable[ConfV0Selection.confV0Type1][1]))
          continue;
        if (!isParticleTOF(posChild, V0ChildTable[ConfV0Selection.confV0Type1][0]) || !isParticleTOF(negChild, V0ChildTable[ConfV0Selection.confV0Type1][1]))
          continue;

        trackHistoV0Type1.fillQABase<false, true>(part, HIST("V0Type1"));
        posChildV0Type1.fillQABase<false, true>(posChild, HIST("posChildV0Type1"));
        negChildV0Type1.fillQABase<false, true>(negChild, HIST("negChildV0Type1"));
        qaRegistry.fill(HIST("V0Type1/hInvMassLambdaVsCent"), multCol, part.mLambda());
        qaRegistry.fill(HIST("V0Type1/hInvMassAntiLambdaVsCent"), multCol, part.mAntiLambda());
      } else {
        if ((posChild.pidCut() & (1u << V0ChildTable[ConfV0Selection.confV0Type1][0])) == 0 || (negChild.pidCut() & (1u << V0ChildTable[ConfV0Selection.confV0Type1][1])) == 0)
          continue;
        if (ConfV0Selection.confUseStrangenessTOF) {
          if (((ConfV0Selection.confV0Type1 == 0) && (part.pidCut() & 3) != 3) || ((ConfV0Selection.confV0Type1 == 1) && (part.pidCut() & 12) != 12) || ((ConfV0Selection.confV0Type1 == 2) && (part.pidCut() & 48) != 48))
            continue;
        } else {
          if ((posChild.pidCut() & (8u << V0ChildTable[ConfV0Selection.confV0Type1][0])) == 0 || (negChild.pidCut() & (8u << V0ChildTable[ConfV0Selection.confV0Type1][1])) == 0)
            continue;
        }
        trackHistoV0Type1.fillQABase<false, false>(part, HIST("V0Type1"));
        posChildV0Type1.fillQABase<false, false>(posChild, HIST("posChildV0Type1"));
        negChildV0Type1.fillQABase<false, false>(negChild, HIST("negChildV0Type1"));
      }
      if constexpr (isMC) {
        effCorrection.fillRecoHist<ParticleNo::ONE, FilteredFDCollisions>(part, kLambda0);
      }
    }

    /// Histogramming same event for second V0 particle
    for (const auto& part : groupPartsTwo) {
      if (!invMLambda(part.mLambda(), part.mAntiLambda(), ConfV0Selection.confV0Type2))
        continue;
      const auto& posChild = parts.iteratorAt(part.globalIndex() - 2 - parts.begin().globalIndex());
      const auto& negChild = parts.iteratorAt(part.globalIndex() - 1 - parts.begin().globalIndex());
      if (posChild.pt() < v0DaughPtLowTable[ConfV0Selection.confV0Type2][0] || negChild.pt() < v0DaughPtLowTable[ConfV0Selection.confV0Type2][1] || posChild.pt() > v0DaughPtHighTable[ConfV0Selection.confV0Type2][0] || negChild.pt() > v0DaughPtHighTable[ConfV0Selection.confV0Type2][1]) {
        continue;
      }

      /// Check daughters of second V0 particle
      if constexpr (std::experimental::is_detected<hasSigma, typename PartType::iterator>::value) {
        if (!isParticleTPC(posChild, V0ChildTable[ConfV0Selection.confV0Type2][0]) || !isParticleTPC(negChild, V0ChildTable[ConfV0Selection.confV0Type2][1]))
          continue;
        if (!isParticleTOF(posChild, V0ChildTable[ConfV0Selection.confV0Type2][0]) || !isParticleTOF(negChild, V0ChildTable[ConfV0Selection.confV0Type2][1]))
          continue;

        trackHistoV0Type2.fillQABase<false, true>(part, HIST("V0Type2"));
        posChildV0Type2.fillQABase<false, true>(posChild, HIST("posChildV0Type2"));
        negChildV0Type2.fillQABase<false, true>(negChild, HIST("negChildV0Type2"));
        qaRegistry.fill(HIST("V0Type2/hInvMassLambdaVsCent"), multCol, part.mLambda());
        qaRegistry.fill(HIST("V0Type2/hInvMassAntiLambdaVsCent"), multCol, part.mAntiLambda());
      } else {
        if ((posChild.pidCut() & (1u << V0ChildTable[ConfV0Selection.confV0Type2][0])) == 0 || (negChild.pidCut() & (1u << V0ChildTable[ConfV0Selection.confV0Type2][1])) == 0)
          continue;
        if (ConfV0Selection.confUseStrangenessTOF) {
          if (((ConfV0Selection.confV0Type2 == 0) && (part.pidCut() & 3) != 3) || ((ConfV0Selection.confV0Type1 == 1) && (part.pidCut() & 12) != 12) || ((ConfV0Selection.confV0Type2 == 2) && (part.pidCut() & 48) != 48))
            continue;
        } else {
          if ((posChild.pidCut() & (8u << V0ChildTable[ConfV0Selection.confV0Type2][0])) == 0 || (negChild.pidCut() & (8u << V0ChildTable[ConfV0Selection.confV0Type2][1])) == 0)
            continue;
        }
        trackHistoV0Type2.fillQABase<false, false>(part, HIST("V0Type2"));
        posChildV0Type2.fillQABase<false, false>(posChild, HIST("posChildV0Type2"));
        negChildV0Type2.fillQABase<false, false>(negChild, HIST("negChildV0Type2"));
      }
      if constexpr (isMC) {
        effCorrection.fillRecoHist<ParticleNo::TWO, FilteredFDCollisions>(part, kLambda0Bar);
      }
    }

    auto pairDuplicateCheckFunc = [&](auto& p1, auto& p2) -> void {
      // V0 inv mass cut for p1
      if (!invMLambda(p1.mLambda(), p1.mAntiLambda(), ConfV0Selection.confV0Type1))
        return;
      // V0 inv mass cut for p2
      if (!invMLambda(p2.mLambda(), p2.mAntiLambda(), ConfV0Selection.confV0Type2))
        return;

      // track cleaning & checking for duplicate pairs
      if (!pairCleanerV0.isCleanPair(p1, p2, parts)) {
        // mark for rejection the cascade that shares a daughter with another cascade and has an invariant mass further from default value
        if (!ConfV0Selection.confV0DuplCosPA) {
          if (std::abs(p1.mLambda() - v0InvMass[ConfV0Selection.confV0Type1]) < std::abs(p2.mLambda() - v0InvMass[ConfV0Selection.confV0Type2])) {
            v0Duplicates.insert(p2.globalIndex());
          } else {
            v0Duplicates.insert(p1.globalIndex());
          }
        } else {
          if (std::abs(p1.tempFitVar() - 1) < std::abs(p2.tempFitVar() - 1)) {
            v0Duplicates.insert(p2.globalIndex());
          } else {
            v0Duplicates.insert(p1.globalIndex());
          }
        }
      }
    };

    auto pairProcessFunc = [&](auto& p1, auto& p2) -> bool {
      if (v0Duplicates.contains(p1.globalIndex()) || v0Duplicates.contains(p2.globalIndex()))
        return false;
      // Lambda invariant mass cut for p1
      if (!invMLambda(p1.mLambda(), p1.mAntiLambda(), ConfV0Selection.confV0Type1))
        return false;
      // Lambda invariant mass cut for p2
      if (!invMLambda(p2.mLambda(), p2.mAntiLambda(), ConfV0Selection.confV0Type2))
        return false;

      const auto& posChild1 = parts.iteratorAt(p1.globalIndex() - 2 - parts.begin().globalIndex());
      const auto& negChild1 = parts.iteratorAt(p1.globalIndex() - 1 - parts.begin().globalIndex());
      if (posChild1.pt() < v0DaughPtLowTable[ConfV0Selection.confV0Type1][0] || negChild1.pt() < v0DaughPtLowTable[ConfV0Selection.confV0Type1][1] || posChild1.pt() > v0DaughPtHighTable[ConfV0Selection.confV0Type1][0] || negChild1.pt() > v0DaughPtHighTable[ConfV0Selection.confV0Type1][1]) {
        return false;
      }

      /// p1 daughters that do not pass this condition are not selected
      if constexpr (std::experimental::is_detected<hasSigma, typename PartType::iterator>::value) {
        if (!isParticleTPC(posChild1, V0ChildTable[ConfV0Selection.confV0Type1][0]) || !isParticleTPC(negChild1, V0ChildTable[ConfV0Selection.confV0Type1][1]))
          return false;
        if (!isParticleTOF(posChild1, V0ChildTable[ConfV0Selection.confV0Type1][0]) || !isParticleTOF(negChild1, V0ChildTable[ConfV0Selection.confV0Type1][1]))
          return false;
      } else {
        if ((posChild1.pidCut() & (1u << V0ChildTable[ConfV0Selection.confV0Type1][0])) == 0 || (negChild1.pidCut() & (1u << V0ChildTable[ConfV0Selection.confV0Type1][1])) == 0)
          return false;
        if (ConfV0Selection.confUseStrangenessTOF) {
          if (((ConfV0Selection.confV0Type1 == 0) && (p1.pidCut() & 3) != 3) || ((ConfV0Selection.confV0Type1 == 1) && (p1.pidCut() & 12) != 12) || ((ConfV0Selection.confV0Type1 == 2) && (p1.pidCut() & 48) != 48))
            return false;
        } else {
          if ((posChild1.pidCut() & (8u << V0ChildTable[ConfV0Selection.confV0Type1][0])) == 0 || (negChild1.pidCut() & (8u << V0ChildTable[ConfV0Selection.confV0Type1][1])) == 0)
            return false;
        }
      }

      const auto& posChild2 = parts.iteratorAt(p2.globalIndex() - 2 - parts.begin().globalIndex());
      const auto& negChild2 = parts.iteratorAt(p2.globalIndex() - 1 - parts.begin().globalIndex());
      if (posChild2.pt() < v0DaughPtLowTable[ConfV0Selection.confV0Type2][0] || negChild2.pt() < v0DaughPtLowTable[ConfV0Selection.confV0Type2][1] || posChild2.pt() > v0DaughPtHighTable[ConfV0Selection.confV0Type2][0] || negChild2.pt() > v0DaughPtHighTable[ConfV0Selection.confV0Type2][1]) {
        return false;
      }

      /// p2 daughters that do not pass this condition are not selected
      if constexpr (std::experimental::is_detected<hasSigma, typename PartType::iterator>::value) {
        if (!isParticleTPC(posChild2, V0ChildTable[ConfV0Selection.confV0Type2][0]) || !isParticleTPC(negChild2, V0ChildTable[ConfV0Selection.confV0Type2][1]))
          return false;
        if (!isParticleTOF(posChild2, V0ChildTable[ConfV0Selection.confV0Type2][0]) || !isParticleTOF(negChild2, V0ChildTable[ConfV0Selection.confV0Type2][1]))
          return false;
      } else {
        if ((posChild2.pidCut() & (1u << V0ChildTable[ConfV0Selection.confV0Type2][0])) == 0 || (negChild2.pidCut() & (1u << V0ChildTable[ConfV0Selection.confV0Type2][1])) == 0)
          return false;
        if (ConfV0Selection.confUseStrangenessTOF) {
          if (((ConfV0Selection.confV0Type2 == 0) && (p2.pidCut() & 3) != 3) || ((ConfV0Selection.confV0Type2 == 1) && (p2.pidCut() & 12) != 12) || ((ConfV0Selection.confV0Type2 == 2) && (p2.pidCut() & 48) != 48))
            return false;
        } else {
          if ((posChild2.pidCut() & (8u << V0ChildTable[ConfV0Selection.confV0Type2][0])) == 0 || (negChild2.pidCut() & (8u << V0ChildTable[ConfV0Selection.confV0Type2][1])) == 0)
            return false;
        }
      }

      if (ConfCPR.confIsCPR.value) {
        if (ConfCPR.confRectV0V0CPR && pairCloseRejectionV0.isClosePair<true>(p1, p2, parts, magFieldTesla, femto_universe_container::EventType::same)) {
          return false;
        } else if (!ConfCPR.confRectV0V0CPR && pairCloseRejectionV0.isClosePair<false>(p1, p2, parts, magFieldTesla, femto_universe_container::EventType::same)) {
          return false;
        }
      }

      if constexpr (std::is_same<PartType, FemtoRecoParticles>::value)
        sameEventCont.setPair<true>(p1, p2, multCol, confUse3D);
      else
        sameEventCont.setPair<false>(p1, p2, multCol, confUse3D);
      return true;
    };

    v0Duplicates.clear();
    for (const auto& [p1, p2] : combinations(CombinationsStrictlyUpperIndexPolicy(groupPartsTwo, groupPartsTwo))) {
      pairDuplicateCheckFunc(p1, p2);
    }
    /// Now build the combinations for V0s
    for (const auto& [p1, p2] : combinations(CombinationsStrictlyUpperIndexPolicy(groupPartsTwo, groupPartsTwo))) {
      if (!pairProcessFunc(p1, p2)) {
        if (pairProcessFunc(p2, p1) && confIsDebug) {
          const auto& posChild1 = parts.iteratorAt(p2.globalIndex() - 2 - parts.begin().globalIndex());
          const auto& negChild1 = parts.iteratorAt(p2.globalIndex() - 1 - parts.begin().globalIndex());
          const auto& posChild2 = parts.iteratorAt(p1.globalIndex() - 2 - parts.begin().globalIndex());
          const auto& negChild2 = parts.iteratorAt(p1.globalIndex() - 1 - parts.begin().globalIndex());
          qaRegistry.fill(HIST("SameEvent/hPtPosDaugh"), posChild1.pt(), posChild2.pt());
          qaRegistry.fill(HIST("SameEvent/hPtNegDaugh"), negChild1.pt(), negChild2.pt());
          qaRegistry.fill(HIST("SameEvent/hDaughMomPart1"), posChild1.pt(), negChild1.pt());
          qaRegistry.fill(HIST("SameEvent/hDaughMomPart2"), posChild2.pt(), negChild2.pt());
        }
      } else if (confIsDebug) {
        const auto& posChild1 = parts.iteratorAt(p1.globalIndex() - 2 - parts.begin().globalIndex());
        const auto& negChild1 = parts.iteratorAt(p1.globalIndex() - 1 - parts.begin().globalIndex());
        const auto& posChild2 = parts.iteratorAt(p2.globalIndex() - 2 - parts.begin().globalIndex());
        const auto& negChild2 = parts.iteratorAt(p2.globalIndex() - 1 - parts.begin().globalIndex());
        qaRegistry.fill(HIST("SameEvent/hPtPosDaugh"), posChild1.pt(), posChild2.pt());
        qaRegistry.fill(HIST("SameEvent/hPtNegDaugh"), negChild1.pt(), negChild2.pt());
        qaRegistry.fill(HIST("SameEvent/hDaughMomPart1"), posChild1.pt(), negChild1.pt());
        qaRegistry.fill(HIST("SameEvent/hDaughMomPart2"), posChild2.pt(), negChild2.pt());
      }
    }
  }

  void processSameEvent(FilteredFDCollision const& col, FemtoFullParticles const& parts)
  {
    doSameEvent(col, parts, partsOneFull, partsTwoFull);
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackV0Extended, processSameEvent, "Enable processing same event for track - V0", false);

  void processSameEventBitmask(FilteredFDCollision const& col, aod::FDParticles const& parts)
  {
    doSameEvent(col, parts, partsOneBasic, partsTwoBasic);
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackV0Extended, processSameEventBitmask, "Enable processing same event for track - V0 using bitmask", false);

  void processSameEventMCReco(FilteredFDCollision const& col, FemtoRecoParticles const& parts, aod::FdMCParticles const& mcparts)
  {
    doSameEvent(col, parts, partsOneMCRecoFull, partsTwoMCRecoFull, mcparts);
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackV0Extended, processSameEventMCReco, "Enable processing same event for track - V0 MC Reco", false);

  /// This function processes the same event for V0 - V0
  void processSameEventV0(FilteredFDCollision const& col, FemtoFullParticles const& parts)
  {
    auto groupPartsTwo = partsTwoFull->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    doSameEventV0<false>(col, parts, groupPartsTwo);
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackV0Extended, processSameEventV0, "Enable processing same event for V0 - V0", false);

  void processSameEventV0Bitmask(FilteredFDCollision const& col, aod::FDParticles const& parts)
  {
    auto groupPartsTwo = partsTwoBasic->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    doSameEventV0<false>(col, parts, groupPartsTwo);
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackV0Extended, processSameEventV0Bitmask, "Enable processing same event for V0 - V0 using bitmask", false);

  void processSameEventV0MCReco(FilteredFDCollision const& col, FemtoRecoParticles const& parts, aod::FdMCParticles const& mcparts)
  {
    auto groupPartsTwo = partsTwoMCRecoFull->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    doSameEventV0<true>(col, parts, groupPartsTwo, mcparts);
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackV0Extended, processSameEventV0MCReco, "Enable processing same event for V0 - V0 MC Reco", false);

  /// This function processes MC same events for Track - V0
  void processMCSameEvent(FilteredFDCollision const& col, FemtoFullParticles const& parts)
  {
    const auto& magFieldTesla = col.magField();

    auto groupPartsOne = partsOneMCFull->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto groupPartsTwo = partsTwoMCFull->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    const int multCol = confUseCent ? col.multV0M() : col.multNtr();

    eventHisto.fillQA(col);

    /// Histogramming same event
    for (const auto& part : groupPartsTwo) {
      int pdgCode = static_cast<int>(part.pidCut());
      if ((ConfV0Selection.confV0Type1 == 0 && pdgCode != kLambda0) || (ConfV0Selection.confV0Type1 == 1 && pdgCode != kLambda0Bar))
        continue;
      trackHistoPartTwo.fillQA<false, true>(part);
    }

    for (const auto& part : groupPartsOne) {
      int pdgCode = static_cast<int>(part.pidCut());
      if (pdgCode != ConfTrkSelection.confTrkPDGCodePartOne)
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
      if (static_cast<int>(p1.pidCut()) != ConfTrkSelection.confTrkPDGCodePartOne)
        continue;
      int pdgCode2 = static_cast<int>(p2.pidCut());
      if ((ConfV0Selection.confV0Type1 == 0 && pdgCode2 != kLambda0) || (ConfV0Selection.confV0Type1 == 1 && pdgCode2 != kLambda0Bar))
        continue;
      // track cleaning
      if (ConfCPR.confIsCPR.value) {
        if (pairCloseRejection.isClosePair(p1, p2, parts, magFieldTesla, femto_universe_container::EventType::same)) {
          continue;
        }
      }
      sameEventCont.setPair<false>(p1, p2, multCol, confUse3D);
    }
  }

  PROCESS_SWITCH(FemtoUniversePairTaskTrackV0Extended, processMCSameEvent, "Enable processing same event for MC truth track - V0", false);

  /// This function processes MC same events for V0 - V0
  void processMCSameEventV0(FilteredFDCollision const& col, FemtoFullParticles const& parts)
  {
    auto groupPartsTwo = partsTwoMCFull->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    const int multCol = confUseCent ? col.multV0M() : col.multNtr();

    eventHisto.fillQA(col);

    /// Histogramming same event
    for (const auto& part : groupPartsTwo) {
      const auto& posChild = parts.iteratorAt(part.index() - 2);
      const auto& negChild = parts.iteratorAt(part.index() - 1);
      int pdgCode = static_cast<int>(part.pidCut());
      if ((ConfV0Selection.confV0Type1 == 0 && pdgCode == kLambda0) || (ConfV0Selection.confV0Type1 == 1 && pdgCode == kLambda0Bar)) {
        trackHistoV0Type1.fillQABase<false, true>(part, HIST("V0Type1"));
        posChildV0Type1.fillQABase<false, true>(posChild, HIST("posChildV0Type1"));
        negChildV0Type1.fillQABase<false, true>(negChild, HIST("negChildV0Type1"));
        qaRegistry.fill(HIST("V0Type1/hInvMassLambdaVsCent"), multCol, part.mLambda());
        qaRegistry.fill(HIST("V0Type1/hInvMassAntiLambdaVsCent"), multCol, part.mAntiLambda());
        effCorrection.fillTruthHist<ParticleNo::ONE, FilteredFDCollisions>(part);
      }
      if ((ConfV0Selection.confV0Type2 == 0 && pdgCode == kLambda0) || (ConfV0Selection.confV0Type2 == 1 && pdgCode == kLambda0Bar)) {
        trackHistoV0Type2.fillQABase<false, true>(part, HIST("V0Type2"));
        posChildV0Type2.fillQABase<false, true>(posChild, HIST("posChildV0Type2"));
        negChildV0Type2.fillQABase<false, true>(negChild, HIST("negChildV0Type2"));
        qaRegistry.fill(HIST("V0Type2/hInvMassLambdaVsCent"), multCol, part.mLambda());
        qaRegistry.fill(HIST("V0Type2/hInvMassAntiLambdaVsCent"), multCol, part.mAntiLambda());
        effCorrection.fillTruthHist<ParticleNo::TWO, FilteredFDCollisions>(part);
      }
    }

    auto pairProcessFunc = [&](auto& p1, auto& p2) -> void {
      int pdgCode1 = static_cast<int>(p1.pidCut());
      if ((ConfV0Selection.confV0Type1 == 0 && pdgCode1 != kLambda0) || (ConfV0Selection.confV0Type1 == 1 && pdgCode1 != kLambda0Bar))
        return;
      int pdgCode2 = static_cast<int>(p2.pidCut());
      if ((ConfV0Selection.confV0Type2 == 0 && pdgCode2 != kLambda0) || (ConfV0Selection.confV0Type2 == 1 && pdgCode2 != kLambda0Bar))
        return;
      sameEventCont.setPair<false>(p1, p2, multCol, confUse3D);
    };
    /// Now build the combinations
    if (ConfV0Selection.confV0Type1 == ConfV0Selection.confV0Type2) {
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
  void doMixedEvent(FilteredFDCollisions const& cols, PartType const& parts, PartitionType& partsOne, PartitionType& partsTwo, [[maybe_unused]] MCParticles mcParts = nullptr)
  {
    ColumnBinningPolicy<aod::collision::PosZ, aod::femtouniversecollision::MultNtr> colBinningMult{{confVtxBins, confMultBins}, true};
    ColumnBinningPolicy<aod::collision::PosZ, aod::femtouniversecollision::MultV0M> colBinningCent{{confVtxBins, confMultBins}, true};

    float v0DaughPtLowTable[3][2] = {{ConfV0Selection.confLPtChildProton, ConfV0Selection.confLPtChildPion}, {ConfV0Selection.confLPtChildPion, ConfV0Selection.confLPtChildProton}, {ConfV0Selection.confLPtChildPion, ConfV0Selection.confLPtChildPion}};
    float v0DaughPtHighTable[3][2] = {{ConfV0Selection.confHPtChildProton, ConfV0Selection.confHPtChildPion}, {ConfV0Selection.confHPtChildPion, ConfV0Selection.confHPtChildProton}, {ConfV0Selection.confHPtChildPion, ConfV0Selection.confHPtChildPion}};

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
        // Lambda invariant mass cut
        if (!invMLambda(p2.mLambda(), p2.mAntiLambda(), ConfV0Selection.confV0Type1))
          continue;
        /// PID using stored binned nsigma
        if constexpr (std::experimental::is_detected<hasSigma, typename PartType::iterator>::value) {
          if (!isParticleCombined(p1, ConfTrkSelection.confTrackChoicePartOne))
            continue;
        } else {
          if ((p1.pidCut() & (64u << ConfTrkSelection.confTrackChoicePartOne)) == 0)
            continue;
        }

        const auto& posChild = parts.iteratorAt(p2.globalIndex() - 2 - parts.begin().globalIndex());
        const auto& negChild = parts.iteratorAt(p2.globalIndex() - 1 - parts.begin().globalIndex());
        if (posChild.pt() < v0DaughPtLowTable[ConfV0Selection.confV0Type1][0] || negChild.pt() < v0DaughPtLowTable[ConfV0Selection.confV0Type1][1] || posChild.pt() > v0DaughPtHighTable[ConfV0Selection.confV0Type1][0] || negChild.pt() > v0DaughPtHighTable[ConfV0Selection.confV0Type1][1]) {
          continue;
        }

        /// Daughters that do not pass this condition are not selected
        if constexpr (std::experimental::is_detected<hasSigma, typename PartType::iterator>::value) {
          if (!isParticleTPC(posChild, V0ChildTable[ConfV0Selection.confV0Type1][0]) || !isParticleTPC(negChild, V0ChildTable[ConfV0Selection.confV0Type1][1]))
            continue;
          if (!isParticleTOF(posChild, V0ChildTable[ConfV0Selection.confV0Type1][0]) || !isParticleTOF(negChild, V0ChildTable[ConfV0Selection.confV0Type1][1]))
            continue;
        } else {
          if ((posChild.pidCut() & (1u << V0ChildTable[ConfV0Selection.confV0Type1][0])) == 0 || (negChild.pidCut() & (1u << V0ChildTable[ConfV0Selection.confV0Type1][1])) == 0)
            continue;
          if (ConfV0Selection.confUseStrangenessTOF) {
            if (((ConfV0Selection.confV0Type1 == 0) && (p2.pidCut() & 3) != 3) || ((ConfV0Selection.confV0Type1 == 1) && (p2.pidCut() & 12) != 12) || ((ConfV0Selection.confV0Type1 == 2) && (p2.pidCut() & 48) != 48))
              continue;
          } else {
            if ((posChild.pidCut() & (8u << V0ChildTable[ConfV0Selection.confV0Type1][0])) == 0 || (negChild.pidCut() & (8u << V0ChildTable[ConfV0Selection.confV0Type1][1])) == 0)
              continue;
          }
        }

        // track cleaning
        if (!pairCleaner.isCleanPair(p1, p2, parts)) {
          continue;
        }
        if (ConfCPR.confIsCPR.value) {
          if (pairCloseRejection.isClosePair(p1, p2, parts, magFieldTesla1, femto_universe_container::EventType::mixed)) {
            continue;
          }
        }
        float weight = 1.0f;
        if (pEffHistp1)
          weight = pEffHistp1.get()->GetBinContent(pEffHistp1->FindBin(p1.pt(), p1.eta())) * pEffHistp2.get()->GetBinContent(pEffHistp2->FindBin(p2.pt(), p2.eta()));

        if constexpr (std::is_same<PartType, FemtoRecoParticles>::value)
          mixedEventCont.setPair<true>(p1, p2, multCol, confUse3D, weight);
        else
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

  /// This function processes the mixed event for V0 - V0
  template <typename PartType, typename PartitionType, typename MCParticles = std::nullptr_t>
  void doMixedEventV0(FilteredFDCollisions const& cols, PartType const& parts, PartitionType& partsTwo, [[maybe_unused]] MCParticles mcParts = nullptr)
  {
    ColumnBinningPolicy<aod::collision::PosZ, aod::femtouniversecollision::MultNtr> colBinningMult{{confVtxBins, confMultBins}, true};
    ColumnBinningPolicy<aod::collision::PosZ, aod::femtouniversecollision::MultV0M> colBinningCent{{confVtxBins, confMultBins}, true};

    float v0DaughPtLowTable[3][2] = {{ConfV0Selection.confLPtChildProton, ConfV0Selection.confLPtChildPion}, {ConfV0Selection.confLPtChildPion, ConfV0Selection.confLPtChildProton}, {ConfV0Selection.confLPtChildPion, ConfV0Selection.confLPtChildPion}};
    float v0DaughPtHighTable[3][2] = {{ConfV0Selection.confHPtChildProton, ConfV0Selection.confHPtChildPion}, {ConfV0Selection.confHPtChildPion, ConfV0Selection.confHPtChildProton}, {ConfV0Selection.confHPtChildPion, ConfV0Selection.confHPtChildPion}};

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
        if (!invMLambda(p1.mLambda(), p1.mAntiLambda(), ConfV0Selection.confV0Type1)) {
          continue;
        }
        // Lambda invariant mass cut for p2
        if (!invMLambda(p2.mLambda(), p2.mAntiLambda(), ConfV0Selection.confV0Type2)) {
          continue;
        }

        const auto& posChild1 = parts.iteratorAt(p1.globalIndex() - 2 - parts.begin().globalIndex());
        const auto& negChild1 = parts.iteratorAt(p1.globalIndex() - 1 - parts.begin().globalIndex());

        if (posChild1.pt() < v0DaughPtLowTable[ConfV0Selection.confV0Type1][0] || negChild1.pt() < v0DaughPtLowTable[ConfV0Selection.confV0Type1][1] || posChild1.pt() > v0DaughPtHighTable[ConfV0Selection.confV0Type1][0] || negChild1.pt() > v0DaughPtHighTable[ConfV0Selection.confV0Type1][1]) {
          continue;
        }

        /// Daughters that do not pass this condition are not selected
        if constexpr (std::experimental::is_detected<hasSigma, typename PartType::iterator>::value) {
          if (!isParticleTPC(posChild1, V0ChildTable[ConfV0Selection.confV0Type1][0]) || !isParticleTPC(negChild1, V0ChildTable[ConfV0Selection.confV0Type1][1]))
            continue;
          if (!isParticleTOF(posChild1, V0ChildTable[ConfV0Selection.confV0Type1][0]) || !isParticleTOF(negChild1, V0ChildTable[ConfV0Selection.confV0Type1][1]))
            continue;
        } else {
          if ((posChild1.pidCut() & (1u << V0ChildTable[ConfV0Selection.confV0Type1][0])) == 0 || (negChild1.pidCut() & (1u << V0ChildTable[ConfV0Selection.confV0Type1][1])) == 0)
            continue;
          if (ConfV0Selection.confUseStrangenessTOF) {
            if (((ConfV0Selection.confV0Type1 == 0) && (p1.pidCut() & 3) != 3) || ((ConfV0Selection.confV0Type1 == 1) && (p1.pidCut() & 12) != 12) || ((ConfV0Selection.confV0Type1 == 2) && (p1.pidCut() & 48) != 48))
              continue;
          } else {
            if ((posChild1.pidCut() & (8u << V0ChildTable[ConfV0Selection.confV0Type1][0])) == 0 || (negChild1.pidCut() & (8u << V0ChildTable[ConfV0Selection.confV0Type1][1])) == 0)
              continue;
          }
        }

        const auto& posChild2 = parts.iteratorAt(p2.globalIndex() - 2 - parts.begin().globalIndex());
        const auto& negChild2 = parts.iteratorAt(p2.globalIndex() - 1 - parts.begin().globalIndex());
        if (posChild2.pt() < v0DaughPtLowTable[ConfV0Selection.confV0Type2][0] || negChild2.pt() < v0DaughPtLowTable[ConfV0Selection.confV0Type2][1] || posChild2.pt() > v0DaughPtHighTable[ConfV0Selection.confV0Type2][0] || negChild2.pt() > v0DaughPtHighTable[ConfV0Selection.confV0Type2][1]) {
          continue;
        }

        /// Daughters that do not pass this condition are not selected
        if constexpr (std::experimental::is_detected<hasSigma, typename PartType::iterator>::value) {
          if (!isParticleTPC(posChild2, V0ChildTable[ConfV0Selection.confV0Type2][0]) || !isParticleTPC(negChild2, V0ChildTable[ConfV0Selection.confV0Type2][1]))
            continue;
          if (!isParticleTOF(posChild2, V0ChildTable[ConfV0Selection.confV0Type2][0]) || !isParticleTOF(negChild2, V0ChildTable[ConfV0Selection.confV0Type2][1]))
            continue;
        } else {
          if ((posChild2.pidCut() & (1u << V0ChildTable[ConfV0Selection.confV0Type2][0])) == 0 || (negChild2.pidCut() & (1u << V0ChildTable[ConfV0Selection.confV0Type2][1])) == 0)
            continue;
          if (ConfV0Selection.confUseStrangenessTOF) {
            if (((ConfV0Selection.confV0Type2 == 0) && (p2.pidCut() & 3) != 3) || ((ConfV0Selection.confV0Type2 == 1) && (p2.pidCut() & 12) != 12) || ((ConfV0Selection.confV0Type2 == 2) && (p2.pidCut() & 48) != 48))
              continue;
          } else {
            if ((posChild2.pidCut() & (8u << V0ChildTable[ConfV0Selection.confV0Type2][0])) == 0 || (negChild2.pidCut() & (8u << V0ChildTable[ConfV0Selection.confV0Type2][1])) == 0)
              continue;
          }
        }

        // track cleaning
        if (!pairCleanerV0.isCleanPair(p1, p2, parts)) {
          continue;
        }
        if (ConfCPR.confIsCPR.value) {
          if (ConfCPR.confRectV0V0CPR && pairCloseRejectionV0.isClosePair<true>(p1, p2, parts, magFieldTesla1, femto_universe_container::EventType::mixed)) {
            continue;
          } else if (!ConfCPR.confRectV0V0CPR && pairCloseRejectionV0.isClosePair<false>(p1, p2, parts, magFieldTesla1, femto_universe_container::EventType::mixed)) {
            continue;
          }
        }

        if constexpr (std::is_same<PartType, FemtoRecoParticles>::value)
          mixedEventCont.setPair<true>(p1, p2, multCol, confUse3D);
        else
          mixedEventCont.setPair<false>(p1, p2, multCol, confUse3D);
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

  void processMixedEvent(FilteredFDCollisions const& cols, FemtoFullParticles const& parts)
  {
    doMixedEvent(cols, parts, partsOneFull, partsTwoFull);
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackV0Extended, processMixedEvent, "Enable processing mixed event for track - V0", false);

  void processMixedEventBitmask(FilteredFDCollisions const& cols, aod::FDParticles const& parts)
  {
    doMixedEvent(cols, parts, partsOneBasic, partsTwoBasic);
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackV0Extended, processMixedEventBitmask, "Enable processing mixed event for track - V0 using bitmask", false);

  void processMixedEventMCReco(FilteredFDCollisions const& cols, FemtoRecoParticles const& parts, aod::FdMCParticles const& mcparts)
  {
    doMixedEvent(cols, parts, partsOneMCRecoFull, partsTwoMCRecoFull, mcparts);
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackV0Extended, processMixedEventMCReco, "Enable processing mixed event for track - V0 for MC Reco", false);

  void processMixedEventV0(FilteredFDCollisions const& cols, FemtoFullParticles const& parts)
  {
    doMixedEventV0(cols, parts, partsTwoFull);
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackV0Extended, processMixedEventV0, "Enable processing mixed events for V0 - V0", false);

  void processMixedEventV0Bitmask(FilteredFDCollisions const& cols, aod::FDParticles const& parts)
  {
    doMixedEventV0(cols, parts, partsTwoBasic);
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackV0Extended, processMixedEventV0Bitmask, "Enable processing mixed events for V0 - V0 using bitmask", false);

  void processMixedEventV0MCReco(FilteredFDCollisions const& cols, FemtoRecoParticles const& parts, aod::FdMCParticles const& mcparts)
  {
    doMixedEventV0(cols, parts, partsTwoMCRecoFull, mcparts);
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackV0Extended, processMixedEventV0MCReco, "Enable processing mixed event for V0 - V0 for MC Reco", false);

  /// This function processes MC mixed events for Track - V0
  void processMCMixedEvent(FilteredFDCollisions const& cols, FemtoFullParticles const& parts)
  {
    ColumnBinningPolicy<aod::collision::PosZ, aod::femtouniversecollision::MultNtr> colBinningMult{{confVtxBins, confMultBins}, true};
    ColumnBinningPolicy<aod::collision::PosZ, aod::femtouniversecollision::MultV0M> colBinningCent{{confVtxBins, confMultBins}, true};

    auto mixedCollProcessFunc = [&](auto& collision1, auto& collision2) -> void {
      const int multCol = confUseCent ? collision1.multV0M() : collision1.multNtr();

      auto groupPartsOne = partsOneMCFull->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);
      auto groupPartsTwo = partsTwoMCFull->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision2.globalIndex(), cache);

      const auto& magFieldTesla1 = collision1.magField();
      const auto& magFieldTesla2 = collision2.magField();

      if (magFieldTesla1 != magFieldTesla2) {
        return;
      }
      for (const auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(groupPartsOne, groupPartsTwo))) {
        if (static_cast<int>(p1.pidCut()) != ConfTrkSelection.confTrkPDGCodePartOne)
          continue;
        int pdgCode2 = static_cast<int>(p2.pidCut());
        if ((ConfV0Selection.confV0Type1 == 0 && pdgCode2 != kLambda0) || (ConfV0Selection.confV0Type1 == 1 && pdgCode2 != kLambda0Bar))
          continue;
        if (ConfCPR.confIsCPR.value) {
          if (pairCloseRejection.isClosePair(p1, p2, parts, magFieldTesla1, femto_universe_container::EventType::mixed)) {
            continue;
          }
        }
        mixedEventCont.setPair<false>(p1, p2, multCol, confUse3D);
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

  PROCESS_SWITCH(FemtoUniversePairTaskTrackV0Extended, processMCMixedEvent, "Enable processing mixed events for MC truth track - V0", false);

  /// This function processes MC mixed events for V0 - V0
  void processMCMixedEventV0(FilteredFDCollisions const& cols, FemtoFullParticles const& /*parts*/)
  {
    ColumnBinningPolicy<aod::collision::PosZ, aod::femtouniversecollision::MultNtr> colBinningMult{{confVtxBins, confMultBins}, true};
    ColumnBinningPolicy<aod::collision::PosZ, aod::femtouniversecollision::MultV0M> colBinningCent{{confVtxBins, confMultBins}, true};

    auto mixedCollProcessFunc = [&](auto& collision1, auto& collision2) -> void {
      const int multCol = confUseCent ? collision1.multV0M() : collision1.multNtr();

      auto groupPartsOne = partsTwoMCFull->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);
      auto groupPartsTwo = partsTwoMCFull->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision2.globalIndex(), cache);

      for (const auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(groupPartsOne, groupPartsTwo))) {
        int pdgCode1 = static_cast<int>(p1.pidCut());
        if ((ConfV0Selection.confV0Type1 == 0 && pdgCode1 != kLambda0) || (ConfV0Selection.confV0Type1 == 1 && pdgCode1 != kLambda0Bar))
          continue;
        int pdgCode2 = static_cast<int>(p2.pidCut());
        if ((ConfV0Selection.confV0Type2 == 0 && pdgCode2 != kLambda0) || (ConfV0Selection.confV0Type2 == 1 && pdgCode2 != kLambda0Bar))
          continue;
        mixedEventCont.setPair<false>(p1, p2, multCol, confUse3D);
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

  PROCESS_SWITCH(FemtoUniversePairTaskTrackV0Extended, processMCMixedEventV0, "Enable processing mixed events for MC truth V0 - V0", false);

  ///--------------------------------------------MC-------------------------------------------------///

  /// This function fills MC truth particles from derived MC table
  void processMCTruth(aod::FDParticles const& parts)
  {
    for (const auto& part : parts) {
      if (part.partType() != uint8_t(aod::femtouniverseparticle::ParticleType::kMCTruthTrack))
        continue;

      int pdgCode = static_cast<int>(part.pidCut());
      const auto& pdgParticle = pdgMC->GetParticle(pdgCode);
      if (!pdgParticle) {
        continue;
      }

      if (pdgCode == kLambda0) {
        registryMCtruth.fill(HIST("plus/MCtruthLambda"), part.pt(), part.eta());
        continue;
      } else if (pdgCode == kLambda0Bar) {
        registryMCtruth.fill(HIST("minus/MCtruthLambda"), part.pt(), part.eta());
        continue;
      }

      if (pdgParticle->Charge() > 0.0) {
        registryMCtruth.fill(HIST("plus/MCtruthAllPt"), part.pt());
      }
      if (pdgCode == kPiPlus) {
        registryMCtruth.fill(HIST("plus/MCtruthPi"), part.pt(), part.eta());
        registryMCtruth.fill(HIST("plus/MCtruthPiPt"), part.pt());
      }
      if (pdgCode == kProton) {
        registryMCtruth.fill(HIST("plus/MCtruthPr"), part.pt(), part.eta());
        registryMCtruth.fill(HIST("plus/MCtruthPrPt"), part.pt());
      }

      if (pdgParticle->Charge() < 0.0) {
        registryMCtruth.fill(HIST("minus/MCtruthAllPt"), part.pt());
      }
      if (pdgCode == kPiMinus) {
        registryMCtruth.fill(HIST("minus/MCtruthPi"), part.pt(), part.eta());
        registryMCtruth.fill(HIST("minus/MCtruthPiPt"), part.pt());
      }
      if (pdgCode == kProtonBar) {
        registryMCtruth.fill(HIST("minus/MCtruthPr"), part.pt(), part.eta());
        registryMCtruth.fill(HIST("minus/MCtruthPrPt"), part.pt());
      }
    }
  }

  PROCESS_SWITCH(FemtoUniversePairTaskTrackV0Extended, processMCTruth, "Process MC truth data", false);

  void processPairFractions(FilteredFDCollisions const& cols, FemtoRecoParticles const& parts, aod::FdMCParticles const& mcparts)
  {
    ColumnBinningPolicy<aod::collision::PosZ, aod::femtouniversecollision::MultNtr> colBinningMult{{confVtxBins, confMultBins}, true};
    ColumnBinningPolicy<aod::collision::PosZ, aod::femtouniversecollision::MultV0M> colBinningCent{{confVtxBins, confMultBins}, true};

    float v0DaughPtLowTable[3][2] = {{ConfV0Selection.confLPtChildProton, ConfV0Selection.confLPtChildPion}, {ConfV0Selection.confLPtChildPion, ConfV0Selection.confLPtChildProton}, {ConfV0Selection.confLPtChildPion, ConfV0Selection.confLPtChildPion}};
    float v0DaughPtHighTable[3][2] = {{ConfV0Selection.confHPtChildProton, ConfV0Selection.confHPtChildPion}, {ConfV0Selection.confHPtChildPion, ConfV0Selection.confHPtChildProton}, {ConfV0Selection.confHPtChildPion, ConfV0Selection.confHPtChildPion}};

    auto mixedCollProcessFunc = [&](auto& collision1, auto& collision2) -> void {
      auto groupPartsOne = partsOneMCRecoFull->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);
      auto groupPartsTwo = partsTwoMCRecoFull->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision2.globalIndex(), cache);

      const auto& magFieldTesla1 = collision1.magField();
      const auto& magFieldTesla2 = collision2.magField();

      if (magFieldTesla1 != magFieldTesla2) {
        return;
      }

      for (const auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(groupPartsOne, groupPartsTwo))) {
        // Lambda invariant mass cut
        if (!invMLambda(p2.mLambda(), p2.mAntiLambda(), ConfV0Selection.confV0Type1))
          continue;
        /// PID using stored binned nsigma
        if (!isParticleCombined(p1, ConfTrkSelection.confTrackChoicePartOne))
          continue;

        const auto& posChild = parts.iteratorAt(p2.globalIndex() - 2 - parts.begin().globalIndex());
        const auto& negChild = parts.iteratorAt(p2.globalIndex() - 1 - parts.begin().globalIndex());
        if (posChild.pt() < v0DaughPtLowTable[ConfV0Selection.confV0Type1][0] || negChild.pt() < v0DaughPtLowTable[ConfV0Selection.confV0Type1][1] || posChild.pt() > v0DaughPtHighTable[ConfV0Selection.confV0Type1][0] || negChild.pt() > v0DaughPtHighTable[ConfV0Selection.confV0Type1][1]) {
          continue;
        }

        /// Daughters that do not pass this condition are not selected
        if (!isParticleTPC(posChild, V0ChildTable[ConfV0Selection.confV0Type1][0]) || !isParticleTPC(negChild, V0ChildTable[ConfV0Selection.confV0Type1][1]))
          continue;
        if (!isParticleTOF(posChild, V0ChildTable[ConfV0Selection.confV0Type1][0]) || !isParticleTOF(negChild, V0ChildTable[ConfV0Selection.confV0Type1][1]))
          continue;

        // track cleaning
        if (!pairCleaner.isCleanPair(p1, p2, parts)) {
          continue;
        }
        if (ConfCPR.confIsCPR.value) {
          if (pairCloseRejection.isClosePair(p1, p2, parts, magFieldTesla1, femto_universe_container::EventType::mixed)) {
            continue;
          }
        }
        registryMCreco.fill(HIST("mothersReco/motherParticle"), p1.motherPDG(), p2.motherPDG());

        auto mcPartId1 = p1.fdMCParticleId();
        if (mcPartId1 == -1)
          continue;
        auto mcPartId2 = p2.fdMCParticleId();
        if (mcPartId2 == -1)
          continue;
        const auto& mcParticle1 = mcparts.iteratorAt(mcPartId1);
        const auto& mcParticle2 = mcparts.iteratorAt(mcPartId2);
        if (mcParticle1.pdgMCTruth() == ConfTrkSelection.confTrkPDGCodePartOne && mcParticle2.pdgMCTruth() == ConfV0Selection.confV0PDGCodePartTwo) {
          registryMCreco.fill(HIST("mothersReco/motherParticlePDGCheck"), p1.motherPDG(), p2.motherPDG());
        }
      }
    };

    for (const auto& [collision1, collision2] : soa::selfCombinations(colBinningCent, confNEventsMix, -1, cols, cols)) {
      mixedCollProcessFunc(collision1, collision2);
    }
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackV0Extended, processPairFractions, "Process MC data to obtain pair fractions", false);

  void processPairFractionsV0(FilteredFDCollisions const& cols, FemtoRecoParticles const& parts, aod::FdMCParticles const& mcparts)
  {
    ColumnBinningPolicy<aod::collision::PosZ, aod::femtouniversecollision::MultNtr> colBinningMult{{confVtxBins, confMultBins}, true};
    ColumnBinningPolicy<aod::collision::PosZ, aod::femtouniversecollision::MultV0M> colBinningCent{{confVtxBins, confMultBins}, true};

    float v0DaughPtLowTable[3][2] = {{ConfV0Selection.confLPtChildProton, ConfV0Selection.confLPtChildPion}, {ConfV0Selection.confLPtChildPion, ConfV0Selection.confLPtChildProton}, {ConfV0Selection.confLPtChildPion, ConfV0Selection.confLPtChildPion}};
    float v0DaughPtHighTable[3][2] = {{ConfV0Selection.confHPtChildProton, ConfV0Selection.confHPtChildPion}, {ConfV0Selection.confHPtChildPion, ConfV0Selection.confHPtChildProton}, {ConfV0Selection.confHPtChildPion, ConfV0Selection.confHPtChildPion}};

    auto mixedCollProcessFunc = [&](auto& collision1, auto& collision2) -> void {
      auto groupPartsOne = partsTwoMCRecoFull->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);
      auto groupPartsTwo = partsTwoMCRecoFull->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision2.globalIndex(), cache);

      const auto& magFieldTesla1 = collision1.magField();
      const auto& magFieldTesla2 = collision2.magField();

      if (magFieldTesla1 != magFieldTesla2) {
        return;
      }

      for (const auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(groupPartsOne, groupPartsTwo))) {
        // Lambda invariant mass cut for p1
        if (!invMLambda(p1.mLambda(), p1.mAntiLambda(), ConfV0Selection.confV0Type1)) {
          continue;
        }
        // Lambda invariant mass cut for p2
        if (!invMLambda(p2.mLambda(), p2.mAntiLambda(), ConfV0Selection.confV0Type2)) {
          continue;
        }

        const auto& posChild1 = parts.iteratorAt(p1.globalIndex() - 2 - parts.begin().globalIndex());
        const auto& negChild1 = parts.iteratorAt(p1.globalIndex() - 1 - parts.begin().globalIndex());
        if (posChild1.pt() < v0DaughPtLowTable[ConfV0Selection.confV0Type1][0] || negChild1.pt() < v0DaughPtLowTable[ConfV0Selection.confV0Type1][1] || posChild1.pt() > v0DaughPtHighTable[ConfV0Selection.confV0Type1][0] || negChild1.pt() > v0DaughPtHighTable[ConfV0Selection.confV0Type1][1]) {
          continue;
        }

        /// Daughters that do not pass this condition are not selected
        if (!isParticleTPC(posChild1, V0ChildTable[ConfV0Selection.confV0Type1][0]) || !isParticleTPC(negChild1, V0ChildTable[ConfV0Selection.confV0Type1][1]))
          continue;
        if (!isParticleTOF(posChild1, V0ChildTable[ConfV0Selection.confV0Type1][0]) || !isParticleTOF(negChild1, V0ChildTable[ConfV0Selection.confV0Type1][1]))
          continue;

        const auto& posChild2 = parts.iteratorAt(p2.globalIndex() - 2 - parts.begin().globalIndex());
        const auto& negChild2 = parts.iteratorAt(p2.globalIndex() - 1 - parts.begin().globalIndex());
        if (posChild2.pt() < v0DaughPtLowTable[ConfV0Selection.confV0Type2][0] || negChild2.pt() < v0DaughPtLowTable[ConfV0Selection.confV0Type2][1] || posChild2.pt() > v0DaughPtHighTable[ConfV0Selection.confV0Type2][0] || negChild2.pt() > v0DaughPtHighTable[ConfV0Selection.confV0Type2][1]) {
          continue;
        }

        /// Daughters that do not pass this condition are not selected
        if (!isParticleTPC(posChild2, V0ChildTable[ConfV0Selection.confV0Type2][0]) || !isParticleTPC(negChild2, V0ChildTable[ConfV0Selection.confV0Type2][1]))
          continue;
        if (!isParticleTOF(negChild2, V0ChildTable[ConfV0Selection.confV0Type1][0]) || !isParticleTOF(negChild2, V0ChildTable[ConfV0Selection.confV0Type1][1]))
          continue;

        // track cleaning
        if (!pairCleanerV0.isCleanPair(p1, p2, parts)) {
          continue;
        }
        if (ConfCPR.confIsCPR.value) {
          if (pairCloseRejectionV0.isClosePair(p1, p2, parts, magFieldTesla1, femto_universe_container::EventType::mixed)) {
            continue;
          }
        }

        registryMCreco.fill(HIST("mothersReco/motherParticle"), p1.motherPDG(), p2.motherPDG());
        auto mcPartId1 = p1.fdMCParticleId();
        if (mcPartId1 == -1)
          continue;
        auto mcPartId2 = p2.fdMCParticleId();
        if (mcPartId2 == -1)
          continue;
        const auto& mcParticle1 = mcparts.iteratorAt(mcPartId1);
        const auto& mcParticle2 = mcparts.iteratorAt(mcPartId2);
        if ((ConfV0Selection.confV0Type1 == 0 && mcParticle1.pdgMCTruth() != kLambda0) || (ConfV0Selection.confV0Type1 == 1 && mcParticle1.pdgMCTruth() != kLambda0Bar))
          continue;
        if ((ConfV0Selection.confV0Type2 == 0 && mcParticle2.pdgMCTruth() != kLambda0) || (ConfV0Selection.confV0Type2 == 1 && mcParticle2.pdgMCTruth() != kLambda0Bar))
          continue;
        registryMCreco.fill(HIST("mothersReco/motherParticlePDGCheck"), p1.motherPDG(), p2.motherPDG());
      }
    };

    if (confUseCent) {
      for (const auto& [collision1, collision2] : soa::selfCombinations(colBinningCent, confNEventsMix, -1, cols, cols)) {
        mixedCollProcessFunc(collision1, collision2);
      }
    } else {
      for (const auto& [collision1, collision2] : soa::selfCombinations(colBinningMult, confNEventsMix, -1, cols, cols)) {
        mixedCollProcessFunc(collision1, collision2);
      }
    }
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackV0Extended, processPairFractionsV0, "Process MC data to obtain pair fractions for V0V0 pairs", false);

  void processPairFractionsMCTruth(FilteredFDCollisions const& cols, FemtoFullParticles const& /*parts*/)
  {
    ColumnBinningPolicy<aod::collision::PosZ, aod::femtouniversecollision::MultNtr> colBinningMult{{confVtxBins, confMultBins}, true};
    ColumnBinningPolicy<aod::collision::PosZ, aod::femtouniversecollision::MultV0M> colBinningCent{{confVtxBins, confMultBins}, true};

    auto mixedCollProcessFunc = [&](auto& collision1, auto& collision2) -> void {
      auto groupPartsOne = partsOneMCFull->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);
      auto groupPartsTwo = partsTwoMCFull->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision2.globalIndex(), cache);

      for (const auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(groupPartsOne, groupPartsTwo))) {
        int pdgCode1 = static_cast<int>(p1.pidCut());
        int pdgCode2 = static_cast<int>(p2.pidCut());

        if (pdgCode1 != ConfTrkSelection.confTrkPDGCodePartOne)
          continue;
        if (pdgCode2 != ConfV0Selection.confV0PDGCodePartTwo)
          continue;

        registryMCtruth.fill(HIST("mothersTruth/motherParticle"), p1.tempFitVar(), p2.tempFitVar());
      }
    };

    if (confUseCent) {
      for (const auto& [collision1, collision2] : soa::selfCombinations(colBinningCent, confNEventsMix, -1, cols, cols)) {
        mixedCollProcessFunc(collision1, collision2);
      }
    } else {
      for (const auto& [collision1, collision2] : soa::selfCombinations(colBinningMult, confNEventsMix, -1, cols, cols)) {
        mixedCollProcessFunc(collision1, collision2);
      }
    }
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackV0Extended, processPairFractionsMCTruth, "Process MC data to obtain pair fractions for MC truth pairs", false);

  void processPairFractionsMCTruthV0(FilteredFDCollisions const& cols, FemtoFullParticles const& /*parts*/)
  {
    ColumnBinningPolicy<aod::collision::PosZ, aod::femtouniversecollision::MultNtr> colBinningMult{{confVtxBins, confMultBins}, true};
    ColumnBinningPolicy<aod::collision::PosZ, aod::femtouniversecollision::MultV0M> colBinningCent{{confVtxBins, confMultBins}, true};

    auto mixedCollProcessFunc = [&](auto& collision1, auto& collision2) -> void {
      auto groupPartsOne = partsTwoMCFull->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);
      auto groupPartsTwo = partsTwoMCFull->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision2.globalIndex(), cache);

      for (const auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(groupPartsOne, groupPartsTwo))) {
        int pdgCode1 = static_cast<int>(p1.pidCut());
        if ((ConfV0Selection.confV0Type1 == 0 && pdgCode1 != kLambda0) || (ConfV0Selection.confV0Type1 == 1 && pdgCode1 != kLambda0Bar))
          continue;
        int pdgCode2 = static_cast<int>(p2.pidCut());
        if ((ConfV0Selection.confV0Type2 == 0 && pdgCode2 != kLambda0) || (ConfV0Selection.confV0Type2 == 1 && pdgCode2 != kLambda0Bar))
          continue;

        registryMCtruth.fill(HIST("mothersTruth/motherParticle"), p1.tempFitVar(), p2.tempFitVar());
      }
    };

    if (confUseCent) {
      for (const auto& [collision1, collision2] : soa::selfCombinations(colBinningCent, confNEventsMix, -1, cols, cols)) {
        mixedCollProcessFunc(collision1, collision2);
      }
    } else {
      for (const auto& [collision1, collision2] : soa::selfCombinations(colBinningMult, confNEventsMix, -1, cols, cols)) {
        mixedCollProcessFunc(collision1, collision2);
      }
    }
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackV0Extended, processPairFractionsMCTruthV0, "Process MC data to obtain pair fractions for V0V0 MC truth pairs", false);

  void processMCReco(FemtoRecoParticles const& parts, aod::FdMCParticles const& mcparts)
  {
    for (const auto& part : parts) {
      auto mcPartId = part.fdMCParticleId();
      if (mcPartId == -1)
        continue; // no MC particle
      const auto& mcpart = mcparts.iteratorAt(mcPartId);
      //
      if (part.partType() == aod::femtouniverseparticle::ParticleType::kV0) {
        if (mcpart.pdgMCTruth() == kLambda0) {
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
        } else if (mcpart.pdgMCTruth() == kLambda0Bar) {
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
          if (mcpart.pdgMCTruth() == kPiPlus && isNSigmaCombined(part.p(), aod::pidtpc_tiny::binning::unPackInTable(part.tpcNSigmaStorePi()), aod::pidtof_tiny::binning::unPackInTable(part.tofNSigmaStorePi()))) {
            registryMCreco.fill(HIST("plus/MCrecoPi"), mcpart.pt(), mcpart.eta());
            registryMCreco.fill(HIST("plus/MCrecoPiPt"), mcpart.pt());
          } else if (mcpart.pdgMCTruth() == kProton && isNSigmaCombined(part.p(), aod::pidtpc_tiny::binning::unPackInTable(part.tpcNSigmaStorePr()), aod::pidtof_tiny::binning::unPackInTable(part.tofNSigmaStorePr()))) {
            registryMCreco.fill(HIST("plus/MCrecoPr"), mcpart.pt(), mcpart.eta());
            registryMCreco.fill(HIST("plus/MCrecoPrPt"), mcpart.pt());
          }
        }

        if (part.sign() < 0) {
          registryMCreco.fill(HIST("minus/MCrecoAllPt"), mcpart.pt());
          if (mcpart.pdgMCTruth() == kPiMinus && isNSigmaCombined(part.p(), aod::pidtpc_tiny::binning::unPackInTable(part.tpcNSigmaStorePi()), aod::pidtof_tiny::binning::unPackInTable(part.tofNSigmaStorePi()))) {
            registryMCreco.fill(HIST("minus/MCrecoPi"), mcpart.pt(), mcpart.eta());
            registryMCreco.fill(HIST("minus/MCrecoPiPt"), mcpart.pt());
          } else if (mcpart.pdgMCTruth() == kProtonBar && isNSigmaCombined(part.p(), aod::pidtpc_tiny::binning::unPackInTable(part.tpcNSigmaStorePr()), aod::pidtof_tiny::binning::unPackInTable(part.tofNSigmaStorePr()))) {
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
