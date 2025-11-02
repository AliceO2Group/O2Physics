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

/// \file femtoUniversePairTaskTrackTrackSpherHarMultKtExtended.cxx
/// \brief Tasks that reads the track tables used for the pairing and builds pairs of two tracks and compute relative pair-momentum in three dimesnions
/// \remark This file is inherited from ~/FemtoUniverse/Tasks/femtoUniversePairTaskTrackTrack3DMultKtExtended.cxx on 17/06/2024
/// \author Pritam Chakraborty, WUT Warsaw, pritam.chakraborty@pw.edu.pl

#include "PWGCF/FemtoUniverse/Core/FemtoUniverseDetaDphiStar.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseEventHisto.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseMath.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniversePairCleaner.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniversePairSHCentMultKt.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseParticleHisto.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseSHContainer.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseTrackSelection.h"
#include "PWGCF/FemtoUniverse/Core/femtoUtils.h"
#include "PWGCF/FemtoUniverse/DataModel/FemtoDerived.h"

#include "Common/DataModel/PIDResponse.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/StepTHn.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/PID.h"

#include "TRandom2.h"

#include <string>
#include <vector>

using namespace o2;
using namespace o2::analysis::femto_universe;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

namespace
{
static constexpr int nPart = 2;
static constexpr int nCuts = 5;
static const std::vector<std::string> partNames{"PartOne", "PartTwo"};
static const std::vector<std::string> cutNames{"MaxPt", "PIDthr", "nSigmaTPC", "nSigmaTPCTOF", "MaxP"};
static const float cutsTable[nPart][nCuts]{
  {4.05f, 1.f, 3.f, 3.f, 100.f},
  {4.05f, 1.f, 3.f, 3.f, 100.f}};
} // namespace

struct femtoUniversePairTaskTrackTrackSpherHarMultKtExtended {

  Service<o2::framework::O2DatabasePDG> pdg;

  /// Particle selection part

  /// Table for both particles
  struct : o2::framework::ConfigurableGroup {
    Configurable<float> ConfNsigmaCombined{"ConfNsigmaCombined", 3.0f, "TPC and TOF Pion Sigma (combined) for momentum > ConfTOFPtMin"};
    Configurable<float> ConfNsigmaTPC{"ConfNsigmaTPC", 3.0f, "TPC Pion Sigma for momentum < ConfTOFPtMin"};
    Configurable<bool> ConfIsElReject{"ConfIsElReject", false, "Is electron rejection activated"};
    Configurable<float> ConfNsigmaTPCElRejectMin{"ConfNsigmaTPCElRejectMin", 2.0f, "TPC Electron SigmaMin for momentum < ConfTOFPtMin"};
    Configurable<float> ConfNsigmaTPCElRejectMax{"ConfNsigmaTPCElRejectMax", 2.0f, "TPC Electron SigmaMax for momentum < ConfTOFPtMin"};
    Configurable<float> ConfTOFPtMin{"ConfTOFPtMin", 0.5f, "Min. Pt for which TOF is required for PID."};
    Configurable<float> ConfEtaMax{"ConfEtaMax", 0.8f, "Higher limit for |Eta| (the same for both particles)"};

    Configurable<LabeledArray<float>> ConfCutTable{"ConfCutTable", {cutsTable[0], nPart, nCuts, partNames, cutNames}, "Particle selections"};
    Configurable<int> ConfNspecies{"ConfNspecies", 2, "Number of particle spieces with PID info"};
    Configurable<bool> ConfIsMC{"ConfIsMC", false, "Enable additional Histogramms in the case of a MonteCarlo Run"};
    Configurable<std::vector<float>> ConfTrkPIDnSigmaMax{"ConfTrkPIDnSigmaMax", std::vector<float>{4.f, 3.f, 2.f}, "This configurable needs to be the same as the one used in the producer task"};
    Configurable<bool> ConfUse3D{"ConfUse3D", false, "Enable three dimensional histogramms (to be used only for analysis with high statistics): k* vs mT vs multiplicity"};
    Configurable<int> ConfPhiBins{"ConfPhiBins", 29, "Number of phi bins in deta dphi"};
    Configurable<int> ConfEtaBins{"ConfEtaBins", 29, "Number of eta bins in deta dphi"};
    Configurable<bool> ConfIsCPR{"ConfIsCPR", true, "Close Pair Rejection"};
    Configurable<bool> ConfCPRPlotPerRadii{"ConfCPRPlotPerRadii", false, "Plot CPR per radii"};
    Configurable<float> ConfCPRdeltaPhiCutMax{"ConfCPRdeltaPhiCutMax", 0.0, "Delta Phi max cut for Close Pair Rejection"};
    Configurable<float> ConfCPRdeltaPhiCutMin{"ConfCPRdeltaPhiCutMin", 0.0, "Delta Phi min cut for Close Pair Rejection"};
    Configurable<float> ConfCPRdeltaEtaCutMax{"ConfCPRdeltaEtaCutMax", 0.0, "Delta Eta max cut for Close Pair Rejection"};
    Configurable<float> ConfCPRdeltaEtaCutMin{"ConfCPRdeltaEtaCutMin", 0.0, "Delta Eta min cut for Close Pair Rejection"};
    Configurable<float> ConfCPRChosenRadii{"ConfCPRChosenRadii", 0.80, "Delta Eta cut for Close Pair Rejection"};
    Configurable<bool> confUseCCImCut{"confUseCCImCut", false, "Fill SH within specific quadrants of qout-qside"};
    Configurable<bool> confUse1stand3rd{"confUse1stand3rd", false, "Use first and third quadrants of qout-qside"};
    Configurable<bool> confUse2ndand4th{"confUse2ndand4th", false, "Use second and fourth quadrants of qout-qside"};
  } twotracksconfigs;

  using FemtoFullParticles = soa::Join<aod::FDParticles, aod::FDExtParticles>;
  // Filters for selecting particles (both p1 and p2)
  Filter trackAdditionalfilter = (nabs(aod::femtouniverseparticle::eta) < twotracksconfigs.ConfEtaMax); // example filtering on configurable
  using FilteredFemtoFullParticles = soa::Filtered<FemtoFullParticles>;
  // using FilteredFemtoFullParticles = FemtoFullParticles; //if no filtering is applied uncomment this option
  using FemtoRecoParticles = soa::Join<aod::FDParticles, aod::FDExtParticles, aod::FDMCLabels>;
  using FilteredFemtoRecoParticles = soa::Filtered<FemtoRecoParticles>;

  SliceCache cache;
  Preslice<FilteredFemtoFullParticles> perCol = aod::femtouniverseparticle::fdCollisionId;

  using FemtoTruthParticles = soa::Join<aod::FDParticles, aod::FDExtParticles, aod::FDMCLabels>;
  Preslice<FemtoTruthParticles> perColMCTruth = aod::femtouniverseparticle::fdCollisionId;

  /// Particle 1
  struct : o2::framework::ConfigurableGroup {
    Configurable<int> ConfPDGCodePartOne{"ConfPDGCodePartOne", 211, "Particle 1 - PDG code"};
    // Configurable<uint32_t> ConfCutPartOne{"ConfCutPartOne", 5542474, "Particle 1 - Selection bit from cutCulator"};
    Configurable<int> ConfPIDPartOne{"ConfPIDPartOne", 2, "Particle 1 - Read from cutCulator"}; // we also need the possibility to specify whether the bit is true/false ->std>>vector<std::pair<int, int>>int>>
    Configurable<float> ConfPtLowPart1{"ConfPtLowPart1", 0.14, "Lower limit for Pt for the first particle"};
    Configurable<float> ConfPtHighPart1{"ConfPtHighPart1", 1.5, "Higher limit for Pt for the first particle"};
    Configurable<int> ConfChargePart1{"ConfChargePart1", 1, "Particle 1 sign"};
  } trackonefilter;

  /// Partition for particle 1
  Partition<FilteredFemtoFullParticles> partsOne = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kTrack)) && aod::femtouniverseparticle::sign == as<int8_t>(trackonefilter.ConfChargePart1) && aod::femtouniverseparticle::pt < trackonefilter.ConfPtHighPart1 && aod::femtouniverseparticle::pt > trackonefilter.ConfPtLowPart1;
  Partition<FilteredFemtoRecoParticles> partsOneMC = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kTrack)) && aod::femtouniverseparticle::sign == as<int8_t>(trackonefilter.ConfChargePart1) && aod::femtouniverseparticle::pt < trackonefilter.ConfPtHighPart1 && aod::femtouniverseparticle::pt > trackonefilter.ConfPtLowPart1;
  Partition<FemtoTruthParticles> partsOneMCTruth = aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kMCTruthTrack);

  //

  /// Histogramming for particle 1
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kTrack, 1> trackHistoPartOne;

  /// Particle 2
  struct : o2::framework::ConfigurableGroup {
    Configurable<int> ConfPDGCodePartTwo{"ConfPDGCodePartTwo", 211, "Particle 2 - PDG code"};
    // Configurable<uint32_t> ConfCutPartTwo{"ConfCutPartTwo", 5542474, "Particle 2 - Selection bit"};
    Configurable<int> ConfPIDPartTwo{"ConfPIDPartTwo", 2, "Particle 2 - Read from cutCulator"}; // we also need the possibility to specify whether the bit is true/false ->std>>vector<std::pair<int, int>>

    Configurable<float> ConfPtLowPart2{"ConfPtLowPart2", 0.14, "Lower limit for Pt for the second particle"};
    Configurable<float> ConfPtHighPart2{"ConfPtHighPart2", 1.5, "Higher limit for Pt for the second particle"};
    Configurable<int> ConfChargePart2{"ConfChargePart2", -1, "Particle 2 sign"};
  } tracktwofilter;

  /// Partition for particle 2
  Partition<FilteredFemtoFullParticles> partsTwo = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kTrack)) && (aod::femtouniverseparticle::sign == as<int8_t>(tracktwofilter.ConfChargePart2)) && aod::femtouniverseparticle::pt < tracktwofilter.ConfPtHighPart2 && aod::femtouniverseparticle::pt > tracktwofilter.ConfPtLowPart2;
  Partition<FilteredFemtoRecoParticles> partsTwoMC = aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kTrack) && (aod::femtouniverseparticle::sign == as<int8_t>(tracktwofilter.ConfChargePart2)) && aod::femtouniverseparticle::pt < tracktwofilter.ConfPtHighPart2 && aod::femtouniverseparticle::pt > tracktwofilter.ConfPtLowPart2;
  Partition<FemtoTruthParticles> partsTwoMCTruth = aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kMCTruthTrack);

  /// Histogramming for particle 2
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kTrack, 2> trackHistoPartTwo;

  /// Histogramming for Event
  FemtoUniverseEventHisto eventHisto;

  /// The configurables need to be passed to an std::vector
  int vPIDPartOne, vPIDPartTwo;
  std::vector<float> kNsigma;

  /// Event part
  Configurable<float> ConfV0MLow{"ConfV0MLow", 0.0, "Lower limit for V0M multiplicity"};
  Configurable<float> ConfV0MHigh{"ConfV0MHigh", 25000.0, "Upper limit for V0M multiplicity"};
  Configurable<int> ConfTPCOccupancyLow{"ConfTPCOccupancyLow", 0, "Lower limit for TPC occupancy"};
  Configurable<int> ConfTPCOccupancyHigh{"ConfTPCOccupancyHigh", 500, "Higher limit for TPC occupancy"};
  Configurable<float> ConfIntRateLow{"ConfIntRateLow", 0.0, "Lower limit for interaction rate"};
  Configurable<float> ConfIntRateHigh{"ConfIntRateHigh", 10000.0, "Higher limit for interaction rate"};
  Configurable<bool> ConfIsCent{"ConfIsCent", true, "Condition to choose centrality of multiplicity for mixing"};

  Filter collfilterFDtable = (o2::aod::femtouniversecollision::multV0M > ConfV0MLow) && (o2::aod::femtouniversecollision::multV0M < ConfV0MHigh);
  Filter collfilterFDExttable = (o2::aod::femtouniversecollision::interactionRate > ConfIntRateLow) && (o2::aod::femtouniversecollision::interactionRate < ConfIntRateHigh) &&
                                (o2::aod::femtouniversecollision::occupancy >= ConfTPCOccupancyLow) && (o2::aod::femtouniversecollision::occupancy < ConfTPCOccupancyHigh);
  using FilteredFDCollisions = soa::Filtered<soa::Join<aod::FdCollisions, aod::FDExtCollisions>>;
  using FilteredFDCollision = FilteredFDCollisions::iterator;
  // Filter trackAdditionalfilter = (nabs(aod::femtouniverseparticle::eta) < twotracksconfigs.ConfEtaMax); // example filtering on configurable

  /// Particle part
  ConfigurableAxis ConfTempFitVarBins{"ConfDTempFitVarBins", {300, -0.15, 0.15}, "binning of the TempFitVar in the pT vs. TempFitVar plot"};
  ConfigurableAxis ConfTempFitVarpTBins{"ConfTempFitVarpTBins", {20, 0.5, 4.05}, "pT binning of the pT vs. TempFitVar plot"};

  /// Correlation part
  ConfigurableAxis ConfMultBinsCent{"ConfMultBinsCent", {VARIABLE_WIDTH, 0.0f, 4.0f, 8.0f, 12.0f, 16.0f, 20.0f, 24.0f, 28.0f, 32.0f, 36.0f, 40.0f, 44.0f, 48.0f, 52.0f, 56.0f, 60.0f, 64.0f, 68.0f, 72.0f, 76.0f, 80.0f, 84.0f, 88.0f, 92.0f, 96.0f, 100.0f, 200.0f, 99999.f}, "Mixing bins - centrality"}; // \todo to be obtained from the hash task
  ConfigurableAxis ConfMultBinsMult{"ConfMultBinsMult", {VARIABLE_WIDTH, 0.0f, 400.0f, 800.0f, 1200.0f, 1600.0f, 2000.0f, 2500.0f, 3000.0f, 3500.0f, 4000.0f, 4500.0f, 5000.0f, 6000.0f, 7000.0f, 8000.0f, 9000.0f, 10000.0f, 11000.0f, 12000.0f, 13000.0f, 14000.0f, 15000.0f, 16000.0f, 17000.0f, 18000.0f, 99999.f}, "Mixing bins - centrality"};
  ConfigurableAxis ConfMultKstarBins{"ConfMultKstarBins", {VARIABLE_WIDTH, 0.0f, 200.0f}, "Bins for kstar analysis in multiplicity or centrality bins (10 is maximum)"};
  ConfigurableAxis ConfKtKstarBins{"ConfKtKstarBins", {VARIABLE_WIDTH, 0.1f, 0.2f, 0.3f, 0.4f}, "Bins for kstar analysis in kT bins"};
  ConfigurableAxis ConfVtxBins{"ConfVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};

  ConfigurableAxis ConfmTBins3D{"ConfmTBins3D", {VARIABLE_WIDTH, 1.02f, 1.14f, 1.20f, 1.26f, 1.38f, 1.56f, 1.86f, 4.50f}, "mT Binning for the 3Dimensional plot: k* vs multiplicity vs mT (set <<twotracksconfigs.ConfUse3D>> to true in order to use)"};
  ConfigurableAxis ConfmultBins3D{"ConfmultBins3D", {VARIABLE_WIDTH, 0.0f, 20.0f, 30.0f, 40.0f, 99999.0f}, "multiplicity Binning for the 3Dimensional plot: k* vs multiplicity vs mT (set <<twotracksconfigs.ConfUse3D>> to true in order to use)"};

  ColumnBinningPolicy<aod::collision::PosZ, aod::femtouniversecollision::MultV0M> colBinningCent{{ConfVtxBins, ConfMultBinsCent}, true};
  ColumnBinningPolicy<aod::collision::PosZ, aod::femtouniversecollision::MultNtr> colBinningNtr{{ConfVtxBins, ConfMultBinsMult}, true};

  ConfigurableAxis ConfkstarBins{"ConfkstarBins", {60, 0.0, 0.3}, "binning kstar"};
  ConfigurableAxis ConfkTBins{"ConfkTBins", {150, 0., 9.}, "binning kT"};
  ConfigurableAxis ConfmTBins{"ConfmTBins", {225, 0., 7.5}, "binning mT"};
  Configurable<bool> ConfIsIden{"ConfIsIden", true, "Choosing identical or non-identical pairs"};
  Configurable<bool> ConfIsLCMS{"ConfIsLCMS", true, "Choosing LCMS or PRF"};
  Configurable<int> ConfNEventsMix{"ConfNEventsMix", 5, "Number of events for mixing"};
  Configurable<int> ConfLMax{"ConfLMax", 2, "Maximum value of l"};
  Configurable<bool> cfgProcessPM{"cfgProcessPM", false, "Process particles of the opposite charge"};
  Configurable<bool> cfgProcessPP{"cfgProcessPP", true, "Process particles of the same, positice charge"};
  Configurable<bool> cfgProcessMM{"cfgProcessMM", true, "Process particles of the same, positice charge"};
  Configurable<bool> cfgProcessMultBins{"cfgProcessMultBins", true, "Process kstar histograms in multiplicity bins (in multiplicity bins)"};
  Configurable<bool> cfgProcessKtBins{"cfgProcessKtBins", true, "Process kstar histograms in kT bins (if cfgProcessMultBins is set false, this will not be processed regardless this Configurable state)"};
  Configurable<bool> cfgProcessKtMt3DCF{"cfgProcessKtMt3DCF", false, "Process 3D histograms in kT and Mult bins"};
  Configurable<bool> ConfIsFillAngqLCMS{"ConfIsFillAngqLCMS", true, "Fill qLCMS vs dEta vs dPhi"};
  Configurable<float> confCPRDistMax{"confCPRDistMax", 0.0, "Max. radial seperation between two closed-pairs"};
  Configurable<float> confCPRFracMax{"confCPRFracMax", 0.0, "Max. allowed fraction bad to all TPC points of radial seperation between two closed-pairs"};
  Configurable<bool> confCPRIsAtITS{"confCPRIsAtITS", false, "Close Pair Rejection at ITS or TPC"};
  Configurable<bool> confCPRDphiAvgOrDist{"confCPRDphiAvgOrDist", true, "Close Pair Rejection by radial or angular seperation"};
  Configurable<bool> confIsCircularCut{"confIsCircularCut", true, "Close Pair Rejection by circular cut"};

  FemtoUniverseSHContainer<femto_universe_sh_container::EventType::same, femto_universe_sh_container::Observable::kstar> sameEventCont;
  FemtoUniverseSHContainer<femto_universe_sh_container::EventType::mixed, femto_universe_sh_container::Observable::kstar> mixedEventCont;

  FemtoUniverseSHContainer<femto_universe_sh_container::EventType::same, femto_universe_sh_container::Observable::kstar> sameEventContPP;
  FemtoUniverseSHContainer<femto_universe_sh_container::EventType::mixed, femto_universe_sh_container::Observable::kstar> mixedEventContPP;

  FemtoUniverseSHContainer<femto_universe_sh_container::EventType::same, femto_universe_sh_container::Observable::kstar> sameEventContMM;
  FemtoUniverseSHContainer<femto_universe_sh_container::EventType::mixed, femto_universe_sh_container::Observable::kstar> mixedEventContMM;

  FemtoUniverseContainer<femto_universe_container::EventType::same, femto_universe_container::Observable::kstar> sameEventCont1D;

  FemtoUniversePairCleaner<aod::femtouniverseparticle::ParticleType::kTrack, aod::femtouniverseparticle::ParticleType::kTrack> pairCleaner;
  FemtoUniverseDetaDphiStar<aod::femtouniverseparticle::ParticleType::kTrack, aod::femtouniverseparticle::ParticleType::kTrack> pairCloseRejection;
  FemtoUniverseTrackSelection trackCuts;

  PairSHCentMultKt<femto_universe_sh_container::EventType::same, femto_universe_sh_container::Observable::kstar> sameEventMultCont;
  PairSHCentMultKt<femto_universe_sh_container::EventType::mixed, femto_universe_sh_container::Observable::kstar> mixedEventMultCont;

  PairSHCentMultKt<femto_universe_sh_container::EventType::same, femto_universe_sh_container::Observable::kstar> sameEventMultContPP;
  PairSHCentMultKt<femto_universe_sh_container::EventType::mixed, femto_universe_sh_container::Observable::kstar> mixedEventMultContPP;

  PairSHCentMultKt<femto_universe_sh_container::EventType::same, femto_universe_sh_container::Observable::kstar> sameEventMultContMM;
  PairSHCentMultKt<femto_universe_sh_container::EventType::mixed, femto_universe_sh_container::Observable::kstar> mixedEventMultContMM;

  float mass1 = -1;
  float mass2 = -1;

  /// Histogram output
  HistogramRegistry qaRegistry{"TrackQA", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry resultRegistry{"Correlations", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry resultRegistry1D{"Correlations1D", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry resultRegistryPM{"CorrelationsPM", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry resultRegistryPP{"CorrelationsPP", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry resultRegistryMM{"CorrelationsMM", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry MixQaRegistry{"MixQaRegistry", {}, OutputObjHandlingPolicy::AnalysisObject};

  HistogramRegistry SameMultRegistryPM{"SameMultRegistryPM", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry MixedMultRegistryPM{"MixedMultRegistryPM", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  HistogramRegistry SameMultRegistryPP{"SameMultRegistryPP", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry MixedMultRegistryPP{"MixedMultRegistryPP", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  HistogramRegistry SameMultRegistryMM{"SameMultRegistryMM", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry MixedMultRegistryMM{"MixedMultRegistryMM", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  TRandom2* randgen;

  // PID for protons
  bool IsProtonNSigma(float mom, float nsigmaTPCPr, float nsigmaTOFPr) // previous version from: https://github.com/alisw/AliPhysics/blob/master/PWGCF/FEMTOSCOPY/AliFemtoUser/AliFemtoMJTrackCut.cxx
  {
    //|nsigma_TPC| < 3 for p < 0.5 GeV/c
    //|nsigma_combined| < 3 for p > 0.5

    // using configurables:
    // ConfTOFPtMin - momentum value when we start using TOF; set to 1000 if TOF not needed
    // ConfNsigmaTPC -> TPC Sigma for momentum < 0.5
    // ConfNsigmaCombined -> TPC and TOF Sigma (combined) for momentum > 0.5

    if (mom < twotracksconfigs.ConfTOFPtMin) {
      if (std::abs(nsigmaTPCPr) < twotracksconfigs.ConfNsigmaTPC) {
        return true;
      } else {
        return false;
      }
    } else {
      if (std::hypot(nsigmaTOFPr, nsigmaTPCPr) < twotracksconfigs.ConfNsigmaCombined) {
        return true;
      } else {
        return false;
      }
    }
    return false;
  }

  bool IsKaonNSigma(float mom, float nsigmaTPCK, float nsigmaTOFK)
  {
    if (mom < 0.3) { // 0.0-0.3
      if (std::abs(nsigmaTPCK) < 3.0) {
        return true;
      } else {
        return false;
      }
    } else if (mom < 0.45) { // 0.30 - 0.45
      if (std::abs(nsigmaTPCK) < 2.0) {
        return true;
      } else {
        return false;
      }
    } else if (mom < 0.55) { // 0.45-0.55
      if (std::abs(nsigmaTPCK) < 1.0) {
        return true;
      } else {
        return false;
      }
    } else if (mom < 1.5) { // 0.55-1.5 (now we use TPC and TOF)
      if ((std::abs(nsigmaTOFK) < 3.0) && (std::abs(nsigmaTPCK) < 3.0)) {
        {
          return true;
        }
      } else {
        return false;
      }
    } else if (mom > 1.5) { // 1.5 -
      if ((std::abs(nsigmaTOFK) < 2.0) && (std::abs(nsigmaTPCK) < 3.0)) {
        return true;
      } else {
        return false;
      }
    } else {
      return false;
    }
  }

  bool IsPionNSigma(float mom, float nsigmaTPCPi, float nsigmaTOFPi, float nsigmaTPCElReject)
  {
    //|nsigma_TPC| < 3 for p < 0.5 GeV/c
    //|nsigma_combined| < 3 for p > 0.5

    // using configurables:
    // ConfTOFPtMin - momentum value when we start using TOF; set to 1000 if TOF not needed
    // ConfNsigmaTPC -> TPC Sigma for momentum < 0.5
    // ConfNsigmaCombined -> TPC and TOF Pion Sigma (combined) for momentum > 0.5
    if (true) {
      if (mom < twotracksconfigs.ConfTOFPtMin) {
        if (twotracksconfigs.ConfIsElReject) {
          if ((std::abs(nsigmaTPCPi) < twotracksconfigs.ConfNsigmaTPC) && (nsigmaTPCElReject < twotracksconfigs.ConfNsigmaTPCElRejectMin || nsigmaTPCElReject > twotracksconfigs.ConfNsigmaTPCElRejectMax)) {
            return true;
          } else {
            return false;
          }
        } else {
          if ((std::abs(nsigmaTPCPi) < twotracksconfigs.ConfNsigmaTPC)) {
            return true;
          } else {
            return false;
          }
        }
      } else {
        if (std::hypot(nsigmaTOFPi, nsigmaTPCPi) < twotracksconfigs.ConfNsigmaCombined) {
          return true;
        } else {
          return false;
        }
      }
    }
    return false;
  }

  bool IsParticleNSigma(int8_t particle_number, float mom, float nsigmaTPCPr, float nsigmaTOFPr, float nsigmaTPCPi, float nsigmaTOFPi, float nsigmaTPCK, float nsigmaTOFK, float nsigmaTPCElReject)
  {
    if (particle_number == 1) {
      switch (trackonefilter.ConfPDGCodePartOne) {
        case 2212:  // Proton
        case -2212: // Antiproton
          return IsProtonNSigma(mom, nsigmaTPCPr, nsigmaTOFPr);
          break;
        case 211:  // Pion+
        case -211: // Pion-
          return IsPionNSigma(mom, nsigmaTPCPi, nsigmaTOFPi, nsigmaTPCElReject);
          break;
        case 321:  // Kaon+
        case -321: // Kaon-
          return IsKaonNSigma(mom, nsigmaTPCK, nsigmaTOFK);
          break;
        default:
          return false;
      }
      return false;
    } else if (particle_number == 2) {
      switch (tracktwofilter.ConfPDGCodePartTwo) {
        case 2212:  // Proton
        case -2212: // Antiproton
          return IsProtonNSigma(mom, nsigmaTPCPr, nsigmaTOFPr);
          break;
        case 211:  // Pion+
        case -211: // Pion-
          return IsPionNSigma(mom, nsigmaTPCPi, nsigmaTOFPi, nsigmaTPCElReject);
          break;
        case 321:  // Kaon+
        case -321: // Kaon-
          return IsKaonNSigma(mom, nsigmaTPCK, nsigmaTOFK);
          break;
        default:
          return false;
      }
      return false;
    } else {
      LOGF(fatal, "Wrong number of particle chosen! It should be 1 or 2. It is -> %d", particle_number);
    }
    return false;
  }

  void init(InitContext&)
  {
    eventHisto.init(&qaRegistry);
    trackHistoPartOne.init(&qaRegistry, ConfTempFitVarpTBins, ConfTempFitVarBins, twotracksconfigs.ConfIsMC, trackonefilter.ConfPDGCodePartOne, true);

    trackHistoPartTwo.init(&qaRegistry, ConfTempFitVarpTBins, ConfTempFitVarBins, twotracksconfigs.ConfIsMC, tracktwofilter.ConfPDGCodePartTwo, true);

    MixQaRegistry.add("MixingQA/hSECollisionBins", ";bin;Entries", kTH1F, {{120, -0.5, 119.5}});
    MixQaRegistry.add("MixingQA/hMECollisionBins", ";bin;Entries", kTH1F, {{120, -0.5, 119.5}});

    mass1 = pdg->Mass(trackonefilter.ConfPDGCodePartOne);
    mass2 = pdg->Mass(tracktwofilter.ConfPDGCodePartTwo);

    if (cfgProcessPM) {
      if (!cfgProcessKtMt3DCF) {
        sameEventCont.init(&resultRegistryPM, ConfkstarBins, ConfLMax);
        mixedEventCont.init(&resultRegistryPM, ConfkstarBins, ConfLMax);
        sameEventCont.setPDGCodes(trackonefilter.ConfPDGCodePartOne, tracktwofilter.ConfPDGCodePartTwo);
        mixedEventCont.setPDGCodes(trackonefilter.ConfPDGCodePartOne, tracktwofilter.ConfPDGCodePartTwo);
      } else {
        sameEventMultCont.init(&SameMultRegistryPM, ConfkstarBins, ConfMultKstarBins, ConfKtKstarBins, ConfLMax);
        mixedEventMultCont.init(&MixedMultRegistryPM, ConfkstarBins, ConfMultKstarBins, ConfKtKstarBins, ConfLMax);
      }
    }

    if (cfgProcessPP) {
      if (!cfgProcessKtMt3DCF) {
        sameEventContPP.init(&resultRegistryPP, ConfkstarBins, ConfLMax);
        mixedEventContPP.init(&resultRegistryPP, ConfkstarBins, ConfLMax);
        sameEventContPP.setPDGCodes(trackonefilter.ConfPDGCodePartOne, tracktwofilter.ConfPDGCodePartTwo);
        mixedEventContPP.setPDGCodes(trackonefilter.ConfPDGCodePartOne, tracktwofilter.ConfPDGCodePartTwo);
      } else {
        sameEventMultContPP.init(&SameMultRegistryPP, ConfkstarBins, ConfMultKstarBins, ConfKtKstarBins, ConfLMax);
        mixedEventMultContPP.init(&MixedMultRegistryPP, ConfkstarBins, ConfMultKstarBins, ConfKtKstarBins, ConfLMax);
      }
      sameEventCont1D.init(&resultRegistry1D, ConfkstarBins, ConfMultBinsCent, ConfkTBins, ConfmTBins, ConfmultBins3D, ConfmTBins3D, twotracksconfigs.ConfEtaBins, twotracksconfigs.ConfPhiBins, twotracksconfigs.ConfIsMC, twotracksconfigs.ConfUse3D);
      sameEventCont1D.setPDGCodes(trackonefilter.ConfPDGCodePartOne, tracktwofilter.ConfPDGCodePartTwo);
    }

    if (cfgProcessMM) {
      if (!cfgProcessKtMt3DCF) {
        sameEventContMM.init(&resultRegistryMM, ConfkstarBins, ConfLMax);
        mixedEventContMM.init(&resultRegistryMM, ConfkstarBins, ConfLMax);
        sameEventContMM.setPDGCodes(trackonefilter.ConfPDGCodePartOne, tracktwofilter.ConfPDGCodePartTwo);
        mixedEventContMM.setPDGCodes(trackonefilter.ConfPDGCodePartOne, tracktwofilter.ConfPDGCodePartTwo);
      } else {
        sameEventMultContMM.init(&SameMultRegistryMM, ConfkstarBins, ConfMultKstarBins, ConfKtKstarBins, ConfLMax);
        mixedEventMultContMM.init(&MixedMultRegistryMM, ConfkstarBins, ConfMultKstarBins, ConfKtKstarBins, ConfLMax);
      }
      sameEventCont1D.init(&resultRegistry1D, ConfkstarBins, ConfMultBinsCent, ConfkTBins, ConfmTBins, ConfmultBins3D, ConfmTBins3D, twotracksconfigs.ConfEtaBins, twotracksconfigs.ConfPhiBins, twotracksconfigs.ConfIsMC, twotracksconfigs.ConfUse3D);
      sameEventCont1D.setPDGCodes(trackonefilter.ConfPDGCodePartOne, tracktwofilter.ConfPDGCodePartTwo);
    }

    pairCleaner.init(&qaRegistry);
    if (twotracksconfigs.ConfIsCPR.value) {
      pairCloseRejection.init(&resultRegistry, &qaRegistry, twotracksconfigs.ConfCPRdeltaPhiCutMin.value, twotracksconfigs.ConfCPRdeltaPhiCutMax.value, twotracksconfigs.ConfCPRdeltaEtaCutMin.value, twotracksconfigs.ConfCPRdeltaEtaCutMax.value, twotracksconfigs.ConfCPRChosenRadii.value, twotracksconfigs.ConfCPRPlotPerRadii.value);
      pairCloseRejection.init_kT(&resultRegistry, ConfKtKstarBins);
    }

    vPIDPartOne = trackonefilter.ConfPIDPartOne.value;
    vPIDPartTwo = tracktwofilter.ConfPIDPartTwo.value;
    kNsigma = twotracksconfigs.ConfTrkPIDnSigmaMax.value;
  }

  template <typename CollisionType>
  void fillCollision(CollisionType col, bool IsCent)
  {
    if (IsCent) {
      MixQaRegistry.fill(HIST("MixingQA/hSECollisionBins"), colBinningCent.getBin({col.posZ(), col.multV0M()}));
    } else {
      MixQaRegistry.fill(HIST("MixingQA/hSECollisionBins"), colBinningNtr.getBin({col.posZ(), col.multNtr()}));
    }
    eventHisto.fillQA(col);
  }

  /// This function processes the same event and takes care of all the histogramming
  /// \todo the trivial loops over the tracks should be factored out since they will be common to all combinations of T-T, T-V0, V0-V0, ...
  /// @tparam PartitionType
  /// @tparam PartType
  /// @tparam isMC: enables Monte Carlo truth specific histograms
  /// @param groupPartsOne partition for the first particle passed by the process function
  /// @param groupPartsTwo partition for the second particle passed by the process function
  /// @param parts femtoUniverseParticles table (in case of Monte Carlo joined with FemtoUniverseMCLabels)
  /// @param magFieldTesla magnetic field of the collision
  /// @param multCol multiplicity of the collision
  template <bool isMC, typename PartitionType, typename PartType>
  void doSameEvent(PartitionType groupPartsOne, PartitionType groupPartsTwo, PartType parts, float magFieldTesla, int multCol, int ContType, bool fillQA)
  {

    /// Histogramming same event
    if ((ContType == 1 || ContType == 2) && fillQA) {
      for (const auto& part : groupPartsOne) {
        if (!IsParticleNSigma((int8_t)1, part.p(), trackCuts.getNsigmaTPC(part, o2::track::PID::Proton), trackCuts.getNsigmaTOF(part, o2::track::PID::Proton), trackCuts.getNsigmaTPC(part, o2::track::PID::Pion), trackCuts.getNsigmaTOF(part, o2::track::PID::Pion), trackCuts.getNsigmaTPC(part, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(part, o2::track::PID::Kaon), trackCuts.getNsigmaTPC(part, o2::track::PID::Electron))) {
          continue;
        }
        trackHistoPartOne.fillQA<isMC, true>(part);
        trackHistoPartOne.fillQAMisIden<isMC, true>(part, trackonefilter.ConfPDGCodePartOne);
      }
    }

    if ((ContType == 1 || ContType == 3) && fillQA) {
      for (const auto& part : groupPartsTwo) {
        if (!IsParticleNSigma((int8_t)2, part.p(), trackCuts.getNsigmaTPC(part, o2::track::PID::Proton), trackCuts.getNsigmaTOF(part, o2::track::PID::Proton), trackCuts.getNsigmaTPC(part, o2::track::PID::Pion), trackCuts.getNsigmaTOF(part, o2::track::PID::Pion), trackCuts.getNsigmaTPC(part, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(part, o2::track::PID::Kaon), trackCuts.getNsigmaTPC(part, o2::track::PID::Electron))) {
          continue;
        }
        trackHistoPartTwo.fillQA<isMC, true>(part);
      }
    }

    if (ContType == 1) {

      /// Now build the combinations for non-identical particle pairs
      for (const auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(groupPartsOne, groupPartsTwo))) {

        if (!IsParticleNSigma((int8_t)1, p1.p(), trackCuts.getNsigmaTPC(p1, o2::track::PID::Proton), trackCuts.getNsigmaTOF(p1, o2::track::PID::Proton), trackCuts.getNsigmaTPC(p1, o2::track::PID::Pion), trackCuts.getNsigmaTOF(p1, o2::track::PID::Pion), trackCuts.getNsigmaTPC(p1, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(p1, o2::track::PID::Kaon), trackCuts.getNsigmaTPC(p1, o2::track::PID::Electron))) {
          continue;
        }

        if (!IsParticleNSigma((int8_t)2, p2.p(), trackCuts.getNsigmaTPC(p2, o2::track::PID::Proton), trackCuts.getNsigmaTOF(p2, o2::track::PID::Proton), trackCuts.getNsigmaTPC(p2, o2::track::PID::Pion), trackCuts.getNsigmaTOF(p2, o2::track::PID::Pion), trackCuts.getNsigmaTPC(p2, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(p2, o2::track::PID::Kaon), trackCuts.getNsigmaTPC(p2, o2::track::PID::Electron))) {
          continue;
        }

        if (twotracksconfigs.ConfIsCPR.value) {
          if (confCPRIsAtITS.value) {
            if (pairCloseRejection.isClosePairAtITS(p1, p2, magFieldTesla, femto_universe_container::EventType::same)) {
              continue;
            }
          } else {
            if (pairCloseRejection.isClosePairFrac(p1, p2, magFieldTesla, femto_universe_container::EventType::same, confCPRDphiAvgOrDist, confCPRDistMax, confCPRFracMax, confIsCircularCut)) {
              continue;
            }
          }
        }

        // track cleaning
        if (!pairCleaner.isCleanPair(p1, p2, parts)) {
          continue;
        }
        float kT = FemtoUniverseMath::getkT(p1, mass1, p2, mass2);
        sameEventMultCont.fillMultNumDen(p1, p2, femto_universe_sh_container::EventType::same, 2, multCol, kT, ConfIsIden);
      }
    } else {
      for (const auto& [p1, p2] : combinations(CombinationsStrictlyUpperIndexPolicy(groupPartsOne, groupPartsOne))) {

        if (!IsParticleNSigma((int8_t)2, p1.p(), trackCuts.getNsigmaTPC(p1, o2::track::PID::Proton), trackCuts.getNsigmaTOF(p1, o2::track::PID::Proton), trackCuts.getNsigmaTPC(p1, o2::track::PID::Pion), trackCuts.getNsigmaTOF(p1, o2::track::PID::Pion), trackCuts.getNsigmaTPC(p1, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(p1, o2::track::PID::Kaon), trackCuts.getNsigmaTPC(p1, o2::track::PID::Electron))) {
          continue;
        }

        if (!IsParticleNSigma((int8_t)2, p2.p(), trackCuts.getNsigmaTPC(p2, o2::track::PID::Proton), trackCuts.getNsigmaTOF(p2, o2::track::PID::Proton), trackCuts.getNsigmaTPC(p2, o2::track::PID::Pion), trackCuts.getNsigmaTOF(p2, o2::track::PID::Pion), trackCuts.getNsigmaTPC(p2, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(p2, o2::track::PID::Kaon), trackCuts.getNsigmaTPC(p2, o2::track::PID::Electron))) {
          continue;
        }

        float kT = FemtoUniverseMath::getkT(p1, mass1, p2, mass2);
        if (twotracksconfigs.ConfIsCPR.value) {
          pairCloseRejection.kTdetadphi(p1, p2, femto_universe_container::EventType::same, kT);
          if (confCPRIsAtITS.value) {
            if (pairCloseRejection.isClosePairAtITS(p1, p2, magFieldTesla, femto_universe_container::EventType::same)) {
              continue;
            }
          } else {
            if (pairCloseRejection.isClosePairFrac(p1, p2, magFieldTesla, femto_universe_container::EventType::same, confCPRDphiAvgOrDist, confCPRDistMax, confCPRFracMax, confIsCircularCut)) {
              continue;
            }
          }
        }

        // track cleaning
        if (!pairCleaner.isCleanPair(p1, p2, parts)) {
          continue;
        }

        double rand;
        rand = randgen->Rndm();

        std::vector<double> f3d;
        double kv;
        float outsideref = 0.0;

        switch (ContType) {
          case 2: {
            if (rand > 0.5) {
              f3d = FemtoUniverseMath::newpairfunc(p1, mass1, p2, mass2, ConfIsIden);
              if (!twotracksconfigs.confUseCCImCut) {
                sameEventMultContPP.fillMultNumDen(p1, p2, femto_universe_sh_container::EventType::same, 2, multCol, kT, ConfIsIden);
              } else {
                if (twotracksconfigs.confUse1stand3rd) {
                  if ((f3d[1] >= outsideref && f3d[2] >= outsideref) || (f3d[1] < outsideref && f3d[2] < outsideref)) {
                    sameEventMultContPP.fillMultNumDen(p1, p2, femto_universe_sh_container::EventType::same, 2, multCol, kT, ConfIsIden);
                  }
                } else if (twotracksconfigs.confUse2ndand4th) {
                  if ((f3d[1] < outsideref && f3d[2] >= outsideref) || (f3d[1] >= outsideref && f3d[2] < outsideref)) {
                    sameEventMultContPP.fillMultNumDen(p1, p2, femto_universe_sh_container::EventType::same, 2, multCol, kT, ConfIsIden);
                  }
                }
              }
              if (twotracksconfigs.ConfIsMC) {
                float weight = 1.0f;
                sameEventCont1D.setPair<isMC>(p1, p2, multCol, twotracksconfigs.ConfUse3D, weight, ConfIsIden);
              }
            } else if (rand <= 0.5) {
              f3d = FemtoUniverseMath::newpairfunc(p2, mass2, p1, mass1, ConfIsIden);
              if (!twotracksconfigs.confUseCCImCut) {
                sameEventMultContPP.fillMultNumDen(p2, p1, femto_universe_sh_container::EventType::same, 2, multCol, kT, ConfIsIden);
              } else {
                if (twotracksconfigs.confUse1stand3rd) {
                  if ((f3d[1] >= outsideref && f3d[2] >= outsideref) || (f3d[1] < outsideref && f3d[2] < outsideref)) {
                    sameEventMultContPP.fillMultNumDen(p2, p1, femto_universe_sh_container::EventType::same, 2, multCol, kT, ConfIsIden);
                  }
                } else if (twotracksconfigs.confUse2ndand4th) {
                  if ((f3d[1] < outsideref && f3d[2] >= outsideref) || (f3d[1] >= outsideref && f3d[2] < outsideref)) {
                    sameEventMultContPP.fillMultNumDen(p2, p1, femto_universe_sh_container::EventType::same, 2, multCol, kT, ConfIsIden);
                  }
                }
              }
              if (twotracksconfigs.ConfIsMC) {
                float weight = 1.0f;
                sameEventCont1D.setPair<isMC>(p2, p1, multCol, twotracksconfigs.ConfUse3D, weight, ConfIsIden);
              }
            }
            if (ConfIsFillAngqLCMS) {
              kv = std::sqrt(f3d[1] * f3d[1] + f3d[2] * f3d[2] + f3d[3] * f3d[3]);
              pairCloseRejection.ClosePairqLCMS(p1, p2, magFieldTesla, femto_universe_container::EventType::same, kv);
            }
            break;
          }

          case 3: {
            if (rand > 0.5) {
              f3d = FemtoUniverseMath::newpairfunc(p1, mass1, p2, mass2, ConfIsIden);
              if (!twotracksconfigs.confUseCCImCut) {
                sameEventMultContMM.fillMultNumDen(p1, p2, femto_universe_sh_container::EventType::same, 2, multCol, kT, ConfIsIden);
              } else {
                if (twotracksconfigs.confUse1stand3rd) {
                  if ((f3d[1] >= outsideref && f3d[2] >= outsideref) || (f3d[1] < outsideref && f3d[2] < outsideref)) {
                    sameEventMultContMM.fillMultNumDen(p1, p2, femto_universe_sh_container::EventType::same, 2, multCol, kT, ConfIsIden);
                  }
                } else if (twotracksconfigs.confUse2ndand4th) {
                  if ((f3d[1] < outsideref && f3d[2] >= outsideref) || (f3d[1] >= outsideref && f3d[2] < outsideref)) {
                    sameEventMultContMM.fillMultNumDen(p1, p2, femto_universe_sh_container::EventType::same, 2, multCol, kT, ConfIsIden);
                  }
                }
              }
              if (twotracksconfigs.ConfIsMC) {
                float weight = 1.0f;
                sameEventCont1D.setPair<isMC>(p1, p2, multCol, twotracksconfigs.ConfUse3D, weight, ConfIsIden);
              }
            } else if (rand <= 0.5) {
              f3d = FemtoUniverseMath::newpairfunc(p2, mass2, p1, mass1, ConfIsIden);
              if (!twotracksconfigs.confUseCCImCut) {
                sameEventMultContMM.fillMultNumDen(p2, p1, femto_universe_sh_container::EventType::same, 2, multCol, kT, ConfIsIden);
              } else {
                if (twotracksconfigs.confUse1stand3rd) {
                  if ((f3d[1] >= outsideref && f3d[2] >= outsideref) || (f3d[1] < outsideref && f3d[2] < outsideref)) {
                    sameEventMultContMM.fillMultNumDen(p2, p1, femto_universe_sh_container::EventType::same, 2, multCol, kT, ConfIsIden);
                  }
                } else if (twotracksconfigs.confUse2ndand4th) {
                  if ((f3d[1] < outsideref && f3d[2] >= outsideref) || (f3d[1] >= outsideref && f3d[2] < outsideref)) {
                    sameEventMultContMM.fillMultNumDen(p2, p1, femto_universe_sh_container::EventType::same, 2, multCol, kT, ConfIsIden);
                  }
                }
              }
              if (twotracksconfigs.ConfIsMC) {
                float weight = 1.0f;
                sameEventCont1D.setPair<isMC>(p2, p1, multCol, twotracksconfigs.ConfUse3D, weight, ConfIsIden);
              }
            }
            if (ConfIsFillAngqLCMS) {
              kv = std::sqrt(f3d[1] * f3d[1] + f3d[2] * f3d[2] + f3d[3] * f3d[3]);
              pairCloseRejection.ClosePairqLCMS(p1, p2, magFieldTesla, femto_universe_container::EventType::same, kv);
            }
            break;
          }
          default:
            break;
        }
      }
    }
  }

  /// process function for to call doSameEvent with Data
  /// \param col subscribe to the collision table (Data)
  /// \param parts subscribe to the femtoUniverseParticleTable
  void processSameEvent(FilteredFDCollision const& col,
                        FilteredFemtoFullParticles const& parts)
  {
    fillCollision(col, ConfIsCent);

    auto thegroupPartsOne = partsOne->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto thegroupPartsTwo = partsTwo->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);

    bool fillQA = true;
    randgen = new TRandom2(0);

    if (ConfIsCent) {
      if (cfgProcessPM) {
        doSameEvent<false>(thegroupPartsOne, thegroupPartsTwo, parts, col.magField(), col.multV0M(), 1, fillQA);
      }
      if (cfgProcessPP) {
        doSameEvent<false>(thegroupPartsOne, thegroupPartsOne, parts, col.magField(), col.multV0M(), 2, fillQA);
      }
      if (cfgProcessMM) {
        doSameEvent<false>(thegroupPartsTwo, thegroupPartsTwo, parts, col.magField(), col.multV0M(), 3, fillQA);
      }
    } else {
      if (cfgProcessPM) {
        doSameEvent<false>(thegroupPartsOne, thegroupPartsTwo, parts, col.magField(), col.multNtr(), 1, fillQA);
      }
      if (cfgProcessPP) {
        doSameEvent<false>(thegroupPartsOne, thegroupPartsOne, parts, col.magField(), col.multNtr(), 2, fillQA);
      }
      if (cfgProcessMM) {
        doSameEvent<false>(thegroupPartsTwo, thegroupPartsTwo, parts, col.magField(), col.multNtr(), 3, fillQA);
      }
    }
    delete randgen;
  }
  PROCESS_SWITCH(femtoUniversePairTaskTrackTrackSpherHarMultKtExtended, processSameEvent, "Enable processing same event", true);

  /// process function for to call doSameEvent with Monte Carlo
  /// \param col subscribe to the collision table (Monte Carlo Reconstructed reconstructed)
  /// \param parts subscribe to joined table FemtoUniverseParticles and FemtoUniverseMCLables to access Monte Carlo truth
  /// \param FemtoUniverseMCParticles subscribe to the Monte Carlo truth table
  void processSameEventMC(o2::aod::FdCollision const& col,
                          FilteredFemtoRecoParticles const& parts,
                          aod::FdMCParticles const&)
  {
    fillCollision(col, ConfIsCent);

    auto thegroupPartsOne = partsOneMC->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto thegroupPartsTwo = partsTwoMC->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    bool fillQA = true;
    randgen = new TRandom2(0);

    if (ConfIsCent) {
      if (cfgProcessPM) {
        doSameEvent<false>(thegroupPartsOne, thegroupPartsTwo, parts, col.magField(), col.multV0M(), 1, fillQA);
      }
      if (cfgProcessPP) {
        doSameEvent<true>(thegroupPartsOne, thegroupPartsOne, parts, col.magField(), col.multV0M(), 2, fillQA);
      }
      if (cfgProcessMM) {
        doSameEvent<true>(thegroupPartsTwo, thegroupPartsTwo, parts, col.magField(), col.multV0M(), 3, fillQA);
      }
    } else {
      if (cfgProcessPM) {
        doSameEvent<false>(thegroupPartsOne, thegroupPartsTwo, parts, col.magField(), col.multNtr(), 1, fillQA);
      }
      if (cfgProcessPP) {
        doSameEvent<true>(thegroupPartsOne, thegroupPartsOne, parts, col.magField(), col.multNtr(), 2, fillQA);
      }
      if (cfgProcessMM) {
        doSameEvent<true>(thegroupPartsTwo, thegroupPartsTwo, parts, col.magField(), col.multNtr(), 3, fillQA);
      }
    }
    delete randgen;
  }
  PROCESS_SWITCH(femtoUniversePairTaskTrackTrackSpherHarMultKtExtended, processSameEventMC, "Enable processing same event for Monte Carlo", false);

  /// This function processes the same event and takes care of all the histogramming
  /// \todo the trivial loops over the tracks should be factored out since they will be common to all combinations of T-T, T-V0, V0-V0, ...
  /// @tparam PartitionType
  /// @tparam PartType
  /// @tparam isMC: enables Monte Carlo truth specific histograms
  /// @param groupPartsOne partition for the first particle passed by the process function
  /// @param groupPartsTwo partition for the second particle passed by the process function
  /// @param parts femtoUniverseParticles table (in case of Monte Carlo joined with FemtoUniverseMCLabels)
  /// @param magFieldTesla magnetic field of the collision
  /// @param multCol multiplicity of the collision
  template <bool isMC, typename PartitionType>
  void doSameEventMCTruth(PartitionType groupPartsOne, PartitionType groupPartsTwo, int multCol, int ContType, bool fillQA)
  {

    randgen = new TRandom2(0);
    /// Histogramming same event
    if ((cfgProcessPM || cfgProcessPP) && fillQA) {
      for (const auto& part : groupPartsOne) {
        if (part.partType() == uint8_t(aod::femtouniverseparticle::ParticleType::kMCTruthTrack)) {
          int pdgCode = static_cast<int>(part.tempFitVar());
          const auto& pdgParticle = pdg->GetParticle(pdgCode);
          if (pdgParticle) {
            trackHistoPartOne.fillQA<isMC, false>(part);
          }
        }
      }
    }

    if ((cfgProcessPM || cfgProcessMM) && fillQA) {
      for (const auto& part : groupPartsTwo) {
        if (part.partType() == uint8_t(aod::femtouniverseparticle::ParticleType::kMCTruthTrack)) {
          int pdgCode = static_cast<int>(part.tempFitVar());
          const auto& pdgParticle = pdg->GetParticle(pdgCode);
          if (pdgParticle) {
            trackHistoPartTwo.fillQA<isMC, false>(part);
          }
        }
      }
    }

    if (cfgProcessPM) {

      /// Now build the combinations for non-identical particle pairs
      for (const auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(groupPartsOne, groupPartsTwo))) {

        int pdgCodePartOne = static_cast<int>(p1.tempFitVar());
        const auto& pdgParticleOne = pdg->GetParticle(pdgCodePartOne);
        int pdgCodePartTwo = static_cast<int>(p2.tempFitVar());
        const auto& pdgParticleTwo = pdg->GetParticle(pdgCodePartTwo);
        if (pdgParticleOne && pdgParticleTwo && (pdgCodePartOne == trackonefilter.ConfPDGCodePartOne) && (pdgCodePartTwo == tracktwofilter.ConfPDGCodePartTwo)) {
          float kT = FemtoUniverseMath::getkT(p1, mass1, p2, mass2);
          sameEventMultCont.fillMultNumDen(p1, p2, femto_universe_sh_container::EventType::same, 2, multCol, kT, ConfIsIden);
        }
      }
    } else {
      for (const auto& [p1, p2] : combinations(CombinationsStrictlyUpperIndexPolicy(groupPartsOne, groupPartsOne))) {

        int pdgCodePartOne = static_cast<int>(p1.tempFitVar());
        const auto& pdgParticleOne = pdg->GetParticle(pdgCodePartOne);
        int pdgCodePartTwo = static_cast<int>(p2.tempFitVar());
        const auto& pdgParticleTwo = pdg->GetParticle(pdgCodePartTwo);

        float kT = FemtoUniverseMath::getkT(p1, mass1, p2, mass2);
        double rand;
        rand = randgen->Rndm();
        std::vector<double> f3d;
        if (pdgParticleOne && pdgParticleTwo && (pdgCodePartOne == trackonefilter.ConfPDGCodePartOne) && (pdgCodePartTwo == tracktwofilter.ConfPDGCodePartTwo)) {

          switch (ContType) {
            case 2: {
              if (rand > 0.5) {
                sameEventMultContPP.fillMultNumDen(p1, p2, femto_universe_sh_container::EventType::same, 2, multCol, kT, ConfIsIden);
              } else if (rand <= 0.5) {
                sameEventMultContPP.fillMultNumDen(p2, p1, femto_universe_sh_container::EventType::same, 2, multCol, kT, ConfIsIden);
              }
              break;
            }

            case 3: {
              if (rand > 0.5) {
                sameEventMultContMM.fillMultNumDen(p1, p2, femto_universe_sh_container::EventType::same, 2, multCol, kT, ConfIsIden);
              } else if (rand <= 0.5) {
                sameEventMultContMM.fillMultNumDen(p2, p1, femto_universe_sh_container::EventType::same, 2, multCol, kT, ConfIsIden);
              }
              break;
            }
            default:
              break;
          }
        }
      }
    }
    delete randgen;
  }

  /// process function for to call doSameEvent with Monte Carlo
  /// \param col subscribe to the collision table (Monte Carlo Reconstructed reconstructed)
  /// \param parts subscribe to joined table FemtoUniverseParticles and FemtoUniverseMCLables to access Monte Carlo truth
  /// \param FemtoUniverseMCParticles subscribe to the Monte Carlo truth table
  void processSameEventMCTruth(o2::aod::FdCollision const& col,
                               FemtoTruthParticles const&)
  {
    fillCollision(col, ConfIsCent);

    auto thegroupPartsOne = partsOneMCTruth->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto thegroupPartsTwo = partsTwoMCTruth->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    bool fillQA = true;

    int pairType = 0;
    if (cfgProcessPM) {
      pairType = 1;
    } else if (cfgProcessPP) {
      pairType = 2;
    } else if (cfgProcessMM) {
      pairType = 3;
    }

    if (ConfIsCent) {
      if (cfgProcessPM) {
        doSameEventMCTruth<false>(thegroupPartsOne, thegroupPartsTwo, col.multV0M(), pairType, fillQA);
      }
      if (cfgProcessPP) {
        doSameEventMCTruth<false>(thegroupPartsOne, thegroupPartsOne, col.multV0M(), pairType, fillQA);
      }
      if (cfgProcessMM) {
        doSameEventMCTruth<false>(thegroupPartsTwo, thegroupPartsTwo, col.multV0M(), pairType, fillQA);
      }
    } else {
      if (cfgProcessPM) {
        doSameEventMCTruth<false>(thegroupPartsOne, thegroupPartsTwo, col.multNtr(), pairType, fillQA);
      }
      if (cfgProcessPP) {
        doSameEventMCTruth<false>(thegroupPartsOne, thegroupPartsOne, col.multNtr(), pairType, fillQA);
      }
      if (cfgProcessMM) {
        doSameEventMCTruth<false>(thegroupPartsTwo, thegroupPartsTwo, col.multNtr(), pairType, fillQA);
      }
    }
  }
  PROCESS_SWITCH(femtoUniversePairTaskTrackTrackSpherHarMultKtExtended, processSameEventMCTruth, "Enable processing same event for MC truth", false);

  /// This function processes the mixed event
  /// \todo the trivial loops over the collisions and tracks should be factored out since they will be common to all combinations of T-T, T-V0, V0-V0, ...
  /// \tparam PartitionType
  /// \tparam PartType
  /// \tparam isMC: enables Monte Carlo truth specific histograms
  /// \param groupPartsOne partition for the first particle passed by the process function
  /// \param groupPartsTwo partition for the second particle passed by the process function
  /// \param parts femtoUniverseParticles table (in case of Monte Carlo joined with FemtoUniverseMCLabels)
  /// \param magFieldTesla magnetic field of the collision
  /// \param multCol multiplicity of the collision
  template <bool isMC, typename PartitionType>
  void doMixedEvent(PartitionType groupPartsOne, PartitionType groupPartsTwo, float magFieldTesla, int multCol, int ContType)
  {

    for (const auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(groupPartsOne, groupPartsTwo))) {

      if (!IsParticleNSigma((int8_t)2, p1.p(), trackCuts.getNsigmaTPC(p1, o2::track::PID::Proton), trackCuts.getNsigmaTOF(p1, o2::track::PID::Proton), trackCuts.getNsigmaTPC(p1, o2::track::PID::Pion), trackCuts.getNsigmaTOF(p1, o2::track::PID::Pion), trackCuts.getNsigmaTPC(p1, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(p1, o2::track::PID::Kaon), trackCuts.getNsigmaTPC(p1, o2::track::PID::Electron))) {
        continue;
      }

      if (!IsParticleNSigma((int8_t)2, p2.p(), trackCuts.getNsigmaTPC(p2, o2::track::PID::Proton), trackCuts.getNsigmaTOF(p2, o2::track::PID::Proton), trackCuts.getNsigmaTPC(p2, o2::track::PID::Pion), trackCuts.getNsigmaTOF(p2, o2::track::PID::Pion), trackCuts.getNsigmaTPC(p2, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(p2, o2::track::PID::Kaon), trackCuts.getNsigmaTPC(p2, o2::track::PID::Electron))) {
        continue;
      }

      float kT = FemtoUniverseMath::getkT(p1, mass1, p2, mass2);
      if (twotracksconfigs.ConfIsCPR.value) {
        pairCloseRejection.kTdetadphi(p1, p2, femto_universe_container::EventType::mixed, kT);
        if (confCPRIsAtITS.value) {
          if (pairCloseRejection.isClosePairAtITS(p1, p2, magFieldTesla, femto_universe_container::EventType::mixed)) {
            continue;
          }
        } else {
          if (pairCloseRejection.isClosePairFrac(p1, p2, magFieldTesla, femto_universe_container::EventType::mixed, confCPRDphiAvgOrDist, confCPRDistMax, confCPRFracMax, confIsCircularCut)) {
            continue;
          }
        }
      }

      double rand;
      rand = randgen->Rndm();

      std::vector<double> f3d;
      double kv;
      float outsideref = 0.0;

      switch (ContType) {
        case 1: {
          if (rand > 0.5) {
            mixedEventMultCont.fillMultNumDen(p1, p2, femto_universe_sh_container::EventType::mixed, 2, multCol, kT, ConfIsIden);
          } else {
            mixedEventMultCont.fillMultNumDen(p2, p1, femto_universe_sh_container::EventType::mixed, 2, multCol, kT, ConfIsIden);
          }
          break;
        }

        case 2: {
          if (rand > 0.5) {
            f3d = FemtoUniverseMath::newpairfunc(p1, mass1, p2, mass2, ConfIsIden);
            if (!twotracksconfigs.confUseCCImCut) {
              mixedEventMultContPP.fillMultNumDen(p1, p2, femto_universe_sh_container::EventType::mixed, 2, multCol, kT, ConfIsIden);
            } else {
              if (twotracksconfigs.confUse1stand3rd) {
                if ((f3d[1] >= outsideref && f3d[2] >= outsideref) || (f3d[1] < outsideref && f3d[2] < outsideref)) {
                  mixedEventMultContPP.fillMultNumDen(p1, p2, femto_universe_sh_container::EventType::mixed, 2, multCol, kT, ConfIsIden);
                }
              } else if (twotracksconfigs.confUse2ndand4th) {
                if ((f3d[1] < outsideref && f3d[2] >= outsideref) || (f3d[1] >= outsideref && f3d[2] < outsideref)) {
                  mixedEventMultContPP.fillMultNumDen(p1, p2, femto_universe_sh_container::EventType::mixed, 2, multCol, kT, ConfIsIden);
                }
              }
            }
          } else {
            f3d = FemtoUniverseMath::newpairfunc(p2, mass2, p1, mass1, ConfIsIden);
            if (!twotracksconfigs.confUseCCImCut) {
              mixedEventMultContPP.fillMultNumDen(p2, p1, femto_universe_sh_container::EventType::mixed, 2, multCol, kT, ConfIsIden);
            } else {
              if (twotracksconfigs.confUse1stand3rd) {
                if ((f3d[1] >= outsideref && f3d[2] >= outsideref) || (f3d[1] < outsideref && f3d[2] < outsideref)) {
                  mixedEventMultContPP.fillMultNumDen(p2, p1, femto_universe_sh_container::EventType::mixed, 2, multCol, kT, ConfIsIden);
                }
              } else if (twotracksconfigs.confUse2ndand4th) {
                if ((f3d[1] < outsideref && f3d[2] >= outsideref) || (f3d[1] >= outsideref && f3d[2] < outsideref)) {
                  mixedEventMultContPP.fillMultNumDen(p2, p1, femto_universe_sh_container::EventType::mixed, 2, multCol, kT, ConfIsIden);
                }
              }
            }
          }
          if (ConfIsFillAngqLCMS) {
            kv = std::sqrt(f3d[1] * f3d[1] + f3d[2] * f3d[2] + f3d[3] * f3d[3]);
            pairCloseRejection.ClosePairqLCMS(p1, p2, magFieldTesla, femto_universe_container::EventType::mixed, kv);
          }
          break;
        }

        case 3: {
          if (rand > 0.5) {
            f3d = FemtoUniverseMath::newpairfunc(p1, mass1, p2, mass2, ConfIsIden);
            if (!twotracksconfigs.confUseCCImCut) {
              mixedEventMultContMM.fillMultNumDen(p1, p2, femto_universe_sh_container::EventType::mixed, 2, multCol, kT, ConfIsIden);
            } else {
              if (twotracksconfigs.confUse1stand3rd) {
                if ((f3d[1] >= outsideref && f3d[2] >= outsideref) || (f3d[1] < outsideref && f3d[2] < outsideref)) {
                  mixedEventMultContMM.fillMultNumDen(p1, p2, femto_universe_sh_container::EventType::mixed, 2, multCol, kT, ConfIsIden);
                }
              } else if (twotracksconfigs.confUse2ndand4th) {
                if ((f3d[1] < outsideref && f3d[2] >= outsideref) || (f3d[1] >= outsideref && f3d[2] < outsideref)) {
                  mixedEventMultContMM.fillMultNumDen(p1, p2, femto_universe_sh_container::EventType::mixed, 2, multCol, kT, ConfIsIden);
                }
              }
            }
          } else {
            f3d = FemtoUniverseMath::newpairfunc(p2, mass2, p1, mass1, ConfIsIden);
            if (!twotracksconfigs.confUseCCImCut) {
              mixedEventMultContMM.fillMultNumDen(p2, p1, femto_universe_sh_container::EventType::mixed, 2, multCol, kT, ConfIsIden);
            } else {
              if (twotracksconfigs.confUse1stand3rd) {
                if ((f3d[1] >= outsideref && f3d[2] >= outsideref) || (f3d[1] < outsideref && f3d[2] < outsideref)) {
                  mixedEventMultContMM.fillMultNumDen(p2, p1, femto_universe_sh_container::EventType::mixed, 2, multCol, kT, ConfIsIden);
                }
              } else if (twotracksconfigs.confUse2ndand4th) {
                if ((f3d[1] < outsideref && f3d[2] >= outsideref) || (f3d[1] >= outsideref && f3d[2] < outsideref)) {
                  mixedEventMultContMM.fillMultNumDen(p2, p1, femto_universe_sh_container::EventType::mixed, 2, multCol, kT, ConfIsIden);
                }
              }
            }
          }
          if (ConfIsFillAngqLCMS) {
            kv = std::sqrt(f3d[1] * f3d[1] + f3d[2] * f3d[2] + f3d[3] * f3d[3]);
            pairCloseRejection.ClosePairqLCMS(p1, p2, magFieldTesla, femto_universe_container::EventType::mixed, kv);
          }
          break;
        }

        default:
          break;
      }
    }
  }

  /// process function for to call doMixedEvent with Data
  /// @param cols subscribe to the collisions table (Data)
  /// @param  subscribe to the femtoUniverseParticleTable
  void processMixedEventCent(FilteredFDCollisions const& cols,
                             FilteredFemtoFullParticles const&)
  {
    randgen = new TRandom2(0);

    for (const auto& [collision1, collision2] : soa::selfCombinations(colBinningCent, ConfNEventsMix, -1, cols, cols)) {

      const int multiplicityCol = collision1.multV0M();
      MixQaRegistry.fill(HIST("MixingQA/hMECollisionBins"), colBinningCent.getBin({collision1.posZ(), multiplicityCol}));

      const auto& magFieldTesla1 = collision1.magField();
      const auto& magFieldTesla2 = collision2.magField();

      if (magFieldTesla1 != magFieldTesla2) {
        continue;
      }

      if (cfgProcessPM) {
        auto groupPartsOne = partsOne->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);
        auto groupPartsTwo = partsTwo->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision2.globalIndex(), cache);
        doMixedEvent<false>(groupPartsOne, groupPartsTwo, magFieldTesla1, multiplicityCol, 1);
      }
      if (cfgProcessPP) {
        auto groupPartsOne = partsOne->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);
        auto groupPartsTwo = partsOne->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision2.globalIndex(), cache);
        doMixedEvent<false>(groupPartsOne, groupPartsTwo, magFieldTesla1, multiplicityCol, 2);
      }
      if (cfgProcessMM) {
        auto groupPartsOne = partsTwo->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);
        auto groupPartsTwo = partsTwo->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision2.globalIndex(), cache);
        doMixedEvent<false>(groupPartsOne, groupPartsTwo, magFieldTesla1, multiplicityCol, 3);
      }
    }
    delete randgen;
  }
  PROCESS_SWITCH(femtoUniversePairTaskTrackTrackSpherHarMultKtExtended, processMixedEventCent, "Enable processing mixed events for centrality", true);

  /// process function for to call doMixedEvent with Data
  /// @param cols subscribe to the collisions table (Data)
  /// @param parts subscribe to the femtoUniverseParticleTable
  void processMixedEventNtr(FilteredFDCollisions const& cols,
                            FilteredFemtoFullParticles const&)
  {
    randgen = new TRandom2(0);

    for (const auto& [collision1, collision2] : soa::selfCombinations(colBinningNtr, ConfNEventsMix, -1, cols, cols)) {

      const int multiplicityCol = collision1.multNtr();
      MixQaRegistry.fill(HIST("MixingQA/hMECollisionBins"), colBinningNtr.getBin({collision1.posZ(), multiplicityCol}));

      const auto& magFieldTesla1 = collision1.magField();
      const auto& magFieldTesla2 = collision2.magField();

      if (magFieldTesla1 != magFieldTesla2) {
        continue;
      }

      if (cfgProcessPM) {
        auto groupPartsOne = partsOne->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);
        auto groupPartsTwo = partsTwo->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision2.globalIndex(), cache);
        doMixedEvent<false>(groupPartsOne, groupPartsTwo, magFieldTesla1, multiplicityCol, 1);
      }
      if (cfgProcessPP) {
        auto groupPartsOne = partsOne->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);
        auto groupPartsTwo = partsOne->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision2.globalIndex(), cache);
        doMixedEvent<false>(groupPartsOne, groupPartsTwo, magFieldTesla1, multiplicityCol, 2);
      }
      if (cfgProcessMM) {
        auto groupPartsOne = partsTwo->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);
        auto groupPartsTwo = partsTwo->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision2.globalIndex(), cache);
        doMixedEvent<false>(groupPartsOne, groupPartsTwo, magFieldTesla1, multiplicityCol, 3);
      }
    }
    delete randgen;
  }
  PROCESS_SWITCH(femtoUniversePairTaskTrackTrackSpherHarMultKtExtended, processMixedEventNtr, "Enable processing mixed events for centrality", false);

  /// brief process function for to call doMixedEvent with Monte Carlo
  /// @param cols subscribe to the collisions table (Monte Carlo Reconstructed reconstructed)
  /// @param parts subscribe to joined table FemtoUniverseParticles and FemtoUniverseMCLables to access Monte Carlo truth
  /// @param FemtoUniverseMCParticles subscribe to the Monte Carlo truth table
  void processMixedEventMCCent(o2::aod::FdCollisions const& cols,
                               soa::Join<FilteredFemtoFullParticles, aod::FDMCLabels> const&,
                               o2::aod::FdMCParticles const&)
  {
    randgen = new TRandom2(0);

    for (const auto& [collision1, collision2] : soa::selfCombinations(colBinningCent, ConfNEventsMix, -1, cols, cols)) {

      const int multiplicityCol = collision1.multV0M();
      MixQaRegistry.fill(HIST("MixingQA/hMECollisionBins"), colBinningCent.getBin({collision1.posZ(), multiplicityCol}));

      const auto& magFieldTesla1 = collision1.magField();
      const auto& magFieldTesla2 = collision2.magField();

      if (magFieldTesla1 != magFieldTesla2) {
        continue;
      }
      /// \todo before mixing we should check whether both collisions contain a pair of particles!
      // if (partsOne.size() == 0 || nPart2Evt1 == 0 || nPart1Evt2 == 0 || partsTwo.size() == 0 ) continue;

      if (cfgProcessPM) {
        auto groupPartsOne = partsOneMC->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);
        auto groupPartsTwo = partsTwoMC->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision2.globalIndex(), cache);
        doMixedEvent<false>(groupPartsOne, groupPartsTwo, magFieldTesla1, multiplicityCol, 1);
      }
      if (cfgProcessPP) {
        auto groupPartsOne = partsOneMC->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);
        auto groupPartsTwo = partsOneMC->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision2.globalIndex(), cache);
        doMixedEvent<false>(groupPartsOne, groupPartsTwo, magFieldTesla1, multiplicityCol, 2);
      }
      if (cfgProcessMM) {
        auto groupPartsOne = partsTwoMC->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);
        auto groupPartsTwo = partsTwoMC->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision2.globalIndex(), cache);
        doMixedEvent<false>(groupPartsOne, groupPartsTwo, magFieldTesla1, multiplicityCol, 3);
      }
    }
    delete randgen;
  }
  PROCESS_SWITCH(femtoUniversePairTaskTrackTrackSpherHarMultKtExtended, processMixedEventMCCent, "Enable processing mixed events MC", false);

  /// brief process function for to call doMixedEvent with Monte Carlo
  /// @param cols subscribe to the collisions table (Monte Carlo Reconstructed reconstructed)
  /// @param parts subscribe to joined table FemtoUniverseParticles and FemtoUniverseMCLables to access Monte Carlo truth
  /// @param FemtoUniverseMCParticles subscribe to the Monte Carlo truth table
  void processMixedEventMCNtr(o2::aod::FdCollisions const& cols,
                              soa::Join<FilteredFemtoFullParticles, aod::FDMCLabels> const&,
                              o2::aod::FdMCParticles const&)
  {
    randgen = new TRandom2(0);

    for (const auto& [collision1, collision2] : soa::selfCombinations(colBinningNtr, ConfNEventsMix, -1, cols, cols)) {

      const int multiplicityCol = collision1.multV0M();
      MixQaRegistry.fill(HIST("MixingQA/hMECollisionBins"), colBinningNtr.getBin({collision1.posZ(), multiplicityCol}));

      const auto& magFieldTesla1 = collision1.magField();
      const auto& magFieldTesla2 = collision2.magField();

      if (magFieldTesla1 != magFieldTesla2) {
        continue;
      }
      /// \todo before mixing we should check whether both collisions contain a pair of particles!
      // if (partsOne.size() == 0 || nPart2Evt1 == 0 || nPart1Evt2 == 0 || partsTwo.size() == 0 ) continue;

      if (cfgProcessPM) {
        auto groupPartsOne = partsOneMC->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);
        auto groupPartsTwo = partsTwoMC->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision2.globalIndex(), cache);
        doMixedEvent<false>(groupPartsOne, groupPartsTwo, magFieldTesla1, multiplicityCol, 1);
      }
      if (cfgProcessPP) {
        auto groupPartsOne = partsOneMC->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);
        auto groupPartsTwo = partsOneMC->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision2.globalIndex(), cache);
        doMixedEvent<false>(groupPartsOne, groupPartsTwo, magFieldTesla1, multiplicityCol, 2);
      }
      if (cfgProcessMM) {
        auto groupPartsOne = partsTwoMC->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);
        auto groupPartsTwo = partsTwoMC->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision2.globalIndex(), cache);
        doMixedEvent<false>(groupPartsOne, groupPartsTwo, magFieldTesla1, multiplicityCol, 3);
      }
    }
    delete randgen;
  }
  PROCESS_SWITCH(femtoUniversePairTaskTrackTrackSpherHarMultKtExtended, processMixedEventMCNtr, "Enable processing mixed events MC", false);

  /// This function processes the mixed event
  /// \todo the trivial loops over the collisions and tracks should be factored out since they will be common to all combinations of T-T, T-V0, V0-V0, ...
  /// \tparam PartitionType
  /// \tparam PartType
  /// \tparam isMC: enables Monte Carlo truth specific histograms
  /// \param groupPartsOne partition for the first particle passed by the process function
  /// \param groupPartsTwo partition for the second particle passed by the process function
  /// \param parts femtoUniverseParticles table (in case of Monte Carlo joined with FemtoUniverseMCLabels)
  /// \param magFieldTesla magnetic field of the collision
  /// \param multCol multiplicity of the collision
  template <bool isMC, typename PartitionType>
  void doMixedEventMCTruth(PartitionType groupPartsOne, PartitionType groupPartsTwo, int multCol, int ContType)
  {
    randgen = new TRandom2(0);
    for (const auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(groupPartsOne, groupPartsTwo))) {

      int pdgCodePartOne = static_cast<int>(p1.tempFitVar());
      const auto& pdgParticleOne = pdg->GetParticle(pdgCodePartOne);
      int pdgCodePartTwo = static_cast<int>(p2.tempFitVar());
      const auto& pdgParticleTwo = pdg->GetParticle(pdgCodePartTwo);

      float kT = FemtoUniverseMath::getkT(p1, mass1, p2, mass2);
      double rand;
      rand = randgen->Rndm();
      if (pdgParticleOne && pdgParticleTwo && (pdgCodePartOne == trackonefilter.ConfPDGCodePartOne) && (pdgCodePartTwo == tracktwofilter.ConfPDGCodePartTwo)) {

        switch (ContType) {
          case 1: {
            if (rand > 0.5) {
              mixedEventMultCont.fillMultNumDen(p1, p2, femto_universe_sh_container::EventType::mixed, 2, multCol, kT, ConfIsIden);
            } else {
              mixedEventMultCont.fillMultNumDen(p2, p1, femto_universe_sh_container::EventType::mixed, 2, multCol, kT, ConfIsIden);
            }
            break;
          }

          case 2: {
            if (rand > 0.5) {
              mixedEventMultContPP.fillMultNumDen(p1, p2, femto_universe_sh_container::EventType::mixed, 2, multCol, kT, ConfIsIden);
            } else {
              mixedEventMultContPP.fillMultNumDen(p2, p1, femto_universe_sh_container::EventType::mixed, 2, multCol, kT, ConfIsIden);
            }
            break;
          }

          case 3: {
            if (rand > 0.5) {
              mixedEventMultContMM.fillMultNumDen(p1, p2, femto_universe_sh_container::EventType::mixed, 2, multCol, kT, ConfIsIden);
            } else {
              mixedEventMultContMM.fillMultNumDen(p2, p1, femto_universe_sh_container::EventType::mixed, 2, multCol, kT, ConfIsIden);
            }
            break;
          }
          default:
            break;
        }
      }
    }
    delete randgen;
  }

  /// process function for to call doMixedEvent with Data
  /// @param cols subscribe to the collisions table (Data)
  /// @param parts subscribe to the femtoUniverseParticleTable
  void processMixedEventNtrMCTruth(o2::aod::FdCollisions const& cols,
                                   FemtoTruthParticles const&)
  {
    int pairType = 0;
    if (cfgProcessPM) {
      pairType = 1;
    } else if (cfgProcessPP) {
      pairType = 2;
    } else if (cfgProcessMM) {
      pairType = 3;
    }

    for (const auto& [collision1, collision2] : soa::selfCombinations(colBinningNtr, ConfNEventsMix, -1, cols, cols)) {

      const int multiplicityCol = collision1.multNtr();
      MixQaRegistry.fill(HIST("MixingQA/hMECollisionBins"), colBinningNtr.getBin({collision1.posZ(), multiplicityCol}));

      if (cfgProcessPM) {
        auto groupPartsOne = partsOneMCTruth->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);
        auto groupPartsTwo = partsTwoMCTruth->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision2.globalIndex(), cache);
        doMixedEventMCTruth<false>(groupPartsOne, groupPartsTwo, multiplicityCol, pairType);
      }
      if (cfgProcessPP) {
        auto groupPartsOne = partsOneMCTruth->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);
        auto groupPartsTwo = partsOneMCTruth->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision2.globalIndex(), cache);
        doMixedEventMCTruth<false>(groupPartsOne, groupPartsTwo, multiplicityCol, pairType);
      }
      if (cfgProcessMM) {
        auto groupPartsOne = partsTwoMCTruth->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);
        auto groupPartsTwo = partsTwoMCTruth->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision2.globalIndex(), cache);
        doMixedEventMCTruth<false>(groupPartsOne, groupPartsTwo, multiplicityCol, pairType);
      }
    }
  }
  PROCESS_SWITCH(femtoUniversePairTaskTrackTrackSpherHarMultKtExtended, processMixedEventNtrMCTruth, "Enable processing MC Truth mixed events for multiplicity", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<femtoUniversePairTaskTrackTrackSpherHarMultKtExtended>(cfgc),
  };
  return workflow;
}
