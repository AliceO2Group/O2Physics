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

/// \file femtoUniversePairTaskTrackTrackExtended.cxx
/// \brief Tasks that reads the track tables used for the pairing and builds pairs of two tracks
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de
/// \author Georgios Mantzaridis, TU München, georgios.mantzaridis@tum.de
/// \author Anton Riedel, TU München, anton.riedel@tum.de
/// \author Zuzanna Chochulska, WUT Warsaw & CTU Prague, zchochul@cern.ch

#include <Framework/O2DatabasePDGPlugin.h>
#include <string>
#include <vector>

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/PID.h"

#include "PWGCF/FemtoUniverse/DataModel/FemtoDerived.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseParticleHisto.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseEventHisto.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniversePairCleaner.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseContainer.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseDetaDphiStar.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseTrackSelection.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseEfficiencyCalculator.h"

using namespace o2;
using namespace o2::analysis::femto_universe;
using namespace o2::analysis::femto_universe::efficiency;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

namespace
{
static constexpr int NPart = 2;
static constexpr int NCuts = 5;
static const std::vector<std::string> partNames{"PartOne", "PartTwo"};
static const std::vector<std::string> cutNames{"MaxPt", "PIDthr", "nSigmaTPC", "nSigmaTPCTOF", "MaxP"};
static const float cutsTable[NPart][NCuts]{
  {4.05f, 1.f, 3.f, 3.f, 100.f},
  {4.05f, 1.f, 3.f, 3.f, 100.f}};
} // namespace

struct FemtoUniversePairTaskTrackTrackExtended {
  Service<o2::framework::O2DatabasePDG> pdgMC;

  /// Particle selection part

  /// Table for both particles
  struct : o2::framework::ConfigurableGroup {
    Configurable<float> confNsigmaCombined{"confNsigmaCombined", 3.0f, "TPC and TOF Pion Sigma (combined) for momentum > confTOFPtMin"};
    Configurable<float> confNsigmaTPC{"confNsigmaTPC", 3.0f, "TPC Pion Sigma for momentum < confTOFPtMin"};
    Configurable<float> confTOFPtMin{"confTOFPtMin", 0.5f, "Min. Pt for which TOF is required for PID."};
    Configurable<float> confEtaMax{"confEtaMax", 0.8f, "Higher limit for |Eta| (the same for both particles)"};

    Configurable<LabeledArray<float>> confCutTable{"confCutTable", {cutsTable[0], NPart, NCuts, partNames, cutNames}, "Particle selections"};
    Configurable<int> confNspecies{"confNspecies", 2, "Number of particle spieces with PID info"};
    Configurable<bool> confIsMC{"confIsMC", false, "Enable additional Histogramms in the case of a MonteCarlo Run"};
    Configurable<std::vector<float>> confTrkPIDnSigmaMax{"confTrkPIDnSigmaMax", std::vector<float>{4.f, 3.f, 2.f}, "This configurable needs to be the same as the one used in the producer task"};
    Configurable<bool> confUse3D{"confUse3D", false, "Enable three dimensional histogramms (to be used only for analysis with high statistics): k* vs mT vs multiplicity"};
  } twotracksconfigs;

  using FemtoFullParticles = soa::Join<aod::FDParticles, aod::FDExtParticles>;
  // Filters for selecting particles (both p1 and p2)
  Filter trackCutFilter = requireGlobalTrackInFilter();                                                 // Global track cuts
  Filter trackAdditionalfilter = (nabs(aod::femtouniverseparticle::eta) < twotracksconfigs.confEtaMax); // example filtering on configurable
  using FilteredFemtoFullParticles = soa::Filtered<FemtoFullParticles>;
  // using FilteredFemtoFullParticles = FemtoFullParticles; //if no filtering is applied uncomment this option

  SliceCache cache;
  Preslice<FilteredFemtoFullParticles> perCol = aod::femtouniverseparticle::fdCollisionId;

  /// Particle 1
  struct : o2::framework::ConfigurableGroup {
    Configurable<int> confPDGCodePartOne{"confPDGCodePartOne", 2212, "Particle 1 - PDG code"};
    Configurable<bool> confIsTrackOneIdentified{"confIsTrackOneIdentified", true, "Enable PID for the track one"};
    // Configurable<uint32_t> confCutPartOne{"confCutPartOne", 5542474, "Particle 1 - Selection bit from cutCulator"};
    Configurable<int> confPIDPartOne{"confPIDPartOne", 2, "Particle 1 - Read from cutCulator"}; // we also need the possibility to specify whether the bit is true/false ->std>>vector<std::pair<int, int>>int>>
    Configurable<float> confPtLowPart1{"confPtLowPart1", 0.5, "Lower limit for Pt for the first particle"};
    Configurable<float> confPtHighPart1{"confPtHighPart1", 1.5, "Higher limit for Pt for the first particle"};
    Configurable<int> confChargePart1{"confChargePart1", 1, "Particle 1 sign"};
  } trackonefilter;

  /// Partition for particle 1
  Partition<FilteredFemtoFullParticles> partsOne = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kTrack)) && aod::femtouniverseparticle::sign == trackonefilter.confChargePart1 && aod::femtouniverseparticle::pt < trackonefilter.confPtHighPart1 && aod::femtouniverseparticle::pt > trackonefilter.confPtLowPart1;
  Partition<soa::Join<FilteredFemtoFullParticles, aod::FDMCLabels>> partsOneMCReco = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kTrack)) && aod::femtouniverseparticle::sign == trackonefilter.confChargePart1 && aod::femtouniverseparticle::pt < trackonefilter.confPtHighPart1 && aod::femtouniverseparticle::pt > trackonefilter.confPtLowPart1;
  // && ((aod::femtouniverseparticle::cut & confCutPartOne) == confCutPartOne);

  Partition<soa::Join<FilteredFemtoFullParticles, aod::FDMCLabels>> partsOneMCTruth =
    aod::femtouniverseparticle::partType == static_cast<uint8_t>(aod::femtouniverseparticle::ParticleType::kMCTruthTrack) &&
    aod::femtouniverseparticle::pt < trackonefilter.confPtHighPart1 &&
    aod::femtouniverseparticle::pt > trackonefilter.confPtLowPart1;

  /// Histogramming for particle 1
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kTrack, 1> trackHistoPartOne;
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kMCTruthTrack, 1> hMCTruth1;

  /// Particle 2
  Configurable<bool> confIsSame{"confIsSame", false, "Pairs of the same particle"};
  struct : o2::framework::ConfigurableGroup {
    Configurable<int> confPDGCodePartTwo{"confPDGCodePartTwo", 2212, "Particle 2 - PDG code"};
    Configurable<bool> confIsTrackTwoIdentified{"confIsTrackTwoIdentified", true, "Enable PID for the track two"};
    Configurable<int> confPIDPartTwo{"confPIDPartTwo", 2, "Particle 2 - Read from cutCulator"}; // we also need the possibility to specify whether the bit is true/false ->std>>vector<std::pair<int, int>>
    Configurable<float> confPtLowPart2{"confPtLowPart2", 0.5, "Lower limit for Pt for the second particle"};
    Configurable<float> confPtHighPart2{"confPtHighPart2", 1.5, "Higher limit for Pt for the second particle"};
    Configurable<int> confChargePart2{"confChargePart2", -1, "Particle 2 sign"};
  } tracktwofilter;
  /// Partition for particle 2
  Partition<FilteredFemtoFullParticles> partsTwo = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kTrack)) && (aod::femtouniverseparticle::sign == tracktwofilter.confChargePart2) && aod::femtouniverseparticle::pt < tracktwofilter.confPtHighPart2 && aod::femtouniverseparticle::pt > tracktwofilter.confPtLowPart2;

  Partition<soa::Join<FilteredFemtoFullParticles, aod::FDMCLabels>> partsTwoMCReco = aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kTrack) && (aod::femtouniverseparticle::sign == tracktwofilter.confChargePart2) && aod::femtouniverseparticle::pt < tracktwofilter.confPtHighPart2 && aod::femtouniverseparticle::pt > tracktwofilter.confPtLowPart2;

  Partition<soa::Join<FilteredFemtoFullParticles, aod::FDMCLabels>> partsTwoMCTruth =
    aod::femtouniverseparticle::partType == static_cast<uint8_t>(aod::femtouniverseparticle::ParticleType::kMCTruthTrack) &&
    aod::femtouniverseparticle::pt < tracktwofilter.confPtHighPart2 &&
    aod::femtouniverseparticle::pt > tracktwofilter.confPtLowPart2;

  /// Histogramming for particle 2
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kTrack, 2> trackHistoPartTwo;
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kMCTruthTrack, 2> hMCTruth2;

  /// Histogramming for Event
  FemtoUniverseEventHisto eventHisto;

  /// The configurables need to be passed to an std::vector
  int vPIDPartOne, vPIDPartTwo;
  std::vector<float> kNsigma;

  /// particle part
  ConfigurableAxis confTempFitVarBins{"confTempFitVarBins", {300, -0.15, 0.15}, "binning of the TempFitVar in the pT vs. TempFitVar plot"};
  ConfigurableAxis confTempFitVarpTBins{"confTempFitVarpTBins", {20, 0.5, 4.05}, "pT binning of the pT vs. TempFitVar plot"};
  ConfigurableAxis confTempFitVarPDGBins{"confTempFitVarPDGBins", {6000, -2300, 2300}, "Binning of the PDG code in the pT vs. TempFitVar plot"};

  /// Correlation part
  ConfigurableAxis confMultBins{"confMultBins", {VARIABLE_WIDTH, 0.0f, 4.0f, 8.0f, 12.0f, 16.0f, 20.0f, 24.0f, 28.0f, 32.0f, 36.0f, 40.0f, 44.0f, 48.0f, 52.0f, 56.0f, 60.0f, 64.0f, 68.0f, 72.0f, 76.0f, 80.0f, 84.0f, 88.0f, 92.0f, 96.0f, 100.0f, 200.0f, 99999.f}, "Mixing bins - multiplicity"}; // \todo to be obtained from the hash task
  // ConfigurableAxis confMultBins{"confMultBins", {VARIABLE_WIDTH, 0.0f, 20.0f, 40.0f, 60.0f, 80.0f, 100.0f, 200.0f, 99999.f}, "Mixing bins - multiplicity"};
  ConfigurableAxis confVtxBins{"confVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};

  ConfigurableAxis confmTBins3D{"confmTBins3D", {VARIABLE_WIDTH, 1.02f, 1.14f, 1.20f, 1.26f, 1.38f, 1.56f, 1.86f, 4.50f}, "mT Binning for the 3Dimensional plot: k* vs multiplicity vs mT (set <<twotracksconfigs.confUse3D>> to true in order to use)"};
  ConfigurableAxis confMultBins3D{"confMultBins3D", {VARIABLE_WIDTH, 0.0f, 20.0f, 30.0f, 40.0f, 99999.0f}, "multiplicity Binning for the 3Dimensional plot: k* vs multiplicity vs mT (set <<twotracksconfigs.confUse3D>> to true in order to use)"};

  ColumnBinningPolicy<aod::collision::PosZ, aod::femtouniversecollision::MultNtr> colBinning{{confVtxBins, confMultBins}, true};

  ConfigurableAxis confkstarBins{"confkstarBins", {1500, 0., 6.}, "binning kstar"};
  ConfigurableAxis confkTBins{"confkTBins", {150, 0., 9.}, "binning kT"};
  ConfigurableAxis confmTBins{"confmTBins", {225, 0., 7.5}, "binning mT"};
  Configurable<int> confNEventsMix{"confNEventsMix", 5, "Number of events for mixing"};
  Configurable<bool> confIsCPR{"confIsCPR", true, "Close Pair Rejection"};
  Configurable<bool> confCPRPlotPerRadii{"confCPRPlotPerRadii", false, "Plot CPR per radii"};
  Configurable<float> confCPRdeltaPhiCutMax{"confCPRdeltaPhiCutMax", 0.0, "Delta Phi max cut for Close Pair Rejection"};
  Configurable<float> confCPRdeltaPhiCutMin{"confCPRdeltaPhiCutMin", 0.0, "Delta Phi min cut for Close Pair Rejection"};
  Configurable<float> confCPRdeltaEtaCutMax{"confCPRdeltaEtaCutMax", 0.0, "Delta Eta max cut for Close Pair Rejection"};
  Configurable<float> confCPRdeltaEtaCutMin{"confCPRdeltaEtaCutMin", 0.0, "Delta Eta min cut for Close Pair Rejection"};
  Configurable<float> confCPRChosenRadii{"confCPRChosenRadii", 0.80, "Delta Eta cut for Close Pair Rejection"};
  Configurable<int> confPhiBins{"confPhiBins", 29, "Number of phi bins in deta dphi"};
  Configurable<int> confEtaBins{"confEtaBins", 29, "Number of eta bins in deta dphi"};

  FemtoUniverseContainer<femto_universe_container::EventType::same, femto_universe_container::Observable::kstar> sameEventCont;
  FemtoUniverseContainer<femto_universe_container::EventType::mixed, femto_universe_container::Observable::kstar> mixedEventCont;
  FemtoUniversePairCleaner<aod::femtouniverseparticle::ParticleType::kTrack, aod::femtouniverseparticle::ParticleType::kTrack> pairCleaner;
  FemtoUniverseDetaDphiStar<aod::femtouniverseparticle::ParticleType::kTrack, aod::femtouniverseparticle::ParticleType::kTrack> pairCloseRejection;
  FemtoUniverseTrackSelection trackCuts;
  /// Histogram output
  HistogramRegistry qaRegistry{"TrackQA", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry resultRegistry{"Correlations", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry mixQaRegistry{"mixQaRegistry", {}, OutputObjHandlingPolicy::AnalysisObject};

  EfficiencyConfigurableGroup effConfGroup;
  EfficiencyCalculator<TH1> efficiencyCalculator{&effConfGroup};

  /// @brief Counter for particle swapping
  int fNeventsProcessed = 0;

  // PID for protons
  bool isProtonNSigma(float mom, float nsigmaTPCPr, float nsigmaTOFPr) // previous version from: https://github.com/alisw/AliPhysics/blob/master/PWGCF/FEMTOSCOPY/AliFemtoUser/AliFemtoMJTrackCut.cxx
  {
    //|nsigma_TPC| < 3 for p < 0.5 GeV/c
    //|nsigma_combined| < 3 for p > 0.5

    // using configurables:
    // confTOFPtMin - momentum value when we start using TOF; set to 1000 if TOF not needed
    // confNsigmaTPC -> TPC Sigma for momentum < confTOFPtMin
    // confNsigmaCombined -> TPC and TOF Sigma (combined) for momentum > confTOFPtMin

    if (mom < twotracksconfigs.confTOFPtMin) {
      if (std::abs(nsigmaTPCPr) < twotracksconfigs.confNsigmaTPC) {
        return true;
      } else {
        return false;
      }
    } else {
      if (std::hypot(nsigmaTOFPr, nsigmaTPCPr) < twotracksconfigs.confNsigmaCombined) {
        return true;
      } else {
        return false;
      }
    }
    return false;
  }

  bool isKaonNSigma(float mom, float nsigmaTPCK, float nsigmaTOFK)
  {
    //|nsigma_TPC| < 3 for p < 0.5 GeV/c
    //|nsigma_combined| < 3 for p > 0.5

    // using configurables:
    // confTOFPtMin - momentum value when we start using TOF; set to 1000 if TOF not needed
    // confNsigmaTPC -> TPC Sigma for momentum < confTOFPtMin
    // confNsigmaCombined -> TPC and TOF Sigma (combined) for momentum > confTOFPtMin
    if (true) {
      if (mom < twotracksconfigs.confTOFPtMin) {
        if (std::abs(nsigmaTPCK) < twotracksconfigs.confNsigmaTPC) {
          return true;
        } else {
          return false;
        }
      } else {
        if (std::hypot(nsigmaTOFK, nsigmaTPCK) < twotracksconfigs.confNsigmaCombined) {
          return true;
        } else {
          return false;
        }
      }
    }
    return false;
  }

  bool isPionNSigma(float mom, float nsigmaTPCPi, float nsigmaTOFPi)
  {
    //|nsigma_TPC| < 3 for p < 0.5 GeV/c
    //|nsigma_combined| < 3 for p > 0.5

    // using configurables:
    // confTOFPtMin - momentum value when we start using TOF; set to 1000 if TOF not needed
    // confNsigmaTPC -> TPC Sigma for momentum < confTOFPtMin
    // confNsigmaCombined -> TPC and TOF Sigma (combined) for momentum > confTOFPtMin
    if (true) {
      if (mom < twotracksconfigs.confTOFPtMin) {
        if (std::abs(nsigmaTPCPi) < twotracksconfigs.confNsigmaTPC) {
          return true;
        } else {
          return false;
        }
      } else {
        if (std::hypot(nsigmaTOFPi, nsigmaTPCPi) < twotracksconfigs.confNsigmaCombined) {
          return true;
        } else {
          return false;
        }
      }
    }
    return false;
  }

  bool isParticleNSigma(int8_t particle_number, float mom, float nsigmaTPCPr, float nsigmaTOFPr, float nsigmaTPCPi, float nsigmaTOFPi, float nsigmaTPCK, float nsigmaTOFK)
  {
    if (particle_number == 1) {
      switch (trackonefilter.confPDGCodePartOne) {
        case 2212:  // Proton
        case -2212: // anty Proton
          return isProtonNSigma(mom, nsigmaTPCPr, nsigmaTOFPr);
          break;
        case 211:  // Pion
        case -211: // Pion-
        case 111:  // Pion 0
          return isPionNSigma(mom, nsigmaTPCPi, nsigmaTOFPi);
          break;
        case 321:  // Kaon+
        case -321: // Kaon-
        case 130:  // Kaon 0 LONG
        case 310:  // Kaon 0 SHORT
          return isKaonNSigma(mom, nsigmaTPCK, nsigmaTOFK);
          break;
        default:
          return false;
      }
      return false;
    } else if (particle_number == 2) {
      switch (tracktwofilter.confPDGCodePartTwo) {
        case 2212:  // Proton
        case -2212: // anty Proton
          return isProtonNSigma(mom, nsigmaTPCPr, nsigmaTOFPr);
          break;
        case 211:  // Pion
        case -211: // Pion-
        case 111:  // Pion 0
          return isPionNSigma(mom, nsigmaTPCPi, nsigmaTOFPi);
          break;
        case 321:  // Kaon+
        case -321: // Kaon-
        case 130:  // Kaon 0 LONG
        case 310:  // Kaon 0 SHORT
          return isKaonNSigma(mom, nsigmaTPCK, nsigmaTOFK);
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

  /// @returns 1 if positive, -1 if negative, 0 if zero
  auto sign(auto number) -> int8_t
  {
    return (number > 0) - (number < 0);
  }

  void init(InitContext&)
  {
    if (twotracksconfigs.confIsMC) {
      hMCTruth1.init(&qaRegistry, confTempFitVarpTBins, confTempFitVarPDGBins, false, trackonefilter.confPDGCodePartOne, false);
      if (!confIsSame) {
        hMCTruth2.init(&qaRegistry, confTempFitVarpTBins, confTempFitVarPDGBins, false, tracktwofilter.confPDGCodePartTwo, false);
      }
    }
    efficiencyCalculator.init();

    eventHisto.init(&qaRegistry);
    trackHistoPartOne.init(&qaRegistry, confTempFitVarpTBins, confTempFitVarBins, twotracksconfigs.confIsMC, trackonefilter.confPDGCodePartOne, true); // last true = isDebug
    if (!confIsSame) {
      trackHistoPartTwo.init(&qaRegistry, confTempFitVarpTBins, confTempFitVarBins, twotracksconfigs.confIsMC, tracktwofilter.confPDGCodePartTwo, true); // last true = isDebug
    }

    mixQaRegistry.add("MixingQA/hSECollisionBins", ";bin;Entries", kTH1F, {{120, -0.5, 119.5}});
    mixQaRegistry.add("MixingQA/hMECollisionBins", ";bin;Entries", kTH1F, {{120, -0.5, 119.5}});

    sameEventCont.init(&resultRegistry, confkstarBins, confMultBins, confkTBins, confmTBins, confMultBins3D, confmTBins3D, confEtaBins, confPhiBins, twotracksconfigs.confIsMC, twotracksconfigs.confUse3D);
    mixedEventCont.init(&resultRegistry, confkstarBins, confMultBins, confkTBins, confmTBins, confMultBins3D, confmTBins3D, confEtaBins, confPhiBins, twotracksconfigs.confIsMC, twotracksconfigs.confUse3D);
    sameEventCont.setPDGCodes(trackonefilter.confPDGCodePartOne, tracktwofilter.confPDGCodePartTwo);
    mixedEventCont.setPDGCodes(trackonefilter.confPDGCodePartOne, tracktwofilter.confPDGCodePartTwo);
    pairCleaner.init(&qaRegistry);
    if (confIsCPR.value) {
      pairCloseRejection.init(&resultRegistry, &qaRegistry, confCPRdeltaPhiCutMin.value, confCPRdeltaPhiCutMax.value, confCPRdeltaEtaCutMin.value, confCPRdeltaEtaCutMax.value, confCPRChosenRadii.value, confCPRPlotPerRadii.value);
    }

    vPIDPartOne = trackonefilter.confPIDPartOne.value;
    vPIDPartTwo = tracktwofilter.confPIDPartTwo.value;
    kNsigma = twotracksconfigs.confTrkPIDnSigmaMax.value;
  }

  template <typename CollisionType>
  void fillCollision(CollisionType col)
  {
    mixQaRegistry.fill(HIST("MixingQA/hSECollisionBins"), colBinning.getBin({col.posZ(), col.multNtr()}));
    eventHisto.fillQA(col);
  }

  template <uint8_t N>
    requires isOneOrTwo<N>
  auto doMCTruth(FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kMCTruthTrack, N> hist, auto parts) -> void
  {
    auto expectedPDG = 0;
    auto expectedCharge = 0.0l;

    if constexpr (N == ParticleNo::ONE) {
      expectedPDG = trackonefilter.confPDGCodePartOne;
      expectedCharge = trackonefilter.confChargePart1;
    } else if constexpr (N == ParticleNo::TWO) {
      expectedPDG = tracktwofilter.confPDGCodePartTwo;
      expectedCharge = tracktwofilter.confChargePart2;
    }

    for (const auto& particle : parts) {
      auto pdgCode = static_cast<int>(particle.pidCut());
      if (pdgCode != expectedPDG) {
        continue;
      }

      const auto& pdgParticle = pdgMC->GetParticle(pdgCode);
      if (!pdgParticle) {
        continue;
      }

      if (sign(pdgParticle->Charge()) == sign(expectedCharge)) {
        hist.template fillQA<false, false>(particle);
      }
    }
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
  void doSameEvent(PartitionType groupPartsOne, PartitionType groupPartsTwo, PartType parts, float magFieldTesla, int multCol)
  {
    // variables for particle swapping
    bool swpart = fNeventsProcessed % 2;
    fNeventsProcessed++;

    /// Histogramming same event
    for (const auto& part : groupPartsOne) {
      // if (part.p() > twotracksconfigs.confCutTable->get("PartOne", "MaxP") || part.pt() > twotracksconfigs.confCutTable->get("PartOne", "MaxPt")) {
      //   continue;
      // }
      // if (!isFullPIDSelected(part.pidCut(),
      //                        part.p(),
      //                        twotracksconfigs.confCutTable->get("PartOne", "PIDthr"),
      //                        vPIDPartOne,
      //                        twotracksconfigs.confNspecies,
      //                        kNsigma,
      //                        twotracksconfigs.confCutTable->get("PartOne", "nSigmaTPC"),
      //                        twotracksconfigs.confCutTable->get("PartOne", "nSigmaTPCTOF"))) {
      //   continue;
      // }

      if (trackonefilter.confIsTrackOneIdentified) {
        if (!isParticleNSigma((int8_t)1, part.p(), trackCuts.getNsigmaTPC(part, o2::track::PID::Proton), trackCuts.getNsigmaTOF(part, o2::track::PID::Proton), trackCuts.getNsigmaTPC(part, o2::track::PID::Pion), trackCuts.getNsigmaTOF(part, o2::track::PID::Pion), trackCuts.getNsigmaTPC(part, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(part, o2::track::PID::Kaon))) {
          continue;
        }
      }

      trackHistoPartOne.fillQA<isMC, true>(part);
    }

    if (!confIsSame) {
      for (const auto& part : groupPartsTwo) {
        // if (part.p() > twotracksconfigs.confCutTable->get("PartTwo", "MaxP") || part.pt() > twotracksconfigs.confCutTable->get("PartTwo", "MaxPt")) {
        //   continue;
        // }
        // if (!isFullPIDSelected(part.pidCut(),
        //                        part.p(),
        //                        twotracksconfigs.confCutTable->get("PartTwo", "PIDthr"),
        //                        vPIDPartTwo,
        //                        twotracksconfigs.confNspecies,
        //                        kNsigma,
        //                        twotracksconfigs.confCutTable->get("PartTwo", "nSigmaTPC"),
        //                        twotracksconfigs.confCutTable->get("PartTwo", "nSigmaTPCTOF"))) {
        //   continue;
        // }

        if (tracktwofilter.confIsTrackTwoIdentified) {
          if (!isParticleNSigma((int8_t)2, part.p(), trackCuts.getNsigmaTPC(part, o2::track::PID::Proton), trackCuts.getNsigmaTOF(part, o2::track::PID::Proton), trackCuts.getNsigmaTPC(part, o2::track::PID::Pion), trackCuts.getNsigmaTOF(part, o2::track::PID::Pion), trackCuts.getNsigmaTPC(part, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(part, o2::track::PID::Kaon))) {
            continue;
          }
        }
        trackHistoPartTwo.fillQA<isMC, true>(part);
      }

      /// Now build the combinations for non-identical particle pairs
      for (const auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(groupPartsOne, groupPartsTwo))) {
        // if (p1.p() > twotracksconfigs.confCutTable->get("PartOne", "MaxP") || p1.pt() > twotracksconfigs.confCutTable->get("PartOne", "MaxPt") || p2.p() > twotracksconfigs.confCutTable->get("PartTwo", "MaxP") || p2.pt() > twotracksconfigs.confCutTable->get("PartTwo", "MaxPt")) {
        //   continue;
        // }
        // if (!isFullPIDSelected(p1.pidCut(),
        //                        p1.p(),
        //                        twotracksconfigs.confCutTable->get("PartOne", "PIDthr"),
        //                        vPIDPartOne,
        //                        twotracksconfigs.confNspecies,
        //                        kNsigma,
        //                        twotracksconfigs.confCutTable->get("PartOne", "nSigmaTPC"),
        //                        twotracksconfigs.confCutTable->get("PartOne", "nSigmaTPCTOF")) ||
        //     !isFullPIDSelected(p2.pidCut(),
        //                        p2.p(),
        //                        twotracksconfigs.confCutTable->get("PartTwo", "PIDthr"),
        //                        vPIDPartTwo,
        //                        twotracksconfigs.confNspecies,
        //                        kNsigma,
        //                        twotracksconfigs.confCutTable->get("PartTwo", "nSigmaTPC"),
        //                        twotracksconfigs.confCutTable->get("PartTwo", "nSigmaTPCTOF"))) {
        //   continue;
        // }

        if (trackonefilter.confIsTrackOneIdentified) {
          if (!isParticleNSigma((int8_t)1, p1.p(), trackCuts.getNsigmaTPC(p1, o2::track::PID::Proton), trackCuts.getNsigmaTOF(p1, o2::track::PID::Proton), trackCuts.getNsigmaTPC(p1, o2::track::PID::Pion), trackCuts.getNsigmaTOF(p1, o2::track::PID::Pion), trackCuts.getNsigmaTPC(p1, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(p1, o2::track::PID::Kaon))) {
            continue;
          }
        }

        if (tracktwofilter.confIsTrackTwoIdentified) {
          if (!isParticleNSigma((int8_t)2, p2.p(), trackCuts.getNsigmaTPC(p2, o2::track::PID::Proton), trackCuts.getNsigmaTOF(p2, o2::track::PID::Proton), trackCuts.getNsigmaTPC(p2, o2::track::PID::Pion), trackCuts.getNsigmaTOF(p2, o2::track::PID::Pion), trackCuts.getNsigmaTPC(p2, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(p2, o2::track::PID::Kaon))) {
            continue;
          }
        }

        if (confIsCPR.value) {
          if (pairCloseRejection.isClosePair(p1, p2, parts, magFieldTesla, femto_universe_container::EventType::same)) {
            continue;
          }
        }

        // track cleaning
        if (!pairCleaner.isCleanPair(p1, p2, parts)) {
          continue;
        }

        float weight = efficiencyCalculator.getWeight(ParticleNo::ONE, p1.pt());
        if (!confIsSame) {
          weight *= efficiencyCalculator.getWeight(ParticleNo::TWO, p2.pt());
        }

        if (swpart)
          sameEventCont.setPair<isMC>(p1, p2, multCol, twotracksconfigs.confUse3D, weight);
        else
          sameEventCont.setPair<isMC>(p2, p1, multCol, twotracksconfigs.confUse3D, weight);

        swpart = !swpart;
      }
    } else {

      /// Now build the combinations for identical particle pairs (different combination policy than for non-identical!)
      for (const auto& [p1, p2] : combinations(CombinationsStrictlyUpperIndexPolicy(groupPartsOne, groupPartsTwo))) {
        // if (p1.p() > twotracksconfigs.confCutTable->get("PartOne", "MaxP") || p1.pt() > twotracksconfigs.confCutTable->get("PartOne", "MaxPt") || p2.p() > twotracksconfigs.confCutTable->get("PartTwo", "MaxP") || p2.pt() > twotracksconfigs.confCutTable->get("PartTwo", "MaxPt")) {
        //   continue;
        // }
        // if (!isFullPIDSelected(p1.pidCut(),
        //                        p1.p(),
        //                        twotracksconfigs.confCutTable->get("PartOne", "PIDthr"),
        //                        vPIDPartOne,
        //                        twotracksconfigs.confNspecies,
        //                        kNsigma,
        //                        twotracksconfigs.confCutTable->get("PartOne", "nSigmaTPC"),
        //                        twotracksconfigs.confCutTable->get("PartOne", "nSigmaTPCTOF")) ||
        //     !isFullPIDSelected(p2.pidCut(),
        //                        p2.p(),
        //                        twotracksconfigs.confCutTable->get("PartTwo", "PIDthr"),
        //                        vPIDPartTwo,
        //                        twotracksconfigs.confNspecies,
        //                        kNsigma,
        //                        twotracksconfigs.confCutTable->get("PartTwo", "nSigmaTPC"),
        //                        twotracksconfigs.confCutTable->get("PartTwo", "nSigmaTPCTOF"))) {
        //   continue;
        // }

        if (trackonefilter.confIsTrackOneIdentified) {
          if (!isParticleNSigma((int8_t)1, p1.p(), trackCuts.getNsigmaTPC(p1, o2::track::PID::Proton), trackCuts.getNsigmaTOF(p1, o2::track::PID::Proton), trackCuts.getNsigmaTPC(p1, o2::track::PID::Pion), trackCuts.getNsigmaTOF(p1, o2::track::PID::Pion), trackCuts.getNsigmaTPC(p1, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(p1, o2::track::PID::Kaon))) {
            continue;
          }
        }

        if (tracktwofilter.confIsTrackTwoIdentified) {
          if (!isParticleNSigma((int8_t)2, p2.p(), trackCuts.getNsigmaTPC(p2, o2::track::PID::Proton), trackCuts.getNsigmaTOF(p2, o2::track::PID::Proton), trackCuts.getNsigmaTPC(p2, o2::track::PID::Pion), trackCuts.getNsigmaTOF(p2, o2::track::PID::Pion), trackCuts.getNsigmaTPC(p2, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(p2, o2::track::PID::Kaon))) {
            continue;
          }
        }

        if (confIsCPR.value) {
          if (pairCloseRejection.isClosePair(p1, p2, parts, magFieldTesla, femto_universe_container::EventType::same)) {
            continue;
          }
        }

        // track cleaning
        if (!pairCleaner.isCleanPair(p1, p2, parts)) {
          continue;
        }

        float weight = efficiencyCalculator.getWeight(ParticleNo::ONE, p1.pt());
        if (!confIsSame) {
          weight *= efficiencyCalculator.getWeight(ParticleNo::TWO, p2.pt());
        }

        sameEventCont.setPair<isMC>(p1, p2, multCol, twotracksconfigs.confUse3D, weight);
      }
    }
  }

  /// process function for to call doSameEvent with Data
  /// \param col subscribe to the collision table (Data)
  /// \param parts subscribe to the femtoUniverseParticleTable
  void processSameEvent(const o2::aod::FdCollision& col,
                        const FilteredFemtoFullParticles& parts)
  {
    fillCollision(col);

    auto thegroupPartsOne = partsOne->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto thegroupPartsTwo = partsTwo->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);

    doSameEvent<false>(thegroupPartsOne, thegroupPartsTwo, parts, col.magField(), col.multNtr());
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackTrackExtended, processSameEvent, "Enable processing same event", true);

  /// process function for to call doSameEvent with Monte Carlo
  /// \param col subscribe to the collision table (Monte Carlo Reconstructed reconstructed)
  /// \param parts subscribe to joined table FemtoUniverseParticles and FemtoUniverseMCLables to access Monte Carlo truth
  /// \param FemtoUniverseMCParticles subscribe to the Monte Carlo truth table
  void processSameEventMC(const o2::aod::FdCollision& col,
                          const soa::Join<FilteredFemtoFullParticles, aod::FDMCLabels>& parts,
                          const o2::aod::FdMCParticles&)
  {
    fillCollision(col);

    auto groupMCTruth1 = partsOneMCTruth->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    doMCTruth<1>(hMCTruth1, groupMCTruth1);

    if (!confIsSame) {
      auto groupMCTruth2 = partsTwoMCTruth->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
      doMCTruth<2>(hMCTruth2, groupMCTruth2);
    }

    auto groupMCReco1 = partsOneMCReco->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto groupMCReco2 = partsTwoMCReco->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);

    doSameEvent<true>(groupMCReco1, groupMCReco2, parts, col.magField(), col.multNtr());
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackTrackExtended, processSameEventMC, "Enable processing same event for Monte Carlo", false);

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
  template <bool isMC, typename PartitionType, typename PartType>
  void doMixedEvent(PartitionType groupPartsOne, PartitionType groupPartsTwo, PartType parts, float magFieldTesla, int multCol)
  {

    // variables for particle swapping
    bool swpart = fNeventsProcessed % 2;
    fNeventsProcessed++;

    for (const auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(groupPartsOne, groupPartsTwo))) {
      // if (p1.p() > twotracksconfigs.confCutTable->get("PartOne", "MaxP") || p1.pt() > twotracksconfigs.confCutTable->get("PartOne", "MaxPt") || p2.p() > twotracksconfigs.confCutTable->get("PartTwo", "MaxP") || p2.pt() > twotracksconfigs.confCutTable->get("PartTwo", "MaxPt")) {
      //   continue;
      // }
      // if (!isFullPIDSelected(p1.pidCut(),
      //                        p1.p(),
      //                        twotracksconfigs.confCutTable->get("PartOne", "PIDthr"),
      //                        vPIDPartOne,
      //                        twotracksconfigs.confNspecies,
      //                        kNsigma,
      //                        twotracksconfigs.confCutTable->get("PartOne", "nSigmaTPC"),
      //                        twotracksconfigs.confCutTable->get("PartOne", "nSigmaTPCTOF")) ||
      //     !isFullPIDSelected(p2.pidCut(),
      //                        p2.p(),
      //                        twotracksconfigs.confCutTable->get("PartTwo", "PIDthr"),
      //                        vPIDPartTwo,
      //                        twotracksconfigs.confNspecies,
      //                        kNsigma,
      //                        twotracksconfigs.confCutTable->get("PartTwo", "nSigmaTPC"),
      //                        twotracksconfigs.confCutTable->get("PartTwo", "nSigmaTPCTOF"))) {
      //   continue;
      // }
      if (trackonefilter.confIsTrackOneIdentified) {
        if (!isParticleNSigma((int8_t)1, p1.p(), trackCuts.getNsigmaTPC(p1, o2::track::PID::Proton), trackCuts.getNsigmaTOF(p1, o2::track::PID::Proton), trackCuts.getNsigmaTPC(p1, o2::track::PID::Pion), trackCuts.getNsigmaTOF(p1, o2::track::PID::Pion), trackCuts.getNsigmaTPC(p1, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(p1, o2::track::PID::Kaon))) {
          continue;
        }
      }

      if (tracktwofilter.confIsTrackTwoIdentified) {
        if (!isParticleNSigma((int8_t)2, p2.p(), trackCuts.getNsigmaTPC(p2, o2::track::PID::Proton), trackCuts.getNsigmaTOF(p2, o2::track::PID::Proton), trackCuts.getNsigmaTPC(p2, o2::track::PID::Pion), trackCuts.getNsigmaTOF(p2, o2::track::PID::Pion), trackCuts.getNsigmaTPC(p2, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(p2, o2::track::PID::Kaon))) {
          continue;
        }
      }

      if (confIsCPR.value) {
        if (pairCloseRejection.isClosePair(p1, p2, parts, magFieldTesla, femto_universe_container::EventType::mixed)) {
          continue;
        }
      }

      float weight = efficiencyCalculator.getWeight(ParticleNo::ONE, p1.pt());
      if (!confIsSame) {
        weight *= efficiencyCalculator.getWeight(ParticleNo::TWO, p2.pt());
      }

      if (swpart)
        mixedEventCont.setPair<isMC>(p1, p2, multCol, twotracksconfigs.confUse3D, weight);
      else
        mixedEventCont.setPair<isMC>(p2, p1, multCol, twotracksconfigs.confUse3D, weight);

      swpart = !swpart;
    }
  }

  /// process function for to call doMixedEvent with Data
  /// @param cols subscribe to the collisions table (Data)
  /// @param parts subscribe to the femtoUniverseParticleTable
  void processMixedEvent(const o2::aod::FdCollisions& cols,
                         const FilteredFemtoFullParticles& parts)
  {
    for (const auto& [collision1, collision2] : soa::selfCombinations(colBinning, 5, -1, cols, cols)) {

      const int multiplicityCol = collision1.multNtr();
      mixQaRegistry.fill(HIST("MixingQA/hMECollisionBins"), colBinning.getBin({collision1.posZ(), multiplicityCol}));

      auto groupPartsOne = partsOne->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);
      auto groupPartsTwo = partsTwo->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision2.globalIndex(), cache);

      const auto& magFieldTesla1 = collision1.magField();
      const auto& magFieldTesla2 = collision2.magField();

      if (magFieldTesla1 != magFieldTesla2) {
        continue;
      }
      /// \todo before mixing we should check whether both collisions contain a pair of particles!
      // if (partsOne.size() == 0 || nPart2Evt1 == 0 || nPart1Evt2 == 0 || partsTwo.size() == 0 ) continue;

      doMixedEvent<false>(groupPartsOne, groupPartsTwo, parts, magFieldTesla1, multiplicityCol);
    }
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackTrackExtended, processMixedEvent, "Enable processing mixed events", true);

  /// brief process function for to call doMixedEvent with Monte Carlo
  /// @param cols subscribe to the collisions table (Monte Carlo Reconstructed reconstructed)
  /// @param parts subscribe to joined table FemtoUniverseParticles and FemtoUniverseMCLables to access Monte Carlo truth
  /// @param FemtoUniverseMCParticles subscribe to the Monte Carlo truth table
  void processMixedEventMC(const o2::aod::FdCollisions& cols,
                           const soa::Join<FilteredFemtoFullParticles, aod::FDMCLabels>& parts,
                           const o2::aod::FdMCParticles&)
  {
    for (const auto& [collision1, collision2] : soa::selfCombinations(colBinning, 5, -1, cols, cols)) {

      const int multiplicityCol = collision1.multNtr();
      mixQaRegistry.fill(HIST("MixingQA/hMECollisionBins"), colBinning.getBin({collision1.posZ(), multiplicityCol}));

      auto groupPartsOne = partsOneMCReco->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);
      auto groupPartsTwo = partsTwoMCReco->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision2.globalIndex(), cache);

      const auto& magFieldTesla1 = collision1.magField();
      const auto& magFieldTesla2 = collision2.magField();

      if (magFieldTesla1 != magFieldTesla2) {
        continue;
      }
      /// \todo before mixing we should check whether both collisions contain a pair of particles!
      // if (partsOne.size() == 0 || nPart2Evt1 == 0 || nPart1Evt2 == 0 || partsTwo.size() == 0 ) continue;

      doMixedEvent<true>(groupPartsOne, groupPartsTwo, parts, magFieldTesla1, multiplicityCol);
    }
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackTrackExtended, processMixedEventMC, "Enable processing mixed events MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<FemtoUniversePairTaskTrackTrackExtended>(cfgc),
  };
  return workflow;
}
