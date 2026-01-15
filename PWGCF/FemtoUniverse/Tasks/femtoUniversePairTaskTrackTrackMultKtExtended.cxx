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

/// \file femtoUniversePairTaskTrackTrackMultKtExtended.cxx
/// \brief Tasks that reads the track tables used for the pairing and builds pairs of two tracks
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de
/// \author Georgios Mantzaridis, TU München, georgios.mantzaridis@tum.de
/// \author Anton Riedel, TU München, anton.riedel@tum.de
/// \author Zuzanna Chochulska, WUT Warsaw & CTU Prague, zchochul@cern.ch
/// \author Alicja Płachta, WUT Warsaw, alicja.plachta.stud@pw.edu.pl

#include "PWGCF/FemtoUniverse/Core/FemtoUniverseDetaDphiStar.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseEventHisto.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseFemtoContainer.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseMath.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniversePairCleaner.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniversePairWithCentMultKt.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseParticleHisto.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseTrackSelection.h"
#include "PWGCF/FemtoUniverse/DataModel/FemtoDerived.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/PID.h"

#include <string>
#include <vector>

using namespace o2;
using namespace o2::analysis::femto_universe;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

namespace
{
static constexpr int Npart = 2;
static constexpr int Ncuts = 5;
static const std::vector<std::string> partNames{"PartOne", "PartTwo"};
static const std::vector<std::string> cutNames{"MaxPt", "PIDthr", "nSigmaTPC", "nSigmaTPCTOF", "MaxP"};
static const float cutsTable[Npart][Ncuts]{{4.05f, 1.f, 3.f, 3.f, 100.f}, {4.05f, 1.f, 3.f, 3.f, 100.f}};
} // namespace

struct FemtoUniversePairTaskTrackTrackMultKtExtended {

  Service<o2::framework::O2DatabasePDG> pdg;

  /// Particle selection part

  /// Table for both particles
  struct : o2::framework::ConfigurableGroup {
    Configurable<bool> isKaonNsigma{"isKaonNsigma", false, "Enable a strict cut selection for K+ and K-"};
    Configurable<float> confNsigmaCombined{"confNsigmaCombined", 3.0f, "TPC and TOF Pion Sigma (combined) for momentum > confTOFpMin"};
    Configurable<float> confNsigmaTPC{"confNsigmaTPC", 3.0f, "TPC Pion Sigma for momentum < confTOFpMin"};
    Configurable<float> confTOFpMin{"confTOFpMin", 0.5f, "Min. momentum for which TOF is required for PID."};
    Configurable<float> confEtaMax{"confEtaMax", 0.8f, "Higher limit for |Eta| (the same for both particles)"};

    Configurable<LabeledArray<float>> confCutTable{"confCutTable", {cutsTable[0], Npart, Ncuts, partNames, cutNames}, "Particle selections"};
    Configurable<int> confNspecies{"confNspecies", 2, "Number of particle spieces with PID info"};
    Configurable<bool> confIsMC{"confIsMC", false, "Enable additional Histogramms in the case of a MonteCarlo Run"};
    Configurable<std::vector<float>> confTrkPIDnSigmaMax{"confTrkPIDnSigmaMax", std::vector<float>{4.f, 3.f, 2.f}, "This configurable needs to be the same as the one used in the producer task"};
    Configurable<bool> confUse3D{"confUse3D", false, "Enable three dimensional histogramms (to be used only for analysis with high statistics): k* vs mT vs multiplicity"};
  } twotracksconfigs;

  using FemtoFullParticles = soa::Join<aod::FDParticles, aod::FDExtParticles>;
  // Filters for selecting particles (both p1 and p2)
  Filter trackAdditionalfilter = (nabs(aod::femtouniverseparticle::eta) < twotracksconfigs.confEtaMax); // example filtering on configurable
  using FilteredFemtoFullParticles = soa::Filtered<FemtoFullParticles>;
  // using FilteredFemtoFullParticles = FemtoFullParticles; //if no filtering is applied uncomment this option

  SliceCache cache;
  Preslice<FilteredFemtoFullParticles> perCol = aod::femtouniverseparticle::fdCollisionId;

  using FemtoTruthParticles = soa::Filtered<soa::Join<aod::FDParticles, aod::FDExtParticles, aod::FDMCLabels>>;
  Preslice<FemtoTruthParticles> perColMCTruth = aod::femtouniverseparticle::fdCollisionId;

  using FemtoRecoParticles = soa::Filtered<soa::Join<aod::FDParticles, aod::FDExtParticles, aod::FDMCLabels, aod::FDExtMCParticles>>;
  Preslice<FemtoRecoParticles> perColMC = aod::femtouniverseparticle::fdCollisionId;

  /// Particle 1
  struct : o2::framework::ConfigurableGroup {
    Configurable<int> confPDGCodePartOne{"confPDGCodePartOne", 211, "Particle 1 -- PDG code"};
    // Configurable<uint32_t> confCutPartOne{"confCutPartOne", 5542474, "Particle 1 -- Selection bit from cutCulator"};
    Configurable<int> confPIDPartOne{"confPIDPartOne", 2, "Particle 1 -- Read from cutCulator"}; // we also need the possibility to specify whether the bit is true/false ->std>>vector<std::pair<int, int>>int>>
    Configurable<float> confpLowPart1{"confpLowPart1", 0.14, "Lower limit for Pt for the first particle"};
    Configurable<float> confPtHighPart1{"confPtHighPart1", 1.5, "Higher limit for Pt for the first particle"};
    Configurable<int> confChargePart1{"confChargePart1", 1, "Particle 1 sign"};
  } trackonefilter;

  /// Partition for particle 1
  Partition<FilteredFemtoFullParticles> partsOne = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kTrack)) && aod::femtouniverseparticle::sign == as<int8_t>(trackonefilter.confChargePart1) && aod::femtouniverseparticle::pt < trackonefilter.confPtHighPart1 && aod::femtouniverseparticle::pt > trackonefilter.confpLowPart1;

  Partition<FemtoRecoParticles> partsOneMC = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kTrack)) && aod::femtouniverseparticle::sign == as<int8_t>(trackonefilter.confChargePart1) && aod::femtouniverseparticle::pt < trackonefilter.confPtHighPart1 && aod::femtouniverseparticle::pt > trackonefilter.confpLowPart1;

  Partition<FemtoTruthParticles> partsOneMCTruth = aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kMCTruthTrack) && aod::femtouniverseparticle::pt < trackonefilter.confPtHighPart1 && aod::femtouniverseparticle::pt > trackonefilter.confpLowPart1;
  ;

  /// Histogramming for particle 1
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kTrack, 1> trackHistoPartOne;

  /// Particle 2
  struct : o2::framework::ConfigurableGroup {
    Configurable<int> confPDGCodePartTwo{"confPDGCodePartTwo", 211, "Particle 2 -- PDG code"};
    // Configurable<uint32_t> confCutPartTwo{"confCutPartTwo", 5542474, "Particle 2 -- Selection bit"};
    Configurable<int> confPIDPartTwo{"confPIDPartTwo", 2, "Particle 2 -- Read from cutCulator"}; // we also need the possibility to specify whether the bit is true/false ->std>>vector<std::pair<int, int>>
    Configurable<float> confpLowPart2{"confpLowPart2", 0.14, "Lower limit for Pt for the second particle"};
    Configurable<float> confPtHighPart2{"confPtHighPart2", 1.5, "Higher limit for Pt for the second particle"};
    Configurable<int> confChargePart2{"confChargePart2", -1, "Particle 2 sign"};
  } tracktwofilter;

  /// Partition for particle 2
  Partition<FilteredFemtoFullParticles> partsTwo = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kTrack)) && (aod::femtouniverseparticle::sign == as<int8_t>(tracktwofilter.confChargePart2)) && aod::femtouniverseparticle::pt < tracktwofilter.confPtHighPart2 && aod::femtouniverseparticle::pt > tracktwofilter.confpLowPart2;

  Partition<FemtoRecoParticles> partsTwoMC = aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kTrack) && (aod::femtouniverseparticle::sign == as<int8_t>(tracktwofilter.confChargePart2)) && aod::femtouniverseparticle::pt < tracktwofilter.confPtHighPart2 && aod::femtouniverseparticle::pt > tracktwofilter.confpLowPart2;

  Partition<FemtoTruthParticles> partsTwoMCTruth = aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kMCTruthTrack) && aod::femtouniverseparticle::pt < tracktwofilter.confPtHighPart2 && aod::femtouniverseparticle::pt > tracktwofilter.confpLowPart2;

  /// Histogramming for particle 2
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kTrack, 2> trackHistoPartTwo;

  /// Histogramming for Event
  FemtoUniverseEventHisto eventHisto;

  /// The configurables need to be passed to an std::vector
  int vPIDPartOne, vPIDPartTwo;
  std::vector<float> kNsigma;

  /// Event part
  struct : o2::framework::ConfigurableGroup {
    Configurable<float> confV0MLow{"confV0MLow", 0.0, "Lower limit for V0M multiplicity"};
    Configurable<float> confV0MHigh{"confV0MHigh", 25000.0, "Upper limit for V0M multiplicity"};
    Configurable<float> confSphericityCutMin{"confSphericityCutMin", 0, "Min. sphericity"};
    Configurable<float> confSphericityCutMax{"confSphericityCutMax", 3, "Max. sphericity"};
  } eventSel;

  Filter collV0Mfilter = ((o2::aod::femtouniversecollision::multV0M > eventSel.confV0MLow) && (o2::aod::femtouniversecollision::multV0M < eventSel.confV0MHigh));
  Filter colSpherfilter = ((o2::aod::femtouniversecollision::sphericity > eventSel.confSphericityCutMin) && (o2::aod::femtouniversecollision::sphericity < eventSel.confSphericityCutMax));
  // Filter trackAdditionalfilter = (nabs(aod::femtouniverseparticle::eta) < twotracksconfigs.confEtaMax); // example filtering on configurable

  /// Particle part
  ConfigurableAxis confTempFitVarBins{"confTempFitVarBins", {300, -0.15, 0.15}, "binning of the TempFitVar in the pT vs. TempFitVar plot"};
  ConfigurableAxis confTempFitVarpTBins{"confTempFitVarpTBins", {20, 0.5, 4.05}, "pT binning of the pT vs. TempFitVar plot"};

  /// Correlation part
  ConfigurableAxis confMultBins{"confMultBins", {VARIABLE_WIDTH, 0.0f, 4.0f, 8.0f, 12.0f, 16.0f, 20.0f, 24.0f, 28.0f, 32.0f, 36.0f, 40.0f, 44.0f, 48.0f, 52.0f, 56.0f, 60.0f, 64.0f, 68.0f, 72.0f, 76.0f, 80.0f, 84.0f, 88.0f, 92.0f, 96.0f, 100.0f, 200.0f, 99999.f}, "Mixing bins - multiplicity or centrality"}; // \todo to be obtained from the hash task
  ConfigurableAxis confMultKstarBins{"confMultKstarBins", {VARIABLE_WIDTH, 0.0f, 13.0f, 20.0f, 30.0f, 40.0f, 50.0f, 100.0f, 99999.f}, "Bins for kstar analysis in multiplicity or centrality bins (10 is maximum)"};
  ConfigurableAxis confKtKstarBins{"confKtKstarBins", {VARIABLE_WIDTH, 0.0f, 0.2f, 0.4f, 0.6f, 0.8f, 1.0f, 2.0f, 99999.f}, "Bins for kstar analysis in kT bins (10 is maximum)"};
  ConfigurableAxis confVtxBins{"confVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};

  ConfigurableAxis confmTBins3D{"confmTBins3D", {VARIABLE_WIDTH, 1.02f, 1.14f, 1.20f, 1.26f, 1.38f, 1.56f, 1.86f, 4.50f}, "mT Binning for the 3Dimensional plot: k* vs multiplicity vs mT (set <<twotracksconfigs.confUse3D>> to true in order to use)"};
  ConfigurableAxis confmultBins3D{"confmultBins3D", {VARIABLE_WIDTH, 0.0f, 20.0f, 30.0f, 40.0f, 99999.0f}, "multiplicity Binning for the 3Dimensional plot: k* vs multiplicity vs mT (set <<twotracksconfigs.confUse3D>> to true in order to use)"};

  ColumnBinningPolicy<aod::collision::PosZ, aod::femtouniversecollision::MultV0M> colBinning{{confVtxBins, confMultBins}, true};

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

  Configurable<bool> isPairIdentical{"isPairIdentical", true, "'true' for identical particles, 'false' for non-identical particles"};
  struct : o2::framework::ConfigurableGroup {
    Configurable<bool> cfgProcessPM{"cfgProcessPM", true, "Process particles of the opposite charge"};
    Configurable<bool> cfgProcessPP{"cfgProcessPP", false, "Process particles of the same, positice charge"};
    Configurable<bool> cfgProcessMM{"cfgProcessMM", false, "Process particles of the same, positice charge"};
  } processPair;
  Configurable<bool> confFillDebug{"confFillDebug", false, "Fill additional histograms: event mixing, sphericity"};
  Configurable<bool> cfgProcessMultBins{"cfgProcessMultBins", false, "Process kstar histograms (in multiplicity bins)"};
  Configurable<bool> cfgProcessKtBins{"cfgProcessKtBins", false, "Process kstar histograms in kT bins (if 'cfgProcessMultBins' is false, it will not be processed regardless of 'cfgProcessKtBins' state)"};
  Configurable<bool> cfgProcessKtMt3DCF{"cfgProcessKtMt3DCF", false, "Process 3D histograms in kT and MultBins"};

  FemtoUniverseFemtoContainer<femto_universe_femto_container::EventType::same, femto_universe_femto_container::Observable::kstar> sameEventCont;
  FemtoUniverseFemtoContainer<femto_universe_femto_container::EventType::mixed, femto_universe_femto_container::Observable::kstar> mixedEventCont;

  FemtoUniverseFemtoContainer<femto_universe_femto_container::EventType::same, femto_universe_femto_container::Observable::kstar> sameEventContPP;
  FemtoUniverseFemtoContainer<femto_universe_femto_container::EventType::mixed, femto_universe_femto_container::Observable::kstar> mixedEventContPP;

  FemtoUniverseFemtoContainer<femto_universe_femto_container::EventType::same, femto_universe_femto_container::Observable::kstar> sameEventContMM;
  FemtoUniverseFemtoContainer<femto_universe_femto_container::EventType::mixed, femto_universe_femto_container::Observable::kstar> mixedEventContMM;

  FemtoUniversePairCleaner<aod::femtouniverseparticle::ParticleType::kTrack, aod::femtouniverseparticle::ParticleType::kTrack> pairCleaner;
  FemtoUniverseDetaDphiStar<aod::femtouniverseparticle::ParticleType::kTrack, aod::femtouniverseparticle::ParticleType::kTrack> pairCloseRejection;
  FemtoUniverseTrackSelection trackCuts;

  FemtoUniversePairWithCentMultKt sameEventMultCont;
  FemtoUniversePairWithCentMultKt mixedEventMultCont;

  FemtoUniversePairWithCentMultKt sameEventMultContPP;
  FemtoUniversePairWithCentMultKt mixedEventMultContPP;

  FemtoUniversePairWithCentMultKt sameEventMultContMM;
  FemtoUniversePairWithCentMultKt mixedEventMultContMM;

  float mass1 = -1;
  float mass2 = -1;

  /// Histogram output
  HistogramRegistry qaRegistry{"TrackQA", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry resultRegistry{"Correlations", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry resultRegistryPM{"CorrelationsPM", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry resultRegistryPP{"CorrelationsPP", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry resultRegistryMM{"CorrelationsMM", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry mixQaRegistry{"mixQaRegistry", {}, OutputObjHandlingPolicy::AnalysisObject};

  HistogramRegistry sameMultRegistryPM{"sameMultRegistryPM", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry mixedMultRegistryPM{"mixedMultRegistryPM", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  HistogramRegistry sameMultRegistryPP{"sameMultRegistryPP", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry mixedMultRegistryPP{"mixedMultRegistryPP", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  HistogramRegistry sameMultRegistryMM{"sameMultRegistryMM", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry mixedMultRegistryMM{"mixedMultRegistryMM", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  HistogramRegistry sphericityRegistry{"SphericityHisto", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  /// TPC Pion/Kaon/Proton Sigma selection (general)
  bool isNSigma(float mom, float nsigmaTPC, float nsigmaTOF)
  {
    // |nsigma_TPC| < 3 for p < 0.5 GeV/c
    // |nsigma_combined| < 3 for p > 0.5

    // using configurables:
    // confTOFpMin -- momentum value when we start using TOF; set to 1000 if TOF not needed
    // confNsigmaTPC -> TPC Sigma for momentum < confTOFpMin
    // confNsigmaCombined -> TPC and TOF Sigma (combined) for momentum > confTOFpMin

    if (mom < twotracksconfigs.confTOFpMin) {
      return std::abs(nsigmaTPC) < twotracksconfigs.confNsigmaTPC;
    } else {
      return std::hypot(nsigmaTOF, nsigmaTPC) < twotracksconfigs.confNsigmaCombined;
    }
  }

  /// TPC Kaon Sigma selection (stricter cuts for K+ and K-) -- based on Run2 results
  bool isKaonNsigma(float mom, float nsigmaTPCK, float nsigmaTOFK)
  {
    double nSigmaTPCRanges[3] = {1, 2, 3};
    double nSigmaTOFRanges[2] = {1.5, 2};
    double momRanges[4] = {0.4, 0.45, 0.8, 1.5};

    if (twotracksconfigs.isKaonNsigma == true) {
      if (mom < momRanges[0]) {
        return std::abs(nsigmaTPCK) < nSigmaTPCRanges[1];
      } else if (mom > momRanges[0] && mom < momRanges[1]) {
        return std::abs(nsigmaTPCK) < nSigmaTPCRanges[0];
      } else if (mom > momRanges[1] && mom < momRanges[2]) {
        return (std::abs(nsigmaTPCK) < nSigmaTPCRanges[2] && std::abs(nsigmaTOFK) < nSigmaTOFRanges[1]);
      } else if (mom > momRanges[2] && mom < momRanges[3]) {
        return (std::abs(nsigmaTPCK) < nSigmaTPCRanges[2] && std::abs(nsigmaTOFK) < nSigmaTOFRanges[0]);
      } else {
        return false;
      }
    } else {
      return isNSigma(mom, nsigmaTPCK, nsigmaTOFK);
    }
  }

  bool isParticleNSigma(int8_t particle_number, float mom, float nsigmaTPCPr, float nsigmaTOFPr, float nsigmaTPCPi, float nsigmaTOFPi, float nsigmaTPCK, float nsigmaTOFK)
  {
    enum ParticleNumber { partOne = 1,
                          partTwo = 2 };
    if (particle_number == partOne) {
      switch (trackonefilter.confPDGCodePartOne) {
        case 2212:  // Proton+
        case -2212: // Proton-
          return isNSigma(mom, nsigmaTPCPr, nsigmaTOFPr);
          break;
        case 211:  // Pion+
        case -211: // Pion-
        case 111:  // Pion 0
          return isNSigma(mom, nsigmaTPCPi, nsigmaTOFPi);
          break;
        case 321:  // Kaon+
        case -321: // Kaon-
          return isKaonNsigma(mom, nsigmaTPCK, nsigmaTOFK);
          break;
        case 130: // Kaon 0 LONG
        case 310: // Kaon 0 SHORT
          return isNSigma(mom, nsigmaTPCK, nsigmaTOFK);
          break;
        default:
          return false;
      }
      return false;
    } else if (particle_number == partTwo) {
      switch (tracktwofilter.confPDGCodePartTwo) {
        case 2212:  // Proton+
        case -2212: // Proton-
          return isNSigma(mom, nsigmaTPCPr, nsigmaTOFPr);
          break;
        case 211:  // Pion+
        case -211: // Pion-
        case 111:  // Pion 0
          return isNSigma(mom, nsigmaTPCPi, nsigmaTOFPi);
          break;
        case 321:  // Kaon+
        case -321: // Kaon-
          return isKaonNsigma(mom, nsigmaTPCK, nsigmaTOFK);
          break;
        case 130: // Kaon 0 LONG
        case 310: // Kaon 0 SHORT
          return isNSigma(mom, nsigmaTPCK, nsigmaTOFK);
          break;
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

    trackHistoPartOne.init(&qaRegistry, confTempFitVarpTBins, confTempFitVarBins, twotracksconfigs.confIsMC, trackonefilter.confPDGCodePartOne, true);

    trackHistoPartTwo.init(&qaRegistry, confTempFitVarpTBins, confTempFitVarBins, twotracksconfigs.confIsMC, tracktwofilter.confPDGCodePartTwo, true);

    if (confFillDebug) {
      sphericityRegistry.add("sphericity", ";Sphericity;Entries", kTH1F, {{150, 0.0, 3, "Sphericity"}});
      mixQaRegistry.add("MixingQA/hSECollisionBins", ";bin;Entries", kTH1F, {{120, -0.5, 119.5}});
      mixQaRegistry.add("MixingQA/hMECollisionBins", ";bin;Entries", kTH1F, {{120, -0.5, 119.5}});
    }

    mass1 = pdg->Mass(trackonefilter.confPDGCodePartOne);
    mass2 = pdg->Mass(tracktwofilter.confPDGCodePartTwo);

    if (processPair.cfgProcessPM) {
      sameEventCont.init(&resultRegistryPM, confkstarBins, confMultBins, confkTBins, confmTBins, confmultBins3D, confmTBins3D, twotracksconfigs.confIsMC, twotracksconfigs.confUse3D);
      mixedEventCont.init(&resultRegistryPM, confkstarBins, confMultBins, confkTBins, confmTBins, confmultBins3D, confmTBins3D, twotracksconfigs.confIsMC, twotracksconfigs.confUse3D);

      sameEventCont.setPDGCodes(trackonefilter.confPDGCodePartOne, tracktwofilter.confPDGCodePartTwo);
      mixedEventCont.setPDGCodes(trackonefilter.confPDGCodePartOne, tracktwofilter.confPDGCodePartTwo);

      if (cfgProcessMultBins) {
        sameEventMultCont.init(&sameMultRegistryPM, confkstarBins, confMultKstarBins, confKtKstarBins, cfgProcessKtBins, cfgProcessKtMt3DCF);
        mixedEventMultCont.init(&mixedMultRegistryPM, confkstarBins, confMultKstarBins, confKtKstarBins, cfgProcessKtBins, cfgProcessKtMt3DCF);
      }
      if (doprocessFractions) {
        mixedMultRegistryPM.add("MCreco/motherParticle", "pair fractions; part1 mother PDG; part2 mother PDG", {HistType::kTH2F, {{8001, -4000, 4000}, {8001, -4000, 4000}}});
        mixedMultRegistryPM.add("MCreco/motherParticlePDGCheck", "pair fractions;part1 mother PDG;part2 mother PDG", {HistType::kTH2F, {{8001, -4000, 4000}, {8001, -4000, 4000}}});
      }
      if (doprocessFractionsMCTruth) {
        mixedMultRegistryPM.add("MCtruth/motherParticle", "pair fractions; part1 mother PDG; part2 mother PDG", {HistType::kTH2F, {{8001, -4000, 4000}, {8001, -4000, 4000}}});
      }
    }

    if (processPair.cfgProcessPP) {
      if (doprocessFractions) {
        mixedMultRegistryPP.add("MCreco/motherParticle", "pair fractions; part1 mother PDG; part2 mother PDG", {HistType::kTH2F, {{8001, -4000, 4000}, {8001, -4000, 4000}}});
        mixedMultRegistryPP.add("MCreco/motherParticlePDGCheck", "pair fractions;part1 mother PDG;part2 mother PDG", {HistType::kTH2F, {{8001, -4000, 4000}, {8001, -4000, 4000}}});
      }
      if (doprocessFractionsMCTruth) {
        mixedMultRegistryPP.add("MCtruth/motherParticle", "pair fractions; part1 mother PDG; part2 mother PDG", {HistType::kTH2F, {{8001, -4000, 4000}, {8001, -4000, 4000}}});
      }
      sameEventContPP.init(&resultRegistryPP, confkstarBins, confMultBins, confkTBins, confmTBins, confmultBins3D, confmTBins3D, twotracksconfigs.confIsMC, twotracksconfigs.confUse3D);
      mixedEventContPP.init(&resultRegistryPP, confkstarBins, confMultBins, confkTBins, confmTBins, confmultBins3D, confmTBins3D, twotracksconfigs.confIsMC, twotracksconfigs.confUse3D);
      sameEventContPP.setPDGCodes(trackonefilter.confPDGCodePartOne, tracktwofilter.confPDGCodePartTwo);
      mixedEventContPP.setPDGCodes(trackonefilter.confPDGCodePartOne, tracktwofilter.confPDGCodePartTwo);

      if (cfgProcessMultBins) {
        sameEventMultContPP.init(&sameMultRegistryPP, confkstarBins, confMultKstarBins, confKtKstarBins, cfgProcessKtBins, cfgProcessKtMt3DCF);
        mixedEventMultContPP.init(&mixedMultRegistryPP, confkstarBins, confMultKstarBins, confKtKstarBins, cfgProcessKtBins, cfgProcessKtMt3DCF);
      }
    }

    if (processPair.cfgProcessMM) {
      if (doprocessFractions) {
        mixedMultRegistryMM.add("MCreco/motherParticle", "pair fractions; part1 mother PDG; part2 mother PDG", {HistType::kTH2F, {{8001, -4000, 4000}, {8001, -4000, 4000}}});
        mixedMultRegistryMM.add("MCreco/motherParticlePDGCheck", "pair fractions;part1 mother PDG;part2 mother PDG", {HistType::kTH2F, {{8001, -4000, 4000}, {8001, -4000, 4000}}});
      }
      if (doprocessFractionsMCTruth) {
        mixedMultRegistryMM.add("MCtruth/motherParticle", "pair fractions; part1 mother PDG; part2 mother PDG", {HistType::kTH2F, {{8001, -4000, 4000}, {8001, -4000, 4000}}});
      }
      sameEventContMM.init(&resultRegistryMM, confkstarBins, confMultBins, confkTBins, confmTBins, confmultBins3D, confmTBins3D, twotracksconfigs.confIsMC, twotracksconfigs.confUse3D);
      mixedEventContMM.init(&resultRegistryMM, confkstarBins, confMultBins, confkTBins, confmTBins, confmultBins3D, confmTBins3D, twotracksconfigs.confIsMC, twotracksconfigs.confUse3D);
      sameEventContMM.setPDGCodes(trackonefilter.confPDGCodePartOne, tracktwofilter.confPDGCodePartTwo);
      mixedEventContMM.setPDGCodes(trackonefilter.confPDGCodePartOne, tracktwofilter.confPDGCodePartTwo);

      if (cfgProcessMultBins) {
        sameEventMultContMM.init(&sameMultRegistryMM, confkstarBins, confMultKstarBins, confKtKstarBins, cfgProcessKtBins, cfgProcessKtMt3DCF);
        mixedEventMultContMM.init(&mixedMultRegistryMM, confkstarBins, confMultKstarBins, confKtKstarBins, cfgProcessKtBins, cfgProcessKtMt3DCF);
      }
    }

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
    if (confFillDebug) {
      mixQaRegistry.fill(HIST("MixingQA/hSECollisionBins"), colBinning.getBin({col.posZ(), col.multV0M()}));
    }
    eventHisto.fillQA(col);
  }

  /// This function processes the same event and takes care of all the histogramming
  /// \todo the trivial loops over the tracks should be factored out since they will be common to all combinations of T-T, T-V0, V0-V0, ...
  /// \tparam PartitionType
  /// \tparam PartType
  /// \tparam isMC: enables Monte Carlo truth specific histograms
  /// \param groupPartsOne partition for the first particle passed by the process function
  /// \param groupPartsTwo partition for the second particle passed by the process function
  /// \param parts femtoUniverseParticles table (in case of Monte Carlo joined with FemtoUniverseMCLabels)
  /// \param magFieldTesla magnetic field of the collision
  /// \param multCol multiplicity of the collision
  /// \param pairType describes charge of correlation pair (plus-minus (1), plus-plus (2), minus-minus (3))
  /// \param fillQA enables filling of QA histograms
  template <bool isMC, typename PartitionType, typename PartType>
  void doSameEvent(PartitionType groupPartsOne, PartitionType groupPartsTwo, PartType parts, float magFieldTesla, int multCol, int pairType, bool fillQA)
  {

    enum PairTypes { type1 = 1,
                     type2 = 2,
                     type3 = 3 };
    /// Histogramming same event
    if ((pairType == type1 || pairType == type2) && fillQA) {
      for (const auto& part : groupPartsOne) {
        if (!isParticleNSigma((int8_t)1, part.p(), trackCuts.getNsigmaTPC(part, o2::track::PID::Proton), trackCuts.getNsigmaTOF(part, o2::track::PID::Proton), trackCuts.getNsigmaTPC(part, o2::track::PID::Pion), trackCuts.getNsigmaTOF(part, o2::track::PID::Pion), trackCuts.getNsigmaTPC(part, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(part, o2::track::PID::Kaon))) {
          continue;
        }
        trackHistoPartOne.fillQA<isMC, true>(part);
      }
    }

    if ((pairType == type1 || pairType == type3) && fillQA) {
      for (const auto& part : groupPartsTwo) {
        if (!isParticleNSigma((int8_t)2, part.p(), trackCuts.getNsigmaTPC(part, o2::track::PID::Proton), trackCuts.getNsigmaTOF(part, o2::track::PID::Proton), trackCuts.getNsigmaTPC(part, o2::track::PID::Pion), trackCuts.getNsigmaTOF(part, o2::track::PID::Pion), trackCuts.getNsigmaTPC(part, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(part, o2::track::PID::Kaon))) {
          continue;
        }
        trackHistoPartTwo.fillQA<isMC, true>(part);
      }
    }

    if (pairType == type1) {

      /// Now build the combinations for non-identical particle pairs
      for (const auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(groupPartsOne, groupPartsTwo))) {

        if (!isParticleNSigma((int8_t)1, p1.p(), trackCuts.getNsigmaTPC(p1, o2::track::PID::Proton), trackCuts.getNsigmaTOF(p1, o2::track::PID::Proton), trackCuts.getNsigmaTPC(p1, o2::track::PID::Pion), trackCuts.getNsigmaTOF(p1, o2::track::PID::Pion), trackCuts.getNsigmaTPC(p1, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(p1, o2::track::PID::Kaon))) {
          continue;
        }

        if (!isParticleNSigma((int8_t)2, p2.p(), trackCuts.getNsigmaTPC(p2, o2::track::PID::Proton), trackCuts.getNsigmaTOF(p2, o2::track::PID::Proton), trackCuts.getNsigmaTPC(p2, o2::track::PID::Pion), trackCuts.getNsigmaTOF(p2, o2::track::PID::Pion), trackCuts.getNsigmaTPC(p2, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(p2, o2::track::PID::Kaon))) {
          continue;
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

        float kstar = FemtoUniverseMath::getkstar(p1, mass1, p2, mass2);
        float kT = FemtoUniverseMath::getkT(p1, mass1, p2, mass2);

        sameEventCont.setPair<isMC>(p1, p2, multCol, twotracksconfigs.confUse3D);
        if (cfgProcessMultBins)
          sameEventMultCont.fill<float>(kstar, multCol, kT);
      }
    } else {
      /// Now build the combinations for identical particles pairs
      for (const auto& [p1, p2] : combinations(CombinationsStrictlyUpperIndexPolicy(groupPartsOne, groupPartsOne))) {

        if (!isParticleNSigma((int8_t)2, p1.p(), trackCuts.getNsigmaTPC(p1, o2::track::PID::Proton), trackCuts.getNsigmaTOF(p1, o2::track::PID::Proton), trackCuts.getNsigmaTPC(p1, o2::track::PID::Pion), trackCuts.getNsigmaTOF(p1, o2::track::PID::Pion), trackCuts.getNsigmaTPC(p1, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(p1, o2::track::PID::Kaon))) {
          continue;
        }

        if (!isParticleNSigma((int8_t)2, p2.p(), trackCuts.getNsigmaTPC(p2, o2::track::PID::Proton), trackCuts.getNsigmaTOF(p2, o2::track::PID::Proton), trackCuts.getNsigmaTPC(p2, o2::track::PID::Pion), trackCuts.getNsigmaTOF(p2, o2::track::PID::Pion), trackCuts.getNsigmaTPC(p2, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(p2, o2::track::PID::Kaon))) {
          continue;
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

        switch (pairType) {
          case 2: {
            if (isPairIdentical == true) {
              float kstar = FemtoUniverseMath::getkstar(p1, mass1, p2, mass1);
              float kT = FemtoUniverseMath::getkT(p1, mass1, p2, mass1);

              sameEventContPP.setPair<isMC>(p1, p2, multCol, twotracksconfigs.confUse3D);
              if (cfgProcessMultBins)
                sameEventMultContPP.fill<float>(kstar, multCol, kT);
            } else {
              float kstar = FemtoUniverseMath::getkstar(p1, mass1, p2, mass2);
              float kT = FemtoUniverseMath::getkT(p1, mass1, p2, mass2);

              sameEventContPP.setPair<isMC>(p1, p2, multCol, twotracksconfigs.confUse3D);
              if (cfgProcessMultBins)
                sameEventMultContPP.fill<float>(kstar, multCol, kT);
            }

            break;
          }

          case 3: {
            if (isPairIdentical == true) {
              float kstar = FemtoUniverseMath::getkstar(p1, mass2, p2, mass2);
              float kT = FemtoUniverseMath::getkT(p1, mass2, p2, mass2);

              sameEventContMM.setPair<isMC>(p1, p2, multCol, twotracksconfigs.confUse3D);
              if (cfgProcessMultBins)
                sameEventMultContMM.fill<float>(kstar, multCol, kT);
            } else {
              float kstar = FemtoUniverseMath::getkstar(p1, mass1, p2, mass2);
              float kT = FemtoUniverseMath::getkT(p1, mass1, p2, mass2);

              sameEventContMM.setPair<isMC>(p1, p2, multCol, twotracksconfigs.confUse3D);
              if (cfgProcessMultBins)
                sameEventMultContMM.fill<float>(kstar, multCol, kT);
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
  void processSameEvent(soa::Filtered<o2::aod::FdCollisions>::iterator const& col,
                        FilteredFemtoFullParticles const& parts)
  {
    fillCollision(col);
    if (confFillDebug) {
      sphericityRegistry.fill(HIST("sphericity"), col.sphericity());
    }

    auto thegroupPartsOne = partsOne->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto thegroupPartsTwo = partsTwo->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);

    bool fillQA = true;

    if (processPair.cfgProcessPM) {
      doSameEvent<false>(thegroupPartsOne, thegroupPartsTwo, parts, col.magField(), col.multV0M(), 1, fillQA);
      fillQA = false;
    }
    if (processPair.cfgProcessPP)
      doSameEvent<false>(thegroupPartsOne, thegroupPartsOne, parts, col.magField(), col.multV0M(), 2, fillQA);
    if (processPair.cfgProcessMM)
      doSameEvent<false>(thegroupPartsTwo, thegroupPartsTwo, parts, col.magField(), col.multV0M(), 3, fillQA);
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackTrackMultKtExtended, processSameEvent, "Enable processing same event", true);

  /// process function for to call doSameEvent with Monte Carlo
  /// \param col subscribe to the collision table (Monte Carlo Reconstructed reconstructed)
  /// \param parts subscribe to joined table FemtoUniverseParticles and FemtoUniverseMCLabels to access Monte Carlo truth
  /// \param FemtoUniverseMCParticles subscribe to the Monte Carlo truth table
  void processSameEventMC(o2::aod::FdCollision const& col,
                          FemtoRecoParticles const& parts,
                          o2::aod::FdMCParticles const&)
  {
    fillCollision(col);

    auto thegroupPartsOne = partsOneMC->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto thegroupPartsTwo = partsTwoMC->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);

    bool fillQA = true;
    if (processPair.cfgProcessPM) {
      doSameEvent<true>(thegroupPartsOne, thegroupPartsTwo, parts, col.magField(), col.multV0M(), 1, fillQA);
      fillQA = false;
    }
    if (processPair.cfgProcessPP)
      doSameEvent<true>(thegroupPartsOne, thegroupPartsOne, parts, col.magField(), col.multV0M(), 2, fillQA);
    if (processPair.cfgProcessMM)
      doSameEvent<true>(thegroupPartsTwo, thegroupPartsTwo, parts, col.magField(), col.multV0M(), 3, fillQA);
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackTrackMultKtExtended, processSameEventMC, "Enable processing same event for Monte Carlo", false);

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
  /// \param pairType describes charge of correlation pair (plus-minus (1), plus-plus (2), minus-minus (3))
  template <bool isMC, typename PartitionType, typename PartType>
  void doMixedEvent(PartitionType groupPartsOne, PartitionType groupPartsTwo, PartType parts, float magFieldTesla, int multCol, int pairType)
  {

    for (const auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(groupPartsOne, groupPartsTwo))) {

      if (!isParticleNSigma((int8_t)2, p1.p(), trackCuts.getNsigmaTPC(p1, o2::track::PID::Proton), trackCuts.getNsigmaTOF(p1, o2::track::PID::Proton), trackCuts.getNsigmaTPC(p1, o2::track::PID::Pion), trackCuts.getNsigmaTOF(p1, o2::track::PID::Pion), trackCuts.getNsigmaTPC(p1, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(p1, o2::track::PID::Kaon))) {
        continue;
      }

      if (!isParticleNSigma((int8_t)2, p2.p(), trackCuts.getNsigmaTPC(p2, o2::track::PID::Proton), trackCuts.getNsigmaTOF(p2, o2::track::PID::Proton), trackCuts.getNsigmaTPC(p2, o2::track::PID::Pion), trackCuts.getNsigmaTOF(p2, o2::track::PID::Pion), trackCuts.getNsigmaTPC(p2, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(p2, o2::track::PID::Kaon))) {
        continue;
      }

      if (confIsCPR.value) {
        if (pairCloseRejection.isClosePair(p1, p2, parts, magFieldTesla, femto_universe_container::EventType::mixed)) {
          continue;
        }
      }

      switch (pairType) {
        case 1: {
          float kstar = FemtoUniverseMath::getkstar(p1, mass1, p2, mass2);
          float kT = FemtoUniverseMath::getkT(p1, mass1, p2, mass2);

          mixedEventCont.setPair<isMC>(p1, p2, multCol, twotracksconfigs.confUse3D);
          if (cfgProcessMultBins)
            mixedEventMultCont.fill<float>(kstar, multCol, kT);

          break;
        }
        case 2: {
          if (isPairIdentical == true) {
            float kstar = FemtoUniverseMath::getkstar(p1, mass1, p2, mass1);
            float kT = FemtoUniverseMath::getkT(p1, mass1, p2, mass1);

            mixedEventContPP.setPair<isMC>(p1, p2, multCol, twotracksconfigs.confUse3D);
            if (cfgProcessMultBins)
              mixedEventMultContPP.fill<float>(kstar, multCol, kT);
          } else {
            float kstar = FemtoUniverseMath::getkstar(p1, mass1, p2, mass2);
            float kT = FemtoUniverseMath::getkT(p1, mass1, p2, mass2);

            mixedEventContPP.setPair<isMC>(p1, p2, multCol, twotracksconfigs.confUse3D);
            if (cfgProcessMultBins)
              mixedEventMultContPP.fill<float>(kstar, multCol, kT);
          }

          break;
        }

        case 3: {
          if (isPairIdentical == true) {
            float kstar = FemtoUniverseMath::getkstar(p1, mass2, p2, mass2);
            float kT = FemtoUniverseMath::getkT(p1, mass2, p2, mass2);

            mixedEventContMM.setPair<isMC>(p1, p2, multCol, twotracksconfigs.confUse3D);
            if (cfgProcessMultBins)
              mixedEventMultContMM.fill<float>(kstar, multCol, kT);
          } else {
            float kstar = FemtoUniverseMath::getkstar(p1, mass1, p2, mass2);
            float kT = FemtoUniverseMath::getkT(p1, mass1, p2, mass2);

            mixedEventContMM.setPair<isMC>(p1, p2, multCol, twotracksconfigs.confUse3D);
            if (cfgProcessMultBins)
              mixedEventMultContMM.fill<float>(kstar, multCol, kT);
          }

          break;
        }
        default:
          break;
      }
    }
  }

  /// process function for to call doMixedEvent with Data
  /// \param cols subscribe to the collisions table (Data)
  /// \param parts subscribe to the femtoUniverseParticleTable
  void processMixedEvent(soa::Filtered<o2::aod::FdCollisions> const& cols,
                         FilteredFemtoFullParticles const& parts)
  {
    for (const auto& [collision1, collision2] : soa::selfCombinations(colBinning, 5, -1, cols, cols)) {

      const int multiplicityCol = collision1.multV0M();
      if (confFillDebug) {
        mixQaRegistry.fill(HIST("MixingQA/hMECollisionBins"), colBinning.getBin({collision1.posZ(), multiplicityCol}));
      }

      const auto& magFieldTesla1 = collision1.magField();
      const auto& magFieldTesla2 = collision2.magField();

      if (magFieldTesla1 != magFieldTesla2) {
        continue;
      }

      if (processPair.cfgProcessPM) {
        auto groupPartsOne = partsOne->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);
        auto groupPartsTwo = partsTwo->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision2.globalIndex(), cache);
        doMixedEvent<false>(groupPartsOne, groupPartsTwo, parts, magFieldTesla1, multiplicityCol, 1);
      }
      if (processPair.cfgProcessPP) {
        auto groupPartsOne = partsOne->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);
        auto groupPartsTwo = partsOne->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision2.globalIndex(), cache);
        doMixedEvent<false>(groupPartsOne, groupPartsTwo, parts, magFieldTesla1, multiplicityCol, 2);
      }
      if (processPair.cfgProcessMM) {
        auto groupPartsOne = partsTwo->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);
        auto groupPartsTwo = partsTwo->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision2.globalIndex(), cache);
        doMixedEvent<false>(groupPartsOne, groupPartsTwo, parts, magFieldTesla1, multiplicityCol, 3);
      }
    }
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackTrackMultKtExtended, processMixedEvent, "Enable processing mixed events", true);

  /// brief process function for to call doMixedEvent with Monte Carlo
  /// \param cols subscribe to the collisions table (Monte Carlo Reconstructed reconstructed)
  /// \param parts subscribe to joined table FemtoUniverseParticles and FemtoUniverseMCLables to access Monte Carlo truth
  /// \param FemtoUniverseMCParticles subscribe to the Monte Carlo truth table
  void processMixedEventMC(o2::aod::FdCollisions const& cols,
                           FemtoRecoParticles const& parts,
                           o2::aod::FdMCParticles const&)
  {
    for (const auto& [collision1, collision2] : soa::selfCombinations(colBinning, 5, -1, cols, cols)) {

      const int multiplicityCol = collision1.multV0M();
      if (confFillDebug) {
        mixQaRegistry.fill(HIST("MixingQA/hMECollisionBins"), colBinning.getBin({collision1.posZ(), multiplicityCol}));
      }
      const auto& magFieldTesla1 = collision1.magField();
      const auto& magFieldTesla2 = collision2.magField();

      if (magFieldTesla1 != magFieldTesla2) {
        continue;
      }
      /// \todo before mixing we should check whether both collisions contain a pair of particles!
      // if (partsOne.size() == 0 || nPart2Evt1 == 0 || nPart1Evt2 == 0 || partsTwo.size() == 0 ) continue;

      if (processPair.cfgProcessPM) {
        auto groupPartsOne = partsOneMC->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);
        auto groupPartsTwo = partsTwoMC->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision2.globalIndex(), cache);
        doMixedEvent<true>(groupPartsOne, groupPartsTwo, parts, magFieldTesla1, multiplicityCol, 1);
      }
      if (processPair.cfgProcessPP) {
        auto groupPartsOne = partsOneMC->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);
        auto groupPartsTwo = partsOneMC->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision2.globalIndex(), cache);
        doMixedEvent<true>(groupPartsOne, groupPartsTwo, parts, magFieldTesla1, multiplicityCol, 2);
      }
      if (processPair.cfgProcessMM) {
        auto groupPartsOne = partsTwoMC->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);
        auto groupPartsTwo = partsTwoMC->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision2.globalIndex(), cache);
        doMixedEvent<true>(groupPartsOne, groupPartsTwo, parts, magFieldTesla1, multiplicityCol, 3);
      }
    }
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackTrackMultKtExtended, processMixedEventMC, "Enable processing mixed events MC reco", false);

  /// brief process function for to call process fractions with Monte Carlo reco
  /// \param cols subscribe to the collisions table (Monte Carlo Reconstructed reconstructed)
  /// \param parts subscribe to joined table FemtoUniverseParticles and FemtoUniverseMCLables to access Monte Carlo truth
  /// \param FemtoUniverseMCParticles subscribe to the Monte Carlo truth table
  void processFractions(o2::aod::FdCollisions const& cols,
                        FemtoRecoParticles const& parts,
                        o2::aod::FdMCParticles const& mparts)
  {

    auto doFractions = [&](auto& p1, auto& p2, auto& magFieldTesla, int partType) -> void {
      if (!isParticleNSigma((int8_t)2, p1.p(), trackCuts.getNsigmaTPC(p1, o2::track::PID::Proton), trackCuts.getNsigmaTOF(p1, o2::track::PID::Proton), trackCuts.getNsigmaTPC(p1, o2::track::PID::Pion), trackCuts.getNsigmaTOF(p1, o2::track::PID::Pion), trackCuts.getNsigmaTPC(p1, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(p1, o2::track::PID::Kaon))) {
        return;
      }

      if (!isParticleNSigma((int8_t)2, p2.p(), trackCuts.getNsigmaTPC(p2, o2::track::PID::Proton), trackCuts.getNsigmaTOF(p2, o2::track::PID::Proton), trackCuts.getNsigmaTPC(p2, o2::track::PID::Pion), trackCuts.getNsigmaTOF(p2, o2::track::PID::Pion), trackCuts.getNsigmaTPC(p2, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(p2, o2::track::PID::Kaon))) {
        return;
      }

      if (confIsCPR.value) {
        if (pairCloseRejection.isClosePair(p1, p2, parts, magFieldTesla, femto_universe_container::EventType::mixed)) {
          return;
        }
      }
      switch (partType) {
        case 1:
          mixedMultRegistryPM.fill(HIST("MCreco/motherParticle"), p1.motherPDG(), p2.motherPDG());
          break;
        case 2:
          mixedMultRegistryPP.fill(HIST("MCreco/motherParticle"), p1.motherPDG(), p2.motherPDG());
          break;
        case 3:
          mixedMultRegistryMM.fill(HIST("MCreco/motherParticle"), p1.motherPDG(), p2.motherPDG());
          break;
      }
      auto mcPartId1 = p1.fdMCParticleId();
      if (mcPartId1 == -1)
        return;
      auto mcPartId2 = p2.fdMCParticleId();
      if (mcPartId2 == -1)
        return;
      const auto& mcParticle1 = mparts.iteratorAt(mcPartId1);
      const auto& mcParticle2 = mparts.iteratorAt(mcPartId2);
      switch (partType) {
        case 1:
          if ((trackonefilter.confPDGCodePartOne != mcParticle1.pdgMCTruth()) || (tracktwofilter.confPDGCodePartTwo != mcParticle2.pdgMCTruth()))
            return;
          mixedMultRegistryPM.fill(HIST("MCreco/motherParticlePDGCheck"), p1.motherPDG(), p2.motherPDG());
          break;
        case 2:
          if ((trackonefilter.confPDGCodePartOne != mcParticle1.pdgMCTruth()) || (trackonefilter.confPDGCodePartOne != mcParticle2.pdgMCTruth()))
            return;
          mixedMultRegistryPP.fill(HIST("MCreco/motherParticlePDGCheck"), p1.motherPDG(), p2.motherPDG());
          break;
        case 3:
          if ((tracktwofilter.confPDGCodePartTwo != mcParticle1.pdgMCTruth()) || (tracktwofilter.confPDGCodePartTwo != mcParticle2.pdgMCTruth()))
            return;
          mixedMultRegistryMM.fill(HIST("MCreco/motherParticlePDGCheck"), p1.motherPDG(), p2.motherPDG());
          break;
      }
    };

    for (const auto& [collision1, collision2] : soa::selfCombinations(colBinning, 5, -1, cols, cols)) {
      const int multiplicityCol = collision1.multV0M();
      if (confFillDebug) {
        mixQaRegistry.fill(HIST("MixingQA/hMECollisionBins"), colBinning.getBin({collision1.posZ(), multiplicityCol}));
      }
      const auto& magFieldTesla1 = collision1.magField();
      const auto& magFieldTesla2 = collision2.magField();

      if (magFieldTesla1 != magFieldTesla2) {
        continue;
      }

      auto groupPartsOne = partsOneMC->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);
      auto groupPartsTwo = partsTwoMC->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision2.globalIndex(), cache);

      if (processPair.cfgProcessPM) {
        for (const auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(groupPartsOne, groupPartsTwo))) {
          doFractions(p1, p2, magFieldTesla1, 1);
        }
      }
      if (processPair.cfgProcessPP) {
        for (const auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(groupPartsOne, groupPartsOne))) {
          doFractions(p1, p2, magFieldTesla1, 2);
        }
      }
      if (processPair.cfgProcessMM) {
        for (const auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(groupPartsTwo, groupPartsTwo))) {
          doFractions(p1, p2, magFieldTesla1, 3);
        }
      }
    }
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackTrackMultKtExtended, processFractions, "Enable processing fractions for MC reco", false);

  /// This function processes fills histograms for fractions analysis
  /// \todo the trivial loops over the collisions and tracks should be factored out since they will be common to all combinations of T-T, T-V0, V0-V0, ...
  /// \tparam PartitionType
  /// \tparam PartType
  /// \tparam isMC: enables Monte Carlo truth specific histograms
  /// \param groupPartsOne partition for the first particle passed by the process function
  /// \param groupPartsTwo partition for the second particle passed by the process function
  template <bool isMC, bool doFractions, typename PartitionType>
  void doFractionsMCTruth(PartitionType groupPartsOne, PartitionType groupPartsTwo, int ContType)
  {

    for (const auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(groupPartsOne, groupPartsTwo))) {

      int pdgCodePartOne = static_cast<int>(p1.pidCut());
      const auto& pdgParticleOne = pdg->GetParticle(pdgCodePartOne);
      int pdgCodePartTwo = static_cast<int>(p2.pidCut());
      const auto& pdgParticleTwo = pdg->GetParticle(pdgCodePartTwo);
      switch (ContType) {
        case 1: {
          if (pdgParticleOne && pdgParticleTwo && (pdgCodePartOne == trackonefilter.confPDGCodePartOne) && (pdgCodePartTwo == tracktwofilter.confPDGCodePartTwo)) {
            continue;
          }
          mixedMultRegistryPM.fill(HIST("MCtruth/motherParticle"), p1.tempFitVar(), p2.tempFitVar());
          break;
        }
        case 2: {
          if (pdgParticleOne && pdgParticleTwo && (pdgCodePartOne == trackonefilter.confPDGCodePartOne) && (pdgCodePartTwo == trackonefilter.confPDGCodePartOne)) {
            continue;
          }
          mixedMultRegistryPP.fill(HIST("MCtruth/motherParticle"), p1.tempFitVar(), p2.tempFitVar());
          break;
        }
        case 3: {
          if (pdgParticleOne && pdgParticleTwo && (pdgCodePartOne == tracktwofilter.confPDGCodePartTwo) && (pdgCodePartTwo == tracktwofilter.confPDGCodePartTwo)) {
            continue;
          }
          mixedMultRegistryMM.fill(HIST("MCtruth/motherParticle"), p1.tempFitVar(), p2.tempFitVar());
        }
        default:
          break;
      }
    }
  }

  /// brief process function for to call doFractionsMCTruth with Monte Carlo
  /// \param cols subscribe to the collisions table (Monte Carlo Reconstructed reconstructed)
  /// \param FemtoTruthParticles subscribe to the Monte Carlo truth particle table
  void processFractionsMCTruth(o2::aod::FdCollisions const& cols,
                               FemtoTruthParticles const&)
  {
    for (const auto& [collision1, collision2] : soa::selfCombinations(colBinning, 5, -1, cols, cols)) {

      const int multiplicityCol = collision1.multV0M();
      if (confFillDebug) {
        mixQaRegistry.fill(HIST("MixingQA/hMECollisionBins"), colBinning.getBin({collision1.posZ(), multiplicityCol}));
      }
      const auto& magFieldTesla1 = collision1.magField();
      const auto& magFieldTesla2 = collision2.magField();

      if (magFieldTesla1 != magFieldTesla2) {
        continue;
      }
      /// \todo before mixing we should check whether both collisions contain a pair of particles!
      // if (partsOne.size() == 0 || nPart2Evt1 == 0 || nPart1Evt2 == 0 || partsTwo.size() == 0 ) continue;
      auto groupPartsOne = partsOneMCTruth->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);
      auto groupPartsTwo = partsTwoMCTruth->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision2.globalIndex(), cache);
      if (processPair.cfgProcessPM) {
        doFractionsMCTruth<true, true>(groupPartsOne, groupPartsTwo, 1);
      }
      if (processPair.cfgProcessPP) {
        doFractionsMCTruth<true, true>(groupPartsOne, groupPartsOne, 2);
      }
      if (processPair.cfgProcessMM) {
        doFractionsMCTruth<true, true>(groupPartsTwo, groupPartsTwo, 3);
      }
    }
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackTrackMultKtExtended, processFractionsMCTruth, "Enable processing residual correlations fractions for MC truth", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<FemtoUniversePairTaskTrackTrackMultKtExtended>(cfgc),
  };
  return workflow;
}
