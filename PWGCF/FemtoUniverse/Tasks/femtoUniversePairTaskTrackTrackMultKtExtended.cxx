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

/// \file femtoUniversePairTaskTrackTrackMultKtExtended.cxx
/// \brief Tasks that reads the track tables used for the pairing and builds pairs of two tracks
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de
/// \author Georgios Mantzaridis, TU München, georgios.mantzaridis@tum.de
/// \author Anton Riedel, TU München, anton.riedel@tum.de
/// \author Zuzanna Chochulska, WUT Warsaw & CTU Prague, zchochul@cern.ch
/// \author Alicja Płachta, WUT Warsaw, alicja.plachta.stud@pw.edu.pl

#include <vector>
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/StepTHn.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "TDatabasePDG.h"
#include "ReconstructionDataFormats/PID.h"
#include "Common/DataModel/PIDResponse.h"

#include "PWGCF/FemtoUniverse/DataModel/FemtoDerived.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseParticleHisto.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseEventHisto.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniversePairCleaner.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseFemtoContainer.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseDetaDphiStar.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUtils.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseMath.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseTrackSelection.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniversePairWithCentMultKt.h"

using namespace o2;
using namespace o2::analysis::femtoUniverse;
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

struct femtoUniversePairTaskTrackTrackMultKtExtended {

  Service<o2::framework::O2DatabasePDG> pdg;

  /// Particle selection part

  /// Table for both particles
  struct : o2::framework::ConfigurableGroup {
    Configurable<bool> IsKaonNsigma{"IsKaonNsigma", false, "Enable a strict cut selection for K+ and K-"};
    Configurable<float> ConfNsigmaCombined{"ConfNsigmaCombined", 3.0f, "TPC and TOF Pion Sigma (combined) for momentum > ConfTOFpMin"};
    Configurable<float> ConfNsigmaTPC{"ConfNsigmaTPC", 3.0f, "TPC Pion Sigma for momentum < ConfTOFpMin"};
    Configurable<float> ConfTOFpMin{"ConfTOFpMin", 0.5f, "Min. momentum for which TOF is required for PID."};
    Configurable<float> ConfEtaMax{"ConfEtaMax", 0.8f, "Higher limit for |Eta| (the same for both particles)"};

    Configurable<LabeledArray<float>> ConfCutTable{"ConfCutTable", {cutsTable[0], nPart, nCuts, partNames, cutNames}, "Particle selections"};
    Configurable<int> ConfNspecies{"ConfNspecies", 2, "Number of particle spieces with PID info"};
    Configurable<bool> ConfIsMC{"ConfIsMC", false, "Enable additional Histogramms in the case of a MonteCarlo Run"};
    Configurable<std::vector<float>> ConfTrkPIDnSigmaMax{"ConfTrkPIDnSigmaMax", std::vector<float>{4.f, 3.f, 2.f}, "This configurable needs to be the same as the one used in the producer task"};
    Configurable<bool> ConfUse3D{"ConfUse3D", false, "Enable three dimensional histogramms (to be used only for analysis with high statistics): k* vs mT vs multiplicity"};
  } twotracksconfigs;

  /// Table for separate deuteron configurables
  struct : o2::framework::ConfigurableGroup {
    Configurable<float> ConfNsigmaTPCDe{"ConfNsigmaTPCDe", 2.0f, "TPC Deuteron Sigma for momentum < ConfTOFpMinDe"};
    Configurable<float> ConfNsigmaTOFDe{"ConfNsigmaTOFDe", 2.0f, "TOF Deuteron Sigma"};
    Configurable<float> ConfTOFpMinDe{"ConfTOFpMinDe", 0.5f, "Min. momentum for deuterons for which TOF is required for PID"};
    Configurable<float> ConfPLowDe{"ConfPLowDe", 0.8f, "Lower limit for momentum for deuterons"};
    Configurable<float> ConfPHighDe{"ConfPHighDe", 1.8f, "Higher limit for momentum for deuterons"};
  } deuteronconfigs;

  /// Table for linear cut for TPC Deuteron Sigma
  struct : o2::framework::ConfigurableGroup {
    Configurable<bool> ConfIsLine{"ConfIsLine", false, "Enable a separation line for clearer TPC Deuteron Sigma"};
    Configurable<float> pLow{"pLow", 0.0f, "Lower limit of momentum for linear cut of TPC Deuteron Sigma"};
    Configurable<float> pHigh{"pHigh", 1.4f, "Higher limit of momentum for linear cut of TPC Deuteron Sigma"};
    Configurable<float> a{"a", -167.0f, "Parameter 'a' of a linear function 'y = a * x + b'"};
    Configurable<float> b{"b", 300.0f, "Parameter 'b' of a linear function 'y = a * x + b'"};
  } lincut;

  /// Table for polynomial 3 cut for TPC Deuteron Sigma
  struct : o2::framework::ConfigurableGroup {
    Configurable<bool> ConfIsPol{"ConfIsPol", false, "Enable a separation polynomial 3 curve for clearer TPC Deuteron Sigma"};
    Configurable<float> A{"A", -52.2f, "Parameter 'A' of a polynomial function 'y = A * x^3 + B * x^2 + C * x + D'"};
    Configurable<float> B{"B", 357.7f, "Parameter 'B' of a polynomial function 'y = A * x^3 + B * x^2 + C * x + D'"};
    Configurable<float> C{"C", -834.7f, "Parameter 'C' of a polynomial function 'y = A * x^3 + B * x^2 + C * x + D'"};
    Configurable<float> D{"D", 705.8f, "Parameter 'D' of a polynomial function 'y = A * x^3 + B * x^2 + C * x + D'"};
  } polcut;

  using FemtoFullParticles = soa::Join<aod::FDParticles, aod::FDExtParticles>;
  // Filters for selecting particles (both p1 and p2)
  Filter trackAdditionalfilter = (nabs(aod::femtouniverseparticle::eta) < twotracksconfigs.ConfEtaMax); // example filtering on configurable
  using FilteredFemtoFullParticles = soa::Filtered<FemtoFullParticles>;
  // using FilteredFemtoFullParticles = FemtoFullParticles; //if no filtering is applied uncomment this option

  SliceCache cache;
  Preslice<FilteredFemtoFullParticles> perCol = aod::femtouniverseparticle::fdCollisionId;

  /// Particle 1
  struct : o2::framework::ConfigurableGroup {
    Configurable<int> ConfPDGCodePartOne{"ConfPDGCodePartOne", 211, "Particle 1 -- PDG code"};
    // Configurable<uint32_t> ConfCutPartOne{"ConfCutPartOne", 5542474, "Particle 1 -- Selection bit from cutCulator"};
    Configurable<int> ConfPIDPartOne{"ConfPIDPartOne", 2, "Particle 1 -- Read from cutCulator"}; // we also need the possibility to specify whether the bit is true/false ->std>>vector<std::pair<int, int>>int>>
    Configurable<float> ConfpLowPart1{"ConfpLowPart1", 0.14, "Lower limit for Pt for the first particle"};
    Configurable<float> ConfPtHighPart1{"ConfPtHighPart1", 1.5, "Higher limit for Pt for the first particle"};
    Configurable<int> ConfChargePart1{"ConfChargePart1", 1, "Particle 1 sign"};
  } trackonefilter;

  /// Partition for particle 1
  Partition<FilteredFemtoFullParticles> partsOne = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kTrack)) && aod::femtouniverseparticle::sign == trackonefilter.ConfChargePart1 && aod::femtouniverseparticle::pt < trackonefilter.ConfPtHighPart1 && aod::femtouniverseparticle::pt > trackonefilter.ConfpLowPart1;

  Partition<soa::Join<FilteredFemtoFullParticles, aod::FDMCLabels>> partsOneMC = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kTrack)) && aod::femtouniverseparticle::sign == trackonefilter.ConfChargePart1 && aod::femtouniverseparticle::pt < trackonefilter.ConfPtHighPart1 && aod::femtouniverseparticle::pt > trackonefilter.ConfpLowPart1;

  /// Histogramming for particle 1
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kTrack, 1> trackHistoPartOne;

  /// Particle 2
  struct : o2::framework::ConfigurableGroup {
    Configurable<int> ConfPDGCodePartTwo{"ConfPDGCodePartTwo", 211, "Particle 2 -- PDG code"};
    // Configurable<uint32_t> ConfCutPartTwo{"ConfCutPartTwo", 5542474, "Particle 2 -- Selection bit"};
    Configurable<int> ConfPIDPartTwo{"ConfPIDPartTwo", 2, "Particle 2 -- Read from cutCulator"}; // we also need the possibility to specify whether the bit is true/false ->std>>vector<std::pair<int, int>>
    Configurable<float> ConfpLowPart2{"ConfpLowPart2", 0.14, "Lower limit for Pt for the second particle"};
    Configurable<float> ConfPtHighPart2{"ConfPtHighPart2", 1.5, "Higher limit for Pt for the second particle"};
    Configurable<int> ConfChargePart2{"ConfChargePart2", -1, "Particle 2 sign"};
  } tracktwofilter;

  /// Partition for particle 2
  Partition<FilteredFemtoFullParticles> partsTwo = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kTrack)) && (aod::femtouniverseparticle::sign == tracktwofilter.ConfChargePart2) && aod::femtouniverseparticle::pt < tracktwofilter.ConfPtHighPart2 && aod::femtouniverseparticle::pt > tracktwofilter.ConfpLowPart2;

  Partition<soa::Join<FilteredFemtoFullParticles, aod::FDMCLabels>> partsTwoMC = aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kTrack) && (aod::femtouniverseparticle::sign == tracktwofilter.ConfChargePart2) && aod::femtouniverseparticle::pt < tracktwofilter.ConfPtHighPart2 && aod::femtouniverseparticle::pt > tracktwofilter.ConfpLowPart2;

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
  Configurable<float> ConfSphericityCutMin{"ConfSphericityCutMin", 0, "Min. sphericity"};
  Configurable<float> ConfSphericityCutMax{"ConfSphericityCutMax", 3, "Max. sphericity"};

  Filter collV0Mfilter = ((o2::aod::femtouniversecollision::multV0M > ConfV0MLow) && (o2::aod::femtouniversecollision::multV0M < ConfV0MHigh));
  Filter colSpherfilter = ((o2::aod::femtouniversecollision::sphericity > ConfSphericityCutMin) && (o2::aod::femtouniversecollision::sphericity < ConfSphericityCutMax));
  // Filter trackAdditionalfilter = (nabs(aod::femtouniverseparticle::eta) < twotracksconfigs.ConfEtaMax); // example filtering on configurable

  /// Particle part
  ConfigurableAxis ConfTempFitVarBins{"ConfDTempFitVarBins", {300, -0.15, 0.15}, "binning of the TempFitVar in the pT vs. TempFitVar plot"};
  ConfigurableAxis ConfTempFitVarpTBins{"ConfTempFitVarpTBins", {20, 0.5, 4.05}, "pT binning of the pT vs. TempFitVar plot"};

  /// Correlation part
  ConfigurableAxis ConfMultBins{"ConfMultBins", {VARIABLE_WIDTH, 0.0f, 4.0f, 8.0f, 12.0f, 16.0f, 20.0f, 24.0f, 28.0f, 32.0f, 36.0f, 40.0f, 44.0f, 48.0f, 52.0f, 56.0f, 60.0f, 64.0f, 68.0f, 72.0f, 76.0f, 80.0f, 84.0f, 88.0f, 92.0f, 96.0f, 100.0f, 200.0f, 99999.f}, "Mixing bins - multiplicity or centrality"}; // \todo to be obtained from the hash task
  ConfigurableAxis ConfMultKstarBins{"ConfMultKstarBins", {VARIABLE_WIDTH, 0.0f, 13.0f, 20.0f, 30.0f, 40.0f, 50.0f, 100.0f, 99999.f}, "Bins for kstar analysis in multiplicity or centrality bins (10 is maximum)"};
  ConfigurableAxis ConfKtKstarBins{"ConfKtKstarBins", {VARIABLE_WIDTH, 0.0f, 0.2f, 0.4f, 0.6f, 0.8f, 1.0f, 2.0f, 99999.f}, "Bins for kstar analysis in kT bins (10 is maximum)"};
  ConfigurableAxis ConfVtxBins{"ConfVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};

  ConfigurableAxis ConfmTBins3D{"ConfmTBins3D", {VARIABLE_WIDTH, 1.02f, 1.14f, 1.20f, 1.26f, 1.38f, 1.56f, 1.86f, 4.50f}, "mT Binning for the 3Dimensional plot: k* vs multiplicity vs mT (set <<twotracksconfigs.ConfUse3D>> to true in order to use)"};
  ConfigurableAxis ConfmultBins3D{"ConfmultBins3D", {VARIABLE_WIDTH, 0.0f, 20.0f, 30.0f, 40.0f, 99999.0f}, "multiplicity Binning for the 3Dimensional plot: k* vs multiplicity vs mT (set <<twotracksconfigs.ConfUse3D>> to true in order to use)"};

  ColumnBinningPolicy<aod::collision::PosZ, aod::femtouniversecollision::MultNtr> colBinning{{ConfVtxBins, ConfMultBins}, true};

  ConfigurableAxis ConfkstarBins{"ConfkstarBins", {1500, 0., 6.}, "binning kstar"};
  ConfigurableAxis ConfkTBins{"ConfkTBins", {150, 0., 9.}, "binning kT"};
  ConfigurableAxis ConfmTBins{"ConfmTBins", {225, 0., 7.5}, "binning mT"};
  Configurable<int> ConfNEventsMix{"ConfNEventsMix", 5, "Number of events for mixing"};
  Configurable<bool> ConfIsCPR{"ConfIsCPR", true, "Close Pair Rejection"};
  Configurable<bool> ConfCPRPlotPerRadii{"ConfCPRPlotPerRadii", false, "Plot CPR per radii"};
  Configurable<float> ConfCPRdeltaPhiCutMax{"ConfCPRdeltaPhiCutMax", 0.0, "Delta Phi max cut for Close Pair Rejection"};
  Configurable<float> ConfCPRdeltaPhiCutMin{"ConfCPRdeltaPhiCutMin", 0.0, "Delta Phi min cut for Close Pair Rejection"};
  Configurable<float> ConfCPRdeltaEtaCutMax{"ConfCPRdeltaEtaCutMax", 0.0, "Delta Eta max cut for Close Pair Rejection"};
  Configurable<float> ConfCPRdeltaEtaCutMin{"ConfCPRdeltaEtaCutMin", 0.0, "Delta Eta min cut for Close Pair Rejection"};
  Configurable<float> ConfCPRChosenRadii{"ConfCPRChosenRadii", 0.80, "Delta Eta cut for Close Pair Rejection"};

  Configurable<bool> IsPairIdentical{"IsPairIdentical", true, "'true' for identical particles, 'false' for non-identical particles"};
  Configurable<bool> cfgProcessPM{"cfgProcessPM", true, "Process differently charged particles (plus-minus)"};
  Configurable<bool> cfgProcessPP{"cfgProcessPP", true, "Process positively charged particles (plus-plus)"};
  Configurable<bool> cfgProcessMM{"cfgProcessMM", true, "Process negatively charged particles (minus-minus)"};
  Configurable<bool> cfgProcessMultBins{"cfgProcessMultBins", true, "Process kstar histograms (in multiplicity bins)"};
  Configurable<bool> cfgProcessKtBins{"cfgProcessKtBins", true, "Process kstar histograms in kT bins (if 'cfgProcessMultBins' is false, it will not be processed regardless of 'cfgProcessKtBins' state)"};
  Configurable<bool> cfgProcessKtMt3DCF{"cfgProcessKtMt3DCF", false, "Process 3D histograms in kT and MultBins"};

  FemtoUniverseFemtoContainer<femtoUniverseFemtoContainer::EventType::same, femtoUniverseFemtoContainer::Observable::kstar> sameEventCont;
  FemtoUniverseFemtoContainer<femtoUniverseFemtoContainer::EventType::mixed, femtoUniverseFemtoContainer::Observable::kstar> mixedEventCont;

  FemtoUniverseFemtoContainer<femtoUniverseFemtoContainer::EventType::same, femtoUniverseFemtoContainer::Observable::kstar> sameEventContPP;
  FemtoUniverseFemtoContainer<femtoUniverseFemtoContainer::EventType::mixed, femtoUniverseFemtoContainer::Observable::kstar> mixedEventContPP;

  FemtoUniverseFemtoContainer<femtoUniverseFemtoContainer::EventType::same, femtoUniverseFemtoContainer::Observable::kstar> sameEventContMM;
  FemtoUniverseFemtoContainer<femtoUniverseFemtoContainer::EventType::mixed, femtoUniverseFemtoContainer::Observable::kstar> mixedEventContMM;

  FemtoUniversePairCleaner<aod::femtouniverseparticle::ParticleType::kTrack, aod::femtouniverseparticle::ParticleType::kTrack> pairCleaner;
  FemtoUniverseDetaDphiStar<aod::femtouniverseparticle::ParticleType::kTrack, aod::femtouniverseparticle::ParticleType::kTrack> pairCloseRejection;
  FemtoUniverseTrackSelection trackCuts;

  PairWithCentMultKt sameEventMultCont;
  PairWithCentMultKt mixedEventMultCont;

  PairWithCentMultKt sameEventMultContPP;
  PairWithCentMultKt mixedEventMultContPP;

  PairWithCentMultKt sameEventMultContMM;
  PairWithCentMultKt mixedEventMultContMM;

  float mass1 = -1;
  float mass2 = -1;

  /// Histogram output
  HistogramRegistry qaRegistry{"TrackQA", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry resultRegistry{"Correlations", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
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

  HistogramRegistry sphericityRegistry{"SphericityHisto", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  /// TPC Pion/Kaon/Proton Sigma selection (general)
  bool IsNSigma(float mom, float nsigmaTPC, float nsigmaTOF)
  {
    // |nsigma_TPC| < 3 for p < 0.5 GeV/c
    // |nsigma_combined| < 3 for p > 0.5

    // using configurables:
    // ConfTOFpMin -- momentum value when we start using TOF; set to 1000 if TOF not needed
    // ConfNsigmaTPC -> TPC Sigma for momentum < ConfTOFpMin
    // ConfNsigmaCombined -> TPC and TOF Sigma (combined) for momentum > ConfTOFpMin

    if (mom < twotracksconfigs.ConfTOFpMin) {
      return TMath::Abs(nsigmaTPC) < twotracksconfigs.ConfNsigmaTPC;
    } else {
      return TMath::Hypot(nsigmaTOF, nsigmaTPC) < twotracksconfigs.ConfNsigmaCombined;
    }
  }

  /// TPC Kaon Sigma selection (stricter cuts for K+ and K-) -- based on Run2 results
  bool IsKaonNSigma(float mom, float nsigmaTPCK, float nsigmaTOFK)
  {
    if (twotracksconfigs.IsKaonNsigma == true) {
      if (mom < 0.4) {
        return TMath::Abs(nsigmaTPCK) < 2;
      } else if (mom > 0.4 && mom < 0.45) {
        return TMath::Abs(nsigmaTPCK) < 1;
      } else if (mom > 0.45 && mom < 0.8) {
        return (TMath::Abs(nsigmaTPCK) < 3 && TMath::Abs(nsigmaTOFK) < 2);
      } else if (mom > 0.8 && mom < 1.5) {
        return (TMath::Abs(nsigmaTPCK) < 3 && TMath::Abs(nsigmaTOFK) < 1.5);
      } else {
        return false;
      }
    } else {
      return IsNSigma(mom, nsigmaTPCK, nsigmaTOFK);
    }
  }

  /// TPC Deuteron Sigma selection
  bool IsDeuteronNSigma(float mom, float nsigmaTPCDe, float nsigmaTOFDe)
  {
    if (mom > deuteronconfigs.ConfPLowDe && mom < deuteronconfigs.ConfPHighDe) {
      if (mom < deuteronconfigs.ConfTOFpMinDe) {
        return TMath::Abs(nsigmaTPCDe) < deuteronconfigs.ConfNsigmaTPCDe;
      } else {
        return (TMath::Abs(nsigmaTOFDe) < deuteronconfigs.ConfNsigmaTOFDe && (TMath::Abs(nsigmaTPCDe) < deuteronconfigs.ConfNsigmaTPCDe));
      }
    } else {
      return false;
    }
  }

  /// Linear cut for clearer TPC Deuteron Sigma
  bool IsDeuteronNSigmaLinearCut(float mom, float nsigmaTPCDe, float nsigmaTOFDe, float tpcSignal)
  {
    if (lincut.ConfIsLine == true) {
      if (mom > lincut.pLow && mom < lincut.pHigh) {
        if (tpcSignal > lincut.a * mom + lincut.b) {
          return IsDeuteronNSigma(mom, nsigmaTPCDe, nsigmaTOFDe);
        } else {
          return false;
        }
      } else {
        return IsDeuteronNSigma(mom, nsigmaTPCDe, nsigmaTOFDe);
      }
    } else {
      return IsDeuteronNSigma(mom, nsigmaTPCDe, nsigmaTOFDe);
    }
  }

  /// Polynomial 3 cut for clearer TPC Deuteron Sigma
  bool IsDeuteronNSigmaPolCut(float mom, float nsigmaTPCDe, float nsigmaTOFDe, float tpcSignal)
  {
    if (polcut.ConfIsPol == true) {
      if (tpcSignal > polcut.A * TMath::Power(mom, 3) + polcut.B * TMath::Power(mom, 2) + polcut.C * mom + polcut.D) {
        return IsDeuteronNSigma(mom, nsigmaTPCDe, nsigmaTOFDe);
      } else {
        return false;
      }
    } else {
      return IsDeuteronNSigma(mom, nsigmaTPCDe, nsigmaTOFDe);
    }
  }

  bool IsParticleNSigma(int8_t particle_number, float mom, float nsigmaTPCPr, float nsigmaTOFPr, float nsigmaTPCPi, float nsigmaTOFPi, float nsigmaTPCK, float nsigmaTOFK, float nsigmaTPCDe, float nsigmaTOFDe, float tpcSignal)
  {
    if (particle_number == 1) {
      switch (trackonefilter.ConfPDGCodePartOne) {
        case 2212:  // Proton+
        case -2212: // Proton-
          return IsNSigma(mom, nsigmaTPCPr, nsigmaTOFPr);
          break;
        case 211:  // Pion+
        case -211: // Pion-
        case 111:  // Pion 0
          return IsNSigma(mom, nsigmaTPCPi, nsigmaTOFPi);
          break;
        case 321:  // Kaon+
        case -321: // Kaon-
          return IsKaonNSigma(mom, nsigmaTPCK, nsigmaTOFK);
          break;
        case 130: // Kaon 0 LONG
        case 310: // Kaon 0 SHORT
          return IsNSigma(mom, nsigmaTPCK, nsigmaTOFK);
          break;
        case 1000010020:  // Deuteron+
        case -1000010020: // Deuteron-
          return IsDeuteronNSigmaPolCut(mom, nsigmaTPCDe, nsigmaTOFDe, tpcSignal);
          break;
        default:
          return false;
      }
      return false;
    } else if (particle_number == 2) {
      switch (tracktwofilter.ConfPDGCodePartTwo) {
        case 2212:  // Proton+
        case -2212: // Proton-
          return IsNSigma(mom, nsigmaTPCPr, nsigmaTOFPr);
          break;
        case 211:  // Pion+
        case -211: // Pion-
        case 111:  // Pion 0
          return IsNSigma(mom, nsigmaTPCPi, nsigmaTOFPi);
          break;
        case 321:  // Kaon+
        case -321: // Kaon-
          return IsKaonNSigma(mom, nsigmaTPCK, nsigmaTOFK);
          break;
        case 130: // Kaon 0 LONG
        case 310: // Kaon 0 SHORT
          return IsNSigma(mom, nsigmaTPCK, nsigmaTOFK);
          break;
        case 1000010020:  // Deuteron+
        case -1000010020: // Deuteron-
          return IsDeuteronNSigmaPolCut(mom, nsigmaTPCDe, nsigmaTOFDe, tpcSignal);
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
    sphericityRegistry.add("sphericity", ";Sphericity;Entries", kTH1F, {{150, 0.0, 3, "Sphericity"}});

    trackHistoPartOne.init(&qaRegistry, ConfTempFitVarpTBins, ConfTempFitVarBins, twotracksconfigs.ConfIsMC, trackonefilter.ConfPDGCodePartOne, true);

    trackHistoPartTwo.init(&qaRegistry, ConfTempFitVarpTBins, ConfTempFitVarBins, twotracksconfigs.ConfIsMC, tracktwofilter.ConfPDGCodePartTwo, true);

    MixQaRegistry.add("MixingQA/hSECollisionBins", ";bin;Entries", kTH1F, {{120, -0.5, 119.5}});
    MixQaRegistry.add("MixingQA/hMECollisionBins", ";bin;Entries", kTH1F, {{120, -0.5, 119.5}});

    mass1 = pdg->Mass(trackonefilter.ConfPDGCodePartOne);
    mass2 = pdg->Mass(tracktwofilter.ConfPDGCodePartTwo);

    if (cfgProcessPM) {
      sameEventCont.init(&resultRegistryPM, ConfkstarBins, ConfMultBins, ConfkTBins, ConfmTBins, ConfmultBins3D, ConfmTBins3D, twotracksconfigs.ConfIsMC, twotracksconfigs.ConfUse3D);
      mixedEventCont.init(&resultRegistryPM, ConfkstarBins, ConfMultBins, ConfkTBins, ConfmTBins, ConfmultBins3D, ConfmTBins3D, twotracksconfigs.ConfIsMC, twotracksconfigs.ConfUse3D);

      sameEventCont.setPDGCodes(trackonefilter.ConfPDGCodePartOne, tracktwofilter.ConfPDGCodePartTwo);
      mixedEventCont.setPDGCodes(trackonefilter.ConfPDGCodePartOne, tracktwofilter.ConfPDGCodePartTwo);

      if (cfgProcessMultBins) {
        sameEventMultCont.init(&SameMultRegistryPM, ConfkstarBins, ConfMultKstarBins, ConfKtKstarBins, cfgProcessKtBins, cfgProcessKtMt3DCF);
        mixedEventMultCont.init(&MixedMultRegistryPM, ConfkstarBins, ConfMultKstarBins, ConfKtKstarBins, cfgProcessKtBins, cfgProcessKtMt3DCF);
      }
    }

    if (cfgProcessPP) {
      sameEventContPP.init(&resultRegistryPP, ConfkstarBins, ConfMultBins, ConfkTBins, ConfmTBins, ConfmultBins3D, ConfmTBins3D, twotracksconfigs.ConfIsMC, twotracksconfigs.ConfUse3D);
      mixedEventContPP.init(&resultRegistryPP, ConfkstarBins, ConfMultBins, ConfkTBins, ConfmTBins, ConfmultBins3D, ConfmTBins3D, twotracksconfigs.ConfIsMC, twotracksconfigs.ConfUse3D);
      sameEventContPP.setPDGCodes(trackonefilter.ConfPDGCodePartOne, tracktwofilter.ConfPDGCodePartTwo);
      mixedEventContPP.setPDGCodes(trackonefilter.ConfPDGCodePartOne, tracktwofilter.ConfPDGCodePartTwo);

      if (cfgProcessMultBins) {
        sameEventMultContPP.init(&SameMultRegistryPP, ConfkstarBins, ConfMultKstarBins, ConfKtKstarBins, cfgProcessKtBins, cfgProcessKtMt3DCF);
        mixedEventMultContPP.init(&MixedMultRegistryPP, ConfkstarBins, ConfMultKstarBins, ConfKtKstarBins, cfgProcessKtBins, cfgProcessKtMt3DCF);
      }
    }

    if (cfgProcessMM) {
      sameEventContMM.init(&resultRegistryMM, ConfkstarBins, ConfMultBins, ConfkTBins, ConfmTBins, ConfmultBins3D, ConfmTBins3D, twotracksconfigs.ConfIsMC, twotracksconfigs.ConfUse3D);
      mixedEventContMM.init(&resultRegistryMM, ConfkstarBins, ConfMultBins, ConfkTBins, ConfmTBins, ConfmultBins3D, ConfmTBins3D, twotracksconfigs.ConfIsMC, twotracksconfigs.ConfUse3D);
      sameEventContMM.setPDGCodes(trackonefilter.ConfPDGCodePartOne, tracktwofilter.ConfPDGCodePartTwo);
      mixedEventContMM.setPDGCodes(trackonefilter.ConfPDGCodePartOne, tracktwofilter.ConfPDGCodePartTwo);

      if (cfgProcessMultBins) {
        sameEventMultContMM.init(&SameMultRegistryMM, ConfkstarBins, ConfMultKstarBins, ConfKtKstarBins, cfgProcessKtBins, cfgProcessKtMt3DCF);
        mixedEventMultContMM.init(&MixedMultRegistryMM, ConfkstarBins, ConfMultKstarBins, ConfKtKstarBins, cfgProcessKtBins, cfgProcessKtMt3DCF);
      }
    }

    pairCleaner.init(&qaRegistry);
    if (ConfIsCPR.value) {
      pairCloseRejection.init(&resultRegistry, &qaRegistry, ConfCPRdeltaPhiCutMin.value, ConfCPRdeltaPhiCutMax.value, ConfCPRdeltaEtaCutMin.value, ConfCPRdeltaEtaCutMax.value, ConfCPRChosenRadii.value, ConfCPRPlotPerRadii.value);
    }

    vPIDPartOne = trackonefilter.ConfPIDPartOne.value;
    vPIDPartTwo = tracktwofilter.ConfPIDPartTwo.value;
    kNsigma = twotracksconfigs.ConfTrkPIDnSigmaMax.value;
  }

  template <typename CollisionType>
  void fillCollision(CollisionType col)
  {
    MixQaRegistry.fill(HIST("MixingQA/hSECollisionBins"), colBinning.getBin({col.posZ(), col.multNtr()}));
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

    /// Histogramming same event
    if ((pairType == 1 || pairType == 2) && fillQA) {
      for (auto& part : groupPartsOne) {
        if (!IsParticleNSigma((int8_t)1, part.p(), trackCuts.getNsigmaTPC(part, o2::track::PID::Proton), trackCuts.getNsigmaTOF(part, o2::track::PID::Proton), trackCuts.getNsigmaTPC(part, o2::track::PID::Pion), trackCuts.getNsigmaTOF(part, o2::track::PID::Pion), trackCuts.getNsigmaTPC(part, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(part, o2::track::PID::Kaon), trackCuts.getNsigmaTPC(part, o2::track::PID::Deuteron), trackCuts.getNsigmaTOF(part, o2::track::PID::Deuteron), part.tpcSignal())) {
          continue;
        }
        trackHistoPartOne.fillQA<isMC, true>(part);
      }
    }

    if ((pairType == 1 || pairType == 3) && fillQA) {
      for (auto& part : groupPartsTwo) {
        if (!IsParticleNSigma((int8_t)2, part.p(), trackCuts.getNsigmaTPC(part, o2::track::PID::Proton), trackCuts.getNsigmaTOF(part, o2::track::PID::Proton), trackCuts.getNsigmaTPC(part, o2::track::PID::Pion), trackCuts.getNsigmaTOF(part, o2::track::PID::Pion), trackCuts.getNsigmaTPC(part, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(part, o2::track::PID::Kaon), trackCuts.getNsigmaTPC(part, o2::track::PID::Deuteron), trackCuts.getNsigmaTOF(part, o2::track::PID::Deuteron), part.tpcSignal())) {
          continue;
        }
        trackHistoPartTwo.fillQA<isMC, true>(part);
      }
    }

    if (pairType == 1) {

      /// Now build the combinations for non-identical particle pairs
      for (auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(groupPartsOne, groupPartsTwo))) {

        if (!IsParticleNSigma((int8_t)1, p1.p(), trackCuts.getNsigmaTPC(p1, o2::track::PID::Proton), trackCuts.getNsigmaTOF(p1, o2::track::PID::Proton), trackCuts.getNsigmaTPC(p1, o2::track::PID::Pion), trackCuts.getNsigmaTOF(p1, o2::track::PID::Pion), trackCuts.getNsigmaTPC(p1, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(p1, o2::track::PID::Kaon), trackCuts.getNsigmaTPC(p1, o2::track::PID::Deuteron), trackCuts.getNsigmaTOF(p1, o2::track::PID::Deuteron), p1.tpcSignal())) {
          continue;
        }

        if (!IsParticleNSigma((int8_t)2, p2.p(), trackCuts.getNsigmaTPC(p2, o2::track::PID::Proton), trackCuts.getNsigmaTOF(p2, o2::track::PID::Proton), trackCuts.getNsigmaTPC(p2, o2::track::PID::Pion), trackCuts.getNsigmaTOF(p2, o2::track::PID::Pion), trackCuts.getNsigmaTPC(p2, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(p2, o2::track::PID::Kaon), trackCuts.getNsigmaTPC(p2, o2::track::PID::Deuteron), trackCuts.getNsigmaTOF(p2, o2::track::PID::Deuteron), p2.tpcSignal())) {
          continue;
        }

        if (ConfIsCPR.value) {
          if (pairCloseRejection.isClosePair(p1, p2, parts, magFieldTesla, femtoUniverseContainer::EventType::same)) {
            continue;
          }
        }

        // track cleaning
        if (!pairCleaner.isCleanPair(p1, p2, parts)) {
          continue;
        }

        float kstar = FemtoUniverseMath::getkstar(p1, mass1, p2, mass2);
        float kT = FemtoUniverseMath::getkT(p1, mass1, p2, mass2);

        sameEventCont.setPair<isMC>(p1, p2, multCol, twotracksconfigs.ConfUse3D);
        if (cfgProcessMultBins)
          sameEventMultCont.fill<float>(kstar, multCol, kT);
      }
    } else {
      /// Now build the combinations for identical particles pairs
      for (auto& [p1, p2] : combinations(CombinationsStrictlyUpperIndexPolicy(groupPartsOne, groupPartsOne))) {

        if (!IsParticleNSigma((int8_t)2, p1.p(), trackCuts.getNsigmaTPC(p1, o2::track::PID::Proton), trackCuts.getNsigmaTOF(p1, o2::track::PID::Proton), trackCuts.getNsigmaTPC(p1, o2::track::PID::Pion), trackCuts.getNsigmaTOF(p1, o2::track::PID::Pion), trackCuts.getNsigmaTPC(p1, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(p1, o2::track::PID::Kaon), trackCuts.getNsigmaTPC(p1, o2::track::PID::Deuteron), trackCuts.getNsigmaTOF(p1, o2::track::PID::Deuteron), p1.tpcSignal())) {
          continue;
        }

        if (!IsParticleNSigma((int8_t)2, p2.p(), trackCuts.getNsigmaTPC(p2, o2::track::PID::Proton), trackCuts.getNsigmaTOF(p2, o2::track::PID::Proton), trackCuts.getNsigmaTPC(p2, o2::track::PID::Pion), trackCuts.getNsigmaTOF(p2, o2::track::PID::Pion), trackCuts.getNsigmaTPC(p2, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(p2, o2::track::PID::Kaon), trackCuts.getNsigmaTPC(p2, o2::track::PID::Deuteron), trackCuts.getNsigmaTOF(p2, o2::track::PID::Deuteron), p2.tpcSignal())) {
          continue;
        }

        if (ConfIsCPR.value) {
          if (pairCloseRejection.isClosePair(p1, p2, parts, magFieldTesla, femtoUniverseContainer::EventType::same)) {
            continue;
          }
        }

        // track cleaning
        if (!pairCleaner.isCleanPair(p1, p2, parts)) {
          continue;
        }

        switch (pairType) {
          case 2: {
            if (IsPairIdentical == true) {
              float kstar = FemtoUniverseMath::getkstar(p1, mass1, p2, mass1);
              float kT = FemtoUniverseMath::getkT(p1, mass1, p2, mass1);

              sameEventContPP.setPair<isMC>(p1, p2, multCol, twotracksconfigs.ConfUse3D);
              if (cfgProcessMultBins)
                sameEventMultContPP.fill<float>(kstar, multCol, kT);
            } else {
              float kstar = FemtoUniverseMath::getkstar(p1, mass1, p2, mass2);
              float kT = FemtoUniverseMath::getkT(p1, mass1, p2, mass2);

              sameEventContPP.setPair<isMC>(p1, p2, multCol, twotracksconfigs.ConfUse3D);
              if (cfgProcessMultBins)
                sameEventMultContPP.fill<float>(kstar, multCol, kT);
            }

            break;
          }

          case 3: {
            if (IsPairIdentical == true) {
              float kstar = FemtoUniverseMath::getkstar(p1, mass2, p2, mass2);
              float kT = FemtoUniverseMath::getkT(p1, mass2, p2, mass2);

              sameEventContMM.setPair<isMC>(p1, p2, multCol, twotracksconfigs.ConfUse3D);
              if (cfgProcessMultBins)
                sameEventMultContMM.fill<float>(kstar, multCol, kT);
            } else {
              float kstar = FemtoUniverseMath::getkstar(p1, mass1, p2, mass2);
              float kT = FemtoUniverseMath::getkT(p1, mass1, p2, mass2);

              sameEventContMM.setPair<isMC>(p1, p2, multCol, twotracksconfigs.ConfUse3D);
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
  void processSameEvent(soa::Filtered<o2::aod::FDCollisions>::iterator& col,
                        FilteredFemtoFullParticles& parts)
  {
    fillCollision(col);
    sphericityRegistry.fill(HIST("sphericity"), col.sphericity());

    auto thegroupPartsOne = partsOne->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto thegroupPartsTwo = partsTwo->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);

    bool fillQA = true;

    if (cfgProcessPM) {
      doSameEvent<false>(thegroupPartsOne, thegroupPartsTwo, parts, col.magField(), col.multNtr(), 1, fillQA);
      fillQA = false;
    }
    if (cfgProcessPP)
      doSameEvent<false>(thegroupPartsOne, thegroupPartsOne, parts, col.magField(), col.multNtr(), 2, fillQA);
    if (cfgProcessMM)
      doSameEvent<false>(thegroupPartsTwo, thegroupPartsTwo, parts, col.magField(), col.multNtr(), 3, fillQA);
  }
  PROCESS_SWITCH(femtoUniversePairTaskTrackTrackMultKtExtended, processSameEvent, "Enable processing same event", true);

  /// process function for to call doSameEvent with Monte Carlo
  /// \param col subscribe to the collision table (Monte Carlo Reconstructed reconstructed)
  /// \param parts subscribe to joined table FemtoUniverseParticles and FemtoUniverseMCLables to access Monte Carlo truth
  /// \param FemtoUniverseMCParticles subscribe to the Monte Carlo truth table
  void processSameEventMC(o2::aod::FDCollision& col,
                          soa::Join<FilteredFemtoFullParticles, aod::FDMCLabels>& parts,
                          o2::aod::FDMCParticles&)
  {
    fillCollision(col);

    auto thegroupPartsOne = partsOneMC->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto thegroupPartsTwo = partsTwoMC->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);

    bool fillQA = true;
    if (cfgProcessPM) {
      doSameEvent<false>(thegroupPartsOne, thegroupPartsTwo, parts, col.magField(), col.multNtr(), 1, fillQA);
      fillQA = false;
    }
    if (cfgProcessPP)
      doSameEvent<true>(thegroupPartsOne, thegroupPartsOne, parts, col.magField(), col.multNtr(), 2, fillQA);
    if (cfgProcessMM)
      doSameEvent<true>(thegroupPartsTwo, thegroupPartsTwo, parts, col.magField(), col.multNtr(), 3, fillQA);
  }
  PROCESS_SWITCH(femtoUniversePairTaskTrackTrackMultKtExtended, processSameEventMC, "Enable processing same event for Monte Carlo", false);

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

    for (auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(groupPartsOne, groupPartsTwo))) {

      if (!IsParticleNSigma((int8_t)2, p1.p(), trackCuts.getNsigmaTPC(p1, o2::track::PID::Proton), trackCuts.getNsigmaTOF(p1, o2::track::PID::Proton), trackCuts.getNsigmaTPC(p1, o2::track::PID::Pion), trackCuts.getNsigmaTOF(p1, o2::track::PID::Pion), trackCuts.getNsigmaTPC(p1, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(p1, o2::track::PID::Kaon), trackCuts.getNsigmaTPC(p1, o2::track::PID::Deuteron), trackCuts.getNsigmaTOF(p1, o2::track::PID::Deuteron), p1.tpcSignal())) {
        continue;
      }

      if (!IsParticleNSigma((int8_t)2, p2.p(), trackCuts.getNsigmaTPC(p2, o2::track::PID::Proton), trackCuts.getNsigmaTOF(p2, o2::track::PID::Proton), trackCuts.getNsigmaTPC(p2, o2::track::PID::Pion), trackCuts.getNsigmaTOF(p2, o2::track::PID::Pion), trackCuts.getNsigmaTPC(p2, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(p2, o2::track::PID::Kaon), trackCuts.getNsigmaTPC(p2, o2::track::PID::Deuteron), trackCuts.getNsigmaTOF(p2, o2::track::PID::Deuteron), p2.tpcSignal())) {
        continue;
      }

      if (ConfIsCPR.value) {
        if (pairCloseRejection.isClosePair(p1, p2, parts, magFieldTesla, femtoUniverseContainer::EventType::mixed)) {
          continue;
        }
      }

      switch (pairType) {
        case 1: {
          float kstar = FemtoUniverseMath::getkstar(p1, mass1, p2, mass2);
          float kT = FemtoUniverseMath::getkT(p1, mass1, p2, mass2);

          mixedEventCont.setPair<isMC>(p1, p2, multCol, twotracksconfigs.ConfUse3D);
          if (cfgProcessMultBins)
            mixedEventMultCont.fill<float>(kstar, multCol, kT);

          break;
        }
        case 2: {
          if (IsPairIdentical == true) {
            float kstar = FemtoUniverseMath::getkstar(p1, mass1, p2, mass1);
            float kT = FemtoUniverseMath::getkT(p1, mass1, p2, mass1);

            mixedEventContPP.setPair<isMC>(p1, p2, multCol, twotracksconfigs.ConfUse3D);
            if (cfgProcessMultBins)
              mixedEventMultContPP.fill<float>(kstar, multCol, kT);
          } else {
            float kstar = FemtoUniverseMath::getkstar(p1, mass1, p2, mass2);
            float kT = FemtoUniverseMath::getkT(p1, mass1, p2, mass2);

            mixedEventContPP.setPair<isMC>(p1, p2, multCol, twotracksconfigs.ConfUse3D);
            if (cfgProcessMultBins)
              mixedEventMultContPP.fill<float>(kstar, multCol, kT);
          }

          break;
        }

        case 3: {
          if (IsPairIdentical == true) {
            float kstar = FemtoUniverseMath::getkstar(p1, mass2, p2, mass2);
            float kT = FemtoUniverseMath::getkT(p1, mass2, p2, mass2);

            mixedEventContMM.setPair<isMC>(p1, p2, multCol, twotracksconfigs.ConfUse3D);
            if (cfgProcessMultBins)
              mixedEventMultContMM.fill<float>(kstar, multCol, kT);
          } else {
            float kstar = FemtoUniverseMath::getkstar(p1, mass1, p2, mass2);
            float kT = FemtoUniverseMath::getkT(p1, mass1, p2, mass2);

            mixedEventContMM.setPair<isMC>(p1, p2, multCol, twotracksconfigs.ConfUse3D);
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
  void processMixedEvent(soa::Filtered<o2::aod::FDCollisions>& cols,
                         FilteredFemtoFullParticles& parts)
  {
    for (auto& [collision1, collision2] : soa::selfCombinations(colBinning, 5, -1, cols, cols)) {

      const int multiplicityCol = collision1.multNtr();
      MixQaRegistry.fill(HIST("MixingQA/hMECollisionBins"), colBinning.getBin({collision1.posZ(), multiplicityCol}));

      const auto& magFieldTesla1 = collision1.magField();
      const auto& magFieldTesla2 = collision2.magField();

      if (magFieldTesla1 != magFieldTesla2) {
        continue;
      }

      if (cfgProcessPM) {
        auto groupPartsOne = partsOne->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);
        auto groupPartsTwo = partsTwo->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision2.globalIndex(), cache);
        doMixedEvent<false>(groupPartsOne, groupPartsTwo, parts, magFieldTesla1, multiplicityCol, 1);
      }
      if (cfgProcessPP) {
        auto groupPartsOne = partsOne->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);
        auto groupPartsTwo = partsOne->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision2.globalIndex(), cache);
        doMixedEvent<false>(groupPartsOne, groupPartsTwo, parts, magFieldTesla1, multiplicityCol, 2);
      }
      if (cfgProcessMM) {
        auto groupPartsOne = partsTwo->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);
        auto groupPartsTwo = partsTwo->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision2.globalIndex(), cache);
        doMixedEvent<false>(groupPartsOne, groupPartsTwo, parts, magFieldTesla1, multiplicityCol, 3);
      }
    }
  }
  PROCESS_SWITCH(femtoUniversePairTaskTrackTrackMultKtExtended, processMixedEvent, "Enable processing mixed events", true);

  /// brief process function for to call doMixedEvent with Monte Carlo
  /// \param cols subscribe to the collisions table (Monte Carlo Reconstructed reconstructed)
  /// \param parts subscribe to joined table FemtoUniverseParticles and FemtoUniverseMCLables to access Monte Carlo truth
  /// \param FemtoUniverseMCParticles subscribe to the Monte Carlo truth table
  void processMixedEventMC(o2::aod::FDCollisions& cols,
                           soa::Join<FilteredFemtoFullParticles, aod::FDMCLabels>& parts,
                           o2::aod::FDMCParticles&)
  {
    for (auto& [collision1, collision2] : soa::selfCombinations(colBinning, 5, -1, cols, cols)) {

      const int multiplicityCol = collision1.multNtr();
      MixQaRegistry.fill(HIST("MixingQA/hMECollisionBins"), colBinning.getBin({collision1.posZ(), multiplicityCol}));

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
        doMixedEvent<true>(groupPartsOne, groupPartsTwo, parts, magFieldTesla1, multiplicityCol, 1);
      }
      if (cfgProcessPP) {
        auto groupPartsOne = partsOneMC->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);
        auto groupPartsTwo = partsOneMC->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision2.globalIndex(), cache);
        doMixedEvent<true>(groupPartsOne, groupPartsTwo, parts, magFieldTesla1, multiplicityCol, 2);
      }
      if (cfgProcessMM) {
        auto groupPartsOne = partsTwoMC->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);
        auto groupPartsTwo = partsTwoMC->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision2.globalIndex(), cache);
        doMixedEvent<true>(groupPartsOne, groupPartsTwo, parts, magFieldTesla1, multiplicityCol, 3);
      }
    }
  }
  PROCESS_SWITCH(femtoUniversePairTaskTrackTrackMultKtExtended, processMixedEventMC, "Enable processing mixed events MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<femtoUniversePairTaskTrackTrackMultKtExtended>(cfgc),
  };
  return workflow;
}
