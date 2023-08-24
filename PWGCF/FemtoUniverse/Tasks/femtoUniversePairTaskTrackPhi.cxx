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

/// \file femtoUniversePairTaskTrackPhi.cxx
/// \brief Tasks that reads the track tables used for the pairing and builds pairs of two tracks
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de
/// \author Georgios Mantzaridis, TU München, georgios.mantzaridis@tum.de
/// \author Anton Riedel, TU München, anton.riedel@tum.de
/// \author Zuzanna Chochulska, WUT Warsaw, zuzanna.chochulska.stud@pw.edu.pl

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
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseAngularContainer.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseDetaDphiStar.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUtils.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseTrackSelection.h"

using namespace o2;
using namespace o2::analysis::femtoUniverse;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

namespace
{
static constexpr int nPart = 2;
static constexpr int nCuts = 5;
static const std::vector<std::string> partNames{"PhiCandidate", "Hadron"};
static const std::vector<std::string> cutNames{"MaxPt", "PIDthr", "nSigmaTPC", "nSigmaTPCTOF", "MaxP"};
static const float cutsTable[nPart][nCuts]{
  {4.05f, 1.f, 3.f, 3.f, 100.f},
  {4.05f, 1.f, 3.f, 3.f, 100.f}};
} // namespace

struct femtoUniversePairTaskTrackPhi {

  using FemtoFullParticles = soa::Join<aod::FDParticles, aod::FDExtParticles>;
  SliceCache cache;
  Preslice<FemtoFullParticles> perCol = aod::femtouniverseparticle::fdCollisionId;

  /// Particle selection part

  /// Table for both particles
  struct : o2::framework::ConfigurableGroup {
    Configurable<float> ConfNsigmaCombinedKaon{"ConfNsigmaCombinedKaon", 3.0, "TPC and TOF Kaon Sigma (combined) for momentum > 0.4"};
    Configurable<float> ConfNsigmaTPCKaon{"ConfNsigmaTPCKaon", 3.0, "TPC Kaon Sigma for momentum < 0.4"};
    Configurable<float> ConfNsigmaCombinedProton{"ConfNsigmaCombinedProton", 3.0, "TPC and TOF Proton Sigma (combined) for momentum > 0.5"};
    Configurable<float> ConfNsigmaTPCProton{"ConfNsigmaTPCProton", 3.0, "TPC Proton Sigma for momentum < 0.5"};
    Configurable<float> ConfNsigmaCombinedPion{"ConfNsigmaCombinedPion", 3.0, "TPC and TOF Pion Sigma (combined) for momentum > 0.5"};
    Configurable<float> ConfNsigmaTPCPion{"ConfNsigmaTPCPion", 3.0, "TPC Pion Sigma for momentum < 0.5"};

    Configurable<LabeledArray<float>> ConfCutTable{"ConfCutTable", {cutsTable[0], nPart, nCuts, partNames, cutNames}, "Particle selections"};
    Configurable<int> ConfNspecies{"ConfNspecies", 2, "Number of particle spieces with PID info"};
    Configurable<bool> ConfIsMC{"ConfIsMC", false, "Enable additional Histogramms in the case of a MonteCarlo Run"};
    Configurable<std::vector<float>> ConfTrkPIDnSigmaMax{"ConfTrkPIDnSigmaMax", std::vector<float>{4.f, 3.f, 2.f}, "This configurable needs to be the same as the one used in the producer task"};
    Configurable<bool> ConfUse3D{"ConfUse3D", false, "Enable three dimensional histogramms (to be used only for analysis with high statistics): k* vs mT vs multiplicity"};
    Configurable<int> ConfPhiBins{"ConfPhiBins", 29, "Number of phi bins in deta dphi"};
    Configurable<int> ConfEtaBins{"ConfEtaBins", 29, "Number of eta bins in deta dphi"};

  } twotracksconfigs;

  /// Particle 1 --- IDENTIFIED HADRON
  Configurable<bool> ConfIsSame{"ConfIsSame", false, "Pairs of the same particle"};
  struct : o2::framework::ConfigurableGroup {
    Configurable<int> ConfTrackChoiceHadron{"ConfTrackChoiceHadron", 0, "Type of particle (track): {0:Proton, 1:Pion, 2:Kaon}"};
    Configurable<int> ConfPDGCodeHadron{"ConfPDGCodeHadron", 2212, "Particle 2 - PDG code"};
    // Configurable<uint32_t> ConfCutHadron{"ConfCutHadron", 5542474, "Particle 2 - Selection bit"};
    Configurable<int> ConfPIDHadron{"ConfPIDHadron", 2, "Particle 2 - Read from cutCulator"}; // we also need the possibility to specify whether the bit is true/false ->std>>vector<std::pair<int, int>>
  } trackHadronfilter;

  /// Partition for particle 1
  Partition<FemtoFullParticles> partsHadron = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kTrack));
  Partition<soa::Join<aod::FDParticles, aod::FDMCLabels>> partsHadronMC = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kTrack));

  /// Histogramming for particle 1
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kTrack, 2> trackHistoPartHadron;

  /// Particle 2 --- PHI
  struct : o2::framework::ConfigurableGroup {
    Configurable<int> ConfPDGCodePhi{"ConfPDGCodePhi", 333, "Phi meson - PDG code"};
    Configurable<uint32_t> ConfCutPhi{"ConfCutPhi", 5542474, "Phi meson - Selection bit from cutCulator"};
    Configurable<int> ConfPIDPhi{"ConfPIDPhi", 2, "Phi meson - Read from cutCulator"};           // we also need the possibility to specify whether the bit is true/false ->std>>vector<std::pair<int, int>>int>>
    Configurable<float> cfgPtLowPhi{"cfgPtLowPhi", 0.5, "Lower limit for Pt for the Phi meson"}; // change according to wrzesa cuts
    Configurable<float> cfgPtHighPhi{"cfgPtHighPhi", 1.5, "Higher limit for Pt for the Phi meson"};
  } trackPhifilter;

  /// Partition for particle 2
  Partition<FemtoFullParticles> partsPhi = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kPhi));
  Partition<soa::Join<aod::FDParticles, aod::FDMCLabels>> partsPhiMC = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kPhi));

  /// Histogramming for particle 2
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kPhi, 0> trackHistoPartPhi;

  /// Histogramming for Event
  FemtoUniverseEventHisto eventHisto;

  /// The configurables need to be passed to an std::vector
  int vPIDPhiCandidate, vPIDHadron;
  std::vector<float> kNsigma;

  /// particle part
  ConfigurableAxis ConfTempFitVarBins{"ConfDTempFitVarBins", {300, -0.15, 0.15}, "binning of the TempFitVar in the pT vs. TempFitVar plot"};
  ConfigurableAxis ConfTempFitVarInvMassBins{"ConfDTempFitVarInvMassBins", {6000, 0.9, 4.0}, "binning of the TempFitVar in the pT vs. TempFitVar plot"};
  ConfigurableAxis ConfTempFitVarpTBins{"ConfTempFitVarpTBins", {20, 0.5, 4.05}, "pT binning of the pT vs. TempFitVar plot"};

  /// Correlation part
  ConfigurableAxis ConfMultBins{"ConfMultBins", {VARIABLE_WIDTH, 0.0f, 4.0f, 8.0f, 12.0f, 16.0f, 20.0f, 24.0f, 28.0f, 32.0f, 36.0f, 40.0f, 44.0f, 48.0f, 52.0f, 56.0f, 60.0f, 64.0f, 68.0f, 72.0f, 76.0f, 80.0f, 84.0f, 88.0f, 92.0f, 96.0f, 100.0f, 200.0f, 99999.f}, "Mixing bins - multiplicity"}; // \todo to be obtained from the hash task
  // ConfigurableAxis ConfMultBins{"CfgMultBins", {VARIABLE_WIDTH, 0.0f, 20.0f, 40.0f, 60.0f, 80.0f, 100.0f, 200.0f, 99999.f}, "Mixing bins - multiplicity"};
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
  Configurable<float> ConfCPRdeltaPhiMax{"ConfCPRdeltaPhiMax", 0.01, "Max. Delta Phi for Close Pair Rejection"};
  Configurable<float> ConfCPRdeltaEtaMax{"ConfCPRdeltaEtaMax", 0.01, "Max. Delta Eta for Close Pair Rejection"};

  FemtoUniverseFemtoContainer<femtoUniverseFemtoContainer::EventType::same, femtoUniverseFemtoContainer::Observable::kstar> sameEventFemtoCont;
  FemtoUniverseFemtoContainer<femtoUniverseFemtoContainer::EventType::mixed, femtoUniverseFemtoContainer::Observable::kstar> mixedEventFemtoCont;
  FemtoUniverseAngularContainer<femtoUniverseAngularContainer::EventType::same, femtoUniverseAngularContainer::Observable::kstar> sameEventAngularCont;
  FemtoUniverseAngularContainer<femtoUniverseAngularContainer::EventType::mixed, femtoUniverseAngularContainer::Observable::kstar> mixedEventAngularCont;
  FemtoUniversePairCleaner<aod::femtouniverseparticle::ParticleType::kTrack, aod::femtouniverseparticle::ParticleType::kPhi> pairCleaner;
  FemtoUniverseDetaDphiStar<aod::femtouniverseparticle::ParticleType::kTrack, aod::femtouniverseparticle::ParticleType::kPhi> pairCloseRejection;
  FemtoUniverseTrackSelection trackCuts;
  /// Histogram output
  HistogramRegistry qaRegistry{"TrackQA", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry resultRegistry{"Correlations", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry MixQaRegistry{"MixQaRegistry", {}, OutputObjHandlingPolicy::AnalysisObject};

  // PID for protons
  bool IsProtonNSigma(float mom, float nsigmaTPCPr, float nsigmaTOFPr) // previous version from: https://github.com/alisw/AliPhysics/blob/master/PWGCF/FEMTOSCOPY/AliFemtoUser/AliFemtoMJTrackCut.cxx
  {
    //|nsigma_TPC| < 3 for p < 0.5 GeV/c
    //|nsigma_combined| < 3 for p > 0.5

    // using configurables:
    // ConfNsigmaTPCProton -> TPC Kaon Sigma for momentum < 0.5
    // ConfNsigmaCombinedProton -> TPC and TOF Kaon Sigma (combined) for momentum > 0.5

    if (mom < 0.5) {
      if (TMath::Abs(nsigmaTPCPr) < twotracksconfigs.ConfNsigmaTPCProton) {
        return true;
      } else {
        return false;
      }
    } else if (mom > 0.4) {
      // if (TMath::Hypot(nsigmaTOFPr, nsigmaTPCPr) < twotracksconfigs.ConfNsigmaCombinedProton) {
      if (TMath::Abs(nsigmaTPCPr) < twotracksconfigs.ConfNsigmaCombinedProton) {
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
      if (TMath::Abs(nsigmaTPCK) < 3.0) {
        return true;
      } else {
        return false;
      }
    } else if (mom < 0.45) { // 0.30 - 0.45
      if (TMath::Abs(nsigmaTPCK) < 2.0) {
        return true;
      } else {
        return false;
      }
    } else if (mom < 0.55) { // 0.45-0.55
      if (TMath::Abs(nsigmaTPCK) < 1.0) {
        return true;
      } else {
        return false;
      }
    } else if (mom < 1.5) { // 0.55-1.5 (now we use TPC and TOF)
      if ((TMath::Abs(nsigmaTOFK) < 3.0) && (TMath::Abs(nsigmaTPCK) < 3.0)) {
        {
          return true;
        }
      } else {
        return false;
      }
    } else if (mom > 1.5) { // 1.5 -
      if ((TMath::Abs(nsigmaTOFK) < 2.0) && (TMath::Abs(nsigmaTPCK) < 3.0)) {
        return true;
      } else {
        return false;
      }
    } else {
      return false;
    }
  }

  bool IsPionNSigma(float mom, float nsigmaTPCPi, float nsigmaTOFPi)
  {
    //|nsigma_TPC| < 3 for p < 0.5 GeV/c
    //|nsigma_combined| < 3 for p > 0.5

    // using configurables:
    // ConfNsigmaTPCPion -> TPC Kaon Sigma for momentum < 0.5
    // ConfNsigmaCombinedPion -> TPC and TOF Pion Sigma (combined) for momentum > 0.5
    if (true) {
      if (mom < 0.5) {
        if (TMath::Abs(nsigmaTPCPi) < twotracksconfigs.ConfNsigmaTPCPion) {
          return true;
        } else {
          return false;
        }
      } else if (mom > 0.5) {
        // if (TMath::Hypot(nsigmaTOFPi, nsigmaTPCPi) < twotracksconfigs.ConfNsigmaCombinedPion) {
        if (TMath::Abs(nsigmaTPCPi) < twotracksconfigs.ConfNsigmaCombinedPion) {
          return true;
        } else {
          return false;
        }
      }
    }
    return false;
  }

  bool IsParticleNSigma(float mom, float nsigmaTPCPr, float nsigmaTOFPr, float nsigmaTPCPi, float nsigmaTOFPi, float nsigmaTPCK, float nsigmaTOFK)
  {
    switch (trackHadronfilter.ConfTrackChoiceHadron) {
      case 0: // Proton
        return IsProtonNSigma(mom, nsigmaTPCPr, nsigmaTOFPr);
        break;
      case 1: // Pion
        return IsPionNSigma(mom, nsigmaTPCPi, nsigmaTOFPi);
        break;
      case 2: // Kaon
        return IsKaonNSigma(mom, nsigmaTPCK, nsigmaTOFK);
        break;
      default:
        return false;
    }
  }

  void init(InitContext&)
  {
    eventHisto.init(&qaRegistry);
    trackHistoPartPhi.init(&qaRegistry, ConfTempFitVarpTBins, ConfTempFitVarInvMassBins, twotracksconfigs.ConfIsMC, trackPhifilter.ConfPDGCodePhi);
    if (!ConfIsSame) {
      trackHistoPartHadron.init(&qaRegistry, ConfTempFitVarpTBins, ConfTempFitVarBins, twotracksconfigs.ConfIsMC, trackHadronfilter.ConfPDGCodeHadron);
    }

    MixQaRegistry.add("MixingQA/hSECollisionBins", ";bin;Entries", kTH1F, {{120, -0.5, 119.5}});
    MixQaRegistry.add("MixingQA/hMECollisionBins", ";bin;Entries", kTH1F, {{120, -0.5, 119.5}});

    sameEventFemtoCont.init(&resultRegistry, ConfkstarBins, ConfMultBins, ConfkTBins, ConfmTBins, ConfmultBins3D, ConfmTBins3D, twotracksconfigs.ConfIsMC, twotracksconfigs.ConfUse3D);
    mixedEventFemtoCont.init(&resultRegistry, ConfkstarBins, ConfMultBins, ConfkTBins, ConfmTBins, ConfmultBins3D, ConfmTBins3D, twotracksconfigs.ConfIsMC, twotracksconfigs.ConfUse3D);
    sameEventAngularCont.init(&resultRegistry, ConfkstarBins, ConfMultBins, ConfkTBins, ConfmTBins, ConfmultBins3D, ConfmTBins3D, twotracksconfigs.ConfEtaBins, twotracksconfigs.ConfPhiBins, twotracksconfigs.ConfIsMC, twotracksconfigs.ConfUse3D);
    mixedEventAngularCont.init(&resultRegistry, ConfkstarBins, ConfMultBins, ConfkTBins, ConfmTBins, ConfmultBins3D, ConfmTBins3D, twotracksconfigs.ConfEtaBins, twotracksconfigs.ConfPhiBins, twotracksconfigs.ConfIsMC, twotracksconfigs.ConfUse3D);

    sameEventFemtoCont.setPDGCodes(trackPhifilter.ConfPDGCodePhi, trackHadronfilter.ConfPDGCodeHadron);
    mixedEventFemtoCont.setPDGCodes(trackPhifilter.ConfPDGCodePhi, trackHadronfilter.ConfPDGCodeHadron);
    sameEventAngularCont.setPDGCodes(trackPhifilter.ConfPDGCodePhi, trackHadronfilter.ConfPDGCodeHadron);
    mixedEventAngularCont.setPDGCodes(trackPhifilter.ConfPDGCodePhi, trackHadronfilter.ConfPDGCodeHadron);

    pairCleaner.init(&qaRegistry);
    if (ConfIsCPR.value) {
      pairCloseRejection.init(&resultRegistry, &qaRegistry, ConfCPRdeltaPhiMax.value, ConfCPRdeltaEtaMax.value, ConfCPRPlotPerRadii.value);
    }

    vPIDPhiCandidate = trackPhifilter.ConfPIDPhi.value;
    vPIDHadron = trackHadronfilter.ConfPIDHadron.value;
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
  /// @tparam PartitionType
  /// @tparam PartType
  /// @tparam isMC: enables Monte Carlo truth specific histograms
  /// @param groupPartsHadron partition for the first particle passed by the process function
  /// @param groupPartsPhi partition for the second particle passed by the process function
  /// @param parts femtoUniverseParticles table (in case of Monte Carlo joined with FemtoUniverseMCLabels)
  /// @param magFieldTesla magnetic field of the collision
  /// @param multCol multiplicity of the collision
  template <bool isMC, typename PartitionType, typename PartType>
  void doSameEvent(PartitionType groupPartsHadron, PartitionType groupPartsPhi, PartType parts, float magFieldTesla, int multCol)
  {

    /// Histogramming same event
    for (auto& phicandidate : groupPartsPhi) {
      trackHistoPartPhi.fillQA<isMC, false>(phicandidate);
    }

    if (!ConfIsSame) {
      for (auto& hadron : groupPartsHadron) {
        // if (hadron.p() > twotracksconfigs.ConfCutTable->get("Hadron", "MaxP") || hadron.pt() > twotracksconfigs.ConfCutTable->get("Hadron", "MaxPt")) {
        //   continue;
        // }
        // if (!isFullPIDSelected(hadron.pidcut(),
        //                        hadron.p(),
        //                        twotracksconfigs.ConfCutTable->get("Hadron", "PIDthr"),
        //                        vPIDHadron,
        //                        twotracksconfigs.ConfNspecies,
        //                        kNsigma,
        //                        twotracksconfigs.ConfCutTable->get("Hadron", "nSigmaTPC"),
        //                        twotracksconfigs.ConfCutTable->get("Hadron", "nSigmaTPCTOF"))) {
        //   continue;
        // }
        if (!IsParticleNSigma(hadron.p(), trackCuts.getNsigmaTPC(hadron, o2::track::PID::Proton), trackCuts.getNsigmaTOF(hadron, o2::track::PID::Proton), trackCuts.getNsigmaTPC(hadron, o2::track::PID::Pion), trackCuts.getNsigmaTOF(hadron, o2::track::PID::Pion), trackCuts.getNsigmaTPC(hadron, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(hadron, o2::track::PID::Kaon))) {
          continue;
        }
        trackHistoPartHadron.fillQA<isMC, false>(hadron);
      }
    }
    /// Now build the combinations
    for (auto& [hadron, phicandidate] : combinations(CombinationsFullIndexPolicy(groupPartsHadron, groupPartsPhi))) {
      // if (hadron.p() > twotracksconfigs.ConfCutTable->get("PhiCandidate", "MaxP") || hadron.pt() > twotracksconfigs.ConfCutTable->get("PhiCandidate", "MaxPt") || phicandidate.p() > twotracksconfigs.ConfCutTable->get("Hadron", "MaxP") || phicandidate.pt() > twotracksconfigs.ConfCutTable->get("Hadron", "MaxPt")) {
      //   continue;
      // }
      // if (!isFullPIDSelected(hadron.pidcut(),
      //                        hadron.p(),
      //                        twotracksconfigs.ConfCutTable->get("PhiCandidate", "PIDthr"),
      //                        vPIDPhiCandidate,
      //                        twotracksconfigs.ConfNspecies,
      //                        kNsigma,
      //                        twotracksconfigs.ConfCutTable->get("PhiCandidate", "nSigmaTPC"),
      //                        twotracksconfigs.ConfCutTable->get("PhiCandidate", "nSigmaTPCTOF")) ||
      //     !isFullPIDSelected(phicandidate.pidcut(),
      //                        phicandidate.p(),
      //                        twotracksconfigs.ConfCutTable->get("Hadron", "PIDthr"),
      //                        vPIDHadron,
      //                        twotracksconfigs.ConfNspecies,
      //                        kNsigma,
      //                        twotracksconfigs.ConfCutTable->get("Hadron", "nSigmaTPC"),
      //                        twotracksconfigs.ConfCutTable->get("Hadron", "nSigmaTPCTOF"))) {
      //   continue;
      // }

      if (!IsParticleNSigma(hadron.p(), trackCuts.getNsigmaTPC(hadron, o2::track::PID::Proton), trackCuts.getNsigmaTOF(hadron, o2::track::PID::Proton), trackCuts.getNsigmaTPC(hadron, o2::track::PID::Pion), trackCuts.getNsigmaTOF(hadron, o2::track::PID::Pion), trackCuts.getNsigmaTPC(hadron, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(hadron, o2::track::PID::Kaon))) {
        continue;
      }

      // Close Pair Rejection
      // if (ConfIsCPR.value) {
      if (pairCloseRejection.isClosePair(hadron, phicandidate, parts, magFieldTesla)) {
        continue;
      }
      //}

      // Track Cleaning
      if (!pairCleaner.isCleanPair(hadron, phicandidate, parts)) {
        continue;
      }
      sameEventFemtoCont.setPair<isMC>(hadron, phicandidate, multCol, twotracksconfigs.ConfUse3D);
      sameEventAngularCont.setPair<isMC>(hadron, phicandidate, multCol, twotracksconfigs.ConfUse3D);
    }
  }

  /// process function for to call doSameEvent with Data
  /// \param col subscribe to the collision table (Data)
  /// \param parts subscribe to the femtoUniverseParticleTable
  void processSameEvent(o2::aod::FDCollision& col,
                        FemtoFullParticles& parts)
  {
    fillCollision(col);

    auto thegroupPartsHadron = partsHadron->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto thegroupPartsPhi = partsPhi->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);

    doSameEvent<false>(thegroupPartsHadron, thegroupPartsPhi, parts, col.magField(), col.multNtr());
  }
  PROCESS_SWITCH(femtoUniversePairTaskTrackPhi, processSameEvent, "Enable processing same event", true);

  /// process function for to call doSameEvent with Monte Carlo
  /// \param col subscribe to the collision table (Monte Carlo Reconstructed reconstructed)
  /// \param parts subscribe to joined table FemtoUniverseParticles and FemtoUniverseMCLables to access Monte Carlo truth
  /// \param FemtoUniverseMCParticles subscribe to the Monte Carlo truth table
  void processSameEventMC(o2::aod::FDCollision& col,
                          soa::Join<o2::aod::FDParticles, o2::aod::FDMCLabels>& parts,
                          o2::aod::FDMCParticles&)
  {
    fillCollision(col);

    auto thegroupPartsPhi = partsPhiMC->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto thegroupPartsHadron = partsHadronMC->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);

    doSameEvent<true>(thegroupPartsHadron, thegroupPartsPhi, parts, col.magField(), col.multNtr());
  }
  PROCESS_SWITCH(femtoUniversePairTaskTrackPhi, processSameEventMC, "Enable processing same event for Monte Carlo", false);

  /// This function processes the mixed event
  /// \todo the trivial loops over the collisions and tracks should be factored out since they will be common to all combinations of T-T, T-V0, V0-V0, ...
  /// \tparam PartitionType
  /// \tparam PartType
  /// \tparam isMC: enables Monte Carlo truth specific histograms
  /// \param groupPartsHadron partition for the identified passed by the process function
  /// \param groupPartsPhi partition for Phi meson passed by the process function
  /// \param parts femtoUniverseParticles table (in case of Monte Carlo joined with FemtoUniverseMCLabels)
  /// \param magFieldTesla magnetic field of the collision
  /// \param multCol multiplicity of the collision
  template <bool isMC, typename PartitionType, typename PartType>
  void doMixedEvent(PartitionType groupPartsHadron, PartitionType groupPartsPhi, PartType parts, float magFieldTesla, int multCol)
  {

    for (auto& [hadron, phicandidate] : combinations(CombinationsFullIndexPolicy(groupPartsHadron, groupPartsPhi))) {
      // if (hadron.p() > twotracksconfigs.ConfCutTable->get("PhiCandidate", "MaxP") || hadron.pt() > twotracksconfigs.ConfCutTable->get("PhiCandidate", "MaxPt") || phicandidate.p() > twotracksconfigs.ConfCutTable->get("Hadron", "MaxP") || phicandidate.pt() > twotracksconfigs.ConfCutTable->get("Hadron", "MaxPt")) {
      //   continue;
      // }
      // if (!isFullPIDSelected(hadron.pidcut(),
      //                        hadron.p(),
      //                        twotracksconfigs.ConfCutTable->get("PhiCandidate", "PIDthr"),
      //                        vPIDPhiCandidate,
      //                        twotracksconfigs.ConfNspecies,
      //                        kNsigma,
      //                        twotracksconfigs.ConfCutTable->get("PhiCandidate", "nSigmaTPC"),
      //                        twotracksconfigs.ConfCutTable->get("PhiCandidate", "nSigmaTPCTOF")) ||
      //     !isFullPIDSelected(phicandidate.pidcut(),
      //                        phicandidate.p(),
      //                        twotracksconfigs.ConfCutTable->get("Hadron", "PIDthr"),
      //                        vPIDHadron,
      //                        twotracksconfigs.ConfNspecies,
      //                        kNsigma,
      //                        twotracksconfigs.ConfCutTable->get("Hadron", "nSigmaTPC"),
      //                        twotracksconfigs.ConfCutTable->get("Hadron", "nSigmaTPCTOF"))) {
      //   continue;
      // }

      if (!IsParticleNSigma(hadron.p(), trackCuts.getNsigmaTPC(hadron, o2::track::PID::Proton), trackCuts.getNsigmaTOF(hadron, o2::track::PID::Proton), trackCuts.getNsigmaTPC(hadron, o2::track::PID::Pion), trackCuts.getNsigmaTOF(hadron, o2::track::PID::Pion), trackCuts.getNsigmaTPC(hadron, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(hadron, o2::track::PID::Kaon))) {
        continue;
      }

      // // if (ConfIsCPR.value) {
      if (pairCloseRejection.isClosePair(hadron, phicandidate, parts, magFieldTesla)) {
        continue;
      }
      // }

      mixedEventFemtoCont.setPair<isMC>(hadron, phicandidate, multCol, twotracksconfigs.ConfUse3D);
      mixedEventAngularCont.setPair<isMC>(hadron, phicandidate, multCol, twotracksconfigs.ConfUse3D);
    }
  }

  /// process function for to call doMixedEvent with Data
  /// @param cols subscribe to the collisions table (Data)
  /// @param parts subscribe to the femtoUniverseParticleTable
  void processMixedEvent(o2::aod::FDCollisions& cols,
                         FemtoFullParticles& parts)
  {
    for (auto& [collision1, collision2] : soa::selfCombinations(colBinning, 5, -1, cols, cols)) {

      const int multiplicityCol = collision1.multNtr();
      MixQaRegistry.fill(HIST("MixingQA/hMECollisionBins"), colBinning.getBin({collision1.posZ(), multiplicityCol}));

      auto groupPartsHadron = partsHadron->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision2.globalIndex(), cache);
      auto groupPartsPhi = partsPhi->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);

      const auto& magFieldTesla1 = collision1.magField();
      const auto& magFieldTesla2 = collision2.magField();

      if (magFieldTesla1 != magFieldTesla2) {
        continue;
      }
      /// \todo before mixing we should check whether both collisions contain a pair of particles!
      // if (partsPhi.size() == 0 || nPart2Evt1 == 0 || nPart1Evt2 == 0 || partsHadron.size() == 0 ) continue;

      doMixedEvent<false>(groupPartsHadron, groupPartsPhi, parts, magFieldTesla1, multiplicityCol);
    }
  }
  PROCESS_SWITCH(femtoUniversePairTaskTrackPhi, processMixedEvent, "Enable processing mixed events", true);

  /// brief process function for to call doMixedEvent with Monte Carlo
  /// @param cols subscribe to the collisions table (Monte Carlo Reconstructed reconstructed)
  /// @param parts subscribe to joined table FemtoUniverseParticles and FemtoUniverseMCLables to access Monte Carlo truth
  /// @param FemtoUniverseMCParticles subscribe to the Monte Carlo truth table
  void processMixedEventMC(o2::aod::FDCollisions& cols,
                           soa::Join<o2::aod::FDParticles, o2::aod::FDMCLabels>& parts,
                           o2::aod::FDMCParticles&)
  {
    for (auto& [collision1, collision2] : soa::selfCombinations(colBinning, 5, -1, cols, cols)) {

      const int multiplicityCol = collision1.multNtr();
      MixQaRegistry.fill(HIST("MixingQA/hMECollisionBins"), colBinning.getBin({collision1.posZ(), multiplicityCol}));

      auto groupPartsHadron = partsHadronMC->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision2.globalIndex(), cache);
      auto groupPartsPhi = partsPhiMC->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);

      const auto& magFieldTesla1 = collision1.magField();
      const auto& magFieldTesla2 = collision2.magField();

      if (magFieldTesla1 != magFieldTesla2) {
        continue;
      }
      /// \todo before mixing we should check whether both collisions contain a pair of particles!
      // if (partsPhi.size() == 0 || nPart2Evt1 == 0 || nPart1Evt2 == 0 || partsHadron.size() == 0 ) continue;

      doMixedEvent<true>(groupPartsHadron, groupPartsPhi, parts, magFieldTesla1, multiplicityCol);
    }
  }
  PROCESS_SWITCH(femtoUniversePairTaskTrackPhi, processMixedEventMC, "Enable processing mixed events MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<femtoUniversePairTaskTrackPhi>(cfgc),
  };
  return workflow;
}
