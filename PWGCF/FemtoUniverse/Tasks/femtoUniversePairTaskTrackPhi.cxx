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
static const std::vector<std::string> partNames{"PhiCandidate", "Track"};
static const std::vector<std::string> cutNames{"MaxPt", "PIDthr", "nSigmaTPC", "nSigmaTPCTOF", "MaxP"};
static const float cutsTable[nPart][nCuts]{
  {4.05f, 1.f, 3.f, 3.f, 100.f},
  {4.05f, 1.f, 3.f, 3.f, 100.f}};
} // namespace

struct femtoUniversePairTaskTrackPhi {

  using FemtoFullParticles = soa::Join<aod::FDParticles, aod::FDExtParticles>;
  SliceCache cache;
  Preslice<FemtoFullParticles> perCol = aod::femtouniverseparticle::fdCollisionId;

  /// Table for both particles
  struct : o2::framework::ConfigurableGroup {
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
  } ConfBothTracks;

  /// Particle 1 --- IDENTIFIED TRACK
  struct : o2::framework::ConfigurableGroup {
    Configurable<bool> ConfIsSame{"ConfIsSame", false, "Pairs of the same particle"};
    Configurable<int> ConfPDGCodeTrack{"ConfPDGCodeTrack", 2212, "Particle 2 - PDG code"};
    // Configurable<uint32_t> ConfCutTrack{"ConfCutTrack", 5542474, "Particle 2 - Selection bit"};
    Configurable<int> ConfPIDTrack{"ConfPIDTrack", 2, "Particle 2 - Read from cutCulator"}; // we also need the possibility to specify whether the bit is true/false ->std>>vector<std::pair<int, int>>
    Configurable<int8_t> ConfTrackSign{"ConfTrackSign", 1, "Track sign"};
    Configurable<bool> ConfIsTrackIdentified{"ConfIsTrackIdentified", true, "Enable PID for the track"};
  } ConfTrack;

  /// Particle 2 --- PHI
  struct : o2::framework::ConfigurableGroup {
    Configurable<int> ConfPDGCodePhi{"ConfPDGCodePhi", 333, "Phi meson - PDG code"};
    // Configurable<uint32_t> ConfCutPhi{"ConfCutPhi", 5542474, "Phi meson - Selection bit from cutCulator"};
    Configurable<int> ConfPIDPhi{"ConfPIDPhi", 2, "Phi meson - Read from cutCulator"}; // we also need the possibility to specify whether the bit is true/false ->std>>vector<std::pair<int, int>>int>>
    // Configurable<float> ConfPtLowPhi{"ConfPtLowPhi", 0.5, "Lower limit for Pt for the Phi meson"}; // change according to wrzesa cuts
    // Configurable<float> ConfPtHighPhi{"ConfPtHighPhi", 1.5, "Higher limit for Pt for the Phi meson"};
  } ConfPhi;

  /// Partitions for particle 1
  Partition<FemtoFullParticles> partsTrack = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kTrack)) && (aod::femtouniverseparticle::sign == int8_t(ConfTrack.ConfTrackSign));
  Partition<soa::Join<aod::FDParticles, aod::FDMCLabels>> partsTrackMC = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kTrack));

  /// Partitions for particle 2
  Partition<FemtoFullParticles> partsPhi = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kPhi));
  Partition<soa::Join<aod::FDParticles, aod::FDMCLabels>> partsPhiMC = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kPhi));

  /// Histogramming for particle 1
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kTrack, 2> trackHistoPartTrack;

  /// Histogramming for particle 2
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kPhi, 0> trackHistoPartPhi;

  /// Histogramming for Event
  FemtoUniverseEventHisto eventHisto;

  /// The configurables need to be passed to an std::vector
  int vPIDPhiCandidate, vPIDTrack;
  std::vector<float> kNsigma;

  /// particle part
  ConfigurableAxis ConfTempFitVarBins{"ConfDTempFitVarBins", {300, -0.15, 0.15}, "binning of the TempFitVar in the pT vs. TempFitVar plot"};
  ConfigurableAxis ConfTempFitVarInvMassBins{"ConfDTempFitVarInvMassBins", {6000, 0.9, 4.0}, "binning of the TempFitVar in the pT vs. TempFitVar plot"};
  ConfigurableAxis ConfTempFitVarpTBins{"ConfTempFitVarpTBins", {20, 0.5, 4.05}, "pT binning of the pT vs. TempFitVar plot"};

  /// Correlation part
  ConfigurableAxis ConfMultBins{"ConfMultBins", {VARIABLE_WIDTH, 0.0f, 4.0f, 8.0f, 12.0f, 16.0f, 20.0f, 24.0f, 28.0f, 32.0f, 36.0f, 40.0f, 44.0f, 48.0f, 52.0f, 56.0f, 60.0f, 64.0f, 68.0f, 72.0f, 76.0f, 80.0f, 84.0f, 88.0f, 92.0f, 96.0f, 100.0f, 200.0f, 99999.f}, "Mixing bins - multiplicity"}; // \todo to be obtained from the hash task
  // ConfigurableAxis ConfMultBins{"CfgMultBins", {VARIABLE_WIDTH, 0.0f, 20.0f, 40.0f, 60.0f, 80.0f, 100.0f, 200.0f, 99999.f}, "Mixing bins - multiplicity"};
  ConfigurableAxis ConfVtxBins{"ConfVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis ConfmTBins3D{"ConfmTBins3D", {VARIABLE_WIDTH, 1.02f, 1.14f, 1.20f, 1.26f, 1.38f, 1.56f, 1.86f, 4.50f}, "mT Binning for the 3Dimensional plot: k* vs multiplicity vs mT (set <<ConfBothTracks.ConfUse3D>> to true in order to use)"};
  ConfigurableAxis ConfmultBins3D{"ConfmultBins3D", {VARIABLE_WIDTH, 0.0f, 20.0f, 30.0f, 40.0f, 99999.0f}, "multiplicity Binning for the 3Dimensional plot: k* vs multiplicity vs mT (set <<ConfBothTracks.ConfUse3D>> to true in order to use)"};

  ColumnBinningPolicy<aod::collision::PosZ, aod::femtouniversecollision::MultNtr> colBinning{{ConfVtxBins, ConfMultBins}, true};

  ConfigurableAxis ConfkstarBins{"ConfkstarBins", {1500, 0., 6.}, "binning kstar"};
  ConfigurableAxis ConfkTBins{"ConfkTBins", {150, 0., 9.}, "binning kT"};
  ConfigurableAxis ConfmTBins{"ConfmTBins", {225, 0., 7.5}, "binning mT"};
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
      if (TMath::Abs(nsigmaTPCPr) < ConfBothTracks.ConfNsigmaTPCProton) {
        return true;
      } else {
        return false;
      }
    } else if (mom > 0.4) {
      if (TMath::Hypot(nsigmaTOFPr, nsigmaTPCPr) < ConfBothTracks.ConfNsigmaCombinedProton) {
        // if (TMath::Abs(nsigmaTPCPr) < ConfBothTracks.ConfNsigmaCombinedProton) {
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
        if (TMath::Abs(nsigmaTPCPi) < ConfBothTracks.ConfNsigmaTPCPion) {
          return true;
        } else {
          return false;
        }
      } else if (mom > 0.5) {
        if (TMath::Hypot(nsigmaTOFPi, nsigmaTPCPi) < ConfBothTracks.ConfNsigmaCombinedPion) {
          // if (TMath::Abs(nsigmaTPCPi) < ConfBothTracks.ConfNsigmaCombinedPion) {
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
    switch (ConfTrack.ConfPDGCodeTrack) {
      case 2212:  // Proton
      case -2212: // anty Proton
        return IsProtonNSigma(mom, nsigmaTPCPr, nsigmaTOFPr);
        break;
      case 211:  // Pion
      case -211: // Pion-
        return IsPionNSigma(mom, nsigmaTPCPi, nsigmaTOFPi);
        break;
      case 321:  // Kaon+
      case -321: // Kaon-
        return IsKaonNSigma(mom, nsigmaTPCK, nsigmaTOFK);
        break;
      default:
        return false;
    }
  }

  void init(InitContext&)
  {
    eventHisto.init(&qaRegistry);
    trackHistoPartPhi.init(&qaRegistry, ConfTempFitVarpTBins, ConfTempFitVarInvMassBins, ConfBothTracks.ConfIsMC, ConfPhi.ConfPDGCodePhi);
    if (!ConfTrack.ConfIsSame) {
      trackHistoPartTrack.init(&qaRegistry, ConfTempFitVarpTBins, ConfTempFitVarBins, ConfBothTracks.ConfIsMC, ConfTrack.ConfPDGCodeTrack);
    }

    MixQaRegistry.add("MixingQA/hSECollisionBins", ";bin;Entries", kTH1F, {{120, -0.5, 119.5}});
    MixQaRegistry.add("MixingQA/hMECollisionBins", ";bin;Entries", kTH1F, {{120, -0.5, 119.5}});

    sameEventFemtoCont.init(&resultRegistry, ConfkstarBins, ConfMultBins, ConfkTBins, ConfmTBins, ConfmultBins3D, ConfmTBins3D, ConfBothTracks.ConfIsMC, ConfBothTracks.ConfUse3D);
    mixedEventFemtoCont.init(&resultRegistry, ConfkstarBins, ConfMultBins, ConfkTBins, ConfmTBins, ConfmultBins3D, ConfmTBins3D, ConfBothTracks.ConfIsMC, ConfBothTracks.ConfUse3D);
    sameEventAngularCont.init(&resultRegistry, ConfkstarBins, ConfMultBins, ConfkTBins, ConfmTBins, ConfmultBins3D, ConfmTBins3D, ConfBothTracks.ConfEtaBins, ConfBothTracks.ConfPhiBins, ConfBothTracks.ConfIsMC, ConfBothTracks.ConfUse3D);
    mixedEventAngularCont.init(&resultRegistry, ConfkstarBins, ConfMultBins, ConfkTBins, ConfmTBins, ConfmultBins3D, ConfmTBins3D, ConfBothTracks.ConfEtaBins, ConfBothTracks.ConfPhiBins, ConfBothTracks.ConfIsMC, ConfBothTracks.ConfUse3D);

    sameEventFemtoCont.setPDGCodes(ConfPhi.ConfPDGCodePhi, ConfTrack.ConfPDGCodeTrack);
    mixedEventFemtoCont.setPDGCodes(ConfPhi.ConfPDGCodePhi, ConfTrack.ConfPDGCodeTrack);
    sameEventAngularCont.setPDGCodes(ConfPhi.ConfPDGCodePhi, ConfTrack.ConfPDGCodeTrack);
    mixedEventAngularCont.setPDGCodes(ConfPhi.ConfPDGCodePhi, ConfTrack.ConfPDGCodeTrack);

    pairCleaner.init(&qaRegistry);
    if (ConfIsCPR.value) {
      pairCloseRejection.init(&resultRegistry, &qaRegistry, ConfCPRdeltaPhiMax.value, ConfCPRdeltaEtaMax.value, ConfCPRPlotPerRadii.value);
    }

    vPIDPhiCandidate = ConfPhi.ConfPIDPhi.value;
    vPIDTrack = ConfTrack.ConfPIDTrack.value;
    kNsigma = ConfBothTracks.ConfTrkPIDnSigmaMax.value;
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
  /// @param groupPartsTrack partition for the first particle passed by the process function
  /// @param groupPartsPhi partition for the second particle passed by the process function
  /// @param parts femtoUniverseParticles table (in case of Monte Carlo joined with FemtoUniverseMCLabels)
  /// @param magFieldTesla magnetic field of the collision
  /// @param multCol multiplicity of the collision
  template <bool isMC, typename PartitionType, typename PartType>
  void doSameEvent(PartitionType groupPartsTrack, PartitionType groupPartsPhi, PartType parts, float magFieldTesla, int multCol)
  {

    /// Histogramming same event
    for (auto& phicandidate : groupPartsPhi) {
      trackHistoPartPhi.fillQA<isMC, false>(phicandidate);
    }

    if (!ConfTrack.ConfIsSame) {
      for (auto& track : groupPartsTrack) {
        // if (track.p() > ConfBothTracks.ConfCutTable->get("Track", "MaxP") || track.pt() > ConfBothTracks.ConfCutTable->get("Track", "MaxPt")) {
        //   continue;
        // }
        // if (!isFullPIDSelected(track.pidcut(),
        //                        track.p(),
        //                        ConfBothTracks.ConfCutTable->get("Track", "PIDthr"),
        //                        vPIDTrack,
        //                        ConfBothTracks.ConfNspecies,
        //                        kNsigma,
        //                        ConfBothTracks.ConfCutTable->get("Track", "nSigmaTPC"),
        //                        ConfBothTracks.ConfCutTable->get("Track", "nSigmaTPCTOF"))) {
        //   continue;
        // }
        if (ConfTrack.ConfIsTrackIdentified) {
          if (!IsParticleNSigma(track.p(), trackCuts.getNsigmaTPC(track, o2::track::PID::Proton), trackCuts.getNsigmaTOF(track, o2::track::PID::Proton), trackCuts.getNsigmaTPC(track, o2::track::PID::Pion), trackCuts.getNsigmaTOF(track, o2::track::PID::Pion), trackCuts.getNsigmaTPC(track, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(track, o2::track::PID::Kaon))) {
            continue;
          }
        }
        trackHistoPartTrack.fillQA<isMC, false>(track);
      }
    }
    /// Now build the combinations
    for (auto& [track, phicandidate] : combinations(CombinationsFullIndexPolicy(groupPartsTrack, groupPartsPhi))) {
      // if (track.p() > ConfBothTracks.ConfCutTable->get("PhiCandidate", "MaxP") || track.pt() > ConfBothTracks.ConfCutTable->get("PhiCandidate", "MaxPt") || phicandidate.p() > ConfBothTracks.ConfCutTable->get("Track", "MaxP") || phicandidate.pt() > ConfBothTracks.ConfCutTable->get("Track", "MaxPt")) {
      //   continue;
      // }
      // if (!isFullPIDSelected(track.pidcut(),
      //                        track.p(),
      //                        ConfBothTracks.ConfCutTable->get("PhiCandidate", "PIDthr"),
      //                        vPIDPhiCandidate,
      //                        ConfBothTracks.ConfNspecies,
      //                        kNsigma,
      //                        ConfBothTracks.ConfCutTable->get("PhiCandidate", "nSigmaTPC"),
      //                        ConfBothTracks.ConfCutTable->get("PhiCandidate", "nSigmaTPCTOF")) ||
      //     !isFullPIDSelected(phicandidate.pidcut(),
      //                        phicandidate.p(),
      //                        ConfBothTracks.ConfCutTable->get("Track", "PIDthr"),
      //                        vPIDTrack,
      //                        ConfBothTracks.ConfNspecies,
      //                        kNsigma,
      //                        ConfBothTracks.ConfCutTable->get("Track", "nSigmaTPC"),
      //                        ConfBothTracks.ConfCutTable->get("Track", "nSigmaTPCTOF"))) {
      //   continue;
      // }
      if (ConfTrack.ConfIsTrackIdentified) {
        if (!IsParticleNSigma(track.p(), trackCuts.getNsigmaTPC(track, o2::track::PID::Proton), trackCuts.getNsigmaTOF(track, o2::track::PID::Proton), trackCuts.getNsigmaTPC(track, o2::track::PID::Pion), trackCuts.getNsigmaTOF(track, o2::track::PID::Pion), trackCuts.getNsigmaTPC(track, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(track, o2::track::PID::Kaon))) {
          continue;
        }
      }
      // // Close Pair Rejection
      if (ConfIsCPR.value) {
        if (pairCloseRejection.isClosePair(track, phicandidate, parts, magFieldTesla)) {
          continue;
        }
      }

      // Track Cleaning
      if (!pairCleaner.isCleanPair(track, phicandidate, parts)) {
        continue;
      }
      sameEventFemtoCont.setPair<isMC>(track, phicandidate, multCol, ConfBothTracks.ConfUse3D);
      sameEventAngularCont.setPair<isMC>(track, phicandidate, multCol, ConfBothTracks.ConfUse3D);
    }
  }

  /// process function for to call doSameEvent with Data
  /// \param col subscribe to the collision table (Data)
  /// \param parts subscribe to the femtoUniverseParticleTable
  void processSameEvent(o2::aod::FDCollision& col,
                        FemtoFullParticles& parts)
  {
    fillCollision(col);

    auto thegroupPartsTrack = partsTrack->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto thegroupPartsPhi = partsPhi->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);

    doSameEvent<false>(thegroupPartsTrack, thegroupPartsPhi, parts, col.magField(), col.multNtr());
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
    auto thegroupPartsTrack = partsTrackMC->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);

    doSameEvent<true>(thegroupPartsTrack, thegroupPartsPhi, parts, col.magField(), col.multNtr());
  }
  PROCESS_SWITCH(femtoUniversePairTaskTrackPhi, processSameEventMC, "Enable processing same event for Monte Carlo", false);

  /// This function processes the mixed event
  /// \todo the trivial loops over the collisions and tracks should be factored out since they will be common to all combinations of T-T, T-V0, V0-V0, ...
  /// \tparam PartitionType
  /// \tparam PartType
  /// \tparam isMC: enables Monte Carlo truth specific histograms
  /// \param groupPartsTrack partition for the identified passed by the process function
  /// \param groupPartsPhi partition for Phi meson passed by the process function
  /// \param parts femtoUniverseParticles table (in case of Monte Carlo joined with FemtoUniverseMCLabels)
  /// \param magFieldTesla magnetic field of the collision
  /// \param multCol multiplicity of the collision
  template <bool isMC, typename PartitionType, typename PartType>
  void doMixedEvent(PartitionType groupPartsTrack, PartitionType groupPartsPhi, PartType parts, float magFieldTesla, int multCol)
  {

    for (auto& [track, phicandidate] : combinations(CombinationsFullIndexPolicy(groupPartsTrack, groupPartsPhi))) {
      // if (track.p() > ConfBothTracks.ConfCutTable->get("PhiCandidate", "MaxP") || track.pt() > ConfBothTracks.ConfCutTable->get("PhiCandidate", "MaxPt") || phicandidate.p() > ConfBothTracks.ConfCutTable->get("Track", "MaxP") || phicandidate.pt() > ConfBothTracks.ConfCutTable->get("Track", "MaxPt")) {
      //   continue;
      // }
      // if (!isFullPIDSelected(track.pidcut(),
      //                        track.p(),
      //                        ConfBothTracks.ConfCutTable->get("PhiCandidate", "PIDthr"),
      //                        vPIDPhiCandidate,
      //                        ConfBothTracks.ConfNspecies,
      //                        kNsigma,
      //                        ConfBothTracks.ConfCutTable->get("PhiCandidate", "nSigmaTPC"),
      //                        ConfBothTracks.ConfCutTable->get("PhiCandidate", "nSigmaTPCTOF")) ||
      //     !isFullPIDSelected(phicandidate.pidcut(),
      //                        phicandidate.p(),
      //                        ConfBothTracks.ConfCutTable->get("Track", "PIDthr"),
      //                        vPIDTrack,
      //                        ConfBothTracks.ConfNspecies,
      //                        kNsigma,
      //                        ConfBothTracks.ConfCutTable->get("Track", "nSigmaTPC"),
      //                        ConfBothTracks.ConfCutTable->get("Track", "nSigmaTPCTOF"))) {
      //   continue;
      // }
      if (ConfTrack.ConfIsTrackIdentified) {
        if (!IsParticleNSigma(track.p(), trackCuts.getNsigmaTPC(track, o2::track::PID::Proton), trackCuts.getNsigmaTOF(track, o2::track::PID::Proton), trackCuts.getNsigmaTPC(track, o2::track::PID::Pion), trackCuts.getNsigmaTOF(track, o2::track::PID::Pion), trackCuts.getNsigmaTPC(track, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(track, o2::track::PID::Kaon))) {
          continue;
        }
      }
      if (ConfIsCPR.value) {
        if (pairCloseRejection.isClosePair(track, phicandidate, parts, magFieldTesla)) {
          continue;
        }
      }

      mixedEventFemtoCont.setPair<isMC>(track, phicandidate, multCol, ConfBothTracks.ConfUse3D);
      mixedEventAngularCont.setPair<isMC>(track, phicandidate, multCol, ConfBothTracks.ConfUse3D);
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

      auto groupPartsTrack = partsTrack->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision2.globalIndex(), cache);
      auto groupPartsPhi = partsPhi->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);

      const auto& magFieldTesla1 = collision1.magField();
      const auto& magFieldTesla2 = collision2.magField();

      if (magFieldTesla1 != magFieldTesla2) {
        continue;
      }
      /// \todo before mixing we should check whether both collisions contain a pair of particles!
      // if (partsPhi.size() == 0 || nPart2Evt1 == 0 || nPart1Evt2 == 0 || partsTrack.size() == 0 ) continue;

      doMixedEvent<false>(groupPartsTrack, groupPartsPhi, parts, magFieldTesla1, multiplicityCol);
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

      auto groupPartsTrack = partsTrackMC->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision2.globalIndex(), cache);
      auto groupPartsPhi = partsPhiMC->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);

      const auto& magFieldTesla1 = collision1.magField();
      const auto& magFieldTesla2 = collision2.magField();

      if (magFieldTesla1 != magFieldTesla2) {
        continue;
      }
      /// \todo before mixing we should check whether both collisions contain a pair of particles!
      // if (partsPhi.size() == 0 || nPart2Evt1 == 0 || nPart1Evt2 == 0 || partsTrack.size() == 0 ) continue;

      doMixedEvent<true>(groupPartsTrack, groupPartsPhi, parts, magFieldTesla1, multiplicityCol);
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
