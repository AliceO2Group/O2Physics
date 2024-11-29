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
/// \author Zuzanna Chochulska, WUT Warsaw & CTU Prague, zchochul@cern.ch

#include <vector>
#include <string>
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
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseAngularContainer.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseDetaDphiStar.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUtils.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseTrackSelection.h"
#include <TFile.h>
#include <TH1.h>
#include "CCDB/BasicCCDBManager.h"

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

  Service<o2::framework::O2DatabasePDG> pdgMC;

  using FemtoFullParticles = soa::Join<aod::FDParticles, aod::FDExtParticles>;
  SliceCache cache;
  Preslice<FemtoFullParticles> perCol = aod::femtouniverseparticle::fdCollisionId;

  using FemtoRecoParticles = soa::Join<aod::FDParticles, aod::FDExtParticles, aod::FDMCLabels>;
  Preslice<FemtoRecoParticles> perColMC = aod::femtouniverseparticle::fdCollisionId;

  // Efficiency
  struct : o2::framework::ConfigurableGroup {
    Configurable<std::string> ConfEfficiencyTrackPath{"ConfEfficiencyTrackPath", "", "Local path to proton efficiency TH2F file"};
    Configurable<std::string> ConfEfficiencyPhiPath{"ConfEfficiencyPhiPath", "", "Local path to Phi efficiency TH2F file"};
    Configurable<int64_t> ConfEfficiencyTrackTimestamp{"ConfEfficiencyTrackTimestamp", 0, "(int64_t) Timestamp for hadron"};
    Configurable<int64_t> ConfEfficiencyPhiTimestamp{"ConfEfficiencyPhiTimestamp", 0, "(int64_t) Timestamp for phi"};
  } ConfEff;

  struct : o2::framework::ConfigurableGroup {
    Configurable<bool> ConfCPRIsEnabled{"ConfCPRIsEnabled", true, "Close Pair Rejection"};
    Configurable<bool> ConfCPRPlotPerRadii{"ConfCPRPlotPerRadii", false, "Plot CPR per radii"};
    Configurable<float> ConfCPRdeltaPhiCutMax{"ConfCPRdeltaPhiCutMax", 0.0, "Delta Phi max cut for Close Pair Rejection"};
    Configurable<float> ConfCPRdeltaPhiCutMin{"ConfCPRdeltaPhiCutMin", 0.0, "Delta Phi min cut for Close Pair Rejection"};
    Configurable<float> ConfCPRdeltaEtaCutMax{"ConfCPRdeltaEtaCutMax", 0.0, "Delta Eta max cut for Close Pair Rejection"};
    Configurable<float> ConfCPRdeltaEtaCutMin{"ConfCPRdeltaEtaCutMin", 0.0, "Delta Eta min cut for Close Pair Rejection"};
    Configurable<float> ConfCPRInvMassCutMin{"ConfCPRInvMassCutMin", 1.014, "Invariant mass (low) cut for Close Pair Rejection"};
    Configurable<float> ConfCPRInvMassCutMax{"ConfCPRInvMassCutMax", 1.026, "Invariant mass (high) cut for Close Pair Rejection"};
    Configurable<float> ConfCPRChosenRadii{"ConfCPRChosenRadii", 0.80, "Delta Eta cut for Close Pair Rejection"};
  } ConfCPR;

  /// Table for both particles
  struct : o2::framework::ConfigurableGroup {
    Configurable<float> ConfPIDProtonNsigmaCombined{"ConfPIDProtonNsigmaCombined", 3.0, "TPC and TOF Proton Sigma (combined) for momentum > 0.5"};
    Configurable<float> ConfPIDProtonNsigmaTPC{"ConfPIDProtonNsigmaTPC", 3.0, "TPC Proton Sigma for momentum < 0.5"};
    Configurable<float> ConfPIDKaonNsigmaReject{"ConfPIDKaonNsigmaReject", 3.0, "Reject if particle could be a Kaon combined nsigma value."};
    Configurable<float> ConfPIDPionNsigmaReject{"ConfPIDPionNsigmaReject", 3.0, "Reject if particle could be a Pion combined nsigma value."};
    Configurable<float> ConfPIDProtonNsigmaReject{"ConfPIDProtonNsigmaReject", 3.0, "Reject if particle could be a Proton combined nsigma value."};
    Configurable<float> ConfPIDPionNsigmaCombined{"ConfPIDPionNsigmaCombined", 3.0, "TPC and TOF Pion Sigma (combined) for momentum > 0.5"};
    Configurable<float> ConfPIDPionNsigmaTPC{"ConfPIDPionNsigmaTPC", 3.0, "TPC Pion Sigma for momentum < 0.5"};

    Configurable<LabeledArray<float>> ConfCutTable{"ConfCutTable", {cutsTable[0], nPart, nCuts, partNames, cutNames}, "Particle selections"};
    Configurable<int> ConfNspecies{"ConfNspecies", 2, "Number of particle spieces with PID info"};
    Configurable<bool> ConfIsMC{"ConfIsMC", false, "Enable additional Histograms in the case of a MonteCarlo Run"};
    Configurable<std::vector<float>> ConfTrkPIDnSigmaMax{"ConfTrkPIDnSigmaMax", std::vector<float>{4.f, 3.f, 2.f}, "This configurable needs to be the same as the one used in the producer task"};
    Configurable<bool> ConfUse3D{"ConfUse3D", false, "Enable three dimensional histogramms (to be used only for analysis with high statistics): k* vs mT vs multiplicity"};
    Configurable<int> ConfBinsPhi{"ConfBinsPhi", 29, "Number of phi bins in deta dphi"};
    Configurable<int> ConfBinsEta{"ConfBinsEta", 29, "Number of eta bins in deta dphi"};
  } ConfBothTracks;

  /// Particle 1 --- IDENTIFIED TRACK
  struct : o2::framework::ConfigurableGroup {
    Configurable<bool> ConfTrackIsSame{"ConfTrackIsSame", false, "Pairs of the same particle"};
    Configurable<int> ConfTrackPDGCode{"ConfTrackPDGCode", 2212, "Particle 2 - PDG code"};
    Configurable<int> ConfTrackPID{"ConfTrackPID", 2, "Particle 2 - Read from cutCulator"}; // we also need the possibility to specify whether the bit is true/false ->std>>vector<std::pair<int, int>>
    Configurable<int> ConfTrackSign{"ConfTrackSign", 1, "Track sign"};
    Configurable<bool> ConfTrackIsIdentified{"ConfTrackIsIdentified", true, "Enable PID for the track"};
    Configurable<bool> ConfTrackIsRejected{"ConfTrackIsRejected", true, "Enable PID rejection for the track other species than the identified one."};
    Configurable<float> ConfTrackPtLowLimit{"ConfTrackPtLowLimit", 0.5, "Lower limit of the Phi pT."};
    Configurable<float> ConfTrackPtHighLimit{"ConfTrackPtHighLimit", 2.5, "Higher limit of the Phi pT."};
  } ConfTrack;

  /// Particle 2 --- PHI
  struct : o2::framework::ConfigurableGroup {
    Configurable<int> ConfPhiPDGCode{"ConfPhiPDGCode", 333, "Phi meson - PDG code"};
    Configurable<int> ConfPhiPID{"ConfPhiPID", 2, "Phi meson - Read from cutCulator"}; // we also need the possibility to specify whether the bit is true/false ->std>>vector<std::pair<int, int>>int>>
  } ConfPhi;

  /// Partitions for particle 1
  Partition<FemtoFullParticles> partsTrack = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kTrack)) && (aod::femtouniverseparticle::sign == ConfTrack.ConfTrackSign.value) && (aod::femtouniverseparticle::pt > ConfTrack.ConfTrackPtLowLimit.value) && (aod::femtouniverseparticle::pt < ConfTrack.ConfTrackPtHighLimit.value);
  Partition<soa::Join<aod::FDParticles, aod::FDMCLabels>> partsTrackMC = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kTrack)) && (aod::femtouniverseparticle::sign == ConfTrack.ConfTrackSign.value) && (aod::femtouniverseparticle::pt > ConfTrack.ConfTrackPtLowLimit.value) && (aod::femtouniverseparticle::pt < ConfTrack.ConfTrackPtHighLimit.value);

  /// Partitions for particle 2
  Partition<FemtoFullParticles> partsPhi = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kPhi));
  Partition<soa::Join<aod::FDParticles, aod::FDMCLabels>> partsPhiMC = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kPhi));

  /// Partitions  for Phi daughters kPhiChild
  Partition<FemtoFullParticles> partsPhiDaugh = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kPhiChild));
  Partition<soa::Join<aod::FDParticles, aod::FDMCLabels>> partsPhiDaughMC = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kPhiChild));

  // Partition for K+K- minv
  Partition<FemtoFullParticles> partsKaons = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kTrack));
  Partition<soa::Join<aod::FDParticles, aod::FDMCLabels>> partsKaonsMC = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kTrack));

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
  ConfigurableAxis ConfBinsTempFitVar{"ConfDTempFitVarBins", {300, -0.15, 0.15}, "binning of the TempFitVar in the pT vs. TempFitVar plot"};
  ConfigurableAxis ConfBinsTempFitVarInvMass{"ConfDTempFitVarInvMassBins", {6000, 0.9, 4.0}, "binning of the TempFitVar in the pT vs. TempFitVar plot"};
  ConfigurableAxis ConfBinsTempFitVarpT{"ConfBinsTempFitVarpT", {20, 0.5, 4.05}, "pT binning of the pT vs. TempFitVar plot"};

  /// Correlation part
  ConfigurableAxis ConfBinsMult{"ConfBinsMult", {VARIABLE_WIDTH, 0.0f, 4.0f, 8.0f, 12.0f, 16.0f, 20.0f, 24.0f, 28.0f, 32.0f, 36.0f, 40.0f, 44.0f, 48.0f, 52.0f, 56.0f, 60.0f, 64.0f, 68.0f, 72.0f, 76.0f, 80.0f, 84.0f, 88.0f, 92.0f, 96.0f, 100.0f, 200.0f, 99999.f}, "Mixing bins - multiplicity"}; // \todo to be obtained from the hash task
  ConfigurableAxis ConfBinsVtx{"ConfBinsVtx", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis ConfBins3DmT{"ConfBins3DmT", {VARIABLE_WIDTH, 1.02f, 1.14f, 1.20f, 1.26f, 1.38f, 1.56f, 1.86f, 4.50f}, "mT Binning for the 3Dimensional plot: k* vs multiplicity vs mT (set <<ConfBothTracks.ConfUse3D>> to true in order to use)"};
  ConfigurableAxis ConfBins3Dmult{"ConfBins3Dmult", {VARIABLE_WIDTH, 0.0f, 20.0f, 30.0f, 40.0f, 99999.0f}, "multiplicity Binning for the 3Dimensional plot: k* vs multiplicity vs mT (set <<ConfBothTracks.ConfUse3D>> to true in order to use)"};

  ColumnBinningPolicy<aod::collision::PosZ, aod::femtouniversecollision::MultNtr> colBinning{{ConfBinsVtx, ConfBinsMult}, true};

  ConfigurableAxis ConfBinskstar{"ConfBinskstar", {1500, 0., 6.}, "binning kstar"};
  ConfigurableAxis ConfBinskT{"ConfBinskT", {150, 0., 9.}, "binning kT"};
  ConfigurableAxis ConfBinsmT{"ConfBinsmT", {225, 0., 7.5}, "binning mT"};

  FemtoUniverseAngularContainer<femtoUniverseAngularContainer::EventType::same, femtoUniverseAngularContainer::Observable::kstar> sameEventAngularCont;
  FemtoUniverseAngularContainer<femtoUniverseAngularContainer::EventType::mixed, femtoUniverseAngularContainer::Observable::kstar> mixedEventAngularCont;
  FemtoUniversePairCleaner<aod::femtouniverseparticle::ParticleType::kTrack, aod::femtouniverseparticle::ParticleType::kPhi> pairCleaner;
  FemtoUniverseDetaDphiStar<aod::femtouniverseparticle::ParticleType::kTrack, aod::femtouniverseparticle::ParticleType::kPhi> pairCloseRejection;
  FemtoUniverseTrackSelection trackCuts;

  /// Histogram output
  HistogramRegistry qaRegistry{"TrackQA", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry resultRegistry{"Correlations", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry MixQaRegistry{"MixQaRegistry", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registryMCtruth{"MCtruthHistos", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry registryMCreco{"MCrecoHistos", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry registryPhiMinvBackground{"registryPhiMinvBackground", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};

  Service<ccdb::BasicCCDBManager> ccdb;
  TH2F* protoneff;
  TH2F* phieff;

  // PID for protons
  bool IsProtonNSigma(float mom, float nsigmaTPCPr, float nsigmaTOFPr) // previous version from: https://github.com/alisw/AliPhysics/blob/master/PWGCF/FEMTOSCOPY/AliFemtoUser/AliFemtoMJTrackCut.cxx
  {
    if (mom < 0.5) {
      if (TMath::Abs(nsigmaTPCPr) < ConfBothTracks.ConfPIDProtonNsigmaTPC.value) {
        return true;
      } else {
        return false;
      }
    } else if (mom > 0.4) {
      if (TMath::Hypot(nsigmaTOFPr, nsigmaTPCPr) < ConfBothTracks.ConfPIDProtonNsigmaCombined.value) {
        return true;
      } else {
        return false;
      }
    }
    return false;
  }

  bool IsProtonRejected(float mom, float nsigmaTPCPi, float nsigmaTOFPi, float nsigmaTPCK, float nsigmaTOFK)
  {
    if (mom < 0.5) {
      return true;
    }
    if (mom > 0.5) {
      if (TMath::Hypot(nsigmaTOFPi, nsigmaTPCPi) < ConfBothTracks.ConfPIDPionNsigmaReject.value) {
        return true;
      } else if (TMath::Hypot(nsigmaTOFK, nsigmaTPCK) < ConfBothTracks.ConfPIDKaonNsigmaReject.value) {
        return true;
      } else {
        return false;
      }
    } else {
      return false;
    }
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

  bool IsKaonRejected(float mom, float nsigmaTPCPr, float nsigmaTOFPr, float nsigmaTPCPi, float nsigmaTOFPi)
  {
    if (mom < 0.5) {
      if (TMath::Abs(nsigmaTPCPi) < ConfBothTracks.ConfPIDPionNsigmaReject.value) {
        return true;
      } else if (TMath::Abs(nsigmaTPCPr) < ConfBothTracks.ConfPIDProtonNsigmaReject.value) {
        return true;
      }
    }
    if (mom > 0.5) {
      if (TMath::Hypot(nsigmaTOFPi, nsigmaTPCPi) < ConfBothTracks.ConfPIDPionNsigmaReject.value) {
        return true;
      } else if (TMath::Hypot(nsigmaTOFPr, nsigmaTPCPr) < ConfBothTracks.ConfPIDProtonNsigmaReject.value) {
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
    if (true) {
      if (mom < 0.5) {
        if (TMath::Abs(nsigmaTPCPi) < ConfBothTracks.ConfPIDPionNsigmaTPC.value) {
          return true;
        } else {
          return false;
        }
      } else if (mom > 0.5) {
        if (TMath::Hypot(nsigmaTOFPi, nsigmaTPCPi) < ConfBothTracks.ConfPIDPionNsigmaCombined.value) {
          return true;
        } else {
          return false;
        }
      }
    }
    return false;
  }

  bool IsPionRejected(float mom, float nsigmaTPCPr, float nsigmaTOFPr, float nsigmaTPCK, float nsigmaTOFK)
  {
    if (mom < 0.5) {
      if (TMath::Abs(nsigmaTPCK) < ConfBothTracks.ConfPIDKaonNsigmaReject.value) {
        return true;
      } else if (TMath::Abs(nsigmaTPCPr) < ConfBothTracks.ConfPIDProtonNsigmaReject.value) {
        return true;
      }
    }
    if (mom > 0.5) {
      if (TMath::Hypot(nsigmaTOFK, nsigmaTPCK) < ConfBothTracks.ConfPIDKaonNsigmaReject.value) {
        return true;
      } else if (TMath::Hypot(nsigmaTOFPr, nsigmaTPCPr) < ConfBothTracks.ConfPIDProtonNsigmaReject.value) {
        return true;
      } else {
        return false;
      }
    } else {
      return false;
    }
  }

  bool IsParticleNSigmaAccepted(float mom, float nsigmaTPCPr, float nsigmaTOFPr, float nsigmaTPCPi, float nsigmaTOFPi, float nsigmaTPCK, float nsigmaTOFK)
  {
    switch (ConfTrack.ConfTrackPDGCode) {
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

  bool IsParticleNSigmaRejected(float mom, float nsigmaTPCPr, float nsigmaTOFPr, float nsigmaTPCPi, float nsigmaTOFPi, float nsigmaTPCK, float nsigmaTOFK)
  {
    switch (ConfTrack.ConfTrackPDGCode) {
      case 2212:  // Proton
      case -2212: // anty Proton
        return IsProtonRejected(mom, nsigmaTPCPi, nsigmaTOFPi, nsigmaTPCK, nsigmaTOFK);
        break;
      case 211:  // Pion
      case -211: // Pion-
        return IsPionRejected(mom, nsigmaTPCPr, nsigmaTOFPr, nsigmaTPCK, nsigmaTOFK);
        break;
      case 321:  // Kaon+
      case -321: // Kaon-
        return IsKaonRejected(mom, nsigmaTPCPr, nsigmaTOFPr, nsigmaTPCPi, nsigmaTOFPi);
        break;
      default:
        return false;
    }
  }

  void init(InitContext&)
  {
    eventHisto.init(&qaRegistry);
    qaRegistry.add("PhiDaugh_pos/nSigmaTPC", "; #it{p} (GeV/#it{c}); n#sigma_{TPC}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("PhiDaugh_pos/nSigmaTOF", "; #it{p} (GeV/#it{c}); n#sigma_{TOF}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("PhiDaugh_pos/pt", "; #it{p_T} (GeV/#it{c}); Counts", kTH1F, {{100, 0, 10}});
    qaRegistry.add("PhiDaugh_pos/eta", "; #it{eta}; Counts", kTH1F, {{200, -1.5, 1.5}});
    qaRegistry.add("PhiDaugh_pos/phi", "; #it{varphi}; Counts", kTH1F, {{200, 0, 2. * M_PI}});
    qaRegistry.add("PhiDaugh_pos/hDCAxy", "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", kTH2F, {{100, 0, 10}, {500, -5, 5}});

    qaRegistry.add("PhiDaugh_neg/nSigmaTPC", "; #it{p} (GeV/#it{c}); n#sigma_{TPC}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("PhiDaugh_neg/nSigmaTOF", "; #it{p} (GeV/#it{c}); n#sigma_{TOF}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("PhiDaugh_neg/pt", "; #it{p_T} (GeV/#it{c}); Counts", kTH1F, {{100, 0, 10}});
    qaRegistry.add("PhiDaugh_neg/eta", "; #it{eta}; Counts", kTH1F, {{200, -1.5, 1.5}});
    qaRegistry.add("PhiDaugh_neg/phi", "; #it{varphi}; Counts", kTH1F, {{200, 0, 2. * M_PI}});
    qaRegistry.add("PhiDaugh_neg/hDCAxy", "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", kTH2F, {{100, 0, 10}, {500, -5, 5}});

    registryPhiMinvBackground.add("InvariantMassKpKp", "; invariant mass K+K+; Counts", kTH1F, {{6000, 0.9, 4.0}});
    registryPhiMinvBackground.add("InvariantMassKmKm", "; invariant mass K-K-; Counts", kTH1F, {{6000, 0.9, 4.0}});

    qaRegistry.add("Hadron/nSigmaTPCPr", "; #it{p} (GeV/#it{c}); n#sigma_{TPCPr}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("Hadron/nSigmaTOFPr", "; #it{p} (GeV/#it{c}); n#sigma_{TOFPr}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("Hadron/nSigmaTPCPi", "; #it{p} (GeV/#it{c}); n#sigma_{TPCPi}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("Hadron/nSigmaTOFPi", "; #it{p} (GeV/#it{c}); n#sigma_{TOFPi}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("Hadron/nSigmaTPCKa", "; #it{p} (GeV/#it{c}); n#sigma_{TPCKa}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("Hadron/nSigmaTOFKa", "; #it{p} (GeV/#it{c}); n#sigma_{TOFKa}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});

    // MC truth
    registryMCtruth.add("MCtruthPhi", "MC truth Phi;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 5}, {400, -1.0, 1.0}}});
    registryMCtruth.add("MCtruthAllPositivePt", "MC truth all positive;#it{p}_{T} (GeV/c); #eta", {HistType::kTH1F, {{500, 0, 5}}});
    registryMCtruth.add("MCtruthAllNegativePt", "MC truth all negative;#it{p}_{T} (GeV/c); #eta", {HistType::kTH1F, {{500, 0, 5}}});
    registryMCtruth.add("MCtruthKp", "MC truth K+;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 5}, {400, -1.0, 1.0}}});
    registryMCtruth.add("MCtruthKm", "MC truth protons;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 5}, {400, -1.0, 1.0}}});
    registryMCtruth.add("MCtruthKpPt", "MC truth kaons positive;#it{p}_{T} (GeV/c)", {HistType::kTH1F, {{500, 0, 5}}});
    registryMCtruth.add("MCtruthKmPt", "MC truth kaons negative;#it{p}_{T} (GeV/c)", {HistType::kTH1F, {{500, 0, 5}}});
    registryMCtruth.add("MCtruthPpos", "MC truth proton;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 5}, {400, -1.0, 1.0}}});
    registryMCtruth.add("MCtruthPneg", "MC truth antiproton;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 5}, {400, -1.0, 1.0}}});

    // MC reco
    registryMCreco.add("MCrecoPhi", "MC reco Phi;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 5}, {400, -1.0, 1.0}}});
    registryMCreco.add("MCrecoAllPositivePt", "MC reco all;#it{p}_{T} (GeV/c); #eta", {HistType::kTH1F, {{500, 0, 5}}});
    registryMCreco.add("MCrecoAllNegativePt", "MC reco all;#it{p}_{T} (GeV/c); #eta", {HistType::kTH1F, {{500, 0, 5}}});
    registryMCreco.add("MCrecoPpos", "MC reco proton;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 5}, {400, -1.0, 1.0}}});
    registryMCreco.add("MCrecoPneg", "MC reco antiproton;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 5}, {400, -1.0, 1.0}}});

    trackHistoPartPhi.init(&qaRegistry, ConfBinsTempFitVarpT, ConfBinsTempFitVarInvMass, ConfBothTracks.ConfIsMC, ConfPhi.ConfPhiPDGCode);
    if (!ConfTrack.ConfTrackIsSame) {
      trackHistoPartTrack.init(&qaRegistry, ConfBinsTempFitVarpT, ConfBinsTempFitVar, ConfBothTracks.ConfIsMC, ConfTrack.ConfTrackPDGCode);
    }

    MixQaRegistry.add("MixingQA/hSECollisionBins", ";bin;Entries", kTH1F, {{120, -0.5, 119.5}});
    MixQaRegistry.add("MixingQA/hMECollisionBins", ";bin;Entries", kTH1F, {{120, -0.5, 119.5}});

    sameEventAngularCont.init(&resultRegistry, ConfBinskstar, ConfBinsMult, ConfBinskT, ConfBinsmT, ConfBins3Dmult, ConfBins3DmT, ConfBothTracks.ConfBinsEta, ConfBothTracks.ConfBinsPhi, ConfBothTracks.ConfIsMC, ConfBothTracks.ConfUse3D);
    mixedEventAngularCont.init(&resultRegistry, ConfBinskstar, ConfBinsMult, ConfBinskT, ConfBinsmT, ConfBins3Dmult, ConfBins3DmT, ConfBothTracks.ConfBinsEta, ConfBothTracks.ConfBinsPhi, ConfBothTracks.ConfIsMC, ConfBothTracks.ConfUse3D);

    sameEventAngularCont.setPDGCodes(ConfPhi.ConfPhiPDGCode, ConfTrack.ConfTrackPDGCode);
    mixedEventAngularCont.setPDGCodes(ConfPhi.ConfPhiPDGCode, ConfTrack.ConfTrackPDGCode);

    pairCleaner.init(&qaRegistry);
    if (ConfCPR.ConfCPRIsEnabled.value) {
      pairCloseRejection.init(&resultRegistry, &qaRegistry, ConfCPR.ConfCPRdeltaPhiCutMin.value, ConfCPR.ConfCPRdeltaPhiCutMax.value, ConfCPR.ConfCPRdeltaEtaCutMin.value, ConfCPR.ConfCPRdeltaEtaCutMax.value, ConfCPR.ConfCPRChosenRadii.value, ConfCPR.ConfCPRPlotPerRadii.value, ConfCPR.ConfCPRInvMassCutMin.value, ConfCPR.ConfCPRInvMassCutMax.value);
    }

    /// Initializing CCDB
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    int64_t now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);

    if (!ConfEff.ConfEfficiencyTrackPath.value.empty()) {
      protoneff = ccdb->getForTimeStamp<TH2F>(ConfEff.ConfEfficiencyTrackPath.value.c_str(), ConfEff.ConfEfficiencyTrackTimestamp.value);
      if (!protoneff || protoneff->IsZombie()) {
        LOGF(fatal, "Could not load efficiency protoneff histogram from %s", ConfEff.ConfEfficiencyTrackPath.value.c_str());
      }
    }
    if (!ConfEff.ConfEfficiencyPhiPath.value.empty()) {
      phieff = ccdb->getForTimeStamp<TH2F>(ConfEff.ConfEfficiencyPhiPath.value.c_str(), ConfEff.ConfEfficiencyPhiTimestamp.value);
      if (!phieff || phieff->IsZombie()) {
        LOGF(fatal, "Could not load efficiency phieff histogram from %s", ConfEff.ConfEfficiencyPhiPath.value.c_str());
      }
    }

    vPIDPhiCandidate = ConfPhi.ConfPhiPID.value;
    vPIDTrack = ConfTrack.ConfTrackPID.value;
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
  void doSameEvent(PartitionType groupPartsTrack, PartitionType groupPartsPhi, PartitionType groupPartsPhiDaugh, PartitionType groupPartsKaons, PartType parts, float magFieldTesla, int multCol)
  {

    /// Histogramming same event
    for (auto& phicandidate : groupPartsPhi) {
      trackHistoPartPhi.fillQA<isMC, false>(phicandidate);
    }

    float tpcNSigma;
    float tofNSigma;
    for (auto& phidaugh : groupPartsPhiDaugh) {
      if (phidaugh.mAntiLambda() == 1) { // workaround
        tpcNSigma = trackCuts.getNsigmaTPC(phidaugh, o2::track::PID::Kaon);
        tofNSigma = trackCuts.getNsigmaTOF(phidaugh, o2::track::PID::Kaon);

        qaRegistry.fill(HIST("PhiDaugh_pos/nSigmaTPC"), phidaugh.p(), tpcNSigma);
        qaRegistry.fill(HIST("PhiDaugh_pos/nSigmaTOF"), phidaugh.p(), tofNSigma);
        qaRegistry.fill(HIST("PhiDaugh_pos/hDCAxy"), phidaugh.p(), phidaugh.tempFitVar());
        qaRegistry.fill(HIST("PhiDaugh_pos/pt"), phidaugh.pt());
        qaRegistry.fill(HIST("PhiDaugh_pos/eta"), phidaugh.eta());
        qaRegistry.fill(HIST("PhiDaugh_pos/phi"), phidaugh.phi());
      } else if (phidaugh.mAntiLambda() == -1) { // workaround
        tpcNSigma = trackCuts.getNsigmaTPC(phidaugh, o2::track::PID::Kaon);
        tofNSigma = trackCuts.getNsigmaTOF(phidaugh, o2::track::PID::Kaon);

        qaRegistry.fill(HIST("PhiDaugh_neg/nSigmaTPC"), phidaugh.p(), tpcNSigma);
        qaRegistry.fill(HIST("PhiDaugh_neg/nSigmaTOF"), phidaugh.p(), tofNSigma);
        qaRegistry.fill(HIST("PhiDaugh_neg/pt"), phidaugh.pt());
        qaRegistry.fill(HIST("PhiDaugh_neg/eta"), phidaugh.eta());
        qaRegistry.fill(HIST("PhiDaugh_neg/phi"), phidaugh.phi());
        qaRegistry.fill(HIST("PhiDaugh_neg/hDCAxy"), phidaugh.p(), phidaugh.tempFitVar());
      }
    }
    float tpcNSigmaPr, tofNSigmaPr, tpcNSigmaPi, tofNSigmaPi, tpcNSigmaKa, tofNSigmaKa;
    if (!ConfTrack.ConfTrackIsSame) {
      for (auto& track : groupPartsTrack) {
        tpcNSigmaPi = trackCuts.getNsigmaTPC(track, o2::track::PID::Pion);
        tofNSigmaPi = trackCuts.getNsigmaTOF(track, o2::track::PID::Pion);
        tpcNSigmaKa = trackCuts.getNsigmaTPC(track, o2::track::PID::Kaon);
        tofNSigmaKa = trackCuts.getNsigmaTOF(track, o2::track::PID::Kaon);
        tpcNSigmaPr = trackCuts.getNsigmaTPC(track, o2::track::PID::Proton);
        tofNSigmaPr = trackCuts.getNsigmaTOF(track, o2::track::PID::Proton);

        if (ConfTrack.ConfTrackIsIdentified) {
          if (!IsParticleNSigmaAccepted(track.p(), tpcNSigmaPr, tofNSigmaPr, tpcNSigmaPi, tofNSigmaPi, tpcNSigmaKa, tofNSigmaKa)) {
            continue;
          }
        }

        if (ConfTrack.ConfTrackIsRejected) {
          if (IsParticleNSigmaRejected(track.p(), tpcNSigmaPr, tofNSigmaPr, tpcNSigmaPi, tofNSigmaPi, tpcNSigmaKa, tofNSigmaKa)) {
            continue;
          }
        }

        trackHistoPartTrack.fillQA<isMC, false>(track);

        qaRegistry.fill(HIST("Hadron/nSigmaTPCPi"), track.p(), tpcNSigmaPi);
        qaRegistry.fill(HIST("Hadron/nSigmaTOFPi"), track.p(), tofNSigmaPi);
        qaRegistry.fill(HIST("Hadron/nSigmaTPCKa"), track.p(), tpcNSigmaKa);
        qaRegistry.fill(HIST("Hadron/nSigmaTOFKa"), track.p(), tofNSigmaKa);
        qaRegistry.fill(HIST("Hadron/nSigmaTPCPr"), track.p(), tpcNSigmaPr);
        qaRegistry.fill(HIST("Hadron/nSigmaTOFPr"), track.p(), tofNSigmaPr);
      }
    }
    /// Now build the combinations
    for (auto& [track, phicandidate] : combinations(CombinationsFullIndexPolicy(groupPartsTrack, groupPartsPhi))) {
      if (ConfTrack.ConfTrackIsIdentified) {
        if (!IsParticleNSigmaAccepted(track.p(), trackCuts.getNsigmaTPC(track, o2::track::PID::Proton), trackCuts.getNsigmaTOF(track, o2::track::PID::Proton), trackCuts.getNsigmaTPC(track, o2::track::PID::Pion), trackCuts.getNsigmaTOF(track, o2::track::PID::Pion), trackCuts.getNsigmaTPC(track, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(track, o2::track::PID::Kaon))) {
          continue;
        }
      }

      if (ConfTrack.ConfTrackIsRejected) {
        if (IsParticleNSigmaRejected(track.p(), trackCuts.getNsigmaTPC(track, o2::track::PID::Proton), trackCuts.getNsigmaTOF(track, o2::track::PID::Proton), trackCuts.getNsigmaTPC(track, o2::track::PID::Pion), trackCuts.getNsigmaTOF(track, o2::track::PID::Pion), trackCuts.getNsigmaTPC(track, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(track, o2::track::PID::Kaon))) {
          continue;
        }
      }

      // // Close Pair Rejection
      if (ConfCPR.ConfCPRIsEnabled.value) {
        if (pairCloseRejection.isClosePair(track, phicandidate, parts, magFieldTesla, femtoUniverseContainer::EventType::same)) {
          continue;
        }
      }

      // Track Cleaning
      if (!pairCleaner.isCleanPair(track, phicandidate, parts)) {
        continue;
      }

      float weight = 1.0f;
      if (phieff) {
        weight = protoneff->GetBinContent(protoneff->FindBin(track.pt(), track.eta())) * phieff->GetBinContent(phieff->FindBin(phicandidate.pt(), phicandidate.eta()));
        sameEventAngularCont.setPair<isMC>(track, phicandidate, multCol, ConfBothTracks.ConfUse3D, weight);
      } else {
        sameEventAngularCont.setPair<isMC>(track, phicandidate, multCol, ConfBothTracks.ConfUse3D, weight);
      }
    }

    TLorentzVector part1Vec;
    TLorentzVector part2Vec;

    float mMassOne = TDatabasePDG::Instance()->GetParticle(321)->Mass();  // FIXME: Get from the PDG service of the common header
    float mMassTwo = TDatabasePDG::Instance()->GetParticle(-321)->Mass(); // FIXME: Get from the PDG service of the common header

    for (auto& [kaon1, kaon2] : combinations(CombinationsStrictlyUpperIndexPolicy(groupPartsKaons, groupPartsKaons))) {
      // empty if statements commented out on 20241114 to get rid of MegaLinter errors
      // if (!IsKaonNSigma(kaon1.p(), trackCuts.getNsigmaTPC(kaon1, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(kaon1, o2::track::PID::Kaon))) {
      // }
      // if (!IsKaonNSigma(kaon2.p(), trackCuts.getNsigmaTPC(kaon2, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(kaon2, o2::track::PID::Kaon))) {
      // }
      if ((kaon1.mAntiLambda() == 1) && (kaon2.mAntiLambda() == 1)) {
        part1Vec.SetPtEtaPhiM(kaon1.pt(), kaon1.eta(), kaon1.phi(), mMassOne);
        part2Vec.SetPtEtaPhiM(kaon2.pt(), kaon2.eta(), kaon2.phi(), mMassOne);
        TLorentzVector sumVec(part1Vec);
        sumVec += part2Vec;
        registryPhiMinvBackground.fill(HIST("InvariantMassKpKp"), sumVec.M());
      }
      if ((kaon1.mAntiLambda() == -1) && (kaon2.mAntiLambda() == -1)) {
        part1Vec.SetPtEtaPhiM(kaon1.pt(), kaon1.eta(), kaon1.phi(), mMassTwo);
        part2Vec.SetPtEtaPhiM(kaon2.pt(), kaon2.eta(), kaon2.phi(), mMassTwo);

        TLorentzVector sumVec(part1Vec);
        sumVec += part2Vec;
        registryPhiMinvBackground.fill(HIST("InvariantMassKmKm"), sumVec.M());
      }
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
    auto thegroupPartsPhiDaugh = partsPhiDaugh->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto thegroupPartsKaons = partsKaons->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);

    doSameEvent<false>(thegroupPartsTrack, thegroupPartsPhi, thegroupPartsPhiDaugh, thegroupPartsKaons, parts, col.magField(), col.multNtr());
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
    auto thegroupPartsPhiDaugh = partsPhiDaughMC->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto thegroupPartsKaons = partsKaonsMC->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);

    doSameEvent<true>(thegroupPartsTrack, thegroupPartsPhi, thegroupPartsPhiDaugh, thegroupPartsKaons, parts, col.magField(), col.multNtr());
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
      if (ConfTrack.ConfTrackIsIdentified) {
        if (!IsParticleNSigmaAccepted(track.p(), trackCuts.getNsigmaTPC(track, o2::track::PID::Proton), trackCuts.getNsigmaTOF(track, o2::track::PID::Proton), trackCuts.getNsigmaTPC(track, o2::track::PID::Pion), trackCuts.getNsigmaTOF(track, o2::track::PID::Pion), trackCuts.getNsigmaTPC(track, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(track, o2::track::PID::Kaon))) {
          continue;
        }
      }

      if (ConfTrack.ConfTrackIsRejected) {
        if (IsParticleNSigmaRejected(track.p(), trackCuts.getNsigmaTPC(track, o2::track::PID::Proton), trackCuts.getNsigmaTOF(track, o2::track::PID::Proton), trackCuts.getNsigmaTPC(track, o2::track::PID::Pion), trackCuts.getNsigmaTOF(track, o2::track::PID::Pion), trackCuts.getNsigmaTPC(track, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(track, o2::track::PID::Kaon))) {
          continue;
        }
      }

      if (ConfCPR.ConfCPRIsEnabled.value) {
        if (pairCloseRejection.isClosePair(track, phicandidate, parts, magFieldTesla, femtoUniverseContainer::EventType::mixed)) {
          continue;
        }
      }

      float weight = 1.0f;
      if (protoneff) {
        weight = protoneff->GetBinContent(protoneff->FindBin(track.pt(), track.eta())) * phieff->GetBinContent(phieff->FindBin(phicandidate.pt(), phicandidate.eta()));
        mixedEventAngularCont.setPair<isMC>(track, phicandidate, multCol, ConfBothTracks.ConfUse3D, weight);
      } else {
        mixedEventAngularCont.setPair<isMC>(track, phicandidate, multCol, ConfBothTracks.ConfUse3D, weight);
      }
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

  ///--------------------------------------------MC-------------------------------------------------///

  /// This function fills MC truth particles from derived MC table
  void processMCTruth(aod::FDParticles const& parts)
  {
    for (auto& part : parts) {
      if (part.partType() != uint8_t(aod::femtouniverseparticle::ParticleType::kMCTruthTrack))
        continue;

      int pdgCode = static_cast<int>(part.pidcut());
      const auto& pdgParticle = pdgMC->GetParticle(pdgCode);
      if (!pdgParticle) {
        continue;
      }

      if (pdgParticle->Charge() > 0.0) {
        registryMCtruth.fill(HIST("MCtruthAllPositivePt"), part.pt());
      }
      if (pdgCode == 321) {
        registryMCtruth.fill(HIST("MCtruthKp"), part.pt(), part.eta());
        registryMCtruth.fill(HIST("MCtruthKpPt"), part.pt());
      }
      if (pdgCode == 333) {
        registryMCtruth.fill(HIST("MCtruthPhi"), part.pt(), part.eta());
        continue;
      }
      if (pdgCode == 2212) {
        registryMCtruth.fill(HIST("MCtruthPpos"), part.pt(), part.eta());
      }

      if (pdgParticle->Charge() < 0.0) {
        registryMCtruth.fill(HIST("MCtruthAllNegativePt"), part.pt());
      }
      if (pdgCode == -321) {
        registryMCtruth.fill(HIST("MCtruthKm"), part.pt(), part.eta());
        registryMCtruth.fill(HIST("MCtruthKmPt"), part.pt());
      }
      if (pdgCode == -2212) {
        registryMCtruth.fill(HIST("MCtruthPneg"), part.pt(), part.eta());
      }
    }
  }
  PROCESS_SWITCH(femtoUniversePairTaskTrackPhi, processMCTruth, "Process MC truth data", false);

  void processMCReco(FemtoRecoParticles const& parts, aod::FDMCParticles const& mcparts)
  {
    for (auto& part : parts) {
      auto mcPartId = part.fdMCParticleId();
      if (mcPartId == -1)
        continue; // no MC particle
      const auto& mcpart = mcparts.iteratorAt(mcPartId);
      if (part.partType() == aod::femtouniverseparticle::ParticleType::kPhi) {
        if ((mcpart.pdgMCTruth() == 333) && (mcpart.partOriginMCTruth() == aod::femtouniverseMCparticle::ParticleOriginMCTruth::kPrimary)) {
          registryMCreco.fill(HIST("MCrecoPhi"), mcpart.pt(), mcpart.eta()); // phi
        }
      } else if (part.partType() == aod::femtouniverseparticle::ParticleType::kTrack) {
        if (part.sign() > 0) {
          registryMCreco.fill(HIST("MCrecoAllPositivePt"), mcpart.pt());
          if (mcpart.pdgMCTruth() == 2212 && IsParticleNSigmaAccepted(part.p(), trackCuts.getNsigmaTPC(part, o2::track::PID::Proton), trackCuts.getNsigmaTOF(part, o2::track::PID::Proton), trackCuts.getNsigmaTPC(part, o2::track::PID::Pion), trackCuts.getNsigmaTOF(part, o2::track::PID::Pion), trackCuts.getNsigmaTPC(part, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(part, o2::track::PID::Kaon))) {
            registryMCreco.fill(HIST("MCrecoPpos"), mcpart.pt(), mcpart.eta());
          }
        }

        if (part.sign() < 0) {
          registryMCreco.fill(HIST("MCrecoAllNegativePt"), mcpart.pt());
          if (mcpart.pdgMCTruth() == -2212 && IsParticleNSigmaAccepted(part.p(), trackCuts.getNsigmaTPC(part, o2::track::PID::Proton), trackCuts.getNsigmaTOF(part, o2::track::PID::Proton), trackCuts.getNsigmaTPC(part, o2::track::PID::Pion), trackCuts.getNsigmaTOF(part, o2::track::PID::Pion), trackCuts.getNsigmaTPC(part, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(part, o2::track::PID::Kaon))) {
            registryMCreco.fill(HIST("MCrecoPneg"), mcpart.pt(), mcpart.eta());
          }
        }
      } // partType
    }
  }

  PROCESS_SWITCH(femtoUniversePairTaskTrackPhi, processMCReco, "Process MC reco data", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<femtoUniversePairTaskTrackPhi>(cfgc),
  };
  return workflow;
}
