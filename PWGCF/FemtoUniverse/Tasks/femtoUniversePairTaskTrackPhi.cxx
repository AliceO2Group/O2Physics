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
/// \brief Tasks that reads the track tables used for the pairing and builds pairs for h-Phi angular correlation analysis
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
#include "Framework/O2DatabasePDGPlugin.h"
#include "ReconstructionDataFormats/PID.h"

#include "PWGCF/FemtoUniverse/DataModel/FemtoDerived.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseParticleHisto.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseEventHisto.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniversePairCleaner.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseAngularContainer.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseDetaDphiStar.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseTrackSelection.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseEfficiencyCalculator.h"
#include <TFile.h>
#include <TH1.h>

using namespace o2;
using namespace o2::analysis::femto_universe;
using namespace o2::analysis::femto_universe::efficiency;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

namespace
{
// static constexpr int NPart = 2;
// static constexpr int NCuts = 5;
static const std::vector<std::string> partNames{"PhiCandidate", "Track"};
static const std::vector<std::string> cutNames{"MaxPt", "PIDthr", "nSigmaTPC", "nSigmaTPCTOF", "MaxP"};
// static const float cutsTable[NPart][NCuts]{ //unused variable
//   {4.05f, 1.f, 3.f, 3.f, 100.f},
//   {4.05f, 1.f, 3.f, 3.f, 100.f}};
} // namespace

struct FemtoUniversePairTaskTrackPhi {

  Service<o2::framework::O2DatabasePDG> pdgMC;

  using FilteredFemtoFullParticles = soa::Join<aod::FDParticles, aod::FDExtParticles>;

  SliceCache cache;
  Preslice<FilteredFemtoFullParticles> perCol = aod::femtouniverseparticle::fdCollisionId;

  using FemtoRecoParticles = soa::Join<aod::FDParticles, aod::FDExtParticles, aod::FDMCLabels>;
  Preslice<FemtoRecoParticles> perColMC = aod::femtouniverseparticle::fdCollisionId;

  struct : o2::framework::ConfigurableGroup {
    Configurable<bool> confCPRIsEnabled{"confCPRIsEnabled", false, "Close Pair Rejection"};
    Configurable<bool> confCPRPlotPerRadii{"confCPRPlotPerRadii", false, "Plot CPR per radii"};
    Configurable<float> confCPRdeltaPhiCutMax{"confCPRdeltaPhiCutMax", 0.0, "Delta Phi max cut for Close Pair Rejection"};
    Configurable<float> confCPRdeltaPhiCutMin{"confCPRdeltaPhiCutMin", 0.0, "Delta Phi min cut for Close Pair Rejection"};
    Configurable<float> confCPRdeltaEtaCutMax{"confCPRdeltaEtaCutMax", 0.0, "Delta Eta max cut for Close Pair Rejection"};
    Configurable<float> confCPRdeltaEtaCutMin{"confCPRdeltaEtaCutMin", 0.0, "Delta Eta min cut for Close Pair Rejection"};
    Configurable<float> confCPRInvMassCutMin{"confCPRInvMassCutMin", 1.014, "Invariant mass (low) cut for Close Pair Rejection"};
    Configurable<float> confCPRInvMassCutMax{"confCPRInvMassCutMax", 1.026, "Invariant mass (high) cut for Close Pair Rejection"};
    Configurable<float> confCPRChosenRadii{"confCPRChosenRadii", 0.80, "Delta Eta cut for Close Pair Rejection"};
  } ConfCPR;

  /// Table for both particles
  struct : o2::framework::ConfigurableGroup {
    Configurable<float> confPIDProtonNsigmaCombined{"confPIDProtonNsigmaCombined", 3.0, "TPC and TOF Proton Sigma (combined) for momentum > 0.5"};
    Configurable<float> confPIDProtonNsigmaTPC{"confPIDProtonNsigmaTPC", 3.0, "TPC Proton Sigma for momentum < 0.5"};
    Configurable<float> confPIDKaonNsigmaReject{"confPIDKaonNsigmaReject", 3.0, "Reject if particle could be a Kaon combined nsigma value."};
    Configurable<float> confPIDPionNsigmaReject{"confPIDPionNsigmaReject", 3.0, "Reject if particle could be a Pion combined nsigma value."};
    Configurable<float> confPIDProtonNsigmaReject{"confPIDProtonNsigmaReject", 3.0, "Reject if particle could be a Proton combined nsigma value."};
    Configurable<float> confPIDPionNsigmaCombined{"confPIDPionNsigmaCombined", 3.0, "TPC and TOF Pion Sigma (combined) for momentum > 0.5"};
    Configurable<float> confPIDPionNsigmaTPC{"confPIDPionNsigmaTPC", 3.0, "TPC Pion Sigma for momentum < 0.5"};
    Configurable<bool> confIsMC{"confIsMC", false, "Enable additional Histograms in the case of a MonteCarlo Run"};
    Configurable<std::vector<float>> confTrkPIDnSigmaMax{"confTrkPIDnSigmaMax", std::vector<float>{4.f, 3.f, 2.f}, "This configurable needs to be the same as the one used in the producer task"};
    Configurable<bool> confUse3D{"confUse3D", false, "Enable three dimensional histogramms (to be used only for analysis with high statistics): k* vs mT vs multiplicity"};
    Configurable<int> confBinsPhi{"confBinsPhi", 29, "Number of phi bins in deta dphi"};
    Configurable<int> confBinsEta{"confBinsEta", 29, "Number of eta bins in deta dphi"};
  } ConfBothTracks;

  /// Particle 1 --- IDENTIFIED TRACK
  struct : o2::framework::ConfigurableGroup {
    Configurable<int> confTrackPDGCode{"confTrackPDGCode", 2212, "Particle 2 - PDG code"};
    Configurable<int> confTrackPID{"confTrackPID", 2, "Particle 2 - Read from cutCulator"}; // we also need the possibility to specify whether the bit is true/false ->std>>vector<std::pair<int, int>>
    Configurable<int> confTrackSign{"confTrackSign", 1, "Track sign"};
    Configurable<bool> confTrackIsIdentified{"confTrackIsIdentified", true, "Enable PID for the track"};
    Configurable<bool> confTrackIsRejected{"confTrackIsRejected", true, "Enable PID rejection for the track other species than the identified one."};
    Configurable<float> confTrackPtPIDLimit{"confTrackPtPIDLimit", 0.5, "Momentum threshold for change of the PID method (from using TPC to TPC and TOF)."};
    Configurable<float> confTrackPtLowLimit{"confTrackPtLowLimit", 0.5, "Lower limit of the hadron pT."};
    Configurable<float> confTrackPtHighLimit{"confTrackPtHighLimit", 2.5, "Higher limit of the hadron pT."};
  } ConfTrack;

  /// Partitions for the track (particle 1)
  Partition<FilteredFemtoFullParticles> partsTrack = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kTrack)) &&
                                                     (aod::femtouniverseparticle::sign == ConfTrack.confTrackSign.value) &&
                                                     (aod::femtouniverseparticle::pt > ConfTrack.confTrackPtLowLimit.value) &&
                                                     (aod::femtouniverseparticle::pt < ConfTrack.confTrackPtHighLimit.value);

  /// Particle 2 --- PHI MESON
  struct : o2::framework::ConfigurableGroup {
    Configurable<float> confPhiPtLowLimit{"confPhiPtLowLimit", 0.8, "Lower limit of the Phi pT."};
    Configurable<float> confPhiPtHighLimit{"confPhiPtHighLimit", 4.0, "Higher limit of the Phi pT."};
  } ConfPhi;

  /// Partitions for the Phi meson (particle 2)
  Partition<FilteredFemtoFullParticles> partsPhi = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kPhi)) &&
                                                   (aod::femtouniverseparticle::pt > ConfPhi.confPhiPtLowLimit.value) &&
                                                   (aod::femtouniverseparticle::pt < ConfPhi.confPhiPtHighLimit.value);

  /// Partitions  for Phi daughters kPhiChild
  Partition<FilteredFemtoFullParticles> partsPhiDaugh = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kPhiChild));
  Partition<soa::Join<FilteredFemtoFullParticles, aod::FDMCLabels>> partsPhiDaughMC = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kPhiChild));

  // Partition for K+K- minv
  Partition<FilteredFemtoFullParticles> partsKaons = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kTrack));
  Partition<soa::Join<FilteredFemtoFullParticles, aod::FDMCLabels>> partsKaonsMC = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kTrack));

  /// Histogramming for particle 1
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kTrack, 2> trackHistoPartTrack;
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kMCTruthTrack, 1> hMCTruth1;

  /// Histogramming for particle 2
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kPhi, 0> trackHistoPartPhi;
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kMCTruthTrack, 2> hMCTruth2;

  /// Histogramming for Event
  FemtoUniverseEventHisto eventHisto;

  /// particle part
  ConfigurableAxis confBinsTempFitVar{"confBinsTempFitVar", {300, -0.15, 0.15}, "binning of the TempFitVar in the pT vs. TempFitVar plot"};
  ConfigurableAxis confBinsTempFitVarInvMass{"confBinsTempFitVarInvMass", {6000, 0.9, 4.0}, "binning of the TempFitVar in the pT vs. TempFitVar plot"};
  ConfigurableAxis confBinsTempFitVarpT{"confBinsTempFitVarpT", {20, 0.5, 4.05}, "pT binning of the pT vs. TempFitVar plot"};
  ConfigurableAxis confBinsTempFitVarPDG{"confBinsTempFitVarPDG", {6000, -2300, 2300}, "Binning of the PDG code in the pT vs. TempFitVar plot"};

  /// Correlation part
  ConfigurableAxis confBinsMult{"confBinsMult", {VARIABLE_WIDTH, 0.0f, 4.0f, 8.0f, 12.0f, 16.0f, 20.0f, 24.0f, 28.0f, 32.0f, 36.0f, 40.0f, 44.0f, 48.0f, 52.0f, 56.0f, 60.0f, 64.0f, 68.0f, 72.0f, 76.0f, 80.0f, 84.0f, 88.0f, 92.0f, 96.0f, 100.0f, 200.0f, 99999.f}, "Mixing bins - multiplicity"}; // \todo to be obtained from the hash task
  ConfigurableAxis confBinsVtx{"confBinsVtx", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis confBins3DmT{"confBins3DmT", {VARIABLE_WIDTH, 1.02f, 1.14f, 1.20f, 1.26f, 1.38f, 1.56f, 1.86f, 4.50f}, "mT Binning for the 3Dimensional plot: k* vs multiplicity vs mT (set <<ConfBothTracks.confUse3D>> to true in order to use)"};
  ConfigurableAxis confBins3Dmult{"confBins3Dmult", {VARIABLE_WIDTH, 0.0f, 20.0f, 30.0f, 40.0f, 99999.0f}, "multiplicity Binning for the 3Dimensional plot: k* vs multiplicity vs mT (set <<ConfBothTracks.confUse3D>> to true in order to use)"};

  ColumnBinningPolicy<aod::collision::PosZ, aod::femtouniversecollision::MultNtr> colBinning{{confBinsVtx, confBinsMult}, true};

  ConfigurableAxis confBinskstar{"confBinskstar", {1500, 0., 6.}, "binning kstar"};
  ConfigurableAxis confBinskT{"confBinskT", {150, 0., 9.}, "binning kT"};
  ConfigurableAxis confBinsmT{"confBinsmT", {225, 0., 7.5}, "binning mT"};

  FemtoUniverseAngularContainer<femto_universe_angular_container::EventType::same, femto_universe_angular_container::Observable::kstar> sameEventAngularCont;
  FemtoUniverseAngularContainer<femto_universe_angular_container::EventType::mixed, femto_universe_angular_container::Observable::kstar> mixedEventAngularCont;
  FemtoUniversePairCleaner<aod::femtouniverseparticle::ParticleType::kTrack, aod::femtouniverseparticle::ParticleType::kPhi> pairCleaner;
  FemtoUniverseDetaDphiStar<aod::femtouniverseparticle::ParticleType::kTrack, aod::femtouniverseparticle::ParticleType::kPhi> pairCloseRejection;
  FemtoUniverseTrackSelection trackCuts;

  /// Histogram output
  HistogramRegistry qaRegistry{"qaRegistry", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry resultRegistry{"Correlations", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry mixQaRegistry{"mixQaRegistry", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registryMCtruth{"registryMCtruth", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry registryMCreco{"registryMCreco", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry registryPhiMinvBackground{"registryPhiMinvBackground", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};

  Configurable<bool> confDoEfficiency{"confDoEfficiency", true, "Do efficiency corrections."};
  EfficiencyConfigurableGroup effConfGroup;
  EfficiencyCalculator efficiencyCalculator{&effConfGroup};

  /// @brief Counter for particle swapping
  int fNeventsProcessed = 0;
  bool swpart = 0;
  float weight = 1;

  // PID for protons
  bool isProtonNSigma(float mom, float nsigmaTPCPr, float nsigmaTOFPr) // previous version from: https://github.com/alisw/AliPhysics/blob/master/PWGCF/FEMTOSCOPY/AliFemtoUser/AliFemtoMJTrackCut.cxx
  {
    if (mom < ConfTrack.confTrackPtPIDLimit.value) {
      if (std::abs(nsigmaTPCPr) < ConfBothTracks.confPIDProtonNsigmaTPC.value) {
        return true;
      } else {
        return false;
      }
    } else if (mom > ConfTrack.confTrackPtPIDLimit.value) {
      if (std::hypot(nsigmaTOFPr, nsigmaTPCPr) < ConfBothTracks.confPIDProtonNsigmaCombined.value) {
        return true;
      } else {
        return false;
      }
    }
    return false;
  }

  bool isProtonRejected(float mom, float nsigmaTPCPi, float nsigmaTOFPi, float nsigmaTPCK, float nsigmaTOFK)
  {
    if (mom < 0.5) {
      return true;
    }
    if (mom > 0.5) {
      if (std::hypot(nsigmaTOFPi, nsigmaTPCPi) < ConfBothTracks.confPIDPionNsigmaReject.value) {
        return true;
      } else if (std::hypot(nsigmaTOFK, nsigmaTPCK) < ConfBothTracks.confPIDKaonNsigmaReject.value) {
        return true;
      } else {
        return false;
      }
    } else {
      return false;
    }
  }

  bool isKaonNSigma(float mom, float nsigmaTPCK, float nsigmaTOFK)
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

  bool isKaonRejected(float mom, float nsigmaTPCPr, float nsigmaTOFPr, float nsigmaTPCPi, float nsigmaTOFPi)
  {
    if (mom < 0.5) {
      if (std::abs(nsigmaTPCPi) < ConfBothTracks.confPIDPionNsigmaReject.value) {
        return true;
      } else if (std::abs(nsigmaTPCPr) < ConfBothTracks.confPIDProtonNsigmaReject.value) {
        return true;
      }
    }
    if (mom > 0.5) {
      if (std::hypot(nsigmaTOFPi, nsigmaTPCPi) < ConfBothTracks.confPIDPionNsigmaReject.value) {
        return true;
      } else if (std::hypot(nsigmaTOFPr, nsigmaTPCPr) < ConfBothTracks.confPIDProtonNsigmaReject.value) {
        return true;
      } else {
        return false;
      }
    } else {
      return false;
    }
  }

  bool isPionNSigma(float mom, float nsigmaTPCPi, float nsigmaTOFPi)
  {
    if (true) {
      if (mom < 0.5) {
        if (std::abs(nsigmaTPCPi) < ConfBothTracks.confPIDPionNsigmaTPC.value) {
          return true;
        } else {
          return false;
        }
      } else if (mom > 0.5) {
        if (std::hypot(nsigmaTOFPi, nsigmaTPCPi) < ConfBothTracks.confPIDPionNsigmaCombined.value) {
          return true;
        } else {
          return false;
        }
      }
    }
    return false;
  }

  bool isPionRejected(float mom, float nsigmaTPCPr, float nsigmaTOFPr, float nsigmaTPCK, float nsigmaTOFK)
  {
    if (mom < 0.5) {
      if (std::abs(nsigmaTPCK) < ConfBothTracks.confPIDKaonNsigmaReject.value) {
        return true;
      } else if (std::abs(nsigmaTPCPr) < ConfBothTracks.confPIDProtonNsigmaReject.value) {
        return true;
      }
    }
    if (mom > 0.5) {
      if (std::hypot(nsigmaTOFK, nsigmaTPCK) < ConfBothTracks.confPIDKaonNsigmaReject.value) {
        return true;
      } else if (std::hypot(nsigmaTOFPr, nsigmaTPCPr) < ConfBothTracks.confPIDProtonNsigmaReject.value) {
        return true;
      } else {
        return false;
      }
    } else {
      return false;
    }
  }

  bool isParticleNSigmaAccepted(float mom, float nsigmaTPCPr, float nsigmaTOFPr, float nsigmaTPCPi, float nsigmaTOFPi, float nsigmaTPCK, float nsigmaTOFK)
  {
    switch (ConfTrack.confTrackPDGCode.value) {
      case 2212:  // Proton
      case -2212: // anty Proton
        return isProtonNSigma(mom, nsigmaTPCPr, nsigmaTOFPr);
        break;
      case 211:  // Pion
      case -211: // Pion-
        return isPionNSigma(mom, nsigmaTPCPi, nsigmaTOFPi);
        break;
      case 321:  // Kaon+
      case -321: // Kaon-
        return isKaonNSigma(mom, nsigmaTPCK, nsigmaTOFK);
        break;
      default:
        return false;
    }
  }

  bool isParticleNSigmaRejected(float mom, float nsigmaTPCPr, float nsigmaTOFPr, float nsigmaTPCPi, float nsigmaTOFPi, float nsigmaTPCK, float nsigmaTOFK)
  {
    switch (ConfTrack.confTrackPDGCode.value) {
      case 2212:  // Proton
      case -2212: // anty Proton
        return isProtonRejected(mom, nsigmaTPCPi, nsigmaTOFPi, nsigmaTPCK, nsigmaTOFK);
        break;
      case 211:  // Pion
      case -211: // Pion-
        return isPionRejected(mom, nsigmaTPCPr, nsigmaTOFPr, nsigmaTPCK, nsigmaTOFK);
        break;
      case 321:  // Kaon+
      case -321: // Kaon-
        return isKaonRejected(mom, nsigmaTPCPr, nsigmaTOFPr, nsigmaTPCPi, nsigmaTOFPi);
        break;
      default:
        return false;
    }
  }

  void init(InitContext&)
  {
    if (effConfGroup.confEfficiencyDoMCTruth) {
      hMCTruth1.init(&qaRegistry, confBinsTempFitVarpT, confBinsTempFitVarPDG, false, ConfTrack.confTrackPDGCode, false);
      hMCTruth2.init(&qaRegistry, confBinsTempFitVarpT, confBinsTempFitVarPDG, false, 333, false);
    }
    efficiencyCalculator.init();

    eventHisto.init(&qaRegistry);
    qaRegistry.add("PhiDaugh_pos/nSigmaTPC", "; #it{p} (GeV/#it{c}); n#sigma_{TPC}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("PhiDaugh_pos/nSigmaTOF", "; #it{p} (GeV/#it{c}); n#sigma_{TOF}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("PhiDaugh_pos/pt", "; #it{p_T} (GeV/#it{c}); Counts", kTH1F, {{100, 0, 10}});
    qaRegistry.add("PhiDaugh_pos/eta", "; #it{eta}; Counts", kTH1F, {{200, -1.5, 1.5}});
    qaRegistry.add("PhiDaugh_pos/phi", "; #it{varphi}; Counts", kTH1F, {{200, 0, o2::constants::math::TwoPI}});
    qaRegistry.add("PhiDaugh_pos/hDCAxy", "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", kTH2F, {{100, 0, 10}, {500, -5, 5}});

    qaRegistry.add("PhiDaugh_neg/nSigmaTPC", "; #it{p} (GeV/#it{c}); n#sigma_{TPC}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("PhiDaugh_neg/nSigmaTOF", "; #it{p} (GeV/#it{c}); n#sigma_{TOF}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("PhiDaugh_neg/pt", "; #it{p_T} (GeV/#it{c}); Counts", kTH1F, {{100, 0, 10}});
    qaRegistry.add("PhiDaugh_neg/eta", "; #it{eta}; Counts", kTH1F, {{200, -1.5, 1.5}});
    qaRegistry.add("PhiDaugh_neg/phi", "; #it{varphi}; Counts", kTH1F, {{200, 0, o2::constants::math::TwoPI}});
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
    registryMCtruth.add("MCtruthAllPositivePt", "MC truth all positive;#it{p}_{T} (GeV/c); #eta", {HistType::kTH1F, {{500, 0, 5}}});
    registryMCtruth.add("MCtruthAllNegativePt", "MC truth all negative;#it{p}_{T} (GeV/c); #eta", {HistType::kTH1F, {{500, 0, 5}}});
    // K+
    registryMCtruth.add("MCtruthKp", "MC truth K+;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 5}, {400, -1.0, 1.0}}});
    registryMCtruth.add("MCtruthKpPt", "MC truth kaons positive;#it{p}_{T} (GeV/c)", {HistType::kTH1F, {{500, 0, 5}}});
    // K-
    registryMCtruth.add("MCtruthKm", "MC truth protons;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 5}, {400, -1.0, 1.0}}});
    registryMCtruth.add("MCtruthKmPt", "MC truth kaons negative;#it{p}_{T} (GeV/c)", {HistType::kTH1F, {{500, 0, 5}}});
    // p
    registryMCtruth.add("MCtruthPpos", "MC truth protons;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 5}, {400, -1.0, 1.0}}});
    registryMCtruth.add("MCtruthPposPt", "MC truth protons;#it{p}_{T} (GeV/c)", {HistType::kTH1F, {{500, 0, 5}}});
    // pbar
    registryMCtruth.add("MCtruthPneg", "MC truth antiprotons;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 5}, {400, -1.0, 1.0}}});
    registryMCtruth.add("MCtruthPnegPt", "MC truth antiproton;#it{p}_{T} (GeV/c)", {HistType::kTH1F, {{500, 0, 5}}});
    // phi
    registryMCtruth.add("MCtruthPhi", "MC truth phi mesons;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 5}, {400, -1.0, 1.0}}});
    registryMCtruth.add("MCtruthPhiPt", "MC truth phi mesons;#it{p}_{T} (GeV/c)", {HistType::kTH1F, {{500, 0, 5}}});

    // MC reco
    registryMCreco.add("MCrecoAllPositivePt", "MC reco all;#it{p}_{T} (GeV/c); #eta", {HistType::kTH1F, {{500, 0, 5}}});
    registryMCreco.add("MCrecoAllNegativePt", "MC reco all;#it{p}_{T} (GeV/c); #eta", {HistType::kTH1F, {{500, 0, 5}}});
    // p
    registryMCreco.add("MCrecoPpos", "MC reco proton;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 5}, {400, -1.0, 1.0}}});
    registryMCreco.add("MCrecoPposPt", "MC reco proton; #it{p_T} (GeV/#it{c}); Counts", kTH1F, {{500, 0, 5}});
    // pbar
    registryMCreco.add("MCrecoPneg", "MC reco antiproton;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 5}, {400, -1.0, 1.0}}});
    registryMCreco.add("MCrecoPnegPt", "MC reco antiproton; #it{p_T} (GeV/#it{c}); Counts", kTH1F, {{500, 0, 5}});
    // phi
    registryMCreco.add("MCrecoPhi", "MC reco Phi;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 5}, {400, -1.0, 1.0}}});
    registryMCreco.add("MCrecoPhiPt", "MC reco Phi; #it{p_T} (GeV/#it{c}); Counts", kTH1F, {{500, 0, 5}});

    trackHistoPartPhi.init(&qaRegistry, confBinsTempFitVarpT, confBinsTempFitVarInvMass, ConfBothTracks.confIsMC.value, 333);
    trackHistoPartTrack.init(&qaRegistry, confBinsTempFitVarpT, confBinsTempFitVar, ConfBothTracks.confIsMC.value, ConfTrack.confTrackPDGCode.value);

    mixQaRegistry.add("MixingQA/hSECollisionBins", ";bin;Entries", kTH1F, {{120, -0.5, 119.5}});
    mixQaRegistry.add("MixingQA/hMECollisionBins", ";bin;Entries", kTH1F, {{120, -0.5, 119.5}});

    sameEventAngularCont.init(&resultRegistry, confBinskstar, confBinsMult, confBinskT, confBinsmT, confBins3Dmult, confBins3DmT, ConfBothTracks.confBinsEta, ConfBothTracks.confBinsPhi, ConfBothTracks.confIsMC, ConfBothTracks.confUse3D);
    mixedEventAngularCont.init(&resultRegistry, confBinskstar, confBinsMult, confBinskT, confBinsmT, confBins3Dmult, confBins3DmT, ConfBothTracks.confBinsEta, ConfBothTracks.confBinsPhi, ConfBothTracks.confIsMC, ConfBothTracks.confUse3D);

    sameEventAngularCont.setPDGCodes(333, ConfTrack.confTrackPDGCode.value);
    mixedEventAngularCont.setPDGCodes(333, ConfTrack.confTrackPDGCode.value);

    pairCleaner.init(&qaRegistry);
    if (ConfCPR.confCPRIsEnabled.value) {
      pairCloseRejection.init(&resultRegistry, &qaRegistry, ConfCPR.confCPRdeltaPhiCutMin.value, ConfCPR.confCPRdeltaPhiCutMax.value, ConfCPR.confCPRdeltaEtaCutMin.value, ConfCPR.confCPRdeltaEtaCutMax.value, ConfCPR.confCPRChosenRadii.value, ConfCPR.confCPRPlotPerRadii.value, ConfCPR.confCPRInvMassCutMin.value, ConfCPR.confCPRInvMassCutMax.value);
    }
  }

  template <typename CollisionType>
  void fillCollision(CollisionType col)
  {
    mixQaRegistry.fill(HIST("MixingQA/hSECollisionBins"), colBinning.getBin({col.posZ(), col.multNtr()}));
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
    // variables for particle swapping
    swpart = fNeventsProcessed % 2;
    fNeventsProcessed++;

    /// Histogramming same event
    for (auto const& phicandidate : groupPartsPhi) {
      trackHistoPartPhi.fillQA<isMC, false>(phicandidate);
    }

    float tpcNSigma;
    float tofNSigma;
    for (auto const& phidaugh : groupPartsPhiDaugh) {
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

    for (auto const& track : groupPartsTrack) {
      tpcNSigmaPi = trackCuts.getNsigmaTPC(track, o2::track::PID::Pion);
      tofNSigmaPi = trackCuts.getNsigmaTOF(track, o2::track::PID::Pion);
      tpcNSigmaKa = trackCuts.getNsigmaTPC(track, o2::track::PID::Kaon);
      tofNSigmaKa = trackCuts.getNsigmaTOF(track, o2::track::PID::Kaon);
      tpcNSigmaPr = trackCuts.getNsigmaTPC(track, o2::track::PID::Proton);
      tofNSigmaPr = trackCuts.getNsigmaTOF(track, o2::track::PID::Proton);

      if (ConfTrack.confTrackIsIdentified.value) {
        if (!isParticleNSigmaAccepted(track.p(), tpcNSigmaPr, tofNSigmaPr, tpcNSigmaPi, tofNSigmaPi, tpcNSigmaKa, tofNSigmaKa)) {
          continue;
        }
      }

      if (ConfTrack.confTrackIsRejected.value) {
        if (isParticleNSigmaRejected(track.p(), tpcNSigmaPr, tofNSigmaPr, tpcNSigmaPi, tofNSigmaPi, tpcNSigmaKa, tofNSigmaKa)) {
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

    /// Now build the combinations
    for (auto const& [track, phicandidate] : combinations(CombinationsFullIndexPolicy(groupPartsTrack, groupPartsPhi))) {
      if (ConfTrack.confTrackIsIdentified.value) {
        if (!isParticleNSigmaAccepted(track.p(), trackCuts.getNsigmaTPC(track, o2::track::PID::Proton), trackCuts.getNsigmaTOF(track, o2::track::PID::Proton), trackCuts.getNsigmaTPC(track, o2::track::PID::Pion), trackCuts.getNsigmaTOF(track, o2::track::PID::Pion), trackCuts.getNsigmaTPC(track, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(track, o2::track::PID::Kaon))) {
          continue;
        }
      }

      if (ConfTrack.confTrackIsRejected.value) {
        if (isParticleNSigmaRejected(track.p(), trackCuts.getNsigmaTPC(track, o2::track::PID::Proton), trackCuts.getNsigmaTOF(track, o2::track::PID::Proton), trackCuts.getNsigmaTPC(track, o2::track::PID::Pion), trackCuts.getNsigmaTOF(track, o2::track::PID::Pion), trackCuts.getNsigmaTPC(track, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(track, o2::track::PID::Kaon))) {
          continue;
        }
      }

      // Close Pair Rejection
      if (ConfCPR.confCPRIsEnabled.value) {
        if (pairCloseRejection.isClosePair(track, phicandidate, parts, magFieldTesla, femto_universe_container::EventType::same)) {
          continue;
        }
      }

      // Track Cleaning
      if (!pairCleaner.isCleanPair(track, phicandidate, parts)) {
        continue;
      }

      weight = 1.0f;
      if (confDoEfficiency.value) {
        weight = efficiencyCalculator.getWeight(ParticleNo::ONE, phicandidate) * efficiencyCalculator.getWeight(ParticleNo::TWO, track);
      }
      if (swpart)
        sameEventAngularCont.setPair<isMC>(track, phicandidate, multCol, ConfBothTracks.confUse3D.value, weight);
      else
        sameEventAngularCont.setPair<isMC>(phicandidate, track, multCol, ConfBothTracks.confUse3D.value, weight);

      swpart = !swpart;
    }

    // Used for better fitting of invariant mass background.

    TLorentzVector part1Vec;
    TLorentzVector part2Vec;

    float mMassOne = o2::constants::physics::MassKPlus;
    float mMassTwo = o2::constants::physics::MassKMinus;

    for (auto const& [kaon1, kaon2] : combinations(CombinationsStrictlyUpperIndexPolicy(groupPartsKaons, groupPartsKaons))) {
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
  void processSameEvent(o2::aod::FdCollision const& col,
                        FilteredFemtoFullParticles const& parts)
  {
    fillCollision(col);

    auto thegroupPartsTrack = partsTrack->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto thegroupPartsPhi = partsPhi->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto thegroupPartsPhiDaugh = partsPhiDaugh->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto thegroupPartsKaons = partsKaons->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);

    doSameEvent<false>(thegroupPartsTrack, thegroupPartsPhi, thegroupPartsPhiDaugh, thegroupPartsKaons, parts, col.magField(), col.multNtr());
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackPhi, processSameEvent, "Enable processing same event", true);

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
    // variables for particle swapping
    bool swpart = fNeventsProcessed % 2;
    fNeventsProcessed++;

    for (auto const& [track, phicandidate] : combinations(CombinationsFullIndexPolicy(groupPartsTrack, groupPartsPhi))) {
      if (ConfTrack.confTrackIsIdentified.value) {
        if (!isParticleNSigmaAccepted(track.p(), trackCuts.getNsigmaTPC(track, o2::track::PID::Proton), trackCuts.getNsigmaTOF(track, o2::track::PID::Proton), trackCuts.getNsigmaTPC(track, o2::track::PID::Pion), trackCuts.getNsigmaTOF(track, o2::track::PID::Pion), trackCuts.getNsigmaTPC(track, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(track, o2::track::PID::Kaon))) {
          continue;
        }
      }

      if (ConfTrack.confTrackIsRejected.value) {
        if (isParticleNSigmaRejected(track.p(), trackCuts.getNsigmaTPC(track, o2::track::PID::Proton), trackCuts.getNsigmaTOF(track, o2::track::PID::Proton), trackCuts.getNsigmaTPC(track, o2::track::PID::Pion), trackCuts.getNsigmaTOF(track, o2::track::PID::Pion), trackCuts.getNsigmaTPC(track, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(track, o2::track::PID::Kaon))) {
          continue;
        }
      }

      if (ConfCPR.confCPRIsEnabled.value) {
        if (pairCloseRejection.isClosePair(track, phicandidate, parts, magFieldTesla, femto_universe_container::EventType::mixed)) {
          continue;
        }
      }

      weight = 1.0f;
      if (confDoEfficiency.value) {
        weight = efficiencyCalculator.getWeight(ParticleNo::ONE, phicandidate) * efficiencyCalculator.getWeight(ParticleNo::TWO, track);
      }

      if (swpart)
        mixedEventAngularCont.setPair<isMC>(track, phicandidate, multCol, ConfBothTracks.confUse3D.value, weight);
      else
        mixedEventAngularCont.setPair<isMC>(phicandidate, track, multCol, ConfBothTracks.confUse3D.value, weight);

      swpart = !swpart;
    }
  }

  /// process function for to call doMixedEvent with Data
  /// @param cols subscribe to the collisions table (Data)
  /// @param parts subscribe to the femtoUniverseParticleTable
  void processMixedEvent(o2::aod::FdCollisions const& cols,
                         FilteredFemtoFullParticles const& parts)
  {
    for (auto const& [collision1, collision2] : soa::selfCombinations(colBinning, 5, -1, cols, cols)) {

      const int multiplicityCol = collision1.multNtr();
      mixQaRegistry.fill(HIST("MixingQA/hMECollisionBins"), colBinning.getBin({collision1.posZ(), multiplicityCol}));

      auto groupPartsTrack = partsTrack->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision2.globalIndex(), cache);
      auto groupPartsPhi = partsPhi->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);

      const auto& magFieldTesla1 = collision1.magField();
      const auto& magFieldTesla2 = collision2.magField();

      if (magFieldTesla1 != magFieldTesla2) {
        continue;
      }
      /// \todo before mixing we should check whether both collisions contain a pair of particles!
      // if (partsPhi.size() == 0 || NPart2Evt1 == 0 || NPart1Evt2 == 0 || partsTrack.size() == 0 ) continue;

      doMixedEvent<false>(groupPartsTrack, groupPartsPhi, parts, magFieldTesla1, multiplicityCol);
    }
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackPhi, processMixedEvent, "Enable processing mixed events", true);

  ///--------------------------------------------MC-------------------------------------------------///

  /// This function fills MC truth particles from derived MC table
  void processMCTruth(aod::FDParticles const& parts)
  {
    for (auto const& part : parts) {
      if (part.partType() != uint8_t(aod::femtouniverseparticle::ParticleType::kMCTruthTrack))
        continue;

      int pdgCode = static_cast<int>(part.pidCut());
      const auto& pdgParticle = pdgMC->GetParticle(pdgCode);
      if (!pdgParticle) {
        continue;
      }

      // charge +
      if (pdgParticle->Charge() > 0.0) {
        registryMCtruth.fill(HIST("MCtruthAllPositivePt"), part.pt());
        if (pdgCode == 2212) {
          registryMCtruth.fill(HIST("MCtruthPpos"), part.pt(), part.eta());
          registryMCtruth.fill(HIST("MCtruthPposPt"), part.pt());
          continue;
        } else if (pdgCode == 321) {
          registryMCtruth.fill(HIST("MCtruthKp"), part.pt(), part.eta());
          registryMCtruth.fill(HIST("MCtruthKpPt"), part.pt());
          continue;
        }
      }
      // charge 0
      if (pdgCode == 333) {
        registryMCtruth.fill(HIST("MCtruthPhi"), part.pt(), part.eta());
        registryMCtruth.fill(HIST("MCtruthPhiPt"), part.pt());
        continue;
      }

      // charge -
      if (pdgParticle->Charge() < 0.0) {
        registryMCtruth.fill(HIST("MCtruthAllNegativePt"), part.pt());

        if (pdgCode == -321) {
          registryMCtruth.fill(HIST("MCtruthKm"), part.pt(), part.eta());
          registryMCtruth.fill(HIST("MCtruthKmPt"), part.pt());
          continue;
        } else if (pdgCode == -2212) {
          registryMCtruth.fill(HIST("MCtruthPneg"), part.pt(), part.eta());
          registryMCtruth.fill(HIST("MCtruthPnegPt"), part.pt());
          continue;
        }
      }
    }
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackPhi, processMCTruth, "Process MC truth data", false);

  void processMCReco(FemtoRecoParticles const& parts, aod::FdMCParticles const& mcparts)
  {
    for (auto const& part : parts) {
      auto mcPartId = part.fdMCParticleId();
      if (mcPartId == -1)
        continue; // no MC particle
      const auto& mcpart = mcparts.iteratorAt(mcPartId);
      if ((part.partType() == aod::femtouniverseparticle::ParticleType::kPhi) && (mcpart.pdgMCTruth() == 333) && (mcpart.partOriginMCTruth() == aod::femtouniverse_mc_particle::ParticleOriginMCTruth::kPrimary)) {
        registryMCreco.fill(HIST("MCrecoPhi"), mcpart.pt(), mcpart.eta()); // phi
        registryMCreco.fill(HIST("MCrecoPhiPt"), mcpart.pt());
      } else if (part.partType() == aod::femtouniverseparticle::ParticleType::kTrack) {
        if (part.sign() > 0) {
          registryMCreco.fill(HIST("MCrecoAllPositivePt"), mcpart.pt());
          if (mcpart.pdgMCTruth() == 2212 && isParticleNSigmaAccepted(part.p(), trackCuts.getNsigmaTPC(part, o2::track::PID::Proton), trackCuts.getNsigmaTOF(part, o2::track::PID::Proton), trackCuts.getNsigmaTPC(part, o2::track::PID::Pion), trackCuts.getNsigmaTOF(part, o2::track::PID::Pion), trackCuts.getNsigmaTPC(part, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(part, o2::track::PID::Kaon))) {
            registryMCreco.fill(HIST("MCrecoPpos"), mcpart.pt(), mcpart.eta());
            registryMCreco.fill(HIST("MCrecoPposPt"), mcpart.pt());
          }
        }

        else if (part.sign() < 0) {
          registryMCreco.fill(HIST("MCrecoAllNegativePt"), mcpart.pt());
          if (mcpart.pdgMCTruth() == -2212 && isParticleNSigmaAccepted(part.p(), trackCuts.getNsigmaTPC(part, o2::track::PID::Proton), trackCuts.getNsigmaTOF(part, o2::track::PID::Proton), trackCuts.getNsigmaTPC(part, o2::track::PID::Pion), trackCuts.getNsigmaTOF(part, o2::track::PID::Pion), trackCuts.getNsigmaTPC(part, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(part, o2::track::PID::Kaon))) {
            registryMCreco.fill(HIST("MCrecoPneg"), mcpart.pt(), mcpart.eta());
            registryMCreco.fill(HIST("MCrecoPnegPt"), mcpart.pt());
          }
        }
      } // partType kTrack
    }
  }

  PROCESS_SWITCH(FemtoUniversePairTaskTrackPhi, processMCReco, "Process MC reco data", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<FemtoUniversePairTaskTrackPhi>(cfgc),
  };
  return workflow;
}
