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
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseContainer.h"
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

  Configurable<float> ConfZVertexCut{"ConfZVertexCut", 10.f, "Event sel: Maximum z-Vertex (cm)"};
  Configurable<int> ConfNEventsMix{"ConfNEventsMix", 5, "Number of events for mixing"};
  Filter collisionFilter = (nabs(aod::collision::posZ) < ConfZVertexCut);
  using FilteredFDCollisions = soa::Filtered<o2::aod::FdCollisions>;
  using FilteredFDCollision = FilteredFDCollisions::iterator;

  Configurable<bool> ConfCPRIsEnabled{"ConfCPRIsEnabled", false, "Close Pair Rejection"};
  Configurable<bool> ConfCPRPlotPerRadii{"ConfCPRPlotPerRadii", false, "Plot CPR per radii"};
  Configurable<float> ConfCPRdeltaPhiCutMax{"ConfCPRdeltaPhiCutMax", 0.0, "Delta Phi max cut for Close Pair Rejection"};
  Configurable<float> ConfCPRdeltaPhiCutMin{"ConfCPRdeltaPhiCutMin", 0.0, "Delta Phi min cut for Close Pair Rejection"};
  Configurable<float> ConfCPRdeltaEtaCutMax{"ConfCPRdeltaEtaCutMax", 0.0, "Delta Eta max cut for Close Pair Rejection"};
  Configurable<float> ConfCPRdeltaEtaCutMin{"ConfCPRdeltaEtaCutMin", 0.0, "Delta Eta min cut for Close Pair Rejection"};
  Configurable<float> ConfCPRInvMassCutMin{"ConfCPRInvMassCutMin", 1.014, "Invariant mass (low) cut for Close Pair Rejection"};
  Configurable<float> ConfCPRInvMassCutMax{"ConfCPRInvMassCutMax", 1.026, "Invariant mass (high) cut for Close Pair Rejection"};
  Configurable<float> ConfCPRChosenRadii{"ConfCPRChosenRadii", 0.80, "Delta Eta cut for Close Pair Rejection"};

  /// Table for both particles
  Configurable<float> ConfPIDProtonNsigmaCombined{"ConfPIDProtonNsigmaCombined", 3.0, "TPC and TOF Proton Sigma (combined) for momentum > 0.5"};
  Configurable<float> ConfPIDProtonNsigmaTPC{"ConfPIDProtonNsigmaTPC", 3.0, "TPC Proton Sigma for momentum < 0.5"};
  Configurable<float> ConfPIDKaonNsigmaReject{"ConfPIDKaonNsigmaReject", 3.0, "Reject if particle could be a Kaon combined nsigma value."};
  Configurable<float> ConfPIDPionNsigmaReject{"ConfPIDPionNsigmaReject", 3.0, "Reject if particle could be a Pion combined nsigma value."};
  Configurable<float> ConfPIDProtonNsigmaReject{"ConfPIDProtonNsigmaReject", 3.0, "Reject if particle could be a Proton combined nsigma value."};
  Configurable<float> ConfPIDPionNsigmaCombined{"ConfPIDPionNsigmaCombined", 3.0, "TPC and TOF Pion Sigma (combined) for momentum > 0.5"};
  Configurable<float> ConfPIDPionNsigmaTPC{"ConfPIDPionNsigmaTPC", 3.0, "TPC Pion Sigma for momentum < 0.5"};
  Configurable<bool> ConfIsMC{"ConfIsMC", false, "Enable additional Histograms in the case of a MonteCarlo Run"};
  Configurable<bool> ConfUse3D{"ConfUse3D", false, "Enable three dimensional histogramms (to be used only for analysis with high statistics): k* vs mT vs multiplicity"};
  Configurable<int> ConfBinsPhi{"ConfBinsPhi", 29, "Number of phi bins in deta dphi"};
  Configurable<int> ConfBinsEta{"ConfBinsEta", 29, "Number of eta bins in deta dphi"};

  /// Particle 1 --- IDENTIFIED TRACK
  Configurable<int> ConfTrackPDGCode{"ConfTrackPDGCode", 2212, "Particle 2 - PDG code"};
  Configurable<int> ConfTrackSign{"ConfTrackSign", 1, "Track sign"};
  Configurable<bool> ConfTrackIsIdentified{"ConfTrackIsIdentified", true, "Enable PID for the track"};
  Configurable<bool> ConfTrackIsRejected{"ConfTrackIsRejected", true, "Enable PID rejection for the track other species than the identified one."};
  Configurable<float> ConfTrackPtPIDLimit{"ConfTrackPtPIDLimit", 0.5, "Momentum threshold for change of the PID method (from using TPC to TPC and TOF)."};
  Configurable<float> ConfTrackPtLow{"ConfTrackPtLow", 0.5, "Lower limit of the hadron pT."};
  Configurable<float> ConfTrackPtHigh{"ConfTrackPtHigh", 2.5, "Higher limit of the hadron pT."};

  /// Partitions for the track (particle 1)
  Partition<FilteredFemtoFullParticles> partsTrack = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kTrack)) &&
                                                     (aod::femtouniverseparticle::sign == ConfTrackSign) &&
                                                     (aod::femtouniverseparticle::pt > ConfTrackPtLow) &&
                                                     (aod::femtouniverseparticle::pt < ConfTrackPtHigh);

  Partition<FemtoRecoParticles> partsTrackMCReco = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kTrack)) &&
                                                   (aod::femtouniverseparticle::sign == ConfTrackSign) &&
                                                   (aod::femtouniverseparticle::pt > ConfTrackPtLow) &&
                                                   (aod::femtouniverseparticle::pt < ConfTrackPtHigh);

  /// Particle 2 --- PHI MESON
  Configurable<float> ConfPhiPtLow{"ConfPhiPtLow", 0.8, "Lower limit of the Phi pT."};
  Configurable<float> ConfPhiPtHigh{"ConfPhiPtHigh", 4.0, "Higher limit of the Phi pT."};
  Configurable<float> confInvMassLowLimitPhi{"confInvMassLowLimitPhi", 1.011, "Lower limit of the Phi invariant mass"}; // change that to do invariant mass cut
  Configurable<float> confInvMassUpLimitPhi{"confInvMassUpLimitPhi", 1.027, "Upper limit of the Phi invariant mass"};

  /// Partitions for the Phi meson (particle 2)
  Partition<FilteredFemtoFullParticles> partsPhi = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kPhi)) &&
                                                   (aod::femtouniverseparticle::pt > ConfPhiPtLow) &&
                                                   (aod::femtouniverseparticle::pt < ConfPhiPtHigh) &&
                                                   (aod::femtouniverseparticle::tempFitVar > confInvMassLowLimitPhi) &&
                                                   (aod::femtouniverseparticle::tempFitVar < confInvMassUpLimitPhi);

  Partition<FemtoRecoParticles> partsPhiMCReco = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kPhi)) &&
                                                 (aod::femtouniverseparticle::pt > ConfPhiPtLow) &&
                                                 (aod::femtouniverseparticle::pt < ConfPhiPtHigh) &&
                                                 (aod::femtouniverseparticle::tempFitVar > confInvMassLowLimitPhi) &&
                                                 (aod::femtouniverseparticle::tempFitVar < confInvMassUpLimitPhi);

  /// Partitions  for Phi daughters kPhiChild
  Partition<FilteredFemtoFullParticles> partsPhiDaugh = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kPhiChild));
  Partition<soa::Join<FilteredFemtoFullParticles, aod::FDMCLabels>> partsPhiDaughMC = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kPhiChild));

  // Partition for K+K- minv
  Partition<FilteredFemtoFullParticles> partsKaons = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kTrack));
  Partition<soa::Join<FilteredFemtoFullParticles, aod::FDMCLabels>> partsKaonsMC = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kTrack));

  /// Histogramming for particle 1
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kTrack, 1> trackHistoPartTrack;
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kTrack, 1> hTrackDCA;

  /// Histogramming for particle 2
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kPhi, 2> trackHistoPartPhi;

  /// Histogramming for Event
  FemtoUniverseEventHisto eventHisto;

  /// particle part
  ConfigurableAxis ConfBinsTempFitVar{"ConfBinsTempFitVar", {300, -0.15, 0.15}, "binning of the TempFitVar in the pT vs. TempFitVar plot"};
  ConfigurableAxis ConfBinsTempFitVarInvMass{"ConfBinsTempFitVarInvMass", {6000, 0.9, 4.0}, "binning of the TempFitVar in the pT vs. TempFitVar plot"};
  ConfigurableAxis ConfBinsTempFitVarpT{"ConfBinsTempFitVarpT", {20, 0.5, 4.05}, "pT binning of the pT vs. TempFitVar plot"};
  ConfigurableAxis ConfBinsTempFitVarPDG{"ConfBinsTempFitVarPDG", {6000, -2300, 2300}, "Binning of the PDG code in the pT vs. TempFitVar plot"};
  ConfigurableAxis ConfBinsTempFitVarDCA{"ConfBinsTempFitVarDCA", {300, -3.0, 3.0}, "binning of the TempFitVar in the pT vs. TempFitVar plot"};

  /// Correlation part
  ConfigurableAxis ConfBinsMult{"ConfBinsMult", {VARIABLE_WIDTH, 0.0f, 4.0f, 8.0f, 12.0f, 16.0f, 20.0f, 24.0f, 28.0f, 32.0f, 36.0f, 40.0f, 44.0f, 48.0f, 52.0f, 56.0f, 60.0f, 64.0f, 68.0f, 72.0f, 76.0f, 80.0f, 84.0f, 88.0f, 92.0f, 96.0f, 100.0f, 200.0f, 99999.f}, "Mixing bins - multiplicity"}; // \todo to be obtained from the hash task
  ConfigurableAxis ConfBinsVtx{"ConfBinsVtx", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis ConfBins3DmT{"ConfBins3DmT", {VARIABLE_WIDTH, 1.02f, 1.14f, 1.20f, 1.26f, 1.38f, 1.56f, 1.86f, 4.50f}, "mT Binning for the 3Dimensional plot: k* vs multiplicity vs mT (set <<ConfUse3D>> to true in order to use)"};
  ConfigurableAxis ConfBins3Dmult{"ConfBins3Dmult", {VARIABLE_WIDTH, 0.0f, 20.0f, 30.0f, 40.0f, 99999.0f}, "multiplicity Binning for the 3Dimensional plot: k* vs multiplicity vs mT (set <<ConfUse3D>> to true in order to use)"};

  ConfigurableAxis ConfBinskstar{"ConfBinskstar", {1500, 0., 6.}, "binning kstar"};
  ConfigurableAxis ConfBinskT{"ConfBinskT", {150, 0., 9.}, "binning kT"};
  ConfigurableAxis ConfBinsmT{"ConfBinsmT", {225, 0., 7.5}, "binning mT"};

  FemtoUniverseContainer<femto_universe_container::EventType::same, femto_universe_container::Observable::kstar> sameEventCont;
  FemtoUniverseContainer<femto_universe_container::EventType::mixed, femto_universe_container::Observable::kstar> mixedEventCont;
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
  HistogramRegistry registryDCA{"registryDCA", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry registryMCpT{"registryMCpT", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};

  EfficiencyConfigurableGroup effConfGroup;
  EfficiencyCalculator<TH1> efficiencyCalculator{&effConfGroup};

  float weight = 1;

  // PID for protons
  bool isProtonNSigma(float mom, float nsigmaTPCPr, float nsigmaTOFPr) // previous version from: https://github.com/alisw/AliPhysics/blob/master/PWGCF/FEMTOSCOPY/AliFemtoUser/AliFemtoMJTrackCut.cxx
  {
    if (mom < ConfTrackPtPIDLimit) {
      if (std::abs(nsigmaTPCPr) < ConfPIDProtonNsigmaTPC) {
        return true;
      } else {
        return false;
      }
    } else if (mom > ConfTrackPtPIDLimit) {
      if (std::hypot(nsigmaTOFPr, nsigmaTPCPr) < ConfPIDProtonNsigmaCombined) {
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
      if (std::hypot(nsigmaTOFPi, nsigmaTPCPi) < ConfPIDPionNsigmaReject) {
        return true;
      } else if (std::hypot(nsigmaTOFK, nsigmaTPCK) < ConfPIDKaonNsigmaReject) {
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
      if (std::abs(nsigmaTPCPi) < ConfPIDPionNsigmaReject) {
        return true;
      } else if (std::abs(nsigmaTPCPr) < ConfPIDProtonNsigmaReject) {
        return true;
      }
    }
    if (mom > 0.5) {
      if (std::hypot(nsigmaTOFPi, nsigmaTPCPi) < ConfPIDPionNsigmaReject) {
        return true;
      } else if (std::hypot(nsigmaTOFPr, nsigmaTPCPr) < ConfPIDProtonNsigmaReject) {
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
        if (std::abs(nsigmaTPCPi) < ConfPIDPionNsigmaTPC) {
          return true;
        } else {
          return false;
        }
      } else if (mom > 0.5) {
        if (std::hypot(nsigmaTOFPi, nsigmaTPCPi) < ConfPIDPionNsigmaCombined) {
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
      if (std::abs(nsigmaTPCK) < ConfPIDKaonNsigmaReject) {
        return true;
      } else if (std::abs(nsigmaTPCPr) < ConfPIDProtonNsigmaReject) {
        return true;
      }
    }
    if (mom > 0.5) {
      if (std::hypot(nsigmaTOFK, nsigmaTPCK) < ConfPIDKaonNsigmaReject) {
        return true;
      } else if (std::hypot(nsigmaTOFPr, nsigmaTPCPr) < ConfPIDProtonNsigmaReject) {
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
    switch (ConfTrackPDGCode) {
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
    switch (ConfTrackPDGCode) {
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
    if (ConfIsMC) {
      hTrackDCA.init(&registryDCA, ConfBinsTempFitVarpT, ConfBinsTempFitVarDCA, true, ConfTrackPDGCode, true);

      registryMCpT.add("MCReco/C_phi_pT", "; #it{p_T} (GeV/#it{c}); Counts", kTH1F, {{100, 0, 10}});
      registryMCpT.add("MCReco/NC_phi_pT", "; #it{p_T} (GeV/#it{c}); Counts", kTH1F, {{100, 0, 10}});

      registryMCpT.add("MCReco/C_p_pT", "; #it{p_T} (GeV/#it{c}); Counts", kTH1F, {{100, 0, 10}});
      registryMCpT.add("MCReco/NC_p_pT", "; #it{p_T} (GeV/#it{c}); Counts", kTH1F, {{100, 0, 10}});
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

    qaRegistry.add("Hadron_pos/nSigmaTPCPr", "; #it{p} (GeV/#it{c}); n#sigma_{TPCPr}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("Hadron_pos/nSigmaTOFPr", "; #it{p} (GeV/#it{c}); n#sigma_{TOFPr}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("Hadron_pos/nSigmaTPCPi", "; #it{p} (GeV/#it{c}); n#sigma_{TPCPi}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("Hadron_pos/nSigmaTOFPi", "; #it{p} (GeV/#it{c}); n#sigma_{TOFPi}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("Hadron_pos/nSigmaTPCKa", "; #it{p} (GeV/#it{c}); n#sigma_{TPCKa}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("Hadron_pos/nSigmaTOFKa", "; #it{p} (GeV/#it{c}); n#sigma_{TOFKa}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});

    qaRegistry.add("Hadron_neg/nSigmaTPCPr", "; #it{p} (GeV/#it{c}); n#sigma_{TPCPr}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("Hadron_neg/nSigmaTOFPr", "; #it{p} (GeV/#it{c}); n#sigma_{TOFPr}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("Hadron_neg/nSigmaTPCPi", "; #it{p} (GeV/#it{c}); n#sigma_{TPCPi}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("Hadron_neg/nSigmaTOFPi", "; #it{p} (GeV/#it{c}); n#sigma_{TOFPi}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("Hadron_neg/nSigmaTPCKa", "; #it{p} (GeV/#it{c}); n#sigma_{TPCKa}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("Hadron_neg/nSigmaTOFKa", "; #it{p} (GeV/#it{c}); n#sigma_{TOFKa}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});

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

    trackHistoPartPhi.init(&qaRegistry, ConfBinsTempFitVarpT, ConfBinsTempFitVarInvMass, ConfIsMC, 333);
    trackHistoPartTrack.init(&qaRegistry, ConfBinsTempFitVarpT, ConfBinsTempFitVar, ConfIsMC, ConfTrackPDGCode);

    mixQaRegistry.add("MixingQA/hSECollisionBins", ";bin;Entries", kTH1F, {{120, -0.5, 119.5}});
    mixQaRegistry.add("MixingQA/hMECollisionBins", ";bin;Entries", kTH1F, {{120, -0.5, 119.5}});

    sameEventCont.init(&resultRegistry, ConfBinskstar, ConfBinsMult, ConfBinskT, ConfBinsmT, ConfBins3Dmult, ConfBins3DmT, ConfBinsEta, ConfBinsPhi, ConfIsMC, ConfUse3D);
    mixedEventCont.init(&resultRegistry, ConfBinskstar, ConfBinsMult, ConfBinskT, ConfBinsmT, ConfBins3Dmult, ConfBins3DmT, ConfBinsEta, ConfBinsPhi, ConfIsMC, ConfUse3D);

    sameEventCont.setPDGCodes(333, ConfTrackPDGCode);
    mixedEventCont.setPDGCodes(333, ConfTrackPDGCode);

    pairCleaner.init(&qaRegistry);
    if (ConfCPRIsEnabled) {
      pairCloseRejection.init(&resultRegistry, &qaRegistry, ConfCPRdeltaPhiCutMin, ConfCPRdeltaPhiCutMax, ConfCPRdeltaEtaCutMin, ConfCPRdeltaEtaCutMax, ConfCPRChosenRadii, ConfCPRPlotPerRadii, ConfCPRInvMassCutMin, ConfCPRInvMassCutMax);
    }
  }

  template <typename PartitionType, typename PartType, typename MCParticles = std::nullptr_t>
  void doSameEvent(FilteredFDCollision const& col, PartType const& parts, PartitionType& groupPartsOne, PartitionType& groupPartsTwo, [[maybe_unused]] MCParticles mcParts = nullptr)
  {
    const auto& magFieldTesla = col.magField();
    const int multCol = col.multNtr();

    eventHisto.fillQA(col);

    for (auto const& phicandidate : groupPartsTwo) {
      // TODO: add phi meson minv cut here
      const auto& posChild = parts.iteratorAt(phicandidate.index() - 2);
      float tpcNSigmaKp = trackCuts.getNsigmaTPC(posChild, o2::track::PID::Kaon);
      float tofNSigmaKp = trackCuts.getNsigmaTOF(posChild, o2::track::PID::Kaon);
      qaRegistry.fill(HIST("PhiDaugh_pos/nSigmaTPC"), posChild.p(), tpcNSigmaKp);
      qaRegistry.fill(HIST("PhiDaugh_pos/nSigmaTOF"), posChild.p(), tofNSigmaKp);
      qaRegistry.fill(HIST("PhiDaugh_pos/hDCAxy"), posChild.p(), posChild.tempFitVar());
      qaRegistry.fill(HIST("PhiDaugh_pos/pt"), posChild.pt());
      qaRegistry.fill(HIST("PhiDaugh_pos/eta"), posChild.eta());
      qaRegistry.fill(HIST("PhiDaugh_pos/phi"), posChild.phi());

      const auto& negChild = parts.iteratorAt(phicandidate.index() - 1);
      float tpcNSigmaKm = trackCuts.getNsigmaTPC(negChild, o2::track::PID::Kaon);
      float tofNSigmaKm = trackCuts.getNsigmaTOF(negChild, o2::track::PID::Kaon);
      qaRegistry.fill(HIST("PhiDaugh_neg/nSigmaTPC"), negChild.p(), tpcNSigmaKm);
      qaRegistry.fill(HIST("PhiDaugh_neg/nSigmaTOF"), negChild.p(), tofNSigmaKm);
      qaRegistry.fill(HIST("PhiDaugh_neg/pt"), negChild.pt());
      qaRegistry.fill(HIST("PhiDaugh_neg/eta"), negChild.eta());
      qaRegistry.fill(HIST("PhiDaugh_neg/phi"), negChild.phi());
      qaRegistry.fill(HIST("PhiDaugh_neg/hDCAxy"), negChild.p(), negChild.tempFitVar());

      trackHistoPartPhi.fillQA<false, false>(phicandidate);
    }

    for (auto const& track : groupPartsOne) {
      float tpcNSigmaPi = trackCuts.getNsigmaTPC(track, o2::track::PID::Pion);
      float tofNSigmaPi = trackCuts.getNsigmaTOF(track, o2::track::PID::Pion);
      float tpcNSigmaKa = trackCuts.getNsigmaTPC(track, o2::track::PID::Kaon);
      float tofNSigmaKa = trackCuts.getNsigmaTOF(track, o2::track::PID::Kaon);
      float tpcNSigmaPr = trackCuts.getNsigmaTPC(track, o2::track::PID::Proton);
      float tofNSigmaPr = trackCuts.getNsigmaTOF(track, o2::track::PID::Proton);

      if (ConfTrackIsIdentified) {
        if (!isParticleNSigmaAccepted(track.p(), tpcNSigmaPr, tofNSigmaPr, tpcNSigmaPi, tofNSigmaPi, tpcNSigmaKa, tofNSigmaKa)) {
          continue;
        }
        if (ConfTrackIsRejected) {
          if (isParticleNSigmaRejected(track.p(), tpcNSigmaPr, tofNSigmaPr, tpcNSigmaPi, tofNSigmaPi, tpcNSigmaKa, tofNSigmaKa)) {
            continue;
          }
        }
      }
      if (track.sign() > 0) {
        qaRegistry.fill(HIST("Hadron_pos/nSigmaTPCPi"), track.p(), tpcNSigmaPi);
        qaRegistry.fill(HIST("Hadron_pos/nSigmaTOFPi"), track.p(), tofNSigmaPi);
        qaRegistry.fill(HIST("Hadron_pos/nSigmaTPCKa"), track.p(), tpcNSigmaKa);
        qaRegistry.fill(HIST("Hadron_pos/nSigmaTOFKa"), track.p(), tofNSigmaKa);
        qaRegistry.fill(HIST("Hadron_pos/nSigmaTPCPr"), track.p(), tpcNSigmaPr);
        qaRegistry.fill(HIST("Hadron_pos/nSigmaTOFPr"), track.p(), tofNSigmaPr);
      } else if (track.sign() < 0) {
        qaRegistry.fill(HIST("Hadron_neg/nSigmaTPCPi"), track.p(), tpcNSigmaPi);
        qaRegistry.fill(HIST("Hadron_neg/nSigmaTOFPi"), track.p(), tofNSigmaPi);
        qaRegistry.fill(HIST("Hadron_neg/nSigmaTPCKa"), track.p(), tpcNSigmaKa);
        qaRegistry.fill(HIST("Hadron_neg/nSigmaTOFKa"), track.p(), tofNSigmaKa);
        qaRegistry.fill(HIST("Hadron_neg/nSigmaTPCPr"), track.p(), tpcNSigmaPr);
        qaRegistry.fill(HIST("Hadron_neg/nSigmaTOFPr"), track.p(), tofNSigmaPr);
      }
      trackHistoPartTrack.fillQA<false, false>(track);
    }

    /// Now build the combinations
    for (auto const& [track, phicandidate] : combinations(CombinationsFullIndexPolicy(groupPartsOne, groupPartsTwo))) {
      // TODO: add phi inv mass cut here

      if (ConfTrackIsIdentified) {
        if (!isParticleNSigmaAccepted(track.p(), trackCuts.getNsigmaTPC(track, o2::track::PID::Proton), trackCuts.getNsigmaTOF(track, o2::track::PID::Proton), trackCuts.getNsigmaTPC(track, o2::track::PID::Pion), trackCuts.getNsigmaTOF(track, o2::track::PID::Pion), trackCuts.getNsigmaTPC(track, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(track, o2::track::PID::Kaon))) {
          continue;
        }
      }

      if (ConfTrackIsRejected) {
        if (isParticleNSigmaRejected(track.p(), trackCuts.getNsigmaTPC(track, o2::track::PID::Proton), trackCuts.getNsigmaTOF(track, o2::track::PID::Proton), trackCuts.getNsigmaTPC(track, o2::track::PID::Pion), trackCuts.getNsigmaTOF(track, o2::track::PID::Pion), trackCuts.getNsigmaTPC(track, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(track, o2::track::PID::Kaon))) {
          continue;
        }
      }

      // Close Pair Rejection
      if (ConfCPRIsEnabled) {
        if (pairCloseRejection.isClosePair(track, phicandidate, parts, magFieldTesla, femto_universe_container::EventType::same)) {
          continue;
        }
      }

      // Track Cleaning
      if (!pairCleaner.isCleanPair(track, phicandidate, parts)) {
        continue;
      }

      weight = efficiencyCalculator.getWeight(ParticleNo::ONE, phicandidate.pt()) * efficiencyCalculator.getWeight(ParticleNo::TWO, track.pt());

      if constexpr (std::is_same<PartType, FemtoRecoParticles>::value)
        sameEventCont.setPair<true>(track, phicandidate, multCol, ConfUse3D, weight);
      else
        sameEventCont.setPair<false>(track, phicandidate, multCol, ConfUse3D, weight);
    }

    // // Used for better fitting of invariant mass background.

    // TLorentzVector part1Vec;
    // TLorentzVector part2Vec;

    // float mMassOne = o2::constants::physics::MassKPlus;
    // float mMassTwo = o2::constants::physics::MassKMinus;

    // for (auto const& [kaon1, kaon2] : combinations(CombinationsStrictlyUpperIndexPolicy(groupPartsKaons, groupPartsKaons))) {
    //   if ((kaon1.mAntiLambda() == 1) && (kaon2.mAntiLambda() == 1)) {
    //     part1Vec.SetPtEtaPhiM(kaon1.pt(), kaon1.eta(), kaon1.phi(), mMassOne);
    //     part2Vec.SetPtEtaPhiM(kaon2.pt(), kaon2.eta(), kaon2.phi(), mMassOne);
    //     TLorentzVector sumVec(part1Vec);
    //     sumVec += part2Vec;
    //     registryPhiMinvBackground.fill(HIST("InvariantMassKpKp"), sumVec.M());
    //   }
    //   if ((kaon1.mAntiLambda() == -1) && (kaon2.mAntiLambda() == -1)) {
    //     part1Vec.SetPtEtaPhiM(kaon1.pt(), kaon1.eta(), kaon1.phi(), mMassTwo);
    //     part2Vec.SetPtEtaPhiM(kaon2.pt(), kaon2.eta(), kaon2.phi(), mMassTwo);

    //     TLorentzVector sumVec(part1Vec);
    //     sumVec += part2Vec;
    //     registryPhiMinvBackground.fill(HIST("InvariantMassKmKm"), sumVec.M());
    //   }
    // }
  }

  void processSameEvent(FilteredFDCollision const& col, FilteredFemtoFullParticles const& parts)
  {
    auto thegroupPartsTrack = partsTrack->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto thegroupPartsPhi = partsPhi->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    // auto thegroupPartsPhiDaugh = partsPhiDaugh->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    // auto thegroupPartsKaons = partsKaons->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);

    doSameEvent(col, parts, thegroupPartsTrack, thegroupPartsPhi);
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackPhi, processSameEvent, "Enable processing same event", true);

  template <typename PartitionType, typename PartType, typename MCParticles = std::nullptr_t>
  void doMixedEvent(FilteredFDCollisions const& cols, PartType const& parts, PartitionType& partitionPhi, PartitionType& partitionTrack, [[maybe_unused]] MCParticles mcParts = nullptr)
  {
    ColumnBinningPolicy<aod::collision::PosZ, aod::femtouniversecollision::MultNtr> colBinning{{ConfBinsVtx, ConfBinsMult}, true};

    auto mixedCollProcessFunc = [&](auto& collision1, auto& collision2) -> void {
      const int multCol = collision1.multNtr();

      auto groupPartsPhi = partitionPhi->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);
      auto groupPartsTrack = partitionTrack->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);

      const auto& magFieldTesla1 = collision1.magField();
      const auto& magFieldTesla2 = collision2.magField();

      if (magFieldTesla1 != magFieldTesla2) {
        return;
      }

      for (const auto& [phicandidate, track] : combinations(CombinationsFullIndexPolicy(groupPartsPhi, groupPartsTrack))) {
        // TODO: move here phi meson mass cut

        if (ConfTrackIsIdentified) {
          if (!isParticleNSigmaAccepted(track.p(), trackCuts.getNsigmaTPC(track, o2::track::PID::Proton), trackCuts.getNsigmaTOF(track, o2::track::PID::Proton), trackCuts.getNsigmaTPC(track, o2::track::PID::Pion), trackCuts.getNsigmaTOF(track, o2::track::PID::Pion), trackCuts.getNsigmaTPC(track, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(track, o2::track::PID::Kaon))) {
            continue;
          }
          if (ConfTrackIsRejected) {
            if (isParticleNSigmaRejected(track.p(), trackCuts.getNsigmaTPC(track, o2::track::PID::Proton), trackCuts.getNsigmaTOF(track, o2::track::PID::Proton), trackCuts.getNsigmaTPC(track, o2::track::PID::Pion), trackCuts.getNsigmaTOF(track, o2::track::PID::Pion), trackCuts.getNsigmaTPC(track, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(track, o2::track::PID::Kaon))) {
              continue;
            }
          }
        }

        if (ConfCPRIsEnabled) {
          if (pairCloseRejection.isClosePair(track, phicandidate, parts, magFieldTesla1, femto_universe_container::EventType::mixed)) {
            continue;
          }
        }

        weight = efficiencyCalculator.getWeight(ParticleNo::ONE, phicandidate.pt()) * efficiencyCalculator.getWeight(ParticleNo::TWO, track.pt());

        if constexpr (std::is_same<PartType, FemtoRecoParticles>::value)
          mixedEventCont.setPair<true>(track, phicandidate, multCol, ConfUse3D, weight);
        else
          mixedEventCont.setPair<false>(track, phicandidate, multCol, ConfUse3D, weight);
      }
    };

    for (const auto& [collision1, collision2] : soa::selfCombinations(colBinning, ConfNEventsMix, -1, cols, cols)) {
      mixedCollProcessFunc(collision1, collision2);
      mixQaRegistry.fill(HIST("MixingQA/hMECollisionBins"), colBinning.getBin({collision1.posZ(), collision1.multNtr()}));
    }
  }

  void processMixedEvent(FilteredFDCollisions const& cols, FilteredFemtoFullParticles const& parts)
  {
    doMixedEvent(cols, parts, partsPhi, partsTrack);
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackPhi, processMixedEvent, "Enable processing mixed events", true);

  ///--------------------------------------------MC-------------------------------------------------///
  void processSameEventMCReco(FilteredFDCollision const& col, FemtoRecoParticles const& parts, aod::FdMCParticles const& mcparts)
  {
    auto thegroupPartsTrack = partsTrackMCReco->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto thegroupPartsPhi = partsPhiMCReco->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    // auto thegroupPartsPhiDaugh = partsPhiDaugh->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    // auto thegroupPartsKaons = partsKaons->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);

    doSameEvent(col, parts, thegroupPartsTrack, thegroupPartsPhi, mcparts);
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackPhi, processSameEventMCReco, "Enable processing same event for MC Reco", true);

  void processMixedEventMCReco(FilteredFDCollisions const& cols, FemtoRecoParticles const& parts, aod::FdMCParticles const& mcparts)
  {
    doMixedEvent(cols, parts, partsPhiMCReco, partsTrackMCReco, mcparts);
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackPhi, processMixedEventMCReco, "Enable processing mixed events for MC Reco", false);

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

      if (mcpart.pdgMCTruth() == ConfTrackPDGCode && (part.pt() > ConfTrackPtLow) && (part.pt() < ConfTrackPtHigh) && isParticleNSigmaAccepted(part.p(), trackCuts.getNsigmaTPC(part, o2::track::PID::Proton), trackCuts.getNsigmaTOF(part, o2::track::PID::Proton), trackCuts.getNsigmaTPC(part, o2::track::PID::Pion), trackCuts.getNsigmaTOF(part, o2::track::PID::Pion), trackCuts.getNsigmaTPC(part, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(part, o2::track::PID::Kaon))) {
        registryMCpT.fill(HIST("MCReco/NC_p_pT"), part.pt());
        float weightTrack = efficiencyCalculator.getWeight(ParticleNo::TWO, part.pt());
        registryMCpT.fill(HIST("MCReco/C_p_pT"), part.pt(), weightTrack);
      }
      if ((mcpart.pdgMCTruth() == 333) && (part.partType() == aod::femtouniverseparticle::ParticleType::kPhi) && (part.pt() > ConfPhiPtLow) && (part.pt() < ConfPhiPtHigh)) {
        registryMCpT.fill(HIST("MCReco/NC_phi_pT"), part.pt());
        float weightPhi = efficiencyCalculator.getWeight(ParticleNo::ONE, part.pt());
        registryMCpT.fill(HIST("MCReco/C_phi_pT"), part.pt(), weightPhi);
      }

      if (isParticleNSigmaAccepted(part.p(), trackCuts.getNsigmaTPC(part, o2::track::PID::Proton), trackCuts.getNsigmaTOF(part, o2::track::PID::Proton), trackCuts.getNsigmaTPC(part, o2::track::PID::Pion), trackCuts.getNsigmaTOF(part, o2::track::PID::Pion), trackCuts.getNsigmaTPC(part, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(part, o2::track::PID::Kaon)))
        hTrackDCA.fillQA<true, true>(part);
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
        } else if (part.sign() < 0) {
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
