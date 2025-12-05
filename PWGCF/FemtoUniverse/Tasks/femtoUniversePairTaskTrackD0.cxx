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

/// \file femtoUniversePairTaskTrackD0.cxx
/// \brief Tasks that reads the track tables and D0/D0bar mesons
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de
/// \author Georgios Mantzaridis, TU München, georgios.mantzaridis@tum.de
/// \author Anton Riedel, TU München, anton.riedel@tum.de
/// \author Zuzanna Chochulska, WUT Warsaw & CTU Prague, zchochul@cern.ch
/// \author Katarzyna Gwiździel, WUT Warsaw, katarzyna.gwizdziel@cern.ch

#include "PWGCF/FemtoUniverse/Core/FemtoUniverseContainer.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseDetaDphiStar.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseEfficiencyCalculator.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseEventHisto.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniversePairCleaner.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseParticleHisto.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseSoftPionRemoval.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseTrackSelection.h"
#include "PWGCF/FemtoUniverse/Core/femtoUtils.h"
#include "PWGCF/FemtoUniverse/DataModel/FemtoDerived.h"
#include "PWGHF/Core/DecayChannels.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

#include "Common/Core/RecoDecay.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/StepTHn.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/PID.h"

#include <string>
#include <vector>

using namespace o2;
using namespace o2::analysis;
using namespace o2::analysis::femto_universe;
using namespace o2::analysis::femto_universe::efficiency;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

namespace
{
static const int8_t genPromptD0 = 1;
static const int8_t genNonPromptD0 = 2;
} // namespace

/// Returns deltaPhi value within the range [-pi/2, 3/2*pi]
///
double getDeltaPhi(double phiD, double phiDbar)
{
  return RecoDecay::constrainAngle(phiDbar - phiD, -o2::constants::math::PIHalf);
}

/// Returns deltaPhi value within the range [0, pi]
///
double wrapDeltaPhi0PI(double phiD, double phiDbar)
{
  double deltaPhi = 0.0;
  deltaPhi = RecoDecay::constrainAngle(phiDbar - phiD, 0.0);
  if (deltaPhi > o2::constants::math::TwoPI) {
    deltaPhi = o2::constants::math::TwoPI - deltaPhi;
  }
  return deltaPhi;
}

struct FemtoUniversePairTaskTrackD0 {

  Service<o2::framework::O2DatabasePDG> pdgMC;

  using FemtoFullParticles = soa::Join<aod::FDParticles, aod::FDExtParticles>;
  SliceCache cache;
  Preslice<FemtoFullParticles> perCol = aod::femtouniverseparticle::fdCollisionId;

  using FemtoMcRecoParticles = soa::Join<aod::FDParticles, aod::FDExtParticles, aod::FDMCLabels>;
  Preslice<FemtoMcRecoParticles> perColMc = aod::femtouniverseparticle::fdCollisionId;

  /// Table for both particles
  struct : o2::framework::ConfigurableGroup {
    Configurable<float> confNsigmaCombinedProton{"confNsigmaCombinedProton", 3.0, "TPC and TOF Proton Sigma (combined)"};
    Configurable<float> confNsigmaTPCProton{"confNsigmaTPCProton", 3.0, "TPC Proton Sigma"};
    Configurable<float> confNsigmaCombinedPion{"confNsigmaCombinedPion", 3.0, "TPC and TOF Pion Sigma (combined)"};
    Configurable<float> confNsigmaTPCPion{"confNsigmaTPCPion", 3.0, "TPC Pion Sigma"};
    Configurable<float> confNsigmaCombinedKaon{"confNsigmaCombinedKaon", 3.0, "TPC and TOF Kaon Sigma (combined)"};
    Configurable<float> confNsigmaTPCKaon{"confNsigmaTPCKaon", 3.0, "TPC Kaon Sigma"};

    Configurable<bool> confIsMC{"confIsMC", false, "Enable additional Histogramms in the case of a MonteCarlo Run"};
    Configurable<std::vector<float>> confTrkPIDnSigmaMax{"confTrkPIDnSigmaMax", std::vector<float>{4.f, 3.f, 2.f}, "This configurable needs to be the same as the one used in the producer task"};
    Configurable<bool> confUse3D{"confUse3D", false, "Enable three dimensional histogramms (to be used only for analysis with high statistics): k* vs mT vs multiplicity"};
    Configurable<int> confPhiBins{"confPhiBins", 29, "Number of phi bins in deta dphi"};
    Configurable<int> confEtaBins{"confEtaBins", 29, "Number of eta bins in deta dphi"};
  } ConfBothTracks;

  /// Particle 1 --- IDENTIFIED TRACK
  struct : o2::framework::ConfigurableGroup {
    Configurable<bool> confIsSame{"confIsSame", false, "Pairs of the same particle"};
    Configurable<int> confPDGCodeTrack{"confPDGCodeTrack", 2212, "Particle 2 - PDG code"};
    Configurable<int> confPIDTrack{"confPIDTrack", 2, "Particle 2 - Read from cutCulator"}; // we also need the possibility to specify whether the bit is true/false ->std>>vector<std::pair<int, int>>
    Configurable<int> confTrackSign{"confTrackSign", 1, "Track sign"};
    Configurable<bool> confIsTrackIdentified{"confIsTrackIdentified", true, "Enable PID for the track"};
    Configurable<float> confTrackLowPtCut{"confTrackLowPtCut", 0.5, "Low pT cut of the track"};
    Configurable<float> confTrackHighPtCut{"confTrackHighPtCut", 2.5, "High pT cut of the track"};
    Configurable<float> confTrackEtaMax{"confTrackEtaMax", 0.8, "max. pseudorapidity of the track"};
    Configurable<float> minPtPidTpcTofProton{"minPtPidTpcTofProton", 0.5, "Momentum threshold for change of the PID method (from using TPC to TPC and TOF)."};
    Configurable<float> minPtPidTpcTofPion{"minPtPidTpcTofPion", 0.5, "Momentum threshold for change of the PID method (from using TPC to TPC and TOF)."};
    Configurable<float> minPtPidTpcTofKaonLF{"minPtPidTpcTofKaonLF", 0.5, "Momentum threshold for change of the PID method (from using TPC to TPC and TOF)."};
  } ConfTrack;

  /// Particle 2 --- D0/D0bar meson
  struct : o2::framework::ConfigurableGroup {
    Configurable<int> confPDGCodeD0{"confPDGCodeD0", 421, "D0 meson - PDG code"};
    Configurable<int> confPDGCodeD0bar{"confPDGCodeD0bar", -421, "D0bar meson - PDG code"};
    Configurable<float> confMinPtD0D0bar{"confMinPtD0D0bar", 1.0, "D0/D0bar sel. - min. pT"};
    Configurable<float> confMaxPtD0D0bar{"confMaxPtD0D0bar", 3.0, "D0/D0bar sel. - max. pT"};
    Configurable<uint8_t> confD0OriginFlag{"confD0OriginFlag", 0, "D0/D0bar origin: 0 - prompt, 1 - non-prompt"};
    Configurable<float> confMinPtD0D0barReco{"confMinPtD0D0barReco", 0.5, "MC Reco - D0/D0bar sel. - min. pT"};
    Configurable<float> confMaxPtD0D0barReco{"confMaxPtD0D0barReco", 24.0, "MC Reco - D0/D0bar sel. - max. pT"};
    Configurable<float> minInvMassD0D0barSignal{"minInvMassD0D0barSignal", 1.81, "Min. inv. mass of D0/D0bar for signal region"};
    Configurable<float> maxInvMassD0D0barSignal{"maxInvMassD0D0barSignal", 1.922, "Max. inv. mass of D0/D0bar for signal region"};
    Configurable<float> minInvMassD0D0barLeftSB{"minInvMassD0D0barLeftSB", 1.65, "Min. inv. mass of D0/D0bar for left SB region"};
    Configurable<float> maxInvMassD0D0barLeftSB{"maxInvMassD0D0barLeftSB", 1.754, "Max. inv. mass of D0/D0bar for left SB region"};
    Configurable<float> minInvMassD0D0barRightSB{"minInvMassD0D0barRightSB", 1.978, "Min. inv. mass of D0/D0bar for right SB region"};
    Configurable<float> maxInvMassD0D0barRightSB{"maxInvMassD0D0barRightSB", 2.09, "Max. inv. mass of D0/D0bar for right SB region"};
  } ConfDmesons;

  struct : o2::framework::ConfigurableGroup {
    Configurable<float> confMaxProbMlClass1Bg{"confMaxProbMlClass1Bg", 0.4, "ML: max prob. that D0/D0bar cand. is from the backgound"};
    Configurable<float> confMinProbMlClass2Prompt{"confMinProbMlClass2Prompt", 0.05, "ML: min prob. that D0/D0bar cand. is prompt"};
    Configurable<float> confMaxProbMlClass3NonPrompt{"confMaxProbMlClass3NonPrompt", 1.0, "ML: max prob. that D0/D0bar cand. is non-prompt"};
    Configurable<float> confClass1BgProbStep{"confClass1BgProbStep", 0.05, "ML: prob. step for score class 1"};
    Configurable<float> confClass1BgProbStart{"confClass1BgProbStart", 0.05, "ML: starting prob. value in optimization for score class 1"};
    Configurable<float> confClass2PromptProbStep{"confClass2PromptProbStep", 0.05, "ML: prob. step for score class 2 - prompt"};
    Configurable<float> confClass2PromptProbStart{"confClass2PromptProbStart", 0.1, "ML: starting prob. value in optimization for score class 2"};
    Configurable<float> confClass3NonPromptProbStep{"confClass3NonPromptProbStep", 0.05, "ML: prob. step for score class 2 - non-prompt"};
    Configurable<float> confClass3NonPromptProbStart{"confClass3NonPromptProbStart", 0.05, "ML: starting prob. value in optimization for score class 3"};
  } ConfMlOpt;

  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_d0_to_pi_k::vecBinsPt}, "pT bin limits"};
  Configurable<uint8_t> confChooseD0trackCorr{"confChooseD0trackCorr", 0, "If 0 correlations with D0s, if 1 with D0bars"};

  // Efficiency
  Configurable<bool> doEfficiencyCorr{"doEfficiencyCorr", false, "Apply efficiency corrections"};

  /// Partitions for particle 1
  Partition<FemtoFullParticles> partsTrack = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kTrack)) && (aod::femtouniverseparticle::sign == int8_t(ConfTrack.confTrackSign)) && (aod::femtouniverseparticle::pt > ConfTrack.confTrackLowPtCut) && (aod::femtouniverseparticle::pt < ConfTrack.confTrackHighPtCut) && (nabs(aod::femtouniverseparticle::eta) < ConfTrack.confTrackEtaMax);
  Partition<FemtoMcRecoParticles> partsTrackMCReco = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kTrack)) && (aod::femtouniverseparticle::sign == int8_t(ConfTrack.confTrackSign)) && (aod::femtouniverseparticle::pt > ConfTrack.confTrackLowPtCut) && (aod::femtouniverseparticle::pt < ConfTrack.confTrackHighPtCut) && (nabs(aod::femtouniverseparticle::eta) < ConfTrack.confTrackEtaMax);
  Partition<FemtoMcRecoParticles> partsTrackMCTruth = (aod::femtouniverseparticle::partType == static_cast<uint8_t>(aod::femtouniverseparticle::ParticleType::kMCTruthTrack)) && (aod::femtouniverseparticle::pidCut == static_cast<uint32_t>(ConfTrack.confPDGCodeTrack)) && (aod::femtouniverseparticle::pt > ConfTrack.confTrackLowPtCut) && (aod::femtouniverseparticle::pt < ConfTrack.confTrackHighPtCut) && (nabs(aod::femtouniverseparticle::eta) < ConfTrack.confTrackEtaMax);

  /// Partitions for particle 2
  /// Partition with all D0/D0bar mesons (which pass one and double mass hypothesis)
  Partition<FemtoFullParticles> partsAllDmesons = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kD0)) && ((aod::femtouniverseparticle::mLambda > 0.0f) || (aod::femtouniverseparticle::mAntiLambda > 0.0f)) && (aod::femtouniverseparticle::decayVtxX < ConfMlOpt.confMaxProbMlClass1Bg) && (aod::femtouniverseparticle::decayVtxY > ConfMlOpt.confMinProbMlClass2Prompt);
  /// Partition with D0 mesons only (one and double mass hypothesis)
  Partition<FemtoFullParticles> partsAllD0s = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kD0)) && (aod::femtouniverseparticle::mLambda > ConfDmesons.minInvMassD0D0barSignal) && (aod::femtouniverseparticle::mLambda < ConfDmesons.maxInvMassD0D0barSignal) && (aod::femtouniverseparticle::pt > ConfDmesons.confMinPtD0D0bar) && (aod::femtouniverseparticle::pt < ConfDmesons.confMaxPtD0D0bar);
  /// Partition with D0 mesons only (one mass hypothesis)
  Partition<FemtoFullParticles> partsD0s = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kD0)) && (aod::femtouniverseparticle::mLambda > ConfDmesons.minInvMassD0D0barSignal) && (aod::femtouniverseparticle::mLambda < ConfDmesons.maxInvMassD0D0barSignal) && (aod::femtouniverseparticle::mAntiLambda < 0.0f) && (aod::femtouniverseparticle::pt > ConfDmesons.confMinPtD0D0bar) && (aod::femtouniverseparticle::pt < ConfDmesons.confMaxPtD0D0bar) && (aod::femtouniverseparticle::decayVtxX < ConfMlOpt.confMaxProbMlClass1Bg) && (aod::femtouniverseparticle::decayVtxY > ConfMlOpt.confMinProbMlClass2Prompt);
  /// Partition with D0s selected from the side-band (SB) regions (candidates with double mass hypothesis included)
  Partition<FemtoFullParticles> partsD0sFromSB = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kD0)) && ((aod::femtouniverseparticle::mLambda > ConfDmesons.minInvMassD0D0barLeftSB && aod::femtouniverseparticle::mLambda < ConfDmesons.maxInvMassD0D0barLeftSB) || (aod::femtouniverseparticle::mLambda > ConfDmesons.minInvMassD0D0barRightSB && aod::femtouniverseparticle::mLambda < ConfDmesons.maxInvMassD0D0barRightSB)) && (aod::femtouniverseparticle::pt > ConfDmesons.confMinPtD0D0bar) && (aod::femtouniverseparticle::pt < ConfDmesons.confMaxPtD0D0bar);
  /// Partition with D0bar mesons only (one and double mass hypothesis)
  Partition<FemtoFullParticles> partsAllD0bars = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kD0)) && (aod::femtouniverseparticle::mAntiLambda > ConfDmesons.minInvMassD0D0barSignal) && (aod::femtouniverseparticle::mAntiLambda < ConfDmesons.maxInvMassD0D0barSignal) && (aod::femtouniverseparticle::pt > ConfDmesons.confMinPtD0D0bar) && (aod::femtouniverseparticle::pt < ConfDmesons.confMaxPtD0D0bar);
  /// Partition with D0bar mesons only (one mass hypothesis)
  Partition<FemtoFullParticles> partsD0bars = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kD0)) && (aod::femtouniverseparticle::mLambda < 0.0f) && (aod::femtouniverseparticle::mAntiLambda > ConfDmesons.minInvMassD0D0barSignal) && (aod::femtouniverseparticle::mAntiLambda < ConfDmesons.maxInvMassD0D0barSignal) && (aod::femtouniverseparticle::pt > ConfDmesons.confMinPtD0D0bar) && (aod::femtouniverseparticle::pt < ConfDmesons.confMaxPtD0D0bar) && (aod::femtouniverseparticle::decayVtxX < ConfMlOpt.confMaxProbMlClass1Bg) && (aod::femtouniverseparticle::decayVtxY > ConfMlOpt.confMinProbMlClass2Prompt);
  /// Partition with D0bars selected from the side-band (SB) regions (candidates with double mass hypothesis included)
  Partition<FemtoFullParticles> partsD0barsFromSB = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kD0)) && ((aod::femtouniverseparticle::mAntiLambda > ConfDmesons.minInvMassD0D0barLeftSB && aod::femtouniverseparticle::mAntiLambda < ConfDmesons.maxInvMassD0D0barLeftSB) || (aod::femtouniverseparticle::mAntiLambda > ConfDmesons.minInvMassD0D0barRightSB && aod::femtouniverseparticle::mAntiLambda < ConfDmesons.maxInvMassD0D0barRightSB)) && (aod::femtouniverseparticle::pt > ConfDmesons.confMinPtD0D0bar) && (aod::femtouniverseparticle::pt < ConfDmesons.confMaxPtD0D0bar);
  /// Partition for MC Reco prompt D0/D0bar mesons
  Partition<FemtoMcRecoParticles> partsPromptD0MCReco = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kD0)) && (aod::femtouniverseparticle::mLambda > ConfDmesons.minInvMassD0D0barSignal) && (aod::femtouniverseparticle::mLambda < ConfDmesons.maxInvMassD0D0barSignal) && (aod::femtouniverseparticle::pt > ConfDmesons.confMinPtD0D0bar) && (aod::femtouniverseparticle::pt < ConfDmesons.confMaxPtD0D0bar) && (aod::femtouniverseparticle::decayVtxX < ConfMlOpt.confMaxProbMlClass1Bg) && (aod::femtouniverseparticle::decayVtxY > ConfMlOpt.confMinProbMlClass2Prompt) && (aod::femtouniverseparticle::tpcNClsFound == ConfDmesons.confD0OriginFlag);
  Partition<FemtoMcRecoParticles> partsPromptD0barMCReco = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kD0)) && (aod::femtouniverseparticle::mAntiLambda > ConfDmesons.minInvMassD0D0barSignal) && (aod::femtouniverseparticle::mAntiLambda < ConfDmesons.maxInvMassD0D0barSignal) && (aod::femtouniverseparticle::pt > ConfDmesons.confMinPtD0D0bar) && (aod::femtouniverseparticle::pt < ConfDmesons.confMaxPtD0D0bar) && (aod::femtouniverseparticle::decayVtxX < ConfMlOpt.confMaxProbMlClass1Bg) && (aod::femtouniverseparticle::decayVtxY > ConfMlOpt.confMinProbMlClass2Prompt) && (aod::femtouniverseparticle::tpcNClsFound == ConfDmesons.confD0OriginFlag);
  /// Partition for MC Reco prompt D0/D0bar mesons (sideband candidates)
  Partition<FemtoMcRecoParticles> partsPromptD0MCRecoSB = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kD0)) && ((aod::femtouniverseparticle::mLambda > ConfDmesons.minInvMassD0D0barLeftSB && aod::femtouniverseparticle::mLambda < ConfDmesons.maxInvMassD0D0barLeftSB) || (aod::femtouniverseparticle::mLambda > ConfDmesons.minInvMassD0D0barRightSB && aod::femtouniverseparticle::mLambda < ConfDmesons.maxInvMassD0D0barRightSB)) && (aod::femtouniverseparticle::pt > ConfDmesons.confMinPtD0D0bar) && (aod::femtouniverseparticle::pt < ConfDmesons.confMaxPtD0D0bar) && (aod::femtouniverseparticle::decayVtxX < ConfMlOpt.confMaxProbMlClass1Bg) && (aod::femtouniverseparticle::decayVtxY > ConfMlOpt.confMinProbMlClass2Prompt) && (aod::femtouniverseparticle::tpcNClsFound == ConfDmesons.confD0OriginFlag);
  Partition<FemtoMcRecoParticles> partsPromptD0barMCRecoSB = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kD0)) && ((aod::femtouniverseparticle::mAntiLambda > ConfDmesons.minInvMassD0D0barLeftSB && aod::femtouniverseparticle::mAntiLambda < ConfDmesons.maxInvMassD0D0barLeftSB) || (aod::femtouniverseparticle::mAntiLambda > ConfDmesons.minInvMassD0D0barRightSB && aod::femtouniverseparticle::mAntiLambda < ConfDmesons.maxInvMassD0D0barRightSB)) && (aod::femtouniverseparticle::pt > ConfDmesons.confMinPtD0D0bar) && (aod::femtouniverseparticle::pt < ConfDmesons.confMaxPtD0D0bar) && (aod::femtouniverseparticle::decayVtxX < ConfMlOpt.confMaxProbMlClass1Bg) && (aod::femtouniverseparticle::decayVtxY > ConfMlOpt.confMinProbMlClass2Prompt) && (aod::femtouniverseparticle::tpcNClsFound == ConfDmesons.confD0OriginFlag);
  // MC Truth D0/D0bar candidates
  Partition<FemtoMcRecoParticles> partsD0MCTruth = (aod::femtouniverseparticle::partType == static_cast<uint8_t>(aod::femtouniverseparticle::ParticleType::kMCTruthTrack)) && (aod::femtouniverseparticle::pidCut == static_cast<uint32_t>(ConfDmesons.confPDGCodeD0)) && (aod::femtouniverseparticle::pt > ConfDmesons.confMinPtD0D0bar) && (aod::femtouniverseparticle::pt < ConfDmesons.confMaxPtD0D0bar);
  Partition<FemtoMcRecoParticles> partsD0barMCTruth = (aod::femtouniverseparticle::partType == static_cast<uint8_t>(aod::femtouniverseparticle::ParticleType::kMCTruthTrack)) && (aod::femtouniverseparticle::pidCut == static_cast<uint32_t>(ConfDmesons.confPDGCodeD0bar)) && (aod::femtouniverseparticle::pt > ConfDmesons.confMinPtD0D0bar) && (aod::femtouniverseparticle::pt < ConfDmesons.confMaxPtD0D0bar);
  /// Partition for D0/D0bar daughters
  Partition<FemtoFullParticles> partsDmesonsChildren = aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kD0Child);

  /// Histogramming for particle 1
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kTrack, 1> trackHistoPartTrack;

  /// Histogramming for particle 2
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kD0, 0> trackHistoPartD0D0bar;

  /// Histogramming for Event
  FemtoUniverseEventHisto eventHisto;

  /// The configurables need to be passed to an std::vector
  int vPIDTrack;
  std::vector<float> kNsigma;

  /// particle part
  ConfigurableAxis confTempFitVarBins{"confTempFitVarBins", {300, -0.15, 0.15}, "binning of the TempFitVar in the pT vs. TempFitVar plot"};
  ConfigurableAxis confTempFitVarInvMassBins{"confTempFitVarInvMassBins", {6000, 0.9, 4.0}, "binning of the TempFitVar in the pT vs. TempFitVar plot"};
  ConfigurableAxis confTempFitVarpTBins{"confTempFitVarpTBins", {20, 0.5, 4.05}, "pT binning of the pT vs. TempFitVar plot"};
  ConfigurableAxis confTempFitVarDCABins{"confTempFitVarDCABins", {500, -1.0, 1.0}, "binning of the DCA histograms"};

  /// Correlation part
  ConfigurableAxis confMultBins{"confMultBins", {VARIABLE_WIDTH, 0.0f, 4.0f, 8.0f, 12.0f, 16.0f, 20.0f, 24.0f, 28.0f, 32.0f, 36.0f, 40.0f, 44.0f, 48.0f, 52.0f, 56.0f, 60.0f, 64.0f, 68.0f, 72.0f, 76.0f, 80.0f, 84.0f, 88.0f, 92.0f, 96.0f, 100.0f, 200.0f, 99999.f}, "Mixing bins - multiplicity"}; // \todo to be obtained from the hash task
  // ConfigurableAxis confMultBins{"confMultBins", {VARIABLE_WIDTH, 0.0f, 20.0f, 40.0f, 60.0f, 80.0f, 100.0f, 200.0f, 99999.f}, "Mixing bins - multiplicity"};
  ConfigurableAxis confVtxBins{"confVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis confmTBins3D{"confmTBins3D", {VARIABLE_WIDTH, 1.02f, 1.14f, 1.20f, 1.26f, 1.38f, 1.56f, 1.86f, 4.50f}, "mT Binning for the 3Dimensional plot: k* vs multiplicity vs mT (set <<ConfBothTracks.confUse3D>> to true in order to use)"};
  ConfigurableAxis confmultBins3D{"confmultBins3D", {VARIABLE_WIDTH, 0.0f, 20.0f, 30.0f, 40.0f, 99999.0f}, "multiplicity Binning for the 3Dimensional plot: k* vs multiplicity vs mT (set <<ConfBothTracks.confUse3D>> to true in order to use)"};

  ColumnBinningPolicy<aod::collision::PosZ, aod::femtouniversecollision::MultNtr> colBinning{{confVtxBins, confMultBins}, true};

  ConfigurableAxis confkstarBins{"confkstarBins", {1500, 0., 6.}, "binning kstar"};
  ConfigurableAxis confkTBins{"confkTBins", {150, 0., 9.}, "binning kT"};
  ConfigurableAxis confmTBins{"confmTBins", {225, 0., 7.5}, "binning mT"};
  ConfigurableAxis confPtBins{"confPtBins", {360, 0., 36.}, "binning pT"};
  ConfigurableAxis confInvMassBins{"confInvMassBins", {500, 0., 5.0}, "binning inv. mass"};

  Configurable<bool> confIsCPR{"confIsCPR", true, "Close Pair Rejection"};
  Configurable<bool> confCPRPlotPerRadii{"confCPRPlotPerRadii", false, "Plot CPR per radii"};
  Configurable<float> confCPRdeltaPhiCutMax{"confCPRdeltaPhiCutMax", 0.0, "Delta Phi max cut for Close Pair Rejection"};
  Configurable<float> confCPRdeltaPhiCutMin{"confCPRdeltaPhiCutMin", 0.0, "Delta Phi min cut for Close Pair Rejection"};
  Configurable<float> confCPRdeltaEtaCutMax{"confCPRdeltaEtaCutMax", 0.0, "Delta Eta max cut for Close Pair Rejection"};
  Configurable<float> confCPRdeltaEtaCutMin{"confCPRdeltaEtaCutMin", 0.0, "Delta Eta min cut for Close Pair Rejection"};
  Configurable<float> confCPRChosenRadii{"confCPRChosenRadii", 0.80, "Delta Eta cut for Close Pair Rejection"};

  Configurable<bool> applyMLOpt{"applyMLOpt", false, "Enable for ML selection optimization"};
  Configurable<bool> confRemoveSoftPions{"confRemoveSoftPions", false, "Enable to remove soft pions from D* decays"};
  Configurable<bool> confSoftPionD0Flag{"confSoftPionD0Flag", false, "Enable soft pion check for D0s"};
  Configurable<bool> confSoftPionD0barFlag{"confSoftPionD0barFlag", false, "Enable soft pion check for D0bars"};
  Configurable<float> sigmaSoftPiInvMass{"sigmaSoftPiInvMass", 0.1, "Sigma value from the inv. mass fit for soft pions"};
  // Event mixing configurables
  Configurable<int> confNEventsMix{"confNEventsMix", 5, "Number of events for mixing"};

  FemtoUniverseContainer<femto_universe_container::EventType::same, femto_universe_container::Observable::kstar> sameEventAngularCont;
  FemtoUniverseContainer<femto_universe_container::EventType::mixed, femto_universe_container::Observable::kstar> mixedEventAngularCont;
  FemtoUniversePairCleaner<aod::femtouniverseparticle::ParticleType::kTrack, aod::femtouniverseparticle::ParticleType::kD0> pairCleaner;
  FemtoUniverseDetaDphiStar<aod::femtouniverseparticle::ParticleType::kTrack, aod::femtouniverseparticle::ParticleType::kD0> pairCloseRejection;
  FemtoUniverseSoftPionRemoval<aod::femtouniverseparticle::ParticleType::kTrack, aod::femtouniverseparticle::ParticleType::kD0> softPionRemoval;
  FemtoUniverseTrackSelection trackCuts;
  // Axes for BDT score classes' histograms
  AxisSpec axisBdtScore{100, 0.f, 1.f};
  AxisSpec axisSelStatus{2, -0.5f, 1.5f};

  /// Histogram output
  HistogramRegistry qaRegistry{"TrackQA", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry resultRegistry{"Correlations", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry mixQaRegistry{"mixQaRegistry", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry mcRecoRegistry{"McRecoHistos", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry mcTruthRegistry{"McGenHistos", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry registryDCA{"registryDCA", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};

  // Efficiency
  EfficiencyConfigurableGroup effConfGroup;
  EfficiencyCalculator<TH1> efficiencyCalculator{&effConfGroup};
  float weight = 1.0;

  HistogramRegistry registry{"registry",
                             {{"hPtD0", "D^{0} cand.;#it{p}_{T} (GeV/#it{c});counts", {HistType::kTH1F, {confPtBins}}},
                              {"hPtD0bar", "#bar{D^{0}};#it{p}_{T} (GeV/#it{c});counts", {HistType::kTH1F, {confPtBins}}},
                              {"hPtD0D0bar", "#bar{D^{0}};#it{p}_{T} (GeV/#it{c});counts", {HistType::kTH1F, {confPtBins}}},
                              {"hPhiD0D0bar", ";#varphi (rad);counts", {HistType::kTH1F, {{80, 0., o2::constants::math::TwoPI}}}},
                              {"hPhiD0", ";#varphi (rad);counts", {HistType::kTH1F, {{80, 0., o2::constants::math::TwoPI}}}},
                              {"hPhiD0bar", ";#varphi (rad);counts", {HistType::kTH1F, {{80, 0., o2::constants::math::TwoPI}}}},
                              {"hEtaD0D0bar", ";#eta ;counts", {HistType::kTH1F, {{200, -1., 1.}}}},
                              {"hEtaD0", ";#eta ;counts", {HistType::kTH1F, {{200, -1., 1.}}}},
                              {"hEtaD0bar", ";#eta ;counts", {HistType::kTH1F, {{200, -1., 1.}}}},
                              {"hDecayLengthD0", ";decay length (cm);counts", {HistType::kTH1F, {{800, 0., 4.}}}},
                              {"hDecayLengthD0bar", ";decay length (cm);counts", {HistType::kTH1F, {{800, 0., 4.}}}},
                              {"hPtDaughters", ";#it{p}_{T} (GeV/#it{c});counts", {HistType::kTH1F, {{300, 0., 12.}}}},
                              {"hSignDaughters", ";sign ;counts", {HistType::kTH1F, {{10, -2.5, 2.5}}}},
                              {"hDCAxyDaughters", "; #it{DCA}_{xy} (cm); counts", {HistType::kTH1F, {{140, 0., 0.14}}}},
                              {"hDCAzDaughters", "; #it{DCA}_{z} (cm); counts", {HistType::kTH1F, {{140, 0., 0.14}}}}}};

  // PID for protons
  bool isProtonNSigma(float mom, float nsigmaTPCPr, float nsigmaTOFPr) // previous version from: https://github.com/alisw/AliPhysics/blob/master/PWGCF/FEMTOSCOPY/AliFemtoUser/AliFemtoMJTrackCut.cxx
  {
    if (mom < ConfTrack.minPtPidTpcTofProton) {
      if (std::abs(nsigmaTPCPr) < ConfBothTracks.confNsigmaTPCProton) {
        return true;
      } else {
        return false;
      }
    } else if (mom > ConfTrack.minPtPidTpcTofProton) {
      if (std::hypot(nsigmaTOFPr, nsigmaTPCPr) < ConfBothTracks.confNsigmaCombinedProton) {
        return true;
      } else {
        return false;
      }
    }
    return false;
  }

  bool isKaonNSigmaLF(float mom, float nsigmaTPCK, float nsigmaTOFK)
  {
    if (mom < ConfTrack.minPtPidTpcTofKaonLF) {
      if (std::abs(nsigmaTPCK) < ConfBothTracks.confNsigmaTPCKaon) {
        return true;
      } else {
        return false;
      }
    } else if (mom >= ConfTrack.minPtPidTpcTofKaonLF) {
      if (std::sqrt(nsigmaTPCK * nsigmaTPCK + nsigmaTOFK * nsigmaTOFK) < ConfBothTracks.confNsigmaCombinedKaon) {
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
    // using configurables:
    // confNsigmaTPCPion -> TPC Pion Sigma for a given momentum
    // confNsigmaCombinedPion -> TPC and TOF Pion Sigma (combined) for a given momentum
    if (true) {
      if (mom < ConfTrack.minPtPidTpcTofPion) {
        if (std::abs(nsigmaTPCPi) < ConfBothTracks.confNsigmaTPCPion) {
          return true;
        } else {
          return false;
        }
      } else if (mom > ConfTrack.minPtPidTpcTofPion) {
        if (std::hypot(nsigmaTOFPi, nsigmaTPCPi) < ConfBothTracks.confNsigmaCombinedPion) {
          return true;
        } else {
          return false;
        }
      }
    }
    return false;
  }

  bool isParticleNSigma(float mom, float nsigmaTPCPr, float nsigmaTOFPr, float nsigmaTPCPi, float nsigmaTOFPi, float nsigmaTPCK, float nsigmaTOFK)
  {
    switch (ConfTrack.confPDGCodeTrack) {
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
        return isKaonNSigmaLF(mom, nsigmaTPCK, nsigmaTOFK);
        break;
      default:
        return false;
    }
  }

  void init(InitContext&)
  {
    efficiencyCalculator.init();

    eventHisto.init(&qaRegistry);
    qaRegistry.add("QA_D0D0barSelection/hInvMassD0", ";#it{M}(K^{-}#pi^{+}) (GeV/#it{c}^{2});counts", kTH1F, {confInvMassBins});
    qaRegistry.add("QA_D0D0barSelection/hPtD0", "D^{0} cand.;#it{p}_{T} (GeV/#it{c});counts", kTH1F, {confPtBins});
    qaRegistry.add("QA_D0D0barSelection/hInvMassD0bar", ";#it{M}(K^{-}#pi^{+}) (GeV/#it{c}^{2});counts", kTH1F, {confInvMassBins});
    qaRegistry.add("QA_D0D0barSelection/hPtD0bar", "#bar{D^{0}} cand.;#it{p}_{T} (GeV/#it{c});counts", kTH1F, {confPtBins});
    qaRegistry.add("QA_D0D0barSelection_SB/hInvMassD0", ";#it{M}(K^{-}#pi^{+}) (GeV/#it{c}^{2});counts", kTH1F, {confInvMassBins});
    qaRegistry.add("QA_D0D0barSelection_SB/hPtD0", "D^{0} cand.;#it{p}_{T} (GeV/#it{c});counts", kTH1F, {confPtBins});

    qaRegistry.add("D0_pos_daugh/nSigmaTPC", "; #it{p} (GeV/#it{c}); n#sigma_{TPC}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("D0_pos_daugh/nSigmaTOF", "; #it{p} (GeV/#it{c}); n#sigma_{TOF}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("D0_pos_daugh/pt", "; #it{p_T} (GeV/#it{c}); Counts", kTH1F, {{100, 0, 10}});
    qaRegistry.add("D0_pos_daugh/eta", "; #it{eta}; Counts", kTH1F, {{200, -1.5, 1.5}});
    qaRegistry.add("D0_pos_daugh/phi", "; #it{varphi}; Counts", kTH1F, {{200, 0, o2::constants::math::TwoPI}});
    qaRegistry.add("D0_pos_daugh/hDCAxy", "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", kTH2F, {{100, 0, 10}, {500, -5, 5}});

    qaRegistry.add("D0_neg_daugh/nSigmaTPC", "; #it{p} (GeV/#it{c}); n#sigma_{TPC}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("D0_neg_daugh/nSigmaTOF", "; #it{p} (GeV/#it{c}); n#sigma_{TOF}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("D0_neg_daugh/pt", "; #it{p_T} (GeV/#it{c}); Counts", kTH1F, {{100, 0, 10}});
    qaRegistry.add("D0_neg_daugh/eta", "; #it{eta}; Counts", kTH1F, {{200, -1.5, 1.5}});
    qaRegistry.add("D0_neg_daugh/phi", "; #it{varphi}; Counts", kTH1F, {{200, 0, o2::constants::math::TwoPI}});
    qaRegistry.add("D0_neg_daugh/hDCAxy", "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", kTH2F, {{100, 0, 10}, {500, -5, 5}});

    qaRegistry.add("D0bar_pos_daugh/nSigmaTPC", "; #it{p} (GeV/#it{c}); n#sigma_{TPC}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("D0bar_pos_daugh/nSigmaTOF", "; #it{p} (GeV/#it{c}); n#sigma_{TOF}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("D0bar_pos_daugh/pt", "; #it{p_T} (GeV/#it{c}); Counts", kTH1F, {{100, 0, 10}});
    qaRegistry.add("D0bar_pos_daugh/eta", "; #it{eta}; Counts", kTH1F, {{200, -1.5, 1.5}});
    qaRegistry.add("D0bar_pos_daugh/phi", "; #it{varphi}; Counts", kTH1F, {{200, 0, o2::constants::math::TwoPI}});
    qaRegistry.add("D0bar_pos_daugh/hDCAxy", "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", kTH2F, {{100, 0, 10}, {500, -5, 5}});

    qaRegistry.add("D0bar_neg_daugh/nSigmaTPC", "; #it{p} (GeV/#it{c}); n#sigma_{TPC}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("D0bar_neg_daugh/nSigmaTOF", "; #it{p} (GeV/#it{c}); n#sigma_{TOF}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("D0bar_neg_daugh/pt", "; #it{p_T} (GeV/#it{c}); Counts", kTH1F, {{100, 0, 10}});
    qaRegistry.add("D0bar_neg_daugh/eta", "; #it{eta}; Counts", kTH1F, {{200, -1.5, 1.5}});
    qaRegistry.add("D0bar_neg_daugh/phi", "; #it{varphi}; Counts", kTH1F, {{200, 0, o2::constants::math::TwoPI}});
    qaRegistry.add("D0bar_neg_daugh/hDCAxy", "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", kTH2F, {{100, 0, 10}, {500, -5, 5}});

    qaRegistry.add("Hadron/nSigmaTPCPr", "; #it{p} (GeV/#it{c}); n#sigma_{TPCPr}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("Hadron/nSigmaTOFPr", "; #it{p} (GeV/#it{c}); n#sigma_{TOFPr}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("Hadron/nSigmaTPCPi", "; #it{p} (GeV/#it{c}); n#sigma_{TPCPi}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("Hadron/nSigmaTOFPi", "; #it{p} (GeV/#it{c}); n#sigma_{TOFPi}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("Hadron/nSigmaTPCKa", "; #it{p} (GeV/#it{c}); n#sigma_{TPCKa}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("Hadron/nSigmaTOFKa", "; #it{p} (GeV/#it{c}); n#sigma_{TOFKa}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});

    // D0/D0bar histograms
    auto vbins = (std::vector<double>)binsPt;
    if (doEfficiencyCorr) {
      registry.add("D0D0bar_oneMassHypo/hMassVsPtEffCorr", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {confInvMassBins, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("D0D0bar_oneMassHypo/hMassVsPtD0EffCorr", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {confInvMassBins, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("D0D0bar_oneMassHypo/hMassVsPtD0barEffCorr", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {confInvMassBins, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("D0D0bar_doubleMassHypo/hMassVsPtEffCorr", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {confInvMassBins, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("D0D0bar_doubleMassHypo/hMassVsPtD0EffCorr", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {confInvMassBins, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("D0D0bar_doubleMassHypo/hMassVsPtD0barEffCorr", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {confInvMassBins, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    } else {
      registry.add("D0D0bar_oneMassHypo/hMassVsPt", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {confInvMassBins, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("D0D0bar_oneMassHypo/hMassVsPtD0", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {confInvMassBins, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("D0D0bar_oneMassHypo/hMassVsPtD0bar", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {confInvMassBins, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("D0D0bar_doubleMassHypo/hMassVsPt", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {confInvMassBins, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("D0D0bar_doubleMassHypo/hMassVsPtD0", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {confInvMassBins, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("D0D0bar_doubleMassHypo/hMassVsPtD0bar", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {confInvMassBins, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    }
    // Histograms for BDT score classes' check
    registry.add("DebugBdt/hBdtScore1VsStatus", ";BDT score;status", {HistType::kTH2F, {axisBdtScore, axisSelStatus}});
    registry.add("DebugBdt/hBdtScore2VsStatus", ";BDT score;status", {HistType::kTH2F, {axisBdtScore, axisSelStatus}});
    registry.add("DebugBdt/hBdtScore3VsStatus", ";BDT score;status", {HistType::kTH2F, {axisBdtScore, axisSelStatus}});
    if (applyMLOpt) {
      registry.add("D0D0bar_MLSel/hMassVsPt1", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {confInvMassBins, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("D0D0bar_MLSel/hMassVsPt2", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {confInvMassBins, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("D0D0bar_MLSel/hMassVsPt3", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {confInvMassBins, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("D0D0bar_MLSel/hMassVsPt4", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {confInvMassBins, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("D0D0bar_MLSel/hMassVsPt5", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {confInvMassBins, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("D0D0bar_MLSel/hMassVsPt6", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {confInvMassBins, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("D0D0bar_MLSel/hMassVsPt7", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {confInvMassBins, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("D0D0bar_MLSel/hMassVsPt8", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {confInvMassBins, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("D0D0bar_MLSel/hMassVsPt9", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {confInvMassBins, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("D0D0bar_MLSel/hMassVsPt10", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {confInvMassBins, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    }
    // MC Reco
    mcRecoRegistry.add("hMcRecD0", "MC Reco all D0s;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {400, -1.0, 1.0}}});
    mcRecoRegistry.add("hMcRecD0Pt", "MC Reco all D0s;#it{p}_{T} (GeV/c); counts", {HistType::kTH1F, {vbins}});
    mcRecoRegistry.add("hMcRecD0Prompt", "MC Reco prompt D0s;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {400, -1.0, 1.0}}});
    mcRecoRegistry.add("hMcRecD0PromptPt", "MC Reco prompt D0s;#it{p}_{T} (GeV/c); counts", {HistType::kTH1F, {vbins}});
    mcRecoRegistry.add("hMcRecD0NonPrompt", "MC Reco non-prompt D0s;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {400, -1.0, 1.0}}});
    mcRecoRegistry.add("hMcRecD0NonPromptPt", "MC Reco non-prompt D0s;#it{p}_{T} (GeV/c); counts", {HistType::kTH1F, {vbins}});
    mcRecoRegistry.add("hMcRecD0Phi", "MC Reco all D0s;#varphi (rad); counts", {HistType::kTH1F, {{80, 0., o2::constants::math::TwoPI}}});
    mcRecoRegistry.add("hMcRecD0PromptPhi", "MC Reco prompt D0s;#varphi (rad); counts", {HistType::kTH1F, {{80, 0., o2::constants::math::TwoPI}}});
    mcRecoRegistry.add("hMcRecD0NonPromptPhi", "MC Reco non-prompt D0s;#varphi (rad); counts", {HistType::kTH1F, {{80, 0., o2::constants::math::TwoPI}}});
    mcRecoRegistry.add("hMcRecD0bar", "MC Reco all D0bars;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {400, -1.0, 1.0}}});
    mcRecoRegistry.add("hMcRecD0barPt", "MC Reco all D0bars;#it{p}_{T} (GeV/c); counts", {HistType::kTH1F, {vbins}});
    mcRecoRegistry.add("hMcRecD0barPrompt", "MC Reco prompt D0bars;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {400, -1.0, 1.0}}});
    mcRecoRegistry.add("hMcRecD0barPromptPt", "MC Reco prompt D0bars;#it{p}_{T} (GeV/c); counts", {HistType::kTH1F, {vbins}});
    mcRecoRegistry.add("hMcRecD0barNonPrompt", "MC Reco non-prompt D0bars;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {400, -1.0, 1.0}}});
    mcRecoRegistry.add("hMcRecD0barNonPromptPt", "MC Reco non-prompt D0bars;#it{p}_{T} (GeV/c); counts", {HistType::kTH1F, {vbins}});
    mcRecoRegistry.add("hMcRecD0barPhi", "MC Reco all D0bars;#varphi (rad); counts", {HistType::kTH1F, {{80, 0., o2::constants::math::TwoPI}}});
    // Inv. mass histograms
    mcRecoRegistry.add("hMassVsPtD0Sig", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {confInvMassBins, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    mcRecoRegistry.add("hMassVsPtD0Refl", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {confInvMassBins, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    mcRecoRegistry.add("hMassVsPtD0Bkg", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {confInvMassBins, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    mcRecoRegistry.add("hMassVsPtD0Prompt", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {confInvMassBins, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    mcRecoRegistry.add("hMassVsPtD0NonPrompt", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {confInvMassBins, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    mcRecoRegistry.add("hMassVsPtD0barSig", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {confInvMassBins, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    mcRecoRegistry.add("hMassVsPtD0barRefl", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {confInvMassBins, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    mcRecoRegistry.add("hMassVsPtD0barBkg", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {confInvMassBins, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    mcRecoRegistry.add("hMassVsPtD0barPrompt", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {confInvMassBins, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    mcRecoRegistry.add("hMassVsPtD0barNonPrompt", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {confInvMassBins, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    // Histograms for identified hadrons
    mcRecoRegistry.add("hMcRecKpPt", "MC Reco K+;#it{p}_{T} (GeV/c); counts", {HistType::kTH1F, {{500, 0, 5}}});
    mcRecoRegistry.add("hMcRecKpPtGenEtaGen", "MC Reco K+;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 5}, {400, -1.0, 1.0}}});
    mcRecoRegistry.add("hMcRecKmPt", "MC Reco K-;#it{p}_{T} (GeV/c); counts", {HistType::kTH1F, {{500, 0, 5}}});
    mcRecoRegistry.add("hMcRecKmPtGenEtaGen", "MC Reco K-;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 5}, {400, -1.0, 1.0}}});
    mcRecoRegistry.add("hMcRecPipPt", "MC Reco #pi+;#it{p}_{T} (GeV/c); counts", {HistType::kTH1F, {{500, 0, 5}}});
    mcRecoRegistry.add("hMcRecPipPtGenEtaGen", "MC Reco #pi+;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 5}, {400, -1.0, 1.0}}});
    mcRecoRegistry.add("hMcRecPimPt", "MC Reco #pi-;#it{p}_{T} (GeV/c); counts", {HistType::kTH1F, {{500, 0, 5}}});
    mcRecoRegistry.add("hMcRecPimPtGenEtaGen", "MC Reco #pi-;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 5}, {400, -1.0, 1.0}}});
    mcRecoRegistry.add("hMcRecPrPt", "MC Reco proton;#it{p}_{T} (GeV/c); counts", {HistType::kTH1F, {{500, 0, 5}}});
    mcRecoRegistry.add("hMcRecPrPtGenEtaGen", "MC Reco proton;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 5}, {400, -1.0, 1.0}}});
    mcRecoRegistry.add("hMcRecAntiPrPt", "MC Reco antiproton;#it{p}_{T} (GeV/c); counts", {HistType::kTH1F, {{500, 0, 5}}});
    mcRecoRegistry.add("hMcRecAntiPrPtGenEtaGen", "MC Reco antiproton;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 5}, {400, -1.0, 1.0}}});

    // MC truth
    mcTruthRegistry.add("hMcGenD0", "MC Truth all D0s;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {400, -1.0, 1.0}}});
    mcTruthRegistry.add("hMcGenD0Pt", "MC Truth all D0s;#it{p}_{T} (GeV/c); counts", {HistType::kTH1F, {vbins}});
    mcTruthRegistry.add("hMcGenD0Prompt", "MC Truth prompt D0s;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {400, -1.0, 1.0}}});
    mcTruthRegistry.add("hMcGenD0PromptPt", "MC Truth prompt D0s;#it{p}_{T} (GeV/c); counts", {HistType::kTH1F, {vbins}});
    mcTruthRegistry.add("hMcGenD0NonPrompt", "MC Truth non-prompt D0s;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {400, -1.0, 1.0}}});
    mcTruthRegistry.add("hMcGenD0NonPromptPt", "MC Truth non-prompt D0s;#it{p}_{T} (GeV/c); counts", {HistType::kTH1F, {vbins}});
    mcTruthRegistry.add("hMcGenD0bar", "MC Truth all D0bars;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {400, -1.0, 1.0}}});
    mcTruthRegistry.add("hMcGenD0barPt", "MC Truth all D0bars;#it{p}_{T} (GeV/c); counts", {HistType::kTH1F, {vbins}});
    mcTruthRegistry.add("hMcGenD0barPrompt", "MC Truth prompt D0bars;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {400, -1.0, 1.0}}});
    mcTruthRegistry.add("hMcGenD0barPromptPt", "MC Truth prompt D0s;#it{p}_{T} (GeV/c); counts", {HistType::kTH1F, {vbins}});
    mcTruthRegistry.add("hMcGenD0barNonPrompt", "MC Truth non-prompt D0bars;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {400, -1.0, 1.0}}});
    mcTruthRegistry.add("hMcGenD0barNonPromptPt", "MC Truth non-prompt D0s;#it{p}_{T} (GeV/c); counts", {HistType::kTH1F, {vbins}});
    mcTruthRegistry.add("hMcGenAllPositivePt", "MC Truth all positive;#it{p}_{T} (GeV/c); counts", {HistType::kTH1F, {{360, 0, 36}}});
    mcTruthRegistry.add("hMcGenAllNegativePt", "MC Truth all negative;#it{p}_{T} (GeV/c); counts", {HistType::kTH1F, {{360, 0, 36}}});
    mcTruthRegistry.add("hMcGenKpPtVsEta", "MC Truth K+;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 5}, {400, -1.0, 1.0}}});
    mcTruthRegistry.add("hMcGenKmPtVsEta", "MC Truth K-;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 5}, {400, -1.0, 1.0}}});
    mcTruthRegistry.add("hMcGenPipPtVsEta", "MC Truth #pi+;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 5}, {400, -1.0, 1.0}}});
    mcTruthRegistry.add("hMcGenPimPtVsEta", "MC Truth #pi-;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 5}, {400, -1.0, 1.0}}});
    mcTruthRegistry.add("hMcGenPrPtVsEta", "MC Truth proton;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 5}, {400, -1.0, 1.0}}});
    mcTruthRegistry.add("hMcGenAntiPrPtVsEta", "MC Truth antiproton;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 5}, {400, -1.0, 1.0}}});
    mcTruthRegistry.add("hMcGenKpPt", "MC Truth K+;#it{p}_{T} (GeV/c); counts", {HistType::kTH1F, {{500, 0, 5}}});
    mcTruthRegistry.add("hMcGenKmPt", "MC Truth K-;#it{p}_{T} (GeV/c); counts", {HistType::kTH1F, {{500, 0, 5}}});
    mcTruthRegistry.add("hMcGenPipPt", "MC Truth #pi+;#it{p}_{T} (GeV/c); counts", {HistType::kTH1F, {{500, 0, 5}}});
    mcTruthRegistry.add("hMcGenPimPt", "MC Truth #pi-;#it{p}_{T} (GeV/c); counts", {HistType::kTH1F, {{500, 0, 5}}});
    mcTruthRegistry.add("hMcGenPrPt", "MC Truth proton;#it{p}_{T} (GeV/c); counts", {HistType::kTH1F, {{500, 0, 5}}});
    mcTruthRegistry.add("hMcGenAntiPrPt", "MC Truth antiproton;#it{p}_{T} (GeV/c); counts", {HistType::kTH1F, {{500, 0, 5}}});

    trackHistoPartTrack.init(&qaRegistry, confTempFitVarpTBins, confTempFitVarBins, ConfBothTracks.confIsMC, ConfTrack.confPDGCodeTrack);
    trackHistoPartD0D0bar.init(&qaRegistry, confTempFitVarpTBins, confTempFitVarInvMassBins, ConfBothTracks.confIsMC, ConfDmesons.confPDGCodeD0);

    mixQaRegistry.add("MixingQA/hSECollisionBins", ";bin;Entries", kTH1F, {{120, -0.5, 119.5}});
    mixQaRegistry.add("MixingQA/hMECollisionBins", ";bin;Entries", kTH1F, {{120, -0.5, 119.5}});

    sameEventAngularCont.init(&resultRegistry, confkstarBins, confMultBins, confkTBins, confmTBins, confmultBins3D, confmTBins3D, ConfBothTracks.confEtaBins, ConfBothTracks.confPhiBins, ConfBothTracks.confIsMC, ConfBothTracks.confUse3D);
    mixedEventAngularCont.init(&resultRegistry, confkstarBins, confMultBins, confkTBins, confmTBins, confmultBins3D, confmTBins3D, ConfBothTracks.confEtaBins, ConfBothTracks.confPhiBins, ConfBothTracks.confIsMC, ConfBothTracks.confUse3D);

    sameEventAngularCont.setPDGCodes(ConfDmesons.confPDGCodeD0, ConfTrack.confPDGCodeTrack);
    mixedEventAngularCont.setPDGCodes(ConfDmesons.confPDGCodeD0, ConfTrack.confPDGCodeTrack);

    softPionRemoval.init(&qaRegistry);
    pairCleaner.init(&qaRegistry);
    if (confIsCPR.value) {
      pairCloseRejection.init(&resultRegistry, &qaRegistry, confCPRdeltaPhiCutMin.value, confCPRdeltaPhiCutMax.value, confCPRdeltaEtaCutMin.value, confCPRdeltaEtaCutMax.value, confCPRChosenRadii.value, confCPRPlotPerRadii.value);
    }

    vPIDTrack = ConfTrack.confPIDTrack.value;
    kNsigma = ConfBothTracks.confTrkPIDnSigmaMax.value;
  }

  template <typename CollisionType>
  void fillCollision(CollisionType col)
  {
    mixQaRegistry.fill(HIST("MixingQA/hSECollisionBins"), colBinning.getBin({col.posZ(), col.multNtr()}));
    eventHisto.fillQA(col);
  }

  void processD0MLOptBg(o2::aod::FdCollision const& col, FemtoFullParticles const&)
  {
    auto groupD0D0barCands = partsAllDmesons->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);

    // loop over selected D0/D0bar candidates
    for (auto const& charmCand : groupD0D0barCands) {
      // D0 candidates
      if (charmCand.mLambda() > 0.0f && charmCand.mAntiLambda() < 0.0f) {
        if (charmCand.tempFitVar() < ConfMlOpt.confClass1BgProbStart)
          registry.fill(HIST("D0D0bar_MLSel/hMassVsPt1"), charmCand.mLambda(), charmCand.pt());
        if (charmCand.tempFitVar() < ConfMlOpt.confClass1BgProbStart + ConfMlOpt.confClass1BgProbStep)
          registry.fill(HIST("D0D0bar_MLSel/hMassVsPt2"), charmCand.mLambda(), charmCand.pt());
        if (charmCand.tempFitVar() < ConfMlOpt.confClass1BgProbStart + 2.0 * ConfMlOpt.confClass1BgProbStep)
          registry.fill(HIST("D0D0bar_MLSel/hMassVsPt3"), charmCand.mLambda(), charmCand.pt());
        if (charmCand.tempFitVar() < ConfMlOpt.confClass1BgProbStart + 3.0 * ConfMlOpt.confClass1BgProbStep)
          registry.fill(HIST("D0D0bar_MLSel/hMassVsPt4"), charmCand.mLambda(), charmCand.pt());
        if (charmCand.tempFitVar() < ConfMlOpt.confClass1BgProbStart + 4.0 * ConfMlOpt.confClass1BgProbStep)
          registry.fill(HIST("D0D0bar_MLSel/hMassVsPt5"), charmCand.mLambda(), charmCand.pt());
        if (charmCand.tempFitVar() < ConfMlOpt.confClass1BgProbStart + 5.0 * ConfMlOpt.confClass1BgProbStep)
          registry.fill(HIST("D0D0bar_MLSel/hMassVsPt6"), charmCand.mLambda(), charmCand.pt());
        if (charmCand.tempFitVar() < ConfMlOpt.confClass1BgProbStart + 6.0 * ConfMlOpt.confClass1BgProbStep)
          registry.fill(HIST("D0D0bar_MLSel/hMassVsPt7"), charmCand.mLambda(), charmCand.pt());
        if (charmCand.tempFitVar() < ConfMlOpt.confClass1BgProbStart + 7.0 * ConfMlOpt.confClass1BgProbStep)
          registry.fill(HIST("D0D0bar_MLSel/hMassVsPt8"), charmCand.mLambda(), charmCand.pt());
        if (charmCand.tempFitVar() < ConfMlOpt.confClass1BgProbStart + 8.0 * ConfMlOpt.confClass1BgProbStep)
          registry.fill(HIST("D0D0bar_MLSel/hMassVsPt9"), charmCand.mLambda(), charmCand.pt());
        if (charmCand.tempFitVar() < ConfMlOpt.confClass1BgProbStart + 9.0 * ConfMlOpt.confClass1BgProbStep)
          registry.fill(HIST("D0D0bar_MLSel/hMassVsPt10"), charmCand.mLambda(), charmCand.pt());
      }
      // DObar candidates
      if (charmCand.mLambda() < 0.0f && charmCand.mAntiLambda() > 0.0f) {
        if (charmCand.tempFitVar() < ConfMlOpt.confClass1BgProbStart)
          registry.fill(HIST("D0D0bar_MLSel/hMassVsPt1"), charmCand.mAntiLambda(), charmCand.pt());
        if (charmCand.tempFitVar() < ConfMlOpt.confClass1BgProbStart + ConfMlOpt.confClass1BgProbStep)
          registry.fill(HIST("D0D0bar_MLSel/hMassVsPt2"), charmCand.mAntiLambda(), charmCand.pt());
        if (charmCand.tempFitVar() < ConfMlOpt.confClass1BgProbStart + 2.0 * ConfMlOpt.confClass1BgProbStep)
          registry.fill(HIST("D0D0bar_MLSel/hMassVsPt3"), charmCand.mAntiLambda(), charmCand.pt());
        if (charmCand.tempFitVar() < ConfMlOpt.confClass1BgProbStart + 3.0 * ConfMlOpt.confClass1BgProbStep)
          registry.fill(HIST("D0D0bar_MLSel/hMassVsPt4"), charmCand.mAntiLambda(), charmCand.pt());
        if (charmCand.tempFitVar() < ConfMlOpt.confClass1BgProbStart + 4.0 * ConfMlOpt.confClass1BgProbStep)
          registry.fill(HIST("D0D0bar_MLSel/hMassVsPt5"), charmCand.mAntiLambda(), charmCand.pt());
        if (charmCand.tempFitVar() < ConfMlOpt.confClass1BgProbStart + 5.0 * ConfMlOpt.confClass1BgProbStep)
          registry.fill(HIST("D0D0bar_MLSel/hMassVsPt6"), charmCand.mAntiLambda(), charmCand.pt());
        if (charmCand.tempFitVar() < ConfMlOpt.confClass1BgProbStart + 6.0 * ConfMlOpt.confClass1BgProbStep)
          registry.fill(HIST("D0D0bar_MLSel/hMassVsPt7"), charmCand.mAntiLambda(), charmCand.pt());
        if (charmCand.tempFitVar() < ConfMlOpt.confClass1BgProbStart + 7.0 * ConfMlOpt.confClass1BgProbStep)
          registry.fill(HIST("D0D0bar_MLSel/hMassVsPt8"), charmCand.mAntiLambda(), charmCand.pt());
        if (charmCand.tempFitVar() < ConfMlOpt.confClass1BgProbStart + 8.0 * ConfMlOpt.confClass1BgProbStep)
          registry.fill(HIST("D0D0bar_MLSel/hMassVsPt9"), charmCand.mAntiLambda(), charmCand.pt());
        if (charmCand.tempFitVar() < ConfMlOpt.confClass1BgProbStart + 9.0 * ConfMlOpt.confClass1BgProbStep)
          registry.fill(HIST("D0D0bar_MLSel/hMassVsPt10"), charmCand.mAntiLambda(), charmCand.pt());
      }
    }
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackD0, processD0MLOptBg, "Enable filling QA plots for ML Bg D0/D0bar selection optimization", false);

  void processD0MLOptBgAndPrompt(o2::aod::FdCollision const& col, FemtoFullParticles const&)
  {
    auto groupD0D0barCands = partsAllDmesons->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);

    // loop over selected D0/D0bar candidates
    for (auto const& charmCand : groupD0D0barCands) {
      // D0 candidates
      if (charmCand.decayVtxY() > ConfMlOpt.confClass2PromptProbStart) {
        if (charmCand.mLambda() > 0.0f && charmCand.mAntiLambda() < 0.0f) {
          if (charmCand.tempFitVar() < ConfMlOpt.confClass1BgProbStart)
            registry.fill(HIST("D0D0bar_MLSel/hMassVsPt1"), charmCand.mLambda(), charmCand.pt());
          if (charmCand.tempFitVar() < ConfMlOpt.confClass1BgProbStart + ConfMlOpt.confClass1BgProbStep)
            registry.fill(HIST("D0D0bar_MLSel/hMassVsPt2"), charmCand.mLambda(), charmCand.pt());
          if (charmCand.tempFitVar() < ConfMlOpt.confClass1BgProbStart + 2.0 * ConfMlOpt.confClass1BgProbStep)
            registry.fill(HIST("D0D0bar_MLSel/hMassVsPt3"), charmCand.mLambda(), charmCand.pt());
          if (charmCand.tempFitVar() < ConfMlOpt.confClass1BgProbStart + 3.0 * ConfMlOpt.confClass1BgProbStep)
            registry.fill(HIST("D0D0bar_MLSel/hMassVsPt4"), charmCand.mLambda(), charmCand.pt());
          if (charmCand.tempFitVar() < ConfMlOpt.confClass1BgProbStart + 4.0 * ConfMlOpt.confClass1BgProbStep)
            registry.fill(HIST("D0D0bar_MLSel/hMassVsPt5"), charmCand.mLambda(), charmCand.pt());
          if (charmCand.tempFitVar() < ConfMlOpt.confClass1BgProbStart + 5.0 * ConfMlOpt.confClass1BgProbStep)
            registry.fill(HIST("D0D0bar_MLSel/hMassVsPt6"), charmCand.mLambda(), charmCand.pt());
          if (charmCand.tempFitVar() < ConfMlOpt.confClass1BgProbStart + 6.0 * ConfMlOpt.confClass1BgProbStep)
            registry.fill(HIST("D0D0bar_MLSel/hMassVsPt7"), charmCand.mLambda(), charmCand.pt());
          if (charmCand.tempFitVar() < ConfMlOpt.confClass1BgProbStart + 7.0 * ConfMlOpt.confClass1BgProbStep)
            registry.fill(HIST("D0D0bar_MLSel/hMassVsPt8"), charmCand.mLambda(), charmCand.pt());
          if (charmCand.tempFitVar() < ConfMlOpt.confClass1BgProbStart + 8.0 * ConfMlOpt.confClass1BgProbStep)
            registry.fill(HIST("D0D0bar_MLSel/hMassVsPt9"), charmCand.mLambda(), charmCand.pt());
          if (charmCand.tempFitVar() < ConfMlOpt.confClass1BgProbStart + 9.0 * ConfMlOpt.confClass1BgProbStep)
            registry.fill(HIST("D0D0bar_MLSel/hMassVsPt10"), charmCand.mLambda(), charmCand.pt());
        }
        // DObar candidates
        if (charmCand.mLambda() < 0.0f && charmCand.mAntiLambda() > 0.0f) {
          if (charmCand.tempFitVar() < ConfMlOpt.confClass1BgProbStart)
            registry.fill(HIST("D0D0bar_MLSel/hMassVsPt1"), charmCand.mAntiLambda(), charmCand.pt());
          if (charmCand.tempFitVar() < ConfMlOpt.confClass1BgProbStart + ConfMlOpt.confClass1BgProbStep)
            registry.fill(HIST("D0D0bar_MLSel/hMassVsPt2"), charmCand.mAntiLambda(), charmCand.pt());
          if (charmCand.tempFitVar() < ConfMlOpt.confClass1BgProbStart + 2.0 * ConfMlOpt.confClass1BgProbStep)
            registry.fill(HIST("D0D0bar_MLSel/hMassVsPt3"), charmCand.mAntiLambda(), charmCand.pt());
          if (charmCand.tempFitVar() < ConfMlOpt.confClass1BgProbStart + 3.0 * ConfMlOpt.confClass1BgProbStep)
            registry.fill(HIST("D0D0bar_MLSel/hMassVsPt4"), charmCand.mAntiLambda(), charmCand.pt());
          if (charmCand.tempFitVar() < ConfMlOpt.confClass1BgProbStart + 4.0 * ConfMlOpt.confClass1BgProbStep)
            registry.fill(HIST("D0D0bar_MLSel/hMassVsPt5"), charmCand.mAntiLambda(), charmCand.pt());
          if (charmCand.tempFitVar() < ConfMlOpt.confClass1BgProbStart + 5.0 * ConfMlOpt.confClass1BgProbStep)
            registry.fill(HIST("D0D0bar_MLSel/hMassVsPt6"), charmCand.mAntiLambda(), charmCand.pt());
          if (charmCand.tempFitVar() < ConfMlOpt.confClass1BgProbStart + 6.0 * ConfMlOpt.confClass1BgProbStep)
            registry.fill(HIST("D0D0bar_MLSel/hMassVsPt7"), charmCand.mAntiLambda(), charmCand.pt());
          if (charmCand.tempFitVar() < ConfMlOpt.confClass1BgProbStart + 7.0 * ConfMlOpt.confClass1BgProbStep)
            registry.fill(HIST("D0D0bar_MLSel/hMassVsPt8"), charmCand.mAntiLambda(), charmCand.pt());
          if (charmCand.tempFitVar() < ConfMlOpt.confClass1BgProbStart + 8.0 * ConfMlOpt.confClass1BgProbStep)
            registry.fill(HIST("D0D0bar_MLSel/hMassVsPt9"), charmCand.mAntiLambda(), charmCand.pt());
          if (charmCand.tempFitVar() < ConfMlOpt.confClass1BgProbStart + 9.0 * ConfMlOpt.confClass1BgProbStep)
            registry.fill(HIST("D0D0bar_MLSel/hMassVsPt10"), charmCand.mAntiLambda(), charmCand.pt());
        }
      }
    }
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackD0, processD0MLOptBgAndPrompt, "Enable filling QA plots for ML Bg and Prompt D0/D0bar selection optimization", false);

  void processQAD0D0barSel(o2::aod::FdCollision const& col, FemtoFullParticles const&)
  {
    auto groupPartsD0s = partsD0s->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto groupPartsD0bars = partsD0bars->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto groupPartsD0sFromSB = partsD0sFromSB->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);

    // loop over selected D0 candidates
    for (auto const& d0cand : groupPartsD0s) {

      qaRegistry.fill(HIST("QA_D0D0barSelection/hInvMassD0"), d0cand.mLambda());
      qaRegistry.fill(HIST("QA_D0D0barSelection/hPtD0"), d0cand.pt());
    }
    // loop over selected D0bar candidates
    for (auto const& d0barcand : groupPartsD0bars) {

      qaRegistry.fill(HIST("QA_D0D0barSelection/hInvMassD0bar"), d0barcand.mAntiLambda());
      qaRegistry.fill(HIST("QA_D0D0barSelection/hPtD0bar"), d0barcand.pt());
    }

    // loop over selected D0 candidates from SB regions
    for (auto const& d0SB : groupPartsD0sFromSB) {

      qaRegistry.fill(HIST("QA_D0D0barSelection_SB/hInvMassD0"), d0SB.mLambda());
      qaRegistry.fill(HIST("QA_D0D0barSelection_SB/hPtD0"), d0SB.pt());
    }
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackD0, processQAD0D0barSel, "Enable filling QA plots for selected D0/D0bar cand.", true);

  void processAllDmesons(o2::aod::FdCollision const& col, FemtoFullParticles const&)
  {
    auto groupPartsAllD0D0barCands = partsAllDmesons->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);

    // loop over D0/D0bar mesons
    for (auto const& d0d0bar : groupPartsAllD0D0barCands) {

      registry.fill(HIST("hPtD0D0bar"), d0d0bar.pt());
      registry.fill(HIST("hPhiD0D0bar"), d0d0bar.phi());
      registry.fill(HIST("hEtaD0D0bar"), d0d0bar.eta());
      // BDT score classes
      registry.fill(HIST("DebugBdt/hBdtScore1VsStatus"), d0d0bar.decayVtxX(), 1);
      registry.fill(HIST("DebugBdt/hBdtScore2VsStatus"), d0d0bar.decayVtxY(), 1);
      registry.fill(HIST("DebugBdt/hBdtScore3VsStatus"), d0d0bar.decayVtxZ(), 1);

      weight = 1.0f;
      if (doEfficiencyCorr) {
        weight = efficiencyCalculator.getWeight(ParticleNo::TWO, d0d0bar.pt());
        if (d0d0bar.mLambda() > 0.0f) {
          registry.fill(HIST("D0D0bar_doubleMassHypo/hMassVsPtEffCorr"), d0d0bar.mLambda(), d0d0bar.pt(), weight);
          registry.fill(HIST("D0D0bar_doubleMassHypo/hMassVsPtD0EffCorr"), d0d0bar.mLambda(), d0d0bar.pt(), weight);
          if (d0d0bar.mAntiLambda() < 0.0f) {
            registry.fill(HIST("D0D0bar_oneMassHypo/hMassVsPtEffCorr"), d0d0bar.mLambda(), d0d0bar.pt(), weight);
            registry.fill(HIST("D0D0bar_oneMassHypo/hMassVsPtD0EffCorr"), d0d0bar.mLambda(), d0d0bar.pt(), weight);
            registry.fill(HIST("hPtD0"), d0d0bar.pt());
            registry.fill(HIST("hPhiD0"), d0d0bar.phi());
            registry.fill(HIST("hEtaD0"), d0d0bar.eta());
          }
        }
        if (d0d0bar.mAntiLambda() > 0.0f) {
          registry.fill(HIST("D0D0bar_doubleMassHypo/hMassVsPtEffCorr"), d0d0bar.mLambda(), d0d0bar.pt(), weight);
          registry.fill(HIST("D0D0bar_doubleMassHypo/hMassVsPtD0barEffCorr"), d0d0bar.mLambda(), d0d0bar.pt(), weight);
          if (d0d0bar.mLambda() < 0.0f) {
            registry.fill(HIST("D0D0bar_oneMassHypo/hMassVsPtEffCorr"), d0d0bar.mAntiLambda(), d0d0bar.pt(), weight);
            registry.fill(HIST("D0D0bar_oneMassHypo/hMassVsPtD0barEffCorr"), d0d0bar.mAntiLambda(), d0d0bar.pt(), weight);
            registry.fill(HIST("hPtD0bar"), d0d0bar.pt());
            registry.fill(HIST("hPhiD0bar"), d0d0bar.phi());
            registry.fill(HIST("hEtaD0bar"), d0d0bar.eta());
          }
        }
      } else {
        if (d0d0bar.mLambda() > 0.0f) {
          registry.fill(HIST("D0D0bar_doubleMassHypo/hMassVsPt"), d0d0bar.mLambda(), d0d0bar.pt());
          registry.fill(HIST("D0D0bar_doubleMassHypo/hMassVsPtD0"), d0d0bar.mLambda(), d0d0bar.pt());
          if (d0d0bar.mAntiLambda() < 0.0f) {
            registry.fill(HIST("D0D0bar_oneMassHypo/hMassVsPt"), d0d0bar.mLambda(), d0d0bar.pt());
            registry.fill(HIST("D0D0bar_oneMassHypo/hMassVsPtD0"), d0d0bar.mLambda(), d0d0bar.pt());
            registry.fill(HIST("hPtD0"), d0d0bar.pt());
            registry.fill(HIST("hPhiD0"), d0d0bar.phi());
            registry.fill(HIST("hEtaD0"), d0d0bar.eta());
          }
        }
        if (d0d0bar.mAntiLambda() > 0.0f) {
          registry.fill(HIST("D0D0bar_doubleMassHypo/hMassVsPt"), d0d0bar.mLambda(), d0d0bar.pt());
          registry.fill(HIST("D0D0bar_doubleMassHypo/hMassVsPtD0bar"), d0d0bar.mLambda(), d0d0bar.pt());
          if (d0d0bar.mLambda() < 0.0f) {
            registry.fill(HIST("D0D0bar_oneMassHypo/hMassVsPt"), d0d0bar.mAntiLambda(), d0d0bar.pt());
            registry.fill(HIST("D0D0bar_oneMassHypo/hMassVsPtD0bar"), d0d0bar.mAntiLambda(), d0d0bar.pt());
            registry.fill(HIST("hPtD0bar"), d0d0bar.pt());
            registry.fill(HIST("hPhiD0bar"), d0d0bar.phi());
            registry.fill(HIST("hEtaD0bar"), d0d0bar.eta());
          }
        }
      }
    }
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackD0, processAllDmesons, "Enable processing over all D meson candidates", false);

  void processD0Children(o2::aod::FdCollision const& col, FemtoFullParticles const&)
  {
    auto groupPartsD0D0barChildren = partsDmesonsChildren->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);

    // loop over D mesons childen
    for (auto const& daughD0D0bar : groupPartsD0D0barChildren) {
      registry.fill(HIST("hPtDaughters"), daughD0D0bar.pt());
      registry.fill(HIST("hSignDaughters"), daughD0D0bar.mLambda());
      // filling QA plots for D0 mesons' positive daughters (K+)
      if (daughD0D0bar.mLambda() == 1 && (daughD0D0bar.mAntiLambda() == 1 || daughD0D0bar.mAntiLambda() == 0)) {
        qaRegistry.fill(HIST("D0_pos_daugh/pt"), daughD0D0bar.pt());
        qaRegistry.fill(HIST("D0_pos_daugh/eta"), daughD0D0bar.eta());
        qaRegistry.fill(HIST("D0_pos_daugh/phi"), daughD0D0bar.phi());
      }
      // filling QA plots for D0 mesons' negative daughters (pi-)
      if (daughD0D0bar.mLambda() == -1 && (daughD0D0bar.mAntiLambda() == 1 || daughD0D0bar.mAntiLambda() == 0)) {
        qaRegistry.fill(HIST("D0_neg_daugh/pt"), daughD0D0bar.pt());
        qaRegistry.fill(HIST("D0_neg_daugh/eta"), daughD0D0bar.eta());
        qaRegistry.fill(HIST("D0_neg_daugh/phi"), daughD0D0bar.phi());
      }
      // filling QA plots for D0bar mesons' positive daughters (pi+)
      if (daughD0D0bar.mLambda() == 1 && (daughD0D0bar.mAntiLambda() == -1 || daughD0D0bar.mAntiLambda() == 0)) {
        qaRegistry.fill(HIST("D0bar_pos_daugh/pt"), daughD0D0bar.pt());
        qaRegistry.fill(HIST("D0bar_pos_daugh/eta"), daughD0D0bar.eta());
        qaRegistry.fill(HIST("D0bar_pos_daugh/phi"), daughD0D0bar.phi());
      }
      // filling QA plots for D0bar mesons' negative daughters (K-)
      if (daughD0D0bar.mLambda() == -1 && (daughD0D0bar.mAntiLambda() == -1 || daughD0D0bar.mAntiLambda() == 0)) {
        qaRegistry.fill(HIST("D0bar_neg_daugh/pt"), daughD0D0bar.pt());
        qaRegistry.fill(HIST("D0bar_neg_daugh/eta"), daughD0D0bar.eta());
        qaRegistry.fill(HIST("D0bar_neg_daugh/phi"), daughD0D0bar.phi());
      }
    }
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackD0, processD0Children, "Enable processing over daughters of D0 meson candidates", false);

  /// This function processes the same event and takes care of all the histogramming
  /// \todo the trivial loops over the tracks should be factored out since they will be common to all combinations of T-T, T-V0, V0-V0, ...
  /// @tparam PartitionType
  /// @tparam PartType
  /// @tparam isMC: enables Monte Carlo truth specific histograms
  /// @param groupPartsTrack partition for the first particle passed by the process function
  /// @param groupPartsD0 partition for the second particle passed by the process function
  /// @param parts femtoUniverseParticles table (in case of Monte Carlo joined with FemtoUniverseMCLabels)
  /// @param magFieldTesla magnetic field of the collision
  /// @param multCol multiplicity of the collision
  template <bool isMC, typename PartitionType, typename PartType>
  void doSameEvent(PartitionType groupPartsTrack, PartitionType groupPartsD0, PartType parts, float magFieldTesla, int multCol)
  {

    /// Histogramming same event
    for (auto const& d0candidate : groupPartsD0) {
      trackHistoPartD0D0bar.fillQA<isMC, false>(d0candidate);
    }

    float tpcNSigmaPr, tofNSigmaPr, tpcNSigmaPi, tofNSigmaPi, tpcNSigmaKa, tofNSigmaKa;

    if (!ConfTrack.confIsSame) {
      for (auto const& track : groupPartsTrack) {
        if (ConfTrack.confIsTrackIdentified) {
          if (!isParticleNSigma(track.p(), trackCuts.getNsigmaTPC(track, o2::track::PID::Proton), trackCuts.getNsigmaTOF(track, o2::track::PID::Proton), trackCuts.getNsigmaTPC(track, o2::track::PID::Pion), trackCuts.getNsigmaTOF(track, o2::track::PID::Pion), trackCuts.getNsigmaTPC(track, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(track, o2::track::PID::Kaon))) {
            continue;
          }
        }
        trackHistoPartTrack.fillQA<isMC, false>(track);

        tpcNSigmaPi = trackCuts.getNsigmaTPC(track, o2::track::PID::Pion);
        tofNSigmaPi = trackCuts.getNsigmaTOF(track, o2::track::PID::Pion);
        tpcNSigmaKa = trackCuts.getNsigmaTPC(track, o2::track::PID::Kaon);
        tofNSigmaKa = trackCuts.getNsigmaTOF(track, o2::track::PID::Kaon);
        tpcNSigmaPr = trackCuts.getNsigmaTPC(track, o2::track::PID::Proton);
        tofNSigmaPr = trackCuts.getNsigmaTOF(track, o2::track::PID::Proton);

        qaRegistry.fill(HIST("Hadron/nSigmaTPCPi"), track.p(), tpcNSigmaPi);
        qaRegistry.fill(HIST("Hadron/nSigmaTOFPi"), track.p(), tofNSigmaPi);
        qaRegistry.fill(HIST("Hadron/nSigmaTPCKa"), track.p(), tpcNSigmaKa);
        qaRegistry.fill(HIST("Hadron/nSigmaTOFKa"), track.p(), tofNSigmaKa);
        qaRegistry.fill(HIST("Hadron/nSigmaTPCPr"), track.p(), tpcNSigmaPr);
        qaRegistry.fill(HIST("Hadron/nSigmaTOFPr"), track.p(), tofNSigmaPr);
      }
    }
    /// Now build the combinations
    for (auto const& [track, d0candidate] : combinations(CombinationsFullIndexPolicy(groupPartsTrack, groupPartsD0))) {
      if (ConfTrack.confIsTrackIdentified) {
        if (!isParticleNSigma(track.p(), trackCuts.getNsigmaTPC(track, o2::track::PID::Proton), trackCuts.getNsigmaTOF(track, o2::track::PID::Proton), trackCuts.getNsigmaTPC(track, o2::track::PID::Pion), trackCuts.getNsigmaTOF(track, o2::track::PID::Pion), trackCuts.getNsigmaTPC(track, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(track, o2::track::PID::Kaon))) {
          continue;
        }
      }
      // Soft Pion Removal
      if (confRemoveSoftPions) {
        if (softPionRemoval.isSoftPion(track, d0candidate, parts, confSoftPionD0Flag, confSoftPionD0barFlag, sigmaSoftPiInvMass)) {
          continue;
        }
      }
      // // Close Pair Rejection
      if (confIsCPR.value) {
        if (pairCloseRejection.isClosePair(track, d0candidate, parts, magFieldTesla, femto_universe_container::EventType::same)) {
          continue;
        }
      }
      // Track Cleaning
      if (!pairCleaner.isCleanPair(track, d0candidate, parts)) {
        continue;
      }
      // Efficiency
      weight = 1.0f;
      if (doEfficiencyCorr) {
        weight = efficiencyCalculator.getWeight(ParticleNo::ONE, track.pt()) * efficiencyCalculator.getWeight(ParticleNo::TWO, d0candidate.pt());
      }
      sameEventAngularCont.setPair<isMC>(track, d0candidate, multCol, ConfBothTracks.confUse3D, weight);
    }
  }

  /// process function for to call doSameEvent with Data
  /// call this process function if you need D0/D0bar candidates which pass only one mass hypothesis
  /// \param col subscribe to the collision table (Data)
  /// \param parts subscribe to the femtoUniverseParticleTable
  void processSameEvent(o2::aod::FdCollision const& col,
                        FemtoFullParticles const& parts)
  {
    fillCollision(col);

    auto theGroupPartsTrack = partsTrack->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto theGroupPartsD0s = partsD0s->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto theGroupPartsD0bars = partsD0bars->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);

    switch (confChooseD0trackCorr) {
      case 0:
        doSameEvent<false>(theGroupPartsTrack, theGroupPartsD0s, parts, col.magField(), col.multNtr());
        break;
      case 1:
        doSameEvent<false>(theGroupPartsTrack, theGroupPartsD0bars, parts, col.magField(), col.multNtr());
        break;
      default:
        break;
    }
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackD0, processSameEvent, "Enable processing same event", true);

  /// process function for to call doSameEvent with Data
  /// call this process function to include candidates which pass as well the selection for both D0 and D0bar candidates
  /// \param col subscribe to the collision table (Data)
  /// \param parts subscribe to the femtoUniverseParticleTable
  void processSameEventWithDoubleHypo(o2::aod::FdCollision const& col,
                                      FemtoFullParticles const& parts)
  {
    fillCollision(col);

    auto theGroupPartsTrack = partsTrack->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto theGroupPartsD0s = partsAllD0s->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto theGroupPartsD0bars = partsAllD0bars->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);

    switch (confChooseD0trackCorr) {
      case 0:
        doSameEvent<false>(theGroupPartsTrack, theGroupPartsD0s, parts, col.magField(), col.multNtr());
        break;
      case 1:
        doSameEvent<false>(theGroupPartsTrack, theGroupPartsD0bars, parts, col.magField(), col.multNtr());
        break;
      default:
        break;
    }
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackD0, processSameEventWithDoubleHypo, "Enable processing same event", false);

  /// process function for to call doSameEvent with Data
  /// call this process to obtain the function for D0/D0bar candidates from side-band regions
  /// \param col subscribe to the collision table (Data)
  /// \param parts subscribe to the femtoUniverseParticleTable
  void processSameEventSB(o2::aod::FdCollision const& col, FemtoFullParticles const& parts)
  {
    fillCollision(col);

    auto groupPartsTrack = partsTrack->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto groupPartsD0sFromSB = partsD0sFromSB->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto groupPartsD0barsFromSB = partsD0barsFromSB->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);

    switch (confChooseD0trackCorr) {
      case 0:
        doSameEvent<false>(groupPartsTrack, groupPartsD0sFromSB, parts, col.magField(), col.multNtr());
        break;
      case 1:
        doSameEvent<false>(groupPartsTrack, groupPartsD0barsFromSB, parts, col.magField(), col.multNtr());
        break;
      default:
        break;
    }
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackD0, processSameEventSB, "Enable processing same event", false);

  /// process function for to call doSameEvent with Monte Carlo
  /// \param col subscribe to the collision table (Monte Carlo Reconstructed)
  /// \param parts subscribe to joined table FemtoUniverseParticles and FemtoUniverseMCLables to access Monte Carlo truth
  /// \param FemtoUniverseMCParticles subscribe to the Monte Carlo truth table
  void processSameEventMCRecoSignal(aod::FdCollision const& col, FemtoMcRecoParticles const& parts, aod::FdMCParticles const&)
  {
    fillCollision(col);

    auto theGroupPartsD0 = partsPromptD0MCReco->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto theGroupPartsD0bar = partsPromptD0barMCReco->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto theGroupPartsTrack = partsTrackMCReco->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);

    switch (confChooseD0trackCorr) {
      case 0:
        doSameEvent<false>(theGroupPartsTrack, theGroupPartsD0, parts, col.magField(), col.multNtr());
        break;
      case 1:
        doSameEvent<false>(theGroupPartsTrack, theGroupPartsD0bar, parts, col.magField(), col.multNtr());
        break;
      default:
        break;
    }
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackD0, processSameEventMCRecoSignal, "MC Reco - enable processing same event for D0/D0bar signal cand.", false);

  /// process function for to call doSameEvent with Monte Carlo
  /// \param col subscribe to the collision table (Monte Carlo Reconstructed)
  /// \param parts subscribe to joined table FemtoUniverseParticles and FemtoUniverseMCLables to access Monte Carlo truth
  /// \param FemtoUniverseMCParticles subscribe to the Monte Carlo truth table
  void processSameEventMCRecoSB(aod::FdCollision const& col, FemtoMcRecoParticles const& parts, aod::FdMCParticles const&)
  {
    fillCollision(col);

    auto theGroupPartsD0 = partsPromptD0MCRecoSB->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto theGroupPartsD0bar = partsPromptD0barMCRecoSB->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto theGroupPartsTrack = partsTrackMCReco->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);

    switch (confChooseD0trackCorr) {
      case 0:
        doSameEvent<false>(theGroupPartsTrack, theGroupPartsD0, parts, col.magField(), col.multNtr());
        break;
      case 1:
        doSameEvent<false>(theGroupPartsTrack, theGroupPartsD0bar, parts, col.magField(), col.multNtr());
        break;
      default:
        break;
    }
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackD0, processSameEventMCRecoSB, "MC Reco - enable processing same event for D0/D0bar side-band cand.", false);

  /// process function for to call doSameEvent with Monte Carlo
  /// \param col subscribe to the collision table
  /// \param parts subscribe to the joined table FemtoUniverseParticles and FemtoUniverseMCLables to access Monte Carlo truth
  /// \param FdMCParticles subscribe to the Monte Carlo truth table
  void processSameEventMcTruth(aod::FdCollision const& col, FemtoMcRecoParticles const& parts, aod::FdMCParticles const&)
  {
    fillCollision(col);
    /// !!! needs to be finished -> add possibility to change between D0 and D0bar
    auto theGroupTruthD0s = partsD0MCTruth->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto theGroupTruthD0bars = partsD0barMCTruth->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto theGroupTruthTracks = partsTrackMCTruth->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);

    switch (confChooseD0trackCorr) {
      case 0:
        doSameEvent<true>(theGroupTruthTracks, theGroupTruthD0s, parts, col.magField(), col.multNtr());
        break;
      case 1:
        doSameEvent<true>(theGroupTruthTracks, theGroupTruthD0bars, parts, col.magField(), col.multNtr());
        break;
      default:
        break;
    }
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackD0, processSameEventMcTruth, "MC Truth - enable processing same event", false);

  /// This function processes the mixed event
  /// \todo the trivial loops over the collisions and tracks should be factored out since they will be common to all combinations of T-T, T-V0, V0-V0, ...
  /// \tparam PartitionType
  /// \tparam PartType
  /// \tparam isMC: enables Monte Carlo truth specific histograms
  /// \param groupPartsTrack partition for the identified passed by the process function
  /// \param groupPartsD0 partition for D0 meson passed by the process function
  /// \param parts femtoUniverseParticles table (in case of Monte Carlo joined with FemtoUniverseMCLabels)
  /// \param magFieldTesla magnetic field of the collision
  /// \param multCol multiplicity of the collision
  template <bool isMC, typename PartitionType, typename PartType>
  void doMixedEvent(PartitionType groupPartsTrack, PartitionType groupPartsD0, PartType parts, float magFieldTesla, int multCol)
  {

    for (auto const& [track, d0candidate] : combinations(CombinationsFullIndexPolicy(groupPartsTrack, groupPartsD0))) {
      if (ConfTrack.confIsTrackIdentified) {
        if (!isParticleNSigma(track.p(), trackCuts.getNsigmaTPC(track, o2::track::PID::Proton), trackCuts.getNsigmaTOF(track, o2::track::PID::Proton), trackCuts.getNsigmaTPC(track, o2::track::PID::Pion), trackCuts.getNsigmaTOF(track, o2::track::PID::Pion), trackCuts.getNsigmaTPC(track, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(track, o2::track::PID::Kaon))) {
          continue;
        }
      }
      // // Soft Pion Removal
      if (confRemoveSoftPions) {
        if (softPionRemoval.isSoftPion(track, d0candidate, parts, confSoftPionD0Flag, confSoftPionD0barFlag, sigmaSoftPiInvMass)) {
          continue;
        }
      }
      // // Close Pair Rejection
      if (confIsCPR.value) {
        if (pairCloseRejection.isClosePair(track, d0candidate, parts, magFieldTesla, femto_universe_container::EventType::mixed)) {
          continue;
        }
      }
      // // Track Cleaning
      if (!pairCleaner.isCleanPair(track, d0candidate, parts)) {
        continue;
      }
      // Efficiency
      weight = 1.0f;
      if (doEfficiencyCorr) {
        weight = efficiencyCalculator.getWeight(ParticleNo::ONE, track.pt()) * efficiencyCalculator.getWeight(ParticleNo::TWO, d0candidate.pt());
      }

      mixedEventAngularCont.setPair<isMC>(track, d0candidate, multCol, ConfBothTracks.confUse3D, weight);
    }
  }

  /// process function for to call doMixedEvent with Data
  /// call this process function if you need D0/D0bar candidates which pass only one mass hypothesis
  /// @param cols subscribe to the collisions table (Data)
  /// @param parts subscribe to the femtoUniverseParticleTable
  void processMixedEvent(o2::aod::FdCollisions const& cols,
                         FemtoFullParticles const& parts)
  {
    for (auto const& [collision1, collision2] : soa::selfCombinations(colBinning, confNEventsMix, -1, cols, cols)) {

      const int multiplicityCol = collision1.multNtr();
      mixQaRegistry.fill(HIST("MixingQA/hMECollisionBins"), colBinning.getBin({collision1.posZ(), multiplicityCol}));

      auto groupPartsTrack = partsTrack->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision2.globalIndex(), cache);
      auto theGroupPartsD0s = partsD0s->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);
      auto theGroupPartsD0bars = partsD0bars->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);

      const auto& magFieldTesla1 = collision1.magField();
      const auto& magFieldTesla2 = collision2.magField();

      if (magFieldTesla1 != magFieldTesla2) {
        continue;
      }
      /// \todo before mixing we should check whether both collisions contain a pair of particles!
      // if (partsD0.size() == 0 || kNPart2Evt1 == 0 || kNPart1Evt2 == 0 || partsTrack.size() == 0 ) continue;

      switch (confChooseD0trackCorr) {
        case 0:
          doMixedEvent<false>(groupPartsTrack, theGroupPartsD0s, parts, magFieldTesla1, multiplicityCol);
          break;
        case 1:
          doMixedEvent<false>(groupPartsTrack, theGroupPartsD0bars, parts, magFieldTesla1, multiplicityCol);
          break;
        default:
          break;
      }
    }
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackD0, processMixedEvent, "Enable processing mixed events", true);

  /// process function for to call doMixedEvent with Data
  /// call this process function to include candidates which pass as well the selection for both D0 and D0bar candidates
  /// @param cols subscribe to the collisions table (Data)
  /// @param parts subscribe to the femtoUniverseParticleTable
  void processMixedEventWithDoubleHypo(o2::aod::FdCollisions const& cols,
                                       FemtoFullParticles const& parts)
  {
    for (auto const& [collision1, collision2] : soa::selfCombinations(colBinning, confNEventsMix, -1, cols, cols)) {

      const int multiplicityCol = collision1.multNtr();
      mixQaRegistry.fill(HIST("MixingQA/hMECollisionBins"), colBinning.getBin({collision1.posZ(), multiplicityCol}));

      auto groupPartsTrack = partsTrack->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision2.globalIndex(), cache);
      auto theGroupPartsD0s = partsAllD0s->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);
      auto theGroupPartsD0bars = partsAllD0bars->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);

      const auto& magFieldTesla1 = collision1.magField();
      const auto& magFieldTesla2 = collision2.magField();

      if (magFieldTesla1 != magFieldTesla2) {
        continue;
      }
      /// \todo before mixing we should check whether both collisions contain a pair of particles!
      // if (partsD0.size() == 0 || kNPart2Evt1 == 0 || kNPart1Evt2 == 0 || partsTrack.size() == 0 ) continue;

      switch (confChooseD0trackCorr) {
        case 0:
          doMixedEvent<false>(groupPartsTrack, theGroupPartsD0s, parts, magFieldTesla1, multiplicityCol);
          break;
        case 1:
          doMixedEvent<false>(groupPartsTrack, theGroupPartsD0bars, parts, magFieldTesla1, multiplicityCol);
          break;
        default:
          break;
      }
    }
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackD0, processMixedEventWithDoubleHypo, "Enable processing mixed events", false);

  /// process function for to call doMixedEvent with Data
  /// call this process to obtain the function for D0/D0bar candidates from side-band regions
  /// @param cols subscribe to the collisions table (Data)
  /// @param parts subscribe to the femtoUniverseParticleTable
  void processMixedEventSB(o2::aod::FdCollisions const& cols, FemtoFullParticles const& parts)
  {
    for (auto const& [collision1, collision2] : soa::selfCombinations(colBinning, confNEventsMix, -1, cols, cols)) {

      const int multiplicityCol = collision1.multNtr();
      mixQaRegistry.fill(HIST("MixingQA/hMECollisionBins"), colBinning.getBin({collision1.posZ(), multiplicityCol}));

      auto groupPartsTrack = partsTrack->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision2.globalIndex(), cache);
      auto groupPartsD0sFromSB = partsD0sFromSB->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);
      auto groupPartsD0barsFromSB = partsD0barsFromSB->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);

      const auto& magFieldTesla1 = collision1.magField();
      const auto& magFieldTesla2 = collision2.magField();

      if (magFieldTesla1 != magFieldTesla2) {
        continue;
      }
      /// \todo before mixing we should check whether both collisions contain a pair of particles!
      // if (partsD0.size() == 0 || kNPart2Evt1 == 0 || kNPart1Evt2 == 0 || partsTrack.size() == 0 ) continue;
      switch (confChooseD0trackCorr) {
        case 0:
          doMixedEvent<false>(groupPartsTrack, groupPartsD0sFromSB, parts, magFieldTesla1, multiplicityCol);
          break;
        case 1:
          doMixedEvent<false>(groupPartsTrack, groupPartsD0barsFromSB, parts, magFieldTesla1, multiplicityCol);
          break;
        default:
          break;
      }
    }
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackD0, processMixedEventSB, "Enable processing mixed events for data", false);

  /// brief process function for to call doMixedEvent with Monte Carlo
  /// @param cols subscribe to the collisions table (Monte Carlo Reconstructed)
  /// @param parts subscribe to joined table FemtoUniverseParticles and FemtoUniverseMCLables to access Monte Carlo truth
  /// @param FemtoUniverseMCParticles subscribe to the Monte Carlo truth table
  void processMixedEventMCRecoSignal(aod::FdCollisions const& cols, FemtoMcRecoParticles const& parts, aod::FdMCParticles const&)
  {
    for (auto const& [collision1, collision2] : soa::selfCombinations(colBinning, confNEventsMix, -1, cols, cols)) {

      const int multiplicityCol = collision1.multNtr();
      mixQaRegistry.fill(HIST("MixingQA/hMECollisionBins"), colBinning.getBin({collision1.posZ(), multiplicityCol}));

      auto groupPartsTrack = partsTrackMCReco->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision2.globalIndex(), cache);
      auto groupPartsD0 = partsPromptD0MCReco->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);
      auto groupPartsD0bar = partsPromptD0barMCReco->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);

      const auto& magFieldTesla1 = collision1.magField();
      const auto& magFieldTesla2 = collision2.magField();

      if (magFieldTesla1 != magFieldTesla2) {
        continue;
      }

      switch (confChooseD0trackCorr) {
        case 0:
          doMixedEvent<false>(groupPartsTrack, groupPartsD0, parts, magFieldTesla1, multiplicityCol);
          break;
        case 1:
          doMixedEvent<false>(groupPartsTrack, groupPartsD0bar, parts, magFieldTesla1, multiplicityCol);
          break;
        default:
          break;
      }
    }
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackD0, processMixedEventMCRecoSignal, "MC Reco - enable processing mixed events for D0/D0bar signal cand.", false);

  /// brief process function for to call doMixedEvent with Monte Carlo Reco
  /// @param cols subscribe to the collisions table (Monte Carlo Reconstructed)
  /// @param parts subscribe to joined table FemtoUniverseParticles and FemtoUniverseMCLables to access Monte Carlo truth
  /// @param FemtoUniverseMCParticles subscribe to the Monte Carlo truth table
  void processMixedEventMCRecoSB(aod::FdCollisions const& cols, FemtoMcRecoParticles const& parts, aod::FdMCParticles const&)
  {
    for (auto const& [collision1, collision2] : soa::selfCombinations(colBinning, confNEventsMix, -1, cols, cols)) {

      const int multiplicityCol = collision1.multNtr();
      mixQaRegistry.fill(HIST("MixingQA/hMECollisionBins"), colBinning.getBin({collision1.posZ(), multiplicityCol}));

      auto groupPartsTrack = partsTrackMCReco->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision2.globalIndex(), cache);
      auto groupPartsD0 = partsPromptD0MCRecoSB->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);
      auto groupPartsD0bar = partsPromptD0barMCRecoSB->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);

      const auto& magFieldTesla1 = collision1.magField();
      const auto& magFieldTesla2 = collision2.magField();

      if (magFieldTesla1 != magFieldTesla2) {
        continue;
      }

      switch (confChooseD0trackCorr) {
        case 0:
          doMixedEvent<false>(groupPartsTrack, groupPartsD0, parts, magFieldTesla1, multiplicityCol);
          break;
        case 1:
          doMixedEvent<false>(groupPartsTrack, groupPartsD0bar, parts, magFieldTesla1, multiplicityCol);
          break;
        default:
          break;
      }
    }
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackD0, processMixedEventMCRecoSB, "MC Reco - enable processing mixed events for D0/D0bar side-band cand.", false);

  /// brief process function for to call doMixedEvent with Monte Carlo Reco
  /// @param cols subscribe to the collisions table (Monte Carlo Reconstructed reconstructed)
  /// @param parts subscribe to joined table FemtoUniverseParticles and FemtoUniverseMCLables to access Monte Carlo truth
  /// @param FemtoUniverseMCParticles subscribe to the Monte Carlo truth table
  void processMixedEventMCTruth(aod::FdCollisions const& cols, FemtoMcRecoParticles const& parts, aod::FdMCParticles const&)
  {
    for (auto const& [collision1, collision2] : soa::selfCombinations(colBinning, confNEventsMix, -1, cols, cols)) {

      const int multiplicityCol = collision1.multNtr();
      mixQaRegistry.fill(HIST("MixingQA/hMECollisionBins"), colBinning.getBin({collision1.posZ(), multiplicityCol}));

      auto groupPartsTrack = partsTrackMCReco->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision2.globalIndex(), cache);
      auto groupPartsD0 = partsD0MCTruth->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);
      auto groupPartsD0bar = partsD0barMCTruth->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);

      const auto& magFieldTesla1 = collision1.magField();
      const auto& magFieldTesla2 = collision2.magField();

      if (magFieldTesla1 != magFieldTesla2) {
        continue;
      }

      switch (confChooseD0trackCorr) {
        case 0:
          doMixedEvent<true>(groupPartsTrack, groupPartsD0, parts, magFieldTesla1, multiplicityCol);
          break;
        case 1:
          doMixedEvent<true>(groupPartsTrack, groupPartsD0bar, parts, magFieldTesla1, multiplicityCol);
          break;
        default:
          break;
      }
    }
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackD0, processMixedEventMCTruth, "MC Truth - enable processing mixed events for gen. D0/D0bar", false);

  void processMcReco(FemtoMcRecoParticles const& recoParts, aod::FdMCParticles const& mcParts)
  {
    for (auto const& part : recoParts) {
      auto mcPartId = part.fdMCParticleId();
      if (mcPartId == -1) {
        continue; // no MC particle
      }
      const auto& mcpart = mcParts.iteratorAt(mcPartId);
      weight = 1.0f;

      // filling the histograms for identified hadrons
      if (part.partType() == aod::femtouniverseparticle::ParticleType::kTrack) {
        if (part.sign() > 0) {
          if ((mcpart.pdgMCTruth() == PDG_t::kProton) && (part.pt() > ConfTrack.confTrackLowPtCut) && (part.pt() < ConfTrack.confTrackHighPtCut) && isParticleNSigma(part.p(), trackCuts.getNsigmaTPC(part, o2::track::PID::Proton), trackCuts.getNsigmaTOF(part, o2::track::PID::Proton), trackCuts.getNsigmaTPC(part, o2::track::PID::Pion), trackCuts.getNsigmaTOF(part, o2::track::PID::Pion), trackCuts.getNsigmaTPC(part, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(part, o2::track::PID::Kaon))) {
            mcRecoRegistry.fill(HIST("hMcRecPrPt"), part.pt());
            if (mcpart.partOriginMCTruth() == aod::femtouniverse_mc_particle::ParticleOriginMCTruth::kPrimary) {
              mcRecoRegistry.fill(HIST("hMcRecPrPtGenEtaGen"), mcpart.pt(), mcpart.eta());
            }
          }
          if ((mcpart.pdgMCTruth() == PDG_t::kPiPlus) && (part.pt() > ConfTrack.confTrackLowPtCut) && (part.pt() < ConfTrack.confTrackHighPtCut) && isParticleNSigma(part.p(), trackCuts.getNsigmaTPC(part, o2::track::PID::Proton), trackCuts.getNsigmaTOF(part, o2::track::PID::Proton), trackCuts.getNsigmaTPC(part, o2::track::PID::Pion), trackCuts.getNsigmaTOF(part, o2::track::PID::Pion), trackCuts.getNsigmaTPC(part, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(part, o2::track::PID::Kaon))) {
            mcRecoRegistry.fill(HIST("hMcRecPipPt"), part.pt());
            if (mcpart.partOriginMCTruth() == aod::femtouniverse_mc_particle::ParticleOriginMCTruth::kPrimary) {
              mcRecoRegistry.fill(HIST("hMcRecPipPtGenEtaGen"), mcpart.pt(), mcpart.eta());
            }
          }
          if ((mcpart.pdgMCTruth() == PDG_t::kKPlus) && (part.pt() > ConfTrack.confTrackLowPtCut) && (part.pt() < ConfTrack.confTrackHighPtCut) && isParticleNSigma(part.p(), trackCuts.getNsigmaTPC(part, o2::track::PID::Proton), trackCuts.getNsigmaTOF(part, o2::track::PID::Proton), trackCuts.getNsigmaTPC(part, o2::track::PID::Pion), trackCuts.getNsigmaTOF(part, o2::track::PID::Pion), trackCuts.getNsigmaTPC(part, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(part, o2::track::PID::Kaon))) {
            mcRecoRegistry.fill(HIST("hMcRecKpPt"), part.pt());
            if (mcpart.partOriginMCTruth() == aod::femtouniverse_mc_particle::ParticleOriginMCTruth::kPrimary) {
              mcRecoRegistry.fill(HIST("hMcRecKpPtGenEtaGen"), mcpart.pt(), mcpart.eta());
            }
          }
        } else if (part.sign() < 0) {
          if ((mcpart.pdgMCTruth() == -PDG_t::kProton) && (part.pt() > ConfTrack.confTrackLowPtCut) && (part.pt() < ConfTrack.confTrackHighPtCut) && isParticleNSigma(part.p(), trackCuts.getNsigmaTPC(part, o2::track::PID::Proton), trackCuts.getNsigmaTOF(part, o2::track::PID::Proton), trackCuts.getNsigmaTPC(part, o2::track::PID::Pion), trackCuts.getNsigmaTOF(part, o2::track::PID::Pion), trackCuts.getNsigmaTPC(part, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(part, o2::track::PID::Kaon))) {
            mcRecoRegistry.fill(HIST("hMcRecAntiPrPt"), part.pt());
            if (mcpart.partOriginMCTruth() == aod::femtouniverse_mc_particle::ParticleOriginMCTruth::kPrimary) {
              mcRecoRegistry.fill(HIST("hMcRecAntiPrPtGenEtaGen"), mcpart.pt(), mcpart.eta());
            }
          }
          if ((mcpart.pdgMCTruth() == -PDG_t::kPiPlus) && (part.pt() > ConfTrack.confTrackLowPtCut) && (part.pt() < ConfTrack.confTrackHighPtCut) && isParticleNSigma(part.p(), trackCuts.getNsigmaTPC(part, o2::track::PID::Proton), trackCuts.getNsigmaTOF(part, o2::track::PID::Proton), trackCuts.getNsigmaTPC(part, o2::track::PID::Pion), trackCuts.getNsigmaTOF(part, o2::track::PID::Pion), trackCuts.getNsigmaTPC(part, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(part, o2::track::PID::Kaon))) {
            mcRecoRegistry.fill(HIST("hMcRecPimPt"), part.pt());
            if (mcpart.partOriginMCTruth() == aod::femtouniverse_mc_particle::ParticleOriginMCTruth::kPrimary) {
              mcRecoRegistry.fill(HIST("hMcRecPimPtGenEtaGen"), mcpart.pt(), mcpart.eta());
            }
          }
          if ((mcpart.pdgMCTruth() == -PDG_t::kKPlus) && (part.pt() > ConfTrack.confTrackLowPtCut) && (part.pt() < ConfTrack.confTrackHighPtCut) && isParticleNSigma(part.p(), trackCuts.getNsigmaTPC(part, o2::track::PID::Proton), trackCuts.getNsigmaTOF(part, o2::track::PID::Proton), trackCuts.getNsigmaTPC(part, o2::track::PID::Pion), trackCuts.getNsigmaTOF(part, o2::track::PID::Pion), trackCuts.getNsigmaTPC(part, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(part, o2::track::PID::Kaon))) {
            mcRecoRegistry.fill(HIST("hMcRecKmPt"), part.pt());
            if (mcpart.partOriginMCTruth() == aod::femtouniverse_mc_particle::ParticleOriginMCTruth::kPrimary) {
              mcRecoRegistry.fill(HIST("hMcRecKmPtGenEtaGen"), mcpart.pt(), mcpart.eta());
            }
          }
        }
      } else if ((part.partType() == aod::femtouniverseparticle::ParticleType::kD0) && (part.pt() > ConfDmesons.confMinPtD0D0barReco) && (part.pt() < ConfDmesons.confMaxPtD0D0barReco)) {
        if (mcpart.pdgMCTruth() == ConfDmesons.confPDGCodeD0) {
          mcRecoRegistry.fill(HIST("hMcRecD0"), part.pt(), part.eta());
          mcRecoRegistry.fill(HIST("hMcRecD0Pt"), part.pt());
          mcRecoRegistry.fill(HIST("hMcRecD0Phi"), part.phi());
          if (part.tpcNClsFound() == 0) { // prompt candidates
            mcRecoRegistry.fill(HIST("hMcRecD0Prompt"), part.pt(), part.eta());
            mcRecoRegistry.fill(HIST("hMcRecD0PromptPt"), part.pt());
            mcRecoRegistry.fill(HIST("hMcRecD0PromptPhi"), part.phi());
          } else if (part.tpcNClsFound() == 1) { // non-prompt candidates
            mcRecoRegistry.fill(HIST("hMcRecD0NonPrompt"), part.pt(), part.eta());
            mcRecoRegistry.fill(HIST("hMcRecD0NonPromptPt"), part.pt());
            mcRecoRegistry.fill(HIST("hMcRecD0NonPromptPhi"), part.phi());
          }
        } else if (mcpart.pdgMCTruth() == ConfDmesons.confPDGCodeD0bar) {
          mcRecoRegistry.fill(HIST("hMcRecD0bar"), part.pt(), part.eta());
          mcRecoRegistry.fill(HIST("hMcRecD0barPt"), part.pt());
          mcRecoRegistry.fill(HIST("hMcRecD0barPhi"), part.phi());
          if (part.tpcNClsFound() == 0) { // prompt candidates
            mcRecoRegistry.fill(HIST("hMcRecD0barPrompt"), part.pt(), part.eta());
            mcRecoRegistry.fill(HIST("hMcRecD0barPromptPt"), part.pt());
          } else if (part.tpcNClsFound() == 1) { // non-prompt candidates
            mcRecoRegistry.fill(HIST("hMcRecD0barNonPrompt"), part.pt(), part.eta());
            mcRecoRegistry.fill(HIST("hMcRecD0barNonPromptPt"), part.pt());
          }
        }
      }
      // if (isParticleNSigma(part.p(), trackCuts.getNsigmaTPC(part, o2::track::PID::Proton), trackCuts.getNsigmaTOF(part, o2::track::PID::Proton), trackCuts.getNsigmaTPC(part, o2::track::PID::Pion), trackCuts.getNsigmaTOF(part, o2::track::PID::Pion), trackCuts.getNsigmaTPC(part, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(part, o2::track::PID::Kaon))) {
      // hTrackDCA.fillQA<true, true>(part);
      //} This part required change
    }
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackD0, processMcReco, "Process MC reco data", false);

  void processMcRecoD0InvMass(FemtoMcRecoParticles const& recoParts)
  {
    for (auto const& part : recoParts) {
      // filling the histograms for identified hadrons
      if ((part.partType() == aod::femtouniverseparticle::ParticleType::kD0) && (part.pt() > ConfDmesons.confMinPtD0D0barReco) && (part.pt() < ConfDmesons.confMaxPtD0D0barReco)) {
        // getting the efficiency value
        if (doEfficiencyCorr) {
          weight = efficiencyCalculator.getWeight(ParticleNo::TWO, part.pt());
        } else {
          weight = 1.0f;
        }
        // filling the inv. mass histograms
        if (part.mLambda() > 0) {
          if (part.sign() == o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK) {
            mcRecoRegistry.fill(HIST("hMassVsPtD0Sig"), part.mLambda(), part.pt(), weight);
          } else if (part.sign() == -o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK) {
            mcRecoRegistry.fill(HIST("hMassVsPtD0Refl"), part.mLambda(), part.pt(), weight);
          } else {
            mcRecoRegistry.fill(HIST("hMassVsPtD0Bkg"), part.mLambda(), part.pt(), weight);
          }
          if (part.tpcNClsFound() == 0) { // prompt candidates
            mcRecoRegistry.fill(HIST("hMassVsPtD0Prompt"), part.mLambda(), part.pt(), weight);
          } else if (part.tpcNClsFound() == 1) { // non-prompt candidates
            mcRecoRegistry.fill(HIST("hMassVsPtD0NonPrompt"), part.mLambda(), part.pt(), weight);
          }
        } else if (part.mAntiLambda() > 0) {
          if (part.sign() == -o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK) {
            mcRecoRegistry.fill(HIST("hMassVsPtD0barSig"), part.mAntiLambda(), part.pt(), weight);
          } else if (part.sign() == o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK) {
            mcRecoRegistry.fill(HIST("hMassVsPtD0barRefl"), part.mAntiLambda(), part.pt(), weight);
          } else {
            mcRecoRegistry.fill(HIST("hMassVsPtD0barBkg"), part.mAntiLambda(), part.pt(), weight);
          }
          if (part.tpcNClsFound() == 0) { // prompt candidates
            mcRecoRegistry.fill(HIST("hMassVsPtD0barPrompt"), part.mAntiLambda(), part.pt(), weight);
          } else if (part.tpcNClsFound() == 1) { // non-prompt candidates
            mcRecoRegistry.fill(HIST("hMassVsPtD0barNonPrompt"), part.mAntiLambda(), part.pt(), weight);
          }
        }
      }
    }
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackD0, processMcRecoD0InvMass, "Process MC reco inv. mass histograms for D0", false);

  void processMcTruth(aod::FDParticles const& parts)
  {
    for (auto const& part : parts) {
      if (part.partType() != uint8_t(aod::femtouniverseparticle::ParticleType::kMCTruthTrack)) {
        continue;
      }

      int pdgCode = static_cast<int>(part.tempFitVar());
      int8_t hfFlagMcGen = static_cast<int8_t>(part.mLambda());
      const auto& pdgParticle = pdgMC->GetParticle(pdgCode);
      if (!pdgParticle) {
        continue;
      }

      if (pdgParticle->Charge() > 0.0) {
        mcTruthRegistry.fill(HIST("hMcGenAllPositivePt"), part.pt());
      }
      if (pdgCode == PDG_t::kPiPlus) {
        mcTruthRegistry.fill(HIST("hMcGenPipPtVsEta"), part.pt(), part.eta());
        mcTruthRegistry.fill(HIST("hMcGenPipPt"), part.pt());
      }
      if (pdgCode == PDG_t::kKPlus) {
        mcTruthRegistry.fill(HIST("hMcGenKpPtVsEta"), part.pt(), part.eta());
        mcTruthRegistry.fill(HIST("hMcGenKpPt"), part.pt());
      }
      if (pdgCode == o2::constants::physics::Pdg::kD0) {
        if (std::abs(hfFlagMcGen) == o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK) {
          mcTruthRegistry.fill(HIST("hMcGenD0"), part.pt(), part.eta());
          mcTruthRegistry.fill(HIST("hMcGenD0Pt"), part.pt());
          if (static_cast<int8_t>(part.mAntiLambda()) == genPromptD0) {
            mcTruthRegistry.fill(HIST("hMcGenD0Prompt"), part.pt(), part.eta());
            mcTruthRegistry.fill(HIST("hMcGenD0PromptPt"), part.pt());
          } else if (static_cast<int8_t>(part.mAntiLambda()) == genNonPromptD0) {
            mcTruthRegistry.fill(HIST("hMcGenD0NonPrompt"), part.pt(), part.eta());
            mcTruthRegistry.fill(HIST("hMcGenD0NonPromptPt"), part.pt());
          }
        }
      }
      if (pdgCode == PDG_t::kProton) {
        mcTruthRegistry.fill(HIST("hMcGenPrPtVsEta"), part.pt(), part.eta());
        mcTruthRegistry.fill(HIST("hMcGenPrPt"), part.pt());
      }

      if (pdgParticle->Charge() < 0.0) {
        mcTruthRegistry.fill(HIST("hMcGenAllNegativePt"), part.pt());
      }
      if (pdgCode == PDG_t::kPiMinus) {
        mcTruthRegistry.fill(HIST("hMcGenPimPtVsEta"), part.pt(), part.eta());
        mcTruthRegistry.fill(HIST("hMcGenPimPt"), part.pt());
      }
      if (pdgCode == PDG_t::kKMinus) {
        mcTruthRegistry.fill(HIST("hMcGenKmPtVsEta"), part.pt(), part.eta());
        mcTruthRegistry.fill(HIST("hMcGenKmPt"), part.pt());
      }
      if (pdgCode == o2::constants::physics::Pdg::kD0Bar) {
        if (std::abs(hfFlagMcGen) == o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK) {
          mcTruthRegistry.fill(HIST("hMcGenD0bar"), part.pt(), part.eta());
          mcTruthRegistry.fill(HIST("hMcGenD0barPt"), part.pt());
          if (static_cast<int8_t>(part.mAntiLambda()) == genPromptD0) {
            mcTruthRegistry.fill(HIST("hMcGenD0barPrompt"), part.pt(), part.eta());
            mcTruthRegistry.fill(HIST("hMcGenD0barPromptPt"), part.pt());
          } else if (static_cast<int8_t>(part.mAntiLambda()) == genNonPromptD0) {
            mcTruthRegistry.fill(HIST("hMcGenD0barNonPrompt"), part.pt(), part.eta());
            mcTruthRegistry.fill(HIST("hMcGenD0barNonPromptPt"), part.pt());
          }
        }
      }
      if (pdgCode == -PDG_t::kProton) {
        mcTruthRegistry.fill(HIST("hMcGenAntiPrPtVsEta"), part.pt(), part.eta());
        mcTruthRegistry.fill(HIST("hMcGenAntiPrPt"), part.pt());
      }
    }
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackD0, processMcTruth, "Process MC truth data", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<FemtoUniversePairTaskTrackD0>(cfgc),
  };
  return workflow;
}
