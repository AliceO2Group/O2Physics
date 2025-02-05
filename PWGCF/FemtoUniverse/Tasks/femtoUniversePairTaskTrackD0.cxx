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

/// \file femtoUniversePairTaskTrackD0.cxx
/// \brief Tasks that reads the track tables and D0/D0bar mesons
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de
/// \author Georgios Mantzaridis, TU München, georgios.mantzaridis@tum.de
/// \author Anton Riedel, TU München, anton.riedel@tum.de
/// \author Zuzanna Chochulska, WUT Warsaw & CTU Prague, zchochul@cern.ch
/// \author Katarzyna Gwiździel, WUT Warsaw, katarzyna.gwizdziel@cern.ch

#include <vector>
#include <string>
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/StepTHn.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "ReconstructionDataFormats/PID.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/Core/RecoDecay.h"

#include "PWGCF/FemtoUniverse/DataModel/FemtoDerived.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseParticleHisto.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseEventHisto.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniversePairCleaner.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseFemtoContainer.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseAngularContainer.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseDetaDphiStar.h"
#include "PWGCF/FemtoUniverse/Core/femtoUtils.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseTrackSelection.h"

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::analysis;
using namespace o2::analysis::femto_universe;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

namespace
{
static constexpr int kNPart = 2;
static constexpr int kNCuts = 5;
static const std::vector<std::string> partNames{"D0", "Track"};
static const std::vector<std::string> cutNames{"MaxPt", "PIDthr", "nSigmaTPC", "nSigmaTPCTOF", "MaxP"};
static const float cutsTable[kNPart][kNCuts]{{4.05f, 1.f, 3.f, 3.f, 100.f}, {4.05f, 1.f, 3.f, 3.f, 100.f}};
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

  using FemtoFullParticles = soa::Join<aod::FDParticles, aod::FDExtParticles>;
  SliceCache cache;
  Preslice<FemtoFullParticles> perCol = aod::femtouniverseparticle::fdCollisionId;

  /// Table for both particles
  struct : o2::framework::ConfigurableGroup {
    Configurable<float> confNsigmaCombinedProton{"confNsigmaCombinedProton", 3.0, "TPC and TOF Proton Sigma (combined) for momentum > 0.5"};
    Configurable<float> confNsigmaTPCProton{"confNsigmaTPCProton", 3.0, "TPC Proton Sigma for momentum < 0.5"};
    Configurable<float> confNsigmaCombinedPion{"confNsigmaCombinedPion", 3.0, "TPC and TOF Pion Sigma (combined) for momentum > 0.5"};
    Configurable<float> confNsigmaTPCPion{"confNsigmaTPCPion", 3.0, "TPC Pion Sigma for momentum < 0.5"};

    Configurable<LabeledArray<float>> confCutTable{"confCutTable", {cutsTable[0], kNPart, kNCuts, partNames, cutNames}, "Particle selections"};
    Configurable<int> confNspecies{"confNspecies", 2, "Number of particle spieces with PID info"};
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
  } ConfTrack;

  /// Particle 2 --- D0/D0bar meson
  struct : o2::framework::ConfigurableGroup {
    Configurable<int> confPDGCodeD0{"confPDGCodeD0", 421, "D0 meson - PDG code"};
    Configurable<int> confPDGCodeD0bar{"confPDGCodeD0bar", -421, "D0bar meson - PDG code"};
    Configurable<float> confMinPtD0D0bar{"confMinPtD0D0bar", 1.0, "D0/D0bar sel. - min. pT"};
    Configurable<float> confMaxPtD0D0bar{"confMaxPtD0D0bar", 3.0, "D0/D0bar sel. - max. pT"};
    Configurable<float> confMinInvMassD0D0bar{"confMinInvMassD0D0bar", 1.65, "D0/D0bar sel. - min. invMass"};
    Configurable<float> confMaxInvMassD0D0bar{"confMaxInvMassD0D0bar", 2.05, "D0/D0bar sel. - max. invMass"};
    Configurable<float> confMaxProbMlClass1Bg{"confMaxProbMlClass1Bg", 0.4, "ML: max prob. that D0/D0bar cand. is from the backgound"};
  } ConfDmesons;

  struct : o2::framework::ConfigurableGroup {
    Configurable<float> confSignalRegionMin{"confSignalRegionMin", 1.810, "Min. inv. mass for D0/D0bar in the signal region"};
    Configurable<float> confSignalRegionMax{"confSignalRegionMax", 1.922, "Max. inv. mass for D0/D0bar in the signal region"};
    Configurable<float> confMinInvMassLeftSB{"confMinInvMassLeftSB", 1.642, "Min. inv. mass for D0/D0bar in the left sideband region"};
    Configurable<float> confMaxInvMassLeftSB{"confMaxInvMassLeftSB", 1.754, "Max. inv. mass for D0/D0bar in the left sideband region"};
    Configurable<float> confMinInvMassRightSB{"confMinInvMassRightSB", 1.978, "Min. inv. mass for D0/D0bar in the right sideband region"};
    Configurable<float> confMaxInvMassRightSB{"confMaxInvMassRightSB", 2.090, "Max. inv. mass for D0/D0bar in the right sideband region"};
  } ConfD0D0barSideBand;

  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_d0_to_pi_k::vecBinsPt}, "pT bin limits"};
  Configurable<uint8_t> confChooseD0trackCorr{"confChooseD0trackCorr", 2, "If 0 - only D0s, 1 - only D0bars, 2 - D0/D0bar (one mass hypo.)"};

  /// Partitions for particle 1
  Partition<FemtoFullParticles> partsTrack = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kTrack)) && (aod::femtouniverseparticle::sign == int8_t(ConfTrack.confTrackSign)) && (aod::femtouniverseparticle::pt > ConfTrack.confTrackLowPtCut) && (aod::femtouniverseparticle::pt < ConfTrack.confTrackHighPtCut);
  Partition<soa::Join<aod::FDParticles, aod::FDMCLabels>> partsTrackMC = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kTrack));

  /// Partitions for particle 2
  /// Partition with all D0/D0bar mesons (which pass double mass hypothesis)
  // Partition<FemtoFullParticles> partsAllDmesons = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kD0)) && (aod::femtouniverseparticle::mLambda > 0.0f) && (aod::femtouniverseparticle::mAntiLambda > 0.0f);
  /// Partition with D0/D0bar candidates, which pass only one mass hypothesis
  Partition<FemtoFullParticles> partsOnlyD0D0bar = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kD0)) && (aod::femtouniverseparticle::mLambda < 0.0f || aod::femtouniverseparticle::mAntiLambda < 0.0f) && (aod::femtouniverseparticle::tempFitVar < ConfDmesons.confMaxProbMlClass1Bg);
  /// Partition with D0 mesons only (one mass hypothesis)
  Partition<FemtoFullParticles> partsD0s = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kD0)) && (aod::femtouniverseparticle::mLambda > ConfDmesons.confMinInvMassD0D0bar) && (aod::femtouniverseparticle::mLambda < ConfDmesons.confMaxInvMassD0D0bar) && (aod::femtouniverseparticle::mAntiLambda < 0.0f) && (aod::femtouniverseparticle::pt > ConfDmesons.confMinPtD0D0bar) && (aod::femtouniverseparticle::pt < ConfDmesons.confMaxPtD0D0bar) && (aod::femtouniverseparticle::tempFitVar < ConfDmesons.confMaxProbMlClass1Bg);
  /// Partition with D0bar mesons only (one mass hypothesis)
  Partition<FemtoFullParticles> partsD0bars = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kD0)) && (aod::femtouniverseparticle::mLambda < 0.0f) && (aod::femtouniverseparticle::mAntiLambda > ConfDmesons.confMinInvMassD0D0bar) && (aod::femtouniverseparticle::mAntiLambda < ConfDmesons.confMaxInvMassD0D0bar) && (aod::femtouniverseparticle::pt > ConfDmesons.confMinPtD0D0bar) && (aod::femtouniverseparticle::pt < ConfDmesons.confMaxPtD0D0bar) && (aod::femtouniverseparticle::tempFitVar < ConfDmesons.confMaxProbMlClass1Bg);
  /// Partition for D0/D0bar mesons from MC
  Partition<soa::Join<aod::FDParticles, aod::FDMCLabels>> partsD0D0barMC = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kD0));

  /// Partition for D0/D0bar daughters
  Partition<FemtoFullParticles> partsDmesonsChildren = aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kD0Child);

  /// Histogramming for particle 1
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kTrack, 2> trackHistoPartTrack;

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
  ConfigurableAxis confInvMassFinerBins{"confInvMassFinerBins", {120, 1.5848, 2.1848}, "finer binning of inv. mass"};

  Configurable<bool> confIsCPR{"confIsCPR", true, "Close Pair Rejection"};
  Configurable<bool> confCPRPlotPerRadii{"confCPRPlotPerRadii", false, "Plot CPR per radii"};
  Configurable<float> confCPRdeltaPhiCutMax{"confCPRdeltaPhiCutMax", 0.0, "Delta Phi max cut for Close Pair Rejection"};
  Configurable<float> confCPRdeltaPhiCutMin{"confCPRdeltaPhiCutMin", 0.0, "Delta Phi min cut for Close Pair Rejection"};
  Configurable<float> confCPRdeltaEtaCutMax{"confCPRdeltaEtaCutMax", 0.0, "Delta Eta max cut for Close Pair Rejection"};
  Configurable<float> confCPRdeltaEtaCutMin{"confCPRdeltaEtaCutMin", 0.0, "Delta Eta min cut for Close Pair Rejection"};
  Configurable<float> confCPRChosenRadii{"confCPRChosenRadii", 0.80, "Delta Eta cut for Close Pair Rejection"};

  Configurable<bool> applyMLOpt{"applyMLOpt", false, "Enable for ML selection optimization"};

  FemtoUniverseFemtoContainer<femto_universe_femto_container::EventType::same, femto_universe_femto_container::Observable::kstar> sameEventFemtoCont;
  FemtoUniverseFemtoContainer<femto_universe_femto_container::EventType::mixed, femto_universe_femto_container::Observable::kstar> mixedEventFemtoCont;
  FemtoUniverseAngularContainer<femto_universe_angular_container::EventType::same, femto_universe_angular_container::Observable::kstar> sameEventAngularCont;
  FemtoUniverseAngularContainer<femto_universe_angular_container::EventType::mixed, femto_universe_angular_container::Observable::kstar> mixedEventAngularCont;
  FemtoUniversePairCleaner<aod::femtouniverseparticle::ParticleType::kTrack, aod::femtouniverseparticle::ParticleType::kD0> pairCleaner;
  FemtoUniverseDetaDphiStar<aod::femtouniverseparticle::ParticleType::kTrack, aod::femtouniverseparticle::ParticleType::kD0> pairCloseRejection;
  FemtoUniverseTrackSelection trackCuts;

  /// Histogram output
  HistogramRegistry qaRegistry{"TrackQA", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry resultRegistry{"Correlations", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry mixQaRegistry{"mixQaRegistry", {}, OutputObjHandlingPolicy::AnalysisObject};

  HistogramRegistry registry{"registry",
                             {{"hInvMassD0", ";#it{M}(K^{-}#pi^{+}) (GeV/#it{c}^{2});counts", {HistType::kTH1F, {confInvMassBins}}},
                              {"hInvMassD0bar", ";#it{M}(#pi^{-}K^{+}) (GeV/#it{c}^{2});counts", {HistType::kTH1F, {confInvMassBins}}},
                              {"hPtDmesonCand", "2-prong candidates;#it{p}_{T} (GeV/#it{c});counts", {HistType::kTH1F, {confPtBins}}},
                              {"hPtD0", "D^{0} cand.;#it{p}_{T} (GeV/#it{c});counts", {HistType::kTH1F, {confPtBins}}},
                              {"hPtD0bar", "#bar{D^{0}};#it{p}_{T} (GeV/#it{c});counts", {HistType::kTH1F, {confPtBins}}},
                              {"hPtD0D0bar", "#bar{D^{0}};#it{p}_{T} (GeV/#it{c});counts", {HistType::kTH1F, {confPtBins}}},
                              {"hPhiDmesonCand", ";#varphi (rad);counts", {HistType::kTH1F, {{80, 0., o2::constants::math::TwoPI}}}},
                              {"hPhiD0", ";#varphi (rad);counts", {HistType::kTH1F, {{80, 0., o2::constants::math::TwoPI}}}},
                              {"hPhiD0bar", ";#varphi (rad);counts", {HistType::kTH1F, {{80, 0., o2::constants::math::TwoPI}}}},
                              {"hEtaDmesonCand", ";#eta ;counts", {HistType::kTH1F, {{200, -1., 1.}}}},
                              {"hEtaD0", ";#eta ;counts", {HistType::kTH1F, {{200, -1., 1.}}}},
                              {"hEtaD0bar", ";#eta ;counts", {HistType::kTH1F, {{200, -1., 1.}}}},
                              {"hDecayLengthD0", ";decay length (cm);counts", {HistType::kTH1F, {{800, 0., 4.}}}},
                              {"hDecayLengthD0bar", ";decay length (cm);counts", {HistType::kTH1F, {{800, 0., 4.}}}},
                              {"hPtDaughters", ";#it{p}_{T} (GeV/#it{c});counts", {HistType::kTH1F, {{300, 0., 12.}}}},
                              {"hSignDaughters", ";sign ;counts", {HistType::kTH1F, {{10, -2.5, 2.5}}}},
                              {"hbetaDaughters", "; p (GeV/#it{c}); TOF #beta", {HistType::kTH2F, {{300, 0., 15.}, {200, 0., 2.}}}},
                              {"hdEdxDaughters", "; p (GeV/#it{c}); TPC dE/dx (KeV/cm)", {HistType::kTH2F, {{300, 0., 15.}, {500, 0., 500.}}}},
                              {"hDCAxyDaughters", "; #it{DCA}_{xy} (cm); counts", {HistType::kTH1F, {{140, 0., 0.14}}}},
                              {"hDCAzDaughters", "; #it{DCA}_{z} (cm); counts", {HistType::kTH1F, {{140, 0., 0.14}}}}}};

  // PID for protons
  bool isProtonNSigma(float mom, float nsigmaTPCPr, float nsigmaTOFPr) // previous version from: https://github.com/alisw/AliPhysics/blob/master/PWGCF/FEMTOSCOPY/AliFemtoUser/AliFemtoMJTrackCut.cxx
  {
    if (mom < 0.5) {
      if (std::abs(nsigmaTPCPr) < ConfBothTracks.confNsigmaTPCProton) {
        return true;
      } else {
        return false;
      }
    } else if (mom > 0.4) {
      if (std::hypot(nsigmaTOFPr, nsigmaTPCPr) < ConfBothTracks.confNsigmaCombinedProton) {
        return true;
      } else {
        return false;
      }
    }
    return false;
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
    } else if (mom > 1.5) {
      if ((std::abs(nsigmaTOFK) < 2.0) && (std::abs(nsigmaTPCK) < 3.0)) {
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
    //|nsigma_TPC| < 3 for p < 0.5 GeV/c
    //|nsigma_combined| < 3 for p > 0.5

    // using configurables:
    // confNsigmaTPCPion -> TPC Kaon Sigma for momentum < 0.5
    // confNsigmaCombinedPion -> TPC and TOF Pion Sigma (combined) for momentum > 0.5
    if (true) {
      if (mom < 0.5) {
        if (std::abs(nsigmaTPCPi) < ConfBothTracks.confNsigmaTPCPion) {
          return true;
        } else {
          return false;
        }
      } else if (mom > 0.5) {
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
        return isKaonNSigma(mom, nsigmaTPCK, nsigmaTOFK);
        break;
      default:
        return false;
    }
  }

  void init(InitContext&)
  {
    eventHisto.init(&qaRegistry);
    qaRegistry.add("QA_D0D0barSelection/hInvMassD0", ";#it{M}(K^{-}#pi^{+}) (GeV/#it{c}^{2});counts", kTH1F, {confInvMassBins});
    qaRegistry.add("QA_D0D0barSelection/hPtD0", "D^{0} cand.;#it{p}_{T} (GeV/#it{c});counts", kTH1F, {confPtBins});
    qaRegistry.add("QA_D0D0barSelection/hInvMassD0bar", ";#it{M}(K^{-}#pi^{+}) (GeV/#it{c}^{2});counts", kTH1F, {confInvMassBins});
    qaRegistry.add("QA_D0D0barSelection/hPtD0bar", "#bar{D^{0}} cand.;#it{p}_{T} (GeV/#it{c});counts", kTH1F, {confPtBins});
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

    trackHistoPartD0D0bar.init(&qaRegistry, confTempFitVarpTBins, confTempFitVarInvMassBins, ConfBothTracks.confIsMC, ConfDmesons.confPDGCodeD0);
    if (!ConfTrack.confIsSame) {
      trackHistoPartTrack.init(&qaRegistry, confTempFitVarpTBins, confTempFitVarBins, ConfBothTracks.confIsMC, ConfTrack.confPDGCodeTrack);
    }

    mixQaRegistry.add("MixingQA/hSECollisionBins", ";bin;Entries", kTH1F, {{120, -0.5, 119.5}});
    mixQaRegistry.add("MixingQA/hMECollisionBins", ";bin;Entries", kTH1F, {{120, -0.5, 119.5}});

    sameEventFemtoCont.init(&resultRegistry, confkstarBins, confMultBins, confkTBins, confmTBins, confmultBins3D, confmTBins3D, ConfBothTracks.confIsMC, ConfBothTracks.confUse3D);
    mixedEventFemtoCont.init(&resultRegistry, confkstarBins, confMultBins, confkTBins, confmTBins, confmultBins3D, confmTBins3D, ConfBothTracks.confIsMC, ConfBothTracks.confUse3D);
    sameEventAngularCont.init(&resultRegistry, confkstarBins, confMultBins, confkTBins, confmTBins, confmultBins3D, confmTBins3D, ConfBothTracks.confEtaBins, ConfBothTracks.confPhiBins, ConfBothTracks.confIsMC, ConfBothTracks.confUse3D);
    mixedEventAngularCont.init(&resultRegistry, confkstarBins, confMultBins, confkTBins, confmTBins, confmultBins3D, confmTBins3D, ConfBothTracks.confEtaBins, ConfBothTracks.confPhiBins, ConfBothTracks.confIsMC, ConfBothTracks.confUse3D);

    sameEventFemtoCont.setPDGCodes(ConfDmesons.confPDGCodeD0, ConfTrack.confPDGCodeTrack);
    mixedEventFemtoCont.setPDGCodes(ConfDmesons.confPDGCodeD0, ConfTrack.confPDGCodeTrack);
    sameEventAngularCont.setPDGCodes(ConfDmesons.confPDGCodeD0, ConfTrack.confPDGCodeTrack);
    mixedEventAngularCont.setPDGCodes(ConfDmesons.confPDGCodeD0, ConfTrack.confPDGCodeTrack);

    pairCleaner.init(&qaRegistry);
    if (confIsCPR.value) {
      pairCloseRejection.init(&resultRegistry, &qaRegistry, confCPRdeltaPhiCutMin.value, confCPRdeltaPhiCutMax.value, confCPRdeltaEtaCutMin.value, confCPRdeltaEtaCutMax.value, confCPRChosenRadii.value, confCPRPlotPerRadii.value);
    }

    vPIDTrack = ConfTrack.confPIDTrack.value;
    kNsigma = ConfBothTracks.confTrkPIDnSigmaMax.value;

    // D0/D0bar histograms
    auto vbins = (std::vector<double>)binsPt;
    registry.add("D0D0bar_oneMassHypo/hMassVsPt", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {confInvMassBins, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("D0D0bar_oneMassHypo/hMassVsPtFinerBinning", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {confInvMassFinerBins, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDeltaPhiSigSig", "SxS correlation;#Delta#varphi (rad);counts", {HistType::kTH1F, {{10, 0.0, o2::constants::math::PI}}});
    registry.add("hDeltaPhiD0BgD0barSig", "B(D0)x S(D0bar) correlation;#Delta#varphi (rad);counts", {HistType::kTH1F, {{10, 0.0, o2::constants::math::PI}}});
    registry.add("hDeltaPhiD0SigD0barBg", "S(D0)x B(D0bar) correlation;#Delta#varphi (rad);counts", {HistType::kTH1F, {{10, 0.0, o2::constants::math::PI}}});
    registry.add("hDeltaPhiBgBg", "BxB correlation;#Delta#varphi (rad);counts", {HistType::kTH1F, {{10, 0.0, o2::constants::math::PI}}});
    registry.add("hPtCand1VsPtCand2", "2-prong candidates;#it{p}_{T} (GeV/#it{c});#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDeltaEtaDeltaPhi", "2-prong candidates;#Delta #eta;#Delta #varphi (rad)", {HistType::kTH2F, {{29, -2., 2.}, {29, 0.0, o2::constants::math::PI}}});
    if (applyMLOpt) {
      registry.add("D0D0bar_MLSel/hMassVsPt001", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {confInvMassBins, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("D0D0bar_MLSel/hMassVsPt0015", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {confInvMassBins, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("D0D0bar_MLSel/hMassVsPt002", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {confInvMassBins, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("D0D0bar_MLSel/hMassVsPt0025", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {confInvMassBins, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("D0D0bar_MLSel/hMassVsPt003", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {confInvMassBins, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("D0D0bar_MLSel/hMassVsPt004", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {confInvMassBins, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("D0D0bar_MLSel/hMassVsPt005", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {confInvMassBins, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("D0D0bar_MLSel/hMassVsPt006", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {confInvMassBins, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("D0D0bar_MLSel/hMassVsPt007", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {confInvMassBins, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("D0D0bar_MLSel/hMassVsPt008", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {confInvMassBins, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("D0D0bar_MLSel/hMassVsPt009", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {confInvMassBins, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("D0D0bar_MLSel/hMassVsPt01", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {confInvMassBins, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("D0D0bar_MLSel/hMassVsPt015", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {confInvMassBins, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("D0D0bar_MLSel/hMassVsPt02", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {confInvMassBins, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("D0D0bar_MLSel/hMassVsPt025", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {confInvMassBins, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
      registry.add("D0D0bar_MLSel/hMassVsPt03", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {confInvMassBins, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    }
  }

  template <typename CollisionType>
  void fillCollision(CollisionType col)
  {
    mixQaRegistry.fill(HIST("MixingQA/hSECollisionBins"), colBinning.getBin({col.posZ(), col.multNtr()}));
    eventHisto.fillQA(col);
  }

  void processD0MLOpt(o2::aod::FdCollision const& col, FemtoFullParticles const&)
  {
    auto groupD0D0barCands = partsOnlyD0D0bar->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);

    // loop over selected D0/D0bar candidates
    for (auto const& charmCand : groupD0D0barCands) {
      // D0 candidates
      if (charmCand.mLambda() > 0.0f && charmCand.mAntiLambda() < 0.0f) {
        if (charmCand.tempFitVar() < 0.01)
          registry.fill(HIST("D0D0bar_MLSel/hMassVsPt001"), charmCand.mLambda(), charmCand.pt());
        if (charmCand.tempFitVar() < 0.015)
          registry.fill(HIST("D0D0bar_MLSel/hMassVsPt0015"), charmCand.mLambda(), charmCand.pt());
        if (charmCand.tempFitVar() < 0.02)
          registry.fill(HIST("D0D0bar_MLSel/hMassVsPt002"), charmCand.mLambda(), charmCand.pt());
        if (charmCand.tempFitVar() < 0.025)
          registry.fill(HIST("D0D0bar_MLSel/hMassVsPt0025"), charmCand.mLambda(), charmCand.pt());
        if (charmCand.tempFitVar() < 0.03)
          registry.fill(HIST("D0D0bar_MLSel/hMassVsPt003"), charmCand.mLambda(), charmCand.pt());
        if (charmCand.tempFitVar() < 0.04)
          registry.fill(HIST("D0D0bar_MLSel/hMassVsPt004"), charmCand.mLambda(), charmCand.pt());
        if (charmCand.tempFitVar() < 0.05)
          registry.fill(HIST("D0D0bar_MLSel/hMassVsPt005"), charmCand.mLambda(), charmCand.pt());
        if (charmCand.tempFitVar() < 0.06)
          registry.fill(HIST("D0D0bar_MLSel/hMassVsPt006"), charmCand.mLambda(), charmCand.pt());
        if (charmCand.tempFitVar() < 0.07)
          registry.fill(HIST("D0D0bar_MLSel/hMassVsPt007"), charmCand.mLambda(), charmCand.pt());
        if (charmCand.tempFitVar() < 0.08)
          registry.fill(HIST("D0D0bar_MLSel/hMassVsPt008"), charmCand.mLambda(), charmCand.pt());
        if (charmCand.tempFitVar() < 0.09)
          registry.fill(HIST("D0D0bar_MLSel/hMassVsPt009"), charmCand.mLambda(), charmCand.pt());
        if (charmCand.tempFitVar() < 0.1)
          registry.fill(HIST("D0D0bar_MLSel/hMassVsPt01"), charmCand.mLambda(), charmCand.pt());
        if (charmCand.tempFitVar() < 0.15)
          registry.fill(HIST("D0D0bar_MLSel/hMassVsPt015"), charmCand.mLambda(), charmCand.pt());
        if (charmCand.tempFitVar() < 0.2)
          registry.fill(HIST("D0D0bar_MLSel/hMassVsPt02"), charmCand.mLambda(), charmCand.pt());
        if (charmCand.tempFitVar() < 0.25)
          registry.fill(HIST("D0D0bar_MLSel/hMassVsPt025"), charmCand.mLambda(), charmCand.pt());
        if (charmCand.tempFitVar() < 0.3)
          registry.fill(HIST("D0D0bar_MLSel/hMassVsPt03"), charmCand.mLambda(), charmCand.pt());
      }
      // DObar candidates
      if (charmCand.mLambda() < 0.0f && charmCand.mAntiLambda() > 0.0f) {
        if (charmCand.tempFitVar() < 0.01)
          registry.fill(HIST("D0D0bar_MLSel/hMassVsPt001"), charmCand.mAntiLambda(), charmCand.pt());
        if (charmCand.tempFitVar() < 0.015)
          registry.fill(HIST("D0D0bar_MLSel/hMassVsPt0015"), charmCand.mAntiLambda(), charmCand.pt());
        if (charmCand.tempFitVar() < 0.02)
          registry.fill(HIST("D0D0bar_MLSel/hMassVsPt002"), charmCand.mAntiLambda(), charmCand.pt());
        if (charmCand.tempFitVar() < 0.025)
          registry.fill(HIST("D0D0bar_MLSel/hMassVsPt0025"), charmCand.mAntiLambda(), charmCand.pt());
        if (charmCand.tempFitVar() < 0.03)
          registry.fill(HIST("D0D0bar_MLSel/hMassVsPt003"), charmCand.mAntiLambda(), charmCand.pt());
        if (charmCand.tempFitVar() < 0.04)
          registry.fill(HIST("D0D0bar_MLSel/hMassVsPt004"), charmCand.mAntiLambda(), charmCand.pt());
        if (charmCand.tempFitVar() < 0.05)
          registry.fill(HIST("D0D0bar_MLSel/hMassVsPt005"), charmCand.mAntiLambda(), charmCand.pt());
        if (charmCand.tempFitVar() < 0.06)
          registry.fill(HIST("D0D0bar_MLSel/hMassVsPt006"), charmCand.mAntiLambda(), charmCand.pt());
        if (charmCand.tempFitVar() < 0.07)
          registry.fill(HIST("D0D0bar_MLSel/hMassVsPt007"), charmCand.mAntiLambda(), charmCand.pt());
        if (charmCand.tempFitVar() < 0.08)
          registry.fill(HIST("D0D0bar_MLSel/hMassVsPt008"), charmCand.mAntiLambda(), charmCand.pt());
        if (charmCand.tempFitVar() < 0.09)
          registry.fill(HIST("D0D0bar_MLSel/hMassVsPt009"), charmCand.mAntiLambda(), charmCand.pt());
        if (charmCand.tempFitVar() < 0.1)
          registry.fill(HIST("D0D0bar_MLSel/hMassVsPt01"), charmCand.mAntiLambda(), charmCand.pt());
        if (charmCand.tempFitVar() < 0.15)
          registry.fill(HIST("D0D0bar_MLSel/hMassVsPt015"), charmCand.mAntiLambda(), charmCand.pt());
        if (charmCand.tempFitVar() < 0.2)
          registry.fill(HIST("D0D0bar_MLSel/hMassVsPt02"), charmCand.mAntiLambda(), charmCand.pt());
        if (charmCand.tempFitVar() < 0.25)
          registry.fill(HIST("D0D0bar_MLSel/hMassVsPt025"), charmCand.mAntiLambda(), charmCand.pt());
        if (charmCand.tempFitVar() < 0.3)
          registry.fill(HIST("D0D0bar_MLSel/hMassVsPt03"), charmCand.mAntiLambda(), charmCand.pt());
      }
    }
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackD0, processD0MLOpt, "Enable filling QA plots for ML D0/D0bar selection optimization", false);

  void processQAD0D0barSel(o2::aod::FdCollision const& col, FemtoFullParticles const&)
  {
    auto groupPartsD0s = partsD0s->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto groupPartsD0bars = partsD0bars->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);

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
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackD0, processQAD0D0barSel, "Enable filling QA plots for selected D0/D0bar cand.", true);

  void processD0mesons(o2::aod::FdCollision const& col, FemtoFullParticles const&)
  {
    auto groupPartsOnlyD0D0bar = partsOnlyD0D0bar->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto groupPartsD0D0barChildren = partsDmesonsChildren->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);

    // loop over D0/D0bar mesons (ONLY)
    for (auto const& d0d0bar : groupPartsOnlyD0D0bar) {

      registry.fill(HIST("hPtD0D0bar"), d0d0bar.pt());

      if (d0d0bar.mLambda() > 0.0f && d0d0bar.mAntiLambda() < 0.0f) {
        registry.fill(HIST("D0D0bar_oneMassHypo/hMassVsPt"), d0d0bar.mLambda(), d0d0bar.pt());
        registry.fill(HIST("hInvMassD0"), d0d0bar.mLambda());
        registry.fill(HIST("hPtD0"), d0d0bar.pt());
        registry.fill(HIST("hPhiD0"), d0d0bar.phi());
        registry.fill(HIST("hEtaD0"), d0d0bar.eta());
      }
      if (d0d0bar.mLambda() < 0.0f && d0d0bar.mAntiLambda() > 0.0f) {
        registry.fill(HIST("D0D0bar_oneMassHypo/hMassVsPt"), d0d0bar.mAntiLambda(), d0d0bar.pt());
        registry.fill(HIST("hInvMassD0bar"), d0d0bar.mAntiLambda());
        registry.fill(HIST("hPtD0bar"), d0d0bar.pt());
        registry.fill(HIST("hPhiD0bar"), d0d0bar.phi());
        registry.fill(HIST("hEtaD0bar"), d0d0bar.eta());
      }
    }

    // loop over D mesons childen
    for (auto const& daughD0D0bar : groupPartsD0D0barChildren) {
      registry.fill(HIST("hPtDaughters"), daughD0D0bar.pt());
      registry.fill(HIST("hSignDaughters"), daughD0D0bar.mLambda());
      // filling QA plots for D0 mesons' positive daughters (K+)
      if (daughD0D0bar.mLambda() == 1 && daughD0D0bar.mAntiLambda() == 1) {
        qaRegistry.fill(HIST("D0_pos_daugh/pt"), daughD0D0bar.pt());
        qaRegistry.fill(HIST("D0_pos_daugh/eta"), daughD0D0bar.eta());
        qaRegistry.fill(HIST("D0_pos_daugh/phi"), daughD0D0bar.phi());
      }
      // filling QA plots for D0 mesons' negative daughters (pi-)
      if (daughD0D0bar.mLambda() == -1 && daughD0D0bar.mAntiLambda() == 1) {
        qaRegistry.fill(HIST("D0_neg_daugh/pt"), daughD0D0bar.pt());
        qaRegistry.fill(HIST("D0_neg_daugh/eta"), daughD0D0bar.eta());
        qaRegistry.fill(HIST("D0_neg_daugh/phi"), daughD0D0bar.phi());
      }
      // filling QA plots for D0bar mesons' positive daughters (pi+)
      if (daughD0D0bar.mLambda() == 1 && daughD0D0bar.mAntiLambda() == -1) {
        qaRegistry.fill(HIST("D0bar_pos_daugh/pt"), daughD0D0bar.pt());
        qaRegistry.fill(HIST("D0bar_pos_daugh/eta"), daughD0D0bar.eta());
        qaRegistry.fill(HIST("D0bar_pos_daugh/phi"), daughD0D0bar.phi());
      }
      // filling QA plots for D0bar mesons' negative daughters (K-)
      if (daughD0D0bar.mLambda() == -1 && daughD0D0bar.mAntiLambda() == -1) {
        qaRegistry.fill(HIST("D0bar_neg_daugh/pt"), daughD0D0bar.pt());
        qaRegistry.fill(HIST("D0bar_neg_daugh/eta"), daughD0D0bar.eta());
        qaRegistry.fill(HIST("D0bar_neg_daugh/phi"), daughD0D0bar.phi());
      }
    }
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackD0, processD0mesons, "Enable processing D0 mesons", true);

  // D0-D0bar pair correlations (side-band methode)
  void processSideBand(o2::aod::FdCollision const& col, FemtoFullParticles const&)
  {
    auto groupPartsOnlyD0D0bar = partsOnlyD0D0bar->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);

    double deltaPhi = 0.0;
    double deltaEta = 0.0;

    // loop over D0/D0bar candidates (ONLY)
    for (auto const& cand1 : groupPartsOnlyD0D0bar) {
      // Check if the first candidate is D0 meson
      if (cand1.mLambda() < 0.0f && cand1.mAntiLambda() > 0.0f) {
        continue;
      }

      for (auto const& cand2 : groupPartsOnlyD0D0bar) {
        // Check if the second candidate is D0bar meson
        if (cand2.mLambda() > 0.0f && cand2.mAntiLambda() < 0.0f) {
          continue;
        }
        // deltaPhi = getDeltaPhi(cand1.phi(), cand2.phi());
        deltaPhi = wrapDeltaPhi0PI(cand1.phi(), cand2.phi());
        deltaEta = cand2.eta() - cand1.eta();

        // General histograms
        registry.fill(HIST("hPtCand1VsPtCand2"), cand1.pt(), cand2.pt());
        registry.fill(HIST("hDeltaEtaDeltaPhi"), deltaEta, deltaPhi);

        // ----------------------------------- Creating D0-D0bar pairs correlations ------------------------------------------------
        if (cand1.mLambda() > ConfD0D0barSideBand.confSignalRegionMin.value && cand1.mLambda() < ConfD0D0barSideBand.confSignalRegionMax.value) {
          // S(D0) x S(D0bar) correlation
          if (cand2.mAntiLambda() > ConfD0D0barSideBand.confSignalRegionMin.value && cand2.mAntiLambda() < ConfD0D0barSideBand.confSignalRegionMax.value) {
            registry.fill(HIST("hDeltaPhiSigSig"), deltaPhi);
          }
          // S(D0) x B(D0bar) correlation
          if ((cand2.mAntiLambda() > ConfD0D0barSideBand.confMinInvMassLeftSB.value && cand2.mAntiLambda() < ConfD0D0barSideBand.confMaxInvMassLeftSB.value) ||
              (cand2.mAntiLambda() > ConfD0D0barSideBand.confMinInvMassRightSB.value && cand2.mAntiLambda() < ConfD0D0barSideBand.confMaxInvMassRightSB.value)) {
            registry.fill(HIST("hDeltaPhiD0SigD0barBg"), deltaPhi);
          }
        }
        if ((cand1.mLambda() > ConfD0D0barSideBand.confMinInvMassLeftSB.value && cand1.mLambda() < ConfD0D0barSideBand.confMaxInvMassLeftSB.value) ||
            (cand1.mLambda() > ConfD0D0barSideBand.confMinInvMassRightSB.value && cand1.mLambda() < ConfD0D0barSideBand.confMaxInvMassRightSB.value)) {
          // B(D0) x S (D0bar) correlation
          if (cand2.mAntiLambda() > ConfD0D0barSideBand.confSignalRegionMin.value && cand2.mAntiLambda() < ConfD0D0barSideBand.confSignalRegionMax.value) {
            registry.fill(HIST("hDeltaPhiD0BgD0barSig"), deltaPhi);
          }
          // B(D0) x B(D0bar) correlation
          if ((cand2.mAntiLambda() > ConfD0D0barSideBand.confMinInvMassLeftSB.value && cand2.mAntiLambda() < ConfD0D0barSideBand.confMaxInvMassLeftSB.value) ||
              (cand2.mAntiLambda() > ConfD0D0barSideBand.confMinInvMassRightSB.value && cand2.mAntiLambda() < ConfD0D0barSideBand.confMaxInvMassRightSB.value)) {
            registry.fill(HIST("hDeltaPhiBgBg"), deltaPhi);
          }
        }
      } // It is the end of the for loop over D0bar mesons
    } // It is the end of the for loop over all candidates
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackD0, processSideBand, "Enable processing side-band methode", false);

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
      sameEventFemtoCont.setPair<isMC>(track, d0candidate, multCol, ConfBothTracks.confUse3D);
      sameEventAngularCont.setPair<isMC>(track, d0candidate, multCol, ConfBothTracks.confUse3D);
    }
  }

  /// process function for to call doSameEvent with Data
  /// \param col subscribe to the collision table (Data)
  /// \param parts subscribe to the femtoUniverseParticleTable
  void processSameEvent(o2::aod::FdCollision const& col,
                        FemtoFullParticles const& parts)
  {
    fillCollision(col);

    auto thegroupPartsTrack = partsTrack->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto thegroupPartsOnlyD0D0bar = partsOnlyD0D0bar->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto theGroupPartsD0s = partsD0s->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto theGroupPartsD0bars = partsD0bars->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);

    switch (confChooseD0trackCorr) {
      case 0:
        doSameEvent<false>(thegroupPartsTrack, theGroupPartsD0s, parts, col.magField(), col.multNtr());
        break;
      case 1:
        doSameEvent<false>(thegroupPartsTrack, theGroupPartsD0bars, parts, col.magField(), col.multNtr());
        break;
      case 2:
        doSameEvent<false>(thegroupPartsTrack, thegroupPartsOnlyD0D0bar, parts, col.magField(), col.multNtr());
        break;
      default:
        break;
    }
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackD0, processSameEvent, "Enable processing same event", true);

  /// process function for to call doSameEvent with Monte Carlo
  /// \param col subscribe to the collision table (Monte Carlo Reconstructed reconstructed)
  /// \param parts subscribe to joined table FemtoUniverseParticles and FemtoUniverseMCLables to access Monte Carlo truth
  /// \param FemtoUniverseMCParticles subscribe to the Monte Carlo truth table
  void processSameEventMC(o2::aod::FdCollision const& col,
                          soa::Join<o2::aod::FDParticles, o2::aod::FDMCLabels> const& parts,
                          o2::aod::FdMCParticles const&)
  {
    fillCollision(col);

    auto thegroupPartsD0 = partsD0D0barMC->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto thegroupPartsTrack = partsTrackMC->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);

    doSameEvent<true>(thegroupPartsTrack, thegroupPartsD0, parts, col.magField(), col.multNtr());
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackD0, processSameEventMC, "Enable processing same event for Monte Carlo", false);

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
      // // Close Pair Rejection
      if (confIsCPR.value) {
        if (pairCloseRejection.isClosePair(track, d0candidate, parts, magFieldTesla, femto_universe_container::EventType::mixed)) {
          continue;
        }
      }

      mixedEventFemtoCont.setPair<isMC>(track, d0candidate, multCol, ConfBothTracks.confUse3D);
      mixedEventAngularCont.setPair<isMC>(track, d0candidate, multCol, ConfBothTracks.confUse3D);
    }
  }

  /// process function for to call doMixedEvent with Data
  /// @param cols subscribe to the collisions table (Data)
  /// @param parts subscribe to the femtoUniverseParticleTable
  void processMixedEvent(o2::aod::FdCollisions const& cols,
                         FemtoFullParticles const& parts)
  {
    for (auto const& [collision1, collision2] : soa::selfCombinations(colBinning, 5, -1, cols, cols)) {

      const int multiplicityCol = collision1.multNtr();
      mixQaRegistry.fill(HIST("MixingQA/hMECollisionBins"), colBinning.getBin({collision1.posZ(), multiplicityCol}));

      auto groupPartsTrack = partsTrack->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision2.globalIndex(), cache);
      auto groupPartsOnlyD0D0bar = partsOnlyD0D0bar->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);
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
        case 2:
          doMixedEvent<false>(groupPartsTrack, groupPartsOnlyD0D0bar, parts, magFieldTesla1, multiplicityCol);
          break;
        default:
          break;
      }
    }
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackD0, processMixedEvent, "Enable processing mixed events", true);

  /// brief process function for to call doMixedEvent with Monte Carlo
  /// @param cols subscribe to the collisions table (Monte Carlo Reconstructed reconstructed)
  /// @param parts subscribe to joined table FemtoUniverseParticles and FemtoUniverseMCLables to access Monte Carlo truth
  /// @param FemtoUniverseMCParticles subscribe to the Monte Carlo truth table
  void processMixedEventMC(o2::aod::FdCollisions const& cols,
                           soa::Join<o2::aod::FDParticles, o2::aod::FDMCLabels> const& parts,
                           o2::aod::FdMCParticles const&)
  {
    for (auto const& [collision1, collision2] : soa::selfCombinations(colBinning, 5, -1, cols, cols)) {

      const int multiplicityCol = collision1.multNtr();
      mixQaRegistry.fill(HIST("MixingQA/hMECollisionBins"), colBinning.getBin({collision1.posZ(), multiplicityCol}));

      auto groupPartsTrack = partsTrackMC->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision2.globalIndex(), cache);
      auto groupPartsD0 = partsD0D0barMC->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);

      const auto& magFieldTesla1 = collision1.magField();
      const auto& magFieldTesla2 = collision2.magField();

      if (magFieldTesla1 != magFieldTesla2) {
        continue;
      }
      /// \todo before mixing we should check whether both collisions contain a pair of particles!
      // if (partsD0.size() == 0 || kNPart2Evt1 == 0 || kNPart1Evt2 == 0 || partsTrack.size() == 0 ) continue;

      doMixedEvent<true>(groupPartsTrack, groupPartsD0, parts, magFieldTesla1, multiplicityCol);
    }
  }
  PROCESS_SWITCH(FemtoUniversePairTaskTrackD0, processMixedEventMC, "Enable processing mixed events MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<FemtoUniversePairTaskTrackD0>(cfgc),
  };
  return workflow;
}
