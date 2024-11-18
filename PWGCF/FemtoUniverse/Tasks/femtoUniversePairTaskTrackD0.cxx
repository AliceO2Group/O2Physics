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

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::analysis;
using namespace o2::analysis::femtoUniverse;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

namespace
{
static constexpr int nPart = 2;
static constexpr int nCuts = 5;
static const std::vector<std::string> partNames{"D0", "Track"};
static const std::vector<std::string> cutNames{"MaxPt", "PIDthr", "nSigmaTPC", "nSigmaTPCTOF", "MaxP"};
static const float cutsTable[nPart][nCuts]{{4.05f, 1.f, 3.f, 3.f, 100.f}, {4.05f, 1.f, 3.f, 3.f, 100.f}};
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
  if (deltaPhi < 0.) {
    deltaPhi = deltaPhi + o2::constants::math::TwoPI;
  }
  if (deltaPhi > o2::constants::math::TwoPI) {
    deltaPhi = o2::constants::math::TwoPI - deltaPhi;
  }
  return deltaPhi;
}

struct femtoUniversePairTaskTrackD0 {

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
    Configurable<int> ConfPIDTrack{"ConfPIDTrack", 2, "Particle 2 - Read from cutCulator"}; // we also need the possibility to specify whether the bit is true/false ->std>>vector<std::pair<int, int>>
    Configurable<int> ConfTrackSign{"ConfTrackSign", 1, "Track sign"};
    Configurable<bool> ConfIsTrackIdentified{"ConfIsTrackIdentified", true, "Enable PID for the track"};
  } ConfTrack;

  /// Particle 2 --- D0/D0bar meson
  struct : o2::framework::ConfigurableGroup {
    Configurable<int> ConfPDGCodeD0{"ConfPDGCodeD0", 421, "D0 meson - PDG code"};
    Configurable<int> ConfPDGCodeD0bar{"ConfPDGCodeD0bar", -421, "D0bar meson - PDG code"};
    Configurable<float> ConfMinPtD0D0bar{"ConfMinPtD0D0bar", 3.0, "D0/D0bar sel. - min. pT"};
    Configurable<float> ConfMaxPtD0D0bar{"ConfMaxPtD0D0bar", 5.0, "D0/D0bar sel. - max. pT"};
    Configurable<float> ConfMinInvMassD0D0bar{"ConfMinInvMassD0D0bar", 1.65, "D0/D0bar sel. - min. invMass"};
    Configurable<float> ConfMaxInvMassD0D0bar{"ConfMaxInvMassD0D0bar", 2.05, "D0/D0bar sel. - max. invMass"};
  } ConfDmesons;

  struct : o2::framework::ConfigurableGroup {
    Configurable<float> ConfSignalRegionMin{"ConfSignalRegionMin", 1.810, "Min. inv. mass for D0/D0bar in the signal region"};
    Configurable<float> ConfSignalRegionMax{"ConfSignalRegionMax", 1.922, "Max. inv. mass for D0/D0bar in the signal region"};
    Configurable<float> ConfMinInvMassLeftSB{"ConfMinInvMassLeftSB", 1.642, "Min. inv. mass for D0/D0bar in the left sideband region"};
    Configurable<float> ConfMaxInvMassLeftSB{"ConfMaxInvMassLeftSB", 1.754, "Max. inv. mass for D0/D0bar in the left sideband region"};
    Configurable<float> ConfMinInvMassRightSB{"ConfMinInvMassRightSB", 1.978, "Min. inv. mass for D0/D0bar in the right sideband region"};
    Configurable<float> ConfMaxInvMassRightSB{"ConfMaxInvMassRightSB", 2.090, "Max. inv. mass for D0/D0bar in the right sideband region"};
  } ConfD0D0barSideBand;

  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_d0_to_pi_k::vecBinsPt}, "pT bin limits"};
  Configurable<uint> ConfChooseD0trackCorr{"ConfChooseD0trackCorr", 3, "If 0 - only D0s, 1 - only D0bars, 2 - D0/D0bar (one mass hypo.), 3 - all D0/D0bar cand."};
  Configurable<bool> ConfUsePtCutForD0D0bar{"ConfUsePtCutForD0D0bar", false, "Include pT cut for D0/D0bar in same and mixed processes."};
  Configurable<bool> ConfUseMassCutForD0D0bar{"ConfUseMassCutForD0D0bar", false, "Switch to save D0/D0bar within declared inv. mass range"};

  /// Partitions for particle 1
  Partition<FemtoFullParticles> partsTrack = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kTrack)) && (aod::femtouniverseparticle::sign == int8_t(ConfTrack.ConfTrackSign));
  Partition<soa::Join<aod::FDParticles, aod::FDMCLabels>> partsTrackMC = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kTrack));

  /// Partitions for particle 2
  /// Partition with all D0/D0bar mesons (which pass double and one mass hypothesis)
  Partition<FemtoFullParticles> partsAllDmesons = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kD0));
  /// Partition with D0/D0bar candidates, which pass only one mass hypothesis
  Partition<FemtoFullParticles> partsOnlyD0D0bar = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kD0)) && (aod::femtouniverseparticle::mLambda < 0.0f || aod::femtouniverseparticle::mAntiLambda < 0.0f);
  /// Partition with D0 mesons only (one mass hypothesis)
  Partition<FemtoFullParticles> partsD0s = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kD0)) && (aod::femtouniverseparticle::mLambda > 0.0f) && (aod::femtouniverseparticle::mAntiLambda < 0.0f);
  /// Partition with D0bar mesons only (one mass hypothesis)
  Partition<FemtoFullParticles> partsD0bars = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kD0)) && (aod::femtouniverseparticle::mLambda < 0.0f) && (aod::femtouniverseparticle::mAntiLambda > 0.0f);
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
  ConfigurableAxis ConfPtBins{"ConfPtBins", {360, 0., 36.}, "binning pT"};
  ConfigurableAxis ConfInvMassBins{"ConfInvMassBins", {500, 0., 5.0}, "binning inv. mass"};
  ConfigurableAxis ConfInvMassFinerBins{"ConfInvMassFinerBins", {120, 1.5848, 2.1848}, "finer binning of inv. mass"};

  Configurable<bool> ConfIsCPR{"ConfIsCPR", true, "Close Pair Rejection"};
  Configurable<bool> ConfCPRPlotPerRadii{"ConfCPRPlotPerRadii", false, "Plot CPR per radii"};
  Configurable<float> ConfCPRdeltaPhiCutMax{"ConfCPRdeltaPhiCutMax", 0.0, "Delta Phi max cut for Close Pair Rejection"};
  Configurable<float> ConfCPRdeltaPhiCutMin{"ConfCPRdeltaPhiCutMin", 0.0, "Delta Phi min cut for Close Pair Rejection"};
  Configurable<float> ConfCPRdeltaEtaCutMax{"ConfCPRdeltaEtaCutMax", 0.0, "Delta Eta max cut for Close Pair Rejection"};
  Configurable<float> ConfCPRdeltaEtaCutMin{"ConfCPRdeltaEtaCutMin", 0.0, "Delta Eta min cut for Close Pair Rejection"};
  Configurable<float> ConfCPRChosenRadii{"ConfCPRChosenRadii", 0.80, "Delta Eta cut for Close Pair Rejection"};

  FemtoUniverseFemtoContainer<femtoUniverseFemtoContainer::EventType::same, femtoUniverseFemtoContainer::Observable::kstar> sameEventFemtoCont;
  FemtoUniverseFemtoContainer<femtoUniverseFemtoContainer::EventType::mixed, femtoUniverseFemtoContainer::Observable::kstar> mixedEventFemtoCont;
  FemtoUniverseAngularContainer<femtoUniverseAngularContainer::EventType::same, femtoUniverseAngularContainer::Observable::kstar> sameEventAngularCont;
  FemtoUniverseAngularContainer<femtoUniverseAngularContainer::EventType::mixed, femtoUniverseAngularContainer::Observable::kstar> mixedEventAngularCont;
  FemtoUniversePairCleaner<aod::femtouniverseparticle::ParticleType::kTrack, aod::femtouniverseparticle::ParticleType::kD0> pairCleaner;
  FemtoUniverseDetaDphiStar<aod::femtouniverseparticle::ParticleType::kTrack, aod::femtouniverseparticle::ParticleType::kD0> pairCloseRejection;
  FemtoUniverseTrackSelection trackCuts;

  /// Histogram output
  HistogramRegistry qaRegistry{"TrackQA", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry resultRegistry{"Correlations", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry MixQaRegistry{"MixQaRegistry", {}, OutputObjHandlingPolicy::AnalysisObject};

  HistogramRegistry registry{"registry",
                             {{"hInvMassD0", ";#it{M}(K^{-}#pi^{+}) (GeV/#it{c}^{2});counts", {HistType::kTH1F, {ConfInvMassBins}}},
                              {"hInvMassD0bar", ";#it{M}(#pi^{-}K^{+}) (GeV/#it{c}^{2});counts", {HistType::kTH1F, {ConfInvMassBins}}},
                              {"hPtDmesonCand", "2-prong candidates;#it{p}_{T} (GeV/#it{c});counts", {HistType::kTH1F, {ConfPtBins}}},
                              {"hPtD0", "D^{0} cand.;#it{p}_{T} (GeV/#it{c});counts", {HistType::kTH1F, {ConfPtBins}}},
                              {"hPtD0bar", "#bar{D^{0}};#it{p}_{T} (GeV/#it{c});counts", {HistType::kTH1F, {ConfPtBins}}},
                              {"hPtD0D0bar", "#bar{D^{0}};#it{p}_{T} (GeV/#it{c});counts", {HistType::kTH1F, {ConfPtBins}}},
                              {"hPhiDmesonCand", ";#varphi (rad);counts", {HistType::kTH1F, {{80, 0., 2. * o2::constants::math::PI}}}},
                              {"hPhiD0", ";#varphi (rad);counts", {HistType::kTH1F, {{80, 0., 2. * o2::constants::math::PI}}}},
                              {"hPhiD0bar", ";#varphi (rad);counts", {HistType::kTH1F, {{80, 0., 2. * o2::constants::math::PI}}}},
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
  bool IsProtonNSigma(float mom, float nsigmaTPCPr, float nsigmaTOFPr) // previous version from: https://github.com/alisw/AliPhysics/blob/master/PWGCF/FEMTOSCOPY/AliFemtoUser/AliFemtoMJTrackCut.cxx
  {
    if (mom < 0.5) {
      if (TMath::Abs(nsigmaTPCPr) < ConfBothTracks.ConfNsigmaTPCProton) {
        return true;
      } else {
        return false;
      }
    } else if (mom > 0.4) {
      if (TMath::Hypot(nsigmaTOFPr, nsigmaTPCPr) < ConfBothTracks.ConfNsigmaCombinedProton) {
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
    } else if (mom > 1.5) {
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
    qaRegistry.add("D0_pos_daugh/nSigmaTPC", "; #it{p} (GeV/#it{c}); n#sigma_{TPC}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("D0_pos_daugh/nSigmaTOF", "; #it{p} (GeV/#it{c}); n#sigma_{TOF}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("D0_pos_daugh/pt", "; #it{p_T} (GeV/#it{c}); Counts", kTH1F, {{100, 0, 10}});
    qaRegistry.add("D0_pos_daugh/eta", "; #it{eta}; Counts", kTH1F, {{200, -1.5, 1.5}});
    qaRegistry.add("D0_pos_daugh/phi", "; #it{varphi}; Counts", kTH1F, {{200, 0, 2. * M_PI}});
    qaRegistry.add("D0_pos_daugh/hDCAxy", "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", kTH2F, {{100, 0, 10}, {500, -5, 5}});

    qaRegistry.add("D0_neg_daugh/nSigmaTPC", "; #it{p} (GeV/#it{c}); n#sigma_{TPC}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("D0_neg_daugh/nSigmaTOF", "; #it{p} (GeV/#it{c}); n#sigma_{TOF}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("D0_neg_daugh/pt", "; #it{p_T} (GeV/#it{c}); Counts", kTH1F, {{100, 0, 10}});
    qaRegistry.add("D0_neg_daugh/eta", "; #it{eta}; Counts", kTH1F, {{200, -1.5, 1.5}});
    qaRegistry.add("D0_neg_daugh/phi", "; #it{varphi}; Counts", kTH1F, {{200, 0, 2. * M_PI}});
    qaRegistry.add("D0_neg_daugh/hDCAxy", "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", kTH2F, {{100, 0, 10}, {500, -5, 5}});

    qaRegistry.add("D0bar_pos_daugh/nSigmaTPC", "; #it{p} (GeV/#it{c}); n#sigma_{TPC}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("D0bar_pos_daugh/nSigmaTOF", "; #it{p} (GeV/#it{c}); n#sigma_{TOF}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("D0bar_pos_daugh/pt", "; #it{p_T} (GeV/#it{c}); Counts", kTH1F, {{100, 0, 10}});
    qaRegistry.add("D0bar_pos_daugh/eta", "; #it{eta}; Counts", kTH1F, {{200, -1.5, 1.5}});
    qaRegistry.add("D0bar_pos_daugh/phi", "; #it{varphi}; Counts", kTH1F, {{200, 0, 2. * M_PI}});
    qaRegistry.add("D0bar_pos_daugh/hDCAxy", "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", kTH2F, {{100, 0, 10}, {500, -5, 5}});

    qaRegistry.add("D0bar_neg_daugh/nSigmaTPC", "; #it{p} (GeV/#it{c}); n#sigma_{TPC}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("D0bar_neg_daugh/nSigmaTOF", "; #it{p} (GeV/#it{c}); n#sigma_{TOF}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("D0bar_neg_daugh/pt", "; #it{p_T} (GeV/#it{c}); Counts", kTH1F, {{100, 0, 10}});
    qaRegistry.add("D0bar_neg_daugh/eta", "; #it{eta}; Counts", kTH1F, {{200, -1.5, 1.5}});
    qaRegistry.add("D0bar_neg_daugh/phi", "; #it{varphi}; Counts", kTH1F, {{200, 0, 2. * M_PI}});
    qaRegistry.add("D0bar_neg_daugh/hDCAxy", "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", kTH2F, {{100, 0, 10}, {500, -5, 5}});

    qaRegistry.add("Hadron/nSigmaTPCPr", "; #it{p} (GeV/#it{c}); n#sigma_{TPCPr}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("Hadron/nSigmaTOFPr", "; #it{p} (GeV/#it{c}); n#sigma_{TOFPr}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("Hadron/nSigmaTPCPi", "; #it{p} (GeV/#it{c}); n#sigma_{TPCPi}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("Hadron/nSigmaTOFPi", "; #it{p} (GeV/#it{c}); n#sigma_{TOFPi}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("Hadron/nSigmaTPCKa", "; #it{p} (GeV/#it{c}); n#sigma_{TPCKa}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});
    qaRegistry.add("Hadron/nSigmaTOFKa", "; #it{p} (GeV/#it{c}); n#sigma_{TOFKa}", kTH2F, {{100, 0, 10}, {200, -4.975, 5.025}});

    trackHistoPartD0D0bar.init(&qaRegistry, ConfTempFitVarpTBins, ConfTempFitVarInvMassBins, ConfBothTracks.ConfIsMC, ConfDmesons.ConfPDGCodeD0);
    if (!ConfTrack.ConfIsSame) {
      trackHistoPartTrack.init(&qaRegistry, ConfTempFitVarpTBins, ConfTempFitVarBins, ConfBothTracks.ConfIsMC, ConfTrack.ConfPDGCodeTrack);
    }

    MixQaRegistry.add("MixingQA/hSECollisionBins", ";bin;Entries", kTH1F, {{120, -0.5, 119.5}});
    MixQaRegistry.add("MixingQA/hMECollisionBins", ";bin;Entries", kTH1F, {{120, -0.5, 119.5}});

    sameEventFemtoCont.init(&resultRegistry, ConfkstarBins, ConfMultBins, ConfkTBins, ConfmTBins, ConfmultBins3D, ConfmTBins3D, ConfBothTracks.ConfIsMC, ConfBothTracks.ConfUse3D);
    mixedEventFemtoCont.init(&resultRegistry, ConfkstarBins, ConfMultBins, ConfkTBins, ConfmTBins, ConfmultBins3D, ConfmTBins3D, ConfBothTracks.ConfIsMC, ConfBothTracks.ConfUse3D);
    sameEventAngularCont.init(&resultRegistry, ConfkstarBins, ConfMultBins, ConfkTBins, ConfmTBins, ConfmultBins3D, ConfmTBins3D, ConfBothTracks.ConfEtaBins, ConfBothTracks.ConfPhiBins, ConfBothTracks.ConfIsMC, ConfBothTracks.ConfUse3D);
    mixedEventAngularCont.init(&resultRegistry, ConfkstarBins, ConfMultBins, ConfkTBins, ConfmTBins, ConfmultBins3D, ConfmTBins3D, ConfBothTracks.ConfEtaBins, ConfBothTracks.ConfPhiBins, ConfBothTracks.ConfIsMC, ConfBothTracks.ConfUse3D);

    sameEventFemtoCont.setPDGCodes(ConfDmesons.ConfPDGCodeD0, ConfTrack.ConfPDGCodeTrack);
    mixedEventFemtoCont.setPDGCodes(ConfDmesons.ConfPDGCodeD0, ConfTrack.ConfPDGCodeTrack);
    sameEventAngularCont.setPDGCodes(ConfDmesons.ConfPDGCodeD0, ConfTrack.ConfPDGCodeTrack);
    mixedEventAngularCont.setPDGCodes(ConfDmesons.ConfPDGCodeD0, ConfTrack.ConfPDGCodeTrack);

    pairCleaner.init(&qaRegistry);
    if (ConfIsCPR.value) {
      pairCloseRejection.init(&resultRegistry, &qaRegistry, ConfCPRdeltaPhiCutMin.value, ConfCPRdeltaPhiCutMax.value, ConfCPRdeltaEtaCutMin.value, ConfCPRdeltaEtaCutMax.value, ConfCPRChosenRadii.value, ConfCPRPlotPerRadii.value);
    }

    vPIDTrack = ConfTrack.ConfPIDTrack.value;
    kNsigma = ConfBothTracks.ConfTrkPIDnSigmaMax.value;

    // D0/D0bar histograms
    auto vbins = (std::vector<double>)binsPt;
    registry.add("hMassVsPt", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {ConfInvMassBins, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassVsPtFinerBinning", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {ConfInvMassFinerBins, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hInvMassVsPtOnlyD0D0bar", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {ConfInvMassBins, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDeltaPhiSigSig", "SxS correlation;#Delta#varphi (rad);counts", {HistType::kTH1F, {{10, 0.0, o2::constants::math::PI}}});
    registry.add("hDeltaPhiD0BgD0barSig", "B(D0)x S(D0bar) correlation;#Delta#varphi (rad);counts", {HistType::kTH1F, {{10, 0.0, o2::constants::math::PI}}});
    registry.add("hDeltaPhiD0SigD0barBg", "S(D0)x B(D0bar) correlation;#Delta#varphi (rad);counts", {HistType::kTH1F, {{10, 0.0, o2::constants::math::PI}}});
    registry.add("hDeltaPhiBgBg", "BxB correlation;#Delta#varphi (rad);counts", {HistType::kTH1F, {{10, 0.0, o2::constants::math::PI}}});
    registry.add("hPtCand1VsPtCand2", "2-prong candidates;#it{p}_{T} (GeV/#it{c});#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{vbins, "#it{p}_{T} (GeV/#it{c})"}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDeltaEtaDeltaPhi", "2-prong candidates;#Delta #eta;#Delta #varphi (rad)", {HistType::kTH2F, {{29, -2., 2.}, {29, 0.0, o2::constants::math::PI}}});
  }

  template <typename CollisionType>
  void fillCollision(CollisionType col)
  {
    MixQaRegistry.fill(HIST("MixingQA/hSECollisionBins"), colBinning.getBin({col.posZ(), col.multNtr()}));
    eventHisto.fillQA(col);
  }

  void processD0mesons(o2::aod::FDCollision& col, FemtoFullParticles&)
  {
    auto groupPartsOnlyD0D0bar = partsOnlyD0D0bar->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto groupPartsAllDmesons = partsAllDmesons->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto groupPartsD0D0barChildren = partsDmesonsChildren->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);

    // loop over all D mesons
    for (auto const& dmeson : groupPartsAllDmesons) {

      if (dmeson.mLambda() > 0.0f) {
        registry.fill(HIST("hMassVsPt"), dmeson.mLambda(), dmeson.pt());
        registry.fill(HIST("hMassVsPtFinerBinning"), dmeson.mLambda(), dmeson.pt());
      }

      if (dmeson.mAntiLambda() > 0.0f) {
        registry.fill(HIST("hMassVsPt"), dmeson.mAntiLambda(), dmeson.pt());
        registry.fill(HIST("hMassVsPtFinerBinning"), dmeson.mAntiLambda(), dmeson.pt());
      }

      registry.fill(HIST("hPtDmesonCand"), dmeson.pt());
      registry.fill(HIST("hPhiDmesonCand"), dmeson.phi());
      registry.fill(HIST("hEtaDmesonCand"), dmeson.eta());
    }

    // loop over D0/D0bar mesons (ONLY)
    for (auto const& d0d0bar : groupPartsOnlyD0D0bar) {

      registry.fill(HIST("hPtD0D0bar"), d0d0bar.pt());

      if (d0d0bar.mLambda() > 0.0f && d0d0bar.mAntiLambda() < 0.0f) {
        registry.fill(HIST("hInvMassVsPtOnlyD0D0bar"), d0d0bar.mLambda(), d0d0bar.pt());
        if (d0d0bar.mLambda() > ConfDmesons.ConfMinInvMassD0D0bar && d0d0bar.mLambda() < ConfDmesons.ConfMaxInvMassD0D0bar) {
          registry.fill(HIST("hInvMassD0"), d0d0bar.mLambda());
        }
        registry.fill(HIST("hPtD0"), d0d0bar.pt());
        registry.fill(HIST("hPhiD0"), d0d0bar.phi());
        registry.fill(HIST("hEtaD0"), d0d0bar.eta());
      }
      if (d0d0bar.mLambda() < 0.0f && d0d0bar.mAntiLambda() > 0.0f) {
        registry.fill(HIST("hInvMassVsPtOnlyD0D0bar"), d0d0bar.mAntiLambda(), d0d0bar.pt());
        if (d0d0bar.mAntiLambda() > ConfDmesons.ConfMinInvMassD0D0bar && d0d0bar.mAntiLambda() < ConfDmesons.ConfMaxInvMassD0D0bar) {
          registry.fill(HIST("hInvMassD0bar"), d0d0bar.mAntiLambda());
        }
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
  PROCESS_SWITCH(femtoUniversePairTaskTrackD0, processD0mesons, "Enable processing D0 mesons", true);

  // D0-D0bar pair correlations (side-band methode)
  void processSideBand(o2::aod::FDCollision& col, FemtoFullParticles&)
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
        if (cand1.mLambda() > ConfD0D0barSideBand.ConfSignalRegionMin.value && cand1.mLambda() < ConfD0D0barSideBand.ConfSignalRegionMax.value) {
          // S(D0) x S(D0bar) correlation
          if (cand2.mAntiLambda() > ConfD0D0barSideBand.ConfSignalRegionMin.value && cand2.mAntiLambda() < ConfD0D0barSideBand.ConfSignalRegionMax.value) {
            registry.fill(HIST("hDeltaPhiSigSig"), deltaPhi);
          }
          // S(D0) x B(D0bar) correlation
          if ((cand2.mAntiLambda() > ConfD0D0barSideBand.ConfMinInvMassLeftSB.value && cand2.mAntiLambda() < ConfD0D0barSideBand.ConfMaxInvMassLeftSB.value) ||
              (cand2.mAntiLambda() > ConfD0D0barSideBand.ConfMinInvMassRightSB.value && cand2.mAntiLambda() < ConfD0D0barSideBand.ConfMaxInvMassRightSB.value)) {
            registry.fill(HIST("hDeltaPhiD0SigD0barBg"), deltaPhi);
          }
        }
        if ((cand1.mLambda() > ConfD0D0barSideBand.ConfMinInvMassLeftSB.value && cand1.mLambda() < ConfD0D0barSideBand.ConfMaxInvMassLeftSB.value) ||
            (cand1.mLambda() > ConfD0D0barSideBand.ConfMinInvMassRightSB.value && cand1.mLambda() < ConfD0D0barSideBand.ConfMaxInvMassRightSB.value)) {
          // B(D0) x S (D0bar) correlation
          if (cand2.mAntiLambda() > ConfD0D0barSideBand.ConfSignalRegionMin.value && cand2.mAntiLambda() < ConfD0D0barSideBand.ConfSignalRegionMax.value) {
            registry.fill(HIST("hDeltaPhiD0BgD0barSig"), deltaPhi);
          }
          // B(D0) x B(D0bar) correlation
          if ((cand2.mAntiLambda() > ConfD0D0barSideBand.ConfMinInvMassLeftSB.value && cand2.mAntiLambda() < ConfD0D0barSideBand.ConfMaxInvMassLeftSB.value) ||
              (cand2.mAntiLambda() > ConfD0D0barSideBand.ConfMinInvMassRightSB.value && cand2.mAntiLambda() < ConfD0D0barSideBand.ConfMaxInvMassRightSB.value)) {
            registry.fill(HIST("hDeltaPhiBgBg"), deltaPhi);
          }
        }
      } // It is the end of the for loop over D0bar mesons
    } // It is the end of the for loop over all candidates
  }
  PROCESS_SWITCH(femtoUniversePairTaskTrackD0, processSideBand, "Enable processing side-band methode", false);

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
    for (auto& d0candidate : groupPartsD0) {
      trackHistoPartD0D0bar.fillQA<isMC, false>(d0candidate);
    }

    float tpcNSigmaPr, tofNSigmaPr, tpcNSigmaPi, tofNSigmaPi, tpcNSigmaKa, tofNSigmaKa;

    if (!ConfTrack.ConfIsSame) {
      for (auto& track : groupPartsTrack) {
        if (ConfTrack.ConfIsTrackIdentified) {
          if (!IsParticleNSigma(track.p(), trackCuts.getNsigmaTPC(track, o2::track::PID::Proton), trackCuts.getNsigmaTOF(track, o2::track::PID::Proton), trackCuts.getNsigmaTPC(track, o2::track::PID::Pion), trackCuts.getNsigmaTOF(track, o2::track::PID::Pion), trackCuts.getNsigmaTPC(track, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(track, o2::track::PID::Kaon))) {
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
    for (auto& [track, d0candidate] : combinations(CombinationsFullIndexPolicy(groupPartsTrack, groupPartsD0))) {
      if (ConfTrack.ConfIsTrackIdentified) {
        if (!IsParticleNSigma(track.p(), trackCuts.getNsigmaTPC(track, o2::track::PID::Proton), trackCuts.getNsigmaTOF(track, o2::track::PID::Proton), trackCuts.getNsigmaTPC(track, o2::track::PID::Pion), trackCuts.getNsigmaTOF(track, o2::track::PID::Pion), trackCuts.getNsigmaTPC(track, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(track, o2::track::PID::Kaon))) {
          continue;
        }
      }
      // // Set pT cut for D0/D0bar candidates
      if (ConfUsePtCutForD0D0bar) {
        if (d0candidate.pt() < ConfDmesons.ConfMinPtD0D0bar && d0candidate.pt() > ConfDmesons.ConfMaxPtD0D0bar) {
          continue;
        }
      }
      // // Set inv. mass cut for D0/D0bar candidates
      if (ConfUseMassCutForD0D0bar) {
        if ((d0candidate.mLambda() < ConfD0D0barSideBand.ConfSignalRegionMin && d0candidate.mLambda() > ConfD0D0barSideBand.ConfSignalRegionMax) || (d0candidate.mAntiLambda() < ConfD0D0barSideBand.ConfSignalRegionMin && d0candidate.mAntiLambda() > ConfD0D0barSideBand.ConfSignalRegionMax)) {
          continue;
        }
      }
      // // Close Pair Rejection
      if (ConfIsCPR.value) {
        if (pairCloseRejection.isClosePair(track, d0candidate, parts, magFieldTesla, femtoUniverseContainer::EventType::same)) {
          continue;
        }
      }

      // Track Cleaning
      if (!pairCleaner.isCleanPair(track, d0candidate, parts)) {
        continue;
      }
      sameEventFemtoCont.setPair<isMC>(track, d0candidate, multCol, ConfBothTracks.ConfUse3D);
      sameEventAngularCont.setPair<isMC>(track, d0candidate, multCol, ConfBothTracks.ConfUse3D);
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
    auto thegroupPartsAllD0D0bar = partsAllDmesons->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto thegroupPartsOnlyD0D0bar = partsOnlyD0D0bar->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto theGroupPartsD0s = partsD0s->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto theGroupPartsD0bars = partsD0bars->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);

    switch (ConfChooseD0trackCorr) {
      case 0:
        doSameEvent<false>(thegroupPartsTrack, theGroupPartsD0s, parts, col.magField(), col.multNtr());
        break;
      case 1:
        doSameEvent<false>(thegroupPartsTrack, theGroupPartsD0bars, parts, col.magField(), col.multNtr());
        break;
      case 2:
        doSameEvent<false>(thegroupPartsTrack, thegroupPartsOnlyD0D0bar, parts, col.magField(), col.multNtr());
        break;
      case 3:
        doSameEvent<false>(thegroupPartsTrack, thegroupPartsAllD0D0bar, parts, col.magField(), col.multNtr());
      default:
        break;
    }
  }
  PROCESS_SWITCH(femtoUniversePairTaskTrackD0, processSameEvent, "Enable processing same event", true);

  /// process function for to call doSameEvent with Monte Carlo
  /// \param col subscribe to the collision table (Monte Carlo Reconstructed reconstructed)
  /// \param parts subscribe to joined table FemtoUniverseParticles and FemtoUniverseMCLables to access Monte Carlo truth
  /// \param FemtoUniverseMCParticles subscribe to the Monte Carlo truth table
  void processSameEventMC(o2::aod::FDCollision& col,
                          soa::Join<o2::aod::FDParticles, o2::aod::FDMCLabels>& parts,
                          o2::aod::FDMCParticles&)
  {
    fillCollision(col);

    auto thegroupPartsD0 = partsD0D0barMC->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    auto thegroupPartsTrack = partsTrackMC->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);

    doSameEvent<true>(thegroupPartsTrack, thegroupPartsD0, parts, col.magField(), col.multNtr());
  }
  PROCESS_SWITCH(femtoUniversePairTaskTrackD0, processSameEventMC, "Enable processing same event for Monte Carlo", false);

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

    for (auto& [track, d0candidate] : combinations(CombinationsFullIndexPolicy(groupPartsTrack, groupPartsD0))) {
      if (ConfTrack.ConfIsTrackIdentified) {
        if (!IsParticleNSigma(track.p(), trackCuts.getNsigmaTPC(track, o2::track::PID::Proton), trackCuts.getNsigmaTOF(track, o2::track::PID::Proton), trackCuts.getNsigmaTPC(track, o2::track::PID::Pion), trackCuts.getNsigmaTOF(track, o2::track::PID::Pion), trackCuts.getNsigmaTPC(track, o2::track::PID::Kaon), trackCuts.getNsigmaTOF(track, o2::track::PID::Kaon))) {
          continue;
        }
      }
      // // Set pT cut for D0/D0bar candidates
      if (ConfUsePtCutForD0D0bar) {
        if (d0candidate.pt() < ConfDmesons.ConfMinPtD0D0bar && d0candidate.pt() > ConfDmesons.ConfMaxPtD0D0bar) {
          continue;
        }
      }
      // // Set inv. mass cut for D0/D0bar candidates
      if (ConfUseMassCutForD0D0bar) {
        if ((d0candidate.mLambda() < ConfD0D0barSideBand.ConfSignalRegionMin && d0candidate.mLambda() > ConfD0D0barSideBand.ConfSignalRegionMax) || (d0candidate.mAntiLambda() < ConfD0D0barSideBand.ConfSignalRegionMin && d0candidate.mAntiLambda() > ConfD0D0barSideBand.ConfSignalRegionMax)) {
          continue;
        }
      }
      // // Close Pair Rejection
      if (ConfIsCPR.value) {
        if (pairCloseRejection.isClosePair(track, d0candidate, parts, magFieldTesla, femtoUniverseContainer::EventType::mixed)) {
          continue;
        }
      }

      mixedEventFemtoCont.setPair<isMC>(track, d0candidate, multCol, ConfBothTracks.ConfUse3D);
      mixedEventAngularCont.setPair<isMC>(track, d0candidate, multCol, ConfBothTracks.ConfUse3D);
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
      auto groupPartsAllD0D0bar = partsAllDmesons->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);
      auto groupPartsOnlyD0D0bar = partsOnlyD0D0bar->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);
      auto theGroupPartsD0s = partsD0s->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);
      auto theGroupPartsD0bars = partsD0bars->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);

      const auto& magFieldTesla1 = collision1.magField();
      const auto& magFieldTesla2 = collision2.magField();

      if (magFieldTesla1 != magFieldTesla2) {
        continue;
      }
      /// \todo before mixing we should check whether both collisions contain a pair of particles!
      // if (partsD0.size() == 0 || nPart2Evt1 == 0 || nPart1Evt2 == 0 || partsTrack.size() == 0 ) continue;

      switch (ConfChooseD0trackCorr) {
        case 0:
          doMixedEvent<false>(groupPartsTrack, theGroupPartsD0s, parts, magFieldTesla1, multiplicityCol);
          break;
        case 1:
          doMixedEvent<false>(groupPartsTrack, theGroupPartsD0bars, parts, magFieldTesla1, multiplicityCol);
          break;
        case 2:
          doMixedEvent<false>(groupPartsTrack, groupPartsOnlyD0D0bar, parts, magFieldTesla1, multiplicityCol);
          break;
        case 3:
          doMixedEvent<false>(groupPartsTrack, groupPartsAllD0D0bar, parts, magFieldTesla1, multiplicityCol);
        default:
          break;
      }
    }
  }
  PROCESS_SWITCH(femtoUniversePairTaskTrackD0, processMixedEvent, "Enable processing mixed events", true);

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
      auto groupPartsD0 = partsD0D0barMC->sliceByCached(aod::femtouniverseparticle::fdCollisionId, collision1.globalIndex(), cache);

      const auto& magFieldTesla1 = collision1.magField();
      const auto& magFieldTesla2 = collision2.magField();

      if (magFieldTesla1 != magFieldTesla2) {
        continue;
      }
      /// \todo before mixing we should check whether both collisions contain a pair of particles!
      // if (partsD0.size() == 0 || nPart2Evt1 == 0 || nPart1Evt2 == 0 || partsTrack.size() == 0 ) continue;

      doMixedEvent<true>(groupPartsTrack, groupPartsD0, parts, magFieldTesla1, multiplicityCol);
    }
  }
  PROCESS_SWITCH(femtoUniversePairTaskTrackD0, processMixedEventMC, "Enable processing mixed events MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<femtoUniversePairTaskTrackD0>(cfgc),
  };
  return workflow;
}
