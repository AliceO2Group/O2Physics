// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file lambda1520_spherocity_analysis.cxx
/// \brief Invariant Mass Reconstruction of Lambda(1520) Resonance
/// \author Yash Patley <yash.patley@cern.ch>

#include <TLorentzVector.h>
#include <TRandom.h>
#include <vector>

#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/runDataProcessing.h"
#include "PWGLF/DataModel/LFResonanceTables.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Common/Core/RecoDecay.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

struct lambdaAnalysis {

  SliceCache cache;
  Preslice<aod::ResoTracks> perRCol = aod::resodaughter::resoCollisionId;
  Preslice<aod::Tracks> perCollision = aod::track::collisionId;

  // Configurables.
  Configurable<int> nBinsPt{"nBinsPt", 500, "N bins in pT histogram"};
  Configurable<int> nBinsInvM{"nBinsInvM", 300, "N bins in InvMass histogram"};
  Configurable<int> nBinsSp{"nBinsSp", 120, "N bins in spherocity histogram"};
  Configurable<bool> doRotate{"doRotate", true, "rotated inv mass spectra"};

  // Tracks
  Configurable<float> cPtMin{"cPtMin", 0.2, "Minimum Track pT"};
  Configurable<float> cEtaCut{"cEtaCut", 0.8, "Pseudorapidity cut"};
  Configurable<float> cDcaz{"cDcazMin", 1., "Minimum DCAz"};
  Configurable<float> cDcaxy{"cDcaxyMin", 0.1, "Minimum DCAxy"};
  Configurable<bool> cPrimaryTrack{"cPrimaryTrack", true, "Primary track selection"};                    // kGoldenChi2 | kDCAxy | kDCAz
  Configurable<bool> cGlobalWoDCATrack{"cGlobalWoDCATrack", true, "Global track selection without DCA"}; // kQualityTracks (kTrackType | kTPCNCls | kTPCCrossedRows | kTPCCrossedRowsOverNCls | kTPCChi2NDF | kTPCRefit | kITSNCls | kITSChi2NDF | kITSRefit | kITSHits) | kInAcceptanceTracks (kPtRange | kEtaRange)
  Configurable<bool> cPVContributor{"cPVContributor", true, "PV contributor track selection"};           // PV Contriuibutor

  // Kinematics cuts
  Configurable<bool> cKinCuts{"cKinCuts", true, "Kinematic Cuts for p-K pair opening angle"};
  Configurable<std::vector<float>> cKinCutsPt{"cKinCutsPt", {0.0, 0.4, 0.6, 0.8, 1.0, 1.4, 1.8, 2.2, 2.6, 3.0, 4.0, 5.0, 6.0, 1e10}, "p_{T} of L* for kinematic cuts"};
  Configurable<std::vector<float>> cKinLowerCutsAlpha{"cKinLowerCutsAlpha", {1.5, 1.0, 0.5, 0.3, 0.2, 0.15, 0.1, 0.08, 0.07, 0.06, 0.04, 0.02}, "Lower cut on Opening angle of p-K of L*"};
  Configurable<std::vector<float>> cKinUpperCutsAlpha{"cKinUpperCutsAlpha", {3.0, 2.0, 1.5, 1.4, 1.0, 0.8, 0.6, 0.5, 0.45, 0.35, 0.3, 0.25, 0.2}, "Upper cut on Opening angle of p-K of L*"};

  // PID Selections
  Configurable<bool> cUseOnlyTOFTrackPr{"cUseOnlyTOFTrackPr", false, "Use only TOF track for PID selection"}; // Use only TOF track for Proton PID selection
  Configurable<bool> cUseOnlyTOFTrackKa{"cUseOnlyTOFTrackKa", false, "Use only TOF track for PID selection"}; // Use only TOF track for Kaon PID selection
  Configurable<bool> cUseTpcOnly{"cUseTpcOnly", false, "Use TPC Only selection"};                             // TPC And TOF tracks
  Configurable<float> cRejNsigmaTpc{"cRejNsigmaTpc", 3.0, "Reject tracks to improve purity of TPC PID"};      // Reject missidentified particles when tpc bands merge
  Configurable<float> cRejNsigmaTof{"cRejNsigmaTof", 3.0, "Reject tracks to improve purity of TOF PID"};      // Reject missidentified particles when tpc bands merge
  // Proton
  Configurable<double> cMaxTPCnSigmaProton{"cMaxTPCnSigmaProton", 4.0, "TPC nSigma cut for Proton"};              // TPC
  Configurable<double> cMaxTOFnSigmaProton{"cMaxTOFnSigmaProton", 3.0, "TOF nSigma cut for Proton"};              // TOF
  Configurable<double> nsigmaCutCombinedProton{"nsigmaCutCombinedProton", 3.0, "Combined nSigma cut for Proton"}; // Combined
  Configurable<std::vector<float>> protonTPCPIDp{"protonTPCPIDp", {0.15, 0.5, 0.62, 0.72, 0.90, 1.0}, "p dependent TPC cuts protons"};
  Configurable<std::vector<float>> protonTPCPIDcut{"protonTPCPIDcut", {5, 3.5, 2.8, 2.5, 2}, "TPC nsigma cuts protons"};
  // Kaon
  Configurable<double> cMaxTPCnSigmaKaon{"cMaxTPCnSigmaKaon", 4.0, "TPC nSigma cut for Kaon"};              // TPC
  Configurable<double> cMaxTOFnSigmaKaon{"cMaxTOFnSigmaKaon", 3.0, "TOF nSigma cut for Kaon"};              // TOF
  Configurable<double> nsigmaCutCombinedKaon{"nsigmaCutCombinedKaon", 3.0, "Combined nSigma cut for Kaon"}; // Combined
  Configurable<std::vector<float>> kaonTPCPIDp{"kaonTPCPIDp", {0.15, 0.3, 0.35, 0.40, 0.54}, "p dependent TPC cuts kaons"};
  Configurable<std::vector<float>> kaonTPCPIDcut{"kaonTPCPIDcut", {5., 4., 3., 2.}, "TPC nsigma cuts kaons"};
  // Event Mixing.
  Configurable<bool> cMixSph{"cMixSph", true, "Include Sph Bins to be mixed"};
  Configurable<int> cNumMixEv{"cNumMixEv", 20, "Number of Events to be mixed"};
  ConfigurableAxis cMixVtxBins{"cMixVtxBins", {VARIABLE_WIDTH, -10.0f, -9.f, -8.f, -7.f, -6.f, -5.f, -4.f, -3.f, -2.f, -1.f, 0.f, 1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis cMixMultBins{"cMixMultBins", {VARIABLE_WIDTH, 0.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.0f}, "Mixing bins - multiplicity"};
  ConfigurableAxis cMixSphBins{"cMixSphBins", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f}, "Mixing bins - spherocity"};

  // Histogram Registry.
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {

    // Define Axis.
    const AxisSpec axisPosZ(240, -12., 12., "V_{z} (cm)");
    const AxisSpec axisSp(nBinsSp, 0., 1., "S_{0}");
    const AxisSpec axisCent(105, 0, 105, "FT0M (%)");
    const AxisSpec axisP_pid(600, 0., 6., "p (GeV/c)");
    const AxisSpec axisPt_pid(600, 0., 6., "p_{T} (GeV/c)");
    const AxisSpec axisPt(nBinsPt, 0., 10., "p_{T} (GeV/c)");
    const AxisSpec axisAlpha(200, 0 - 0.15, 2 * TMath::Pi() + 0.15, "#alpha");
    const AxisSpec axisEta(40, -1, 1, "#eta");
    const AxisSpec axisPhi(128, -0.05, 6.35, "#phi (rad)");
    const AxisSpec axisDCAz(500, -0.5, 0.5, {"DCA_{z} (cm)"});
    const AxisSpec axisDCAxy(240, -0.12, 0.12, {"DCA_{xy} (cm)"});
    const AxisSpec axisTPCNCls(200, 0, 200, {"TPCNCls"});
    const AxisSpec axisTPCNsigma(401, -10.025, 10.025, {"n#sigma^{TPC}"});
    const AxisSpec axisTOFNsigma(401, -10.025, 10.025, {"n#sigma^{TOF}"});
    const AxisSpec axisdEdx(380, 10, 200, {"#frac{dE}{dx}"});
    const AxisSpec axisInvM(nBinsInvM, 1.425, 2.025, {"M_{inv} (GeV/#it{c}^{2})"});

    // Create Histograms.
    // Event
    histos.add("Event/h1d_ft0m_mult_percentile", "FT0M (%)", kTH1F, {axisCent});
    histos.add("Event/h1d_spherocity", "Event Spherocity", kTH1F, {axisSp});
    histos.add("Event/h2d_sph_vs_multpercentile", "Spherocity vs FT0M(%)", kTH2F, {axisCent, axisSp});

    // QA Before
    histos.add("QAbefore/Proton/h2d_pr_nsigma_tpc_p", "n#sigma^{TPC} Protons", kTH2F, {axisP_pid, axisTPCNsigma});
    histos.add("QAbefore/Proton/h2d_pr_nsigma_tof_p", "n#sigma^{TOF} Protons", kTH2F, {axisP_pid, axisTOFNsigma});
    histos.add("QAbefore/Proton/h2d_pr_nsigma_tof_vs_tpc", "n#sigma^{TPC} vs n#sigma^{TOF} Protons", kTH2F, {axisTPCNsigma, axisTOFNsigma});
    histos.add("QAbefore/Kaon/h2d_ka_nsigma_tpc_p", "n#sigma^{TPC} Kaons", kTH2F, {axisP_pid, axisTPCNsigma});
    histos.add("QAbefore/Kaon/h2d_ka_nsigma_tof_p", "n#sigma^{TOF} Kaons", kTH2F, {axisP_pid, axisTOFNsigma});
    histos.add("QAbefore/Kaon/h2d_ka_nsigma_tof_vs_tpc", "n#sigma^{TPC} vs n#sigma^{TOF} Kaons", kTH2F, {axisTPCNsigma, axisTOFNsigma});

    // QA After
    histos.add("QAafter/Proton/h2d_pr_dca_z", "Protons", kTH2F, {axisPt_pid, axisDCAz});
    histos.add("QAafter/Proton/h2d_pr_dca_xy", "Protons", kTH2F, {axisPt_pid, axisDCAxy});
    histos.add("QAafter/Proton/h2d_pr_dEdx_p", "TPC Signal Protons", kTH2F, {axisP_pid, axisdEdx});
    histos.add("QAafter/Proton/h2d_pr_nsigma_tpc_pt", "Protons", kTH2F, {axisPt_pid, axisTPCNsigma});
    histos.add("QAafter/Proton/h2d_pr_nsigma_tpc_p", "Protons", kTH2F, {axisP_pid, axisTPCNsigma});
    histos.add("QAafter/Proton/h2d_pr_nsigma_tof_pt", "Protons", kTH2F, {axisPt_pid, axisTOFNsigma});
    histos.add("QAafter/Proton/h2d_pr_nsigma_tof_p", "Protons", kTH2F, {axisP_pid, axisTOFNsigma});
    histos.add("QAafter/Proton/h2d_pr_nsigma_tof_vs_tpc", "n#sigma(TOF) vs n#sigma(TPC) Protons", kTH2F, {axisTPCNsigma, axisTOFNsigma});
    histos.add("QAafter/Kaon/h2d_ka_dca_z", "Kaons", kTH2F, {axisPt_pid, axisDCAz});
    histos.add("QAafter/Kaon/h2d_ka_dca_xy", "Kaons", kTH2F, {axisPt_pid, axisDCAxy});
    histos.add("QAafter/Kaon/h2d_ka_dEdx_p", "TPC Signal Kaon", kTH2F, {axisP_pid, axisdEdx});
    histos.add("QAafter/Kaon/h2d_ka_nsigma_tpc_pt", "Kaons", kTH2F, {axisPt_pid, axisTPCNsigma});
    histos.add("QAafter/Kaon/h2d_ka_nsigma_tpc_p", "Kaons", kTH2F, {axisP_pid, axisTPCNsigma});
    histos.add("QAafter/Kaon/h2d_ka_nsigma_tof_pt", "Kaons", kTH2F, {axisPt_pid, axisTOFNsigma});
    histos.add("QAafter/Kaon/h2d_ka_nsigma_tof_p", "Kaons", kTH2F, {axisP_pid, axisTOFNsigma});
    histos.add("QAafter/Kaon/h2d_ka_nsigma_tof_vs_tpc", "n#sigma(TOF) vs n#sigma(TPC) Kaons", kTH2F, {axisTPCNsigma, axisTOFNsigma});

    // QA checks for protons and kaons
    histos.add("QAChecks/h1d_pr_pt", "p_{T}-spectra Protons", kTH1F, {axisPt_pid});
    histos.add("QAChecks/h1d_ka_pt", "p_{T}-spectra Kaons", kTH1F, {axisPt_pid});
    histos.add("QAChecks/h2d_lstar_alpha_vs_pt", "#alpha_{oa} vs p_{T} Before Cuts", kTH2F, {axisPt, axisAlpha});
    histos.add("QAChecks/h2d_sel_kincuts_lstar_alpha_vs_pt", "#alpha_{oa} vs p_{T} After Cuts", kTH2F, {axisPt, axisAlpha});

    // Analysis
    // Lambda Invariant Mass
    histos.add("Analysis/h1d_lstar_invm_US", "#Lambda(1520) M_{inv}", kTH1D, {axisInvM});
    histos.add("Analysis/h1d_lstar_invm_LS", "Like Signs M_{inv} p/#bar{p} K^{#mp}", kTH1D, {axisInvM});
    histos.add("Analysis/h1d_lstar_invm_PP", "Like Signs M_{inv} p K^{+}", kTH1D, {axisInvM});
    histos.add("Analysis/h1d_lstar_invm_MM", "Like Signs M_{inv} #bar{p} K^{-}", kTH1D, {axisInvM});
    histos.add("Analysis/h1d_lstar_invm_rot", "Rotated Spectra", kTH1D, {axisInvM});
    histos.add("Analysis/h1d_lstar_invm_US_mix", "Mixed Events M_{inv}", kTH1D, {axisInvM});
    histos.add("Analysis/h1d_lstar_invm_LS_mix", "Mixed Events M_{inv}", kTH1D, {axisInvM});
    histos.add("Analysis/h4d_lstar_invm_US", "THn #Lambda(1520)", kTHnSparseD, {axisInvM, axisPt, axisSp, axisCent});
    histos.add("Analysis/h4d_lstar_invm_LS", "THn Like Signs p/#bar{p} K^{#mp}", kTHnSparseD, {axisInvM, axisPt, axisSp, axisCent});
    histos.add("Analysis/h4d_lstar_invm_PP", "THn Like Signs p K^{+}", kTHnSparseD, {axisInvM, axisPt, axisSp, axisCent});
    histos.add("Analysis/h4d_lstar_invm_MM", "THn Like Signs #bar{p} K^{-}", kTHnSparseD, {axisInvM, axisPt, axisSp, axisCent});
    histos.add("Analysis/h4d_lstar_invm_rot", "THn Rotated", kTHnSparseD, {axisInvM, axisPt, axisSp, axisCent});
    histos.add("Analysis/h4d_lstar_invm_US_mix", "THn Mixed Events", kTHnSparseD, {axisInvM, axisPt, axisSp, axisCent});
    histos.add("Analysis/h4d_lstar_invm_LS_mix", "THn Mixed Events", kTHnSparseD, {axisInvM, axisPt, axisSp, axisCent});

    // MC
    if (doprocessMC) {
      histos.add("Event/h1d_rec_sph", "Reconstructed S_{0}", kTH1F, {axisSp});

      // histos.add("MCTruth/h1d_gen_sph", "Reconstructed S_{0}", kTH1F, {axisSp});
      histos.add("MCTruth/h1d_gen_posZ", "Generated PosZ", kTH1F, {axisPosZ});
      histos.add("MCTruth/h1d_ch_gen_phi", "Generated #phi distribution", kTH1F, {axisPhi});
      histos.add("MCTruth/h1d_pr_gen_eta", "Generated #eta Protons", kTH1F, {axisEta});
      histos.add("MCTruth/h1d_ka_gen_eta", "Generated #eta Kaons", kTH1F, {axisEta});

      // QAChecks
      histos.add("QAChecks/h1d_pr_rec_pt", "Reconstructed p_{T}-spectra Protons", kTH1F, {axisPt_pid});
      histos.add("QAChecks/h1d_ka_rec_pt", "Recondstucted p_{T}-spectra Kaons", kTH1F, {axisPt_pid});
      histos.add("QAChecks/h1d_pr_gen_pt", "Generated p_{T}-spectra Protons", kTH1F, {axisPt_pid});
      histos.add("QAChecks/h1d_ka_gen_pt", "Generated p_{T}-spectra Kaons", kTH1F, {axisPt_pid});

      // lstar
      histos.add("Analysis/h1d_gen_lstar", "Generated #Lambda(1520) p_{T}", kTH1D, {axisPt});
      histos.add("Analysis/h1d_gen_lstar_anti", "Generated #bar{#Lambda}(1520) p_{T}", kTH1D, {axisPt});
      histos.add("Analysis/h1d_rec_lstar", "Reconstructed #Lambda(1520) p_{T}", kTH1D, {axisPt});
      histos.add("Analysis/h1d_rec_lstar_anti", "Reconstructed #bar{#Lambda}(1520) p_{T}", kTH1D, {axisPt});
      histos.add("Analysis/h1d_rec_invm_lstar", "Recostructed #Lambda(1520)", kTH1D, {axisInvM});
      histos.add("Analysis/h1d_rec_invm_lstar_anti", "Recostructed #bar{#Lambda}(1520)", kTH1D, {axisInvM});
    }
  }

  template <typename T>
  bool selTracks(T const& track)
  {

    if (track.pt() < cPtMin)
      return false;

    if (std::abs(track.eta()) > cEtaCut)
      return false;

    if (std::abs(track.dcaZ()) > cDcaz)
      return false;

    if (std::abs(track.dcaXY()) > cDcaxy)
      return false;

    if (cPrimaryTrack && !track.isPrimaryTrack())
      return false;

    if (cGlobalWoDCATrack && !track.isGlobalTrackWoDCA())
      return false;

    if (cPVContributor && !track.isPVContributor())
      return false;

    return true;
  }
  // PID selection tools
  template <typename T>
  bool selectionPIDProton(const T& candidate, float p)
  {
    bool tpcPIDPassed{false}, tofPIDPassed{false};
    auto tpcPIDp = static_cast<std::vector<float>>(protonTPCPIDp);
    auto tpcPIDcut = static_cast<std::vector<float>>(protonTPCPIDcut);
    int nitr = static_cast<int>(tpcPIDp.size());

    float tpcNsigmaPi = std::abs(candidate.tpcNSigmaPi());
    float tpcNsigmaKa = std::abs(candidate.tpcNSigmaKa());
    float tpcNsigmaPr = std::abs(candidate.tpcNSigmaPr());
    float tofNsigmaPi = std::abs(candidate.tofNSigmaPi());
    float tofNsigmaKa = std::abs(candidate.tofNSigmaKa());
    float tofNsigmaPr = std::abs(candidate.tofNSigmaPr());

    float tpcTofNsigmaPi = tpcNsigmaPi * tpcNsigmaPi + tofNsigmaPi * tofNsigmaPi;
    float tpcTofNsigmaKa = tpcNsigmaKa * tpcNsigmaKa + tofNsigmaKa * tofNsigmaKa;
    float tpcTofNsigmaPr = tpcNsigmaPr * tpcNsigmaPr + tofNsigmaPr * tofNsigmaPr;
    float combinedCut = nsigmaCutCombinedProton * nsigmaCutCombinedProton;
    float combinedRejCut = cRejNsigmaTof * cRejNsigmaTpc;

    if (!cUseTpcOnly && candidate.hasTOF()) {
      if (tofNsigmaPr < cMaxTOFnSigmaProton && tofNsigmaPi > cRejNsigmaTof && tofNsigmaKa > cRejNsigmaTof) {
        tofPIDPassed = true;
      }
      // square cut
      if ((nsigmaCutCombinedProton < 0) && (tpcNsigmaPr < cMaxTPCnSigmaProton)) {
        tpcPIDPassed = true;
      }
      // circular cut
      if ((nsigmaCutCombinedProton > 0) && (tpcTofNsigmaPr < combinedCut && tpcTofNsigmaPi > combinedRejCut && tpcTofNsigmaKa > combinedRejCut)) {
        tofPIDPassed = true;
        tpcPIDPassed = true;
      }
    } else {
      tofPIDPassed = true;
      if (cUseTpcOnly) {
        if (tpcNsigmaPr < cMaxTPCnSigmaProton && tpcNsigmaPi > cRejNsigmaTpc && tpcNsigmaKa > cRejNsigmaTpc) {
          tpcPIDPassed = true;
        }
      } else {
        for (int i = 0; i < nitr - 1; ++i) {
          if (p >= tpcPIDp[i] && p < tpcPIDp[i + 1] && (tpcNsigmaPr < tpcPIDcut[i] && tpcNsigmaPi > cRejNsigmaTpc && tpcNsigmaKa > cRejNsigmaTpc)) {
            tpcPIDPassed = true;
          }
        }
        if (tpcPIDPassed && ((tpcNsigmaPr > tpcNsigmaKa) || (tpcNsigmaPr > tpcNsigmaPi))) {
          tpcPIDPassed = false;
        }
      }
    }
    if (tpcPIDPassed && tofPIDPassed) {
      return true;
    }
    return false;
  }
  template <typename T>
  bool selectionPIDKaon(const T& candidate, float p)
  {
    bool tpcPIDPassed{false}, tofPIDPassed{false};
    auto tpcPIDp = static_cast<std::vector<float>>(kaonTPCPIDp);
    auto tpcPIDcut = static_cast<std::vector<float>>(kaonTPCPIDcut);
    int nitr = static_cast<int>(tpcPIDp.size());

    float tpcNsigmaPi = std::abs(candidate.tpcNSigmaPi());
    float tpcNsigmaKa = std::abs(candidate.tpcNSigmaKa());
    float tpcNsigmaPr = std::abs(candidate.tpcNSigmaPr());
    float tofNsigmaPi = std::abs(candidate.tofNSigmaPi());
    float tofNsigmaKa = std::abs(candidate.tofNSigmaKa());
    float tofNsigmaPr = std::abs(candidate.tofNSigmaPr());

    float tpcTofNsigmaPi = tpcNsigmaPi * tpcNsigmaPi + tofNsigmaPi * tofNsigmaPi;
    float tpcTofNsigmaKa = tpcNsigmaKa * tpcNsigmaKa + tofNsigmaKa * tofNsigmaKa;
    float tpcTofNsigmaPr = tpcNsigmaPr * tpcNsigmaPr + tofNsigmaPr * tofNsigmaPr;
    float combinedCut = nsigmaCutCombinedKaon * nsigmaCutCombinedKaon;
    float combinedRejCut = cRejNsigmaTpc * cRejNsigmaTof;

    if (!cUseTpcOnly && candidate.hasTOF()) {
      if (tofNsigmaKa < cMaxTOFnSigmaKaon && tofNsigmaPi > cRejNsigmaTof && tofNsigmaPr > cRejNsigmaTof) {
        tofPIDPassed = true;
      }
      // square cut
      if ((nsigmaCutCombinedKaon < 0) && (tpcNsigmaKa < cMaxTPCnSigmaKaon)) {
        tpcPIDPassed = true;
      }
      // circular
      if ((nsigmaCutCombinedKaon > 0) && (tpcTofNsigmaKa < combinedCut && tpcTofNsigmaPi > combinedRejCut && tpcTofNsigmaPr > combinedRejCut)) {
        tofPIDPassed = true;
        tpcPIDPassed = true;
      }
    } else {
      tofPIDPassed = true;
      if (cUseTpcOnly) {
        if (tpcNsigmaKa < cMaxTPCnSigmaKaon && tpcNsigmaPi > cRejNsigmaTpc && tpcNsigmaPr > cRejNsigmaTpc) {
          tpcPIDPassed = true;
        }
      } else {
        for (int i = 0; i < nitr - 1; ++i) {
          if (p >= tpcPIDp[i] && p < tpcPIDp[i + 1] && (tpcNsigmaKa < tpcPIDcut[i] && tpcNsigmaPi > cRejNsigmaTpc && tpcNsigmaPr > cRejNsigmaTpc)) {
            tpcPIDPassed = true;
          }
        }
        if (tpcPIDPassed && ((tpcNsigmaKa > tpcNsigmaPr) || (tpcNsigmaKa > tpcNsigmaPi))) {
          tpcPIDPassed = false;
        }
      }
    }
    if (tpcPIDPassed && tofPIDPassed) {
      return true;
    }
    return false;
  }

  // kinematic cuts method
  template <typename trackType, typename T>
  bool kinCuts(trackType trkPr, trackType trkKa, T p, float& alpha)
  {
    // initialize
    std::vector<float> kinCutsPt = static_cast<std::vector<float>>(cKinCutsPt);
    std::vector<float> kinLowerCutsAlpha = static_cast<std::vector<float>>(cKinLowerCutsAlpha);
    std::vector<float> kinUpperCutsAlpha = static_cast<std::vector<float>>(cKinUpperCutsAlpha);
    int kinCutsSize = static_cast<int>(kinUpperCutsAlpha.size());

    TVector3 v1(trkPr.px(), trkPr.py(), trkPr.pz());
    TVector3 v2(trkKa.px(), trkKa.py(), trkKa.pz());
    alpha = v1.Angle(v2);

    for (int i = 0; i < kinCutsSize; ++i) {
      if ((p.Pt() > kinCutsPt[i] && p.Pt() <= kinCutsPt[i + 1]) && (alpha < kinLowerCutsAlpha[i] || alpha > kinUpperCutsAlpha[i])) {
        return false;
      }
    }

    return true;
  }

  template <typename T>
  void fillQAHistos(T const& track)
  {

    // get total momentum
    float p = RecoDecay::p(track.px(), track.py(), track.pz());

    // fill before QA first
    histos.fill(HIST("QAbefore/Proton/h2d_pr_nsigma_tpc_p"), p, track.tpcNSigmaPr());
    if (track.hasTOF()) {
      histos.fill(HIST("QAbefore/Proton/h2d_pr_nsigma_tof_p"), p, track.tofNSigmaPr());
      histos.fill(HIST("QAbefore/Proton/h2d_pr_nsigma_tof_vs_tpc"), track.tpcNSigmaPr(), track.tofNSigmaPr());
    }
    histos.fill(HIST("QAbefore/Kaon/h2d_ka_nsigma_tpc_p"), p, track.tpcNSigmaKa());
    if (track.hasTOF()) {
      histos.fill(HIST("QAbefore/Kaon/h2d_ka_nsigma_tof_p"), p, track.tofNSigmaKa());
      histos.fill(HIST("QAbefore/Kaon/h2d_ka_nsigma_tof_vs_tpc"), track.tpcNSigmaKa(), track.tofNSigmaKa());
    }

    // select particle (Proton/Kaon)
    // Proton
    if (selectionPIDProton(track, p)) {
      histos.fill(HIST("QAChecks/h1d_ka_pt"), track.pt());
      histos.fill(HIST("QAafter/Proton/h2d_pr_dca_z"), track.pt(), track.dcaZ());
      histos.fill(HIST("QAafter/Proton/h2d_pr_dca_xy"), track.pt(), track.dcaXY());
      histos.fill(HIST("QAafter/Proton/h2d_pr_dEdx_p"), p, track.tpcSignal());
      histos.fill(HIST("QAafter/Proton/h2d_pr_nsigma_tpc_p"), p, track.tpcNSigmaPr());
      histos.fill(HIST("QAafter/Proton/h2d_pr_nsigma_tpc_pt"), track.pt(), track.tpcNSigmaPr());
      if (!cUseTpcOnly && track.hasTOF()) {
        histos.fill(HIST("QAafter/Proton/h2d_pr_nsigma_tof_p"), p, track.tofNSigmaPr());
        histos.fill(HIST("QAafter/Proton/h2d_pr_nsigma_tof_pt"), track.pt(), track.tofNSigmaPr());
        histos.fill(HIST("QAafter/Proton/h2d_pr_nsigma_tof_vs_tpc"), track.tpcNSigmaPr(), track.tofNSigmaPr());
      }
    }
    // Kaon
    if (selectionPIDKaon(track, p)) {
      histos.fill(HIST("QAChecks/h1d_ka_pt"), track.pt());
      histos.fill(HIST("QAafter/Kaon/h2d_ka_dca_z"), track.pt(), track.dcaZ());
      histos.fill(HIST("QAafter/Kaon/h2d_ka_dca_xy"), track.pt(), track.dcaXY());
      histos.fill(HIST("QAafter/Kaon/h2d_ka_dEdx_p"), p, track.tpcSignal());
      histos.fill(HIST("QAafter/Kaon/h2d_ka_nsigma_tpc_p"), p, track.tpcNSigmaKa());
      histos.fill(HIST("QAafter/Kaon/h2d_ka_nsigma_tpc_pt"), track.pt(), track.tpcNSigmaKa());
      if (!cUseTpcOnly && track.hasTOF()) {
        histos.fill(HIST("QAafter/Kaon/h2d_ka_nsigma_tof_p"), p, track.tofNSigmaKa());
        histos.fill(HIST("QAafter/Kaon/h2d_ka_nsigma_tof_pt"), track.pt(), track.tofNSigmaKa());
        histos.fill(HIST("QAafter/Kaon/h2d_ka_nsigma_tof_vs_tpc"), track.tpcNSigmaKa(), track.tofNSigmaKa());
      }
    }
  }

  template <bool mix, bool mc, typename trackType>
  void fillInvMassHistos(trackType const& trk1, trackType const& trk2, float const& sph, float const& mult)
  {
    TLorentzVector p1, p2, p;
    TRandom* rn = new TRandom();
    float p_ptot = 0., k_ptot = 0.;

    for (auto const& [trkPr, trkKa] : soa::combinations(soa::CombinationsFullIndexPolicy(trk1, trk2))) {
      // Do not analyse same index tracks.
      if (trkPr.index() == trkKa.index())
        continue;

      // pT, DCA, Global Tracks and PVcontrib selection.
      if (!selTracks(trkPr) || !selTracks(trkKa))
        continue;

      p_ptot = RecoDecay::p(trkPr.px(), trkPr.py(), trkPr.pz());
      k_ptot = RecoDecay::p(trkKa.px(), trkKa.py(), trkKa.pz());

      // Apply PID Selection
      if (cUseOnlyTOFTrackPr && !trkPr.hasTOF())
        continue;
      if (cUseOnlyTOFTrackKa && !trkKa.hasTOF())
        continue;
      if (!selectionPIDProton(trkPr, p_ptot) || !selectionPIDKaon(trkKa, k_ptot))
        continue;

      // Invariant mass reconstruction.
      p1.SetXYZM(trkPr.px(), trkPr.py(), trkPr.pz(), MassProton);
      p2.SetXYZM(trkKa.px(), trkKa.py(), trkKa.pz(), MassKaonCharged);
      p = p1 + p2;

      // rapidity cut
      if (std::abs(p.Rapidity()) > 0.5)
        continue;

      // Get the opening angle b/w proton and kaon and kinematic cut flag.
      float alpha = 0;
      bool kinCutFlag{true};
      if (cKinCuts) {
        kinCutFlag = kinCuts(trkPr, trkKa, p, alpha);
        if constexpr (!mix && !mc) {
          histos.fill(HIST("QAChecks/h2d_lstar_alpha_vs_pt"), p.Pt(), alpha);
        }
      }

      // apply kincuts
      if (cKinCuts && !kinCutFlag) {
        continue;
      }

      // Fill Invariant Mass Histograms.
      if constexpr (!mix && !mc) {
        histos.fill(HIST("QAChecks/h2d_sel_kincuts_lstar_alpha_vs_pt"), p.Pt(), alpha);
        if (trkPr.sign() * trkKa.sign() < 0) {
          histos.fill(HIST("Analysis/h1d_lstar_invm_US"), p.M());
          histos.fill(HIST("Analysis/h4d_lstar_invm_US"), p.M(), p.Pt(), sph, mult);
          if (doRotate) {
            float theta = rn->Uniform(0.01, 0.1);
            p1.RotateZ(theta);
            p = p1 + p2;
            if (std::abs(p.Rapidity()) < 0.5) {
              histos.fill(HIST("Analysis/h1d_lstar_invm_rot"), p.M());
              histos.fill(HIST("Analysis/h4d_lstar_invm_rot"), p.M(), p.Pt(), sph, mult);
            }
          }
        } else {
          if (trkPr.sign() == 1) {
            histos.fill(HIST("Analysis/h1d_lstar_invm_LS"), p.M());
            histos.fill(HIST("Analysis/h4d_lstar_invm_LS"), p.M(), p.Pt(), sph, mult);
            histos.fill(HIST("Analysis/h1d_lstar_invm_PP"), p.M());
            histos.fill(HIST("Analysis/h4d_lstar_invm_PP"), p.M(), p.Pt(), sph, mult);
          } else {
            histos.fill(HIST("Analysis/h1d_lstar_invm_LS"), p.M());
            histos.fill(HIST("Analysis/h4d_lstar_invm_LS"), p.M(), p.Pt(), sph, mult);
            histos.fill(HIST("Analysis/h1d_lstar_invm_MM"), p.M());
            histos.fill(HIST("Analysis/h4d_lstar_invm_MM"), p.M(), p.Pt(), sph, mult);
          }
        }
      }

      if constexpr (mix) {
        if (trkPr.sign() * trkKa.sign() < 0) {
          histos.fill(HIST("Analysis/h1d_lstar_invm_US_mix"), p.M());
          histos.fill(HIST("Analysis/h4d_lstar_invm_US_mix"), p.M(), p.Pt(), sph, mult);
        } else {
          histos.fill(HIST("Analysis/h1d_lstar_invm_LS_mix"), p.M());
          histos.fill(HIST("Analysis/h4d_lstar_invm_LS_mix"), p.M(), p.Pt(), sph, mult);
        }
      }

      if constexpr (mc) {
        if (std::abs(trkPr.pdgCode()) != 2212 || std::abs(trkKa.pdgCode()) != 321)
          continue;

        if (trkPr.motherId() != trkKa.motherId())
          continue;

        if (std::abs(trkPr.motherPDG()) != 102134) // L* pdg_code = 102134
          continue;

        // MC histograms
        if (trkPr.motherPDG() > 0) {
          histos.fill(HIST("Analysis/h1d_rec_lstar"), p.Pt());
          histos.fill(HIST("Analysis/h1d_rec_invm_lstar"), p.M());
        } else {
          histos.fill(HIST("Analysis/h1d_rec_lstar_anti"), p.Pt());
          histos.fill(HIST("Analysis/h1d_rec_invm_lstar_anti"), p.M());
        }
      }
    }
  }

  using resoCols = soa::Join<aod::ResoCollisions, aod::ResoSpheroCollisions>;
  using resoTracks = aod::ResoTracks;

  void processData(resoCols::iterator const& collision, resoTracks const& tracks)
  {

    histos.fill(HIST("Event/h1d_ft0m_mult_percentile"), collision.cent());
    histos.fill(HIST("Event/h1d_spherocity"), collision.spherocity());
    histos.fill(HIST("Event/h2d_sph_vs_multpercentile"), collision.cent(), collision.spherocity());

    for (auto const& track : tracks) {
      if (!selTracks(track))
        continue;

      // QA histos
      fillQAHistos(track);
    }

    // get invariant mass histograms
    fillInvMassHistos<false, false>(tracks, tracks, collision.spherocity(), collision.cent());
  }

  PROCESS_SWITCH(lambdaAnalysis, processData, "Process for Same Event Data", true);

  void processMC(resoCols::iterator const& collision,
                 soa::Join<aod::ResoTracks, aod::ResoMCTracks> const& tracks)
  {
    histos.fill(HIST("Event/h1d_rec_sph"), collision.spherocity());

    // get MC reco pT-spectra
    for (auto const& track : tracks) {

      // get the reconstructed level pT spectra of protons and kaons
      if (!selTracks(track))
        continue;

      // QA histos
      fillQAHistos(track);

      float p = TMath::Sqrt(track.px() * track.px() + track.py() * track.py() + track.pz() * track.pz());

      if (selectionPIDKaon(track, p) && std::abs(track.pdgCode()) == 321) {
        histos.fill(HIST("QAChecks/h1d_ka_rec_pt"), track.pt());
      }

      if (selectionPIDProton(track, p) && std::abs(track.pdgCode()) == 2212) {
        histos.fill(HIST("QAChecks/h1d_pr_rec_pt"), track.pt());
      }
    }

    // get invariant mass histograms
    fillInvMassHistos<false, true>(tracks, tracks, collision.spherocity(), collision.cent());
  }
  PROCESS_SWITCH(lambdaAnalysis, processMC, "Process Event for MC", false);

  void processMCTrueDaughters(aod::McCollisions::iterator const& McCollision, aod::McParticles const& McParts)
  {

    // IP range selection
    if (std::abs(McCollision.posZ()) > 10.)
      return;

    histos.fill(HIST("MCTruth/h1d_gen_posZ"), McCollision.posZ());

    for (auto const& part : McParts) {

      // kinematic acceptance of particles
      if (part.pt() < cPtMin || std::abs(part.eta()) > cEtaCut || !part.isPhysicalPrimary())
        continue;

      histos.fill(HIST("MCTruth/h1d_ch_gen_phi"), part.phi());

      if (std::abs(part.pdgCode()) == 321) {
        histos.fill(HIST("MCTruth/h1d_ka_gen_eta"), part.eta());
        histos.fill(HIST("QAChecks/h1d_ka_gen_pt"), part.pt());
      }

      if (std::abs(part.pdgCode()) == 2212) {
        histos.fill(HIST("MCTruth/h1d_pr_gen_eta"), part.eta());
        histos.fill(HIST("QAChecks/h1d_pr_gen_pt"), part.pt());
      }
    }
  }
  PROCESS_SWITCH(lambdaAnalysis, processMCTrueDaughters, "Process Event for MC truth of protons and kaons", false);

  void processMCTrue(aod::ResoMCParents const& resoParents)
  {

    for (auto const& part : resoParents) {

      if (abs(part.pdgCode()) != 102134) // // L* pdg_code = 102134
        continue;
      if (abs(part.y()) > 0.5) { // rapidity cut
        continue;
      }

      bool pass1 = false;
      bool pass2 = false;

      if (abs(part.daughterPDG1()) == 2212 || abs(part.daughterPDG2()) == 2212) { // At least one decay to Proton
        pass1 = true;
      }
      if (abs(part.daughterPDG1()) == 321 || abs(part.daughterPDG2()) == 321) { // At least one decay to Kaon
        pass2 = true;
      }

      if (!pass1 || !pass2) // If we have both decay products
        continue;

      if (part.pdgCode() > 0)
        histos.fill(HIST("Analysis/h1d_gen_lstar"), part.pt());
      else
        histos.fill(HIST("Analysis/h1d_gen_lstar_anti"), part.pt());
    }
  }
  PROCESS_SWITCH(lambdaAnalysis, processMCTrue, "Process Event for MC", false);

  // Processing Event Mixing
  using BinningType1 = ColumnBinningPolicy<aod::collision::PosZ, aod::resocollision::Cent, aod::resocollision::Spherocity>;
  using BinningType2 = ColumnBinningPolicy<aod::collision::PosZ, aod::resocollision::Cent>;
  void processMix(resoCols& collisions, resoTracks const& tracks)
  {

    LOGF(debug, "Event Mixing Started");
    BinningType1 binningPositions1{{cMixVtxBins, cMixMultBins, cMixSphBins}, true};
    BinningType2 binningPositions2{{cMixVtxBins, cMixMultBins}, true};
    auto tracksTuple = std::make_tuple(tracks);
    if (cMixSph) {
      SameKindPair<resoCols, resoTracks, BinningType1> pairs{binningPositions1, cNumMixEv, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip
      for (auto& [c1, t1, c2, t2] : pairs) {
        fillInvMassHistos<true, false>(t1, t2, c1.spherocity(), c1.cent());
      }
    } else {
      SameKindPair<resoCols, resoTracks, BinningType2> pairs{binningPositions2, cNumMixEv, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip
      for (auto& [c1, t1, c2, t2] : pairs) {
        fillInvMassHistos<true, false>(t1, t2, c1.spherocity(), c1.cent());
      }
    }
  }

  PROCESS_SWITCH(lambdaAnalysis, processMix, "Process for Mixed Events", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<lambdaAnalysis>(cfgc)};
}
