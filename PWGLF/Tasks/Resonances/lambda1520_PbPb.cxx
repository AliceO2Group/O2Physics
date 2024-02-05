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

/// \file
/// ///        Invariant Mass Reconstruction of Lambda(1520) Resonance.
///
/// \author Yash Patley <yash.patley@cern.ch>
///  \author Nasir Mehdi Malik

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <TLorentzVector.h>
#include <TRandom.h>
#include <cmath>

#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/runDataProcessing.h"
#include "PWGLF/DataModel/LFResonanceTables.h"
#include "CommonConstants/PhysicsConstants.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

struct lambdaAnalysis_pb {

  SliceCache cache;
  Preslice<aod::ResoTracks> perRCol = aod::resodaughter::resoCollisionId;
  Preslice<aod::Tracks> perCollision = aod::track::collisionId;

  // Configurables.
  Configurable<int> nBinsPt{"nBinsPt", 100, "N bins in pT histogram"};
  Configurable<int> nBinsInvM{"nBinsInvM", 120, "N bins in InvMass histogram"};
  Configurable<bool> doRotate{"doRotate", true, "rotated inv mass spectra"};

  // Tracks
  Configurable<float> cPtMin{"cPtMin", 0.15, "Minimum Track pT"};
  Configurable<float> cEtaCut{"cEtaCut", 0.8, "Pseudorapidity cut"};
  Configurable<float> cDcaz{"cDcazMin", 1., "Minimum DCAz"};
  Configurable<float> cDcaxy{"cDcaxyMin", 0.1, "Minimum DCAxy"};
  Configurable<bool> cKinCuts{"cKinCuts", false, "Kinematic Cuts for p-K pair opening angle"};
  Configurable<bool> cITSRefit{"cITSRefit", false, "ITS TPC refit"};
  Configurable<bool> cTPCRefit{"cTPCRefit", false, "ITS TPC refit"};
  Configurable<bool> cPrimaryTrack{"cPrimaryTrack", true, "Primary track selection"};                    // kGoldenChi2 | kDCAxy | kDCAz
  Configurable<bool> cGlobalWoDCATrack{"cGlobalWoDCATrack", true, "Global track selection without DCA"}; // kQualityTracks (kTrackType | kTPCNCls | kTPCCrossedRows | kTPCCrossedRowsOverNCls | kTPCChi2NDF | kTPCRefit | kITSNCls | kITSChi2NDF | kITSRefit | kITSHits) | kInAcceptanceTracks (kPtRange | kEtaRange)
  Configurable<bool> cPVContributor{"cPVContributor", true, "PV contributor track selection"};           // PV Contriuibutor

  // PID Selections
  Configurable<bool> cUseOnlyTOFTrackPr{"cUseOnlyTOFTrackPr", false, "Use only TOF track for PID selection"}; // Use only TOF track for Proton PID selection
  Configurable<bool> cUseOnlyTOFTrackKa{"cUseOnlyTOFTrackKa", false, "Use only TOF track for PID selection"}; // Use only TOF track for Kaon PID selection
  Configurable<bool> cUseTpcOnly{"cUseTpcOnly", false, "Use TPC Only selection"};                             // TPC And TOF tracks
  Configurable<float> cRejNsigmaTpc{"cRejNsigmaTpc", 3.0, "Reject tracks to improve purity of TPC PID"};      // Reject missidentified particles when tpc bands merge
  Configurable<float> cRejNsigmaTof{"cRejNsigmaTof", 3.0, "Reject tracks to improve purity of TOF PID"};      // Reject missidentified particles when tpc bands merge
  // Proton
  Configurable<double> cMaxTPCnSigmaProton{"cMaxTPCnSigmaProton", 3.0, "TPC nSigma cut for Proton"};              // TPC
  Configurable<double> cMaxTOFnSigmaProton{"cMaxTOFnSigmaProton", 3.0, "TOF nSigma cut for Proton"};              // TOF
  Configurable<double> nsigmaCutCombinedProton{"nsigmaCutCombinedProton", 3.0, "Combined nSigma cut for Proton"}; // Combined
  Configurable<std::vector<float>> protonTPCPIDp{"protonTPCPIDp", {0, 0.5, 0.7, 0.8}, "p dependent TPC cuts protons"};
  Configurable<std::vector<float>> protonTPCPIDcut{"protonTPCPIDcut", {5., 3.5, 2.5}, "TPC nsigma cuts protons"};
  // Kaon
  Configurable<double> cMaxTPCnSigmaKaon{"cMaxTPCnSigmaKaon", 3.0, "TPC nSigma cut for Kaon"};              // TPC
  Configurable<double> cMaxTOFnSigmaKaon{"cMaxTOFnSigmaKaon", 3.0, "TOF nSigma cut for Kaon"};              // TOF
  Configurable<double> nsigmaCutCombinedKaon{"nsigmaCutCombinedKaon", 3.0, "Combined nSigma cut for Kaon"}; // Combined
  Configurable<std::vector<float>> kaonTPCPIDp{"kaonTPCPIDp", {0., 0.25, 0.3, 0.45}, "pT dependent TPC cuts kaons"};
  Configurable<std::vector<float>> kaonTPCPIDcut{"kaonTPCPIDcut", {6, 3.5, 2.5}, "TPC nsigma cuts kaons"};
  // Event Mixing.
  Configurable<int> cNumMixEv{"cNumMixEv", 20, "Number of Events to be mixed"};

  ConfigurableAxis cMixVtxBins{"cMixVtxBins", {VARIABLE_WIDTH, -10.0f, -9.f, -8.f, -7.f, -6.f, -5.f, -4.f, -3.f, -2.f, -1.f, 0.f, 1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis cMixMultBins{"cMixMultBins", {VARIABLE_WIDTH, 0.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.0f, 200.0f}, "Mixing bins - multiplicity"};

  // Histogram Registry.
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {

    // Define Axis.
    const AxisSpec axisCent(110, 0, 110, "FT0 (%)");
    const AxisSpec axisP_pid(200, 0., 10., "p (GeV/c)");
    const AxisSpec axisPt_pid(200, 0., 10., "p_{T} (GeV/c)");
    const AxisSpec axisPt(nBinsPt, 0., 10., "p_{T} (GeV/c)");
    const AxisSpec axisEta(40, -1, 1, "#eta");
    const AxisSpec axisDCAz(500, -0.5, 0.5, {"DCA_{z} (cm)"});
    const AxisSpec axisDCAxy(240, -0.12, 0.12, {"DCA_{xy} (cm)"});
    const AxisSpec axisTPCNCls(200, 0, 200, {"TPCNCls"});
    const AxisSpec axisTPCNsigma(401, -10.025, 10.025, {"n#sigma^{TPC}"});
    const AxisSpec axisTOFNsigma(401, -10.025, 10.025, {"n#sigma^{TOF}"});
    const AxisSpec axisdEdx(380, 10, 200, {"#frac{dE}{dx}"});
    const AxisSpec axisInvM(nBinsInvM, 1.44, 2.04, {"M_{inv} (GeV/c^{2})"});

    // Create Histograms.
    // Event
    histos.add("Event/h1d_ft0_mult_percentile", "FT0 (%)", kTH1F, {axisCent});

    histos.add("Event/mixing_vzVsmultpercentile", "FT0(%)", kTH1F, {axisCent});

    // QA Before
    histos.add("QAbefore/Proton/h2d_pr_nsigma_tpc_p", "n#sigma^{TPC} Protons", kTH2F, {axisP_pid, axisTPCNsigma});
    histos.add("QAbefore/Proton/h2d_pr_nsigma_tof_p", "n#sigma^{TOF} Protons", kTH2F, {axisP_pid, axisTOFNsigma});
    histos.add("QAbefore/Proton/h2d_pr_nsigma_tof_vs_tpc", "n#sigma^{TPC} vs n#sigma^{TOF} Protons", kTH2F, {axisTPCNsigma, axisTOFNsigma});
    histos.add("QAbefore/Kaon/h2d_ka_nsigma_tpc_p", "n#sigma^{TPC} Kaons", kTH2F, {axisP_pid, axisTPCNsigma});
    histos.add("QAbefore/Kaon/h2d_ka_nsigma_tof_p", "n#sigma^{TOF} Kaons", kTH2F, {axisP_pid, axisTOFNsigma});
    histos.add("QAbefore/Kaon/h2d_ka_nsigma_tof_vs_tpc", "n#sigma^{TPC} vs n#sigma^{TOF} Kaons", kTH2F, {axisTPCNsigma, axisTOFNsigma});

    // QA After
    histos.add("QAafter/Proton/h1d_pr_pt", "p_{T}-spectra Protons", kTH1F, {axisPt_pid});
    histos.add("QAafter/Proton/h2d_pr_dca_z", "dca_{z} Protons", kTH2F, {axisPt_pid, axisDCAz});
    histos.add("QAafter/Proton/h2d_pr_dca_xy", "dca_{xy} Protons", kTH2F, {axisPt_pid, axisDCAxy});
    histos.add("QAafter/Proton/h2d_pr_dEdx_p", "TPC Signal Protons", kTH2F, {axisP_pid, axisdEdx});
    histos.add("QAafter/Proton/h2d_pr_nsigma_tpc_pt", " Protons", kTH2F, {axisPt_pid, axisTPCNsigma});
    histos.add("QAafter/Proton/h2d_Prpi_nsigma_tpc_p", " Protons pion", kTH2F, {axisPt_pid, axisTPCNsigma});
    histos.add("QAafter/Proton/h2d_Prka_nsigma_tpc_p", " Protons kaon", kTH2F, {axisPt_pid, axisTPCNsigma});
    histos.add("QAafter/Proton/h2d_pr_nsigma_tpc_p", " Protons", kTH2F, {axisP_pid, axisTPCNsigma});
    histos.add("QAafter/Proton/h2d_pr_nsigma_tof_pt", " Protons", kTH2F, {axisPt_pid, axisTOFNsigma});
    histos.add("QAafter/Proton/h2d_pr_nsigma_tof_p", " Protons", kTH2F, {axisP_pid, axisTOFNsigma});
    histos.add("QAafter/Proton/h2d_Prpi_nsigma_tof_p", " Protons pion", kTH2F, {axisP_pid, axisTOFNsigma});
    histos.add("QAafter/Proton/h2d_Prka_nsigma_tof_p", " Protons kaon", kTH2F, {axisP_pid, axisTOFNsigma});
    histos.add("QAafter/Proton/h2d_pr_nsigma_tof_vs_tpc", "n#sigma(TOF) vs n#sigma(TPC) Protons", kTH2F, {axisTPCNsigma, axisTOFNsigma});
    histos.add("QAafter/Kaon/h1d_ka_pt", "p_{T}-spectra Kaons", kTH1F, {axisPt_pid});
    histos.add("QAafter/Kaon/h2d_ka_dca_z", "dca_{z} Kaons", kTH2F, {axisPt_pid, axisDCAz});
    histos.add("QAafter/Kaon/h2d_ka_dca_xy", "dca_{xy} Kaons", kTH2F, {axisPt_pid, axisDCAxy});
    histos.add("QAafter/Kaon/h2d_ka_dEdx_p", "TPC Signal Kaon", kTH2F, {axisP_pid, axisdEdx});
    histos.add("QAafter/Kaon/h2d_Kapi_nsigma_tpc_p", " Kaons pion", kTH2F, {axisPt_pid, axisTPCNsigma});
    histos.add("QAafter/Kaon/h2d_Kapr_nsigma_tpc_p", " Kaons proton", kTH2F, {axisP_pid, axisTPCNsigma});
    histos.add("QAafter/Kaon/h2d_ka_nsigma_tpc_pt", " Kaons", kTH2F, {axisPt_pid, axisTPCNsigma});
    histos.add("QAafter/Kaon/h2d_ka_nsigma_tpc_p", " Kaons", kTH2F, {axisP_pid, axisTPCNsigma});
    histos.add("QAafter/Kaon/h2d_ka_nsigma_tof_pt", " Kaons", kTH2F, {axisPt_pid, axisTOFNsigma});
    histos.add("QAafter/Kaon/h2d_ka_nsigma_tof_p", " Kaons", kTH2F, {axisP_pid, axisTOFNsigma});
    histos.add("QAafter/Kaon/h2d_Kapi_nsigma_tof_p", " Kaons pion", kTH2F, {axisP_pid, axisTOFNsigma});
    histos.add("QAafter/Kaon/h2d_Kapr_nsigma_tof_p", " Kaons proton", kTH2F, {axisP_pid, axisTOFNsigma});
    histos.add("QAafter/Kaon/h2d_ka_nsigma_tof_vs_tpc", "n#sigma(TOF) vs n#sigma(TPC) Kaons", kTH2F, {axisTPCNsigma, axisTOFNsigma});

    // QA checks for protons and kaons
    histos.add("QAChecks/h1d_pr_pt", "p_{T}-spectra Protons", kTH1F, {axisPt_pid});
    histos.add("QAChecks/h1d_ka_pt", "p_{T}-spectra Kaons", kTH1F, {axisPt_pid});
    histos.add("QAChecks/h1d_pr_rec_pt", "Reconstructed p_{T}-spectra Protons", kTH1F, {axisPt_pid});
    histos.add("QAChecks/h1d_ka_rec_pt", "Recondstucted p_{T}-spectra Kaons", kTH1F, {axisPt_pid});
    histos.add("QAChecks/h1d_pr_gen_pt", "Generated p_{T}-spectra Protons", kTH1F, {axisPt_pid});
    histos.add("QAChecks/h1d_ka_gen_pt", "Generated p_{T}-spectra Kaons", kTH1F, {axisPt_pid});

    // Analysis
    // Lambda Invariant Mass
    histos.add("Analysis/h1d_lstar_invm_US", "#Lambda(1520) M_{inv}", kTH1D, {axisInvM});
    histos.add("Analysis/h1d_lstar_invm_PP", "Like Signs M_{inv} p K^{+}", kTH1D, {axisInvM});
    histos.add("Analysis/h1d_lstar_invm_MM", "Like Signs M_{inv} #bar{p} K^{-}", kTH1D, {axisInvM});
    histos.add("Analysis/h1d_lstar_invm_rot", "Rotated Spectra", kTH1D, {axisInvM});
    histos.add("Analysis/h1d_lstar_invm_US_mix", "Mixed Events M_{inv}", kTH1D, {axisInvM});
    histos.add("Analysis/h1d_lstar_invm_LS_mix", "Mixed Events M_{inv}", kTH1D, {axisInvM});
    histos.add("Analysis/h4d_lstar_invm_US", "THn #Lambda(1520)", kTHnSparseD, {axisInvM, axisPt, axisCent});
    histos.add("Analysis/h4d_lstar_invm_PP", "THn Like Signs p K^{+}", kTHnSparseD, {axisInvM, axisPt, axisCent});
    histos.add("Analysis/h4d_lstar_invm_MM", "THn Like Signs #bar{p} K^{-}", kTHnSparseD, {axisInvM, axisPt, axisCent});
    histos.add("Analysis/h4d_lstar_invm_rot", "THn Rotated", kTHnSparseD, {axisInvM, axisPt, axisCent});
    histos.add("Analysis/h4d_lstar_invm_US_mix", "THn Mixed Events", kTHnSparseD, {axisInvM, axisPt, axisCent});
    histos.add("Analysis/h4d_lstar_invm_LS_mix", "THn Mixed Events", kTHnSparseD, {axisInvM, axisPt, axisCent});

    // MC
    if (doprocessMC) {

      histos.add("Event/h1d_rec_cent", "Reconstructed  FT0(%)", kTH2F, {axisCent});
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
    // float tofNsigmaEl = std::abs();

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
      }
    }
    if (tpcPIDPassed && tofPIDPassed) {
      return true;
    }
    return false;
  }

  template <bool mix, bool mc, typename trackType>
  void fillDataHistos(trackType const& trk1, trackType const& trk2, float const& mult)
  {
    TLorentzVector p1, p2, p;
    TRandom* rn = new TRandom();
    float p_ptot = 0., k_ptot = 0.;

    for (auto const& [trkPr, trkKa] : soa::combinations(soa::CombinationsFullIndexPolicy(trk1, trk2))) {
      // Do not analyse same index tracks.
      if (trkPr.index() == trkKa.index() && !mix)
        continue;

      // pT, DCA, Global Tracks and PVcontrib selection.
      if (!selTracks(trkPr) || !selTracks(trkKa))
        continue;

      if (cITSRefit && !trkPr.passedITSRefit())
        continue;
      if (cITSRefit && !trkKa.passedITSRefit())
        continue;
      if (cTPCRefit && !trkPr.passedTPCRefit())
        continue;
      if (cTPCRefit && !trkKa.passedTPCRefit())
        continue;
      auto _pxPr = trkPr.px();
      auto _pyPr = trkPr.py();
      auto _pzPr = trkPr.pz();
      auto _pxKa = trkKa.px();
      auto _pyKa = trkKa.py();
      auto _pzKa = trkKa.pz();

      p_ptot = TMath::Sqrt(_pxPr * _pxPr + _pyPr * _pyPr + _pzPr * _pzPr);
      k_ptot = TMath::Sqrt(_pxKa * _pxKa + _pyKa * _pyKa + _pzKa * _pzKa);

      // Fill QA before track selection.
      if (!mix) {
        auto _tpcnsigmaPr = trkPr.tpcNSigmaPr();

        histos.fill(HIST("QAbefore/Proton/h2d_pr_nsigma_tpc_p"), p_ptot, _tpcnsigmaPr);
        // histos.fill(HIST("QAbefore/Proton/h2d_prel_nsigma_tpc_p"), p_ptot, trkPr.tpcNSigmaEl());
        if (trkPr.hasTOF()) {
          auto _tofnsigmaPr = trkPr.tofNSigmaPr();
          histos.fill(HIST("QAbefore/Proton/h2d_pr_nsigma_tof_p"), p_ptot, _tofnsigmaPr);
          histos.fill(HIST("QAbefore/Proton/h2d_pr_nsigma_tof_vs_tpc"), _tpcnsigmaPr, _tofnsigmaPr);
        }
        auto _tpcnsigmaKa = trkKa.tpcNSigmaKa();
        histos.fill(HIST("QAbefore/Kaon/h2d_ka_nsigma_tpc_p"), k_ptot, _tpcnsigmaKa);
        if (trkKa.hasTOF()) {
          auto _tofnsigmaKa = trkKa.tofNSigmaKa();
          histos.fill(HIST("QAbefore/Kaon/h2d_ka_nsigma_tof_p"), k_ptot, _tofnsigmaKa);
          histos.fill(HIST("QAbefore/Kaon/h2d_ka_nsigma_tof_vs_tpc"), _tpcnsigmaKa, _tofnsigmaKa);
        }
      }

      // Apply PID Selection
      if (cUseOnlyTOFTrackPr && !trkPr.hasTOF())
        continue;
      if (cUseOnlyTOFTrackKa && !trkKa.hasTOF())
        continue;
      if (!selectionPIDProton(trkPr, p_ptot) || !selectionPIDKaon(trkKa, k_ptot))
        continue;

      // Fill QA after track selection.
      if constexpr (!mix) {
        auto _ptPr = trkPr.pt();
        auto _tpcnsigmaPr = trkPr.tpcNSigmaPr();

        // Proton
        histos.fill(HIST("QAafter/Proton/h1d_pr_pt"), _ptPr);
        histos.fill(HIST("QAafter/Proton/h2d_pr_dca_z"), _ptPr, trkPr.dcaZ());
        histos.fill(HIST("QAafter/Proton/h2d_pr_dca_xy"), _ptPr, trkPr.dcaXY());
        histos.fill(HIST("QAafter/Proton/h2d_pr_dEdx_p"), p_ptot, trkPr.tpcSignal());
        histos.fill(HIST("QAafter/Proton/h2d_Prpi_nsigma_tpc_p"), p_ptot, trkPr.tpcNSigmaPi());
        histos.fill(HIST("QAafter/Proton/h2d_Prka_nsigma_tpc_p"), p_ptot, trkPr.tpcNSigmaKa());
        histos.fill(HIST("QAafter/Proton/h2d_pr_nsigma_tpc_p"), p_ptot, _tpcnsigmaPr);
        histos.fill(HIST("QAafter/Proton/h2d_pr_nsigma_tpc_pt"), _ptPr, _tpcnsigmaPr);
        if (!cUseTpcOnly && trkPr.hasTOF()) {
          auto _tofnsigmaPr = trkPr.tofNSigmaPr();
          histos.fill(HIST("QAafter/Proton/h2d_pr_nsigma_tof_p"), p_ptot, _tofnsigmaPr);
          histos.fill(HIST("QAafter/Proton/h2d_pr_nsigma_tof_pt"), _ptPr, _tofnsigmaPr);
          histos.fill(HIST("QAafter/Proton/h2d_Prpi_nsigma_tof_p"), p_ptot, trkPr.tofNSigmaPi());
          histos.fill(HIST("QAafter/Proton/h2d_Prka_nsigma_tof_p"), p_ptot, trkPr.tofNSigmaKa());
          histos.fill(HIST("QAafter/Proton/h2d_pr_nsigma_tof_vs_tpc"), _tpcnsigmaPr, _tofnsigmaPr);
        }
        auto _ptKa = trkKa.pt();
        auto _tpcnsigmaKa = trkKa.tpcNSigmaKa();

        // Kaon
        histos.fill(HIST("QAafter/Kaon/h1d_ka_pt"), _ptKa);
        histos.fill(HIST("QAafter/Kaon/h2d_ka_dca_z"), _ptKa, trkKa.dcaZ());
        histos.fill(HIST("QAafter/Kaon/h2d_ka_dca_xy"), _ptKa, trkKa.dcaXY());
        histos.fill(HIST("QAafter/Kaon/h2d_ka_dEdx_p"), k_ptot, trkKa.tpcSignal());
        histos.fill(HIST("QAafter/Kaon/h2d_Kapi_nsigma_tpc_p"), k_ptot, trkKa.tpcNSigmaPi());
        histos.fill(HIST("QAafter/Kaon/h2d_Kapr_nsigma_tpc_p"), k_ptot, trkKa.tpcNSigmaPr());
        histos.fill(HIST("QAafter/Kaon/h2d_ka_nsigma_tpc_p"), k_ptot, _tpcnsigmaKa);
        histos.fill(HIST("QAafter/Kaon/h2d_ka_nsigma_tpc_pt"), _ptKa, _tpcnsigmaKa);
        if (!cUseTpcOnly && trkKa.hasTOF()) {
          auto _tofnsigmaKa = trkKa.tofNSigmaKa();
          histos.fill(HIST("QAafter/Kaon/h2d_ka_nsigma_tof_p"), k_ptot, _tofnsigmaKa);
          histos.fill(HIST("QAafter/Kaon/h2d_ka_nsigma_tof_pt"), _ptKa, _tofnsigmaKa);
          histos.fill(HIST("QAafter/Kaon/h2d_Kapi_nsigma_tof_p"), k_ptot, trkKa.tofNSigmaPi());
          histos.fill(HIST("QAafter/Kaon/h2d_Kapr_nsigma_tof_p"), k_ptot, trkPr.tofNSigmaPr());
          histos.fill(HIST("QAafter/Kaon/h2d_ka_nsigma_tof_vs_tpc"), _tpcnsigmaKa, _tofnsigmaKa);
        }
      }

      // Invariant mass reconstruction.
      p1.SetXYZM(trkPr.px(), trkPr.py(), trkPr.pz(), MassProton);
      p2.SetXYZM(trkKa.px(), trkKa.py(), trkKa.pz(), MassKaonCharged);
      p = p1 + p2;

      if (std::abs(p.Rapidity()) > 0.5)
        continue;

      // Apply kinematic cuts.
      if (cKinCuts) {
        TVector3 v1(_pxPr, _pyPr, _pzPr);
        TVector3 v2(_pxKa, _pyKa, _pzKa);
        float alpha = v1.Angle(v2);
        if (alpha > 1.4 && alpha < 2.4)
          continue;
      }

      // Fill Invariant Mass Histograms.
      if constexpr (!mix && !mc) {
        if (trkPr.sign() * trkKa.sign() < 0) {
          histos.fill(HIST("Analysis/h1d_lstar_invm_US"), p.M());
          histos.fill(HIST("Analysis/h4d_lstar_invm_US"), p.M(), p.Pt(), mult);
          if (doRotate) {
            float theta = rn->Uniform(1.56, 1.58);
            p1.RotateZ(theta);
            p = p1 + p2;
            if (std::abs(p.Rapidity()) < 0.5) {
              histos.fill(HIST("Analysis/h1d_lstar_invm_rot"), p.M());
              histos.fill(HIST("Analysis/h4d_lstar_invm_rot"), p.M(), p.Pt(), mult);
            }
          }
        } else {
          if (trkPr.sign() == 1) {
            histos.fill(HIST("Analysis/h1d_lstar_invm_PP"), p.M());
            histos.fill(HIST("Analysis/h4d_lstar_invm_PP"), p.M(), p.Pt(), mult);
          } else {
            histos.fill(HIST("Analysis/h1d_lstar_invm_MM"), p.M());
            histos.fill(HIST("Analysis/h4d_lstar_invm_MM"), p.M(), p.Pt(), mult);
          }
        }
      }

      if constexpr (mix) {
        if (trkPr.sign() * trkKa.sign() < 0) {
          histos.fill(HIST("Analysis/h1d_lstar_invm_US_mix"), p.M());
          histos.fill(HIST("Analysis/h4d_lstar_invm_US_mix"), p.M(), p.Pt(), mult);
        } else {
          histos.fill(HIST("Analysis/h1d_lstar_invm_LS_mix"), p.M());
          histos.fill(HIST("Analysis/h4d_lstar_invm_LS_mix"), p.M(), p.Pt(), mult);
        }
      }

      if constexpr (mc) {
        if (std::abs(trkPr.pdgCode()) != 2212 || std::abs(trkKa.pdgCode()) != 321)
          continue;

        if (trkPr.motherId() != trkKa.motherId())
          continue;

        if (std::abs(trkPr.motherPDG()) != 3124) // L* pdg_code = 3124
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

  using resoCols = aod::ResoCollisions;
  using resoTracks = aod::ResoTracks;

  void processData(resoCols::iterator const& collision, resoTracks const& tracks)
  {

    histos.fill(HIST("Event/h1d_ft0_mult_percentile"), collision.cent());
    fillDataHistos<false, false>(tracks, tracks, collision.cent());

    // get proton and kaon pT-spectra
    for (auto const& track : tracks) {
      if (!selTracks(track))
        continue;

      float p = TMath::Sqrt(track.px() * track.px() + track.py() * track.py() + track.pz() * track.pz());

      if (selectionPIDKaon(track, p)) {
        histos.fill(HIST("QAChecks/h1d_ka_pt"), track.pt());
      }

      if (selectionPIDProton(track, p)) {
        histos.fill(HIST("QAChecks/h1d_pr_pt"), track.pt());
      }
    }
  }

  PROCESS_SWITCH(lambdaAnalysis_pb, processData, "Process for Same Event Data", true);

  void processMC(resoCols::iterator const& collision,
                 soa::Join<aod::ResoTracks, aod::ResoMCTracks> const& tracks)
  {

    fillDataHistos<false, true>(tracks, tracks, collision.cent());

    // get MC pT-spectra
    for (auto const& track : tracks) {

      // get the generated level pT spectra of protons and kaons
      if (std::abs(track.pdgCode()) == 321)
        histos.fill(HIST("QAChecks/h1d_ka_gen_pt"), track.pt());

      if (std::abs(track.pdgCode()) == 2212)
        histos.fill(HIST("QAChecks/h1d_pr_gen_pt"), track.pt());

      // get the reconstructed level pT spectra of protons and kaons
      if (!selTracks(track))
        continue;

      float p = TMath::Sqrt(track.px() * track.px() + track.py() * track.py() + track.pz() * track.pz());

      if (selectionPIDKaon(track, p) && std::abs(track.pdgCode()) == 321) {
        histos.fill(HIST("QAChecks/h1d_ka_rec_pt"), track.pt());
      }

      if (selectionPIDProton(track, p) && std::abs(track.pdgCode()) == 2212) {
        histos.fill(HIST("QAChecks/h1d_pr_rec_pt"), track.pt());
      }
    }
  }
  PROCESS_SWITCH(lambdaAnalysis_pb, processMC, "Process Event for MC", false);

  void processMCTrue(aod::ResoMCParents const& resoParents)
  {

    for (auto const& part : resoParents) {

      if (abs(part.pdgCode()) != 3124) // // L* pdg_code = 3124
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
  PROCESS_SWITCH(lambdaAnalysis_pb, processMCTrue, "Process Event for MC", false);

  using BinningType2 = ColumnBinningPolicy<aod::collision::PosZ, aod::resocollision::Cent>;

  void processMix(resoCols& collisions, resoTracks const& tracks)
  {

    LOGF(debug, "Event Mixing Started");

    BinningType2 binningPositions2{{cMixVtxBins, cMixMultBins}, true};
    auto tracksTuple = std::make_tuple(tracks);

    SameKindPair<resoCols, resoTracks, BinningType2> pairs{binningPositions2, cNumMixEv, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip
    for (auto& [c1, t1, c2, t2] : pairs) {
      histos.fill(HIST("Event/mixing_vzVsmultpercentile"), c1.cent());
      fillDataHistos<true, false>(t1, t2, c1.cent());
    }
  }

  PROCESS_SWITCH(lambdaAnalysis_pb, processMix, "Process for Mixed Events", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<lambdaAnalysis_pb>(cfgc)};
}
