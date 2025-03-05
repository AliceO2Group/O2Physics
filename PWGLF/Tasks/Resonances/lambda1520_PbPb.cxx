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
/// \author Nasir Mehdi Malik

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <TLorentzVector.h>
#include <TRandom.h>
#include <fairlogger/Logger.h>

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

  Configurable<bool> ConfEvtOccupancyInTimeRange{"ConfEvtOccupancyInTimeRange", false, "occupancy selection true or false"};
  Configurable<int> nBinsPt{"nBinsPt", 100, "N bins in pT histogram"};
  Configurable<int> nBinsInvM{"nBinsInvM", 120, "N bins in InvMass histogram"};
  Configurable<int> lambda1520id{"lambda1520id", 3124, "pdg"};
  Configurable<bool> doRotate{"doRotate", true, "rotated inv mass spectra"};

  // Tracks
  Configurable<float> cPtMin{"cPtMin", 0.15, "Minimum Track pT"};
  Configurable<float> cPMin{"cPMin", 0., "Minimum Track p"};
  Configurable<float> cEtaCut{"cEtaCut", 0.8, "Pseudorapidity cut"};
  Configurable<float> cDcaz{"cDcazMin", 1., "Minimum DCAz"};
  Configurable<float> cDcaxy{"cDcaxyMin", 0.1, "Minimum DCAxy"};
  Configurable<bool> isonlyQC{"isonlyQC", false, "only QC"};
  Configurable<bool> isDeepAngle{"isDeepAngle", false, "Deep Angle cut"};
  Configurable<double> cfgDeepAngle{"cfgDeepAngle", 0.04, "Deep Angle cut value"};
  Configurable<bool> cKinCuts{"cKinCuts", false, "Kinematic Cuts for p-K pair opening angle"};
  Configurable<bool> cPrimaryTrack{"cPrimaryTrack", true, "Primary track selection"};                    // kGoldenChi2 | kDCAxy | kDCAz
  Configurable<bool> cGlobalWoDCATrack{"cGlobalWoDCATrack", true, "Global track selection without DCA"}; // kQualityTracks (kTrackType | kTPCNCls | kTPCCrossedRows | kTPCCrossedRowsOverNCls | kTPCChi2NDF | kTPCRefit | kITSNCls | kITSChi2NDF | kITSRefit | kITSHits) | kInAcceptanceTracks (kPtRange | kEtaRange)
  Configurable<bool> cPVContributor{"cPVContributor", true, "PV contributor track selection"};           // PV Contriuibutor

  // PID Selections
  Configurable<bool> cUseOnlyTOFTrackPr{"cUseOnlyTOFTrackPr", false, "Use only TOF track for PID selection"}; // Use only TOF track for Proton PID selection
  Configurable<bool> cUseOnlyTOFTrackKa{"cUseOnlyTOFTrackKa", false, "Use only TOF track for PID selection"}; // Use only TOF track for Kaon PID selection
  Configurable<bool> cUseTpcOnly{"cUseTpcOnly", false, "Use TPC Only selection"};                             // TPC And TOF tracks
  Configurable<float> cRejNsigmaTpc{"cRejNsigmaTpc", 3.0, "Reject tracks to improve purity of TPC PID"};      // Reject missidentified particles when tpc bands merge
  Configurable<float> cRejNsigmaTpcPi{"cRejNsigmaTpcPi", 3.0, "Reject tracks to improve purity of TPC PID"};  // TPC And TOF tracks
                                                                                                              // Configurable<float> cRejNsigmaTpcPr{"cRejNsigmaTpcPr", 3.0, "Reject tracks to improve purity of TPC PID"};
  Configurable<float> cRejNsigmaTpcKa{"cRejNsigmaTpcKa", 3.0, "Reject tracks to improve purity of TPC PID"};
  Configurable<float> cRejNsigmakTpcPi{"cRejNsigmakTpcPi", 3.0, "Reject tracks to improve purity of TPC PID"};
  Configurable<float> cRejNsigmakTpcPr{"cRejNsigmakTpcPr", 3.0, "Reject tracks to improve purity of TPC PID"};
  Configurable<float> minnsigmatpcKa{"minnsigmatpcKa", -6.0, "Reject tracks to improve purity of TPC PID"};
  Configurable<float> minnsigmatpcPr{"minnsigmatpcPr", -6.0, "Reject tracks to improve purity of TPC PID"};
  Configurable<float> minnsigmatofKa{"minnsigmatofKa", -6.0, "Reject tracks to improve purity of TofPID"};
  Configurable<float> minnsigmatofPr{"minnsigmatofPr", -6.0, "Reject tracks to improve purity of Tof PID"};
  Configurable<float> minnsigmatpctofKa{"minnsigmatpctofKa", -6.0, "Reject tracks to improve purity of TPC PID"};
  Configurable<float> minnsigmatpctofPr{"minnsigmatpctofPr", -6.0, "Reject tracks to improve purity of TPC PID"};
  // Configurable<float> cRejNsigmaTpcPr{"cRejNsigmaTpcPr", 3.0, "Reject tracks to improve purity of TPC PID"};
  Configurable<float> cRejNsigmaTpcVeto{"cRejNsigmaTpcVeto", 3.0, "Reject tracks to improve purity of TPC PID"}; // Reject missidentified particles when tpc bands merge
  Configurable<float> cRejNsigmaTof{"cRejNsigmaTof", 3.0, "Reject tracks to improve purity of TOF PID"};         // Reject missidentified particles when tpc bands merge
  // Proton
  Configurable<double> cMaxTPCnSigmaProton{"cMaxTPCnSigmaProton", 3.0, "TPC nSigma cut for Proton"}; // TPC
  //  Configurable<double> cMaxTOFnSigmaProton{"cMaxTOFnSigmaProton", 3.0, "TOF nSigma cut for Proton"};              // TOF
  Configurable<double> nsigmaCutCombinedProton{"nsigmaCutCombinedProton", 3.0, "Combined nSigma cut for Proton"}; // Combined
  Configurable<std::vector<float>> protonTPCPIDp{"protonTPCPIDp", {0, 0.5, 0.7, 0.8}, "p dependent TPC cuts protons"};
  Configurable<std::vector<float>> protonTPCPIDcut{"protonTPCPIDcut", {5., 3.5, 2.5}, "TPC nsigma cuts protons"};
  Configurable<std::vector<float>> protonTOFPIDp{"protonTOFPIDp", {0., 999.}, "p dependent TOF cuts protons"};
  Configurable<std::vector<float>> protonTOFPIDcut{"protonTOFPIDcut", {3.0}, "TOF nsigma cuts protons"};
  // Kaon
  Configurable<double> cMaxTPCnSigmaKaon{"cMaxTPCnSigmaKaon", 3.0, "TPC nSigma cut for Kaon"}; // TPC
  // Configurable<double> cMaxTOFnSigmaKaon{"cMaxTOFnSigmaKaon", 3.0, "TOF nSigma cut for Kaon"};              // TOF
  Configurable<double> nsigmaCutCombinedKaon{"nsigmaCutCombinedKaon", 3.0, "Combined nSigma cut for Kaon"}; // Combined
  Configurable<std::vector<float>> kaonTPCPIDp{"kaonTPCPIDp", {0., 0.25, 0.3, 0.45}, "pT dependent TPC cuts kaons"};
  Configurable<std::vector<float>> kaonTPCPIDcut{"kaonTPCPIDcut", {6, 3.5, 2.5}, "TPC nsigma cuts kaons"};
  Configurable<std::vector<float>> kaonTOFPIDp{"kaonTOFPIDp", {0., 999.}, "p dependent TOF cuts kaons"};
  Configurable<std::vector<float>> kaonTOFPIDcut{"kaonTOFPIDcut", {3.0}, "TOF nsigma cuts kaons"};
  // Event Mixing.
  Configurable<int> cNumMixEv{"cNumMixEv", 20, "Number of Events to be mixed"};

  ConfigurableAxis cMixVtxBins{"cMixVtxBins", {VARIABLE_WIDTH, -10.0f, -9.f, -8.f, -7.f, -6.f, -5.f, -4.f, -3.f, -2.f, -1.f, 0.f, 1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis cMixMultBins{"cMixMultBins", {VARIABLE_WIDTH, 0.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.0f, 200.0f}, "Mixing bins - multiplicity"};
  ConfigurableAxis cMixEPAngle{"cMixEPAngle", {VARIABLE_WIDTH, -1.5708f, -1.25664f, -0.942478f, -0.628319f, 0.f, 0.628319f, 0.942478f, 1.25664f, 1.5708f}, "event plane"};
  ConfigurableAxis occupancy_bins{"occupancy_bins", {VARIABLE_WIDTH, 0.0, 100, 500, 600, 1000, 1100, 1500, 1600, 2000, 2100, 2500, 2600, 3000, 3100, 3500, 3600, 4000, 4100, 4500, 4600, 5000, 5100, 9999}, "Binning of the occupancy axis"};
  // Histogram Registry.
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {
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
    const AxisSpec axisVz(120, -12, 12, {"vz"});
    const AxisSpec axisEP(120, -3.14, 3.14, {"#theta"});
    const AxisSpec axisInvM(nBinsInvM, 1.44, 2.04, {"M_{inv} (GeV/c^{2})"});
    AxisSpec axisOccupancy = {occupancy_bins, "Occupancy [-40,100]"};

    histos.add("Event/h1d_ft0_mult_percentile", "FT0 (%)", kTH2F, {axisCent, axisOccupancy});
    if (doprocessMix || doprocessMixDF || doprocessMixepDF) {
      histos.add("Event/mixing_vzVsmultpercentile", "FT0(%)", kTH3F, {axisCent, axisVz, axisEP});
    }
    // QA Before
    histos.add("QAbefore/Proton/h2d_pr_nsigma_tpc_p", "n#sigma^{TPC} Protons", kTH2F, {axisP_pid, axisTPCNsigma});
    histos.add("QAbefore/Proton/h2d_pr_nsigma_tof_p", "n#sigma^{TOF} Protons", kTH2F, {axisP_pid, axisTOFNsigma});
    histos.add("QAbefore/Proton/h2d_pr_nsigma_tof_vs_tpc", "n#sigma^{TPC} vs n#sigma^{TOF} Protons", kTH2F, {axisTPCNsigma, axisTOFNsigma});
    histos.add("QAbefore/Kaon/h2d_ka_nsigma_tpc_p", "n#sigma^{TPC} Kaons", kTH2F, {axisP_pid, axisTPCNsigma});
    histos.add("QAbefore/Kaon/h2d_ka_nsigma_tof_p", "n#sigma^{TOF} Kaons", kTH2F, {axisP_pid, axisTOFNsigma});
    histos.add("QAbefore/Kaon/h2d_ka_nsigma_tof_vs_tpc", "n#sigma^{TPC} vs n#sigma^{TOF} Kaons", kTH2F, {axisTPCNsigma, axisTOFNsigma});

    // QA After
    histos.add("QAafter/Proton/hd_pr_pt", "p_{T}-spectra Protons", kTH2F, {axisPt_pid, axisCent});
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
    histos.add("QAafter/Kaon/hd_ka_pt", "p_{T}-spectra Kaons", kTH2F, {axisPt_pid, axisCent});
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

    // Analysis
    // Lambda Invariant Mass
    if (!doprocessMC) {
      histos.add("Analysis/h4d_lstar_invm_US_PM", "THn #Lambda(1520)", kTHnSparseF, {axisInvM, axisPt, axisCent, axisOccupancy});
      histos.add("Analysis/h4d_lstar_invm_US_MP", "THn #bar #Lambda(1520)", kTHnSparseF, {axisInvM, axisPt, axisCent, axisOccupancy});
      histos.add("Analysis/h4d_lstar_invm_PP", "THn Like Signs p K^{+}", kTHnSparseF, {axisInvM, axisPt, axisCent, axisOccupancy});
      histos.add("Analysis/h4d_lstar_invm_MM", "THn Like Signs #bar{p} K^{-}", kTHnSparseF, {axisInvM, axisPt, axisCent, axisOccupancy});
      histos.add("Analysis/h4d_lstar_invm_rot", "THn Rotated", kTHnSparseF, {axisInvM, axisPt, axisCent, axisOccupancy});
      histos.add("Analysis/h4d_lstar_invm_US_PM_mix", "THn Mixed Events", kTHnSparseF, {axisInvM, axisPt, axisCent, axisOccupancy});
      histos.add("Analysis/h4d_lstar_invm_US_MP_mix", "THn anti Mixed Events", kTHnSparseF, {axisInvM, axisPt, axisCent, axisOccupancy});
      histos.add("Analysis/h4d_lstar_invm_LS_PP_mix", "THn Mixed Events PP", kTHnSparseF, {axisInvM, axisPt, axisCent, axisOccupancy});
      histos.add("Analysis/h4d_lstar_invm_LS_MM_mix", "THn Mixed Events MM", kTHnSparseF, {axisInvM, axisPt, axisCent, axisOccupancy});
    }
    // MC
    if (doprocessMC) {

      histos.add("QAChecks/h1d_pr_rec_pt", "Reconstructed p_{T}-spectra Protons", kTH1F, {axisPt_pid});
      histos.add("QAChecks/h1d_ka_rec_pt", "Recondstucted p_{T}-spectra Kaons", kTH1F, {axisPt_pid});
      histos.add("QAChecks/h1d_pr_gen_pt", "Generated p_{T}-spectra Protons", kTH1F, {axisPt_pid});
      histos.add("QAChecks/h1d_ka_gen_pt", "Generated p_{T}-spectra Kaons", kTH1F, {axisPt_pid});

      histos.add("Analysis/h3d_gen_lstar_PM", "Generated #Lambda(1520) p_{T}", kTHnSparseF, {axisInvM, axisPt, axisCent});
      histos.add("Analysis/h3d_gen_lstar_MP", "Generated #bar{#Lambda}(1520) p_{T}", kTHnSparseF, {axisInvM, axisPt, axisCent});
      histos.add("Analysis/h3d_rec_lstar_PM", "Reconstructed #Lambda(1520) p_{T}", kTHnSparseF, {axisInvM, axisPt, axisCent});
      histos.add("Analysis/h3d_rec_lstar_MP", "Reconstructed #bar{#Lambda}(1520) p_{T}", kTHnSparseF, {axisInvM, axisPt, axisCent});
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
    auto tofPIDp = static_cast<std::vector<float>>(protonTOFPIDp);
    auto tofPIDcut = static_cast<std::vector<float>>(protonTOFPIDcut);
    int nitr = static_cast<int>(tpcPIDp.size());
    int nitrtof = static_cast<int>(tofPIDp.size());

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
      if (candidate.tofNSigmaPr() < minnsigmatofPr)
        return false;
      if (nsigmaCutCombinedProton < 0 && p >= cPMin) {

        for (int i = 0; i < nitrtof - 1; ++i) {
          if (p >= tofPIDp[i] && p < tofPIDp[i + 1] && (tofNsigmaPr < tofPIDcut[i] && tofNsigmaPi > cRejNsigmaTof && tofNsigmaKa > cRejNsigmaTof))
            tofPIDPassed = true;
        }
        if (candidate.tpcNSigmaPr() < minnsigmatpctofPr)
          return false;
        if (tpcNsigmaPr < cMaxTPCnSigmaProton && tpcNsigmaPi > cRejNsigmaTpcVeto && tpcNsigmaKa > cRejNsigmaTpcVeto)
          tpcPIDPassed = true;
      }

      // circular cut
      if ((nsigmaCutCombinedProton > 0) && p >= cPMin && (tpcTofNsigmaPr < combinedCut && tpcTofNsigmaPi > combinedRejCut && tpcTofNsigmaKa > combinedRejCut)) {
        tofPIDPassed = true;
        tpcPIDPassed = true;
      }

      if (p < cPMin && tpcNsigmaPr < cMaxTPCnSigmaProton) {

        tofPIDPassed = true;
        tpcPIDPassed = true;
      }
    } else {
      tofPIDPassed = true;
      if (candidate.tpcNSigmaPr() < minnsigmatpcPr)
        return false;
      for (int i = 0; i < nitr - 1; ++i) {
        if (p >= tpcPIDp[i] && p < tpcPIDp[i + 1] && (tpcNsigmaPr < tpcPIDcut[i] && tpcNsigmaPi > cRejNsigmaTpcPi && tpcNsigmaKa > cRejNsigmaTpcKa)) {
          tpcPIDPassed = true;
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
    auto tofPIDp = static_cast<std::vector<float>>(kaonTOFPIDp);
    auto tofPIDcut = static_cast<std::vector<float>>(kaonTOFPIDcut);
    int nitr = static_cast<int>(tpcPIDp.size());
    int nitrtof = static_cast<int>(tofPIDp.size());

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
      if (candidate.tofNSigmaKa() < minnsigmatofKa)
        return false;
      if (nsigmaCutCombinedKaon < 0 && p >= cPMin) {

        for (int i = 0; i < nitrtof - 1; ++i) {
          if (p >= tofPIDp[i] && p < tofPIDp[i + 1] && (tofNsigmaKa < tofPIDcut[i] && tofNsigmaPi > cRejNsigmaTof && tofNsigmaPr > cRejNsigmaTof))
            tofPIDPassed = true;
        }
        if (candidate.tpcNSigmaKa() < minnsigmatpctofKa)
          return false;
        if (tpcNsigmaKa < cMaxTPCnSigmaKaon && tpcNsigmaPi > cRejNsigmaTpcVeto && tpcNsigmaPr > cRejNsigmaTpcVeto)
          tpcPIDPassed = true;
      }

      // circular
      if ((nsigmaCutCombinedKaon > 0) && p >= cPMin && (tpcTofNsigmaKa < combinedCut && tpcTofNsigmaPi > combinedRejCut && tpcTofNsigmaPr > combinedRejCut)) {
        tofPIDPassed = true;
        tpcPIDPassed = true;
      }

      if (p < cPMin && tpcNsigmaKa < cMaxTPCnSigmaKaon) {

        tofPIDPassed = true;
        tpcPIDPassed = true;
      }

    } else {
      tofPIDPassed = true;
      if (candidate.tpcNSigmaKa() < minnsigmatpcKa)
        return false;
      for (int i = 0; i < nitr - 1; ++i) {
        if (p >= tpcPIDp[i] && p < tpcPIDp[i + 1] && (tpcNsigmaKa < tpcPIDcut[i] && tpcNsigmaPi > cRejNsigmakTpcPi && tpcNsigmaPr > cRejNsigmakTpcPr)) {
          tpcPIDPassed = true;
        }
      }
    }
    if (tpcPIDPassed && tofPIDPassed) {
      return true;
    }
    return false;
  }

  template <bool mix, bool mc, typename trackType>
  void fillDataHistos(trackType const& trk1, trackType const& trk2, float mult, int occup = 100)
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
      //  LOGF(info, "eork 4 %d, %d  %d  ",trkPr.index(),trk1.size(),trkPr.size());
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
      if (isDeepAngle && TMath::ACos((trkPr.pt() * trkKa.pt() + _pzPr * _pzKa) / (p_ptot * k_ptot)) < cfgDeepAngle)
        continue;

      // Fill QA after track selection.
      if constexpr (!mix) {
        auto _ptPr = trkPr.pt();
        auto _tpcnsigmaPr = trkPr.tpcNSigmaPr();

        // Proton
        histos.fill(HIST("QAafter/Proton/hd_pr_pt"), _ptPr, mult);
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
        histos.fill(HIST("QAafter/Kaon/hd_ka_pt"), _ptKa, mult);
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
          histos.fill(HIST("QAafter/Kaon/h2d_Kapr_nsigma_tof_p"), k_ptot, trkKa.tofNSigmaPr());
          histos.fill(HIST("QAafter/Kaon/h2d_ka_nsigma_tof_vs_tpc"), _tpcnsigmaKa, _tofnsigmaKa);
        }
      }

      if (isonlyQC)
        continue;
      // Invariant mass reconstruction.
      p1.SetXYZM(_pxPr, _pyPr, _pzPr, MassProton);
      p2.SetXYZM(_pxKa, _pyKa, _pzKa, MassKaonCharged);
      p = p1 + p2;

      if (std::abs(p.Rapidity()) > 0.5)
        continue;

      auto _M = p.M();
      auto _pt = p.Pt();

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
          if (trkPr.sign() > 0)
            histos.fill(HIST("Analysis/h4d_lstar_invm_US_PM"), _M, _pt, mult, occup);
          else
            histos.fill(HIST("Analysis/h4d_lstar_invm_US_MP"), _M, _pt, mult, occup);
          if (doRotate) {
            float theta = rn->Uniform(1.56, 1.58);
            p1.RotateZ(theta);
            p = p1 + p2;
            if (std::abs(p.Rapidity()) < 0.5) {
              histos.fill(HIST("Analysis/h4d_lstar_invm_rot"), p.M(), p.Pt(), mult, occup);
            }
          }
        } else {
          if (trkPr.sign() > 0) {
            histos.fill(HIST("Analysis/h4d_lstar_invm_PP"), _M, _pt, mult, occup);
          } else {
            histos.fill(HIST("Analysis/h4d_lstar_invm_MM"), _M, _pt, mult, occup);
          }
        }
      }

      if constexpr (mix) {
        if (trkPr.sign() * trkKa.sign() < 0) {
          if (trkPr.sign() > 0)
            histos.fill(HIST("Analysis/h4d_lstar_invm_US_PM_mix"), _M, _pt, mult, occup);
          else
            histos.fill(HIST("Analysis/h4d_lstar_invm_US_MP_mix"), _M, _pt, mult, occup);
        } else {
          if (trkPr.sign() > 0)
            histos.fill(HIST("Analysis/h4d_lstar_invm_LS_PP_mix"), _M, _pt, mult, occup);
          else
            histos.fill(HIST("Analysis/h4d_lstar_invm_LS_MM_mix"), _M, _pt, mult, occup);
        }
      }

      if constexpr (mc) {
        if (trkPr.sign() * trkKa.sign() < 0) {
          if (std::abs(trkPr.pdgCode()) != 2212 || std::abs(trkKa.pdgCode()) != 321)
            continue;

          if (trkPr.motherId() != trkKa.motherId())
            continue;

          if (std::abs(trkPr.motherPDG()) != lambda1520id) // L* pdg_code = 3124
            continue;

          // MC histograms
          if (trkPr.motherPDG() > 0) {
            histos.fill(HIST("Analysis/h3d_rec_lstar_PM"), _M, _pt, mult);
          } else {
            histos.fill(HIST("Analysis/h3d_rec_lstar_MP"), _M, _pt, mult);
          }
        }
      }
    }
  }

  using resoCols = aod::ResoCollisions;
  using resoTracks = aod::ResoTracks;

  void processData(resoCols::iterator const& collision, resoTracks const& tracks)
  {

    // LOGF(info, " collisions: Index = %d %d", collision.globalIndex(),tracks.size());
    histos.fill(HIST("Event/h1d_ft0_mult_percentile"), collision.cent(), 100);
    fillDataHistos<false, false>(tracks, tracks, collision.cent());
  }

  PROCESS_SWITCH(lambdaAnalysis_pb, processData, "Process for Same Event Data", true);

  void processMC(resoCols::iterator const& collision,
                 soa::Join<aod::ResoTracks, aod::ResoMCTracks> const& tracks, aod::ResoMCParents const& resoParents)
  {

    auto mult = collision.cent();
    histos.fill(HIST("Event/h1d_ft0_mult_percentile"), mult);
    fillDataHistos<false, true>(tracks, tracks, mult);

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

    for (auto const& part : resoParents) {
      if (abs(part.pdgCode()) != lambda1520id) // // L* pdg_code = 3124
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

      TLorentzVector p4;
      p4.SetPxPyPzE(part.px(), part.py(), part.pz(), part.e());
      auto mass = p4.M();
      if (part.pdgCode() > 0)
        histos.fill(HIST("Analysis/h3d_gen_lstar_PM"), mass, part.pt(), mult);
      else
        histos.fill(HIST("Analysis/h3d_gen_lstar_MP"), mass, part.pt(), mult);
    }
  }
  PROCESS_SWITCH(lambdaAnalysis_pb, processMC, "Process Event for MC", false);

  using BinningType2 = ColumnBinningPolicy<aod::collision::PosZ, aod::resocollision::Cent>;

  void processMix(resoCols& collisions, resoTracks const& tracks)
  {

    LOGF(debug, "Event Mixing Started");

    BinningType2 binningPositions2{{cMixVtxBins, cMixMultBins}, true};
    auto tracksTuple = std::make_tuple(tracks);

    SameKindPair<resoCols, resoTracks, BinningType2> pairs{binningPositions2, cNumMixEv, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip
    for (auto& [c1, t1, c2, t2] : pairs) {
      // LOGF(info, "processMCMixedDerived: Mixed collisions : %d (%.3f, %.3f,%d), %d (%.3f, %.3f,%d)",c1.globalIndex(), c1.posZ(), c1.cent(),c1.mult(), c2.globalIndex(), c2.posZ(), c2.cent(),c2.mult());
      histos.fill(HIST("Event/mixing_vzVsmultpercentile"), c1.cent(), c1.posZ(), c1.evtPl());
      fillDataHistos<true, false>(t1, t2, c1.cent());
    }
  }

  PROCESS_SWITCH(lambdaAnalysis_pb, processMix, "Process for Mixed Events", false);

  Preslice<aod::ResoTrackDFs> perRColdf = aod::resodaughter::resoCollisionDFId;

  using resoColDFs = aod::ResoCollisionDFs;
  using resoTrackDFs = aod::ResoTrackDFs;

  void processDatadf(resoColDFs::iterator const& collision, resoTrackDFs const& tracks)
  {

    if (doprocessData)
      LOG(error) << "Disable processData() first!";
    auto _occup = 100;
    if (ConfEvtOccupancyInTimeRange)
      _occup = collision.trackOccupancyInTimeRange();

    // LOGF(info, "inside df collisions: Index = %d %d", collision.globalIndex(),tracks.size());
    histos.fill(HIST("Event/h1d_ft0_mult_percentile"), collision.cent(), _occup);
    fillDataHistos<false, false>(tracks, tracks, collision.cent(), _occup);
  }

  PROCESS_SWITCH(lambdaAnalysis_pb, processDatadf, "Process for data merged DF", false);

  using BinningTypeDF = ColumnBinningPolicy<aod::collision::PosZ, aod::resocollision::Cent>;
  void processMixDF(resoColDFs& collisions, resoTrackDFs const& tracks)
  {
    if (doprocessMix)
      LOG(fatal) << "Disable processMix() first!";
    LOGF(debug, "Event Mixing Started");

    BinningTypeDF binningPositions2{{cMixVtxBins, cMixMultBins}, true};
    auto tracksTuple = std::make_tuple(tracks);

    SameKindPair<resoColDFs, resoTrackDFs, BinningTypeDF> pairs{binningPositions2, cNumMixEv, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip
    for (auto& [c1, t1, c2, t2] : pairs) {
      auto _occup = 100;
      if (ConfEvtOccupancyInTimeRange)
        _occup = c1.trackOccupancyInTimeRange();

      // LOGF(info, "processMCMixedDerived: Mixed collisions : %d (%.3f, %.3f,%d), %d (%.3f, %.3f,%d)",c1.globalIndex(), c1.posZ(), c1.cent(),c1.mult(), c2.globalIndex(), c2.posZ(), c2.cent(),c2.mult());
      histos.fill(HIST("Event/mixing_vzVsmultpercentile"), c1.cent(), c1.posZ(), c1.evtPl());
      fillDataHistos<true, false>(t1, t2, c1.cent(), _occup);
    }
  }

  PROCESS_SWITCH(lambdaAnalysis_pb, processMixDF, "Process for merged DF  Mixed Events", false);

  using BinningTypeEP = ColumnBinningPolicy<aod::collision::PosZ, aod::resocollision::Cent, aod::resocollision::EvtPl>;
  void processMixepDF(resoColDFs& collisions, resoTrackDFs const& tracks)
  {
    if (doprocessMix || doprocessMixDF)
      LOG(fatal) << "Disable processMix() or processMixDF() first!";
    LOGF(debug, "Event Mixing Started");
    BinningTypeEP binningPositions2{{cMixVtxBins, cMixMultBins, cMixEPAngle}, true};
    auto tracksTuple = std::make_tuple(tracks);

    SameKindPair<resoColDFs, resoTrackDFs, BinningTypeEP> pairs{binningPositions2, cNumMixEv, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip
    for (auto& [c1, t1, c2, t2] : pairs) {
      //  LOGF(info, "processMCMixedDerived: Mixed collisions : %d (%.3f, %.3f,%.3f), %d (%.3f, %.3f, %.3f)",c1.globalIndex(), c1.posZ(), c1.cent(),c1.evtPl(), c2.globalIndex(), c2.posZ(), c2.cent(),c2.evtPl());
      histos.fill(HIST("Event/mixing_vzVsmultpercentile"), c1.cent(), c1.posZ(), c1.evtPl());
      fillDataHistos<true, false>(t1, t2, c1.cent());
    }
  }

  PROCESS_SWITCH(lambdaAnalysis_pb, processMixepDF, "Process for merged DF  Mixed Events", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<lambdaAnalysis_pb>(cfgc)};
}
