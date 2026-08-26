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

/// \file lambda1520pbpb.cxx
/// \brief Invariant-mass reconstruction and cosThetaStar polarization/spin-alignment
///        analysis (Helicity / Collins-Soper / Production / Beam-axis / Random frames)
///        of the Lambda(1520) resonance in Pb-Pb collisions.
///
/// \author Dukhishyam Mallick
/// \author Nasir Mehdi Malik <nasir.mehdi.malik@cern.ch>

#include "PWGLF/DataModel/LFResonanceTables.h"

#include "Common/Core/RecoDecay.h"

#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/GroupedCombinations.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TPDGCode.h>

#include <array>
#include <cmath>
#include <cstdlib>
#include <random>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

struct Lambda1520pbpb {
  SliceCache cache;
  Preslice<aod::ResoTracks> perRCol = aod::resodaughter::resoCollisionId;
  Preslice<aod::Tracks> perCollision = aod::track::collisionId;
  // Configurables.
  aod::ResoMCParents const* mResoParents = nullptr;
  Configurable<bool> cfgEvtOccupancyInTimeRange{"cfgEvtOccupancyInTimeRange", false, "occupancy selection true or false"};
  Configurable<int> nBinsPt{"nBinsPt", 100, "N bins in pT histogram"};
  Configurable<int> nBinsInvM{"nBinsInvM", 120, "N bins in InvMass histogram"};
  Configurable<int> lambda1520id{"lambda1520id", 3124, "pdg"};
  Configurable<bool> doRotate{"doRotate", true, "rotated inv mass spectra"};
  // Tracks
  Configurable<float> cPtMin{"cPtMin", 0.15, "Minimum Track pT"};
  Configurable<float> cPMin{"cPMin", 0., "Minimum Track p"};
  Configurable<float> cEtaCut{"cEtaCut", 0.8, "Pseudorapidity cut"};
  Configurable<float> cDcazMin{"cDcazMin", 1., "Minimum DCAz"};
  Configurable<float> cfgRapidityShift{"cfgRapidityShift", 0., " rapidity shift"};
  Configurable<float> cfgRapidityCut{"cfgRapidityCut", 0.5, "Rapidity window"};
  // TPC crossed rows (absolute)
  Configurable<int> cfgMinCrossedRows{"cfgMinCrossedRows", 70, "min TPC crossed rows"};
  Configurable<bool> cfgUseCrossedRows{"cfgUseCrossedRows", false, "apply crossed rows cut"};

  Configurable<int> cfgMinTPCcls{"cfgMinTPCcls", 70, "min TPC clusters found"};
  Configurable<bool> cfgUseTPCcls{"cfgUseTPCcls", false, "apply TPC clusters cut"};

  Configurable<std::vector<float>> cDcaPtBinsPr{"cDcaPtBinsPr", {0.0f, 0.5f, 1.0f, 2.0f, 3.0f, 5.0f, 1000.0f}, "Proton pT bin edges for DCAxy cut"};
  Configurable<std::vector<float>> cDcaXYBinsPr{"cDcaXYBinsPr", {0.020f, 0.015f, 0.010f, 0.007f, 0.005f, 0.004f}, "Proton max |DCAxy| per pT bin (cm)"};

  // Kaon DCAxy — pT binned
  Configurable<std::vector<float>> cDcaPtBinsKa{"cDcaPtBinsKa", {0.0f, 0.3f, 0.6f, 1.0f, 2.0f, 1000.0f}, "Kaon pT bin edges for DCAxy cut"};
  Configurable<std::vector<float>> cDcaXYBinsKa{"cDcaXYBinsKa", {0.025f, 0.018f, 0.012f, 0.008f, 0.004f}, "Kaon max |DCAxy| per pT bin (cm)"};

  Configurable<bool> isonlyQC{"isonlyQC", false, "only QC"};
  Configurable<bool> isDeepAngle{"isDeepAngle", false, "Deep Angle cut"};
  Configurable<double> cfgDeepAngle{"cfgDeepAngle", 0.04, "Deep Angle cut value"};
  Configurable<bool> cKinCuts{"cKinCuts", false, "Kinematic Cuts for p-K pair opening angle"};
  Configurable<bool> cPrimaryTrack{"cPrimaryTrack", true, "Primary track selection"};
  Configurable<bool> cGlobalWoDCATrack{"cGlobalWoDCATrack", true, "Global track selection without DCA"};
  Configurable<bool> cPVContributor{"cPVContributor", true, "PV contributor track selection"};

  // PID Selections
  Configurable<bool> cUseOnlyTOFTrackPr{"cUseOnlyTOFTrackPr", false, "Use only TOF track for PID selection"};
  Configurable<bool> cUseOnlyTOFTrackKa{"cUseOnlyTOFTrackKa", false, "Use only TOF track for PID selection"};
  Configurable<bool> cUseTpcOnly{"cUseTpcOnly", false, "Use TPC Only selection"};
  Configurable<float> cRejNsigmaTpc{"cRejNsigmaTpc", 3.0, "Reject tracks to improve purity of TPC PID"};
  Configurable<float> cRejNsigmaTpcPi{"cRejNsigmaTpcPi", 3.0, "Reject tracks to improve purity of TPC PID"};
  Configurable<float> cRejNsigmaTpcKa{"cRejNsigmaTpcKa", 3.0, "Reject tracks to improve purity of TPC PID"};
  Configurable<float> cRejNsigmakTpcPi{"cRejNsigmakTpcPi", 3.0, "Reject tracks to improve purity of TPC PID"};
  Configurable<float> cRejNsigmakTpcPr{"cRejNsigmakTpcPr", 3.0, "Reject tracks to improve purity of TPC PID"};
  Configurable<float> minnsigmatpcKa{"minnsigmatpcKa", -6.0, "Reject tracks to improve purity of TPC PID"};
  Configurable<float> minnsigmatpcPr{"minnsigmatpcPr", -6.0, "Reject tracks to improve purity of TPC PID"};
  Configurable<float> minnsigmatofKa{"minnsigmatofKa", -6.0, "Reject tracks to improve purity of TofPID"};
  Configurable<float> minnsigmatofPr{"minnsigmatofPr", -6.0, "Reject tracks to improve purity of Tof PID"};
  Configurable<float> minnsigmatpctofKa{"minnsigmatpctofKa", -6.0, "Reject tracks to improve purity of TPC PID"};
  Configurable<float> minnsigmatpctofPr{"minnsigmatpctofPr", -6.0, "Reject tracks to improve purity of TPC PID"};
  Configurable<float> cRejNsigmaTpcVeto{"cRejNsigmaTpcVeto", 3.0, "Reject tracks to improve purity of TPC PID"};
  Configurable<float> cRejNsigmaTof{"cRejNsigmaTof", 3.0, "Reject tracks to improve purity of TOF PID"};

  // Proton
  Configurable<double> cMaxTPCnSigmaProton{"cMaxTPCnSigmaProton", 3.0, "TPC nSigma cut for Proton"};
  Configurable<double> nsigmaCutCombinedProton{"nsigmaCutCombinedProton", 3.0, "Combined nSigma cut for Proton"};
  Configurable<std::vector<float>> protonTPCPIDp{"protonTPCPIDp", {0, 0.5, 0.7, 0.8}, "p dependent TPC cuts protons"};
  Configurable<std::vector<float>> protonTPCPIDcut{"protonTPCPIDcut", {5., 3.5, 2.5}, "TPC nsigma cuts protons"};
  Configurable<std::vector<float>> protonTOFPIDp{"protonTOFPIDp", {0., 999.}, "p dependent TOF cuts protons"};
  Configurable<std::vector<float>> protonTOFPIDcut{"protonTOFPIDcut", {3.0}, "TOF nsigma cuts protons"};

  // Kaon
  Configurable<double> cMaxTPCnSigmaKaon{"cMaxTPCnSigmaKaon", 3.0, "TPC nSigma cut for Kaon"};
  Configurable<double> nsigmaCutCombinedKaon{"nsigmaCutCombinedKaon", 3.0, "Combined nSigma cut for Kaon"};
  Configurable<std::vector<float>> kaonTPCPIDp{"kaonTPCPIDp", {0., 0.25, 0.3, 0.45}, "pT dependent TPC cuts kaons"};
  Configurable<std::vector<float>> kaonTPCPIDcut{"kaonTPCPIDcut", {6, 3.5, 2.5}, "TPC nsigma cuts kaons"};
  Configurable<std::vector<float>> kaonTOFPIDp{"kaonTOFPIDp", {0., 999.}, "p dependent TOF cuts kaons"};
  Configurable<std::vector<float>> kaonTOFPIDcut{"kaonTOFPIDcut", {3.0}, "TOF nsigma cuts kaons"};

  // Event Mixing.
  Configurable<int> cNumMixEv{"cNumMixEv", 20, "Number of Events to be mixed"};
  ConfigurableAxis cDcaZBins{"cDcaZBins", {VARIABLE_WIDTH, -1.2f, -1.0f, -0.9f, -0.8f, -0.7f, -0.6f, -0.5f, -0.4f, -0.3f, -0.2f, -0.1f, 0.f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.2f}, "DCA - z-vertex"};
  ConfigurableAxis cMixVtxBins{"cMixVtxBins", {VARIABLE_WIDTH, -10.0f, -9.f, -8.f, -7.f, -6.f, -5.f, -4.f, -3.f, -2.f, -1.f, 0.f, 1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis cMixMultBins{"cMixMultBins", {VARIABLE_WIDTH, 0.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.0f, 200.0f}, "Mixing bins - multiplicity"};
  ConfigurableAxis cMixEPAngle{"cMixEPAngle", {VARIABLE_WIDTH, -1.5708f, -1.25664f, -0.942478f, -0.628319f, 0.f, 0.628319f, 0.942478f, 1.25664f, 1.5708f}, "event plane"};
  ConfigurableAxis occupancyBins{"occupancyBins", {VARIABLE_WIDTH, 0.0, 100, 500, 600, 1000, 1100, 1500, 1600, 2000, 2100, 2500, 2600, 3000, 3100, 3500, 3600, 4000, 4100, 4500, 4600, 5000, 5100, 9999}, "Binning of the occupancy axis"};
  Configurable<int> cNofRotations{"cNofRotations", 10, "Number of rotations for rotational background"};
  Configurable<float> rotationalcut{"rotationalcut", 6.f, "Rotational background angle window: PI/rotationalcut"};

  // ── MC Event Selection Configurables ─────────────────────────────────────
  Configurable<bool> cEvtMCAfterAllCuts{"cEvtMCAfterAllCuts", false, "MC event sel: isInAfterAllCuts"};
  Configurable<bool> cEvtMCINELgt0{"cEvtMCINELgt0", false, "MC event sel: isINELgt0"};
  Configurable<bool> cEvtMCSel8{"cEvtMCSel8", false, "MC event sel: isInSel8"};
  Configurable<bool> cEvtMCVtxIn10{"cEvtMCVtxIn10", false, "MC event sel: isVtxIn10"};
  Configurable<bool> cEvtMCTriggerTVX{"cEvtMCTriggerTVX", false, "MC event sel: isTriggerTVX"};
  Configurable<bool> cEvtMCRecINELgt0{"cEvtMCRecINELgt0", false, "MC event sel: isRecINELgt0"};

  // ── Polarization / angular-distribution (cosThetaStar) Configurables ─────
  Configurable<int> cPolarizationRefDaughter{"cPolarizationRefDaughter", 0, "Reference daughter for cosThetaStar: 0 = proton, 1 = kaon"};
  Configurable<bool> cActivateHelicityFrame{"cActivateHelicityFrame", false, "Compute cosThetaStar in the Helicity frame"};
  Configurable<bool> cActivateCollinsSoperFrame{"cActivateCollinsSoperFrame", false, "Compute cosThetaStar in the Collins-Soper frame"};
  Configurable<bool> cActivateProductionFrame{"cActivateProductionFrame", false, "Compute cosThetaStar wrt the production-plane normal"};
  Configurable<bool> cActivateBeamAxisFrame{"cActivateBeamAxisFrame", false, "Compute cosThetaStar in the beam-axis (Gottfried-Jackson) frame"};
  Configurable<bool> cActivateRandomFrame{"cActivateRandomFrame", false, "Compute cosThetaStar wrt a random axis (null/control test, should be flat)"};
  Configurable<bool> doPolarizationRot{"doPolarizationRot", false, "Also fill rotated-background polarization histograms"};
  Configurable<bool> doPolarizationMix{"doPolarizationMix", false, "Also fill mixed-event-background polarization histograms"};
  Configurable<float> cfgBeamEnergy{"cfgBeamEnergy", 13600.0, "sqrt(s) beam energy in GeV (used for Collins-Soper/beam-axis/production frames)"};
  ConfigurableAxis cCosThetaBins{"cCosThetaBins", {20, -1.0f, 1.0f}, "cosThetaStar binning"};

  float beam1E = 0.f;
  std::array<float, 3> beam1P{};
  float beam2E = 0.f;
  std::array<float, 3> beam2P{};

  std::mt19937 randGen{std::random_device{}()};
  std::uniform_real_distribution<double> randUniform01{0., 1.};

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {
    const AxisSpec axisCent(110, 0, 110, "FT0 (%)");
    const AxisSpec axisPpid(200, 0., 10., "p (GeV/c)");
    const AxisSpec axisPtpid(200, 0., 10., "p_{T} (GeV/c)");
    const AxisSpec axisPt(nBinsPt, 0., 10., "p_{T} (GeV/c)");
    const AxisSpec axisEta(40, -1, 1, "#eta");
    const AxisSpec axisDCAxy(240, -0.12, 0.12, {"DCA_{xy} (cm)"});
    const AxisSpec axisTPCNsigma(401, -10.025, 10.025, {"n#sigma^{TPC}"});
    const AxisSpec axisTOFNsigma(401, -10.025, 10.025, {"n#sigma^{TOF}"});
    const AxisSpec axisdEdx(380, 10, 200, {"#frac{dE}{dx}"});
    const AxisSpec axisVz(120, -12, 12, {"vz"});
    const AxisSpec axisEP(120, -3.14, 3.14, {"#theta"});
    const AxisSpec axisInvM(nBinsInvM, 1.4, 2.0, {"M_{inv} (GeV/c^{2})"});
    AxisSpec axisOccupancy = {occupancyBins, "Occupancy [-40,100]"};
    AxisSpec axisDCAz = {cDcaZBins, "DCA_{z} (cm)"};

    double beamMomentum = std::sqrt(cfgBeamEnergy * cfgBeamEnergy / 4. - MassProton * MassProton);
    beam1E = static_cast<float>(cfgBeamEnergy * 0.5);
    beam1P = {0.f, 0.f, static_cast<float>(-beamMomentum)};
    beam2E = static_cast<float>(cfgBeamEnergy * 0.5);
    beam2P = {0.f, 0.f, static_cast<float>(beamMomentum)};

    histos.add("Event/h1d_ft0_mult_percentile", "FT0 (%)", kTH2F, {axisCent, axisOccupancy});
    histos.add("Event/h_ft0_vz", "Collision Vertex Z position", kTH1F, {{100, -15., 15.}});
    if (doprocessMix || doprocessMixDF || doprocessMixepDF) {
      histos.add("Event/mixing_vzVsmultpercentile", "FT0(%)", kTH3F, {axisCent, axisVz, axisEP});
    }

    // QA Before
    histos.add("QAbefore/hEta_rec", "Reco dN/d#eta; #eta; dN/d#eta", kTH1F, {{50, -1.0, 1.0}});
    histos.add("QAbefore/hPt_rec", "Reco pT; p_{T} (GeV/c); Tracks", kTH1F, {axisPpid});
    histos.add("QAbefore/hPhi_rec", "Reco #varphi; #varphi (rad); Tracks", kTH1F, {{72, 0, 6.2832}});
    histos.add("QAbefore/hEtaPhi_rec", "Reco #eta vs #varphi; #eta; #varphi", kTH2F, {axisEta, {72, 0, 6.2832}});
    histos.add("QAbefore/Proton/h2d_pr_nsigma_tpc_p", "n#sigma^{TPC} Protons", kTH2F, {axisPpid, axisTPCNsigma});
    histos.add("QAbefore/Proton/h2d_pr_nsigma_tof_p", "n#sigma^{TOF} Protons", kTH2F, {axisPpid, axisTOFNsigma});
    histos.add("QAbefore/Proton/h2d_pr_nsigma_tof_vs_tpc", "n#sigma^{TPC} vs n#sigma^{TOF} Protons", kTH2F, {axisTPCNsigma, axisTOFNsigma});
    histos.add("QAbefore/Kaon/h2d_ka_nsigma_tpc_p", "n#sigma^{TPC} Kaons", kTH2F, {axisPpid, axisTPCNsigma});
    histos.add("QAbefore/Kaon/h2d_ka_nsigma_tof_p", "n#sigma^{TOF} Kaons", kTH2F, {axisPpid, axisTOFNsigma});
    histos.add("QAbefore/Kaon/h2d_ka_nsigma_tof_vs_tpc", "n#sigma^{TPC} vs n#sigma^{TOF} Kaons", kTH2F, {axisTPCNsigma, axisTOFNsigma});

    // ========================================================================
    // QA After - Protons
    // ========================================================================
    // 1. Base QA (Filled for all tracks)
    histos.add("QAafter/Proton/h2d_pr_dca_z", "dca_{z} Protons", kTH2F, {axisPtpid, axisDCAz});
    histos.add("QAafter/Proton/h2d_pr_dca_xy", "dca_{xy} Protons", kTH2F, {axisPtpid, axisDCAxy});
    histos.add("QAafter/Proton/h2d_pr_dEdx_p", "TPC Signal Protons", kTH2F, {axisPpid, axisdEdx});
    histos.add("QAafter/Proton/hTPCNClsCrossedRowsVsPt", "TPC Crossed Rows vs pT;p_{T} (GeV/c);N_{cls,crossed};Counts", kTH2F, {axisPtpid, {200, 0, 200}});
    histos.add("QAafter/Proton/hTPCNClsFoundVsPt", "TPC Found Clusters vs pT;p_{T} (GeV/c);N_{cls,found};Counts", kTH2F, {axisPtpid, {200, 0, 200}});
    histos.add("QAafter/Proton/hFracTPConly_Proton", "Proton candidates using TPC-only PID;p (GeV/c);TPC-only (0=no,1=yes)", kTH2F, {axisPpid, {2, 0, 2}});

    // 2. TPC-Only Split
    histos.add("QAafter/Proton/hd_pr_pt_TPConly", "p_{T}-spectra Protons, TPC-only PID", kTH2F, {axisPtpid, axisCent});
    histos.add("QAafter/Proton/h2d_pr_nsigma_tpc_pt_TPConly", "Protons, TPC-only PID", kTH2F, {axisPtpid, axisTPCNsigma});
    histos.add("QAafter/Proton/h2d_Prpi_nsigma_tpc_p_TPConly", "Protons pion, TPC-only PID", kTH2F, {axisPtpid, axisTPCNsigma});
    histos.add("QAafter/Proton/h2d_Prka_nsigma_tpc_p_TPConly", "Protons kaon, TPC-only PID", kTH2F, {axisPtpid, axisTPCNsigma});
    histos.add("QAafter/Proton/h2d_pr_nsigma_tpc_p_TPConly", "Protons, TPC-only PID", kTH2F, {axisPpid, axisTPCNsigma});

    // 3. hasTOF Split
    histos.add("QAafter/Proton/hd_pr_pt_hasTOF", "p_{T}-spectra Protons, hasTOF PID", kTH2F, {axisPtpid, axisCent});
    histos.add("QAafter/Proton/h2d_pr_nsigma_tpc_pt_hasTOF", "Protons, hasTOF PID", kTH2F, {axisPtpid, axisTPCNsigma});
    histos.add("QAafter/Proton/h2d_Prpi_nsigma_tpc_p_hasTOF", "Protons pion, hasTOF PID", kTH2F, {axisPtpid, axisTPCNsigma});
    histos.add("QAafter/Proton/h2d_Prka_nsigma_tpc_p_hasTOF", "Protons kaon, hasTOF PID", kTH2F, {axisPtpid, axisTPCNsigma});
    histos.add("QAafter/Proton/h2d_pr_nsigma_tpc_p_hasTOF", "Protons, hasTOF PID", kTH2F, {axisPpid, axisTPCNsigma});

    // 4. TOF-Specific (Only possible if hasTOF)
    histos.add("QAafter/Proton/h2d_pr_nsigma_tof_pt", " Protons", kTH2F, {axisPtpid, axisTOFNsigma});
    histos.add("QAafter/Proton/h2d_pr_nsigma_tof_p", " Protons", kTH2F, {axisPpid, axisTOFNsigma});
    histos.add("QAafter/Proton/h2d_Prpi_nsigma_tof_p", " Protons pion", kTH2F, {axisPpid, axisTOFNsigma});
    histos.add("QAafter/Proton/h2d_Prka_nsigma_tof_p", " Protons kaon", kTH2F, {axisPpid, axisTOFNsigma});
    histos.add("QAafter/Proton/h2d_pr_nsigma_tof_vs_tpc", "n#sigma(TOF) vs n#sigma(TPC) Protons", kTH2F, {axisTPCNsigma, axisTOFNsigma});

    // ========================================================================
    // QA After - Kaons
    // ========================================================================
    // 1. Base QA (Filled for all tracks)
    histos.add("QAafter/Kaon/h2d_ka_dca_z", "dca_{z} Kaons", kTH2F, {axisPtpid, axisDCAz});
    histos.add("QAafter/Kaon/h2d_ka_dca_xy", "dca_{xy} Kaons", kTH2F, {axisPtpid, axisDCAxy});
    histos.add("QAafter/Kaon/h2d_ka_dEdx_p", "TPC Signal Kaon", kTH2F, {axisPpid, axisdEdx});
    histos.add("QAafter/Kaon/hTPCNClsCrossedRowsVsPt", "TPC Crossed Rows vs pT;p_{T} (GeV/c);N_{cls,crossed};Counts", kTH2F, {axisPtpid, {200, 0, 200}});
    histos.add("QAafter/Kaon/hTPCNClsFoundVsPt", "TPC Found Clusters vs pT;p_{T} (GeV/c);N_{cls,found};Counts", kTH2F, {axisPtpid, {200, 0, 200}});
    histos.add("QAafter/Kaon/hFracTPConly_Kaon", "Kaon candidates using TPC-only PID;p (GeV/c);TPC-only (0=no,1=yes)", kTH2F, {axisPpid, {2, 0, 2}});

    // 2. TPC-Only Split
    histos.add("QAafter/Kaon/hd_ka_pt_TPConly", "p_{T}-spectra Kaons, TPC-only PID", kTH2F, {axisPtpid, axisCent});
    histos.add("QAafter/Kaon/h2d_ka_nsigma_tpc_pt_TPConly", "Kaons, TPC-only PID", kTH2F, {axisPtpid, axisTPCNsigma});
    histos.add("QAafter/Kaon/h2d_Kapi_nsigma_tpc_p_TPConly", "Kaons pion, TPC-only PID", kTH2F, {axisPtpid, axisTPCNsigma});
    histos.add("QAafter/Kaon/h2d_Kapr_nsigma_tpc_p_TPConly", "Kaons proton, TPC-only PID", kTH2F, {axisPpid, axisTPCNsigma});
    histos.add("QAafter/Kaon/h2d_ka_nsigma_tpc_p_TPConly", "Kaons, TPC-only PID", kTH2F, {axisPpid, axisTPCNsigma});

    // 3. hasTOF Split
    histos.add("QAafter/Kaon/hd_ka_pt_hasTOF", "p_{T}-spectra Kaons, hasTOF PID", kTH2F, {axisPtpid, axisCent});
    histos.add("QAafter/Kaon/h2d_ka_nsigma_tpc_pt_hasTOF", "Kaons, hasTOF PID", kTH2F, {axisPtpid, axisTPCNsigma});
    histos.add("QAafter/Kaon/h2d_Kapi_nsigma_tpc_p_hasTOF", "Kaons pion, hasTOF PID", kTH2F, {axisPtpid, axisTPCNsigma});
    histos.add("QAafter/Kaon/h2d_Kapr_nsigma_tpc_p_hasTOF", "Kaons proton, hasTOF PID", kTH2F, {axisPpid, axisTPCNsigma});
    histos.add("QAafter/Kaon/h2d_ka_nsigma_tpc_p_hasTOF", "Kaons, hasTOF PID", kTH2F, {axisPpid, axisTPCNsigma});

    // 4. TOF-Specific (Only possible if hasTOF)
    histos.add("QAafter/Kaon/h2d_ka_nsigma_tof_pt", " Kaons", kTH2F, {axisPtpid, axisTOFNsigma});
    histos.add("QAafter/Kaon/h2d_ka_nsigma_tof_p", " Kaons", kTH2F, {axisPpid, axisTOFNsigma});
    histos.add("QAafter/Kaon/h2d_Kapi_nsigma_tof_p", " Kaons pion", kTH2F, {axisPpid, axisTOFNsigma});
    histos.add("QAafter/Kaon/h2d_Kapr_nsigma_tof_p", " Kaons proton", kTH2F, {axisPpid, axisTOFNsigma});
    histos.add("QAafter/Kaon/h2d_ka_nsigma_tof_vs_tpc", "n#sigma(TOF) vs n#sigma(TPC) Kaons", kTH2F, {axisTPCNsigma, axisTOFNsigma});

    // Analysis
    if (!doprocessMC) {
      int nActiveFrames = static_cast<int>(cActivateHelicityFrame) +
                          static_cast<int>(cActivateCollinsSoperFrame) +
                          static_cast<int>(cActivateProductionFrame) +
                          static_cast<int>(cActivateBeamAxisFrame) +
                          static_cast<int>(cActivateRandomFrame);

      bool polActive = nActiveFrames > 0;

      if (nActiveFrames > 1) {
        LOG(fatal) << "Multiple polarization frames activated, but axisPolFrame was removed! "
                   << "Please activate only ONE frame at a time to prevent data mixing in the 4D sparse.";
      }

      if (!polActive) {
        histos.add("Analysis/h4d_lstar_invm_US_PM", "THn #Lambda(1520)", kTHnSparseF, {axisInvM, axisPt, axisCent, axisOccupancy});
        histos.add("Analysis/h4d_lstar_invm_US_MP", "THn #bar #Lambda(1520)", kTHnSparseF, {axisInvM, axisPt, axisCent, axisOccupancy});
        histos.add("Analysis/h4d_lstar_invm_PP", "THn Like Signs p K^{+}", kTHnSparseF, {axisInvM, axisPt, axisCent, axisOccupancy});
        histos.add("Analysis/h4d_lstar_invm_MM", "THn Like Signs #bar{p} K^{-}", kTHnSparseF, {axisInvM, axisPt, axisCent, axisOccupancy});

        if (doRotate) {
          histos.add("Analysis/h4d_lstar_invm_rot_PM", "THn Rotated", kTHnSparseF, {axisInvM, axisPt, axisCent, axisOccupancy});
          histos.add("Analysis/h4d_lstar_invm_rot_MP", "THn Rotated", kTHnSparseF, {axisInvM, axisPt, axisCent, axisOccupancy});
        }

        if (doprocessMix || doprocessMixDF || doprocessMixepDF) {
          histos.add("Analysis/h4d_lstar_invm_US_PM_mix", "THn Mixed Events", kTHnSparseF, {axisInvM, axisPt, axisCent, axisOccupancy});
          histos.add("Analysis/h4d_lstar_invm_US_MP_mix", "THn anti Mixed Events", kTHnSparseF, {axisInvM, axisPt, axisCent, axisOccupancy});
          histos.add("Analysis/h4d_lstar_invm_LS_PP_mix", "THn Mixed Events PP", kTHnSparseF, {axisInvM, axisPt, axisCent, axisOccupancy});
          histos.add("Analysis/h4d_lstar_invm_LS_MM_mix", "THn Mixed Events MM", kTHnSparseF, {axisInvM, axisPt, axisCent, axisOccupancy});
        }
      }

      // ── 4D Polarization / cosThetaStar histograms ──────────────────────────
      if (polActive) {
        const AxisSpec axisCosTheta{cCosThetaBins, "cos#theta*"};

        histos.add("Analysis/h4d_pol_US_PM", "cos#theta* #Lambda(1520)", kTHnSparseF, {axisInvM, axisPt, axisCent, axisCosTheta});
        histos.add("Analysis/h4d_pol_US_MP", "cos#theta* #bar{#Lambda}(1520)", kTHnSparseF, {axisInvM, axisPt, axisCent, axisCosTheta});
        histos.add("Analysis/h4d_pol_LS_PP", "cos#theta* Like Signs p K^{+}", kTHnSparseF, {axisInvM, axisPt, axisCent, axisCosTheta});
        histos.add("Analysis/h4d_pol_LS_MM", "cos#theta* Like Signs #bar{p} K^{-}", kTHnSparseF, {axisInvM, axisPt, axisCent, axisCosTheta});

        if (doPolarizationRot) {
          histos.add("Analysis/h4d_pol_rot_PM", "cos#theta* Rotated bkg", kTHnSparseF, {axisInvM, axisPt, axisCent, axisCosTheta});
          histos.add("Analysis/h4d_pol_rot_MP", "cos#theta* Rotated bkg", kTHnSparseF, {axisInvM, axisPt, axisCent, axisCosTheta});
        }

        if (doPolarizationMix && (doprocessMix || doprocessMixDF || doprocessMixepDF)) {
          histos.add("Analysis/h4d_pol_mix_PM", "cos#theta* Mixed-event bkg", kTHnSparseF, {axisInvM, axisPt, axisCent, axisCosTheta});
          histos.add("Analysis/h4d_pol_mix_MP", "cos#theta* Mixed-event bkg", kTHnSparseF, {axisInvM, axisPt, axisCent, axisCosTheta});
          histos.add("Analysis/h4d_pol_LS_PP_mix", "cos#theta* Mixed PP", kTHnSparseF, {axisInvM, axisPt, axisCent, axisCosTheta});
          histos.add("Analysis/h4d_pol_LS_MM_mix", "cos#theta* Mixed MM", kTHnSparseF, {axisInvM, axisPt, axisCent, axisCosTheta});
        }
      }
    }

    // MC
    if (doprocessMC) {
      histos.add("Event/hMCEventCutflow", "MC Event Cutflow", kTH1F, {{7, 0, 7}});
      histos.add("QAChecks/h1d_pr_rec_pt", "Reconstructed p_{T}-spectra Protons", kTH1F, {axisPtpid});
      histos.add("QAChecks/h1d_ka_rec_pt", "Recondstucted p_{T}-spectra Kaons", kTH1F, {axisPtpid});
      histos.add("QAChecks/h1d_pr_gen_pt", "Generated p_{T}-spectra Protons", kTH1F, {axisPtpid});
      histos.add("QAChecks/h1d_ka_gen_pt", "Generated p_{T}-spectra Kaons", kTH1F, {axisPtpid});

      histos.add("Analysis/h3d_gen_lstar_PM", "Generated #Lambda(1520) p_{T}", kTHnSparseF, {axisInvM, axisPt, axisCent});
      histos.add("Analysis/h3d_gen_lstar_MP", "Generated #bar{#Lambda}(1520) p_{T}", kTHnSparseF, {axisInvM, axisPt, axisCent});
      histos.add("Analysis/h3d_rec_lstar_PM", "Reconstructed #Lambda(1520) p_{T}", kTHnSparseF, {axisInvM, axisPt, axisCent});
      histos.add("Analysis/h3d_rec_lstar_MP", "Reconstructed #bar{#Lambda}(1520) p_{T}", kTHnSparseF, {axisInvM, axisPt, axisCent});
      histos.add("Analysis/h3d_reso_lstar_PM", "Resolution #Lambda(1520) p_{T}", kTHnSparseF, {{200, -0.05, 0.05}, axisPt, axisCent});
      histos.add("Analysis/h3d_reso_lstar_MP", "Resolution #bar{#Lambda}(1520) p_{T}", kTHnSparseF, {{200, -0.05, 0.05}, axisPt, axisCent});
    }

    if (doprocessMCGen) {
      histos.add("SignalLoss/hMCEventCutflow", "MC Event Cutflow", kTH1F, {{7, 0, 7}});
      histos.add("SignalLoss/hGen_mT_scaled_Proton", "mT Scaled #Lambda(1520) from Proton", kTHnSparseF, {axisPt, axisCent});
      histos.add("SignalLoss/hGen_mT_scaled_AntiProton", "mT Scaled #bar{#Lambda}(1520) from AntiProton", kTHnSparseF, {axisPt, axisCent});
      histos.add("SignalLoss/hGen_mT_scaled_Lambda0", "mT Scaled #Lambda(1520) from Lambda0", kTHnSparseF, {axisPt, axisCent});
      histos.add("SignalLoss/hGen_mT_scaled_AntiLambda0", "mT Scaled #bar{#Lambda}(1520) from AntiLambda0", kTHnSparseF, {axisPt, axisCent});
      histos.add("SignalLoss/hGen_mT_scaled_XiMinus", "mT Scaled #Lambda(1520) from Xi-", kTHnSparseF, {axisPt, axisCent});
      histos.add("SignalLoss/hGen_mT_scaled_XiPlus", "mT Scaled #bar{#Lambda}(1520) from Xi+", kTHnSparseF, {axisPt, axisCent});
      histos.add("SignalLoss/hGen_mT_scaled_Xi0", "mT Scaled #Lambda(1520) from Xi0", kTHnSparseF, {axisPt, axisCent});
      histos.add("SignalLoss/hGen_mT_scaled_AntiXi0", "mT Scaled #bar{#Lambda}(1520) from AntiXi0", kTHnSparseF, {axisPt, axisCent});
      histos.add("SignalLoss/hGen_mT_scaled_OmegaMinus", "mT Scaled #Lambda(1520) from Omega-", kTHnSparseF, {axisPt, axisCent});
      histos.add("SignalLoss/hGen_mT_scaled_OmegaPlus", "mT Scaled #bar{#Lambda}(1520) from Omega+", kTHnSparseF, {axisPt, axisCent});
    }
  }

  template <typename T>
  bool selTracks(T const& track)
  {
    if (track.pt() < cPtMin)
      return false;
    if (std::abs(track.eta()) > cEtaCut)
      return false;
    if (cPrimaryTrack && !track.isPrimaryTrack())
      return false;
    if (cGlobalWoDCATrack && !track.isGlobalTrackWoDCA())
      return false;
    if (cPVContributor && !track.isPVContributor())
      return false;
    if (cfgUseCrossedRows && track.tpcNClsCrossedRows() < cfgMinCrossedRows)
      return false;
    if (cfgUseTPCcls && track.tpcNClsFound() < cfgMinTPCcls)
      return false;
    return true;
  }

  template <typename T>
  bool dcaSelectionProton(T const& track, float p)
  {
    auto ptBinsPr = static_cast<std::vector<float>>(cDcaPtBinsPr);
    auto dcaXYPr = static_cast<std::vector<float>>(cDcaXYBinsPr);
    int nBinsPr = static_cast<int>(ptBinsPr.size()) - 1;

    bool dcaXYPassed = false;
    for (int i = 0; i < nBinsPr; i++) {
      if (p >= ptBinsPr[i] && p < ptBinsPr[i + 1] &&
          std::abs(track.dcaXY()) < dcaXYPr[i])
        dcaXYPassed = true;
    }
    if (!dcaXYPassed)
      return false;
    if (std::abs(track.dcaZ()) > cDcazMin)
      return false;
    return true;
  }

  template <typename T>
  bool dcaSelectionKaon(T const& track, float p)
  {
    auto ptBinsKa = static_cast<std::vector<float>>(cDcaPtBinsKa);
    auto dcaXYKa = static_cast<std::vector<float>>(cDcaXYBinsKa);
    int nBinsKa = static_cast<int>(ptBinsKa.size()) - 1;

    bool dcaXYPassed = false;
    for (int i = 0; i < nBinsKa; i++) {
      if (p >= ptBinsKa[i] && p < ptBinsKa[i + 1] &&
          std::abs(track.dcaXY()) < dcaXYKa[i])
        dcaXYPassed = true;
    }
    if (!dcaXYPassed)
      return false;
    if (std::abs(track.dcaZ()) > cDcazMin)
      return false;
    return true;
  }

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
    if (tpcPIDPassed && tofPIDPassed)
      return true;
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
    if (tpcPIDPassed && tofPIDPassed)
      return true;
    return false;
  }

  static constexpr float VecUnitEpsilon = 1e-9f;
  static constexpr float MomentumEpsilon = 1e-6f;

  static float vecDot(const std::array<float, 3>& a, const std::array<float, 3>& b)
  {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  }

  static std::array<float, 3> vecSub(const std::array<float, 3>& a, const std::array<float, 3>& b)
  {
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
  }

  static std::array<float, 3> vecCross(const std::array<float, 3>& a, const std::array<float, 3>& b)
  {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
  }

  static std::array<float, 3> vecUnit(const std::array<float, 3>& a)
  {
    float mag = std::sqrt(vecDot(a, a));
    if (mag < VecUnitEpsilon)
      return {0.f, 0.f, 0.f};
    return {a[0] / mag, a[1] / mag, a[2] / mag};
  }

  static std::array<float, 3> boostToRestFrame(float motherE, const std::array<float, 3>& motherP,
                                               float dE, const std::array<float, 3>& dP)
  {
    float pMag2 = vecDot(motherP, motherP);
    float pMag = std::sqrt(pMag2);
    if (pMag < MomentumEpsilon)
      return dP;

    std::array<float, 3> nHat = {motherP[0] / pMag, motherP[1] / pMag, motherP[2] / pMag};
    float beta = pMag / motherE;
    float motherMass2 = motherE * motherE - pMag2;
    float gamma = motherMass2 > 0.f ? motherE / std::sqrt(motherMass2) : 1.f;

    float pParallel = vecDot(dP, nHat);
    std::array<float, 3> pPerp = {dP[0] - pParallel * nHat[0], dP[1] - pParallel * nHat[1], dP[2] - pParallel * nHat[2]};
    float pParallelPrime = gamma * (pParallel - beta * dE);

    return {pPerp[0] + pParallelPrime * nHat[0], pPerp[1] + pParallelPrime * nHat[1], pPerp[2] + pParallelPrime * nHat[2]};
  }

  enum class PolBkgMode { Signal,
                          Rotated,
                          Mixed,
                          LikeSign,
                          LikeSignMixed };

  template <PolBkgMode Mode>
  void fillPolarization(float candMass, float candPt,
                        const std::array<float, 3>& motherP, float motherE,
                        const std::array<float, 3>& protonP, float protonE,
                        const std::array<float, 3>& kaonP, float kaonE,
                        float mult, bool protonIsPositive)
  {
    if (!cActivateHelicityFrame && !cActivateCollinsSoperFrame && !cActivateProductionFrame &&
        !cActivateBeamAxisFrame && !cActivateRandomFrame)
      return;

    std::array<float, 3> protonCM = boostToRestFrame(motherE, motherP, protonE, protonP);
    std::array<float, 3> kaonCM = boostToRestFrame(motherE, motherP, kaonE, kaonP);
    std::array<float, 3> refCM = vecUnit(cPolarizationRefDaughter == 0 ? protonCM : kaonCM);

    auto fillFrame = [&](float cosTheta) {
      if constexpr (Mode == PolBkgMode::Rotated) {
        if (protonIsPositive)
          histos.fill(HIST("Analysis/h4d_pol_rot_PM"), candMass, candPt, mult, cosTheta);
        else
          histos.fill(HIST("Analysis/h4d_pol_rot_MP"), candMass, candPt, mult, cosTheta);
      } else if constexpr (Mode == PolBkgMode::Mixed) {
        if (protonIsPositive)
          histos.fill(HIST("Analysis/h4d_pol_mix_PM"), candMass, candPt, mult, cosTheta);
        else
          histos.fill(HIST("Analysis/h4d_pol_mix_MP"), candMass, candPt, mult, cosTheta);
      } else if constexpr (Mode == PolBkgMode::LikeSign) {
        if (protonIsPositive)
          histos.fill(HIST("Analysis/h4d_pol_LS_PP"), candMass, candPt, mult, cosTheta);
        else
          histos.fill(HIST("Analysis/h4d_pol_LS_MM"), candMass, candPt, mult, cosTheta);
      } else if constexpr (Mode == PolBkgMode::LikeSignMixed) {
        if (protonIsPositive)
          histos.fill(HIST("Analysis/h4d_pol_LS_PP_mix"), candMass, candPt, mult, cosTheta);
        else
          histos.fill(HIST("Analysis/h4d_pol_LS_MM_mix"), candMass, candPt, mult, cosTheta);
      } else {
        if (protonIsPositive)
          histos.fill(HIST("Analysis/h4d_pol_US_PM"), candMass, candPt, mult, cosTheta);
        else
          histos.fill(HIST("Analysis/h4d_pol_US_MP"), candMass, candPt, mult, cosTheta);
      }
    };

    if (cActivateHelicityFrame) {
      std::array<float, 3> zAxisHE = vecUnit(motherP);
      fillFrame(vecDot(refCM, zAxisHE));
    }

    if (cActivateCollinsSoperFrame || cActivateProductionFrame || cActivateBeamAxisFrame) {
      std::array<float, 3> beam1CM = vecUnit(boostToRestFrame(motherE, motherP, beam1E, beam1P));
      std::array<float, 3> beam2CM = vecUnit(boostToRestFrame(motherE, motherP, beam2E, beam2P));

      if (cActivateCollinsSoperFrame) {
        std::array<float, 3> zAxisCS = vecUnit(vecSub(beam1CM, beam2CM));
        fillFrame(vecDot(refCM, zAxisCS));
      }
      if (cActivateProductionFrame) {
        std::array<float, 3> motherDir = vecUnit(motherP);
        std::array<float, 3> yAxisPR = vecUnit(vecCross(beam1CM, motherDir));
        fillFrame(vecDot(refCM, yAxisPR));
      }
      if (cActivateBeamAxisFrame) {
        fillFrame(vecDot(refCM, beam1CM));
      }
    }

    if (cActivateRandomFrame) {
      double randCos = randUniform01(randGen) * 2. - 1.;
      double randPhi = randUniform01(randGen) * o2::constants::math::TwoPI;
      double randSin = std::sqrt(1. - randCos * randCos);
      std::array<float, 3> zAxisRD = {static_cast<float>(randSin * std::cos(randPhi)),
                                      static_cast<float>(randSin * std::sin(randPhi)),
                                      static_cast<float>(randCos)};
      fillFrame(vecDot(refCM, zAxisRD));
    }
  }

  template <bool mix, bool mc, typename trackType>
  void fillDataHistos(trackType const& trk1, trackType const& trk2, float mult, int occup = 100)
  {
    float prptot = 0., kaptot = 0.;
    bool polActive = static_cast<bool>(cActivateHelicityFrame) ||
                     static_cast<bool>(cActivateCollinsSoperFrame) ||
                     static_cast<bool>(cActivateProductionFrame) ||
                     static_cast<bool>(cActivateBeamAxisFrame) ||
                     static_cast<bool>(cActivateRandomFrame);

    for (auto const& [trkPr, trkKa] : soa::combinations(soa::CombinationsFullIndexPolicy(trk1, trk2))) {
      if (trkPr.index() == trkKa.index())
        continue;
      if (!selTracks(trkPr) || !selTracks(trkKa))
        continue;

      auto pxPr = trkPr.px();
      auto pyPr = trkPr.py();
      auto pzPr = trkPr.pz();
      auto pxKa = trkKa.px();
      auto pyKa = trkKa.py();
      auto pzKa = trkKa.pz();

      prptot = std::sqrt(pxPr * pxPr + pyPr * pyPr + pzPr * pzPr);
      kaptot = std::sqrt(pxKa * pxKa + pyKa * pyKa + pzKa * pzKa);

      if (!mix) {
        auto tpcNsigmaPr = trkPr.tpcNSigmaPr();
        histos.fill(HIST("QAbefore/Proton/h2d_pr_nsigma_tpc_p"), prptot, tpcNsigmaPr);
        if (trkPr.hasTOF()) {
          auto tofNsigmaPr = trkPr.tofNSigmaPr();
          histos.fill(HIST("QAbefore/Proton/h2d_pr_nsigma_tof_p"), prptot, tofNsigmaPr);
          histos.fill(HIST("QAbefore/Proton/h2d_pr_nsigma_tof_vs_tpc"), tpcNsigmaPr, tofNsigmaPr);
        }
        auto tpcNsigmaKa = trkKa.tpcNSigmaKa();
        histos.fill(HIST("QAbefore/Kaon/h2d_ka_nsigma_tpc_p"), kaptot, tpcNsigmaKa);
        if (trkKa.hasTOF()) {
          auto tofNsigmaKa = trkKa.tofNSigmaKa();
          histos.fill(HIST("QAbefore/Kaon/h2d_ka_nsigma_tof_p"), kaptot, tofNsigmaKa);
          histos.fill(HIST("QAbefore/Kaon/h2d_ka_nsigma_tof_vs_tpc"), tpcNsigmaKa, tofNsigmaKa);
        }
      }

      if (cUseOnlyTOFTrackPr && !trkPr.hasTOF())
        continue;
      if (cUseOnlyTOFTrackKa && !trkKa.hasTOF())
        continue;
      if (!selectionPIDProton(trkPr, prptot) || !selectionPIDKaon(trkKa, kaptot))
        continue;
      if (!dcaSelectionProton(trkPr, prptot) || !dcaSelectionKaon(trkKa, kaptot))
        continue;

      if (isDeepAngle && std::acos((trkPr.pt() * trkKa.pt() + pzPr * pzKa) / (prptot * kaptot)) < cfgDeepAngle)
        continue;

      if constexpr (!mix) {
        // ====================================================================
        // Proton QA Fills
        // ====================================================================
        auto ptPr = trkPr.pt();
        auto tpcNsigmaPr = trkPr.tpcNSigmaPr();

        // 1. Base QA (Filled for all)
        histos.fill(HIST("QAafter/Proton/h2d_pr_dca_z"), ptPr, trkPr.dcaZ());
        histos.fill(HIST("QAafter/Proton/h2d_pr_dca_xy"), ptPr, trkPr.dcaXY());
        histos.fill(HIST("QAafter/Proton/h2d_pr_dEdx_p"), prptot, trkPr.tpcSignal());
        histos.fill(HIST("QAafter/Proton/hTPCNClsCrossedRowsVsPt"), ptPr, trkPr.tpcNClsCrossedRows());
        histos.fill(HIST("QAafter/Proton/hTPCNClsFoundVsPt"), ptPr, trkPr.tpcNClsFound());

        bool prIsTpcOnly = cUseTpcOnly || !trkPr.hasTOF();
        histos.fill(HIST("QAafter/Proton/hFracTPConly_Proton"), prptot, prIsTpcOnly ? 1 : 0);

        // 2 & 3. Mutually Exclusive TPC vs hasTOF Split
        if (prIsTpcOnly) {
          histos.fill(HIST("QAafter/Proton/hd_pr_pt_TPConly"), ptPr, mult);
          histos.fill(HIST("QAafter/Proton/h2d_pr_nsigma_tpc_pt_TPConly"), ptPr, tpcNsigmaPr);
          histos.fill(HIST("QAafter/Proton/h2d_Prpi_nsigma_tpc_p_TPConly"), ptPr, trkPr.tpcNSigmaPi());
          histos.fill(HIST("QAafter/Proton/h2d_Prka_nsigma_tpc_p_TPConly"), ptPr, trkPr.tpcNSigmaKa());
          histos.fill(HIST("QAafter/Proton/h2d_pr_nsigma_tpc_p_TPConly"), prptot, tpcNsigmaPr);
        } else {
          histos.fill(HIST("QAafter/Proton/hd_pr_pt_hasTOF"), ptPr, mult);
          histos.fill(HIST("QAafter/Proton/h2d_pr_nsigma_tpc_pt_hasTOF"), ptPr, tpcNsigmaPr);
          histos.fill(HIST("QAafter/Proton/h2d_Prpi_nsigma_tpc_p_hasTOF"), ptPr, trkPr.tpcNSigmaPi());
          histos.fill(HIST("QAafter/Proton/h2d_Prka_nsigma_tpc_p_hasTOF"), ptPr, trkPr.tpcNSigmaKa());
          histos.fill(HIST("QAafter/Proton/h2d_pr_nsigma_tpc_p_hasTOF"), prptot, tpcNsigmaPr);

          // 4. TOF-Specific (Only possible if hasTOF)
          auto tofNsigmaPr = trkPr.tofNSigmaPr();
          histos.fill(HIST("QAafter/Proton/h2d_pr_nsigma_tof_pt"), ptPr, tofNsigmaPr);
          histos.fill(HIST("QAafter/Proton/h2d_pr_nsigma_tof_p"), prptot, tofNsigmaPr);
          histos.fill(HIST("QAafter/Proton/h2d_Prpi_nsigma_tof_p"), prptot, trkPr.tofNSigmaPi());
          histos.fill(HIST("QAafter/Proton/h2d_Prka_nsigma_tof_p"), prptot, trkPr.tofNSigmaKa());
          histos.fill(HIST("QAafter/Proton/h2d_pr_nsigma_tof_vs_tpc"), tpcNsigmaPr, tofNsigmaPr);
        }

        // ====================================================================
        // Kaon QA Fills
        // ====================================================================
        auto ptKa = trkKa.pt();
        auto tpcNsigmaKa = trkKa.tpcNSigmaKa();

        // 1. Base QA (Filled for all)
        histos.fill(HIST("QAafter/Kaon/h2d_ka_dca_z"), ptKa, trkKa.dcaZ());
        histos.fill(HIST("QAafter/Kaon/h2d_ka_dca_xy"), ptKa, trkKa.dcaXY());
        histos.fill(HIST("QAafter/Kaon/h2d_ka_dEdx_p"), kaptot, trkKa.tpcSignal());
        histos.fill(HIST("QAafter/Kaon/hTPCNClsCrossedRowsVsPt"), ptKa, trkKa.tpcNClsCrossedRows());
        histos.fill(HIST("QAafter/Kaon/hTPCNClsFoundVsPt"), ptKa, trkKa.tpcNClsFound());

        bool kaIsTpcOnly = cUseTpcOnly || !trkKa.hasTOF();
        histos.fill(HIST("QAafter/Kaon/hFracTPConly_Kaon"), kaptot, kaIsTpcOnly ? 1 : 0);

        // 2 & 3. Mutually Exclusive TPC vs hasTOF Split
        if (kaIsTpcOnly) {
          histos.fill(HIST("QAafter/Kaon/hd_ka_pt_TPConly"), ptKa, mult);
          histos.fill(HIST("QAafter/Kaon/h2d_ka_nsigma_tpc_pt_TPConly"), ptKa, tpcNsigmaKa);
          histos.fill(HIST("QAafter/Kaon/h2d_Kapi_nsigma_tpc_p_TPConly"), ptKa, trkKa.tpcNSigmaPi());
          histos.fill(HIST("QAafter/Kaon/h2d_Kapr_nsigma_tpc_p_TPConly"), kaptot, trkKa.tpcNSigmaPr());
          histos.fill(HIST("QAafter/Kaon/h2d_ka_nsigma_tpc_p_TPConly"), kaptot, tpcNsigmaKa);
        } else {
          histos.fill(HIST("QAafter/Kaon/hd_ka_pt_hasTOF"), ptKa, mult);
          histos.fill(HIST("QAafter/Kaon/h2d_ka_nsigma_tpc_pt_hasTOF"), ptKa, tpcNsigmaKa);
          histos.fill(HIST("QAafter/Kaon/h2d_Kapi_nsigma_tpc_p_hasTOF"), ptKa, trkKa.tpcNSigmaPi());
          histos.fill(HIST("QAafter/Kaon/h2d_Kapr_nsigma_tpc_p_hasTOF"), kaptot, trkKa.tpcNSigmaPr());
          histos.fill(HIST("QAafter/Kaon/h2d_ka_nsigma_tpc_p_hasTOF"), kaptot, tpcNsigmaKa);

          // 4. TOF-Specific (Only possible if hasTOF)
          auto tofNsigmaKa = trkKa.tofNSigmaKa();
          histos.fill(HIST("QAafter/Kaon/h2d_ka_nsigma_tof_pt"), ptKa, tofNsigmaKa);
          histos.fill(HIST("QAafter/Kaon/h2d_ka_nsigma_tof_p"), kaptot, tofNsigmaKa);
          histos.fill(HIST("QAafter/Kaon/h2d_Kapi_nsigma_tof_p"), kaptot, trkKa.tofNSigmaPi());
          histos.fill(HIST("QAafter/Kaon/h2d_Kapr_nsigma_tof_p"), kaptot, trkKa.tofNSigmaPr());
          histos.fill(HIST("QAafter/Kaon/h2d_ka_nsigma_tof_vs_tpc"), tpcNsigmaKa, tofNsigmaKa);
        }
      }

      if (isonlyQC)
        continue;

      std::array<float, 3> pvec0 = {pxPr, pyPr, pzPr};
      std::array<float, 3> pvec1 = {pxKa, pyKa, pzKa};
      std::array<std::array<float, 3>, 2> arrMomrec = {pvec0, pvec1};
      float candMass = RecoDecay::m(arrMomrec, std::array{MassProton, MassKaonCharged});
      float candPt = RecoDecay::pt(std::array{pxPr + pxKa, pyPr + pyKa});
      float candY = std::abs(RecoDecay::y(std::array{pxPr + pxKa, pyPr + pyKa, pzPr + pzKa}, candMass));
      float candYShift = candY - cfgRapidityShift;

      if (std::abs(candYShift) > cfgRapidityCut)
        continue;

      if constexpr (!mix && !mc) {
        if (trkPr.sign() * trkKa.sign() < 0) {
          if (!polActive) {
            if (trkPr.sign() > 0)
              histos.fill(HIST("Analysis/h4d_lstar_invm_US_PM"), candMass, candPt, mult, occup);
            else
              histos.fill(HIST("Analysis/h4d_lstar_invm_US_MP"), candMass, candPt, mult, occup);
          }

          if (polActive) {
            float eProton = std::sqrt(pxPr * pxPr + pyPr * pyPr + pzPr * pzPr + MassProton * MassProton);
            float eKaon = std::sqrt(pxKa * pxKa + pyKa * pyKa + pzKa * pzKa + MassKaonCharged * MassKaonCharged);
            std::array<float, 3> protonP = {pxPr, pyPr, pzPr};
            std::array<float, 3> kaonP = {pxKa, pyKa, pzKa};
            std::array<float, 3> motherP = {pxPr + pxKa, pyPr + pyKa, pzPr + pzKa};
            float motherE = eProton + eKaon;
            fillPolarization<PolBkgMode::Signal>(candMass, candPt, motherP, motherE, protonP, eProton, kaonP, eKaon, mult, trkPr.sign() > 0);

            if (doRotate) {
              for (int i = 0; i < cNofRotations; i++) {
                float delta = o2::constants::math::PI / rotationalcut;
                float theta2 = (cNofRotations == 1) ? o2::constants::math::PI : (o2::constants::math::PI - delta) + i * (2.f * delta / (cNofRotations - 1));
                float phiRot = RecoDecay::constrainAngle(trkKa.phi() + theta2, 0.f);
                float pxKaRot = trkKa.pt() * std::cos(phiRot);
                float pyKaRot = trkKa.pt() * std::sin(phiRot);

                std::array<float, 3> pvec0rot = {pxPr, pyPr, pzPr};
                std::array<float, 3> pvec1rot = {pxKaRot, pyKaRot, pzKa};
                std::array<std::array<float, 3>, 2> arrMomRot = {pvec0rot, pvec1rot};

                float candMassRot = RecoDecay::m(arrMomRot, std::array{MassProton, MassKaonCharged});
                float candPtRot = RecoDecay::pt(std::array{pxPr + pxKaRot, pyPr + pyKaRot});
                float candYRot = std::abs(RecoDecay::y(std::array{pxPr + pxKaRot, pyPr + pyKaRot, pzPr + pzKa}, MassLambda1520));
                float candYShiftRot = candYRot - cfgRapidityShift;

                if (std::abs(candYShiftRot) > cfgRapidityCut)
                  continue;

                if (doPolarizationRot) {
                  float eKaonRot = std::sqrt(pxKaRot * pxKaRot + pyKaRot * pyKaRot + pzKa * pzKa + MassKaonCharged * MassKaonCharged);
                  std::array<float, 3> kaonRotP = {pxKaRot, pyKaRot, pzKa};
                  std::array<float, 3> motherRotP = {pxPr + pxKaRot, pyPr + pyKaRot, pzPr + pzKa};
                  float motherRotE = eProton + eKaonRot;
                  fillPolarization<PolBkgMode::Rotated>(candMassRot, candPtRot, motherRotP, motherRotE, protonP, eProton, kaonRotP, eKaonRot, mult, trkPr.sign() > 0);
                }
              }
            }
          } else if (doRotate) {
            for (int i = 0; i < cNofRotations; i++) {
              float delta = o2::constants::math::PI / rotationalcut;
              float theta2 = (cNofRotations == 1) ? o2::constants::math::PI : (o2::constants::math::PI - delta) + i * (2.f * delta / (cNofRotations - 1));
              float phiRot = RecoDecay::constrainAngle(trkKa.phi() + theta2, 0.f);
              float pxKaRot = trkKa.pt() * std::cos(phiRot);
              float pyKaRot = trkKa.pt() * std::sin(phiRot);
              std::array<float, 3> pvec0rot = {pxPr, pyPr, pzPr};
              std::array<float, 3> pvec1rot = {pxKaRot, pyKaRot, pzKa};
              std::array<std::array<float, 3>, 2> arrMomRot = {pvec0rot, pvec1rot};
              float candMassRot = RecoDecay::m(arrMomRot, std::array{MassProton, MassKaonCharged});
              float candPtRot = RecoDecay::pt(std::array{pxPr + pxKaRot, pyPr + pyKaRot});
              float candYRot = std::abs(RecoDecay::y(std::array{pxPr + pxKaRot, pyPr + pyKaRot, pzPr + pzKa}, MassLambda1520));

              if (std::abs(candYRot - cfgRapidityShift) > cfgRapidityCut)
                continue;

              if (trkPr.sign() > 0)
                histos.fill(HIST("Analysis/h4d_lstar_invm_rot_PM"), candMassRot, candPtRot, mult, occup);
              else
                histos.fill(HIST("Analysis/h4d_lstar_invm_rot_MP"), candMassRot, candPtRot, mult, occup);
            }
          }
        } else {
          if (!polActive) {
            if (trkPr.sign() > 0)
              histos.fill(HIST("Analysis/h4d_lstar_invm_PP"), candMass, candPt, mult, occup);
            else
              histos.fill(HIST("Analysis/h4d_lstar_invm_MM"), candMass, candPt, mult, occup);
          } else {
            float eProton = std::sqrt(pxPr * pxPr + pyPr * pyPr + pzPr * pzPr + MassProton * MassProton);
            float eKaon = std::sqrt(pxKa * pxKa + pyKa * pyKa + pzKa * pzKa + MassKaonCharged * MassKaonCharged);
            std::array<float, 3> protonP = {pxPr, pyPr, pzPr};
            std::array<float, 3> kaonP = {pxKa, pyKa, pzKa};
            std::array<float, 3> motherP = {pxPr + pxKa, pyPr + pyKa, pzPr + pzKa};
            float motherE = eProton + eKaon;
            fillPolarization<PolBkgMode::LikeSign>(candMass, candPt, motherP, motherE, protonP, eProton, kaonP, eKaon, mult, trkPr.sign() > 0);
          }
        }
      }

      if constexpr (mix) {
        if (trkPr.sign() * trkKa.sign() < 0) {
          if (!polActive) {
            if (trkPr.sign() > 0)
              histos.fill(HIST("Analysis/h4d_lstar_invm_US_PM_mix"), candMass, candPt, mult, occup);
            else
              histos.fill(HIST("Analysis/h4d_lstar_invm_US_MP_mix"), candMass, candPt, mult, occup);
          }
          if (doPolarizationMix && polActive) {
            float eProton = std::sqrt(pxPr * pxPr + pyPr * pyPr + pzPr * pzPr + MassProton * MassProton);
            float eKaon = std::sqrt(pxKa * pxKa + pyKa * pyKa + pzKa * pzKa + MassKaonCharged * MassKaonCharged);
            std::array<float, 3> protonP = {pxPr, pyPr, pzPr};
            std::array<float, 3> kaonP = {pxKa, pyKa, pzKa};
            std::array<float, 3> motherP = {pxPr + pxKa, pyPr + pyKa, pzPr + pzKa};
            float motherE = eProton + eKaon;
            fillPolarization<PolBkgMode::Mixed>(candMass, candPt, motherP, motherE, protonP, eProton, kaonP, eKaon, mult, trkPr.sign() > 0);
          }
        } else {
          if (!polActive) {
            if (trkPr.sign() > 0)
              histos.fill(HIST("Analysis/h4d_lstar_invm_LS_PP_mix"), candMass, candPt, mult, occup);
            else
              histos.fill(HIST("Analysis/h4d_lstar_invm_LS_MM_mix"), candMass, candPt, mult, occup);
          }
          if (doPolarizationMix && polActive) {
            float eProton = std::sqrt(pxPr * pxPr + pyPr * pyPr + pzPr * pzPr + MassProton * MassProton);
            float eKaon = std::sqrt(pxKa * pxKa + pyKa * pyKa + pzKa * pzKa + MassKaonCharged * MassKaonCharged);
            std::array<float, 3> protonP = {pxPr, pyPr, pzPr};
            std::array<float, 3> kaonP = {pxKa, pyKa, pzKa};
            std::array<float, 3> motherP = {pxPr + pxKa, pyPr + pyKa, pzPr + pzKa};
            float motherE = eProton + eKaon;
            fillPolarization<PolBkgMode::LikeSignMixed>(candMass, candPt, motherP, motherE, protonP, eProton, kaonP, eKaon, mult, trkPr.sign() > 0);
          }
        }
      }
      if constexpr (mc) {
        if (trkPr.sign() * trkKa.sign() < 0) {
          if (std::abs(trkPr.pdgCode()) != kProton || std::abs(trkKa.pdgCode()) != kKPlus)
            continue;
          if (trkPr.motherId() != trkKa.motherId())
            continue;
          if (trkPr.motherPDG() != trkKa.motherPDG())
            continue;
          if (trkPr.pdgCode() == 0 || trkKa.pdgCode() == 0)
            continue;
          if (trkPr.motherPDG() == -1 || trkKa.motherPDG() == -1)
            continue;
          if (std::abs(trkPr.motherPDG()) != lambda1520id)
            continue;

          float massParent = 0.;
          for (auto const& resoParent : *mResoParents) {
            if (resoParent.mcParticleId() == trkPr.motherId()) {
              std::array<float, 3> pvecParent = {resoParent.px(), resoParent.py(), resoParent.pz()};
              massParent = RecoDecay::m(pvecParent, resoParent.e());
              break;
            }
          }

          float candMassRes = candMass - massParent;
          if (trkPr.motherPDG() > 0) {
            histos.fill(HIST("Analysis/h3d_rec_lstar_PM"), candMass, candPt, mult);
            histos.fill(HIST("Analysis/h3d_reso_lstar_PM"), candMassRes, candPt, mult);
          } else {
            histos.fill(HIST("Analysis/h3d_rec_lstar_MP"), candMass, candPt, mult);
            histos.fill(HIST("Analysis/h3d_reso_lstar_MP"), candMassRes, candPt, mult);
          }
        }
      }
    }
  }

  using ResoCols = soa::Join<aod::ResoCollisions, aod::ResoEvtPlCollisions>;
  using ResoMCCols = soa::Join<aod::ResoCollisions, aod::ResoMCCollisions>;
  using ResoTracks = aod::ResoTracks;

  void processData(ResoCols::iterator const& collision, ResoTracks const& tracks)
  {
    if (cEvtMCRecINELgt0 && !collision.isRecINELgt0())
      return;
    histos.fill(HIST("Event/h1d_ft0_mult_percentile"), collision.cent(), 100);
    histos.fill(HIST("Event/h_ft0_vz"), collision.posZ());
    fillDataHistos<false, false>(tracks, tracks, collision.cent());
  }

  PROCESS_SWITCH(Lambda1520pbpb, processData, "Process for Same Event Data", true);

  void processMC(ResoMCCols::iterator const& collision, soa::Join<aod::ResoTracks, aod::ResoMCTracks> const& tracks, aod::ResoMCParents const& resoParents)
  {
    histos.fill(HIST("Event/hMCEventCutflow"), 0);
    if (cEvtMCTriggerTVX && !collision.isTriggerTVX())
      return;
    histos.fill(HIST("Event/hMCEventCutflow"), 1);
    if (cEvtMCVtxIn10 && !collision.isVtxIn10())
      return;
    histos.fill(HIST("Event/hMCEventCutflow"), 2);
    if (cEvtMCINELgt0 && !collision.isINELgt0())
      return;
    histos.fill(HIST("Event/hMCEventCutflow"), 3);
    if (cEvtMCSel8 && !collision.isInSel8())
      return;
    histos.fill(HIST("Event/hMCEventCutflow"), 4);
    if (cEvtMCRecINELgt0 && !collision.isRecINELgt0())
      return;
    histos.fill(HIST("Event/hMCEventCutflow"), 5);
    if (cEvtMCAfterAllCuts && !collision.isInAfterAllCuts())
      return;
    histos.fill(HIST("Event/hMCEventCutflow"), 6);

    auto mult = collision.cent();
    histos.fill(HIST("Event/h1d_ft0_mult_percentile"), mult, 100);
    histos.fill(HIST("Event/h_ft0_vz"), collision.posZ());
    mResoParents = &resoParents;
    fillDataHistos<false, true>(tracks, tracks, mult);

    for (auto const& track : tracks) {
      histos.fill(HIST("QAbefore/hEta_rec"), track.eta());
      histos.fill(HIST("QAbefore/hPt_rec"), track.pt());
      histos.fill(HIST("QAbefore/hPhi_rec"), track.phi());
      histos.fill(HIST("QAbefore/hEtaPhi_rec"), track.eta(), track.phi());

      if (std::abs(track.pdgCode()) == kKPlus)
        histos.fill(HIST("QAChecks/h1d_ka_gen_pt"), track.pt());
      if (std::abs(track.pdgCode()) == kProton)
        histos.fill(HIST("QAChecks/h1d_pr_gen_pt"), track.pt());

      if (!selTracks(track))
        continue;
      float p = std::sqrt(track.px() * track.px() + track.py() * track.py() + track.pz() * track.pz());

      if (selectionPIDKaon(track, p)) {
        if (std::abs(track.pdgCode()) == kKPlus)
          histos.fill(HIST("QAChecks/h1d_ka_rec_pt"), track.pt());
      }
      if (selectionPIDProton(track, p)) {
        if (std::abs(track.pdgCode()) == kProton)
          histos.fill(HIST("QAChecks/h1d_pr_rec_pt"), track.pt());
      }
    }

    for (auto const& part : resoParents) {
      if (std::abs(part.pdgCode()) != lambda1520id)
        continue;
      float yshift = std::abs(part.y()) - cfgRapidityShift;
      if (std::abs(yshift) > cfgRapidityCut)
        continue;

      bool pass1 = false, pass2 = false;
      if (std::abs(part.daughterPDG1()) == kProton || std::abs(part.daughterPDG2()) == kProton)
        pass1 = true;
      if (std::abs(part.daughterPDG1()) == kKPlus || std::abs(part.daughterPDG2()) == kKPlus)
        pass2 = true;

      if (!pass1 || !pass2)
        continue;
      std::array<float, 3> pvec = {part.px(), part.py(), part.pz()};
      float mass = RecoDecay::m(pvec, part.e());
      if (part.pdgCode() > 0)
        histos.fill(HIST("Analysis/h3d_gen_lstar_PM"), mass, part.pt(), mult);
      else
        histos.fill(HIST("Analysis/h3d_gen_lstar_MP"), mass, part.pt(), mult);
    }
  }
  PROCESS_SWITCH(Lambda1520pbpb, processMC, "Process Event for MC", false);

  void processMCGen(ResoMCCols::iterator const& collision, aod::ResoMCParents const& resoParents)
  {
    float centrality = collision.cent();
    histos.fill(HIST("SignalLoss/hMCEventCutflow"), 0);
    if (cEvtMCTriggerTVX && !collision.isTriggerTVX())
      return;
    histos.fill(HIST("SignalLoss/hMCEventCutflow"), 1);
    if (cEvtMCVtxIn10 && !collision.isVtxIn10())
      return;
    histos.fill(HIST("SignalLoss/hMCEventCutflow"), 2);
    if (cEvtMCINELgt0 && !collision.isINELgt0())
      return;
    histos.fill(HIST("SignalLoss/hMCEventCutflow"), 3);
    if (cEvtMCSel8 && !collision.isInSel8())
      return;
    histos.fill(HIST("SignalLoss/hMCEventCutflow"), 4);
    if (cEvtMCRecINELgt0 && !collision.isRecINELgt0())
      return;
    histos.fill(HIST("SignalLoss/hMCEventCutflow"), 5);
    if (cEvtMCAfterAllCuts && !collision.isInAfterAllCuts())
      return;
    histos.fill(HIST("SignalLoss/hMCEventCutflow"), 6);

    for (auto const& part : resoParents) {
      float yshift = std::abs(part.y()) - cfgRapidityShift;
      if (std::abs(yshift) > cfgRapidityCut)
        continue;

      int pdg = part.pdgCode();
      float ptRef = part.pt();
      double ptSq = -1.0;

      std::array<float, 3> pvec = {part.px(), part.py(), part.pz()};
      float mass = RecoDecay::m(pvec, part.e());

      if (pdg == kProton) {
        ptSq = (ptRef * ptRef) + (mass * mass) - (o2::constants::physics::MassLambda1520 * o2::constants::physics::MassLambda1520);
        if (ptSq > 0)
          histos.fill(HIST("SignalLoss/hGen_mT_scaled_Proton"), std::sqrt(ptSq), centrality);
      } else if (pdg == -kProton) {
        ptSq = (ptRef * ptRef) + (mass * mass) - (o2::constants::physics::MassLambda1520 * o2::constants::physics::MassLambda1520);
        if (ptSq > 0)
          histos.fill(HIST("SignalLoss/hGen_mT_scaled_AntiProton"), std::sqrt(ptSq), centrality);
      } else if (pdg == kLambda0) {
        ptSq = (ptRef * ptRef) + (mass * mass) - (o2::constants::physics::MassLambda1520 * o2::constants::physics::MassLambda1520);
        if (ptSq > 0)
          histos.fill(HIST("SignalLoss/hGen_mT_scaled_Lambda0"), std::sqrt(ptSq), centrality);
      } else if (pdg == -kLambda0) {
        ptSq = (ptRef * ptRef) + (mass * mass) - (o2::constants::physics::MassLambda1520 * o2::constants::physics::MassLambda1520);
        if (ptSq > 0)
          histos.fill(HIST("SignalLoss/hGen_mT_scaled_AntiLambda0"), std::sqrt(ptSq), centrality);
      } else if (pdg == kXiMinus) {
        ptSq = (ptRef * ptRef) + (mass * mass) - (o2::constants::physics::MassLambda1520 * o2::constants::physics::MassLambda1520);
        if (ptSq > 0)
          histos.fill(HIST("SignalLoss/hGen_mT_scaled_XiMinus"), std::sqrt(ptSq), centrality);
      } else if (pdg == -kXiMinus) {
        ptSq = (ptRef * ptRef) + (mass * mass) - (o2::constants::physics::MassLambda1520 * o2::constants::physics::MassLambda1520);
        if (ptSq > 0)
          histos.fill(HIST("SignalLoss/hGen_mT_scaled_XiPlus"), std::sqrt(ptSq), centrality);
      } else if (pdg == kXi0) {
        ptSq = (ptRef * ptRef) + (mass * mass) - (o2::constants::physics::MassLambda1520 * o2::constants::physics::MassLambda1520);
        if (ptSq > 0)
          histos.fill(HIST("SignalLoss/hGen_mT_scaled_Xi0"), std::sqrt(ptSq), centrality);
      } else if (pdg == -kXi0) {
        ptSq = (ptRef * ptRef) + (mass * mass) - (o2::constants::physics::MassLambda1520 * o2::constants::physics::MassLambda1520);
        if (ptSq > 0)
          histos.fill(HIST("SignalLoss/hGen_mT_scaled_AntiXi0"), std::sqrt(ptSq), centrality);
      } else if (pdg == kOmegaMinus) {
        ptSq = (ptRef * ptRef) + (mass * mass) - (o2::constants::physics::MassLambda1520 * o2::constants::physics::MassLambda1520);
        if (ptSq > 0)
          histos.fill(HIST("SignalLoss/hGen_mT_scaled_OmegaMinus"), std::sqrt(ptSq), centrality);
      } else if (pdg == -kOmegaMinus) {
        ptSq = (ptRef * ptRef) + (mass * mass) - (o2::constants::physics::MassLambda1520 * o2::constants::physics::MassLambda1520);
        if (ptSq > 0)
          histos.fill(HIST("SignalLoss/hGen_mT_scaled_OmegaPlus"), std::sqrt(ptSq), centrality);
      }
    }
  }

  PROCESS_SWITCH(Lambda1520pbpb, processMCGen, "Process Event for MC", false);

  using BinningType2 = ColumnBinningPolicy<aod::collision::PosZ, aod::resocollision::Cent>;

  void processMix(ResoCols const& collisions, ResoTracks const& tracks)
  {
    LOGF(debug, "Event Mixing Started");
    BinningType2 binningPositions2{{cMixVtxBins, cMixMultBins}, true};
    auto tracksTuple = std::make_tuple(tracks);

    SameKindPair<ResoCols, ResoTracks, BinningType2> pairs{binningPositions2, cNumMixEv, -1, collisions, tracksTuple, &cache};
    for (auto const& [c1, t1, c2, t2] : pairs) {
      if (cEvtMCRecINELgt0 && !c1.isRecINELgt0())
        return;
      histos.fill(HIST("Event/mixing_vzVsmultpercentile"), c1.cent(), c1.posZ(), c1.evtPl());
      fillDataHistos<true, false>(t1, t2, c1.cent());
    }
  }

  PROCESS_SWITCH(Lambda1520pbpb, processMix, "Process for Mixed Events", true);

  Preslice<aod::ResoTrackDFs> perRColdf = aod::resodaughter::resoCollisionDFId;
  using ResoColDFs = aod::ResoCollisionDFs;
  using ResoTrackDFs = aod::ResoTrackDFs;

  void processDatadf(ResoColDFs::iterator const& collision, ResoTrackDFs const& tracks)
  {
    if (doprocessData)
      LOG(error) << "Disable processData() first!";
    if (cEvtMCRecINELgt0 && !collision.isRecINELgt0())
      return;

    auto occup = 100;
    if (cfgEvtOccupancyInTimeRange)
      occup = collision.trackOccupancyInTimeRange();

    histos.fill(HIST("Event/h1d_ft0_mult_percentile"), collision.cent(), occup);
    fillDataHistos<false, false>(tracks, tracks, collision.cent(), occup);
  }

  PROCESS_SWITCH(Lambda1520pbpb, processDatadf, "Process for data merged DF", false);

  using BinningTypeDF = ColumnBinningPolicy<aod::collision::PosZ, aod::resocollision::Cent>;
  void processMixDF(ResoColDFs const& collisions, ResoTrackDFs const& tracks)
  {
    if (doprocessMix)
      LOG(fatal) << "Disable processMix() first!";
    LOGF(debug, "Event Mixing Started");

    BinningTypeDF binningPositions2{{cMixVtxBins, cMixMultBins}, true};
    auto tracksTuple = std::make_tuple(tracks);

    SameKindPair<ResoColDFs, ResoTrackDFs, BinningTypeDF> pairs{binningPositions2, cNumMixEv, -1, collisions, tracksTuple, &cache};
    for (auto const& [c1, t1, c2, t2] : pairs) {
      if (cEvtMCRecINELgt0 && !c1.isRecINELgt0())
        return;

      auto occup = 100;
      if (cfgEvtOccupancyInTimeRange)
        occup = c1.trackOccupancyInTimeRange();

      histos.fill(HIST("Event/mixing_vzVsmultpercentile"), c1.cent(), c1.posZ(), c1.evtPl());
      fillDataHistos<true, false>(t1, t2, c1.cent(), occup);
    }
  }

  PROCESS_SWITCH(Lambda1520pbpb, processMixDF, "Process for merged DF  Mixed Events", false);

  using BinningTypeEP = ColumnBinningPolicy<aod::collision::PosZ, aod::resocollision::Cent, aod::resocollision::EvtPl>;
  void processMixepDF(ResoColDFs const& collisions, ResoTrackDFs const& tracks)
  {
    if (doprocessMix || doprocessMixDF)
      LOG(fatal) << "Disable processMix() or processMixDF() first!";
    LOGF(debug, "Event Mixing Started");
    BinningTypeEP binningPositions2{{cMixVtxBins, cMixMultBins, cMixEPAngle}, true};
    auto tracksTuple = std::make_tuple(tracks);

    SameKindPair<ResoColDFs, ResoTrackDFs, BinningTypeEP> pairs{binningPositions2, cNumMixEv, -1, collisions, tracksTuple, &cache};
    for (auto const& [c1, t1, c2, t2] : pairs) {
      histos.fill(HIST("Event/mixing_vzVsmultpercentile"), c1.cent(), c1.posZ(), c1.evtPl());
      fillDataHistos<true, false>(t1, t2, c1.cent());
    }
  }

  PROCESS_SWITCH(Lambda1520pbpb, processMixepDF, "Process for merged DF  Mixed Events", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<Lambda1520pbpb>(cfgc)};
}
