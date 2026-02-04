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

/// \file k892_spherocity_analysis.cxx
/// \brief Invariant Mass Reconstruction of K*(892) Resonance
/// \author Sayan Dhani <sayan.dhani@cern.ch>

#include "PWGLF/DataModel/LFResonanceTables.h"

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include <TLorentzVector.h>
#include <TRandom.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

struct k892Analysis {

  SliceCache cache;
  Preslice<aod::ResoTracks> perRCol = aod::resodaughter::resoCollisionId;
  Preslice<aod::Tracks> perCollision = aod::track::collisionId;

  // Configurables.
  Configurable<int> nBinsPt{"nBinsPt", 100, "N bins in pT histogram"};
  Configurable<int> nBinsInvM{"nBinsInvM", 120, "N bins in InvMass histogram"};
  Configurable<int> nBinsSp{"nBinsSp", 120, "N bins in spherocity histogram"};
  Configurable<bool> doRotate{"doRotate", true, "rotated inv mass spectra"};

  // Tracks
  Configurable<float> cPtMin{"cPtMin", 0.15, "Minimum Track pT"};
  Configurable<float> cEtaCut{"cEtaCut", 0.8, "Pseudorapidity cut"};
  Configurable<float> cDcaz{"cDcazMin", 1., "Minimum DCAz"};
  Configurable<float> cDcaxy{"cDcaxyMin", 0.1, "Minimum DCAxy"};
  Configurable<float> cPIDprecut{"cPIDprecut", 5, "Preselection PID TPC TOF cut"};
  Configurable<bool> cKinCuts{"cKinCuts", false, "Kinematic Cuts for p-K pair opening angle"};
  Configurable<bool> cPrimaryTrack{"cPrimaryTrack", true, "Primary track selection"};                    // kGoldenChi2 | kDCAxy | kDCAz
  Configurable<bool> cGlobalWoDCATrack{"cGlobalWoDCATrack", true, "Global track selection without DCA"}; // kQualityTracks (kTrackType | kTPCNCls | kTPCCrossedRows | kTPCCrossedRowsOverNCls | kTPCChi2NDF | kTPCRefit | kITSNCls | kITSChi2NDF | kITSRefit | kITSHits) | kInAcceptanceTracks (kPtRange | kEtaRange)
  Configurable<bool> cPVContributor{"cPVContributor", true, "PV contributor track selection"};           // PV Contriuibutor

  // PID Selections
  Configurable<bool> cUseOnlyTOFTrackPi{"cUseOnlyTOFTrackPi", false, "Use only TOF track for PID selection"}; // Use only TOF track for Pion PID selection
  Configurable<bool> cUseOnlyTOFTrackKa{"cUseOnlyTOFTrackKa", false, "Use only TOF track for PID selection"}; // Use only TOF track for Kaon PID selection
  Configurable<bool> cUseTpcOnly{"cUseTpcOnly", false, "Use TPC Only selection"};                             // TPC And TOF tracks
  Configurable<float> cRejNsigma{"cRejNsigma", 1.0, "Reject tracks to improve purity of PID"};                // Reject missidentified particles when tpc bands merge
  // Pion
  Configurable<double> cMaxTPCnSigmaPion{"cMaxTPCnSigmaPion", 3.0, "TPC nSigma cut for Pion"};              // TPC
  Configurable<double> cMaxTOFnSigmaPion{"cMaxTOFnSigmaPion", 3.0, "TOF nSigma cut for Pion"};              // TOF
  Configurable<double> nsigmaCutCombinedPion{"nsigmaCutCombinedPion", 3.0, "Combined nSigma cut for Pion"}; // Combined
  Configurable<std::vector<float>> pionTPCPIDp{"pionTPCPIDp", {0, 0.5, 0.7, 0.8}, "pT dependent TPC cuts pions"};
  Configurable<std::vector<float>> pionTPCPIDcut{"pionsTPCPIDcut", {5., 3.5, 2.5}, "TPC nsigma cuts pions"};
  // Kaon
  Configurable<double> cMaxTPCnSigmaKaon{"cMaxTPCnSigmaKaon", 3.0, "TPC nSigma cut for Kaon"};              // TPC
  Configurable<double> cMaxTOFnSigmaKaon{"cMaxTOFnSigmaKaon", 3.0, "TOF nSigma cut for Kaon"};              // TOF
  Configurable<double> nsigmaCutCombinedKaon{"nsigmaCutCombinedKaon", 3.0, "Combined nSigma cut for Kaon"}; // Combined
  Configurable<std::vector<float>> kaonTPCPIDp{"kaonTPCPIDp", {0., 0.25, 0.3, 0.45}, "pT dependent TPC cuts kaons"};
  Configurable<std::vector<float>> kaonTPCPIDcut{"kaonTPCPIDcut", {6, 3.5, 2.5}, "TPC nsigma cuts kaons"};
  // Event Mixing.
  Configurable<bool> cMixSph{"cMixSph", true, "Include Sph Bins to be mixed"};
  Configurable<int> cNumMixEv{"cNumMixEv", 20, "Number of Events to be mixed"};
  ConfigurableAxis cMixVtxBins{"cMixVtxBins", {VARIABLE_WIDTH, -10.0f, -9.f, -8.f, -7.f, -6.f, -5.f, -4.f, -3.f, -2.f, -1.f, 0.f, 1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis cMixMultBins{"cMixMultBins", {VARIABLE_WIDTH, 0.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.0f, 200.0f}, "Mixing bins - multiplicity"};
  ConfigurableAxis cMixSphBins{"cMixSphBins", {VARIABLE_WIDTH, 0.0f, 0.2f, 0.4f, 0.6f, 0.8f, 1.0f}, "Mixing bins - spherocity"};

  // Histogram Registry.
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {

    // Define Axis.
    const AxisSpec axisSp(nBinsSp, 0., 1., "S_{0}");
    const AxisSpec axisCent(105, 0, 105, "FT0M (%)");
    const AxisSpec axisP_pid(400, 0., 4., "p (GeV/c)");
    const AxisSpec axisPt_pid(400, 0., 4., "p_{T} (GeV/c)");
    const AxisSpec axisPt(nBinsPt, 0., 10., "p_{T} (GeV/c)");
    const AxisSpec axisEta(40, -1, 1, "#eta");
    const AxisSpec axisDCAz(500, -0.5, 0.5, {"DCA_{z} (cm)"});
    const AxisSpec axisDCAxy(240, -0.12, 0.12, {"DCA_{xy} (cm)"});
    const AxisSpec axisTPCNCls(200, 0, 200, {"TPCNCls"});
    const AxisSpec axisTPCNsigma(120, -6, 6, {"n#sigma^{TPC}"});
    const AxisSpec axisTOFNsigma(120, -6, 6, {"n#sigma^{TOF}"});
    const AxisSpec axisdEdx(380, 10, 200, {"#frac{dE}{dx}"});
    const AxisSpec axisInvM(nBinsInvM, 0.64, 1.04, {"M_{inv} (GeV/c^{2})"});

    // Create Histograms.
    // Event
    histos.add("Event/h1d_ft0m_mult_percentile", "FT0M (%)", kTH1F, {axisCent});
    histos.add("Event/h1d_spherocity", "Event Spherocity", kTH1F, {axisSp});
    histos.add("Event/h2d_sph_vs_multpercentile", "Spherocity vs FT0M(%)", kTH2F, {axisCent, axisSp});

    // QA Before
    histos.add("QAbefore/Pion/h2d_pi_nsigma_tpc_p", "n#sigma^{TPC} Pions", kTH2F, {axisP_pid, axisTPCNsigma});
    histos.add("QAbefore/Pion/h2d_pi_nsigma_tof_p", "n#sigma^{TOF} Pions", kTH2F, {axisP_pid, axisTOFNsigma});
    histos.add("QAbefore/Pion/h2d_pi_nsigma_tof_vs_tpc", "n#sigma^{TPC} vs n#sigma^{TOF} Pions", kTH2F, {axisTPCNsigma, axisTOFNsigma});
    histos.add("QAbefore/Kaon/h2d_ka_nsigma_tpc_p", "n#sigma^{TPC} Kaons", kTH2F, {axisP_pid, axisTPCNsigma});
    histos.add("QAbefore/Kaon/h2d_ka_nsigma_tof_p", "n#sigma^{TOF} Kaons", kTH2F, {axisP_pid, axisTOFNsigma});
    histos.add("QAbefore/Kaon/h2d_ka_nsigma_tof_vs_tpc", "n#sigma^{TPC} vs n#sigma^{TOF} Kaons", kTH2F, {axisTPCNsigma, axisTOFNsigma});

    // QA After
    histos.add("QAafter/Pion/h1d_pi_pt", "p_{T}-spectra Pions", kTH1F, {axisPt});
    histos.add("QAafter/Pion/h2d_pi_dca_z", "dca_{z} Pions", kTH2F, {axisPt_pid, axisDCAz});
    histos.add("QAafter/Pion/h2d_pi_dca_xy", "dca_{xy} Pions", kTH2F, {axisPt_pid, axisDCAxy});
    histos.add("QAafter/Pion/h2d_pi_dEdx_p", "TPC Signal Pions", kTH2F, {axisP_pid, axisdEdx});
    histos.add("QAafter/Pion/h2d_pi_nsigma_tpc_pt", " Pions", kTH2F, {axisPt_pid, axisTPCNsigma});
    histos.add("QAafter/Pion/h2d_pi_nsigma_tpc_p", " Pions", kTH2F, {axisP_pid, axisTPCNsigma});
    histos.add("QAafter/Pion/h2d_pi_nsigma_tof_pt", " Pions", kTH2F, {axisPt_pid, axisTOFNsigma});
    histos.add("QAafter/Pion/h2d_pi_nsigma_tof_p", " Pions", kTH2F, {axisP_pid, axisTOFNsigma});
    histos.add("QAafter/Pion/h2d_pi_nsigma_tof_vs_tpc", "n#sigma(TOF) vs n#sigma(TPC) Pions", kTH2F, {axisTPCNsigma, axisTOFNsigma});
    histos.add("QAafter/Kaon/h1d_ka_pt", "p_{T}-spectra Kaons", kTH1F, {axisPt});
    histos.add("QAafter/Kaon/h2d_ka_dca_z", "dca_{z} Kaons", kTH2F, {axisPt_pid, axisDCAz});
    histos.add("QAafter/Kaon/h2d_ka_dca_xy", "dca_{xy} Kaons", kTH2F, {axisPt_pid, axisDCAxy});
    histos.add("QAafter/Kaon/h2d_ka_dEdx_p", "TPC Signal Kaon", kTH2F, {axisP_pid, axisdEdx});
    histos.add("QAafter/Kaon/h2d_ka_nsigma_tpc_pt", " Kaons", kTH2F, {axisPt_pid, axisTPCNsigma});
    histos.add("QAafter/Kaon/h2d_ka_nsigma_tpc_p", " Kaons", kTH2F, {axisP_pid, axisTPCNsigma});
    histos.add("QAafter/Kaon/h2d_ka_nsigma_tof_pt", " Kaons", kTH2F, {axisPt_pid, axisTOFNsigma});
    histos.add("QAafter/Kaon/h2d_ka_nsigma_tof_p", " Kaons", kTH2F, {axisP_pid, axisTOFNsigma});
    histos.add("QAafter/Kaon/h2d_ka_nsigma_tof_vs_tpc", "n#sigma(TOF) vs n#sigma(TPC) Kaons", kTH2F, {axisTPCNsigma, axisTOFNsigma});

    // Analysis
    // k892 Invariant Mass
    histos.add("Analysis/h1d_kstar_invm_US", "K(892) M_{inv}", kTH1D, {axisInvM});
    histos.add("Analysis/h1d_kstar_invm_PP", "Like Signs M_{inv} #{pi}^{+} K^{+}", kTH1D, {axisInvM});
    histos.add("Analysis/h1d_kstar_invm_MM", "Like Signs M_{inv} #{pi}^{-} K^{-}", kTH1D, {axisInvM});
    histos.add("Analysis/h1d_kstar_invm_rot", "Rotated Spectra", kTH1D, {axisInvM});
    histos.add("Analysis/h1d_kstar_invm_US_mix", "Mixed Events M_{inv}", kTH1D, {axisInvM});
    histos.add("Analysis/h1d_kstar_invm_LS_mix", "Mixed Events M_{inv}", kTH1D, {axisInvM});
    histos.add("Analysis/h4d_kstar_invm_US", "THn #K(892)", kTHnSparseD, {axisInvM, axisPt, axisSp, axisCent});
    histos.add("Analysis/h4d_kstar_invm_PP", "THn Like Signs #{pi}^{+} K^{+}", kTHnSparseD, {axisInvM, axisPt, axisSp, axisCent});
    histos.add("Analysis/h4d_kstar_invm_MM", "THn Like Signs #{pi}^{-} K^{-}", kTHnSparseD, {axisInvM, axisPt, axisSp, axisCent});
    histos.add("Analysis/h4d_kstar_invm_rot", "THn Rotated", kTHnSparseD, {axisInvM, axisPt, axisSp, axisCent});
    histos.add("Analysis/h4d_kstar_invm_US_mix", "THn Mixed Events", kTHnSparseD, {axisInvM, axisPt, axisSp, axisCent});
    histos.add("Analysis/h4d_kstar_invm_LS_mix", "THn Mixed Events", kTHnSparseD, {axisInvM, axisPt, axisSp, axisCent});

    // MC}
    if (doprocessMC) {
      histos.add("Event/h1d_rec_sph", "Reconstructed S_{0}", kTH1F, {axisSp});
      histos.add("Event/h1d_rec_sph_vs_cent", "Reconstructed S_{0} vs FT0M(%)", kTH2F, {axisCent, axisSp});
      histos.add("Analysis/h1d_gen_kstar", "Generated K(892) p_{T}", kTH1D, {axisPt});
      histos.add("Analysis/h1d_gen_kstar_anti", "Generated #bar{K}(892) p_{T}", kTH1D, {axisPt});
      histos.add("Analysis/h1d_rec_kstar", "Reconstructed K(892) p_{T}", kTH1D, {axisPt});
      histos.add("Analysis/h1d_rec_kstar_anti", "Reconstructed #bar{K}(892) p_{T}", kTH1D, {axisPt});
      histos.add("Analysis/h1d_rec_invm_kstar", "Recostructed K(892)", kTH1D, {axisInvM});
      histos.add("Analysis/h1d_rec_invm_kstar_anti", "Recostructed #bar{K}(892)", kTH1D, {axisInvM});
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
  bool selectionPIDPion(const T& candidate, float p)
  {
    bool tpcPIDPassed{false}, tofPIDPassed{false};
    auto tpcPIDp = static_cast<std::vector<float>>(pionTPCPIDp);
    auto tpcPIDcut = static_cast<std::vector<float>>(pionTPCPIDcut);
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
    float combinedCut = nsigmaCutCombinedPion * nsigmaCutCombinedPion;
    float combinedRejCut = cRejNsigma * cRejNsigma;

    if (!cUseTpcOnly && candidate.hasTOF()) {
      if (tofNsigmaPi < cMaxTOFnSigmaPion && tofNsigmaPr > cRejNsigma && tofNsigmaKa > cRejNsigma) {
        tofPIDPassed = true;
      }
      // square cut
      if ((nsigmaCutCombinedPion < 0) && (tpcNsigmaPi < cMaxTPCnSigmaPion)) {
        tpcPIDPassed = true;
      }
      // circular cut
      if ((nsigmaCutCombinedPion > 0) && (tpcTofNsigmaPi < combinedCut && tpcTofNsigmaPr > combinedRejCut && tpcTofNsigmaKa > combinedRejCut)) {
        tofPIDPassed = true;
        tpcPIDPassed = true;
      }
    } else {
      tofPIDPassed = true;
      if (cUseTpcOnly) {
        if (tpcNsigmaPi < cMaxTPCnSigmaPion && tpcNsigmaPr > cRejNsigma && tpcNsigmaKa > cRejNsigma) {
          tpcPIDPassed = true;
        }
      } else {
        for (int i = 0; i < nitr - 1; ++i) {
          if (p >= tpcPIDp[i] && p < tpcPIDp[i + 1] && (tpcNsigmaPi < tpcPIDcut[i] && tpcNsigmaPr > cRejNsigma && tpcNsigmaKa > cRejNsigma)) {
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

    float tpcTofNsigmaPi = tpcNsigmaPi * tpcNsigmaPi + tofNsigmaPi * tofNsigmaPi;
    float tpcTofNsigmaKa = tpcNsigmaKa * tpcNsigmaKa + tofNsigmaKa * tofNsigmaKa;
    float tpcTofNsigmaPr = tpcNsigmaPr * tpcNsigmaPr + tofNsigmaPr * tofNsigmaPr;
    float combinedCut = nsigmaCutCombinedKaon * nsigmaCutCombinedKaon;
    float combinedRejCut = cRejNsigma * cRejNsigma;

    if (!cUseTpcOnly && candidate.hasTOF()) {
      if (tofNsigmaKa < cMaxTOFnSigmaKaon && tofNsigmaPi > cRejNsigma && tofNsigmaPr > cRejNsigma) {
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
        if (tpcNsigmaKa < cMaxTPCnSigmaKaon && tpcNsigmaPi > cRejNsigma && tpcNsigmaPr > cRejNsigma) {
          tpcPIDPassed = true;
        }
      } else {
        for (int i = 0; i < nitr - 1; ++i) {
          if (p >= tpcPIDp[i] && p < tpcPIDp[i + 1] && (tpcNsigmaKa < tpcPIDcut[i] && tpcNsigmaPi > cRejNsigma && tpcNsigmaPr > cRejNsigma)) {
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
  void fillDataHistos(trackType const& trk1, trackType const& trk2, float const& sph, float const& mult)
  {
    TLorentzVector p1, p2, p;
    TRandom* rn = new TRandom();
    float p_ptot = 0., k_ptot = 0.;

    for (auto const& [trkPi, trkKa] : soa::combinations(soa::CombinationsFullIndexPolicy(trk1, trk2))) {
      // Do not analyse same index tracks.
      if (trkPi.index() == trkKa.index())
        continue;

      // pT, DCA, Global Tracks and PVcontrib selection.
      if (!selTracks(trkPi) || !selTracks(trkKa))
        continue;

      p_ptot = TMath::Sqrt(trkPi.px() * trkPi.px() + trkPi.py() * trkPi.py() + trkPi.pz() * trkPi.pz());
      k_ptot = TMath::Sqrt(trkKa.px() * trkKa.px() + trkKa.py() * trkKa.py() + trkKa.pz() * trkKa.pz());

      // Fill QA before track selection.
      if (!mix) {
        histos.fill(HIST("QAbefore/Pion/h2d_pi_nsigma_tpc_p"), p_ptot, trkPi.tpcNSigmaPi());
        if (trkPi.hasTOF()) {
          histos.fill(HIST("QAbefore/Pion/h2d_pi_nsigma_tof_p"), p_ptot, trkPi.tofNSigmaPi());
          histos.fill(HIST("QAbefore/Pion/h2d_pi_nsigma_tof_vs_tpc"), trkPi.tpcNSigmaPi(), trkPi.tofNSigmaPi());
        }
        histos.fill(HIST("QAbefore/Kaon/h2d_ka_nsigma_tpc_p"), k_ptot, trkKa.tpcNSigmaKa());
        if (trkKa.hasTOF()) {
          histos.fill(HIST("QAbefore/Kaon/h2d_ka_nsigma_tof_p"), k_ptot, trkKa.tofNSigmaKa());
          histos.fill(HIST("QAbefore/Kaon/h2d_ka_nsigma_tof_vs_tpc"), trkKa.tpcNSigmaKa(), trkKa.tofNSigmaKa());
        }
      }

      // Apply PID Selection
      if (cUseOnlyTOFTrackPi && !trkPi.hasTOF())
        continue;
      if (cUseOnlyTOFTrackKa && !trkKa.hasTOF())
        continue;
      if (!selectionPIDPion(trkPi, p_ptot) || !selectionPIDKaon(trkKa, k_ptot))
        continue;

      // Fill QA after track selection.
      if constexpr (!mix) {
        // Pion
        histos.fill(HIST("QAafter/Pion/h1d_pi_pt"), trkPi.pt());
        histos.fill(HIST("QAafter/Pion/h2d_pi_dca_z"), trkPi.pt(), trkPi.dcaZ());
        histos.fill(HIST("QAafter/Pion/h2d_pi_dca_xy"), trkPi.pt(), trkPi.dcaXY());
        histos.fill(HIST("QAafter/Pion/h2d_pi_dEdx_p"), p_ptot, trkPi.tpcSignal());
        histos.fill(HIST("QAafter/Pion/h2d_pi_nsigma_tpc_p"), p_ptot, trkPi.tpcNSigmaPi());
        histos.fill(HIST("QAafter/Pion/h2d_pi_nsigma_tpc_pt"), trkPi.pt(), trkPi.tpcNSigmaPi());
        if (!cUseTpcOnly && trkPi.hasTOF()) {
          histos.fill(HIST("QAafter/Pion/h2d_pi_nsigma_tof_p"), p_ptot, trkPi.tofNSigmaPi());
          histos.fill(HIST("QAafter/Pion/h2d_pi_nsigma_tof_pt"), trkPi.pt(), trkPi.tofNSigmaPi());
          histos.fill(HIST("QAafter/Pion/h2d_pi_nsigma_tof_vs_tpc"), trkPi.tpcNSigmaPi(), trkPi.tofNSigmaPi());
        }
        // Kaon
        histos.fill(HIST("QAafter/Kaon/h1d_ka_pt"), trkKa.pt());
        histos.fill(HIST("QAafter/Kaon/h2d_ka_dca_z"), trkKa.pt(), trkKa.dcaZ());
        histos.fill(HIST("QAafter/Kaon/h2d_ka_dca_xy"), trkKa.pt(), trkKa.dcaXY());
        histos.fill(HIST("QAafter/Kaon/h2d_ka_dEdx_p"), k_ptot, trkKa.tpcSignal());
        histos.fill(HIST("QAafter/Kaon/h2d_ka_nsigma_tpc_p"), k_ptot, trkKa.tpcNSigmaKa());
        histos.fill(HIST("QAafter/Kaon/h2d_ka_nsigma_tpc_pt"), trkKa.pt(), trkKa.tpcNSigmaKa());
        if (!cUseTpcOnly && trkKa.hasTOF()) {
          histos.fill(HIST("QAafter/Kaon/h2d_ka_nsigma_tof_p"), k_ptot, trkKa.tofNSigmaKa());
          histos.fill(HIST("QAafter/Kaon/h2d_ka_nsigma_tof_pt"), trkKa.pt(), trkKa.tofNSigmaKa());
          histos.fill(HIST("QAafter/Kaon/h2d_ka_nsigma_tof_vs_tpc"), trkKa.tpcNSigmaKa(), trkKa.tofNSigmaKa());
        }
      }

      // Invariant mass reconstruction.
      p1.SetXYZM(trkPi.px(), trkPi.py(), trkPi.pz(), MassPionCharged);
      p2.SetXYZM(trkKa.px(), trkKa.py(), trkKa.pz(), MassKaonCharged);
      p = p1 + p2;

      if (std::abs(p.Rapidity()) > 0.5)
        continue;

      // Apply kinematic cuts.
      if (cKinCuts) {
        TVector3 v1(trkPi.px(), trkPi.py(), trkPi.pz());
        TVector3 v2(trkKa.px(), trkKa.py(), trkKa.pz());
        float alpha = v1.Angle(v2);
        if (alpha > 1.4 && alpha < 2.4)
          continue;
      }

      // Fill Invariant Mass Histograms.
      if constexpr (!mix && !mc) {
        if (trkPi.sign() * trkKa.sign() < 0) {
          histos.fill(HIST("Analysis/h1d_kstar_invm_US"), p.M());
          histos.fill(HIST("Analysis/h4d_kstar_invm_US"), p.M(), p.Pt(), sph, mult);
          if (doRotate) {
            float theta = rn->Uniform(1.56, 1.58);
            p1.RotateZ(theta);
            p = p1 + p2;
            if (std::abs(p.Rapidity()) < 0.5) {
              histos.fill(HIST("Analysis/h1d_kstar_invm_rot"), p.M());
              histos.fill(HIST("Analysis/h4d_kstar_invm_rot"), p.M(), p.Pt(), sph, mult);
            }
          }
        } else {
          if (trkPi.sign() == 1) {
            histos.fill(HIST("Analysis/h1d_kstar_invm_PP"), p.M());
            histos.fill(HIST("Analysis/h4d_kstar_invm_PP"), p.M(), p.Pt(), sph, mult);
          } else {
            histos.fill(HIST("Analysis/h1d_kstar_invm_MM"), p.M());
            histos.fill(HIST("Analysis/h4d_kstar_invm_MM"), p.M(), p.Pt(), sph, mult);
          }
        }
      }

      if constexpr (mix) {
        if (trkPi.sign() * trkKa.sign() < 0) {
          histos.fill(HIST("Analysis/h1d_kstar_invm_US_mix"), p.M());
          histos.fill(HIST("Analysis/h4d_kstar_invm_US_mix"), p.M(), p.Pt(), sph, mult);
        } else {
          histos.fill(HIST("Analysis/h1d_kstar_invm_LS_mix"), p.M());
          histos.fill(HIST("Analysis/h4d_kstar_invm_LS_mix"), p.M(), p.Pt(), sph, mult);
        }
      }

      if constexpr (mc) {
        if (abs(trkPi.pdgCode()) != 211 || abs(trkKa.pdgCode()) != 321)
          continue;

        if (trkPi.motherId() != trkKa.motherId())
          continue;

        if (abs(trkPi.motherPDG()) != 313) // L* pdg_code = 313
          continue;

        // MC histograms
        if (trkPi.motherPDG() > 0) {
          histos.fill(HIST("Analysis/h1d_rec_kstar"), p.Pt());
          histos.fill(HIST("Analysis/h1d_rec_invm_kstar"), p.M());
        } else {
          histos.fill(HIST("Analysis/h1d_rec_kstar_anti"), p.Pt());
          histos.fill(HIST("Analysis/h1d_rec_invm_kstar_anti"), p.M());
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

    fillDataHistos<false, false>(tracks, tracks, collision.spherocity(), collision.cent());
  }

  PROCESS_SWITCH(k892Analysis, processData, "Process for Same Event Data", true);

  void processMC(resoCols::iterator const& collision,
                 soa::Join<aod::ResoTracks, aod::ResoMCTracks> const& tracks)
  {
    histos.fill(HIST("Event/h1d_rec_sph"), collision.spherocity());
    fillDataHistos<false, true>(tracks, tracks, collision.spherocity(), collision.cent());
  }
  PROCESS_SWITCH(k892Analysis, processMC, "Process Event for MC", false);

  void processMCTrue(aod::ResoMCParents const& resoParents)
  {

    for (auto const& part : resoParents) {

      if (abs(part.pdgCode()) != 313) // // L* pdg_code = 313
        continue;
      if (abs(part.y()) > 0.5) { // rapidity cut
        continue;
      }

      bool pass1 = false;
      bool pass2 = false;

      if (abs(part.daughterPDG1()) == 211 || abs(part.daughterPDG2()) == 211) { // At least one decay to Pion
        pass1 = true;
      }
      if (abs(part.daughterPDG1()) == 321 || abs(part.daughterPDG2()) == 321) { // At least one decay to Kaon
        pass2 = true;
      }

      if (!pass1 || !pass2) // If we have both decay products
        continue;

      if (part.pdgCode() > 0)
        histos.fill(HIST("Analysis/h1d_gen_kstar"), part.pt());
      else
        histos.fill(HIST("Analysis/h1d_gen_kstar_anti"), part.pt());
    }
  }
  PROCESS_SWITCH(k892Analysis, processMCTrue, "Process Event for MC", false);

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
        fillDataHistos<true, false>(t1, t2, c1.spherocity(), c1.cent());
      }
    } else {
      SameKindPair<resoCols, resoTracks, BinningType2> pairs{binningPositions2, cNumMixEv, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip
      for (auto& [c1, t1, c2, t2] : pairs) {
        fillDataHistos<true, false>(t1, t2, c1.spherocity(), c1.cent());
      }
    }
  }

  PROCESS_SWITCH(k892Analysis, processMix, "Process for Mixed Events", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<k892Analysis>(cfgc)};
}
