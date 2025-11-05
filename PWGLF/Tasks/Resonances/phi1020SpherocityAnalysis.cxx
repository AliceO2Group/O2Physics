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

/// \file phi1020SpherocityAnalysis.cxx
/// \brief Invariant Mass Reconstruction of phi(1020) Resonance in K+K- channel and Spherocity dependence study in pp collision at 13.6 TeV.
/// \author Swadhin Behera <swadhin.behera@cern.ch> , Balaram Singh <balaram.singh@cern.ch>

#include "PWGLF/DataModel/LFResonanceTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include <array>
#include <random>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

struct Phi1020SpherocityAnalysis {

  // Named constants to avoid magic numbers
  static constexpr float kRapidityCut = 0.5F;
  static constexpr float kRotateMin = 0.01F;
  static constexpr float kRotateMax = 0.1F;
  static constexpr float kMaxMcPosZ = 10.0F;

  SliceCache cache;
  Preslice<aod::ResoTracks> perRCol = aod::resodaughter::resoCollisionId;
  Preslice<aod::Tracks> perCollision = aod::track::collisionId;

  // Configurables.
  Configurable<int> nBinsPt{"nBinsPt", 200, "N bins in pT histogram"};
  Configurable<int> nBinsInvM{"nBinsInvM", 200, "N bins in InvMass histogram"};
  Configurable<int> nBinsSp{"nBinsSp", 120, "N bins in spherocity histogram"};
  Configurable<bool> doRotate{"doRotate", true, "rotated inv mass spectra"};

  // Tracks
  Configurable<float> ptMin{"ptMin", 0.15, "Minimum Track pT"};
  Configurable<float> etaCut{"etaCut", 0.8, "Pseudorapidity cut"};
  Configurable<float> dcazMin{"dcazMin", 1., "Minimum DCAz"};
  Configurable<float> dcaxyMin{"dcaxyMin", 0.1, "Minimum DCAxy"};
  Configurable<bool> primaryTrack{"primaryTrack", true, "Primary track selection"};
  Configurable<bool> globalWoDCATrack{"globalWoDCATrack", true, "Global track selection without DCA"};
  Configurable<bool> pvContributor{"pvContributor", true, "PV contributor track selection"};

  // PID Selections (Kaons)
  Configurable<bool> useOnlyTOFTrackK{"useOnlyTOFTrackK", false, "Use only TOF track for PID selection (Kaon)"};
  Configurable<bool> useTpcOnly{"useTpcOnly", false, "Use TPC Only selection"};
  Configurable<float> rejNsigma{"rejNsigma", 1.0, "Reject tracks to improve purity of PID"};
  Configurable<double> maxTPCnSigmaKaon{"maxTPCnSigmaKaon", 3.0, "TPC nSigma cut for Kaon"};
  Configurable<double> maxTOFnSigmaKaon{"maxTOFnSigmaKaon", 3.0, "TOF nSigma cut for Kaon"};
  Configurable<double> nsigmaCutCombinedKaon{"nsigmaCutCombinedKaon", 3.0, "Combined nSigma cut for Kaon"};
  Configurable<std::vector<float>> kaonTPCPIDp{"kaonTPCPIDp", {0.0, 0.25, 0.5, 1.0}, "pT dependent TPC cuts kaons"};
  Configurable<std::vector<float>> kaonTPCPIDcut{"kaonTPCPIDcut", {6., 3.5, 2.5}, "TPC nsigma cuts kaons"};

  // Event Mixing.
  Configurable<bool> mixSph{"mixSph", true, "Include Sph Bins to be mixed"};
  Configurable<int> numMixEv{"numMixEv", 20, "Number of Events to be mixed"};
  ConfigurableAxis mixVtxBins{"mixVtxBins", {VARIABLE_WIDTH, -10.0f, -9.f, -8.f, -7.f, -6.f, -5.f, -4.f, -3.f, -2.f, -1.f, 0.f, 1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis mixMultBins{"mixMultBins", {VARIABLE_WIDTH, 0.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.0f, 200.0f}, "Mixing bins - multiplicity"};
  ConfigurableAxis mixSphBins{"mixSphBins", {VARIABLE_WIDTH, 0.0f, 0.2f, 0.4f, 0.6f, 0.8f, 1.0f}, "Mixing bins - spherocity"};

  // Histogram Registry.
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {
    const AxisSpec axisSp(nBinsSp, 0., 1., "S_{0}");
    const AxisSpec axisCent(105, 0, 105, "FT0M (%)");
    const AxisSpec axisPpid(400, 0., 4., "p (GeV/c)");
    const AxisSpec axisPtPid(400, 0., 4., "p_{T} (GeV/c)");
    const AxisSpec axisPt(nBinsPt, 0., 10., "p_{T} (GeV/c)");
    const AxisSpec axisEta(40, -1, 1, "#eta");
    const AxisSpec axisDCAz(500, -0.5, 0.5, {"DCA_{z} (cm)"});
    const AxisSpec axisDCAxy(240, -0.12, 0.12, {"DCA_{xy} (cm)"});
    const AxisSpec axisTPCNCls(200, 0, 200, {"TPCNCls"});
    const AxisSpec axisTPCNsigma(120, -6, 6, {"n#sigma^{TPC}"});
    const AxisSpec axisTOFNsigma(120, -6, 6, {"n#sigma^{TOF}"});
    const AxisSpec axisdEdx(380, 10, 200, {"#frac{dE}{dx}"});
    const AxisSpec axisInvM(nBinsInvM, 0.98, 1.06, {"M_{inv} (GeV/#it{c}^{2})"});

    // Event
    histos.add("Event/h1d_ft0m_mult_percentile", "FT0M (%)", kTH1F, {axisCent});
    histos.add("Event/h1d_spherocity", "Event Spherocity", kTH1F, {axisSp});
    histos.add("Event/h2d_sph_vs_multpercentile", "Spherocity vs FT0M(%)", kTH2F, {axisCent, axisSp});

    // QA Before
    histos.add("QAbefore/Kaon/h2d_ka_nsigma_tpc_p", "n#sigma^{TPC} Kaons", kTH2F, {axisPpid, axisTPCNsigma});
    histos.add("QAbefore/Kaon/h2d_ka_nsigma_tof_p", "n#sigma^{TOF} Kaons", kTH2F, {axisPpid, axisTOFNsigma});
    histos.add("QAbefore/Kaon/h2d_ka_nsigma_tof_vs_tpc", "n#sigma^{TPC} vs n#sigma^{TOF} Kaons", kTH2F, {axisTPCNsigma, axisTOFNsigma});

    // QA After
    histos.add("QAafter/Kaon/h1d_ka_pt", "p_{T}-spectra Kaons", kTH1F, {axisPt});
    histos.add("QAafter/Kaon/h2d_ka_dca_z", "dca_{z} Kaons", kTH2F, {axisPtPid, axisDCAz});
    histos.add("QAafter/Kaon/h2d_ka_dca_xy", "dca_{xy} Kaons", kTH2F, {axisPtPid, axisDCAxy});
    histos.add("QAafter/Kaon/h2d_ka_dEdx_p", "TPC Signal Kaon", kTH2F, {axisPpid, axisdEdx});
    histos.add("QAafter/Kaon/h2d_ka_nsigma_tpc_pt", " Kaons", kTH2F, {axisPtPid, axisTPCNsigma});
    histos.add("QAafter/Kaon/h2d_ka_nsigma_tpc_p", " Kaons", kTH2F, {axisPpid, axisTPCNsigma});
    histos.add("QAafter/Kaon/h2d_ka_nsigma_tof_pt", " Kaons", kTH2F, {axisPtPid, axisTOFNsigma});
    histos.add("QAafter/Kaon/h2d_ka_nsigma_tof_p", " Kaons", kTH2F, {axisPpid, axisTOFNsigma});
    histos.add("QAafter/Kaon/h2d_ka_nsigma_tof_vs_tpc", "n#sigma(TOF) vs n#sigma(TPC) Kaons", kTH2F, {axisTPCNsigma, axisTOFNsigma});

    // Analysis - phi(1020) -> K+ K-
    histos.add("Analysis/h1d_phi_invm_US", "phi(1020) M_{inv} unlike-sign", kTH1D, {axisInvM});
    histos.add("Analysis/h1d_phi_invm_PP", "Like Signs M_{inv} K^{+}K^{+}", kTH1D, {axisInvM});
    histos.add("Analysis/h1d_phi_invm_MM", "Like Signs M_{inv} K^{-}K^{-}", kTH1D, {axisInvM});
    histos.add("Analysis/h1d_phi_invm_rot", "Rotated Spectra", kTH1D, {axisInvM});
    histos.add("Analysis/h1d_phi_invm_US_mix", "Mixed Events M_{inv}", kTH1D, {axisInvM});
    histos.add("Analysis/h1d_phi_invm_LS_mix", "Mixed Events M_{inv}", kTH1D, {axisInvM});
    histos.add("Analysis/h4d_phi_invm_US", "THn #phi(1020)", kTHnSparseD, {axisInvM, axisPt, axisSp, axisCent});
    histos.add("Analysis/h4d_phi_invm_PP", "THn Like Signs K^{+}K^{+}", kTHnSparseD, {axisInvM, axisPt, axisSp, axisCent});
    histos.add("Analysis/h4d_phi_invm_MM", "THn Like Signs K^{-}K^{-}", kTHnSparseD, {axisInvM, axisPt, axisSp, axisCent});
    histos.add("Analysis/h4d_phi_invm_rot", "THn Rotated", kTHnSparseD, {axisInvM, axisPt, axisSp, axisCent});
    histos.add("Analysis/h4d_phi_invm_US_mix", "THn Mixed Events", kTHnSparseD, {axisInvM, axisPt, axisSp, axisCent});
    histos.add("Analysis/h4d_phi_invm_LS_mix", "THn Mixed Events", kTHnSparseD, {axisInvM, axisPt, axisSp, axisCent});

    // MC
    if (doprocessMC) {
      histos.add("Event/h1d_rec_sph", "Reconstructed S_{0}", kTH1F, {axisSp});
      histos.add("Analysis/h1d_gen_phi", "Generated phi(1020) p_{T}", kTH1D, {axisPt});
      histos.add("Analysis/h1d_rec_phi", "Reconstructed phi(1020) p_{T}", kTH1D, {axisPt});
      histos.add("Analysis/h1d_rec_invm_phi", "Recostructed phi(1020)", kTH1D, {axisInvM});
      // MC truth histograms for daughters
      histos.add("MCTruth/h1d_gen_posZ", "Generated PosZ", kTH1F, {{240, -12., 12.}});
      histos.add("MCTruth/h1d_ch_gen_phi", "Generated #phi distribution", kTH1F, {{128, -0.05, 6.35}});
      histos.add("MCTruth/h1d_ka_gen_eta", "Generated #eta Kaons", kTH1F, {{40, -1, 1}});
      histos.add("QAChecks/h1d_ka_gen_pt", "Generated p_{T}-spectra Kaons", kTH1F, {{400, 0., 4.}});
    }
  }

  template <typename T>
  bool selTracks(T const& track)
  {
    if (track.pt() < ptMin)
      return false;
    if (std::abs(track.eta()) > etaCut)
      return false;
    if (std::abs(track.dcaZ()) > dcazMin)
      return false;
    if (std::abs(track.dcaXY()) > dcaxyMin)
      return false;
    if (primaryTrack && !track.isPrimaryTrack())
      return false;
    if (globalWoDCATrack && !track.isGlobalTrackWoDCA())
      return false;
    if (pvContributor && !track.isPVContributor())
      return false;
    return true;
  }

  // PID selection tools (Kaon focused)
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
    float combinedRejCut = rejNsigma * rejNsigma;

    if (!useTpcOnly && candidate.hasTOF()) {
      if (tofNsigmaKa < maxTOFnSigmaKaon && tofNsigmaPi > rejNsigma && tofNsigmaPr > rejNsigma) {
        tofPIDPassed = true;
      }
      // square cut
      if ((nsigmaCutCombinedKaon < 0) && (tpcNsigmaKa < maxTPCnSigmaKaon)) {
        tpcPIDPassed = true;
      }
      // circular
      if ((nsigmaCutCombinedKaon > 0) && (tpcTofNsigmaKa < combinedCut && tpcTofNsigmaPi > combinedRejCut && tpcTofNsigmaPr > combinedRejCut)) {
        tofPIDPassed = true;
        tpcPIDPassed = true;
      }
    } else {
      tofPIDPassed = true;
      if (useTpcOnly) {
        if (tpcNsigmaKa < maxTPCnSigmaKaon && tpcNsigmaPi > rejNsigma && tpcNsigmaPr > rejNsigma) {
          tpcPIDPassed = true;
        }
      } else {
        for (int i = 0; i < nitr - 1; ++i) {
          if (p >= tpcPIDp[i] && p < tpcPIDp[i + 1] && (tpcNsigmaKa < tpcPIDcut[i] && tpcNsigmaPi > rejNsigma && tpcNsigmaPr > rejNsigma)) {
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
  void fillQAHistos(T const& track)
  {
    // get total momentum
    float p = RecoDecay::p(track.px(), track.py(), track.pz());

    // fill before QA
    histos.fill(HIST("QAbefore/Kaon/h2d_ka_nsigma_tpc_p"), p, track.tpcNSigmaKa());
    if (track.hasTOF()) {
      histos.fill(HIST("QAbefore/Kaon/h2d_ka_nsigma_tof_p"), p, track.tofNSigmaKa());
      histos.fill(HIST("QAbefore/Kaon/h2d_ka_nsigma_tof_vs_tpc"), track.tpcNSigmaKa(), track.tofNSigmaKa());
    }

    // fill after QA if it passes PID
    if (selectionPIDKaon(track, p)) {
      histos.fill(HIST("QAafter/Kaon/h1d_ka_pt"), track.pt());
      histos.fill(HIST("QAafter/Kaon/h2d_ka_dca_z"), track.pt(), track.dcaZ());
      histos.fill(HIST("QAafter/Kaon/h2d_ka_dca_xy"), track.pt(), track.dcaXY());
      histos.fill(HIST("QAafter/Kaon/h2d_ka_dEdx_p"), p, track.tpcSignal());
      histos.fill(HIST("QAafter/Kaon/h2d_ka_nsigma_tpc_p"), p, track.tpcNSigmaKa());
      histos.fill(HIST("QAafter/Kaon/h2d_ka_nsigma_tpc_pt"), track.pt(), track.tpcNSigmaKa());
      if (track.hasTOF() && !useTpcOnly) {
        histos.fill(HIST("QAafter/Kaon/h2d_ka_nsigma_tof_p"), p, track.tofNSigmaKa());
        histos.fill(HIST("QAafter/Kaon/h2d_ka_nsigma_tof_pt"), track.pt(), track.tofNSigmaKa());
        histos.fill(HIST("QAafter/Kaon/h2d_ka_nsigma_tof_vs_tpc"), track.tpcNSigmaKa(), track.tofNSigmaKa());
      }
    }
  }

  template <bool mix, bool mc, typename trackType>
  void fillInvMassHistos(trackType const& trk1, trackType const& trk2, float const& sph, float const& mult)
  {
    // use std random for rotation
    static thread_local std::mt19937 rng{std::random_device{}()};
    std::uniform_real_distribution<float> distRotate(kRotateMin, kRotateMax);
    float k1Ptot = 0., k2Ptot = 0.;

    for (auto const& pair : soa::combinations(soa::CombinationsFullIndexPolicy(trk1, trk2))) {
      auto const& trkK1 = std::get<0>(pair);
      auto const& trkK2 = std::get<1>(pair);
      // Do not analyse same index tracks.
      if (trkK1.index() == trkK2.index())
        continue;

      // pT, DCA, Global Tracks and PVcontrib selection.
      if (!selTracks(trkK1) || !selTracks(trkK2))
        continue;

      k1Ptot = RecoDecay::p(trkK1.px(), trkK1.py(), trkK1.pz());
      k2Ptot = RecoDecay::p(trkK2.px(), trkK2.py(), trkK2.pz());

      // Apply PID Selection
      if (useOnlyTOFTrackK && !trkK1.hasTOF())
        continue;
      if (useOnlyTOFTrackK && !trkK2.hasTOF())
        continue;
      if (!selectionPIDKaon(trkK1, k1Ptot) || !selectionPIDKaon(trkK2, k2Ptot))
        continue;

      // Invariant mass reconstruction using RecoDecay utilities.
      std::array<float, 3> mom1{trkK1.px(), trkK1.py(), trkK1.pz()};
      std::array<float, 3> mom2{trkK2.px(), trkK2.py(), trkK2.pz()};
      auto momTot = RecoDecay::pVec(mom1, mom2);
      double mass = RecoDecay::m(std::array{mom1, mom2}, std::array{MassKaonCharged, MassKaonCharged});
      double rap = RecoDecay::y(momTot, mass);
      double ptTot = RecoDecay::pt(momTot);

      // rapidity cut
      if (std::abs(rap) > kRapidityCut)
        continue;

      // Fill Invariant Mass Histograms.
      if constexpr (!mix && !mc) {
        if (trkK1.sign() * trkK2.sign() < 0) {
          histos.fill(HIST("Analysis/h1d_phi_invm_US"), mass);
          histos.fill(HIST("Analysis/h4d_phi_invm_US"), mass, ptTot, sph, mult);
          if (doRotate) {
            float theta = distRotate(rng);
            // rotate first daughter's transverse momentum
            float px1 = mom1[0];
            float py1 = mom1[1];
            float px1r = px1 * std::cos(theta) - py1 * std::sin(theta);
            float py1r = px1 * std::sin(theta) + py1 * std::cos(theta);
            std::array<float, 3> mom1r{px1r, py1r, mom1[2]};
            auto momTotR = RecoDecay::pVec(mom1r, mom2);
            double massR = RecoDecay::m(std::array{mom1r, mom2}, std::array{MassKaonCharged, MassKaonCharged});
            double rapR = RecoDecay::y(momTotR, massR);
            double ptR = RecoDecay::pt(momTotR);
            if (std::abs(rapR) < kRapidityCut) {
              histos.fill(HIST("Analysis/h1d_phi_invm_rot"), massR);
              histos.fill(HIST("Analysis/h4d_phi_invm_rot"), massR, ptR, sph, mult);
            }
          }
        } else {
          if (trkK1.sign() == 1) {
            histos.fill(HIST("Analysis/h1d_phi_invm_PP"), mass);
            histos.fill(HIST("Analysis/h4d_phi_invm_PP"), mass, ptTot, sph, mult);
          } else {
            histos.fill(HIST("Analysis/h1d_phi_invm_MM"), mass);
            histos.fill(HIST("Analysis/h4d_phi_invm_MM"), mass, ptTot, sph, mult);
          }
        }
      }

      if constexpr (mix) {
        if (trkK1.sign() * trkK2.sign() < 0) {
          histos.fill(HIST("Analysis/h1d_phi_invm_US_mix"), mass);
          histos.fill(HIST("Analysis/h4d_phi_invm_US_mix"), mass, ptTot, sph, mult);
        } else {
          histos.fill(HIST("Analysis/h1d_phi_invm_LS_mix"), mass);
          histos.fill(HIST("Analysis/h4d_phi_invm_LS_mix"), mass, ptTot, sph, mult);
        }
      }

      if constexpr (mc) {
        if (std::abs(trkK1.pdgCode()) != static_cast<int>(PDG_t::kKPlus) || std::abs(trkK2.pdgCode()) != static_cast<int>(PDG_t::kKPlus))
          continue;

        if (trkK1.motherId() != trkK2.motherId())
          continue;

        if (std::abs(trkK1.motherPDG()) != static_cast<int>(Pdg::kPhi)) // phi pdg
          continue;

        // MC histograms
        histos.fill(HIST("Analysis/h1d_rec_phi"), ptTot);
        histos.fill(HIST("Analysis/h1d_rec_invm_phi"), mass);
      }
    }
  }

  using ResoCols = soa::Join<aod::ResoCollisions, aod::ResoSpheroCollisions>;
  using ResoTracks = aod::ResoTracks;

  void processData(ResoCols::iterator const& collision, ResoTracks const& tracks)
  {
    histos.fill(HIST("Event/h1d_ft0m_mult_percentile"), collision.cent());
    histos.fill(HIST("Event/h1d_spherocity"), collision.spherocity());
    histos.fill(HIST("Event/h2d_sph_vs_multpercentile"), collision.cent(), collision.spherocity());

    // QA per track
    for (auto const& track : tracks) {
      if (!selTracks(track))
        continue;
      fillQAHistos(track);
    }

    // get invariant mass histograms
    fillInvMassHistos<false, false>(tracks, tracks, collision.spherocity(), collision.cent());
  }

  PROCESS_SWITCH(Phi1020SpherocityAnalysis, processData, "Process for Same Event Data", true);

  void processMC(ResoCols::iterator const& collision,
                 soa::Join<aod::ResoTracks, aod::ResoMCTracks> const& tracks)
  {
    histos.fill(HIST("Event/h1d_rec_sph"), collision.spherocity());

    // get MC reco pT-spectra and QA
    for (auto const& track : tracks) {
      if (!selTracks(track))
        continue;
      fillQAHistos(track);
    }

    // get invariant mass histograms
    fillInvMassHistos<false, true>(tracks, tracks, collision.spherocity(), collision.cent());
  }
  PROCESS_SWITCH(Phi1020SpherocityAnalysis, processMC, "Process Event for MC", false);

  void processMCTrue(aod::ResoMCParents const& resoParents)
  {
    for (auto const& part : resoParents) {
      if (std::abs(part.pdgCode()) != static_cast<int>(Pdg::kPhi)) // phi pdg
        continue;
      if (std::abs(part.y()) > kRapidityCut)
        continue;

      bool passKa1 = false;
      bool passKa2 = false;
      if (std::abs(part.daughterPDG1()) == static_cast<int>(PDG_t::kKPlus) || std::abs(part.daughterPDG2()) == static_cast<int>(PDG_t::kKPlus))
        passKa1 = true;
      if (std::abs(part.daughterPDG1()) == static_cast<int>(PDG_t::kKPlus) || std::abs(part.daughterPDG2()) == static_cast<int>(PDG_t::kKPlus))
        passKa2 = true;

      if (!passKa1 || !passKa2)
        continue;

      histos.fill(HIST("Analysis/h1d_gen_phi"), part.pt());
    }
  }
  PROCESS_SWITCH(Phi1020SpherocityAnalysis, processMCTrue, "Process Event for MC", false);

  void processMCTrueDaughters(aod::McCollisions::iterator const& McCollision, aod::McParticles const& McParts)
  {
    // IP range selection
    if (std::abs(McCollision.posZ()) > kMaxMcPosZ)
      return;

    histos.fill(HIST("MCTruth/h1d_gen_posZ"), McCollision.posZ());

    for (auto const& part : McParts) {
      // kinematic acceptance of particles
      if (part.pt() < ptMin || std::abs(part.eta()) > etaCut || !part.isPhysicalPrimary())
        continue;

      histos.fill(HIST("MCTruth/h1d_ch_gen_phi"), part.phi());

      if (std::abs(part.pdgCode()) == static_cast<int>(PDG_t::kKPlus)) {
        histos.fill(HIST("MCTruth/h1d_ka_gen_eta"), part.eta());
        histos.fill(HIST("QAChecks/h1d_ka_gen_pt"), part.pt());
      }
    }
  }
  PROCESS_SWITCH(Phi1020SpherocityAnalysis, processMCTrueDaughters, "Process Event for MC truth of kaons", false);

  // Processing Event Mixing
  using BinningType1 = ColumnBinningPolicy<aod::collision::PosZ, aod::resocollision::Cent, aod::resocollision::Spherocity>;
  using BinningType2 = ColumnBinningPolicy<aod::collision::PosZ, aod::resocollision::Cent>;
  void processMix(ResoCols const& collisions, ResoTracks const& tracks)
  {
    LOGF(debug, "Event Mixing Started");
    BinningType1 binningPositions1{{mixVtxBins, mixMultBins, mixSphBins}, true};
    BinningType2 binningPositions2{{mixVtxBins, mixMultBins}, true};
    auto tracksTuple = std::make_tuple(tracks);
    if (mixSph) {
      SameKindPair<ResoCols, ResoTracks, BinningType1> pairs{binningPositions1, numMixEv, -1, collisions, tracksTuple, &cache};
      for (auto const& item : pairs) {
        auto const& c1 = std::get<0>(item);
        auto const& t1 = std::get<1>(item);
        auto const& t2 = std::get<3>(item);
        fillInvMassHistos<true, false>(t1, t2, c1.spherocity(), c1.cent());
      }
    } else {
      SameKindPair<ResoCols, ResoTracks, BinningType2> pairs{binningPositions2, numMixEv, -1, collisions, tracksTuple, &cache};
      for (auto const& item : pairs) {
        auto const& c1 = std::get<0>(item);
        auto const& t1 = std::get<1>(item);
        auto const& t2 = std::get<3>(item);
        fillInvMassHistos<true, false>(t1, t2, c1.spherocity(), c1.cent());
      }
    }
  }
  PROCESS_SWITCH(Phi1020SpherocityAnalysis, processMix, "Process for Mixed Events", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<Phi1020SpherocityAnalysis>(cfgc)};
}
