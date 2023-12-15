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
/// \brief Produce Spherocity table.
///        Invariant Mass Reconstruction of Lambda(1520) Resonance.
///
/// \author Yash Patley <yash.patley@cern.ch>

#include <TLorentzVector.h>
#include <TRandom.h>

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

struct lambdaAnalysis {

  SliceCache cache;
  Preslice<aod::ResoTracks> perRCol = aod::resodaughter::resoCollisionId;
  Preslice<aod::Tracks> perCollision = aod::track::collisionId;

  // Configurables.
  Configurable<int> nBinsPt{"nBinsPt", 100, "N bins in pT histogram"};
  Configurable<int> nBinsInvM{"nBinsInvM", 120, "N bins in InvMass histogram"};
  Configurable<int> nBinsSp{"nBinsSp", 100, "N bins in spherocity histogram"};
  Configurable<bool> doRotate{"doRotate", true, "rotated inv mass spectra"};

  // Tracks
  Configurable<float> cfgPtMin{"ptMin", 0.15, "Minimum Track pT"};
  Configurable<float> cfgEtaCut{"etaCut", 0.8, "Pseudorapidity cut"};
  Configurable<float> cfgDcaz{"dcazMin", 1., "Minimum DCAz"};
  Configurable<float> cfgDcaxy{"dcaxyMin", 0.1, "Minimum DCAxy"};
  Configurable<float> cfgPIDprecut{"cfgPIDprecut", 5, "Preselection PID TPC TOF cut"};
  Configurable<bool> cfgKinCuts{"cfgKinCuts", false, "Kinematic Cuts for p-K pair opening angle"};
  Configurable<bool> cfgPrimaryTrack{"cfgPrimaryTrack", true, "Primary track selection"};                    // kGoldenChi2 | kDCAxy | kDCAz
  Configurable<bool> cfgGlobalWoDCATrack{"cfgGlobalWoDCATrack", true, "Global track selection without DCA"}; // kQualityTracks (kTrackType | kTPCNCls | kTPCCrossedRows | kTPCCrossedRowsOverNCls | kTPCChi2NDF | kTPCRefit | kITSNCls | kITSChi2NDF | kITSRefit | kITSHits) | kInAcceptanceTracks (kPtRange | kEtaRange)
  Configurable<bool> cfgPVContributor{"cfgPVContributor", true, "PV contributor track selection"};           // PV Contriuibutor

  // PID Selections
  Configurable<bool> cUseOnlyTOFTrackPr{"cUseOnlyTOFTrackPr", false, "Use only TOF track for PID selection"}; // Use only TOF track for Proton PID selection
  Configurable<bool> cUseOnlyTOFTrackKa{"cUseOnlyTOFTrackKa", false, "Use only TOF track for PID selection"}; // Use only TOF track for Kaon PID selection
  Configurable<bool> cUseTpcAndTof{"cUseTpcAndTof", true, "Use TPC and TOF PID selection"};                   // TPC And TOF tracks
  Configurable<bool> cUseTpcOnly{"cUseTpcOnly", false, "Use TPC Only tracks (No TOF Veto)"};                  // TPC only selection
  Configurable<float> cRejNsigma{"cRejNsigma", 3.0, "Reject tracks to improve purity of PID"};
  Configurable<float> cMaxPtpc{"cMaxPtpc", 1.8, "Maximum p for tracks without TOF"};
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
  Configurable<bool> doSphMix{"doSphMix", true, "Include Sph Bins to be mixed"};
  Configurable<int> nMix{"nMix", 10, "Number of Events to be mixed"};
  ConfigurableAxis cfgVtxBins{"cfgVtxBins", {VARIABLE_WIDTH, -10.0f, -9.f, -8.f, -7.f, -6.f, -5.f, -4.f, -3.f, -2.f, -1.f, 0.f, 1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis cfgMultBins{"cfgMultBins", {VARIABLE_WIDTH, 0.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.0f, 200.0f}, "Mixing bins - multiplicity"};
  ConfigurableAxis cfgSphBins{"cfgSphBins", {VARIABLE_WIDTH, 0.0f, 0.2f, 0.4f, 0.6f, 0.8f, 1.0f}, "Mixing bins - spherocity"};

  // Histogram Registry.
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {

    // Define Axis.
    const AxisSpec axisSp(nBinsSp, 0., 1., "S_{0}");
    const AxisSpec axisCent(105, 0, 105, "FT0M (%)");
    const AxisSpec axisP_pid(400, 0., 4., "p (GeV/c)");
    const AxisSpec axisPt(nBinsPt, 0., 10., "p_{T} (GeV/c)");
    const AxisSpec axisEta(40, -1, 1, "#eta");
    const AxisSpec axisDCAz(500, -0.5, 0.5, {"DCA_{z} (cm)"});
    const AxisSpec axisDCAxy(240, -0.12, 0.12, {"DCA_{xy} (cm)"});
    const AxisSpec axisTPCNCls(200, 0, 200, {"TPCNCls"});
    const AxisSpec axisTPCNsigma(120, -6, 6, {"n#sigma^{TPC}"});
    const AxisSpec axisTOFNsigma(120, -6, 6, {"n#sigma^{TOF}"});
    const AxisSpec axisInvM(nBinsInvM, 1.44, 2.44, {"M_{inv} (GeV/c^{2})"});

    // Create Histograms.
    // Event
    histos.add("Event/hCent", "FT0M (%)", kTH1F, {axisCent});
    histos.add("Event/hSph", "Event Spherocity", kTH1F, {axisSp});
    histos.add("Event/hSpCent", "Spherocity vs FT0M(%)", kTH2F, {axisCent, axisSp});

    // QA Before
    histos.add("QAbefore/Proton/hTPCNsigma", "n#sigma^{TPC} Protons", kTH2F, {axisP_pid, axisTPCNsigma});
    histos.add("QAbefore/Proton/hTOFNsigma", "n#sigma^{TOF} Protons", kTH2F, {axisP_pid, axisTOFNsigma});
    histos.add("QAbefore/Proton/hTpcTofNsigma", "n#sigma^{TPC} vs n#sigma^{TOF} Protons", kTH2F, {axisTPCNsigma, axisTOFNsigma});
    histos.add("QAbefore/Kaon/hTPCNsigma", "n#sigma^{TPC} Kaons", kTH2F, {axisP_pid, axisTPCNsigma});
    histos.add("QAbefore/Kaon/hTOFNsigma", "n#sigma^{TOF} Kaons", kTH2F, {axisP_pid, axisTOFNsigma});
    histos.add("QAbefore/Kaon/hTpcTofNsigma", "n#sigma^{TPC} vs n#sigma^{TOF} Kaons", kTH2F, {axisTPCNsigma, axisTOFNsigma});

    // QA After
    histos.add("QAafter/Proton/hPt", "p_{T}-spectra Protons", kTH1F, {axisPt});
    histos.add("QAafter/Proton/hDcaZ", "dca_{z} Protons", kTH2F, {axisP_pid, axisDCAz});
    histos.add("QAafter/Proton/hDcaXY", "dca_{xy} Protons", kTH2F, {axisP_pid, axisDCAxy});
    histos.add("QAafter/Proton/hTPCNsigmaFull", "n#sigma(TPC) Protons", kTH2F, {axisP_pid, axisTPCNsigma});
    histos.add("QAafter/Proton/hTPCNsigma", "n#sigma(TPC) Protons", kTH2F, {axisP_pid, axisTPCNsigma});
    histos.add("QAafter/Proton/hTPCNsigmaTOF", "n#sigma(TPC) Protons", kTH2F, {axisP_pid, axisTPCNsigma});
    histos.add("QAafter/Proton/hTOFNsigma", "n#sigma(TOF) Protons", kTH2F, {axisP_pid, axisTOFNsigma});
    histos.add("QAafter/Proton/hTpcTofNsigma", "n#sigma(TOF) vs n#sigma(TPC) Protons", kTH2F, {axisTPCNsigma, axisTOFNsigma});
    histos.add("QAafter/Kaon/hPt", "p_{T}-spectra Kaons", kTH1F, {axisPt});
    histos.add("QAafter/Kaon/hDcaZ", "dca_{z} Kaons", kTH2F, {axisP_pid, axisDCAz});
    histos.add("QAafter/Kaon/hDcaXY", "dca_{xy} Kaons", kTH2F, {axisP_pid, axisDCAxy});
    histos.add("QAafter/Kaon/hTPCNsigmaFull", "n#sigma(TPC) Kaons", kTH2F, {axisP_pid, axisTPCNsigma});
    histos.add("QAafter/Kaon/hTPCNsigma", "n#sigma(TPC) Kaons", kTH2F, {axisP_pid, axisTPCNsigma});
    histos.add("QAafter/Kaon/hTPCNsigmaTOF", "n#sigma(TPC) Kaons", kTH2F, {axisP_pid, axisTPCNsigma});
    histos.add("QAafter/Kaon/hTOFNsigma", "n#sigma(TOF) Kaons", kTH2F, {axisP_pid, axisTOFNsigma});
    histos.add("QAafter/Kaon/hTpcTofNsigma", "n#sigma(TOF) vs n#sigma(TPC) Kaons", kTH2F, {axisTPCNsigma, axisTOFNsigma});

    // Analysis
    // Lambda Invariant Mass
    histos.add("Analysis/hInvMass", "#Lambda(1520) M_{inv}", kTH1D, {axisInvM});
    histos.add("Analysis/hInvMassLS1", "Like Signs M_{inv} p K^{+}", kTH1D, {axisInvM});
    histos.add("Analysis/hInvMassLS2", "Like Signs M_{inv} #bar{p} K^{-}", kTH1D, {axisInvM});
    histos.add("Analysis/hInvMassR", "Rotated Spectra", kTH1D, {axisInvM});
    histos.add("Analysis/hInvMassMix", "Mixed Events M_{inv}", kTH1D, {axisInvM});
    histos.add("Analysis/hInvMassMixLS", "Mixed Events M_{inv}", kTH1D, {axisInvM});
    histos.add("Analysis/h4InvMass", "THn #Lambda(1520)", kTHnSparseD, {axisInvM, axisPt, axisSp, axisCent});
    histos.add("Analysis/h4InvMassLS1", "THn Like Signs p K^{+}", kTHnSparseD, {axisInvM, axisPt, axisSp, axisCent});
    histos.add("Analysis/h4InvMassLS2", "THn Like Signs #bar{p} K^{-}", kTHnSparseD, {axisInvM, axisPt, axisSp, axisCent});
    histos.add("Analysis/h4InvMassR", "THn Rotated", kTHnSparseD, {axisInvM, axisPt, axisSp, axisCent});
    histos.add("Analysis/h4InvMassMix", "THn Mixed Events", kTHnSparseD, {axisInvM, axisPt, axisSp, axisCent});
    histos.add("Analysis/h4InvMassMixLS", "THn Mixed Events", kTHnSparseD, {axisInvM, axisPt, axisSp, axisCent});

    // MC
    if (doprocessMC) {
      histos.add("Event/hSphRec", "Reconstructed S_{0}", kTH1F, {axisSp});
      histos.add("Event/hSpCentRec", "Reconstructed S_{0} vs FT0M(%)", kTH2F, {axisCent, axisSp});
      histos.add("QAMCTrue/DcaZ_pr", "dca_{z}^{MC} Protons", kTH2F, {axisP_pid, axisDCAz});
      histos.add("QAMCTrue/DcaZ_ka", "dca_{z}^{MC} Kaons", kTH2F, {axisP_pid, axisDCAz});
      histos.add("QAMCTrue/DcaXY_pr", "dca_{xy}^{MC} Protons", kTH2F, {axisP_pid, axisDCAxy});
      histos.add("QAMCTrue/DcaXY_ka", "dca_{xy}^{MC} Kaons", kTH2F, {axisP_pid, axisDCAxy});
      histos.add("Analysis/hLambdaGen", "Generated #Lambda(1520) p_{T}", kTH1D, {axisPt});
      histos.add("Analysis/hLambdaGenAnti", "Generated #bar{#Lambda}(1520) p_{T}", kTH1D, {axisPt});
      histos.add("Analysis/hLambdaRec", "Reconstructed #Lambda(1520) p_{T}", kTH1D, {axisPt});
      histos.add("Analysis/hLambdaRecAnti", "Reconstructed #bar{#Lambda}(1520) p_{T}", kTH1D, {axisPt});
      histos.add("Analysis/hInvMassLambdaRec", "Recostructed #Lambda(1520)", kTH1D, {axisInvM});
      histos.add("Analysis/hInvMassLambdaRecAnti", "Recostructed #bar{#Lambda}(1520)", kTH1D, {axisInvM});
      histos.add("Analysis/h4InvMassLambdaRec", "Recostructed #Lambda(1520)", kTHnSparseD, {axisInvM, axisPt, axisSp, axisCent});
      histos.add("Analysis/h4InvMassLambdaRecAnti", "Recostructed #bar{#Lambda}(1520)", kTHnSparseD, {axisInvM, axisPt, axisSp, axisCent});
    }
  }

  template <typename T>
  bool selTracks(T const& track)
  {

    if (track.pt() < cfgPtMin)
      return false;

    if (std::abs(track.dcaZ()) > cfgDcaz)
      return false;

    if (std::abs(track.dcaXY()) > cfgDcaxy)
      return false;

    if (cfgPrimaryTrack && !track.isPrimaryTrack())
      return false;

    if (cfgGlobalWoDCATrack && !track.isGlobalTrackWoDCA())
      return false;

    if (cfgPVContributor && !track.isPVContributor())
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
    float combinedRejCut = cRejNsigma * cRejNsigma;

    if (!cUseTpcOnly && candidate.hasTOF()) {
      if (tofNsigmaPr < cMaxTOFnSigmaProton && tofNsigmaPi > cRejNsigma && tofNsigmaKa > cRejNsigma) {
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
        if (tpcNsigmaPr < cMaxTPCnSigmaProton && tpcNsigmaPi > cRejNsigma && tpcNsigmaKa > cRejNsigma) {
          tpcPIDPassed = true;
        }
      } else {
        for (int i = 0; i < nitr - 1; ++i) {
          if (p >= tpcPIDp[i] && p < tpcPIDp[i + 1] && (tpcNsigmaPr < tpcPIDcut[i] && tpcNsigmaPi > cRejNsigma && tpcNsigmaKa > cRejNsigma)) {
            tpcPIDPassed = true;
          }
          if (!cUseTpcAndTof && p >= tpcPIDp[nitr - 1] && p < cMaxPtpc && (tpcNsigmaPr < tpcPIDcut[nitr - 2] && tpcNsigmaPi > cRejNsigma && tpcNsigmaKa > cRejNsigma)) {
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
          if (!cUseTpcAndTof && p >= tpcPIDp[nitr - 1] && p < cMaxPtpc && (tpcNsigmaKa < tpcPIDcut[nitr - 2] && tpcNsigmaPi > cRejNsigma && tpcNsigmaPr > cRejNsigma)) {
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

    for (auto const& [trkPr, trkKa] : soa::combinations(soa::CombinationsFullIndexPolicy(trk1, trk2))) {
      // Do not analyse same index tracks.
      if (trkPr.index() == trkKa.index())
        continue;

      // pT, DCA, Global Tracks and PVcontrib selection.
      if (!selTracks(trkPr) || !selTracks(trkKa))
        continue;

      // PID preselection
      if ((std::abs(trkPr.tpcNSigmaPr()) > cfgPIDprecut) || (std::abs(trkKa.tpcNSigmaKa()) > cfgPIDprecut))
        continue;

      if ((trkPr.hasTOF() && std::abs(trkPr.tofNSigmaPr()) > cfgPIDprecut) || (trkKa.hasTOF() && std::abs(trkKa.tofNSigmaKa()) > cfgPIDprecut))
        continue;

      p_ptot = TMath::Sqrt(trkPr.px() * trkPr.px() + trkPr.py() * trkPr.py() + trkPr.py() * trkPr.py());
      k_ptot = TMath::Sqrt(trkKa.px() * trkKa.px() + trkKa.py() * trkKa.py() + trkKa.py() * trkKa.py());

      // Fill QA before track selection.
      if (!mix) {
        histos.fill(HIST("QAbefore/Proton/hTPCNsigma"), p_ptot, trkPr.tpcNSigmaPr());
        if (trkPr.hasTOF()) {
          histos.fill(HIST("QAbefore/Proton/hTOFNsigma"), p_ptot, trkPr.tofNSigmaPr());
          histos.fill(HIST("QAbefore/Proton/hTpcTofNsigma"), trkPr.tpcNSigmaPr(), trkPr.tofNSigmaPr());
        }
        histos.fill(HIST("QAbefore/Kaon/hTPCNsigma"), k_ptot, trkKa.tpcNSigmaKa());
        if (trkKa.hasTOF()) {
          histos.fill(HIST("QAbefore/Kaon/hTOFNsigma"), k_ptot, trkKa.tofNSigmaKa());
          histos.fill(HIST("QAbefore/Kaon/hTpcTofNsigma"), trkKa.tpcNSigmaKa(), trkKa.tofNSigmaKa());
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
      if (!mix) {
        histos.fill(HIST("QAafter/Proton/hPt"), trkPr.pt());
        histos.fill(HIST("QAafter/Proton/hDcaZ"), p_ptot, trkPr.dcaZ());
        histos.fill(HIST("QAafter/Proton/hDcaXY"), p_ptot, trkPr.dcaXY());
        histos.fill(HIST("QAafter/Proton/hTPCNsigmaFull"), p_ptot, trkPr.tpcNSigmaPr());
        if (!cUseTpcOnly && trkPr.hasTOF()) {
          histos.fill(HIST("QAafter/Proton/hTPCNsigmaTOF"), p_ptot, trkPr.tpcNSigmaPr());
          histos.fill(HIST("QAafter/Proton/hTOFNsigma"), p_ptot, trkPr.tofNSigmaPr());
          histos.fill(HIST("QAafter/Proton/hTpcTofNsigma"), trkPr.tpcNSigmaPr(), trkPr.tofNSigmaPr());
        } else {
          histos.fill(HIST("QAafter/Proton/hTPCNsigma"), p_ptot, trkPr.tpcNSigmaPr());
        }
        histos.fill(HIST("QAafter/Kaon/hPt"), trkKa.pt());
        histos.fill(HIST("QAafter/Kaon/hDcaZ"), k_ptot, trkKa.dcaZ());
        histos.fill(HIST("QAafter/Kaon/hDcaXY"), k_ptot, trkKa.dcaXY());
        histos.fill(HIST("QAafter/Kaon/hTPCNsigmaFull"), k_ptot, trkKa.tpcNSigmaKa());
        if (!cUseTpcOnly && trkKa.hasTOF()) {
          histos.fill(HIST("QAafter/Kaon/hTPCNsigmaTOF"), k_ptot, trkKa.tpcNSigmaKa());
          histos.fill(HIST("QAafter/Kaon/hTOFNsigma"), k_ptot, trkKa.tofNSigmaKa());
          histos.fill(HIST("QAafter/Kaon/hTpcTofNsigma"), trkKa.tpcNSigmaKa(), trkKa.tofNSigmaKa());
        } else {
          histos.fill(HIST("QAafter/Kaon/hTPCNsigma"), k_ptot, trkKa.tpcNSigmaKa());
        }
      }

      // Invariant mass reconstruction.
      p1.SetXYZM(trkPr.px(), trkPr.py(), trkPr.pz(), MassProton);
      p2.SetXYZM(trkKa.px(), trkKa.py(), trkKa.pz(), MassKaonCharged);
      p = p1 + p2;

      if (std::abs(p.Rapidity()) > 0.5)
        continue;

      // Apply kinematic cuts.
      if (cfgKinCuts) {
        TVector3 v1(trkPr.px(), trkPr.py(), trkPr.pz());
        TVector3 v2(trkKa.px(), trkKa.py(), trkKa.pz());
        float alpha = v1.Angle(v2);
        if (alpha > 1.4 && alpha < 2.4)
          continue;
      }

      // Fill Invariant Mass Histograms.
      if constexpr (!mix && !mc) {
        if (trkPr.sign() * trkKa.sign() < 0) {
          histos.fill(HIST("Analysis/hInvMass"), p.M());
          histos.fill(HIST("Analysis/h4InvMass"), p.M(), p.Pt(), sph, mult);
          if (doRotate) {
            float theta = rn->Uniform(3.12, 3.16);
            p1.RotateZ(theta);
            p = p1 + p2;
            if (std::abs(p.Rapidity()) < 0.5) {
              histos.fill(HIST("Analysis/hInvMassR"), p.M());
              histos.fill(HIST("Analysis/h4InvMassR"), p.M(), p.Pt(), sph, mult);
            }
          }
        } else {
          if (trkPr.sign() == 1) {
            histos.fill(HIST("Analysis/hInvMassLS1"), p.M());
            histos.fill(HIST("Analysis/h4InvMassLS1"), p.M(), p.Pt(), sph, mult);
          } else {
            histos.fill(HIST("Analysis/hInvMassLS2"), p.M());
            histos.fill(HIST("Analysis/h4InvMassLS2"), p.M(), p.Pt(), sph, mult);
          }
        }
      }

      if constexpr (mc) {
        if (abs(trkPr.pdgCode()) != 2212 || abs(trkKa.pdgCode()) != 321)
          continue;

        if (trkPr.motherId() != trkKa.motherId())
          continue;

        if (abs(trkPr.motherPDG()) != 3124)
          continue;

        // Track selection check.
        histos.fill(HIST("QAMCTrue/DcaXY_pr"), p_ptot, trkPr.dcaXY());
        histos.fill(HIST("QAMCTrue/DcaXY_ka"), k_ptot, trkKa.dcaXY());
        histos.fill(HIST("QAMCTrue/DcaZ_pr"), p_ptot, trkPr.dcaZ());
        histos.fill(HIST("QAMCTrue/DcaZ_ka"), k_ptot, trkKa.dcaZ());

        // MC histograms
        if (trkPr.motherPDG() > 0) {
          histos.fill(HIST("Analysis/hLambdaRec"), p.Pt());
          histos.fill(HIST("Analysis/hInvMassLambdaRec"), p.M());
          histos.fill(HIST("Analysis/h4InvMassLambdaRec"), p.M(), p.Pt(), sph, mult);
        } else {
          histos.fill(HIST("Analysis/hLambdaRecAnti"), p.Pt());
          histos.fill(HIST("Analysis/hInvMassLambdaRecAnti"), p.M());
          histos.fill(HIST("Analysis/h4InvMassLambdaRecAnti"), p.M(), p.Pt(), sph, mult);
        }
      }

      if constexpr (mix) {
        if (trkPr.sign() * trkKa.sign() < 0) {
          histos.fill(HIST("Analysis/hInvMassMix"), p.M());
          histos.fill(HIST("Analysis/h4InvMassMix"), p.M(), p.Pt(), sph, mult);
        } else {
          histos.fill(HIST("Analysis/hInvMassMixLS"), p.M());
          histos.fill(HIST("Analysis/h4InvMassMixLS"), p.M(), p.Pt(), sph, mult);
        }
      }
    }
  }

  using resoCols = aod::ResoCollisions;
  using resoTracks = aod::ResoTracks;

  void processData(resoCols::iterator const& collision, resoTracks const& tracks)
  {

    histos.fill(HIST("Event/hCent"), collision.cent());
    histos.fill(HIST("Event/hSph"), collision.spherocity());
    histos.fill(HIST("Event/hSpCent"), collision.cent(), collision.spherocity());

    fillDataHistos<false, false>(tracks, tracks, collision.spherocity(), collision.cent());
  }

  PROCESS_SWITCH(lambdaAnalysis, processData, "Process for Same Event Data", true);

  void processMC(resoCols::iterator const& collision,
                 soa::Join<aod::ResoTracks, aod::ResoMCTracks> const& tracks)
  {

    histos.fill(HIST("Event/hSphRec"), collision.spherocity());
    histos.fill(HIST("Event/hSpCentRec"), collision.cent(), collision.spherocity());
    fillDataHistos<false, true>(tracks, tracks, collision.spherocity(), collision.cent());
  }
  PROCESS_SWITCH(lambdaAnalysis, processMC, "Process Event for MC", false);

  void processMCTrue(aod::ResoMCParents const& resoParents)
  {

    for (auto const& part : resoParents) {

      if (abs(part.pdgCode()) != 3124) // #Lambda(1520) PDGCode = 3124
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
        histos.fill(HIST("Analysis/hLambdaGen"), part.pt());
      else
        histos.fill(HIST("Analysis/hLambdaGenAnti"), part.pt());
    }
  }
  PROCESS_SWITCH(lambdaAnalysis, processMCTrue, "Process Event for MC", false);

  // Processing Event Mixing
  using BinningType1 = ColumnBinningPolicy<aod::collision::PosZ, aod::resocollision::Cent, aod::resocollision::Spherocity>;
  using BinningType2 = ColumnBinningPolicy<aod::collision::PosZ, aod::resocollision::Cent>;
  void processMix(resoCols& collisions, resoTracks const& tracks)
  {

    LOGF(debug, "Event Mixing Started");
    BinningType1 binningPositions1{{cfgVtxBins, cfgMultBins, cfgSphBins}, true};
    BinningType2 binningPositions2{{cfgVtxBins, cfgMultBins}, true};
    auto tracksTuple = std::make_tuple(tracks);
    if (doSphMix) {
      SameKindPair<resoCols, resoTracks, BinningType1> pairs{binningPositions1, nMix, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip
      for (auto& [c1, t1, c2, t2] : pairs) {
        fillDataHistos<true, false>(t1, t2, c1.spherocity(), c1.cent());
      }
    } else {
      SameKindPair<resoCols, resoTracks, BinningType2> pairs{binningPositions2, nMix, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip
      for (auto& [c1, t1, c2, t2] : pairs) {
        fillDataHistos<true, false>(t1, t2, c1.spherocity(), c1.cent());
      }
    }
  }

  PROCESS_SWITCH(lambdaAnalysis, processMix, "Process for Mixed Events", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<lambdaAnalysis>(cfgc)};
}
