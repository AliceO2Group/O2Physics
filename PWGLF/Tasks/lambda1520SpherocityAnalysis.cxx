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
  // Proton
  Configurable<double> cMaxTPCnSigmaProton{"cMaxTPCnSigmaProton", 3.0, "TPC nSigma cut for Proton"};              // TPC
  Configurable<double> cMaxTOFnSigmaProton{"cMaxTOFnSigmaProton", 3.0, "TOF nSigma cut for Proton"};              // TOF
  Configurable<double> nsigmaCutCombinedProton{"nsigmaCutCombinedProton", 3.0, "Combined nSigma cut for Proton"}; // Combined
  Configurable<std::vector<float>> protonTPCPIDpt{"protonTPCPIDpt", {0, 0.5, 0.7, 0.8}, "pT dependent TPC cuts protons"};
  Configurable<std::vector<float>> protonTPCPIDcut{"protonTPCPIDcut", {5., 3.5, 2.5}, "TPC nsigma cuts protons"};
  // Kaon
  Configurable<double> cMaxTPCnSigmaKaon{"cMaxTPCnSigmaKaon", 3.0, "TPC nSigma cut for Kaon"};              // TPC
  Configurable<double> cMaxTOFnSigmaKaon{"cMaxTOFnSigmaKaon", 3.0, "TOF nSigma cut for Kaon"};              // TOF
  Configurable<double> nsigmaCutCombinedKaon{"nsigmaCutCombinedKaon", 3.0, "Combined nSigma cut for Kaon"}; // Combined
  Configurable<std::vector<float>> kaonTPCPIDpt{"kaonTPCPIDpt", {0., 0.25, 0.3, 0.45}, "pT dependent TPC cuts kaons"};
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
    const AxisSpec axisPtQA(40, 0., 2., "p_{T} (GeV/c)");
    const AxisSpec axisPt(nBinsPt, 0., 10., "p_{T} (GeV/c)");
    const AxisSpec axisEta(40, -1, 1, "#eta");
    const AxisSpec axisDCAz(500, -0.5, 0.5, {"DCA_{z} (cm)"});
    const AxisSpec axisDCAxy(240, -0.12, 0.12, {"DCA_{xy} (cm)"});
    const AxisSpec axisTPCNCls(200, 0, 200, {"TPCNCls"});
    const AxisSpec axisTPCNsigma(120, -6, 6, {"n#sigma^{TPC}"});
    const AxisSpec axisTOFNsigma(120, -6, 6, {"n#sigma^{TOF}"});
    const AxisSpec axisInvM(nBinsInvM, 1.44, 2.04, {"M_{inv} (GeV/c^{2})"});

    // Create Histograms.
    // Event
    histos.add("Event/hCent", "FT0M (%)", kTH1F, {axisCent});
    histos.add("Event/hSph", "Event Spherocity", kTH1F, {axisSp});
    histos.add("Event/hSpCent", "Spherocity vs FT0M(%)", kTH2F, {axisCent, axisSp});

    // QA Before
    histos.add("QAbefore/Proton/hTPCNsigma", "n#sigma^{TPC} Protons", kTH2F, {axisPtQA, axisTPCNsigma});
    histos.add("QAbefore/Proton/hTOFNsigma", "n#sigma^{TOF} Protons", kTH2F, {axisPtQA, axisTOFNsigma});
    histos.add("QAbefore/Proton/hTpcTofNsigma", "n#sigma^{TPC} vs n#sigma^{TOF} Protons", kTH2F, {axisTPCNsigma, axisTOFNsigma});
    histos.add("QAbefore/Kaon/hTPCNsigma", "n#sigma^{TPC} Kaons", kTH2F, {axisPtQA, axisTPCNsigma});
    histos.add("QAbefore/Kaon/hTOFNsigma", "n#sigma^{TOF} Kaons", kTH2F, {axisPtQA, axisTOFNsigma});
    histos.add("QAbefore/Kaon/hTpcTofNsigma", "n#sigma^{TPC} vs n#sigma^{TOF} Kaons", kTH2F, {axisTPCNsigma, axisTOFNsigma});

    // QA After
    histos.add("QAafter/Proton/hPt", "p_{T}-spectra Protons", kTH1F, {axisPt});
    histos.add("QAafter/Proton/hDcaZ", "dca_{z} Protons", kTH2F, {axisPtQA, axisDCAz});
    histos.add("QAafter/Proton/hDcaXY", "dca_{xy} Protons", kTH2F, {axisPtQA, axisDCAxy});
    histos.add("QAafter/Proton/hTPCNsigma", "n#sigma^{TPC} only Protons", kTH2F, {axisPtQA, axisTPCNsigma});
    histos.add("QAafter/Proton/hTPCNsigmaTOF", "n#sigma^{TPC} Protons", kTH2F, {axisPtQA, axisTPCNsigma});
    histos.add("QAafter/Proton/hTOFNsigma", "n#sigma^{TOF} Protons", kTH2F, {axisPtQA, axisTOFNsigma});
    histos.add("QAafter/Proton/hTpcTofNsigma", "n#sigma^{TPC} vs n#sigma^{TOF} Protons", kTH2F, {axisTPCNsigma, axisTOFNsigma});
    histos.add("QAafter/Kaon/hPt", "p_{T}-spectra Kaons", kTH1F, {axisPt});
    histos.add("QAafter/Kaon/hDcaZ", "dca_{z} Kaons", kTH2F, {axisPtQA, axisDCAz});
    histos.add("QAafter/Kaon/hDcaXY", "dca_{xy} Kaons", kTH2F, {axisPtQA, axisDCAxy});
    histos.add("QAafter/Kaon/hTPCNsigma", "n#sigma^{TPC} only Kaons", kTH2F, {axisPtQA, axisTPCNsigma});
    histos.add("QAafter/Kaon/hTPCNsigmaTOF", "n#sigma^{TPC} Kaons", kTH2F, {axisPtQA, axisTPCNsigma});
    histos.add("QAafter/Kaon/hTOFNsigma", "n#sigma^{TOF} Kaons", kTH2F, {axisPtQA, axisTOFNsigma});
    histos.add("QAafter/Kaon/hTpcTofNsigma", "n#sigma^{TPC} vs n#sigma^{TOF} Kaons", kTH2F, {axisTPCNsigma, axisTOFNsigma});

    // Analysis
    // Lambda Invariant Mass
    histos.add("Analysis/hInvMass", "#Lambda(1520) M_{inv}", kTH1D, {axisInvM});
    histos.add("Analysis/hInvMassLS1", "Like Signs M_{inv} p K^{+}", kTH1D, {axisInvM});
    histos.add("Analysis/hInvMassLS2", "Like Signs M_{inv} #bar{p} K^{-}", kTH1D, {axisInvM});
    histos.add("Analysis/hInvMassR", "Rotated Spectra", kTH1D, {axisInvM});
    histos.add("Analysis/hInvMassMix", "Mixed Events M_{inv}", kTH1D, {axisInvM});
    histos.add("Analysis/h4InvMass", "THn #Lambda(1520)", kTHnSparseD, {axisInvM, axisPt, axisSp, axisCent});
    histos.add("Analysis/h4InvMassLS1", "THn Like Signs p K^{+}", kTHnSparseD, {axisInvM, axisPt, axisSp, axisCent});
    histos.add("Analysis/h4InvMassLS2", "THn Like Signs #bar{p} K^{-}", kTHnSparseD, {axisInvM, axisPt, axisSp, axisCent});
    histos.add("Analysis/h4InvMassR", "THn Rotated", kTHnSparseD, {axisInvM, axisPt, axisSp, axisCent});
    histos.add("Analysis/h4InvMassMix", "THn Mixed Events", kTHnSparseD, {axisInvM, axisPt, axisSp, axisCent});

    // MC
    if (doprocessMC) {
      histos.add("Event/hSphRec", "Reconstructed S_{0}", kTH1F, {axisSp});
      histos.add("Event/hSpCentRec", "Reconstructed S_{0} vs FT0M(%)", kTH2F, {axisCent, axisSp});
      histos.add("QAMCTrue/DcaZ_pr", "dca_{z}^{MC} Protons", kTH2F, {axisPtQA, axisDCAz});
      histos.add("QAMCTrue/DcaZ_ka", "dca_{z}^{MC} Kaons", kTH2F, {axisPtQA, axisDCAz});
      histos.add("QAMCTrue/DcaXY_pr", "dca_{xy}^{MC} Protons", kTH2F, {axisPtQA, axisDCAxy});
      histos.add("QAMCTrue/DcaXY_ka", "dca_{xy}^{MC} Kaons", kTH2F, {axisPtQA, axisDCAxy});
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
  bool selectionPIDProton(const T& candidate)
  {
    bool tpcPIDPassed{false}, tofPIDPassed{false};
    auto tpcPIDpt = static_cast<std::vector<float>>(protonTPCPIDpt);
    auto tpcPIDcut = static_cast<std::vector<float>>(protonTPCPIDcut);
    int nitr = static_cast<int>(tpcPIDpt.size());

    if (candidate.hasTOF()) {
      if (std::abs(candidate.tofNSigmaPr()) < cMaxTOFnSigmaProton) {
        tofPIDPassed = true;
      }
      // square cut
      if ((nsigmaCutCombinedProton < 0) && (std::abs(candidate.tpcNSigmaPr()) < cMaxTPCnSigmaProton)) {
        tpcPIDPassed = true;
      }
      // circular cut
      if ((nsigmaCutCombinedProton > 0) && (candidate.tpcNSigmaPr() * candidate.tpcNSigmaPr() + candidate.tofNSigmaPr() * candidate.tofNSigmaPr() < nsigmaCutCombinedProton * nsigmaCutCombinedProton)) {
        tofPIDPassed = true;
        tpcPIDPassed = true;
      }
    } else {
      tofPIDPassed = true;
      for (int i = 0; i < nitr - 1; ++i) {
        if (candidate.pt() >= tpcPIDpt[i] && candidate.pt() < tpcPIDpt[i + 1] && std::abs(candidate.tpcNSigmaPr()) < tpcPIDcut[i]) {
          tpcPIDPassed = true;
        }
        if (!cUseTpcAndTof && candidate.pt() >= tpcPIDpt[nitr - 1] && std::abs(candidate.tpcNSigmaPr()) < tpcPIDcut[nitr - 2]) {
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
  bool selectionPIDKaon(const T& candidate)
  {
    bool tpcPIDPassed{false}, tofPIDPassed{false};
    auto tpcPIDpt = static_cast<std::vector<float>>(kaonTPCPIDpt);
    auto tpcPIDcut = static_cast<std::vector<float>>(kaonTPCPIDcut);
    int nitr = static_cast<int>(tpcPIDpt.size());

    if (candidate.hasTOF()) {
      if (std::abs(candidate.tofNSigmaKa()) < cMaxTOFnSigmaKaon) {
        tofPIDPassed = true;
      }
      // square cut
      if ((nsigmaCutCombinedKaon < 0) && (std::abs(candidate.tpcNSigmaKa()) < cMaxTPCnSigmaKaon)) {
        tpcPIDPassed = true;
      }
      // circular
      if ((nsigmaCutCombinedKaon > 0) && (candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa() + candidate.tofNSigmaKa() * candidate.tofNSigmaKa() < nsigmaCutCombinedKaon * nsigmaCutCombinedKaon)) {
        tofPIDPassed = true;
        tpcPIDPassed = true;
      }
    } else {
      tofPIDPassed = true;
      for (int i = 0; i < nitr - 1; ++i) {
        if (candidate.pt() >= tpcPIDpt[i] && candidate.pt() < tpcPIDpt[i + 1] && std::abs(candidate.tpcNSigmaKa()) < tpcPIDcut[i]) {
          tpcPIDPassed = true;
        }
        if (!cUseTpcAndTof && candidate.pt() >= tpcPIDpt[nitr - 1] && std::abs(candidate.tpcNSigmaKa()) < tpcPIDcut[nitr - 2]) {
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
  void fillDataHistos(trackType const& trk1, trackType const& trk2, float const& sph, float const& mult)
  {
    TLorentzVector p1, p2, p;
    TRandom* rn = new TRandom();

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

      // Fill QA before track selection.
      if (!mix) {
        histos.fill(HIST("QAbefore/Proton/hTPCNsigma"), trkPr.pt(), trkPr.tpcNSigmaPr());
        if (trkPr.hasTOF()) {
          histos.fill(HIST("QAbefore/Proton/hTOFNsigma"), trkPr.pt(), trkPr.tofNSigmaPr());
          histos.fill(HIST("QAbefore/Proton/hTpcTofNsigma"), trkPr.tpcNSigmaPr(), trkPr.tofNSigmaPr());
        }
        histos.fill(HIST("QAbefore/Kaon/hTPCNsigma"), trkKa.pt(), trkKa.tpcNSigmaKa());
        if (trkKa.hasTOF()) {
          histos.fill(HIST("QAbefore/Kaon/hTOFNsigma"), trkKa.pt(), trkKa.tofNSigmaKa());
          histos.fill(HIST("QAbefore/Kaon/hTpcTofNsigma"), trkKa.tpcNSigmaKa(), trkKa.tofNSigmaKa());
        }
      }

      // Apply PID Selection
      if (cUseOnlyTOFTrackPr && !trkPr.hasTOF())
        continue;
      if (cUseOnlyTOFTrackKa && !trkKa.hasTOF())
        continue;
      if (!selectionPIDProton(trkPr) || !selectionPIDKaon(trkKa))
        continue;

      // Fill QA after track selection.
      if (!mix) {
        histos.fill(HIST("QAafter/Proton/hPt"), trkPr.pt());
        histos.fill(HIST("QAafter/Proton/hDcaZ"), trkPr.pt(), trkPr.dcaZ());
        histos.fill(HIST("QAafter/Proton/hDcaXY"), trkPr.pt(), trkPr.dcaXY());
        if (trkPr.hasTOF()) {
          histos.fill(HIST("QAafter/Proton/hTPCNsigmaTOF"), trkPr.pt(), trkPr.tpcNSigmaPr());
          histos.fill(HIST("QAafter/Proton/hTOFNsigma"), trkPr.pt(), trkPr.tofNSigmaPr());
          histos.fill(HIST("QAafter/Proton/hTpcTofNsigma"), trkPr.tpcNSigmaPr(), trkPr.tofNSigmaPr());
        } else {
          histos.fill(HIST("QAafter/Proton/hTPCNsigma"), trkPr.pt(), trkPr.tpcNSigmaPr());
        }
        histos.fill(HIST("QAafter/Kaon/hPt"), trkKa.pt());
        histos.fill(HIST("QAafter/Kaon/hDcaZ"), trkKa.pt(), trkKa.dcaZ());
        histos.fill(HIST("QAafter/Kaon/hDcaXY"), trkKa.pt(), trkKa.dcaXY());
        if (trkKa.hasTOF()) {
          histos.fill(HIST("QAafter/Kaon/hTPCNsigmaTOF"), trkKa.pt(), trkKa.tpcNSigmaKa());
          histos.fill(HIST("QAafter/Kaon/hTOFNsigma"), trkKa.pt(), trkKa.tofNSigmaKa());
          histos.fill(HIST("QAafter/Kaon/hTpcTofNsigma"), trkKa.tpcNSigmaKa(), trkKa.tofNSigmaKa());
        } else {
          histos.fill(HIST("QAafter/Kaon/hTPCNsigma"), trkKa.pt(), trkKa.tpcNSigmaKa());
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
        histos.fill(HIST("QAMCTrue/DcaXY_pr"), trkPr.pt(), trkPr.dcaXY());
        histos.fill(HIST("QAMCTrue/DcaXY_ka"), trkKa.pt(), trkKa.dcaXY());
        histos.fill(HIST("QAMCTrue/DcaZ_pr"), trkPr.pt(), trkPr.dcaZ());
        histos.fill(HIST("QAMCTrue/DcaZ_ka"), trkKa.pt(), trkKa.dcaZ());

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
