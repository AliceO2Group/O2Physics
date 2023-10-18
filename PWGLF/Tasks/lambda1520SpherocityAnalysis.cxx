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

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// PDG p -> 2212, K -> 321
const float massProton = TDatabasePDG::Instance()->GetParticle(2212)->Mass();
const float massKaon = TDatabasePDG::Instance()->GetParticle(321)->Mass();

struct lambdaAnalysis {

  // Configurables.
  Configurable<int> nBinsPt{"nBinsPt", 200, "N bins in pT histogram"};
  Configurable<int> nBinsInvM{"nBinsInvM", 400, "N bins in InvMass histograms"};
  Configurable<bool> doRotate{"doRotate", true, "rotated inv mass spectra"};

  // Tracks
  Configurable<float> cfgPtMin{"ptMin", 0.15, "Minimum Track pT"};
  Configurable<float> cfgPtMax{"ptMax", 36., "Maximum Track pT"};
  Configurable<float> cfgEtaCut{"etaCut", 0.8, "Pseudorapidity cut"};
  Configurable<float> cfgDcaz{"dcazMin", 1., "Minimum DCAz"};
  Configurable<float> cfgDcaxy{"dcaxyMin", 0.1, "Minimum DCAxy"};
  Configurable<float> cfgPIDprecut{"cfgPIDprecut", 6, "Preselection PID TPC TOF cut"};
  Configurable<bool> cfgGlobalTrackWoDCA{"cfgGlobalTrackWoDCA", true, "Global Track Selection"};
  Configurable<bool> cfgPVContributor{"cfgPVContributor", true, "PV Contributor Track Selection"};

  // TPC TOF Protons
  Configurable<float> tpcProtonMaxPt{"tpcProtonMaxPt", 1.0, "max pT for tpc protons"};
  Configurable<float> tpcNSigmaProton{"tpcNSigmaProton", 3, "nsigma tpc for Proton when Tof signal is present"};
  Configurable<std::vector<float>> protonTPCPIDpt{"protonTPCPIDpt", {0, 0.5, 0.7, 0.8, 1.0}, "pT dependent TPC cuts protons"};
  Configurable<std::vector<float>> protonTPCPIDcut{"protonTPCPIDcut", {5., 3.5, 2.5, 1.5}, "TPC cuts protons"};
  Configurable<std::vector<float>> protonTOFPIDpt{"protonTOFPIDpt", {36.}, "pT dependent TOF cuts protons"};
  Configurable<std::vector<float>> protonTOFPIDcut{"protonTOFPIDCut", {3}, "TOF cuts protons"};

  // TPC TOF Kaons
  Configurable<float> tpcKaonMaxPt{"tpcKaonMaxPt", 0.5, "max pT for tpc kaons"};
  Configurable<float> tpcNSigmaKaon{"tpcNSigmaKaon", 3, "nsigma tpc for Kaon when Tof signal is present"};
  Configurable<std::vector<float>> kaonTPCPIDpt{"kaonTPCPIDpt", {0, 0.2, 0.3, 0.4, 0.50}, "pT dependent TPC cuts kaons"};
  Configurable<std::vector<float>> kaonTPCPIDcut{"kaonTPCPIDcut", {5, 3.5, 2., 1.5}, "TPC cuts kaons"};
  Configurable<std::vector<float>> kaonTOFPIDpt{"kaonTOFPIDpt", {36.}, "pT dependent TOF cuts kaons"};
  Configurable<std::vector<float>> kaonTOFPIDcut{"kaonTOFPIDcut", {3}, "TOF cuts kaons"};

  // Event Mixing.
  Configurable<bool> doSphMix{"doSphMix", true, "Include Sph Bins to be mixed"};
  Configurable<int> nMix{"nMix", 10, "Number of Events to be mixed"};
  ConfigurableAxis cfgVtxBins{"cfgVtxBins", {VARIABLE_WIDTH, -10.0f, -9.f, -8.f, -7.f, -6.f, -5.f, -4.f, -3.f, -2.f, -1.f, 0.f, 1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis cfgMultBins{"cfgMultBins", {VARIABLE_WIDTH, 0.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.0f, 200.0f}, "Mixing bins - multiplicity"};
  ConfigurableAxis cfgSphBins{"cfgSphBins", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f}, "Mixing bins - spherocity"};

  // Histogram Registry.
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {

    // Define Axis.
    const AxisSpec axisSp(1000, 0., 1., "S_{0}");
    const AxisSpec axisCent(105, 0, 105, "FT0M (%)");
    const AxisSpec axisPtQA(200, 0., 2., "p_{T} (GeV/c)");
    const AxisSpec axisPt(nBinsPt, 0., 10., "p_{T} (GeV/c)");
    const AxisSpec axisEta(200, -1, 1, "#eta");
    const AxisSpec axisDCAz(500, -0.5, 0.5, {"DCA_{z} (cm)"});
    const AxisSpec axisDCAxy(240, -0.12, 0.12, {"DCA_{xy} (cm)"});
    const AxisSpec axisTPCNCls(200, 0, 200, {"TPCNCls"});
    const AxisSpec axisTPCNsigma(140, -10, 10, {"n#sigma^{TPC}"});
    const AxisSpec axisTOFNsigma(140, -10, 10, {"n#sigma^{TOF}"});
    const AxisSpec axisInvM(nBinsInvM, 1.4, 3.4, {"M_{inv} (GeV/c^{2})"});

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
    histos.add("QAafter/Proton/hDcaZ", "dca_{z} Protons", kTH2F, {axisPtQA, axisDCAz});
    histos.add("QAafter/Proton/hDcaXY", "dca_{xy} Protons", kTH2F, {axisPtQA, axisDCAxy});
    histos.add("QAafter/Proton/hTPCNsigma", "n#sigma^{TPC} Protons", kTH2F, {axisPtQA, axisTPCNsigma});
    histos.add("QAafter/Proton/hTOFNsigma", "n#sigma^{TOF} Protons", kTH2F, {axisPtQA, axisTOFNsigma});
    histos.add("QAafter/Proton/hTpcTofNsigma", "n#sigma^{TPC} vs n#sigma^{TOF} Protons", kTH2F, {axisTPCNsigma, axisTOFNsigma});
    histos.add("QAafter/Kaon/hDcaZ", "dca_{z} Kaons", kTH2F, {axisPtQA, axisDCAz});
    histos.add("QAafter/Kaon/hDcaXY", "dca_{xy} Kaons", kTH2F, {axisPtQA, axisDCAxy});
    histos.add("QAafter/Kaon/hTPCNsigma", "n#sigma^{TPC} Kaons", kTH2F, {axisPtQA, axisTPCNsigma});
    histos.add("QAafter/Kaon/hTOFNsigma", "n#sigma^{TOF} Kaons", kTH2F, {axisPtQA, axisTOFNsigma});
    histos.add("QAafter/Kaon/hTpcTofNsigma", "n#sigma^{TPC} vs n#sigma^{TOF} Kaons", kTH2F, {axisTPCNsigma, axisTOFNsigma});

    // Analysis
    // Proton Kaon
    histos.add("Analysis/hPtProton", "Protons p_{T}", kTH1F, {axisPt});
    histos.add("Analysis/hEtaProton", "Protons #eta", kTH1F, {axisEta});
    histos.add("Analysis/hPtKaon", "Kaons p_{T}", kTH1F, {axisPt});
    histos.add("Analysis/hEtaKaon", "Kaons #eta", kTH1F, {axisEta});

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

    if (track.pt() < cfgPtMin || track.pt() > cfgPtMax)
      return false;

    if (std::abs(track.dcaZ()) > cfgDcaz)
      return false;

    if (std::abs(track.dcaXY()) > cfgDcaxy)
      return false;

    if (cfgGlobalTrackWoDCA && !track.isGlobalTrackWoDCA())
      return false;

    if (cfgPVContributor && !track.isPVContributor())
      return false;

    return true;
  }

  template <bool mix, bool mc, typename trackType>
  void fillDataHistos(trackType const& trk1, trackType const& trk2, float const& sph, float const& mult)
  {

    bool isTrk1Proton{true}, isTrk2Kaon{true}, trk1HasTOF{false}, trk2HasTOF{false}, selTrk1{true}, selTrk2{true};

    auto prTpcPIDpt = static_cast<std::vector<float>>(protonTPCPIDpt);
    auto prTpcPIDcut = static_cast<std::vector<float>>(protonTPCPIDcut);
    auto prTofPIDpt = static_cast<std::vector<float>>(protonTOFPIDpt);
    auto prTofPIDcut = static_cast<std::vector<float>>(protonTOFPIDcut);
    auto kaTpcPIDpt = static_cast<std::vector<float>>(kaonTPCPIDpt);
    auto kaTpcPIDcut = static_cast<std::vector<float>>(kaonTPCPIDcut);
    auto kaTofPIDpt = static_cast<std::vector<float>>(kaonTOFPIDpt);
    auto kaTofPIDcut = static_cast<std::vector<float>>(kaonTOFPIDcut);

    TLorentzVector p1, p2, p;
    TRandom* rn = new TRandom();

    for (auto const& [trkPr, trkKa] : soa::combinations(soa::CombinationsFullIndexPolicy(trk1, trk2))) {
      // Do not analyse same index tracks.
      if (trkPr.index() == trkKa.index())
        continue;

      selTrk1 = true;
      selTrk2 = true;
      isTrk1Proton = true;
      isTrk2Kaon = true;
      trk1HasTOF = false;
      trk2HasTOF = false;

      // pT, DCA, Global Tracks and PVcontrib selection.
      selTrk1 = selTracks(trkPr);
      selTrk2 = selTracks(trkKa);

      if (!selTrk1 || !selTrk2)
        continue;

      // Protons.
      if (trkPr.pt() < tpcProtonMaxPt) {
        // TPC only
        for (int i = 1; i < static_cast<int>(prTpcPIDpt.size()); ++i) {
          if (trkPr.pt() >= prTpcPIDpt[i - 1] && trkPr.pt() < prTpcPIDpt[i]) {
            if (std::abs(trkPr.tpcNSigmaPr()) >= prTpcPIDcut[i - 1]) {
              isTrk1Proton = false;
            }
          }
        }
      } else {
        // TPC + TOF
        if (trkPr.hasTOF()) {
          trk1HasTOF = true;
          for (int i = 0; i < static_cast<int>(prTofPIDpt.size()); ++i) {
            if (trkPr.pt() < prTofPIDpt[i]) {
              if (std::abs(trkPr.tofNSigmaPr()) >= prTofPIDcut[i] || std::abs(trkPr.tpcNSigmaPr()) >= tpcNSigmaProton) {
                isTrk1Proton = false;
                trk1HasTOF = false;
              }
            }
          }
        } else {
          isTrk1Proton = false;
        }
      }

      // Kaons
      if (trkKa.pt() < tpcKaonMaxPt) {
        // TPC only
        for (int i = 1; i < static_cast<int>(kaTpcPIDpt.size()); ++i) {
          if (trkKa.pt() >= kaTpcPIDpt[i - 1] && trkKa.pt() < kaTpcPIDpt[i]) {
            if (std::abs(trkKa.tpcNSigmaKa()) >= kaTpcPIDcut[i - 1]) {
              isTrk2Kaon = false;
            }
          }
        }
      } else {
        // TPC + TOF
        if (trkKa.hasTOF()) {
          trk2HasTOF = true;
          for (int i = 0; i < static_cast<int>(kaTofPIDpt.size()); ++i) {
            if (trkKa.pt() < kaTofPIDpt[i]) {
              if (std::abs(trkKa.tofNSigmaKa()) >= kaTofPIDcut[i] || std::abs(trkKa.tpcNSigmaKa()) >= tpcNSigmaKaon) {
                isTrk2Kaon = false;
                trk2HasTOF = false;
              }
            }
          }
        } else {
          isTrk2Kaon = false;
        }
      }

      // Fill QA before track selection.
      if (!mix & !mc) {
        if (std::abs(trkPr.tpcNSigmaPr()) < cfgPIDprecut) {
          histos.fill(HIST("QAbefore/Proton/hTPCNsigma"), trkPr.pt(), trkPr.tpcNSigmaPr());
          if (std::abs(trkPr.tofNSigmaPr()) < cfgPIDprecut) {
            histos.fill(HIST("QAbefore/Proton/hTOFNsigma"), trkPr.pt(), trkPr.tofNSigmaPr());
            histos.fill(HIST("QAbefore/Proton/hTpcTofNsigma"), trkPr.tpcNSigmaPr(), trkPr.tofNSigmaPr());
          }
        }
        if (std::abs(trkKa.tpcNSigmaKa()) < cfgPIDprecut) {
          histos.fill(HIST("QAbefore/Kaon/hTPCNsigma"), trkKa.pt(), trkKa.tpcNSigmaKa());
          if (std::abs(trkKa.tofNSigmaKa()) < cfgPIDprecut) {
            histos.fill(HIST("QAbefore/Kaon/hTOFNsigma"), trkKa.pt(), trkKa.tofNSigmaKa());
            histos.fill(HIST("QAbefore/Kaon/hTpcTofNsigma"), trkKa.tpcNSigmaKa(), trkKa.tofNSigmaKa());
          }
        }
      }

      // Apply PID Selection.
      if (!isTrk1Proton || !isTrk2Kaon)
        continue;

      // Fill QA after track selection.
      if (!mix && !mc) {
        histos.fill(HIST("QAafter/Proton/hDcaZ"), trkPr.pt(), trkPr.dcaZ());
        histos.fill(HIST("QAafter/Proton/hDcaXY"), trkPr.pt(), trkPr.dcaXY());
        histos.fill(HIST("QAafter/Proton/hTPCNsigma"), trkPr.pt(), trkPr.tpcNSigmaPr());
        if (trk1HasTOF) {
          histos.fill(HIST("QAafter/Proton/hTOFNsigma"), trkPr.pt(), trkPr.tofNSigmaPr());
          histos.fill(HIST("QAafter/Proton/hTpcTofNsigma"), trkPr.tpcNSigmaPr(), trkPr.tofNSigmaPr());
        }
        histos.fill(HIST("QAafter/Kaon/hDcaZ"), trkKa.pt(), trkKa.dcaZ());
        histos.fill(HIST("QAafter/Kaon/hDcaXY"), trkKa.pt(), trkKa.dcaXY());
        histos.fill(HIST("QAafter/Kaon/hTPCNsigma"), trkKa.pt(), trkKa.tpcNSigmaKa());
        if (trk2HasTOF) {
          histos.fill(HIST("QAafter/Kaon/hTOFNsigma"), trkKa.pt(), trkKa.tofNSigmaKa());
          histos.fill(HIST("QAafter/Kaon/hTpcTofNsigma"), trkKa.tpcNSigmaKa(), trkKa.tofNSigmaKa());
        }
      }

      // Invariant mass reconstruction.
      p1.SetXYZM(trkPr.px(), trkPr.py(), trkPr.pz(), massProton);
      p2.SetXYZM(trkKa.px(), trkKa.py(), trkKa.pz(), massKaon);
      p = p1 + p2;

      if (std::abs(p.Rapidity()) > 0.5)
        continue;

      // Fill Invariant Mass Histograms.
      if constexpr (!mix && !mc) {

        // Protons and Kaons pT spectra and Pseudorapidity distribution.
        histos.fill(HIST("Analysis/hPtProton"), trkPr.pt());
        histos.fill(HIST("Analysis/hEtaProton"), trkPr.eta());
        histos.fill(HIST("Analysis/hPtKaon"), trkKa.pt());
        histos.fill(HIST("Analysis/hEtaKaon"), trkKa.eta());

        if (trkPr.sign() * trkKa.sign() < 0) {
          histos.fill(HIST("Analysis/hInvMass"), p.M());
          histos.fill(HIST("Analysis/h4InvMass"), p.M(), p.Pt(), sph, mult);
          if (doRotate) {
            float theta = rn->Uniform(2.0, 2.5);
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
          histos.fill(HIST("Analysis/h4InvMassLambdaAnti"), p.M(), p.Pt(), sph, mult);
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

    histos.fill(HIST("Event/hCent"), collision.multV0M());
    histos.fill(HIST("Event/hSph"), collision.spherocity());
    histos.fill(HIST("Event/hSpCent"), collision.multV0M(), collision.spherocity());

    fillDataHistos<false, false>(tracks, tracks, collision.spherocity(), collision.multV0M());
  }

  PROCESS_SWITCH(lambdaAnalysis, processData, "Process for Same Event Data", true);

  void processMC(resoCols::iterator const& collision,
                 soa::Join<aod::ResoTracks, aod::ResoMCTracks> const& tracks, aod::McParticles const& mcParticles)
  {
    fillDataHistos<false, true>(tracks, tracks, collision.spherocity(), collision.multV0M());
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
  SliceCache cache;

  using BinningType1 = ColumnBinningPolicy<aod::collision::PosZ, aod::resocollision::MultV0M, aod::resocollision::Spherocity>;
  using BinningType2 = ColumnBinningPolicy<aod::collision::PosZ, aod::resocollision::MultV0M>;
  void processMix(resoCols& collisions, resoTracks const& tracks)
  {

    LOGF(debug, "Event Mixing Started");
    BinningType1 binningPositions1{{cfgVtxBins, cfgMultBins, cfgSphBins}, true};
    BinningType2 binningPositions2{{cfgVtxBins, cfgMultBins}, true};
    auto tracksTuple = std::make_tuple(tracks);
    if (doSphMix) {
      SameKindPair<resoCols, resoTracks, BinningType1> pairs{binningPositions1, nMix, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip
      for (auto& [c1, t1, c2, t2] : pairs) {
        fillDataHistos<true, false>(t1, t2, c1.spherocity(), c1.multV0M());
      }
    } else {
      SameKindPair<resoCols, resoTracks, BinningType2> pairs{binningPositions2, nMix, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip
      for (auto& [c1, t1, c2, t2] : pairs) {
        fillDataHistos<true, false>(t1, t2, c1.spherocity(), c1.multV0M());
      }
    }
  }

  PROCESS_SWITCH(lambdaAnalysis, processMix, "Process for Mixed Events", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<lambdaAnalysis>(cfgc)};
}
