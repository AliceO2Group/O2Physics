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
  Configurable<int> nBinsPt{"nBinsPt", 100, "N bins in pT histogram"};
  Configurable<int> nBinsInvM{"nBinsInvM", 500, "N bins in InvMass histograms"};
  Configurable<int> nTracksSph{"nTracksSph", 3, "min #tracks for spherocity distribution"};

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
  Configurable<float> tpcProtonMaxPt{"tpcProtonMaxPt", 1.2, "max pT for tpc protons"};
  Configurable<float> tpcNSigmaProton{"tpcNSigmaProton", 3, "nsigma tpc for Proton when Tof signal is present"};
  Configurable<std::vector<float>> protonTPCPIDpt{"protonTPCPIDpt", {0, 0.5, 0.7, 1.2}, "pT dependent TPC cuts protons"};
  Configurable<std::vector<float>> protonTPCPIDcut{"protonTPCPIDcut", {4, 3, 2}, "TPC cuts protons"};
  Configurable<std::vector<float>> protonTOFPIDpt{"protonTOFPIDpt", {36.}, "pT dependent TOF cuts protons"};
  Configurable<std::vector<float>> protonTOFPIDcut{"protonTOFPIDCut", {3}, "TOF cuts protons"};

  // TPC TOF Kaons
  Configurable<float> tpcKaonMaxPt{"tpcKaonMaxPt", 0.6, "max pT for tpc kaons"};
  Configurable<float> tpcNSigmaKaon{"tpcNSigmaKaon", 3, "nsigma tpc for Kaon when Tof signal is present"};
  Configurable<std::vector<float>> kaonTPCPIDpt{"kaonTPCPIDpt", {0, 0.25, 0.4, 0.6}, "pT dependent TPC cuts kaons"};
  Configurable<std::vector<float>> kaonTPCPIDcut{"kaonTPCPIDcut", {4, 3, 2}, "TPC cuts kaons"};
  Configurable<std::vector<float>> kaonTOFPIDpt{"kaonTOFPIDpt", {36.}, "pT dependent TOF cuts kaons"};
  Configurable<std::vector<float>> kaonTOFPIDcut{"kaonTOFPIDcut", {3}, "TOF cuts kaons"};

  // Event Mixing.
  Configurable<int> nMix{"nMix", 5, "Number of Events to be mixed"};
  ConfigurableAxis cfgVtxBins{"cfgVtxBins", {VARIABLE_WIDTH, -10.0f, -9.f, -8.f, -7.f, -6.f, -5.f, -4.f, -3.f, -2.f, -1.f, 0.f, 1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis cfgMultBins{"cfgMultBins", {VARIABLE_WIDTH, 0.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.0f, 200.0f}, "Mixing bins - multiplicity"};

  // Histogram Registry.
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {

    // Define Axis.
    const AxisSpec axisSp(100, 0., 1., "S_{0}");
    const AxisSpec axisCent(105, 0, 105, "V0M (%)");
    const AxisSpec axisPtQA(200, 0., 2., "p_{T} (GeV/c)");
    const AxisSpec axisPt(nBinsPt, 0., 10., "p_{T} (GeV/c)");
    const AxisSpec axisEta(200, -1, 1, "#eta");
    const AxisSpec axisDCAz(500, -0.5, 0.5, {"DCA_{z} (#cm)"});
    const AxisSpec axisDCAxy(240, -0.12, 0.12, {"DCA_{xy} (#cm)"});
    const AxisSpec axisTPCNCls(200, 0, 200, {"TPCNCls"});
    const AxisSpec axisTPCNsigma(28, -7, 7, {"TPC N^{Sigma}"});
    const AxisSpec axisTOFNsigma(28, -7, 7, {"TOF N^{Sigma}"});
    const AxisSpec axisInvM(nBinsInvM, 1.4, 2.4, {"M_{inv} (GeV/c^{2})"});

    // Create Histograms.
    histos.add("Event/hCent", "V0M (%)", kTH1F, {axisCent});
    histos.add("Event/hSph", "Event Spherocity", kTH1F, {axisSp});
    histos.add("Event/hSpCent", "Spherocity vs V0M(%)", kTH2F, {axisCent, axisSp});
    histos.add("QAbefore/Proton/hTPCNsigma", "n#sigma^{TPC} Protons", kTH2F, {axisPtQA, axisTPCNsigma});
    histos.add("QAbefore/Proton/hTOFNsigma", "n#sigma^{TOF} Protons", kTH2F, {axisPtQA, axisTOFNsigma});
    histos.add("QAbefore/Proton/hTpcTofNsigma", "n#sigma^{TPC} vs n#sigma^{TOF} Protons", kTH2F, {axisTPCNsigma, axisTOFNsigma});
    histos.add("QAbefore/Kaon/hTPCNsigma", "n#sigma^{TPC} Kaons", kTH2F, {axisPtQA, axisTPCNsigma});
    histos.add("QAbefore/Kaon/hTOFNsigma", "n#sigma^{TOF} Kaons", kTH2F, {axisPtQA, axisTOFNsigma});
    histos.add("QAbefore/Kaon/hTpcTofNsigma", "n#sigma^{TPC} vs n#sigma^{TOF} Kaons", kTH2F, {axisTPCNsigma, axisTOFNsigma});
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
    histos.add("Analysis/hPtProton", "Protons p_{T}", kTH1F, {axisPt});
    histos.add("Analysis/hEtaProton", "Protons #eta", kTH1F, {axisEta});
    histos.add("Analysis/hPtKaon", "Kaons p_{T}", kTH1F, {axisPt});
    histos.add("Analysis/hEtaKaon", "Kaons #eta", kTH1F, {axisEta});
    histos.add("Analysis/hInvMass", "#Lambda(1520) M_{inv}", kTH1F, {axisInvM});
    histos.add("Analysis/hInvMassLS", "Like Signs M_{inv}", kTH1F, {axisInvM});
    histos.add("Analysis/hInvMassMix", "Mixed Events M_{inv}", kTH1F, {axisInvM});
    histos.add("Analysis/h4InvMass", "THn #Lambda(1520)", kTHnSparseF, {axisInvM, axisPt, axisSp, axisCent});
    histos.add("Analysis/h4InvMassLS", "THn Like Signs", kTHnSparseF, {axisInvM, axisPt, axisSp, axisCent});
    histos.add("Analysis/h4InvMassMix", "THn Mixed Events", kTHnSparseF, {axisInvM, axisPt, axisSp, axisCent});
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
        if ((trkPr.tofPIDselectionFlag() & aod::resodaughter::kHasTOF) == aod::resodaughter::kHasTOF) {
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
        if ((trkKa.tofPIDselectionFlag() & aod::resodaughter::kHasTOF) == aod::resodaughter::kHasTOF) {
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
      if (!mix) {
        if (std::abs(trkPr.tpcNSigmaPr()) < 6) {
          histos.fill(HIST("QAbefore/Proton/hTPCNsigma"), trkPr.pt(), trkPr.tpcNSigmaPr());
          if (std::abs(trkPr.tofNSigmaPr()) < 6) {
            histos.fill(HIST("QAbefore/Proton/hTOFNsigma"), trkPr.pt(), trkPr.tofNSigmaPr());
            histos.fill(HIST("QAbefore/Proton/hTpcTofNsigma"), trkPr.tpcNSigmaPr(), trkPr.tofNSigmaPr());
          }
        }
        if (std::abs(trkKa.tpcNSigmaKa()) < 6) {
          histos.fill(HIST("QAbefore/Kaon/hTPCNsigma"), trkKa.pt(), trkKa.tpcNSigmaKa());
          if (std::abs(trkKa.tofNSigmaKa()) < 6) {
            histos.fill(HIST("QAbefore/Kaon/hTOFNsigma"), trkKa.pt(), trkKa.tofNSigmaKa());
            histos.fill(HIST("QAbefore/Kaon/hTpcTofNsigma"), trkKa.tpcNSigmaKa(), trkKa.tofNSigmaKa());
          }
        }
      }

      // Apply PID Selection.
      if (!isTrk1Proton || !isTrk2Kaon)
        continue;

      // Fill QA after track selection.
      if (!mix) {
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

      // Protons and Kaons pT spectra and Pseudorapidity distribution.
      histos.fill(HIST("Analysis/hPtProton"), trkPr.pt());
      histos.fill(HIST("Analysis/hEtaProton"), trkPr.eta());
      histos.fill(HIST("Analysis/hPtKaon"), trkKa.pt());
      histos.fill(HIST("Analysis/hEtaKaon"), trkKa.eta());

      // Invariant mass reconstruction.
      p1.SetXYZM(trkPr.px(), trkPr.py(), trkPr.pz(), massProton);
      p2.SetXYZM(trkKa.px(), trkKa.py(), trkKa.pz(), massKaon);
      p = p1 + p2;

      if (std::abs(p.Rapidity()) > 0.5)
        continue;

      // Fill Invariant Mass Histograms.
      if (trkPr.sign() * trkKa.sign() < 0) {
        if (!mix) {
          histos.fill(HIST("Analysis/hInvMass"), p.M());
          histos.fill(HIST("Analysis/h4InvMass"), p.M(), p.Pt(), sph, mult);
        } else {
          histos.fill(HIST("Analysis/hInvMassMix"), p.M());
          histos.fill(HIST("Analysis/h4InvMassMix"), p.M(), p.Pt(), sph, mult);
        }
      }

      if (trkPr.sign() * trkKa.sign() > 0 && !mix) {
        histos.fill(HIST("Analysis/hInvMassLS"), p.M());
        histos.fill(HIST("Analysis/h4InvMassLS"), p.M(), p.Pt(), sph, mult);
      }
    }
  }

  using resoCols = aod::ResoCollisions;
  using resoTracks = aod::ResoTracks;

  void processData(resoCols::iterator const& collision, resoTracks const& tracks)
  {

    histos.fill(HIST("Event/hCent"), collision.multV0M());

    if (tracks.size() >= nTracksSph) {
      histos.fill(HIST("Event/hSph"), collision.spherocity());
      histos.fill(HIST("Event/hSpCent"), collision.multV0M(), collision.spherocity());
    }

    fillDataHistos<false, false>(tracks, tracks, collision.spherocity(), collision.multV0M());
  }

  PROCESS_SWITCH(lambdaAnalysis, processData, "Process for Same Event Data", true);

  // Processing Event Mixing
  SliceCache cache;
  using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::resocollision::MultV0M>;
  BinningType binningPositions{{cfgVtxBins, cfgMultBins}, true};

  void processMix(resoCols& collisions, resoTracks const& tracks)
  {

    LOGF(debug, "Event Mixing Started");
    auto tracksTuple = std::make_tuple(tracks);
    SameKindPair<resoCols, resoTracks, BinningType> pairs{binningPositions, nMix, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip
    for (auto& [c1, t1, c2, t2] : pairs) {
      fillDataHistos<true, false>(t1, t2, c1.spherocity(), c1.multV0M());
    }
  }

  PROCESS_SWITCH(lambdaAnalysis, processMix, "Process for Mixed Events", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<lambdaAnalysis>(cfgc)};
}
