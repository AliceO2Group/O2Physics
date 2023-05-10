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

/// \file taskXic.cxx
/// \brief Ξc± analysis task
/// \note Inspired from taskLc.cxx and SigmaC.cxx
///
/// \author Mattia Faggin <mattia.faggin@cern.ch>, University and INFN PADOVA
/// \author Anton Alkin <anton.alkin@cern.ch>, CERN
/// \author Jinjoo Seo <jin.joo.seo@cern.ch>, Inha University

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::hf_cand_3prong;
using namespace o2::analysis::hf_cuts_xic_to_p_k_pi;

#include "Framework/runDataProcessing.h"

/// Ξc± analysis task

struct HfTaskXic {

  Configurable<int> selectionFlagXic{"selectionFlagXic", 1, "Selection Flag for Xic"};
  Configurable<double> yCandMax{"yCandMax", -1., "max. cand. rapidity"};
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_xic_to_p_k_pi::vecBinsPt}, "pT bin limits"};

  Filter filterSelectCandidates = (aod::hf_sel_candidate_xic::isSelXicToPKPi >= selectionFlagXic || aod::hf_sel_candidate_xic::isSelXicToPiKP >= selectionFlagXic);

  Partition<soa::Join<aod::HfCand3Prong, aod::HfSelXicToPKPi, aod::HfCand3ProngMcRec>> selectedMCXicCandidates = (aod::hf_sel_candidate_xic::isSelXicToPKPi >= selectionFlagXic || aod::hf_sel_candidate_xic::isSelXicToPiKP >= selectionFlagXic);

  HistogramRegistry registry{
    "registry", // histo not in pt bins
    {
      {"Data/hPt", "3-prong candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},      // pt Xic
      {"Data/hEta", "3-prong candidates;candidate #it{eta};entries", {HistType::kTH1F, {{100, -5., 5.}}}},                    // eta Xic
      {"Data/hPhi", "3-prong candidates;candidate #varphi;entries", {HistType::kTH1F, {{72, 0., 2. * constants::math::PI}}}}, // phi Xic
      {"Data/hMass", "3-prong candidates; inv. mass (p K #pi) (GeV/#it{c}^{2})", {HistType::kTH1F, {{600, 2.18, 2.58}}}},     // mass Xic
      {"Data/hMultiplicity", "multiplicity;multiplicity;entries", {HistType::kTH1F, {{1000, 0., 1000.}}}},

      {"MC/reconstructed/signal/hMassRecCand", "3-prong candidates (matched, prompt); inv. mass (p K #pi) (GeV/#it{c}^{2})", {HistType::kTH1F, {{600, 2.18, 2.58}}}}, // mass Xic
      {"MC/reconstructed/signal/hPtRecCand", "3-prong candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},                    // pt Xic
      {"MC/reconstructed/signal/hPtRecSig", "3-prong candidates (matched);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}}},
      {"MC/reconstructed/background/hPtRecBg", "3-prong candidates (unmatched);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}}},
      {"MC/generated/signal/hPtGen", "MC particles (matched);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}}},
      {"MC/generated/signal/hPtGenSig", "3-prong candidates (matched);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}}} ///
    }};

  void init(o2::framework::InitContext&)
  {

    AxisSpec axisPidP = {100, 0.f, 10.0f, "#it{p} (GeV/#it{c})"};
    AxisSpec axisNSigmaPr = {100, -6.f, 6.f, "n#it{#sigma}_{p}"};
    AxisSpec axisNSigmaPi = {100, -6.f, 6.f, "n#it{#sigma}_{#pi}"};
    AxisSpec axisNSigmaKa = {100, -6.f, 6.f, "n#it{#sigma}_{K}"};

    auto vbins = (std::vector<double>)binsPt; // histo in pt bins
    registry.add("Data/hMassVsPt", "3-prong candidates;inv. mass (p K #pi) (GeV/#it{c}^{2});;entries", {HistType::kTH2F, {{500, 2., 3.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/hDecLengthVsPt", "3-prong candidates;decay length (cm);;entries", {HistType::kTH2F, {{800, 0., 4.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/hd0Prong0VsPt", "3-prong candidates;prong 0 DCAxy to prim. vertex (cm);;entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/hd0Prong1VsPt", "3-prong candidates;prong 1 DCAxy to prim. vertex (cm);;entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/hd0Prong2VsPt", "3-prong candidates;prong 1 DCAxy to prim. vertex (cm);;entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/hCtVsPt", "3-prong candidates;proper lifetime (#Xi_{c}) * #it{c} (cm);;entries", {HistType::kTH2F, {{120, -20., 100.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/hCPAVsPt", "3-prong candidates;cosine of pointing angle;;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/hEtaVsPt", "3-prong candidates;candidate #it{#eta};;entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/hSelectionStatusVsPt", "3-prong candidates;selection status;;entries", {HistType::kTH2F, {{5, -0.5, 4.5}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/hImpParErr0VsPt", "3-prong candidates;impact parameter error prong 0 (cm);;entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/hImpParErr1VsPt", "3-prong candidates;impact parameter error prong 1 (cm);;entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/hImpParErr2VsPt", "3-prong candidates;impact parameter error prong 2 (cm);;entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/hDecLenErrVsPt", "3-prong candidates;decay length error (cm);;entries", {HistType::kTH2F, {{100, 0., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/hPtProng0VsPt", "3-prong candidates;prong 0 #it{p}_{T} (GeV/#it{c});;entries", {HistType::kTH2F, {{100, 0., 10.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/hPtProng1VsPt", "3-prong candidates;prong 1 #it{p}_{T} (GeV/#it{c});;entries", {HistType::kTH2F, {{100, 0., 10.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/hPtProng2VsPt", "3-prong candidates;prong 2 #it{p}_{T} (GeV/#it{c});;entries", {HistType::kTH2F, {{100, 0., 10.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/hChi2PCAVsPt", "3-prong candidates;prong Chi2PCA to sec. vertex (cm);;entries", {HistType::kTH2F, {{100, 0, 0.5}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    // TPC nSigma histograms
    registry.add("Data/hPVsTPCNSigmaPr_Prong0", "3-prong candidates;#it{p} bachelor (GeV/#it{c}) ;n#it{#sigma}_{p} TPC", {HistType::kTH2F, {axisPidP, axisNSigmaPr}});
    registry.add("Data/hPVsTPCNSigmaPi_Prong0", "3-prong candidates;#it{p} bachelor (GeV/#it{c}) ;n#it{#sigma}_{#pi} TPC", {HistType::kTH2F, {axisPidP, axisNSigmaPi}});
    registry.add("Data/hPVsTPCNSigmaKa_Prong0", "3-prong candidates;#it{p} bachelor (GeV/#it{c}) ;n#it{#sigma}_{K} TPC", {HistType::kTH2F, {axisPidP, axisNSigmaKa}});

    registry.add("Data/hPVsTPCNSigmaPr_Prong1", "3-prong candidates;#it{p} bachelor (GeV/#it{c}) ;n#it{#sigma}_{p} TPC", {HistType::kTH2F, {axisPidP, axisNSigmaPr}});
    registry.add("Data/hPVsTPCNSigmaPi_Prong1", "3-prong candidates;#it{p} bachelor (GeV/#it{c}) ;n#it{#sigma}_{#pi} TPC", {HistType::kTH2F, {axisPidP, axisNSigmaPi}});
    registry.add("Data/hPVsTPCNSigmaKa_Prong1", "3-prong candidates;#it{p} bachelor (GeV/#it{c}) ;n#it{#sigma}_{K} TPC", {HistType::kTH2F, {axisPidP, axisNSigmaKa}});

    registry.add("Data/hPVsTPCNSigmaPr_Prong2", "3-prong candidates;#it{p} bachelor (GeV/#it{c}) ;n#it{#sigma}_{p} TPC", {HistType::kTH2F, {axisPidP, axisNSigmaPr}});
    registry.add("Data/hPVsTPCNSigmaPi_Prong2", "3-prong candidates;#it{p} bachelor (GeV/#it{c}) ;n#it{#sigma}_{#pi} TPC", {HistType::kTH2F, {axisPidP, axisNSigmaPi}});
    registry.add("Data/hPVsTPCNSigmaKa_Prong2", "3-prong candidates;#it{p} bachelor (GeV/#it{c}) ;n#it{#sigma}_{K} TPC", {HistType::kTH2F, {axisPidP, axisNSigmaKa}});

    // TOF nSigma histograms
    registry.add("Data/hPVsTOFNSigmaPr_Prong0", "3-prong candidates;#it{p} bachelor (GeV/#it{c}) ;n#it{#sigma}_{p} TOF", {HistType::kTH2F, {axisPidP, axisNSigmaPr}});
    registry.add("Data/hPVsTOFNSigmaPi_Prong0", "3-prong candidates;#it{p} bachelor (GeV/#it{c}) ;n#it{#sigma}_{#pi} TOF", {HistType::kTH2F, {axisPidP, axisNSigmaPi}});
    registry.add("Data/hPVsTOFNSigmaKa_Prong0", "3-prong candidates;#it{p} bachelor (GeV/#it{c}) ;n#it{#sigma}_{K} TOF", {HistType::kTH2F, {axisPidP, axisNSigmaKa}});

    registry.add("Data/hPVsTOFNSigmaPr_Prong1", "3-prong candidates;#it{p} bachelor (GeV/#it{c}) ;n#it{#sigma}_{p} TOF", {HistType::kTH2F, {axisPidP, axisNSigmaPr}});
    registry.add("Data/hPVsTOFNSigmaPi_Prong1", "3-prong candidates;#it{p} bachelor (GeV/#it{c}) ;n#it{#sigma}_{#pi} TOF", {HistType::kTH2F, {axisPidP, axisNSigmaPi}});
    registry.add("Data/hPVsTOFNSigmaKa_Prong1", "3-prong candidates;#it{p} bachelor (GeV/#it{c}) ;n#it{#sigma}_{K} TOF", {HistType::kTH2F, {axisPidP, axisNSigmaKa}});

    registry.add("Data/hPVsTOFNSigmaPr_Prong2", "3-prong candidates;#it{p} bachelor (GeV/#it{c}) ;n#it{#sigma}_{p} TOF", {HistType::kTH2F, {axisPidP, axisNSigmaPr}});
    registry.add("Data/hPVsTOFNSigmaPi_Prong2", "3-prong candidates;#it{p} bachelor (GeV/#it{c}) ;n#it{#sigma}_{#pi} TOF", {HistType::kTH2F, {axisPidP, axisNSigmaPi}});
    registry.add("Data/hPVsTOFNSigmaKa_Prong2", "3-prong candidates;#it{p} bachelor (GeV/#it{c}) ;n#it{#sigma}_{K} TOF", {HistType::kTH2F, {axisPidP, axisNSigmaKa}});

    registry.add("MC/reconstructed/signal/hMassSig", "Invariant mass (matched);m (p K #pi) (GeV/#it{c}^{2});;entries", {HistType::kTH2F, {{500, 1.6, 3.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/background/hMassBg", "Invariant mass (unmatched);m (p K #pi) (GeV/#it{c}^{2});;entries", {HistType::kTH2F, {{500, 1.6, 3.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/generated/signal/hEtaGen", "3-prong candidates;candidate #it{#eta};;entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    registry.add("MC/reconstructed/signal/hDecLengthRecSig", "3-prong candidates;decay length (cm);;entries", {HistType::kTH2F, {{200, 0., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hd0Prong0RecSig", "3-prong candidates;prong 0 DCAxy to prim. vertex (cm);;entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hd0Prong1RecSig", "3-prong candidates;prong 1 DCAxy to prim. vertex (cm);;entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hd0Prong2RecSig", "3-prong candidates;prong 1 DCAxy to prim. vertex (cm);;entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hCtRecSig", "3-prong candidates;proper lifetime (#Xi_{c}) * #it{c} (cm);;entries", {HistType::kTH2F, {{120, -20., 100.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hCPARecSig", "3-prong candidates;cosine of pointing angle;;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hEtaRecSig", "3-prong candidates;candidate #it{#eta};;entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hPtProng0RecSig", "3-prong candidates;prong 0 #it{p}_{T} (GeV/#it{c});;entries", {HistType::kTH2F, {{100, 0., 10.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hPtProng1RecSig", "3-prong candidates;prong 1 #it{p}_{T} (GeV/#it{c});;entries", {HistType::kTH2F, {{100, 0., 10.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hPtProng2RecSig", "3-prong candidates;prong 2 #it{p}_{T} (GeV/#it{c});;entries", {HistType::kTH2F, {{100, 0., 10.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    registry.add("MC/reconstructed/background/hDecLengthRecBg", "3-prong candidates;decay length (cm);;entries", {HistType::kTH2F, {{200, 0., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/background/hd0Prong0RecBg", "3-prong candidates;prong 0 DCAxy to prim. vertex (cm);;entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/background/hd0Prong1RecBg", "3-prong candidates;prong 1 DCAxy to prim. vertex (cm);;entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/background/hd0Prong2RecBg", "3-prong candidates;prong 1 DCAxy to prim. vertex (cm);;entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/background/hCtRecBg", "3-prong candidates;proper lifetime (#Xi_{c}) * #it{c} (cm);;entries", {HistType::kTH2F, {{120, -20., 100.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/background/hCPARecBg", "3-prong candidates;cosine of pointing angle;;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/background/hEtaRecBg", "3-prong candidates;candidate #it{#eta};;entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/background/hPtProng0RecBg", "3-prong candidates;prong 0 #it{p}_{T} (GeV/#it{c});;entries", {HistType::kTH2F, {{100, 0., 10.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/background/hPtProng1RecBg", "3-prong candidates;prong 1 #it{p}_{T} (GeV/#it{c});;entries", {HistType::kTH2F, {{100, 0., 10.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/background/hPtProng2RecBg", "3-prong candidates;prong 2 #it{p}_{T} (GeV/#it{c});;entries", {HistType::kTH2F, {{100, 0., 10.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
  }

  void process(const o2::aod::Collision& collision, const soa::Join<aod::Tracks, aod::TracksDCA>& tracks, soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelXicToPKPi>> const& candidates, aod::BigTracksPID const&)
  {
    int nTracks = 0;

    if (collision.numContrib() > 1) {
      for (auto const& track : tracks) {
        if (std::abs(track.eta()) > 4.0) { // hard coded cut TO BE REFINED
          continue;
        }
        if (std::abs(track.dcaXY()) > 0.0025 || std::abs(track.dcaZ()) > 0.0025) { // hardcoded cut TO BE REFINED
          continue;
        }
        nTracks++;
      }
    }

    registry.fill(HIST("Data/hMultiplicity"), nTracks); // filling the histo for multiplicity
    for (auto const& candidate : candidates) {
      auto candidatePt = candidate.pt();
      if (!(candidate.hfflag() & 1 << DecayType::XicToPKPi)) {
        continue;
      }

      if (yCandMax >= 0. && std::abs(yXic(candidate)) > yCandMax) {
        continue;
      }

      if (candidate.isSelXicToPKPi() >= selectionFlagXic) { // pKpi
        registry.fill(HIST("Data/hMassVsPt"), invMassXicToPKPi(candidate), candidatePt);
        registry.fill(HIST("Data/hMass"), invMassXicToPKPi(candidate));
      }
      if (candidate.isSelXicToPiKP() >= selectionFlagXic) { // piKp
        registry.fill(HIST("Data/hMassVsPt"), invMassXicToPiKP(candidate), candidatePt);
        registry.fill(HIST("Data/hMass"), invMassXicToPiKP(candidate));
      }
      registry.fill(HIST("Data/hPt"), candidatePt);
      registry.fill(HIST("Data/hEta"), candidate.eta());
      registry.fill(HIST("Data/hPhi"), candidate.phi());
      registry.fill(HIST("Data/hPtProng0VsPt"), candidate.ptProng0(), candidatePt);
      registry.fill(HIST("Data/hPtProng1VsPt"), candidate.ptProng1(), candidatePt);
      registry.fill(HIST("Data/hPtProng2VsPt"), candidate.ptProng2(), candidatePt);
      registry.fill(HIST("Data/hDecLengthVsPt"), candidate.decayLength(), candidatePt);
      registry.fill(HIST("Data/hd0Prong0VsPt"), candidate.impactParameter0(), candidatePt);
      registry.fill(HIST("Data/hd0Prong1VsPt"), candidate.impactParameter1(), candidatePt);
      registry.fill(HIST("Data/hd0Prong2VsPt"), candidate.impactParameter2(), candidatePt);
      registry.fill(HIST("Data/hCtVsPt"), ctXic(candidate), candidatePt);
      registry.fill(HIST("Data/hCPAVsPt"), candidate.cpa(), candidatePt);
      registry.fill(HIST("Data/hEtaVsPt"), candidate.eta(), candidatePt);
      registry.fill(HIST("Data/hSelectionStatusVsPt"), candidate.isSelXicToPKPi(), candidatePt);
      registry.fill(HIST("Data/hSelectionStatusVsPt"), candidate.isSelXicToPiKP(), candidatePt);
      registry.fill(HIST("Data/hImpParErr0VsPt"), candidate.errorImpactParameter0(), candidatePt);
      registry.fill(HIST("Data/hImpParErr1VsPt"), candidate.errorImpactParameter1(), candidatePt);
      registry.fill(HIST("Data/hImpParErr2VsPt"), candidate.errorImpactParameter2(), candidatePt);
      registry.fill(HIST("Data/hDecLenErrVsPt"), candidate.errorDecayLength(), candidatePt);
      registry.fill(HIST("Data/hChi2PCAVsPt"), candidate.chi2PCA(), candidatePt);

      const auto& trackProng0 = candidate.prong0_as<aod::BigTracksPID>(); // bachelor track
      const auto& trackProng1 = candidate.prong1_as<aod::BigTracksPID>(); // bachelor track
      const auto& trackProng2 = candidate.prong2_as<aod::BigTracksPID>(); // bachelor track

      auto prong0_P = trackProng0.p();
      auto prong1_P = trackProng1.p();
      auto prong2_P = trackProng2.p();

      // TPC nSigma histograms
      registry.fill(HIST("Data/hPVsTPCNSigmaPr_Prong0"), prong0_P, trackProng0.tpcNSigmaPr());
      registry.fill(HIST("Data/hPVsTPCNSigmaPi_Prong0"), prong0_P, trackProng0.tpcNSigmaPi());
      registry.fill(HIST("Data/hPVsTPCNSigmaKa_Prong0"), prong0_P, trackProng0.tpcNSigmaKa());

      registry.fill(HIST("Data/hPVsTPCNSigmaPr_Prong1"), prong1_P, trackProng1.tpcNSigmaPr());
      registry.fill(HIST("Data/hPVsTPCNSigmaPi_Prong1"), prong1_P, trackProng1.tpcNSigmaPi());
      registry.fill(HIST("Data/hPVsTPCNSigmaKa_Prong1"), prong1_P, trackProng1.tpcNSigmaKa());

      registry.fill(HIST("Data/hPVsTPCNSigmaPr_Prong2"), prong2_P, trackProng2.tpcNSigmaPr());
      registry.fill(HIST("Data/hPVsTPCNSigmaPi_Prong2"), prong2_P, trackProng2.tpcNSigmaPi());
      registry.fill(HIST("Data/hPVsTPCNSigmaKa_Prong2"), prong2_P, trackProng2.tpcNSigmaKa());

      // TOF nSigma histograms
      registry.fill(HIST("Data/hPVsTOFNSigmaPr_Prong0"), prong0_P, trackProng0.tofNSigmaPr());
      registry.fill(HIST("Data/hPVsTOFNSigmaPi_Prong0"), prong0_P, trackProng0.tofNSigmaPi());
      registry.fill(HIST("Data/hPVsTOFNSigmaKa_Prong0"), prong0_P, trackProng0.tofNSigmaKa());

      registry.fill(HIST("Data/hPVsTOFNSigmaPr_Prong1"), prong1_P, trackProng1.tofNSigmaPr());
      registry.fill(HIST("Data/hPVsTOFNSigmaPi_Prong1"), prong1_P, trackProng1.tofNSigmaPi());
      registry.fill(HIST("Data/hPVsTOFNSigmaKa_Prong1"), prong1_P, trackProng1.tofNSigmaKa());

      registry.fill(HIST("Data/hPVsTOFNSigmaPr_Prong2"), prong2_P, trackProng2.tofNSigmaPr());
      registry.fill(HIST("Data/hPVsTOFNSigmaPi_Prong2"), prong2_P, trackProng2.tofNSigmaPi());
      registry.fill(HIST("Data/hPVsTOFNSigmaKa_Prong2"), prong2_P, trackProng2.tofNSigmaKa());
    }
  }
  // Fill MC histograms
  void processMc(soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelXicToPKPi, aod::HfCand3ProngMcRec>> const& candidates,
                 soa::Join<aod::McParticles, aod::HfCand3ProngMcGen> const& particlesMC, aod::BigTracksMC const& /*tracks*/)
  {

    // MC rec.
    for (auto& candidate : candidates) {
      // Selected Xic
      if (!(candidate.hfflag() & 1 << DecayType::XicToPKPi)) {
        continue;
      } // rapidity selection
      if (yCandMax >= 0. && std::abs(yXic(candidate)) > yCandMax) {
        continue;
      }

      auto mass = 0.;
      auto massC = 0.;
      auto candidatePt = candidate.pt();

      if (candidate.isSelXicToPKPi() >= selectionFlagXic) {
        mass = invMassXicToPKPi(candidate);
      }
      if (candidate.isSelXicToPiKP() >= selectionFlagXic) {
        massC = invMassXicToPiKP(candidate); // mass conjiugated
      }

      if (std::abs(candidate.flagMcMatchRec()) == 1 << DecayType::XicToPKPi) {
        // Signal
        auto indexMother = RecoDecay::getMother(particlesMC, candidate.prong0_as<aod::BigTracksMC>().mcParticle_as<soa::Join<aod::McParticles, aod::HfCand3ProngMcGen>>(), pdg::Code::kXiCPlus, true);
        auto particleMother = particlesMC.rawIteratorAt(indexMother);

        registry.fill(HIST("MC/generated/signal/hPtGenSig"), particleMother.pt()); // gen. level pT
        registry.fill(HIST("MC/reconstructed/signal/hPtRecSig"), candidatePt);     // rec. level pT

        if (mass != 0.) {
          registry.fill(HIST("MC/reconstructed/signal/hMassSig"), mass, candidatePt);
        }
        if (massC != 0.) {
          registry.fill(HIST("MC/reconstructed/signal/hMassSig"), massC, candidatePt);
        }

        registry.fill(HIST("MC/reconstructed/signal/hDecLengthRecSig"), candidate.decayLength(), candidatePt);
        registry.fill(HIST("MC/reconstructed/signal/hPtProng0RecSig"), candidate.ptProng0(), candidatePt);
        registry.fill(HIST("MC/reconstructed/signal/hPtProng1RecSig"), candidate.ptProng1(), candidatePt);
        registry.fill(HIST("MC/reconstructed/signal/hPtProng2RecSig"), candidate.ptProng2(), candidatePt);
        registry.fill(HIST("MC/reconstructed/signal/hd0Prong0RecSig"), candidate.impactParameter0(), candidatePt);
        registry.fill(HIST("MC/reconstructed/signal/hd0Prong1RecSig"), candidate.impactParameter1(), candidatePt);
        registry.fill(HIST("MC/reconstructed/signal/hd0Prong2RecSig"), candidate.impactParameter2(), candidatePt);
        registry.fill(HIST("MC/reconstructed/signal/hCtRecSig"), ctXic(candidate), candidatePt);
        registry.fill(HIST("MC/reconstructed/signal/hCPARecSig"), candidate.cpa(), candidatePt);
        registry.fill(HIST("MC/reconstructed/signal/hEtaRecSig"), candidate.eta(), candidatePt);
      } else {
        // Background
        registry.fill(HIST("MC/reconstructed/background/hPtRecBg"), candidatePt);

        if (mass != 0.) {
          registry.fill(HIST("MC/reconstructed/background/hMassBg"), mass, candidatePt);
        }
        if (massC != 0.) {
          registry.fill(HIST("MC/reconstructed/background/hMassBg"), massC, candidatePt);
        }

        registry.fill(HIST("MC/reconstructed/background/hDecLengthRecBg"), candidate.decayLength(), candidatePt);
        registry.fill(HIST("MC/reconstructed/background/hPtProng0RecBg"), candidate.ptProng0(), candidatePt);
        registry.fill(HIST("MC/reconstructed/background/hPtProng1RecBg"), candidate.ptProng1(), candidatePt);
        registry.fill(HIST("MC/reconstructed/background/hPtProng2RecBg"), candidate.ptProng2(), candidatePt);
        registry.fill(HIST("MC/reconstructed/background/hd0Prong0RecBg"), candidate.impactParameter0(), candidatePt);
        registry.fill(HIST("MC/reconstructed/background/hd0Prong1RecBg"), candidate.impactParameter1(), candidatePt);
        registry.fill(HIST("MC/reconstructed/background/hd0Prong2RecBg"), candidate.impactParameter2(), candidatePt);
        registry.fill(HIST("MC/reconstructed/background/hCtRecBg"), ctXic(candidate), candidatePt);
        registry.fill(HIST("MC/reconstructed/background/hCPARecBg"), candidate.cpa(), candidatePt);
        registry.fill(HIST("MC/reconstructed/background/hEtaRecBg"), candidate.eta(), candidatePt);
      }
    }
    // MC gen.
    for (auto& particle : particlesMC) {
      if (std::abs(particle.flagMcMatchGen()) == 1 << DecayType::XicToPKPi) {
        if (yCandMax >= 0. && std::abs(RecoDecay::y(array{particle.px(), particle.py(), particle.pz()}, RecoDecay::getMassPDG(particle.pdgCode()))) > yCandMax) {
          continue;
        }
        registry.fill(HIST("MC/generated/signal/hPtGen"), particle.pt());
        registry.fill(HIST("MC/generated/signal/hEtaGen"), particle.eta(), particle.pt());
      }
    }
  }

  PROCESS_SWITCH(HfTaskXic, processMc, "Process MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskXic>(cfgc)};
}
