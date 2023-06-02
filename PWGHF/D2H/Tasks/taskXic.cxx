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
/// \note Inspired from taskLc.cxx
///
/// \author Mattia Faggin <mattia.faggin@cern.ch>, University and INFN PADOVA
/// \author Anton Alkin <anton.alkin@cern.ch>, CERN
/// \author Jinjoo Seo <jin.joo.seo@cern.ch>, Inha University
/// \author Himanshu Sharma <himanshu.sharma@cern.ch>, University and INFN Padova
/// \author Cristina Terrevoli <cristina.terrevoli@cern.ch>, INFN Bari

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::hf_cand_3prong;
using namespace o2::analysis::hf_cuts_xic_to_p_k_pi;

/// Ξc± analysis task
struct HfTaskXic {
  Configurable<int> selectionFlagXic{"selectionFlagXic", 1, "Selection Flag for Xic"};
  Configurable<double> yCandMax{"yCandMax", -1., "max. cand. rapidity"};
  Configurable<float> etaTrackMax{"etaTrackMax", 4.0, "max. track eta"};
  Configurable<float> dcaXYTrackMax{"dcaXYTrackMax", 0.0025, "max. DCAxy for track"};
  Configurable<float> dcaZTrackMax{"dcaZTrackMax", 0.0025, "max. DCAz for track"};
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_xic_to_p_k_pi::vecBinsPt}, "pT bin limits"};

  Filter filterSelectCandidates = (aod::hf_sel_candidate_xic::isSelXicToPKPi >= selectionFlagXic || aod::hf_sel_candidate_xic::isSelXicToPiKP >= selectionFlagXic);

  Partition<soa::Join<aod::HfCand3Prong, aod::HfSelXicToPKPi, aod::HfCand3ProngMcRec>> selectedMCXicCandidates = (aod::hf_sel_candidate_xic::isSelXicToPKPi >= selectionFlagXic || aod::hf_sel_candidate_xic::isSelXicToPiKP >= selectionFlagXic);

  HistogramRegistry registry{
    "registry", // histo not in pt bins
    {
      {"Data/hPt", "3-prong candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},    // pt Xic
      {"Data/hEta", "3-prong candidates;candidate #it{eta};entries", {HistType::kTH1F, {{100, -5., 5.}}}},                  // eta Xic
      {"Data/hPhi", "3-prong candidates;candidate #varphi;entries", {HistType::kTH1F, {{72, 0., constants::math::TwoPI}}}}, // phi Xic
      {"Data/hMass", "3-prong candidates; inv. mass (p K #pi) (GeV/#it{c}^{2})", {HistType::kTH1F, {{600, 2.18, 2.58}}}},   // mass Xic
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

    AxisSpec axisPPid = {100, 0.f, 10.0f, "#it{p} (GeV/#it{c})"};
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
    registry.add("Data/hPVsTPCNSigmaPr_Prong0", "3-prong candidates;#it{p} bachelor (GeV/#it{c}) ;n#it{#sigma}_{p} TPC", {HistType::kTH2F, {axisPPid, axisNSigmaPr}});
    registry.add("Data/hPVsTPCNSigmaPi_Prong0", "3-prong candidates;#it{p} bachelor (GeV/#it{c}) ;n#it{#sigma}_{#pi} TPC", {HistType::kTH2F, {axisPPid, axisNSigmaPi}});
    registry.add("Data/hPVsTPCNSigmaKa_Prong0", "3-prong candidates;#it{p} bachelor (GeV/#it{c}) ;n#it{#sigma}_{K} TPC", {HistType::kTH2F, {axisPPid, axisNSigmaKa}});

    registry.add("Data/hPVsTPCNSigmaPr_Prong1", "3-prong candidates;#it{p} bachelor (GeV/#it{c}) ;n#it{#sigma}_{p} TPC", {HistType::kTH2F, {axisPPid, axisNSigmaPr}});
    registry.add("Data/hPVsTPCNSigmaPi_Prong1", "3-prong candidates;#it{p} bachelor (GeV/#it{c}) ;n#it{#sigma}_{#pi} TPC", {HistType::kTH2F, {axisPPid, axisNSigmaPi}});
    registry.add("Data/hPVsTPCNSigmaKa_Prong1", "3-prong candidates;#it{p} bachelor (GeV/#it{c}) ;n#it{#sigma}_{K} TPC", {HistType::kTH2F, {axisPPid, axisNSigmaKa}});

    registry.add("Data/hPVsTPCNSigmaPr_Prong2", "3-prong candidates;#it{p} bachelor (GeV/#it{c}) ;n#it{#sigma}_{p} TPC", {HistType::kTH2F, {axisPPid, axisNSigmaPr}});
    registry.add("Data/hPVsTPCNSigmaPi_Prong2", "3-prong candidates;#it{p} bachelor (GeV/#it{c}) ;n#it{#sigma}_{#pi} TPC", {HistType::kTH2F, {axisPPid, axisNSigmaPi}});
    registry.add("Data/hPVsTPCNSigmaKa_Prong2", "3-prong candidates;#it{p} bachelor (GeV/#it{c}) ;n#it{#sigma}_{K} TPC", {HistType::kTH2F, {axisPPid, axisNSigmaKa}});

    // TOF nSigma histograms
    registry.add("Data/hPVsTOFNSigmaPr_Prong0", "3-prong candidates;#it{p} bachelor (GeV/#it{c}) ;n#it{#sigma}_{p} TOF", {HistType::kTH2F, {axisPPid, axisNSigmaPr}});
    registry.add("Data/hPVsTOFNSigmaPi_Prong0", "3-prong candidates;#it{p} bachelor (GeV/#it{c}) ;n#it{#sigma}_{#pi} TOF", {HistType::kTH2F, {axisPPid, axisNSigmaPi}});
    registry.add("Data/hPVsTOFNSigmaKa_Prong0", "3-prong candidates;#it{p} bachelor (GeV/#it{c}) ;n#it{#sigma}_{K} TOF", {HistType::kTH2F, {axisPPid, axisNSigmaKa}});

    registry.add("Data/hPVsTOFNSigmaPr_Prong1", "3-prong candidates;#it{p} bachelor (GeV/#it{c}) ;n#it{#sigma}_{p} TOF", {HistType::kTH2F, {axisPPid, axisNSigmaPr}});
    registry.add("Data/hPVsTOFNSigmaPi_Prong1", "3-prong candidates;#it{p} bachelor (GeV/#it{c}) ;n#it{#sigma}_{#pi} TOF", {HistType::kTH2F, {axisPPid, axisNSigmaPi}});
    registry.add("Data/hPVsTOFNSigmaKa_Prong1", "3-prong candidates;#it{p} bachelor (GeV/#it{c}) ;n#it{#sigma}_{K} TOF", {HistType::kTH2F, {axisPPid, axisNSigmaKa}});

    registry.add("Data/hPVsTOFNSigmaPr_Prong2", "3-prong candidates;#it{p} bachelor (GeV/#it{c}) ;n#it{#sigma}_{p} TOF", {HistType::kTH2F, {axisPPid, axisNSigmaPr}});
    registry.add("Data/hPVsTOFNSigmaPi_Prong2", "3-prong candidates;#it{p} bachelor (GeV/#it{c}) ;n#it{#sigma}_{#pi} TOF", {HistType::kTH2F, {axisPPid, axisNSigmaPi}});
    registry.add("Data/hPVsTOFNSigmaKa_Prong2", "3-prong candidates;#it{p} bachelor (GeV/#it{c}) ;n#it{#sigma}_{K} TOF", {HistType::kTH2F, {axisPPid, axisNSigmaKa}});

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

  void process(aod::Collision const& collision, soa::Join<aod::Tracks, aod::TracksDCA> const& tracks, soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelXicToPKPi>> const& candidates, aod::BigTracksPID const&)
  {
    int nTracks = 0;

    if (collision.numContrib() > 1) {
      for (auto const& track : tracks) {
        if (std::abs(track.eta()) > etaTrackMax) { // hard coded cut TO BE REFINED
          continue;
        }
        if (std::abs(track.dcaXY()) > dcaXYTrackMax || std::abs(track.dcaZ()) > dcaZTrackMax) { // hardcoded cut TO BE REFINED
          continue;
        }
        nTracks++;
      }
    }

    registry.fill(HIST("Data/hMultiplicity"), nTracks); // filling the histo for multiplicity
    for (auto const& candidate : candidates) {
      auto ptCandidate = candidate.pt();
      if (!(candidate.hfflag() & 1 << DecayType::XicToPKPi)) {
        continue;
      }

      if (yCandMax >= 0. && std::abs(yXic(candidate)) > yCandMax) {
        continue;
      }

      if (candidate.isSelXicToPKPi() >= selectionFlagXic) { // pKpi
        registry.fill(HIST("Data/hMassVsPt"), invMassXicToPKPi(candidate), ptCandidate);
        registry.fill(HIST("Data/hMass"), invMassXicToPKPi(candidate));
      }
      if (candidate.isSelXicToPiKP() >= selectionFlagXic) { // piKp
        registry.fill(HIST("Data/hMassVsPt"), invMassXicToPiKP(candidate), ptCandidate);
        registry.fill(HIST("Data/hMass"), invMassXicToPiKP(candidate));
      }
      registry.fill(HIST("Data/hPt"), ptCandidate);
      registry.fill(HIST("Data/hEta"), candidate.eta());
      registry.fill(HIST("Data/hPhi"), candidate.phi());
      registry.fill(HIST("Data/hPtProng0VsPt"), candidate.ptProng0(), ptCandidate);
      registry.fill(HIST("Data/hPtProng1VsPt"), candidate.ptProng1(), ptCandidate);
      registry.fill(HIST("Data/hPtProng2VsPt"), candidate.ptProng2(), ptCandidate);
      registry.fill(HIST("Data/hDecLengthVsPt"), candidate.decayLength(), ptCandidate);
      registry.fill(HIST("Data/hd0Prong0VsPt"), candidate.impactParameter0(), ptCandidate);
      registry.fill(HIST("Data/hd0Prong1VsPt"), candidate.impactParameter1(), ptCandidate);
      registry.fill(HIST("Data/hd0Prong2VsPt"), candidate.impactParameter2(), ptCandidate);
      registry.fill(HIST("Data/hCtVsPt"), ctXic(candidate), ptCandidate);
      registry.fill(HIST("Data/hCPAVsPt"), candidate.cpa(), ptCandidate);
      registry.fill(HIST("Data/hEtaVsPt"), candidate.eta(), ptCandidate);
      registry.fill(HIST("Data/hSelectionStatusVsPt"), candidate.isSelXicToPKPi(), ptCandidate);
      registry.fill(HIST("Data/hSelectionStatusVsPt"), candidate.isSelXicToPiKP(), ptCandidate);
      registry.fill(HIST("Data/hImpParErr0VsPt"), candidate.errorImpactParameter0(), ptCandidate);
      registry.fill(HIST("Data/hImpParErr1VsPt"), candidate.errorImpactParameter1(), ptCandidate);
      registry.fill(HIST("Data/hImpParErr2VsPt"), candidate.errorImpactParameter2(), ptCandidate);
      registry.fill(HIST("Data/hDecLenErrVsPt"), candidate.errorDecayLength(), ptCandidate);
      registry.fill(HIST("Data/hChi2PCAVsPt"), candidate.chi2PCA(), ptCandidate);

      const auto& trackProng0 = candidate.prong0_as<aod::BigTracksPID>(); // bachelor track
      const auto& trackProng1 = candidate.prong1_as<aod::BigTracksPID>(); // bachelor track
      const auto& trackProng2 = candidate.prong2_as<aod::BigTracksPID>(); // bachelor track

      auto momentumProng0 = trackProng0.p();
      auto momentumProng1 = trackProng1.p();
      auto momentumProng2 = trackProng2.p();

      // TPC nSigma histograms
      registry.fill(HIST("Data/hPVsTPCNSigmaPr_Prong0"), momentumProng0, trackProng0.tpcNSigmaPr());
      registry.fill(HIST("Data/hPVsTPCNSigmaPi_Prong0"), momentumProng0, trackProng0.tpcNSigmaPi());
      registry.fill(HIST("Data/hPVsTPCNSigmaKa_Prong0"), momentumProng0, trackProng0.tpcNSigmaKa());

      registry.fill(HIST("Data/hPVsTPCNSigmaPr_Prong1"), momentumProng1, trackProng1.tpcNSigmaPr());
      registry.fill(HIST("Data/hPVsTPCNSigmaPi_Prong1"), momentumProng1, trackProng1.tpcNSigmaPi());
      registry.fill(HIST("Data/hPVsTPCNSigmaKa_Prong1"), momentumProng1, trackProng1.tpcNSigmaKa());

      registry.fill(HIST("Data/hPVsTPCNSigmaPr_Prong2"), momentumProng2, trackProng2.tpcNSigmaPr());
      registry.fill(HIST("Data/hPVsTPCNSigmaPi_Prong2"), momentumProng2, trackProng2.tpcNSigmaPi());
      registry.fill(HIST("Data/hPVsTPCNSigmaKa_Prong2"), momentumProng2, trackProng2.tpcNSigmaKa());

      // TOF nSigma histograms
      registry.fill(HIST("Data/hPVsTOFNSigmaPr_Prong0"), momentumProng0, trackProng0.tofNSigmaPr());
      registry.fill(HIST("Data/hPVsTOFNSigmaPi_Prong0"), momentumProng0, trackProng0.tofNSigmaPi());
      registry.fill(HIST("Data/hPVsTOFNSigmaKa_Prong0"), momentumProng0, trackProng0.tofNSigmaKa());

      registry.fill(HIST("Data/hPVsTOFNSigmaPr_Prong1"), momentumProng1, trackProng1.tofNSigmaPr());
      registry.fill(HIST("Data/hPVsTOFNSigmaPi_Prong1"), momentumProng1, trackProng1.tofNSigmaPi());
      registry.fill(HIST("Data/hPVsTOFNSigmaKa_Prong1"), momentumProng1, trackProng1.tofNSigmaKa());

      registry.fill(HIST("Data/hPVsTOFNSigmaPr_Prong2"), momentumProng2, trackProng2.tofNSigmaPr());
      registry.fill(HIST("Data/hPVsTOFNSigmaPi_Prong2"), momentumProng2, trackProng2.tofNSigmaPi());
      registry.fill(HIST("Data/hPVsTOFNSigmaKa_Prong2"), momentumProng2, trackProng2.tofNSigmaKa());
    }
  }
  // Fill MC histograms
  void processMc(soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelXicToPKPi, aod::HfCand3ProngMcRec>> const& candidates,
                 soa::Join<aod::McParticles, aod::HfCand3ProngMcGen> const& particlesMC, aod::BigTracksMC const& /*tracks*/)
  {

    // MC rec.
    for (auto const& candidate : candidates) {
      // Selected Xic
      if (!(candidate.hfflag() & 1 << DecayType::XicToPKPi)) {
        continue;
      } // rapidity selection
      if (yCandMax >= 0. && std::abs(yXic(candidate)) > yCandMax) {
        continue;
      }

      auto massXicToPKPi = 0.;
      auto massXicToPiKP = 0.;
      auto ptCandidate = candidate.pt();

      if (candidate.isSelXicToPKPi() >= selectionFlagXic) {
        massXicToPKPi = invMassXicToPKPi(candidate);
      }
      if (candidate.isSelXicToPiKP() >= selectionFlagXic) {
        massXicToPiKP = invMassXicToPiKP(candidate); // mass conjugate
      }

      if (std::abs(candidate.flagMcMatchRec()) == 1 << DecayType::XicToPKPi) {
        // Signal
        auto indexMother = RecoDecay::getMother(particlesMC, candidate.prong0_as<aod::BigTracksMC>().mcParticle_as<soa::Join<aod::McParticles, aod::HfCand3ProngMcGen>>(), pdg::Code::kXiCPlus, true);
        auto particleMother = particlesMC.rawIteratorAt(indexMother);

        registry.fill(HIST("MC/generated/signal/hPtGenSig"), particleMother.pt()); // gen. level pT
        registry.fill(HIST("MC/reconstructed/signal/hPtRecSig"), ptCandidate);     // rec. level pT

        if (massXicToPKPi != 0.) {
          registry.fill(HIST("MC/reconstructed/signal/hMassSig"), massXicToPKPi, ptCandidate);
        }
        if (massXicToPiKP != 0.) {
          registry.fill(HIST("MC/reconstructed/signal/hMassSig"), massXicToPiKP, ptCandidate);
        }

        registry.fill(HIST("MC/reconstructed/signal/hDecLengthRecSig"), candidate.decayLength(), ptCandidate);
        registry.fill(HIST("MC/reconstructed/signal/hPtProng0RecSig"), candidate.ptProng0(), ptCandidate);
        registry.fill(HIST("MC/reconstructed/signal/hPtProng1RecSig"), candidate.ptProng1(), ptCandidate);
        registry.fill(HIST("MC/reconstructed/signal/hPtProng2RecSig"), candidate.ptProng2(), ptCandidate);
        registry.fill(HIST("MC/reconstructed/signal/hd0Prong0RecSig"), candidate.impactParameter0(), ptCandidate);
        registry.fill(HIST("MC/reconstructed/signal/hd0Prong1RecSig"), candidate.impactParameter1(), ptCandidate);
        registry.fill(HIST("MC/reconstructed/signal/hd0Prong2RecSig"), candidate.impactParameter2(), ptCandidate);
        registry.fill(HIST("MC/reconstructed/signal/hCtRecSig"), ctXic(candidate), ptCandidate);
        registry.fill(HIST("MC/reconstructed/signal/hCPARecSig"), candidate.cpa(), ptCandidate);
        registry.fill(HIST("MC/reconstructed/signal/hEtaRecSig"), candidate.eta(), ptCandidate);
      } else {
        // Background
        registry.fill(HIST("MC/reconstructed/background/hPtRecBg"), ptCandidate);

        if (massXicToPKPi != 0.) {
          registry.fill(HIST("MC/reconstructed/background/hMassBg"), massXicToPKPi, ptCandidate);
        }
        if (massXicToPiKP != 0.) {
          registry.fill(HIST("MC/reconstructed/background/hMassBg"), massXicToPiKP, ptCandidate);
        }

        registry.fill(HIST("MC/reconstructed/background/hDecLengthRecBg"), candidate.decayLength(), ptCandidate);
        registry.fill(HIST("MC/reconstructed/background/hPtProng0RecBg"), candidate.ptProng0(), ptCandidate);
        registry.fill(HIST("MC/reconstructed/background/hPtProng1RecBg"), candidate.ptProng1(), ptCandidate);
        registry.fill(HIST("MC/reconstructed/background/hPtProng2RecBg"), candidate.ptProng2(), ptCandidate);
        registry.fill(HIST("MC/reconstructed/background/hd0Prong0RecBg"), candidate.impactParameter0(), ptCandidate);
        registry.fill(HIST("MC/reconstructed/background/hd0Prong1RecBg"), candidate.impactParameter1(), ptCandidate);
        registry.fill(HIST("MC/reconstructed/background/hd0Prong2RecBg"), candidate.impactParameter2(), ptCandidate);
        registry.fill(HIST("MC/reconstructed/background/hCtRecBg"), ctXic(candidate), ptCandidate);
        registry.fill(HIST("MC/reconstructed/background/hCPARecBg"), candidate.cpa(), ptCandidate);
        registry.fill(HIST("MC/reconstructed/background/hEtaRecBg"), candidate.eta(), ptCandidate);
      }
    }
    // MC gen.
    for (auto const& particle : particlesMC) {
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
