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
#include "TTree.h"
#include "Common/DataModel/PIDResponse.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::hf_cand_3prong;
using namespace o2::analysis::hf_cuts_xic_to_p_k_pi;

#include "Framework/runDataProcessing.h"

/// Ξc± analysis task

struct XicTTree {
  float fMass;
  float fPt;
  float fEta;
  float fCpa;
  float fCpaXY;
  float fDecayLength;
  float fDecayLengthErr;
  float fChi2PCA;

  float fPidTpcNsigPrProng0;
  float fPidTpcNsigPiProng0;
  float fPidTpcNsigKProng0;
  float fPidTpcNsigPrProng1;
  float fPidTpcNsigPiProng1;
  float fPidTpcNsigKProng1;
  float fPidTpcNsigPrProng2;
  float fPidTpcNsigPiProng2;
  float fPidTpcNsigKProng2;

  float fPidTofNsigPrProng0;
  float fPidTofNsigPiProng0;
  float fPidTofNsigKProng0;
  float fPidTofNsigPrProng1;
  float fPidTofNsigPiProng1;
  float fPidTofNsigKProng1;
  float fPidTofNsigPrProng2;
  float fPidTofNsigPiProng2;
  float fPidTofNsigKProng2;
};

struct HfTaskXic {

  OutputObj<TTree> tree{"XicTTree"};
  XicTTree xictreeObj;

  Configurable<int> selectionFlagXic{"selectionFlagXic", 1, "Selection Flag for Xic"};
  Configurable<double> yCandMax{"yCandMax", -1., "max. cand. rapidity"};
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_xic_to_p_k_pi::vecBinsPt}, "pT bin limits"};
  Configurable<bool> fConfigWriteTTree{"cfgWriteTTree", false, "write tree output"};
  Configurable<bool> fFillTOFnSigma{"fFillTOFnSigma", false, "write TOF nSigma info to tree"};
  Configurable<bool> fFillTPCnSigma{"fFillTPCnSigma", false, "write TPC nSigma info to tree"};

  Filter filterSelectCandidates = (aod::hf_sel_candidate_xic::isSelXicToPKPi >= selectionFlagXic || aod::hf_sel_candidate_xic::isSelXicToPiKP >= selectionFlagXic);

  Partition<soa::Join<aod::HfCand3Prong, aod::HfSelXicToPKPi, aod::HfCand3ProngMcRec>> selectedMCXicCandidates = (aod::hf_sel_candidate_xic::isSelXicToPKPi >= selectionFlagXic || aod::hf_sel_candidate_xic::isSelXicToPiKP >= selectionFlagXic);

  HistogramRegistry registry{
    "registry", // histo not in pt bins
    {
      {"Data/hPt", "3-prong candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},                                       // pt Xic
      {"MC/reconstructed/signal/hPtRecCand", "3-prong candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},             // pt Xic
      {"Data/hEta", "3-prong candidates;candidate #it{eta};entries", {HistType::kTH1F, {{100, -5., 5.}}}},                                                     // eta Xic
      {"Data/hPhi", "3-prong candidates;candidate #varphi;entries", {HistType::kTH1F, {{72, 0., 2. * constants::math::PI}}}},                                  // phi Xic
      {"Data/hMass", "3-prong candidates; inv. mass (p K #pi) (GeV/#it{c}^{2})", {HistType::kTH1F, {{600, 2.18, 2.58}}}},                                      // mass Xic
      {"MC/reconstructed/hMassRecCand", "3-prong candidates (matched, prompt); inv. mass (p K #pi) (GeV/#it{c}^{2})", {HistType::kTH1F, {{600, 2.18, 2.58}}}}, // mass Xic
      {"MC/reconstructed/hPtRecSig", "3-prong candidates (matched);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}}},
      {"MC/reconstructed/hPtRecBg", "3-prong candidates (unmatched);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}}},
      {"MC/generated/hPtGen", "MC particles (matched);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}}},
      {"MC/generated/hPtGenSig", "3-prong candidates (matched);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}}}, ///
      {"Data/hMultiplicity", "multiplicity;multiplicity;entries", {HistType::kTH1F, {{1000, 0., 1000.}}}}}};

  void init(o2::framework::InitContext&)
  {
    if (fConfigWriteTTree) {
      SetTree();
    }

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
    registry.add("MC/reconstructed/hMassBg", "Invariant mass (unmatched);m (p K #pi) (GeV/#it{c}^{2});;entries", {HistType::kTH2F, {{500, 1.6, 3.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/generated/hEtaGen", "3-prong candidates;candidate #it{#eta};;entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    registry.add("MC/reconstructed/hDecLengthRecSig", "3-prong candidates;decay length (cm);;entries", {HistType::kTH2F, {{200, 0., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/hd0Prong0RecSig", "3-prong candidates;prong 0 DCAxy to prim. vertex (cm);;entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/hd0Prong1RecSig", "3-prong candidates;prong 1 DCAxy to prim. vertex (cm);;entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/hd0Prong2RecSig", "3-prong candidates;prong 1 DCAxy to prim. vertex (cm);;entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/hCtRecSig", "3-prong candidates;proper lifetime (#Xi_{c}) * #it{c} (cm);;entries", {HistType::kTH2F, {{120, -20., 100.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/hCPARecSig", "3-prong candidates;cosine of pointing angle;;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/hEtaRecSig", "3-prong candidates;candidate #it{#eta};;entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/hPtProng0RecSig", "3-prong candidates;prong 0 #it{p}_{T} (GeV/#it{c});;entries", {HistType::kTH2F, {{100, 0., 10.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/hPtProng1RecSig", "3-prong candidates;prong 1 #it{p}_{T} (GeV/#it{c});;entries", {HistType::kTH2F, {{100, 0., 10.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/hPtProng2RecSig", "3-prong candidates;prong 2 #it{p}_{T} (GeV/#it{c});;entries", {HistType::kTH2F, {{100, 0., 10.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    registry.add("MC/reconstructed/hDecLengthRecBg", "3-prong candidates;decay length (cm);;entries", {HistType::kTH2F, {{200, 0., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/hd0Prong0RecBg", "3-prong candidates;prong 0 DCAxy to prim. vertex (cm);;entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/hd0Prong1RecBg", "3-prong candidates;prong 1 DCAxy to prim. vertex (cm);;entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/hd0Prong2RecBg", "3-prong candidates;prong 1 DCAxy to prim. vertex (cm);;entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/hCtRecBg", "3-prong candidates;proper lifetime (#Xi_{c}) * #it{c} (cm);;entries", {HistType::kTH2F, {{120, -20., 100.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/hCPARecBg", "3-prong candidates;cosine of pointing angle;;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/hEtaRecBg", "3-prong candidates;candidate #it{#eta};;entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/hPtProng0RecBg", "3-prong candidates;prong 0 #it{p}_{T} (GeV/#it{c});;entries", {HistType::kTH2F, {{100, 0., 10.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/hPtProng1RecBg", "3-prong candidates;prong 1 #it{p}_{T} (GeV/#it{c});;entries", {HistType::kTH2F, {{100, 0., 10.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/hPtProng2RecBg", "3-prong candidates;prong 2 #it{p}_{T} (GeV/#it{c});;entries", {HistType::kTH2F, {{100, 0., 10.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
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
      if (!(candidate.hfflag() & 1 << DecayType::XicToPKPi)) {
        continue;
      }

      if (yCandMax >= 0. && std::abs(yXic(candidate)) > yCandMax) {
        continue;
      }

      if (candidate.isSelXicToPKPi() >= selectionFlagXic) { // pKpi
        registry.fill(HIST("Data/hMassVsPt"), invMassXicToPKPi(candidate), candidate.pt());
        registry.fill(HIST("Data/hMass"), invMassXicToPKPi(candidate));
        xictreeObj.fMass = invMassXicToPKPi(candidate);
      }
      if (candidate.isSelXicToPiKP() >= selectionFlagXic) { // piKp
        registry.fill(HIST("Data/hMassVsPt"), invMassXicToPiKP(candidate), candidate.pt());
        registry.fill(HIST("Data/hMass"), invMassXicToPiKP(candidate));
        xictreeObj.fMass = invMassXicToPiKP(candidate);
      }
      registry.fill(HIST("Data/hPt"), candidate.pt());
      registry.fill(HIST("Data/hEta"), candidate.eta());
      registry.fill(HIST("Data/hPhi"), candidate.phi());
      registry.fill(HIST("Data/hPtProng0VsPt"), candidate.ptProng0(), candidate.pt());
      registry.fill(HIST("Data/hPtProng1VsPt"), candidate.ptProng1(), candidate.pt());
      registry.fill(HIST("Data/hPtProng2VsPt"), candidate.ptProng2(), candidate.pt());
      registry.fill(HIST("Data/hDecLengthVsPt"), candidate.decayLength(), candidate.pt());
      registry.fill(HIST("Data/hd0Prong0VsPt"), candidate.impactParameter0(), candidate.pt());
      registry.fill(HIST("Data/hd0Prong1VsPt"), candidate.impactParameter1(), candidate.pt());
      registry.fill(HIST("Data/hd0Prong2VsPt"), candidate.impactParameter2(), candidate.pt());
      registry.fill(HIST("Data/hCtVsPt"), ctXic(candidate), candidate.pt());
      registry.fill(HIST("Data/hCPAVsPt"), candidate.cpa(), candidate.pt());
      registry.fill(HIST("Data/hEtaVsPt"), candidate.eta(), candidate.pt());
      registry.fill(HIST("Data/hSelectionStatusVsPt"), candidate.isSelXicToPKPi(), candidate.pt());
      registry.fill(HIST("Data/hSelectionStatusVsPt"), candidate.isSelXicToPiKP(), candidate.pt());
      registry.fill(HIST("Data/hImpParErr0VsPt"), candidate.errorImpactParameter0(), candidate.pt());
      registry.fill(HIST("Data/hImpParErr1VsPt"), candidate.errorImpactParameter1(), candidate.pt());
      registry.fill(HIST("Data/hImpParErr2VsPt"), candidate.errorImpactParameter2(), candidate.pt());
      registry.fill(HIST("Data/hDecLenErrVsPt"), candidate.errorDecayLength(), candidate.pt());
      registry.fill(HIST("Data/hChi2PCAVsPt"), candidate.chi2PCA(), candidate.pt());

      const auto& trackProng0 = candidate.prong0_as<aod::BigTracksPID>(); // bachelor track
      const auto& trackProng1 = candidate.prong1_as<aod::BigTracksPID>(); // bachelor track
      const auto& trackProng2 = candidate.prong2_as<aod::BigTracksPID>(); // bachelor track

      // TPC nSigma histograms
      registry.fill(HIST("Data/hPVsTPCNSigmaPr_Prong0"), trackProng0.p(), trackProng0.tpcNSigmaPr());
      registry.fill(HIST("Data/hPVsTPCNSigmaPi_Prong0"), trackProng0.p(), trackProng0.tpcNSigmaPi());
      registry.fill(HIST("Data/hPVsTPCNSigmaKa_Prong0"), trackProng0.p(), trackProng0.tpcNSigmaKa());

      registry.fill(HIST("Data/hPVsTPCNSigmaPr_Prong1"), trackProng1.p(), trackProng1.tpcNSigmaPr());
      registry.fill(HIST("Data/hPVsTPCNSigmaPi_Prong1"), trackProng1.p(), trackProng1.tpcNSigmaPi());
      registry.fill(HIST("Data/hPVsTPCNSigmaKa_Prong1"), trackProng1.p(), trackProng1.tpcNSigmaKa());

      registry.fill(HIST("Data/hPVsTPCNSigmaPr_Prong2"), trackProng2.p(), trackProng2.tpcNSigmaPr());
      registry.fill(HIST("Data/hPVsTPCNSigmaPi_Prong2"), trackProng2.p(), trackProng2.tpcNSigmaPi());
      registry.fill(HIST("Data/hPVsTPCNSigmaKa_Prong2"), trackProng2.p(), trackProng2.tpcNSigmaKa());

      // TOF nSigma histograms
      registry.fill(HIST("Data/hPVsTOFNSigmaPr_Prong0"), trackProng0.p(), trackProng0.tofNSigmaPr());
      registry.fill(HIST("Data/hPVsTOFNSigmaPi_Prong0"), trackProng0.p(), trackProng0.tofNSigmaPi());
      registry.fill(HIST("Data/hPVsTOFNSigmaKa_Prong0"), trackProng0.p(), trackProng0.tofNSigmaKa());

      registry.fill(HIST("Data/hPVsTOFNSigmaPr_Prong1"), trackProng1.p(), trackProng1.tofNSigmaPr());
      registry.fill(HIST("Data/hPVsTOFNSigmaPi_Prong1"), trackProng1.p(), trackProng1.tofNSigmaPi());
      registry.fill(HIST("Data/hPVsTOFNSigmaKa_Prong1"), trackProng1.p(), trackProng1.tofNSigmaKa());

      registry.fill(HIST("Data/hPVsTOFNSigmaPr_Prong2"), trackProng2.p(), trackProng2.tofNSigmaPr());
      registry.fill(HIST("Data/hPVsTOFNSigmaPi_Prong2"), trackProng2.p(), trackProng2.tofNSigmaPi());
      registry.fill(HIST("Data/hPVsTOFNSigmaKa_Prong2"), trackProng2.p(), trackProng2.tofNSigmaKa());

      // fill TTree varaibles
      xictreeObj.fPt = candidate.pt();
      xictreeObj.fCpa = candidate.cpa();
      xictreeObj.fCpaXY = candidate.cpaXY();
      xictreeObj.fDecayLength = candidate.decayLength();
      xictreeObj.fDecayLengthErr = candidate.errorDecayLength();
      xictreeObj.fChi2PCA = candidate.chi2PCA();

      if (fFillTOFnSigma) {
        xictreeObj.fPidTofNsigPrProng0 = trackProng0.tofNSigmaPr();
        xictreeObj.fPidTofNsigPiProng0 = trackProng0.tofNSigmaPi();
        xictreeObj.fPidTofNsigKProng0 = trackProng0.tofNSigmaKa();

        xictreeObj.fPidTofNsigPrProng1 = trackProng1.tofNSigmaPr();
        xictreeObj.fPidTofNsigPiProng1 = trackProng1.tofNSigmaPi();
        xictreeObj.fPidTofNsigKProng1 = trackProng1.tofNSigmaKa();

        xictreeObj.fPidTofNsigPrProng2 = trackProng2.tofNSigmaPr();
        xictreeObj.fPidTofNsigPiProng2 = trackProng2.tofNSigmaPi();
        xictreeObj.fPidTofNsigKProng2 = trackProng2.tofNSigmaKa();
      }

      if (fFillTPCnSigma) {
        xictreeObj.fPidTpcNsigPrProng0 = trackProng0.tpcNSigmaPr();
        xictreeObj.fPidTpcNsigPiProng0 = trackProng0.tpcNSigmaPi();
        xictreeObj.fPidTpcNsigKProng0 = trackProng0.tpcNSigmaKa();

        xictreeObj.fPidTpcNsigPrProng1 = trackProng1.tpcNSigmaPr();
        xictreeObj.fPidTpcNsigPiProng1 = trackProng1.tpcNSigmaPi();
        xictreeObj.fPidTpcNsigKProng1 = trackProng1.tpcNSigmaKa();

        xictreeObj.fPidTpcNsigPrProng2 = trackProng2.tpcNSigmaPr();
        xictreeObj.fPidTpcNsigPiProng2 = trackProng2.tpcNSigmaPi();
        xictreeObj.fPidTpcNsigKProng2 = trackProng2.tpcNSigmaKa();
      }

      // fill TTree
      if (fConfigWriteTTree) { // many parameters not set for photons: d1DCA,fd2DCA, fdectyp,fdau3pdg,fwMultpT,fwMultpT2,fwMultmT2
        tree->Fill();
      }
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
        registry.fill(HIST("MC/reconstructed/signal/hPtRecSig"), candidate.pt());  // rec. level pT

        if (mass != 0.) {
          registry.fill(HIST("MC/reconstructed/signal/hMassSig"), mass, candidate.pt());
        }
        if (massC != 0.) {
          registry.fill(HIST("MC/reconstructed/signal/hMassSig"), massC, candidate.pt());
        }

        registry.fill(HIST("MC/reconstructed/signal/hDecLengthRecSig"), candidate.decayLength(), candidate.pt());
        registry.fill(HIST("MC/reconstructed/signal/hPtProng0RecSig"), candidate.ptProng0(), candidate.pt());
        registry.fill(HIST("MC/reconstructed/signal/hPtProng1RecSig"), candidate.ptProng1(), candidate.pt());
        registry.fill(HIST("MC/reconstructed/signal/hPtProng2RecSig"), candidate.ptProng2(), candidate.pt());
        registry.fill(HIST("MC/reconstructed/signal/hd0Prong0RecSig"), candidate.impactParameter0(), candidate.pt());
        registry.fill(HIST("MC/reconstructed/signal/hd0Prong1RecSig"), candidate.impactParameter1(), candidate.pt());
        registry.fill(HIST("MC/reconstructed/signal/hd0Prong2RecSig"), candidate.impactParameter2(), candidate.pt());
        registry.fill(HIST("MC/reconstructed/signal/hCtRecSig"), ctXic(candidate), candidate.pt());
        registry.fill(HIST("MC/reconstructed/signal/hCPARecSig"), candidate.cpa(), candidate.pt());
        registry.fill(HIST("MC/reconstructed/signal/hEtaRecSig"), candidate.eta(), candidate.pt());
      } else {
        // Background
        registry.fill(HIST("MC/reconstructed/bakground/hPtRecBg"), candidate.pt());

        if (mass != 0.) {
          registry.fill(HIST("MC/reconstructed/bakground/hMassBg"), mass, candidate.pt());
        }
        if (massC != 0.) {
          registry.fill(HIST("MC/reconstructed/bakground/hMassBg"), massC, candidate.pt());
        }

        registry.fill(HIST("MC/reconstructed/bakground/hDecLengthRecBg"), candidate.decayLength(), candidate.pt());
        registry.fill(HIST("MC/reconstructed/bakground/hPtProng0RecBg"), candidate.ptProng0(), candidate.pt());
        registry.fill(HIST("MC/reconstructed/bakground/hPtProng1RecBg"), candidate.ptProng1(), candidate.pt());
        registry.fill(HIST("MC/reconstructed/bakground/hPtProng2RecBg"), candidate.ptProng2(), candidate.pt());
        registry.fill(HIST("MC/reconstructed/bakground/hd0Prong0RecBg"), candidate.impactParameter0(), candidate.pt());
        registry.fill(HIST("MC/reconstructed/bakground/hd0Prong1RecBg"), candidate.impactParameter1(), candidate.pt());
        registry.fill(HIST("MC/reconstructed/bakground/hd0Prong2RecBg"), candidate.impactParameter2(), candidate.pt());
        registry.fill(HIST("MC/reconstructed/bakground/hCtRecBg"), ctXic(candidate), candidate.pt());
        registry.fill(HIST("MC/reconstructed/bakground/hCPARecBg"), candidate.cpa(), candidate.pt());
        registry.fill(HIST("MC/reconstructed/bakground/hEtaRecBg"), candidate.eta(), candidate.pt());
      }
    }
    // MC gen.
    for (auto& particle : particlesMC) {
      if (std::abs(particle.flagMcMatchGen()) == 1 << DecayType::XicToPKPi) {
        if (yCandMax >= 0. && std::abs(RecoDecay::y(array{particle.px(), particle.py(), particle.pz()}, RecoDecay::getMassPDG(particle.pdgCode()))) > yCandMax) {
          continue;
        }
        registry.fill(HIST("MC/genearated/signal/hPtGen"), particle.pt());
        registry.fill(HIST("MC/genearated/signal/hEtaGen"), particle.eta(), particle.pt());
      }
    }
  }

  void SetTree()
  {
    tree.setObject(new TTree("XicTTree", "XicTTree"));

    tree->Branch("fPt", &xictreeObj.fPt, "fPt/F");
    tree->Branch("fMass", &xictreeObj.fMass, "fMass/F");
    tree->Branch("fDecayLength", &xictreeObj.fDecayLength, "fDecayLength/F");
    tree->Branch("fDecayLengthErr", &xictreeObj.fDecayLengthErr, "fDecayLengthErr/F");
    tree->Branch("fCpa", &xictreeObj.fCpa, "fCpa/F");
    tree->Branch("fCpaXY", &xictreeObj.fCpaXY, "fCpaXY/F");
    tree->Branch("fChi2PCA", &xictreeObj.fChi2PCA, "fChi2PCA/F");

    //---tpc nSigma for all prongs and Pi, K, Pr hypothesis
    tree->Branch("fPidTpcNsigPrProng0", &xictreeObj.fPidTpcNsigPrProng0, "fPidTpcNsigPrProng0/F");
    tree->Branch("fPidTpcNsigKProng0", &xictreeObj.fPidTpcNsigKProng0, "fPidTpcNsigKProng0/F");
    tree->Branch("fPidTpcNsigPiProng0", &xictreeObj.fPidTpcNsigPiProng0, "fPidTpcNsigPiProng0/F");

    tree->Branch("fPidTpcNsigPrProng1", &xictreeObj.fPidTpcNsigPrProng1, "fPidTpcNsigPrProng1/F");
    tree->Branch("fPidTpcNsigKProng1", &xictreeObj.fPidTpcNsigKProng1, "fPidTpcNsigKProng1/F");
    tree->Branch("fPidTpcNsigPiProng1", &xictreeObj.fPidTpcNsigPiProng1, "fPidTpcNsigPiProng1/F");

    tree->Branch("fPidTpcNsigPrProng2", &xictreeObj.fPidTpcNsigPrProng2, "fPidTpcNsigPrProng2/F");
    tree->Branch("fPidTpcNsigKProng2", &xictreeObj.fPidTpcNsigKProng2, "fPidTpcNsigKProng2/F");
    tree->Branch("fPidTpcNsigPiProng2", &xictreeObj.fPidTpcNsigPiProng2, "fPidTpcNsigPiProng2/F");

    //---tof nSigma for all prongs and Pi, K, Pr hypothesis
    tree->Branch("fPidTofNsigPrProng0", &xictreeObj.fPidTofNsigPrProng0, "fPidTofNsigPrProng0/F");
    tree->Branch("fPidTofNsigKProng0", &xictreeObj.fPidTofNsigKProng0, "fPidTofNsigKProng0/F");
    tree->Branch("fPidTofNsigPiProng0", &xictreeObj.fPidTofNsigPiProng0, "fPidTofNsigPiProng0/F");

    tree->Branch("fPidTofNsigPrProng1", &xictreeObj.fPidTofNsigPrProng1, "fPidTofNsigPrProng1/F");
    tree->Branch("fPidTofNsigKProng1", &xictreeObj.fPidTofNsigKProng1, "fPidTofNsigKProng1/F");
    tree->Branch("fPidTofNsigPiProng1", &xictreeObj.fPidTofNsigPiProng1, "fPidTofNsigPiProng1/F");

    tree->Branch("fPidTofNsigPrProng2", &xictreeObj.fPidTofNsigPrProng2, "fPidTofNsigPrProng2/F");
    tree->Branch("fPidTofNsigKProng2", &xictreeObj.fPidTofNsigKProng2, "fPidTofNsigKProng2/F");
    tree->Branch("fPidTofNsigPiProng2", &xictreeObj.fPidTofNsigPiProng2, "fPidTofNsigPiProng2/F");
  }

  PROCESS_SWITCH(HfTaskXic, processMc, "Process MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskXic>(cfgc)};
}
