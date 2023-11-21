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

/// \file taskD0.cxx
/// \brief D0 analysis task
///
/// \author Gian Michele Innocenti <gian.michele.innocenti@cern.ch>, CERN
/// \author Vít Kučera <vit.kucera@cern.ch>, CERN

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::framework::expressions;

/// D0 analysis task
struct HfTaskD0 {
  Configurable<int> selectionFlagD0{"selectionFlagD0", 1, "Selection Flag for D0"};
  Configurable<int> selectionFlagD0bar{"selectionFlagD0bar", 1, "Selection Flag for D0bar"};
  Configurable<double> yCandGenMax{"yCandGenMax", 0.5, "max. gen particle rapidity"};
  Configurable<double> yCandRecoMax{"yCandRecoMax", 0.8, "max. cand. rapidity"};
  Configurable<int> selectionFlagHf{"selectionFlagHf", 1, "Selection Flag for HF flagged candidates"};
  Configurable<int> selectionTopol{"selectionTopol", 1, "Selection Flag for topologically selected candidates"};
  Configurable<int> selectionCand{"selectionCand", 1, "Selection Flag for conj. topol. selected candidates"};
  Configurable<int> selectionPid{"selectionPid", 1, "Selection Flag for reco PID candidates"};
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_d0_to_pi_k::vecBinsPt}, "pT bin limits"};

  HfHelper hfHelper;

  using D0Candidates = soa::Join<aod::HfCand2Prong, aod::HfSelD0>;
  using D0CandidatesMc = soa::Join<D0Candidates, aod::HfCand2ProngMcRec>;
  using D0CandidatesKF = soa::Join<D0Candidates, aod::HfCand2ProngKF>;
  using D0CandidatesMcKF = soa::Join<D0CandidatesKF, aod::HfCand2ProngMcRec>;

  Partition<D0Candidates> selectedD0Candidates = aod::hf_sel_candidate_d0::isSelD0 >= selectionFlagD0 || aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlagD0bar;
  Partition<D0CandidatesKF> selectedD0CandidatesKF = aod::hf_sel_candidate_d0::isSelD0 >= selectionFlagD0 || aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlagD0bar;
  Partition<D0CandidatesMc> selectedD0CandidatesMc = aod::hf_sel_candidate_d0::isRecoHfFlag >= selectionFlagHf;
  Partition<D0CandidatesMcKF> selectedD0CandidatesMcKF = aod::hf_sel_candidate_d0::isRecoHfFlag >= selectionFlagHf;

  HistogramRegistry registry{
    "registry",
    {{"hPtCand", "2-prong candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"hPtProng0", "2-prong candidates;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"hPtProng1", "2-prong candidates;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"hPtRecSig", "2-prong candidates (matched);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"hPtRecSigPrompt", "2-prong candidates (matched, prompt);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"hPtRecSigNonPrompt", "2-prong candidates (matched, non-prompt);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"hPtRecBg", "2-prong candidates (unmatched);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"hPtGen", "MC particles (matched);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"hPtGenPrompt", "MC particles (matched, prompt);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"hPtGenNonPrompt", "MC particles (matched, non-prompt);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"hYGenPrompt", "MC particles (matched, prompt);#it{y}^{gen.};entries", {HistType::kTH1F, {{300, -1.5, 1.5}}}},
     {"hYGenNonPrompt", "MC particles (matched, non-prompt);#it{y}^{gen.};entries", {HistType::kTH1F, {{300, -1.5, 1.5}}}},
     {"hPtGenSig", "2-prong candidates (matched);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"hCPARecSig", "2-prong candidates (matched);cosine of pointing angle;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}}},
     {"hCPARecBg", "2-prong candidates (unmatched);cosine of pointing angle;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}}},
     {"hEtaRecSig", "2-prong candidates (matched);#it{#eta};entries", {HistType::kTH1F, {{100, -5., 5.}}}},
     {"hEtaRecBg", "2-prong candidates (unmatched);#it{#eta};entries", {HistType::kTH1F, {{100, -5., 5.}}}},
     {"hEtaGen", "MC particles (matched);#it{#eta};entries", {HistType::kTH1F, {{100, -5., 5.}}}},
     {"hPtGenVsPtRecSig", "2-prong candidates (matched);#it{p}_{T}^{gen.} (GeV/#it{c});#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH2F, {{360, 0., 36.}, {360, 0., 36.}}}},
     {"hYGenVsYRecSig", "2-prong candidates (matched);#it{y}^{gen.} ;#it{y}^{rec.} ;entries", {HistType::kTH2F, {{300, -1.5, 1.5}, {300, -1.5, 1.5}}}},
     {"hPtProng0Sig", "prong0 pt (matched); #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {10, -5., 5.}}}},
     {"hPtProng1Sig", "prong1 pt (matched); #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {10, -5., 5.}}}},
     {"hDecLengthSig", "2-prong candidates (matched);decay length (cm); #it{y}", {HistType::kTH2F, {{200, 0., 2.}, {10, -5., 5.}}}},
     {"hDecLengthXYSig", "2-prong candidates (matched);decay length xy (cm); #it{y}", {HistType::kTH2F, {{200, 0., 2.}, {10, -5., 5.}}}},
     {"hNormalisedDecLengthSig", "2-prong candidates (matched);normalised decay length (cm); #it{y}", {HistType::kTH2F, {{200, 0., 10.}, {10, -5., 5.}}}},
     {"hNormalisedDecLengthXYSig", "2-prong candidates (matched);normalised decay length xy (cm); #it{y}", {HistType::kTH2F, {{200, 0., 10.}, {10, -5., 5.}}}},
     {"hd0Prong0Sig", "2-prong candidates (matched);prong 0 DCAxy to prim. vertex (cm); #it{y}", {HistType::kTH2F, {{100, -1., 1.}, {10, -5., 5.}}}},
     {"hd0Prong1Sig", "2-prong candidates (matched);prong 1 DCAxy to prim. vertex (cm); #it{y}", {HistType::kTH2F, {{100, -1., 1.}, {10, -5., 5.}}}},
     {"hd0d0Sig", "2-prong candidates (matched);product of DCAxy to prim. vertex (cm^{2}); #it{y}", {HistType::kTH2F, {{500, -1., 1.}, {10, -5., 5.}}}},
     {"hCTSSig", "2-prong candidates (matched);cos #it{#theta}* (D^{0}); #it{y}", {HistType::kTH2F, {{110, -1.1, 1.1}, {10, -5., 5.}}}},
     {"hCtSig", "2-prong candidates (matched);proper lifetime (D^{0}) * #it{c} (cm); #it{y}", {HistType::kTH2F, {{120, -20., 100.}, {10, -5., 5.}}}},
     {"hCPASig", "2-prong candidates (matched);cosine of pointing angle; #it{y}", {HistType::kTH2F, {{440, -1.1, 1.1}, {10, -5., 5.}}}},
     {"hCPAxySig", "2-prong candidates (matched);cosine of pointing angle xy; #it{y}", {HistType::kTH2F, {{440, -1.1, 1.1}, {10, -5., 5.}}}},
     {"hPtProng0Bkg", "prong0 pt (matched); #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {10, -5., 5.}}}},
     {"hPtProng1Bkg", "prong1 pt (matched); #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {10, -5., 5.}}}},
     {"hDecLengthBkg", "2-prong candidates (checked);decay length (cm); #it{y}", {HistType::kTH2F, {{200, 0., 2.}, {10, -5., 5.}}}},
     {"hDecLengthXYBkg", "2-prong candidates (checked);decay length xy (cm); #it{y}", {HistType::kTH2F, {{200, 0., 2.}, {10, -5., 5.}}}},
     {"hNormalisedDecLengthBkg", "2-prong candidates (checked);normalised decay length (cm); #it{y}", {HistType::kTH2F, {{200, 0., 10.}, {10, -5., 5.}}}},
     {"hNormalisedDecLengthXYBkg", "2-prong candidates (checked);normalised decay length xy (cm); #it{y}", {HistType::kTH2F, {{200, 0., 10.}, {10, -5., 5.}}}},
     {"hd0Prong0Bkg", "2-prong candidates (checked);prong 0 DCAxy to prim. vertex (cm); #it{y}", {HistType::kTH2F, {{100, -1., 1.}, {10, -5., 5.}}}},
     {"hd0Prong1Bkg", "2-prong candidates (checked);prong 1 DCAxy to prim. vertex (cm); #it{y}", {HistType::kTH2F, {{100, -1., 1.}, {10, -5., 5.}}}},
     {"hd0d0Bkg", "2-prong candidates (checked);product of DCAxy to prim. vertex (cm^{2}); #it{y}", {HistType::kTH2F, {{500, -1., 1.}, {10, -5., 5.}}}},
     {"hCTSBkg", "2-prong candidates (checked);cos #it{#theta}* (D^{0}); #it{y}", {HistType::kTH2F, {{110, -1.1, 1.1}, {10, -5., 5.}}}},
     {"hCtBkg", "2-prong candidates (checked);proper lifetime (D^{0}) * #it{c} (cm); #it{y}", {HistType::kTH2F, {{120, -20., 100.}, {10, -5., 5.}}}},
     {"hCPABkg", "2-prong candidates (checked);cosine of pointing angle; #it{y}", {HistType::kTH2F, {{440, -1.1, 1.1}, {10, -5., 5.}}}},
     {"hCPAxyBkg", "2-prong candidates (checked);cosine of pointing angle xy; #it{y}", {HistType::kTH2F, {{440, -1.1, 1.1}, {10, -5., 5.}}}},
     {"hPtVsYRecSig_RecoPID", "2-prong candidates (RecoPID - matched);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hPtVsYRecSigPromptRecoPID", "2-prong candidates (RecoPID - matched, prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hPtVsYRecSigNonPromptRecoPID", "2-prong candidates (RecoPID - matched, non-prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hPtVsYRecSigRecoCand", "2-prong candidates (RecoCand - matched);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hPtVsYRecSigPromptRecoCand", "2-prong candidates (RecoCand - matched, prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hPtVsYRecSigNonPromptRecoCand", "2-prong candidates (RecoCand - matched, non-prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hPtVsYRecSigRecoTopol", "2-prong candidates (RecoTopol - matched);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hPtVsYRecSigPromptRecoTopol", "2-prong candidates (RecoTopol - matched, prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hPtVsYRecSigNonPromptRecoTopol", "2-prong candidates (RecoTopol - matched, non-prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hPtVsYRecSigRecoHFFlag", "2-prong candidates (RecoHFFlag - matched);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hPtVsYRecSigPromptRecoHFFlag", "2-prong candidates (RecoHFFlag - matched, prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hPtVsYRecSigNonPromptRecoHFFlag", "2-prong candidates (RecoHFFlag - matched, non-prompt);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hPtVsYRecSigReco", "2-prong candidates (Reco - matched);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hPtVsYRecSigPromptReco", "2-prong candidates (Reco - matched);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hPtVsYRecSigNonPromptReco", "2-prong candidates (Reco - matched);#it{p}_{T}^{rec.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hPtVsYGen", "2-prong candidates (matched);#it{p}_{T}^{gen.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hPtVsYGenPrompt", "2-prong candidates (matched, prompt);#it{p}_{T}^{gen.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hPtVsYGenNonPrompt", "2-prong candidates (matched, non-prompt);#it{p}_{T}^{gen.}; #it{y}", {HistType::kTH2F, {{360, 0., 36.}, {100, -5., 5.}}}},
     {"hMassVsPtGenVsPtRecSig", "2-prong candidates (matched);#it{m}_{inv} (GeV/#it{c}^{2});#it{p}_{T}^{gen.} (GeV/#it{c});#it{p}_{T}^{rec.} (GeV/#it{c})", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {150, 0., 30.}}}},
     {"hMassSigD0", "2-prong candidates (matched);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {20, -5., 5.}}}},
     {"hMassBkgD0", "2-prong candidates (checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {20, -5., 5.}}}},
     {"hMassReflBkgD0", "2-prong candidates (matched);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {20, -5., 5.}}}},
     {"hMassSigBkgD0", "2-prong candidates (not checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {20, -5., 5.}}}},
     {"hMassSigD0bar", "2-prong candidates (matched);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {20, -5., 5.}}}},
     {"hMassBkgD0bar", "2-prong candidates (checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {20, -5., 5.}}}},
     {"hMassReflBkgD0bar", "2-prong candidates (matched);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {20, -5., 5.}}}},
     {"hMassSigBkgD0bar", "2-prong candidates (not checked);#it{m}_{inv} (GeV/#it{c}^{2}); #it{p}_{T}; #it{y}", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {150, 0., 30.}, {20, -5., 5.}}}}}};

  void init(InitContext&)
  {
    std::array<bool, 4> doprocess{doprocessDataWithDCAFitterN, doprocessDataWithKFParticle, doprocessMcWithDCAFitterN, doprocessMcWithKFParticle};
    if ((std::accumulate(doprocess.begin(), doprocess.end(), 0)) == 0) {
      LOGP(fatal, "At least one process function should be enabled at a time.");
    }
    if ((doprocessDataWithDCAFitterN || doprocessMcWithDCAFitterN) && (doprocessDataWithKFParticle || doprocessMcWithKFParticle)) {
      LOGP(fatal, "DCAFitterN and KFParticle can not be enabled at a time.");
    }

    auto vbins = (std::vector<double>)binsPt;
    registry.add("hMass", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{500, 0., 5.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassVsPhi", "2-prong candidates vs phi;inv. mass (#pi K) (GeV/#it{c}^{2});phi (rad);entries", {HistType::kTH3F, {{120, 1.5848, 2.1848}, {vbins, "#it{p}_{T} (GeV/#it{c})"}, {32, 0, o2::constants::math::TwoPI}}});
    registry.add("hDecLength", "2-prong candidates;decay length (cm);entries", {HistType::kTH2F, {{800, 0., 4.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLengthxy", "2-prong candidates;decay length xy (cm);entries", {HistType::kTH2F, {{800, 0., 4.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLenErr", "2-prong candidates;decay length error (cm);entries", {HistType::kTH2F, {{800, 0., 0.2}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLenXYErr", "2-prong candidates;decay length xy error (cm);entries", {HistType::kTH2F, {{800, 0., 0.2}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hNormalisedDecLength", "2-prong candidates;decay length (cm);entries", {HistType::kTH2F, {{800, 0., 40.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hNormalisedDecLengthxy", "2-prong candidates;decay length xy (cm);entries", {HistType::kTH2F, {{800, 0., 40.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong0", "2-prong candidates;prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong1", "2-prong candidates;prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0ErrProng0", "2-prong candidates;prong 0 DCAxy to prim. vertex error (cm);entries", {HistType::kTH2F, {{800, 0., 0.2}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0ErrProng1", "2-prong candidates;prong 1 DCAxy to prim. vertex error (cm);entries", {HistType::kTH2F, {{800, 0., 0.2}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0d0", "2-prong candidates;product of DCAxy to prim. vertex (cm^{2});entries", {HistType::kTH2F, {{500, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCTS", "2-prong candidates;cos #it{#theta}* (D^{0});entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCt", "2-prong candidates;proper lifetime (D^{0}) * #it{c} (cm);entries", {HistType::kTH2F, {{120, -20., 100.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPA", "2-prong candidates;cosine of pointing angle;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hEta", "2-prong candidates;candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hSelectionStatus", "2-prong candidates;selection status;entries", {HistType::kTH2F, {{5, -0.5, 4.5}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassFinerBinning", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{120, 1.5848, 2.1848}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLengthFinerBinning", "2-prong candidates;decay length (cm);entries", {HistType::kTH2F, {{400, 0., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLengthxyFinerBinning", "2-prong candidates;decay length xy (cm);entries", {HistType::kTH2F, {{400, 0., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong0FinerBinning", "2-prong candidates;prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{500, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong1FinerBinning", "2-prong candidates;prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{500, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0d0FinerBinning", "2-prong candidates;product of DCAxy to prim. vertex (cm^{2});entries", {HistType::kTH2F, {{500, -0.1, 0.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCTSFinerBinning", "2-prong candidates;cos #it{#theta}* (D^{0});entries", {HistType::kTH2F, {{200, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCtFinerBinning", "2-prong candidates;proper lifetime (D^{0}) * #it{c} (cm);entries", {HistType::kTH2F, {{500, -0., 100.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPAFinerBinning", "2-prong candidates;cosine of pointing angle;entries", {HistType::kTH2F, {{200, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPAXYFinerBinning", "2-prong candidates;cosine of pointing angle xy;entries", {HistType::kTH2F, {{200, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong0VsPtSig", "2-prong candidates;prong 0 DCAxy to prim. vertex (cm) vs #it{p}_{T} for signal;entries", {HistType::kTH2F, {{500, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong1VsPtSig", "2-prong candidates;prong 1 DCAxy to prim. vertex (cm) vs #it{p}_{T} for signal;entries", {HistType::kTH2F, {{500, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0d0VsPtSig", "2-prong candidates;product of DCAxy to prim. vertex (cm^{2}) vs #it{p}_{T} for signal;entries", {HistType::kTH2F, {{500, -0.1, 0.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCTSVsPtSig", "2-prong candidates;cos #it{#theta}* (D^{0}) vs #it{p}_{T} for signal;entries", {HistType::kTH2F, {{200, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPAVsPtSig", "2-prong candidates;cosine of pointing angle vs #it{p}_{T} for signal;entries", {HistType::kTH2F, {{200, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPAXYVsPtSig", "2-prong candidates;cosine of pointing angle xy vs #it{p}_{T} for signal;entries", {HistType::kTH2F, {{200, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hNormalisedDecLengthxyVsPtSig", "2-prong candidates;decay length xy (cm) vs #it{p}_{T} for signal;entries", {HistType::kTH2F, {{800, 0., 40.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLengthVsPtSig", "2-prong candidates;decay length (cm) vs #it{p}_{T} for signal;entries", {HistType::kTH2F, {{800, 0., 4.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLengthxyVsPtSig", "2-prong candidates;decay length xy (cm) vs #it{p}_{T} for signal;entries", {HistType::kTH2F, {{800, 0., 4.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
  }

  template <int reconstructionType, typename CandType>
  void processData(CandType const& candidates)
  {
    for (const auto& candidate : candidates) {
      if (!(candidate.hfflag() & 1 << aod::hf_cand_2prong::DecayType::D0ToPiK)) {
        continue;
      }
      if (yCandRecoMax >= 0. && std::abs(hfHelper.yD0(candidate)) > yCandRecoMax) {
        continue;
      }

      float massD0, massD0bar;
      if constexpr (reconstructionType == aod::hf_cand::VertexerType::KfParticle) {
        massD0 = candidate.kfGeoMassD0();
        massD0bar = candidate.kfGeoMassD0bar();
      } else {
        massD0 = hfHelper.invMassD0ToPiK(candidate);
        massD0bar = hfHelper.invMassD0barToKPi(candidate);
      }
      auto ptCandidate = candidate.pt();

      if (candidate.isSelD0() >= selectionFlagD0) {
        registry.fill(HIST("hMass"), massD0, ptCandidate);
        registry.fill(HIST("hMassFinerBinning"), massD0, ptCandidate);
        registry.fill(HIST("hMassVsPhi"), massD0, ptCandidate, candidate.phi());
      }
      if (candidate.isSelD0bar() >= selectionFlagD0bar) {
        registry.fill(HIST("hMass"), massD0bar, ptCandidate);
        registry.fill(HIST("hMassFinerBinning"), massD0bar, ptCandidate);
        registry.fill(HIST("hMassVsPhi"), massD0bar, ptCandidate, candidate.phi());
      }
      registry.fill(HIST("hPtCand"), ptCandidate);
      registry.fill(HIST("hPtProng0"), candidate.ptProng0());
      registry.fill(HIST("hPtProng1"), candidate.ptProng1());
      registry.fill(HIST("hDecLength"), candidate.decayLength(), ptCandidate);
      registry.fill(HIST("hDecLengthxy"), candidate.decayLengthXY(), ptCandidate);
      registry.fill(HIST("hDecLenErr"), candidate.errorDecayLength(), ptCandidate);
      registry.fill(HIST("hDecLenXYErr"), candidate.errorDecayLengthXY(), ptCandidate);
      registry.fill(HIST("hNormalisedDecLength"), candidate.decayLengthNormalised(), ptCandidate);
      registry.fill(HIST("hNormalisedDecLengthxy"), candidate.decayLengthXYNormalised(), ptCandidate);
      registry.fill(HIST("hd0Prong0"), candidate.impactParameter0(), ptCandidate);
      registry.fill(HIST("hd0Prong1"), candidate.impactParameter1(), ptCandidate);
      registry.fill(HIST("hd0ErrProng0"), candidate.errorImpactParameter0(), ptCandidate);
      registry.fill(HIST("hd0ErrProng1"), candidate.errorImpactParameter1(), ptCandidate);
      registry.fill(HIST("hd0d0"), candidate.impactParameterProduct(), ptCandidate);
      registry.fill(HIST("hCTS"), hfHelper.cosThetaStarD0(candidate), ptCandidate);
      registry.fill(HIST("hCt"), hfHelper.ctD0(candidate), ptCandidate);
      registry.fill(HIST("hCPA"), candidate.cpa(), ptCandidate);
      registry.fill(HIST("hEta"), candidate.eta(), ptCandidate);
      registry.fill(HIST("hSelectionStatus"), candidate.isSelD0() + (candidate.isSelD0bar() * 2), ptCandidate);
      registry.fill(HIST("hDecLengthFinerBinning"), candidate.decayLength(), ptCandidate);
      registry.fill(HIST("hDecLengthxyFinerBinning"), candidate.decayLengthXY(), ptCandidate);
      registry.fill(HIST("hd0Prong0FinerBinning"), candidate.impactParameter0(), ptCandidate);
      registry.fill(HIST("hd0Prong1FinerBinning"), candidate.impactParameter1(), ptCandidate);
      registry.fill(HIST("hd0d0FinerBinning"), candidate.impactParameterProduct(), ptCandidate);
      registry.fill(HIST("hCTSFinerBinning"), hfHelper.cosThetaStarD0(candidate), ptCandidate);
      registry.fill(HIST("hCtFinerBinning"), hfHelper.ctD0(candidate), ptCandidate);
      registry.fill(HIST("hCPAFinerBinning"), candidate.cpa(), ptCandidate);
      registry.fill(HIST("hCPAXYFinerBinning"), candidate.cpaXY(), ptCandidate);
    }
  }
  void processDataWithDCAFitterN(D0Candidates const&)
  {
    processData<aod::hf_cand::VertexerType::DCAFitter>(selectedD0Candidates);
  }
  PROCESS_SWITCH(HfTaskD0, processDataWithDCAFitterN, "process taskD0 with DCAFitterN", true);

  void processDataWithKFParticle(D0CandidatesKF const&)
  {
    processData<aod::hf_cand::VertexerType::KfParticle>(selectedD0CandidatesKF);
  }
  PROCESS_SWITCH(HfTaskD0, processDataWithKFParticle, "process taskD0 with KFParticle", false);

  template <int reconstructionType, typename CandType>
  void processMc(CandType const& candidates,
                 soa::Join<aod::McParticles, aod::HfCand2ProngMcGen> const& mcParticles,
                 aod::TracksWMc const& tracks)
  {
    // MC rec.
    for (const auto& candidate : candidates) {
      if (!(candidate.hfflag() & 1 << aod::hf_cand_2prong::DecayType::D0ToPiK)) {
        continue;
      }
      if (yCandRecoMax >= 0. && std::abs(hfHelper.yD0(candidate)) > yCandRecoMax) {
        continue;
      }
      float massD0, massD0bar;
      if constexpr (reconstructionType == aod::hf_cand::VertexerType::KfParticle) {
        massD0 = candidate.kfGeoMassD0();
        massD0bar = candidate.kfGeoMassD0bar();
      } else {
        massD0 = hfHelper.invMassD0ToPiK(candidate);
        massD0bar = hfHelper.invMassD0barToKPi(candidate);
      }
      if (std::abs(candidate.flagMcMatchRec()) == 1 << aod::hf_cand_2prong::DecayType::D0ToPiK) {
        // Get the corresponding MC particle.
        auto indexMother = RecoDecay::getMother(mcParticles, candidate.template prong0_as<aod::TracksWMc>().template mcParticle_as<soa::Join<aod::McParticles, aod::HfCand2ProngMcGen>>(), pdg::Code::kD0, true);
        auto particleMother = mcParticles.rawIteratorAt(indexMother);
        auto ptGen = particleMother.pt();                                                                                                                     // gen. level pT
        auto yGen = RecoDecay::y(std::array{particleMother.px(), particleMother.py(), particleMother.pz()}, o2::analysis::pdg::MassD0);                       // gen. level y
        registry.fill(HIST("hPtGenSig"), ptGen);                                                                                                              // gen. level pT
        auto ptRec = candidate.pt();
        auto yRec = hfHelper.yD0(candidate);
        if (candidate.isRecoHfFlag() >= selectionFlagHf) {
          registry.fill(HIST("hPtVsYRecSigRecoHFFlag"), ptRec, yRec);
          registry.fill(HIST("hPtGenVsPtRecSig"), ptGen, ptRec);
          registry.fill(HIST("hYGenVsYRecSig"), yGen, yRec);
          if (candidate.isSelD0() >= selectionFlagD0) {
            registry.fill(HIST("hMassVsPtGenVsPtRecSig"), massD0, ptGen, ptRec);
          }
          if (candidate.isSelD0bar() >= selectionFlagD0bar) {
            registry.fill(HIST("hMassVsPtGenVsPtRecSig"), massD0bar, ptGen, ptRec);
          }
        }
        if (candidate.isRecoTopol() >= selectionTopol) {
          registry.fill(HIST("hPtVsYRecSigRecoTopol"), ptRec, yRec);
        }
        if (candidate.isRecoCand() >= selectionCand) {
          registry.fill(HIST("hPtVsYRecSigRecoCand"), ptRec, yRec);
        }
        if (candidate.isRecoPid() >= selectionPid) {
          registry.fill(HIST("hPtVsYRecSig_RecoPID"), ptRec, yRec);
        }
        if (candidate.isSelD0() >= selectionFlagD0 || candidate.isSelD0bar() >= selectionFlagD0bar) {
          registry.fill(HIST("hPtVsYRecSigReco"), ptRec, yRec); // rec. level pT
          registry.fill(HIST("hPtRecSig"), ptRec);
        }
        if (candidate.originMcRec() == RecoDecay::OriginType::Prompt) {
          if (candidate.isRecoHfFlag() >= selectionFlagHf) {
            registry.fill(HIST("hPtVsYRecSigPromptRecoHFFlag"), ptRec, yRec);
          }
          if (candidate.isRecoTopol() >= selectionTopol) {
            registry.fill(HIST("hPtVsYRecSigPromptRecoTopol"), ptRec, yRec);
          }
          if (candidate.isRecoCand() >= selectionCand) {
            registry.fill(HIST("hPtVsYRecSigPromptRecoCand"), ptRec, yRec);
          }
          if (candidate.isRecoPid() >= selectionPid) {
            registry.fill(HIST("hPtVsYRecSigPromptRecoPID"), ptRec, yRec);
          }
          if (candidate.isSelD0() >= selectionFlagD0 || candidate.isSelD0bar() >= selectionFlagD0bar) {
            registry.fill(HIST("hPtVsYRecSigPromptReco"), ptRec, yRec); // rec. level pT, prompt
            registry.fill(HIST("hPtRecSigPrompt"), ptRec);
          }
        } else {
          if (candidate.isRecoHfFlag() >= selectionFlagHf) {
            registry.fill(HIST("hPtVsYRecSigNonPromptRecoHFFlag"), ptRec, yRec);
          }
          if (candidate.isRecoTopol() >= selectionTopol) {
            registry.fill(HIST("hPtVsYRecSigNonPromptRecoTopol"), ptRec, yRec);
          }
          if (candidate.isRecoCand() >= selectionCand) {
            registry.fill(HIST("hPtVsYRecSigNonPromptRecoCand"), ptRec, yRec);
          }
          if (candidate.isRecoPid() >= selectionPid) {
            registry.fill(HIST("hPtVsYRecSigNonPromptRecoPID"), ptRec, yRec);
          }
          if (candidate.isSelD0() >= selectionFlagD0 || candidate.isSelD0bar() >= selectionFlagD0bar) {
            registry.fill(HIST("hPtVsYRecSigNonPromptReco"), ptRec, yRec); // rec. level pT, non-prompt
            registry.fill(HIST("hPtRecSigNonPrompt"), ptRec);
          }
        }
        registry.fill(HIST("hCPARecSig"), candidate.cpa());
        registry.fill(HIST("hEtaRecSig"), candidate.eta());
      } else {
        registry.fill(HIST("hPtRecBg"), candidate.pt());
        registry.fill(HIST("hCPARecBg"), candidate.cpa());
        registry.fill(HIST("hEtaRecBg"), candidate.eta());
      }
      auto ptCandidate = candidate.pt();
      auto ptProng0 = candidate.ptProng0();
      auto ptProng1 = candidate.ptProng1();
      auto rapidityCandidate = hfHelper.yD0(candidate);
      auto declengthCandidate = candidate.decayLength();
      auto declengthxyCandidate = candidate.decayLengthXY();
      auto normaliseddeclengthCandidate = candidate.decayLengthNormalised();
      auto normaliseddeclengthxyCandidate = candidate.decayLengthXYNormalised();
      auto d0Prong0 = candidate.impactParameter0();
      auto d0Prong1 = candidate.impactParameter1();
      auto d0d0Candidate = candidate.impactParameterProduct();
      auto ctsCandidate = hfHelper.cosThetaStarD0(candidate);
      auto ctCandidate = hfHelper.ctD0(candidate);
      auto cpaCandidate = candidate.cpa();
      auto cpaxyCandidate = candidate.cpaXY();
      if (candidate.isSelD0() >= selectionFlagD0) {
        registry.fill(HIST("hMassSigBkgD0"), massD0, ptCandidate, rapidityCandidate);
        if (candidate.flagMcMatchRec() == (1 << aod::hf_cand_2prong::DecayType::D0ToPiK)) {
          registry.fill(HIST("hPtProng0Sig"), ptProng0, rapidityCandidate);
          registry.fill(HIST("hPtProng1Sig"), ptProng1, rapidityCandidate);
          registry.fill(HIST("hDecLengthSig"), declengthCandidate, rapidityCandidate);
          registry.fill(HIST("hDecLengthXYSig"), declengthxyCandidate, rapidityCandidate);
          registry.fill(HIST("hNormalisedDecLengthSig"), normaliseddeclengthCandidate, rapidityCandidate);
          registry.fill(HIST("hNormalisedDecLengthXYSig"), normaliseddeclengthxyCandidate, rapidityCandidate);
          registry.fill(HIST("hd0Prong0Sig"), d0Prong0, rapidityCandidate);
          registry.fill(HIST("hd0Prong1Sig"), d0Prong1, rapidityCandidate);
          registry.fill(HIST("hd0d0Sig"), d0d0Candidate, rapidityCandidate);
          registry.fill(HIST("hCTSSig"), ctsCandidate, rapidityCandidate);
          registry.fill(HIST("hCtSig"), ctCandidate, rapidityCandidate);
          registry.fill(HIST("hCPASig"), cpaCandidate, rapidityCandidate);
          registry.fill(HIST("hCPAxySig"), cpaxyCandidate, rapidityCandidate);
          registry.fill(HIST("hd0Prong0VsPtSig"), d0Prong0, ptCandidate);
          registry.fill(HIST("hd0Prong1VsPtSig"), d0Prong1, ptCandidate);
          registry.fill(HIST("hd0d0VsPtSig"), d0d0Candidate, ptCandidate);
          registry.fill(HIST("hCTSVsPtSig"), ctsCandidate, ptCandidate);
          registry.fill(HIST("hCPAVsPtSig"), cpaCandidate, ptCandidate);
          registry.fill(HIST("hCPAXYVsPtSig"), cpaxyCandidate, ptCandidate);
          registry.fill(HIST("hNormalisedDecLengthxyVsPtSig"), normaliseddeclengthxyCandidate, ptCandidate);
          registry.fill(HIST("hDecLengthVsPtSig"), declengthCandidate, ptCandidate);
          registry.fill(HIST("hDecLengthxyVsPtSig"), declengthxyCandidate, ptCandidate);
          registry.fill(HIST("hMassSigD0"), massD0, ptCandidate, rapidityCandidate);
        } else {
          registry.fill(HIST("hPtProng0Bkg"), ptProng0, rapidityCandidate);
          registry.fill(HIST("hPtProng1Bkg"), ptProng1, rapidityCandidate);
          registry.fill(HIST("hDecLengthBkg"), declengthCandidate, rapidityCandidate);
          registry.fill(HIST("hDecLengthXYBkg"), declengthxyCandidate, rapidityCandidate);
          registry.fill(HIST("hNormalisedDecLengthBkg"), normaliseddeclengthCandidate, rapidityCandidate);
          registry.fill(HIST("hNormalisedDecLengthXYBkg"), normaliseddeclengthxyCandidate, rapidityCandidate);
          registry.fill(HIST("hd0Prong0Bkg"), d0Prong0, rapidityCandidate);
          registry.fill(HIST("hd0Prong1Bkg"), d0Prong1, rapidityCandidate);
          registry.fill(HIST("hd0d0Bkg"), d0d0Candidate, rapidityCandidate);
          registry.fill(HIST("hCTSBkg"), ctsCandidate, rapidityCandidate);
          registry.fill(HIST("hCtBkg"), ctCandidate, rapidityCandidate);
          registry.fill(HIST("hCPABkg"), cpaCandidate, rapidityCandidate);
          registry.fill(HIST("hCPAxyBkg"), cpaxyCandidate, rapidityCandidate);
          registry.fill(HIST("hMassBkgD0"), massD0, ptCandidate, rapidityCandidate);
          if (candidate.flagMcMatchRec() == -(1 << aod::hf_cand_2prong::DecayType::D0ToPiK)) {
            registry.fill(HIST("hMassReflBkgD0"), massD0, ptCandidate, rapidityCandidate);
          }
        }
      }
      if (candidate.isSelD0bar() >= selectionFlagD0bar) {
        registry.fill(HIST("hMassSigBkgD0bar"), massD0bar, ptCandidate, rapidityCandidate);
        if (candidate.flagMcMatchRec() == -(1 << aod::hf_cand_2prong::DecayType::D0ToPiK)) {
          registry.fill(HIST("hMassSigD0bar"), massD0bar, ptCandidate, rapidityCandidate);
        } else {
          registry.fill(HIST("hMassBkgD0bar"), massD0bar, ptCandidate, rapidityCandidate);
          if (candidate.flagMcMatchRec() == (1 << aod::hf_cand_2prong::DecayType::D0ToPiK)) {
            registry.fill(HIST("hMassReflBkgD0bar"), massD0bar, ptCandidate, rapidityCandidate);
          }
        }
      }
    }
    // MC gen.
    for (const auto& particle : mcParticles) {
      if (std::abs(particle.flagMcMatchGen()) == 1 << aod::hf_cand_2prong::DecayType::D0ToPiK) {
        if (yCandGenMax >= 0. && std::abs(RecoDecay::y(std::array{particle.px(), particle.py(), particle.pz()}, o2::analysis::pdg::MassD0)) > yCandGenMax) {
          continue;
        }
        auto ptGen = particle.pt();
        auto yGen = RecoDecay::y(std::array{particle.px(), particle.py(), particle.pz()}, o2::analysis::pdg::MassD0);
        registry.fill(HIST("hPtGen"), ptGen);
        registry.fill(HIST("hPtVsYGen"), ptGen, yGen);
        if (particle.originMcGen() == RecoDecay::OriginType::Prompt) {
          registry.fill(HIST("hPtGenPrompt"), ptGen);
          registry.fill(HIST("hYGenPrompt"), yGen);
          registry.fill(HIST("hPtVsYGenPrompt"), ptGen, yGen);
        } else {
          registry.fill(HIST("hPtGenNonPrompt"), ptGen);
          registry.fill(HIST("hYGenNonPrompt"), yGen);
          registry.fill(HIST("hPtVsYGenNonPrompt"), ptGen, yGen);
        }
        registry.fill(HIST("hEtaGen"), particle.eta());
      }
    }
  }

  void processMcWithDCAFitterN(D0CandidatesMc const&,
                               soa::Join<aod::McParticles, aod::HfCand2ProngMcGen> const& mcParticles,
                               aod::TracksWMc const& tracks)
  {
    processMc<aod::hf_cand::VertexerType::DCAFitter>(selectedD0CandidatesMc, mcParticles, tracks);
  }
  PROCESS_SWITCH(HfTaskD0, processMcWithDCAFitterN, "Process MC with DCAFitterN", false);

  void processMcWithKFParticle(D0CandidatesMcKF const&,
                               soa::Join<aod::McParticles, aod::HfCand2ProngMcGen> const& mcParticles,
                               aod::TracksWMc const& tracks)
  {
    processMc<aod::hf_cand::VertexerType::KfParticle>(selectedD0CandidatesMcKF, mcParticles, tracks);
  }
  PROCESS_SWITCH(HfTaskD0, processMcWithKFParticle, "Process MC with KFParticle", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskD0>(cfgc)};
}
