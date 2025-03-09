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

/// \file taskLc.cxx
/// \brief Λc± → p± K∓ π± analysis task
/// \note Extended from taskD0
///
/// \author Luigi Dello Stritto <luigi.dello.stritto@cern.ch>, University and INFN SALERNO
/// \author Vít Kučera <vit.kucera@cern.ch>, CERN
/// \author Annalena Kalteyer <annalena.sophie.kalteyer@cern.ch>, GSI Darmstadt
/// \author Biao Zhang <biao.zhang@cern.ch>, Heidelberg University

#include <vector> // std::vector

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/CentralityEstimation.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/Utils/utilsEvSelHf.h"

using namespace o2;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::hf_centrality;
using namespace o2::hf_occupancy;

/// Λc± → p± K∓ π± analysis task
struct HfTaskLc {
  Configurable<int> selectionFlagLc{"selectionFlagLc", 1, "Selection Flag for Lc"};
  Configurable<double> yCandGenMax{"yCandGenMax", 0.5, "max. gen particle rapidity"};
  Configurable<double> yCandRecoMax{"yCandRecoMax", 0.8, "max. cand. rapidity"};
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_lc_to_p_k_pi::vecBinsPt}, "pT bin limits"};
  // ThnSparse for ML outputScores and Vars
  Configurable<bool> fillTHn{"fillTHn", false, "fill THn"};
  Configurable<bool> storeOccupancy{"storeOccupancy", true, "Flag to store occupancy information"};
  Configurable<int> occEstimator{"occEstimator", 2, "Occupancy estimation (None: 0, ITS: 1, FT0C: 2)"};

  HfHelper hfHelper;
  SliceCache cache;

  using Collisions = soa::Join<aod::Collisions, aod::EvSels>;
  using CollisionsMc = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels>;
  using CollisionsWithFT0C = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs>;
  using CollisionsMcWithFT0C = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::CentFT0Cs>;
  using CollisionsWithFT0M = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms>;
  using CollisionsMcWithFT0M = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::CentFT0Ms>;

  using LcCandidates = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelLc>>;
  using LcCandidatesMl = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelLc, aod::HfMlLcToPKPi>>;

  using LcCandidatesMc = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelLc, aod::HfCand3ProngMcRec>>;
  using LcCandidatesMlMc = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelLc, aod::HfMlLcToPKPi, aod::HfCand3ProngMcRec>>;
  using McParticles3ProngMatched = soa::Join<aod::McParticles, aod::HfCand3ProngMcGen>;
  Filter filterSelectCandidates = aod::hf_sel_candidate_lc::isSelLcToPKPi >= selectionFlagLc || aod::hf_sel_candidate_lc::isSelLcToPiKP >= selectionFlagLc;
  Preslice<aod::HfCand3Prong> candLcPerCollision = aod::hf_cand::collisionId;
  PresliceUnsorted<aod::McCollisionLabels> colPerMcCollision = aod::mcparticle::mcCollisionId;

  ConfigurableAxis thnConfigAxisPt{"thnConfigAxisPt", {72, 0, 36}, ""};
  ConfigurableAxis thnConfigAxisMass{"thnConfigAxisMass", {300, 1.98, 2.58}, ""};
  ConfigurableAxis thnConfigAxisPtProng{"thnConfigAxisPtProng", {100, 0, 20}, ""};
  ConfigurableAxis thnConfigAxisCentrality{"thnConfigAxisCentrality", {100, 0, 100}, ""};
  ConfigurableAxis thnConfigAxisChi2PCA{"thnConfigAxisChi2PCA", {100, 0, 20}, ""};
  ConfigurableAxis thnConfigAxisDecLength{"thnConfigAxisDecLength", {10, 0, 0.05}, ""};
  ConfigurableAxis thnConfigAxisCPA{"thnConfigAxisCPA", {20, 0.8, 1}, ""};
  ConfigurableAxis thnConfigAxisBdtScoreBkg{"thnConfigAxisBdtScoreBkg", {1000, 0., 1.}, ""};
  ConfigurableAxis thnConfigAxisBdtScoreSignal{"thnConfigAxisBdtScoreSignal", {100, 0., 1.}, ""};
  ConfigurableAxis thnConfigAxisCanType{"thnConfigAxisCanType", {5, 0., 5.}, ""};
  ConfigurableAxis thnAxisRapidity{"thnAxisRapidity", {20, -1, 1}, "Cand. rapidity bins"};
  ConfigurableAxis thnConfigAxisGenPtB{"thnConfigAxisGenPtB", {1000, 0, 100}, "Gen Pt B"};
  ConfigurableAxis thnConfigAxisNumPvContr{"thnConfigAxisNumPvContr", {200, -0.5, 199.5}, "Number of PV contributors"};
  ConfigurableAxis thnConfigAxisOccupancy{"thnConfigAxisOccupancy", {14, 0, 14000}, "axis for centrality"};

  HistogramRegistry registry{
    "registry",
    {/// mass candidate
     {"Data/hMass", "3-prong candidates;inv. mass (p K #pi) (GeV/#it{c}^{2})", {HistType::kTH1F, {{600, 1.98, 2.58}}}},
     {"MC/reconstructed/signal/hMassRecSig", "3-prong candidates (matched);inv. mass (p K #pi) (GeV/#it{c}^{2})", {HistType::kTH1F, {{600, 1.98, 2.58}}}},
     {"MC/reconstructed/prompt/hMassRecSigPrompt", "3-prong candidates (matched, prompt);inv. mass (p K #pi) (GeV/#it{c}^{2})", {HistType::kTH1F, {{600, 1.98, 2.58}}}},
     {"MC/reconstructed/nonprompt/hMassRecSigNonPrompt", "3-prong candidates (matched, non-prompt);inv. mass (p K #pi) (GeV/#it{c}^{2})", {HistType::kTH1F, {{600, 1.98, 2.58}}}},
     /// pT
     {"Data/hPt", "3-prong candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"MC/reconstructed/signal/hPtRecSig", "3-prong candidates (matched);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"MC/reconstructed/prompt/hPtRecSigPrompt", "3-prong candidates (matched, prompt);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"MC/reconstructed/nonprompt/hPtRecSigNonPrompt", "3-prong candidates (matched, non-prompt);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"MC/generated/signal/hPtGen", "MC particles (matched);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"MC/generated/prompt/hPtGenPrompt", "MC particles (matched, prompt);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"MC/generated/nonprompt/hPtGenNonPrompt", "MC particles (matched, non-prompt);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"MC/generated/signal/hPtGenSig", "3-prong candidates (matched);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"Data/hPtProng0", "3-prong candidates;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"MC/reconstructed/signal/hPtRecProng0Sig", "3-prong candidates (matched);prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"MC/reconstructed/prompt/hPtRecProng0SigPrompt", "3-prong candidates (matched, prompt);prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"MC/reconstructed/nonprompt/hPtRecProng0SigNonPrompt", "3-prong candidates (matched, non-prompt);prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"Data/hPtProng1", "3-prong candidates;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"MC/reconstructed/signal/hPtRecProng1Sig", "3-prong candidates (matched);prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"MC/reconstructed/prompt/hPtRecProng1SigPrompt", "3-prong candidates (matched, prompt);prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"MC/reconstructed/nonprompt/hPtRecProng1SigNonPrompt", "3-prong candidates (matched, non-prompt);prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"Data/hPtProng2", "3-prong candidates;prong 2 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"MC/reconstructed/signal/hPtRecProng2Sig", "3-prong candidates (matched);prong 2 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"MC/reconstructed/prompt/hPtRecProng2SigPrompt", "3-prong candidates (matched, prompt);prong 2 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     {"MC/reconstructed/nonprompt/hPtRecProng2SigNonPrompt", "3-prong candidates (matched, non-prompt);prong 2 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}},
     /// DCAxy to prim. vertex prongs
     {"Data/hd0Prong0", "3-prong candidates;prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{600, -0.4, 0.4}}}},
     {"MC/reconstructed/signal/hd0RecProng0Sig", "3-prong candidates (matched);prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{600, -0.4, 0.4}}}},
     {"MC/reconstructed/prompt/hd0RecProng0SigPrompt", "3-prong candidates (matched, prompt);prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{600, -0.4, 0.4}}}},
     {"MC/reconstructed/nonprompt/hd0RecProng0SigNonPrompt", "3-prong candidates (matched, non-prompt);prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{600, -0.4, 0.4}}}},
     {"Data/hd0Prong1", "3-prong candidates;prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{600, -0.4, 0.4}}}},
     {"MC/reconstructed/signal/hd0RecProng1Sig", "3-prong candidates (matched);prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{600, -0.4, 0.4}}}},
     {"MC/reconstructed/prompt/hd0RecProng1SigPrompt", "3-prong candidates (matched, prompt);prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{600, -0.4, 0.4}}}},
     {"MC/reconstructed/nonprompt/hd0RecProng1SigNonPrompt", "3-prong candidates (matched, non-prompt);prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{600, -0.4, 0.4}}}},
     {"Data/hd0Prong2", "3-prong candidates;prong 2 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{600, -0.4, 0.4}}}},
     {"MC/reconstructed/signal/hd0RecProng2Sig", "3-prong candidates (matched);prong 2 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{600, -0.4, 0.4}}}},
     {"MC/reconstructed/prompt/hd0RecProng2SigPrompt", "3-prong candidates (matched, prompt);prong 2 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{600, -0.4, 0.4}}}},
     {"MC/reconstructed/nonprompt/hd0RecProng2SigNonPrompt", "3-prong candidates (matched, non-prompt);prong 2 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {{600, -0.4, 0.4}}}},
     /// decay length candidate
     {"Data/hDecLength", "3-prong candidates;decay length (cm);entries", {HistType::kTH1F, {{400, 0., 1.}}}},
     {"MC/reconstructed/signal/hDecLengthRecSig", "3-prong candidates (matched);decay length (cm);entries", {HistType::kTH1F, {{400, 0., 1.}}}},
     {"MC/reconstructed/prompt/hDecLengthRecSigPrompt", "3-prong candidates (matched, prompt);decay length (cm);entries", {HistType::kTH1F, {{400, 0., 1.}}}},
     {"MC/reconstructed/nonprompt/hDecLengthRecSigNonPrompt", "3-prong candidates (matched, non-prompt);decay length (cm);entries", {HistType::kTH1F, {{400, 0., 1.}}}},
     /// decay length xy candidate
     {"Data/hDecLengthxy", "3-prong candidates;decay length xy (cm);entries", {HistType::kTH1F, {{400, 0., 1.}}}},
     {"MC/reconstructed/signal/hDecLengthxyRecSig", "3-prong candidates (matched);decay length xy (cm);entries", {HistType::kTH1F, {{400, 0., 1.}}}},
     {"MC/reconstructed/prompt/hDecLengthxyRecSigPrompt", "3-prong candidates (matched, prompt);decay length xy (cm);entries", {HistType::kTH1F, {{400, 0., 1.}}}},
     {"MC/reconstructed/nonprompt/hDecLengthxyRecSigNonPrompt", "3-prong candidates (matched, non-prompt);decay length xy (cm);entries", {HistType::kTH1F, {{400, 0., 1.}}}},
     /// proper lifetime
     {"Data/hCt", "3-prong candidates;proper lifetime (#Lambda_{c}) * #it{c} (cm);entries", {HistType::kTH1F, {{100, 0., 0.2}}}},
     {"MC/reconstructed/signal/hCtRecSig", "3-prong candidates (matched);proper lifetime (#Lambda_{c}) * #it{c} (cm);entries", {HistType::kTH1F, {{100, 0., 0.2}}}},
     {"MC/reconstructed/prompt/hCtRecSigPrompt", "3-prong candidates (matched, prompt);proper lifetime (#Lambda_{c}) * #it{c} (cm);entries", {HistType::kTH1F, {{100, 0., 0.2}}}},
     {"MC/reconstructed/nonprompt/hCtRecSigNonPrompt", "3-prong candidates (matched, non-prompt);proper lifetime (#Lambda_{c}) * #it{c} (cm);entries", {HistType::kTH1F, {{100, 0., 0.2}}}},
     /// cosine of pointing angle
     {"Data/hCPA", "3-prong candidates;cosine of pointing angle;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}}},
     {"MC/reconstructed/signal/hCPARecSig", "3-prong candidates;cosine of pointing angle;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}}},
     {"MC/reconstructed/prompt/hCPARecSigPrompt", "3-prong candidates;cosine of pointing angle;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}}},
     {"MC/reconstructed/nonprompt/hCPARecSigNonPrompt", "3-prong candidates;cosine of pointing angle;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}}},
     /// cosine of pointing angle xy
     {"Data/hCPAxy", "3-prong candidates;cosine of pointing angle xy;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}}},
     {"MC/reconstructed/signal/hCPAxyRecSig", "3-prong candidates;cosine of pointing angle xy;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}}},
     {"MC/reconstructed/prompt/hCPAxyRecSigPrompt", "3-prong candidates;cosine of pointing angle xy;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}}},
     {"MC/reconstructed/nonprompt/hCPAxyRecSigNonPrompt", "3-prong candidates;cosine of pointing angle xy;entries", {HistType::kTH1F, {{110, -1.1, 1.1}}}},
     /// Chi 2 PCA to sec. vertex
     {"Data/hDca2", "3-prong candidates;prong Chi2PCA to sec. vertex (cm);entries", {HistType::kTH1F, {{400, 0., 20.}}}},
     {"MC/reconstructed/signal/hDca2RecSig", "3-prong candidates (matched);prong Chi2PCA to sec. vertex (cm);entries", {HistType::kTH1F, {{400, 0., 20.}}}},
     {"MC/reconstructed/prompt/hDca2RecSigPrompt", "3-prong candidates (matched);prong Chi2PCA to sec. vertex (cm);entries", {HistType::kTH1F, {{400, 0., 20.}}}},
     {"MC/reconstructed/nonprompt/hDca2RecSigNonPrompt", "3-prong candidates (matched);prong Chi2PCA to sec. vertex (cm);entries", {HistType::kTH1F, {{400, 0., 20.}}}},
     /// eta
     {"Data/hEta", "3-prong candidates;#it{#eta};entries", {HistType::kTH1F, {{100, -2., 2.}}}},
     {"MC/reconstructed/signal/hEtaRecSig", "3-prong candidates (matched);#it{#eta};entries", {HistType::kTH1F, {{100, -2., 2.}}}},
     {"MC/reconstructed/prompt/hEtaRecSigPrompt", "3-prong candidates (matched, prompt);#it{#eta};entries", {HistType::kTH1F, {{100, -2., 2.}}}},
     {"MC/reconstructed/nonprompt/hEtaRecSigNonPrompt", "3-prong candidates (matched, non-prompt);#it{#eta};entries", {HistType::kTH1F, {{100, -2., 2.}}}},
     {"MC/generated/signal/hEtaGen", "MC particles (matched);#it{#eta};entries", {HistType::kTH1F, {{100, -2., 2.}}}},
     {"MC/generated/prompt/hEtaGenPrompt", "MC particles (matched, prompt);#it{#eta};entries", {HistType::kTH1F, {{100, -2., 2.}}}},
     {"MC/generated/nonprompt/hEtaGenNonPrompt", "MC particles (matched, non-prompt);#it{#eta};entries", {HistType::kTH1F, {{100, -2., 2.}}}},
     {"MC/generated/signal/hYGen", "MC particles (matched);#it{y};entries", {HistType::kTH1F, {{100, -2., 2.}}}},
     {"MC/generated/prompt/hYGenPrompt", "MC particles (matched, prompt);#it{y};entries", {HistType::kTH1F, {{100, -2., 2.}}}},
     {"MC/generated/nonprompt/hYGenNonPrompt", "MC particles (matched, non-prompt);#it{y};entries", {HistType::kTH1F, {{100, -2., 2.}}}},
     /// phi
     {"Data/hPhi", "3-prong candidates;#it{#Phi};entries", {HistType::kTH1F, {{100, 0., 6.3}}}},
     {"MC/reconstructed/signal/hPhiRecSig", "3-prong candidates (matched);#it{#Phi};entries", {HistType::kTH1F, {{100, 0., 6.3}}}},
     {"MC/reconstructed/prompt/hPhiRecSigPrompt", "3-prong candidates (matched, prompt);#it{#Phi};entries", {HistType::kTH1F, {{100, 0., 6.3}}}},
     {"MC/reconstructed/nonprompt/hPhiRecSigNonPrompt", "3-prong candidates (matched, non-prompt);#it{#Phi};entries", {HistType::kTH1F, {{100, 0., 6.3}}}},
     {"MC/generated/signal/hPhiGen", "MC particles (matched);#it{#Phi};entries", {HistType::kTH1F, {{100, 0., 6.3}}}},
     {"MC/generated/prompt/hPhiGenPrompt", "MC particles (matched, prompt);#it{#Phi};entries", {HistType::kTH1F, {{100, 0., 6.3}}}},
     {"MC/generated/nonprompt/hPhiGenNonPrompt", "MC particles (matched, non-prompt);#it{#Phi};entries", {HistType::kTH1F, {{100, 0., 6.3}}}}}};

  void init(InitContext&)
  {
    std::array<bool, 12> doprocess{doprocessDataStd, doprocessDataStdWithFT0C, doprocessDataStdWithFT0M, doprocessDataWithMl, doprocessDataWithMlWithFT0C, doprocessDataWithMlWithFT0M, doprocessMcStd, doprocessMcStdWithFT0C, doprocessMcStdWithFT0M, doprocessMcWithMl, doprocessMcWithMlWithFT0C, doprocessMcWithMlWithFT0M};
    if ((std::accumulate(doprocess.begin(), doprocess.end(), 0)) != 1) {
      LOGP(fatal, "no or more than one process function enabled! Please check your configuration!");
    }

    auto vbins = (std::vector<double>)binsPt;
    /// mass candidate
    registry.add("Data/hMassVsPtVsNPvContributors", "3-prong candidates;inv. mass (p K #pi) (GeV/#it{c}^{2}); p_{T}; Number of PV contributors", {HistType::kTH3F, {{600, 1.98, 2.58}, {vbins, "#it{p}_{T} (GeV/#it{c})"}, {5000, 0., 10000.}}});
    registry.add("Data/hMassVsPt", "3-prong candidates;inv. mass (p K #pi) (GeV/#it{c}^{2}); p_{T}", {HistType::kTH2F, {{600, 1.98, 2.58}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hMassVsPtRecSig", "3-prong candidates (matched);inv. mass (p K #pi) (GeV/#it{c}^{2}); p_{T}", {HistType::kTH2F, {{600, 1.98, 2.58}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/prompt/hMassVsPtRecSigPrompt", "3-prong candidates (matched, prompt);inv. mass (p K #pi) (GeV/#it{c}^{2}); p_{T}", {HistType::kTH2F, {{600, 1.98, 2.58}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/nonprompt/hMassVsPtRecSigNonPrompt", "3-prong candidates (matched, non-prompt);inv. mass (p K #pi) (GeV/#it{c}^{2}); p_{T}", {HistType::kTH2F, {{600, 1.98, 2.58}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    /// DCAxy to prim. vertex prongs
    registry.add("Data/hd0VsPtProng0", "3-prong candidates;prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{600, -0.4, 0.4}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hd0VsPtRecProng0Sig", "3-prong candidates (matched);prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{600, -0.4, 0.4}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/prompt/hd0VsPtRecProng0SigPrompt", "3-prong candidates (matched, prompt);prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{600, -0.4, 0.4}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/nonprompt/hd0VsPtRecProng0SigNonPrompt", "3-prong candidates (matched, non-prompt);prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{600, -0.4, 0.4}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/hd0VsPtProng1", "3-prong candidates;prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{600, -0.4, 0.4}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hd0VsPtRecProng1Sig", "3-prong candidates (matched);prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{600, -0.4, 0.4}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/prompt/hd0VsPtRecProng1SigPrompt", "3-prong candidates (matched, prompt);prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{600, -0.4, 0.4}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/nonprompt/hd0VsPtRecProng1SigNonPrompt", "3-prong candidates (matched, non-prompt);prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{600, -0.4, 0.4}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/hd0VsPtProng2", "3-prong candidates;prong 2 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{600, -0.4, 0.4}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hd0VsPtRecProng2Sig", "3-prong candidates (matched);prong 2 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{600, -0.4, 0.4}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/prompt/hd0VsPtRecProng2SigPrompt", "3-prong candidates (matched, prompt);prong 2 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{600, -0.4, 0.4}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/nonprompt/hd0VsPtRecProng2SigNonPrompt", "3-prong candidates (matched, non-prompt);prong 2 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{600, -0.4, 0.4}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    /// decay length candidate
    registry.add("Data/hDecLengthVsPt", "3-prong candidates;decay length (cm);entries", {HistType::kTH2F, {{400, 0., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hDecLengthVsPtRecSig", "3-prong candidates (matched);decay length (cm);entries", {HistType::kTH2F, {{400, 0., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/prompt/hDecLengthVsPtRecSigPrompt", "3-prong candidates (matched, prompt);decay length (cm);entries", {HistType::kTH2F, {{400, 0., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/nonprompt/hDecLengthVsPtRecSigNonPrompt", "3-prong candidates (matched, non-prompt);decay length (cm);entries", {HistType::kTH2F, {{400, 0., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    /// decay length xy candidate
    registry.add("Data/hDecLengthxyVsPt", "3-prong candidates;decay length xy(cm);entries", {HistType::kTH2F, {{400, 0., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hDecLengthxyVsPtRecSig", "3-prong candidates (matched);decay length xy(cm);entries", {HistType::kTH2F, {{400, 0., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/prompt/hDecLengthxyVsPtRecSigPrompt", "3-prong candidates (matched, prompt);decay length xy(cm);entries", {HistType::kTH2F, {{400, 0., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/nonprompt/hDecLengthxyVsPtRecSigNonPrompt", "3-prong candidates (matched, non-prompt);decay length xy(cm);entries", {HistType::kTH2F, {{400, 0., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    /// proper lifetime
    registry.add("Data/hCtVsPt", "3-prong candidates;proper lifetime (#Lambda_{c}) * #it{c} (cm);entries", {HistType::kTH2F, {{100, 0., 0.2}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hCtVsPtRecSig", "3-prong candidates (matched);proper lifetime (#Lambda_{c}) * #it{c} (cm);entries", {HistType::kTH2F, {{100, 0., 0.2}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/prompt/hCtVsPtRecSigPrompt", "3-prong candidates (matched, prompt);proper lifetime (#Lambda_{c}) * #it{c} (cm);entries", {HistType::kTH2F, {{100, 0., 0.2}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/nonprompt/hCtVsPtRecSigNonPrompt", "3-prong candidates (matched, non-prompt);proper lifetime (#Lambda_{c}) * #it{c} (cm);entries", {HistType::kTH2F, {{100, 0., 0.2}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    /// cosine of pointing angle
    registry.add("Data/hCPAVsPt", "3-prong candidates;cosine of pointing angle;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hCPAVsPtRecSig", "3-prong candidates (matched);cosine of pointing angle;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/prompt/hCPAVsPtRecSigPrompt", "3-prong candidates (matched, prompt);cosine of pointing angle;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/nonprompt/hCPAVsPtRecSigNonPrompt", "3-prong candidates (matched, non-prompt);cosine of pointing angle;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    /// cosine of pointing angle xy
    registry.add("Data/hCPAxyVsPt", "3-prong candidates;cosine of pointing angle xy;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hCPAxyVsPtRecSig", "3-prong candidates (matched);cosine of pointing angle xy;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/prompt/hCPAxyVsPtRecSigPrompt", "3-prong candidates (matched, prompt);cosine of pointing angle xy;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/nonprompt/hCPAxyVsPtRecSigNonPrompt", "3-prong candidates (matched, non-prompt);cosine of pointing angle xy;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    /// Chi 2 PCA to sec. vertex
    registry.add("Data/hDca2VsPt", "3-prong candidates;prong Chi2PCA to sec. vertex (cm);entries", {HistType::kTH2F, {{400, 0., 20.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hDca2VsPtRecSig", "3-prong candidates (matched);prong Chi2PCA to sec. vertex (cm);entries", {HistType::kTH2F, {{400, 0., 20.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/prompt/hDca2VsPtRecSigPrompt", "3-prong candidates (matched, prompt);prong Chi2PCA to sec. vertex (cm);entries", {HistType::kTH2F, {{400, 0., 20.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/nonprompt/hDca2VsPtRecSigNonPrompt", "3-prong candidates (matched, non-prompt);prong Chi2PCA to sec. vertex (cm);entries", {HistType::kTH2F, {{400, 0., 20.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    /// eta
    registry.add("Data/hEtaVsPt", "3-prong candidates;candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hEtaVsPtRecSig", "3-prong candidates (matched);candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/prompt/hEtaVsPtRecSigPrompt", "3-prong candidates (matched, prompt);candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/nonprompt/hEtaVsPtRecSigNonPrompt", "3-prong candidates (matched, non-prompt);candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/generated/signal/hEtaVsPtGenSig", "3-prong candidates (matched);candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/generated/prompt/hEtaVsPtGenSigPrompt", "3-prong candidates (matched, prompt);candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/generated/nonprompt/hEtaVsPtGenSigNonPrompt", "3-prong candidates (matched, non-prompt);candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    /// y
    registry.add("MC/generated/signal/hYVsPtGenSig", "3-prong candidates (matched);candidate #it{y};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/generated/prompt/hYVsPtGenSigPrompt", "3-prong candidates (matched, prompt);candidate #it{y};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/generated/nonprompt/hYVsPtGenSigNonPrompt", "3-prong candidates (matched, non-prompt);candidate #it{y};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    /// phi
    registry.add("Data/hPhiVsPt", "3-prong candidates;candidate #it{#Phi};entries", {HistType::kTH2F, {{100, 0., 6.3}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hPhiVsPtRecSig", "3-prong candidates (matched);candidate #it{#Phi};entries", {HistType::kTH2F, {{100, 0., 6.3}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/prompt/hPhiVsPtRecSigPrompt", "3-prong candidates (matched, prompt);candidate #it{#Phi};entries", {HistType::kTH2F, {{100, 0., 6.3}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/nonprompt/hPhiVsPtRecSigNonPrompt", "3-prong candidates (matched, non-prompt);candidate #it{#Phi};entries", {HistType::kTH2F, {{100, 0., 6.3}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/generated/signal/hPhiVsPtGenSig", "3-prong candidates (matched);candidate #it{#Phi};entries", {HistType::kTH2F, {{100, 0., 6.3}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/generated/prompt/hPhiVsPtGenSigPrompt", "3-prong candidates (matched, prompt);candidate #it{#Phi};entries", {HistType::kTH2F, {{100, 0., 6.3}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/generated/nonprompt/hPhiVsPtGenSigNonPrompt", "3-prong candidates (matched, non-prompt);candidate #it{#Phi};entries", {HistType::kTH2F, {{100, 0., 6.3}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    /// selection status
    registry.add("hSelectionStatus", "3-prong candidates;selection status;entries", {HistType::kTH2F, {{5, -0.5, 4.5}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    /// impact parameter error
    registry.add("Data/hImpParErrProng0", "3-prong candidates;prong 0 impact parameter error (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/hImpParErrProng1", "3-prong candidates;prong 1 impact parameter error (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/hImpParErrProng2", "3-prong candidates;prong 2 impact parameter error (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hImpParErrProng0Sig", "3-prong candidates (matched);prong 0 impact parameter error (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/prompt/hImpParErrProng0SigPrompt", "3-prong candidates (matched, prompt);prong 0 impact parameter error (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/nonprompt/hImpParErrProng0SigNonPrompt", "3-prong candidates (matched, non-prompt);prong 0 impact parameter error (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hImpParErrProng1Sig", "3-prong candidates (matched);prong 1 impact parameter error (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/prompt/hImpParErrProng1SigPrompt", "3-prong candidates (matched, prompt);prong 1 impact parameter error (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/nonprompt/hImpParErrProng1SigNonPrompt", "3-prong candidates (matched, non-prompt);prong 1 impact parameter error (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hImpParErrProng2Sig", "3-prong candidates (matched);prong 2 impact parameter error (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/prompt/hImpParErrProng2SigPrompt", "3-prong candidates (matched, prompt);prong 2 impact parameter error (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/nonprompt/hImpParErrProng2SigNonPrompt", "3-prong candidates (matched, non-prompt);prong 2 impact parameter error (cm);entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    /// decay length error
    registry.add("Data/hDecLenErr", "3-prong candidates;decay length error (cm);entries", {HistType::kTH2F, {{100, 0., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hDecLenErrSig", "3-prong candidates (matched);decay length error (cm);entries", {HistType::kTH2F, {{100, 0., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/prompt/hDecLenErrSigPrompt", "3-prong candidates (matched, prompt);decay length error (cm);entries", {HistType::kTH2F, {{100, 0., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/nonprompt/hDecLenErrSigNonPrompt", "3-prong candidates (matched, non-prompt);decay length error (cm);entries", {HistType::kTH2F, {{100, 0., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    if (fillTHn) {
      const AxisSpec thnAxisMass{thnConfigAxisMass, "inv. mass (p K #pi) (GeV/#it{c}^{2})"};
      const AxisSpec thnAxisPt{thnConfigAxisPt, "#it{p}_{T}(#Lambda_{c}^{+}) (GeV/#it{c})"};
      const AxisSpec thnAxisPtProng0{thnConfigAxisPtProng, "#it{p}_{T}(prong0) (GeV/#it{c})"};
      const AxisSpec thnAxisPtProng1{thnConfigAxisPtProng, "#it{p}_{T}(prong1) (GeV/#it{c})"};
      const AxisSpec thnAxisPtProng2{thnConfigAxisPtProng, "#it{p}_{T}(prong2) (GeV/#it{c})"};
      const AxisSpec thnAxisCentrality{thnConfigAxisCentrality, "centrality (FT0C)"};
      const AxisSpec thnAxisChi2PCA{thnConfigAxisChi2PCA, "Chi2PCA to sec. vertex (cm)"};
      const AxisSpec thnAxisDecLength{thnConfigAxisDecLength, "decay length (cm)"};
      const AxisSpec thnAxisCPA{thnConfigAxisCPA, "cosine of pointing angle"};
      const AxisSpec thnAxisBdtScoreLcBkg{thnConfigAxisBdtScoreBkg, "BDT bkg score (Lc)"};
      const AxisSpec thnAxisBdtScoreLcPrompt{thnConfigAxisBdtScoreSignal, "BDT prompt score (Lc)"};
      const AxisSpec thnAxisBdtScoreLcNonPrompt{thnConfigAxisBdtScoreSignal, "BDT non-prompt score (Lc)"};
      const AxisSpec thnAxisCanType{thnConfigAxisCanType, "candidates type"};
      const AxisSpec thnAxisY{thnAxisRapidity, "rapidity"};
      const AxisSpec thnAxisPtB{thnConfigAxisGenPtB, "#it{p}_{T}^{B} (GeV/#it{c})"};
      const AxisSpec thnAxisTracklets{thnConfigAxisNumPvContr, "Number of PV contributors"};
      const AxisSpec thnAxisOccupancy{thnConfigAxisOccupancy, "Occupancy"};

      bool isDataWithMl = doprocessDataWithMl || doprocessDataWithMlWithFT0C || doprocessDataWithMlWithFT0M;
      bool isMcWithMl = doprocessMcWithMl || doprocessMcWithMlWithFT0C || doprocessMcWithMlWithFT0M;
      bool isDataStd = doprocessDataStd || doprocessDataStdWithFT0C || doprocessDataStdWithFT0M;
      bool isMcStd = doprocessMcStd || doprocessMcStdWithFT0C || doprocessMcStdWithFT0M;

      std::vector<AxisSpec> axesStd, axesWithBdt, axesGen;

      if (isDataStd) {
        axesStd = {thnAxisMass, thnAxisPt, thnAxisCentrality, thnAxisPtProng0, thnAxisPtProng1, thnAxisPtProng2, thnAxisChi2PCA, thnAxisDecLength, thnAxisCPA, thnAxisTracklets};
      }
      if (isMcStd) {
        axesStd = {thnAxisMass, thnAxisPt, thnAxisCentrality, thnAxisPtProng0, thnAxisPtProng1, thnAxisPtProng2, thnAxisChi2PCA, thnAxisDecLength, thnAxisCPA, thnAxisTracklets, thnAxisPtB, thnAxisCanType};
      }
      if (isMcStd || isMcWithMl) {
        axesGen = {thnAxisPt, thnAxisCentrality, thnAxisY, thnAxisTracklets, thnAxisPtB, thnAxisCanType};
      }
      if (isDataWithMl) {
        axesWithBdt = {thnAxisMass, thnAxisPt, thnAxisCentrality, thnAxisBdtScoreLcBkg, thnAxisBdtScoreLcPrompt, thnAxisBdtScoreLcNonPrompt, thnAxisTracklets};
      }
      if (isMcWithMl) {
        axesWithBdt = {thnAxisMass, thnAxisPt, thnAxisCentrality, thnAxisBdtScoreLcBkg, thnAxisBdtScoreLcPrompt, thnAxisBdtScoreLcNonPrompt, thnAxisTracklets, thnAxisPtB, thnAxisCanType};
      }

      if (storeOccupancy) {
        if (!axesWithBdt.empty())
          axesWithBdt.push_back(thnAxisOccupancy);
        if (!axesStd.empty())
          axesStd.push_back(thnAxisOccupancy);
        if (!axesGen.empty())
          axesGen.push_back(thnAxisOccupancy);
      }

      if (isDataWithMl) {
        registry.add("hnLcVarsWithBdt", "THn for Lambdac candidates with BDT scores for data with ML", HistType::kTHnSparseF, axesWithBdt);
      } else if (isMcWithMl) {
        registry.add("hnLcVarsWithBdt", "THn for Lambdac candidates with BDT scores for mc with ML", HistType::kTHnSparseF, axesWithBdt);
        registry.add("hnLcVarsGen", "THn for Generated Lambdac", HistType::kTHnSparseF, axesGen);
      } else if (isDataStd) {
        registry.add("hnLcVars", "THn for Reconstructed Lambdac candidates for data without ML", HistType::kTHnSparseF, axesStd);
      } else {
        registry.add("hnLcVars", "THn for Reconstructed Lambdac candidates for mc without ML", HistType::kTHnSparseF, axesStd);
        registry.add("hnLcVarsGen", "THn for Generated Lambdac", HistType::kTHnSparseF, axesGen);
      }
    }
  }

  /// Evaluate centrality/multiplicity percentile (centrality estimator is automatically selected based on the used table)
  /// \param candidate is candidate
  /// \return centrality/multiplicity percentile of the collision
  template <typename Coll>
  float evaluateCentralityColl(const Coll& collision)
  {
    return o2::hf_centrality::getCentralityColl<Coll>(collision);
  }

  /// Fill MC histograms at reconstruction level
  /// \tparam fillMl switch to fill ML histograms
  template <bool fillMl, typename CollType, typename CandLcMcRec, typename CandLcMcGen>
  void fillHistosMcRec(CollType const& collision, CandLcMcRec const& candidates, CandLcMcGen const& mcParticles)
  {

    auto thisCollId = collision.globalIndex();
    auto groupedLcCandidates = candidates.sliceBy(candLcPerCollision, thisCollId);

    for (const auto& candidate : groupedLcCandidates) {
      /// Select Lc
      if (!(candidate.hfflag() & 1 << aod::hf_cand_3prong::DecayType::LcToPKPi)) {
        continue;
      }
      /// rapidity selection
      if (yCandRecoMax >= 0. && std::abs(hfHelper.yLc(candidate)) > yCandRecoMax) {
        continue;
      }

      if (std::abs(candidate.flagMcMatchRec()) == 1 << aod::hf_cand_3prong::DecayType::LcToPKPi) {
        // Get the corresponding MC particle.
        auto mcParticleProng0 = candidate.template prong0_as<aod::TracksWMc>().template mcParticle_as<soa::Join<aod::McParticles, aod::HfCand3ProngMcGen>>();
        auto pdgCodeProng0 = std::abs(mcParticleProng0.pdgCode());
        auto indexMother = RecoDecay::getMother(mcParticles, mcParticleProng0, o2::constants::physics::Pdg::kLambdaCPlus, true);
        auto particleMother = mcParticles.rawIteratorAt(indexMother);
        registry.fill(HIST("MC/generated/signal/hPtGenSig"), particleMother.pt()); // gen. level pT

        auto pt = candidate.pt();
        auto ptProng0 = candidate.ptProng0();
        auto ptProng1 = candidate.ptProng1();
        auto ptProng2 = candidate.ptProng2();
        auto decayLength = candidate.decayLength();
        auto decayLengthXY = candidate.decayLengthXY();
        auto chi2PCA = candidate.chi2PCA();
        auto cpa = candidate.cpa();
        auto cpaXY = candidate.cpaXY();
        auto originType = candidate.originMcRec();
        auto numPvContributors = collision.numContrib();
        auto ptRecB = candidate.ptBhadMotherPart();

        /// MC reconstructed signal
        if ((candidate.isSelLcToPKPi() >= selectionFlagLc) && pdgCodeProng0 == kProton) {
          registry.fill(HIST("MC/reconstructed/signal/hMassRecSig"), hfHelper.invMassLcToPKPi(candidate));
          registry.fill(HIST("MC/reconstructed/signal/hMassVsPtRecSig"), hfHelper.invMassLcToPKPi(candidate), pt);
        }
        if ((candidate.isSelLcToPiKP() >= selectionFlagLc) && pdgCodeProng0 == kPiPlus) {
          registry.fill(HIST("MC/reconstructed/signal/hMassRecSig"), hfHelper.invMassLcToPiKP(candidate));
          registry.fill(HIST("MC/reconstructed/signal/hMassVsPtRecSig"), hfHelper.invMassLcToPiKP(candidate), pt);
        }
        registry.fill(HIST("MC/reconstructed/signal/hPtRecSig"), pt);
        registry.fill(HIST("MC/reconstructed/signal/hPtRecProng0Sig"), ptProng0);
        registry.fill(HIST("MC/reconstructed/signal/hPtRecProng1Sig"), ptProng1);
        registry.fill(HIST("MC/reconstructed/signal/hPtRecProng2Sig"), ptProng2);

        registry.fill(HIST("MC/reconstructed/signal/hd0RecProng0Sig"), candidate.impactParameter0());
        registry.fill(HIST("MC/reconstructed/signal/hd0RecProng1Sig"), candidate.impactParameter1());
        registry.fill(HIST("MC/reconstructed/signal/hd0RecProng2Sig"), candidate.impactParameter2());
        registry.fill(HIST("MC/reconstructed/signal/hd0VsPtRecProng0Sig"), candidate.impactParameter0(), pt);
        registry.fill(HIST("MC/reconstructed/signal/hd0VsPtRecProng1Sig"), candidate.impactParameter1(), pt);
        registry.fill(HIST("MC/reconstructed/signal/hd0VsPtRecProng2Sig"), candidate.impactParameter2(), pt);
        registry.fill(HIST("MC/reconstructed/signal/hDecLengthRecSig"), decayLength);
        registry.fill(HIST("MC/reconstructed/signal/hDecLengthVsPtRecSig"), decayLength, pt);
        registry.fill(HIST("MC/reconstructed/signal/hDecLengthxyRecSig"), decayLengthXY);
        registry.fill(HIST("MC/reconstructed/signal/hDecLengthxyVsPtRecSig"), decayLengthXY, pt);
        registry.fill(HIST("MC/reconstructed/signal/hCtRecSig"), hfHelper.ctLc(candidate));
        registry.fill(HIST("MC/reconstructed/signal/hCtVsPtRecSig"), hfHelper.ctLc(candidate), pt);
        registry.fill(HIST("MC/reconstructed/signal/hCPARecSig"), cpa);
        registry.fill(HIST("MC/reconstructed/signal/hCPAVsPtRecSig"), cpa, pt);
        registry.fill(HIST("MC/reconstructed/signal/hCPAxyRecSig"), cpaXY);
        registry.fill(HIST("MC/reconstructed/signal/hCPAxyVsPtRecSig"), cpaXY, pt);
        registry.fill(HIST("MC/reconstructed/signal/hDca2RecSig"), chi2PCA);
        registry.fill(HIST("MC/reconstructed/signal/hDca2VsPtRecSig"), chi2PCA, pt);
        registry.fill(HIST("MC/reconstructed/signal/hEtaRecSig"), candidate.eta());
        registry.fill(HIST("MC/reconstructed/signal/hEtaVsPtRecSig"), candidate.eta(), pt);
        registry.fill(HIST("MC/reconstructed/signal/hPhiRecSig"), candidate.phi());
        registry.fill(HIST("MC/reconstructed/signal/hPhiVsPtRecSig"), candidate.phi(), pt);
        registry.fill(HIST("MC/reconstructed/signal/hImpParErrProng0Sig"), candidate.errorImpactParameter0(), pt);
        registry.fill(HIST("MC/reconstructed/signal/hImpParErrProng1Sig"), candidate.errorImpactParameter1(), pt);
        registry.fill(HIST("MC/reconstructed/signal/hImpParErrProng2Sig"), candidate.errorImpactParameter2(), pt);
        registry.fill(HIST("MC/reconstructed/signal/hDecLenErrSig"), candidate.errorDecayLength(), pt);

        /// reconstructed signal prompt
        if (candidate.originMcRec() == RecoDecay::OriginType::Prompt) {
          if ((candidate.isSelLcToPKPi() >= selectionFlagLc) && pdgCodeProng0 == kProton) {
            registry.fill(HIST("MC/reconstructed/prompt/hMassRecSigPrompt"), hfHelper.invMassLcToPKPi(candidate));
            registry.fill(HIST("MC/reconstructed/prompt/hMassVsPtRecSigPrompt"), hfHelper.invMassLcToPKPi(candidate), pt);
          }
          if ((candidate.isSelLcToPiKP() >= selectionFlagLc) && pdgCodeProng0 == kPiPlus) {
            registry.fill(HIST("MC/reconstructed/prompt/hMassRecSigPrompt"), hfHelper.invMassLcToPiKP(candidate));
            registry.fill(HIST("MC/reconstructed/prompt/hMassVsPtRecSigPrompt"), hfHelper.invMassLcToPiKP(candidate), pt);
          }
          registry.fill(HIST("MC/reconstructed/prompt/hPtRecSigPrompt"), pt);
          registry.fill(HIST("MC/reconstructed/prompt/hPtRecProng0SigPrompt"), ptProng0);
          registry.fill(HIST("MC/reconstructed/prompt/hPtRecProng1SigPrompt"), ptProng1);
          registry.fill(HIST("MC/reconstructed/prompt/hPtRecProng2SigPrompt"), ptProng2);
          registry.fill(HIST("MC/reconstructed/prompt/hd0RecProng0SigPrompt"), candidate.impactParameter0());
          registry.fill(HIST("MC/reconstructed/prompt/hd0RecProng1SigPrompt"), candidate.impactParameter1());
          registry.fill(HIST("MC/reconstructed/prompt/hd0RecProng2SigPrompt"), candidate.impactParameter2());
          registry.fill(HIST("MC/reconstructed/prompt/hd0VsPtRecProng0SigPrompt"), candidate.impactParameter0(), pt);
          registry.fill(HIST("MC/reconstructed/prompt/hd0VsPtRecProng1SigPrompt"), candidate.impactParameter1(), pt);
          registry.fill(HIST("MC/reconstructed/prompt/hd0VsPtRecProng2SigPrompt"), candidate.impactParameter2(), pt);
          registry.fill(HIST("MC/reconstructed/prompt/hDecLengthRecSigPrompt"), decayLength);
          registry.fill(HIST("MC/reconstructed/prompt/hDecLengthVsPtRecSigPrompt"), decayLength, pt);
          registry.fill(HIST("MC/reconstructed/prompt/hDecLengthxyRecSigPrompt"), decayLengthXY);
          registry.fill(HIST("MC/reconstructed/prompt/hDecLengthxyVsPtRecSigPrompt"), decayLengthXY, pt);
          registry.fill(HIST("MC/reconstructed/prompt/hCtRecSigPrompt"), hfHelper.ctLc(candidate));
          registry.fill(HIST("MC/reconstructed/prompt/hCtVsPtRecSigPrompt"), hfHelper.ctLc(candidate), pt);
          registry.fill(HIST("MC/reconstructed/prompt/hCPARecSigPrompt"), cpa);
          registry.fill(HIST("MC/reconstructed/prompt/hCPAVsPtRecSigPrompt"), cpa, pt);
          registry.fill(HIST("MC/reconstructed/prompt/hCPAxyRecSigPrompt"), cpaXY);
          registry.fill(HIST("MC/reconstructed/prompt/hCPAxyVsPtRecSigPrompt"), cpaXY, pt);
          registry.fill(HIST("MC/reconstructed/prompt/hDca2RecSigPrompt"), chi2PCA);
          registry.fill(HIST("MC/reconstructed/prompt/hDca2VsPtRecSigPrompt"), chi2PCA, pt);
          registry.fill(HIST("MC/reconstructed/prompt/hEtaRecSigPrompt"), candidate.eta());
          registry.fill(HIST("MC/reconstructed/prompt/hEtaVsPtRecSigPrompt"), candidate.eta(), pt);
          registry.fill(HIST("MC/reconstructed/prompt/hPhiRecSigPrompt"), candidate.phi());
          registry.fill(HIST("MC/reconstructed/prompt/hPhiVsPtRecSigPrompt"), candidate.phi(), pt);
          registry.fill(HIST("MC/reconstructed/prompt/hImpParErrProng0SigPrompt"), candidate.errorImpactParameter0(), pt);
          registry.fill(HIST("MC/reconstructed/prompt/hImpParErrProng1SigPrompt"), candidate.errorImpactParameter1(), pt);
          registry.fill(HIST("MC/reconstructed/prompt/hImpParErrProng2SigPrompt"), candidate.errorImpactParameter2(), pt);
          registry.fill(HIST("MC/reconstructed/prompt/hDecLenErrSigPrompt"), candidate.errorDecayLength(), pt);
        } else {
          if ((candidate.isSelLcToPKPi() >= selectionFlagLc) && pdgCodeProng0 == kProton) {
            registry.fill(HIST("MC/reconstructed/nonprompt/hMassRecSigNonPrompt"), hfHelper.invMassLcToPKPi(candidate));
            registry.fill(HIST("MC/reconstructed/nonprompt/hMassVsPtRecSigNonPrompt"), hfHelper.invMassLcToPKPi(candidate), pt);
          }
          if ((candidate.isSelLcToPiKP() >= selectionFlagLc) && pdgCodeProng0 == kPiPlus) {
            registry.fill(HIST("MC/reconstructed/nonprompt/hMassRecSigNonPrompt"), hfHelper.invMassLcToPiKP(candidate));
            registry.fill(HIST("MC/reconstructed/nonprompt/hMassVsPtRecSigNonPrompt"), hfHelper.invMassLcToPiKP(candidate), pt);
          }
          registry.fill(HIST("MC/reconstructed/nonprompt/hPtRecSigNonPrompt"), pt);
          registry.fill(HIST("MC/reconstructed/nonprompt/hPtRecProng0SigNonPrompt"), ptProng0);
          registry.fill(HIST("MC/reconstructed/nonprompt/hPtRecProng1SigNonPrompt"), ptProng1);
          registry.fill(HIST("MC/reconstructed/nonprompt/hPtRecProng2SigNonPrompt"), ptProng2);
          registry.fill(HIST("MC/reconstructed/nonprompt/hd0RecProng0SigNonPrompt"), candidate.impactParameter0());
          registry.fill(HIST("MC/reconstructed/nonprompt/hd0RecProng1SigNonPrompt"), candidate.impactParameter1());
          registry.fill(HIST("MC/reconstructed/nonprompt/hd0RecProng2SigNonPrompt"), candidate.impactParameter2());
          registry.fill(HIST("MC/reconstructed/nonprompt/hd0VsPtRecProng0SigNonPrompt"), candidate.impactParameter0(), pt);
          registry.fill(HIST("MC/reconstructed/nonprompt/hd0VsPtRecProng1SigNonPrompt"), candidate.impactParameter1(), pt);
          registry.fill(HIST("MC/reconstructed/nonprompt/hd0VsPtRecProng2SigNonPrompt"), candidate.impactParameter2(), pt);
          registry.fill(HIST("MC/reconstructed/nonprompt/hDecLengthRecSigNonPrompt"), decayLength);
          registry.fill(HIST("MC/reconstructed/nonprompt/hDecLengthVsPtRecSigNonPrompt"), decayLength, pt);
          registry.fill(HIST("MC/reconstructed/nonprompt/hDecLengthxyRecSigNonPrompt"), decayLengthXY);
          registry.fill(HIST("MC/reconstructed/nonprompt/hDecLengthxyVsPtRecSigNonPrompt"), decayLengthXY, pt);
          registry.fill(HIST("MC/reconstructed/nonprompt/hCtRecSigNonPrompt"), hfHelper.ctLc(candidate));
          registry.fill(HIST("MC/reconstructed/nonprompt/hCtVsPtRecSigNonPrompt"), hfHelper.ctLc(candidate), pt);
          registry.fill(HIST("MC/reconstructed/nonprompt/hCPARecSigNonPrompt"), cpa);
          registry.fill(HIST("MC/reconstructed/nonprompt/hCPAVsPtRecSigNonPrompt"), cpa, pt);
          registry.fill(HIST("MC/reconstructed/nonprompt/hCPAxyRecSigNonPrompt"), cpaXY);
          registry.fill(HIST("MC/reconstructed/nonprompt/hCPAxyVsPtRecSigNonPrompt"), cpaXY, pt);
          registry.fill(HIST("MC/reconstructed/nonprompt/hDca2RecSigNonPrompt"), chi2PCA);
          registry.fill(HIST("MC/reconstructed/nonprompt/hDca2VsPtRecSigNonPrompt"), chi2PCA, pt);
          registry.fill(HIST("MC/reconstructed/nonprompt/hEtaRecSigNonPrompt"), candidate.eta());
          registry.fill(HIST("MC/reconstructed/nonprompt/hEtaVsPtRecSigNonPrompt"), candidate.eta(), pt);
          registry.fill(HIST("MC/reconstructed/nonprompt/hPhiRecSigNonPrompt"), candidate.phi());
          registry.fill(HIST("MC/reconstructed/nonprompt/hPhiVsPtRecSigNonPrompt"), candidate.phi(), pt);
          registry.fill(HIST("MC/reconstructed/nonprompt/hImpParErrProng0SigNonPrompt"), candidate.errorImpactParameter0(), pt);
          registry.fill(HIST("MC/reconstructed/nonprompt/hImpParErrProng1SigNonPrompt"), candidate.errorImpactParameter1(), pt);
          registry.fill(HIST("MC/reconstructed/nonprompt/hImpParErrProng2SigNonPrompt"), candidate.errorImpactParameter2(), pt);
          registry.fill(HIST("MC/reconstructed/nonprompt/hDecLenErrSigNonPrompt"), candidate.errorDecayLength(), pt);
        }
        if (fillTHn) {
          float cent = evaluateCentralityColl(collision);
          float occ{-1.};
          if (storeOccupancy && occEstimator != o2::hf_occupancy::OccupancyEstimator::None) {
            occ = o2::hf_occupancy::getOccupancyColl(collision, occEstimator);
          }
          double massLc(-1);
          double outputBkg(-1), outputPrompt(-1), outputFD(-1);
          if ((candidate.isSelLcToPKPi() >= selectionFlagLc) && pdgCodeProng0 == kProton) {
            massLc = hfHelper.invMassLcToPKPi(candidate);

            if constexpr (fillMl) {
              if (candidate.mlProbLcToPKPi().size() == 3) {
                outputBkg = candidate.mlProbLcToPKPi()[0];    /// bkg score
                outputPrompt = candidate.mlProbLcToPKPi()[1]; /// prompt score
                outputFD = candidate.mlProbLcToPKPi()[2];     /// non-prompt score
              }
              /// Fill the ML outputScores and variables of candidate
              if (storeOccupancy && occEstimator != o2::hf_occupancy::OccupancyEstimator::None) {
                registry.get<THnSparse>(HIST("hnLcVarsWithBdt"))->Fill(massLc, pt, cent, outputBkg, outputPrompt, outputFD, numPvContributors, ptRecB, originType, occ);
              } else {
                registry.get<THnSparse>(HIST("hnLcVarsWithBdt"))->Fill(massLc, pt, cent, outputBkg, outputPrompt, outputFD, numPvContributors, ptRecB, originType);
              }

            } else {

              if (storeOccupancy && occEstimator != o2::hf_occupancy::OccupancyEstimator::None) {
                registry.get<THnSparse>(HIST("hnLcVars"))->Fill(massLc, pt, cent, ptProng0, ptProng1, ptProng2, chi2PCA, decayLength, cpa, numPvContributors, ptRecB, originType, occ);

              } else {

                registry.get<THnSparse>(HIST("hnLcVars"))->Fill(massLc, pt, cent, ptProng0, ptProng1, ptProng2, chi2PCA, decayLength, cpa, numPvContributors, ptRecB, originType);
              }
            }
          }
          if ((candidate.isSelLcToPiKP() >= selectionFlagLc) && pdgCodeProng0 == kPiPlus) {
            massLc = hfHelper.invMassLcToPiKP(candidate);

            if constexpr (fillMl) {
              if (candidate.mlProbLcToPiKP().size() == 3) {
                outputBkg = candidate.mlProbLcToPiKP()[0];    /// bkg score
                outputPrompt = candidate.mlProbLcToPiKP()[1]; /// prompt score
                outputFD = candidate.mlProbLcToPiKP()[2];     /// non-prompt score
              }
              /// Fill the ML outputScores and variables of candidate (todo: add multiplicity)
              if (storeOccupancy && occEstimator != o2::hf_occupancy::OccupancyEstimator::None) {
                registry.get<THnSparse>(HIST("hnLcVarsWithBdt"))->Fill(massLc, pt, cent, outputBkg, outputPrompt, outputFD, numPvContributors, ptRecB, originType, occ);

              } else {
                registry.get<THnSparse>(HIST("hnLcVarsWithBdt"))->Fill(massLc, pt, cent, outputBkg, outputPrompt, outputFD, numPvContributors, ptRecB, originType);
              }
            } else {
              if (storeOccupancy && occEstimator != o2::hf_occupancy::OccupancyEstimator::None) {
                registry.get<THnSparse>(HIST("hnLcVars"))->Fill(massLc, pt, cent, ptProng0, ptProng1, ptProng2, chi2PCA, decayLength, cpa, numPvContributors, ptRecB, originType, occ);
              } else {
                registry.get<THnSparse>(HIST("hnLcVars"))->Fill(massLc, pt, cent, ptProng0, ptProng1, ptProng2, chi2PCA, decayLength, cpa, numPvContributors, ptRecB, originType);
              }
            }
          }
        }
      }
    }
  }

  /// Fill MC histograms at generated level
  /// \tparam fillMl switch to fill ML histograms
  template <typename CandLcMcGen, typename Coll>
  void fillHistosMcGen(CandLcMcGen const& mcParticles, Coll const& recoCollisions)
  {
    // MC gen.
    for (const auto& particle : mcParticles) {
      if (std::abs(particle.flagMcMatchGen()) == 1 << aod::hf_cand_3prong::DecayType::LcToPKPi) {
        auto yGen = RecoDecay::y(particle.pVector(), o2::constants::physics::MassLambdaCPlus);
        if (yCandGenMax >= 0. && std::abs(yGen) > yCandGenMax) {
          continue;
        }
        auto ptGen = particle.pt();
        auto originType = particle.originMcGen();
        float ptGenB = -1.;
        unsigned int numPvContributors = 0;
        const auto& recoCollsPerMcColl = recoCollisions.sliceBy(colPerMcCollision, particle.mcCollision().globalIndex());
        for (const auto& recCol : recoCollsPerMcColl) {
          numPvContributors = recCol.numContrib() > numPvContributors ? recCol.numContrib() : numPvContributors;
        }
        float cent = o2::hf_centrality::getCentralityGenColl(recoCollsPerMcColl);
        float occ{-1.};
        if (storeOccupancy && occEstimator != o2::hf_occupancy::OccupancyEstimator::None) {
          occ = o2::hf_occupancy::getOccupancyGenColl(recoCollsPerMcColl, occEstimator);
        }

        registry.fill(HIST("MC/generated/signal/hPtGen"), ptGen);
        registry.fill(HIST("MC/generated/signal/hEtaGen"), particle.eta());
        registry.fill(HIST("MC/generated/signal/hYGen"), yGen);
        registry.fill(HIST("MC/generated/signal/hPhiGen"), particle.phi());
        registry.fill(HIST("MC/generated/signal/hEtaVsPtGenSig"), particle.eta(), ptGen);
        registry.fill(HIST("MC/generated/signal/hYVsPtGenSig"), yGen, ptGen);
        registry.fill(HIST("MC/generated/signal/hPhiVsPtGenSig"), particle.phi(), ptGen);

        if (particle.originMcGen() == RecoDecay::OriginType::Prompt) {
          if (fillTHn) {
            if (storeOccupancy && occEstimator != o2::hf_occupancy::OccupancyEstimator::None) {
              registry.get<THnSparse>(HIST("hnLcVarsGen"))->Fill(ptGen, cent, yGen, numPvContributors, ptGenB, originType, occ);

            } else {
              registry.get<THnSparse>(HIST("hnLcVarsGen"))->Fill(ptGen, cent, yGen, numPvContributors, ptGenB, originType);
            }
          }
          registry.fill(HIST("MC/generated/prompt/hPtGenPrompt"), ptGen);
          registry.fill(HIST("MC/generated/prompt/hEtaGenPrompt"), particle.eta());
          registry.fill(HIST("MC/generated/prompt/hYGenPrompt"), yGen);
          registry.fill(HIST("MC/generated/prompt/hPhiGenPrompt"), particle.phi());
          registry.fill(HIST("MC/generated/prompt/hEtaVsPtGenSigPrompt"), particle.eta(), ptGen);
          registry.fill(HIST("MC/generated/prompt/hYVsPtGenSigPrompt"), yGen, ptGen);
          registry.fill(HIST("MC/generated/prompt/hPhiVsPtGenSigPrompt"), particle.phi(), ptGen);
        }
        if (particle.originMcGen() == RecoDecay::OriginType::NonPrompt) {
          ptGenB = mcParticles.rawIteratorAt(particle.idxBhadMotherPart()).pt();
          if (fillTHn) {
            if (storeOccupancy && occEstimator != o2::hf_occupancy::OccupancyEstimator::None) {
              registry.get<THnSparse>(HIST("hnLcVarsGen"))->Fill(ptGen, cent, yGen, numPvContributors, ptGenB, originType, occ);

            } else {
              registry.get<THnSparse>(HIST("hnLcVarsGen"))->Fill(ptGen, cent, yGen, numPvContributors, ptGenB, originType);
            }
          }
          registry.fill(HIST("MC/generated/nonprompt/hPtGenNonPrompt"), ptGen);
          registry.fill(HIST("MC/generated/nonprompt/hEtaGenNonPrompt"), particle.eta());
          registry.fill(HIST("MC/generated/nonprompt/hYGenNonPrompt"), yGen);
          registry.fill(HIST("MC/generated/nonprompt/hPhiGenNonPrompt"), particle.phi());
          registry.fill(HIST("MC/generated/nonprompt/hEtaVsPtGenSigNonPrompt"), particle.eta(), ptGen);
          registry.fill(HIST("MC/generated/nonprompt/hYVsPtGenSigNonPrompt"), yGen, ptGen);
          registry.fill(HIST("MC/generated/nonprompt/hPhiVsPtGenSigNonPrompt"), particle.phi(), ptGen);
        }
      }
    }
  }

  /// Fill histograms for real data
  /// \tparam fillMl switch to fill ML histograms
  template <bool fillMl, typename CollType, typename CandType>
  void fillHistosData(CollType const& collision, CandType const& candidates)
  {
    auto thisCollId = collision.globalIndex();
    auto groupedLcCandidates = candidates.sliceBy(candLcPerCollision, thisCollId);
    auto numPvContributors = collision.numContrib();

    for (const auto& candidate : groupedLcCandidates) {
      if (!(candidate.hfflag() & 1 << aod::hf_cand_3prong::DecayType::LcToPKPi)) {
        continue;
      }
      if (yCandRecoMax >= 0. && std::abs(hfHelper.yLc(candidate)) > yCandRecoMax) {
        continue;
      }
      auto pt = candidate.pt();
      auto ptProng0 = candidate.ptProng0();
      auto ptProng1 = candidate.ptProng1();
      auto ptProng2 = candidate.ptProng2();
      auto decayLength = candidate.decayLength();
      auto decayLengthXY = candidate.decayLengthXY();
      auto chi2PCA = candidate.chi2PCA();
      auto cpa = candidate.cpa();
      auto cpaXY = candidate.cpaXY();

      if (candidate.isSelLcToPKPi() >= selectionFlagLc) {
        registry.fill(HIST("Data/hMass"), hfHelper.invMassLcToPKPi(candidate));
        registry.fill(HIST("Data/hMassVsPtVsNPvContributors"), hfHelper.invMassLcToPKPi(candidate), pt, numPvContributors);
        registry.fill(HIST("Data/hMassVsPt"), hfHelper.invMassLcToPKPi(candidate), pt);
      }
      if (candidate.isSelLcToPiKP() >= selectionFlagLc) {
        registry.fill(HIST("Data/hMass"), hfHelper.invMassLcToPiKP(candidate));
        registry.fill(HIST("Data/hMassVsPtVsNPvContributors"), hfHelper.invMassLcToPiKP(candidate), pt, numPvContributors);
        registry.fill(HIST("Data/hMassVsPt"), hfHelper.invMassLcToPiKP(candidate), pt);
      }
      registry.fill(HIST("Data/hPt"), pt);
      registry.fill(HIST("Data/hPtProng0"), ptProng0);
      registry.fill(HIST("Data/hPtProng1"), ptProng1);
      registry.fill(HIST("Data/hPtProng2"), ptProng2);
      registry.fill(HIST("Data/hd0Prong0"), candidate.impactParameter0());
      registry.fill(HIST("Data/hd0Prong1"), candidate.impactParameter1());
      registry.fill(HIST("Data/hd0Prong2"), candidate.impactParameter2());
      registry.fill(HIST("Data/hd0VsPtProng0"), candidate.impactParameter0(), pt);
      registry.fill(HIST("Data/hd0VsPtProng1"), candidate.impactParameter1(), pt);
      registry.fill(HIST("Data/hd0VsPtProng2"), candidate.impactParameter2(), pt);
      registry.fill(HIST("Data/hDecLength"), decayLength);
      registry.fill(HIST("Data/hDecLengthVsPt"), decayLength, pt);
      registry.fill(HIST("Data/hDecLengthxy"), decayLengthXY);
      registry.fill(HIST("Data/hDecLengthxyVsPt"), decayLengthXY, pt);
      registry.fill(HIST("Data/hCt"), hfHelper.ctLc(candidate));
      registry.fill(HIST("Data/hCtVsPt"), hfHelper.ctLc(candidate), pt);
      registry.fill(HIST("Data/hCPA"), cpa);
      registry.fill(HIST("Data/hCPAVsPt"), cpa, pt);
      registry.fill(HIST("Data/hCPAxy"), cpaXY);
      registry.fill(HIST("Data/hCPAxyVsPt"), cpaXY, pt);
      registry.fill(HIST("Data/hDca2"), chi2PCA);
      registry.fill(HIST("Data/hDca2VsPt"), chi2PCA, pt);
      registry.fill(HIST("Data/hEta"), candidate.eta());
      registry.fill(HIST("Data/hEtaVsPt"), candidate.eta(), pt);
      registry.fill(HIST("Data/hPhi"), candidate.phi());
      registry.fill(HIST("Data/hPhiVsPt"), candidate.phi(), pt);
      registry.fill(HIST("hSelectionStatus"), candidate.isSelLcToPKPi(), pt);
      registry.fill(HIST("hSelectionStatus"), candidate.isSelLcToPiKP(), pt);
      registry.fill(HIST("Data/hImpParErrProng0"), candidate.errorImpactParameter0(), pt);
      registry.fill(HIST("Data/hImpParErrProng1"), candidate.errorImpactParameter1(), pt);
      registry.fill(HIST("Data/hImpParErrProng2"), candidate.errorImpactParameter2(), pt);
      registry.fill(HIST("Data/hDecLenErr"), candidate.errorDecayLength(), pt);

      if (fillTHn) {
        float cent = evaluateCentralityColl(collision);
        float occ{-1.};
        if (storeOccupancy && occEstimator != o2::hf_occupancy::OccupancyEstimator::None) {
          occ = o2::hf_occupancy::getOccupancyColl(collision, occEstimator);
        }
        double massLc(-1);
        double outputBkg(-1), outputPrompt(-1), outputFD(-1);
        if (candidate.isSelLcToPKPi() >= selectionFlagLc) {
          massLc = hfHelper.invMassLcToPKPi(candidate);

          if constexpr (fillMl) {
            if (candidate.mlProbLcToPKPi().size() == 3) {
              outputBkg = candidate.mlProbLcToPKPi()[0];    /// bkg score
              outputPrompt = candidate.mlProbLcToPKPi()[1]; /// prompt score
              outputFD = candidate.mlProbLcToPKPi()[2];     /// non-prompt score
            }
            /// Fill the ML outputScores and variables of candidate
            if (storeOccupancy && occEstimator != o2::hf_occupancy::OccupancyEstimator::None) {
              registry.get<THnSparse>(HIST("hnLcVarsWithBdt"))->Fill(massLc, pt, cent, outputBkg, outputPrompt, outputFD, numPvContributors, occ);

            } else {
              registry.get<THnSparse>(HIST("hnLcVarsWithBdt"))->Fill(massLc, pt, cent, outputBkg, outputPrompt, outputFD, numPvContributors);
            }
          } else {
            if (storeOccupancy && occEstimator != o2::hf_occupancy::OccupancyEstimator::None) {

              registry.get<THnSparse>(HIST("hnLcVars"))->Fill(massLc, pt, cent, ptProng0, ptProng1, ptProng2, chi2PCA, decayLength, cpa, numPvContributors, occ);
            } else {

              registry.get<THnSparse>(HIST("hnLcVars"))->Fill(massLc, pt, cent, ptProng0, ptProng1, ptProng2, chi2PCA, decayLength, cpa, numPvContributors);
            }
          }
        }
        if (candidate.isSelLcToPiKP() >= selectionFlagLc) {
          massLc = hfHelper.invMassLcToPiKP(candidate);

          if constexpr (fillMl) {
            if (candidate.mlProbLcToPiKP().size() == 3) {
              outputBkg = candidate.mlProbLcToPiKP()[0];    /// bkg score
              outputPrompt = candidate.mlProbLcToPiKP()[1]; /// prompt score
              outputFD = candidate.mlProbLcToPiKP()[2];     /// non-prompt score
            }
            /// Fill the ML outputScores and variables of candidate
            if (storeOccupancy && occEstimator != o2::hf_occupancy::OccupancyEstimator::None) {
              registry.get<THnSparse>(HIST("hnLcVarsWithBdt"))->Fill(massLc, pt, cent, outputBkg, outputPrompt, outputFD, numPvContributors, occ);
            } else {
              registry.get<THnSparse>(HIST("hnLcVarsWithBdt"))->Fill(massLc, pt, cent, outputBkg, outputPrompt, outputFD, numPvContributors);
            }
          } else {
            if (storeOccupancy && occEstimator != o2::hf_occupancy::OccupancyEstimator::None) {
              registry.get<THnSparse>(HIST("hnLcVars"))->Fill(massLc, pt, cent, ptProng0, ptProng1, ptProng2, chi2PCA, decayLength, cpa, numPvContributors, occ);
            } else {
              registry.get<THnSparse>(HIST("hnLcVars"))->Fill(massLc, pt, cent, ptProng0, ptProng1, ptProng2, chi2PCA, decayLength, cpa, numPvContributors);
            }
          }
        }
      }
    }
  }
  /// Run the analysis on real data
  /// \tparam fillMl switch to fill ML histograms
  template <bool fillMl, typename CollType, typename CandType>
  void runAnalysisPerCollisionData(CollType const& collisions,
                                   CandType const& candidates)
  {

    for (const auto& collision : collisions) {
      fillHistosData<fillMl>(collision, candidates);
    }
  }

  /// Run the analysis on MC data
  /// \tparam fillMl switch to fill ML histograms
  template <bool fillMl, typename CollType, typename CandType, typename CandLcMcGen>
  void runAnalysisPerCollisionMc(CollType const& collisions,
                                 CandType const& candidates,
                                 CandLcMcGen const& mcParticles)
  {
    for (const auto& collision : collisions) {
      // MC Rec.
      fillHistosMcRec<fillMl>(collision, candidates, mcParticles);
    }
    // MC gen.
    fillHistosMcGen(mcParticles, collisions);
  }

  void processDataStd(Collisions const& collisions,
                      LcCandidates const& selectedLcCandidates,
                      aod::Tracks const&)
  {
    runAnalysisPerCollisionData<false>(collisions, selectedLcCandidates);
  }
  PROCESS_SWITCH(HfTaskLc, processDataStd, "Process Data with the standard method", true);

  void processDataWithMl(Collisions const& collisions,
                         LcCandidatesMl const& selectedLcCandidatesMl,
                         aod::Tracks const&)
  {
    runAnalysisPerCollisionData<true>(collisions, selectedLcCandidatesMl);
  }
  PROCESS_SWITCH(HfTaskLc, processDataWithMl, "Process real data with the ML method and without centrality", false);

  void processDataStdWithFT0C(CollisionsWithFT0C const& collisions,
                              LcCandidates const& selectedLcCandidates,
                              aod::Tracks const&)
  {
    runAnalysisPerCollisionData<false>(collisions, selectedLcCandidates);
  }
  PROCESS_SWITCH(HfTaskLc, processDataStdWithFT0C, "Process real data with the standard method and with FT0C centrality", false);

  void processDataWithMlWithFT0C(CollisionsWithFT0C const& collisions,
                                 LcCandidatesMl const& selectedLcCandidatesMl,
                                 aod::Tracks const&)
  {
    runAnalysisPerCollisionData<true>(collisions, selectedLcCandidatesMl);
  }
  PROCESS_SWITCH(HfTaskLc, processDataWithMlWithFT0C, "Process real data with the ML method and with FT0C centrality", false);

  void processDataStdWithFT0M(CollisionsWithFT0M const& collisions,
                              LcCandidates const& selectedLcCandidates,
                              aod::Tracks const&)
  {
    runAnalysisPerCollisionData<false>(collisions, selectedLcCandidates);
  }
  PROCESS_SWITCH(HfTaskLc, processDataStdWithFT0M, "Process real data with the standard method and with FT0M centrality", false);

  void processDataWithMlWithFT0M(CollisionsWithFT0M const& collisions,
                                 LcCandidatesMl const& selectedLcCandidatesMl,
                                 aod::Tracks const&)
  {
    runAnalysisPerCollisionData<true>(collisions, selectedLcCandidatesMl);
  }
  PROCESS_SWITCH(HfTaskLc, processDataWithMlWithFT0M, "Process real data with the ML method and with FT0M centrality", false);

  void processMcStd(CollisionsMc const& collisions,
                    LcCandidatesMc const& selectedLcCandidatesMc,
                    McParticles3ProngMatched const& mcParticles,
                    aod::McCollisions const&,
                    aod::TracksWMc const&)
  {
    runAnalysisPerCollisionMc<false>(collisions, selectedLcCandidatesMc, mcParticles);
  }
  PROCESS_SWITCH(HfTaskLc, processMcStd, "Process MC with the standard method", false);

  void processMcWithMl(CollisionsMc const& collisions,
                       LcCandidatesMlMc const& selectedLcCandidatesMlMc,
                       McParticles3ProngMatched const& mcParticles,
                       aod::McCollisions const&,
                       aod::TracksWMc const&)
  {
    runAnalysisPerCollisionMc<true>(collisions, selectedLcCandidatesMlMc, mcParticles);
  }
  PROCESS_SWITCH(HfTaskLc, processMcWithMl, "Process Mc with the ML method and without centrality", false);

  void processMcStdWithFT0C(CollisionsMcWithFT0C const& collisions,
                            LcCandidatesMc const& selectedLcCandidatesMc,
                            McParticles3ProngMatched const& mcParticles,
                            aod::McCollisions const&,
                            aod::TracksWMc const&)
  {
    runAnalysisPerCollisionMc<false>(collisions, selectedLcCandidatesMc, mcParticles);
  }
  PROCESS_SWITCH(HfTaskLc, processMcStdWithFT0C, "Process MC with the standard method with FT0C centrality", false);

  void processMcWithMlWithFT0C(CollisionsMcWithFT0C const& collisions,
                               LcCandidatesMlMc const& selectedLcCandidatesMlMc,
                               McParticles3ProngMatched const& mcParticles,
                               aod::McCollisions const&,
                               aod::TracksWMc const&)
  {
    runAnalysisPerCollisionMc<true>(collisions, selectedLcCandidatesMlMc, mcParticles);
  }
  PROCESS_SWITCH(HfTaskLc, processMcWithMlWithFT0C, "Process Mc with the ML method with FT0C centrality", false);

  void processMcStdWithFT0M(CollisionsMcWithFT0M const& collisions,
                            LcCandidatesMc const& selectedLcCandidatesMc,
                            McParticles3ProngMatched const& mcParticles,
                            aod::McCollisions const&,
                            aod::TracksWMc const&)
  {
    runAnalysisPerCollisionMc<false>(collisions, selectedLcCandidatesMc, mcParticles);
  }
  PROCESS_SWITCH(HfTaskLc, processMcStdWithFT0M, "Process MC with the standard method with FT0M centrality", false);

  void processMcWithMlWithFT0M(CollisionsMcWithFT0M const& collisions,
                               LcCandidatesMlMc const& selectedLcCandidatesMlMc,
                               McParticles3ProngMatched const& mcParticles,
                               aod::McCollisions const&,
                               aod::TracksWMc const&)
  {
    runAnalysisPerCollisionMc<true>(collisions, selectedLcCandidatesMlMc, mcParticles);
  }
  PROCESS_SWITCH(HfTaskLc, processMcWithMlWithFT0M, "Process Mc with the ML method with FT0M centrality", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskLc>(cfgc)};
}
