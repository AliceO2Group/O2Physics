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
/// \author Himanshu Sharma <himanshu.sharma@cern.ch>, University and INFN Padova
/// \author Cristina Terrevoli <cristina.terrevoli@cern.ch>, INFN Bari

#include "PWGHF/Core/DecayChannels.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/AliasTables.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/DataModel/TrackIndexSkimmingTables.h"

#include "Common/Core/RecoDecay.h"

#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/runDataProcessing.h>

#include <THnSparse.h>
#include <TPDGCode.h>

#include <array>
#include <numeric>
#include <vector>

using namespace o2;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::framework::expressions;

/// Ξc± analysis task

struct HfTaskXic {
  Configurable<int> selectionFlagXic{"selectionFlagXic", 1, "Selection Flag for Xic"};
  Configurable<double> yCandGenMax{"yCandGenMax", 0.5, "max. gen particle rapidity"};
  Configurable<double> yCandRecoMax{"yCandRecoMax", 0.8, "max. cand. rapidity"};
  Configurable<float> etaTrackMax{"etaTrackMax", 4.0, "max. track eta"};
  Configurable<float> dcaXYTrackMax{"dcaXYTrackMax", 0.0025, "max. DCAxy for track"};
  Configurable<float> dcaZTrackMax{"dcaZTrackMax", 0.0025, "max. DCAz for track"};
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_xic_to_p_k_pi::vecBinsPt}, "pT bin limits"};

  Configurable<bool> enableTHn{"enableTHn", false, "enable THn for Xic"};
  HfHelper hfHelper;
  Service<o2::framework::O2DatabasePDG> pdg;

  Filter filterSelectCandidates = (aod::hf_sel_candidate_xic::isSelXicToPKPi >= selectionFlagXic || aod::hf_sel_candidate_xic::isSelXicToPiKP >= selectionFlagXic);

  // THnSparse for ML outputScores and Vars
  ConfigurableAxis thnConfigAxisPt{"thnConfigAxisPt", {36, 0, 36}, ""};
  ConfigurableAxis thnConfigAxisMass{"thnConfigAxisMass", {300, 1.98, 2.58}, ""};
  ConfigurableAxis thnConfigAxisPtProng{"thnConfigAxisPtProng", {100, 0, 20}, ""};
  ConfigurableAxis thnConfigAxisChi2PCA{"thnConfigAxisChi2PCA", {100, 0, 20}, ""};
  ConfigurableAxis thnConfigAxisDecLength{"thnConfigAxisDecLength", {10, 0, 0.05}, ""};
  ConfigurableAxis thnConfigAxisDecLengthXY{"thnConfigAxisDecLengthXY", {10, 0, 0.05}, ""};
  ConfigurableAxis thnConfigAxisCPA{"thnConfigAxisCPA", {20, 0.8, 1}, ""};
  ConfigurableAxis thnConfigAxisBdtScoreBkg{"thnConfigAxisBdtScoreBkg", {100, 0., 1.}, ""};
  ConfigurableAxis thnConfigAxisBdtScoreSignal{"thnConfigAxisBdtScoreSignal", {100, 0., 1.}, ""};
  ConfigurableAxis thnConfigAxisYMC{"thnConfigAxisYMC", {100, -2., 2.}, ""};
  //

  float etaMaxAcceptance = 0.8;
  float ptMinAcceptance = 0.1;

  HistogramRegistry registry{
    "registry", // histo not in pt bins
    {
      {"Data/hPt", "3-prong candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1D, {{360, 0., 36.}}}},    // pt Xic
      {"Data/hEta", "3-prong candidates;candidate #it{eta};entries", {HistType::kTH1D, {{100, -2., 2.}}}},                  // eta Xic
      {"Data/hPhi", "3-prong candidates;candidate #varphi;entries", {HistType::kTH1D, {{72, 0., constants::math::TwoPI}}}}, // phi Xic
      {"Data/hMass", "3-prong candidates; inv. mass (p K #pi) (GeV/#it{c}^{2})", {HistType::kTH1D, {{600, 1.98, 2.58}}}},   // mass Xic
      {"Data/hMultiplicity", "multiplicity;multiplicity;entries", {HistType::kTH1D, {{1000, 0., 1000.}}}},
      {"MC/generated/signal/hPtGenSig", "3-prong candidates (matched);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1D, {{360, 0., 36.}}}},
      {"MC/generated/signal/hPtGen", "MC particles (matched);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1D, {{360, 0., 36.}}}},
      {"MC/generated/prompt/hPtGenPrompt", "MC particles (matched, prompt);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1D, {{360, 0., 36.}}}},
      {"MC/generated/nonprompt/hPtGenNonPrompt", "MC particles (matched, non-prompt);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1D, {{360, 0., 36.}}}},
      {"MC/generated/signal/hEtaGen", "MC particles (matched);#it{#eta};entries", {HistType::kTH1D, {{100, -2., 2.}}}},
      {"MC/generated/prompt/hEtaGenPrompt", "MC particles (matched, prompt);#it{#eta};entries", {HistType::kTH1D, {{100, -2., 2.}}}},
      {"MC/generated/nonprompt/hEtaGenNonPrompt", "MC particles (matched, non-prompt);#it{#eta};entries", {HistType::kTH1D, {{100, -2., 2.}}}},
      {"MC/generated/signal/hYGen", "MC particles (matched);#it{y};entries", {HistType::kTH1D, {{100, -2., 2.}}}},
      {"MC/generated/prompt/hYGenPrompt", "MC particles (matched, prompt);#it{y};entries", {HistType::kTH1D, {{100, -2., 2.}}}},
      {"MC/generated/nonprompt/hYGenNonPrompt", "MC particles (matched, non-prompt);#it{y};entries", {HistType::kTH1D, {{100, -2., 2.}}}},
      {"MC/reconstructed/signal/hMassRecSig", "3-prong candidates (matched);inv. mass (p K #pi) (GeV/#it{c}^{2})", {HistType::kTH1D, {{600, 2.18, 2.58}}}},
      {"MC/reconstructed/prompt/hMassRecSigPrompt", "3-prong candidates (matched, prompt);inv. mass (p K #pi) (GeV/#it{c}^{2})", {HistType::kTH1D, {{600, 2.18, 2.58}}}},
      {"MC/reconstructed/nonprompt/hMassRecSigNonPrompt", "3-prong candidates (matched, non-prompt);inv. mass (p K #pi) (GeV/#it{c}^{2})", {HistType::kTH1D, {{600, 2.18, 2.58}}}},
      {"MC/reconstructed/signal/hPtRecSig", "3-prong candidates (matched);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1D, {{360, 0., 36.}}}},
      {"MC/reconstructed/prompt/hPtRecSigPrompt", "3-prong candidates (matched, prompt);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1D, {{360, 0., 36.}}}},
      {"MC/reconstructed/nonprompt/hPtRecSigNonPrompt", "3-prong candidates (matched, non-prompt);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1D, {{360, 0., 36.}}}},
      {"MC/reconstructed/background/hPtRecBg", "3-prong candidates (unmatched);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1D, {{100, 0., 10.}}}}}};

  void init(InitContext&)
  {
    std::array<bool, 5> doprocess{doprocessDataStd, doprocessDataWithMl, doprocessMcStd, doprocessMcWithMl};
    if ((std::accumulate(doprocess.begin(), doprocess.end(), 0)) != 1) {
      LOGP(fatal, "no or more than one process function enabled! Please check your configuration!");
    }

    AxisSpec const axisPPid = {100, 0.f, 10.0f, "#it{p} (GeV/#it{c})"};
    AxisSpec const axisNSigmaPr = {100, -6.f, 6.f, "n#it{#sigma}_{p}"};
    AxisSpec const axisNSigmaPi = {100, -6.f, 6.f, "n#it{#sigma}_{#pi}"};
    AxisSpec const axisNSigmaKa = {100, -6.f, 6.f, "n#it{#sigma}_{K}"};

    auto vbins = (std::vector<double>)binsPt; // histo in pt bins
    registry.add("Data/hMassVsPt", "3-prong candidates;inv. mass (p K #pi) (GeV/#it{c}^{2});;entries", {HistType::kTH2F, {{500, 2., 3.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/hDecLength", "3-prong candidates;decay length (cm);;entries", {HistType::kTH2F, {{200, 0., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/hDecLengthXY", "3-prong candidates;decay length XY (cm);;entries", {HistType::kTH2F, {{200, 0., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/hNormalisedDecayLengthXY", "3-prong candidates;norm. decay length XY (cm);;entries", {HistType::kTH2F, {{80, 0., 80.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/hImpParXY", "3-prong candidates;candidate DCAxy to prim. vertex (cm);;entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/hd0Prong0", "3-prong candidates;prong 0 DCAxy to prim. vertex (cm);;entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/hd0Prong1", "3-prong candidates;prong 1 DCAxy to prim. vertex (cm);;entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/hd0Prong2", "3-prong candidates;prong 1 DCAxy to prim. vertex (cm);;entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/hCt", "3-prong candidates;proper lifetime (#Xi_{c}) * #it{c} (cm);;entries", {HistType::kTH2F, {{120, -20., 100.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/hCPA", "3-prong candidates;cosine of pointing angle;;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/hCPAXY", "3-prong candidates;cosine of pointing angle XY;;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/hEtaVsPt", "3-prong candidates;candidate #it{#eta};;entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/hSelectionStatus", "3-prong candidates;selection status;;entries", {HistType::kTH2F, {{5, -0.5, 4.5}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/hImpParErr0", "3-prong candidates;impact parameter error prong 0 (cm);;entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/hImpParErr1", "3-prong candidates;impact parameter error prong 1 (cm);;entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/hImpParErr2", "3-prong candidates;impact parameter error prong 2 (cm);;entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/hDecLenErr", "3-prong candidates;decay length error (cm);;entries", {HistType::kTH2F, {{100, 0., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/hPtProng0", "3-prong candidates;prong 0 #it{p}_{T} (GeV/#it{c});;entries", {HistType::kTH2F, {{100, 0., 10.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/hPtProng1", "3-prong candidates;prong 1 #it{p}_{T} (GeV/#it{c});;entries", {HistType::kTH2F, {{100, 0., 10.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/hPtProng2", "3-prong candidates;prong 2 #it{p}_{T} (GeV/#it{c});;entries", {HistType::kTH2F, {{100, 0., 10.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("Data/hChi2PCA", "3-prong candidates;prong Chi2PCA to sec. vertex (cm);;entries", {HistType::kTH2F, {{100, 0, 120}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

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

    // MC generated
    registry.add("MC/generated/hYGenWithProngsInAcceptance", "MC particles (generated-daughters in acceptance);candidate #it{y}^{gen};entries", {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/generated/hEtaGenWithProngsInAcceptance", "MC particles (generated-daughters in acceptance);candidate #it{#eta}^{gen};entries", {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});

    // MC reconstructed:signal
    registry.add("MC/reconstructed/signal/hMassSig", "Invariant mass (matched);m (p K #pi) (GeV/#it{c}^{2});;entries", {HistType::kTH2F, {{500, 1.6, 3.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hDecLengthRecSig", "3-prong candidates;decay length (cm);;entries", {HistType::kTH2F, {{200, 0., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hDecLengthXYRecSig", "3-prong candidates;decay length XY (cm);;entries", {HistType::kTH2F, {{200, 0., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hNormalisedDecayLengthXYRecSig", "3-prong candidates;norm. decay length XY (cm);;entries", {HistType::kTH2F, {{80, 0., 80.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hImpParXYRecSig", "3-prong candidates;candidate DCAxy to prim. vertex (cm);;entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hd0Prong0RecSig", "3-prong candidates;prong 0 DCAxy to prim. vertex (cm);;entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hd0Prong1RecSig", "3-prong candidates;prong 1 DCAxy to prim. vertex (cm);;entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hd0Prong2RecSig", "3-prong candidates;prong 1 DCAxy to prim. vertex (cm);;entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hCtRecSig", "3-prong candidates;proper lifetime (#Xi_{c}) * #it{c} (cm);;entries", {HistType::kTH2F, {{120, -20., 100.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hCPARecSig", "3-prong candidates;cosine of pointing angle;;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hCPAXYRecSig", "3-prong candidates;cosine of pointing angle XY;;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hEtaRecSig", "3-prong candidates;candidate #it{#eta};;entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hImpParErrSig", "3-prong candidates;impact parameter error (cm);;entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hDecLenErrSig", "3-prong candidates;decay length error (cm);;entries", {HistType::kTH2F, {{100, 0., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hPtProng0RecSig", "3-prong candidates;prong 0 #it{p}_{T} (GeV/#it{c});;entries", {HistType::kTH2F, {{100, 0., 10.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hPtProng1RecSig", "3-prong candidates;prong 1 #it{p}_{T} (GeV/#it{c});;entries", {HistType::kTH2F, {{100, 0., 10.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hPtProng2RecSig", "3-prong candidates;prong 2 #it{p}_{T} (GeV/#it{c});;entries", {HistType::kTH2F, {{100, 0., 10.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/signal/hChi2PCARecSig", "3-prong candidates;prong Chi2PCA to sec.  vertex (cm);; entries", {HistType::kTH2F, {{100, 0, 120}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    registry.add("MC/reconstructed/signal/hMassVsPtRecSig", "3-prong candidates (matched);inv. mass (p K #pi) (GeV/#it{c}^{2}); p_{T}", {HistType::kTH2F, {{600, 2.18, 2.58}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/prompt/hMassVsPtRecSigPrompt", "3-prong candidates (matched, prompt);inv. mass (p K #pi) (GeV/#it{c}^{2}); p_{T}", {HistType::kTH2F, {{600, 2.18, 2.58}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/nonprompt/hMassVsPtRecSigNonPrompt", "3-prong candidates (matched, non-prompt);inv. mass (p K #pi) (GeV/#it{c}^{2}); p_{T}", {HistType::kTH2F, {{600, 2.18, 2.58}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    registry.add("MC/reconstructed/signal/hEtaVsPtRecSig", "3-prong candidates (matched);candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/prompt/hEtaVsPtRecSigPrompt", "3-prong candidates (matched, prompt);candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/nonprompt/hEtaVsPtRecSigNonPrompt", "3-prong candidates (matched, non-prompt);candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    registry.add("MC/generated/signal/hEtaVsPtGenSig", "3-prong candidates (matched);candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/generated/prompt/hEtaVsPtGenSigPrompt", "3-prong candidates (matched, prompt);candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/generated/nonprompt/hEtaVsPtGenSigNonPrompt", "3-prong candidates (matched, non-prompt);candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    registry.add("MC/generated/signal/hYVsPtGenSig", "3-prong candidates (matched);candidate #it{y};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/generated/prompt/hYVsPtGenSigPrompt", "3-prong candidates (matched, prompt);candidate #it{y};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/generated/nonprompt/hYVsPtGenSigNonPrompt", "3-prong candidates (matched, non-prompt);candidate #it{y};entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    // background
    registry.add("MC/reconstructed/background/hMassBg", "Invariant mass (unmatched);m (p K #pi) (GeV/#it{c}^{2});;entries", {HistType::kTH2F, {{500, 1.6, 3.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/background/hDecLengthRecBg", "3-prong candidates;decay length (cm);;entries", {HistType::kTH2F, {{200, 0., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/background/hDecLengthXYRecBg", "3-prong candidates;decay length xy(cm);;entries", {HistType::kTH2F, {{200, 0., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/background/hNormalisedDecayLengthXYRecBg", "3-prong candidates;norm. decay length XY (cm);;entries", {HistType::kTH2F, {{80, 0., 80.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/background/hImpParXYRecBg", "3-prong candidates;candidate DCAxy to prim. vertex (cm);;entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/background/hd0Prong0RecBg", "3-prong candidates;prong 0 DCAxy to prim. vertex (cm);;entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/background/hd0Prong1RecBg", "3-prong candidates;prong 1 DCAxy to prim. vertex (cm);;entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/background/hd0Prong2RecBg", "3-prong candidates;prong 1 DCAxy to prim. vertex (cm);;entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/background/hCtRecBg", "3-prong candidates;proper lifetime (#Xi_{c}) * #it{c} (cm);;entries", {HistType::kTH2F, {{120, -20., 100.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/background/hCPARecBg", "3-prong candidates;cosine of pointing angle;;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/background/hCPAXYRecBg", "3-prong candidates;cosine of pointing angle XY;;entries", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/background/hEtaRecBg", "3-prong candidates;candidate #it{#eta};;entries", {HistType::kTH2F, {{100, -2., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/background/hImpParErrBg", "3-prong candidates;impact parameter error (cm);;entries", {HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/background/hDecLenErrBg", "3-prong candidates;decay length error (cm);;entries", {HistType::kTH2F, {{100, 0., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/background/hPtProng0RecBg", "3-prong candidates;prong 0 #it{p}_{T} (GeV/#it{c});;entries", {HistType::kTH2F, {{100, 0., 10.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/background/hPtProng1RecBg", "3-prong candidates;prong 1 #it{p}_{T} (GeV/#it{c});;entries", {HistType::kTH2F, {{100, 0., 10.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/background/hPtProng2RecBg", "3-prong candidates;prong 2 #it{p}_{T} (GeV/#it{c});;entries", {HistType::kTH2F, {{100, 0., 10.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/background/hChi2PCARecBg", "3-prong candidates;prong    Chi2PCA to sec.  vertex (cm);; entries", {HistType::kTH2F, {{100, 0, 120}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    // THn for candidates Xic+ cut variation
    if (enableTHn) {
      const AxisSpec thnAxisMass{thnConfigAxisMass, "inv. mass (p K #pi) (GeV/#it{c}^{2})"};
      const AxisSpec thnAxisPt{thnConfigAxisPt, "#it{p}_{T}(#Xi_{c}^{+}) (GeV/#it{c})"};
      const AxisSpec thnAxisPtMC{thnConfigAxisPt, "#it{p}_{T}(#Xi_{c}^{+} MC) (GeV/#it{c})"};
      const AxisSpec thnAxisYMC{thnConfigAxisYMC, "#it{y}_{MC}(#Xi_{c}^{+} MC)"};
      const AxisSpec thnAxisChi2PCA{thnConfigAxisChi2PCA, "Chi2PCA to sec. vertex (cm)"};
      const AxisSpec thnAxisDecLength{thnConfigAxisDecLength, "decay length (cm)"};
      const AxisSpec thnAxisDecLengthXY{thnConfigAxisDecLengthXY, "decay length (cm)"};
      const AxisSpec thnAxisCPA{thnConfigAxisCPA, "cosine of pointing angle"};
      const AxisSpec thnAxisBdtScoreXicBkg{thnConfigAxisBdtScoreBkg, "BDT bkg score (Xic)"};
      const AxisSpec thnAxisBdtScoreXicPrompt{thnConfigAxisBdtScoreSignal, "BDT prompt score (Xic)"};
      const AxisSpec thnAxisBdtScoreXicNonPrompt{thnConfigAxisBdtScoreSignal, "BDT non-prompt score (Xic)"};
      const AxisSpec thnAxisMcOrigin{3, -0.5, 2.5, "MC origin"};
      const AxisSpec thnAxisMCAllProngAccepted{2, -0.5, 1.5, "All MC prongs accepted"};

      if (doprocessDataWithMl || doprocessMcWithMl) { // with ML
        registry.add("hnXicVarsWithBdt", "THn for Xic candidates with BDT scores", HistType::kTHnSparseF, {thnAxisMass, thnAxisPt, thnAxisBdtScoreXicBkg, thnAxisBdtScoreXicPrompt, thnAxisBdtScoreXicNonPrompt, thnAxisMcOrigin});
      } else {
        registry.add("hnXicVars", "THn for Xic candidates", HistType::kTHnSparseF, {thnAxisMass, thnAxisPt, thnAxisChi2PCA, thnAxisDecLength, thnAxisDecLengthXY, thnAxisCPA, thnAxisMcOrigin});
      }
    }
  } // end init

  /// Selection of Xic daughters in geometrical acceptanc
  /// \param etaProng is the pseudorapidity of Xic prong
  /// \param ptProng is the pT of Xic prong
  /// \return true if prong is in geometrical acceptance
  template <typename T = float>
  bool isProngInAcceptance(const T& etaProng, const T& ptProng)
  {
    return std::abs(etaProng) <= etaMaxAcceptance && ptProng >= ptMinAcceptance;
  }
  template <bool UseMl, typename Cands>
  void analysisData(aod::Collision const& collision,
                    Cands const& candidates,
                    aod::TracksWDca const& tracks)
  {
    int nTracks = 0;

    if (collision.numContrib() > 1) {
      for (const auto& track : tracks) {
        if (std::abs(track.eta()) > etaTrackMax) {
          continue;
        }
        if (std::abs(track.dcaXY()) > dcaXYTrackMax || std::abs(track.dcaZ()) > dcaZTrackMax) {
          continue;
        }
        nTracks++;
      }
    }

    registry.fill(HIST("Data/hMultiplicity"), nTracks); // filling the histo for multiplicity

    for (const auto& candidate : candidates) {
      if (!(candidate.hfflag() & 1 << aod::hf_cand_3prong::DecayType::XicToPKPi)) {
        continue;
      }
      if (yCandRecoMax >= 0. && std::abs(hfHelper.yXic(candidate)) > yCandRecoMax) {
        continue;
      }
      auto ptCandidate = candidate.pt();

      if (candidate.isSelXicToPKPi() >= selectionFlagXic) { // pKpi
        registry.fill(HIST("Data/hMassVsPt"), hfHelper.invMassXicToPKPi(candidate), ptCandidate);
        registry.fill(HIST("Data/hMass"), hfHelper.invMassXicToPKPi(candidate));
      }
      if (candidate.isSelXicToPiKP() >= selectionFlagXic) { // piKp
        registry.fill(HIST("Data/hMassVsPt"), hfHelper.invMassXicToPiKP(candidate), ptCandidate);
        registry.fill(HIST("Data/hMass"), hfHelper.invMassXicToPiKP(candidate));
      }
      registry.fill(HIST("Data/hPt"), ptCandidate);
      registry.fill(HIST("Data/hEta"), candidate.eta());
      registry.fill(HIST("Data/hPhi"), candidate.phi());
      registry.fill(HIST("Data/hPtProng0"), candidate.ptProng0(), ptCandidate);
      registry.fill(HIST("Data/hPtProng1"), candidate.ptProng1(), ptCandidate);
      registry.fill(HIST("Data/hPtProng2"), candidate.ptProng2(), ptCandidate);
      registry.fill(HIST("Data/hDecLength"), candidate.decayLength(), ptCandidate);
      registry.fill(HIST("Data/hDecLengthXY"), candidate.decayLengthXY(), ptCandidate);
      registry.fill(HIST("Data/hNormalisedDecayLengthXY"), candidate.decayLengthXYNormalised(), ptCandidate);
      registry.fill(HIST("Data/hImpParXY"), candidate.impactParameterXY(), ptCandidate);
      registry.fill(HIST("Data/hd0Prong0"), candidate.impactParameter0(), ptCandidate);
      registry.fill(HIST("Data/hd0Prong1"), candidate.impactParameter1(), ptCandidate);
      registry.fill(HIST("Data/hd0Prong2"), candidate.impactParameter2(), ptCandidate);
      registry.fill(HIST("Data/hCt"), hfHelper.ctXic(candidate), ptCandidate);
      registry.fill(HIST("Data/hCPA"), candidate.cpa(), ptCandidate);
      registry.fill(HIST("Data/hCPAXY"), candidate.cpaXY(), ptCandidate);
      registry.fill(HIST("Data/hEtaVsPt"), candidate.eta(), ptCandidate);
      registry.fill(HIST("Data/hSelectionStatus"), candidate.isSelXicToPKPi(), ptCandidate);
      registry.fill(HIST("Data/hImpParErr0"), candidate.errorImpactParameter0(), ptCandidate);
      registry.fill(HIST("Data/hImpParErr1"), candidate.errorImpactParameter1(), ptCandidate);
      registry.fill(HIST("Data/hImpParErr2"), candidate.errorImpactParameter2(), ptCandidate);
      registry.fill(HIST("Data/hDecLenErr"), candidate.errorDecayLength(), ptCandidate);
      registry.fill(HIST("Data/hChi2PCA"), candidate.chi2PCA(), ptCandidate);

      // PID histos
      auto trackProng0 = candidate.template prong0_as<aod::TracksWDca>();
      auto trackProng1 = candidate.template prong1_as<aod::TracksWDca>();
      auto trackProng2 = candidate.template prong2_as<aod::TracksWDca>();

      auto momentumProng0 = trackProng0.p();
      auto momentumProng1 = trackProng1.p();
      auto momentumProng2 = trackProng2.p();

      // TPC nSigma histograms
      registry.fill(HIST("Data/hPVsTPCNSigmaPr_Prong0"), momentumProng0, candidate.nSigTpcPr0());
      registry.fill(HIST("Data/hPVsTPCNSigmaPi_Prong0"), momentumProng0, candidate.nSigTpcPi0());
      registry.fill(HIST("Data/hPVsTPCNSigmaKa_Prong0"), momentumProng0, candidate.nSigTpcKa0());

      registry.fill(HIST("Data/hPVsTPCNSigmaPr_Prong1"), momentumProng1, candidate.nSigTpcPr1());
      registry.fill(HIST("Data/hPVsTPCNSigmaPi_Prong1"), momentumProng1, candidate.nSigTpcPi1());
      registry.fill(HIST("Data/hPVsTPCNSigmaKa_Prong1"), momentumProng1, candidate.nSigTpcKa1());

      registry.fill(HIST("Data/hPVsTPCNSigmaPr_Prong2"), momentumProng2, candidate.nSigTpcPr2());
      registry.fill(HIST("Data/hPVsTPCNSigmaPi_Prong2"), momentumProng2, candidate.nSigTpcPi2());
      registry.fill(HIST("Data/hPVsTPCNSigmaKa_Prong2"), momentumProng2, candidate.nSigTpcKa2());

      // TOF nSigma histograms
      registry.fill(HIST("Data/hPVsTOFNSigmaPr_Prong0"), momentumProng0, candidate.nSigTofPr0());
      registry.fill(HIST("Data/hPVsTOFNSigmaPi_Prong0"), momentumProng0, candidate.nSigTofPi0());
      registry.fill(HIST("Data/hPVsTOFNSigmaKa_Prong0"), momentumProng0, candidate.nSigTofKa0());

      registry.fill(HIST("Data/hPVsTOFNSigmaPr_Prong1"), momentumProng1, candidate.nSigTofPr1());
      registry.fill(HIST("Data/hPVsTOFNSigmaPi_Prong1"), momentumProng1, candidate.nSigTofPi1());
      registry.fill(HIST("Data/hPVsTOFNSigmaKa_Prong1"), momentumProng1, candidate.nSigTofKa1());

      registry.fill(HIST("Data/hPVsTOFNSigmaPr_Prong2"), momentumProng2, candidate.nSigTofPr2());
      registry.fill(HIST("Data/hPVsTOFNSigmaPi_Prong2"), momentumProng2, candidate.nSigTofPi2());
      registry.fill(HIST("Data/hPVsTOFNSigmaKa_Prong2"), momentumProng2, candidate.nSigTofKa2());

      // THnSparse
      if (enableTHn) {
        double massXic(-1);
        double outputBkg(-1), outputPrompt(-1), outputFD(-1);
        const int ternaryCl = 3;
        if (candidate.isSelXicToPKPi() >= selectionFlagXic) {
          massXic = hfHelper.invMassXicToPKPi(candidate);
          if constexpr (UseMl) {
            if (candidate.mlProbXicToPKPi().size() == ternaryCl) {
              outputBkg = candidate.mlProbXicToPKPi()[0];    /// bkg score
              outputPrompt = candidate.mlProbXicToPKPi()[1]; /// prompt score
              outputFD = candidate.mlProbXicToPKPi()[2];     /// non-prompt score
            }
            /// Fill the ML outputScores and variables of candidate Xic
            registry.get<THnSparse>(HIST("hnXicVarsWithBdt"))->Fill(massXic, ptCandidate, outputBkg, outputPrompt, outputFD, false);
          } else {
            registry.get<THnSparse>(HIST("hnXicVars"))->Fill(massXic, ptCandidate, candidate.chi2PCA(), candidate.decayLength(), candidate.decayLengthXY(), candidate.cpa(), false);
          }
        }
        if (candidate.isSelXicToPiKP() >= selectionFlagXic) {
          massXic = hfHelper.invMassXicToPiKP(candidate);
          if constexpr (UseMl) {
            if (candidate.mlProbXicToPiKP().size() == ternaryCl) {
              outputBkg = candidate.mlProbXicToPiKP()[0];    /// bkg score
              outputPrompt = candidate.mlProbXicToPiKP()[1]; /// prompt score
              outputFD = candidate.mlProbXicToPiKP()[2];     /// non-prompt score
            }
            /// Fill the ML outputScores and variables of candidate
            registry.get<THnSparse>(HIST("hnXicVarsWithBdt"))->Fill(massXic, ptCandidate, outputBkg, outputPrompt, outputFD, false);
          } else {
            registry.get<THnSparse>(HIST("hnXicVars"))->Fill(massXic, ptCandidate, candidate.chi2PCA(), candidate.decayLength(), candidate.decayLengthXY(), candidate.cpa(), false);
          }
        }
      } // thn for Xic
    } // loop candidates
  } // end process data

  void processDataStd(aod::Collision const& collision,
                      soa::Filtered<soa::Join<aod::HfCand3ProngWPidPiKaPr, aod::HfSelXicToPKPi>> const& candidates,
                      aod::TracksWDca const& tracks)
  {
    analysisData<false>(collision, candidates, tracks);
  }
  PROCESS_SWITCH(HfTaskXic, processDataStd, "Process Data with the standard method", true);

  void processDataWithMl(aod::Collision const& collision,
                         soa::Filtered<soa::Join<aod::HfCand3ProngWPidPiKaPr, aod::HfSelXicToPKPi, aod::HfMlXicToPKPi>> const& candidatesMl, aod::TracksWDca const& tracks)
  {
    analysisData<true>(collision, candidatesMl, tracks);
  }
  PROCESS_SWITCH(HfTaskXic, processDataWithMl, "Process Data with the ML method", false);

  // Fill MC histograms
  template <bool UseMl, typename Cands>
  void analysisMc(Cands const& candidates,
                  soa::Join<aod::McParticles, aod::HfCand3ProngMcGen> const& mcParticles,
                  aod::TracksWMc const&)
  {
    // MC rec.
    for (const auto& candidate : candidates) {
      auto ptCandidate = candidate.pt();
      // Selected Xic
      if (!(candidate.hfflag() & 1 << aod::hf_cand_3prong::DecayType::XicToPKPi)) {
        continue;
      }
      if (yCandRecoMax >= 0. && std::abs(hfHelper.yXic(candidate)) > yCandRecoMax) {
        continue;
      }

      auto massXicToPKPi = 0.;
      auto massXicToPiKP = 0.;

      if (candidate.isSelXicToPKPi() >= selectionFlagXic) {
        massXicToPKPi = hfHelper.invMassXicToPKPi(candidate);
      }
      if (candidate.isSelXicToPiKP() >= selectionFlagXic) {
        massXicToPiKP = hfHelper.invMassXicToPiKP(candidate); // mass conjugate
      }

      if (std::abs(candidate.flagMcMatchRec()) == hf_decay::hf_cand_3prong::DecayChannelMain::XicToPKPi) {
        // Get the corresponding MC particle.
        auto mcParticleProng0 = candidate.template prong0_as<aod::TracksWMc>().template mcParticle_as<soa::Join<aod::McParticles, aod::HfCand3ProngMcGen>>();
        auto pdgCodeProng0 = std::abs(mcParticleProng0.pdgCode());
        // Signal
        registry.fill(HIST("MC/reconstructed/signal/hPtRecSig"), ptCandidate); // rec. level pT

        if (massXicToPKPi != 0.) {
          registry.fill(HIST("MC/reconstructed/signal/hMassVsPtRecSig"), massXicToPKPi, ptCandidate);
          registry.fill(HIST("MC/reconstructed/signal/hMassRecSig"), massXicToPKPi);
        }
        if (massXicToPiKP != 0.) {
          registry.fill(HIST("MC/reconstructed/signal/hMassVsPtRecSig"), massXicToPiKP, ptCandidate);
          registry.fill(HIST("MC/reconstructed/signal/hMassRecSig"), massXicToPiKP);
        }
        registry.fill(HIST("MC/reconstructed/signal/hDecLengthRecSig"), candidate.decayLength(), ptCandidate);
        registry.fill(HIST("MC/reconstructed/signal/hDecLengthXYRecSig"), candidate.decayLengthXY(), ptCandidate);
        registry.fill(HIST("MC/reconstructed/signal/hNormalisedDecayLengthXYRecSig"), candidate.decayLengthXYNormalised(), ptCandidate);
        registry.fill(HIST("MC/reconstructed/signal/hPtProng0RecSig"), candidate.ptProng0(), ptCandidate);
        registry.fill(HIST("MC/reconstructed/signal/hPtProng1RecSig"), candidate.ptProng1(), ptCandidate);
        registry.fill(HIST("MC/reconstructed/signal/hPtProng2RecSig"), candidate.ptProng2(), ptCandidate);
        registry.fill(HIST("MC/reconstructed/signal/hImpParXYRecSig"), candidate.impactParameterXY(), ptCandidate);
        registry.fill(HIST("MC/reconstructed/signal/hd0Prong0RecSig"), candidate.impactParameter0(), ptCandidate);
        registry.fill(HIST("MC/reconstructed/signal/hd0Prong1RecSig"), candidate.impactParameter1(), ptCandidate);
        registry.fill(HIST("MC/reconstructed/signal/hd0Prong2RecSig"), candidate.impactParameter2(), ptCandidate);
        registry.fill(HIST("MC/reconstructed/signal/hCtRecSig"), hfHelper.ctXic(candidate), ptCandidate);
        registry.fill(HIST("MC/reconstructed/signal/hCPARecSig"), candidate.cpa(), ptCandidate);
        registry.fill(HIST("MC/reconstructed/signal/hCPAXYRecSig"), candidate.cpaXY(), ptCandidate);
        registry.fill(HIST("MC/reconstructed/signal/hEtaRecSig"), candidate.eta(), ptCandidate);
        registry.fill(HIST("MC/reconstructed/signal/hImpParErrSig"), candidate.errorImpactParameter0(), ptCandidate);
        registry.fill(HIST("MC/reconstructed/signal/hDecLenErrSig"), candidate.errorDecayLength(), ptCandidate);
        registry.fill(HIST("MC/reconstructed/signal/hChi2PCARecSig"), candidate.chi2PCA(), ptCandidate);
        registry.fill(HIST("MC/reconstructed/signal/hEtaVsPtRecSig"), candidate.eta(), ptCandidate);

        /// reconstructed signal prompt
        int const origin = candidate.originMcRec();
        if (origin == RecoDecay::OriginType::Prompt) {
          if ((candidate.isSelXicToPKPi() >= selectionFlagXic) && pdgCodeProng0 == kProton) {
            registry.fill(HIST("MC/reconstructed/prompt/hMassRecSigPrompt"), massXicToPKPi);
            registry.fill(HIST("MC/reconstructed/prompt/hMassVsPtRecSigPrompt"), massXicToPKPi, ptCandidate);
          }
          if ((candidate.isSelXicToPiKP() >= selectionFlagXic) && pdgCodeProng0 == kPiPlus) {
            registry.fill(HIST("MC/reconstructed/prompt/hMassRecSigPrompt"), massXicToPiKP);
            registry.fill(HIST("MC/reconstructed/prompt/hMassVsPtRecSigPrompt"), massXicToPiKP, ptCandidate);
          }
          registry.fill(HIST("MC/reconstructed/prompt/hEtaVsPtRecSigPrompt"), candidate.eta(), ptCandidate);
          registry.fill(HIST("MC/reconstructed/prompt/hPtRecSigPrompt"), ptCandidate);
        } else {
          if ((candidate.isSelXicToPKPi() >= selectionFlagXic) && pdgCodeProng0 == kProton) {
            registry.fill(HIST("MC/reconstructed/nonprompt/hMassRecSigNonPrompt"), massXicToPKPi);
            registry.fill(HIST("MC/reconstructed/nonprompt/hMassVsPtRecSigNonPrompt"), massXicToPKPi, ptCandidate);
          }
          if ((candidate.isSelXicToPiKP() >= selectionFlagXic) && pdgCodeProng0 == kPiPlus) {
            registry.fill(HIST("MC/reconstructed/nonprompt/hMassRecSigNonPrompt"), massXicToPiKP);
            registry.fill(HIST("MC/reconstructed/nonprompt/hMassVsPtRecSigNonPrompt"), massXicToPiKP, ptCandidate);
          }
          registry.fill(HIST("MC/reconstructed/nonprompt/hEtaVsPtRecSigNonPrompt"), candidate.eta(), ptCandidate);
          registry.fill(HIST("MC/reconstructed/nonprompt/hPtRecSigNonPrompt"), ptCandidate);
        }

        if (enableTHn) {
          double outputBkg(-1), outputPrompt(-1), outputFD(-1);
          const int ternaryCl = 3;
          if ((candidate.isSelXicToPKPi() >= selectionFlagXic) && pdgCodeProng0 == kProton) {
            if constexpr (UseMl) {
              if (candidate.mlProbXicToPKPi().size() == ternaryCl) {
                outputBkg = candidate.mlProbXicToPKPi()[0];    /// bkg score
                outputPrompt = candidate.mlProbXicToPKPi()[1]; /// prompt score
                outputFD = candidate.mlProbXicToPKPi()[2];     /// non-prompt score
              }
              /// Fill the ML outputScores and variables of candidate (todo: add multiplicity)
              registry.get<THnSparse>(HIST("hnXicVarsWithBdt"))->Fill(massXicToPKPi, ptCandidate, outputBkg, outputPrompt, outputFD, origin);
            } else {
              registry.get<THnSparse>(HIST("hnXicVars"))->Fill(massXicToPKPi, ptCandidate, candidate.chi2PCA(), candidate.decayLength(), candidate.decayLengthXY(), candidate.cpa(), origin);
            }
          }
          if ((candidate.isSelXicToPiKP() >= selectionFlagXic) && pdgCodeProng0 == kPiPlus) {
            if constexpr (UseMl) {
              if (candidate.mlProbXicToPiKP().size() == ternaryCl) {
                outputBkg = candidate.mlProbXicToPiKP()[0];    /// bkg score
                outputPrompt = candidate.mlProbXicToPiKP()[1]; /// prompt score
                outputFD = candidate.mlProbXicToPiKP()[2];     /// non-prompt score
              }
              /// Fill the ML outputScores and variables of candidate (todo: add multiplicity)
              // add here the pT_Mother, y_Mother, level (reco, Gen, Gen + Acc)
              registry.get<THnSparse>(HIST("hnXicVarsWithBdt"))->Fill(massXicToPiKP, ptCandidate, outputBkg, outputPrompt, outputFD, origin);
            } else {
              registry.get<THnSparse>(HIST("hnXicVars"))->Fill(massXicToPiKP, ptCandidate, candidate.chi2PCA(), candidate.decayLength(), candidate.decayLengthXY(), candidate.cpa(), origin);
            }
          }
        } // enable THn
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
        registry.fill(HIST("MC/reconstructed/background/hDecLengthXYRecBg"), candidate.decayLengthXY(), ptCandidate);
        registry.fill(HIST("MC/reconstructed/background/hNormalisedDecayLengthXYRecBg"), candidate.decayLengthXYNormalised(), ptCandidate);
        registry.fill(HIST("MC/reconstructed/background/hPtProng0RecBg"), candidate.ptProng0(), ptCandidate);
        registry.fill(HIST("MC/reconstructed/background/hPtProng1RecBg"), candidate.ptProng1(), ptCandidate);
        registry.fill(HIST("MC/reconstructed/background/hPtProng2RecBg"), candidate.ptProng2(), ptCandidate);
        registry.fill(HIST("MC/reconstructed/background/hImpParXYRecBg"), candidate.impactParameterXY(), ptCandidate);
        registry.fill(HIST("MC/reconstructed/background/hd0Prong0RecBg"), candidate.impactParameter0(), ptCandidate);
        registry.fill(HIST("MC/reconstructed/background/hd0Prong1RecBg"), candidate.impactParameter1(), ptCandidate);
        registry.fill(HIST("MC/reconstructed/background/hd0Prong2RecBg"), candidate.impactParameter2(), ptCandidate);
        registry.fill(HIST("MC/reconstructed/background/hCtRecBg"), hfHelper.ctXic(candidate), ptCandidate);
        registry.fill(HIST("MC/reconstructed/background/hCPARecBg"), candidate.cpa(), ptCandidate);
        registry.fill(HIST("MC/reconstructed/background/hCPAXYRecBg"), candidate.cpaXY(), ptCandidate);
        registry.fill(HIST("MC/reconstructed/background/hEtaRecBg"), candidate.eta(), ptCandidate);
        registry.fill(HIST("MC/reconstructed/background/hImpParErrBg"), candidate.errorImpactParameter0(), ptCandidate);
        registry.fill(HIST("MC/reconstructed/background/hDecLenErrBg"), candidate.errorDecayLength(), ptCandidate);
        registry.fill(HIST("MC/reconstructed/background/hChi2PCARecBg"), candidate.chi2PCA(), ptCandidate);
      } // Xic background
    } // candidate loop

    // MC gen.
    for (const auto& particle : mcParticles) {
      if (std::abs(particle.flagMcMatchGen()) == hf_decay::hf_cand_3prong::DecayChannelMain::XicToPKPi) {
        auto yGen = RecoDecay::y(particle.pVector(), o2::constants::physics::MassXiCPlus);
        if (yCandGenMax >= 0. && std::abs(yGen) > yCandGenMax) {
          continue;
        }
        auto ptGen = particle.pt();
        registry.fill(HIST("MC/generated/signal/hPtGen"), ptGen);
        registry.fill(HIST("MC/generated/signal/hEtaGen"), particle.eta());
        registry.fill(HIST("MC/generated/signal/hYGen"), yGen);
        registry.fill(HIST("MC/generated/signal/hEtaVsPtGenSig"), particle.eta(), ptGen);
        registry.fill(HIST("MC/generated/signal/hYVsPtGenSig"), yGen, ptGen);

        if (particle.originMcGen() == RecoDecay::OriginType::Prompt) {
          registry.fill(HIST("MC/generated/prompt/hPtGenPrompt"), ptGen);
          registry.fill(HIST("MC/generated/prompt/hEtaGenPrompt"), particle.eta());
          registry.fill(HIST("MC/generated/prompt/hYGenPrompt"), yGen);
          registry.fill(HIST("MC/generated/prompt/hEtaVsPtGenSigPrompt"), particle.eta(), ptGen);
          registry.fill(HIST("MC/generated/prompt/hYVsPtGenSigPrompt"), yGen, ptGen);
        }
        if (particle.originMcGen() == RecoDecay::OriginType::NonPrompt) {
          registry.fill(HIST("MC/generated/nonprompt/hPtGenNonPrompt"), ptGen);
          registry.fill(HIST("MC/generated/nonprompt/hEtaGenNonPrompt"), particle.eta());
          registry.fill(HIST("MC/generated/nonprompt/hYGenNonPrompt"), yGen);
          registry.fill(HIST("MC/generated/nonprompt/hEtaVsPtGenSigNonPrompt"), particle.eta(), ptGen);
          registry.fill(HIST("MC/generated/nonprompt/hYVsPtGenSigNonPrompt"), yGen, ptGen);
        }
      }
    }
  }
  void processMcStd(soa::Filtered<soa::Join<aod::HfCand3ProngWPidPiKaPr, aod::HfSelXicToPKPi, aod::HfCand3ProngMcRec>> const& selectedCandidatesMc,
                    soa::Join<aod::McParticles, aod::HfCand3ProngMcGen> const& mcParticles,
                    aod::TracksWMc const& tracksWithMc)
  {
    analysisMc<false>(selectedCandidatesMc, mcParticles, tracksWithMc);
  }
  PROCESS_SWITCH(HfTaskXic, processMcStd, "Process MC with the standard method", false);

  void processMcWithMl(soa::Filtered<soa::Join<aod::HfCand3ProngWPidPiKaPr, aod::HfSelXicToPKPi, aod::HfMlXicToPKPi, aod::HfCand3ProngMcRec>> const& selectedCandidatesMlMc,
                       soa::Join<aod::McParticles, aod::HfCand3ProngMcGen> const& mcParticles,
                       aod::TracksWMc const& tracksWithMc)
  {
    analysisMc<true>(selectedCandidatesMlMc, mcParticles, tracksWithMc);
  }
  PROCESS_SWITCH(HfTaskXic, processMcWithMl, "Process Mc with the ML method", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskXic>(cfgc)};
}
