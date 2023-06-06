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

  float etaMaxAcceptance = 0.8;
  float ptMinAcceptance = 0.1;

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

      {"MC/reconstructed/signal/hPtRecSig", "3-prong candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}}}, // pt Xic
      {"MC/reconstructed/background/hPtRecBg", "3-prong candidates (unmatched);#it{p}_{T}^{rec.} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}}},
      {"MC/generated/hPtGen", "MC particles (matched);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}}},
      {"MC/generated/hPtGenWithProngsInAcceptance", "MC particles (generated-daughters in acceptance); #it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}}},
      {"MC/generated/hEtaGen", "MC particles; #it{eta}^{gen} ;#it{p}_{T}^{gen.};entries", {HistType::kTH1F, {{100, -2., 2.}}}},
      {"MC/generated/hYGen", "MC particles;  #it{y}^{gen} ;#it{p}_{T}^{gen.} ;entries", {HistType::kTH1F, {{100, -2., 2.}}}},
      // add generated in acceptance!!

    }};

  void init(o2::framework::InitContext&)
  {

    AxisSpec axisPPid = {100, 0.f, 10.0f, "#it{p} (GeV/#it{c})"};
    AxisSpec axisNSigmaPr = {100, -6.f, 6.f, "n#it{#sigma}_{p}"};
    AxisSpec axisNSigmaPi = {100, -6.f, 6.f, "n#it{#sigma}_{#pi}"};
    AxisSpec axisNSigmaKa = {100, -6.f, 6.f, "n#it{#sigma}_{K}"};

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

    // background
    registry.add("MC/reconstructed/background/hMassBg", "Invariant mass (unmatched);m (p K #pi) (GeV/#it{c}^{2});;entries", {HistType::kTH2F, {{500, 1.6, 3.1}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/background/hDecLengthRecBg", "3-prong candidates;decay length (cm);;entries", {HistType::kTH2F, {{200, 0., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/background/hDecLengthXYRecBg", "3-prong candidates;decay length xy(cm);;entries", {HistType::kTH2F, {{200, 0., 2.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/background/hNormalisedDecayLengthXYRecBg", "3-prong candidates;norm. decay length XY (cm);;entries", {HistType::kTH2F, {{80, 0., 80.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("MC/reconstructed/background/hImpParXYRecBg", "3-prong candidates;candidate DCAxy to prim. vertex (cm);;entries",{HistType::kTH2F, {{100, -1., 1.}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
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
  }

  /// Selection of Xic daughters in geometrical acceptanc
  /// \param etaProng is the pseudorapidity of Xic prong
  /// \param ptProng is the pT of Xic prong
  /// \return true if prong is in geometrical acceptance
  template <typename T = float>
  bool isProngInAcceptance(const T& etaProng, const T& ptProng)
  {
    return std::abs(etaProng) <= etaMaxAcceptance && ptProng >= ptMinAcceptance;
  }

  void process(aod::Collision const& collision, soa::Join<aod::Tracks, aod::TracksDCA> const& tracks, soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelXicToPKPi>> const& candidates, aod::BigTracksPID const&)
  {
    int nTracks = 0;

    if (collision.numContrib() > 1) {
      for (auto const& track : tracks) {
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
      registry.fill(HIST("Data/hCt"), ctXic(candidate), ptCandidate);
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
        registry.fill(HIST("MC/reconstructed/signal/hPtRecSig"), ptCandidate);     // rec. level pT

        if (massXicToPKPi != 0.) {
          registry.fill(HIST("MC/reconstructed/signal/hMassSig"), massXicToPKPi, ptCandidate);
        }
        if (massXicToPiKP != 0.) {
          registry.fill(HIST("MC/reconstructed/signal/hMassSig"), massXicToPiKP, ptCandidate);
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
        registry.fill(HIST("MC/reconstructed/signal/hCtRecSig"), ctXic(candidate), ptCandidate);
        registry.fill(HIST("MC/reconstructed/signal/hCPARecSig"), candidate.cpa(), ptCandidate);
        registry.fill(HIST("MC/reconstructed/signal/hCPAXYRecSig"), candidate.cpaXY(), ptCandidate);
        registry.fill(HIST("MC/reconstructed/signal/hEtaRecSig"), candidate.eta(), ptCandidate);
        registry.fill(HIST("MC/reconstructed/signal/hImpParErrSig"), candidate.errorImpactParameter0(), ptCandidate);
        registry.fill(HIST("MC/reconstructed/signal/hDecLenErrSig"), candidate.errorDecayLength(), ptCandidate);
        registry.fill(HIST("MC/reconstructed/signal/hChi2PCARecSig"), candidate.chi2PCA(), ptCandidate);
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
        registry.fill(HIST("MC/reconstructed/background/hCtRecBg"), ctXic(candidate), ptCandidate);
        registry.fill(HIST("MC/reconstructed/background/hCPARecBg"), candidate.cpa(), ptCandidate);
        registry.fill(HIST("MC/reconstructed/background/hCPAXYRecBg"), candidate.cpaXY(), ptCandidate);
        registry.fill(HIST("MC/reconstructed/background/hEtaRecBg"), candidate.eta(), ptCandidate);
        registry.fill(HIST("MC/reconstructed/background/hImpParErrBg"), candidate.errorImpactParameter0(), ptCandidate);
        registry.fill(HIST("MC/reconstructed/background/hDecLenErrBg"), candidate.errorDecayLength(), ptCandidate);
        registry.fill(HIST("MC/reconstructed/background/hChi2PCARecBg"), candidate.chi2PCA(), ptCandidate);
      }
    }
    // MC gen.
    for (auto const& particle : particlesMC) {
      if (std::abs(particle.flagMcMatchGen()) == 1 << DecayType::XicToPKPi) {
        auto yParticle = RecoDecay::y(array{particle.px(), particle.py(), particle.pz()}, RecoDecay::getMassPDG(particle.pdgCode()));
        if (yCandMax >= 0. && std::abs(yParticle) > yCandMax) {
          continue;
        }
        auto ptParticle = particle.pt();
        std::array<float, 3> ptProngs;
        std::array<float, 3> yProngs;
        std::array<float, 3> etaProngs;
        int counter = 0;
        for (auto const& daught : particle.daughters_as<aod::McParticles>()) {
          ptProngs[counter] = daught.pt();
          etaProngs[counter] = daught.eta();
          yProngs[counter] = RecoDecay::y(array{daught.px(), daught.py(), daught.pz()}, RecoDecay::getMassPDG(daught.pdgCode()));
          counter++;
        }

        registry.fill(HIST("MC/generated/hPtGen"), ptParticle);
        registry.fill(HIST("MC/generated/hEtaGen"), particle.eta(), ptParticle);
        registry.fill(HIST("MC/generated/hYGen"), yParticle, ptParticle);
        // reject Xic daughters that are not in geometrical acceptance
        if (!isProngInAcceptance(etaProngs[0], ptProngs[0]) || !isProngInAcceptance(etaProngs[1], ptProngs[1]) || !isProngInAcceptance(etaProngs[2], ptProngs[2])) {
          continue;
        }
        registry.fill(HIST("MC/generated/hPtGenWithProngsInAcceptance"), ptParticle);
        registry.fill(HIST("MC/generated/hYGenWithProngsInAcceptance"), yParticle, ptParticle);
        registry.fill(HIST("MC/generated/hEtaGenWithProngsInAcceptance"), particle.eta(), ptParticle);
      }
    }
  }

  PROCESS_SWITCH(HfTaskXic, processMc, "Process MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskXic>(cfgc)};
}
