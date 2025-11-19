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

/// \file hfTask3Prong.cxx
/// \brief 3-prong candidates analysis task for ALICE 3 simulation studies
/// \author Marcello Di Costanzo <marcello.di.costanzo@cern.ch>, Polytechnic University of Turin and INFN Turin

#include "ALICE3/DataModel/A3DecayFinderTables.h"
#include "ALICE3/DataModel/OTFPIDTrk.h"
#include "ALICE3/DataModel/OTFRICH.h"
#include "ALICE3/DataModel/OTFTOF.h"
#include "ALICE3/DataModel/RICH.h"
#include "ALICE3/Utils/utilsHfAlice3.h"
#include "ALICE3/Utils/utilsSelectionsAlice3.h"
#include "Common/Core/RecoDecay.h"

#include <CCDB/BasicCCDBManager.h>
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
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <THnSparse.h>
#include <TPDGCode.h>

#include <array>
#include <cmath>
#include <numeric>
#include <string>
#include <string_view>
#include <vector>

using namespace o2;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::framework::expressions;

/// Λc± → p± K∓ π± analysis task
struct Alice3HfTask3Prong {
  Configurable<double> yCandGenMax{"yCandGenMax", 0.8, "max. gen particle rapidity"};
  Configurable<double> yCandRecoMax{"yCandRecoMax", 0.8, "max. cand. rapidity"};
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_3prongs_alice3::vecBinsPt}, "pT bin limits"};
  Configurable<bool> fillThn{"fillThn", false, "fill Thn"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> ccdbPathGrp{"ccdbPathGrp", "GLO/GRP/GRP", "Path of the grp file (Run 2)"};
  Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object (Run 3)"};

  HfHelperAlice3 hfHelper;
  SliceCache cache;
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  int selectedPdg{-1};

  using Cands3PReco = soa::Filtered<soa::Join<aod::Alice3Cand3Ps, aod::Alice3Sel3Ps, aod::Alice3McRecFlags>>;
  using Cands3PRecoWMl = soa::Filtered<soa::Join<aod::Alice3Cand3Ps, aod::Alice3Sel3Ps, aod::Alice3Ml3Ps, aod::Alice3McRecFlags>>;
  using Cands3PGen = soa::Join<aod::McParticles, aod::Alice3McGenFlags>;

  Filter filterSelectCandidates = (aod::a3_hf_sel_3prong::isSelMassHypo0 == true || aod::a3_hf_sel_3prong::isSelMassHypo1 == true);

  Partition<Cands3PGen> candsGenLcs = nabs(aod::a3_mc_truth::flagMcGen) == static_cast<int>(CharmHadAlice3::Lc);

  ConfigurableAxis thnConfigAxisPt{"thnConfigAxisPt", {72, 0, 36}, ""};
  ConfigurableAxis thnConfigAxisMass{"thnConfigAxisMass", {300, 1.98, 2.58}, ""};
  ConfigurableAxis thnConfigAxisBdtScoreBkg{"thnConfigAxisBdtScoreBkg", {1000, 0., 1.}, ""};
  ConfigurableAxis thnConfigAxisBdtScoreSignal{"thnConfigAxisBdtScoreSignal", {100, 0., 1.}, ""};
  ConfigurableAxis thnConfigAxisCanType{"thnConfigAxisCanType", {5, 0., 5.}, ""};
  ConfigurableAxis thnAxisRapidity{"thnAxisRapidity", {20, -1, 1}, "Cand. rapidity bins"};
  ConfigurableAxis thnConfigAxisGenPtB{"thnConfigAxisGenPtB", {1000, 0, 100}, "Gen Pt B"};

  HistogramRegistry registry{"registry", {}};

  // Names of folders and suffixes for MC signal histograms
  constexpr static std::string_view SignalFolders[] = {"signal", "prompt", "nonprompt"};
  constexpr static std::string_view SignalSuffixes[] = {"", "Prompt", "NonPrompt"};

  enum SignalClasses : int {
    Signal = 0,
    Prompt,
    NonPrompt
  };

  void init(InitContext&)
  {
    const std::array<bool, 2> doprocess{doprocessLc, doprocessLcWMl};
    if ((std::accumulate(doprocess.begin(), doprocess.end(), 0)) != 1) {
      LOGP(fatal, "no or more than one process function enabled! Please check your configuration!");
    }

    if (doprocessLc || doprocessLcWMl) {
      selectedPdg = CharmHadAlice3::Lc;
    }

    auto addHistogramsRec = [&](const std::string& histoName, const std::string& xAxisTitle, const std::string& yAxisTitle, const HistogramConfigSpec& configSpec) {
      registry.add(("MC/rec/signal/" + histoName + "RecSig").c_str(), ("3-prong cands (matched);" + xAxisTitle + ";" + yAxisTitle).c_str(), configSpec);
      registry.add(("MC/rec/prompt/" + histoName + "RecSigPrompt").c_str(), ("3-prong cands (matched, prompt);" + xAxisTitle + ";" + yAxisTitle).c_str(), configSpec);
      registry.add(("MC/rec/nonprompt/" + histoName + "RecSigNonPrompt").c_str(), ("3-prong cands (matched, non-prompt);" + xAxisTitle + ";" + yAxisTitle).c_str(), configSpec);
    };

    auto addHistogramsGen = [&](const std::string& histoName, const std::string& xAxisTitle, const std::string& yAxisTitle, const HistogramConfigSpec& configSpec) {
      registry.add(("MC/gen/signal/" + histoName + "Gen").c_str(), ("MC particles (matched);" + xAxisTitle + ";" + yAxisTitle).c_str(), configSpec);
      registry.add(("MC/gen/prompt/" + histoName + "GenPrompt").c_str(), ("MC particles (matched, prompt);" + xAxisTitle + ";" + yAxisTitle).c_str(), configSpec);
      registry.add(("MC/gen/nonprompt/" + histoName + "GenNonPrompt").c_str(), ("MC particles (matched, non-prompt);" + xAxisTitle + ";" + yAxisTitle).c_str(), configSpec);
    };

    auto vbins = (std::vector<double>)binsPt;

    /// Reconstructed Histograms
    addHistogramsRec("hMass", "inv. mass (p K #pi) (GeV/#it{c}^{2})", "", {HistType::kTH1F, {{600, 1.98, 2.58}}});
    addHistogramsRec("hPt", "#it{p}_{T}^{rec.} (GeV/#it{c})", "entries", {HistType::kTH1F, {{360, 0., 36.}}});
    addHistogramsRec("hPhi", "#it{#Phi}", "entries", {HistType::kTH1F, {{100, 0., 6.3}}});
    addHistogramsRec("hPtProng0", "prong 0 #it{p}_{T} (GeV/#it{c})", "entries", {HistType::kTH1F, {{360, 0., 36.}}});
    addHistogramsRec("hPtProng1", "prong 1 #it{p}_{T} (GeV/#it{c})", "entries", {HistType::kTH1F, {{360, 0., 36.}}});
    addHistogramsRec("hPtProng2", "prong 2 #it{p}_{T} (GeV/#it{c})", "entries", {HistType::kTH1F, {{360, 0., 36.}}});
    addHistogramsRec("hd0Prong0", "prong 0 DCAxy to prim. vertex (cm)", "entries", {HistType::kTH1F, {{600, -0.4, 0.4}}});
    addHistogramsRec("hd0Prong1", "prong 1 DCAxy to prim. vertex (cm)", "entries", {HistType::kTH1F, {{600, -0.4, 0.4}}});
    addHistogramsRec("hd0Prong2", "prong 2 DCAxy to prim. vertex (cm)", "entries", {HistType::kTH1F, {{600, -0.4, 0.4}}});
    addHistogramsRec("hDecLength", "decay length (cm)", "entries", {HistType::kTH1F, {{400, 0., 1.}}});
    addHistogramsRec("hDecLengthxy", "decay length xy (cm)", "entries", {HistType::kTH1F, {{400, 0., 1.}}});
    addHistogramsRec("hCPA", "cosine of pointing angle", "entries", {HistType::kTH1F, {{110, -1.1, 1.1}}});
    addHistogramsRec("hCPAxy", "cosine of pointing angle xy", "entries", {HistType::kTH1F, {{110, -1.1, 1.1}}});
    addHistogramsRec("hDca2", "prong Chi2PCA to sec. vertex (cm)", "entries", {HistType::kTH1F, {{400, 0., 20.}}});
    addHistogramsRec("hEta", "#it{#eta}", "entries", {HistType::kTH1F, {{100, -2., 2.}}});
    addHistogramsRec("hMassVsPt", "inv. mass (p K #pi) (GeV/#it{c}^{2})", "#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{600, 1.98, 2.58}, {vbins}}});
    addHistogramsRec("hd0VsPtProng0", "prong 0 DCAxy to prim. vertex (cm)", "#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{600, -0.4, 0.4}, {vbins}}});
    addHistogramsRec("hd0VsPtProng1", "prong 1 DCAxy to prim. vertex (cm)", "#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{600, -0.4, 0.4}, {vbins}}});
    addHistogramsRec("hd0VsPtProng2", "prong 2 DCAxy to prim. vertex (cm)", "#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{600, -0.4, 0.4}, {vbins}}});
    addHistogramsRec("hDecLengthVsPt", "decay length (cm)", "#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{400, 0., 1.}, {vbins}}});
    addHistogramsRec("hDecLengthxyVsPt", "decay length xy (cm)", "#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{400, 0., 1.}, {vbins}}});
    addHistogramsRec("hCPAVsPt", "cosine of pointing angle", "#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins}}});
    addHistogramsRec("hCPAxyVsPt", "cosine of pointing angle xy", "#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins}}});
    addHistogramsRec("hDca2VsPt", "prong Chi2PCA to sec. vertex (cm)", "#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{400, 0., 20.}, {vbins}}});
    addHistogramsRec("hEtaVsPt", "candidate #it{#eta}", "#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, -2., 2.}, {vbins}}});
    addHistogramsRec("hPhiVsPt", "candidate #it{#Phi}", "#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, 0., 6.3}, {vbins}}});
    addHistogramsRec("hImpParErrProng0VsPt", "prong 0 impact parameter error (cm)", "#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, -1., 1.}, {vbins}}});
    addHistogramsRec("hImpParErrProng1VsPt", "prong 1 impact parameter error (cm)", "#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, -1., 1.}, {vbins}}});
    addHistogramsRec("hImpParErrProng2VsPt", "prong 2 impact parameter error (cm)", "#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, -1., 1.}, {vbins}}});
    addHistogramsRec("hDecLenErrVsPt", "decay length error (cm)", "#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, 0., 1.}, {vbins}}});

    /// Generated Histograms
    addHistogramsGen("hPt", "#it{p}_{T}^{gen.} (GeV/#it{c})", "entries", {HistType::kTH1F, {{360, 0., 36.}}});
    addHistogramsGen("hEta", "#it{#eta}", "entries", {HistType::kTH1F, {{100, -2., 2.}}});
    addHistogramsGen("hPhi", "#it{#Phi}", "entries", {HistType::kTH1F, {{100, 0., 6.3}}});
    addHistogramsGen("hY", "#it{y}", "entries", {HistType::kTH1F, {{100, -2., 2.}}});
    addHistogramsGen("hEtaVsPt", "#it{#eta}", "#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, -2., 2.}, {vbins}}});
    addHistogramsGen("hYVsPt", "#it{y}", "#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, -2., 2.}, {vbins}}});
    addHistogramsGen("hPhiVsPt", "#it{#Phi}", "#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, 0., 6.3}, {vbins}}});

    /// selection status
    registry.add("hSelectionStatus", "3-prong cands;selection status;entries", {HistType::kTH2F, {{5, -0.5, 4.5}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

    if (fillThn) {
      const AxisSpec thnAxisMass{thnConfigAxisMass, "inv. mass (p K #pi) (GeV/#it{c}^{2})"};
      const AxisSpec thnAxisPt{thnConfigAxisPt, "#it{p}_{T}(#Lambda_{c}^{+}) (GeV/#it{c})"};
      const AxisSpec thnAxisScoreBkg{thnConfigAxisBdtScoreBkg, "BDT bkg score"};
      const AxisSpec thnAxisScorePrompt{thnConfigAxisBdtScoreSignal, "BDT prompt score"};
      const AxisSpec thnAxisScoreNonPrompt{thnConfigAxisBdtScoreSignal, "BDT non-prompt score"};
      const AxisSpec thnAxisCanType{thnConfigAxisCanType, "candidates type"};
      const AxisSpec thnAxisY{thnAxisRapidity, "rapidity"};
      const AxisSpec thnAxisPtB{thnConfigAxisGenPtB, "#it{p}_{T}^{B} (GeV/#it{c})"};

      std::vector<AxisSpec> axesWithBdt = {thnAxisMass, thnAxisPt, thnAxisScoreBkg, thnAxisScorePrompt, thnAxisScoreNonPrompt, thnAxisPtB, thnAxisCanType};
      registry.add("hSparseRec", "Thn for reco cands", HistType::kTHnSparseF, axesWithBdt);
      std::vector<AxisSpec> axesGen = {thnAxisPt, thnAxisY, thnAxisPtB, thnAxisCanType};
      registry.add("hSparseGen", "Thn for gen cands", HistType::kTHnSparseF, axesGen);
    }

    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
  }

  /// Helper function for filling MC reconstructed histograms for prompt, nonpromt and common (signal)
  /// \param candidate is a reconstructed candidate
  /// \tparam SignalType is an enum defining which histogram in which folder (signal, prompt or nonpromt) to fill
  template <CharmHadAlice3 CharmHad, int SignalType, typename CandidateType>
  void fillHistogramsRecSig(CandidateType const& candidate, float mass)
  {
    registry.fill(HIST("MC/rec/") + HIST(SignalFolders[SignalType]) + HIST("/hMassRecSig") + HIST(SignalSuffixes[SignalType]), mass);
    registry.fill(HIST("MC/rec/") + HIST(SignalFolders[SignalType]) + HIST("/hMassVsPtRecSig") + HIST(SignalSuffixes[SignalType]), mass, candidate.pt());
    registry.fill(HIST("MC/rec/") + HIST(SignalFolders[SignalType]) + HIST("/hPtRecSig") + HIST(SignalSuffixes[SignalType]), candidate.pt());
    registry.fill(HIST("MC/rec/") + HIST(SignalFolders[SignalType]) + HIST("/hPtProng0RecSig") + HIST(SignalSuffixes[SignalType]), candidate.ptProng0());
    registry.fill(HIST("MC/rec/") + HIST(SignalFolders[SignalType]) + HIST("/hPtProng1RecSig") + HIST(SignalSuffixes[SignalType]), candidate.ptProng1());
    registry.fill(HIST("MC/rec/") + HIST(SignalFolders[SignalType]) + HIST("/hPtProng2RecSig") + HIST(SignalSuffixes[SignalType]), candidate.ptProng2());

    registry.fill(HIST("MC/rec/") + HIST(SignalFolders[SignalType]) + HIST("/hd0Prong0RecSig") + HIST(SignalSuffixes[SignalType]), candidate.impactParameterY0());
    registry.fill(HIST("MC/rec/") + HIST(SignalFolders[SignalType]) + HIST("/hd0Prong1RecSig") + HIST(SignalSuffixes[SignalType]), candidate.impactParameterY1());
    registry.fill(HIST("MC/rec/") + HIST(SignalFolders[SignalType]) + HIST("/hd0Prong2RecSig") + HIST(SignalSuffixes[SignalType]), candidate.impactParameterY2());
    registry.fill(HIST("MC/rec/") + HIST(SignalFolders[SignalType]) + HIST("/hd0VsPtProng0RecSig") + HIST(SignalSuffixes[SignalType]), candidate.impactParameterY0(), candidate.pt());
    registry.fill(HIST("MC/rec/") + HIST(SignalFolders[SignalType]) + HIST("/hd0VsPtProng1RecSig") + HIST(SignalSuffixes[SignalType]), candidate.impactParameterY1(), candidate.pt());
    registry.fill(HIST("MC/rec/") + HIST(SignalFolders[SignalType]) + HIST("/hd0VsPtProng2RecSig") + HIST(SignalSuffixes[SignalType]), candidate.impactParameterY2(), candidate.pt());
    registry.fill(HIST("MC/rec/") + HIST(SignalFolders[SignalType]) + HIST("/hDecLengthRecSig") + HIST(SignalSuffixes[SignalType]), candidate.decayLength());
    registry.fill(HIST("MC/rec/") + HIST(SignalFolders[SignalType]) + HIST("/hDecLengthVsPtRecSig") + HIST(SignalSuffixes[SignalType]), candidate.decayLength(), candidate.pt());
    registry.fill(HIST("MC/rec/") + HIST(SignalFolders[SignalType]) + HIST("/hDecLengthxyRecSig") + HIST(SignalSuffixes[SignalType]), candidate.decayLengthXY());
    registry.fill(HIST("MC/rec/") + HIST(SignalFolders[SignalType]) + HIST("/hDecLengthxyVsPtRecSig") + HIST(SignalSuffixes[SignalType]), candidate.decayLengthXY(), candidate.pt());
    registry.fill(HIST("MC/rec/") + HIST(SignalFolders[SignalType]) + HIST("/hCPARecSig") + HIST(SignalSuffixes[SignalType]), candidate.cpa());
    registry.fill(HIST("MC/rec/") + HIST(SignalFolders[SignalType]) + HIST("/hCPAVsPtRecSig") + HIST(SignalSuffixes[SignalType]), candidate.cpa(), candidate.pt());
    registry.fill(HIST("MC/rec/") + HIST(SignalFolders[SignalType]) + HIST("/hCPAxyRecSig") + HIST(SignalSuffixes[SignalType]), candidate.cpaXY());
    registry.fill(HIST("MC/rec/") + HIST(SignalFolders[SignalType]) + HIST("/hCPAxyVsPtRecSig") + HIST(SignalSuffixes[SignalType]), candidate.cpaXY(), candidate.pt());
    registry.fill(HIST("MC/rec/") + HIST(SignalFolders[SignalType]) + HIST("/hDca2RecSig") + HIST(SignalSuffixes[SignalType]), candidate.chi2PCA());
    registry.fill(HIST("MC/rec/") + HIST(SignalFolders[SignalType]) + HIST("/hDca2VsPtRecSig") + HIST(SignalSuffixes[SignalType]), candidate.chi2PCA(), candidate.pt());
    registry.fill(HIST("MC/rec/") + HIST(SignalFolders[SignalType]) + HIST("/hEtaRecSig") + HIST(SignalSuffixes[SignalType]), candidate.eta());
    registry.fill(HIST("MC/rec/") + HIST(SignalFolders[SignalType]) + HIST("/hEtaVsPtRecSig") + HIST(SignalSuffixes[SignalType]), candidate.eta(), candidate.pt());
    registry.fill(HIST("MC/rec/") + HIST(SignalFolders[SignalType]) + HIST("/hPhiRecSig") + HIST(SignalSuffixes[SignalType]), candidate.phi());
    registry.fill(HIST("MC/rec/") + HIST(SignalFolders[SignalType]) + HIST("/hPhiVsPtRecSig") + HIST(SignalSuffixes[SignalType]), candidate.phi(), candidate.pt());
    registry.fill(HIST("MC/rec/") + HIST(SignalFolders[SignalType]) + HIST("/hImpParErrProng0VsPtRecSig") + HIST(SignalSuffixes[SignalType]), candidate.errorImpactParameterY0(), candidate.pt());
    registry.fill(HIST("MC/rec/") + HIST(SignalFolders[SignalType]) + HIST("/hImpParErrProng1VsPtRecSig") + HIST(SignalSuffixes[SignalType]), candidate.errorImpactParameterY1(), candidate.pt());
    registry.fill(HIST("MC/rec/") + HIST(SignalFolders[SignalType]) + HIST("/hImpParErrProng2VsPtRecSig") + HIST(SignalSuffixes[SignalType]), candidate.errorImpactParameterY2(), candidate.pt());
    registry.fill(HIST("MC/rec/") + HIST(SignalFolders[SignalType]) + HIST("/hDecLenErrVsPtRecSig") + HIST(SignalSuffixes[SignalType]), candidate.errorDecayLength(), candidate.pt());
  }

  /// Fill MC histograms at reconstruction level
  /// \tparam CharmHad is the charm hadron species
  /// \tparam SaveMl indicates whether ML scores are saved in the THnSparse
  /// \tparam CandsRec is the type of the reconstructed candidates collection
  /// \param candidates is the collection of reconstructed candidates
  template <CharmHadAlice3 CharmHad, bool SaveMl, typename CandsRec>
  void fillHistosMcRec(CandsRec const& candidates)
  {
    for (const auto& candidate : candidates) {
      /// rapidity selection
      if (yCandRecoMax >= 0. && std::abs(hfHelper.getCandY<CharmHad>(candidate)) > yCandRecoMax) {
        continue;
      }

      if (candidate.flagMcRec() != 0) {
        // Get the corresponding MC particle.

        const auto pt = candidate.pt();
        const auto originType = candidate.originMcRec();

        if (fillThn) {
          if (candidate.isSelMassHypo0()) {
            registry.fill(HIST("hSelectionStatus"), 0., pt);
            double mass = hfHelper.getCandMass<CharmHad, false>(candidate);
            /// Fill histograms
            fillHistogramsRecSig<CharmHad, Signal>(candidate, mass);
            if (originType == RecoDecay::OriginType::Prompt) {
              fillHistogramsRecSig<CharmHad, Prompt>(candidate, mass);
            } else if (originType == RecoDecay::OriginType::NonPrompt) {
              fillHistogramsRecSig<CharmHad, NonPrompt>(candidate, mass);
            }
            std::vector<double> valuesToFill{mass, pt};
            if constexpr (SaveMl) {
              LOGP(fatal, "Trying to access ML scores, but SaveMl is false!");
              valuesToFill.push_back(candidate.mlScore0());
              valuesToFill.push_back(candidate.mlScore1());
              valuesToFill.push_back(candidate.mlScore2());
            }
            valuesToFill.push_back(static_cast<double>(originType));
            registry.get<THnSparse>(HIST("hSparseRec"))->Fill(valuesToFill.data());
          }
          if (candidate.isSelMassHypo1()) {
            registry.fill(HIST("hSelectionStatus"), 1., pt);
            double mass = hfHelper.getCandMass<CharmHad, true>(candidate);
            /// Fill histograms
            fillHistogramsRecSig<CharmHad, Signal>(candidate, mass);
            if (originType == RecoDecay::OriginType::Prompt) {
              fillHistogramsRecSig<CharmHad, Prompt>(candidate, mass);
            } else if (originType == RecoDecay::OriginType::NonPrompt) {
              fillHistogramsRecSig<CharmHad, NonPrompt>(candidate, mass);
            }
            std::vector<double> valuesToFill{mass, pt};
            if constexpr (SaveMl) {
              LOGP(fatal, "Trying to access ML scores, but SaveMl is false!");
              valuesToFill.push_back(candidate.mlScore0());
              valuesToFill.push_back(candidate.mlScore1());
              valuesToFill.push_back(candidate.mlScore2());
            }
            valuesToFill.push_back(static_cast<double>(originType));
            registry.get<THnSparse>(HIST("hSparseRec"))->Fill(valuesToFill.data());
          }
        }
      }
    }
  }

  /// Helper function for filling MC generated histograms for prompt, nonpromt and common (signal)
  /// \tparam CharmHad is the charm hadron species
  /// \tparam SignalType is an enum defining which histogram in which folder (signal, prompt or nonpromt) to fill
  /// \tparam ParticleType is the type of the generated particle
  /// \param particle is a generated particle
  template <CharmHadAlice3 CharmHad, int SignalType, typename ParticleType>
  void fillHistogramsGen(ParticleType const& particle)
  {
    LOG(debug) << "Filling generated histograms for signal type " << SignalType;
    registry.fill(HIST("MC/gen/") + HIST(SignalFolders[SignalType]) + HIST("/hPtGen") + HIST(SignalSuffixes[SignalType]), particle.pt());
    registry.fill(HIST("MC/gen/") + HIST(SignalFolders[SignalType]) + HIST("/hEtaGen") + HIST(SignalSuffixes[SignalType]), particle.eta());
    registry.fill(HIST("MC/gen/") + HIST(SignalFolders[SignalType]) + HIST("/hYGen") + HIST(SignalSuffixes[SignalType]), hfHelper.getCandY<CharmHad>(particle));
    registry.fill(HIST("MC/gen/") + HIST(SignalFolders[SignalType]) + HIST("/hPhiGen") + HIST(SignalSuffixes[SignalType]), particle.phi());
    registry.fill(HIST("MC/gen/") + HIST(SignalFolders[SignalType]) + HIST("/hEtaVsPtGen") + HIST(SignalSuffixes[SignalType]), particle.eta(), particle.pt());
    registry.fill(HIST("MC/gen/") + HIST(SignalFolders[SignalType]) + HIST("/hYVsPtGen") + HIST(SignalSuffixes[SignalType]), hfHelper.getCandY<CharmHad>(particle), particle.pt());
    registry.fill(HIST("MC/gen/") + HIST(SignalFolders[SignalType]) + HIST("/hPhiVsPtGen") + HIST(SignalSuffixes[SignalType]), particle.phi(), particle.pt());
  }

  /// Fill MC histograms at generated level
  /// \tparam CharmHad is the charm hadron species
  /// \tparam CandsGen is the type of the generated candidates collection
  /// \param mcParticles is the collection of generated particles
  template <CharmHadAlice3 CharmHad, typename CandsGen>
  void fillHistosMcGen(CandsGen const& mcParticles)
  {
    // MC gen.
    for (const auto& particle : mcParticles) {
      if (std::abs(particle.flagMcGen()) == selectedPdg) {
        double yGen = hfHelper.getCandY<CharmHad>(particle);
        if (yCandGenMax >= 0. && std::abs(yGen) > yCandGenMax) {
          continue;
        }
        const auto ptGen = particle.pt();
        const auto originType = particle.originMcGen();

        fillHistogramsGen<CharmHad, Signal>(particle);

        float ptGenB = -1.f;
        if (originType == RecoDecay::OriginType::Prompt) {
          fillHistogramsGen<CharmHad, Prompt>(particle);
        } else if (particle.originMcGen() == RecoDecay::OriginType::NonPrompt) {
          ptGenB = particle.bHadMotherPtGen();
          fillHistogramsGen<CharmHad, NonPrompt>(particle);
        }

        if (fillThn) {
          std::vector<double> valuesToFill{ptGen, yGen, ptGenB, static_cast<double>(originType)};
          registry.get<THnSparse>(HIST("hSparseGen"))->Fill(valuesToFill.data());
        }
      }
    }
  }

  void processLc(Cands3PReco const& candsLc,
                 Cands3PGen const&)
  {
    fillHistosMcRec<CharmHadAlice3::Lc, false>(candsLc);
    fillHistosMcGen<CharmHadAlice3::Lc>(candsGenLcs);
  }
  PROCESS_SWITCH(Alice3HfTask3Prong, processLc, "Process Lc w/o ML sels", true);

  void processLcWMl(Cands3PRecoWMl const& candsLcWMl,
                    Cands3PGen const&)
  {
    fillHistosMcRec<CharmHadAlice3::Lc, true>(candsLcWMl);
    fillHistosMcGen<CharmHadAlice3::Lc>(candsGenLcs);
  }
  PROCESS_SWITCH(Alice3HfTask3Prong, processLcWMl, "Process Lc with ML sels", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<Alice3HfTask3Prong>(cfgc)};
}
