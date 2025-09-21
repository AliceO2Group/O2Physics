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

#include "PWGHF/Core/CentralityEstimation.h"
#include "PWGHF/Core/DecayChannels.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/Utils/utilsEvSelHf.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/runDataProcessing.h>

#include <THnSparse.h>
#include <TPDGCode.h>

#include <array>
#include <numeric>
#include <string>
#include <vector> // std::vector

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
  Configurable<bool> storeProperLifetime{"storeProperLifetime", false, "Flag to store proper lifetime"};

  constexpr static float CtToProperLifetimePs = 1.f / o2::constants::physics::LightSpeedCm2PS;
  constexpr static float NanoToPico = 1000.f;

  enum MlClasses : int {
    MlClassBackground = 0,
    MlClassPrompt,
    MlClassNonPrompt,
    NumberOfMlClasses
  };

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
  ConfigurableAxis thnConfigAxisProperLifetime{"thnConfigAxisProperLifetime", {200, 0, 2}, "Proper lifetime, ps"};

  HistogramRegistry registry{"registry", {}};

  void init(InitContext&)
  {
    std::array<bool, 12> doprocess{doprocessDataStd, doprocessDataStdWithFT0C, doprocessDataStdWithFT0M, doprocessDataWithMl, doprocessDataWithMlWithFT0C, doprocessDataWithMlWithFT0M, doprocessMcStd, doprocessMcStdWithFT0C, doprocessMcStdWithFT0M, doprocessMcWithMl, doprocessMcWithMlWithFT0C, doprocessMcWithMlWithFT0M};
    if ((std::accumulate(doprocess.begin(), doprocess.end(), 0)) != 1) {
      LOGP(fatal, "no or more than one process function enabled! Please check your configuration!");
    }

    auto addHistogramsRec = [&](const std::string& histoName, const std::string& xAxisTitle, const std::string& yAxisTitle, const HistogramConfigSpec& configSpec) {
      registry.add(("Data/" + histoName).c_str(), ("3-prong candidates;" + xAxisTitle + ";" + yAxisTitle).c_str(), configSpec);
      registry.add(("MC/reconstructed/signal/" + histoName + "RecSig").c_str(), ("3-prong candidates (matched);" + xAxisTitle + ";" + yAxisTitle).c_str(), configSpec);
      registry.add(("MC/reconstructed/prompt/" + histoName + "RecSigPrompt").c_str(), ("3-prong candidates (matched, prompt);" + xAxisTitle + ";" + yAxisTitle).c_str(), configSpec);
      registry.add(("MC/reconstructed/nonprompt/" + histoName + "RecSigNonPrompt").c_str(), ("3-prong candidates (matched, non-prompt);" + xAxisTitle + ";" + yAxisTitle).c_str(), configSpec);
    };

    auto addHistogramsGen = [&](const std::string& histoName, const std::string& xAxisTitle, const std::string& yAxisTitle, const HistogramConfigSpec& configSpec) {
      registry.add(("MC/generated/signal/" + histoName + "Gen").c_str(), ("MC particles (matched);" + xAxisTitle + ";" + yAxisTitle).c_str(), configSpec);
      registry.add(("MC/generated/prompt/" + histoName + "GenPrompt").c_str(), ("MC particles (matched, prompt);" + xAxisTitle + ";" + yAxisTitle).c_str(), configSpec);
      registry.add(("MC/generated/nonprompt/" + histoName + "GenNonPrompt").c_str(), ("MC particles (matched, non-prompt);" + xAxisTitle + ";" + yAxisTitle).c_str(), configSpec);
    };

    /// mass candidate
    addHistogramsRec("hMass", "inv. mass (p K #pi) (GeV/#it{c}^{2})", "", {HistType::kTH1F, {{600, 1.98, 2.58}}});
    /// pT
    addHistogramsRec("hPt", "#it{p}_{T}^{rec.} (GeV/#it{c})", "entries", {HistType::kTH1F, {{360, 0., 36.}}});
    addHistogramsGen("hPt", "#it{p}_{T}^{gen.} (GeV/#it{c})", "entries", {HistType::kTH1F, {{360, 0., 36.}}});
    registry.add("MC/generated/signal/hPtGenSig", "3-prong candidates (matched);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{360, 0., 36.}}});
    addHistogramsRec("hPtProng0", "prong 0 #it{p}_{T} (GeV/#it{c})", "entries", {HistType::kTH1F, {{360, 0., 36.}}});
    addHistogramsRec("hPtProng1", "prong 1 #it{p}_{T} (GeV/#it{c})", "entries", {HistType::kTH1F, {{360, 0., 36.}}});
    addHistogramsRec("hPtProng2", "prong 2 #it{p}_{T} (GeV/#it{c})", "entries", {HistType::kTH1F, {{360, 0., 36.}}});
    /// DCAxy to prim. vertex prongs
    addHistogramsRec("hd0Prong0", "prong 0 DCAxy to prim. vertex (cm)", "entries", {HistType::kTH1F, {{600, -0.4, 0.4}}});
    addHistogramsRec("hd0Prong1", "prong 1 DCAxy to prim. vertex (cm)", "entries", {HistType::kTH1F, {{600, -0.4, 0.4}}});
    addHistogramsRec("hd0Prong2", "prong 2 DCAxy to prim. vertex (cm)", "entries", {HistType::kTH1F, {{600, -0.4, 0.4}}});
    /// decay length candidate
    addHistogramsRec("hDecLength", "decay length (cm)", "entries", {HistType::kTH1F, {{400, 0., 1.}}});
    /// decay length xy candidate
    addHistogramsRec("hDecLengthxy", "decay length xy (cm)", "entries", {HistType::kTH1F, {{400, 0., 1.}}});
    /// proper lifetime
    addHistogramsRec("hCt", "proper lifetime (#Lambda_{c}) * #it{c} (cm)", "entries", {HistType::kTH1F, {{100, 0., 0.2}}});
    /// cosine of pointing angle
    addHistogramsRec("hCPA", "cosine of pointing angle", "entries", {HistType::kTH1F, {{110, -1.1, 1.1}}});
    /// cosine of pointing angle xy
    addHistogramsRec("hCPAxy", "cosine of pointing angle xy", "entries", {HistType::kTH1F, {{110, -1.1, 1.1}}});
    /// Chi 2 PCA to sec. vertex
    addHistogramsRec("hDca2", "prong Chi2PCA to sec. vertex (cm)", "entries", {HistType::kTH1F, {{400, 0., 20.}}});
    /// eta
    addHistogramsRec("hEta", "#it{#eta}", "entries", {HistType::kTH1F, {{100, -2., 2.}}});
    addHistogramsGen("hEta", "#it{#eta}", "entries", {HistType::kTH1F, {{100, -2., 2.}}});
    addHistogramsGen("hY", "#it{y}", "entries", {HistType::kTH1F, {{100, -2., 2.}}});
    /// phi
    addHistogramsRec("hPhi", "#it{#Phi}", "entries", {HistType::kTH1F, {{100, 0., 6.3}}});
    addHistogramsGen("hPhi", "#it{#Phi}", "entries", {HistType::kTH1F, {{100, 0., 6.3}}});

    auto vbins = (std::vector<double>)binsPt;
    /// mass candidate
    registry.add("Data/hMassVsPtVsNPvContributors", "3-prong candidates;inv. mass (p K #pi) (GeV/#it{c}^{2}); p_{T}; Number of PV contributors", {HistType::kTH3F, {{600, 1.98, 2.58}, {vbins, "#it{p}_{T} (GeV/#it{c})"}, {5000, 0., 10000.}}});
    addHistogramsRec("hMassVsPt", "inv. mass (p K #pi) (GeV/#it{c}^{2})", "#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{600, 1.98, 2.58}, {vbins}}});
    /// DCAxy to prim. vertex prongs
    addHistogramsRec("hd0VsPtProng0", "prong 0 DCAxy to prim. vertex (cm)", "#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{600, -0.4, 0.4}, {vbins}}});
    addHistogramsRec("hd0VsPtProng1", "prong 1 DCAxy to prim. vertex (cm)", "#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{600, -0.4, 0.4}, {vbins}}});
    addHistogramsRec("hd0VsPtProng2", "prong 2 DCAxy to prim. vertex (cm)", "#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{600, -0.4, 0.4}, {vbins}}});

    /// decay length candidate
    addHistogramsRec("hDecLengthVsPt", "decay length (cm)", "#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{400, 0., 1.}, {vbins}}});

    /// decay length xy candidate
    addHistogramsRec("hDecLengthxyVsPt", "decay length xy (cm)", "#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{400, 0., 1.}, {vbins}}});

    /// proper lifetime
    addHistogramsRec("hCtVsPt", "proper lifetime (#Lambda_{c}) * #it{c} (cm)", "#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, 0., 0.2}, {vbins}}});

    /// cosine of pointing angle
    addHistogramsRec("hCPAVsPt", "cosine of pointing angle", "#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins}}});

    /// cosine of pointing angle xy
    addHistogramsRec("hCPAxyVsPt", "cosine of pointing angle xy", "#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{110, -1.1, 1.1}, {vbins}}});

    /// Chi 2 PCA to sec. vertex
    addHistogramsRec("hDca2VsPt", "prong Chi2PCA to sec. vertex (cm)", "#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{400, 0., 20.}, {vbins}}});
    /// eta
    addHistogramsRec("hEtaVsPt", "candidate #it{#eta}", "#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, -2., 2.}, {vbins}}});
    addHistogramsGen("hEtaVsPt", "#it{#eta}", "#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, -2., 2.}, {vbins}}});

    /// y
    addHistogramsGen("hYVsPt", "#it{y}", "#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, -2., 2.}, {vbins}}});

    /// phi
    addHistogramsRec("hPhiVsPt", "candidate #it{#Phi}", "#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, 0., 6.3}, {vbins}}});
    addHistogramsGen("hPhiVsPt", "#it{#Phi}", "#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, 0., 6.3}, {vbins}}});

    /// selection status
    registry.add("hSelectionStatus", "3-prong candidates;selection status;entries", {HistType::kTH2F, {{5, -0.5, 4.5}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    /// impact parameter error
    addHistogramsRec("hImpParErrProng0VsPt", "prong 0 impact parameter error (cm)", "#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, -1., 1.}, {vbins}}});
    addHistogramsRec("hImpParErrProng1VsPt", "prong 1 impact parameter error (cm)", "#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, -1., 1.}, {vbins}}});
    addHistogramsRec("hImpParErrProng2VsPt", "prong 2 impact parameter error (cm)", "#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, -1., 1.}, {vbins}}});
    /// decay length error
    addHistogramsRec("hDecLenErrVsPt", "decay length error (cm)", "#it{p}_{T} (GeV/#it{c})", {HistType::kTH2F, {{100, 0., 1.}, {vbins}}});

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
      const AxisSpec thnAxisProperLifetime{thnConfigAxisProperLifetime, "T_{proper} (ps)"};

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
        for (const auto& axes : std::array<std::vector<AxisSpec>*, 3>{&axesWithBdt, &axesStd, &axesGen}) {
          if (!axes->empty()) {
            axes->push_back(thnAxisOccupancy);
          }
        }
      }
      if (storeProperLifetime) {
        for (const auto& axes : std::array<std::vector<AxisSpec>*, 3>{&axesWithBdt, &axesStd, &axesGen}) {
          if (!axes->empty()) {
            axes->push_back(thnAxisProperLifetime);
          }
        }
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

      if (std::abs(candidate.flagMcMatchRec()) == hf_decay::hf_cand_3prong::DecayChannelMain::LcToPKPi) {
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

        auto fillHistogramsRecSig = [&](const std::string& signalFolder, const std::string& signalSuffix) {
          if ((candidate.isSelLcToPKPi() >= selectionFlagLc) && pdgCodeProng0 == kProton) {
            registry.fill(HIST(("MC/reconstructed/" + signalFolder + "/hMassRecSig" + signalSuffix)).c_str(), hfHelper.invMassLcToPKPi(candidate));
            registry.fill(HIST(("MC/reconstructed/" + signalFolder + "/hMassVsPtRecSig" + signalSuffix)).c_str(), hfHelper.invMassLcToPKPi(candidate), pt);
          }
          if ((candidate.isSelLcToPiKP() >= selectionFlagLc) && pdgCodeProng0 == kPiPlus) {
            registry.fill(HIST(("MC/reconstructed/" + signalFolder + "/hMassRecSig" + signalSuffix)).c_str(), hfHelper.invMassLcToPiKP(candidate));
            registry.fill(HIST(("MC/reconstructed/" + signalFolder + "/hMassVsPtRecSig" + signalSuffix)).c_str(), hfHelper.invMassLcToPiKP(candidate), pt);
          }
          registry.fill(HIST(("MC/reconstructed/" + signalFolder + "/hPtRecSig" + signalSuffix)).c_str(), pt);
          registry.fill(HIST(("MC/reconstructed/" + signalFolder + "/hPtProng0RecSig" + signalSuffix)).c_str(), ptProng0);
          registry.fill(HIST(("MC/reconstructed/" + signalFolder + "/hPtProng1RecSig" + signalSuffix)).c_str(), ptProng1);
          registry.fill(HIST(("MC/reconstructed/" + signalFolder + "/hPtProng2RecSig" + signalSuffix)).c_str(), ptProng2);
          registry.fill(HIST(("MC/reconstructed/" + signalFolder + "/hd0Prong0RecSig" + signalSuffix)).c_str(), candidate.impactParameter0());
          registry.fill(HIST(("MC/reconstructed/" + signalFolder + "/hd0Prong1RecSig" + signalSuffix)).c_str(), candidate.impactParameter1());
          registry.fill(HIST(("MC/reconstructed/" + signalFolder + "/hd0Prong2RecSig" + signalSuffix)).c_str(), candidate.impactParameter2());
          registry.fill(HIST(("MC/reconstructed/" + signalFolder + "/hd0VsPtProng0RecSig" + signalSuffix)).c_str(), candidate.impactParameter0(), pt);
          registry.fill(HIST(("MC/reconstructed/" + signalFolder + "/hd0VsPtProng1RecSig" + signalSuffix)).c_str(), candidate.impactParameter1(), pt);
          registry.fill(HIST(("MC/reconstructed/" + signalFolder + "/hd0VsPtProng2RecSig" + signalSuffix)).c_str(), candidate.impactParameter2(), pt);
          registry.fill(HIST(("MC/reconstructed/" + signalFolder + "/hDecLengthRecSig" + signalSuffix)).c_str(), decayLength);
          registry.fill(HIST(("MC/reconstructed/" + signalFolder + "/hDecLengthVsPtRecSig" + signalSuffix)).c_str(), decayLength, pt);
          registry.fill(HIST(("MC/reconstructed/" + signalFolder + "/hDecLengthxyRecSig" + signalSuffix)).c_str(), decayLengthXY);
          registry.fill(HIST(("MC/reconstructed/" + signalFolder + "/hDecLengthxyVsPtRecSig" + signalSuffix)).c_str(), decayLengthXY, pt);
          registry.fill(HIST(("MC/reconstructed/" + signalFolder + "/hCtRecSig" + signalSuffix)).c_str(), hfHelper.ctLc(candidate));
          registry.fill(HIST(("MC/reconstructed/" + signalFolder + "/hCtVsPtRecSig" + signalSuffix)).c_str(), hfHelper.ctLc(candidate), pt);
          registry.fill(HIST(("MC/reconstructed/" + signalFolder + "/hCPARecSig" + signalSuffix)).c_str(), cpa);
          registry.fill(HIST(("MC/reconstructed/" + signalFolder + "/hCPAVsPtRecSig" + signalSuffix)).c_str(), cpa, pt);
          registry.fill(HIST(("MC/reconstructed/" + signalFolder + "/hCPAxyRecSig" + signalSuffix)).c_str(), cpaXY);
          registry.fill(HIST(("MC/reconstructed/" + signalFolder + "/hCPAxyVsPtRecSig" + signalSuffix)).c_str(), cpaXY, pt);
          registry.fill(HIST(("MC/reconstructed/" + signalFolder + "/hDca2RecSig" + signalSuffix)).c_str(), chi2PCA);
          registry.fill(HIST(("MC/reconstructed/" + signalFolder + "/hDca2VsPtRecSig" + signalSuffix)).c_str(), chi2PCA, pt);
          registry.fill(HIST(("MC/reconstructed/" + signalFolder + "/hEtaRecSig" + signalSuffix)).c_str(), candidate.eta());
          registry.fill(HIST(("MC/reconstructed/" + signalFolder + "/hEtaVsPtRecSig" + signalSuffix)).c_str(), candidate.eta(), pt);
          registry.fill(HIST(("MC/reconstructed/" + signalFolder + "/hPhiRecSig" + signalSuffix)).c_str(), candidate.phi());
          registry.fill(HIST(("MC/reconstructed/" + signalFolder + "/hPhiVsPtRecSig" + signalSuffix)).c_str(), candidate.phi(), pt);
          registry.fill(HIST(("MC/reconstructed/" + signalFolder + "/hImpParErrProng0VsPtRecSig" + signalSuffix)).c_str(), candidate.errorImpactParameter0(), pt);
          registry.fill(HIST(("MC/reconstructed/" + signalFolder + "/hImpParErrProng1VsPtRecSig" + signalSuffix)).c_str(), candidate.errorImpactParameter1(), pt);
          registry.fill(HIST(("MC/reconstructed/" + signalFolder + "/hImpParErrProng2VsPtRecSig" + signalSuffix)).c_str(), candidate.errorImpactParameter2(), pt);
          registry.fill(HIST(("MC/reconstructed/" + signalFolder + "/hDecLenErrVsPtRecSig" + signalSuffix)).c_str(), candidate.errorDecayLength(), pt);
        };

        /// MC reconstructed signal
        fillHistogramsRecSig("signal", "");
        /// reconstructed signal prompt
        if (candidate.originMcRec() == RecoDecay::OriginType::Prompt) {
          fillHistogramsRecSig("prompt", "Prompt");
          /// reconstructed signal nonprompt
        } else if (candidate.originMcRec() == RecoDecay::OriginType::NonPrompt) {
          fillHistogramsRecSig("nonprompt", "NonPrompt");
        }

        if (fillTHn) {
          float cent = evaluateCentralityColl(collision);
          float occ{-1.};
          if (storeOccupancy && occEstimator != o2::hf_occupancy::OccupancyEstimator::None) {
            occ = o2::hf_occupancy::getOccupancyColl(collision, occEstimator);
          }
          double massLc(-1);
          double outputBkg(-1), outputPrompt(-1), outputFD(-1);
          const float properLifetime = hfHelper.ctLc(candidate) * CtToProperLifetimePs;
          if ((candidate.isSelLcToPKPi() >= selectionFlagLc) && pdgCodeProng0 == kProton) {
            massLc = hfHelper.invMassLcToPKPi(candidate);

            if constexpr (fillMl) {
              if (candidate.mlProbLcToPKPi().size() == NumberOfMlClasses) {
                outputBkg = candidate.mlProbLcToPKPi()[MlClassBackground]; /// bkg score
                outputPrompt = candidate.mlProbLcToPKPi()[MlClassPrompt];  /// prompt score
                outputFD = candidate.mlProbLcToPKPi()[MlClassNonPrompt];   /// non-prompt score
              }
              /// Fill the ML outputScores and variables of candidate
              std::vector<double> valuesToFill{massLc, pt, cent, outputBkg, outputPrompt, outputFD, static_cast<double>(numPvContributors), ptRecB, static_cast<double>(originType)};
              if (storeOccupancy && occEstimator != o2::hf_occupancy::OccupancyEstimator::None) {
                valuesToFill.push_back(occ);
              }
              if (storeProperLifetime) {
                valuesToFill.push_back(properLifetime);
              }
              registry.get<THnSparse>(HIST("hnLcVarsWithBdt"))->Fill(valuesToFill.data());
            } else {
              std::vector<double> valuesToFill{massLc, pt, cent, ptProng0, ptProng1, ptProng2, chi2PCA, decayLength, cpa, static_cast<double>(numPvContributors), ptRecB, static_cast<double>(originType)};
              if (storeOccupancy && occEstimator != o2::hf_occupancy::OccupancyEstimator::None) {
                valuesToFill.push_back(occ);
              }
              if (storeProperLifetime) {
                valuesToFill.push_back(properLifetime);
              }
              registry.get<THnSparse>(HIST("hnLcVars"))->Fill(valuesToFill.data());
            }
          }
          if ((candidate.isSelLcToPiKP() >= selectionFlagLc) && pdgCodeProng0 == kPiPlus) {
            massLc = hfHelper.invMassLcToPiKP(candidate);

            if constexpr (fillMl) {
              if (candidate.mlProbLcToPiKP().size() == NumberOfMlClasses) {
                outputBkg = candidate.mlProbLcToPiKP()[MlClassBackground]; /// bkg score
                outputPrompt = candidate.mlProbLcToPiKP()[MlClassPrompt];  /// prompt score
                outputFD = candidate.mlProbLcToPiKP()[MlClassNonPrompt];   /// non-prompt score
              }
              /// Fill the ML outputScores and variables of candidate (todo: add multiplicity)
              std::vector<double> valuesToFill{massLc, pt, cent, outputBkg, outputPrompt, outputFD, static_cast<double>(numPvContributors), ptRecB, static_cast<double>(originType)};
              if (storeOccupancy && occEstimator != o2::hf_occupancy::OccupancyEstimator::None) {
                valuesToFill.push_back(occ);
              }
              if (storeProperLifetime) {
                valuesToFill.push_back(properLifetime);
              }
              registry.get<THnSparse>(HIST("hnLcVarsWithBdt"))->Fill(valuesToFill.data());
            } else {
              std::vector<double> valuesToFill{massLc, pt, cent, ptProng0, ptProng1, ptProng2, chi2PCA, decayLength, cpa, static_cast<double>(numPvContributors), ptRecB, static_cast<double>(originType)};
              if (storeOccupancy && occEstimator != o2::hf_occupancy::OccupancyEstimator::None) {
                valuesToFill.push_back(occ);
              }
              if (storeProperLifetime) {
                valuesToFill.push_back(properLifetime);
              }
              registry.get<THnSparse>(HIST("hnLcVars"))->Fill(valuesToFill.data());
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
      if (std::abs(particle.flagMcMatchGen()) == hf_decay::hf_cand_3prong::DecayChannelMain::LcToPKPi) {
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

        const auto mcDaughter0 = particle.template daughters_as<soa::Join<aod::McParticles, aod::HfCand3ProngMcGen>>().begin();
        const float p2m = particle.p() / o2::constants::physics::MassLambdaCPlus;
        const float gamma = std::sqrt(1 + p2m * p2m);                       // mother's particle Lorentz factor
        const float properLifetime = mcDaughter0.vt() * NanoToPico / gamma; // from ns to ps * from lab time to proper time

        registry.fill(HIST("MC/generated/signal/hPtGen"), ptGen);
        registry.fill(HIST("MC/generated/signal/hEtaGen"), particle.eta());
        registry.fill(HIST("MC/generated/signal/hYGen"), yGen);
        registry.fill(HIST("MC/generated/signal/hPhiGen"), particle.phi());
        registry.fill(HIST("MC/generated/signal/hEtaVsPtGen"), particle.eta(), ptGen);
        registry.fill(HIST("MC/generated/signal/hYVsPtGen"), yGen, ptGen);
        registry.fill(HIST("MC/generated/signal/hPhiVsPtGen"), particle.phi(), ptGen);

        if (particle.originMcGen() == RecoDecay::OriginType::Prompt) {
          if (fillTHn) {
            std::vector<double> valuesToFill{ptGen, cent, yGen, static_cast<double>(numPvContributors), ptGenB, static_cast<double>(originType)};
            if (storeOccupancy && occEstimator != o2::hf_occupancy::OccupancyEstimator::None) {
              valuesToFill.push_back(occ);
            }
            if (storeProperLifetime) {
              valuesToFill.push_back(properLifetime);
            }
            registry.get<THnSparse>(HIST("hnLcVarsGen"))->Fill(valuesToFill.data());
          }
          registry.fill(HIST("MC/generated/prompt/hPtGenPrompt"), ptGen);
          registry.fill(HIST("MC/generated/prompt/hEtaGenPrompt"), particle.eta());
          registry.fill(HIST("MC/generated/prompt/hYGenPrompt"), yGen);
          registry.fill(HIST("MC/generated/prompt/hPhiGenPrompt"), particle.phi());
          registry.fill(HIST("MC/generated/prompt/hEtaVsPtGenPrompt"), particle.eta(), ptGen);
          registry.fill(HIST("MC/generated/prompt/hYVsPtGenPrompt"), yGen, ptGen);
          registry.fill(HIST("MC/generated/prompt/hPhiVsPtGenPrompt"), particle.phi(), ptGen);
        }
        if (particle.originMcGen() == RecoDecay::OriginType::NonPrompt) {
          ptGenB = mcParticles.rawIteratorAt(particle.idxBhadMotherPart()).pt();
          if (fillTHn) {
            std::vector<double> valuesToFill{ptGen, cent, yGen, static_cast<double>(numPvContributors), ptGenB, static_cast<double>(originType)};
            if (storeOccupancy && occEstimator != o2::hf_occupancy::OccupancyEstimator::None) {
              valuesToFill.push_back(occ);
            }
            if (storeProperLifetime) {
              valuesToFill.push_back(properLifetime);
            }
            registry.get<THnSparse>(HIST("hnLcVarsGen"))->Fill(valuesToFill.data());
          }
          registry.fill(HIST("MC/generated/nonprompt/hPtGenNonPrompt"), ptGen);
          registry.fill(HIST("MC/generated/nonprompt/hEtaGenNonPrompt"), particle.eta());
          registry.fill(HIST("MC/generated/nonprompt/hYGenNonPrompt"), yGen);
          registry.fill(HIST("MC/generated/nonprompt/hPhiGenNonPrompt"), particle.phi());
          registry.fill(HIST("MC/generated/nonprompt/hEtaVsPtGenNonPrompt"), particle.eta(), ptGen);
          registry.fill(HIST("MC/generated/nonprompt/hYVsPtGenNonPrompt"), yGen, ptGen);
          registry.fill(HIST("MC/generated/nonprompt/hPhiVsPtGenNonPrompt"), particle.phi(), ptGen);
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
      registry.fill(HIST("Data/hImpParErrProng0VsPt"), candidate.errorImpactParameter0(), pt);
      registry.fill(HIST("Data/hImpParErrProng1VsPt"), candidate.errorImpactParameter1(), pt);
      registry.fill(HIST("Data/hImpParErrProng2VsPt"), candidate.errorImpactParameter2(), pt);
      registry.fill(HIST("Data/hDecLenErrVsPt"), candidate.errorDecayLength(), pt);

      if (fillTHn) {
        float cent = evaluateCentralityColl(collision);
        float occ{-1.};
        if (storeOccupancy && occEstimator != o2::hf_occupancy::OccupancyEstimator::None) {
          occ = o2::hf_occupancy::getOccupancyColl(collision, occEstimator);
        }
        double massLc(-1);
        double outputBkg(-1), outputPrompt(-1), outputFD(-1);
        const float properLifetime = hfHelper.ctLc(candidate) * CtToProperLifetimePs;
        if (candidate.isSelLcToPKPi() >= selectionFlagLc) {
          massLc = hfHelper.invMassLcToPKPi(candidate);

          if constexpr (fillMl) {
            if (candidate.mlProbLcToPKPi().size() == NumberOfMlClasses) {
              outputBkg = candidate.mlProbLcToPKPi()[MlClassBackground]; /// bkg score
              outputPrompt = candidate.mlProbLcToPKPi()[MlClassPrompt];  /// prompt score
              outputFD = candidate.mlProbLcToPKPi()[MlClassNonPrompt];   /// non-prompt score
            }
            /// Fill the ML outputScores and variables of candidate
            std::vector<double> valuesToFill{massLc, pt, cent, outputBkg, outputPrompt, outputFD, static_cast<double>(numPvContributors)};
            if (storeOccupancy && occEstimator != o2::hf_occupancy::OccupancyEstimator::None) {
              valuesToFill.push_back(occ);
            }
            if (storeProperLifetime) {
              valuesToFill.push_back(properLifetime);
            }
            registry.get<THnSparse>(HIST("hnLcVarsWithBdt"))->Fill(valuesToFill.data());
          } else {
            std::vector<double> valuesToFill{massLc, pt, cent, ptProng0, ptProng1, ptProng2, chi2PCA, decayLength, cpa, static_cast<double>(numPvContributors)};
            if (storeOccupancy && occEstimator != o2::hf_occupancy::OccupancyEstimator::None) {
              valuesToFill.push_back(occ);
            }
            if (storeProperLifetime) {
              valuesToFill.push_back(properLifetime);
            }
            registry.get<THnSparse>(HIST("hnLcVars"))->Fill(valuesToFill.data());
          }
        }
        if (candidate.isSelLcToPiKP() >= selectionFlagLc) {
          massLc = hfHelper.invMassLcToPiKP(candidate);

          if constexpr (fillMl) {
            if (candidate.mlProbLcToPiKP().size() == NumberOfMlClasses) {
              outputBkg = candidate.mlProbLcToPiKP()[MlClassBackground]; /// bkg score
              outputPrompt = candidate.mlProbLcToPiKP()[MlClassPrompt];  /// prompt score
              outputFD = candidate.mlProbLcToPiKP()[MlClassNonPrompt];   /// non-prompt score
            }
            /// Fill the ML outputScores and variables of candidate
            std::vector<double> valuesToFill{massLc, pt, cent, outputBkg, outputPrompt, outputFD, static_cast<double>(numPvContributors)};
            if (storeOccupancy && occEstimator != o2::hf_occupancy::OccupancyEstimator::None) {
              valuesToFill.push_back(occ);
            }
            if (storeProperLifetime) {
              valuesToFill.push_back(properLifetime);
            }
            registry.get<THnSparse>(HIST("hnLcVarsWithBdt"))->Fill(valuesToFill.data());
          } else {
            std::vector<double> valuesToFill{massLc, pt, cent, ptProng0, ptProng1, ptProng2, chi2PCA, decayLength, cpa, static_cast<double>(numPvContributors)};
            if (storeOccupancy && occEstimator != o2::hf_occupancy::OccupancyEstimator::None) {
              valuesToFill.push_back(occ);
            }
            if (storeProperLifetime) {
              valuesToFill.push_back(properLifetime);
            }
            registry.get<THnSparse>(HIST("hnLcVars"))->Fill(valuesToFill.data());
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
