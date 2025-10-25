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

/// \file candidateCreatorSigmac0plusplusCascade.cxx
/// \brief Σc0,++ → Λc(→ K0sP) + π-,+ candidate builder
/// \note Here the Lc from the cascade channel is obtained using the task taskLcToK0sP.cxx
/// \author Rutuparna Rath <rrath@cern.ch>, INFN BOLOGNA and GSI Darmstadt
/// In collaboration with Andrea Alici <aalici@cern.ch>, INFN BOLOGNA

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/AliasTables.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>

#include <cstdint>
#include <vector>

using namespace o2;
using namespace o2::analysis;
using namespace o2::constants::physics;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct HfCandidateCreatorSigmac0plusplusCascade {

  /// Table with Σc0,++ info
  Produces<aod::HfCandScCasBase> rowScCandidates;

  Configurable<float> trkMinPt{"trkMinPt", 0.15, "track min pT"};
  Configurable<float> trkMaxEta{"trkMaxEta", 0.8, "track max Eta"};
  Configurable<float> maxDCAxyToPVcut{"maxDCAxyToPVcut", 2.0, "Track DCAxy cut to PV Maximum"};
  Configurable<float> maxDCAzToPVcut{"maxDCAzToPVcut", 2.0, "Track DCAz cut to PV Maximum"};
  Configurable<float> nTpcNClsFound{"nTpcNClsFound", 120, "nFindable TPC Clusters"};
  Configurable<float> nTPCCrossedRows{"nTPCCrossedRows", 70, "nCrossed TPC Rows"};
  Configurable<float> nTPCChi2{"nTPCChi2", 4.0, "nTPC Chi2 per Cluster"};
  Configurable<float> nITSChi2{"nITSChi2", 36.0, "nITS Chi2 per Cluster"};
  Configurable<float> tpcnSigmaPi{"tpcnSigmaPi", 3.0, "TPC nSigma selection"};

  /// Selection of candidates Λc+
  Configurable<int> selectionFlagLc{"selectionFlagLc", 1, "Selection Flag for Lc"};
  Configurable<double> yCandLcMax{"yCandLcMax", -1., "max. candLc. Lc rapidity"};
  Configurable<int> selectionFlagLcToK0sP{"selectionFlagLcToK0sP", 1, "Selection Flag for Lc"};
  Configurable<int> selectionFlagLcbarToK0sP{"selectionFlagLcbarToK0sP", 1, "Selection Flag for Lcbar"};
  Configurable<double> cutsMassLcMax{"cutsMassLcMax", 0.08, "Lc candidate mass selection"};
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_sigmac_to_p_k_pi::vecBinsPt}, "pT bin limits"};

  /// Selections on candidate soft π-,+
  Configurable<bool> applyGlobalTrkWoDcaCutsSoftPi{"applyGlobalTrkWoDcaCutsSoftPi", false, "Switch on the application of the global-track w/o dca cuts for soft pion BEFORE ALL OTHER CUSTOM CUTS"};
  Configurable<float> softPiEtaMax{"softPiEtaMax", 0.9f, "Soft pion max value for pseudorapidity (abs vale)"};
  Configurable<float> softPiChi2Max{"softPiChi2Max", 36.f, "Soft pion max value for chi2 ITS"};
  Configurable<int> softPiItsHitMap{"softPiItsHitMap", 127, "Soft pion ITS hitmap"};
  Configurable<int> softPiItsHitsMin{"softPiItsHitsMin", 1, "Minimum number of ITS layers crossed by the soft pion among those in \"softPiItsHitMap\""};
  Configurable<float> softPiDcaXYMax{"softPiDcaXYMax", 0.065, "Soft pion max dcaXY (cm)"};
  Configurable<float> softPiDcaZMax{"softPiDcaZMax", 0.065, "Soft pion max dcaZ (cm)"};
  Configurable<bool> addQA{"addQA", true, "Switch for the qa PLOTS"};

  HfHelper hfHelper;

  using TracksWithPID = soa::Join<aod::TracksWDcaExtra, aod::TracksPidPi, aod::TracksPidKa>;

  /// Filter the candidate Λc+ used for the Σc0,++ creation
  Filter filterSelectCandidateLc = (aod::hf_sel_candidate_lc_to_k0s_p::isSelLcToK0sP >= selectionFlagLcToK0sP ||
                                    aod::hf_sel_candidate_lc_to_k0s_p::isSelLcToK0sP >= selectionFlagLcbarToK0sP);

  // slice by hand the assoc. track with the  Λc+ collisionId
  Preslice<TracksWithPID> trackIndicesPerCollision = aod::track::collisionId;

  HistogramRegistry registry;

  void init(InitContext&)
  {
    // axes
    AxisSpec const axisBinsPt = {binsPt, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec const axisPt = {300, 0.0f, 30.0f, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec const axisEta = {500, -2.0f, 2.0f, "#it{#eta}"};
    AxisSpec const axisPhi = {100, 0.f, 6.3f, "#it{#phi}"};
    AxisSpec const axisMassCand = {600, 1.98f, 2.58f, "inv. mass (p K_{S}^{0}) (GeV/#it{c}^{2})"};
    AxisSpec const axisd0 = {500, -0.5f, 0.5f, "DCAxy (cm)"};
    AxisSpec const axisd0V0Daughters = {1000, -5.0f, 5.0f, "DCAxy (cm)"};
    AxisSpec const axisV0CPA = {500, 0.98f, 1.0001f, "v0 cos pointing angle"};
    AxisSpec const axisV0Radius = {1000, 0.f, 40.f, "V0 radius (cm)"};
    AxisSpec const axisV0DCADaughters = {200, 0.f, 2.f, "DCA (cm)"};
    AxisSpec const axisMassK0Short = {500, 0.4f, 0.6f, "#it{m}(K_{S}^{0}) (GeV/#it{c}^{2})"};
    AxisSpec const axisMassLambda = {500, 1.0f, 1.2f, "#it{m}(#Lambda) (GeV/#it{c}^{2})"};
    AxisSpec const axisMassGamma = {500, 0.0f, 0.4f, "#it{m}(#gamma) (GeV/#it{c}^{2})"};
    AxisSpec const axisCPACand = {110, -1.1f, 1.1f, "candiate cos pointing angle"};
    AxisSpec const axisDecLength = {200, 0.f, 2.0f, "decay length (cm)"};
    AxisSpec const axisProperLifetime = {100, 0.f, 0.2f, "#it{c#tau} (cm)"};
    AxisSpec const axisProperLifetimeV0 = {1000, 0.f, 80.f, "#it{c#tau} (cm)"};
    AxisSpec const axisNSigma = {100, -6.f, 6.f, "n#it{#sigma}_{p}"};
    AxisSpec const axisPidP = {100, 0.f, 10.0f, "#it{p} (GeV/#it{c})"};

    auto h = registry.add<TH1>("candidateStat", "", kTH1D, {{3, 0.5, 3.5}});
    h->GetXaxis()->SetBinLabel(1, "Lc candidates");
    h->GetXaxis()->SetBinLabel(2, "soft #pi (before cuts)");
    h->GetXaxis()->SetBinLabel(3, "soft #pi (after track cuts)");
    // data
    if (addQA) {
      registry.add("lc/hPtCand", "cascade candidates;candidateLc #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPt}});
      registry.add("lc/hEtaCand", "cascade candidates;candidateLc #it{#eta};entries", {HistType::kTH1F, {axisEta}});
      registry.add("lc/hEtaCandVsPtCand", "cascade candidates;candidateLc #it{#eta};p_{T}", {HistType::kTH2F, {axisEta, axisBinsPt}});
      registry.add("lc/hPhiCand", "cascade candidates;candidateLc #it{#phi};entries", {HistType::kTH1F, {axisPhi}});
      registry.add("lc/hPhiCandVsPtCand", "cascade candidates;candidateLc #it{#phi};p_{T}", {HistType::kTH2F, {axisPhi, axisBinsPt}});
      registry.add("lc/hMass", "cascade candidates;inv. mass (p K_{S}^{0}) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {axisMassCand}});
      registry.add("lc/hMassVsPtCand", "cascade candidates;inv. mass (p K_{S}^{0}) (GeV/#it{c}^{2});p_{T}", {HistType::kTH2F, {axisMassCand, axisBinsPt}});
      registry.add("lc/hPtBach", "cascade candidates;bachelor #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPt}});
      registry.add("lc/hPtBachVsPtCand", "cascade candidates;bachelor #it{p}_{T} (GeV/#it{c});p_{T}", {HistType::kTH2F, {axisPt, axisBinsPt}});
      registry.add("lc/hPtV0", "cascade candidates;v0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPt}});
      registry.add("lc/hPtV0VsPtCand", "cascade candidates;v0 #it{p}_{T} (GeV/#it{c});p_{T}", {HistType::kTH2F, {axisPt, axisBinsPt}});
      registry.add("lc/hd0Bach", "cascade candidates;bachelor DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {axisd0}});
      registry.add("lc/hd0BachVsPtCand", "cascade candidates;bachelor DCAxy to prim. vertex (cm);p_{T}", {HistType::kTH2F, {axisd0, axisBinsPt}});
      registry.add("lc/hd0V0", "cascade candidates;V0 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {axisd0}});
      registry.add("lc/hd0V0VsPtCand", "cascade candidates;V0 DCAxy to prim. vertex (cm);p_{T}", {HistType::kTH2F, {axisd0, axisBinsPt}});
      registry.add("lc/hd0V0pos", "cascade candidates;pos daugh v0 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {axisd0V0Daughters}});
      registry.add("lc/hd0V0posVsPtCand", "cascade candidates;pos daugh v0 DCAxy to prim. vertex (cm);p_{T}", {HistType::kTH2F, {axisd0V0Daughters, axisBinsPt}});
      registry.add("lc/hd0V0neg", "cascade candidates;neg daugh v0 DCAxy to prim. vertex (cm);entries", {HistType::kTH1F, {axisd0V0Daughters}});
      registry.add("lc/hd0V0negVsPtCand", "cascade candidates;neg daugh v0 DCAxy to prim. vertex (cm);p_{T}", {HistType::kTH2F, {axisd0V0Daughters, axisBinsPt}});
      registry.add("lc/hPtV0pos", "cascade candidates;pos daugh v0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPt}});
      registry.add("lc/hPtV0posVsPtCand", "cascade candidates;pos daugh v0 #it{p}_{T} (GeV/#it{c});p_{T}", {HistType::kTH2F, {axisPt, axisBinsPt}});
      registry.add("lc/hPtV0neg", "cascade candidates;neg daugh v0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {axisPt}});
      registry.add("lc/hPtV0negVsPtCand", "cascade candidates;neg daugh v0 #it{p}_{T} (GeV/#it{c});p_{T}", {HistType::kTH2F, {axisPt, axisBinsPt}});
      registry.add("lc/hV0CPA", "cascade candidates;v0 cosine of pointing angle;entries", {HistType::kTH1F, {axisV0CPA}});
      registry.add("lc/hV0CPAVsPtCand", "cascade candidates;v0 cosine of pointing angle;p_{T}", {HistType::kTH2F, {axisV0CPA, axisBinsPt}});
      registry.add("lc/hV0Radius", "cascade candidates;v0 radius (cm);entries", {HistType::kTH1F, {axisV0Radius}});
      registry.add("lc/hV0RadiusVsPtCand", "cascade candidates;v0 radius (cm);p_{T}", {HistType::kTH2F, {axisV0Radius, axisBinsPt}});
      registry.add("lc/hV0DCADaughters", "cascade candidates;v0 dca daughters (cm);entries", {HistType::kTH1F, {axisV0DCADaughters}});
      registry.add("lc/hV0DCADaughtersVsPtCand", "cascade candidates;v0 dca daughters (cm);p_{T}", {HistType::kTH2F, {axisV0DCADaughters, axisBinsPt}});
      registry.add("lc/hV0MK0Short", "cascade candidates;v0 mass K0s (GeV/#it{c}^{2});entries", {HistType::kTH1F, {axisMassK0Short}});
      registry.add("lc/hV0MK0ShortVsPtCand", "cascade candidates;v0 mass K0s (GeV/#it{c}^{2});p_{T}", {HistType::kTH2F, {axisMassK0Short, axisBinsPt}});
      registry.add("lc/hV0MLambda", "cascade candidates;v0 mass Lambda (GeV/#it{c}^{2});entries", {HistType::kTH1F, {axisMassLambda}});
      registry.add("lc/hV0MLambdaVsPtCand", "cascade candidates;v0 mass Lambda (GeV/#it{c}^{2});p_{T}", {HistType::kTH2F, {axisMassLambda, axisBinsPt}});
      registry.add("lc/hV0MAntiLambda", "cascade candidates;v0 mass AntiLambda (GeV/#it{c}^{2});entries", {HistType::kTH1F, {axisMassLambda}});
      registry.add("lc/hV0MAntiLambdaVsPtCand", "cascade candidates;v0 mass AntiLambda (GeV/#it{c}^{2});p_{T}", {HistType::kTH2F, {axisMassLambda, axisBinsPt}});
      registry.add("lc/hV0MGamma", "cascade candidates;v0 mass Gamma (GeV/#it{c}^{2});entries", {HistType::kTH1F, {axisMassGamma}});
      registry.add("lc/hV0MGammaVsPtCand", "cascade candidates;v0 mass Gamma (GeV/#it{c}^{2});p_{T}", {HistType::kTH2F, {axisMassGamma, axisBinsPt}});
      registry.add("lc/hCtV0K0Short", "cascade candidates;proper lifetime (V0) * #it{c} (cm);entries", {HistType::kTH1F, {axisProperLifetimeV0}});
      registry.add("lc/hCtV0K0ShortVsPtCand", "cascade candidates;proper lifetime (V0) * #it{c} (cm);p_{T}", {HistType::kTH2F, {axisProperLifetimeV0, axisBinsPt}});
      registry.add("lc/hCtV0Lambda", "cascade candidates;proper lifetime (V0) * #it{c} (cm);entries", {HistType::kTH1F, {axisProperLifetimeV0}});
      registry.add("lc/hCtV0LambdaVsPtCand", "cascade candidates;proper lifetime (V0) * #it{c} (cm);p_{T}", {HistType::kTH2F, {axisProperLifetimeV0, axisBinsPt}});
      registry.add("lc/hCPACand", "cascade candidates;cosine pointing angle;entries", {HistType::kTH1F, {axisCPACand}});
      registry.add("lc/hCPACandVsPtCand", "cascade candidates;cosine pointing angle;p_{T}", {HistType::kTH2F, {axisCPACand, axisBinsPt}});
      registry.add("lc/hCPAxyCand", "cascade candidates;cosine pointing angle xy;entries", {HistType::kTH1F, {axisCPACand}});
      registry.add("lc/hCPAxyCandVsPtCand", "cascade candidates;cosine pointing angle xy;p_{T}", {HistType::kTH2F, {axisCPACand, axisBinsPt}});
      registry.add("lc/hDecLengthCand", "cascade candidates;decay length (cm);entries", {HistType::kTH1F, {axisDecLength}});
      registry.add("lc/hDecLengthCandVsPtCand", "cascade candidates;decay length (cm);p_{T}", {HistType::kTH2F, {axisDecLength, axisBinsPt}});
      registry.add("lc/hDecLengthXYCand", "cascade candidates;decay length xy (cm);entries", {HistType::kTH1F, {axisDecLength}});
      registry.add("lc/hDecLengthXYCandVsPtCand", "cascade candidates;decay length xy (cm);p_{T}", {HistType::kTH2F, {axisDecLength, axisBinsPt}});
      registry.add("lc/hCtCand", "cascade candidates;proper lifetime (#Lambda_{c}) * #it{c} (cm);entries", {HistType::kTH1F, {axisProperLifetime}});
      registry.add("lc/hCtCandVsPtCand", "cascade candidates;proper lifetime (#Lambda_{c}) * #it{c} (cm);p_{T}", {HistType::kTH2F, {axisProperLifetime, axisBinsPt}});

      // soft pion
      registry.add("pion/data/hPtSoftPi", "#pi candidates; #it{p}_{T}(#pi) (GeV/#it{c}); counts;", {HistType::kTH1F, {axisPt}});
      registry.add("pion/data/hEtaSoftPi", "#pi candidates; #eta ; counts;", {HistType::kTH1F, {axisEta}});
      registry.add("pion/data/hPhiSoftPi", "#pi candidates; #Phi ; counts;", {HistType::kTH1F, {axisPhi}});
    }
  }

  template <typename TrackType>
  bool isTrackSelected(const TrackType& track)
  {
    if (track.pt() < trkMinPt) {
      return false;
    }
    if (std::abs(track.eta()) > trkMaxEta) {
      return false;
    }
    if (std::abs(track.dcaXY()) > maxDCAxyToPVcut) {
      return false;
    }
    if (std::abs(track.dcaZ()) > maxDCAzToPVcut) {
      return false;
    }
    if (track.tpcNClsFound() < nTpcNClsFound) {
      return false;
    }
    if (track.tpcNClsCrossedRows() < nTPCCrossedRows) {
      return false;
    }
    if (track.tpcChi2NCl() > nTPCChi2) {
      return false;
    }
    if (track.itsChi2NCl() > nITSChi2) {
      return false;
    }
    if (track.tpcNSigmaPi() > tpcnSigmaPi) {
      return false;
    }

    return true;
  }

  /// @param tracks are the tracks (with dcaXY, dcaZ information) → soft-pion candidate tracks
  /// @param candidatesLc are 2-prong candidates satisfying the analysis selections for Λc+ → Ks0P (and charge conj.)
  void processData(soa::Filtered<soa::Join<aod::HfCandCascExt,
                                           aod::HfSelLcToK0sP>> const& candidates,
                   aod::Collisions const&,
                   TracksWithPID const& tracks,
                   aod::V0s const&)
  {
    for (const auto& candidateLc : candidates) {
      /// slice the tracks based on the collisionId
      const auto& tracksInThisCollision = tracks.sliceBy(trackIndicesPerCollision, candidateLc.collision().globalIndex());
      const auto& bachProton = candidateLc.prong0_as<TracksWithPID>();
      int chargeLc = bachProton.sign(); // Lc charge depends on its bach charge (here it is proton)

      auto ptCand = candidateLc.pt();
      auto eta = candidateLc.eta();
      auto phi = candidateLc.phi();
      auto invMassLcToK0sP = hfHelper.invMassLcToK0sP(candidateLc);
      auto ptProng0 = candidateLc.ptProng0();
      auto ptProng1 = candidateLc.ptProng1();
      auto impactParameter0 = candidateLc.impactParameter0();
      auto impactParameter1 = candidateLc.impactParameter1();
      auto dcaPosToPV = candidateLc.dcapostopv();
      auto dcaNegToPV = candidateLc.dcanegtopv();
      auto ptV0Pos = candidateLc.ptV0Pos();
      auto ptV0Neg = candidateLc.ptV0Neg();
      auto v0CosPA = candidateLc.v0cosPA();
      auto v0Radius = candidateLc.v0radius();
      auto dcaV0Daughters = candidateLc.dcaV0daughters();
      auto mK0Short = candidateLc.mK0Short();
      auto mLambda = candidateLc.mLambda();
      auto mAntiLambda = candidateLc.mAntiLambda();
      auto mGamma = candidateLc.mGamma();
      auto ctV0K0Short = hfHelper.ctV0K0s(candidateLc);
      auto ctV0Lambda = hfHelper.ctV0Lambda(candidateLc);
      auto cpa = candidateLc.cpa();
      auto cpaXY = candidateLc.cpaXY();
      auto decayLength = candidateLc.decayLength();
      auto decayLengthXY = candidateLc.decayLengthXY();
      auto ctLc = hfHelper.ctLc(candidateLc);
      if (addQA) {
        registry.fill(HIST("lc/hPtCand"), ptCand);
        registry.fill(HIST("lc/hEtaCand"), eta);
        registry.fill(HIST("lc/hEtaCandVsPtCand"), eta, ptCand);
        registry.fill(HIST("lc/hPhiCand"), phi);
        registry.fill(HIST("lc/hPhiCandVsPtCand"), phi, ptCand);
        registry.fill(HIST("lc/hMass"), invMassLcToK0sP);
        registry.fill(HIST("lc/hMassVsPtCand"), invMassLcToK0sP, ptCand);
        registry.fill(HIST("lc/hPtBach"), ptProng0);
        registry.fill(HIST("lc/hPtBachVsPtCand"), ptProng0, ptCand);
        registry.fill(HIST("lc/hPtV0"), ptProng1);
        registry.fill(HIST("lc/hPtV0VsPtCand"), ptProng1, ptCand);
        registry.fill(HIST("lc/hd0Bach"), impactParameter0);
        registry.fill(HIST("lc/hd0BachVsPtCand"), impactParameter0, ptCand);
        registry.fill(HIST("lc/hd0V0"), impactParameter1);
        registry.fill(HIST("lc/hd0V0VsPtCand"), impactParameter1, ptCand);
        registry.fill(HIST("lc/hd0V0pos"), dcaPosToPV);
        registry.fill(HIST("lc/hd0V0posVsPtCand"), dcaPosToPV, ptCand);
        registry.fill(HIST("lc/hd0V0neg"), dcaNegToPV);
        registry.fill(HIST("lc/hd0V0negVsPtCand"), dcaNegToPV, ptCand);
        registry.fill(HIST("lc/hPtV0pos"), ptV0Pos);
        registry.fill(HIST("lc/hPtV0posVsPtCand"), ptV0Pos, ptCand);
        registry.fill(HIST("lc/hPtV0neg"), ptV0Neg);
        registry.fill(HIST("lc/hPtV0negVsPtCand"), ptV0Neg, ptCand);
        registry.fill(HIST("lc/hV0CPA"), v0CosPA);
        registry.fill(HIST("lc/hV0CPAVsPtCand"), v0CosPA, ptCand);
        registry.fill(HIST("lc/hV0Radius"), v0Radius);
        registry.fill(HIST("lc/hV0RadiusVsPtCand"), v0Radius, ptCand);
        registry.fill(HIST("lc/hV0DCADaughters"), dcaV0Daughters);
        registry.fill(HIST("lc/hV0DCADaughtersVsPtCand"), dcaV0Daughters, ptCand);
        registry.fill(HIST("lc/hV0MK0Short"), mK0Short);
        registry.fill(HIST("lc/hV0MK0ShortVsPtCand"), mK0Short, ptCand);
        registry.fill(HIST("lc/hV0MLambda"), mLambda);
        registry.fill(HIST("lc/hV0MLambdaVsPtCand"), mLambda, ptCand);
        registry.fill(HIST("lc/hV0MAntiLambda"), mAntiLambda);
        registry.fill(HIST("lc/hV0MAntiLambdaVsPtCand"), mAntiLambda, ptCand);
        registry.fill(HIST("lc/hV0MGamma"), mGamma);
        registry.fill(HIST("lc/hV0MGammaVsPtCand"), mGamma, ptCand);
        registry.fill(HIST("lc/hCtV0K0Short"), ctV0K0Short);
        registry.fill(HIST("lc/hCtV0K0ShortVsPtCand"), ctV0K0Short, ptCand);
        registry.fill(HIST("lc/hCtV0Lambda"), ctV0Lambda);
        registry.fill(HIST("lc/hCtV0LambdaVsPtCand"), ctV0Lambda, ptCand);
        registry.fill(HIST("lc/hCPACand"), cpa);
        registry.fill(HIST("lc/hCPACandVsPtCand"), cpa, ptCand);
        registry.fill(HIST("lc/hCPAxyCand"), cpaXY);
        registry.fill(HIST("lc/hCPAxyCandVsPtCand"), cpaXY, ptCand);
        registry.fill(HIST("lc/hDecLengthCand"), decayLength);
        registry.fill(HIST("lc/hDecLengthCandVsPtCand"), decayLength, ptCand);
        registry.fill(HIST("lc/hDecLengthXYCand"), decayLengthXY);
        registry.fill(HIST("lc/hDecLengthXYCandVsPtCand"), decayLengthXY, ptCand);
        registry.fill(HIST("lc/hCtCand"), ctLc);
        registry.fill(HIST("lc/hCtCandVsPtCand"), ctLc, ptCand);
      }
      if (std::abs(invMassLcToK0sP - MassLambdaCPlus) > cutsMassLcMax) {
        continue;
      }
      registry.fill(HIST("candidateStat"), 1);
      auto k0Short = candidateLc.v0_as<o2::aod::V0s>(); // get the soft pions for the given collId
      auto pos = k0Short.template posTrack_as<TracksWithPID>();
      auto neg = k0Short.template negTrack_as<TracksWithPID>();
      for (const auto& trackSoftPi : tracksInThisCollision) {
        int chargeSoftPi = trackSoftPi.sign();
        if (chargeSoftPi == pos.sign() && trackSoftPi.globalIndex() == pos.globalIndex()) {
          continue;
        }
        if (chargeSoftPi == neg.sign() && trackSoftPi.globalIndex() == neg.globalIndex()) {
          continue;
        }
        if (chargeSoftPi == bachProton.sign() && trackSoftPi.globalIndex() == bachProton.globalIndex()) {
          continue;
        }
        registry.fill(HIST("candidateStat"), 2);
        if (!isTrackSelected(trackSoftPi)) {
          continue;
        }
        registry.fill(HIST("candidateStat"), 3);

        /// fill histograms for softpion
        if (addQA) {
          registry.fill(HIST("pion/data/hPtSoftPi"), trackSoftPi.pt());
          registry.fill(HIST("pion/data/hEtaSoftPi"), trackSoftPi.eta());
          registry.fill(HIST("pion/data/hPhiSoftPi"), trackSoftPi.phi()); // π ← Σc0
        }
        /// determine the Σc candidate charge
        int8_t chargeSigmac = chargeLc + chargeSoftPi;

        /// fill the Σc0,++ candidate table
        rowScCandidates(/* general columns */
                        candidateLc.collisionId(),
                        /* 2-prong specific columns */
                        candidateLc.px(), candidateLc.py(), candidateLc.pz(), // Lc info
                        trackSoftPi.px(), trackSoftPi.py(), trackSoftPi.pz(), // soft pion info
                        candidateLc.collision().globalIndex(), trackSoftPi.globalIndex(),
                        chargeLc,
                        chargeSoftPi,
                        // candLc.hfflag(),
                        /* Σc0,++ specific columns */
                        chargeSigmac);
      }
    } // SC candidate
  }
  PROCESS_SWITCH(HfCandidateCreatorSigmac0plusplusCascade, processData, "Process Data", true);
};
struct HfCandidateCreatorSigmac0plusplusCascadeExpressions {
  Spawns<aod::HfCandScCasExt> candidatesSigmac;
  void processMc(aod::Tracks const&) {}
  PROCESS_SWITCH(HfCandidateCreatorSigmac0plusplusCascadeExpressions, processMc, "Process MC tracks", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateCreatorSigmac0plusplusCascade>(cfgc),
    adaptAnalysisTask<HfCandidateCreatorSigmac0plusplusCascadeExpressions>(cfgc)};
}
