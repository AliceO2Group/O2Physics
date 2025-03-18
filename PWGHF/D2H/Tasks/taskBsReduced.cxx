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

/// \file taskBsReduced.cxx
/// \brief Bs → Ds- π+ → (K- K+ π-) π+ analysis task
///
/// \author Fabio Catalano <fabio.catalano@cern.ch>, CERN

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "Common/Core/RecoDecay.h"

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/D2H/DataModel/ReducedDataModel.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{
namespace hf_cand_bs_lite
{
DECLARE_SOA_COLUMN(PtProng0, ptProng0, float);                               //! Transverse momentum of prong0 (GeV/c)
DECLARE_SOA_COLUMN(PtProng1, ptProng1, float);                               //! Transverse momentum of prong1 (GeV/c)
DECLARE_SOA_COLUMN(MProng0, mProng0, float);                                 //! Invariant mass of prong0 (GeV/c)
DECLARE_SOA_COLUMN(M, m, float);                                             //! Invariant mass of candidate (GeV/c2)
DECLARE_SOA_COLUMN(Pt, pt, float);                                           //! Transverse momentum of candidate (GeV/c)
DECLARE_SOA_COLUMN(PtGen, ptGen, float);                                     //! Transverse momentum of candidate (GeV/c)
DECLARE_SOA_COLUMN(P, p, float);                                             //! Momentum of candidate (GeV/c)
DECLARE_SOA_COLUMN(Y, y, float);                                             //! Rapidity of candidate
DECLARE_SOA_COLUMN(Eta, eta, float);                                         //! Pseudorapidity of candidate
DECLARE_SOA_COLUMN(Phi, phi, float);                                         //! Azimuth angle of candidate
DECLARE_SOA_COLUMN(E, e, float);                                             //! Energy of candidate (GeV)
DECLARE_SOA_COLUMN(NSigTpcPi1, nSigTpcPi1, float);                           //! TPC Nsigma separation for prong1 with pion mass hypothesis
DECLARE_SOA_COLUMN(NSigTofPi1, nSigTofPi1, float);                           //! TOF Nsigma separation for prong1 with pion mass hypothesis
DECLARE_SOA_COLUMN(DecayLength, decayLength, float);                         //! Decay length of candidate (cm)
DECLARE_SOA_COLUMN(DecayLengthXY, decayLengthXY, float);                     //! Transverse decay length of candidate (cm)
DECLARE_SOA_COLUMN(DecayLengthNormalised, decayLengthNormalised, float);     //! Normalised decay length of candidate
DECLARE_SOA_COLUMN(DecayLengthXYNormalised, decayLengthXYNormalised, float); //! Normalised transverse decay length of candidate
DECLARE_SOA_COLUMN(ImpactParameterProduct, impactParameterProduct, float);   //! Impact parameter product of candidate
DECLARE_SOA_COLUMN(Cpa, cpa, float);                                         //! Cosine pointing angle of candidate
DECLARE_SOA_COLUMN(CpaXY, cpaXY, float);                                     //! Cosine pointing angle of candidate in transverse plane
DECLARE_SOA_COLUMN(MaxNormalisedDeltaIP, maxNormalisedDeltaIP, float);       //! Maximum normalized difference between measured and expected impact parameter of candidate prongs
DECLARE_SOA_COLUMN(MlScoreSig, mlScoreSig, float);                           //! ML score for signal class
DECLARE_SOA_COLUMN(FlagWrongCollision, flagWrongCollision, int8_t);          //! Flag for association with wrong collision
} // namespace hf_cand_bs_lite

DECLARE_SOA_TABLE(HfRedCandBsLites, "AOD", "HFREDCANDBSLITE", //! Table with some Bs properties
                  hf_cand::Chi2PCA,
                  hf_cand_bs_lite::DecayLength,
                  hf_cand_bs_lite::DecayLengthXY,
                  hf_cand_bs_lite::DecayLengthNormalised,
                  hf_cand_bs_lite::DecayLengthXYNormalised,
                  hf_cand_bs_lite::MProng0,
                  hf_cand_bs_lite::PtProng0,
                  hf_cand_bs_lite::PtProng1,
                  hf_cand::ImpactParameter0,
                  hf_cand::ImpactParameter1,
                  hf_cand_bs_lite::ImpactParameterProduct,
                  hf_cand_bs_lite::NSigTpcPi1,
                  hf_cand_bs_lite::NSigTofPi1,
                  hf_cand_bs_reduced::Prong0MlScoreBkg,
                  hf_cand_bs_reduced::Prong0MlScorePrompt,
                  hf_cand_bs_reduced::Prong0MlScoreNonprompt,
                  hf_cand_bs_lite::MlScoreSig,
                  hf_sel_candidate_bs::IsSelBsToDsPi,
                  hf_cand_bs_lite::M,
                  hf_cand_bs_lite::Pt,
                  hf_cand_bs_lite::Cpa,
                  hf_cand_bs_lite::CpaXY,
                  hf_cand_bs_lite::MaxNormalisedDeltaIP,
                  hf_cand_bs_lite::Eta,
                  hf_cand_bs_lite::Phi,
                  hf_cand_bs_lite::Y,
                  hf_cand_3prong::FlagMcMatchRec,
                  hf_cand_3prong::OriginMcRec,
                  hf_cand_bs_lite::FlagWrongCollision,
                  hf_cand_bs_lite::PtGen);

DECLARE_SOA_TABLE(HfRedBsMcCheck, "AOD", "HFREDBSMCCHECK", //! Table with MC decay type check
                  hf_cand_3prong::FlagMcMatchRec,
                  hf_cand_bs_lite::FlagWrongCollision,
                  hf_cand_bs_lite::MProng0,
                  hf_cand_bs_lite::PtProng0,
                  hf_cand_bs_lite::M,
                  hf_cand_bs_lite::Pt,
                  hf_cand_bs_lite::MlScoreSig,
                  hf_bs_mc::PdgCodeBeautyMother,
                  hf_bs_mc::PdgCodeCharmMother,
                  hf_bs_mc::PdgCodeProng0,
                  hf_bs_mc::PdgCodeProng1,
                  hf_bs_mc::PdgCodeProng2,
                  hf_bs_mc::PdgCodeProng3);
} // namespace o2::aod

/// Bs analysis task
struct HfTaskBsReduced {
  Produces<aod::HfRedCandBsLites> hfRedCandBsLite;
  Produces<aod::HfRedBsMcCheck> hfRedBsMcCheck;

  Configurable<int> selectionFlagBs{"selectionFlagBs", 1, "Selection Flag for Bs"};
  Configurable<float> yCandGenMax{"yCandGenMax", 0.5, "max. gen particle rapidity"};
  Configurable<float> yCandRecoMax{"yCandRecoMax", 0.8, "max. cand. rapidity"};
  Configurable<float> etaTrackMax{"etaTrackMax", 0.8, "max. track pseudo-rapidity for acceptance calculation"};
  Configurable<float> ptTrackMin{"ptTrackMin", 0.1, "min. track transverse momentum for acceptance calculation"};
  Configurable<bool> fillHistograms{"fillHistograms", true, "Flag to enable histogram filling"};
  Configurable<bool> fillSparses{"fillSparses", false, "Flag to enable sparse filling"};
  Configurable<bool> fillTree{"fillTree", false, "Flag to enable tree filling"};
  Configurable<bool> fillBackground{"fillBackground", false, "Flag to enable filling of background histograms/sparses/tree (only MC)"};
  Configurable<float> downSampleBkgFactor{"downSampleBkgFactor", 1., "Fraction of background_{s} candidates to keep for ML trainings"};
  Configurable<float> ptMaxForDownSample{"ptMaxForDownSample", 10., "Maximum pt for the application of the downsampling factor"};

  HfHelper hfHelper;

  Filter filterSelectCandidates = (aod::hf_sel_candidate_bs::isSelBsToDsPi >= selectionFlagBs);

  HistogramRegistry registry{"registry"};

  using TracksPion = soa::Join<HfRedTracks, HfRedTracksPid>;

  void init(InitContext&)
  {
    std::array<bool, 3> processFuncData{doprocessData, doprocessDataWithDmesMl, doprocessDataWithBsMl};
    if ((std::accumulate(processFuncData.begin(), processFuncData.end(), 0)) > 1) {
      LOGP(fatal, "Only one process function for data can be enabled at a time.");
    }
    std::array<bool, 6> processFuncMc{doprocessMc, doprocessMcWithDecayTypeCheck, doprocessMcWithDmesMl, doprocessMcWithDmesMlAndDecayTypeCheck, doprocessMcWithBsMl, doprocessMcWithBsMlAndDecayTypeCheck};
    if ((std::accumulate(processFuncMc.begin(), processFuncMc.end(), 0)) > 1) {
      LOGP(fatal, "Only one process function for MC can be enabled at a time.");
    }

    const AxisSpec axisMlScore{100, 0.f, 1.f};
    const AxisSpec axisMassBs{300, 4.5f, 6.0f};
    const AxisSpec axisMassDs{300, 1.75f, 2.05f};
    const AxisSpec axisDecayLength{200, 0.f, 0.4f};
    const AxisSpec axisNormDecayLength{100, 0.f, 50.f};
    const AxisSpec axisDca{100, -0.05f, 0.05f};
    const AxisSpec axisCosp{110, 0.f, 1.1f};
    const AxisSpec axisEta{30, -1.5f, 1.5f};
    const AxisSpec axisError{100, 0.f, 1.f};
    const AxisSpec axisImpParProd{100, -1.e-3, 1.e-3};
    const AxisSpec axisPtBs{100, 0.f, 50.f};
    const AxisSpec axisPtDs{100, 0.f, 50.f};
    const AxisSpec axisPtPi{100, 0.f, 10.f};

    if (doprocessData || doprocessDataWithDmesMl || doprocessDataWithBsMl) {
      if (fillHistograms) {
        registry.add("hMass", "B^{0}_{s} candidates;#it{p}_{T}(B^{0}_{s}) (GeV/#it{c});#it{M} (D_{s}#pi) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {axisPtBs, axisMassBs}});
        registry.add("hDecLength", "B^{0}_{s} candidates;#it{p}_{T}(B^{0}_{s}) (GeV/#it{c});B^{0}_{s} candidate decay length (cm);entries", {HistType::kTH2F, {axisPtBs, axisDecayLength}});
        registry.add("hDecLengthXy", "B^{0}_{s} candidates;#it{p}_{T}(B^{0}_{s}) (GeV/#it{c});B^{0}_{s} candidate decay length XY (cm);entries", {HistType::kTH2F, {axisPtBs, axisDecayLength}});
        registry.add("hNormDecLengthXy", "B^{0}_{s} candidates;#it{p}_{T}(B^{0}_{s}) (GeV/#it{c});B^{0}_{s} candidate norm. decay length XY (cm);entries", {HistType::kTH2F, {axisPtBs, axisNormDecayLength}});
        registry.add("hDcaProng0", "B^{0}_{s} candidates;#it{p}_{T}(B^{0}_{s}) (GeV/#it{c});prong 0 (D_{s}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {axisPtBs, axisDca}});
        registry.add("hDcaProng1", "B^{0}_{s} candidates;#it{p}_{T}(B^{0}_{s}) (GeV/#it{c});prong 1 (#pi) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {axisPtBs, axisDca}});
        registry.add("hPtProng0", "B^{0}_{s} candidates;#it{p}_{T}(B^{0}_{s}) (GeV/#it{c});#it{p}_{T}(D_{s}) (GeV/#it{c});entries", {HistType::kTH2F, {axisPtBs, axisPtDs}});
        registry.add("hPtProng1", "B^{0}_{s} candidates;#it{p}_{T}(B^{0}_{s}) (GeV/#it{c});#it{p}_{T}(#pi) (GeV/#it{c});entries", {HistType::kTH2F, {axisPtBs, axisPtPi}});
        registry.add("hCosp", "B^{0}_{s} candidates;#it{p}_{T}(B^{0}_{s}) (GeV/#it{c});B^{0}_{s} candidate cos(#vartheta_{P});entries", {HistType::kTH2F, {axisPtBs, axisCosp}});
        registry.add("hCospXy", "B^{0}_{s} candidates;#it{p}_{T}(B^{0}_{s}) (GeV/#it{c});B^{0}_{s} candidate cos(#vartheta_{P}^{XY});entries", {HistType::kTH2F, {axisPtBs, axisCosp}});
        registry.add("hEta", "B^{0}_{s} candidates;#it{p}_{T}(B^{0}_{s}) (GeV/#it{c});B^{0}_{s} candidate #it{#eta};entries", {HistType::kTH2F, {axisPtBs, axisEta}});
        registry.add("hRapidity", "B^{0}_{s} candidates;#it{p}_{T}(B^{0}_{s}) (GeV/#it{c});B^{0}_{s} candidate #it{y};entries", {HistType::kTH2F, {axisPtBs, axisEta}});
        registry.add("hImpParProd", "B^{0}_{s} candidates;#it{p}_{T}(B^{0}_{s}) (GeV/#it{c});B^{0}_{s} candidate impact parameter product;entries", {HistType::kTH2F, {axisPtBs, axisImpParProd}});
        registry.add("hInvMassD", "B^{0}_{s} candidates;#it{p}_{T}(D_{s}) (GeV/#it{c});prong0, #it{M}(KK#pi) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {axisPtDs, axisMassDs}});
        registry.add("hDecLengthD", "B^{0}_{s} candidates;#it{p}_{T}(D_{s}) (GeV/#it{c});D_{s} candidate decay length (cm);entries", {HistType::kTH2F, {axisPtDs, axisDecayLength}});
        registry.add("hDecLengthXyD", "B^{0}_{s} candidates;#it{p}_{T}(D_{s}) (GeV/#it{c});D_{s} candidate decay length XY (cm);entries", {HistType::kTH2F, {axisPtDs, axisDecayLength}});
        registry.add("hCospD", "B^{0}_{s} candidates;#it{p}_{T}(D_{s}) (GeV/#it{c});D_{s} candidate cos(#vartheta_{P});entries", {HistType::kTH2F, {axisPtDs, axisCosp}});
        registry.add("hCospXyD", "B^{0}_{s} candidates;#it{p}_{T}(D_{s}) (GeV/#it{c});D_{s} candidate cos(#vartheta_{P}^{XY});entries", {HistType::kTH2F, {axisPtDs, axisCosp}});

        // ML scores of Ds- daughter
        if (doprocessDataWithDmesMl) {
          registry.add("hMlScoreBkgDs", "B^{0}_{s} candidates;#it{p}_{T}(D_{s}) (GeV/#it{c});prong0, Ds ML background score;entries", {HistType::kTH2F, {axisPtDs, axisMlScore}});
          registry.add("hMlScorePromptDs", "B^{0}_{s} candidates;#it{p}_{T}(D_{s}) (GeV/#it{c});prong0, Ds ML prompt score;entries", {HistType::kTH2F, {axisPtDs, axisMlScore}});
          registry.add("hMlScoreNonPromptDs", "B^{0}_{s} candidates;#it{p}_{T}(D_{s}) (GeV/#it{c});prong0, Ds ML nonprompt score;entries", {HistType::kTH2F, {axisPtDs, axisMlScore}});
        }

        // ML scores of Bs candidate
        if (doprocessDataWithBsMl) {
          registry.add("hMlScoreSigBs", "B^{0}_{s} candidates;#it{p}_{T}(B^{0}_{s}) (GeV/#it{c});prong0, B^{0}_{s} ML signal score;entries", {HistType::kTH2F, {axisPtBs, axisMlScore}});
        }
      }
      if (fillSparses) {
        if (!(doprocessDataWithDmesMl || doprocessDataWithBsMl)) {
          registry.add("hMassPtCutVars", "B^{0}_{s} candidates;#it{M} (D_{s}#pi) (GeV/#it{c}^{2});#it{p}_{T}(B^{0}_{s}) (GeV/#it{c});B^{0}_{s} candidate decay length (cm);B^{0}_{s} candidate norm. decay length XY (cm);B^{0}_{s} candidate impact parameter product (cm);B^{0}_{s} candidate cos(#vartheta_{P});#it{M} (KK#pi) (GeV/#it{c}^{2});#it{p}_{T}(D_{s}) (GeV/#it{c});D_{s} candidate decay length (cm);D_{s} candidate cos(#vartheta_{P})", {HistType::kTHnSparseF, {axisMassBs, axisPtBs, axisDecayLength, axisNormDecayLength, axisImpParProd, axisCosp, axisMassDs, axisPtDs, axisDecayLength, axisCosp}});
        } else {
          registry.add("hMassPtCutVars", "B^{0}_{s} candidates;#it{M} (D_{s}#pi) (GeV/#it{c}^{2});#it{p}_{T}(B^{0}_{s}) (GeV/#it{c});B^{0}_{s} candidate decay length (cm);B^{0}_{s} candidate norm. decay length XY (cm);B^{0}_{s} candidate impact parameter product (cm);B^{0}_{s} candidate cos(#vartheta_{P});#it{M} (KK#pi) (GeV/#it{c}^{2});#it{p}_{T}(D_{s}) (GeV/#it{c});D_{s} candidate ML score bkg;D_{s} candidate ML score nonprompt", {HistType::kTHnSparseF, {axisMassBs, axisPtBs, axisDecayLength, axisNormDecayLength, axisImpParProd, axisCosp, axisMassDs, axisPtDs, axisMlScore, axisMlScore}});
        }
      }
    }

    if (doprocessMc || doprocessMcWithDecayTypeCheck || doprocessMcWithDmesMl || doprocessMcWithDmesMlAndDecayTypeCheck || doprocessMcWithBsMl || doprocessMcWithBsMlAndDecayTypeCheck) {
      if (fillHistograms) {
        // gen histos
        registry.add("hEtaGen", "B^{0}_{s} particles (generated);#it{p}_{T}^{gen}(B^{0}_{s}) (GeV/#it{c});#it{#eta}^{gen}(B^{0}_{s});entries", {HistType::kTH2F, {axisPtBs, axisEta}});
        registry.add("hYGen", "B^{0}_{s} particles (generated);#it{p}_{T}^{gen}(B^{0}_{s}) (GeV/#it{c});#it{y}^{gen}(B^{0}_{s});entries", {HistType::kTH2F, {axisPtBs, axisEta}});
        registry.add("hYGenWithProngsInAcceptance", "MC particles (generated-daughters in acceptance);#it{p}_{T}^{gen}(B^{0}_{s}) (GeV/#it{c});#it{y}^{gen}(B^{0}_{s});entries", {HistType::kTH2F, {axisPtBs, axisEta}});
        registry.add("hPtProng0Gen", "B^{0}_{s} particles (generated);#it{p}_{T}^{gen}(B^{0}_{s}) (GeV/#it{c});#it{p}_{T}^{gen}(D_{s}) (GeV/#it{c});entries", {HistType::kTH2F, {axisPtBs, axisPtDs}});
        registry.add("hPtProng1Gen", "B^{0}_{s} particles (generated);#it{p}_{T}^{gen}(B^{0}_{s}) (GeV/#it{c});#it{p}_{T}^{gen}(#pi) (GeV/#it{c});entries", {HistType::kTH2F, {axisPtBs, axisPtPi}});
        registry.add("hYProng0Gen", "B^{0}_{s} particles (generated);#it{p}_{T}^{gen}(B^{0}_{s}) (GeV/#it{c});#it{y}^{gen}(D_{s});entries", {HistType::kTH2F, {axisPtBs, axisEta}});
        registry.add("hYProng1Gen", "B^{0}_{s} particles (generated);#it{p}_{T}^{gen}(B^{0}_{s}) (GeV/#it{c});#it{y}^{gen}(#pi);entries", {HistType::kTH2F, {axisPtBs, axisEta}});
        registry.add("hEtaProng0Gen", "B^{0}_{s} particles (generated);#it{p}_{T}^{gen}(B^{0}_{s}) (GeV/#it{c});#it{#eta}^{gen}(D_{s});entries", {HistType::kTH2F, {axisPtBs, axisEta}});
        registry.add("hEtaProng1Gen", "B^{0}_{s} particles (generated);#it{p}_{T}^{gen}(B^{0}_{s}) (GeV/#it{c});#it{#eta}^{gen}(#pi);entries", {HistType::kTH2F, {axisPtBs, axisEta}});

        // reco histos
        // signal
        registry.add("hMassRecSig", "B^{0}_{s} candidates (matched);#it{p}_{T}(B^{0}_{s}) (GeV/#it{c});#it{M} (D_{s}#pi) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {axisPtBs, axisMassBs}});
        registry.add("hDecLengthRecSig", "B^{0}_{s} candidates (matched);#it{p}_{T}(B^{0}_{s}) (GeV/#it{c});B^{0}_{s} candidate decay length (cm);entries", {HistType::kTH2F, {axisPtBs, axisDecayLength}});
        registry.add("hDecLengthXyRecSig", "B^{0}_{s} candidates (matched);#it{p}_{T}(B^{0}_{s}) (GeV/#it{c});B^{0}_{s} candidate decay length XY (cm);entries", {HistType::kTH2F, {axisPtBs, axisDecayLength}});
        registry.add("hNormDecLengthXyRecSig", "B^{0}_{s} candidates (matched);#it{p}_{T}(B^{0}_{s}) (GeV/#it{c});B^{0}_{s} candidate norm. decay length XY (cm);entries", {HistType::kTH2F, {axisPtBs, axisNormDecayLength}});
        registry.add("hDcaProng0RecSig", "B^{0}_{s} candidates (matched);#it{p}_{T}(B^{0}_{s}) (GeV/#it{c});prong 0 (D_{s}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {axisPtBs, axisDca}});
        registry.add("hDcaProng1RecSig", "B^{0}_{s} candidates (matched);#it{p}_{T}(B^{0}_{s}) (GeV/#it{c});prong 1 (#pi) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {axisPtBs, axisDca}});
        registry.add("hPtProng0RecSig", "B^{0}_{s} candidates (matched);#it{p}_{T}(B^{0}_{s}) (GeV/#it{c});#it{p}_{T}(D_{s}) (GeV/#it{c});entries", {HistType::kTH2F, {axisPtBs, axisPtDs}});
        registry.add("hPtProng1RecSig", "B^{0}_{s} candidates (matched);#it{p}_{T}(B^{0}_{s}) (GeV/#it{c});#it{p}_{T}(#pi) (GeV/#it{c});entries", {HistType::kTH2F, {axisPtBs, axisPtPi}});
        registry.add("hCospRecSig", "B^{0}_{s} candidates (matched);#it{p}_{T}(B^{0}_{s}) (GeV/#it{c});B^{0}_{s} candidate cos(#vartheta_{P});entries", {HistType::kTH2F, {axisPtBs, axisCosp}});
        registry.add("hCospXyRecSig", "B^{0}_{s} candidates (matched);#it{p}_{T}(B^{0}_{s}) (GeV/#it{c});B^{0}_{s} candidate cos(#vartheta_{P}^{XY});entries", {HistType::kTH2F, {axisPtBs, axisCosp}});
        registry.add("hEtaRecSig", "B^{0}_{s} candidates (matched);#it{p}_{T}(B^{0}_{s}) (GeV/#it{c});B^{0}_{s} candidate #it{#eta};entries", {HistType::kTH2F, {axisPtBs, axisEta}});
        registry.add("hRapidityRecSig", "B^{0}_{s} candidates (matched);#it{p}_{T}(B^{0}_{s}) (GeV/#it{c});B^{0}_{s} candidate #it{y};entries", {HistType::kTH2F, {axisPtBs, axisEta}});
        registry.add("hImpParProdRecSig", "B^{0}_{s} candidates (matched);#it{p}_{T}(B^{0}_{s}) (GeV/#it{c});B^{0}_{s} candidate impact parameter product;entries", {HistType::kTH2F, {axisPtBs, axisImpParProd}});
        registry.add("hInvMassDRecSig", "B^{0}_{s} candidates (matched);#it{p}_{T}(D_{s}) (GeV/#it{c});prong0, #it{M}(KK#pi) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {axisPtDs, axisMassDs}});
        registry.add("hDecLengthDRecSig", "B^{0}_{s} candidates (matched);#it{p}_{T}(D_{s}) (GeV/#it{c});D_{s} candidate decay length (cm);entries", {HistType::kTH2F, {axisPtDs, axisDecayLength}});
        registry.add("hDecLengthXyDRecSig", "B^{0}_{s} candidates (matched);#it{p}_{T}(D_{s}) (GeV/#it{c});D_{s} candidate decay length XY (cm);entries", {HistType::kTH2F, {axisPtDs, axisDecayLength}});
        registry.add("hCospDRecSig", "B^{0}_{s} candidates (matched);#it{p}_{T}(D_{s}) (GeV/#it{c});D_{s} candidate cos(#vartheta_{P});entries", {HistType::kTH2F, {axisPtDs, axisCosp}});
        registry.add("hCospXyDRecSig", "B^{0}_{s} candidates (matched);#it{p}_{T}(D_{s}) (GeV/#it{c});D_{s} candidate cos(#vartheta_{P}^{XY});entries", {HistType::kTH2F, {axisPtDs, axisCosp}});
        // background
        if (fillBackground) {
          registry.add("hMassRecBg", "B^{0}_{s} candidates (unmatched);#it{p}_{T}(B^{0}_{s}) (GeV/#it{c});#it{M} (D_{s}#pi) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {axisPtBs, axisMassBs}});
          registry.add("hDecLengthRecBg", "B^{0}_{s} candidates (unmatched);#it{p}_{T}(B^{0}_{s}) (GeV/#it{c});B^{0}_{s} candidate decay length (cm);entries", {HistType::kTH2F, {axisPtBs, axisDecayLength}});
          registry.add("hDecLengthXyRecBg", "B^{0}_{s} candidates (unmatched);#it{p}_{T}(B^{0}_{s}) (GeV/#it{c});B^{0}_{s} candidate decay length XY (cm);entries", {HistType::kTH2F, {axisPtBs, axisDecayLength}});
          registry.add("hNormDecLengthXyRecBg", "B^{0}_{s} candidates (unmatched);#it{p}_{T}(B^{0}_{s}) (GeV/#it{c});B^{0}_{s} candidate norm. decay length XY (cm);entries", {HistType::kTH2F, {axisPtBs, axisNormDecayLength}});
          registry.add("hDcaProng0RecBg", "B^{0}_{s} candidates (unmatched);#it{p}_{T}(B^{0}_{s}) (GeV/#it{c});prong 0 (D_{s}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {axisPtBs, axisDca}});
          registry.add("hDcaProng1RecBg", "B^{0}_{s} candidates (unmatched);#it{p}_{T}(B^{0}_{s}) (GeV/#it{c});prong 1 (#pi) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {axisPtBs, axisDca}});
          registry.add("hPtProng0RecBg", "B^{0}_{s} candidates (unmatched);#it{p}_{T}(B^{0}_{s}) (GeV/#it{c});#it{p}_{T}(D_{s}) (GeV/#it{c});entries", {HistType::kTH2F, {axisPtBs, axisPtDs}});
          registry.add("hPtProng1RecBg", "B^{0}_{s} candidates (unmatched);#it{p}_{T}(B^{0}_{s}) (GeV/#it{c});#it{p}_{T}(#pi) (GeV/#it{c});entries", {HistType::kTH2F, {axisPtBs, axisPtPi}});
          registry.add("hCospRecBg", "B^{0}_{s} candidates (unmatched);#it{p}_{T}(B^{0}_{s}) (GeV/#it{c});B^{0}_{s} candidate cos(#vartheta_{P});entries", {HistType::kTH2F, {axisPtBs, axisCosp}});
          registry.add("hCospXyRecBg", "B^{0}_{s} candidates (unmatched);#it{p}_{T}(B^{0}_{s}) (GeV/#it{c});B^{0}_{s} candidate cos(#vartheta_{P}^{XY});entries", {HistType::kTH2F, {axisPtBs, axisCosp}});
          registry.add("hEtaRecBg", "B^{0}_{s} candidates (unmatched);#it{p}_{T}(B^{0}_{s}) (GeV/#it{c});B^{0}_{s} candidate #it{#eta};entries", {HistType::kTH2F, {axisPtBs, axisEta}});
          registry.add("hRapidityRecBg", "B^{0}_{s} candidates (unmatched);#it{p}_{T}(B^{0}_{s}) (GeV/#it{c});B^{0}_{s} candidate #it{y};entries", {HistType::kTH2F, {axisPtBs, axisEta}});
          registry.add("hImpParProdRecBg", "B^{0}_{s} candidates (unmatched);#it{p}_{T}(B^{0}_{s}) (GeV/#it{c});B^{0}_{s} candidate impact parameter product;entries", {HistType::kTH2F, {axisPtBs, axisImpParProd}});
          registry.add("hInvMassDRecBg", "B^{0}_{s} candidates (unmatched);#it{p}_{T}(D_{s}) (GeV/#it{c});prong0, #it{M}(KK#pi) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {axisPtDs, axisMassDs}});
          registry.add("hDecLengthDRecBg", "B^{0}_{s} candidates (unmatched);#it{p}_{T}(D_{s}) (GeV/#it{c});D_{s} candidate decay length (cm);entries", {HistType::kTH2F, {axisPtDs, axisDecayLength}});
          registry.add("hDecLengthXyDRecBg", "B^{0}_{s} candidates (unmatched);#it{p}_{T}(D_{s}) (GeV/#it{c});D_{s} candidate decay length XY (cm);entries", {HistType::kTH2F, {axisPtDs, axisDecayLength}});
          registry.add("hCospDRecBg", "B^{0}_{s} candidates (unmatched);#it{p}_{T}(D_{s}) (GeV/#it{c});D_{s} candidate cos(#vartheta_{P});entries", {HistType::kTH2F, {axisPtDs, axisCosp}});
          registry.add("hCospXyDRecBg", "B^{0}_{s} candidates (unmatched);#it{p}_{T}(D_{s}) (GeV/#it{c});D_{s} candidate cos(#vartheta_{P}^{XY});entries", {HistType::kTH2F, {axisPtDs, axisCosp}});
        }
        // MC checks
        if (doprocessMcWithDecayTypeCheck || doprocessMcWithBsMlAndDecayTypeCheck || doprocessMcWithDmesMlAndDecayTypeCheck) {
          constexpr uint8_t kNBinsDecayTypeMc = hf_cand_bs::DecayTypeMc::NDecayTypeMc;
          TString labels[kNBinsDecayTypeMc];
          labels[hf_cand_bs::DecayTypeMc::BsToDsPiToPhiPiPiToKKPiPi] = "B^{0}_{s} #rightarrow (D_{s} #rightarrow #Phi#pi #rightarrow KK#pi) #pi";
          labels[hf_cand_bs::DecayTypeMc::BsToDsPiToK0starKPiToKKPiPi] = "B^{0}_{s} #rightarrow (D_{s} #rightarrow K^{0*}K #rightarrow KK#pi) #pi";
          labels[hf_cand_bs::DecayTypeMc::B0ToDsPiToPhiPiPiToKKPiPi] = "B^{0} #rightarrow (D_{s} #rightarrow #Phi#pi #rightarrow KK#pi) #pi";
          labels[hf_cand_bs::DecayTypeMc::B0ToDsPiToK0starKPiToKKPiPi] = "B^{0} #rightarrow (D_{s} #rightarrow K^{0*}K #rightarrow KK#pi) #pi";
          labels[hf_cand_bs::DecayTypeMc::PartlyRecoDecay] = "Partly reconstructed decay channel";
          labels[hf_cand_bs::DecayTypeMc::OtherDecay] = "Other decays";
          static const AxisSpec axisDecayType = {kNBinsDecayTypeMc, 0.5, kNBinsDecayTypeMc + 0.5, ""};
          registry.add("hDecayTypeMc", "DecayType", {HistType::kTH3F, {axisDecayType, axisMassBs, axisPtBs}});
          for (uint8_t iBin = 0; iBin < kNBinsDecayTypeMc; ++iBin) {
            registry.get<TH3>(HIST("hDecayTypeMc"))->GetXaxis()->SetBinLabel(iBin + 1, labels[iBin]);
          }
        }
        // ML scores of Ds- daughter
        if (doprocessMcWithDmesMl || doprocessMcWithDmesMlAndDecayTypeCheck) {
          // signal
          registry.add("hMlScoreBkgDsRecSig", "B^{0}_{s} candidates (matched);#it{p}_{T}(D_{s}) (GeV/#it{c});prong0, Ds ML background score;entries", {HistType::kTH2F, {axisPtDs, axisMlScore}});
          registry.add("hMlScorePromptDsRecSig", "B^{0}_{s} candidates (matched);#it{p}_{T}(D_{s}) (GeV/#it{c});prong0, Ds ML prompt score;entries", {HistType::kTH2F, {axisPtDs, axisMlScore}});
          registry.add("hMlScoreNonPromptDsRecSig", "B^{0}_{s} candidates (matched);#it{p}_{T}(D_{s}) (GeV/#it{c});prong0, Ds ML nonprompt score;entries", {HistType::kTH2F, {axisPtDs, axisMlScore}});
          // background
          registry.add("hMlScoreBkgDsRecBg", "B^{0}_{s} candidates (unmatched);#it{p}_{T}(D_{s}) (GeV/#it{c});prong0, Ds ML background score;entries", {HistType::kTH2F, {axisPtDs, axisMlScore}});
          registry.add("hMlScorePromptDsRecBg", "B^{0}_{s} candidates (unmatched);#it{p}_{T}(D_{s}) (GeV/#it{c});prong0, Ds ML prompt score;entries", {HistType::kTH2F, {axisPtDs, axisMlScore}});
          registry.add("hMlScoreNonPromptDsRecBg", "B^{0}_{s} candidates (unmatched);#it{p}_{T}(D_{s}) (GeV/#it{c});prong0, Ds ML nonprompt score;entries", {HistType::kTH2F, {axisPtDs, axisMlScore}});
        }
        // ML scores of Bs candidate
        if (doprocessMcWithBsMl || doprocessMcWithBsMlAndDecayTypeCheck) {
          // signal
          registry.add("hMlScoreSigBsRecSig", "B^{0}_{s} candidates (matched);#it{p}_{T}(B^{0}_{s}) (GeV/#it{c});prong0, B^{0}_{s} ML signal score;entries", {HistType::kTH2F, {axisPtBs, axisMlScore}});
          // background
          registry.add("hMlScoreSigBsRecBg", "B^{0}_{s} candidates (unmatched);#it{p}_{T}(B^{0}_{s}) (GeV/#it{c});prong0, B^{0}_{s} ML signal score;entries", {HistType::kTH2F, {axisPtBs, axisMlScore}});
        }
      }
      if (fillSparses) {
        // gen sparses
        registry.add("hPtYGenSig", "B^{0}_{s} particles (generated);#it{p}_{T}(B^{0}_{s}) (GeV/#it{c});#it{y}(B^{0}_{s})", {HistType::kTHnSparseF, {axisPtBs, axisEta}});
        registry.add("hPtYWithProngsInAccepanceGenSig", "B^{0}_{s} particles (generated-daughters in acceptance);#it{p}_{T}(B^{0}_{s}) (GeV/#it{c});#it{y}(B^{0}_{s})", {HistType::kTHnSparseF, {axisPtBs, axisEta}});

        // reco sparses
        if (!(doprocessDataWithDmesMl || doprocessDataWithBsMl)) {
          registry.add("hMassPtCutVarsRecSig", "B^{0}_{s} candidates (matched);#it{M} (D_{s}#pi) (GeV/#it{c}^{2});#it{p}_{T}(B^{0}_{s}) (GeV/#it{c});B^{0}_{s} candidate decay length (cm);B^{0}_{s} candidate norm. decay length XY (cm);B^{0}_{s} candidate impact parameter product (cm);B^{0}_{s} candidate cos(#vartheta_{P});#it{M} (KK#pi) (GeV/#it{c}^{2});#it{p}_{T}(D_{s}) (GeV/#it{c});D_{s} candidate decay length (cm);D_{s} candidate cos(#vartheta_{P})", {HistType::kTHnSparseF, {axisMassBs, axisPtBs, axisDecayLength, axisNormDecayLength, axisImpParProd, axisCosp, axisMassDs, axisPtDs, axisDecayLength, axisCosp}});
          if (fillBackground) {
            registry.add("hMassPtCutVarsRecBg", "B^{0}_{s} candidates (unmatched);#it{M} (D_{s}#pi) (GeV/#it{c}^{2});#it{p}_{T}(B^{0}_{s}) (GeV/#it{c});B^{0}_{s} candidate decay length (cm);B^{0}_{s} candidate norm. decay length XY (cm);B^{0}_{s} candidate impact parameter product (cm);B^{0}_{s} candidate cos(#vartheta_{P});#it{M} (KK#pi) (GeV/#it{c}^{2});#it{p}_{T}(D_{s}) (GeV/#it{c});D_{s} candidate decay length (cm);D_{s} candidate cos(#vartheta_{P})", {HistType::kTHnSparseF, {axisMassBs, axisPtBs, axisDecayLength, axisNormDecayLength, axisImpParProd, axisCosp, axisMassDs, axisPtDs, axisDecayLength, axisCosp}});
          }
        } else {
          registry.add("hMassPtCutVarsRecSig", "B^{0}_{s} candidates (matched);#it{M} (D_{s}#pi) (GeV/#it{c}^{2});#it{p}_{T}(B^{0}_{s}) (GeV/#it{c});B^{0}_{s} candidate decay length (cm);B^{0}_{s} candidate norm. decay length XY (cm);B^{0}_{s} candidate impact parameter product (cm);B^{0}_{s} candidate cos(#vartheta_{P});#it{M} (KK#pi) (GeV/#it{c}^{2});#it{p}_{T}(D_{s}) (GeV/#it{c});D_{s} candidate ML score bkg;D_{s} candidate ML score nonprompt", {HistType::kTHnSparseF, {axisMassBs, axisPtBs, axisDecayLength, axisNormDecayLength, axisImpParProd, axisCosp, axisMassDs, axisPtDs, axisMlScore, axisMlScore}});
          if (fillBackground) {
            registry.add("hMassPtCutVarsRecBg", "B^{0}_{s} candidates (unmatched);#it{M} (D_{s}#pi) (GeV/#it{c}^{2});#it{p}_{T}(B^{0}_{s}) (GeV/#it{c});B^{0}_{s} candidate decay length (cm);B^{0}_{s} candidate norm. decay length XY (cm);B^{0}_{s} candidate impact parameter product (cm);B^{0}_{s} candidate cos(#vartheta_{P});#it{M} (KK#pi) (GeV/#it{c}^{2});#it{p}_{T}(D_{s}) (GeV/#it{c});D_{s} candidate ML score bkg;D_{s} candidate ML score nonprompt", {HistType::kTHnSparseF, {axisMassBs, axisPtBs, axisDecayLength, axisNormDecayLength, axisImpParProd, axisCosp, axisMassDs, axisPtDs, axisMlScore, axisMlScore}});
          }
        }
      }
    }
  }

  /// Selection of Bs daughter in geometrical acceptance
  /// \param etaProng is the pseudorapidity of Bs prong
  /// \param ptProng is the pT of Bs prong
  /// \return true if prong is in geometrical acceptance
  template <typename T = float>
  bool isProngInAcceptance(const T& etaProng, const T& ptProng)
  {
    return std::abs(etaProng) <= etaTrackMax && ptProng >= ptTrackMin;
  }

  /// Fill candidate information at reconstruction level
  /// \param doMc is the flag to enable the filling with MC information
  /// \param withDecayTypeCheck is the flag to enable MC with decay type check
  /// \param withDmesMl is the flag to enable the filling with ML scores for the Ds- daughter
  /// \param withBsMl is the flag to enable the filling with ML scores for the Bs candidate
  /// \param candidate is the Bs candidate
  /// \param candidatesD is the table with Ds- candidates
  template <bool doMc, bool withDecayTypeCheck, bool withDmesMl, bool withBsMl, typename Cand>
  void fillCand(Cand const& candidate,
                aod::HfRed3Prongs const&)
  {
    auto ptCandBs = candidate.pt();
    auto invMassBs = hfHelper.invMassBsToDsPi(candidate);
    auto candDs = candidate.template prong0_as<aod::HfRed3Prongs>();
    auto ptDs = candidate.ptProng0();
    auto invMassDs = candDs.invMassHypo0() > 0 ? candDs.invMassHypo0() : candDs.invMassHypo1();
    // TODO: here we are assuming that only one of the two hypotheses is filled, to be checked
    std::array<float, 3> posPv{candidate.posX(), candidate.posY(), candidate.posZ()};
    std::array<float, 3> posSvDs{candDs.xSecondaryVertex(), candDs.ySecondaryVertex(), candDs.zSecondaryVertex()};
    std::array<float, 3> momDs{candDs.pVector()};
    auto cospDs = RecoDecay::cpa(posPv, posSvDs, momDs);
    auto cospXyDs = RecoDecay::cpaXY(posPv, posSvDs, momDs);
    auto decLenDs = RecoDecay::distance(posPv, posSvDs);
    auto decLenXyDs = RecoDecay::distanceXY(posPv, posSvDs);

    int8_t flagMcMatchRec = 0;
    int8_t flagWrongCollision = 0;
    bool isSignal = false;
    if constexpr (doMc) {
      flagMcMatchRec = candidate.flagMcMatchRec();
      flagWrongCollision = candidate.flagWrongCollision();
      isSignal = TESTBIT(std::abs(flagMcMatchRec), hf_cand_bs::DecayTypeMc::BsToDsPiToPhiPiPiToKKPiPi);
    }

    if (fillHistograms) {
      if constexpr (doMc) {
        if (isSignal) {
          registry.fill(HIST("hMassRecSig"), ptCandBs, invMassBs);
          registry.fill(HIST("hPtProng0RecSig"), ptCandBs, candidate.ptProng0());
          registry.fill(HIST("hPtProng1RecSig"), ptCandBs, candidate.ptProng1());
          registry.fill(HIST("hImpParProdRecSig"), ptCandBs, candidate.impactParameterProduct());
          registry.fill(HIST("hDecLengthRecSig"), ptCandBs, candidate.decayLength());
          registry.fill(HIST("hDecLengthXyRecSig"), ptCandBs, candidate.decayLengthXY());
          registry.fill(HIST("hNormDecLengthXyRecSig"), ptCandBs, candidate.decayLengthXY() / candidate.errorDecayLengthXY());
          registry.fill(HIST("hDcaProng0RecSig"), ptCandBs, candidate.impactParameter0());
          registry.fill(HIST("hDcaProng1RecSig"), ptCandBs, candidate.impactParameter1());
          registry.fill(HIST("hCospRecSig"), ptCandBs, candidate.cpa());
          registry.fill(HIST("hCospXyRecSig"), ptCandBs, candidate.cpaXY());
          registry.fill(HIST("hEtaRecSig"), ptCandBs, candidate.eta());
          registry.fill(HIST("hRapidityRecSig"), ptCandBs, hfHelper.yBs(candidate));
          registry.fill(HIST("hInvMassDRecSig"), ptDs, invMassDs);
          registry.fill(HIST("hDecLengthDRecSig"), ptDs, decLenDs);
          registry.fill(HIST("hDecLengthXyDRecSig"), ptDs, decLenXyDs);
          registry.fill(HIST("hCospDRecSig"), ptDs, cospDs);
          registry.fill(HIST("hCospXyDRecSig"), ptDs, cospXyDs);
          if constexpr (withDecayTypeCheck) {
            registry.fill(HIST("hDecayTypeMc"), 1 + hf_cand_bs::DecayTypeMc::BsToDsPiToPhiPiPiToKKPiPi, invMassBs, ptCandBs);
          }
          if constexpr (withDmesMl) {
            registry.fill(HIST("hMlScoreBkgDsRecSig"), ptDs, candidate.prong0MlScoreBkg());
            registry.fill(HIST("hMlScorePromptDsRecSig"), ptDs, candidate.prong0MlScorePrompt());
            registry.fill(HIST("hMlScoreNonPromptDsRecSig"), ptDs, candidate.prong0MlScoreNonprompt());
          }
          if constexpr (withBsMl) {
            registry.fill(HIST("hMlScoreSigBsRecSig"), ptCandBs, candidate.mlProbBsToDsPi()[1]);
          }
        } else if (fillBackground) {
          registry.fill(HIST("hMassRecBg"), ptCandBs, invMassBs);
          registry.fill(HIST("hPtProng0RecBg"), ptCandBs, candidate.ptProng0());
          registry.fill(HIST("hPtProng1RecBg"), ptCandBs, candidate.ptProng1());
          registry.fill(HIST("hImpParProdRecBg"), ptCandBs, candidate.impactParameterProduct());
          registry.fill(HIST("hDecLengthRecBg"), ptCandBs, candidate.decayLength());
          registry.fill(HIST("hDecLengthXyRecBg"), ptCandBs, candidate.decayLengthXY());
          registry.fill(HIST("hNormDecLengthXyRecBg"), ptCandBs, candidate.decayLengthXY() / candidate.errorDecayLengthXY());
          registry.fill(HIST("hDcaProng0RecBg"), ptCandBs, candidate.impactParameter0());
          registry.fill(HIST("hDcaProng1RecBg"), ptCandBs, candidate.impactParameter1());
          registry.fill(HIST("hCospRecBg"), ptCandBs, candidate.cpa());
          registry.fill(HIST("hCospXyRecBg"), ptCandBs, candidate.cpaXY());
          registry.fill(HIST("hEtaRecBg"), ptCandBs, candidate.eta());
          registry.fill(HIST("hRapidityRecBg"), ptCandBs, hfHelper.yBs(candidate));
          registry.fill(HIST("hInvMassDRecBg"), ptDs, invMassDs);
          registry.fill(HIST("hDecLengthDRecBg"), ptDs, decLenDs);
          registry.fill(HIST("hDecLengthXyDRecBg"), ptDs, decLenXyDs);
          registry.fill(HIST("hCospDRecBg"), ptDs, cospDs);
          registry.fill(HIST("hCospXyDRecBg"), ptDs, cospXyDs);
          if constexpr (withDmesMl) {
            registry.fill(HIST("hMlScoreBkgDsRecBg"), ptDs, candidate.prong0MlScoreBkg());
            registry.fill(HIST("hMlScorePromptDsRecBg"), ptDs, candidate.prong0MlScorePrompt());
            registry.fill(HIST("hMlScoreNonPromptDsRecBg"), ptDs, candidate.prong0MlScoreNonprompt());
          }
          if constexpr (withBsMl) {
            registry.fill(HIST("hMlScoreSigBsRecBg"), ptCandBs, candidate.mlProbBsToDsPi()[1]);
          }
        } else if constexpr (withDecayTypeCheck) {
          for (uint8_t iFlag = 1; iFlag < hf_cand_bs::DecayTypeMc::NDecayTypeMc; ++iFlag) {
            if (TESTBIT(flagMcMatchRec, iFlag)) {
              registry.fill(HIST("hDecayTypeMc"), 1 + iFlag, invMassBs, ptCandBs);
            }
          }
        }
      } else {
        registry.fill(HIST("hMass"), ptCandBs, invMassBs);
        registry.fill(HIST("hPtProng0"), ptCandBs, candidate.ptProng0());
        registry.fill(HIST("hPtProng1"), ptCandBs, candidate.ptProng1());
        registry.fill(HIST("hImpParProd"), ptCandBs, candidate.impactParameterProduct());
        registry.fill(HIST("hDecLength"), ptCandBs, candidate.decayLength());
        registry.fill(HIST("hDecLengthXy"), ptCandBs, candidate.decayLengthXY());
        registry.fill(HIST("hNormDecLengthXy"), ptCandBs, candidate.decayLengthXY() / candidate.errorDecayLengthXY());
        registry.fill(HIST("hDcaProng0"), ptCandBs, candidate.impactParameter0());
        registry.fill(HIST("hDcaProng1"), ptCandBs, candidate.impactParameter1());
        registry.fill(HIST("hCosp"), ptCandBs, candidate.cpa());
        registry.fill(HIST("hCospXy"), ptCandBs, candidate.cpaXY());
        registry.fill(HIST("hEta"), ptCandBs, candidate.eta());
        registry.fill(HIST("hRapidity"), ptCandBs, hfHelper.yBs(candidate));
        registry.fill(HIST("hInvMassD"), ptDs, invMassDs);
        registry.fill(HIST("hDecLengthD"), ptDs, decLenDs);
        registry.fill(HIST("hDecLengthXyD"), ptDs, decLenXyDs);
        registry.fill(HIST("hCospD"), ptDs, cospDs);
        registry.fill(HIST("hCospXyD"), ptDs, cospXyDs);

        if constexpr (withDmesMl) {
          registry.fill(HIST("hMlScoreBkgDs"), ptDs, candidate.prong0MlScoreBkg());
          registry.fill(HIST("hMlScorePromptDs"), ptDs, candidate.prong0MlScorePrompt());
          registry.fill(HIST("hMlScoreNonPromptDs"), ptDs, candidate.prong0MlScoreNonprompt());
        }
        if constexpr (withBsMl) {
          registry.fill(HIST("hMlScoreSigBs"), ptCandBs, candidate.mlProbBsToDsPi()[1]);
        }
      }
    }
    if (fillSparses) {
      if constexpr (withDmesMl) {
        if (isSignal) {
          if constexpr (withDmesMl) {
            registry.fill(HIST("hMassPtCutVarsRecSig"), invMassBs, ptCandBs, candidate.decayLength(), candidate.decayLengthXY() / candidate.errorDecayLengthXY(), candidate.impactParameterProduct(), candidate.cpa(), invMassDs, ptDs, candidate.prong0MlScoreBkg(), candidate.prong0MlScoreNonprompt());
          } else {
            registry.fill(HIST("hMassPtCutVarsRecSig"), invMassBs, ptCandBs, candidate.decayLength(), candidate.decayLengthXY() / candidate.errorDecayLengthXY(), candidate.impactParameterProduct(), candidate.cpa(), invMassDs, ptDs, decLenDs, cospDs);
          }
        } else if (fillBackground) {
          if constexpr (withDmesMl) {
            registry.fill(HIST("hMassPtCutVarsRecBg"), invMassBs, ptCandBs, candidate.decayLength(), candidate.decayLengthXY() / candidate.errorDecayLengthXY(), candidate.impactParameterProduct(), candidate.cpa(), invMassDs, ptDs, candidate.prong0MlScoreBkg(), candidate.prong0MlScoreNonprompt());
          } else {
            registry.fill(HIST("hMassPtCutVarsRecBg"), invMassBs, ptCandBs, candidate.decayLength(), candidate.decayLengthXY() / candidate.errorDecayLengthXY(), candidate.impactParameterProduct(), candidate.cpa(), invMassDs, ptDs, decLenDs, cospDs);
          }
        }
      } else {
        if constexpr (withDmesMl) {
          registry.fill(HIST("hMassPtCutVars"), invMassBs, ptCandBs, candidate.decayLength(), candidate.decayLengthXY() / candidate.errorDecayLengthXY(), candidate.impactParameterProduct(), candidate.cpa(), invMassDs, ptDs, candidate.prong0MlScoreBkg(), candidate.prong0MlScoreNonprompt());
        } else {
          registry.fill(HIST("hMassPtCutVars"), invMassBs, ptCandBs, candidate.decayLength(), candidate.decayLengthXY() / candidate.errorDecayLengthXY(), candidate.impactParameterProduct(), candidate.cpa(), invMassDs, ptDs, decLenDs, cospDs);
        }
      }
    }
    if (fillTree) {
      float pseudoRndm = ptDs * 1000. - static_cast<int64_t>(ptDs * 1000);
      if (flagMcMatchRec != 0 || (((doMc && fillBackground) || !doMc) && (ptCandBs >= ptMaxForDownSample || pseudoRndm < downSampleBkgFactor))) {
        float prong0MlScoreBkg = -1.;
        float prong0MlScorePrompt = -1.;
        float prong0MlScoreNonprompt = -1.;
        float candidateMlScoreSig = -1;
        if constexpr (withDmesMl) {
          prong0MlScoreBkg = candidate.prong0MlScoreBkg();
          prong0MlScorePrompt = candidate.prong0MlScorePrompt();
          prong0MlScoreNonprompt = candidate.prong0MlScoreNonprompt();
        }
        if constexpr (withBsMl) {
          candidateMlScoreSig = candidate.mlProbBsToDsPi()[1];
        }
        auto prong1 = candidate.template prong1_as<TracksPion>();

        float ptMother = -1.;
        if constexpr (doMc) {
          ptMother = candidate.ptMother();
        }

        hfRedCandBsLite(
          candidate.chi2PCA(),
          candidate.decayLength(),
          candidate.decayLengthXY(),
          candidate.decayLengthNormalised(),
          candidate.decayLengthXYNormalised(),
          invMassDs,
          ptDs,
          candidate.ptProng1(),
          candidate.impactParameter0(),
          candidate.impactParameter1(),
          candidate.impactParameterProduct(),
          prong1.tpcNSigmaPi(),
          prong1.tofNSigmaPi(),
          prong0MlScoreBkg,
          prong0MlScorePrompt,
          prong0MlScoreNonprompt,
          candidateMlScoreSig,
          candidate.isSelBsToDsPi(),
          invMassBs,
          ptCandBs,
          candidate.cpa(),
          candidate.cpaXY(),
          candidate.maxNormalisedDeltaIP(),
          candidate.eta(),
          candidate.phi(),
          hfHelper.yBs(candidate),
          flagMcMatchRec,
          isSignal,
          flagWrongCollision,
          ptMother);

        if constexpr (withDecayTypeCheck) {
          float candidateMlScoreSig = -1;
          if constexpr (withBsMl) {
            candidateMlScoreSig = candidate.mlProbBsToDsPi()[1];
          }
          hfRedBsMcCheck(
            flagMcMatchRec,
            flagWrongCollision,
            invMassDs,
            ptDs,
            invMassBs,
            ptCandBs,
            candidateMlScoreSig,
            candidate.pdgCodeBeautyMother(),
            candidate.pdgCodeCharmMother(),
            candidate.pdgCodeProng0(),
            candidate.pdgCodeProng1(),
            candidate.pdgCodeProng2(),
            candidate.pdgCodeProng3());
        }
      }
    }
  }

  /// Fill particle histograms (gen MC truth)
  void fillCandMcGen(aod::HfMcGenRedBss::iterator const& particle)
  {
    // keep only generated Bs with the analysis decay channel
    if (!TESTBIT(std::abs(particle.flagMcMatchGen()), hf_cand_bs::DecayTypeMc::BsToDsPiToPhiPiPiToKKPiPi)) {
      return;
    }
    auto ptParticle = particle.ptTrack();
    auto yParticle = particle.yTrack();
    auto etaParticle = particle.etaTrack();
    if (yCandGenMax >= 0. && std::abs(yParticle) > yCandGenMax) {
      return;
    }
    std::array<float, 2> ptProngs = {particle.ptProng0(), particle.ptProng1()};
    std::array<float, 2> yProngs = {particle.yProng0(), particle.yProng1()};
    std::array<float, 2> etaProngs = {particle.etaProng0(), particle.etaProng1()};
    bool prongsInAcc = isProngInAcceptance(etaProngs[0], ptProngs[0]) && isProngInAcceptance(etaProngs[1], ptProngs[1]);

    if (fillHistograms) {
      registry.fill(HIST("hPtProng0Gen"), ptParticle, ptProngs[0]);
      registry.fill(HIST("hPtProng1Gen"), ptParticle, ptProngs[1]);
      registry.fill(HIST("hYProng0Gen"), ptParticle, yProngs[0]);
      registry.fill(HIST("hYProng1Gen"), ptParticle, yProngs[1]);
      registry.fill(HIST("hEtaProng0Gen"), ptParticle, etaProngs[0]);
      registry.fill(HIST("hEtaProng1Gen"), ptParticle, etaProngs[1]);

      registry.fill(HIST("hYGen"), ptParticle, yParticle);
      registry.fill(HIST("hEtaGen"), ptParticle, etaParticle);

      // generated Bs with daughters in geometrical acceptance
      if (prongsInAcc) {
        registry.fill(HIST("hYGenWithProngsInAcceptance"), ptParticle, yParticle);
      }
    }
    if (fillSparses) {
      registry.fill(HIST("hPtYGenSig"), ptParticle, yParticle);
      if (prongsInAcc) {
        registry.fill(HIST("hPtYWithProngsInAccepanceGenSig"), ptParticle, yParticle);
      }
    }
  }

  // Process functions
  void processData(soa::Filtered<soa::Join<aod::HfRedCandBs, aod::HfSelBsToDsPi>> const& candidates,
                   aod::HfRed3Prongs const& candidatesD,
                   TracksPion const&)
  {
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(hfHelper.yBs(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCand<false, false, false, false>(candidate, candidatesD);
    } // candidate loop
  } // processData
  PROCESS_SWITCH(HfTaskBsReduced, processData, "Process data without ML scores for Bs and D daughter", true);

  void processDataWithDmesMl(soa::Filtered<soa::Join<aod::HfRedCandBs, aod::HfRedBsDsMls, aod::HfSelBsToDsPi>> const& candidates,
                             aod::HfRed3Prongs const& candidatesD,
                             TracksPion const&)
  {
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(hfHelper.yBs(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCand<false, false, true, false>(candidate, candidatesD);
    } // candidate loop
  } // processDataWithDmesMl
  PROCESS_SWITCH(HfTaskBsReduced, processDataWithDmesMl, "Process data with(out) ML scores for D daughter (Bs)", false);

  void processDataWithBsMl(soa::Filtered<soa::Join<aod::HfRedCandBs, aod::HfMlBsToDsPi, aod::HfSelBsToDsPi>> const& candidates,
                           aod::HfRed3Prongs const& candidatesD,
                           TracksPion const&)
  {
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(hfHelper.yBs(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCand<false, false, false, true>(candidate, candidatesD);
    } // candidate loop
  } // processDataWithBsMl
  PROCESS_SWITCH(HfTaskBsReduced, processDataWithBsMl, "Process data with(out) ML scores for Bs (D daughter)", false);

  void processMc(soa::Filtered<soa::Join<aod::HfRedCandBs, aod::HfSelBsToDsPi, aod::HfMcRecRedBss>> const& candidates,
                 aod::HfMcGenRedBss const& mcParticles,
                 aod::HfRed3Prongs const& candidatesD,
                 TracksPion const&)
  {
    // MC rec
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(hfHelper.yBs(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCand<true, false, false, false>(candidate, candidatesD);
    } // rec

    // MC gen. level
    for (const auto& particle : mcParticles) {
      fillCandMcGen(particle);
    } // gen
  } // processMc
  PROCESS_SWITCH(HfTaskBsReduced, processMc, "Process MC without ML scores for Bs and D daughter", false);

  void processMcWithDecayTypeCheck(soa::Filtered<soa::Join<aod::HfRedCandBs, aod::HfSelBsToDsPi, aod::HfMcRecRedBss, aod::HfMcCheckBss>> const& candidates,
                                   aod::HfMcGenRedBss const& mcParticles,
                                   aod::HfRed3Prongs const& candidatesD,
                                   TracksPion const&)
  {
    // MC rec
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(hfHelper.yBs(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCand<true, true, false, false>(candidate, candidatesD);
    } // rec

    // MC gen. level
    for (const auto& particle : mcParticles) {
      fillCandMcGen(particle);
    } // gen
  } // processMc
  PROCESS_SWITCH(HfTaskBsReduced, processMcWithDecayTypeCheck, "Process MC with decay type check and without ML scores for Bs and D daughter", false);

  void processMcWithDmesMl(soa::Filtered<soa::Join<aod::HfRedCandBs, aod::HfRedBsDsMls, aod::HfSelBsToDsPi, aod::HfMcRecRedBss>> const& candidates,
                           aod::HfMcGenRedBss const& mcParticles,
                           aod::HfRed3Prongs const& candidatesD,
                           TracksPion const&)
  {
    // MC rec
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(hfHelper.yBs(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCand<true, false, true, false>(candidate, candidatesD);
    } // rec

    // MC gen. level
    for (const auto& particle : mcParticles) {
      fillCandMcGen(particle);
    } // gen
  } // processMcWithDmesMl
  PROCESS_SWITCH(HfTaskBsReduced, processMcWithDmesMl, "Process MC with(out) ML scores for D daughter (Bs)", false);

  void processMcWithDmesMlAndDecayTypeCheck(soa::Filtered<soa::Join<aod::HfRedCandBs, aod::HfRedBsDsMls, aod::HfSelBsToDsPi, aod::HfMcRecRedBss, aod::HfMcCheckBss>> const& candidates,
                                            aod::HfMcGenRedBss const& mcParticles,
                                            aod::HfRed3Prongs const& candidatesD,
                                            TracksPion const&)
  {
    // MC rec
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(hfHelper.yBs(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCand<true, true, true, false>(candidate, candidatesD);
    } // rec

    // MC gen. level
    for (const auto& particle : mcParticles) {
      fillCandMcGen(particle);
    } // gen
  } // processMc
  PROCESS_SWITCH(HfTaskBsReduced, processMcWithDmesMlAndDecayTypeCheck, "Process MC with decay type check and with(out) ML scores for Bs (D daughter)", false);

  void processMcWithBsMl(soa::Filtered<soa::Join<aod::HfRedCandBs, aod::HfMlBsToDsPi, aod::HfSelBsToDsPi, aod::HfMcRecRedBss>> const& candidates,
                         aod::HfMcGenRedBss const& mcParticles,
                         aod::HfRed3Prongs const& candidatesD,
                         TracksPion const&)
  {
    // MC rec
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(hfHelper.yBs(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCand<true, false, false, true>(candidate, candidatesD);
    } // rec

    // MC gen. level
    for (const auto& particle : mcParticles) {
      fillCandMcGen(particle);
    } // gen
  } // processMcWithBsMl
  PROCESS_SWITCH(HfTaskBsReduced, processMcWithBsMl, "Process MC with(out) ML scores for Bs (D daughter)", false);

  void processMcWithBsMlAndDecayTypeCheck(soa::Filtered<soa::Join<aod::HfRedCandBs, aod::HfMlBsToDsPi, aod::HfSelBsToDsPi, aod::HfMcRecRedBss, aod::HfMcCheckBss>> const& candidates,
                                          aod::HfMcGenRedBss const& mcParticles,
                                          aod::HfRed3Prongs const& candidatesD,
                                          TracksPion const&)
  {
    // MC rec
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(hfHelper.yBs(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCand<true, true, false, true>(candidate, candidatesD);
    } // rec

    // MC gen. level
    for (const auto& particle : mcParticles) {
      fillCandMcGen(particle);
    } // gen
  } // processMc
  PROCESS_SWITCH(HfTaskBsReduced, processMcWithBsMlAndDecayTypeCheck, "Process MC with decay type check and with(out) ML scores for B0 (D daughter)", false);
}; // struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskBsReduced>(cfgc)};
}
