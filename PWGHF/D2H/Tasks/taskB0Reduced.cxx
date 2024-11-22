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

/// \file taskB0Reduced.cxx
/// \brief B0 → D- π+ → (π- K+ π-) π+ analysis task
///
/// \author Alexandre Bigot <alexandre.bigot@cern.ch>, IPHC Strasbourg
/// \author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN

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
namespace hf_cand_b0_lite
{
DECLARE_SOA_COLUMN(PtD, ptD, float);                                                     //! Transverse momentum of D-meson daughter candidate (GeV/c)
DECLARE_SOA_COLUMN(PtBach, ptBach, float);                                               //! Transverse momentum of bachelor pion (GeV/c)
DECLARE_SOA_COLUMN(AbsEtaBach, absEtaBach, float);                                       //! Absolute pseudorapidity of bachelor pion
DECLARE_SOA_COLUMN(ItsNClsBach, itsNClsBach, int);                                       //! Number of ITS clusters of bachelor pion
DECLARE_SOA_COLUMN(TpcNClsCrossedRowsBach, tpcNClsCrossedRowsBach, int);                 //! Number of TPC crossed rows of prongs of bachelor pion
DECLARE_SOA_COLUMN(TpcChi2NClBach, tpcChi2NClBach, float);                               //! Maximum TPC chi2 of prongs of D-meson daughter candidate
DECLARE_SOA_COLUMN(PtDmesProngMin, ptProngDmesMin, float);                               //! Minimum pT of prongs of D-meson daughter candidate (GeV/c)
DECLARE_SOA_COLUMN(AbsEtaDmesProngMin, absEtaProngDmesMin, float);                       //! Minimum absolute pseudorapidity of prongs of D-meson daughter candidate
DECLARE_SOA_COLUMN(ItsNClsDmesProngMin, itsNClsDmesProngMin, int);                       //! Minimum number of ITS clusters of prongs of D-meson daughter candidate
DECLARE_SOA_COLUMN(TpcNClsCrossedRowsDmesProngMin, tpcNClsCrossedRowsDmesProngMin, int); //! Minimum number of TPC crossed rows of prongs of D-meson daughter candidate
DECLARE_SOA_COLUMN(TpcChi2NClDmesProngMax, tpcChi2NClDmesProngMax, float);               //! Maximum TPC chi2 of prongs of D-meson daughter candidate
DECLARE_SOA_COLUMN(MD, mD, float);                                                       //! Invariant mass of D-meson daughter candidates (GeV/c)
DECLARE_SOA_COLUMN(M, m, float);                                                         //! Invariant mass of candidate (GeV/c2)
DECLARE_SOA_COLUMN(Pt, pt, float);                                                       //! Transverse momentum of candidate (GeV/c)
DECLARE_SOA_COLUMN(PtGen, ptGen, float);                                                 //! Transverse momentum of candidate (GeV/c)
DECLARE_SOA_COLUMN(P, p, float);                                                         //! Momentum of candidate (GeV/c)
DECLARE_SOA_COLUMN(Y, y, float);                                                         //! Rapidity of candidate
DECLARE_SOA_COLUMN(Eta, eta, float);                                                     //! Pseudorapidity of candidate
DECLARE_SOA_COLUMN(Phi, phi, float);                                                     //! Azimuth angle of candidate
DECLARE_SOA_COLUMN(E, e, float);                                                         //! Energy of candidate (GeV)
DECLARE_SOA_COLUMN(NSigTpcPiBachelor, nSigTpcPiBachelor, float);                         //! TPC Nsigma separation for bachelor with pion mass hypothesis
DECLARE_SOA_COLUMN(NSigTofPiBachelor, nSigTofPiBachelor, float);                         //! TOF Nsigma separation for bachelor with pion mass hypothesis
DECLARE_SOA_COLUMN(NSigTpcTofPiBachelor, nSigTpcTofPiBachelor, float);                   //! Combined TPC and TOF Nsigma separation for bachelor with pion mass hypothesis
DECLARE_SOA_COLUMN(NSigTpcPiDmesProng0, nSigTpcPiDmesProng0, float);                     //! TPC Nsigma separation for D-meson prong0 with pion mass hypothesis
DECLARE_SOA_COLUMN(NSigTofPiDmesProng0, nSigTofPiDmesProng0, float);                     //! TOF Nsigma separation for D-meson prong0 with pion mass hypothesis
DECLARE_SOA_COLUMN(NSigTpcTofPiDmesProng0, nSigTpcTofPiDmesProng0, float);               //! Combined TPC and TOF Nsigma separation for D-meson prong0 with pion mass hypothesis
DECLARE_SOA_COLUMN(NSigTpcKaDmesProng1, nSigTpcKaDmesProng1, float);                     //! TPC Nsigma separation for D-meson prong1 with kaon mass hypothesis
DECLARE_SOA_COLUMN(NSigTofKaDmesProng1, nSigTofKaDmesProng1, float);                     //! TOF Nsigma separation for D-meson prong1 with kaon mass hypothesis
DECLARE_SOA_COLUMN(NSigTpcTofKaDmesProng1, nSigTpcTofKaDmesProng1, float);               //! Combined TPC and TOF Nsigma separation for D-meson prong1 with kaon mass hypothesis
DECLARE_SOA_COLUMN(NSigTpcPiDmesProng2, nSigTpcPiDmesProng2, float);                     //! TPC Nsigma separation for D-meson prong2 with pion mass hypothesis
DECLARE_SOA_COLUMN(NSigTofPiDmesProng2, nSigTofPiDmesProng2, float);                     //! TOF Nsigma separation for D-meson prong2 with pion mass hypothesis
DECLARE_SOA_COLUMN(NSigTpcTofPiDmesProng2, nSigTpcTofPiDmesProng2, float);               //! Combined TPC and TOF Nsigma separation for D-meson prong0 with pion mass hypothesis
DECLARE_SOA_COLUMN(DecayLength, decayLength, float);                                     //! Decay length of candidate (cm)
DECLARE_SOA_COLUMN(DecayLengthXY, decayLengthXY, float);                                 //! Transverse decay length of candidate (cm)
DECLARE_SOA_COLUMN(DecayLengthNormalised, decayLengthNormalised, float);                 //! Normalised decay length of candidate
DECLARE_SOA_COLUMN(DecayLengthXYNormalised, decayLengthXYNormalised, float);             //! Normalised transverse decay length of candidate
DECLARE_SOA_COLUMN(DecayLengthD, decayLengthD, float);                                   //! Decay length of D-meson daughter candidate (cm)
DECLARE_SOA_COLUMN(DecayLengthXYD, decayLengthXYD, float);                               //! Transverse decay length of D-meson daughter candidate (cm)
DECLARE_SOA_COLUMN(ImpactParameterD, impactParameterD, float);                           //! Impact parameter product of D-meson daughter candidate
DECLARE_SOA_COLUMN(ImpactParameterBach, impactParameterBach, float);                     //! Impact parameter product of bachelor pion
DECLARE_SOA_COLUMN(ImpactParameterProduct, impactParameterProduct, float);               //! Impact parameter product of daughters
DECLARE_SOA_COLUMN(Cpa, cpa, float);                                                     //! Cosine pointing angle of candidate
DECLARE_SOA_COLUMN(CpaXY, cpaXY, float);                                                 //! Cosine pointing angle of candidate in transverse plane
DECLARE_SOA_COLUMN(MaxNormalisedDeltaIP, maxNormalisedDeltaIP, float);                   //! Maximum normalized difference between measured and expected impact parameter of candidate prongs
DECLARE_SOA_COLUMN(MlScoreSig, mlScoreSig, float);                                       //! ML score for signal class
DECLARE_SOA_COLUMN(FlagWrongCollision, flagWrongCollision, int8_t);                      //! Flag for association with wrong collision
} // namespace hf_cand_b0_lite

DECLARE_SOA_TABLE(HfRedCandB0Lites, "AOD", "HFREDCANDB0LITE", //! Table with some B0 properties
                                                              // B meson features                  hf_cand_b0_lite::M,
                  hf_cand_b0_lite::M,
                  hf_cand_b0_lite::Pt,
                  hf_cand_b0_lite::Eta,
                  hf_cand_b0_lite::Phi,
                  hf_cand_b0_lite::Y,
                  hf_cand_b0_lite::Cpa,
                  hf_cand_b0_lite::CpaXY,
                  hf_cand::Chi2PCA,
                  hf_cand_b0_lite::DecayLength,
                  hf_cand_b0_lite::DecayLengthXY,
                  hf_cand_b0_lite::DecayLengthNormalised,
                  hf_cand_b0_lite::DecayLengthXYNormalised,
                  hf_cand_b0_lite::ImpactParameterProduct,
                  hf_cand_b0_lite::MaxNormalisedDeltaIP,
                  hf_cand_b0_lite::MlScoreSig,
                  hf_sel_candidate_b0::IsSelB0ToDPi,
                  // D meson features
                  hf_cand_b0_lite::MD,
                  hf_cand_b0_lite::PtD,
                  hf_cand_b0_lite::DecayLengthD,
                  hf_cand_b0_lite::DecayLengthXYD,
                  hf_cand_b0_lite::ImpactParameterD,
                  hf_cand_b0_lite::PtDmesProngMin,
                  hf_cand_b0_lite::AbsEtaDmesProngMin,
                  hf_cand_b0_lite::ItsNClsDmesProngMin,
                  hf_cand_b0_lite::TpcNClsCrossedRowsDmesProngMin,
                  hf_cand_b0_lite::TpcChi2NClDmesProngMax,
                  hf_cand_b0_lite::NSigTpcPiDmesProng0,
                  hf_cand_b0_lite::NSigTofPiDmesProng0,
                  hf_cand_b0_lite::NSigTpcTofPiDmesProng0,
                  hf_cand_b0_lite::NSigTpcKaDmesProng1,
                  hf_cand_b0_lite::NSigTofKaDmesProng1,
                  hf_cand_b0_lite::NSigTpcTofKaDmesProng1,
                  hf_cand_b0_lite::NSigTpcPiDmesProng2,
                  hf_cand_b0_lite::NSigTofPiDmesProng2,
                  hf_cand_b0_lite::NSigTpcTofPiDmesProng2,
                  hf_cand_b0_reduced::Prong0MlScoreBkg,
                  hf_cand_b0_reduced::Prong0MlScorePrompt,
                  hf_cand_b0_reduced::Prong0MlScoreNonprompt,
                  // pion features
                  hf_cand_b0_lite::PtBach,
                  hf_cand_b0_lite::AbsEtaBach,
                  hf_cand_b0_lite::ItsNClsBach,
                  hf_cand_b0_lite::TpcNClsCrossedRowsBach,
                  hf_cand_b0_lite::TpcChi2NClBach,
                  hf_cand_b0_lite::ImpactParameterBach,
                  hf_cand_b0_lite::NSigTpcPiBachelor,
                  hf_cand_b0_lite::NSigTofPiBachelor,
                  hf_cand_b0_lite::NSigTpcTofPiBachelor,
                  // MC truth
                  hf_cand_3prong::FlagMcMatchRec,
                  hf_cand_3prong::OriginMcRec,
                  hf_cand_b0_lite::FlagWrongCollision,
                  hf_cand_b0_lite::PtGen);

DECLARE_SOA_TABLE(HfRedB0McCheck, "AOD", "HFREDB0MCCHECK", //! Table with MC decay type check
                  hf_cand_3prong::FlagMcMatchRec,
                  hf_cand_b0_lite::FlagWrongCollision,
                  hf_cand_b0_lite::MD,
                  hf_cand_b0_lite::PtD,
                  hf_cand_b0_lite::M,
                  hf_cand_b0_lite::Pt,
                  hf_cand_b0_lite::MlScoreSig,
                  hf_b0_mc::PdgCodeBeautyMother,
                  hf_b0_mc::PdgCodeCharmMother,
                  hf_b0_mc::PdgCodeProng0,
                  hf_b0_mc::PdgCodeProng1,
                  hf_b0_mc::PdgCodeProng2,
                  hf_b0_mc::PdgCodeProng3);
} // namespace o2::aod

/// B0 analysis task
struct HfTaskB0Reduced {
  Produces<aod::HfRedCandB0Lites> hfRedCandB0Lite;
  Produces<aod::HfRedB0McCheck> hfRedB0McCheck;

  Configurable<int> selectionFlagB0{"selectionFlagB0", 1, "Selection Flag for B0"};
  Configurable<float> yCandGenMax{"yCandGenMax", 0.5, "max. gen particle rapidity"};
  Configurable<float> yCandRecoMax{"yCandRecoMax", 0.8, "max. cand. rapidity"};
  Configurable<float> etaTrackMax{"etaTrackMax", 0.8, "max. track pseudo-rapidity for acceptance calculation"};
  Configurable<float> ptTrackMin{"ptTrackMin", 0.1, "min. track transverse momentum for acceptance calculation"};
  Configurable<bool> fillHistograms{"fillHistograms", true, "Flag to enable histogram filling"};
  Configurable<bool> fillSparses{"fillSparses", false, "Flag to enable sparse filling"};
  Configurable<bool> fillTree{"fillTree", false, "Flag to enable tree filling"};
  Configurable<bool> fillBackground{"fillBackground", false, "Flag to enable filling of background histograms/sparses/tree (only MC)"};
  Configurable<float> downSampleBkgFactor{"downSampleBkgFactor", 1., "Fraction of background candidates to keep for ML trainings"};
  Configurable<float> ptMaxForDownSample{"ptMaxForDownSample", 10., "Maximum pt for the application of the downsampling factor"};

  HfHelper hfHelper;

  Filter filterSelectCandidates = (aod::hf_sel_candidate_b0::isSelB0ToDPi >= selectionFlagB0);

  HistogramRegistry registry{"registry"};

  using TracksPion = soa::Join<HfRedTracks, HfRedTracksPid>;
  using CandsDplus = soa::Join<HfRed3Prongs, HfRedPidDau0s, HfRedPidDau1s, HfRedPidDau2s>;

  void init(InitContext&)
  {
    std::array<bool, 3> processFuncData{doprocessData, doprocessDataWithDmesMl, doprocessDataWithB0Ml};
    if ((std::accumulate(processFuncData.begin(), processFuncData.end(), 0)) > 1) {
      LOGP(fatal, "Only one process function for data can be enabled at a time.");
    }
    std::array<bool, 6> processFuncMc{doprocessMc, doprocessMcWithDecayTypeCheck, doprocessMcWithDmesMl, doprocessMcWithDmesMlAndDecayTypeCheck, doprocessMcWithB0Ml, doprocessMcWithB0MlAndDecayTypeCheck};
    if ((std::accumulate(processFuncMc.begin(), processFuncMc.end(), 0)) > 1) {
      LOGP(fatal, "Only one process function for MC can be enabled at a time.");
    }

    const AxisSpec axisMlScore{100, 0.f, 1.f};
    const AxisSpec axisMassB0{300, 4.5f, 6.0f};
    const AxisSpec axisMassDminus{300, 1.75f, 2.05f};
    const AxisSpec axisDecayLength{200, 0.f, 0.4f};
    const AxisSpec axisNormDecayLength{100, 0.f, 50.f};
    const AxisSpec axisDca{100, -0.05f, 0.05f};
    const AxisSpec axisCosp{110, 0.f, 1.1f};
    const AxisSpec axisEta{30, -1.5f, 1.5f};
    const AxisSpec axisError{100, 0.f, 1.f};
    const AxisSpec axisImpParProd{100, -1.e-3, 1.e-3};
    const AxisSpec axisPtB0{100, 0.f, 50.f};
    const AxisSpec axisPtDminus{100, 0.f, 50.f};
    const AxisSpec axisPtPi{100, 0.f, 10.f};

    if (doprocessData || doprocessDataWithDmesMl || doprocessDataWithB0Ml) {
      if (fillHistograms) {
        registry.add("hMass", "B^{0} candidates;#it{p}_{T}(B^{0}) (GeV/#it{c});#it{M} (D#pi) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {axisPtB0, axisMassB0}});
        registry.add("hDecLength", "B^{0} candidates;#it{p}_{T}(B^{0}) (GeV/#it{c});B^{0} candidate decay length (cm);entries", {HistType::kTH2F, {axisPtB0, axisDecayLength}});
        registry.add("hDecLengthXy", "B^{0} candidates;#it{p}_{T}(B^{0}) (GeV/#it{c});decay length XY (cm);entries", {HistType::kTH2F, {axisPtB0, axisDecayLength}});
        registry.add("hNormDecLengthXy", "B^{0} candidates;#it{p}_{T}(B^{0}) (GeV/#it{c});B^{0} candidate norm. decay length XY (cm);entries", {HistType::kTH2F, {axisPtB0, axisNormDecayLength}});
        registry.add("hDcaProng0", "B^{0} candidates;#it{p}_{T}(B^{0}) (GeV/#it{c});prong 0 (D^{#minus}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {axisPtB0, axisDca}});
        registry.add("hDcaProng1", "B^{0} candidates;#it{p}_{T}(B^{0}) (GeV/#it{c});prong 1 (#pi^{#plus}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {axisPtB0, axisDca}});
        registry.add("hPtProng0", "B^{0} candidates;#it{p}_{T}(B^{0}) (GeV/#it{c});#it{p}_{T}(D^{#minus}) (GeV/#it{c});entries", {HistType::kTH2F, {axisPtB0, axisPtDminus}});
        registry.add("hPtProng1", "B^{0} candidates;#it{p}_{T}(B^{0}) (GeV/#it{c});#it{p}_{T}(#pi^{#plus}) (GeV/#it{c});entries", {HistType::kTH2F, {axisPtB0, axisPtPi}});
        registry.add("hCosp", "B^{0} candidates;#it{p}_{T}(B^{0}) (GeV/#it{c});B^{0} candidate cos(#vartheta_{P});entries", {HistType::kTH2F, {axisPtB0, axisCosp}});
        registry.add("hCospXy", "B^{0} candidates;#it{p}_{T}(B^{0}) (GeV/#it{c});B^{0} candidate cos(#vartheta_{P}^{XY});entries", {HistType::kTH2F, {axisPtB0, axisCosp}});
        registry.add("hEta", "B^{0} candidates;#it{p}_{T}(B^{0}) (GeV/#it{c});B^{0} candidate #it{#eta};entries", {HistType::kTH2F, {axisPtB0, axisEta}});
        registry.add("hRapidity", "B^{0} candidates;#it{p}_{T}(B^{0}) (GeV/#it{c});B^{0} candidate #it{y};entries", {HistType::kTH2F, {axisPtB0, axisEta}});
        registry.add("hImpParProd", "B^{0} candidates;#it{p}_{T}(B^{0}) (GeV/#it{c});B^{0} candidate impact parameter product;entries", {HistType::kTH2F, {axisPtB0, axisImpParProd}});
        registry.add("hInvMassD", "B^{0} candidates;#it{p}_{T}(D^{#minus}) (GeV/#it{c});prong0, #it{M}(K#pi) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {axisPtDminus, axisMassDminus}});
        registry.add("hDecLengthD", "B^{0} candidates;#it{p}_{T}(D^{#minus}) (GeV/#it{c});D^{#minus} candidate decay length (cm);entries", {HistType::kTH2F, {axisPtDminus, axisDecayLength}});
        registry.add("hDecLengthXyD", "B^{0} candidates;#it{p}_{T}(D^{#minus}) (GeV/#it{c});decay length XY (cm);entries", {HistType::kTH2F, {axisPtDminus, axisDecayLength}});
        registry.add("hCospD", "B^{0} candidates;#it{p}_{T}(D^{#minus}) (GeV/#it{c});D^{#minus} candidate cos(#vartheta_{P});entries", {HistType::kTH2F, {axisPtDminus, axisCosp}});
        registry.add("hCospXyD", "B^{0} candidates;#it{p}_{T}(D^{#minus}) (GeV/#it{c});D^{#minus} candidate cos(#vartheta_{P}^{XY});entries", {HistType::kTH2F, {axisPtDminus, axisCosp}});

        // ML scores of D- daughter
        if (doprocessDataWithDmesMl) {
          registry.add("hMlScoreBkgD", "B^{0} candidates;#it{p}_{T}(D^{#minus}) (GeV/#it{c});prong0, D^{#minus} ML background score;entries", {HistType::kTH2F, {axisPtDminus, axisMlScore}});
          registry.add("hMlScorePromptD", "B^{0} candidates;#it{p}_{T}(D^{#minus}) (GeV/#it{c});prong0, D^{#minus} ML prompt score;entries", {HistType::kTH2F, {axisPtDminus, axisMlScore}});
          registry.add("hMlScoreNonPromptD", "B^{0} candidates;#it{p}_{T}(D^{#minus}) (GeV/#it{c});prong0, D^{#minus} ML nonprompt score;entries", {HistType::kTH2F, {axisPtDminus, axisMlScore}});
        }

        // ML scores of B0 candidate
        if (doprocessDataWithB0Ml) {
          registry.add("hMlScoreSigB0", "B^{0} candidates;#it{p}_{T}(B^{0}) (GeV/#it{c});prong0, B^{0} ML signal score;entries", {HistType::kTH2F, {axisPtB0, axisMlScore}});
        }
      }
      if (fillSparses) {
        if (!(doprocessDataWithDmesMl || doprocessDataWithB0Ml)) {
          registry.add("hMassPtCutVars", "B^{0} candidates;#it{M} (D#pi) (GeV/#it{c}^{2});#it{p}_{T}(B^{0}) (GeV/#it{c});B^{0} candidate decay length (cm);B^{0} candidate norm. decay length XY (cm);B^{0} candidate impact parameter product (cm);B^{0} candidate cos(#vartheta_{P});#it{M} (K#pi) (GeV/#it{c}^{2});#it{p}_{T}(D^{#minus}) (GeV/#it{c});D^{#minus} candidate decay length (cm);D^{#minus} candidate cos(#vartheta_{P})", {HistType::kTHnSparseF, {axisMassB0, axisPtB0, axisDecayLength, axisNormDecayLength, axisImpParProd, axisCosp, axisMassDminus, axisPtDminus, axisDecayLength, axisCosp}});
        } else {
          registry.add("hMassPtCutVars", "B^{0} candidates;#it{M} (D#pi) (GeV/#it{c}^{2});#it{p}_{T}(B^{0}) (GeV/#it{c});B^{0} candidate decay length (cm);B^{0} candidate norm. decay length XY (cm);B^{0} candidate impact parameter product (cm);B^{0} candidate cos(#vartheta_{P});#it{M} (K#pi) (GeV/#it{c}^{2});#it{p}_{T}(D^{#minus}) (GeV/#it{c});D^{#minus} candidate ML score bkg;D^{#minus} candidate ML score nonprompt", {HistType::kTHnSparseF, {axisMassB0, axisPtB0, axisDecayLength, axisNormDecayLength, axisImpParProd, axisCosp, axisMassDminus, axisPtDminus, axisMlScore, axisMlScore}});
        }
      }
    }

    if (doprocessMc || doprocessMcWithDecayTypeCheck || doprocessMcWithDmesMl || doprocessMcWithDmesMlAndDecayTypeCheck || doprocessMcWithB0Ml || doprocessMcWithB0MlAndDecayTypeCheck) {
      if (fillHistograms) {
        // gen histos
        registry.add("hEtaGen", "B^{0} particles (generated);#it{p}_{T}^{gen}(B^{0}) (GeV/#it{c});#it{#eta}^{gen}(B^{0});entries", {HistType::kTH2F, {axisPtB0, axisEta}});
        registry.add("hYGen", "B^{0} particles (generated);#it{p}_{T}^{gen}(B^{0}) (GeV/#it{c});#it{y}^{gen}(B^{0});entries", {HistType::kTH2F, {axisPtB0, axisEta}});
        registry.add("hYGenWithProngsInAcceptance", "MC particles (generated-daughters in acceptance);#it{p}_{T}^{gen}(B^{0}) (GeV/#it{c});#it{y}^{gen}(B^{0});entries", {HistType::kTH2F, {axisPtB0, axisEta}});
        registry.add("hPtProng0Gen", "B^{0} particles (generated);#it{p}_{T}^{gen}(B^{0}) (GeV/#it{c});#it{p}_{T}^{gen}(D^{#minus}) (GeV/#it{c});entries", {HistType::kTH2F, {axisPtB0, axisPtDminus}});
        registry.add("hPtProng1Gen", "B^{0} particles (generated);#it{p}_{T}^{gen}(B^{0}) (GeV/#it{c});#it{p}_{T}^{gen}(#pi^{#plus}) (GeV/#it{c});entries", {HistType::kTH2F, {axisPtB0, axisPtPi}});
        registry.add("hYProng0Gen", "B^{0} particles (generated);#it{p}_{T}^{gen}(B^{0}) (GeV/#it{c});#it{y}^{gen}(D^{#minus});entries", {HistType::kTH2F, {axisPtB0, axisEta}});
        registry.add("hYProng1Gen", "B^{0} particles (generated);#it{p}_{T}^{gen}(B^{0}) (GeV/#it{c});#it{y}^{gen}(#pi^{#plus});entries", {HistType::kTH2F, {axisPtB0, axisEta}});
        registry.add("hEtaProng0Gen", "B^{0} particles (generated);#it{p}_{T}^{gen}(B^{0}) (GeV/#it{c});#it{#eta}^{gen}(D^{#minus});entries", {HistType::kTH2F, {axisPtB0, axisEta}});
        registry.add("hEtaProng1Gen", "B^{0} particles (generated);#it{p}_{T}^{gen}(B^{0}) (GeV/#it{c});#it{#eta}^{gen}(#pi^{#plus});entries", {HistType::kTH2F, {axisPtB0, axisEta}});

        // reco histos
        // signal
        registry.add("hMassRecSig", "B^{0} candidates (matched);#it{p}_{T}(B^{0}) (GeV/#it{c});#it{M} (D#pi) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {axisPtB0, axisMassB0}});
        registry.add("hDecLengthRecSig", "B^{0} candidates (matched);#it{p}_{T}(B^{0}) (GeV/#it{c});B^{0} candidate decay length (cm);entries", {HistType::kTH2F, {axisPtB0, axisDecayLength}});
        registry.add("hDecLengthXyRecSig", "B^{0} candidates (matched);#it{p}_{T}(B^{0}) (GeV/#it{c});decay length XY (cm);entries", {HistType::kTH2F, {axisPtB0, axisDecayLength}});
        registry.add("hNormDecLengthXyRecSig", "B^{0} candidates (matched);#it{p}_{T}(B^{0}) (GeV/#it{c});B^{0} candidate norm. decay length XY (cm);entries", {HistType::kTH2F, {axisPtB0, axisNormDecayLength}});
        registry.add("hDcaProng0RecSig", "B^{0} candidates (matched);#it{p}_{T}(B^{0}) (GeV/#it{c});prong 0 (D^{#minus}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {axisPtB0, axisDca}});
        registry.add("hDcaProng1RecSig", "B^{0} candidates (matched);#it{p}_{T}(B^{0}) (GeV/#it{c});prong 1 (#pi^{#plus}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {axisPtB0, axisDca}});
        registry.add("hPtProng0RecSig", "B^{0} candidates (matched);#it{p}_{T}(B^{0}) (GeV/#it{c});#it{p}_{T}(D^{#minus}) (GeV/#it{c});entries", {HistType::kTH2F, {axisPtB0, axisPtDminus}});
        registry.add("hPtProng1RecSig", "B^{0} candidates (matched);#it{p}_{T}(B^{0}) (GeV/#it{c});#it{p}_{T}(#pi^{#plus}) (GeV/#it{c});entries", {HistType::kTH2F, {axisPtB0, axisPtPi}});
        registry.add("hCospRecSig", "B^{0} candidates (matched);#it{p}_{T}(B^{0}) (GeV/#it{c});B^{0} candidate cos(#vartheta_{P});entries", {HistType::kTH2F, {axisPtB0, axisCosp}});
        registry.add("hCospXyRecSig", "B^{0} candidates (matched);#it{p}_{T}(B^{0}) (GeV/#it{c});B^{0} candidate cos(#vartheta_{P}^{XY});entries", {HistType::kTH2F, {axisPtB0, axisCosp}});
        registry.add("hEtaRecSig", "B^{0} candidates (matched);#it{p}_{T}(B^{0}) (GeV/#it{c});B^{0} candidate #it{#eta};entries", {HistType::kTH2F, {axisPtB0, axisEta}});
        registry.add("hRapidityRecSig", "B^{0} candidates (matched);#it{p}_{T}(B^{0}) (GeV/#it{c});B^{0} candidate #it{y};entries", {HistType::kTH2F, {axisPtB0, axisEta}});
        registry.add("hImpParProdRecSig", "B^{0} candidates (matched);#it{p}_{T}(B^{0}) (GeV/#it{c});B^{0} candidate impact parameter product;entries", {HistType::kTH2F, {axisPtB0, axisImpParProd}});
        registry.add("hInvMassDRecSig", "B^{0} candidates (matched);#it{p}_{T}(D^{#minus}) (GeV/#it{c});prong0, #it{M}(K#pi) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {axisPtDminus, axisMassDminus}});
        registry.add("hDecLengthDRecSig", "B^{0} candidates (matched);#it{p}_{T}(D^{#minus}) (GeV/#it{c});D^{#minus} candidate decay length (cm);entries", {HistType::kTH2F, {axisPtDminus, axisDecayLength}});
        registry.add("hDecLengthXyDRecSig", "B^{0} candidates (matched);#it{p}_{T}(D^{#minus}) (GeV/#it{c});decay length XY (cm);entries", {HistType::kTH2F, {axisPtDminus, axisDecayLength}});
        registry.add("hCospDRecSig", "B^{0} candidates (matched);#it{p}_{T}(D^{#minus}) (GeV/#it{c});D^{#minus} candidate cos(#vartheta_{P});entries", {HistType::kTH2F, {axisPtDminus, axisCosp}});
        registry.add("hCospXyDRecSig", "B^{0} candidates (matched);#it{p}_{T}(D^{#minus}) (GeV/#it{c});D^{#minus} candidate cos(#vartheta_{P}^{XY});entries", {HistType::kTH2F, {axisPtDminus, axisCosp}});
        // background
        if (fillBackground) {
          registry.add("hMassRecBg", "B^{0} candidates (unmatched);#it{p}_{T}(B^{0}) (GeV/#it{c});#it{M} (D#pi) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {axisPtB0, axisMassB0}});
          registry.add("hDecLengthRecBg", "B^{0} candidates (unmatched);#it{p}_{T}(B^{0}) (GeV/#it{c});B^{0} candidate decay length (cm);entries", {HistType::kTH2F, {axisPtB0, axisDecayLength}});
          registry.add("hDecLengthXyRecBg", "B^{0} candidates (unmatched);#it{p}_{T}(B^{0}) (GeV/#it{c});decay length XY (cm);entries", {HistType::kTH2F, {axisPtB0, axisDecayLength}});
          registry.add("hNormDecLengthXyRecBg", "B^{0} candidates (unmatched);#it{p}_{T}(B^{0}) (GeV/#it{c});B^{0} candidate norm. decay length XY (cm);entries", {HistType::kTH2F, {axisPtB0, axisNormDecayLength}});
          registry.add("hDcaProng0RecBg", "B^{0} candidates (unmatched);#it{p}_{T}(B^{0}) (GeV/#it{c});prong 0 (D^{#minus}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {axisPtB0, axisDca}});
          registry.add("hDcaProng1RecBg", "B^{0} candidates (unmatched);#it{p}_{T}(B^{0}) (GeV/#it{c});prong 1 (#pi^{#plus}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {axisPtB0, axisDca}});
          registry.add("hPtProng0RecBg", "B^{0} candidates (unmatched);#it{p}_{T}(B^{0}) (GeV/#it{c});#it{p}_{T}(D^{#minus}) (GeV/#it{c});entries", {HistType::kTH2F, {axisPtB0, axisPtDminus}});
          registry.add("hPtProng1RecBg", "B^{0} candidates (unmatched);#it{p}_{T}(B^{0}) (GeV/#it{c});#it{p}_{T}(#pi^{#plus}) (GeV/#it{c});entries", {HistType::kTH2F, {axisPtB0, axisPtPi}});
          registry.add("hCospRecBg", "B^{0} candidates (unmatched);#it{p}_{T}(B^{0}) (GeV/#it{c});B^{0} candidate cos(#vartheta_{P});entries", {HistType::kTH2F, {axisPtB0, axisCosp}});
          registry.add("hCospXyRecBg", "B^{0} candidates (unmatched);#it{p}_{T}(B^{0}) (GeV/#it{c});B^{0} candidate cos(#vartheta_{P}^{XY});entries", {HistType::kTH2F, {axisPtB0, axisCosp}});
          registry.add("hEtaRecBg", "B^{0} candidates (unmatched);#it{p}_{T}(B^{0}) (GeV/#it{c});B^{0} candidate #it{#eta};entries", {HistType::kTH2F, {axisPtB0, axisEta}});
          registry.add("hRapidityRecBg", "B^{0} candidates (unmatched);#it{p}_{T}(B^{0}) (GeV/#it{c});B^{0} candidate #it{y};entries", {HistType::kTH2F, {axisPtB0, axisEta}});
          registry.add("hImpParProdRecBg", "B^{0} candidates (unmatched);#it{p}_{T}(B^{0}) (GeV/#it{c});B^{0} candidate impact parameter product;entries", {HistType::kTH2F, {axisPtB0, axisImpParProd}});
          registry.add("hInvMassDRecBg", "B^{0} candidates (unmatched);#it{p}_{T}(D^{#minus}) (GeV/#it{c});prong0, #it{M}(K#pi) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {axisPtDminus, axisMassDminus}});
          registry.add("hDecLengthDRecBg", "B^{0} candidates (unmatched);#it{p}_{T}(D^{#minus}) (GeV/#it{c});D^{#minus} candidate decay length (cm);entries", {HistType::kTH2F, {axisPtDminus, axisDecayLength}});
          registry.add("hDecLengthXyDRecBg", "B^{0} candidates (unmatched);#it{p}_{T}(D^{#minus}) (GeV/#it{c});decay length XY (cm);entries", {HistType::kTH2F, {axisPtDminus, axisDecayLength}});
          registry.add("hCospDRecBg", "B^{0} candidates (unmatched);#it{p}_{T}(D^{#minus}) (GeV/#it{c});D^{#minus} candidate cos(#vartheta_{P});entries", {HistType::kTH2F, {axisPtDminus, axisCosp}});
          registry.add("hCospXyDRecBg", "B^{0} candidates (unmatched);#it{p}_{T}(D^{#minus}) (GeV/#it{c});D^{#minus} candidate cos(#vartheta_{P}^{XY});entries", {HistType::kTH2F, {axisPtDminus, axisCosp}});
        }
        // MC checks
        if (doprocessMcWithDecayTypeCheck || doprocessMcWithB0MlAndDecayTypeCheck || doprocessMcWithDmesMlAndDecayTypeCheck) {
          constexpr uint8_t kNBinsDecayTypeMc = hf_cand_b0::DecayTypeMc::NDecayTypeMc;
          TString labels[kNBinsDecayTypeMc];
          labels[hf_cand_b0::DecayTypeMc::B0ToDplusPiToPiKPiPi] = "B^{0} #rightarrow (D^{#minus} #rightarrow #pi^{#minus} K^{#plus} #pi^{#minus}) #pi^{#plus}";
          labels[hf_cand_b0::DecayTypeMc::B0ToDsPiToKKPiPi] = "B^{0} #rightarrow (D^{#minus}_{s} #rightarrow K^{#minus} K^{#plus} #pi^{#minus}) #pi^{#plus}";
          labels[hf_cand_b0::DecayTypeMc::PartlyRecoDecay] = "Partly reconstructed decay channel";
          labels[hf_cand_b0::DecayTypeMc::OtherDecay] = "Other decays";
          static const AxisSpec axisDecayType = {kNBinsDecayTypeMc, 0.5, kNBinsDecayTypeMc + 0.5, ""};
          registry.add("hDecayTypeMc", "DecayType", {HistType::kTH3F, {axisDecayType, axisMassB0, axisPtB0}});
          for (uint8_t iBin = 0; iBin < kNBinsDecayTypeMc; ++iBin) {
            registry.get<TH3>(HIST("hDecayTypeMc"))->GetXaxis()->SetBinLabel(iBin + 1, labels[iBin]);
          }
        }
        // ML scores of D- daughter
        if (doprocessMcWithDmesMl || doprocessMcWithDmesMlAndDecayTypeCheck) {
          // signal
          registry.add("hMlScoreBkgDRecSig", "B^{0} candidates (matched);#it{p}_{T}(D^{#minus}) (GeV/#it{c});prong0, D^{#minus} ML background score;entries", {HistType::kTH2F, {axisPtDminus, axisMlScore}});
          registry.add("hMlScorePromptDRecSig", "B^{0} candidates (matched);#it{p}_{T}(D^{#minus}) (GeV/#it{c});prong0, D^{#minus} ML prompt score;entries", {HistType::kTH2F, {axisPtDminus, axisMlScore}});
          registry.add("hMlScoreNonPromptDRecSig", "B^{0} candidates (matched);#it{p}_{T}(D^{#minus}) (GeV/#it{c});prong0, D^{#minus} ML nonprompt score;entries", {HistType::kTH2F, {axisPtDminus, axisMlScore}});
          // background
          registry.add("hMlScoreBkgDRecBg", "B^{0} candidates (unmatched);#it{p}_{T}(D^{#minus}) (GeV/#it{c});prong0, D^{#minus} ML background score;entries", {HistType::kTH2F, {axisPtDminus, axisMlScore}});
          registry.add("hMlScorePromptDRecBg", "B^{0} candidates (unmatched);#it{p}_{T}(D^{#minus}) (GeV/#it{c});prong0, D^{#minus} ML prompt score;entries", {HistType::kTH2F, {axisPtDminus, axisMlScore}});
          registry.add("hMlScoreNonPromptDRecBg", "B^{0} candidates (unmatched);#it{p}_{T}(D^{#minus}) (GeV/#it{c});prong0, D^{#minus} ML nonprompt score;entries", {HistType::kTH2F, {axisPtDminus, axisMlScore}});
        }
        // ML scores of B0 candidate
        if (doprocessMcWithB0Ml || doprocessMcWithB0MlAndDecayTypeCheck) {
          // signal
          registry.add("hMlScoreSigB0RecSig", "B^{0} candidates (matched);#it{p}_{T}(B^{0}) (GeV/#it{c});prong0, B^{0} ML signal score;entries", {HistType::kTH2F, {axisPtB0, axisMlScore}});
          // background
          registry.add("hMlScoreSigB0RecBg", "B^{0} candidates (unmatched);#it{p}_{T}(B^{0}) (GeV/#it{c});prong0, B^{0} ML signal score;entries", {HistType::kTH2F, {axisPtB0, axisMlScore}});
        }
      }
      if (fillSparses) {
        // gen sparses
        registry.add("hPtYGenSig", "B^{0} particles (generated);#it{p}_{T}(B^{0}) (GeV/#it{c});#it{y}(B^{0})", {HistType::kTHnSparseF, {axisPtB0, axisEta}});
        registry.add("hPtYWithProngsInAccepanceGenSig", "B^{0} particles (generated-daughters in acceptance);#it{p}_{T}(B^{0}) (GeV/#it{c});#it{y}(B^{0})", {HistType::kTHnSparseF, {axisPtB0, axisEta}});

        // reco sparses
        if (!(doprocessDataWithDmesMl || doprocessDataWithB0Ml)) {
          registry.add("hMassPtCutVarsRecSig", "B^{0} candidates (matched);#it{M} (D#pi) (GeV/#it{c}^{2});#it{p}_{T}(B^{0}) (GeV/#it{c});B^{0} candidate decay length (cm);B^{0} candidate norm. decay length XY (cm);B^{0} candidate impact parameter product (cm);B^{0} candidate cos(#vartheta_{P});#it{M} (K#pi) (GeV/#it{c}^{2});#it{p}_{T}(D^{#minus}) (GeV/#it{c});D^{#minus} candidate decay length (cm);D^{#minus} candidate cos(#vartheta_{P})", {HistType::kTHnSparseF, {axisMassB0, axisPtB0, axisDecayLength, axisNormDecayLength, axisImpParProd, axisCosp, axisMassDminus, axisPtDminus, axisDecayLength, axisCosp}});
          if (fillBackground) {
            registry.add("hMassPtCutVarsRecBg", "B^{0} candidates (unmatched);#it{M} (D#pi) (GeV/#it{c}^{2});#it{p}_{T}(B^{0}) (GeV/#it{c});B^{0} candidate decay length (cm);B^{0} candidate norm. decay length XY (cm);B^{0} candidate impact parameter product (cm);B^{0} candidate cos(#vartheta_{P});#it{M} (K#pi) (GeV/#it{c}^{2});#it{p}_{T}(D^{#minus}) (GeV/#it{c});D^{#minus} candidate decay length (cm);D^{#minus} candidate cos(#vartheta_{P})", {HistType::kTHnSparseF, {axisMassB0, axisPtB0, axisDecayLength, axisNormDecayLength, axisImpParProd, axisCosp, axisMassDminus, axisPtDminus, axisDecayLength, axisCosp}});
          }
        } else {
          registry.add("hMassPtCutVarsRecSig", "B^{0} candidates (matched);#it{M} (D#pi) (GeV/#it{c}^{2});#it{p}_{T}(B^{0}) (GeV/#it{c});B^{0} candidate decay length (cm);B^{0} candidate norm. decay length XY (cm);B^{0} candidate impact parameter product (cm);B^{0} candidate cos(#vartheta_{P});#it{M} (K#pi) (GeV/#it{c}^{2});#it{p}_{T}(D^{#minus}) (GeV/#it{c});D^{#minus} candidate ML score bkg;D^{#minus} candidate ML score nonprompt", {HistType::kTHnSparseF, {axisMassB0, axisPtB0, axisDecayLength, axisNormDecayLength, axisImpParProd, axisCosp, axisMassDminus, axisPtDminus, axisMlScore, axisMlScore}});
          if (fillBackground) {
            registry.add("hMassPtCutVarsRecBg", "B^{0} candidates (unmatched);#it{M} (D#pi) (GeV/#it{c}^{2});#it{p}_{T}(B^{0}) (GeV/#it{c});B^{0} candidate decay length (cm);B^{0} candidate norm. decay length XY (cm);B^{0} candidate impact parameter product (cm);B^{0} candidate cos(#vartheta_{P});#it{M} (K#pi) (GeV/#it{c}^{2});#it{p}_{T}(D^{#minus}) (GeV/#it{c});D^{#minus} candidate ML score bkg;D^{#minus} candidate ML score nonprompt", {HistType::kTHnSparseF, {axisMassB0, axisPtB0, axisDecayLength, axisNormDecayLength, axisImpParProd, axisCosp, axisMassDminus, axisPtDminus, axisMlScore, axisMlScore}});
          }
        }
      }
    }
  }

  /// Selection of B0 daughter in geometrical acceptance
  /// \param etaProng is the pseudorapidity of B0 prong
  /// \param ptProng is the pT of B0 prong
  /// \return true if prong is in geometrical acceptance
  template <typename T = float>
  bool isProngInAcceptance(const T& etaProng, const T& ptProng)
  {
    return std::abs(etaProng) <= etaTrackMax && ptProng >= ptTrackMin;
  }

  /// Fill candidate information at reconstruction level
  /// \param doMc is the flag to enable the filling with MC information
  /// \param withDecayTypeCheck is the flag to enable MC with decay type check
  /// \param withDmesMl is the flag to enable the filling with ML scores for the D- daughter
  /// \param withB0Ml is the flag to enable the filling with ML scores for the B0 candidate
  /// \param candidate is the B0 candidate
  /// \param candidatesD is the table with D- candidates
  template <bool doMc, bool withDecayTypeCheck, bool withDmesMl, bool withB0Ml, typename Cand, typename CandsDmes>
  void fillCand(Cand const& candidate,
                CandsDmes const&)
  {
    auto ptCandB0 = candidate.pt();
    auto invMassB0 = hfHelper.invMassB0ToDPi(candidate);
    auto candD = candidate.template prong0_as<CandsDmes>();
    auto ptD = candidate.ptProng0();
    auto invMassD = candD.invMassHypo0();
    std::array<float, 3> posPv{candidate.posX(), candidate.posY(), candidate.posZ()};
    std::array<float, 3> posSvD{candD.xSecondaryVertex(), candD.ySecondaryVertex(), candD.zSecondaryVertex()};
    std::array<float, 3> momD{candD.pVector()};
    auto cospD = RecoDecay::cpa(posPv, posSvD, momD);
    auto cospXyD = RecoDecay::cpaXY(posPv, posSvD, momD);
    auto decLenD = RecoDecay::distance(posPv, posSvD);
    auto decLenXyD = RecoDecay::distanceXY(posPv, posSvD);

    int8_t flagMcMatchRec = 0;
    int8_t flagWrongCollision = 0;
    bool isSignal = false;
    if constexpr (doMc) {
      flagMcMatchRec = candidate.flagMcMatchRec();
      flagWrongCollision = candidate.flagWrongCollision();
      isSignal = TESTBIT(std::abs(flagMcMatchRec), hf_cand_b0::DecayTypeMc::B0ToDplusPiToPiKPiPi);
    }

    if (fillHistograms) {
      if constexpr (doMc) {
        if (isSignal) {
          registry.fill(HIST("hMassRecSig"), ptCandB0, hfHelper.invMassB0ToDPi(candidate));
          registry.fill(HIST("hPtProng0RecSig"), ptCandB0, candidate.ptProng0());
          registry.fill(HIST("hPtProng1RecSig"), ptCandB0, candidate.ptProng1());
          registry.fill(HIST("hImpParProdRecSig"), ptCandB0, candidate.impactParameterProduct());
          registry.fill(HIST("hDecLengthRecSig"), ptCandB0, candidate.decayLength());
          registry.fill(HIST("hDecLengthXyRecSig"), ptCandB0, candidate.decayLengthXY());
          registry.fill(HIST("hNormDecLengthXyRecSig"), ptCandB0, candidate.decayLengthXY() / candidate.errorDecayLengthXY());
          registry.fill(HIST("hDcaProng0RecSig"), ptCandB0, candidate.impactParameter0());
          registry.fill(HIST("hDcaProng1RecSig"), ptCandB0, candidate.impactParameter1());
          registry.fill(HIST("hCospRecSig"), ptCandB0, candidate.cpa());
          registry.fill(HIST("hCospXyRecSig"), ptCandB0, candidate.cpaXY());
          registry.fill(HIST("hEtaRecSig"), ptCandB0, candidate.eta());
          registry.fill(HIST("hRapidityRecSig"), ptCandB0, hfHelper.yB0(candidate));
          registry.fill(HIST("hInvMassDRecSig"), ptD, invMassD);
          registry.fill(HIST("hDecLengthDRecSig"), ptD, decLenD);
          registry.fill(HIST("hDecLengthXyDRecSig"), ptD, decLenXyD);
          registry.fill(HIST("hCospDRecSig"), ptD, cospD);
          registry.fill(HIST("hCospXyDRecSig"), ptD, cospXyD);
          if constexpr (withDecayTypeCheck) {
            registry.fill(HIST("hDecayTypeMc"), 1 + hf_cand_b0::DecayTypeMc::B0ToDplusPiToPiKPiPi, invMassB0, ptCandB0);
          }
          if constexpr (withDmesMl) {
            registry.fill(HIST("hMlScoreBkgDRecSig"), ptD, candidate.prong0MlScoreBkg());
            registry.fill(HIST("hMlScorePromptDRecSig"), ptD, candidate.prong0MlScorePrompt());
            registry.fill(HIST("hMlScoreNonPromptDRecSig"), ptD, candidate.prong0MlScoreNonprompt());
          }
          if constexpr (withB0Ml) {
            registry.fill(HIST("hMlScoreSigB0RecSig"), ptCandB0, candidate.mlProbB0ToDPi());
          }
        } else if (fillBackground) {
          registry.fill(HIST("hMassRecBg"), ptCandB0, hfHelper.invMassB0ToDPi(candidate));
          registry.fill(HIST("hPtProng0RecBg"), ptCandB0, candidate.ptProng0());
          registry.fill(HIST("hPtProng1RecBg"), ptCandB0, candidate.ptProng1());
          registry.fill(HIST("hImpParProdRecBg"), ptCandB0, candidate.impactParameterProduct());
          registry.fill(HIST("hDecLengthRecBg"), ptCandB0, candidate.decayLength());
          registry.fill(HIST("hDecLengthXyRecBg"), ptCandB0, candidate.decayLengthXY());
          registry.fill(HIST("hNormDecLengthXyRecBg"), ptCandB0, candidate.decayLengthXY() / candidate.errorDecayLengthXY());
          registry.fill(HIST("hDcaProng0RecBg"), ptCandB0, candidate.impactParameter0());
          registry.fill(HIST("hDcaProng1RecBg"), ptCandB0, candidate.impactParameter1());
          registry.fill(HIST("hCospRecBg"), ptCandB0, candidate.cpa());
          registry.fill(HIST("hCospXyRecBg"), ptCandB0, candidate.cpaXY());
          registry.fill(HIST("hEtaRecBg"), ptCandB0, candidate.eta());
          registry.fill(HIST("hRapidityRecBg"), ptCandB0, hfHelper.yB0(candidate));
          registry.fill(HIST("hInvMassDRecBg"), ptD, invMassD);
          registry.fill(HIST("hDecLengthDRecBg"), ptD, decLenD);
          registry.fill(HIST("hDecLengthXyDRecBg"), ptD, decLenXyD);
          registry.fill(HIST("hCospDRecBg"), ptD, cospD);
          registry.fill(HIST("hCospXyDRecBg"), ptD, cospXyD);
          if constexpr (withDmesMl) {
            registry.fill(HIST("hMlScoreBkgDRecBg"), ptD, candidate.prong0MlScoreBkg());
            registry.fill(HIST("hMlScorePromptDRecBg"), ptD, candidate.prong0MlScorePrompt());
            registry.fill(HIST("hMlScoreNonPromptDRecBg"), ptD, candidate.prong0MlScoreNonprompt());
          }
          if constexpr (withB0Ml) {
            registry.fill(HIST("hMlScoreSigB0RecBg"), ptCandB0, candidate.mlProbB0ToDPi());
          }
        } else if constexpr (withDecayTypeCheck) {
          if (TESTBIT(flagMcMatchRec, hf_cand_b0::DecayTypeMc::B0ToDsPiToKKPiPi)) { // B0 → Ds- π+ → (K- K+ π-) π+
            registry.fill(HIST("hDecayTypeMc"), 1 + hf_cand_b0::DecayTypeMc::B0ToDsPiToKKPiPi, invMassB0, ptCandB0);
          } else if (TESTBIT(flagMcMatchRec, hf_cand_b0::DecayTypeMc::PartlyRecoDecay)) { // Partly reconstructed decay channel
            registry.fill(HIST("hDecayTypeMc"), 1 + hf_cand_b0::DecayTypeMc::PartlyRecoDecay, invMassB0, ptCandB0);
          } else {
            registry.fill(HIST("hDecayTypeMc"), 1 + hf_cand_b0::DecayTypeMc::OtherDecay, invMassB0, ptCandB0);
          }
        }
      } else {
        registry.fill(HIST("hMass"), ptCandB0, invMassB0);
        registry.fill(HIST("hPtProng0"), ptCandB0, candidate.ptProng0());
        registry.fill(HIST("hPtProng1"), ptCandB0, candidate.ptProng1());
        registry.fill(HIST("hImpParProd"), ptCandB0, candidate.impactParameterProduct());
        registry.fill(HIST("hDecLength"), ptCandB0, candidate.decayLength());
        registry.fill(HIST("hDecLengthXy"), ptCandB0, candidate.decayLengthXY());
        registry.fill(HIST("hNormDecLengthXy"), ptCandB0, candidate.decayLengthXY() / candidate.errorDecayLengthXY());
        registry.fill(HIST("hDcaProng0"), ptCandB0, candidate.impactParameter0());
        registry.fill(HIST("hDcaProng1"), ptCandB0, candidate.impactParameter1());
        registry.fill(HIST("hCosp"), ptCandB0, candidate.cpa());
        registry.fill(HIST("hCospXy"), ptCandB0, candidate.cpaXY());
        registry.fill(HIST("hEta"), ptCandB0, candidate.eta());
        registry.fill(HIST("hRapidity"), ptCandB0, hfHelper.yB0(candidate));
        registry.fill(HIST("hInvMassD"), ptD, invMassD);
        registry.fill(HIST("hDecLengthD"), ptD, decLenD);
        registry.fill(HIST("hDecLengthXyD"), ptD, decLenXyD);
        registry.fill(HIST("hCospD"), ptD, cospD);
        registry.fill(HIST("hCospXyD"), ptD, cospXyD);

        if constexpr (withDmesMl) {
          registry.fill(HIST("hMlScoreBkgD"), ptD, candidate.prong0MlScoreBkg());
          registry.fill(HIST("hMlScorePromptD"), ptD, candidate.prong0MlScorePrompt());
          registry.fill(HIST("hMlScoreNonPromptD"), ptD, candidate.prong0MlScoreNonprompt());
        }
        if constexpr (withB0Ml) {
          registry.fill(HIST("hMlScoreSigB0"), ptCandB0, candidate.mlProbB0ToDPi());
        }
      }
    }
    if (fillSparses) {
      if constexpr (withDmesMl) {
        if (isSignal) {
          if constexpr (withDmesMl) {
            registry.fill(HIST("hMassPtCutVarsRecSig"), invMassB0, ptCandB0, candidate.decayLength(), candidate.decayLengthXY() / candidate.errorDecayLengthXY(), candidate.impactParameterProduct(), candidate.cpa(), invMassD, ptD, candidate.prong0MlScoreBkg(), candidate.prong0MlScoreNonprompt());
          } else {
            registry.fill(HIST("hMassPtCutVarsRecSig"), invMassB0, ptCandB0, candidate.decayLength(), candidate.decayLengthXY() / candidate.errorDecayLengthXY(), candidate.impactParameterProduct(), candidate.cpa(), invMassD, ptD, decLenD, cospD);
          }
        } else if (fillBackground) {
          if constexpr (withDmesMl) {
            registry.fill(HIST("hMassPtCutVarsRecBg"), invMassB0, ptCandB0, candidate.decayLength(), candidate.decayLengthXY() / candidate.errorDecayLengthXY(), candidate.impactParameterProduct(), candidate.cpa(), invMassD, ptD, candidate.prong0MlScoreBkg(), candidate.prong0MlScoreNonprompt());
          } else {
            registry.fill(HIST("hMassPtCutVarsRecBg"), invMassB0, ptCandB0, candidate.decayLength(), candidate.decayLengthXY() / candidate.errorDecayLengthXY(), candidate.impactParameterProduct(), candidate.cpa(), invMassD, ptD, decLenD, cospD);
          }
        }
      } else {
        if constexpr (withDmesMl) {
          registry.fill(HIST("hMassPtCutVars"), invMassB0, ptCandB0, candidate.decayLength(), candidate.decayLengthXY() / candidate.errorDecayLengthXY(), candidate.impactParameterProduct(), candidate.cpa(), invMassD, ptD, candidate.prong0MlScoreBkg(), candidate.prong0MlScoreNonprompt());
        } else {
          registry.fill(HIST("hMassPtCutVars"), invMassB0, ptCandB0, candidate.decayLength(), candidate.decayLengthXY() / candidate.errorDecayLengthXY(), candidate.impactParameterProduct(), candidate.cpa(), invMassD, ptD, decLenD, cospD);
        }
      }
    }
    if (fillTree) {
      float pseudoRndm = ptD * 1000. - static_cast<int64_t>(ptD * 1000);
      if (flagMcMatchRec != 0 || (((doMc && fillBackground) || !doMc) && (ptCandB0 >= ptMaxForDownSample || pseudoRndm < downSampleBkgFactor))) {
        float prong0MlScoreBkg = -1.;
        float prong0MlScorePrompt = -1.;
        float prong0MlScoreNonprompt = -1.;
        float candidateMlScoreSig = -1;
        if constexpr (withDmesMl) {
          prong0MlScoreBkg = candidate.prong0MlScoreBkg();
          prong0MlScorePrompt = candidate.prong0MlScorePrompt();
          prong0MlScoreNonprompt = candidate.prong0MlScoreNonprompt();
        }
        if constexpr (withB0Ml) {
          candidateMlScoreSig = candidate.mlProbB0ToDPi();
        }
        auto prong1 = candidate.template prong1_as<TracksPion>();

        float ptMother = -1.;
        if constexpr (doMc) {
          ptMother = candidate.ptMother();
        }

        hfRedCandB0Lite(
          // B-meson features
          invMassB0,
          ptCandB0,
          candidate.eta(),
          candidate.phi(),
          hfHelper.yB0(candidate),
          candidate.cpa(),
          candidate.cpaXY(),
          candidate.chi2PCA(),
          candidate.decayLength(),
          candidate.decayLengthXY(),
          candidate.decayLengthNormalised(),
          candidate.decayLengthXYNormalised(),
          candidate.impactParameterProduct(),
          candidate.maxNormalisedDeltaIP(),
          candidateMlScoreSig,
          candidate.isSelB0ToDPi(),
          // D-meson features
          invMassD,
          ptD,
          decLenD,
          decLenXyD,
          candidate.impactParameter0(),
          candD.ptProngMin(),
          candD.absEtaProngMin(),
          candD.itsNClsProngMin(),
          candD.tpcNClsCrossedRowsProngMin(),
          candD.tpcChi2NClProngMax(),
          candD.tpcNSigmaPiProng0(),
          candD.tofNSigmaPiProng0(),
          candD.tpcTofNSigmaPiProng0(),
          candD.tpcNSigmaKaProng1(),
          candD.tofNSigmaKaProng1(),
          candD.tpcTofNSigmaKaProng1(),
          candD.tpcNSigmaPiProng2(),
          candD.tofNSigmaPiProng2(),
          candD.tpcTofNSigmaPiProng2(),
          prong0MlScoreBkg,
          prong0MlScorePrompt,
          prong0MlScoreNonprompt,
          // pion features
          candidate.ptProng1(),
          std::abs(RecoDecay::eta(prong1.pVector())),
          prong1.itsNCls(),
          prong1.tpcNClsCrossedRows(),
          prong1.tpcChi2NCl(),
          candidate.impactParameter1(),
          prong1.tpcNSigmaPi(),
          prong1.tofNSigmaPi(),
          prong1.tpcTofNSigmaPi(),
          // MC truth
          flagMcMatchRec,
          isSignal,
          flagWrongCollision,
          ptMother);

        if constexpr (withDecayTypeCheck) {
          hfRedB0McCheck(
            flagMcMatchRec,
            flagWrongCollision,
            invMassD,
            ptD,
            invMassB0,
            ptCandB0,
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
  void fillCandMcGen(aod::HfMcGenRedB0s::iterator const& particle)
  {
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

      // generated B0 with daughters in geometrical acceptance
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
  void processData(soa::Filtered<soa::Join<aod::HfRedCandB0, aod::HfSelB0ToDPi>> const& candidates,
                   CandsDplus const& candidatesD,
                   TracksPion const&)
  {
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(hfHelper.yB0(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCand<false, false, false, false>(candidate, candidatesD);
    } // candidate loop
  } // processData
  PROCESS_SWITCH(HfTaskB0Reduced, processData, "Process data without ML scores for B0 and D daughter", true);

  void processDataWithDmesMl(soa::Filtered<soa::Join<aod::HfRedCandB0, aod::HfRedB0DpMls, aod::HfSelB0ToDPi>> const& candidates,
                             CandsDplus const& candidatesD,
                             TracksPion const&)
  {
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(hfHelper.yB0(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCand<false, false, true, false>(candidate, candidatesD);
    } // candidate loop
  } // processDataWithDmesMl
  PROCESS_SWITCH(HfTaskB0Reduced, processDataWithDmesMl, "Process data with(out) ML scores for D daughter (B0)", false);

  void processDataWithB0Ml(soa::Filtered<soa::Join<aod::HfRedCandB0, aod::HfMlB0ToDPi, aod::HfSelB0ToDPi>> const& candidates,
                           CandsDplus const& candidatesD,
                           TracksPion const&)
  {
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(hfHelper.yB0(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCand<false, false, false, true>(candidate, candidatesD);
    } // candidate loop
  } // processDataWithB0Ml
  PROCESS_SWITCH(HfTaskB0Reduced, processDataWithB0Ml, "Process data with(out) ML scores for B0 (D daughter)", false);

  void processMc(soa::Filtered<soa::Join<aod::HfRedCandB0, aod::HfSelB0ToDPi, aod::HfMcRecRedB0s>> const& candidates,
                 aod::HfMcGenRedB0s const& mcParticles,
                 CandsDplus const& candidatesD,
                 TracksPion const&)
  {
    // MC rec
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(hfHelper.yB0(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCand<true, false, false, false>(candidate, candidatesD);
    } // rec

    // MC gen. level
    for (const auto& particle : mcParticles) {
      fillCandMcGen(particle);
    } // gen
  } // processMc
  PROCESS_SWITCH(HfTaskB0Reduced, processMc, "Process MC without ML scores for B0 and D daughter", false);

  void processMcWithDecayTypeCheck(soa::Filtered<soa::Join<aod::HfRedCandB0, aod::HfSelB0ToDPi, aod::HfMcRecRedB0s, aod::HfMcCheckB0s>> const& candidates,
                                   aod::HfMcGenRedB0s const& mcParticles,
                                   CandsDplus const& candidatesD,
                                   TracksPion const&)
  {
    // MC rec
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(hfHelper.yB0(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCand<true, true, false, false>(candidate, candidatesD);
    } // rec

    // MC gen. level
    for (const auto& particle : mcParticles) {
      fillCandMcGen(particle);
    } // gen
  } // processMc
  PROCESS_SWITCH(HfTaskB0Reduced, processMcWithDecayTypeCheck, "Process MC with decay type check and without ML scores for B0 and D daughter", false);

  void processMcWithDmesMl(soa::Filtered<soa::Join<aod::HfRedCandB0, aod::HfRedB0DpMls, aod::HfSelB0ToDPi, aod::HfMcRecRedB0s>> const& candidates,
                           aod::HfMcGenRedB0s const& mcParticles,
                           CandsDplus const& candidatesD,
                           TracksPion const&)
  {
    // MC rec
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(hfHelper.yB0(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCand<true, false, true, false>(candidate, candidatesD);
    } // rec

    // MC gen. level
    for (const auto& particle : mcParticles) {
      fillCandMcGen(particle);
    } // gen
  } // processMcWithDmesMl
  PROCESS_SWITCH(HfTaskB0Reduced, processMcWithDmesMl, "Process MC with(out) ML scores for D daughter (B0)", false);

  void processMcWithDmesMlAndDecayTypeCheck(soa::Filtered<soa::Join<aod::HfRedCandB0, aod::HfRedB0DpMls, aod::HfSelB0ToDPi, aod::HfMcRecRedB0s, aod::HfMcCheckB0s>> const& candidates,
                                            aod::HfMcGenRedB0s const& mcParticles,
                                            CandsDplus const& candidatesD,
                                            TracksPion const&)
  {
    // MC rec
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(hfHelper.yB0(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCand<true, true, true, false>(candidate, candidatesD);
    } // rec

    // MC gen. level
    for (const auto& particle : mcParticles) {
      fillCandMcGen(particle);
    } // gen
  } // processMc
  PROCESS_SWITCH(HfTaskB0Reduced, processMcWithDmesMlAndDecayTypeCheck, "Process MC with decay type check and with(out) ML scores for B0 (D daughter)", false);

  void processMcWithB0Ml(soa::Filtered<soa::Join<aod::HfRedCandB0, aod::HfMlB0ToDPi, aod::HfSelB0ToDPi, aod::HfMcRecRedB0s>> const& candidates,
                         aod::HfMcGenRedB0s const& mcParticles,
                         CandsDplus const& candidatesD,
                         TracksPion const&)
  {
    // MC rec
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(hfHelper.yB0(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCand<true, false, false, true>(candidate, candidatesD);
    } // rec

    // MC gen. level
    for (const auto& particle : mcParticles) {
      fillCandMcGen(particle);
    } // gen
  } // processMcWithB0Ml
  PROCESS_SWITCH(HfTaskB0Reduced, processMcWithB0Ml, "Process MC with(out) ML scores for B0 (D daughter)", false);

  void processMcWithB0MlAndDecayTypeCheck(soa::Filtered<soa::Join<aod::HfRedCandB0, aod::HfMlB0ToDPi, aod::HfSelB0ToDPi, aod::HfMcRecRedB0s, aod::HfMcCheckB0s>> const& candidates,
                                          aod::HfMcGenRedB0s const& mcParticles,
                                          CandsDplus const& candidatesD,
                                          TracksPion const&)
  {
    // MC rec
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(hfHelper.yB0(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCand<true, true, false, true>(candidate, candidatesD);
    } // rec

    // MC gen. level
    for (const auto& particle : mcParticles) {
      fillCandMcGen(particle);
    } // gen
  } // processMc
  PROCESS_SWITCH(HfTaskB0Reduced, processMcWithB0MlAndDecayTypeCheck, "Process MC with decay type check and with(out) ML scores for B0 (D daughter)", false);
}; // struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskB0Reduced>(cfgc)};
}
