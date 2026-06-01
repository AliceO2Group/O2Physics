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

/// \file taskLbReduced.cxx
/// \brief Lb → Lc+ π- → (pK- π+) π- analysis task
///
/// \author Biao Zhang <biao.zhang@cern.ch>, Heidelberg University

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/D2H/DataModel/ReducedDataModel.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

#include "Common/Core/RecoDecay.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/runDataProcessing.h>

#include <TH3.h>
#include <TString.h>

#include <Rtypes.h>

#include <array>
#include <cstdint>
#include <numeric>

using namespace o2;
using namespace o2::aod;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{
namespace hf_cand_lb_lite
{
DECLARE_SOA_COLUMN(PtLc, ptLc, float);                                               //! Transverse momentum of Lc-baryon daughter candidate (GeV/c)
DECLARE_SOA_COLUMN(PtBach, ptBach, float);                                           //! Transverse momentum of bachelor pion (GeV/c)
DECLARE_SOA_COLUMN(AbsEtaBach, absEtaBach, float);                                   //! Absolute pseudorapidity of bachelor pion
DECLARE_SOA_COLUMN(ItsNClsBach, itsNClsBach, int);                                   //! Number of ITS clusters of bachelor pion
DECLARE_SOA_COLUMN(TpcNClsCrossedRowsBach, tpcNClsCrossedRowsBach, int);             //! Number of TPC crossed rows of prongs of bachelor pion
DECLARE_SOA_COLUMN(TpcChi2NClBach, tpcChi2NClBach, float);                           //! Maximum TPC chi2 of prongs of Lc-baryon daughter candidate
DECLARE_SOA_COLUMN(PtLcProngMin, ptLcProngMin, float);                               //! Minimum pT of prongs of Lc-baryon daughter candidate (GeV/c)
DECLARE_SOA_COLUMN(EtaLcProngMin, etaLcProngMin, float);                             //! Minimum absolute pseudorapidity of prongs of Lc-baryon daughter candidate
DECLARE_SOA_COLUMN(ItsNClsLcProngMin, itsNClsLcProngMin, int);                       //! Minimum number of ITS clusters of prongs of Lc-baryon daughter candidate
DECLARE_SOA_COLUMN(TpcNClsCrossedRowsLcProngMin, tpcNClsCrossedRowsLcProngMin, int); //! Minimum number of TPC crossed rows of prongs of Lc-baryon daughter candidate
DECLARE_SOA_COLUMN(TpcChi2NClLcProngMax, tpcChi2NClLcProngMax, float);               //! Maximum TPC chi2 of prongs of Lc-baryon daughter candidate
DECLARE_SOA_COLUMN(MLc, mLc, float);                                                 //! Invariant mass of Lc-baryon daughter candidates (GeV/c)
DECLARE_SOA_COLUMN(M, m, float);                                                     //! Invariant mass of candidate (GeV/c2)
DECLARE_SOA_COLUMN(Pt, pt, float);                                                   //! Transverse momentum of candidate (GeV/c)
DECLARE_SOA_COLUMN(PtGen, ptGen, float);                                             //! Transverse momentum of candidate (GeV/c)
DECLARE_SOA_COLUMN(P, p, float);                                                     //! Momentum of candidate (GeV/c)
DECLARE_SOA_COLUMN(Y, y, float);                                                     //! Rapidity of candidate
DECLARE_SOA_COLUMN(Eta, eta, float);                                                 //! Pseudorapidity of candidate
DECLARE_SOA_COLUMN(Phi, phi, float);                                                 //! Azimuth angle of candidate
DECLARE_SOA_COLUMN(E, e, float);                                                     //! Energy of candidate (GeV)
DECLARE_SOA_COLUMN(NSigTpcPiBachelor, nSigTpcPiBachelor, float);                     //! TPC Nsigma separation for bachelor with pion mass hypothesis
DECLARE_SOA_COLUMN(NSigTofPiBachelor, nSigTofPiBachelor, float);                     //! TOF Nsigma separation for bachelor with pion mass hypothesis
DECLARE_SOA_COLUMN(NSigTpcTofPiBachelor, nSigTpcTofPiBachelor, float);               //! Combined TPC and TOF Nsigma separation for bachelor with pion mass hypothesis
DECLARE_SOA_COLUMN(NSigTpcPrLcProng0, nSigTpcPrLcProng0, float);                     //! TPC Nsigma separation for Lc-baryon prong0 with proton mass hypothesis
DECLARE_SOA_COLUMN(NSigTofPrLcProng0, nSigTofPrLcProng0, float);                     //! TOF Nsigma separation for Lc-baryon prong0 with proton mass hypothesis
DECLARE_SOA_COLUMN(NSigTpcTofPrLcProng0, nSigTpcTofPrLcProng0, float);               //! Combined TPC and TOF Nsigma separation for Lc-baryon prong0 with proton mass hypothesis
DECLARE_SOA_COLUMN(NSigTpcKaLcProng1, nSigTpcKaLcProng1, float);                     //! TPC Nsigma separation for Lc-baryon prong1 with kaon mass hypothesis
DECLARE_SOA_COLUMN(NSigTofKaLcProng1, nSigTofKaLcProng1, float);                     //! TOF Nsigma separation for Lc-baryon prong1 with kaon mass hypothesis
DECLARE_SOA_COLUMN(NSigTpcTofKaLcProng1, nSigTpcTofKaLcProng1, float);               //! Combined TPC and TOF Nsigma separation for Lc-baryon prong1 with kaon mass hypothesis
DECLARE_SOA_COLUMN(NSigTpcPiLcProng2, nSigTpcPiLcProng2, float);                     //! TPC Nsigma separation for Lc-baryon prong2 with pion mass hypothesis
DECLARE_SOA_COLUMN(NSigTofPiLcProng2, nSigTofPiLcProng2, float);                     //! TOF Nsigma separation for Lc-baryon prong2 with pion mass hypothesis
DECLARE_SOA_COLUMN(NSigTpcTofPiLcProng2, nSigTpcTofPiLcProng2, float);               //! Combined TPC and TOF Nsigma separation for Lc-baryon prong0 with pion mass hypothesis
DECLARE_SOA_COLUMN(DecayLength, decayLength, float);                                 //! Decay length of candidate (cm)
DECLARE_SOA_COLUMN(DecayLengthXY, decayLengthXY, float);                             //! Transverse decay length of candidate (cm)
DECLARE_SOA_COLUMN(DecayLengthNormalised, decayLengthNormalised, float);             //! Normalised decay length of candidate
DECLARE_SOA_COLUMN(DecayLengthXYNormalised, decayLengthXYNormalised, float);         //! Normalised transverse decay length of candidate
DECLARE_SOA_COLUMN(DecayLengthLc, decayLengthLc, float);                             //! Decay length of Lc-baryon daughter candidate (cm)
DECLARE_SOA_COLUMN(DecayLengthXYLc, decayLengthXYLc, float);                         //! Transverse decay length of Lc-baryon daughter candidate (cm)
DECLARE_SOA_COLUMN(ImpactParameterLc, impactParameterLc, float);                     //! Impact parameter product of Lc-baryon daughter candidate
DECLARE_SOA_COLUMN(ImpactParameterBach, impactParameterBach, float);                 //! Impact parameter product of bachelor pion
DECLARE_SOA_COLUMN(ImpactParameterProduct, impactParameterProduct, float);           //! Impact parameter product of daughters
DECLARE_SOA_COLUMN(Cpa, cpa, float);                                                 //! Cosine pointing angle of candidate
DECLARE_SOA_COLUMN(CpaXY, cpaXY, float);                                             //! Cosine pointing angle of candidate in transverse plane
DECLARE_SOA_COLUMN(MaxNormalisedDeltaIP, maxNormalisedDeltaIP, float);               //! Maximum normalized difference between measured and expected impact parameter of candidate prongs
DECLARE_SOA_COLUMN(MlScoreSig, mlScoreSig, float);                                   //! ML score for signal class
DECLARE_SOA_COLUMN(FlagWrongCollision, flagWrongCollision, int8_t);                  //! Flag for association with wrong collision
} // namespace hf_cand_lb_lite

DECLARE_SOA_TABLE(HfRedCandLbLites, "AOD", "HFREDCANDLBLITE", //! Table with some Lb properties
                                                              // B meson features                  hf_cand_lb_lite::M,
                  hf_cand_lb_lite::M,
                  hf_cand_lb_lite::Pt,
                  hf_cand_lb_lite::Eta,
                  hf_cand_lb_lite::Phi,
                  hf_cand_lb_lite::Y,
                  hf_cand_lb_lite::Cpa,
                  hf_cand_lb_lite::CpaXY,
                  hf_cand::Chi2PCA,
                  hf_cand_lb_lite::DecayLength,
                  hf_cand_lb_lite::DecayLengthXY,
                  hf_cand_lb_lite::DecayLengthNormalised,
                  hf_cand_lb_lite::DecayLengthXYNormalised,
                  hf_cand_lb_lite::ImpactParameterProduct,
                  hf_cand_lb_lite::MaxNormalisedDeltaIP,
                  hf_cand_lb_lite::MlScoreSig,
                  hf_sel_candidate_lb::IsSelLbToLcPi,
                  // Lc baryon features
                  hf_cand_lb_lite::MLc,
                  hf_cand_lb_lite::PtLc,
                  hf_cand_lb_lite::DecayLengthLc,
                  hf_cand_lb_lite::DecayLengthXYLc,
                  hf_cand_lb_lite::ImpactParameterLc,
                  hf_cand_lb_lite::PtLcProngMin,
                  hf_cand_lb_lite::EtaLcProngMin,
                  hf_cand_lb_lite::ItsNClsLcProngMin,
                  hf_cand_lb_lite::TpcNClsCrossedRowsLcProngMin,
                  hf_cand_lb_lite::TpcChi2NClLcProngMax,
                  hf_cand_lb_lite::NSigTpcPrLcProng0,
                  hf_cand_lb_lite::NSigTofPrLcProng0,
                  hf_cand_lb_lite::NSigTpcTofPrLcProng0,
                  hf_cand_lb_lite::NSigTpcKaLcProng1,
                  hf_cand_lb_lite::NSigTofKaLcProng1,
                  hf_cand_lb_lite::NSigTpcTofKaLcProng1,
                  hf_cand_lb_lite::NSigTpcPiLcProng2,
                  hf_cand_lb_lite::NSigTofPiLcProng2,
                  hf_cand_lb_lite::NSigTpcTofPiLcProng2,
                  hf_cand_lb_reduced::Prong0MlScoreBkg,
                  hf_cand_lb_reduced::Prong0MlScorePrompt,
                  hf_cand_lb_reduced::Prong0MlScoreNonprompt,
                  // pion features
                  hf_cand_lb_lite::PtBach,
                  hf_cand_lb_lite::AbsEtaBach,
                  hf_cand_lb_lite::ItsNClsBach,
                  hf_cand_lb_lite::TpcNClsCrossedRowsBach,
                  hf_cand_lb_lite::TpcChi2NClBach,
                  hf_cand_lb_lite::ImpactParameterBach,
                  hf_cand_lb_lite::NSigTpcPiBachelor,
                  hf_cand_lb_lite::NSigTofPiBachelor,
                  hf_cand_lb_lite::NSigTpcTofPiBachelor,
                  // MC truth
                  hf_cand_mc_flag::FlagMcMatchRec,
                  hf_cand_mc_flag::OriginMcRec,
                  hf_cand_lb_lite::FlagWrongCollision,
                  hf_cand_lb_lite::PtGen);

DECLARE_SOA_TABLE(HfRedLbMcCheck, "AOD", "HFREDLBMCCHECK", //! Table with MC decay type check
                  hf_cand_mc_flag::FlagMcMatchRec,
                  hf_cand_lb_lite::FlagWrongCollision,
                  hf_cand_lb_lite::MLc,
                  hf_cand_lb_lite::PtLc,
                  hf_cand_lb_lite::M,
                  hf_cand_lb_lite::Pt,
                  hf_cand_lb_lite::MlScoreSig,
                  hf_lb_mc::PdgCodeBeautyMother,
                  hf_lb_mc::PdgCodeCharmMother,
                  hf_lb_mc::PdgCodeProng0,
                  hf_lb_mc::PdgCodeProng1,
                  hf_lb_mc::PdgCodeProng2,
                  hf_lb_mc::PdgCodeProng3);
} // namespace o2::aod

/// Lb analysis task
struct HfTaskLbReduced {
  Produces<aod::HfRedCandLbLites> hfRedCandLbLite;
  Produces<aod::HfRedLbMcCheck> hfRedLbMcCheck;

  Configurable<int> selectionFlagLb{"selectionFlagLb", 1, "Selection Flag for Lb"};
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

  using TracksPion = soa::Join<HfRedTracks, HfRedTracksPid>;
  using CandsLc = soa::Join<HfRed3Prongs, HfRedPidDau0s, HfRedPidDau1s, HfRedPidDau2s>;
  Filter filterSelectCandidates = (aod::hf_sel_candidate_lb::isSelLbToLcPi >= selectionFlagLb);
  HistogramRegistry registry{"registry"};

  void init(InitContext&)
  {
    std::array<bool, 3> processFuncData{doprocessData, doprocessDataWithLcMl, doprocessDataWithLbMl};
    if ((std::accumulate(processFuncData.begin(), processFuncData.end(), 0)) > 1) {
      LOGP(fatal, "Only one process function for data can be enabled at a time.");
    }
    std::array<bool, 6> processFuncMc{doprocessMc, doprocessMcWithDecayTypeCheck, doprocessMcWithLcMl, doprocessMcWithLcMlAndDecayTypeCheck, doprocessMcWithLbMl, doprocessMcWithLbMlAndDecayTypeCheck};
    if ((std::accumulate(processFuncMc.begin(), processFuncMc.end(), 0)) > 1) {
      LOGP(fatal, "Only one process function for MC can be enabled at a time.");
    }

    const AxisSpec axisMlScore{100, 0.f, 1.f};
    const AxisSpec axisMassLb{300, 4.5f, 6.0f};
    const AxisSpec axisMassLc{300, 2.15f, 2.45f};
    const AxisSpec axisDecayLength{200, 0.f, 0.4f};
    const AxisSpec axisNormDecayLength{100, 0.f, 50.f};
    const AxisSpec axisDca{100, -0.05f, 0.05f};
    const AxisSpec axisCosp{110, 0.f, 1.1f};
    const AxisSpec axisEta{30, -1.5f, 1.5f};
    const AxisSpec axisError{100, 0.f, 1.f};
    const AxisSpec axisImpParProd{100, -1.e-3, 1.e-3};
    const AxisSpec axisPtLb{100, 0.f, 50.f};
    const AxisSpec axisPtLc{100, 0.f, 50.f};
    const AxisSpec axisPtPi{100, 0.f, 10.f};

    if (doprocessData || doprocessDataWithLcMl || doprocessDataWithLbMl) {
      if (fillHistograms) {
        registry.add("hMass", "#Lambda_{b}^{0} candidates;#it{p}_{T}(#Lambda_{b}^{0}) (GeV/#it{c});#it{M} (D#pi) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {axisPtLb, axisMassLb}});
        registry.add("hDecLength", "#Lambda_{b}^{0} candidates;#it{p}_{T}(#Lambda_{b}^{0}) (GeV/#it{c});#Lambda_{b}^{0} candidate decay length (cm);entries", {HistType::kTH2F, {axisPtLb, axisDecayLength}});
        registry.add("hDecLengthXy", "#Lambda_{b}^{0} candidates;#it{p}_{T}(#Lambda_{b}^{0}) (GeV/#it{c});decay length XY (cm);entries", {HistType::kTH2F, {axisPtLb, axisDecayLength}});
        registry.add("hNormDecLengthXy", "#Lambda_{b}^{0} candidates;#it{p}_{T}(#Lambda_{b}^{0}) (GeV/#it{c});#Lambda_{b}^{0} candidate norm. decay length XY (cm);entries", {HistType::kTH2F, {axisPtLb, axisNormDecayLength}});
        registry.add("hDcaProng0", "#Lambda_{b}^{0} candidates;#it{p}_{T}(#Lambda_{b}^{0}) (GeV/#it{c});prong 0 (#Lambda_{c}^{+}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {axisPtLb, axisDca}});
        registry.add("hDcaProng1", "#Lambda_{b}^{0} candidates;#it{p}_{T}(#Lambda_{b}^{0}) (GeV/#it{c});prong 1 (#pi^{#plus}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {axisPtLb, axisDca}});
        registry.add("hPtProng0", "#Lambda_{b}^{0} candidates;#it{p}_{T}(#Lambda_{b}^{0}) (GeV/#it{c});#it{p}_{T}(#Lambda_{c}^{#plus}) (GeV/#it{c});entries", {HistType::kTH2F, {axisPtLb, axisPtLc}});
        registry.add("hPtProng1", "#Lambda_{b}^{0} candidates;#it{p}_{T}(#Lambda_{b}^{0}) (GeV/#it{c});#it{p}_{T}(#pi^{#plus}) (GeV/#it{c});entries", {HistType::kTH2F, {axisPtLb, axisPtPi}});
        registry.add("hCosp", "#Lambda_{b}^{0} candidates;#it{p}_{T}(#Lambda_{b}^{0}) (GeV/#it{c});#Lambda_{b}^{0} candidate cos(#vartheta_{P});entries", {HistType::kTH2F, {axisPtLb, axisCosp}});
        registry.add("hCospXy", "#Lambda_{b}^{0} candidates;#it{p}_{T}(#Lambda_{b}^{0}) (GeV/#it{c});#Lambda_{b}^{0} candidate cos(#vartheta_{P}^{XY});entries", {HistType::kTH2F, {axisPtLb, axisCosp}});
        registry.add("hEta", "#Lambda_{b}^{0} candidates;#it{p}_{T}(#Lambda_{b}^{0}) (GeV/#it{c});#Lambda_{b}^{0} candidate #it{#eta};entries", {HistType::kTH2F, {axisPtLb, axisEta}});
        registry.add("hRapidity", "#Lambda_{b}^{0} candidates;#it{p}_{T}(#Lambda_{b}^{0}) (GeV/#it{c});#Lambda_{b}^{0} candidate #it{y};entries", {HistType::kTH2F, {axisPtLb, axisEta}});
        registry.add("hImpParProd", "#Lambda_{b}^{0} candidates;#it{p}_{T}(#Lambda_{b}^{0}) (GeV/#it{c});#Lambda_{b}^{0} candidate impact parameter product;entries", {HistType::kTH2F, {axisPtLb, axisImpParProd}});
        registry.add("hinvMassLc", "#Lambda_{b}^{0} candidates;#it{p}_{T}(#Lambda_{c}^{#plus}) (GeV/#it{c});prong0, #it{M}(pK#pi) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {axisPtLc, axisMassLc}});
        registry.add("hDecLengthLc", "#Lambda_{b}^{0} candidates;#it{p}_{T}(#Lambda_{c}^{#plus}) (GeV/#it{c});#Lambda_{c}^{+} candidate decay length (cm);entries", {HistType::kTH2F, {axisPtLc, axisDecayLength}});
        registry.add("hDecLengthXyLc", "#Lambda_{b}^{0} candidates;#it{p}_{T}(#Lambda_{c}^{#plus}) (GeV/#it{c});decay length XY (cm);entries", {HistType::kTH2F, {axisPtLc, axisDecayLength}});
        registry.add("hCospLc", "#Lambda_{b}^{0} candidates;#it{p}_{T}(#Lambda_{c}^{#plus}) (GeV/#it{c});#Lambda_{c}^{+} candidate cos(#vartheta_{P});entries", {HistType::kTH2F, {axisPtLc, axisCosp}});
        registry.add("hCospXyLc", "#Lambda_{b}^{0} candidates;#it{p}_{T}(#Lambda_{c}^{#plus}) (GeV/#it{c});#Lambda_{c}^{+} candidate cos(#vartheta_{P}^{XY});entries", {HistType::kTH2F, {axisPtLc, axisCosp}});

        // ML scores of Lc daughter
        if (doprocessDataWithLcMl) {
          registry.add("hMlScoreBkgLc", "#Lambda_{b}^{0} candidates;#it{p}_{T}(#Lambda_{c}^{#plus}) (GeV/#it{c});prong0, #Lambda_{c}^{+} ML background score;entries", {HistType::kTH2F, {axisPtLc, axisMlScore}});
          registry.add("hMlScorePromptLc", "#Lambda_{b}^{0} candidates;#it{p}_{T}(#Lambda_{c}^{#plus}) (GeV/#it{c});prong0, #Lambda_{c}^{+} ML prompt score;entries", {HistType::kTH2F, {axisPtLc, axisMlScore}});
          registry.add("hMlScoreNonPromptLc", "#Lambda_{b}^{0} candidates;#it{p}_{T}(#Lambda_{c}^{#plus}) (GeV/#it{c});prong0, #Lambda_{c}^{+} ML nonprompt score;entries", {HistType::kTH2F, {axisPtLc, axisMlScore}});
        }

        // ML scores of Lb candidate
        if (doprocessDataWithLbMl) {
          registry.add("hMlScoreSigLb", "#Lambda_{b}^{0} candidates;#it{p}_{T}(#Lambda_{b}^{0}) (GeV/#it{c});prong0, #Lambda_{b}^{0} ML signal score;entries", {HistType::kTH2F, {axisPtLb, axisMlScore}});
        }
      }
      if (fillSparses) {
        if (!(doprocessDataWithLcMl || doprocessDataWithLbMl)) {
          registry.add("hMassPtCutVars", "#Lambda_{b}^{0} candidates;#it{M} (D#pi) (GeV/#it{c}^{2});#it{p}_{T}(#Lambda_{b}^{0}) (GeV/#it{c});#Lambda_{b}^{0} candidate decay length (cm);#Lambda_{b}^{0} candidate norm. decay length XY (cm);#Lambda_{b}^{0} candidate impact parameter product (cm);#Lambda_{b}^{0} candidate cos(#vartheta_{P});#it{M} (pK#pi) (GeV/#it{c}^{2});#it{p}_{T}(#Lambda_{c}^{#plus}) (GeV/#it{c});#Lambda_{c}^{+} candidate decay length (cm);#Lambda_{c}^{+} candidate cos(#vartheta_{P})", {HistType::kTHnSparseF, {axisMassLb, axisPtLb, axisDecayLength, axisNormDecayLength, axisImpParProd, axisCosp, axisMassLc, axisPtLc, axisDecayLength, axisCosp}});
        } else {
          registry.add("hMassPtCutVars", "#Lambda_{b}^{0} candidates;#it{M} (D#pi) (GeV/#it{c}^{2});#it{p}_{T}(#Lambda_{b}^{0}) (GeV/#it{c});#Lambda_{b}^{0} candidate decay length (cm);#Lambda_{b}^{0} candidate norm. decay length XY (cm);#Lambda_{b}^{0} candidate impact parameter product (cm);#Lambda_{b}^{0} candidate cos(#vartheta_{P});#it{M} (pK#pi) (GeV/#it{c}^{2});#it{p}_{T}(#Lambda_{c}^{#plus}) (GeV/#it{c});#Lambda_{c}^{+} candidate ML score bkg;#Lambda_{c}^{+} candidate ML score nonprompt", {HistType::kTHnSparseF, {axisMassLb, axisPtLb, axisDecayLength, axisNormDecayLength, axisImpParProd, axisCosp, axisMassLc, axisPtLc, axisMlScore, axisMlScore}});
        }
      }
    }

    if (doprocessMc || doprocessMcWithDecayTypeCheck || doprocessMcWithLcMl || doprocessMcWithLcMlAndDecayTypeCheck || doprocessMcWithLbMl || doprocessMcWithLbMlAndDecayTypeCheck) {
      if (fillHistograms) {
        // gen histos
        registry.add("hEtaGen", "#Lambda_{b}^{0} particles (generated);#it{p}_{T}^{gen}(#Lambda_{b}^{0}) (GeV/#it{c});#it{#eta}^{gen}(#Lambda_{b}^{0});entries", {HistType::kTH2F, {axisPtLb, axisEta}});
        registry.add("hYGen", "#Lambda_{b}^{0} particles (generated);#it{p}_{T}^{gen}(#Lambda_{b}^{0}) (GeV/#it{c});#it{y}^{gen}(#Lambda_{b}^{0});entries", {HistType::kTH2F, {axisPtLb, axisEta}});
        registry.add("hYGenWithProngsInAcceptance", "MC particles (generated-daughters in acceptance);#it{p}_{T}^{gen}(#Lambda_{b}^{0}) (GeV/#it{c});#it{y}^{gen}(#Lambda_{b}^{0});entries", {HistType::kTH2F, {axisPtLb, axisEta}});
        registry.add("hPtProng0Gen", "#Lambda_{b}^{0} particles (generated);#it{p}_{T}^{gen}(#Lambda_{b}^{0}) (GeV/#it{c});#it{p}_{T}^{gen}(#Lambda_{c}^{+}) (GeV/#it{c});entries", {HistType::kTH2F, {axisPtLb, axisPtLc}});
        registry.add("hPtProng1Gen", "#Lambda_{b}^{0} particles (generated);#it{p}_{T}^{gen}(#Lambda_{b}^{0}) (GeV/#it{c});#it{p}_{T}^{gen}(#pi^{#plus}) (GeV/#it{c});entries", {HistType::kTH2F, {axisPtLb, axisPtPi}});
        registry.add("hYProng0Gen", "#Lambda_{b}^{0} particles (generated);#it{p}_{T}^{gen}(#Lambda_{b}^{0}) (GeV/#it{c});#it{y}^{gen}(#Lambda_{c}^{+});entries", {HistType::kTH2F, {axisPtLb, axisEta}});
        registry.add("hYProng1Gen", "#Lambda_{b}^{0} particles (generated);#it{p}_{T}^{gen}(#Lambda_{b}^{0}) (GeV/#it{c});#it{y}^{gen}(#pi^{#plus});entries", {HistType::kTH2F, {axisPtLb, axisEta}});
        registry.add("hEtaProng0Gen", "#Lambda_{b}^{0} particles (generated);#it{p}_{T}^{gen}(#Lambda_{b}^{0}) (GeV/#it{c});#it{#eta}^{gen}(#Lambda_{c}^{+});entries", {HistType::kTH2F, {axisPtLb, axisEta}});
        registry.add("hEtaProng1Gen", "#Lambda_{b}^{0} particles (generated);#it{p}_{T}^{gen}(#Lambda_{b}^{0}) (GeV/#it{c});#it{#eta}^{gen}(#pi^{#plus});entries", {HistType::kTH2F, {axisPtLb, axisEta}});

        // reco histos
        // signal
        registry.add("hMassRecSig", "#Lambda_{b}^{0} candidates (matched);#it{p}_{T}(#Lambda_{b}^{0}) (GeV/#it{c});#it{M} (D#pi) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {axisPtLb, axisMassLb}});
        registry.add("hDecLengthRecSig", "#Lambda_{b}^{0} candidates (matched);#it{p}_{T}(#Lambda_{b}^{0}) (GeV/#it{c});#Lambda_{b}^{0} candidate decay length (cm);entries", {HistType::kTH2F, {axisPtLb, axisDecayLength}});
        registry.add("hDecLengthXyRecSig", "#Lambda_{b}^{0} candidates (matched);#it{p}_{T}(#Lambda_{b}^{0}) (GeV/#it{c});decay length XY (cm);entries", {HistType::kTH2F, {axisPtLb, axisDecayLength}});
        registry.add("hNormDecLengthXyRecSig", "#Lambda_{b}^{0} candidates (matched);#it{p}_{T}(#Lambda_{b}^{0}) (GeV/#it{c});#Lambda_{b}^{0} candidate norm. decay length XY (cm);entries", {HistType::kTH2F, {axisPtLb, axisNormDecayLength}});
        registry.add("hDcaProng0RecSig", "#Lambda_{b}^{0} candidates (matched);#it{p}_{T}(#Lambda_{b}^{0}) (GeV/#it{c});prong 0 (#Lambda_{c}^{+}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {axisPtLb, axisDca}});
        registry.add("hDcaProng1RecSig", "#Lambda_{b}^{0} candidates (matched);#it{p}_{T}(#Lambda_{b}^{0}) (GeV/#it{c});prong 1 (#pi^{#plus}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {axisPtLb, axisDca}});
        registry.add("hPtProng0RecSig", "#Lambda_{b}^{0} candidates (matched);#it{p}_{T}(#Lambda_{b}^{0}) (GeV/#it{c});#it{p}_{T}(#Lambda_{c}^{#plus}) (GeV/#it{c});entries", {HistType::kTH2F, {axisPtLb, axisPtLc}});
        registry.add("hPtProng1RecSig", "#Lambda_{b}^{0} candidates (matched);#it{p}_{T}(#Lambda_{b}^{0}) (GeV/#it{c});#it{p}_{T}(#pi^{#plus}) (GeV/#it{c});entries", {HistType::kTH2F, {axisPtLb, axisPtPi}});
        registry.add("hCospRecSig", "#Lambda_{b}^{0} candidates (matched);#it{p}_{T}(#Lambda_{b}^{0}) (GeV/#it{c});#Lambda_{b}^{0} candidate cos(#vartheta_{P});entries", {HistType::kTH2F, {axisPtLb, axisCosp}});
        registry.add("hCospXyRecSig", "#Lambda_{b}^{0} candidates (matched);#it{p}_{T}(#Lambda_{b}^{0}) (GeV/#it{c});#Lambda_{b}^{0} candidate cos(#vartheta_{P}^{XY});entries", {HistType::kTH2F, {axisPtLb, axisCosp}});
        registry.add("hEtaRecSig", "#Lambda_{b}^{0} candidates (matched);#it{p}_{T}(#Lambda_{b}^{0}) (GeV/#it{c});#Lambda_{b}^{0} candidate #it{#eta};entries", {HistType::kTH2F, {axisPtLb, axisEta}});
        registry.add("hRapidityRecSig", "#Lambda_{b}^{0} candidates (matched);#it{p}_{T}(#Lambda_{b}^{0}) (GeV/#it{c});#Lambda_{b}^{0} candidate #it{y};entries", {HistType::kTH2F, {axisPtLb, axisEta}});
        registry.add("hImpParProdRecSig", "#Lambda_{b}^{0} candidates (matched);#it{p}_{T}(#Lambda_{b}^{0}) (GeV/#it{c});#Lambda_{b}^{0} candidate impact parameter product;entries", {HistType::kTH2F, {axisPtLb, axisImpParProd}});
        registry.add("hinvMassLcRecSig", "#Lambda_{b}^{0} candidates (matched);#it{p}_{T}(#Lambda_{c}^{#plus}) (GeV/#it{c});prong0, #it{M}(pK#pi) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {axisPtLc, axisMassLc}});
        registry.add("hDecLengthLcRecSig", "#Lambda_{b}^{0} candidates (matched);#it{p}_{T}(#Lambda_{c}^{#plus}) (GeV/#it{c});#Lambda_{c}^{+} candidate decay length (cm);entries", {HistType::kTH2F, {axisPtLc, axisDecayLength}});
        registry.add("hDecLengthXyLcRecSig", "#Lambda_{b}^{0} candidates (matched);#it{p}_{T}(#Lambda_{c}^{#plus}) (GeV/#it{c});decay length XY (cm);entries", {HistType::kTH2F, {axisPtLc, axisDecayLength}});
        registry.add("hCospLcRecSig", "#Lambda_{b}^{0} candidates (matched);#it{p}_{T}(#Lambda_{c}^{#plus}) (GeV/#it{c});#Lambda_{c}^{+} candidate cos(#vartheta_{P});entries", {HistType::kTH2F, {axisPtLc, axisCosp}});
        registry.add("hCospXyLcRecSig", "#Lambda_{b}^{0} candidates (matched);#it{p}_{T}(#Lambda_{c}^{#plus}) (GeV/#it{c});#Lambda_{c}^{+} candidate cos(#vartheta_{P}^{XY});entries", {HistType::kTH2F, {axisPtLc, axisCosp}});
        // background
        if (fillBackground) {
          registry.add("hMassRecBg", "#Lambda_{b}^{0} candidates (unmatched);#it{p}_{T}(#Lambda_{b}^{0}) (GeV/#it{c});#it{M} (D#pi) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {axisPtLb, axisMassLb}});
          registry.add("hDecLengthRecBg", "#Lambda_{b}^{0} candidates (unmatched);#it{p}_{T}(#Lambda_{b}^{0}) (GeV/#it{c});#Lambda_{b}^{0} candidate decay length (cm);entries", {HistType::kTH2F, {axisPtLb, axisDecayLength}});
          registry.add("hDecLengthXyRecBg", "#Lambda_{b}^{0} candidates (unmatched);#it{p}_{T}(#Lambda_{b}^{0}) (GeV/#it{c});decay length XY (cm);entries", {HistType::kTH2F, {axisPtLb, axisDecayLength}});
          registry.add("hNormDecLengthXyRecBg", "#Lambda_{b}^{0} candidates (unmatched);#it{p}_{T}(#Lambda_{b}^{0}) (GeV/#it{c});#Lambda_{b}^{0} candidate norm. decay length XY (cm);entries", {HistType::kTH2F, {axisPtLb, axisNormDecayLength}});
          registry.add("hDcaProng0RecBg", "#Lambda_{b}^{0} candidates (unmatched);#it{p}_{T}(#Lambda_{b}^{0}) (GeV/#it{c});prong 0 (#Lambda_{c}^{+}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {axisPtLb, axisDca}});
          registry.add("hDcaProng1RecBg", "#Lambda_{b}^{0} candidates (unmatched);#it{p}_{T}(#Lambda_{b}^{0}) (GeV/#it{c});prong 1 (#pi^{#plus}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {axisPtLb, axisDca}});
          registry.add("hPtProng0RecBg", "#Lambda_{b}^{0} candidates (unmatched);#it{p}_{T}(#Lambda_{b}^{0}) (GeV/#it{c});#it{p}_{T}(#Lambda_{c}^{#plus}) (GeV/#it{c});entries", {HistType::kTH2F, {axisPtLb, axisPtLc}});
          registry.add("hPtProng1RecBg", "#Lambda_{b}^{0} candidates (unmatched);#it{p}_{T}(#Lambda_{b}^{0}) (GeV/#it{c});#it{p}_{T}(#pi^{#plus}) (GeV/#it{c});entries", {HistType::kTH2F, {axisPtLb, axisPtPi}});
          registry.add("hCospRecBg", "#Lambda_{b}^{0} candidates (unmatched);#it{p}_{T}(#Lambda_{b}^{0}) (GeV/#it{c});#Lambda_{b}^{0} candidate cos(#vartheta_{P});entries", {HistType::kTH2F, {axisPtLb, axisCosp}});
          registry.add("hCospXyRecBg", "#Lambda_{b}^{0} candidates (unmatched);#it{p}_{T}(#Lambda_{b}^{0}) (GeV/#it{c});#Lambda_{b}^{0} candidate cos(#vartheta_{P}^{XY});entries", {HistType::kTH2F, {axisPtLb, axisCosp}});
          registry.add("hEtaRecBg", "#Lambda_{b}^{0} candidates (unmatched);#it{p}_{T}(#Lambda_{b}^{0}) (GeV/#it{c});#Lambda_{b}^{0} candidate #it{#eta};entries", {HistType::kTH2F, {axisPtLb, axisEta}});
          registry.add("hRapidityRecBg", "#Lambda_{b}^{0} candidates (unmatched);#it{p}_{T}(#Lambda_{b}^{0}) (GeV/#it{c});#Lambda_{b}^{0} candidate #it{y};entries", {HistType::kTH2F, {axisPtLb, axisEta}});
          registry.add("hImpParProdRecBg", "#Lambda_{b}^{0} candidates (unmatched);#it{p}_{T}(#Lambda_{b}^{0}) (GeV/#it{c});#Lambda_{b}^{0} candidate impact parameter product;entries", {HistType::kTH2F, {axisPtLb, axisImpParProd}});
          registry.add("hinvMassLcRecBg", "#Lambda_{b}^{0} candidates (unmatched);#it{p}_{T}(#Lambda_{c}^{#plus}) (GeV/#it{c});prong0, #it{M}(pK#pi) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {axisPtLc, axisMassLc}});
          registry.add("hDecLengthLcRecBg", "#Lambda_{b}^{0} candidates (unmatched);#it{p}_{T}(#Lambda_{c}^{#plus}) (GeV/#it{c});#Lambda_{c}^{+} candidate decay length (cm);entries", {HistType::kTH2F, {axisPtLc, axisDecayLength}});
          registry.add("hDecLengthXyLcRecBg", "#Lambda_{b}^{0} candidates (unmatched);#it{p}_{T}(#Lambda_{c}^{#plus}) (GeV/#it{c});decay length XY (cm);entries", {HistType::kTH2F, {axisPtLc, axisDecayLength}});
          registry.add("hCospLcRecBg", "#Lambda_{b}^{0} candidates (unmatched);#it{p}_{T}(#Lambda_{c}^{#plus}) (GeV/#it{c});#Lambda_{c}^{+} candidate cos(#vartheta_{P});entries", {HistType::kTH2F, {axisPtLc, axisCosp}});
          registry.add("hCospXyLcRecBg", "#Lambda_{b}^{0} candidates (unmatched);#it{p}_{T}(#Lambda_{c}^{#plus}) (GeV/#it{c});#Lambda_{c}^{+} candidate cos(#vartheta_{P}^{XY});entries", {HistType::kTH2F, {axisPtLc, axisCosp}});
        }
        // MC checks
        if (doprocessMcWithDecayTypeCheck || doprocessMcWithLbMlAndDecayTypeCheck || doprocessMcWithLcMlAndDecayTypeCheck) {
          constexpr uint8_t kNBinsDecayTypeMc = hf_cand_lb::DecayTypeMc::NDecayTypeMc;
          TString labels[kNBinsDecayTypeMc];
          labels[hf_cand_lb::DecayTypeMc::LbToLcPiToPKPiPi] = "#Lambda_{b}^{0} #rightarrow (#Lambda_{c}^{#plus} #rightarrow p K^{#minus} #pi^{#plus}) #pi^{#minus}";
          labels[hf_cand_lb::DecayTypeMc::B0ToDplusPiToPiKPiPi] = "B^{0} #rightarrow (D^{#minus} #rightarrow #pi^{#minus} K^{#plus} #pi^{#minus}) #pi^{#plus}";
          labels[hf_cand_lb::DecayTypeMc::LbToLcKToPKPiK] = "#Lambda_{b}^{0} #rightarrow (#Lambda_{c}^{#plus} #rightarrow p K^{#minus} #pi^{#plus}) K^{#minus}";
          labels[hf_cand_lb::DecayTypeMc::PartlyRecoDecay] = "Partly reconstructed decay channel";
          labels[hf_cand_lb::DecayTypeMc::OtherDecay] = "Other decays";
          static const AxisSpec axisDecayType = {kNBinsDecayTypeMc, 0.5, kNBinsDecayTypeMc + 0.5, ""};
          registry.add("hDecayTypeMc", "DecayType", {HistType::kTH3F, {axisDecayType, axisMassLb, axisPtLb}});
          for (uint8_t iBin = 0; iBin < kNBinsDecayTypeMc; ++iBin) {
            registry.get<TH3>(HIST("hDecayTypeMc"))->GetXaxis()->SetBinLabel(iBin + 1, labels[iBin]);
          }
        }
        // ML scores of Lc daughter
        if (doprocessMcWithLcMl || doprocessMcWithLcMlAndDecayTypeCheck) {
          // signal
          registry.add("hMlScoreBkgLcRecSig", "#Lambda_{b}^{0} candidates (matched);#it{p}_{T}(#Lambda_{c}^{#plus}) (GeV/#it{c});prong0, #Lambda_{c}^{+} ML background score;entries", {HistType::kTH2F, {axisPtLc, axisMlScore}});
          registry.add("hMlScorePromptLcRecSig", "#Lambda_{b}^{0} candidates (matched);#it{p}_{T}(#Lambda_{c}^{#plus}) (GeV/#it{c});prong0, #Lambda_{c}^{+} ML prompt score;entries", {HistType::kTH2F, {axisPtLc, axisMlScore}});
          registry.add("hMlScoreNonPromptLcRecSig", "#Lambda_{b}^{0} candidates (matched);#it{p}_{T}(#Lambda_{c}^{#plus}) (GeV/#it{c});prong0, #Lambda_{c}^{+} ML nonprompt score;entries", {HistType::kTH2F, {axisPtLc, axisMlScore}});
          // background
          registry.add("hMlScoreBkgLcRecBg", "#Lambda_{b}^{0} candidates (unmatched);#it{p}_{T}(#Lambda_{c}^{#plus}) (GeV/#it{c});prong0, #Lambda_{c}^{+} ML background score;entries", {HistType::kTH2F, {axisPtLc, axisMlScore}});
          registry.add("hMlScorePromptLcRecBg", "#Lambda_{b}^{0} candidates (unmatched);#it{p}_{T}(#Lambda_{c}^{#plus}) (GeV/#it{c});prong0, #Lambda_{c}^{+} ML prompt score;entries", {HistType::kTH2F, {axisPtLc, axisMlScore}});
          registry.add("hMlScoreNonPromptLcRecBg", "#Lambda_{b}^{0} candidates (unmatched);#it{p}_{T}(#Lambda_{c}^{#plus}) (GeV/#it{c});prong0, #Lambda_{c}^{+} ML nonprompt score;entries", {HistType::kTH2F, {axisPtLc, axisMlScore}});
        }
        // ML scores of Lb candidate
        if (doprocessMcWithLbMl || doprocessMcWithLbMlAndDecayTypeCheck) {
          // signal
          registry.add("hMlScoreSigLbRecSig", "#Lambda_{b}^{0} candidates (matched);#it{p}_{T}(#Lambda_{b}^{0}) (GeV/#it{c});prong0, #Lambda_{b}^{0} ML signal score;entries", {HistType::kTH2F, {axisPtLb, axisMlScore}});
          // background
          registry.add("hMlScoreSigLbRecBg", "#Lambda_{b}^{0} candidates (unmatched);#it{p}_{T}(#Lambda_{b}^{0}) (GeV/#it{c});prong0, #Lambda_{b}^{0} ML signal score;entries", {HistType::kTH2F, {axisPtLb, axisMlScore}});
        }
      }
      if (fillSparses) {
        // gen sparses
        registry.add("hPtYGenSig", "#Lambda_{b}^{0} particles (generated);#it{p}_{T}(#Lambda_{b}^{0}) (GeV/#it{c});#it{y}(#Lambda_{b}^{0})", {HistType::kTHnSparseF, {axisPtLb, axisEta}});
        registry.add("hPtYWithProngsInAccepanceGenSig", "#Lambda_{b}^{0} particles (generated-daughters in acceptance);#it{p}_{T}(#Lambda_{b}^{0}) (GeV/#it{c});#it{y}(#Lambda_{b}^{0})", {HistType::kTHnSparseF, {axisPtLb, axisEta}});

        // reco sparses
        if (!(doprocessDataWithLcMl || doprocessDataWithLbMl)) {
          registry.add("hMassPtCutVarsRecSig", "#Lambda_{b}^{0} candidates (matched);#it{M} (D#pi) (GeV/#it{c}^{2});#it{p}_{T}(#Lambda_{b}^{0}) (GeV/#it{c});#Lambda_{b}^{0} candidate decay length (cm);#Lambda_{b}^{0} candidate norm. decay length XY (cm);#Lambda_{b}^{0} candidate impact parameter product (cm);#Lambda_{b}^{0} candidate cos(#vartheta_{P});#it{M} (pK#pi) (GeV/#it{c}^{2});#it{p}_{T}(#Lambda_{c}^{#plus}) (GeV/#it{c});#Lambda_{c}^{+} candidate decay length (cm);#Lambda_{c}^{+} candidate cos(#vartheta_{P})", {HistType::kTHnSparseF, {axisMassLb, axisPtLb, axisDecayLength, axisNormDecayLength, axisImpParProd, axisCosp, axisMassLc, axisPtLc, axisDecayLength, axisCosp}});
          if (fillBackground) {
            registry.add("hMassPtCutVarsRecBg", "#Lambda_{b}^{0} candidates (unmatched);#it{M} (D#pi) (GeV/#it{c}^{2});#it{p}_{T}(#Lambda_{b}^{0}) (GeV/#it{c});#Lambda_{b}^{0} candidate decay length (cm);#Lambda_{b}^{0} candidate norm. decay length XY (cm);#Lambda_{b}^{0} candidate impact parameter product (cm);#Lambda_{b}^{0} candidate cos(#vartheta_{P});#it{M} (pK#pi) (GeV/#it{c}^{2});#it{p}_{T}(#Lambda_{c}^{#plus}) (GeV/#it{c});#Lambda_{c}^{+} candidate decay length (cm);#Lambda_{c}^{+} candidate cos(#vartheta_{P})", {HistType::kTHnSparseF, {axisMassLb, axisPtLb, axisDecayLength, axisNormDecayLength, axisImpParProd, axisCosp, axisMassLc, axisPtLc, axisDecayLength, axisCosp}});
          }
        } else {
          registry.add("hMassPtCutVarsRecSig", "#Lambda_{b}^{0} candidates (matched);#it{M} (D#pi) (GeV/#it{c}^{2});#it{p}_{T}(#Lambda_{b}^{0}) (GeV/#it{c});#Lambda_{b}^{0} candidate decay length (cm);#Lambda_{b}^{0} candidate norm. decay length XY (cm);#Lambda_{b}^{0} candidate impact parameter product (cm);#Lambda_{b}^{0} candidate cos(#vartheta_{P});#it{M} (pK#pi) (GeV/#it{c}^{2});#it{p}_{T}(#Lambda_{c}^{#plus}) (GeV/#it{c});#Lambda_{c}^{+} candidate ML score bkg;#Lambda_{c}^{+} candidate ML score nonprompt", {HistType::kTHnSparseF, {axisMassLb, axisPtLb, axisDecayLength, axisNormDecayLength, axisImpParProd, axisCosp, axisMassLc, axisPtLc, axisMlScore, axisMlScore}});
          if (fillBackground) {
            registry.add("hMassPtCutVarsRecBg", "#Lambda_{b}^{0} candidates (unmatched);#it{M} (D#pi) (GeV/#it{c}^{2});#it{p}_{T}(#Lambda_{b}^{0}) (GeV/#it{c});#Lambda_{b}^{0} candidate decay length (cm);#Lambda_{b}^{0} candidate norm. decay length XY (cm);#Lambda_{b}^{0} candidate impact parameter product (cm);#Lambda_{b}^{0} candidate cos(#vartheta_{P});#it{M} (pK#pi) (GeV/#it{c}^{2});#it{p}_{T}(#Lambda_{c}^{#plus}) (GeV/#it{c});#Lambda_{c}^{+} candidate ML score bkg;#Lambda_{c}^{+} candidate ML score nonprompt", {HistType::kTHnSparseF, {axisMassLb, axisPtLb, axisDecayLength, axisNormDecayLength, axisImpParProd, axisCosp, axisMassLc, axisPtLc, axisMlScore, axisMlScore}});
          }
        }
      }
    }
  }

  /// Selection of Lb daughter in geometrical acceptance
  /// \param etaProng is the pseudorapidity of Lb prong
  /// \param ptProng is the pT of Lb prong
  /// \return true if prong is in geometrical acceptance
  template <typename T = float>
  bool isProngInAcceptance(const T& etaProng, const T& ptProng)
  {
    return std::abs(etaProng) <= etaTrackMax && ptProng >= ptTrackMin;
  }

  /// Fill candidate information at reconstruction level
  /// \param doMc is the flag to enable the filling with MC information
  /// \param withDecayTypeCheck is the flag to enable MC with decay type check
  /// \param withLcMl is the flag to enable the filling with ML scores for the Lc daughter
  /// \param withLbMl is the flag to enable the filling with ML scores for the Lb candidate
  /// \param candidate is the Lb candidate
  /// \param candidatesLc is the table with Lc candidates
  template <bool DoMc, bool WithDecayTypeCheck, bool WithLcMl, bool WithLbMl, typename Cand, typename CandsLc>
  void fillCand(Cand const& candidate,
                CandsLc const&)
  {
    auto ptCandLb = candidate.pt();
    auto invMassLb = HfHelper::invMassLbToLcPi(candidate);
    auto candLc = candidate.template prong0_as<CandsLc>();
    auto ptLc = candidate.ptProng0();
    auto invMassLc = candLc.invMassHypo0() > 0 ? candLc.invMassHypo0() : candLc.invMassHypo1();
    // TODO: here we are assuming that only one of the two hypotheses is filled, to be checked
    std::array<float, 3> const posPv{candidate.posX(), candidate.posY(), candidate.posZ()};
    std::array<float, 3> const posSvLc{candLc.xSecondaryVertex(), candLc.ySecondaryVertex(), candLc.zSecondaryVertex()};
    std::array<float, 3> const momLc{candLc.pVector()};
    auto cospLc = RecoDecay::cpa(posPv, posSvLc, momLc);
    auto cospXyLc = RecoDecay::cpaXY(posPv, posSvLc, momLc);
    auto decLenLc = RecoDecay::distance(posPv, posSvLc);
    auto decLenXyLc = RecoDecay::distanceXY(posPv, posSvLc);

    int8_t flagMcMatchRec = 0;
    int8_t flagWrongCollision = 0;
    bool isSignal = false;
    if constexpr (DoMc) {
      flagMcMatchRec = candidate.flagMcMatchRec();
      flagWrongCollision = candidate.flagWrongCollision();
      isSignal = TESTBIT(std::abs(flagMcMatchRec), hf_cand_lb::DecayTypeMc::LbToLcPiToPKPiPi);
    }

    if (fillHistograms) {
      if constexpr (DoMc) {
        if (isSignal) {
          registry.fill(HIST("hMassRecSig"), ptCandLb, HfHelper::invMassLbToLcPi(candidate));
          registry.fill(HIST("hPtProng0RecSig"), ptCandLb, candidate.ptProng0());
          registry.fill(HIST("hPtProng1RecSig"), ptCandLb, candidate.ptProng1());
          registry.fill(HIST("hImpParProdRecSig"), ptCandLb, candidate.impactParameterProduct());
          registry.fill(HIST("hDecLengthRecSig"), ptCandLb, candidate.decayLength());
          registry.fill(HIST("hDecLengthXyRecSig"), ptCandLb, candidate.decayLengthXY());
          registry.fill(HIST("hNormDecLengthXyRecSig"), ptCandLb, candidate.decayLengthXY() / candidate.errorDecayLengthXY());
          registry.fill(HIST("hDcaProng0RecSig"), ptCandLb, candidate.impactParameter0());
          registry.fill(HIST("hDcaProng1RecSig"), ptCandLb, candidate.impactParameter1());
          registry.fill(HIST("hCospRecSig"), ptCandLb, candidate.cpa());
          registry.fill(HIST("hCospXyRecSig"), ptCandLb, candidate.cpaXY());
          registry.fill(HIST("hEtaRecSig"), ptCandLb, candidate.eta());
          registry.fill(HIST("hRapidityRecSig"), ptCandLb, HfHelper::yLb(candidate));
          registry.fill(HIST("hinvMassLcRecSig"), ptLc, invMassLc);
          registry.fill(HIST("hDecLengthLcRecSig"), ptLc, decLenLc);
          registry.fill(HIST("hDecLengthXyLcRecSig"), ptLc, decLenXyLc);
          registry.fill(HIST("hCospLcRecSig"), ptLc, cospLc);
          registry.fill(HIST("hCospXyLcRecSig"), ptLc, cospXyLc);
          if constexpr (WithDecayTypeCheck) {
            registry.fill(HIST("hDecayTypeMc"), 1 + hf_cand_lb::DecayTypeMc::LbToLcPiToPKPiPi, invMassLb, ptCandLb);
          }
          if constexpr (WithLcMl) {
            registry.fill(HIST("hMlScoreBkgLcRecSig"), ptLc, candidate.prong0MlScoreBkg());
            registry.fill(HIST("hMlScorePromptLcRecSig"), ptLc, candidate.prong0MlScorePrompt());
            registry.fill(HIST("hMlScoreNonPromptLcRecSig"), ptLc, candidate.prong0MlScoreNonprompt());
          }
          if constexpr (WithLbMl) {
            registry.fill(HIST("hMlScoreSigLbRecSig"), ptCandLb, candidate.mlProbLbToLcPi());
          }
        } else if (fillBackground) {
          registry.fill(HIST("hMassRecBg"), ptCandLb, HfHelper::invMassLbToLcPi(candidate));
          registry.fill(HIST("hPtProng0RecBg"), ptCandLb, candidate.ptProng0());
          registry.fill(HIST("hPtProng1RecBg"), ptCandLb, candidate.ptProng1());
          registry.fill(HIST("hImpParProdRecBg"), ptCandLb, candidate.impactParameterProduct());
          registry.fill(HIST("hDecLengthRecBg"), ptCandLb, candidate.decayLength());
          registry.fill(HIST("hDecLengthXyRecBg"), ptCandLb, candidate.decayLengthXY());
          registry.fill(HIST("hNormDecLengthXyRecBg"), ptCandLb, candidate.decayLengthXY() / candidate.errorDecayLengthXY());
          registry.fill(HIST("hDcaProng0RecBg"), ptCandLb, candidate.impactParameter0());
          registry.fill(HIST("hDcaProng1RecBg"), ptCandLb, candidate.impactParameter1());
          registry.fill(HIST("hCospRecBg"), ptCandLb, candidate.cpa());
          registry.fill(HIST("hCospXyRecBg"), ptCandLb, candidate.cpaXY());
          registry.fill(HIST("hEtaRecBg"), ptCandLb, candidate.eta());
          registry.fill(HIST("hRapidityRecBg"), ptCandLb, HfHelper::yLb(candidate));
          registry.fill(HIST("hinvMassLcRecBg"), ptLc, invMassLc);
          registry.fill(HIST("hDecLengthLcRecBg"), ptLc, decLenLc);
          registry.fill(HIST("hDecLengthXyLcRecBg"), ptLc, decLenXyLc);
          registry.fill(HIST("hCospLcRecBg"), ptLc, cospLc);
          registry.fill(HIST("hCospXyLcRecBg"), ptLc, cospXyLc);
          if constexpr (WithLcMl) {
            registry.fill(HIST("hMlScoreBkgLcRecBg"), ptLc, candidate.prong0MlScoreBkg());
            registry.fill(HIST("hMlScorePromptLcRecBg"), ptLc, candidate.prong0MlScorePrompt());
            registry.fill(HIST("hMlScoreNonPromptLcRecBg"), ptLc, candidate.prong0MlScoreNonprompt());
          }
          if constexpr (WithLbMl) {
            registry.fill(HIST("hMlScoreSigLbRecBg"), ptCandLb, candidate.mlProbLbToLcPi());
          }
        } else if constexpr (WithDecayTypeCheck) {
          if (TESTBIT(flagMcMatchRec, hf_cand_lb::DecayTypeMc::LbToLcKToPKPiK)) { // Lb → Lc+ K- → (pK-π+) K-
            registry.fill(HIST("hDecayTypeMc"), 1 + hf_cand_lb::DecayTypeMc::LbToLcKToPKPiK, invMassLb, ptCandLb);
          } else if (TESTBIT(flagMcMatchRec, hf_cand_lb::DecayTypeMc::B0ToDplusPiToPiKPiPi)) { // // B0 → D- π+ → (π- K+ π-) π+
            registry.fill(HIST("hDecayTypeMc"), 1 + hf_cand_lb::DecayTypeMc::B0ToDplusPiToPiKPiPi, invMassLb, ptCandLb);
          } else if (TESTBIT(flagMcMatchRec, hf_cand_lb::DecayTypeMc::PartlyRecoDecay)) { // Partly reconstructed decay channel
            registry.fill(HIST("hDecayTypeMc"), 1 + hf_cand_lb::DecayTypeMc::PartlyRecoDecay, invMassLb, ptCandLb);
          } else {
            registry.fill(HIST("hDecayTypeMc"), 1 + hf_cand_lb::DecayTypeMc::OtherDecay, invMassLb, ptCandLb);
          }
        }
      } else {
        registry.fill(HIST("hMass"), ptCandLb, invMassLb);
        registry.fill(HIST("hPtProng0"), ptCandLb, candidate.ptProng0());
        registry.fill(HIST("hPtProng1"), ptCandLb, candidate.ptProng1());
        registry.fill(HIST("hImpParProd"), ptCandLb, candidate.impactParameterProduct());
        registry.fill(HIST("hDecLength"), ptCandLb, candidate.decayLength());
        registry.fill(HIST("hDecLengthXy"), ptCandLb, candidate.decayLengthXY());
        registry.fill(HIST("hNormDecLengthXy"), ptCandLb, candidate.decayLengthXY() / candidate.errorDecayLengthXY());
        registry.fill(HIST("hDcaProng0"), ptCandLb, candidate.impactParameter0());
        registry.fill(HIST("hDcaProng1"), ptCandLb, candidate.impactParameter1());
        registry.fill(HIST("hCosp"), ptCandLb, candidate.cpa());
        registry.fill(HIST("hCospXy"), ptCandLb, candidate.cpaXY());
        registry.fill(HIST("hEta"), ptCandLb, candidate.eta());
        registry.fill(HIST("hRapidity"), ptCandLb, HfHelper::yLb(candidate));
        registry.fill(HIST("hinvMassLc"), ptLc, invMassLc);
        registry.fill(HIST("hDecLengthLc"), ptLc, decLenLc);
        registry.fill(HIST("hDecLengthXyLc"), ptLc, decLenXyLc);
        registry.fill(HIST("hCospLc"), ptLc, cospLc);
        registry.fill(HIST("hCospXyLc"), ptLc, cospXyLc);

        if constexpr (WithLcMl) {
          registry.fill(HIST("hMlScoreBkgLc"), ptLc, candidate.prong0MlScoreBkg());
          registry.fill(HIST("hMlScorePromptLc"), ptLc, candidate.prong0MlScorePrompt());
          registry.fill(HIST("hMlScoreNonPromptLc"), ptLc, candidate.prong0MlScoreNonprompt());
        }
        if constexpr (WithLbMl) {
          registry.fill(HIST("hMlScoreSigLb"), ptCandLb, candidate.mlProbLbToLcPi());
        }
      }
    }
    if (fillSparses) {
      if constexpr (WithLcMl) {
        if (isSignal) {
          if constexpr (WithLcMl) {
            registry.fill(HIST("hMassPtCutVarsRecSig"), invMassLb, ptCandLb, candidate.decayLength(), candidate.decayLengthXY() / candidate.errorDecayLengthXY(), candidate.impactParameterProduct(), candidate.cpa(), invMassLc, ptLc, candidate.prong0MlScoreBkg(), candidate.prong0MlScoreNonprompt());
          } else {
            registry.fill(HIST("hMassPtCutVarsRecSig"), invMassLb, ptCandLb, candidate.decayLength(), candidate.decayLengthXY() / candidate.errorDecayLengthXY(), candidate.impactParameterProduct(), candidate.cpa(), invMassLc, ptLc, decLenLc, cospLc);
          }
        } else if (fillBackground) {
          if constexpr (WithLcMl) {
            registry.fill(HIST("hMassPtCutVarsRecBg"), invMassLb, ptCandLb, candidate.decayLength(), candidate.decayLengthXY() / candidate.errorDecayLengthXY(), candidate.impactParameterProduct(), candidate.cpa(), invMassLc, ptLc, candidate.prong0MlScoreBkg(), candidate.prong0MlScoreNonprompt());
          } else {
            registry.fill(HIST("hMassPtCutVarsRecBg"), invMassLb, ptCandLb, candidate.decayLength(), candidate.decayLengthXY() / candidate.errorDecayLengthXY(), candidate.impactParameterProduct(), candidate.cpa(), invMassLc, ptLc, decLenLc, cospLc);
          }
        }
      } else {
        if constexpr (WithLcMl) {
          registry.fill(HIST("hMassPtCutVars"), invMassLb, ptCandLb, candidate.decayLength(), candidate.decayLengthXY() / candidate.errorDecayLengthXY(), candidate.impactParameterProduct(), candidate.cpa(), invMassLc, ptLc, candidate.prong0MlScoreBkg(), candidate.prong0MlScoreNonprompt());
        } else {
          registry.fill(HIST("hMassPtCutVars"), invMassLb, ptCandLb, candidate.decayLength(), candidate.decayLengthXY() / candidate.errorDecayLengthXY(), candidate.impactParameterProduct(), candidate.cpa(), invMassLc, ptLc, decLenLc, cospLc);
        }
      }
    }
    if (fillTree) {
      float const pseudoRndm = ptLc * 1000. - static_cast<int64_t>(ptLc * 1000);
      if (flagMcMatchRec != 0 || (((DoMc && fillBackground) || !DoMc) && (ptCandLb >= ptMaxForDownSample || pseudoRndm < downSampleBkgFactor))) {
        float prong0MlScoreBkg = -1.;
        float prong0MlScorePrompt = -1.;
        float prong0MlScoreNonprompt = -1.;
        float candidateMlScoreSig = -1;
        if constexpr (WithLcMl) {
          prong0MlScoreBkg = candidate.prong0MlScoreBkg();
          prong0MlScorePrompt = candidate.prong0MlScorePrompt();
          prong0MlScoreNonprompt = candidate.prong0MlScoreNonprompt();
        }
        if constexpr (WithLbMl) {
          candidateMlScoreSig = candidate.mlProbLbToLcPi();
        }
        auto prong1 = candidate.template prong1_as<TracksPion>();

        float ptMother = -1.;
        if constexpr (DoMc) {
          ptMother = candidate.ptMother();
        }

        hfRedCandLbLite(
          // Lb features
          invMassLb,
          ptCandLb,
          candidate.eta(),
          candidate.phi(),
          HfHelper::yLb(candidate),
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
          candidate.isSelLbToLcPi(),
          // Lc-baryon features
          invMassLc,
          ptLc,
          decLenLc,
          decLenXyLc,
          candidate.impactParameter0(),
          candLc.ptProngMin(),
          candLc.absEtaProngMin(),
          candLc.itsNClsProngMin(),
          candLc.tpcNClsCrossedRowsProngMin(),
          candLc.tpcChi2NClProngMax(),
          candLc.tpcNSigmaPrProng0(),
          candLc.tofNSigmaPrProng0(),
          candLc.tpcTofNSigmaPrProng0(),
          candLc.tpcNSigmaKaProng1(),
          candLc.tofNSigmaKaProng1(),
          candLc.tpcTofNSigmaKaProng1(),
          candLc.tpcNSigmaPiProng2(),
          candLc.tofNSigmaPiProng2(),
          candLc.tpcTofNSigmaPiProng2(),
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

        if constexpr (WithDecayTypeCheck) {
          hfRedLbMcCheck(
            flagMcMatchRec,
            flagWrongCollision,
            invMassLc,
            ptLc,
            invMassLb,
            ptCandLb,
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
  void fillCandMcGen(aod::HfMcGenRedLbs::iterator const& particle)
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
    bool const prongsInAcc = isProngInAcceptance(etaProngs[0], ptProngs[0]) && isProngInAcceptance(etaProngs[1], ptProngs[1]);

    if (fillHistograms) {
      registry.fill(HIST("hPtProng0Gen"), ptParticle, ptProngs[0]);
      registry.fill(HIST("hPtProng1Gen"), ptParticle, ptProngs[1]);
      registry.fill(HIST("hYProng0Gen"), ptParticle, yProngs[0]);
      registry.fill(HIST("hYProng1Gen"), ptParticle, yProngs[1]);
      registry.fill(HIST("hEtaProng0Gen"), ptParticle, etaProngs[0]);
      registry.fill(HIST("hEtaProng1Gen"), ptParticle, etaProngs[1]);

      registry.fill(HIST("hYGen"), ptParticle, yParticle);
      registry.fill(HIST("hEtaGen"), ptParticle, etaParticle);

      // generated Lb with daughters in geometrical acceptance
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
  void processData(soa::Filtered<soa::Join<aod::HfRedCandLb, aod::HfSelLbToLcPi>> const& candidates,
                   CandsLc const& candidatesLc,
                   TracksPion const&)
  {
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(HfHelper::yLb(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCand<false, false, false, false>(candidate, candidatesLc);
    } // candidate loop
  } // processData
  PROCESS_SWITCH(HfTaskLbReduced, processData, "Process data without ML scores for Lb and Lc daughter", true);

  void processDataWithLcMl(soa::Filtered<soa::Join<aod::HfRedCandLb, aod::HfRedLbLcMls, aod::HfSelLbToLcPi>> const& candidates,
                           CandsLc const& candidatesLc,
                           TracksPion const&)
  {
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(HfHelper::yLb(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCand<false, false, true, false>(candidate, candidatesLc);
    } // candidate loop
  } // processDataWithLcMl
  PROCESS_SWITCH(HfTaskLbReduced, processDataWithLcMl, "Process data with(out) ML scores for Lc daughter (Lb)", false);

  void processDataWithLbMl(soa::Filtered<soa::Join<aod::HfRedCandLb, aod::HfMlLbToLcPi, aod::HfSelLbToLcPi>> const& candidates,
                           CandsLc const& candidatesLc,
                           TracksPion const&)
  {
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(HfHelper::yLb(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCand<false, false, false, true>(candidate, candidatesLc);
    } // candidate loop
  } // processDataWithLbMl
  PROCESS_SWITCH(HfTaskLbReduced, processDataWithLbMl, "Process data with(out) ML scores for Lb (Lc daughter)", false);

  void processMc(soa::Filtered<soa::Join<aod::HfRedCandLb, aod::HfSelLbToLcPi, aod::HfMcRecRedLbs>> const& candidates,
                 aod::HfMcGenRedLbs const& mcParticles,
                 CandsLc const& candidatesLc,
                 TracksPion const&)
  {
    // MC rec
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(HfHelper::yLb(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCand<true, false, false, false>(candidate, candidatesLc);
    } // rec

    // MC gen. level
    for (const auto& particle : mcParticles) {
      fillCandMcGen(particle);
    } // gen
  } // processMc
  PROCESS_SWITCH(HfTaskLbReduced, processMc, "Process MC without ML scores for Lb and Lc daughter", false);

  void processMcWithDecayTypeCheck(soa::Filtered<soa::Join<aod::HfRedCandLb, aod::HfSelLbToLcPi, aod::HfMcRecRedLbs, aod::HfMcCheckLbs>> const& candidates,
                                   aod::HfMcGenRedLbs const& mcParticles,
                                   CandsLc const& candidatesLc,
                                   TracksPion const&)
  {
    // MC rec
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(HfHelper::yLb(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCand<true, true, false, false>(candidate, candidatesLc);
    } // rec

    // MC gen. level
    for (const auto& particle : mcParticles) {
      fillCandMcGen(particle);
    } // gen
  } // processMc
  PROCESS_SWITCH(HfTaskLbReduced, processMcWithDecayTypeCheck, "Process MC with decay type check and without ML scores for Lb and D daughter", false);

  void processMcWithLcMl(soa::Filtered<soa::Join<aod::HfRedCandLb, aod::HfRedLbLcMls, aod::HfSelLbToLcPi, aod::HfMcRecRedLbs>> const& candidates,
                         aod::HfMcGenRedLbs const& mcParticles,
                         CandsLc const& candidatesLc,
                         TracksPion const&)
  {
    // MC rec
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(HfHelper::yLb(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCand<true, false, true, false>(candidate, candidatesLc);
    } // rec

    // MC gen. level
    for (const auto& particle : mcParticles) {
      fillCandMcGen(particle);
    } // gen
  } // processMcWithLcMl
  PROCESS_SWITCH(HfTaskLbReduced, processMcWithLcMl, "Process MC with(out) ML scores for Lc daughter (Lb)", false);

  void processMcWithLcMlAndDecayTypeCheck(soa::Filtered<soa::Join<aod::HfRedCandLb, aod::HfRedLbLcMls, aod::HfSelLbToLcPi, aod::HfMcRecRedLbs, aod::HfMcCheckLbs>> const& candidates,
                                          aod::HfMcGenRedLbs const& mcParticles,
                                          CandsLc const& candidatesLc,
                                          TracksPion const&)
  {
    // MC rec
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(HfHelper::yLb(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCand<true, true, true, false>(candidate, candidatesLc);
    } // rec

    // MC gen. level
    for (const auto& particle : mcParticles) {
      fillCandMcGen(particle);
    } // gen
  } // processMc
  PROCESS_SWITCH(HfTaskLbReduced, processMcWithLcMlAndDecayTypeCheck, "Process MC with decay type check and with(out) ML scores for Lb (Lc daughter)", false);

  void processMcWithLbMl(soa::Filtered<soa::Join<aod::HfRedCandLb, aod::HfMlLbToLcPi, aod::HfSelLbToLcPi, aod::HfMcRecRedLbs>> const& candidates,
                         aod::HfMcGenRedLbs const& mcParticles,
                         CandsLc const& candidatesLc,
                         TracksPion const&)
  {
    // MC rec
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(HfHelper::yLb(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCand<true, false, false, true>(candidate, candidatesLc);
    } // rec

    // MC gen. level
    for (const auto& particle : mcParticles) {
      fillCandMcGen(particle);
    } // gen
  } // processMcWithLbMl
  PROCESS_SWITCH(HfTaskLbReduced, processMcWithLbMl, "Process MC with(out) ML scores for Lb (Lc daughter)", false);

  void processMcWithLbMlAndDecayTypeCheck(soa::Filtered<soa::Join<aod::HfRedCandLb, aod::HfMlLbToLcPi, aod::HfSelLbToLcPi, aod::HfMcRecRedLbs, aod::HfMcCheckLbs>> const& candidates,
                                          aod::HfMcGenRedLbs const& mcParticles,
                                          CandsLc const& candidatesLc,
                                          TracksPion const&)
  {
    // MC rec
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(HfHelper::yLb(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCand<true, true, false, true>(candidate, candidatesLc);
    } // rec

    // MC gen. level
    for (const auto& particle : mcParticles) {
      fillCandMcGen(particle);
    } // gen
  } // processMc
  PROCESS_SWITCH(HfTaskLbReduced, processMcWithLbMlAndDecayTypeCheck, "Process MC with decay type check and with(out) ML scores for Lb (Lc daughter)", false);
}; // struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskLbReduced>(cfgc)};
}
