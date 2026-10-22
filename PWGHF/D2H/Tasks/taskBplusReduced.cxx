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

/// \file taskBplusReduced.cxx
/// \brief B+ → D0bar π+ → (π+ K-) π+ analysis task
///
/// \author Antonio Palasciano <antonio.palasciano@cern.ch>, Università degli Studi di Bari & INFN, Sezione di Bari

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/D2H/DataModel/ReducedDataModel.h"
#include "PWGHF/D2H/Utils/utilsFlow.h"
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
#include <cmath>
#include <cstdint>
#include <numeric>

using namespace o2;
using namespace o2::aod;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::analysis::hf_flow_utils;

namespace o2::aod
{
namespace hf_cand_bplus_lite
{
DECLARE_SOA_COLUMN(PtD, ptD, float);                                                     //! Transverse momentum of D-meson daughter candidate (GeV/c)
DECLARE_SOA_COLUMN(PtBach, ptBach, float);                                               //! Transverse momentum of bachelor pion (GeV/c)
DECLARE_SOA_COLUMN(AbsEtaBach, absEtaBach, float);                                       //! Absolute pseudorapidity of bachelor pion
DECLARE_SOA_COLUMN(ItsNClsBach, itsNClsBach, int);                                       //! Number of ITS clusters of bachelor pion
DECLARE_SOA_COLUMN(TpcNClsCrossedRowsBach, tpcNClsCrossedRowsBach, int);                 //! Number of TPC crossed rows of prongs of bachelor pion
DECLARE_SOA_COLUMN(TpcChi2NClBach, tpcChi2NClBach, float);                               //! Maximum TPC chi2 of prongs of D0-meson daughter candidate
DECLARE_SOA_COLUMN(PtDmesProngMin, ptDmesProngMin, float);                               //! Minimum pT of prongs of D-meson daughter candidate (GeV/c)
DECLARE_SOA_COLUMN(AbsEtaDmesProngMin, absEtaDmesProngMin, float);                       //! Minimum absolute pseudorapidity of prongs of D-meson daughter candidate
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
DECLARE_SOA_COLUMN(DecayLength, decayLength, float);                                     //! Decay length of candidate (cm)
DECLARE_SOA_COLUMN(DecayLengthXY, decayLengthXY, float);                                 //! Transverse decay length of candidate (cm)
DECLARE_SOA_COLUMN(DecayLengthNormalised, decayLengthNormalised, float);                 //! Normalised decay length of candidate
DECLARE_SOA_COLUMN(DecayLengthXYNormalised, decayLengthXYNormalised, float);             //! Normalised transverse decay length of candidate
DECLARE_SOA_COLUMN(DecayLengthD, decayLengthD, float);                                   //! Decay length of D-meson daughter candidate (cm)
DECLARE_SOA_COLUMN(DecayLengthXYD, decayLengthXYD, float);                               //! Transverse decay length of D-meson daughter candidate (cm)
DECLARE_SOA_COLUMN(DecayLengthDToB, decayLengthDToB, float);                             //! Decay length of D-meson daughter candidate from B decay vertex (cm)
DECLARE_SOA_COLUMN(DecayLengthXYDToB, decayLengthXYDToB, float);                         //! Transverse decay length of D-meson daughter candidate from B decay vertex (cm)
DECLARE_SOA_COLUMN(ImpactParameterD, impactParameterD, float);                           //! Impact parameter product of D-meson daughter candidate
DECLARE_SOA_COLUMN(ImpactParameterBach, impactParameterBach, float);                     //! Impact parameter product of bachelor pion
DECLARE_SOA_COLUMN(ImpactParameterProduct, impactParameterProduct, float);               //! Impact parameter product of daughters
DECLARE_SOA_COLUMN(Cpa, cpa, float);                                                     //! Cosine pointing angle of candidate
DECLARE_SOA_COLUMN(CpaXY, cpaXY, float);                                                 //! Cosine pointing angle of candidate in transverse plane
DECLARE_SOA_COLUMN(CpaD, cpaD, float);                                                   //! Cosine pointing angle of D candidate to PV
DECLARE_SOA_COLUMN(CpaDXY, cpaDXY, float);                                               //! Cosine pointing angle of D candidate to PV in transverse plane
DECLARE_SOA_COLUMN(CpaDToB, cpaDToB, float);                                             //! Cosine pointing angle of D candidate to B decay vertex
DECLARE_SOA_COLUMN(CpaDToBXY, cpaDToBXY, float);                                         //! Cosine pointing angle of D candidate to B decay vertex in transverse plane
DECLARE_SOA_COLUMN(ScalarProd, scalarProd, float);                                       //! Scalar product of B candidate
DECLARE_SOA_COLUMN(Centrality, centrality, float);                                       //! Centrality of the collision
DECLARE_SOA_COLUMN(MaxNormalisedDeltaIP, maxNormalisedDeltaIP, float);                   //! Maximum normalized difference between measured and expected impact parameter of candidate prongs
DECLARE_SOA_COLUMN(MlScoreSig, mlScoreSig, float);                                       //! ML score for signal class
DECLARE_SOA_COLUMN(FlagWrongCollision, flagWrongCollision, int8_t);                      //! Flag for association with wrong collision
} // namespace hf_cand_bplus_lite

DECLARE_SOA_TABLE(HfRedCandBpLites, "AOD", "HFREDCANDBPLITE", //! Table with some B+ properties
                                                              // B meson features
                  hf_cand_bplus_lite::M,
                  hf_cand_bplus_lite::Pt,
                  hf_cand_bplus_lite::Eta,
                  hf_cand_bplus_lite::Phi,
                  hf_cand_bplus_lite::Y,
                  hf_cand_bplus_lite::Cpa,
                  hf_cand_bplus_lite::CpaXY,
                  hf_cand::Chi2PCA,
                  hf_cand_bplus_lite::DecayLength,
                  hf_cand_bplus_lite::DecayLengthXY,
                  hf_cand_bplus_lite::DecayLengthNormalised,
                  hf_cand_bplus_lite::DecayLengthXYNormalised,
                  hf_cand_bplus_lite::ImpactParameterProduct,
                  hf_cand_bplus_lite::MaxNormalisedDeltaIP,
                  hf_cand_bplus_lite::MlScoreSig,
                  hf_sel_candidate_bplus::IsSelBplusToD0Pi,
                  hf_cand_bplus_lite::ScalarProd,
                  hf_cand_bplus_lite::Centrality,
                  // D meson features
                  hf_cand_bplus_lite::MD,
                  hf_cand_bplus_lite::PtD,
                  hf_cand_bplus_lite::DecayLengthD,
                  hf_cand_bplus_lite::DecayLengthXYD,
                  hf_cand_bplus_lite::DecayLengthDToB,
                  hf_cand_bplus_lite::DecayLengthXYDToB,
                  hf_cand_bplus_lite::CpaD,
                  hf_cand_bplus_lite::CpaDXY,
                  hf_cand_bplus_lite::CpaDToB,
                  hf_cand_bplus_lite::CpaDToBXY,
                  hf_cand_bplus_lite::ImpactParameterD,
                  hf_cand_bplus_lite::PtDmesProngMin,
                  hf_cand_bplus_lite::AbsEtaDmesProngMin,
                  hf_cand_bplus_lite::ItsNClsDmesProngMin,
                  hf_cand_bplus_lite::TpcNClsCrossedRowsDmesProngMin,
                  hf_cand_bplus_lite::TpcChi2NClDmesProngMax,
                  hf_cand_bplus_lite::NSigTpcPiDmesProng0,
                  hf_cand_bplus_lite::NSigTofPiDmesProng0,
                  hf_cand_bplus_lite::NSigTpcTofPiDmesProng0,
                  hf_cand_bplus_lite::NSigTpcKaDmesProng1,
                  hf_cand_bplus_lite::NSigTofKaDmesProng1,
                  hf_cand_bplus_lite::NSigTpcTofKaDmesProng1,
                  hf_cand_bplus_reduced::Prong0MlScoreBkg,
                  hf_cand_bplus_reduced::Prong0MlScorePrompt,
                  hf_cand_bplus_reduced::Prong0MlScoreNonprompt,
                  // pion features
                  hf_cand_bplus_lite::PtBach,
                  hf_cand_bplus_lite::AbsEtaBach,
                  hf_cand_bplus_lite::ItsNClsBach,
                  hf_cand_bplus_lite::TpcNClsCrossedRowsBach,
                  hf_cand_bplus_lite::TpcChi2NClBach,
                  hf_cand_bplus_lite::ImpactParameterBach,
                  hf_cand_bplus_lite::NSigTpcPiBachelor,
                  hf_cand_bplus_lite::NSigTofPiBachelor,
                  hf_cand_bplus_lite::NSigTpcTofPiBachelor,
                  // MC truth
                  hf_cand_mc_flag::FlagMcMatchRec,
                  hf_cand_mc_flag::OriginMcRec,
                  hf_cand_bplus_lite::FlagWrongCollision,
                  hf_cand_bplus_lite::PtGen);

DECLARE_SOA_TABLE(HfRedBpMcCheck, "AOD", "HFREDBPMCCHECK", //! Table with MC decay type check
                  hf_cand_mc_flag::FlagMcMatchRec,
                  hf_cand_bplus_lite::FlagWrongCollision,
                  hf_cand_bplus_lite::MD,
                  hf_cand_bplus_lite::PtD,
                  hf_cand_bplus_lite::M,
                  hf_cand_bplus_lite::Pt,
                  hf_cand_bplus_lite::MlScoreSig,
                  hf_bplus_mc::PdgCodeBeautyMother,
                  hf_bplus_mc::PdgCodeCharmMother,
                  hf_bplus_mc::PdgCodeProng0,
                  hf_bplus_mc::PdgCodeProng1,
                  hf_bplus_mc::PdgCodeProng2);
} // namespace o2::aod

// string definitions, used for histogram axis labels
const TString stringPt = "#it{p}_{T}(B^{+}) (GeV/#it{c});";
const TString stringPtGen = "#it{p}_{T}^{gen}(B^{+}) (GeV/#it{c});";
const TString stringPtD = "#it{p}_{T}(D^{0}) (GeV/#it{c});";
const TString bPlusCandTitle = "B^{+} candidates;";
const TString entries = "entries";
const TString bPlusCandMatch = "B^{+} candidates (matched);";
const TString bPlusCandUnmatch = "B^{+} candidates (unmatched);";
const TString bPlusPartGen = "B^{+} particles (generated);";

/// B+ analysis task
struct HfTaskBplusReduced {
  Produces<aod::HfRedCandBpLites> hfRedCandBpLite;
  Produces<aod::HfRedBpMcCheck> hfRedBpMcCheck;

  Configurable<int> selectionFlagBplus{"selectionFlagBplus", 1, "Selection Flag for Bplus"};
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
  Configurable<int> centralityEstimator{"centralityEstimator", 1, "Centrality estimator (FT0M: 0, FT0C: 1)"};
  Configurable<int> harmonic{"harmonic", 2, "harmonic number"};
  Configurable<int> qVecDetector{"qVecDetector", 3, "Detector for Q vector estimation (FT0M: 1, FT0A: 2, FT0C: 3, TPC Pos: 4, TPC Neg: 5, TPC Tot: 6)"};

  using CollisionsWithCents = soa::Join<aod::HfRedCollisions, aod::HfRedCollCents>;
  using CollisionsWithCentsQvecs = soa::Join<aod::HfRedCollisions, aod::HfRedCollCents, aod::HfRedQvectors>;
  using TracksPion = soa::Join<HfRedTracks, HfRedTracksPid>;
  using CandsD0 = soa::Join<HfRed2Prongs, HfRedPidDau0s, HfRedPidDau1s>;

  Filter filterSelectCandidates = (aod::hf_sel_candidate_bplus::isSelBplusToD0Pi >= selectionFlagBplus);

  HistogramRegistry registry{"registry"};

  void init(InitContext&)
  {
    std::array<bool, 5> processFuncData{doprocessData, doprocessDataWithDmesMl, doprocessDataWithDmesMlScalarProd, doprocessDataWithBplusMl, doprocessDataWithBplusMlScalarProd};
    if ((std::accumulate(processFuncData.begin(), processFuncData.end(), 0)) > 1) {
      LOGP(fatal, "Only one process function for data can be enabled at a time.");
    }
    std::array<bool, 6> processFuncMc{doprocessMc, doprocessMcWithDecayTypeCheck, doprocessMcWithDmesMl, doprocessMcWithDmesMlAndDecayTypeCheck, doprocessMcWithBplusMl, doprocessMcWithBplusMlAndDecayTypeCheck};
    if ((std::accumulate(processFuncMc.begin(), processFuncMc.end(), 0)) > 1) {
      LOGP(fatal, "Only one process function for MC can be enabled at a time.");
    }

    const AxisSpec axisMlScore{100, 0.f, 1.f};
    const AxisSpec axisMassBplus{300, 4.5f, 6.0f};
    const AxisSpec axisMassD0{300, 1.75f, 2.05f};
    const AxisSpec axisDecayLength{200, 0.f, 0.4f};
    const AxisSpec axisNormDecayLength{100, 0.f, 50.f};
    const AxisSpec axisDca{100, -0.05f, 0.05f};
    const AxisSpec axisCosp{110, 0.f, 1.1f};
    const AxisSpec axisEta{30, -1.5f, 1.5f};
    const AxisSpec axisImpParProd{100, -1.e-3, 1.e-3};
    const AxisSpec axisPtBplus{100, 0.f, 50.f};
    const AxisSpec axisPtD0{100, 0.f, 50.f};
    const AxisSpec axisPtPi{100, 0.f, 10.f};

    if (doprocessData || doprocessDataWithDmesMl || doprocessDataWithDmesMlScalarProd || doprocessDataWithBplusMl || doprocessDataWithBplusMlScalarProd) {
      if (fillHistograms) {
        registry.add("hMass", bPlusCandTitle + stringPt + "#it{M} (D^{0}#pi) (GeV/#it{c}^{2});" + entries, {HistType::kTH2F, {axisPtBplus, axisMassBplus}});
        registry.add("hDecLength", bPlusCandTitle + stringPt + "B^{+} candidate decay length (cm);" + entries, {HistType::kTH2F, {axisPtBplus, axisDecayLength}});
        registry.add("hDecLengthXy", bPlusCandTitle + stringPt + "decay length XY (cm);" + entries, {HistType::kTH2F, {axisPtBplus, axisDecayLength}});
        registry.add("hNormDecLengthXy", bPlusCandTitle + stringPt + "B^{+} candidate norm. decay length XY (cm);" + entries, {HistType::kTH2F, {axisPtBplus, axisNormDecayLength}});
        registry.add("hDcaProng0", bPlusCandTitle + stringPt + "prong 0 (D^{0}) DCAxy to prim. vertex (cm);" + entries, {HistType::kTH2F, {axisPtBplus, axisDca}});
        registry.add("hDcaProng1", bPlusCandTitle + stringPt + "prong 1 (#pi^{#plus}) DCAxy to prim. vertex (cm);" + entries, {HistType::kTH2F, {axisPtBplus, axisDca}});
        registry.add("hCosp", bPlusCandTitle + stringPt + "B^{+} candidate cos(#vartheta_{P});" + entries, {HistType::kTH2F, {axisPtBplus, axisCosp}});
        registry.add("hCospXy", bPlusCandTitle + stringPt + "B^{+} candidate cos(#vartheta_{P}^{XY});" + entries, {HistType::kTH2F, {axisPtBplus, axisCosp}});
        registry.add("hEta", bPlusCandTitle + stringPt + "B^{+} candidate #it{#eta};" + entries, {HistType::kTH2F, {axisPtBplus, axisEta}});
        registry.add("hRapidity", bPlusCandTitle + stringPt + "B^{+} candidate #it{y};" + entries, {HistType::kTH2F, {axisPtBplus, axisEta}});
        registry.add("hImpParProd", bPlusCandTitle + stringPt + "B^{+} candidate impact parameter product;" + entries, {HistType::kTH2F, {axisPtBplus, axisImpParProd}});
        registry.add("hPtProng0", bPlusCandTitle + stringPt + stringPtD + entries, {HistType::kTH2F, {axisPtBplus, axisPtD0}});
        registry.add("hPtProng1", bPlusCandTitle + stringPt + "#it{p}_{T}(#pi^{#plus}) (GeV/#it{c});" + entries, {HistType::kTH2F, {axisPtBplus, axisPtPi}});
        registry.add("hInvMassD", bPlusCandTitle + stringPtD + "prong0, #it{M}(K#pi) (GeV/#it{c}^{2});" + entries, {HistType::kTH2F, {axisPtD0, axisMassD0}});
        registry.add("hDecLengthD", bPlusCandTitle + stringPtD + "D^{0} candidate decay length (cm);" + entries, {HistType::kTH2F, {axisPtD0, axisDecayLength}});
        registry.add("hDecLengthXyD", bPlusCandTitle + stringPtD + "decay length XY (cm);" + entries, {HistType::kTH2F, {axisPtD0, axisDecayLength}});
        registry.add("hCospD", bPlusCandTitle + stringPtD + "D^{0} candidate cos(#vartheta_{P});" + entries, {HistType::kTH2F, {axisPtD0, axisCosp}});
        registry.add("hCospXyD", bPlusCandTitle + stringPtD + "D^{0} candidate cos(#vartheta_{P}^{XY});" + entries, {HistType::kTH2F, {axisPtD0, axisCosp}});

        // ML scores of D0 daughter
        if (doprocessDataWithDmesMl || doprocessDataWithDmesMlScalarProd) {
          registry.add("hMlScoreBkgD", bPlusCandTitle + stringPtD + "prong0, D^{0} ML background score;" + entries, {HistType::kTH2F, {axisPtD0, axisMlScore}});
          registry.add("hMlScorePromptD", bPlusCandTitle + stringPtD + "prong0, D^{0} ML prompt score;" + entries, {HistType::kTH2F, {axisPtD0, axisMlScore}});
          registry.add("hMlScoreNonPromptD", bPlusCandTitle + stringPtD + "prong0, D^{0} ML nonprompt score;" + entries, {HistType::kTH2F, {axisPtD0, axisMlScore}});
        }

        // ML scores of B+ candidate
        if (doprocessDataWithBplusMl || doprocessDataWithBplusMlScalarProd) {
          registry.add("hMlScoreSigBplus", bPlusCandTitle + stringPt + "prong0, B^{+} ML signal score;" + entries, {HistType::kTH2F, {axisPtBplus, axisMlScore}});
        }
      }
      if (fillSparses) {
        if (!(doprocessDataWithDmesMl || doprocessDataWithDmesMlScalarProd)) {
          registry.add("hMassPtCutVars", bPlusCandTitle + "#it{M} (D^{0}#pi) (GeV/#it{c}^{2});" + stringPt + "B^{+} candidate decay length (cm);B^{+} candidate norm. decay length XY (cm);B^{+} candidate impact parameter product (cm);B^{+} candidate cos(#vartheta_{P});#it{M} (K#pi) (GeV/#it{c}^{2});" + stringPtD + "D^{0} candidate decay length (cm);D^{0} candidate cos(#vartheta_{P})", {HistType::kTHnSparseF, {axisMassBplus, axisPtBplus, axisDecayLength, axisNormDecayLength, axisImpParProd, axisCosp, axisMassD0, axisPtD0, axisDecayLength, axisCosp}});
        } else {
          registry.add("hMassPtCutVars", bPlusCandTitle + "#it{M} (D^{0}#pi) (GeV/#it{c}^{2});" + stringPt + "B^{+} candidate decay length (cm);B^{+} candidate norm. decay length XY (cm);B^{+} candidate impact parameter product (cm);B^{+} candidate cos(#vartheta_{P});#it{M} (K#pi) (GeV/#it{c}^{2});" + stringPtD + "D^{0} candidate ML score bkg;D^{0} candidate ML score nonprompt", {HistType::kTHnSparseF, {axisMassBplus, axisPtBplus, axisDecayLength, axisNormDecayLength, axisImpParProd, axisCosp, axisMassD0, axisPtD0, axisMlScore, axisMlScore}});
        }
      }
    }

    if (doprocessMc || doprocessMcWithDecayTypeCheck || doprocessMcWithDmesMl || doprocessMcWithDmesMlAndDecayTypeCheck || doprocessMcWithBplusMl || doprocessMcWithBplusMlAndDecayTypeCheck) {
      if (fillHistograms) {
        // gen histos
        registry.add("hEtaGen", bPlusPartGen + stringPtGen + "#it{#eta}^{gen}(B^{+});" + entries, {HistType::kTH2F, {axisPtBplus, axisEta}});
        registry.add("hYGen", bPlusPartGen + stringPtGen + "#it{y}^{gen}(B^{+});" + entries, {HistType::kTH2F, {axisPtBplus, axisEta}});
        registry.add("hYGenWithProngsInAcceptance", "MC particles (generated-daughters in acceptance);" + stringPtGen + "#it{y}^{gen}(B^{+});" + entries, {HistType::kTH2F, {axisPtBplus, axisEta}});
        registry.add("hPtProng0Gen", bPlusPartGen + stringPtGen + "#it{p}_{T}^{gen}(D^{0}) (GeV/#it{c});" + entries, {HistType::kTH2F, {axisPtBplus, axisPtD0}});
        registry.add("hPtProng1Gen", bPlusPartGen + stringPtGen + "#it{p}_{T}^{gen}(#pi^{#plus}) (GeV/#it{c});" + entries, {HistType::kTH2F, {axisPtBplus, axisPtPi}});
        registry.add("hYProng0Gen", bPlusPartGen + stringPtGen + "#it{y}^{gen}(D^{0});" + entries, {HistType::kTH2F, {axisPtBplus, axisEta}});
        registry.add("hYProng1Gen", bPlusPartGen + stringPtGen + "#it{y}^{gen}(#pi^{#plus});" + entries, {HistType::kTH2F, {axisPtBplus, axisEta}});
        registry.add("hEtaProng0Gen", bPlusPartGen + stringPtGen + "#it{#eta}^{gen}(D^{0});" + entries, {HistType::kTH2F, {axisPtBplus, axisEta}});
        registry.add("hEtaProng1Gen", bPlusPartGen + stringPtGen + "#it{#eta}^{gen}(#pi^{#plus});" + entries, {HistType::kTH2F, {axisPtBplus, axisEta}});

        // reco histos
        // signal
        registry.add("hMassRecSig", bPlusCandMatch + stringPt + "#it{M} (D^{0}#pi) (GeV/#it{c}^{2});" + entries, {HistType::kTH2F, {axisPtBplus, axisMassBplus}});
        registry.add("hDecLengthRecSig", bPlusCandMatch + stringPt + "B^{+} candidate decay length (cm);" + entries, {HistType::kTH2F, {axisPtBplus, axisDecayLength}});
        registry.add("hDecLengthXyRecSig", bPlusCandMatch + stringPt + "decay length XY (cm);" + entries, {HistType::kTH2F, {axisPtBplus, axisDecayLength}});
        registry.add("hNormDecLengthXyRecSig", bPlusCandMatch + stringPt + "B^{+} candidate norm. decay length XY (cm);" + entries, {HistType::kTH2F, {axisPtBplus, axisNormDecayLength}});
        registry.add("hDcaProng0RecSig", bPlusCandMatch + stringPt + "prong 0 (D^{0}) DCAxy to prim. vertex (cm);" + entries, {HistType::kTH2F, {axisPtBplus, axisDca}});
        registry.add("hDcaProng1RecSig", bPlusCandMatch + stringPt + "prong 1 (#pi^{#plus}) DCAxy to prim. vertex (cm);" + entries, {HistType::kTH2F, {axisPtBplus, axisDca}});
        registry.add("hCospRecSig", bPlusCandMatch + stringPt + "B^{+} candidate cos(#vartheta_{P});" + entries, {HistType::kTH2F, {axisPtBplus, axisCosp}});
        registry.add("hCospXyRecSig", bPlusCandMatch + stringPt + "B^{+} candidate cos(#vartheta_{P}^{XY});" + entries, {HistType::kTH2F, {axisPtBplus, axisCosp}});
        registry.add("hEtaRecSig", bPlusCandMatch + stringPt + "B^{+} candidate #it{#eta};" + entries, {HistType::kTH2F, {axisPtBplus, axisEta}});
        registry.add("hRapidityRecSig", bPlusCandMatch + stringPt + "B^{+} candidate #it{y};" + entries, {HistType::kTH2F, {axisPtBplus, axisEta}});
        registry.add("hImpParProdRecSig", bPlusCandMatch + stringPt + "B^{+} candidate impact parameter product;" + entries, {HistType::kTH2F, {axisPtBplus, axisImpParProd}});
        registry.add("hPtProng0RecSig", bPlusCandMatch + stringPt + stringPtD + entries, {HistType::kTH2F, {axisPtBplus, axisPtD0}});
        registry.add("hPtProng1RecSig", bPlusCandMatch + stringPt + "#it{p}_{T}(#pi^{#plus}) (GeV/#it{c});" + entries, {HistType::kTH2F, {axisPtBplus, axisPtPi}});
        registry.add("hInvMassDRecSig", bPlusCandMatch + stringPtD + "prong0, #it{M}(K#pi) (GeV/#it{c}^{2});" + entries, {HistType::kTH2F, {axisPtD0, axisMassD0}});
        registry.add("hDecLengthDRecSig", bPlusCandMatch + stringPtD + "D^{0} candidate decay length (cm);" + entries, {HistType::kTH2F, {axisPtD0, axisDecayLength}});
        registry.add("hDecLengthXyDRecSig", bPlusCandMatch + stringPtD + "decay length XY (cm);" + entries, {HistType::kTH2F, {axisPtD0, axisDecayLength}});
        registry.add("hCospDRecSig", bPlusCandMatch + stringPtD + "D^{0} candidate cos(#vartheta_{P});" + entries, {HistType::kTH2F, {axisPtD0, axisCosp}});
        registry.add("hCospXyDRecSig", bPlusCandMatch + stringPtD + "D^{0} candidate cos(#vartheta_{P}^{XY});" + entries, {HistType::kTH2F, {axisPtD0, axisCosp}});

        // background
        if (fillBackground) {
          registry.add("hMassRecBg", bPlusCandUnmatch + stringPt + "#it{M} (D^{0}#pi) (GeV/#it{c}^{2});" + entries, {HistType::kTH2F, {axisPtBplus, axisMassBplus}});
          registry.add("hDecLengthRecBg", bPlusCandUnmatch + stringPt + "B^{+} candidate decay length (cm);" + entries, {HistType::kTH2F, {axisPtBplus, axisDecayLength}});
          registry.add("hDecLengthXyRecBg", bPlusCandUnmatch + stringPt + "decay length XY (cm);" + entries, {HistType::kTH2F, {axisPtBplus, axisDecayLength}});
          registry.add("hNormDecLengthXyRecBg", bPlusCandUnmatch + stringPt + "B^{+} candidate norm. decay length XY (cm);" + entries, {HistType::kTH2F, {axisPtBplus, axisNormDecayLength}});
          registry.add("hDcaProng0RecBg", bPlusCandUnmatch + stringPt + "prong 0 (D^{0}) DCAxy to prim. vertex (cm);" + entries, {HistType::kTH2F, {axisPtBplus, axisDca}});
          registry.add("hDcaProng1RecBg", bPlusCandUnmatch + stringPt + "prong 1 (#pi^{#plus}) DCAxy to prim. vertex (cm);" + entries, {HistType::kTH2F, {axisPtBplus, axisDca}});
          registry.add("hCospRecBg", bPlusCandUnmatch + stringPt + "B^{+} candidate cos(#vartheta_{P});" + entries, {HistType::kTH2F, {axisPtBplus, axisCosp}});
          registry.add("hCospXyRecBg", bPlusCandUnmatch + stringPt + "B^{+} candidate cos(#vartheta_{P}^{XY});" + entries, {HistType::kTH2F, {axisPtBplus, axisCosp}});
          registry.add("hEtaRecBg", bPlusCandUnmatch + stringPt + "B^{+} candidate #it{#eta};" + entries, {HistType::kTH2F, {axisPtBplus, axisEta}});
          registry.add("hRapidityRecBg", bPlusCandUnmatch + stringPt + "B^{+} candidate #it{y};" + entries, {HistType::kTH2F, {axisPtBplus, axisEta}});
          registry.add("hImpParProdRecBg", bPlusCandUnmatch + stringPt + "B^{+} candidate impact parameter product;" + entries, {HistType::kTH2F, {axisPtBplus, axisImpParProd}});
          registry.add("hPtProng0RecBg", bPlusCandUnmatch + stringPt + stringPtD + entries, {HistType::kTH2F, {axisPtBplus, axisPtD0}});
          registry.add("hPtProng1RecBg", bPlusCandUnmatch + stringPt + "#it{p}_{T}(#pi^{#plus}) (GeV/#it{c});" + entries, {HistType::kTH2F, {axisPtBplus, axisPtPi}});
          registry.add("hInvMassDRecBg", bPlusCandUnmatch + stringPtD + "prong0, #it{M}(K#pi) (GeV/#it{c}^{2});" + entries, {HistType::kTH2F, {axisPtD0, axisMassD0}});
          registry.add("hDecLengthDRecBg", bPlusCandUnmatch + stringPtD + "D^{0} candidate decay length (cm);" + entries, {HistType::kTH2F, {axisPtD0, axisDecayLength}});
          registry.add("hDecLengthXyDRecBg", bPlusCandUnmatch + stringPtD + "decay length XY (cm);" + entries, {HistType::kTH2F, {axisPtD0, axisDecayLength}});
          registry.add("hCospDRecBg", bPlusCandUnmatch + stringPtD + "D^{0} candidate cos(#vartheta_{P});" + entries, {HistType::kTH2F, {axisPtD0, axisCosp}});
          registry.add("hCospXyDRecBg", bPlusCandUnmatch + stringPtD + "D^{0} candidate cos(#vartheta_{P}^{XY});" + entries, {HistType::kTH2F, {axisPtD0, axisCosp}});
        }
        // MC checks
        if (doprocessMcWithDecayTypeCheck || doprocessMcWithDmesMlAndDecayTypeCheck || doprocessMcWithBplusMlAndDecayTypeCheck) {
          constexpr uint8_t kNBinsDecayTypeMc = hf_cand_bplus::DecayTypeMc::NDecayTypeMc;
          TString labels[kNBinsDecayTypeMc];
          labels[hf_cand_bplus::DecayTypeMc::BplusToD0PiToKPiPi] = "B^{+} #rightarrow (#overline{D^{0}} #rightarrow K^{#plus} #pi^{#minus}) #pi^{#plus}";
          labels[hf_cand_bplus::DecayTypeMc::BplusToD0KToKPiK] = "B^{+} #rightarrow (#overline{D^{0}} #rightarrow K^{#plus} #pi^{#minus}) K^{#plus}";
          labels[hf_cand_bplus::DecayTypeMc::PartlyRecoDecay] = "Partly reconstructed decay channel";
          labels[hf_cand_bplus::DecayTypeMc::OtherDecay] = "Other decays";
          static const AxisSpec axisDecayType = {kNBinsDecayTypeMc, 0.5, kNBinsDecayTypeMc + 0.5, ""};
          registry.add("hDecayTypeMc", "DecayType", {HistType::kTH3F, {axisDecayType, axisMassBplus, axisPtBplus}});
          for (uint8_t iBin = 0; iBin < kNBinsDecayTypeMc; ++iBin) {
            registry.get<TH3>(HIST("hDecayTypeMc"))->GetXaxis()->SetBinLabel(iBin + 1, labels[iBin]);
          }
        }
        // ML scores of D0 daughter
        if (doprocessMcWithDmesMl || doprocessMcWithDmesMlAndDecayTypeCheck) {
          // signal
          registry.add("hMlScoreBkgDRecSig", bPlusCandMatch + stringPtD + "prong0, D^{0} ML background score;" + entries, {HistType::kTH2F, {axisPtD0, axisMlScore}});
          registry.add("hMlScorePromptDRecSig", bPlusCandMatch + stringPtD + "prong0, D^{0} ML prompt score;" + entries, {HistType::kTH2F, {axisPtD0, axisMlScore}});
          registry.add("hMlScoreNonPromptDRecSig", bPlusCandMatch + stringPtD + "prong0, D^{0} ML nonprompt score;" + entries, {HistType::kTH2F, {axisPtD0, axisMlScore}});
          // background
          registry.add("hMlScoreBkgDRecBg", bPlusCandUnmatch + stringPtD + "prong0, D^{0} ML background score;" + entries, {HistType::kTH2F, {axisPtD0, axisMlScore}});
          registry.add("hMlScorePromptDRecBg", bPlusCandUnmatch + stringPtD + "prong0, D^{0} ML prompt score;" + entries, {HistType::kTH2F, {axisPtD0, axisMlScore}});
          registry.add("hMlScoreNonPromptDRecBg", bPlusCandUnmatch + stringPtD + "prong0, D^{0} ML nonprompt score;" + entries, {HistType::kTH2F, {axisPtD0, axisMlScore}});
        }
        // ML scores of B+ candidate
        if (doprocessMcWithBplusMl || doprocessMcWithBplusMlAndDecayTypeCheck) {
          // signal
          registry.add("hMlScoreSigBplusRecSig", bPlusCandMatch + stringPt + "prong0, B^{+} ML signal score;" + entries, {HistType::kTH2F, {axisPtBplus, axisMlScore}});
          // background
          registry.add("hMlScoreSigBplusRecBg", bPlusCandUnmatch + stringPt + "prong0, B^{+} ML signal score;" + entries, {HistType::kTH2F, {axisPtBplus, axisMlScore}});
        }
      }
      if (fillSparses) {
        // gen sparses
        registry.add("hPtYGenSig", bPlusPartGen + stringPt + "#it{y}(B^{+})", {HistType::kTHnSparseF, {axisPtBplus, axisEta}});
        registry.add("hPtYWithProngsInAccepanceGenSig", "B^{+} particles (generated-daughters in acceptance);" + stringPt + "#it{y}(B^{+})", {HistType::kTHnSparseF, {axisPtBplus, axisEta}});

        // reco sparses
        if (!(doprocessMcWithDmesMl || doprocessMcWithDmesMlAndDecayTypeCheck)) {
          registry.add("hMassPtCutVarsRecSig", bPlusCandMatch + "#it{M} (D^{0}#pi) (GeV/#it{c}^{2});" + stringPt + "B^{+} candidate decay length (cm);B^{+} candidate norm. decay length XY (cm);B^{+} candidate impact parameter product (cm);B^{+} candidate cos(#vartheta_{P});#it{M} (K#pi) (GeV/#it{c}^{2});" + stringPtD + "D^{0} candidate decay length (cm);D^{0} candidate cos(#vartheta_{P})", {HistType::kTHnSparseF, {axisMassBplus, axisPtBplus, axisDecayLength, axisNormDecayLength, axisImpParProd, axisCosp, axisMassD0, axisPtD0, axisDecayLength, axisCosp}});
          if (fillBackground) {
            registry.add("hMassPtCutVarsRecBg", bPlusCandUnmatch + "#it{M} (D^{0}#pi) (GeV/#it{c}^{2});" + stringPt + "B^{+} candidate decay length (cm);B^{+} candidate norm. decay length XY (cm);B^{+} candidate impact parameter product (cm);B^{+} candidate cos(#vartheta_{P});#it{M} (K#pi) (GeV/#it{c}^{2});" + stringPtD + "D^{0} candidate decay length (cm);D^{0} candidate cos(#vartheta_{P})", {HistType::kTHnSparseF, {axisMassBplus, axisPtBplus, axisDecayLength, axisNormDecayLength, axisImpParProd, axisCosp, axisMassD0, axisPtD0, axisDecayLength, axisCosp}});
          }
        } else {
          registry.add("hMassPtCutVarsRecSig", bPlusCandMatch + "#it{M} (D^{0}#pi) (GeV/#it{c}^{2});" + stringPt + "B^{+} candidate decay length (cm);B^{+} candidate norm. decay length XY (cm);B^{+} candidate impact parameter product (cm);B^{+} candidate cos(#vartheta_{P});#it{M} (K#pi) (GeV/#it{c}^{2});" + stringPtD + "D^{0} candidate ML score bkg;D^{0} candidate ML score nonprompt", {HistType::kTHnSparseF, {axisMassBplus, axisPtBplus, axisDecayLength, axisNormDecayLength, axisImpParProd, axisCosp, axisMassD0, axisPtD0, axisMlScore, axisMlScore}});
          if (fillBackground) {
            registry.add("hMassPtCutVarsRecBg", bPlusCandUnmatch + "#it{M} (D^{0}#pi) (GeV/#it{c}^{2});" + stringPt + "B^{+} candidate decay length (cm);B^{+} candidate norm. decay length XY (cm);B^{+} candidate impact parameter product (cm);B^{+} candidate cos(#vartheta_{P});#it{M} (K#pi) (GeV/#it{c}^{2});" + stringPtD + "D^{0} candidate ML score bkg;D^{0} candidate ML score nonprompt", {HistType::kTHnSparseF, {axisMassBplus, axisPtBplus, axisDecayLength, axisNormDecayLength, axisImpParProd, axisCosp, axisMassD0, axisPtD0, axisMlScore, axisMlScore}});
          }
        }
      }
    }
  }

  /// Selection of B+ daughter in geometrical acceptance
  /// \param etaProng is the pseudorapidity of B+ prong
  /// \param ptProng is the pT of B+ prong
  /// \return true if prong is in geometrical acceptance
  template <typename T = float>
  bool isProngInAcceptance(const T& etaProng, const T& ptProng)
  {
    return std::abs(etaProng) <= etaTrackMax && ptProng >= ptTrackMin;
  }

  /// Fill candidate information at reconstruction level
  /// \param doMc is the flag to enable the filling with MC information
  /// \param withDecayTypeCheck is the flag to enable MC with decay type check
  /// \param withDmesMl is the flag to enable the filling with ML scores for the D0 daughter
  /// \param withBplusMl is the flag to enable the filling with ML scores for the B+ candidate
  /// \param withScalarProd is the flag to enable the filling with scalar product
  /// \param collision is the collision
  /// \param candidate is the B+ candidate
  /// \param candidatesD is the table with D0 candidates
  template <bool DoMc, bool WithDecayTypeCheck, bool WithDmesMl, bool WithBplusMl, bool WithScalarProd, typename Collision, typename Cand, typename CandsDmes>
  void fillCand(Collision const& collision,
                Cand const& candidate,
                CandsDmes const&)
  {
    auto ptCandBplus = candidate.pt();
    auto invMassBplus = HfHelper::invMassBplusToD0Pi(candidate);
    auto candD0 = candidate.template prong0_as<CandsDmes>();
    auto prong1 = candidate.template prong1_as<TracksPion>();
    auto ptD0 = candidate.ptProng0();
    auto invMassD0 = (prong1.signed1Pt() < 0) ? candD0.invMassHypo0() : candD0.invMassHypo1();
    std::array<float, 3> const posPv{candidate.posX(), candidate.posY(), candidate.posZ()};
    std::array<float, 3> const posSvD{candD0.xSecondaryVertex(), candD0.ySecondaryVertex(), candD0.zSecondaryVertex()};
    std::array<float, 3> const momD{candD0.pVector()};
    auto cpaD0 = RecoDecay::cpa(posPv, posSvD, momD);
    auto cpaXyD0 = RecoDecay::cpaXY(posPv, posSvD, momD);
    auto decLenD0 = RecoDecay::distance(posPv, posSvD);
    auto decLenXyD0 = RecoDecay::distanceXY(posPv, posSvD);

    int8_t flagMcMatchRec = 0;
    int8_t flagWrongCollision = 0;
    bool isSignal = false;
    if constexpr (DoMc) {
      flagMcMatchRec = candidate.flagMcMatchRec();
      flagWrongCollision = candidate.flagWrongCollision();
      isSignal = TESTBIT(std::abs(flagMcMatchRec), hf_cand_bplus::DecayTypeMc::BplusToD0PiToKPiPi);
    }

    if (fillHistograms) {
      if constexpr (DoMc) {
        if (isSignal) {
          registry.fill(HIST("hMassRecSig"), ptCandBplus, invMassBplus);
          registry.fill(HIST("hPtProng0RecSig"), ptCandBplus, candidate.ptProng0());
          registry.fill(HIST("hPtProng1RecSig"), ptCandBplus, candidate.ptProng1());
          registry.fill(HIST("hImpParProdRecSig"), ptCandBplus, candidate.impactParameterProduct());
          registry.fill(HIST("hDecLengthRecSig"), ptCandBplus, candidate.decayLength());
          registry.fill(HIST("hDecLengthXyRecSig"), ptCandBplus, candidate.decayLengthXY());
          registry.fill(HIST("hNormDecLengthXyRecSig"), ptCandBplus, candidate.decayLengthXY() / candidate.errorDecayLengthXY());
          registry.fill(HIST("hDcaProng0RecSig"), ptCandBplus, candidate.impactParameter0());
          registry.fill(HIST("hDcaProng1RecSig"), ptCandBplus, candidate.impactParameter1());
          registry.fill(HIST("hCospRecSig"), ptCandBplus, candidate.cpa());
          registry.fill(HIST("hCospXyRecSig"), ptCandBplus, candidate.cpaXY());
          registry.fill(HIST("hEtaRecSig"), ptCandBplus, candidate.eta());
          registry.fill(HIST("hRapidityRecSig"), ptCandBplus, HfHelper::yBplus(candidate));
          registry.fill(HIST("hInvMassDRecSig"), ptD0, invMassD0);
          registry.fill(HIST("hDecLengthDRecSig"), ptD0, decLenD0);
          registry.fill(HIST("hDecLengthXyDRecSig"), ptD0, decLenXyD0);
          registry.fill(HIST("hCospDRecSig"), ptD0, cpaD0);
          registry.fill(HIST("hCospXyDRecSig"), ptD0, cpaXyD0);
          if constexpr (WithDecayTypeCheck) {
            registry.fill(HIST("hDecayTypeMc"), 1 + hf_cand_bplus::DecayTypeMc::BplusToD0PiToKPiPi, invMassBplus, ptCandBplus);
          }
          if constexpr (WithDmesMl) {
            registry.fill(HIST("hMlScoreBkgDRecSig"), ptD0, candidate.prong0MlScoreBkg());
            registry.fill(HIST("hMlScorePromptDRecSig"), ptD0, candidate.prong0MlScorePrompt());
            registry.fill(HIST("hMlScoreNonPromptDRecSig"), ptD0, candidate.prong0MlScoreNonprompt());
          }
          if constexpr (WithBplusMl) {
            registry.fill(HIST("hMlScoreSigBplusRecSig"), ptCandBplus, candidate.mlProbBplusToD0Pi());
          }
        } else if (fillBackground) {
          registry.fill(HIST("hMassRecBg"), ptCandBplus, invMassBplus);
          registry.fill(HIST("hPtProng0RecBg"), ptCandBplus, candidate.ptProng0());
          registry.fill(HIST("hPtProng1RecBg"), ptCandBplus, candidate.ptProng1());
          registry.fill(HIST("hImpParProdRecBg"), ptCandBplus, candidate.impactParameterProduct());
          registry.fill(HIST("hDecLengthRecBg"), ptCandBplus, candidate.decayLength());
          registry.fill(HIST("hDecLengthXyRecBg"), ptCandBplus, candidate.decayLengthXY());
          registry.fill(HIST("hNormDecLengthXyRecBg"), ptCandBplus, candidate.decayLengthXY() / candidate.errorDecayLengthXY());
          registry.fill(HIST("hDcaProng0RecBg"), ptCandBplus, candidate.impactParameter0());
          registry.fill(HIST("hDcaProng1RecBg"), ptCandBplus, candidate.impactParameter1());
          registry.fill(HIST("hCospRecBg"), ptCandBplus, candidate.cpa());
          registry.fill(HIST("hCospXyRecBg"), ptCandBplus, candidate.cpaXY());
          registry.fill(HIST("hEtaRecBg"), ptCandBplus, candidate.eta());
          registry.fill(HIST("hRapidityRecBg"), ptCandBplus, HfHelper::yBplus(candidate));
          registry.fill(HIST("hInvMassDRecBg"), ptD0, invMassD0);
          registry.fill(HIST("hDecLengthDRecBg"), ptD0, decLenD0);
          registry.fill(HIST("hDecLengthXyDRecBg"), ptD0, decLenXyD0);
          registry.fill(HIST("hCospDRecBg"), ptD0, cpaD0);
          registry.fill(HIST("hCospXyDRecBg"), ptD0, cpaXyD0);
          if constexpr (WithDmesMl) {
            registry.fill(HIST("hMlScoreBkgDRecBg"), ptD0, candidate.prong0MlScoreBkg());
            registry.fill(HIST("hMlScorePromptDRecBg"), ptD0, candidate.prong0MlScorePrompt());
            registry.fill(HIST("hMlScoreNonPromptDRecBg"), ptD0, candidate.prong0MlScoreNonprompt());
          }
          if constexpr (WithBplusMl) {
            registry.fill(HIST("hMlScoreSigBplusRecBg"), ptCandBplus, candidate.mlProbBplusToD0Pi());
          }
        } else if constexpr (WithDecayTypeCheck) {
          if (TESTBIT(flagMcMatchRec, hf_cand_bplus::DecayTypeMc::BplusToD0KToKPiK)) { // B+ → D0bar K+ → (K+ π-) K+
            registry.fill(HIST("hDecayTypeMc"), 1 + hf_cand_bplus::DecayTypeMc::BplusToD0KToKPiK, invMassBplus, ptCandBplus);
          } else if (TESTBIT(flagMcMatchRec, hf_cand_bplus::DecayTypeMc::PartlyRecoDecay)) { // Partly reconstructed decay channel
            registry.fill(HIST("hDecayTypeMc"), 1 + hf_cand_bplus::DecayTypeMc::PartlyRecoDecay, invMassBplus, ptCandBplus);
          } else {
            registry.fill(HIST("hDecayTypeMc"), 1 + hf_cand_bplus::DecayTypeMc::OtherDecay, invMassBplus, ptCandBplus);
          }
        }
      } else {
        registry.fill(HIST("hMass"), ptCandBplus, invMassBplus);
        registry.fill(HIST("hPtProng0"), ptCandBplus, candidate.ptProng0());
        registry.fill(HIST("hPtProng1"), ptCandBplus, candidate.ptProng1());
        registry.fill(HIST("hImpParProd"), ptCandBplus, candidate.impactParameterProduct());
        registry.fill(HIST("hDecLength"), ptCandBplus, candidate.decayLength());
        registry.fill(HIST("hDecLengthXy"), ptCandBplus, candidate.decayLengthXY());
        registry.fill(HIST("hNormDecLengthXy"), ptCandBplus, candidate.decayLengthXY() / candidate.errorDecayLengthXY());
        registry.fill(HIST("hDcaProng0"), ptCandBplus, candidate.impactParameter0());
        registry.fill(HIST("hDcaProng1"), ptCandBplus, candidate.impactParameter1());
        registry.fill(HIST("hCosp"), ptCandBplus, candidate.cpa());
        registry.fill(HIST("hCospXy"), ptCandBplus, candidate.cpaXY());
        registry.fill(HIST("hEta"), ptCandBplus, candidate.eta());
        registry.fill(HIST("hRapidity"), ptCandBplus, HfHelper::yBplus(candidate));
        registry.fill(HIST("hInvMassD"), ptD0, invMassD0);
        registry.fill(HIST("hDecLengthD"), ptD0, decLenD0);
        registry.fill(HIST("hDecLengthXyD"), ptD0, decLenXyD0);
        registry.fill(HIST("hCospD"), ptD0, cpaD0);
        registry.fill(HIST("hCospXyD"), ptD0, cpaXyD0);

        if constexpr (WithDmesMl) {
          registry.fill(HIST("hMlScoreBkgD"), ptD0, candidate.prong0MlScoreBkg());
          registry.fill(HIST("hMlScorePromptD"), ptD0, candidate.prong0MlScorePrompt());
          registry.fill(HIST("hMlScoreNonPromptD"), ptD0, candidate.prong0MlScoreNonprompt());
        }
        if constexpr (WithBplusMl) {
          registry.fill(HIST("hMlScoreSigBplus"), ptCandBplus, candidate.mlProbBplusToD0Pi());
        }
      }
    }
    if (fillSparses) {
      if constexpr (DoMc) {
        if (isSignal) {
          if constexpr (WithDmesMl) {
            registry.fill(HIST("hMassPtCutVarsRecSig"), invMassBplus, ptCandBplus, candidate.decayLength(), candidate.decayLengthXY() / candidate.errorDecayLengthXY(), candidate.impactParameterProduct(), candidate.cpa(), invMassD0, ptD0, candidate.prong0MlScoreBkg(), candidate.prong0MlScoreNonprompt());
          } else {
            registry.fill(HIST("hMassPtCutVarsRecSig"), invMassBplus, ptCandBplus, candidate.decayLength(), candidate.decayLengthXY() / candidate.errorDecayLengthXY(), candidate.impactParameterProduct(), candidate.cpa(), invMassD0, ptD0, decLenD0, cpaD0);
          }
        } else if (fillBackground) {
          if constexpr (WithDmesMl) {
            registry.fill(HIST("hMassPtCutVarsRecBg"), invMassBplus, ptCandBplus, candidate.decayLength(), candidate.decayLengthXY() / candidate.errorDecayLengthXY(), candidate.impactParameterProduct(), candidate.cpa(), invMassD0, ptD0, candidate.prong0MlScoreBkg(), candidate.prong0MlScoreNonprompt());
          } else {
            registry.fill(HIST("hMassPtCutVarsRecBg"), invMassBplus, ptCandBplus, candidate.decayLength(), candidate.decayLengthXY() / candidate.errorDecayLengthXY(), candidate.impactParameterProduct(), candidate.cpa(), invMassD0, ptD0, decLenD0, cpaD0);
          }
        }
      } else {
        if constexpr (WithDmesMl) {
          registry.fill(HIST("hMassPtCutVars"), invMassBplus, ptCandBplus, candidate.decayLength(), candidate.decayLengthXY() / candidate.errorDecayLengthXY(), candidate.impactParameterProduct(), candidate.cpa(), invMassD0, ptD0, candidate.prong0MlScoreBkg(), candidate.prong0MlScoreNonprompt());
        } else {
          registry.fill(HIST("hMassPtCutVars"), invMassBplus, ptCandBplus, candidate.decayLength(), candidate.decayLengthXY() / candidate.errorDecayLengthXY(), candidate.impactParameterProduct(), candidate.cpa(), invMassD0, ptD0, decLenD0, cpaD0);
        }
      }
    }
    if (fillTree) {
      float const pseudoRndm = ptD0 * 1000. - static_cast<int64_t>(ptD0 * 1000);
      if (flagMcMatchRec != 0 || (((DoMc && fillBackground) || !DoMc) && (ptCandBplus >= ptMaxForDownSample || pseudoRndm < downSampleBkgFactor))) {
        float prong0MlScoreBkg = -1.;
        float prong0MlScorePrompt = -1.;
        float prong0MlScoreNonprompt = -1.;
        float candidateMlScoreSig = -1;
        float scalarProd = -1;
        if constexpr (WithDmesMl) {
          prong0MlScoreBkg = candidate.prong0MlScoreBkg();
          prong0MlScorePrompt = candidate.prong0MlScorePrompt();
          prong0MlScoreNonprompt = candidate.prong0MlScoreNonprompt();
        }
        if constexpr (WithBplusMl) {
          candidateMlScoreSig = candidate.mlProbBplusToD0Pi();
        }
        if constexpr (WithScalarProd) {
          auto qVecs = getQvec(collision, qVecDetector.value);
          float xQVec = qVecs[0];
          float yQVec = qVecs[1];
          scalarProd = std::cos(harmonic * candidate.phi()) * xQVec + std::sin(harmonic * candidate.phi()) * yQVec;
        }
        float tpcNSigmaPi, tofNSigmaPi, tpcTofNSigmaPi, tpcNSigmaKa, tofNSigmaKa, tpcTofNSigmaKa;
        if (prong1.signed1Pt() < 0) {
          tpcNSigmaPi = candD0.tpcNSigmaPiProng1();
          tofNSigmaPi = candD0.tofNSigmaPiProng1();
          tpcTofNSigmaPi = candD0.tpcTofNSigmaPiProng1();
          tpcNSigmaKa = candD0.tpcNSigmaKaProng0();
          tofNSigmaKa = candD0.tofNSigmaKaProng0();
          tpcTofNSigmaKa = candD0.tpcTofNSigmaKaProng0();
        } else {
          tpcNSigmaPi = candD0.tpcNSigmaPiProng0();
          tofNSigmaPi = candD0.tofNSigmaPiProng0();
          tpcTofNSigmaPi = candD0.tpcTofNSigmaPiProng0();
          tpcNSigmaKa = candD0.tpcNSigmaKaProng1();
          tofNSigmaKa = candD0.tofNSigmaKaProng1();
          tpcTofNSigmaKa = candD0.tpcTofNSigmaKaProng1();
        }

        float ptMother = -1.;
        if constexpr (DoMc) {
          ptMother = candidate.ptMother();
        }

        hfRedCandBpLite(
          // B-meson features
          invMassBplus,
          ptCandBplus,
          candidate.eta(),
          candidate.phi(),
          HfHelper::yBplus(candidate),
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
          candidate.isSelBplusToD0Pi(),
          scalarProd,
          centralityEstimator == 0 ? collision.centFT0M() : collision.centFT0C(),
          // D-meson features
          invMassD0,
          ptD0,
          candidate.decayLengthDToPv(),
          candidate.decayLengthDToPvXY(),
          candidate.decayLengthDToB(),
          candidate.decayLengthDToBXY(),
          candidate.cpaDToPv(),
          candidate.cpaDToPvXY(),
          candidate.cpaDToB(),
          candidate.cpaDToBXY(),
          candidate.impactParameter0(),
          candD0.ptProngMin(),
          candD0.absEtaProngMin(),
          candD0.itsNClsProngMin(),
          candD0.tpcNClsCrossedRowsProngMin(),
          candD0.tpcChi2NClProngMax(),
          tpcNSigmaPi,
          tofNSigmaPi,
          tpcTofNSigmaPi,
          tpcNSigmaKa,
          tofNSigmaKa,
          tpcTofNSigmaKa,
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
          hfRedBpMcCheck(
            flagMcMatchRec,
            flagWrongCollision,
            invMassD0,
            ptD0,
            invMassBplus,
            ptCandBplus,
            candidateMlScoreSig,
            candidate.pdgCodeBeautyMother(),
            candidate.pdgCodeCharmMother(),
            candidate.pdgCodeProng0(),
            candidate.pdgCodeProng1(),
            candidate.pdgCodeProng2());
        }
      }
    }
  }

  /// Fill particle histograms (gen MC truth)
  void fillCandMcGen(aod::HfMcGenRedBps::iterator const& particle)
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

      // generated B+ with daughters in geometrical acceptance
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
  void processData(CollisionsWithCents const& collisions,
                   soa::Filtered<soa::Join<aod::HfRedCandBplus, aod::HfSelBplusToD0Pi>> const& candidates,
                   CandsD0 const& candidatesD,
                   TracksPion const&)
  {
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(HfHelper::yBplus(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCand<false, false, false, false, false>(collisions.rawIteratorAt(candidate.collisionId()), candidate, candidatesD);
    } // candidate loop
  } // processData
  PROCESS_SWITCH(HfTaskBplusReduced, processData, "Process data without ML scores for B+ and D0 daughter", true);

  void processDataWithDmesMl(CollisionsWithCents const& collisions,
                             soa::Filtered<soa::Join<aod::HfRedCandBplus, aod::HfRedBplusD0Mls, aod::HfSelBplusToD0Pi>> const& candidates,
                             CandsD0 const& candidatesD,
                             TracksPion const&)
  {
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(HfHelper::yBplus(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCand<false, false, true, false, false>(collisions.rawIteratorAt(candidate.collisionId()), candidate, candidatesD);
    } // candidate loop
  } // processDataWithDmesMl
  PROCESS_SWITCH(HfTaskBplusReduced, processDataWithDmesMl, "Process data with(out) ML scores for D0 daughter (B+)", false);

  void processDataWithDmesMlScalarProd(CollisionsWithCentsQvecs const& collisions,
                                       soa::Filtered<soa::Join<aod::HfRedCandBplus, aod::HfRedBplusD0Mls, aod::HfSelBplusToD0Pi>> const& candidates,
                                       CandsD0 const& candidatesD,
                                       TracksPion const&)
  {
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(HfHelper::yBplus(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCand<false, false, true, false, true>(collisions.rawIteratorAt(candidate.collisionId()), candidate, candidatesD);
    } // candidate loop
  } // processDataWithDmesMlScalarProd
  PROCESS_SWITCH(HfTaskBplusReduced, processDataWithDmesMlScalarProd, "Process data with(out) ML scores for D0 daughter (B+) and scalar products", false);

  void processDataWithBplusMl(CollisionsWithCents const& collisions,
                              soa::Filtered<soa::Join<aod::HfRedCandBplus, aod::HfMlBplusToD0Pi, aod::HfSelBplusToD0Pi>> const& candidates,
                              CandsD0 const& candidatesD,
                              TracksPion const&)
  {
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(HfHelper::yBplus(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCand<false, false, false, true, false>(collisions.rawIteratorAt(candidate.collisionId()), candidate, candidatesD);
    } // candidate loop
  } // processDataWithBplusMl
  PROCESS_SWITCH(HfTaskBplusReduced, processDataWithBplusMl, "Process data with(out) ML scores for B+ (D0 daughter)", false);

  void processDataWithBplusMlScalarProd(CollisionsWithCentsQvecs const& collisions,
                                        soa::Filtered<soa::Join<aod::HfRedCandBplus, aod::HfMlBplusToD0Pi, aod::HfSelBplusToD0Pi>> const& candidates,
                                        CandsD0 const& candidatesD,
                                        TracksPion const&)
  {
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(HfHelper::yBplus(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCand<false, false, false, true, true>(collisions.rawIteratorAt(candidate.collisionId()), candidate, candidatesD);
    } // candidate loop
  } // processDataWithBplusMlScalarProd
  PROCESS_SWITCH(HfTaskBplusReduced, processDataWithBplusMlScalarProd, "Process data with(out) ML scores for B+ (D0 daughter) and scalar products", false);

  void processMc(CollisionsWithCents const& collisions,
                 soa::Filtered<soa::Join<aod::HfRedCandBplus, aod::HfSelBplusToD0Pi, aod::HfMcRecRedBps>> const& candidates,
                 aod::HfMcGenRedBps const& mcParticles,
                 CandsD0 const& candidatesD,
                 TracksPion const&)
  {
    // MC rec
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(HfHelper::yBplus(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCand<true, false, false, false, false>(collisions.rawIteratorAt(candidate.collisionId()), candidate, candidatesD);
    } // rec

    // MC gen. level
    for (const auto& particle : mcParticles) {
      fillCandMcGen(particle);
    } // gen
  } // processMc
  PROCESS_SWITCH(HfTaskBplusReduced, processMc, "Process MC without ML scores for B+ and D0 daughter", false);

  void processMcWithDecayTypeCheck(CollisionsWithCents const& collisions,
                                   soa::Filtered<soa::Join<aod::HfRedCandBplus, aod::HfSelBplusToD0Pi, aod::HfMcRecRedBps, aod::HfMcCheckBps>> const& candidates,
                                   aod::HfMcGenRedBps const& mcParticles,
                                   CandsD0 const& candidatesD,
                                   TracksPion const&)
  {
    // MC rec
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(HfHelper::yBplus(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCand<true, true, false, false, false>(collisions.rawIteratorAt(candidate.collisionId()), candidate, candidatesD);
    } // rec

    // MC gen. level
    for (const auto& particle : mcParticles) {
      fillCandMcGen(particle);
    } // gen
  } // processMcWithDecayTypeCheck
  PROCESS_SWITCH(HfTaskBplusReduced, processMcWithDecayTypeCheck, "Process MC with decay type check and without ML scores for B+ and D0 daughter", false);

  void processMcWithDmesMl(CollisionsWithCents const& collisions,
                           soa::Filtered<soa::Join<aod::HfRedCandBplus, aod::HfRedBplusD0Mls, aod::HfSelBplusToD0Pi, aod::HfMcRecRedBps>> const& candidates,
                           aod::HfMcGenRedBps const& mcParticles,
                           CandsD0 const& candidatesD,
                           TracksPion const&)
  {
    // MC rec
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(HfHelper::yBplus(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCand<true, false, true, false, false>(collisions.rawIteratorAt(candidate.collisionId()), candidate, candidatesD);
    } // rec

    // MC gen. level
    for (const auto& particle : mcParticles) {
      fillCandMcGen(particle);
    } // gen
  } // processMcWithDmesMl
  PROCESS_SWITCH(HfTaskBplusReduced, processMcWithDmesMl, "Process MC with(out) ML scores for D0 daughter (B+)", false);

  void processMcWithDmesMlAndDecayTypeCheck(CollisionsWithCents const& collisions,
                                            soa::Filtered<soa::Join<aod::HfRedCandBplus, aod::HfRedBplusD0Mls, aod::HfSelBplusToD0Pi, aod::HfMcRecRedBps, aod::HfMcCheckBps>> const& candidates,
                                            aod::HfMcGenRedBps const& mcParticles,
                                            CandsD0 const& candidatesD,
                                            TracksPion const&)
  {
    // MC rec
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(HfHelper::yBplus(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCand<true, true, true, false, false>(collisions.rawIteratorAt(candidate.collisionId()), candidate, candidatesD);
    } // rec

    // MC gen. level
    for (const auto& particle : mcParticles) {
      fillCandMcGen(particle);
    } // gen
  } // processMcWithDmesMlAndDecayTypeCheck
  PROCESS_SWITCH(HfTaskBplusReduced, processMcWithDmesMlAndDecayTypeCheck, "Process MC with decay type check and with(out) ML scores for B+ (D0 daughter)", false);

  void processMcWithBplusMl(CollisionsWithCents const& collisions,
                            soa::Filtered<soa::Join<aod::HfRedCandBplus, aod::HfMlBplusToD0Pi, aod::HfSelBplusToD0Pi, aod::HfMcRecRedBps>> const& candidates,
                            aod::HfMcGenRedBps const& mcParticles,
                            CandsD0 const& candidatesD,
                            TracksPion const&)
  {
    // MC rec
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(HfHelper::yBplus(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCand<true, false, false, true, false>(collisions.rawIteratorAt(candidate.collisionId()), candidate, candidatesD);
    } // rec

    // MC gen. level
    for (const auto& particle : mcParticles) {
      fillCandMcGen(particle);
    } // gen
  } // processMcWithBplusMl
  PROCESS_SWITCH(HfTaskBplusReduced, processMcWithBplusMl, "Process MC with(out) ML scores for B+ (D0 daughter)", false);

  void processMcWithBplusMlAndDecayTypeCheck(CollisionsWithCents const& collisions,
                                             soa::Filtered<soa::Join<aod::HfRedCandBplus, aod::HfMlBplusToD0Pi, aod::HfSelBplusToD0Pi, aod::HfMcRecRedBps, aod::HfMcCheckBps>> const& candidates,
                                             aod::HfMcGenRedBps const& mcParticles,
                                             CandsD0 const& candidatesD,
                                             TracksPion const&)
  {
    // MC rec
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(HfHelper::yBplus(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCand<true, true, false, true, false>(collisions.rawIteratorAt(candidate.collisionId()), candidate, candidatesD);
    } // rec

    // MC gen. level
    for (const auto& particle : mcParticles) {
      fillCandMcGen(particle);
    } // gen
  } // processMcWithBplusMlAndDecayTypeCheck
  PROCESS_SWITCH(HfTaskBplusReduced, processMcWithBplusMlAndDecayTypeCheck, "Process MC with decay type check and with(out) ML scores for B+ (D0 daughter)", false);
}; // struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskBplusReduced>(cfgc)};
}
