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
#include "PWGHF/Core/SelectorCuts.h"
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
#include <vector>

using namespace o2;
using namespace o2::aod;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::framework::expressions;

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
DECLARE_SOA_COLUMN(ImpactParameterD, impactParameterD, float);                           //! Impact parameter product of D-meson daughter candidate
DECLARE_SOA_COLUMN(ImpactParameterBach, impactParameterBach, float);                     //! Impact parameter product of bachelor pion
DECLARE_SOA_COLUMN(ImpactParameterProduct, impactParameterProduct, float);               //! Impact parameter product of daughters
DECLARE_SOA_COLUMN(Cpa, cpa, float);                                                     //! Cosine pointing angle of candidate
DECLARE_SOA_COLUMN(CpaXY, cpaXY, float);                                                 //! Cosine pointing angle of candidate in transverse plane
DECLARE_SOA_COLUMN(CpaD, cpaD, float);                                                   //! Cosine pointing angle of D-meson daughter candidate
DECLARE_SOA_COLUMN(CpaXYD, cpaXYD, float);                                               //! Cosine pointing angle in transverse plane of D-meson daughter candidate
DECLARE_SOA_COLUMN(MaxNormalisedDeltaIP, maxNormalisedDeltaIP, float);                   //! Maximum normalized difference between measured and expected impact parameter of candidate prongs
DECLARE_SOA_COLUMN(MlScoreSig, mlScoreSig, float);                                       //! ML score for signal class
DECLARE_SOA_COLUMN(FlagWrongCollision, flagWrongCollision, int8_t);                      //! Flag for association with wrong collision
} // namespace hf_cand_bplus_lite

DECLARE_SOA_TABLE(HfRedCandBpLites, "AOD", "HFREDCANDBPLITE", //! Table with some B+ properties
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
                  // D meson features
                  hf_cand_bplus_lite::MD,
                  hf_cand_bplus_lite::PtD,
                  hf_cand_bplus_lite::DecayLengthD,
                  hf_cand_bplus_lite::DecayLengthXYD,
                  hf_cand_bplus_lite::ImpactParameterD,
                  hf_cand_bplus_lite::CpaD,
                  hf_cand_bplus_lite::CpaXYD,
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
const TString stringPt = "#it{p}_{T} (GeV/#it{c})";
const TString stringPtD = "#it{p}_{T}(D0) (GeV/#it{c});";
const TString bPlusCandTitle = "B+ candidates;";
const TString entries = "entries";
const TString bPlusCandMatch = "B+ candidates (matched);";
const TString bPlusCandUnmatch = "B+ candidates (unmatched);";
const TString mcParticleMatched = "MC particles (matched);";

/// B+ analysis task
struct HfTaskBplusReduced {
  Produces<aod::HfRedCandBpLites> hfRedCandBpLite;
  Produces<aod::HfRedBpMcCheck> hfRedBpMcCheck;

  Configurable<int> selectionFlagBplus{"selectionFlagBplus", 1, "Selection Flag for Bplus"};
  Configurable<double> yCandGenMax{"yCandGenMax", 0.5, "max. gen particle rapidity"};
  Configurable<double> yCandRecoMax{"yCandRecoMax", 0.8, "max. cand. rapidity"};
  Configurable<float> etaTrackMax{"etaTrackMax", 0.8, "max. track pseudo-rapidity"};
  Configurable<float> ptTrackMin{"ptTrackMin", 0.1, "min. track transverse momentum"};
  Configurable<bool> fillHistograms{"fillHistograms", true, "Flag to enable histogram filling"};
  Configurable<bool> fillSparses{"fillSparses", false, "Flag to enable sparse filling"};
  Configurable<bool> fillTree{"fillTree", false, "Flag to enable tree filling"};
  Configurable<bool> fillBackground{"fillBackground", false, "Flag to enable filling of background histograms/sparses/tree (only MC)"};
  Configurable<float> downSampleBkgFactor{"downSampleBkgFactor", 1., "Fraction of background candidates to keep for ML trainings"};
  Configurable<float> ptMaxForDownSample{"ptMaxForDownSample", 10., "Maximum pt for the application of the downsampling factor"};
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_bplus_to_d0_pi::vecBinsPt}, "pT bin limits"};

  using TracksPion = soa::Join<HfRedTracks, HfRedTracksPid>;
  using CandsD0 = soa::Join<HfRed2Prongs, HfRedPidDau0s, HfRedPidDau1s>;

  Filter filterSelectCandidates = (aod::hf_sel_candidate_bplus::isSelBplusToD0Pi >= selectionFlagBplus);

  HistogramRegistry registry{"registry"};

  void init(InitContext&)
  {
    std::array<bool, 3> processFuncData{doprocessData, doprocessDataWithDmesMl, doprocessDataWithBplusMl};
    if ((std::accumulate(processFuncData.begin(), processFuncData.end(), 0)) > 1) {
      LOGP(fatal, "Only one process function for data can be enabled at a time.");
    }
    std::array<bool, 6> processFuncMc{doprocessMc, doprocessMcWithDecayTypeCheck, doprocessMcWithDmesMl, doprocessMcWithDmesMlAndDecayTypeCheck, doprocessMcWithBplusMl, doprocessMcWithBplusMlAndDecayTypeCheck};
    if ((std::accumulate(processFuncMc.begin(), processFuncMc.end(), 0)) > 1) {
      LOGP(fatal, "Only one process function for MC can be enabled at a time.");
    }

    const AxisSpec axisMlScore{100, 0.f, 1.f};
    const AxisSpec axisMassBplus{150, 4.5, 6.0};
    const AxisSpec axisMassD0{300, 1.75f, 2.05f};
    const AxisSpec axisCpa{120, -1.1, 1.1};
    const AxisSpec axisCpaD{101, 0.9, 1.01};
    const AxisSpec axisPtProng{100, 0., 10.};
    const AxisSpec axisD0Prong{200, -0.05, 0.05};
    const AxisSpec axisImpParProd{200, -0.001, 0.001};
    const AxisSpec axisDecLength{100, 0., 0.5};
    const AxisSpec axisNormDecLength{40, 0., 20};
    const AxisSpec axisEta{100, -2., 2.};
    const AxisSpec axisRapidity{100, -2., 2.};
    const AxisSpec axisPtD0{100, 0., 50.};
    const AxisSpec axisPtB{(std::vector<double>)binsPt, "#it{p}_{T}^{B^{+}} (GeV/#it{c})"};
    const AxisSpec axisPtPi{100, 0.f, 10.f};

    if (doprocessData || doprocessDataWithDmesMl || doprocessDataWithBplusMl) {
      if (fillHistograms) {
        registry.add("hMass", bPlusCandTitle + "inv. mass #bar{D^{0}}#pi^{+} (GeV/#it{c}^{2});" + stringPt, {HistType::kTH2F, {axisMassBplus, axisPtB}});
        registry.add("hPtCand", bPlusCandTitle + "candidate #it{p}_{T} (GeV/#it{c});" + entries, {HistType::kTH1F, {{1000, 0., 50.}}});
        registry.add("hPtProng0", bPlusCandTitle + "prong 0 #it{p}_{T}^{gen} (GeV/#it{c});" + stringPt, {HistType::kTH2F, {axisPtProng, axisPtB}});
        registry.add("hPtProng1", bPlusCandTitle + "prong 1 #it{p}_{T}^{gen} (GeV/#it{c});" + stringPt, {HistType::kTH2F, {axisPtProng, axisPtB}});
        registry.add("hDecLength", bPlusCandTitle + "decay length (cm);" + stringPt, {HistType::kTH2F, {axisDecLength, axisPtB}});
        registry.add("hDecLengthXy", bPlusCandTitle + "decay length xy (cm);" + stringPt, {HistType::kTH2F, {axisDecLength, axisPtB}});
        registry.add("hNormDecLengthXy", bPlusCandTitle + "Norm. decay length xy (cm);" + stringPt, {HistType::kTH2F, {axisNormDecLength, axisPtB}});
        registry.add("hd0Prong0", bPlusCandTitle + "prong 0 DCAxy to prim. vertex (cm);" + stringPt, {HistType::kTH2F, {axisD0Prong, axisPtB}});
        registry.add("hd0Prong1", bPlusCandTitle + "prong 1 DCAxy to prim. vertex (cm);" + stringPt, {HistType::kTH2F, {axisD0Prong, axisPtB}});
        registry.add("hCpa", bPlusCandTitle + "candidate cosine of pointing angle;" + stringPt, {HistType::kTH2F, {axisCpa, axisPtB}});
        registry.add("hCpaXy", bPlusCandTitle + "candidate cosine of pointing angle xy;" + stringPt, {HistType::kTH2F, {axisCpa, axisPtB}});
        registry.add("hEta", bPlusCandTitle + "candidate #it{#eta};" + stringPt, {HistType::kTH2F, {axisEta, axisPtB}});
        registry.add("hRapidity", bPlusCandTitle + "candidate #it{y};" + stringPt, {HistType::kTH2F, {axisRapidity, axisPtB}});
        registry.add("hd0d0", bPlusCandTitle + "candidate product of DCAxy to prim. vertex (cm^{2});" + stringPt, {HistType::kTH2F, {axisImpParProd, axisPtB}});
        registry.add("hInvMassD0", bPlusCandTitle + "prong0, D0 inv. mass (GeV/#it{c}^{2});" + stringPt, {HistType::kTH2F, {axisMassD0, axisPtD0}});
        registry.add("hDecLengthD0", bPlusCandTitle + "D^{0} candidate decay length (cm);#it{p}_{T}(D^{0}) (GeV/#it{c});entries", {HistType::kTH2F, {axisDecLength, axisPtD0}});
        registry.add("hDecLengthXyD0", bPlusCandTitle + "decay length XY (cm);#it{p}_{T}(D^{0}) (GeV/#it{c});entries", {HistType::kTH2F, {axisDecLength, axisPtD0}});
        registry.add("hCpaD0", bPlusCandTitle + "D^{0} candidate cos(#vartheta_{P});#it{p}_{T}(D^{0}) (GeV/#it{c});entries", {HistType::kTH2F, {axisCpaD, axisPtD0}});
        registry.add("hCpaXyD0", bPlusCandTitle + "D^{0} candidate cos(#vartheta_{P}^{XY});#it{p}_{T}(D^{0}) (GeV/#it{c});entries", {HistType::kTH2F, {axisCpaD, axisPtD0}});

        // ML scores of D0 daughter
        if (doprocessDataWithDmesMl) {
          registry.add("hMlScoreBkgD", bPlusCandTitle + "#it{p}_{T}(D0) (GeV/#it{c});prong0, D0 ML background score;entries", {HistType::kTH2F, {axisPtD0, axisMlScore}});
          registry.add("hMlScorePromptD", bPlusCandTitle + "#it{p}_{T}(D0) (GeV/#it{c});prong0, D0 ML prompt score;entries", {HistType::kTH2F, {axisPtD0, axisMlScore}});
          registry.add("hMlScoreNonPromptD", bPlusCandTitle + "#it{p}_{T}(D0) (GeV/#it{c});prong0, D0 ML nonprompt score;entries", {HistType::kTH2F, {axisPtD0, axisMlScore}});
        }

        // ML scores of B+ candidate
        if (doprocessDataWithBplusMl) {
          registry.add("hMlScoreSigBplus", bPlusCandTitle + "#it{p}_{T} (GeV/#it{c});prong0, B^{+} ML signal score;entries", {HistType::kTH2F, {axisPtB, axisMlScore}});
        }
      }
      if (fillSparses) {
        if (!(doprocessDataWithDmesMl || doprocessDataWithBplusMl)) {
          registry.add("hMassPtCutVars", bPlusCandTitle + "#it{M} (D0#pi) (GeV/#it{c}^{2});#it{p}_{T}(B^{+}) (GeV/#it{c});B^{+} candidate decay length (cm);B^{+} candidate norm. decay length XY (cm);B^{+} candidate impact parameter product (cm);B^{+} candidate cos(#vartheta_{P});#it{M} (K#pi) (GeV/#it{c}^{2});#it{p}_{T}(D0) (GeV/#it{c});D0 candidate decay length (cm);D0 candidate cos(#vartheta_{P})", {HistType::kTHnSparseF, {axisMassBplus, axisPtB, axisDecLength, axisNormDecLength, axisImpParProd, axisCpa, axisMassD0, axisPtD0, axisDecLength, axisCpa}});
        } else {
          registry.add("hMassPtCutVars", bPlusCandTitle + "#it{M} (D0#pi) (GeV/#it{c}^{2});#it{p}_{T}(B^{+}) (GeV/#it{c});B^{+} candidate decay length (cm);B^{+} candidate norm. decay length XY (cm);B^{+} candidate impact parameter product (cm);B^{+} candidate cos(#vartheta_{P});#it{M} (K#pi) (GeV/#it{c}^{2});#it{p}_{T}(D0) (GeV/#it{c});D0 candidate ML score bkg;D0 candidate ML score nonprompt", {HistType::kTHnSparseF, {axisMassBplus, axisPtB, axisDecLength, axisNormDecLength, axisImpParProd, axisCpa, axisMassD0, axisPtD0, axisMlScore, axisMlScore}});
        }
      }
    }

    // histograms processMC
    if (doprocessMc || doprocessMcWithDecayTypeCheck || doprocessMcWithDmesMl || doprocessMcWithDmesMlAndDecayTypeCheck || doprocessMcWithBplusMl || doprocessMcWithBplusMlAndDecayTypeCheck) {
      if (fillHistograms) {
        //  Gen Level
        registry.add("hEtaGen", mcParticleMatched + "candidate #it{#eta}^{gen};" + stringPt, {HistType::kTH2F, {axisEta, axisPtB}});
        registry.add("hYGen", mcParticleMatched + "candidate #it{y}^{gen};" + stringPt, {HistType::kTH2F, {axisRapidity, axisPtB}});
        registry.add("hPtProng0Gen", mcParticleMatched + "prong 0 #it{p}_{T}^{gen} (GeV/#it{c});" + stringPt, {HistType::kTH2F, {axisPtProng, axisPtB}});
        registry.add("hPtProng1Gen", mcParticleMatched + "prong 1 #it{p}_{T}^{gen} (GeV/#it{c});" + stringPt, {HistType::kTH2F, {axisPtProng, axisPtB}});
        registry.add("hYProng0Gen", mcParticleMatched + "prong 0 #it{y}^{gen};" + stringPt, {HistType::kTH2F, {axisRapidity, axisPtB}});
        registry.add("hYProng1Gen", mcParticleMatched + "prong 1 #it{y}^{gen};" + stringPt, {HistType::kTH2F, {axisRapidity, axisPtB}});
        registry.add("hEtaProng0Gen", mcParticleMatched + "prong 0 #it{#eta}^{gen};" + stringPt, {HistType::kTH2F, {axisEta, axisPtB}});
        registry.add("hEtaProng1Gen", mcParticleMatched + "prong 1 #it{#eta}^{gen};" + stringPt, {HistType::kTH2F, {axisEta, axisPtB}});
        registry.add("hPtProngsVsPtBGen", mcParticleMatched + "prong 0 #it{p}_{T}^{gen} (GeV/#it{c});prong 1 #it{p}_{T}^{gen} (GeV/#it{c});" + stringPt, {HistType::kTH3F, {axisPtProng, axisPtProng, axisPtB}});
        registry.add("hYProngsVsBplusGen", mcParticleMatched + "prong 0 #it{y}^{gen};prong 1 #it{y}^{gen};" + stringPt, {HistType::kTH3F, {axisRapidity, axisRapidity, axisPtB}});
        registry.add("hEtaProngsVsBplusGen", mcParticleMatched + "prong 0 #it{#eta}^{gen};prong 1 #it{#eta}^{gen};" + stringPt, {HistType::kTH3F, {axisEta, axisEta, axisPtB}});
        registry.add("hYGenWithProngsInAcceptance", "MC particles (generated-daughters in acceptance);B^{+} candidate #it{y}^{gen};entries", {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
        // Reco Level - Signal
        registry.add("hPtRecSig", bPlusCandMatch + "candidate #it{p}_{T} (GeV/#it{c});" + entries, {HistType::kTH1F, {{300, 0., 30.}}});
        registry.add("hCpaRecSig", bPlusCandMatch + "candidate cosine of pointing angle;" + stringPt, {HistType::kTH2F, {axisCpa, axisPtB}});
        registry.add("hCpaXyRecSig", bPlusCandMatch + "candidate CpaXy;" + stringPt, {HistType::kTH2F, {axisCpa, axisPtB}});
        registry.add("hEtaRecSig", bPlusCandMatch + "candidate #it{#eta};" + stringPt, {HistType::kTH2F, {axisEta, axisPtB}});
        registry.add("hRapidityRecSig", bPlusCandMatch + "candidate #it{y};" + stringPt, {HistType::kTH2F, {axisRapidity, axisPtB}});
        registry.add("hPtProng0RecSig", bPlusCandMatch + "prong 0 #it{p}_{T} (GeV/#it{c});" + stringPt, {HistType::kTH2F, {axisPtProng, axisPtB}});
        registry.add("hPtProng1RecSig", bPlusCandMatch + "prong 1 #it{p}_{T} (GeV/#it{c});" + stringPt, {HistType::kTH2F, {axisPtProng, axisPtB}});
        registry.add("hMassRecSig", bPlusCandMatch + "inv. mass #bar{D^{0}}#pi^{+} (GeV/#it{c}^{2});" + stringPt, {HistType::kTH2F, {axisMassBplus, axisPtB}});
        registry.add("hd0Prong0RecSig", bPlusCandMatch + "prong 0 DCAxy to prim. vertex (cm);" + stringPt, {HistType::kTH2F, {axisD0Prong, axisPtB}});
        registry.add("hd0Prong1RecSig", bPlusCandMatch + "prong 1 DCAxy to prim. vertex (cm);" + stringPt, {HistType::kTH2F, {axisD0Prong, axisPtB}});
        registry.add("hDecLengthRecSig", bPlusCandMatch + "candidate decay length (cm);" + stringPt, {HistType::kTH2F, {axisDecLength, axisPtB}});
        registry.add("hDecLengthXyRecSig", bPlusCandMatch + "candidate decay length xy (cm);" + stringPt, {HistType::kTH2F, {axisDecLength, axisPtB}});
        registry.add("hNormDecLengthXyRecSig", bPlusCandMatch + "candidate normalized decay length (cm);" + stringPt, {HistType::kTH2F, {axisNormDecLength, axisPtB}});
        registry.add("hd0d0RecSig", bPlusCandMatch + "product of DCAxy to prim. vertex (cm^{2});" + stringPt, {HistType::kTH2F, {axisImpParProd, axisPtB}});
        registry.add("hCpaD0RecSig", bPlusCandMatch + "prong0 (D^{0}) cosine of pointing angle;" + stringPtD + entries, {HistType::kTH2F, {{220, 0., 1.1}, {120, 0., 60.}}});
        registry.add("hCpaXyD0RecSig", bPlusCandMatch + "prong0 (D^{0}) cosine of pointing angle;" + stringPtD + entries, {HistType::kTH2F, {{220, 0., 1.1}, {120, 0., 60.}}});
        registry.add("hDecLengthD0RecSig", bPlusCandMatch + "prong0 D^{0} decay length (cm);" + stringPtD + entries, {HistType::kTH2F, {{100, 0., 0.5}, {120, 0., 60.}}});
        registry.add("hDecLengthXyD0RecSig", bPlusCandMatch + "prong0 D^{0} decay length XY (cm);" + stringPtD + entries, {HistType::kTH2F, {{100, 0., 0.5}, {120, 0., 60.}}});
        // background
        if (fillBackground) {
          registry.add("hPtRecBg", bPlusCandUnmatch + "candidate #it{p}_{T} (GeV/#it{c});" + entries, {HistType::kTH1F, {{300, 0., 30.}}});
          registry.add("hCpaRecBg", bPlusCandUnmatch + "candidate cosine of pointing angle;" + stringPt, {HistType::kTH2F, {axisCpa, axisPtB}});
          registry.add("hCpaXyRecBg", bPlusCandUnmatch + "candidate CpaXy;" + stringPt, {HistType::kTH2F, {axisCpa, axisPtB}});
          registry.add("hEtaRecBg", bPlusCandUnmatch + "candidate #it{#eta};" + stringPt, {HistType::kTH2F, {axisEta, axisPtB}});
          registry.add("hRapidityRecBg", bPlusCandUnmatch + "candidate #it{#y};" + stringPt, {HistType::kTH2F, {axisRapidity, axisPtB}});
          registry.add("hPtProngsVsBplusRecBg", bPlusCandUnmatch + "prong 0 #it{p}_{T} (GeV/#it{c});prong 1 #it{p}_{T} (GeV/#it{c});" + stringPt, {HistType::kTH3F, {axisPtProng, axisPtProng, axisPtB}});
          registry.add("hPtProng0RecBg", bPlusCandUnmatch + "prong 0 #it{p}_{T} (GeV/#it{c});" + stringPt, {HistType::kTH2F, {axisPtProng, axisPtB}});
          registry.add("hPtProng1RecBg", bPlusCandUnmatch + "prong 1 #it{p}_{T} (GeV/#it{c});" + stringPt, {HistType::kTH2F, {axisPtProng, axisPtB}});
          registry.add("hMassRecBg", bPlusCandUnmatch + "inv. mass #bar{D^{0}}#pi^{+} (GeV/#it{c}^{2});" + stringPt, {HistType::kTH2F, {axisMassBplus, axisPtB}});
          registry.add("hd0Prong0RecBg", bPlusCandUnmatch + "prong 0 DCAxy to prim. vertex (cm);" + stringPt, {HistType::kTH2F, {axisD0Prong, axisPtB}});
          registry.add("hd0Prong1RecBg", bPlusCandUnmatch + "prong 1 DCAxy to prim. vertex (cm);" + stringPt, {HistType::kTH2F, {axisD0Prong, axisPtB}});
          registry.add("hDecLengthRecBg", bPlusCandUnmatch + "candidate decay length (cm);" + stringPt, {HistType::kTH2F, {axisDecLength, axisPtB}});
          registry.add("hDecLengthXyRecBg", bPlusCandUnmatch + "candidate decay length xy(cm);" + stringPt, {HistType::kTH2F, {axisDecLength, axisPtB}});
          registry.add("hNormDecLengthXyRecBg", bPlusCandUnmatch + "candidate normalized decay length (cm);" + stringPt, {HistType::kTH2F, {axisNormDecLength, axisPtB}});
          registry.add("hd0d0RecBg", bPlusCandUnmatch + "product of DCAxy to prim. vertex (cm^{2});" + stringPt, {HistType::kTH2F, {axisImpParProd, axisPtB}});
          registry.add("hCpaD0RecBg", bPlusCandUnmatch + "prong0 (D^{0}) cosine of pointing angle;" + stringPtD + entries, {HistType::kTH2F, {{220, 0., 1.1}, {120, 0., 60.}}});
          registry.add("hCpaXyD0RecBg", bPlusCandUnmatch + "prong0 (D^{0}) cosine of pointing angle;" + stringPtD + entries, {HistType::kTH2F, {{220, 0., 1.1}, {120, 0., 60.}}});
          registry.add("hDecLengthD0RecBg", bPlusCandUnmatch + "prong0 D^{0} decay length (cm);" + stringPtD + entries, {HistType::kTH2F, {{100, 0., 0.5}, {120, 0., 60.}}});
          registry.add("hDecLengthXyD0RecBg", bPlusCandUnmatch + "prong0 D^{0} decay length XY (cm);" + stringPtD + entries, {HistType::kTH2F, {{100, 0., 0.5}, {120, 0., 60.}}});
        }
        // MC checks
        if (doprocessMcWithDecayTypeCheck || doprocessMcWithDmesMlAndDecayTypeCheck || doprocessMcWithBplusMlAndDecayTypeCheck) {
          constexpr uint8_t kNBinsDecayTypeMc = hf_cand_bplus::DecayTypeMc::NDecayTypeMc;
          TString labels[kNBinsDecayTypeMc];
          labels[hf_cand_bplus::DecayTypeMc::BplusToD0PiToKPiPi] = "B^{+} #rightarrow (#overline{D^{0}} #rightarrow K^{#plus} #pi^{#minus}) #pi^{#plus}";
          labels[hf_cand_bplus::DecayTypeMc::BplusToD0KToKPiK] = "B^{+} #rightarrow (#overline{D^{0}} #rightarrow K^{#plus} #pi^{#minus}) #K^{#plus}";
          labels[hf_cand_bplus::DecayTypeMc::PartlyRecoDecay] = "Partly reconstructed decay channel";
          labels[hf_cand_bplus::DecayTypeMc::OtherDecay] = "Other decays";
          static const AxisSpec axisDecayType = {kNBinsDecayTypeMc, 0.5, kNBinsDecayTypeMc + 0.5, ""};
          registry.add("hDecayTypeMc", "DecayType", {HistType::kTH3F, {axisDecayType, axisMassBplus, axisPtB}});
          for (uint8_t iBin = 0; iBin < kNBinsDecayTypeMc; ++iBin) {
            registry.get<TH3>(HIST("hDecayTypeMc"))->GetXaxis()->SetBinLabel(iBin + 1, labels[iBin]);
          }
        }
        // ML scores of D0 daughter
        if (doprocessMcWithDmesMl || doprocessMcWithDmesMlAndDecayTypeCheck) {
          // signal
          registry.add("hMlScoreBkgDRecSig", bPlusCandMatch + stringPtD + "prong0, D0 ML background score;entries", {HistType::kTH2F, {axisPtD0, axisMlScore}});
          registry.add("hMlScorePromptDRecSig", bPlusCandMatch + stringPtD + "prong0, D0 ML prompt score;entries", {HistType::kTH2F, {axisPtD0, axisMlScore}});
          registry.add("hMlScoreNonPromptDRecSig", bPlusCandMatch + stringPtD + "prong0, D0 ML nonprompt score;entries", {HistType::kTH2F, {axisPtD0, axisMlScore}});
          // background
          registry.add("hMlScoreBkgDRecBg", bPlusCandUnmatch + stringPtD + "prong0, D0 ML background score;entries", {HistType::kTH2F, {axisPtD0, axisMlScore}});
          registry.add("hMlScorePromptDRecBg", bPlusCandUnmatch + stringPtD + "prong0, D0 ML prompt score;entries", {HistType::kTH2F, {axisPtD0, axisMlScore}});
          registry.add("hMlScoreNonPromptDRecBg", bPlusCandUnmatch + stringPtD + "prong0, D0 ML nonprompt score;entries", {HistType::kTH2F, {axisPtD0, axisMlScore}});
        }
        // ML scores of B+ candidate
        if (doprocessMcWithBplusMl || doprocessMcWithBplusMlAndDecayTypeCheck) {
          // signal
          registry.add("hMlScoreSigBplusRecSig", bPlusCandMatch + "#it{p}_{T}(B^{+}) (GeV/#it{c});prong0, B^{+} ML signal score;entries", {HistType::kTH2F, {axisPtB, axisMlScore}});
          // background
          registry.add("hMlScoreSigBplusRecBg", bPlusCandUnmatch + "#it{p}_{T}(B^{+}) (GeV/#it{c});prong0, B^{+} ML signal score;entries", {HistType::kTH2F, {axisPtB, axisMlScore}});
        }
      }
      if (fillSparses) {
        // gen sparses
        registry.add("hPtYGenSig", "B^{+} particles (generated);#it{p}_{T}(B^{+}) (GeV/#it{c});#it{y}(B^{+})", {HistType::kTHnSparseF, {axisPtB, axisEta}});
        registry.add("hPtYWithProngsInAccepanceGenSig", "B^{+} particles (generated-daughters in acceptance);#it{p}_{T}(B^{+}) (GeV/#it{c});#it{y}(B^{+})", {HistType::kTHnSparseF, {axisPtB, axisEta}});
        // reco sparses
        if (!doprocessDataWithDmesMl) {
          registry.add("hMassPtCutVarsRecSig", bPlusCandMatch + "#it{M} (D^{0}#pi) (GeV/#it{c}^{2});#it{p}_{T}(B^{+}) (GeV/#it{c});B^{+} candidate decay length (cm);B^{+} candidate norm. decay length XY (cm);B^{+} candidate impact parameter product (cm);B^{+} candidate cos(#vartheta_{P});#it{M} (K#pi) (GeV/#it{c}^{2});#it{p}_{T}(D0) (GeV/#it{c});D0 candidate decay length (cm);D0 candidate cos(#vartheta_{P})", {HistType::kTHnSparseF, {axisMassBplus, axisPtB, axisDecLength, axisNormDecLength, axisImpParProd, axisCpa, axisMassD0, axisPtD0, axisDecLength, axisCpa}});
          if (fillBackground) {
            registry.add("hMassPtCutVarsRecBg", bPlusCandUnmatch + "#it{M} (D^{0}#pi) (GeV/#it{c}^{2});#it{p}_{T}(B^{+}) (GeV/#it{c});B^{+} candidate decay length (cm);B^{+} candidate norm. decay length XY (cm);B^{+} candidate impact parameter product (cm);B^{+} candidate cos(#vartheta_{P});#it{M} (K#pi) (GeV/#it{c}^{2});#it{p}_{T}(D0) (GeV/#it{c});D0 candidate decay length (cm);D0 candidate cos(#vartheta_{P})", {HistType::kTHnSparseF, {axisMassBplus, axisPtB, axisDecLength, axisNormDecLength, axisImpParProd, axisCpa, axisMassD0, axisPtD0, axisDecLength, axisCpa}});
          }
        } else {
          registry.add("hMassPtCutVarsRecSig", bPlusCandMatch + "#it{M} (D^{0}#pi) (GeV/#it{c}^{2});#it{p}_{T}(B^{+}) (GeV/#it{c});B^{+} candidate decay length (cm);B^{+} candidate norm. decay length XY (cm);B^{+} candidate impact parameter product (cm);B^{+} candidate cos(#vartheta_{P});#it{M} (K#pi) (GeV/#it{c}^{2});#it{p}_{T}(D0) (GeV/#it{c});D0 candidate ML score bkg;D0 candidate ML score nonprompt", {HistType::kTHnSparseF, {axisMassBplus, axisPtB, axisDecLength, axisNormDecLength, axisImpParProd, axisCpa, axisMassD0, axisPtD0, axisMlScore, axisMlScore}});
          if (fillBackground) {
            registry.add("hMassPtCutVarsRecBg", bPlusCandUnmatch + "#it{M} (D^{0}#pi) (GeV/#it{c}^{2});#it{p}_{T}(B^{+}) (GeV/#it{c});B^{+} candidate decay length (cm);B^{+} candidate norm. decay length XY (cm);B^{+} candidate impact parameter product (cm);B^{+} candidate cos(#vartheta_{P});#it{M} (K#pi) (GeV/#it{c}^{2});#it{p}_{T}(D0) (GeV/#it{c});D0 candidate ML score bkg;D0 candidate ML score nonprompt", {HistType::kTHnSparseF, {axisMassBplus, axisPtB, axisDecLength, axisNormDecLength, axisImpParProd, axisCpa, axisMassD0, axisPtD0, axisMlScore, axisMlScore}});
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
  /// \param candidate is the B+ candidate
  /// \param candidatesD is the table with D0 candidates
  template <bool DoMc, bool WithDecayTypeCheck, bool WithDmesMl, bool WithBplusMl, typename Cand, typename CandsDmes>
  void fillCand(Cand const& candidate,
                CandsDmes const& /*candidatesD*/,
                TracksPion const&)
  {
    auto ptCandBplus = candidate.pt();
    auto invMassBplus = HfHelper::invMassBplusToD0Pi(candidate);
    auto candD0 = candidate.template prong0_as<CandsDmes>();
    auto candPi = candidate.template prong1_as<TracksPion>();
    auto ptD0 = candidate.ptProng0();
    auto invMassD0 = (candPi.signed1Pt() < 0) ? candD0.invMassHypo0() : candD0.invMassHypo1();
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
          registry.fill(HIST("hMassRecSig"), invMassBplus, ptCandBplus);
          registry.fill(HIST("hPtRecSig"), ptCandBplus);
          registry.fill(HIST("hPtProng0RecSig"), ptD0, ptCandBplus);
          registry.fill(HIST("hPtProng1RecSig"), candidate.ptProng1(), ptCandBplus);
          registry.fill(HIST("hCpaRecSig"), candidate.cpa(), ptCandBplus);
          registry.fill(HIST("hCpaXyRecSig"), candidate.cpaXY(), ptCandBplus);
          registry.fill(HIST("hEtaRecSig"), candidate.eta(), ptCandBplus);
          registry.fill(HIST("hRapidityRecSig"), HfHelper::yBplus(candidate), ptCandBplus);
          registry.fill(HIST("hDecLengthRecSig"), candidate.decayLength(), ptCandBplus);
          registry.fill(HIST("hDecLengthXyRecSig"), candidate.decayLengthXY(), ptCandBplus);
          registry.fill(HIST("hNormDecLengthXyRecSig"), candidate.decayLengthXYNormalised(), ptCandBplus);
          registry.fill(HIST("hd0Prong0RecSig"), candidate.impactParameter0(), ptCandBplus);
          registry.fill(HIST("hd0Prong1RecSig"), candidate.impactParameter1(), ptCandBplus);
          registry.fill(HIST("hd0d0RecSig"), candidate.impactParameterProduct(), ptCandBplus);
          registry.fill(HIST("hPtProng0RecSig"), candidate.ptProng0(), ptCandBplus);
          registry.fill(HIST("hPtProng1RecSig"), candidate.ptProng1(), ptCandBplus);
          registry.fill(HIST("hDecLengthD0RecSig"), decLenD0, candidate.ptProng0());
          registry.fill(HIST("hDecLengthXyD0RecSig"), decLenXyD0, ptD0);
          registry.fill(HIST("hCpaD0RecSig"), cpaD0, ptD0);
          registry.fill(HIST("hCpaXyD0RecSig"), cpaXyD0, ptD0);
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
          registry.fill(HIST("hMassRecBg"), invMassBplus, ptCandBplus);
          registry.fill(HIST("hPtRecBg"), ptCandBplus);
          registry.fill(HIST("hPtProng0RecBg"), candidate.ptProng0(), ptCandBplus);
          registry.fill(HIST("hPtProng1RecBg"), candidate.ptProng1(), ptCandBplus);
          registry.fill(HIST("hCpaRecBg"), candidate.cpa(), ptCandBplus);
          registry.fill(HIST("hCpaXyRecBg"), candidate.cpaXY(), ptCandBplus);
          registry.fill(HIST("hEtaRecBg"), candidate.eta(), ptCandBplus);
          registry.fill(HIST("hRapidityRecBg"), HfHelper::yBplus(candidate), ptCandBplus);
          registry.fill(HIST("hDecLengthRecBg"), candidate.decayLength(), ptCandBplus);
          registry.fill(HIST("hDecLengthXyRecBg"), candidate.decayLengthXY(), ptCandBplus);
          registry.fill(HIST("hNormDecLengthXyRecBg"), candidate.decayLengthXYNormalised(), ptCandBplus);
          registry.fill(HIST("hd0Prong0RecBg"), candidate.impactParameter0(), ptCandBplus);
          registry.fill(HIST("hd0Prong1RecBg"), candidate.impactParameter1(), ptCandBplus);
          registry.fill(HIST("hd0d0RecBg"), candidate.impactParameterProduct(), ptCandBplus);
          registry.fill(HIST("hDecLengthD0RecBg"), decLenD0, candidate.ptProng0());
          registry.fill(HIST("hInvMassDRecBg"), invMassD0, ptD0);
          registry.fill(HIST("hDecLengthDRecBg"), decLenD0, ptD0);
          registry.fill(HIST("hDecLengthXyDRecBg"), decLenXyD0, ptD0);
          registry.fill(HIST("hCpaDRecBg"), cpaD0, ptD0);
          registry.fill(HIST("hCpaXyDRecBg"), cpaXyD0, ptD0);
          if constexpr (WithDmesMl) {
            registry.fill(HIST("hMlScoreBkgDRecBg"), ptD0, candidate.prong0MlScoreBkg());
            registry.fill(HIST("hMlScorePromptDRecBg"), ptD0, candidate.prong0MlScorePrompt());
            registry.fill(HIST("hMlScoreNonPromptDRecBg"), ptD0, candidate.prong0MlScoreNonprompt());
          }
          if constexpr (WithBplusMl) {
            registry.fill(HIST("hMlScoreSigBplusRecBg"), ptCandBplus, candidate.mlProbBplusToD0Pi());
          }
        } else if constexpr (WithDecayTypeCheck) {
          if (TESTBIT(flagMcMatchRec, hf_cand_bplus::DecayTypeMc::BplusToD0KToKPiK)) { // Partly reconstructed decay channel
            registry.fill(HIST("hDecayTypeMc"), 1 + hf_cand_bplus::DecayTypeMc::BplusToD0KToKPiK, invMassBplus, ptCandBplus);
          } else if (TESTBIT(flagMcMatchRec, hf_cand_bplus::DecayTypeMc::PartlyRecoDecay)) { // Partly reconstructed decay channel
            registry.fill(HIST("hDecayTypeMc"), 1 + hf_cand_bplus::DecayTypeMc::PartlyRecoDecay, invMassBplus, ptCandBplus);
          } else {
            registry.fill(HIST("hDecayTypeMc"), 1 + hf_cand_bplus::DecayTypeMc::OtherDecay, invMassBplus, ptCandBplus);
          }
        }
      } else {
        registry.fill(HIST("hMass"), invMassBplus, ptCandBplus);
        registry.fill(HIST("hPtCand"), ptCandBplus);
        registry.fill(HIST("hPtProng0"), ptD0, ptCandBplus);
        registry.fill(HIST("hPtProng1"), candidate.ptProng1(), ptCandBplus);
        registry.fill(HIST("hd0d0"), candidate.impactParameterProduct(), ptCandBplus);
        registry.fill(HIST("hDecLength"), candidate.decayLength(), ptCandBplus);
        registry.fill(HIST("hDecLengthXy"), candidate.decayLengthXY(), ptCandBplus);
        registry.fill(HIST("hNormDecLengthXy"), candidate.decayLengthXY() / candidate.errorDecayLengthXY(), ptCandBplus);
        registry.fill(HIST("hd0Prong0"), candidate.impactParameter0(), ptCandBplus);
        registry.fill(HIST("hd0Prong1"), candidate.impactParameter1(), ptCandBplus);
        registry.fill(HIST("hCpa"), candidate.cpa(), ptCandBplus);
        registry.fill(HIST("hCpaXy"), candidate.cpaXY(), ptCandBplus);
        registry.fill(HIST("hEta"), candidate.eta(), ptCandBplus);
        registry.fill(HIST("hRapidity"), HfHelper::yBplus(candidate), ptCandBplus);
        registry.fill(HIST("hInvMassD0"), invMassD0, ptCandBplus);
        registry.fill(HIST("hDecLengthD0"), decLenD0, ptD0);
        registry.fill(HIST("hDecLengthXyD0"), decLenXyD0, ptD0);
        registry.fill(HIST("hCpaD0"), cpaD0, ptD0);
        registry.fill(HIST("hCpaXyD0"), cpaXyD0, ptD0);
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
      if (ptCandBplus >= ptMaxForDownSample || pseudoRndm < downSampleBkgFactor) {
        float prong0MlScoreBkg = -1.;
        float prong0MlScorePrompt = -1.;
        float prong0MlScoreNonprompt = -1.;
        float candidateMlScoreSig = -1;
        if constexpr (WithDmesMl) {
          prong0MlScoreBkg = candidate.prong0MlScoreBkg();
          prong0MlScorePrompt = candidate.prong0MlScorePrompt();
          prong0MlScoreNonprompt = candidate.prong0MlScoreNonprompt();
        }
        if constexpr (WithBplusMl) {
          candidateMlScoreSig = candidate.mlProbBplusToD0Pi();
        }
        auto prong1 = candidate.template prong1_as<TracksPion>();
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
          // B+ - meson features
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
          // D-meson features
          invMassD0,
          ptD0,
          decLenD0,
          decLenXyD0,
          candidate.impactParameter0(),
          cpaD0,
          cpaXyD0,
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
  void processData(soa::Filtered<soa::Join<aod::HfRedCandBplus, aod::HfSelBplusToD0Pi>> const& candidates,
                   CandsD0 const& candidatesD,
                   TracksPion const& pionTracks)
  {
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(HfHelper::yBplus(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCand<false, false, false, false>(candidate, candidatesD, pionTracks);
    } // candidate loop
  } // processData
  PROCESS_SWITCH(HfTaskBplusReduced, processData, "Process data without ML scores for D0 daughter", true);

  void processDataWithDmesMl(soa::Filtered<soa::Join<aod::HfRedCandBplus, aod::HfRedBplusD0Mls, aod::HfSelBplusToD0Pi>> const& candidates,
                             CandsD0 const& candidatesD,
                             TracksPion const& pionTracks)
  {
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(HfHelper::yBplus(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCand<false, false, true, false>(candidate, candidatesD, pionTracks);
    } // candidate loop
  } // processDataWithDmesMl
  PROCESS_SWITCH(HfTaskBplusReduced, processDataWithDmesMl, "Process data with ML scores for D0 daughter", false);

  void processDataWithBplusMl(soa::Filtered<soa::Join<aod::HfRedCandBplus, aod::HfMlBplusToD0Pi, aod::HfSelBplusToD0Pi>> const& candidates,
                              CandsD0 const& candidatesD,
                              TracksPion const& pionTracks)
  {
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(HfHelper::yBplus(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCand<false, false, false, true>(candidate, candidatesD, pionTracks);
    } // candidate loop
  } // processDataWithBplusMl
  PROCESS_SWITCH(HfTaskBplusReduced, processDataWithBplusMl, "Process data with(out) ML scores for B+ (D0 daughter)", false);

  void processMc(soa::Filtered<soa::Join<aod::HfRedCandBplus, aod::HfSelBplusToD0Pi, aod::HfMcRecRedBps>> const& candidates,
                 aod::HfMcGenRedBps const& mcParticles,
                 CandsD0 const& candidatesD,
                 TracksPion const& pionTracks)
  {
    // MC rec
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(HfHelper::yBplus(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCand<true, false, false, false>(candidate, candidatesD, pionTracks);
    } // rec

    // MC gen. level
    for (const auto& particle : mcParticles) {
      fillCandMcGen(particle);
    } // gen
  } // processMc
  PROCESS_SWITCH(HfTaskBplusReduced, processMc, "Process MC without ML scores for B+ and D0 daughter", false);

  void processMcWithDecayTypeCheck(soa::Filtered<soa::Join<aod::HfRedCandBplus, aod::HfSelBplusToD0Pi, aod::HfMcRecRedBps, aod::HfMcCheckBps>> const& candidates,
                                   aod::HfMcGenRedBps const& mcParticles,
                                   CandsD0 const& candidatesD,
                                   TracksPion const& pionTracks)
  {
    // MC rec
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(HfHelper::yBplus(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCand<true, true, false, false>(candidate, candidatesD, pionTracks);
    } // rec

    // MC gen. level
    for (const auto& particle : mcParticles) {
      fillCandMcGen(particle);
    } // gen
  } // processMc
  PROCESS_SWITCH(HfTaskBplusReduced, processMcWithDecayTypeCheck, "Process MC with decay type check and without ML scores for B+ and D daughter", false);

  void processMcWithDmesMl(soa::Filtered<soa::Join<aod::HfRedCandBplus, aod::HfRedBplusD0Mls, aod::HfSelBplusToD0Pi, aod::HfMcRecRedBps>> const& candidates,
                           aod::HfMcGenRedBps const& mcParticles,
                           CandsD0 const& candidatesD,
                           TracksPion const& pionTracks)
  {
    // MC rec
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(HfHelper::yBplus(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCand<true, false, true, false>(candidate, candidatesD, pionTracks);
    } // rec

    // MC gen. level
    for (const auto& particle : mcParticles) {
      fillCandMcGen(particle);
    } // gen
  } // processMcWithDmesMl
  PROCESS_SWITCH(HfTaskBplusReduced, processMcWithDmesMl, "Process MC with(out) ML scores for D0 daughter (B+)", false);

  void processMcWithDmesMlAndDecayTypeCheck(soa::Filtered<soa::Join<aod::HfRedCandBplus, aod::HfRedBplusD0Mls, aod::HfSelBplusToD0Pi, aod::HfMcRecRedBps, aod::HfMcCheckBps>> const& candidates,
                                            aod::HfMcGenRedBps const& mcParticles,
                                            CandsD0 const& candidatesD,
                                            TracksPion const& pionTracks)
  {
    // MC rec
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(HfHelper::yBplus(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCand<true, true, true, false>(candidate, candidatesD, pionTracks);
    } // rec

    // MC gen. level
    for (const auto& particle : mcParticles) {
      fillCandMcGen(particle);
    } // gen
  } // processMc
  PROCESS_SWITCH(HfTaskBplusReduced, processMcWithDmesMlAndDecayTypeCheck, "Process MC with decay type check and with(out) ML scores for B+ (D0 daughter)", false);

  void processMcWithBplusMl(soa::Filtered<soa::Join<aod::HfRedCandBplus, aod::HfMlBplusToD0Pi, aod::HfSelBplusToD0Pi, aod::HfMcRecRedBps>> const& candidates,
                            aod::HfMcGenRedBps const& mcParticles,
                            CandsD0 const& candidatesD,
                            TracksPion const& pionTracks)
  {
    // MC rec
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(HfHelper::yBplus(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCand<true, false, false, true>(candidate, candidatesD, pionTracks);
    } // rec

    // MC gen. level
    for (const auto& particle : mcParticles) {
      fillCandMcGen(particle);
    } // gen
  } // processMcWithBplusMl
  PROCESS_SWITCH(HfTaskBplusReduced, processMcWithBplusMl, "Process MC with(out) ML scores for B+ (D0 daughter)", false);

  void processMcWithBplusMlAndDecayTypeCheck(soa::Filtered<soa::Join<aod::HfRedCandBplus, aod::HfMlBplusToD0Pi, aod::HfSelBplusToD0Pi, aod::HfMcRecRedBps, aod::HfMcCheckBps>> const& candidates,
                                             aod::HfMcGenRedBps const& mcParticles,
                                             CandsD0 const& candidatesD,
                                             TracksPion const& pionTracks)
  {
    // MC rec
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(HfHelper::yBplus(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCand<true, true, false, true>(candidate, candidatesD, pionTracks);
    } // rec

    // MC gen. level
    for (const auto& particle : mcParticles) {
      fillCandMcGen(particle);
    } // gen
  } // processMc
  PROCESS_SWITCH(HfTaskBplusReduced, processMcWithBplusMlAndDecayTypeCheck, "Process MC with decay type check and with(out) ML scores for B+ (D0 daughter)", false);
}; // struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskBplusReduced>(cfgc)};
}
