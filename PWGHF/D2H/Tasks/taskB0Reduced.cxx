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
#include <string>

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
                                                              // B meson features
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

  using TracksBachPions = soa::Join<HfRedTracks, HfRedTracksPid>;
  using CandsDplus = soa::Join<HfRed3Prongs, HfRedPidDau0s, HfRedPidDau1s, HfRedPidDau2s>;
  using CandsDstar = soa::Join<HfRed2Prongs, HfRedPidDau0s, HfRedPidDau1s, HfRedSoftPiPid>;
  using TracksSoftPions = soa::Join<aod::HfRedSoftPiBases, aod::HfRedSoftPiCov, aod::HfRedSoftPiPid>;

  Filter filterSelectCandidates = (aod::hf_sel_candidate_b0::isSelB0ToDPi >= selectionFlagB0);

  HistogramRegistry registry{"registry"};

  void init(InitContext&)
  {
    std::array<bool, 6> processFuncData{doprocessDataDplusPi, doprocessDataDplusPiWithDmesMl, doprocessDataDplusPiWithB0Ml,
                                        doprocessDataDstarPi, doprocessDataDstarPiWithDmesMl};
    if ((std::accumulate(processFuncData.begin(), processFuncData.end(), 0)) > 1) {
      LOGP(fatal, "Only one process function for data can be enabled at a time.");
    }
    std::array<bool, 8> processFuncMc{doprocessMcDplusPi, doprocessMcDplusPiWithDecayTypeCheck, doprocessMcDplusPiWithDmesMl, doprocessMcDplusPiWithDmesMlAndDecayTypeCheck, doprocessMcDplusPiWithB0Ml, doprocessMcDplusPiWithB0MlAndDecayTypeCheck,
                                      doprocessMcDstarPi, doprocessMcDstarPiWithDmesMl};
    if ((std::accumulate(processFuncMc.begin(), processFuncMc.end(), 0)) > 1) {
      LOGP(fatal, "Only one process function for MC can be enabled at a time.");
    }

    const AxisSpec axisMlScore{100, 0.f, 1.f};
    const AxisSpec axisMassB0{300, 4.5f, 6.0f};
    const AxisSpec axisMassDminus{300, 1.75f, 2.05f};
    const AxisSpec axisMassDeltaMassDStar{300, 0.05f, 0.3f};
    const AxisSpec axisDecayLength{200, 0.f, 0.4f};
    const AxisSpec axisNormDecayLength{100, 0.f, 50.f};
    const AxisSpec axisDca{100, -0.05f, 0.05f};
    const AxisSpec axisCosp{110, 0.f, 1.1f};
    const AxisSpec axisEta{30, -1.5f, 1.5f};
    const AxisSpec axisError{100, 0.f, 1.f};
    const AxisSpec axisImpParProd{100, -1.e-3, 1.e-3};
    const AxisSpec axisImpParProngSqSum{100, 0, 1.e-3};
    const AxisSpec axisPtB0{100, 0.f, 50.f};
    const AxisSpec axisPtDminus{100, 0.f, 50.f};
    const AxisSpec axisPtPi{100, 0.f, 10.f};
    const AxisSpec axisPtSoftPi{100, 0.f, 1.f};

    std::array<bool, 9> processFuncDplusPi = {doprocessDataDplusPi, doprocessDataDplusPiWithDmesMl, doprocessDataDplusPiWithB0Ml,
                                              doprocessMcDplusPi, doprocessMcDplusPiWithDecayTypeCheck, doprocessMcDplusPiWithDmesMl,
                                              doprocessMcDplusPiWithDmesMlAndDecayTypeCheck, doprocessMcDplusPiWithB0Ml,
                                              doprocessMcDplusPiWithB0MlAndDecayTypeCheck};
    const AxisSpec axisMass = ((std::accumulate(processFuncDplusPi.begin(), processFuncDplusPi.end(), 0)) > 0) ? axisMassDminus : axisMassDeltaMassDStar;
    std::string dMesSpecie;
    if ((std::accumulate(processFuncDplusPi.begin(), processFuncDplusPi.end(), 0)) > 0) {
      dMesSpecie += "D^{#minus}";
    } else {
      dMesSpecie += "D^{0}#pi^{#minus}";
    }

    if (doprocessDataDplusPi || doprocessDataDplusPiWithDmesMl || doprocessDataDplusPiWithB0Ml || doprocessDataDstarPi || doprocessDataDstarPiWithDmesMl) {
      if (fillHistograms) {
        registry.add("hMass", "B^{0} candidates;#it{p}_{T}(B^{0}) (GeV/#it{c});#it{M} (D#pi) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {axisPtB0, axisMassB0}});
        registry.add("hDecLength", "B^{0} candidates;#it{p}_{T}(B^{0}) (GeV/#it{c});B^{0} candidate decay length (cm);entries", {HistType::kTH2F, {axisPtB0, axisDecayLength}});
        registry.add("hDecLengthXy", "B^{0} candidates;#it{p}_{T}(B^{0}) (GeV/#it{c});decay length XY (cm);entries", {HistType::kTH2F, {axisPtB0, axisDecayLength}});
        registry.add("hNormDecLengthXy", "B^{0} candidates;#it{p}_{T}(B^{0}) (GeV/#it{c});B^{0} candidate norm. decay length XY (cm);entries", {HistType::kTH2F, {axisPtB0, axisNormDecayLength}});
        registry.add("hDcaProng0", Form("B^{0} candidates;#it{p}_{T}(B^{0}) (GeV/#it{c});prong 0 (%s) DCAxy to prim. vertex (cm);entries", dMesSpecie.c_str()), {HistType::kTH2F, {axisPtB0, axisDca}});
        registry.add("hDcaProng1", "B^{0} candidates;#it{p}_{T}(B^{0}) (GeV/#it{c});prong 1 (#pi^{#plus}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {axisPtB0, axisDca}});
        registry.add("hCosp", "B^{0} candidates;#it{p}_{T}(B^{0}) (GeV/#it{c});B^{0} candidate cos(#vartheta_{P});entries", {HistType::kTH2F, {axisPtB0, axisCosp}});
        registry.add("hCospXy", "B^{0} candidates;#it{p}_{T}(B^{0}) (GeV/#it{c});B^{0} candidate cos(#vartheta_{P}^{XY});entries", {HistType::kTH2F, {axisPtB0, axisCosp}});
        registry.add("hEta", "B^{0} candidates;#it{p}_{T}(B^{0}) (GeV/#it{c});B^{0} candidate #it{#eta};entries", {HistType::kTH2F, {axisPtB0, axisEta}});
        registry.add("hRapidity", "B^{0} candidates;#it{p}_{T}(B^{0}) (GeV/#it{c});B^{0} candidate #it{y};entries", {HistType::kTH2F, {axisPtB0, axisEta}});
        registry.add("hInvMassD", Form("B^{0} candidates;#it{p}_{T}(%s) (GeV/#it{c});prong0, #it{M}(K#pi) (GeV/#it{c}^{2});entries", dMesSpecie.c_str()), {HistType::kTH2F, {axisPtDminus, axisMass}});
        registry.add("hDecLengthD", Form("B^{0} candidates;#it{p}_{T}(%s) (GeV/#it{c});%s candidate decay length (cm);entries", dMesSpecie.c_str(), dMesSpecie.c_str()), {HistType::kTH2F, {axisPtDminus, axisDecayLength}});
        registry.add("hDecLengthXyD", Form("B^{0} candidates;#it{p}_{T}(%s) (GeV/#it{c});decay length XY (cm);entries", dMesSpecie.c_str()), {HistType::kTH2F, {axisPtDminus, axisDecayLength}});
        registry.add("hCospD", Form("B^{0} candidates;#it{p}_{T}(%s) (GeV/#it{c});%s candidate cos(#vartheta_{P});entries", dMesSpecie.c_str(), dMesSpecie.c_str()), {HistType::kTH2F, {axisPtDminus, axisCosp}});
        registry.add("hCospXyD", Form("B^{0} candidates;#it{p}_{T}(%s) (GeV/#it{c});%s candidate cos(#vartheta_{P}^{XY});entries", dMesSpecie.c_str(), dMesSpecie.c_str()), {HistType::kTH2F, {axisPtDminus, axisCosp}});

        if ((std::accumulate(processFuncDplusPi.begin(), processFuncDplusPi.end(), 0)) > 0) {
          registry.add("hImpParProd", "B^{0} candidates;#it{p}_{T}(B^{0}) (GeV/#it{c});B^{0} candidate impact parameter product;entries", {HistType::kTH2F, {axisPtB0, axisImpParProd}});
          registry.add("hPtProng0", Form("B^{0} candidates;#it{p}_{T}(B^{0}) (GeV/#it{c});#it{p}_{T}(%s) (GeV/#it{c});entries", dMesSpecie.c_str()), {HistType::kTH2F, {axisPtB0, axisPtDminus}});
          registry.add("hPtProng1", "B^{0} candidates;#it{p}_{T}(B^{0}) (GeV/#it{c});#it{p}_{T}(#pi^{#plus}) (GeV/#it{c});entries", {HistType::kTH2F, {axisPtB0, axisPtPi}});
        } else {
          registry.add("hImpParProngSqSum", "B^{0} candidates;#it{p}_{T}(B^{0}) (GeV/#it{c});B^{0} candidate impact parameter sum;entries", {HistType::kTH2F, {axisPtB0, axisImpParProngSqSum}});
          registry.add("hPtProngD0", Form("B^{0} candidates;#it{p}_{T}(B^{0}) (GeV/#it{c});#it{p}_{T}(%s) (GeV/#it{c});entries", dMesSpecie.c_str()), {HistType::kTH2F, {axisPtB0, axisPtDminus}});
          registry.add("hPtProngSoftPi", "B^{0} candidates;#it{p}_{T}(B^{0}) (GeV/#it{c});#it{p}_{T}(#pi^{#plus}) (GeV/#it{c});entries", {HistType::kTH2F, {axisPtB0, axisPtSoftPi}});
          registry.add("hPtProngBachPi", "B^{0} candidates;#it{p}_{T}(B^{0}) (GeV/#it{c});#it{p}_{T}(#pi^{#plus}) (GeV/#it{c});entries", {HistType::kTH2F, {axisPtB0, axisPtPi}});
          registry.add("hDcaProng2", "B^{0} candidates;#it{p}_{T}(B^{0}) (GeV/#it{c});prong 2 (#pi^{#plus}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {axisPtB0, axisDca}});
        }
        // ML scores of D- daughter
        if (doprocessDataDplusPiWithDmesMl || doprocessDataDstarPiWithDmesMl) {
          registry.add("hMlScoreBkgD", Form("B^{0} candidates;#it{p}_{T}(%s) (GeV/#it{c});prong0, %s ML background score;entries", dMesSpecie.c_str(), dMesSpecie.c_str()), {HistType::kTH2F, {axisPtDminus, axisMlScore}});
          registry.add("hMlScorePromptD", Form("B^{0} candidates;#it{p}_{T}(%s) (GeV/#it{c});prong0, %s ML prompt score;entries", dMesSpecie.c_str(), dMesSpecie.c_str()), {HistType::kTH2F, {axisPtDminus, axisMlScore}});
          registry.add("hMlScoreNonPromptD", Form("B^{0} candidates;#it{p}_{T}(%s) (GeV/#it{c});prong0, %s ML nonprompt score;entries", dMesSpecie.c_str(), dMesSpecie.c_str()), {HistType::kTH2F, {axisPtDminus, axisMlScore}});
        }

        // ML scores of B0 candidate
        if (doprocessDataDplusPiWithB0Ml) {
          registry.add("hMlScoreSigB0", "B^{0} candidates;#it{p}_{T}(B^{0}) (GeV/#it{c});prong0, B^{0} ML signal score;entries", {HistType::kTH2F, {axisPtB0, axisMlScore}});
        }
      }
      if (fillSparses) {
        if (!(doprocessDataDplusPiWithDmesMl || doprocessDataDplusPiWithB0Ml || doprocessDataDstarPiWithDmesMl)) {
          if ((std::accumulate(processFuncDplusPi.begin(), processFuncDplusPi.end(), 0)) > 0) {
            registry.add("hMassPtCutVars", "B^{0} candidates;#it{M} (D#pi) (GeV/#it{c}^{2});#it{p}_{T}(B^{0}) (GeV/#it{c});B^{0} candidate decay length (cm);B^{0} candidate norm. decay length XY (cm);B^{0} candidate impact parameter product (cm);B^{0} candidate cos(#vartheta_{P});#it{M} (K#pi) (GeV/#it{c}^{2});#it{p}_{T}(%s) (GeV/#it{c});%s candidate decay length (cm);%s candidate cos(#vartheta_{P})", {HistType::kTHnSparseF, {axisMassB0, axisPtB0, axisDecayLength, axisNormDecayLength, axisImpParProd, axisCosp, axisMass, axisPtDminus, axisDecayLength, axisCosp}});
          } else {
            registry.add("hMassPtCutVars", "B^{0} candidates;#it{M} (D#pi) (GeV/#it{c}^{2});#it{p}_{T}(B^{0}) (GeV/#it{c});B^{0} candidate decay length (cm);B^{0} candidate norm. decay length XY (cm);B^{0} candidate impact parameter product (cm);B^{0} candidate cos(#vartheta_{P});#it{M} (K#pi) (GeV/#it{c}^{2});#it{p}_{T}(%s) (GeV/#it{c});%s candidate decay length (cm);%s candidate cos(#vartheta_{P})", {HistType::kTHnSparseF, {axisMassB0, axisPtB0, axisDecayLength, axisNormDecayLength, axisImpParProngSqSum, axisCosp, axisMass, axisPtDminus, axisDecayLength, axisCosp}});
          }
        } else {
          if ((std::accumulate(processFuncDplusPi.begin(), processFuncDplusPi.end(), 0)) > 0) {
            registry.add("hMassPtCutVars", "B^{0} candidates;#it{M} (D#pi) (GeV/#it{c}^{2});#it{p}_{T}(B^{0}) (GeV/#it{c});B^{0} candidate decay length (cm);B^{0} candidate norm. decay length XY (cm);B^{0} candidate impact parameter product (cm);B^{0} candidate cos(#vartheta_{P});#it{M} (K#pi) (GeV/#it{c}^{2});#it{p}_{T}(%s) (GeV/#it{c});%s candidate ML score bkg;%s candidate ML score nonprompt", {HistType::kTHnSparseF, {axisMassB0, axisPtB0, axisDecayLength, axisNormDecayLength, axisImpParProd, axisCosp, axisMass, axisPtDminus, axisMlScore, axisMlScore}});
          } else {
            registry.add("hMassPtCutVars", "B^{0} candidates;#it{M} (D#pi) (GeV/#it{c}^{2});#it{p}_{T}(B^{0}) (GeV/#it{c});B^{0} candidate decay length (cm);B^{0} candidate norm. decay length XY (cm);B^{0} candidate impact parameter product (cm);B^{0} candidate cos(#vartheta_{P});#it{M} (K#pi) (GeV/#it{c}^{2});#it{p}_{T}(%s) (GeV/#it{c});%s candidate ML score bkg;%s candidate ML score nonprompt", {HistType::kTHnSparseF, {axisMassB0, axisPtB0, axisDecayLength, axisNormDecayLength, axisImpParProngSqSum, axisCosp, axisMass, axisPtDminus, axisMlScore, axisMlScore}});
          }
        }
      }
    }

    if (doprocessMcDplusPi || doprocessMcDplusPiWithDecayTypeCheck || doprocessMcDplusPiWithDmesMl || doprocessMcDplusPiWithDmesMlAndDecayTypeCheck || doprocessMcDplusPiWithB0Ml || doprocessMcDplusPiWithB0MlAndDecayTypeCheck ||
        doprocessMcDstarPi || doprocessMcDstarPiWithDmesMl) {
      if (fillHistograms) {
        // gen histos
        registry.add("hEtaGen", "B^{0} particles (generated);#it{p}_{T}^{gen}(B^{0}) (GeV/#it{c});#it{#eta}^{gen}(B^{0});entries", {HistType::kTH2F, {axisPtB0, axisEta}});
        registry.add("hYGen", "B^{0} particles (generated);#it{p}_{T}^{gen}(B^{0}) (GeV/#it{c});#it{y}^{gen}(B^{0});entries", {HistType::kTH2F, {axisPtB0, axisEta}});
        registry.add("hYGenWithProngsInAcceptance", "MC particles (generated-daughters in acceptance);#it{p}_{T}^{gen}(B^{0}) (GeV/#it{c});#it{y}^{gen}(B^{0});entries", {HistType::kTH2F, {axisPtB0, axisEta}});
        registry.add("hPtProng0Gen", Form("B^{0} particles (generated);#it{p}_{T}^{gen}(B^{0}) (GeV/#it{c});#it{p}_{T}^{gen}(%s) (GeV/#it{c});entries", dMesSpecie.c_str()), {HistType::kTH2F, {axisPtB0, axisPtDminus}});
        registry.add("hPtProng1Gen", "B^{0} particles (generated);#it{p}_{T}^{gen}(B^{0}) (GeV/#it{c});#it{p}_{T}^{gen}(#pi^{#plus}) (GeV/#it{c});entries", {HistType::kTH2F, {axisPtB0, axisPtPi}});
        registry.add("hYProng0Gen", Form("B^{0} particles (generated);#it{p}_{T}^{gen}(B^{0}) (GeV/#it{c});#it{y}^{gen}(%s);entries", dMesSpecie.c_str()), {HistType::kTH2F, {axisPtB0, axisEta}});
        registry.add("hYProng1Gen", "B^{0} particles (generated);#it{p}_{T}^{gen}(B^{0}) (GeV/#it{c});#it{y}^{gen}(#pi^{#plus});entries", {HistType::kTH2F, {axisPtB0, axisEta}});
        registry.add("hEtaProng0Gen", Form("B^{0} particles (generated);#it{p}_{T}^{gen}(B^{0}) (GeV/#it{c});#it{#eta}^{gen}(%s);entries", dMesSpecie.c_str()), {HistType::kTH2F, {axisPtB0, axisEta}});
        registry.add("hEtaProng1Gen", "B^{0} particles (generated);#it{p}_{T}^{gen}(B^{0}) (GeV/#it{c});#it{#eta}^{gen}(#pi^{#plus});entries", {HistType::kTH2F, {axisPtB0, axisEta}});

        // reco histos
        // signal
        registry.add("hMassRecSig", "B^{0} candidates (matched);#it{p}_{T}(B^{0}) (GeV/#it{c});#it{M} (D#pi) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {axisPtB0, axisMassB0}});
        registry.add("hDecLengthRecSig", "B^{0} candidates (matched);#it{p}_{T}(B^{0}) (GeV/#it{c});B^{0} candidate decay length (cm);entries", {HistType::kTH2F, {axisPtB0, axisDecayLength}});
        registry.add("hDecLengthXyRecSig", "B^{0} candidates (matched);#it{p}_{T}(B^{0}) (GeV/#it{c});decay length XY (cm);entries", {HistType::kTH2F, {axisPtB0, axisDecayLength}});
        registry.add("hNormDecLengthXyRecSig", "B^{0} candidates (matched);#it{p}_{T}(B^{0}) (GeV/#it{c});B^{0} candidate norm. decay length XY (cm);entries", {HistType::kTH2F, {axisPtB0, axisNormDecayLength}});
        registry.add("hDcaProng0RecSig", Form("B^{0} candidates (matched);#it{p}_{T}(B^{0}) (GeV/#it{c});prong 0 (%s) DCAxy to prim. vertex (cm);entries", dMesSpecie.c_str()), {HistType::kTH2F, {axisPtB0, axisDca}});
        registry.add("hDcaProng1RecSig", "B^{0} candidates (matched);#it{p}_{T}(B^{0}) (GeV/#it{c});prong 1 (#pi^{#plus}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {axisPtB0, axisDca}});
        registry.add("hCospRecSig", "B^{0} candidates (matched);#it{p}_{T}(B^{0}) (GeV/#it{c});B^{0} candidate cos(#vartheta_{P});entries", {HistType::kTH2F, {axisPtB0, axisCosp}});
        registry.add("hCospXyRecSig", "B^{0} candidates (matched);#it{p}_{T}(B^{0}) (GeV/#it{c});B^{0} candidate cos(#vartheta_{P}^{XY});entries", {HistType::kTH2F, {axisPtB0, axisCosp}});
        registry.add("hEtaRecSig", "B^{0} candidates (matched);#it{p}_{T}(B^{0}) (GeV/#it{c});B^{0} candidate #it{#eta};entries", {HistType::kTH2F, {axisPtB0, axisEta}});
        registry.add("hRapidityRecSig", "B^{0} candidates (matched);#it{p}_{T}(B^{0}) (GeV/#it{c});B^{0} candidate #it{y};entries", {HistType::kTH2F, {axisPtB0, axisEta}});
        registry.add("hInvMassDRecSig", Form("B^{0} candidates (matched);#it{p}_{T}(%s) (GeV/#it{c});prong0, #it{M}(K#pi) (GeV/#it{c}^{2});entries", dMesSpecie.c_str()), {HistType::kTH2F, {axisPtDminus, axisMass}});
        registry.add("hDecLengthDRecSig", Form("B^{0} candidates (matched);#it{p}_{T}(%s) (GeV/#it{c});%s candidate decay length (cm);entries", dMesSpecie.c_str(), dMesSpecie.c_str()), {HistType::kTH2F, {axisPtDminus, axisDecayLength}});
        registry.add("hDecLengthXyDRecSig", Form("B^{0} candidates (matched);#it{p}_{T}(%s) (GeV/#it{c});decay length XY (cm);entries", dMesSpecie.c_str()), {HistType::kTH2F, {axisPtDminus, axisDecayLength}});
        registry.add("hCospDRecSig", Form("B^{0} candidates (matched);#it{p}_{T}(%s) (GeV/#it{c});%s candidate cos(#vartheta_{P});entries", dMesSpecie.c_str(), dMesSpecie.c_str()), {HistType::kTH2F, {axisPtDminus, axisCosp}});
        registry.add("hCospXyDRecSig", Form("B^{0} candidates (matched);#it{p}_{T}(%s) (GeV/#it{c});%s candidate cos(#vartheta_{P}^{XY});entries", dMesSpecie.c_str(), dMesSpecie.c_str()), {HistType::kTH2F, {axisPtDminus, axisCosp}});

        if ((std::accumulate(processFuncDplusPi.begin(), processFuncDplusPi.end(), 0)) > 0) {
          registry.add("hPtProng0RecSig", Form("B^{0} candidates (matched);#it{p}_{T}(B^{0}) (GeV/#it{c});#it{p}_{T}(%s) (GeV/#it{c});entries", dMesSpecie.c_str()), {HistType::kTH2F, {axisPtB0, axisPtDminus}});
          registry.add("hPtProng1RecSig", "B^{0} candidates (matched);#it{p}_{T}(B^{0}) (GeV/#it{c});#it{p}_{T}(#pi^{#plus}) (GeV/#it{c});entries", {HistType::kTH2F, {axisPtB0, axisPtPi}});
          registry.add("hImpParProdRecSig", "B^{0} candidates (matched);#it{p}_{T}(B^{0}) (GeV/#it{c});B^{0} candidate impact parameter product;entries", {HistType::kTH2F, {axisPtB0, axisImpParProd}});
        } else {
          registry.add("hPtProngD0RecSig", Form("B^{0} candidates (matched);#it{p}_{T}(B^{0}) (GeV/#it{c});#it{p}_{T}(%s) (GeV/#it{c});entries", dMesSpecie.c_str()), {HistType::kTH2F, {axisPtB0, axisPtDminus}});
          registry.add("hPtProngSoftPiRecSig", "B^{0} candidates (matched);#it{p}_{T}(B^{0}) (GeV/#it{c});#it{p}_{T}(#pi^{#plus}) (GeV/#it{c});entries", {HistType::kTH2F, {axisPtB0, axisPtSoftPi}});
          registry.add("hPtProngBachPiRecSig", "B^{0} candidates (matched);#it{p}_{T}(B^{0}) (GeV/#it{c});#it{p}_{T}(#pi^{#plus}) (GeV/#it{c});entries", {HistType::kTH2F, {axisPtB0, axisPtPi}});
          registry.add("hDcaProng2RecSig", "B^{0} candidates (matched);#it{p}_{T}(B^{0}) (GeV/#it{c});prong 2 (#pi^{#minus}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {axisPtB0, axisDca}});
          registry.add("hImpParProngSqSumRecSig", "B^{0} candidates (matched);#it{p}_{T}(B^{0}) (GeV/#it{c});B^{0} candidate impact parameter product;entries", {HistType::kTH2F, {axisPtB0, axisImpParProngSqSum}});
        }

        // background
        if (fillBackground) {
          registry.add("hMassRecBg", "B^{0} candidates (unmatched);#it{p}_{T}(B^{0}) (GeV/#it{c});#it{M} (D#pi) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {axisPtB0, axisMassB0}});
          registry.add("hDecLengthRecBg", "B^{0} candidates (unmatched);#it{p}_{T}(B^{0}) (GeV/#it{c});B^{0} candidate decay length (cm);entries", {HistType::kTH2F, {axisPtB0, axisDecayLength}});
          registry.add("hDecLengthXyRecBg", "B^{0} candidates (unmatched);#it{p}_{T}(B^{0}) (GeV/#it{c});decay length XY (cm);entries", {HistType::kTH2F, {axisPtB0, axisDecayLength}});
          registry.add("hNormDecLengthXyRecBg", "B^{0} candidates (unmatched);#it{p}_{T}(B^{0}) (GeV/#it{c});B^{0} candidate norm. decay length XY (cm);entries", {HistType::kTH2F, {axisPtB0, axisNormDecayLength}});
          registry.add("hDcaProng0RecBg", Form("B^{0} candidates (unmatched);#it{p}_{T}(B^{0}) (GeV/#it{c});prong 0 (%s) DCAxy to prim. vertex (cm);entries", dMesSpecie.c_str()), {HistType::kTH2F, {axisPtB0, axisDca}});
          registry.add("hDcaProng1RecBg", "B^{0} candidates (unmatched);#it{p}_{T}(B^{0}) (GeV/#it{c});prong 1 (#pi^{#plus}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {axisPtB0, axisDca}});
          registry.add("hCospRecBg", "B^{0} candidates (unmatched);#it{p}_{T}(B^{0}) (GeV/#it{c});B^{0} candidate cos(#vartheta_{P});entries", {HistType::kTH2F, {axisPtB0, axisCosp}});
          registry.add("hCospXyRecBg", "B^{0} candidates (unmatched);#it{p}_{T}(B^{0}) (GeV/#it{c});B^{0} candidate cos(#vartheta_{P}^{XY});entries", {HistType::kTH2F, {axisPtB0, axisCosp}});
          registry.add("hEtaRecBg", "B^{0} candidates (unmatched);#it{p}_{T}(B^{0}) (GeV/#it{c});B^{0} candidate #it{#eta};entries", {HistType::kTH2F, {axisPtB0, axisEta}});
          registry.add("hRapidityRecBg", "B^{0} candidates (unmatched);#it{p}_{T}(B^{0}) (GeV/#it{c});B^{0} candidate #it{y};entries", {HistType::kTH2F, {axisPtB0, axisEta}});
          registry.add("hInvMassDRecBg", Form("B^{0} candidates (unmatched);#it{p}_{T}(%s) (GeV/#it{c});prong0, #it{M}(K#pi) (GeV/#it{c}^{2});entries", dMesSpecie.c_str()), {HistType::kTH2F, {axisPtDminus, axisMass}});
          registry.add("hDecLengthDRecBg", Form("B^{0} candidates (unmatched);#it{p}_{T}(%s) (GeV/#it{c});%s candidate decay length (cm);entries", dMesSpecie.c_str(), dMesSpecie.c_str()), {HistType::kTH2F, {axisPtDminus, axisDecayLength}});
          registry.add("hDecLengthXyDRecBg", Form("B^{0} candidates (unmatched);#it{p}_{T}(%s) (GeV/#it{c});decay length XY (cm);entries", dMesSpecie.c_str()), {HistType::kTH2F, {axisPtDminus, axisDecayLength}});
          registry.add("hCospDRecBg", Form("B^{0} candidates (unmatched);#it{p}_{T}(%s) (GeV/#it{c});%s candidate cos(#vartheta_{P});entries", dMesSpecie.c_str(), dMesSpecie.c_str()), {HistType::kTH2F, {axisPtDminus, axisCosp}});
          registry.add("hCospXyDRecBg", Form("B^{0} candidates (unmatched);#it{p}_{T}(%s) (GeV/#it{c});%s candidate cos(#vartheta_{P}^{XY});entries", dMesSpecie.c_str(), dMesSpecie.c_str()), {HistType::kTH2F, {axisPtDminus, axisCosp}});

          if ((std::accumulate(processFuncDplusPi.begin(), processFuncDplusPi.end(), 0)) > 0) {
            registry.add("hPtProng0RecBg", Form("B^{0} candidates (unmatched);#it{p}_{T}(B^{0}) (GeV/#it{c});#it{p}_{T}(%s) (GeV/#it{c});entries", dMesSpecie.c_str()), {HistType::kTH2F, {axisPtB0, axisPtDminus}});
            registry.add("hPtProng1RecBg", "B^{0} candidates (unmatched);#it{p}_{T}(B^{0}) (GeV/#it{c});#it{p}_{T}(#pi^{#plus}) (GeV/#it{c});entries", {HistType::kTH2F, {axisPtB0, axisPtPi}});
            registry.add("hImpParProdRecBg", "B^{0} candidates (unmatched);#it{p}_{T}(B^{0}) (GeV/#it{c});B^{0} candidate impact parameter product;entries", {HistType::kTH2F, {axisPtB0, axisImpParProd}});
          } else {
            registry.add("hPtProngD0RecBg", Form("B^{0} candidates (unmatched);#it{p}_{T}(B^{0}) (GeV/#it{c});#it{p}_{T}(%s) (GeV/#it{c});entries", dMesSpecie.c_str()), {HistType::kTH2F, {axisPtB0, axisPtDminus}});
            registry.add("hPtProngSoftPiRecBg", "B^{0} candidates (unmatched);#it{p}_{T}(B^{0}) (GeV/#it{c});#it{p}_{T}(#pi^{#plus}) (GeV/#it{c});entries", {HistType::kTH2F, {axisPtB0, axisPtSoftPi}});
            registry.add("hPtProngBachPiRecBg", "B^{0} candidates (unmatched);#it{p}_{T}(B^{0}) (GeV/#it{c});#it{p}_{T}(#pi^{#plus}) (GeV/#it{c});entries", {HistType::kTH2F, {axisPtB0, axisPtPi}});
            registry.add("hImpParProngSqSumRecBg", "B^{0} candidates (unmatched);#it{p}_{T}(B^{0}) (GeV/#it{c});B^{0} candidate impact parameter product;entries", {HistType::kTH2F, {axisPtB0, axisImpParProngSqSum}});
            registry.add("hDcaProng2RecBg", "B^{0} candidates (unmatched);#it{p}_{T}(B^{0}) (GeV/#it{c});prong 2 (#pi^{#plus}) DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {axisPtB0, axisDca}});
          }
        }
        // MC checks
        if (doprocessMcDplusPiWithDecayTypeCheck || doprocessMcDplusPiWithB0MlAndDecayTypeCheck || doprocessMcDplusPiWithDmesMlAndDecayTypeCheck) {
          constexpr uint8_t kNBinsDecayTypeMc = hf_cand_b0::DecayTypeMc::NDecayTypeMc;
          TString labels[kNBinsDecayTypeMc];
          labels[hf_cand_b0::DecayTypeMc::B0ToDplusPiToPiKPiPi] = Form("B^{0} #rightarrow (%s #rightarrow #pi^{#minus} K^{#plus} #pi^{#minus}) #pi^{#plus}", dMesSpecie.c_str());
          labels[hf_cand_b0::DecayTypeMc::B0ToDsPiToKKPiPi] = Form("B^{0} #rightarrow (%s_{s} #rightarrow K^{#minus} K^{#plus} #pi^{#minus}) #pi^{#plus}", dMesSpecie.c_str());
          labels[hf_cand_b0::DecayTypeMc::BsToDsPiToKKPiPi] = Form("B_{s}^{0} #rightarrow (%s_{s} #rightarrow K^{#minus} K^{#plus} #pi^{#minus}) #pi^{#plus}", dMesSpecie.c_str());
          labels[hf_cand_b0::DecayTypeMc::B0ToDplusKToPiKPiK] = Form("B^{0} #rightarrow (%s #rightarrow #pi^{#minus} K^{#plus} #pi^{#minus}) K^{#plus}", dMesSpecie.c_str());
          labels[hf_cand_b0::DecayTypeMc::PartlyRecoDecay] = "Partly reconstructed decay channel";
          labels[hf_cand_b0::DecayTypeMc::OtherDecay] = "Other decays";
          static const AxisSpec axisDecayType = {kNBinsDecayTypeMc, 0.5, kNBinsDecayTypeMc + 0.5, ""};
          registry.add("hDecayTypeMc", "DecayType", {HistType::kTH3F, {axisDecayType, axisMassB0, axisPtB0}});
          for (uint8_t iBin = 0; iBin < kNBinsDecayTypeMc; ++iBin) {
            registry.get<TH3>(HIST("hDecayTypeMc"))->GetXaxis()->SetBinLabel(iBin + 1, labels[iBin]);
          }
        }
        // ML scores of D- daughter
        if (doprocessMcDplusPiWithDmesMl || doprocessMcDplusPiWithDmesMlAndDecayTypeCheck || doprocessMcDstarPiWithDmesMl) {
          // signal
          registry.add("hMlScoreBkgDRecSig", Form("B^{0} candidates (matched);#it{p}_{T}(%s) (GeV/#it{c});prong0, %s ML background score;entries", dMesSpecie.c_str(), dMesSpecie.c_str()), {HistType::kTH2F, {axisPtDminus, axisMlScore}});
          registry.add("hMlScorePromptDRecSig", Form("B^{0} candidates (matched);#it{p}_{T}(%s) (GeV/#it{c});prong0, %s ML prompt score;entries", dMesSpecie.c_str(), dMesSpecie.c_str()), {HistType::kTH2F, {axisPtDminus, axisMlScore}});
          registry.add("hMlScoreNonPromptDRecSig", Form("B^{0} candidates (matched);#it{p}_{T}(%s) (GeV/#it{c});prong0, %s ML nonprompt score;entries", dMesSpecie.c_str(), dMesSpecie.c_str()), {HistType::kTH2F, {axisPtDminus, axisMlScore}});
          // background
          registry.add("hMlScoreBkgDRecBg", Form("B^{0} candidates (unmatched);#it{p}_{T}(%s) (GeV/#it{c});prong0, %s ML background score;entries", dMesSpecie.c_str(), dMesSpecie.c_str()), {HistType::kTH2F, {axisPtDminus, axisMlScore}});
          registry.add("hMlScorePromptDRecBg", Form("B^{0} candidates (unmatched);#it{p}_{T}(%s) (GeV/#it{c});prong0, %s ML prompt score;entries", dMesSpecie.c_str(), dMesSpecie.c_str()), {HistType::kTH2F, {axisPtDminus, axisMlScore}});
          registry.add("hMlScoreNonPromptDRecBg", Form("B^{0} candidates (unmatched);#it{p}_{T}(%s) (GeV/#it{c});prong0, %s ML nonprompt score;entries", dMesSpecie.c_str(), dMesSpecie.c_str()), {HistType::kTH2F, {axisPtDminus, axisMlScore}});
        }
        // ML scores of B0 candidate
        if (doprocessMcDplusPiWithB0Ml || doprocessMcDplusPiWithB0MlAndDecayTypeCheck) {
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
        if (!(doprocessDataDplusPiWithDmesMl || doprocessDataDplusPiWithB0Ml || doprocessDataDstarPiWithDmesMl)) {
          registry.add("hMassPtCutVarsRecSig", Form("B^{0} candidates (matched);#it{M} (D#pi) (GeV/#it{c}^{2});#it{p}_{T}(B^{0}) (GeV/#it{c});B^{0} candidate decay length (cm);B^{0} candidate norm. decay length XY (cm);B^{0} candidate impact parameter product (cm);B^{0} candidate cos(#vartheta_{P});#it{M} (K#pi) (GeV/#it{c}^{2});#it{p}_{T}(%s) (GeV/#it{c});%s candidate decay length (cm);%s candidate cos(#vartheta_{P})", dMesSpecie.c_str(), dMesSpecie.c_str(), dMesSpecie.c_str()), {HistType::kTHnSparseF, {axisMassB0, axisPtB0, axisDecayLength, axisNormDecayLength, axisImpParProd, axisCosp, axisMass, axisPtDminus, axisDecayLength, axisCosp}});
          if (fillBackground) {
            registry.add("hMassPtCutVarsRecBg", Form("B^{0} candidates (unmatched);#it{M} (D#pi) (GeV/#it{c}^{2});#it{p}_{T}(B^{0}) (GeV/#it{c});B^{0} candidate decay length (cm);B^{0} candidate norm. decay length XY (cm);B^{0} candidate impact parameter product (cm);B^{0} candidate cos(#vartheta_{P});#it{M} (K#pi) (GeV/#it{c}^{2});#it{p}_{T}(%s) (GeV/#it{c});%s candidate decay length (cm);%s candidate cos(#vartheta_{P})", dMesSpecie.c_str(), dMesSpecie.c_str(), dMesSpecie.c_str()), {HistType::kTHnSparseF, {axisMassB0, axisPtB0, axisDecayLength, axisNormDecayLength, axisImpParProd, axisCosp, axisMass, axisPtDminus, axisDecayLength, axisCosp}});
          }
        } else {
          registry.add("hMassPtCutVarsRecSig", Form("B^{0} candidates (matched);#it{M} (D#pi) (GeV/#it{c}^{2});#it{p}_{T}(B^{0}) (GeV/#it{c});B^{0} candidate decay length (cm);B^{0} candidate norm. decay length XY (cm);B^{0} candidate impact parameter product (cm);B^{0} candidate cos(#vartheta_{P});#it{M} (K#pi) (GeV/#it{c}^{2});#it{p}_{T}(%s) (GeV/#it{c});%s candidate ML score bkg;%s candidate ML score nonprompt", dMesSpecie.c_str(), dMesSpecie.c_str(), dMesSpecie.c_str()), {HistType::kTHnSparseF, {axisMassB0, axisPtB0, axisDecayLength, axisNormDecayLength, axisImpParProd, axisCosp, axisMass, axisPtDminus, axisMlScore, axisMlScore}});
          if (fillBackground) {
            registry.add("hMassPtCutVarsRecBg", Form("B^{0} candidates (unmatched);#it{M} (D#pi) (GeV/#it{c}^{2});#it{p}_{T}(B^{0}) (GeV/#it{c});B^{0} candidate decay length (cm);B^{0} candidate norm. decay length XY (cm);B^{0} candidate impact parameter product (cm);B^{0} candidate cos(#vartheta_{P});#it{M} (K#pi) (GeV/#it{c}^{2});#it{p}_{T}(%s) (GeV/#it{c});%s candidate ML score bkg;%s candidate ML score nonprompt", dMesSpecie.c_str(), dMesSpecie.c_str(), dMesSpecie.c_str()), {HistType::kTHnSparseF, {axisMassB0, axisPtB0, axisDecayLength, axisNormDecayLength, axisImpParProd, axisCosp, axisMass, axisPtDminus, axisMlScore, axisMlScore}});
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
  template <bool DoMc, bool WithDecayTypeCheck, bool WithDmesMl, typename Cand, typename CandsDmes, typename SoftPions>
  void fillCandDStarPi(Cand const& candidate,
                       SoftPions const& softPions,
                       CandsDmes const&)
  {
    auto ptCandB0 = candidate.pt();
    auto invMassB0 = HfHelper::invMassB0ToDPi(candidate);
    auto candD = candidate.template prongD0_as<CandsDmes>();
    auto ptD = candidate.ptProng0();
    auto invMassD = candD.invMassHypo0();
    auto softPi = softPions.rawIteratorAt(candD.globalIndex());
    std::array<float, 3> const posPv{candidate.posX(), candidate.posY(), candidate.posZ()};
    std::array<float, 3> const posSvD{candD.xSecondaryVertex(), candD.ySecondaryVertex(), candD.zSecondaryVertex()};
    std::array<float, 3> const momD{candD.pVector()};
    auto cospD = RecoDecay::cpa(posPv, posSvD, momD);
    auto cospXyD = RecoDecay::cpaXY(posPv, posSvD, momD);
    auto decLenD = RecoDecay::distance(posPv, posSvD);
    auto decLenXyD = RecoDecay::distanceXY(posPv, posSvD);

    int8_t flagMcMatchRec = 0;
    int8_t flagWrongCollision = 0;
    bool isSignal = false;
    if constexpr (DoMc) {
      flagMcMatchRec = candidate.flagMcMatchRec();
      flagWrongCollision = candidate.flagWrongCollision();
      isSignal = TESTBIT(std::abs(flagMcMatchRec), hf_cand_b0::DecayTypeMc::B0ToDstarPiToD0PiPiToKPiPiPi);
    }

    if (fillHistograms) {
      if constexpr (DoMc) {
        if (isSignal) {
          registry.fill(HIST("hMassRecSig"), ptCandB0, HfHelper::invMassB0ToDPi(candidate));
          registry.fill(HIST("hPtProngD0RecSig"), ptCandB0, candidate.ptProng0());
          registry.fill(HIST("hPtProngSoftPiRecSig"), ptCandB0, candidate.ptProng1());
          registry.fill(HIST("hPtProngBachPiRecSig"), ptCandB0, candidate.ptProng2());
          registry.fill(HIST("hImpParProngSqSumRecSig"), ptCandB0, candidate.impactParameterProngSqSum());
          registry.fill(HIST("hDecLengthRecSig"), ptCandB0, candidate.decayLength());
          registry.fill(HIST("hDecLengthXyRecSig"), ptCandB0, candidate.decayLengthXY());
          registry.fill(HIST("hNormDecLengthXyRecSig"), ptCandB0, candidate.decayLengthXY() / candidate.errorDecayLengthXY());
          registry.fill(HIST("hDcaProng0RecSig"), ptCandB0, candidate.impactParameter0());
          registry.fill(HIST("hDcaProng1RecSig"), ptCandB0, candidate.impactParameter1());
          registry.fill(HIST("hDcaProng2RecSig"), ptCandB0, candidate.impactParameter2());
          registry.fill(HIST("hCospRecSig"), ptCandB0, candidate.cpa());
          registry.fill(HIST("hCospXyRecSig"), ptCandB0, candidate.cpaXY());
          registry.fill(HIST("hEtaRecSig"), ptCandB0, candidate.eta());
          registry.fill(HIST("hRapidityRecSig"), ptCandB0, HfHelper::yB0(candidate));
          registry.fill(HIST("hInvMassDRecSig"), ptD, invMassD);
          registry.fill(HIST("hDecLengthDRecSig"), ptD, decLenD);
          registry.fill(HIST("hDecLengthXyDRecSig"), ptD, decLenXyD);
          registry.fill(HIST("hCospDRecSig"), ptD, cospD);
          registry.fill(HIST("hCospXyDRecSig"), ptD, cospXyD);
          if constexpr (WithDmesMl) {
            registry.fill(HIST("hMlScoreBkgDRecSig"), ptD, candidate.prong0MlScoreBkg());
            registry.fill(HIST("hMlScorePromptDRecSig"), ptD, candidate.prong0MlScorePrompt());
            registry.fill(HIST("hMlScoreNonPromptDRecSig"), ptD, candidate.prong0MlScoreNonprompt());
          }
        } else if (fillBackground) {
          registry.fill(HIST("hMassRecBg"), ptCandB0, HfHelper::invMassB0ToDPi(candidate));
          registry.fill(HIST("hPtProngD0RecBg"), ptCandB0, candidate.ptProng0());
          registry.fill(HIST("hPtProngSoftPiRecBg"), ptCandB0, candidate.ptProng1());
          registry.fill(HIST("hPtProngBachPiRecBg"), ptCandB0, candidate.ptProng2());
          registry.fill(HIST("hImpParProngSqSumRecBg"), ptCandB0, candidate.impactParameterProngSqSum());
          registry.fill(HIST("hDecLengthRecBg"), ptCandB0, candidate.decayLength());
          registry.fill(HIST("hDecLengthXyRecBg"), ptCandB0, candidate.decayLengthXY());
          registry.fill(HIST("hNormDecLengthXyRecBg"), ptCandB0, candidate.decayLengthXY() / candidate.errorDecayLengthXY());
          registry.fill(HIST("hDcaProng0RecBg"), ptCandB0, candidate.impactParameter0());
          registry.fill(HIST("hDcaProng1RecBg"), ptCandB0, candidate.impactParameter1());
          registry.fill(HIST("hDcaProng2RecBg"), ptCandB0, candidate.impactParameter2());
          registry.fill(HIST("hCospRecBg"), ptCandB0, candidate.cpa());
          registry.fill(HIST("hCospXyRecBg"), ptCandB0, candidate.cpaXY());
          registry.fill(HIST("hEtaRecBg"), ptCandB0, candidate.eta());
          registry.fill(HIST("hRapidityRecBg"), ptCandB0, HfHelper::yB0(candidate));
          registry.fill(HIST("hInvMassDRecBg"), ptD, invMassD);
          registry.fill(HIST("hDecLengthDRecBg"), ptD, decLenD);
          registry.fill(HIST("hDecLengthXyDRecBg"), ptD, decLenXyD);
          registry.fill(HIST("hCospDRecBg"), ptD, cospD);
          registry.fill(HIST("hCospXyDRecBg"), ptD, cospXyD);
          if constexpr (WithDmesMl) {
            registry.fill(HIST("hMlScoreBkgDRecBg"), ptD, candidate.prong0MlScoreBkg());
            registry.fill(HIST("hMlScorePromptDRecBg"), ptD, candidate.prong0MlScorePrompt());
            registry.fill(HIST("hMlScoreNonPromptDRecBg"), ptD, candidate.prong0MlScoreNonprompt());
          }
        }
      } else {
        registry.fill(HIST("hMass"), ptCandB0, invMassB0);
        registry.fill(HIST("hPtProngD0"), ptCandB0, candidate.ptProng0());
        registry.fill(HIST("hPtProngSoftPi"), ptCandB0, candidate.ptProng1());
        registry.fill(HIST("hPtProngBachPi"), ptCandB0, candidate.ptProng2());
        registry.fill(HIST("hDecLength"), ptCandB0, candidate.decayLength());
        registry.fill(HIST("hDecLengthXy"), ptCandB0, candidate.decayLengthXY());
        registry.fill(HIST("hNormDecLengthXy"), ptCandB0, candidate.decayLengthXY() / candidate.errorDecayLengthXY());
        registry.fill(HIST("hImpParProngSqSum"), ptCandB0, candidate.impactParameterProngSqSum());
        registry.fill(HIST("hDcaProng0"), ptCandB0, candidate.impactParameter0());
        registry.fill(HIST("hDcaProng1"), ptCandB0, candidate.impactParameter1());
        registry.fill(HIST("hDcaProng2"), ptCandB0, candidate.impactParameter2());
        registry.fill(HIST("hCosp"), ptCandB0, candidate.cpa());
        registry.fill(HIST("hCospXy"), ptCandB0, candidate.cpaXY());
        registry.fill(HIST("hEta"), ptCandB0, candidate.eta());
        registry.fill(HIST("hRapidity"), ptCandB0, HfHelper::yB0(candidate));
        registry.fill(HIST("hInvMassD"), ptD, invMassD);
        registry.fill(HIST("hDecLengthD"), ptD, decLenD);
        registry.fill(HIST("hDecLengthXyD"), ptD, decLenXyD);
        registry.fill(HIST("hCospD"), ptD, cospD);
        registry.fill(HIST("hCospXyD"), ptD, cospXyD);

        if constexpr (WithDmesMl) {
          registry.fill(HIST("hMlScoreBkgD"), ptD, candidate.prong0MlScoreBkg());
          registry.fill(HIST("hMlScorePromptD"), ptD, candidate.prong0MlScorePrompt());
          registry.fill(HIST("hMlScoreNonPromptD"), ptD, candidate.prong0MlScoreNonprompt());
        }
      }
    }
    if (fillSparses) {
      if constexpr (DoMc) {
        if (isSignal) {
          if constexpr (WithDmesMl) {
            registry.fill(HIST("hMassPtCutVarsRecSig"), invMassB0, ptCandB0, candidate.decayLength(), candidate.decayLengthXY() / candidate.errorDecayLengthXY(), candidate.impactParameterProngSqSum(), candidate.cpa(), invMassD, ptD, candidate.prong0MlScoreBkg(), candidate.prong0MlScoreNonprompt());
          } else {
            registry.fill(HIST("hMassPtCutVarsRecSig"), invMassB0, ptCandB0, candidate.decayLength(), candidate.decayLengthXY() / candidate.errorDecayLengthXY(), candidate.impactParameterProngSqSum(), candidate.cpa(), invMassD, ptD, decLenD, cospD);
          }
        } else if (fillBackground) {
          if constexpr (WithDmesMl) {
            registry.fill(HIST("hMassPtCutVarsRecBg"), invMassB0, ptCandB0, candidate.decayLength(), candidate.decayLengthXY() / candidate.errorDecayLengthXY(), candidate.impactParameterProngSqSum(), candidate.cpa(), invMassD, ptD, candidate.prong0MlScoreBkg(), candidate.prong0MlScoreNonprompt());
          } else {
            registry.fill(HIST("hMassPtCutVarsRecBg"), invMassB0, ptCandB0, candidate.decayLength(), candidate.decayLengthXY() / candidate.errorDecayLengthXY(), candidate.impactParameterProngSqSum(), candidate.cpa(), invMassD, ptD, decLenD, cospD);
          }
        }
      } else {
        if constexpr (WithDmesMl) {
          registry.fill(HIST("hMassPtCutVars"), invMassB0, ptCandB0, candidate.decayLength(), candidate.decayLengthXY() / candidate.errorDecayLengthXY(), candidate.impactParameterProngSqSum(), candidate.cpa(), invMassD, ptD, candidate.prong0MlScoreBkg(), candidate.prong0MlScoreNonprompt());
        } else {
          registry.fill(HIST("hMassPtCutVars"), invMassB0, ptCandB0, candidate.decayLength(), candidate.decayLengthXY() / candidate.errorDecayLengthXY(), candidate.impactParameterProngSqSum(), candidate.cpa(), invMassD, ptD, decLenD, cospD);
        }
      }
    }
    if (fillTree) {
      float const pseudoRndm = ptD * 1000. - static_cast<int64_t>(ptD * 1000);
      if (flagMcMatchRec != 0 || (((DoMc && fillBackground) || !DoMc) && (ptCandB0 >= ptMaxForDownSample || pseudoRndm < downSampleBkgFactor))) {
        float prong0MlScoreBkg = -1.;
        float prong0MlScorePrompt = -1.;
        float prong0MlScoreNonprompt = -1.;
        float const candidateMlScoreSig = -1;
        if constexpr (WithDmesMl) {
          prong0MlScoreBkg = candidate.prong0MlScoreBkg();
          prong0MlScorePrompt = candidate.prong0MlScorePrompt();
          prong0MlScoreNonprompt = candidate.prong0MlScoreNonprompt();
        }
        auto prongBachPi = candidate.template prongBachPi_as<TracksBachPions>();

        float ptMother = -1.;
        if constexpr (DoMc) {
          ptMother = candidate.ptMother();
        }

        hfRedCandB0Lite(
          // B-meson features
          invMassB0,
          ptCandB0,
          candidate.eta(),
          candidate.phi(),
          HfHelper::yB0(candidate),
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
          softPi.tpcNSigmaPiSoftPi(),
          softPi.tofNSigmaPiSoftPi(),
          softPi.tpcTofNSigmaPiSoftPi(),
          prong0MlScoreBkg,
          prong0MlScorePrompt,
          prong0MlScoreNonprompt,
          // pion features
          candidate.ptProng2(),
          std::abs(RecoDecay::eta(prongBachPi.pVector())),
          prongBachPi.itsNCls(),
          prongBachPi.tpcNClsCrossedRows(),
          prongBachPi.tpcChi2NCl(),
          candidate.impactParameter2(),
          prongBachPi.tpcNSigmaPi(),
          prongBachPi.tofNSigmaPi(),
          prongBachPi.tpcTofNSigmaPi(),
          // MC truth
          flagMcMatchRec,
          isSignal,
          flagWrongCollision,
          ptMother);

        if constexpr (WithDecayTypeCheck) {
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

  /// Fill candidate information at reconstruction level
  /// \param doMc is the flag to enable the filling with MC information
  /// \param withDecayTypeCheck is the flag to enable MC with decay type check
  /// \param withDmesMl is the flag to enable the filling with ML scores for the D- daughter
  /// \param withB0Ml is the flag to enable the filling with ML scores for the B0 candidate
  /// \param candidate is the B0 candidate
  /// \param candidatesD is the table with D- candidates
  template <bool DoMc, bool WithDecayTypeCheck, bool WithDmesMl, bool WithB0Ml, typename Cand, typename CandsDmes>
  void fillCand(Cand const& candidate,
                CandsDmes const&)
  {
    auto ptCandB0 = candidate.pt();
    auto invMassB0 = HfHelper::invMassB0ToDPi(candidate);
    auto candD = candidate.template prong0_as<CandsDmes>();
    auto ptD = candidate.ptProng0();
    auto invMassD = candD.invMassHypo0();
    std::array<float, 3> const posPv{candidate.posX(), candidate.posY(), candidate.posZ()};
    std::array<float, 3> const posSvD{candD.xSecondaryVertex(), candD.ySecondaryVertex(), candD.zSecondaryVertex()};
    std::array<float, 3> const momD{candD.pVector()};
    auto cospD = RecoDecay::cpa(posPv, posSvD, momD);
    auto cospXyD = RecoDecay::cpaXY(posPv, posSvD, momD);
    auto decLenD = RecoDecay::distance(posPv, posSvD);
    auto decLenXyD = RecoDecay::distanceXY(posPv, posSvD);

    int8_t flagMcMatchRec = 0;
    int8_t flagWrongCollision = 0;
    bool isSignal = false;
    if constexpr (DoMc) {
      flagMcMatchRec = candidate.flagMcMatchRec();
      flagWrongCollision = candidate.flagWrongCollision();
      isSignal = TESTBIT(std::abs(flagMcMatchRec), hf_cand_b0::DecayTypeMc::B0ToDplusPiToPiKPiPi);
    }

    if (fillHistograms) {
      if constexpr (DoMc) {
        if (isSignal) {
          registry.fill(HIST("hMassRecSig"), ptCandB0, HfHelper::invMassB0ToDPi(candidate));
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
          registry.fill(HIST("hRapidityRecSig"), ptCandB0, HfHelper::yB0(candidate));
          registry.fill(HIST("hInvMassDRecSig"), ptD, invMassD);
          registry.fill(HIST("hDecLengthDRecSig"), ptD, decLenD);
          registry.fill(HIST("hDecLengthXyDRecSig"), ptD, decLenXyD);
          registry.fill(HIST("hCospDRecSig"), ptD, cospD);
          registry.fill(HIST("hCospXyDRecSig"), ptD, cospXyD);
          if constexpr (WithDecayTypeCheck) {
            registry.fill(HIST("hDecayTypeMc"), 1 + hf_cand_b0::DecayTypeMc::B0ToDplusPiToPiKPiPi, invMassB0, ptCandB0);
          }
          if constexpr (WithDmesMl) {
            registry.fill(HIST("hMlScoreBkgDRecSig"), ptD, candidate.prong0MlScoreBkg());
            registry.fill(HIST("hMlScorePromptDRecSig"), ptD, candidate.prong0MlScorePrompt());
            registry.fill(HIST("hMlScoreNonPromptDRecSig"), ptD, candidate.prong0MlScoreNonprompt());
          }
          if constexpr (WithB0Ml) {
            registry.fill(HIST("hMlScoreSigB0RecSig"), ptCandB0, candidate.mlProbB0ToDPi());
          }
        } else if (fillBackground) {
          registry.fill(HIST("hMassRecBg"), ptCandB0, HfHelper::invMassB0ToDPi(candidate));
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
          registry.fill(HIST("hRapidityRecBg"), ptCandB0, HfHelper::yB0(candidate));
          registry.fill(HIST("hInvMassDRecBg"), ptD, invMassD);
          registry.fill(HIST("hDecLengthDRecBg"), ptD, decLenD);
          registry.fill(HIST("hDecLengthXyDRecBg"), ptD, decLenXyD);
          registry.fill(HIST("hCospDRecBg"), ptD, cospD);
          registry.fill(HIST("hCospXyDRecBg"), ptD, cospXyD);
          if constexpr (WithDmesMl) {
            registry.fill(HIST("hMlScoreBkgDRecBg"), ptD, candidate.prong0MlScoreBkg());
            registry.fill(HIST("hMlScorePromptDRecBg"), ptD, candidate.prong0MlScorePrompt());
            registry.fill(HIST("hMlScoreNonPromptDRecBg"), ptD, candidate.prong0MlScoreNonprompt());
          }
          if constexpr (WithB0Ml) {
            registry.fill(HIST("hMlScoreSigB0RecBg"), ptCandB0, candidate.mlProbB0ToDPi());
          }
        } else if constexpr (WithDecayTypeCheck) {
          if (TESTBIT(flagMcMatchRec, hf_cand_b0::DecayTypeMc::B0ToDsPiToKKPiPi)) { // B0 → Ds- π+ → (K- K+ π-) π+
            registry.fill(HIST("hDecayTypeMc"), 1 + hf_cand_b0::DecayTypeMc::B0ToDsPiToKKPiPi, invMassB0, ptCandB0);
          } else if (TESTBIT(flagMcMatchRec, hf_cand_b0::DecayTypeMc::BsToDsPiToKKPiPi)) { // B0s → Ds- π+ → (K- K+ π-) π+
            registry.fill(HIST("hDecayTypeMc"), 1 + hf_cand_b0::DecayTypeMc::BsToDsPiToKKPiPi, invMassB0, ptCandB0);
          } else if (TESTBIT(flagMcMatchRec, hf_cand_b0::DecayTypeMc::B0ToDplusKToPiKPiK)) { // B0 → D- K+ → (π- K+ π-) K+
            registry.fill(HIST("hDecayTypeMc"), 1 + hf_cand_b0::DecayTypeMc::B0ToDplusKToPiKPiK, invMassB0, ptCandB0);
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
        registry.fill(HIST("hRapidity"), ptCandB0, HfHelper::yB0(candidate));
        registry.fill(HIST("hInvMassD"), ptD, invMassD);
        registry.fill(HIST("hDecLengthD"), ptD, decLenD);
        registry.fill(HIST("hDecLengthXyD"), ptD, decLenXyD);
        registry.fill(HIST("hCospD"), ptD, cospD);
        registry.fill(HIST("hCospXyD"), ptD, cospXyD);

        if constexpr (WithDmesMl) {
          registry.fill(HIST("hMlScoreBkgD"), ptD, candidate.prong0MlScoreBkg());
          registry.fill(HIST("hMlScorePromptD"), ptD, candidate.prong0MlScorePrompt());
          registry.fill(HIST("hMlScoreNonPromptD"), ptD, candidate.prong0MlScoreNonprompt());
        }
        if constexpr (WithB0Ml) {
          registry.fill(HIST("hMlScoreSigB0"), ptCandB0, candidate.mlProbB0ToDPi());
        }
      }
    }
    if (fillSparses) {
      if constexpr (DoMc) {
        if (isSignal) {
          if constexpr (WithDmesMl) {
            registry.fill(HIST("hMassPtCutVarsRecSig"), invMassB0, ptCandB0, candidate.decayLength(), candidate.decayLengthXY() / candidate.errorDecayLengthXY(), candidate.impactParameterProduct(), candidate.cpa(), invMassD, ptD, candidate.prong0MlScoreBkg(), candidate.prong0MlScoreNonprompt());
          } else {
            registry.fill(HIST("hMassPtCutVarsRecSig"), invMassB0, ptCandB0, candidate.decayLength(), candidate.decayLengthXY() / candidate.errorDecayLengthXY(), candidate.impactParameterProduct(), candidate.cpa(), invMassD, ptD, decLenD, cospD);
          }
        } else if (fillBackground) {
          if constexpr (WithDmesMl) {
            registry.fill(HIST("hMassPtCutVarsRecBg"), invMassB0, ptCandB0, candidate.decayLength(), candidate.decayLengthXY() / candidate.errorDecayLengthXY(), candidate.impactParameterProduct(), candidate.cpa(), invMassD, ptD, candidate.prong0MlScoreBkg(), candidate.prong0MlScoreNonprompt());
          } else {
            registry.fill(HIST("hMassPtCutVarsRecBg"), invMassB0, ptCandB0, candidate.decayLength(), candidate.decayLengthXY() / candidate.errorDecayLengthXY(), candidate.impactParameterProduct(), candidate.cpa(), invMassD, ptD, decLenD, cospD);
          }
        }
      } else {
        if constexpr (WithDmesMl) {
          registry.fill(HIST("hMassPtCutVars"), invMassB0, ptCandB0, candidate.decayLength(), candidate.decayLengthXY() / candidate.errorDecayLengthXY(), candidate.impactParameterProduct(), candidate.cpa(), invMassD, ptD, candidate.prong0MlScoreBkg(), candidate.prong0MlScoreNonprompt());
        } else {
          registry.fill(HIST("hMassPtCutVars"), invMassB0, ptCandB0, candidate.decayLength(), candidate.decayLengthXY() / candidate.errorDecayLengthXY(), candidate.impactParameterProduct(), candidate.cpa(), invMassD, ptD, decLenD, cospD);
        }
      }
    }
    if (fillTree) {
      float const pseudoRndm = ptD * 1000. - static_cast<int64_t>(ptD * 1000);
      if (flagMcMatchRec != 0 || (((DoMc && fillBackground) || !DoMc) && (ptCandB0 >= ptMaxForDownSample || pseudoRndm < downSampleBkgFactor))) {
        float prong0MlScoreBkg = -1.;
        float prong0MlScorePrompt = -1.;
        float prong0MlScoreNonprompt = -1.;
        float candidateMlScoreSig = -1;
        if constexpr (WithDmesMl) {
          prong0MlScoreBkg = candidate.prong0MlScoreBkg();
          prong0MlScorePrompt = candidate.prong0MlScorePrompt();
          prong0MlScoreNonprompt = candidate.prong0MlScoreNonprompt();
        }
        if constexpr (WithB0Ml) {
          candidateMlScoreSig = candidate.mlProbB0ToDPi();
        }
        auto prong1 = candidate.template prong1_as<TracksBachPions>();

        float ptMother = -1.;
        if constexpr (DoMc) {
          ptMother = candidate.ptMother();
        }

        hfRedCandB0Lite(
          // B-meson features
          invMassB0,
          ptCandB0,
          candidate.eta(),
          candidate.phi(),
          HfHelper::yB0(candidate),
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

        if constexpr (WithDecayTypeCheck) {
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
  void processDataDplusPi(soa::Filtered<soa::Join<aod::HfRedCandB0, aod::HfSelB0ToDPi>> const& candidates,
                          CandsDplus const& candidatesD,
                          TracksBachPions const&)
  {
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(HfHelper::yB0(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCand<false, false, false, false>(candidate, candidatesD);
    } // candidate loop
  } // processDataDplusPi
  PROCESS_SWITCH(HfTaskB0Reduced, processDataDplusPi, "Process data without ML scores for B0 and Dplus daughter", true);

  void processDataDplusPiWithDmesMl(soa::Filtered<soa::Join<aod::HfRedCandB0, aod::HfRedB0DpMls, aod::HfSelB0ToDPi>> const& candidates,
                                    CandsDplus const& candidatesD,
                                    TracksBachPions const&)
  {
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(HfHelper::yB0(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCand<false, false, true, false>(candidate, candidatesD);
    } // candidate loop
  } // processDataDplusPiWithDmesMl
  PROCESS_SWITCH(HfTaskB0Reduced, processDataDplusPiWithDmesMl, "Process data with(out) ML scores for Dplus daughter (B0)", false);

  void processDataDplusPiWithB0Ml(soa::Filtered<soa::Join<aod::HfRedCandB0, aod::HfMlB0ToDPi, aod::HfSelB0ToDPi>> const& candidates,
                                  CandsDplus const& candidatesD,
                                  TracksBachPions const&)
  {
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(HfHelper::yB0(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCand<false, false, false, true>(candidate, candidatesD);
    } // candidate loop
  } // processDataDplusPiWithB0Ml
  PROCESS_SWITCH(HfTaskB0Reduced, processDataDplusPiWithB0Ml, "Process data with(out) ML scores for B0 (Dplus daughter)", false);

  // Process functions
  void processDataDstarPi(soa::Filtered<soa::Join<aod::HfRedCandB0DStar, aod::HfSelB0ToDPi>> const& candidates,
                          CandsDstar const& candidatesD,
                          TracksSoftPions const& softPions,
                          TracksBachPions const&)
  {
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(HfHelper::yB0(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCandDStarPi<false, false, false>(candidate, softPions, candidatesD);
    } // candidate loop
  } // processDataDstarPi
  PROCESS_SWITCH(HfTaskB0Reduced, processDataDstarPi, "Process data without ML scores for B0 and Dstar daughter", false);

  void processDataDstarPiWithDmesMl(soa::Filtered<soa::Join<aod::HfRedCandB0DStar, aod::HfRedB0DpMls, aod::HfSelB0ToDPi>> const& candidates,
                                    CandsDstar const& candidatesD,
                                    TracksSoftPions const& softPions,
                                    TracksBachPions const&)
  {
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(HfHelper::yB0(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCandDStarPi<false, false, true>(candidate, softPions, candidatesD);
    } // candidate loop
  } // processDataDstarPiWithDmesMl
  PROCESS_SWITCH(HfTaskB0Reduced, processDataDstarPiWithDmesMl, "Process data with(out) ML scores for Dstar daughter (B0)", false);

  void processMcDplusPi(soa::Filtered<soa::Join<aod::HfRedCandB0, aod::HfSelB0ToDPi, aod::HfMcRecRedB0s>> const& candidates,
                        aod::HfMcGenRedB0s const& mcParticles,
                        CandsDplus const& candidatesD,
                        TracksBachPions const&)
  {
    // MC rec
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(HfHelper::yB0(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCand<true, false, false, false>(candidate, candidatesD);
    } // rec

    // MC gen. level
    for (const auto& particle : mcParticles) {
      fillCandMcGen(particle);
    } // gen
  } // processMcDplusPi
  PROCESS_SWITCH(HfTaskB0Reduced, processMcDplusPi, "Process MC without ML scores for B0 and Dplus daughter", false);

  void processMcDplusPiWithDecayTypeCheck(soa::Filtered<soa::Join<aod::HfRedCandB0, aod::HfSelB0ToDPi, aod::HfMcRecRedB0s, aod::HfMcCheckB0s>> const& candidates,
                                          aod::HfMcGenRedB0s const& mcParticles,
                                          CandsDplus const& candidatesD,
                                          TracksBachPions const&)
  {
    // MC rec
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(HfHelper::yB0(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCand<true, true, false, false>(candidate, candidatesD);
    } // rec

    // MC gen. level
    for (const auto& particle : mcParticles) {
      fillCandMcGen(particle);
    } // gen
  } // processMcDplusPi
  PROCESS_SWITCH(HfTaskB0Reduced, processMcDplusPiWithDecayTypeCheck, "Process MC with decay type check and without ML scores for B0 and Dplus daughter", false);

  void processMcDplusPiWithDmesMl(soa::Filtered<soa::Join<aod::HfRedCandB0, aod::HfRedB0DpMls, aod::HfSelB0ToDPi, aod::HfMcRecRedB0s>> const& candidates,
                                  aod::HfMcGenRedB0s const& mcParticles,
                                  CandsDplus const& candidatesD,
                                  TracksBachPions const&)
  {
    // MC rec
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(HfHelper::yB0(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCand<true, false, true, false>(candidate, candidatesD);
    } // rec

    // MC gen. level
    for (const auto& particle : mcParticles) {
      fillCandMcGen(particle);
    } // gen
  } // processMcDplusPiWithDmesMl
  PROCESS_SWITCH(HfTaskB0Reduced, processMcDplusPiWithDmesMl, "Process MC with(out) ML scores for Dplus daughter (B0)", false);

  void processMcDplusPiWithDmesMlAndDecayTypeCheck(soa::Filtered<soa::Join<aod::HfRedCandB0, aod::HfRedB0DpMls, aod::HfSelB0ToDPi, aod::HfMcRecRedB0s, aod::HfMcCheckB0s>> const& candidates,
                                                   aod::HfMcGenRedB0s const& mcParticles,
                                                   CandsDplus const& candidatesD,
                                                   TracksBachPions const&)
  {
    // MC rec
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(HfHelper::yB0(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCand<true, true, true, false>(candidate, candidatesD);
    } // rec

    // MC gen. level
    for (const auto& particle : mcParticles) {
      fillCandMcGen(particle);
    } // gen
  } // processMcDplusPi
  PROCESS_SWITCH(HfTaskB0Reduced, processMcDplusPiWithDmesMlAndDecayTypeCheck, "Process MC with decay type check and with(out) ML scores for B0 (Dplus daughter)", false);

  void processMcDplusPiWithB0Ml(soa::Filtered<soa::Join<aod::HfRedCandB0, aod::HfMlB0ToDPi, aod::HfSelB0ToDPi, aod::HfMcRecRedB0s>> const& candidates,
                                aod::HfMcGenRedB0s const& mcParticles,
                                CandsDplus const& candidatesD,
                                TracksBachPions const&)
  {
    // MC rec
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(HfHelper::yB0(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCand<true, false, false, true>(candidate, candidatesD);
    } // rec

    // MC gen. level
    for (const auto& particle : mcParticles) {
      fillCandMcGen(particle);
    } // gen
  } // processMcDplusPiWithB0Ml
  PROCESS_SWITCH(HfTaskB0Reduced, processMcDplusPiWithB0Ml, "Process MC with(out) ML scores for B0 (Dplus daughter)", false);

  void processMcDplusPiWithB0MlAndDecayTypeCheck(soa::Filtered<soa::Join<aod::HfRedCandB0, aod::HfMlB0ToDPi, aod::HfSelB0ToDPi, aod::HfMcRecRedB0s, aod::HfMcCheckB0s>> const& candidates,
                                                 aod::HfMcGenRedB0s const& mcParticles,
                                                 CandsDplus const& candidatesD,
                                                 TracksBachPions const&)
  {
    // MC rec
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(HfHelper::yB0(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCand<true, true, false, true>(candidate, candidatesD);
    } // rec

    // MC gen. level
    for (const auto& particle : mcParticles) {
      fillCandMcGen(particle);
    } // gen
  } // processMcDplusPi
  PROCESS_SWITCH(HfTaskB0Reduced, processMcDplusPiWithB0MlAndDecayTypeCheck, "Process MC with decay type check and with(out) ML scores for B0 (Dplus daughter)", false);

  void processMcDstarPi(soa::Filtered<soa::Join<aod::HfRedCandB0DStar, aod::HfSelB0ToDPi, aod::HfMcRecRedB0s>> const& candidates,
                        aod::HfMcGenRedB0s const& mcParticles,
                        CandsDstar const& candidatesD,
                        TracksSoftPions const& softPions,
                        TracksBachPions const&)
  {
    // MC rec
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(HfHelper::yB0(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCandDStarPi<true, false, false>(candidate, softPions, candidatesD);
    } // rec

    // MC gen. level
    for (const auto& particle : mcParticles) {
      fillCandMcGen(particle);
    } // gen
  } // processMcDstarPi
  PROCESS_SWITCH(HfTaskB0Reduced, processMcDstarPi, "Process MC without ML scores for B0 and Dstar daughter", false);

  void processMcDstarPiWithDmesMl(soa::Filtered<soa::Join<aod::HfRedCandB0DStar, aod::HfRedB0DpMls, aod::HfSelB0ToDPi, aod::HfMcRecRedB0s>> const& candidates,
                                  aod::HfMcGenRedB0s const& mcParticles,
                                  CandsDstar const& candidatesD,
                                  TracksSoftPions const& softPions,
                                  TracksBachPions const&)
  {
    // MC rec
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(HfHelper::yB0(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCandDStarPi<true, false, true>(candidate, softPions, candidatesD);
    } // rec

    // MC gen. level
    for (const auto& particle : mcParticles) {
      fillCandMcGen(particle);
    } // gen
  } // processMcDstarPiWithDmesMl
  PROCESS_SWITCH(HfTaskB0Reduced, processMcDstarPiWithDmesMl, "Process MC with(out) ML scores for Dstar daughter (B0)", false);

}; // struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskB0Reduced>(cfgc)};
}
