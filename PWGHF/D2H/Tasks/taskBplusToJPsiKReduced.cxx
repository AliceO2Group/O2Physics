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

/// \file taskBplusToJPsiKReduced.cxx
/// \brief B+ → JPsi K+ → (µ+ µ-) K+ analysis task
///
/// \author Fabrizio Chinu <fabrizio.chinu@cern.ch>, Università degli Studi and INFN Torino
/// \author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN

#include <vector>
#include <string>

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "Common/Core/RecoDecay.h"

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/HfMlResponseBplusToJPsiKReduced.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/D2H/DataModel/ReducedDataModel.h"
#include "PWGHF/Utils/utilsPid.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::pid_tpc_tof_utils;

namespace o2::aod
{
namespace hf_cand_bplustojpsik_lite
{
DECLARE_SOA_COLUMN(PtJPsi, ptJPsi, float); //! Transverse momentum of JPsi daughter candidate (GeV/c)
DECLARE_SOA_COLUMN(PtBach, ptBach, float); //! Transverse momentum of bachelor kaon (GeV/c)
// DECLARE_SOA_COLUMN(AbsEtaBach, absEtaBach, float);                                       //! Absolute pseudorapidity of bachelor kaon
// DECLARE_SOA_COLUMN(ItsNClsBach, itsNClsBach, int);                                       //! Number of ITS clusters of bachelor kaon
// DECLARE_SOA_COLUMN(TpcNClsCrossedRowsBach, tpcNClsCrossedRowsBach, int);                 //! Number of TPC crossed rows of prongs of bachelor kaon
// DECLARE_SOA_COLUMN(TpcChi2NClBach, tpcChi2NClBach, float);                               //! Maximum TPC chi2 of prongs of JPsi-meson daughter candidate
// DECLARE_SOA_COLUMN(PtJPsiProngMin, ptJPsiProngMin, float);                               //! Minimum pT of prongs of JPsi daughter candidate (GeV/c)
// DECLARE_SOA_COLUMN(AbsEtaJPsiProngMin, absEtaJPsiProngMin, float);                       //! Minimum absolute pseudorapidity of prongs of JPsi daughter candidate
// DECLARE_SOA_COLUMN(ItsNClsJPsiProngMin, itsNClsJPsiProngMin, int);                       //! Minimum number of ITS clusters of prongs of JPsi daughter candidate
// DECLARE_SOA_COLUMN(TpcNClsCrossedRowsJPsiProngMin, tpcNClsCrossedRowsJPsiProngMin, int); //! Minimum number of TPC crossed rows of prongs of JPsi daughter candidate
// DECLARE_SOA_COLUMN(TpcChi2NClJPsiProngMax, tpcChi2NClJPsiProngMax, float);               //! Maximum TPC chi2 of prongs of JPsi daughter candidate
DECLARE_SOA_COLUMN(MJPsi, mJPsi, float);                                           //! Invariant mass of JPsi daughter candidates (GeV/c)
DECLARE_SOA_COLUMN(M, m, float);                                                   //! Invariant mass of candidate (GeV/c2)
DECLARE_SOA_COLUMN(Pt, pt, float);                                                 //! Transverse momentum of candidate (GeV/c)
DECLARE_SOA_COLUMN(PtGen, ptGen, float);                                           //! Transverse momentum of candidate (GeV/c)
DECLARE_SOA_COLUMN(P, p, float);                                                   //! Momentum of candidate (GeV/c)
DECLARE_SOA_COLUMN(Y, y, float);                                                   //! Rapidity of candidate
DECLARE_SOA_COLUMN(Eta, eta, float);                                               //! Pseudorapidity of candidate
DECLARE_SOA_COLUMN(Phi, phi, float);                                               //! Azimuth angle of candidate
DECLARE_SOA_COLUMN(E, e, float);                                                   //! Energy of candidate (GeV)
DECLARE_SOA_COLUMN(NSigTpcKaBachelor, nSigTpcKaBachelor, float);                   //! TPC Nsigma separation for bachelor with kaon mass hypothesis
DECLARE_SOA_COLUMN(NSigTofKaBachelor, nSigTofKaBachelor, float);                   //! TOF Nsigma separation for bachelor with kaon mass hypothesis
DECLARE_SOA_COLUMN(NSigTpcTofKaBachelor, nSigTpcTofKaBachelor, float);             //! Combined TPC and TOF Nsigma separation for bachelor with kaon mass hypothesis
DECLARE_SOA_COLUMN(NSigTpcMuJPsiDauPos, nSigTpcMuJPsiDauPos, float);               //! TPC Nsigma separation for JPsi DauPos with muon mass hypothesis
DECLARE_SOA_COLUMN(NSigTofMuJPsiDauPos, nSigTofMuJPsiDauPos, float);               //! TOF Nsigma separation for JPsi DauPos with muon mass hypothesis
DECLARE_SOA_COLUMN(NSigTpcTofMuJPsiDauPos, nSigTpcTofMuJPsiDauPos, float);         //! Combined TPC and TOF Nsigma separation for JPsi prong0 with muon mass hypothesis
DECLARE_SOA_COLUMN(NSigTpcMuJPsiDauNeg, nSigTpcMuJPsiDauNeg, float);               //! TPC Nsigma separation for JPsi DauNeg with muon mass hypothesis
DECLARE_SOA_COLUMN(NSigTofMuJPsiDauNeg, nSigTofMuJPsiDauNeg, float);               //! TOF Nsigma separation for JPsi DauNeg with muon mass hypothesis
DECLARE_SOA_COLUMN(NSigTpcTofMuJPsiDauNeg, nSigTpcTofMuJPsiDauNeg, float);         //! Combined TPC and TOF Nsigma separation for JPsi prong1 with muon mass hypothesis
DECLARE_SOA_COLUMN(DecayLength, decayLength, float);                               //! Decay length of candidate (cm)
DECLARE_SOA_COLUMN(DecayLengthXY, decayLengthXY, float);                           //! Transverse decay length of candidate (cm)
DECLARE_SOA_COLUMN(DecayLengthNormalised, decayLengthNormalised, float);           //! Normalised decay length of candidate
DECLARE_SOA_COLUMN(DecayLengthXYNormalised, decayLengthXYNormalised, float);       //! Normalised transverse decay length of candidate
DECLARE_SOA_COLUMN(ImpactParameterProduct, impactParameterProduct, float);         //! Impact parameter product of B daughters
DECLARE_SOA_COLUMN(ImpactParameterProductJPsi, impactParameterProductJPsi, float); //! Impact parameter product of JPsi daughters
DECLARE_SOA_COLUMN(ImpactParameterJPsiDauPos, impactParameterJPsiDauPos, float);   //! Impact parameter of JPsi daughter candidate
DECLARE_SOA_COLUMN(ImpactParameterJPsiDauNeg, impactParameterJPsiDauNeg, float);   //! Impact parameter of JPsi daughter candidate
DECLARE_SOA_COLUMN(ImpactParameterLfTrack0, impactParameterLfTrack0, float);       //! Impact parameter of Phi daughter candidate
DECLARE_SOA_COLUMN(Cpa, cpa, float);                                               //! Cosine pointing angle of candidate
DECLARE_SOA_COLUMN(CpaXY, cpaXY, float);                                           //! Cosine pointing angle of candidate in transverse plane
DECLARE_SOA_COLUMN(CpaJPsi, cpaJPsi, float);                                       //! Cosine pointing angle of JPsi daughter candidate
DECLARE_SOA_COLUMN(CpaXYJPsi, cpaXYJPsi, float);                                   //! Cosine pointing angle in transverse plane of JPsi daughter candidate
DECLARE_SOA_COLUMN(MaxNormalisedDeltaIP, maxNormalisedDeltaIP, float);             //! Maximum normalized difference between measured and expected impact parameter of candidate prongs
DECLARE_SOA_COLUMN(MlScoreSig, mlScoreSig, float);                                 //! ML score for signal class
DECLARE_SOA_COLUMN(FlagWrongCollision, flagWrongCollision, int8_t);                //! Flag for association with wrong collision
} // namespace hf_cand_bplustojpsik_lite

DECLARE_SOA_TABLE(HfRedCandBpLites, "AOD", "HFREDCANDBPLITE", //! Table with some B+ properties
                  hf_cand_bplustojpsik_lite::M,
                  hf_cand_bplustojpsik_lite::Pt,
                  hf_cand_bplustojpsik_lite::Eta,
                  hf_cand_bplustojpsik_lite::Phi,
                  hf_cand_bplustojpsik_lite::Y,
                  hf_cand_bplustojpsik_lite::Cpa,
                  hf_cand_bplustojpsik_lite::CpaXY,
                  hf_cand::Chi2PCA,
                  hf_cand_bplustojpsik_lite::DecayLength,
                  hf_cand_bplustojpsik_lite::DecayLengthXY,
                  hf_cand_bplustojpsik_lite::DecayLengthNormalised,
                  hf_cand_bplustojpsik_lite::DecayLengthXYNormalised,
                  hf_cand_bplustojpsik_lite::ImpactParameterProduct,
                  hf_cand_bplustojpsik_lite::ImpactParameterProductJPsi,
                  hf_cand_bplustojpsik_lite::MaxNormalisedDeltaIP,
                  hf_cand_bplustojpsik_lite::MlScoreSig,
                  // hf_sel_candidate_bplus::IsSelBplusToJPsiPi,
                  //  JPsi meson features
                  hf_cand_bplustojpsik_lite::MJPsi,
                  hf_cand_bplustojpsik_lite::PtJPsi,
                  hf_cand_bplustojpsik_lite::ImpactParameterJPsiDauPos,
                  hf_cand_bplustojpsik_lite::ImpactParameterJPsiDauNeg,
                  hf_cand_bplustojpsik_lite::ImpactParameterLfTrack0,
                  // hf_cand_bplustojpsik_lite::PtJPsiProngMin,
                  // hf_cand_bplustojpsik_lite::AbsEtaJPsiProngMin,
                  // hf_cand_bplustojpsik_lite::ItsNClsJPsiProngMin,
                  // hf_cand_bplustojpsik_lite::TpcNClsCrossedRowsJPsiProngMin,
                  // hf_cand_bplustojpsik_lite::TpcChi2NClJPsiProngMax,
                  // kaon features
                  hf_cand_bplustojpsik_lite::PtBach,
                  // hf_cand_bplustojpsik_lite::AbsEtaBach,
                  // hf_cand_bplustojpsik_lite::ItsNClsBach,
                  // hf_cand_bplustojpsik_lite::TpcNClsCrossedRowsBach,
                  // hf_cand_bplustojpsik_lite::TpcChi2NClBach,
                  hf_cand_bplustojpsik_lite::NSigTpcKaBachelor,
                  hf_cand_bplustojpsik_lite::NSigTofKaBachelor,
                  hf_cand_bplustojpsik_lite::NSigTpcTofKaBachelor,
                  // MC truth
                  hf_cand_2prong::FlagMcMatchRec,
                  hf_cand_2prong::OriginMcRec,
                  hf_cand_bplustojpsik_lite::FlagWrongCollision,
                  hf_cand_bplustojpsik_lite::PtGen);

// DECLARE_SOA_TABLE(HfRedBpMcCheck, "AOD", "HFREDBPMCCHECK", //! Table with MC decay type check
//                   hf_cand_2prong::FlagMcMatchRec,
//                   hf_cand_bplustojpsik_lite::FlagWrongCollision,
//                   hf_cand_bplustojpsik_lite::MJPsi,
//                   hf_cand_bplustojpsik_lite::PtJPsi,
//                   hf_cand_bplustojpsik_lite::M,
//                   hf_cand_bplustojpsik_lite::Pt,
//                   // hf_cand_bplustojpsik_lite::MlScoreSig,
//                   hf_bplus_mc::PdgCodeBeautyMother,
//                   hf_bplus_mc::PdgCodeCharmMother,
//                   hf_bplus_mc::PdgCodeDauPos,
//                   hf_bplus_mc::PdgCodeDauNeg,
//                   hf_bplus_mc::PdgCodeProng2);
} // namespace o2::aod

// string definitions, used for histogram axis labels
const TString stringPt = "#it{p}_{T} (GeV/#it{c})";
const TString stringPtJPsi = "#it{p}_{T}(JPsi) (GeV/#it{c});";
const TString bPlusCandTitle = "B+ candidates;";
const TString entries = "entries";
const TString bPlusCandMatch = "B+ candidates (matched);";
const TString bPlusCandUnmatch = "B+ candidates (unmatched);";
const TString mcParticleMatched = "MC particles (matched);";

/// B+ analysis task
struct HfTaskBplusToJPsiKReduced {
  Produces<aod::HfRedCandBpLites> hfRedCandBpLite;
  // Produces<aod::HfRedBpMcCheck> hfRedBpMcCheck;

  Configurable<uint8_t> selectionFlagBplus{"selectionFlagBplus", 1, "Selection Flag for Bplus"};
  Configurable<double> yCandGenMax{"yCandGenMax", 0.5, "max. gen particle rapidity"};
  Configurable<double> yCandRecoMax{"yCandRecoMax", 0.8, "max. cand. rapidity"};
  Configurable<float> etaTrackMax{"etaTrackMax", 0.8, "max. track pseudo-rapidity"};
  Configurable<float> ptTrackMin{"ptTrackMin", 0.1, "min. track transverse momentum"};
  Configurable<bool> fillBackground{"fillBackground", false, "Flag to enable filling of background histograms/sparses/tree (only MC)"};
  Configurable<float> downSampleBkgFactor{"downSampleBkgFactor", 1., "Fraction of background candidates to keep for ML trainings"};
  Configurable<float> ptMaxForDownSample{"ptMaxForDownSample", 10., "Maximum pt for the application of the downsampling factor"};
  // topological cuts
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_bplus_to_jpsi_k::vecBinsPt}, "pT bin limits"};
  Configurable<LabeledArray<double>> cuts{"cuts", {hf_cuts_bplus_to_jpsi_k::Cuts[0], hf_cuts_bplus_to_jpsi_k::NBinsPt, hf_cuts_bplus_to_jpsi_k::NCutVars, hf_cuts_bplus_to_jpsi_k::labelsPt, hf_cuts_bplus_to_jpsi_k::labelsCutVar}, "B+ candidate selection per pT bin"};
  // Enable PID
  Configurable<int> kaonPidMethod{"kaonPidMethod", PidMethod::TpcOrTof, "PID selection method for the bachelor kaon (PidMethod::NoPid: none, PidMethod::TpcOrTof: TPC or TOF, PidMethod::TpcAndTof: TPC and TOF)"};
  Configurable<bool> acceptPIDNotApplicable{"acceptPIDNotApplicable", true, "Switch to accept Status::NotApplicable [(NotApplicable for one detector) and (NotApplicable or Conditional for the other)] in PID selection"};
  // TPC PID
  Configurable<double> ptPidTpcMin{"ptPidTpcMin", 0.15, "Lower bound of track pT for TPC PID"};
  Configurable<double> ptPidTpcMax{"ptPidTpcMax", 20., "Upper bound of track pT for TPC PID"};
  Configurable<double> nSigmaTpcMax{"nSigmaTpcMax", 5., "Nsigma cut on TPC only"};
  Configurable<double> nSigmaTpcCombinedMax{"nSigmaTpcCombinedMax", 5., "Nsigma cut on TPC combined with TOF"};
  // TOF PID
  Configurable<double> ptPidTofMin{"ptPidTofMin", 0.15, "Lower bound of track pT for TOF PID"};
  Configurable<double> ptPidTofMax{"ptPidTofMax", 20., "Upper bound of track pT for TOF PID"};
  Configurable<double> nSigmaTofMax{"nSigmaTofMax", 5., "Nsigma cut on TOF only"};
  Configurable<double> nSigmaTofCombinedMax{"nSigmaTofCombinedMax", 5., "Nsigma cut on TOF combined with TPC"};
  // B+ ML inference
  Configurable<std::vector<double>> binsPtBpMl{"binsPtBpMl", std::vector<double>{hf_cuts_ml::vecBinsPt}, "pT bin limits for ML application"};
  Configurable<std::vector<int>> cutDirBpMl{"cutDirBpMl", std::vector<int>{hf_cuts_ml::vecCutDir}, "Whether to reject score values greater or smaller than the threshold"};
  Configurable<LabeledArray<double>> cutsBpMl{"cutsBpMl", {hf_cuts_ml::Cuts[0], hf_cuts_ml::NBinsPt, hf_cuts_ml::NCutScores, hf_cuts_ml::labelsPt, hf_cuts_ml::labelsCutScore}, "ML selections per pT bin"};
  Configurable<int> nClassesBpMl{"nClassesBpMl", static_cast<int>(hf_cuts_ml::NCutScores), "Number of classes in ML model"};
  Configurable<std::vector<std::string>> namesInputFeatures{"namesInputFeatures", std::vector<std::string>{"feature1", "feature2"}, "Names of ML model input features"};
  // CCDB configuration
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::vector<std::string>> modelPathsCCDB{"modelPathsCCDB", std::vector<std::string>{"path_ccdb/BDT_BPLUS/"}, "Paths of models on CCDB"};
  Configurable<std::vector<std::string>> onnxFileNames{"onnxFileNames", std::vector<std::string>{"ModelHandler_onnx_BPLUSToJPSIK.onnx"}, "ONNX file names for each pT bin (if not from CCDB full path)"};
  Configurable<int64_t> timestampCCDB{"timestampCCDB", -1, "timestamp of the ONNX file for ML model used to query in CCDB"};
  Configurable<bool> loadModelsFromCCDB{"loadModelsFromCCDB", false, "Flag to enable or disable the loading of models from CCDB"};

  HfHelper hfHelper;
  TrackSelectorKa selectorKaon;
  o2::analysis::HfMlResponseBplusToJPsiKReduced<float> hfMlResponse;
  o2::ccdb::CcdbApi ccdbApi;

  using TracksKaon = soa::Join<HfRedTracks, HfRedTracksPid>;
  std::vector<float> outputMl = {};

  // Filter filterSelectCandidates = (aod::hf_sel_candidate_bplus::isSelBplusToJPsiPi >= selectionFlagBplus);

  HistogramRegistry registry{"registry"};

  void init(InitContext&)
  {
    std::array<bool, 2> processFuncData{doprocessData, doprocessDataWithBplusMl};
    if ((std::accumulate(processFuncData.begin(), processFuncData.end(), 0)) > 1) {
      LOGP(fatal, "Only one process function for data can be enabled at a time.");
    }
    std::array<bool, 2> processFuncMc{doprocessMc, doprocessMcWithBplusMl};
    if ((std::accumulate(processFuncMc.begin(), processFuncMc.end(), 0)) > 1) {
      LOGP(fatal, "Only one process function for MC can be enabled at a time.");
    }

    if (kaonPidMethod < 0 || kaonPidMethod >= PidMethod::NPidMethods) {
      LOGP(fatal, "Invalid PID option in configurable, please set 0 (no PID), 1 (TPC or TOF), or 2 (TPC and TOF)");
    }

    if (kaonPidMethod != PidMethod::NoPid) {
      selectorKaon.setRangePtTpc(ptPidTpcMin, ptPidTpcMax);
      selectorKaon.setRangeNSigmaTpc(-nSigmaTpcMax, nSigmaTpcMax);
      selectorKaon.setRangeNSigmaTpcCondTof(-nSigmaTpcCombinedMax, nSigmaTpcCombinedMax);
      selectorKaon.setRangePtTof(ptPidTofMin, ptPidTofMax);
      selectorKaon.setRangeNSigmaTof(-nSigmaTofMax, nSigmaTofMax);
      selectorKaon.setRangeNSigmaTofCondTpc(-nSigmaTofCombinedMax, nSigmaTofCombinedMax);
    }

    const AxisSpec axisMassBplus{150, 4.5, 6.0};
    const AxisSpec axisMassJPsi{600, 2.8f, 3.4f};
    const AxisSpec axisPtProng{100, 0., 10.};
    const AxisSpec axisImpactPar{200, -0.05, 0.05};
    const AxisSpec axisPtJPsi{100, 0., 50.};
    const AxisSpec axisRapidity{100, -2., 2.};
    const AxisSpec axisPtB{(std::vector<double>)binsPt, "#it{p}_{T}^{B^{+}} (GeV/#it{c})"};
    const AxisSpec axisPtKa{100, 0.f, 10.f};

    registry.add("hMass", bPlusCandTitle + "inv. mass J/#Psi K^{+} (GeV/#it{c}^{2});" + stringPt, {HistType::kTH2F, {axisMassBplus, axisPtB}});
    registry.add("hMassJPsi", bPlusCandTitle + "inv. mass #mu^{+}#mu^{#minus} (GeV/#it{c}^{2});" + stringPt, {HistType::kTH2F, {axisMassJPsi, axisPtJPsi}});
    registry.add("hd0K", bPlusCandTitle + "Kaon DCAxy to prim. vertex (cm);" + stringPt, {HistType::kTH2F, {axisImpactPar, axisPtKa}});

    // histograms processMC
    if (doprocessMc || doprocessMcWithBplusMl) {
      registry.add("hPtJPsiGen", mcParticleMatched + "J/#Psi #it{p}_{T}^{gen} (GeV/#it{c}); B^{+} " + stringPt, {HistType::kTH2F, {axisPtProng, axisPtB}});
      registry.add("hPtKGen", mcParticleMatched + "Kaon #it{p}_{T}^{gen} (GeV/#it{c}); B^{+} " + stringPt, {HistType::kTH2F, {axisPtProng, axisPtB}});
      registry.add("hYGenWithProngsInAcceptance", mcParticleMatched + "Kaon #it{p}_{T}^{gen} (GeV/#it{c}); B^{+} " + stringPt, {HistType::kTH2F, {axisPtProng, axisRapidity}});
      registry.add("hMassRecSig", bPlusCandMatch + "inv. mass J/#Psi K^{+} (GeV/#it{c}^{2}); B^{+} " + stringPt, {HistType::kTH2F, {axisMassBplus, axisPtB}});
      registry.add("hMassJPsiRecSig", bPlusCandMatch + "inv. mass #mu^{+}#mu^{#minus} (GeV/#it{c}^{2}); J/#Psi " + stringPt, {HistType::kTH2F, {axisMassJPsi, axisPtJPsi}});
      registry.add("hd0KRecSig", bPlusCandMatch + "Kaon DCAxy to prim. vertex (cm); K^{+} " + stringPt, {HistType::kTH2F, {axisImpactPar, axisPtKa}});
      registry.add("hMassRecBg", bPlusCandUnmatch + "inv. mass J/#Psi K^{+} (GeV/#it{c}^{2}); B^{+} " + stringPt, {HistType::kTH2F, {axisMassBplus, axisPtB}});
      registry.add("hMassJPsiRecBg", bPlusCandUnmatch + "inv. mass #mu^{+}#mu^{#minus} (GeV/#it{c}^{2}); J/#Psi " + stringPt, {HistType::kTH2F, {axisMassJPsi, axisPtJPsi}});
      registry.add("hd0KRecBg", bPlusCandMatch + "Kaon DCAxy to prim. vertex (cm); K^{+} " + stringPt, {HistType::kTH2F, {axisImpactPar, axisPtKa}});
    }

    if (doprocessDataWithBplusMl || doprocessMcWithBplusMl) {
      hfMlResponse.configure(binsPtBpMl, cutsBpMl, cutDirBpMl, nClassesBpMl);
      if (loadModelsFromCCDB) {
        ccdbApi.init(ccdbUrl);
        hfMlResponse.setModelPathsCCDB(onnxFileNames, ccdbApi, modelPathsCCDB, timestampCCDB);
      } else {
        hfMlResponse.setModelPathsLocal(onnxFileNames);
      }
      hfMlResponse.cacheInputFeaturesIndices(namesInputFeatures);
      hfMlResponse.init();
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
  /// \param withBplusMl is the flag to enable the filling with ML scores for the B+ candidate
  /// \param candidate is the B+ candidate
  /// \param candidatesJPsi is the table with JPsi candidates
  template <bool doMc, bool withBplusMl, typename Cand>
  void fillCand(Cand const& candidate,
                aod::HfRedJPsis const& /*candidatesJPsi*/,
                aod::HfRedBach0Tracks const&)
  {
    auto ptCandBplus = candidate.pt();
    auto invMassBplus = hfHelper.invMassBplusToJPsiK(candidate);
    auto candJPsi = candidate.template jPsi_as<aod::HfRedJPsis>();
    auto candKa = candidate.template bachKa_as<aod::HfRedBach0Tracks>();
    auto ptJPsi = candidate.ptProng0();
    auto invMassJPsi = candJPsi.m();
    uint8_t statusBplus = 0;

    int8_t flagMcMatchRec = 0;
    int8_t flagWrongCollision = 0;
    bool isSignal = false;
    if constexpr (doMc) {
      flagMcMatchRec = candidate.flagMcMatchRec();
      flagWrongCollision = candidate.flagWrongCollision();
      isSignal = TESTBIT(std::abs(flagMcMatchRec), static_cast<int>(hf_cand_bplus::DecayTypeBToJPsiMc::BplusToJPsiKToMuMuK));
    }

    SETBIT(statusBplus, SelectionStep::RecoSkims);
    if (hfHelper.selectionBplusToJPsiKTopol(candidate, cuts, binsPt)) {
      SETBIT(statusBplus, SelectionStep::RecoTopol);
    } else if (selectionFlagBplus >= BIT(SelectionStep::RecoTopol) * 2 - 1) {
      return;
    }
    // track-level PID selection
    // auto trackKa = candidate.template prong1_as<TracksPion>();
    if (kaonPidMethod == PidMethod::TpcOrTof || kaonPidMethod == PidMethod::TpcAndTof) {
      int pidTrackKa{TrackSelectorPID::Status::NotApplicable};
      if (kaonPidMethod == PidMethod::TpcOrTof) {
        pidTrackKa = selectorKaon.statusTpcOrTof(candKa);
      } else if (kaonPidMethod == PidMethod::TpcAndTof) {
        pidTrackKa = selectorKaon.statusTpcAndTof(candKa);
      }
      if (hfHelper.selectionBplusToJPsiKPid(pidTrackKa, acceptPIDNotApplicable.value)) {
        // LOGF(info, "B+ candidate selection failed at PID selection");
        SETBIT(statusBplus, SelectionStep::RecoPID);
      } else if (selectionFlagBplus >= BIT(SelectionStep::RecoPID) * 2 - 1) {
        return;
      }
    }

    float candidateMlScoreSig = -1;
    if constexpr (withBplusMl) {
      // B+ ML selections
      std::vector<float> inputFeatures = hfMlResponse.getInputFeatures(candidate, candKa);
      if (hfMlResponse.isSelectedMl(inputFeatures, ptCandBplus, outputMl)) {
        SETBIT(statusBplus, SelectionStep::RecoMl);
      } else if (selectionFlagBplus >= BIT(SelectionStep::RecoMl) * 2 - 1) {
        return;
      }
      candidateMlScoreSig = outputMl[1];
    }

    registry.fill(HIST("hMass"), invMassBplus, ptCandBplus);
    registry.fill(HIST("hMassJPsi"), invMassJPsi, candidate.ptProng0());
    registry.fill(HIST("hd0K"), candidate.impactParameter1(), candidate.ptProng1());
    if constexpr (doMc) {
      if (isSignal) {
        registry.fill(HIST("hMassRecSig"), invMassBplus, ptCandBplus);
        registry.fill(HIST("hMassJPsiRecSig"), invMassJPsi, candidate.ptProng0());
        registry.fill(HIST("hd0KRecSig"), candidate.impactParameter1(), candidate.ptProng1());
      } else if (fillBackground) {
        registry.fill(HIST("hMassRecBg"), invMassBplus, ptCandBplus);
        registry.fill(HIST("hMassJPsiRecBg"), invMassJPsi, candidate.ptProng0());
        registry.fill(HIST("hd0KRecBg"), candidate.impactParameter1(), candidate.ptProng1());
      }
    }

    float pseudoRndm = ptJPsi * 1000. - static_cast<int64_t>(ptJPsi * 1000);
    if (ptCandBplus >= ptMaxForDownSample || pseudoRndm < downSampleBkgFactor) {
      float ptMother = -1.;
      if constexpr (doMc) {
        ptMother = candidate.ptMother();
      }

      hfRedCandBpLite(
        // B+ - meson features
        invMassBplus,
        ptCandBplus,
        candidate.eta(),
        candidate.phi(),
        hfHelper.yBplus(candidate),
        candidate.cpa(),
        candidate.cpaXY(),
        candidate.chi2PCA(),
        candidate.decayLength(),
        candidate.decayLengthXY(),
        candidate.decayLengthNormalised(),
        candidate.decayLengthXYNormalised(),
        candidate.impactParameterProduct(),
        candidate.impactParameterProductJPsi(),
        candidate.maxNormalisedDeltaIP(),
        candidateMlScoreSig,
        // J/Psi features
        invMassJPsi,
        ptJPsi,
        candidate.impactParameter0(),
        candidate.impactParameter1(),
        candidate.impactParameter2(),
        // candJPsi.ptProngMin(),
        // candJPsi.absEtaProngMin(),
        // candJPsi.itsNClsProngMin(),
        // candJPsi.tpcNClsCrossedRowsProngMin(),
        // candJPsi.tpcChi2NClProngMax(),
        // kaon features
        candidate.ptProng1(),
        // std::abs(RecoDecay::eta(candKa.pVector())),
        // candKa.itsNCls(),
        // candKa.tpcNClsCrossedRows(),
        // candKa.tpcChi2NCl(),
        candKa.tpcNSigmaKa(),
        candKa.tofNSigmaKa(),
        candKa.tpcTofNSigmaKa(),
        // MC truth
        flagMcMatchRec,
        isSignal,
        flagWrongCollision,
        ptMother);
    }
  }

  /// Fill particle histograms (gen MC truth)
  void fillCandMcGen(aod::HfMcGenRedBps::iterator const& particle)
  {
    auto ptParticle = particle.ptTrack();
    auto yParticle = particle.yTrack();
    if (yCandGenMax >= 0. && std::abs(yParticle) > yCandGenMax) {
      return;
    }
    std::array<float, 2> ptProngs = {particle.ptProng0(), particle.ptProng1()};
    std::array<float, 2> etaProngs = {particle.etaProng0(), particle.etaProng1()};
    bool prongsInAcc = isProngInAcceptance(etaProngs[0], ptProngs[0]) && isProngInAcceptance(etaProngs[1], ptProngs[1]);

    registry.fill(HIST("hPtJPsiGen"), ptProngs[0], ptParticle);
    registry.fill(HIST("hPtKGen"), ptProngs[1], ptParticle);

    // generated B+ with daughters in geometrical acceptance
    if (prongsInAcc) {
      registry.fill(HIST("hYGenWithProngsInAcceptance"), ptParticle, yParticle);
    }
  }

  // Process functions
  void processData(aod::HfRedCandBplusToJPsiK const& candidates,
                   aod::HfRedJPsis const& candidatesJPsi,
                   aod::HfRedBach0Tracks const& kaonTracks)
  {
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(hfHelper.yBplus(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCand<false, false>(candidate, candidatesJPsi, kaonTracks);
    } // candidate loop
  } // processData
  PROCESS_SWITCH(HfTaskBplusToJPsiKReduced, processData, "Process data without ML for B+", true);

  void processDataWithBplusMl(aod::HfRedCandBplusToJPsiK const& candidates,
                              aod::HfRedJPsis const& candidatesJPsi,
                              aod::HfRedBach0Tracks const& kaonTracks)
  {
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(hfHelper.yBplus(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCand<false, true>(candidate, candidatesJPsi, kaonTracks);
    } // candidate loop
  } // processDataWithBplusMl
  PROCESS_SWITCH(HfTaskBplusToJPsiKReduced, processDataWithBplusMl, "Process data with ML for B+", false);

  void processMc(soa::Join<aod::HfRedCandBplusToJPsiK, aod::HfMcRecRedBps> const& candidates,
                 aod::HfMcGenRedBps const& mcParticles,
                 aod::HfRedJPsis const& candidatesJPsi,
                 aod::HfRedBach0Tracks const& kaonTracks)
  {
    // MC rec
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(hfHelper.yBplus(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCand<true, false>(candidate, candidatesJPsi, kaonTracks);
    } // rec

    // MC gen. level
    for (const auto& particle : mcParticles) {
      fillCandMcGen(particle);
    } // gen
  } // processMc
  PROCESS_SWITCH(HfTaskBplusToJPsiKReduced, processMc, "Process MC without ML for B+", false);

  void processMcWithBplusMl(soa::Join<aod::HfRedCandBplusToJPsiK, aod::HfMcRecRedBps> const& candidates,
                            aod::HfMcGenRedBps const& mcParticles,
                            aod::HfRedJPsis const& candidatesJPsi,
                            aod::HfRedBach0Tracks const& kaonTracks)
  {
    // MC rec
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(hfHelper.yBplus(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCand<true, true>(candidate, candidatesJPsi, kaonTracks);
    } // rec

    // MC gen. level
    for (const auto& particle : mcParticles) {
      fillCandMcGen(particle);
    } // gen
  } // processMcWithBplusMl
  PROCESS_SWITCH(HfTaskBplusToJPsiKReduced, processMcWithBplusMl, "Process MC with ML for B+", false);

}; // struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskBplusToJPsiKReduced>(cfgc)};
}
