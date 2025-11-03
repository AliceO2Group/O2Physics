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

/// \file taskBsToJpsiPhiReduced.cxx
/// \brief Bs → Jpsi phi → (µ+ µ-) (K+K-) analysis task
///
/// \author Fabrizio Chinu <fabrizio.chinu@cern.ch>, Università degli Studi and INFN Torino
/// \author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN

#include "PWGHF/Core/DecayChannels.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/HfMlResponseBsToJpsiPhiReduced.h"
#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/D2H/DataModel/ReducedDataModel.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/Utils/utilsPid.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelectorPID.h"

#include <CCDB/CcdbApi.h>
#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Array2D.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/WorkflowSpec.h>
#include <Framework/runDataProcessing.h>

#include <TString.h>

#include <Rtypes.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::aod;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::pid_tpc_tof_utils;

namespace o2::aod
{
namespace hf_cand_bstojpsiphi_lite
{
DECLARE_SOA_COLUMN(PtJpsi, ptJpsi, float);                                           //! Transverse momentum of Jpsi daughter candidate (GeV/c)
DECLARE_SOA_COLUMN(PtBach0, ptBach0, float);                                         //! Transverse momentum of bachelor kaon(<- phi) (GeV/c)
DECLARE_SOA_COLUMN(PtBach1, ptBach1, float);                                         //! Transverse momentum of bachelor kaon(<- phi) (GeV/c)
DECLARE_SOA_COLUMN(ItsNClsJpsiDauPos, itsNClsJpsiDauPos, int);                       //! Number of clusters in ITS
DECLARE_SOA_COLUMN(TpcNClsCrossedRowsJpsiDauPos, tpcNClsCrossedRowsJpsiDauPos, int); //! Number of TPC crossed rows
DECLARE_SOA_COLUMN(ItsChi2NClJpsiDauPos, itsChi2NClJpsiDauPos, float);               //! ITS chi2 / Number of clusters
DECLARE_SOA_COLUMN(TpcChi2NClJpsiDauPos, tpcChi2NClJpsiDauPos, float);               //! TPC chi2 / Number of clusters
DECLARE_SOA_COLUMN(AbsEtaJpsiDauPos, absEtaJpsiDauPos, float);                       //! |eta|
DECLARE_SOA_COLUMN(ItsNClsJpsiDauNeg, itsNClsJpsiDauNeg, int);                       //! Number of clusters in ITS
DECLARE_SOA_COLUMN(TpcNClsCrossedRowsJpsiDauNeg, tpcNClsCrossedRowsJpsiDauNeg, int); //! Number of TPC crossed rows
DECLARE_SOA_COLUMN(ItsChi2NClJpsiDauNeg, itsChi2NClJpsiDauNeg, float);               //! ITS chi2 / Number of clusters
DECLARE_SOA_COLUMN(TpcChi2NClJpsiDauNeg, tpcChi2NClJpsiDauNeg, float);               //! TPC chi2 / Number of clusters
DECLARE_SOA_COLUMN(AbsEtaJpsiDauNeg, absEtaJpsiDauNeg, float);                       //! |eta|
DECLARE_SOA_COLUMN(ItsNClsLfTrack0, itsNClsLfTrack0, int);                           //! Number of clusters in ITS
DECLARE_SOA_COLUMN(TpcNClsCrossedRowsLfTrack0, tpcNClsCrossedRowsLfTrack0, int);     //! Number of TPC crossed rows
DECLARE_SOA_COLUMN(ItsChi2NClLfTrack0, itsChi2NClLfTrack0, float);                   //! ITS chi2 / Number of clusters
DECLARE_SOA_COLUMN(TpcChi2NClLfTrack0, tpcChi2NClLfTrack0, float);                   //! TPC chi2 / Number of clusters
DECLARE_SOA_COLUMN(AbsEtaLfTrack0, absEtaLfTrack0, float);                           //! |eta|
DECLARE_SOA_COLUMN(ItsNClsLfTrack1, itsNClsLfTrack1, int);                           //! Number of clusters in ITS
DECLARE_SOA_COLUMN(TpcNClsCrossedRowsLfTrack1, tpcNClsCrossedRowsLfTrack1, int);     //! Number of TPC crossed rows
DECLARE_SOA_COLUMN(ItsChi2NClLfTrack1, itsChi2NClLfTrack1, float);                   //! ITS chi2 / Number of clusters
DECLARE_SOA_COLUMN(TpcChi2NClLfTrack1, tpcChi2NClLfTrack1, float);                   //! TPC chi2 / Number of clusters
DECLARE_SOA_COLUMN(AbsEtaLfTrack1, absEtaLfTrack1, float);                           //! |eta|
DECLARE_SOA_COLUMN(MJpsi, mJpsi, float);                                             //! Invariant mass of Jpsi daughter candidates (GeV/c)
DECLARE_SOA_COLUMN(MPhi, mPhi, float);                                               //! Invariant mass of phi daughter candidates (GeV/c)
DECLARE_SOA_COLUMN(M, m, float);                                                     //! Invariant mass of candidate (GeV/c2)
DECLARE_SOA_COLUMN(Pt, pt, float);                                                   //! Transverse momentum of candidate (GeV/c)
DECLARE_SOA_COLUMN(PtGen, ptGen, float);                                             //! Transverse momentum of candidate (GeV/c)
DECLARE_SOA_COLUMN(P, p, float);                                                     //! Momentum of candidate (GeV/c)
DECLARE_SOA_COLUMN(Y, y, float);                                                     //! Rapidity of candidate
DECLARE_SOA_COLUMN(Eta, eta, float);                                                 //! Pseudorapidity of candidate
DECLARE_SOA_COLUMN(Phi, phi, float);                                                 //! Azimuth angle of candidate
DECLARE_SOA_COLUMN(E, e, float);                                                     //! Energy of candidate (GeV)
DECLARE_SOA_COLUMN(NSigTpcKaBachelor0, nSigTpcKaBachelor0, float);                   //! TPC Nsigma separation for bachelor 0 with kaon mass hypothesis
DECLARE_SOA_COLUMN(NSigTofKaBachelor0, nSigTofKaBachelor0, float);                   //! TOF Nsigma separation for bachelor 0 with kaon mass hypothesis
DECLARE_SOA_COLUMN(NSigTpcTofKaBachelor0, nSigTpcTofKaBachelor0, float);             //! Combined TPC and TOF Nsigma separation for bachelor 0 with kaon mass hypothesis
DECLARE_SOA_COLUMN(NSigTpcKaBachelor1, nSigTpcKaBachelor1, float);                   //! TPC Nsigma separation for bachelor 1 with kaon mass hypothesis
DECLARE_SOA_COLUMN(NSigTofKaBachelor1, nSigTofKaBachelor1, float);                   //! TOF Nsigma separation for bachelor 1 with kaon mass hypothesis
DECLARE_SOA_COLUMN(NSigTpcTofKaBachelor1, nSigTpcTofKaBachelor1, float);             //! Combined TPC and TOF Nsigma separation for bachelor 1 with kaon mass hypothesis
DECLARE_SOA_COLUMN(NSigTpcMuJpsiDauPos, nSigTpcMuJpsiDauPos, float);                 //! TPC Nsigma separation for Jpsi DauPos with muon mass hypothesis
DECLARE_SOA_COLUMN(NSigTofMuJpsiDauPos, nSigTofMuJpsiDauPos, float);                 //! TOF Nsigma separation for Jpsi DauPos with muon mass hypothesis
DECLARE_SOA_COLUMN(NSigTpcTofMuJpsiDauPos, nSigTpcTofMuJpsiDauPos, float);           //! Combined TPC and TOF Nsigma separation for Jpsi prong0 with muon mass hypothesis
DECLARE_SOA_COLUMN(NSigTpcMuJpsiDauNeg, nSigTpcMuJpsiDauNeg, float);                 //! TPC Nsigma separation for Jpsi DauNeg with muon mass hypothesis
DECLARE_SOA_COLUMN(NSigTofMuJpsiDauNeg, nSigTofMuJpsiDauNeg, float);                 //! TOF Nsigma separation for Jpsi DauNeg with muon mass hypothesis
DECLARE_SOA_COLUMN(NSigTpcTofMuJpsiDauNeg, nSigTpcTofMuJpsiDauNeg, float);           //! Combined TPC and TOF Nsigma separation for Jpsi prong1 with muon mass hypothesis
DECLARE_SOA_COLUMN(DecayLength, decayLength, float);                                 //! Decay length of candidate (cm)
DECLARE_SOA_COLUMN(DecayLengthXY, decayLengthXY, float);                             //! Transverse decay length of candidate (cm)
DECLARE_SOA_COLUMN(DecayLengthNormalised, decayLengthNormalised, float);             //! Normalised decay length of candidate
DECLARE_SOA_COLUMN(DecayLengthXYNormalised, decayLengthXYNormalised, float);         //! Normalised transverse decay length of candidate
DECLARE_SOA_COLUMN(CtXY, ctXY, float);                                               //! Pseudo-proper decay length of candidate
DECLARE_SOA_COLUMN(ImpactParameterProduct, impactParameterProduct, float);           //! Impact parameter product of B daughters
DECLARE_SOA_COLUMN(ImpactParameterProductJpsi, impactParameterProductJpsi, float);   //! Impact parameter product of Jpsi daughters
DECLARE_SOA_COLUMN(ImpactParameterProductPhi, impactParameterProductPhi, float);     //! Impact parameter product of Phi daughters
DECLARE_SOA_COLUMN(ImpactParameterJpsiDauPos, impactParameterJpsiDauPos, float);     //! Impact parameter of Jpsi daughter candidate
DECLARE_SOA_COLUMN(ImpactParameterJpsiDauNeg, impactParameterJpsiDauNeg, float);     //! Impact parameter of Jpsi daughter candidate
DECLARE_SOA_COLUMN(ImpactParameterLfTrack0, impactParameterLfTrack0, float);         //! Impact parameter of Phi daughter candidate
DECLARE_SOA_COLUMN(ImpactParameterLfTrack1, impactParameterLfTrack1, float);         //! Impact parameter of Phi daughter candidate
DECLARE_SOA_COLUMN(Cpa, cpa, float);                                                 //! Cosine pointing angle of candidate
DECLARE_SOA_COLUMN(CpaXY, cpaXY, float);                                             //! Cosine pointing angle of candidate in transverse plane
DECLARE_SOA_COLUMN(CpaJpsi, cpaJpsi, float);                                         //! Cosine pointing angle of Jpsi daughter candidate
DECLARE_SOA_COLUMN(CpaXYJpsi, cpaXYJpsi, float);                                     //! Cosine pointing angle in transverse plane of Jpsi daughter candidate
DECLARE_SOA_COLUMN(MaxNormalisedDeltaIP, maxNormalisedDeltaIP, float);               //! Maximum normalized difference between measured and expected impact parameter of candidate prongs
DECLARE_SOA_COLUMN(MlScoreSig, mlScoreSig, float);                                   //! ML score for signal class
DECLARE_SOA_COLUMN(FlagWrongCollision, flagWrongCollision, int8_t);                  //! Flag for association with wrong collision
} // namespace hf_cand_bstojpsiphi_lite

DECLARE_SOA_TABLE(HfRedCandBsLites, "AOD", "HFREDCANDBSLITE", //! Table with some Bs properties
                  hf_cand_bstojpsiphi_lite::M,
                  hf_cand_bstojpsiphi_lite::Pt,
                  hf_cand_bstojpsiphi_lite::Eta,
                  hf_cand_bstojpsiphi_lite::Phi,
                  hf_cand_bstojpsiphi_lite::Y,
                  hf_cand_bstojpsiphi_lite::Cpa,
                  hf_cand_bstojpsiphi_lite::CpaXY,
                  hf_cand::Chi2PCA,
                  hf_cand_bstojpsiphi_lite::DecayLength,
                  hf_cand_bstojpsiphi_lite::DecayLengthXY,
                  hf_cand_bstojpsiphi_lite::DecayLengthNormalised,
                  hf_cand_bstojpsiphi_lite::DecayLengthXYNormalised,
                  hf_cand_bstojpsiphi_lite::CtXY,
                  hf_cand_bstojpsiphi_lite::ImpactParameterProduct,
                  hf_cand_bstojpsiphi_lite::ImpactParameterProductJpsi,
                  hf_cand_bstojpsiphi_lite::ImpactParameterProductPhi,
                  hf_cand_bstojpsiphi_lite::MaxNormalisedDeltaIP,
                  hf_cand_bstojpsiphi_lite::MlScoreSig,
                  // hf_sel_candidate_bplus::IsSelBsToJpsiPi,
                  //  Jpsi meson features
                  hf_cand_bstojpsiphi_lite::MJpsi,
                  hf_cand_bstojpsiphi_lite::PtJpsi,
                  hf_cand_bstojpsiphi_lite::MPhi,
                  hf_cand_bstojpsiphi_lite::ImpactParameterJpsiDauPos,
                  hf_cand_bstojpsiphi_lite::ImpactParameterJpsiDauNeg,
                  hf_cand_bstojpsiphi_lite::ImpactParameterLfTrack0,
                  hf_cand_bstojpsiphi_lite::ImpactParameterLfTrack1,
                  // Jpsi daughter features
                  hf_cand_bstojpsiphi_lite::ItsNClsJpsiDauPos,
                  hf_cand_bstojpsiphi_lite::TpcNClsCrossedRowsJpsiDauPos,
                  hf_cand_bstojpsiphi_lite::ItsChi2NClJpsiDauPos,
                  hf_cand_bstojpsiphi_lite::TpcChi2NClJpsiDauPos,
                  hf_cand_bstojpsiphi_lite::AbsEtaJpsiDauPos,
                  hf_cand_bstojpsiphi_lite::ItsNClsJpsiDauNeg,
                  hf_cand_bstojpsiphi_lite::TpcNClsCrossedRowsJpsiDauNeg,
                  hf_cand_bstojpsiphi_lite::ItsChi2NClJpsiDauNeg,
                  hf_cand_bstojpsiphi_lite::TpcChi2NClJpsiDauNeg,
                  hf_cand_bstojpsiphi_lite::AbsEtaJpsiDauNeg,
                  // kaon features
                  hf_cand_bstojpsiphi_lite::PtBach0,
                  hf_cand_bstojpsiphi_lite::ItsNClsLfTrack0,
                  hf_cand_bstojpsiphi_lite::TpcNClsCrossedRowsLfTrack0,
                  hf_cand_bstojpsiphi_lite::ItsChi2NClLfTrack0,
                  hf_cand_bstojpsiphi_lite::TpcChi2NClLfTrack0,
                  hf_cand_bstojpsiphi_lite::AbsEtaLfTrack0,
                  hf_cand_bstojpsiphi_lite::NSigTpcKaBachelor0,
                  hf_cand_bstojpsiphi_lite::NSigTofKaBachelor0,
                  hf_cand_bstojpsiphi_lite::NSigTpcTofKaBachelor0,
                  hf_cand_bstojpsiphi_lite::PtBach1,
                  hf_cand_bstojpsiphi_lite::ItsNClsLfTrack1,
                  hf_cand_bstojpsiphi_lite::TpcNClsCrossedRowsLfTrack1,
                  hf_cand_bstojpsiphi_lite::ItsChi2NClLfTrack1,
                  hf_cand_bstojpsiphi_lite::TpcChi2NClLfTrack1,
                  hf_cand_bstojpsiphi_lite::AbsEtaLfTrack1,
                  hf_cand_bstojpsiphi_lite::NSigTpcKaBachelor1,
                  hf_cand_bstojpsiphi_lite::NSigTofKaBachelor1,
                  hf_cand_bstojpsiphi_lite::NSigTpcTofKaBachelor1,
                  // MC truth
                  hf_cand_bs::FlagMcMatchRec,
                  hf_cand_bs::FlagMcDecayChanRec,
                  hf_cand_bs::OriginMcRec,
                  hf_cand_bstojpsiphi_lite::FlagWrongCollision,
                  hf_cand_bstojpsiphi_lite::PtGen);

// DECLARE_SOA_TABLE(HfRedBsMcCheck, "AOD", "HFREDBPMCCHECK", //! Table with MC decay type check
//                   hf_cand_2prong::FlagMcMatchRec,
//                   hf_cand_bstojpsiphi_lite::FlagWrongCollision,
//                   hf_cand_bstojpsiphi_lite::MJpsi,
//                   hf_cand_bstojpsiphi_lite::PtJpsi,
//                   hf_cand_bstojpsiphi_lite::M,
//                   hf_cand_bstojpsiphi_lite::Pt,
//                   // hf_cand_bstojpsiphi_lite::MlScoreSig,
//                   hf_bplus_mc::PdgCodeBeautyMother,
//                   hf_bplus_mc::PdgCodeCharmMother,
//                   hf_bplus_mc::PdgCodeDauPos,
//                   hf_bplus_mc::PdgCodeDauNeg,
//                   hf_bplus_mc::PdgCodeProng2);
} // namespace o2::aod

// string definitions, used for histogram axis labels
const TString stringPt = "#it{p}_{T} (GeV/#it{c})";
const TString stringPtJpsi = "#it{p}_{T}(Jpsi) (GeV/#it{c});";
const TString bSCandTitle = "B_{s}^{0} candidates;";
const TString entries = "entries";
const TString bSCandMatch = "B_{s}^{0} candidates (matched);";
const TString bSCandUnmatch = "B_{s}^{0} candidates (unmatched);";
const TString mcParticleMatched = "MC particles (matched);";

/// Bs analysis task
struct HfTaskBsToJpsiPhiReduced {
  Produces<aod::HfRedCandBsLites> hfRedCandBsLite;
  // Produces<aod::HfRedBsMcCheck> hfRedBsMcCheck;

  Configurable<uint8_t> selectionFlagBs{"selectionFlagBs", 1, "Selection Flag for Bs"};
  Configurable<double> yCandGenMax{"yCandGenMax", 0.5, "max. gen particle rapidity"};
  Configurable<double> yCandRecoMax{"yCandRecoMax", 0.8, "max. cand. rapidity"};
  Configurable<float> etaTrackMax{"etaTrackMax", 0.8, "max. track pseudo-rapidity"};
  Configurable<float> ptTrackMin{"ptTrackMin", 0.1, "min. track transverse momentum"};
  Configurable<bool> fillBackground{"fillBackground", false, "Flag to enable filling of background histograms/sparses/tree (only MC)"};
  Configurable<float> downSampleBkgFactor{"downSampleBkgFactor", 1., "Fraction of background candidates to keep for ML trainings"};
  Configurable<float> ptMaxForDownSample{"ptMaxForDownSample", 10., "Maximum pt for the application of the downsampling factor"};
  // topological cuts
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_bs_to_jpsi_phi::vecBinsPt}, "pT bin limits"};
  Configurable<LabeledArray<double>> cuts{"cuts", {hf_cuts_bs_to_jpsi_phi::Cuts[0], hf_cuts_bs_to_jpsi_phi::NBinsPt, hf_cuts_bs_to_jpsi_phi::NCutVars, hf_cuts_bs_to_jpsi_phi::labelsPt, hf_cuts_bs_to_jpsi_phi::labelsCutVar}, "Bs candidate selection per pT bin"};
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
  // Bs ML inference
  Configurable<std::vector<double>> binsPtBsMl{"binsPtBsMl", std::vector<double>{hf_cuts_ml::vecBinsPt}, "pT bin limits for ML application"};
  Configurable<std::vector<int>> cutDirBsMl{"cutDirBsMl", std::vector<int>{hf_cuts_ml::vecCutDir}, "Whether to reject score values greater or smaller than the threshold"};
  Configurable<LabeledArray<double>> cutsBsMl{"cutsBsMl", {hf_cuts_ml::Cuts[0], hf_cuts_ml::NBinsPt, hf_cuts_ml::NCutScores, hf_cuts_ml::labelsPt, hf_cuts_ml::labelsCutScore}, "ML selections per pT bin"};
  Configurable<int> nClassesBsMl{"nClassesBsMl", static_cast<int>(hf_cuts_ml::NCutScores), "Number of classes in ML model"};
  Configurable<std::vector<std::string>> namesInputFeatures{"namesInputFeatures", std::vector<std::string>{"feature1", "feature2"}, "Names of ML model input features"};
  // CCDB configuration
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::vector<std::string>> modelPathsCCDB{"modelPathsCCDB", std::vector<std::string>{"path_ccdb/BDT_BS/"}, "Paths of models on CCDB"};
  Configurable<std::vector<std::string>> onnxFileNames{"onnxFileNames", std::vector<std::string>{"ModelHandler_onnx_BSToJpsiPhi.onnx"}, "ONNX file names for each pT bin (if not from CCDB full path)"};
  Configurable<int64_t> timestampCCDB{"timestampCCDB", -1, "timestamp of the ONNX file for ML model used to query in CCDB"};
  Configurable<bool> loadModelsFromCCDB{"loadModelsFromCCDB", false, "Flag to enable or disable the loading of models from CCDB"};

  TrackSelectorKa selectorKaon;
  o2::analysis::HfMlResponseBsToJpsiPhiReduced<float> hfMlResponse;
  o2::ccdb::CcdbApi ccdbApi;

  using TracksKaon = soa::Join<HfRedTracks, HfRedTracksPid>;
  std::vector<float> outputMl;

  // Filter filterSelectCandidates = (aod::hf_sel_candidate_bplus::isSelBsToJpsiPi >= selectionFlagBs);

  HistogramRegistry registry{"registry"};

  void init(InitContext&)
  {
    std::array<bool, 2> processFuncData{doprocessData, doprocessDataWithBsMl};
    if ((std::accumulate(processFuncData.begin(), processFuncData.end(), 0)) > 1) {
      LOGP(fatal, "Only one process function for data can be enabled at a time.");
    }
    std::array<bool, 2> processFuncMc{doprocessMc, doprocessMcWithBsMl};
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

    const AxisSpec axisMassBs{150, 4.5, 6.0};
    const AxisSpec axisMassJpsi{600, 2.8f, 3.4f};
    const AxisSpec axisMassPhi{200, 0.9f, 1.1f};
    const AxisSpec axisPtProng{100, 0., 10.};
    const AxisSpec axisImpactPar{200, -0.05, 0.05};
    const AxisSpec axisPtJpsi{100, 0., 50.};
    const AxisSpec axisRapidity{100, -2., 2.};
    const AxisSpec axisPtB{(std::vector<double>)binsPt, "#it{p}_{T}^{B_{s}^{0}} (GeV/#it{c})"};
    const AxisSpec axisPtKa{100, 0.f, 10.f};
    const AxisSpec axisPtPhi{100, 0.f, 10.f};

    registry.add("hMass", bSCandTitle + "inv. mass J/#Psi K^{+} (GeV/#it{c}^{2});" + stringPt, {HistType::kTH2F, {axisMassBs, axisPtB}});
    registry.add("hMassJpsi", bSCandTitle + "inv. mass #mu^{+}#mu^{#minus} (GeV/#it{c}^{2});" + stringPt, {HistType::kTH2F, {axisMassJpsi, axisPtJpsi}});
    registry.add("hMassPhi", bSCandTitle + "inv. mass K^{+}K^{#minus} (GeV/#it{c}^{2});" + stringPt, {HistType::kTH2F, {axisMassPhi, axisPtPhi}});
    registry.add("hd0K", bSCandTitle + "Kaon DCAxy to prim. vertex (cm);" + stringPt, {HistType::kTH2F, {axisImpactPar, axisPtKa}});

    // histograms processMC
    if (doprocessMc || doprocessMcWithBsMl) {
      registry.add("hPtJpsiGen", mcParticleMatched + "J/#Psi #it{p}_{T}^{gen} (GeV/#it{c}); B_{s}^{0} " + stringPt, {HistType::kTH2F, {axisPtProng, axisPtB}});
      registry.add("hPtPhiGen", mcParticleMatched + "#phi #it{p}_{T}^{gen} (GeV/#it{c}); B_{s}^{0} " + stringPt, {HistType::kTH2F, {axisPtProng, axisPtB}});
      registry.add("hPtKGen", mcParticleMatched + "Kaon #it{p}_{T}^{gen} (GeV/#it{c}); B_{s}^{0} " + stringPt, {HistType::kTH2F, {axisPtProng, axisPtB}});
      registry.add("hYGenWithProngsInAcceptance", mcParticleMatched + "B_{s}^{0} #it{p}_{T}^{gen} (GeV/#it{c}); B_{s}^{0} #it{y}", {HistType::kTH2F, {axisPtProng, axisRapidity}});
      registry.add("hMassRecSig", bSCandMatch + "inv. mass J/#Psi K^{+} (GeV/#it{c}^{2}); B_{s}^{0} " + stringPt, {HistType::kTH2F, {axisMassBs, axisPtB}});
      registry.add("hMassJpsiRecSig", bSCandMatch + "inv. mass #mu^{+}#mu^{#minus} (GeV/#it{c}^{2}); J/#Psi " + stringPt, {HistType::kTH2F, {axisMassJpsi, axisPtJpsi}});
      registry.add("hMassPhiRecSig", bSCandMatch + "inv. mass K^{+}K^{#minus} (GeV/#it{c}^{2}); #phi " + stringPt, {HistType::kTH2F, {axisMassPhi, axisPtPhi}});
      registry.add("hd0KRecSig", bSCandMatch + "Kaon DCAxy to prim. vertex (cm); K^{+} " + stringPt, {HistType::kTH2F, {axisImpactPar, axisPtKa}});
      registry.add("hMassRecBg", bSCandUnmatch + "inv. mass J/#Psi K^{+} (GeV/#it{c}^{2}); B_{s}^{0} " + stringPt, {HistType::kTH2F, {axisMassBs, axisPtB}});
      registry.add("hMassJpsiRecBg", bSCandUnmatch + "inv. mass #mu^{+}#mu^{#minus} (GeV/#it{c}^{2}); J/#Psi " + stringPt, {HistType::kTH2F, {axisMassJpsi, axisPtJpsi}});
      registry.add("hMassPhiRecBg", bSCandMatch + "inv. mass K^{+}K^{#minus} (GeV/#it{c}^{2}); #phi " + stringPt, {HistType::kTH2F, {axisMassPhi, axisPtPhi}});
      registry.add("hd0KRecBg", bSCandMatch + "Kaon DCAxy to prim. vertex (cm); K^{+} " + stringPt, {HistType::kTH2F, {axisImpactPar, axisPtKa}});
    }

    if (doprocessDataWithBsMl || doprocessMcWithBsMl) {
      hfMlResponse.configure(binsPtBsMl, cutsBsMl, cutDirBsMl, nClassesBsMl);
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

  /// Selection of Bs daughter in geometrical acceptance
  /// \param etaProng is the pseudorapidity of Bs prong
  /// \param ptProng is the pT of Bs prong
  /// \return true if prong is in geometrical acceptance
  template <typename T = float>
  bool isProngInAcceptance(const T& etaProng, const T& ptProng)
  {
    return std::abs(etaProng) <= etaTrackMax && ptProng >= ptTrackMin;
  }

  /// Calculate pseudorapidity from track tan(lambda)
  /// \param tgl is the track tangent of the dip angle
  /// \return pseudorapidity
  float absEta(float tgl)
  {
    return std::abs(std::log(std::tan(o2::constants::math::PIQuarter - 0.5f * std::atan(tgl))));
  }

  /// Fill candidate information at reconstruction level
  /// \param doMc is the flag to enable the filling with MC information
  /// \param withBsMl is the flag to enable the filling with ML scores for the Bs candidate
  /// \param candidate is the Bs candidate
  /// \param candidatesJpsi is the table with Jpsi candidates
  template <bool DoMc, bool WithBsMl, typename Cand>
  void fillCand(Cand const& candidate,
                aod::HfRedJpsis const& /*candidatesJpsi*/,
                aod::HfRedBach0Tracks const&,
                aod::HfRedBach1Tracks const&)
  {
    auto ptCandBs = candidate.pt();
    auto invMassBs = HfHelper::invMassBsToJpsiPhi(candidate);
    auto candJpsi = candidate.template jpsi_as<aod::HfRedJpsis>();
    auto candKa0 = candidate.template prong0Phi_as<aod::HfRedBach0Tracks>();
    auto candKa1 = candidate.template prong1Phi_as<aod::HfRedBach1Tracks>();
    std::array<float, 3> const pVecKa0 = {candKa0.px(), candKa0.py(), candKa0.pz()};
    std::array<float, 3> const pVecKa1 = {candKa1.px(), candKa1.py(), candKa1.pz()};
    auto ptJpsi = candidate.ptProng0();
    auto invMassJpsi = candJpsi.m();
    auto invMassPhi = RecoDecay::m(std::array{pVecKa0, pVecKa1}, std::array{o2::constants::physics::MassKPlus, o2::constants::physics::MassKPlus});
    uint8_t statusBs = 0;

    int8_t flagMcMatchRec{0}, flagMcDecayChanRec{0}, flagWrongCollision{0};
    bool isSignal = false;
    if constexpr (DoMc) {
      flagMcMatchRec = candidate.flagMcMatchRec();
      flagMcDecayChanRec = candidate.flagMcDecayChanRec();
      flagWrongCollision = candidate.flagWrongCollision();
      isSignal = flagMcMatchRec == o2::hf_decay::hf_cand_beauty::BsToJpsiKK &&
                 flagMcDecayChanRec == o2::hf_decay::hf_cand_beauty::BsToJpsiPhi;
    }

    SETBIT(statusBs, SelectionStep::RecoSkims);
    if (HfHelper::selectionBsToJpsiPhiTopol(candidate, candKa0, candKa1, cuts, binsPt)) {
      SETBIT(statusBs, SelectionStep::RecoTopol);
    } else if (selectionFlagBs >= BIT(SelectionStep::RecoTopol) * 2 - 1) {
      return;
    }
    // track-level PID selection
    // auto trackKa = candidate.template prong1_as<TracksPion>();
    if (kaonPidMethod == PidMethod::TpcOrTof || kaonPidMethod == PidMethod::TpcAndTof) {
      int pidTrackKa0{TrackSelectorPID::Status::NotApplicable};
      int pidTrackKa1{TrackSelectorPID::Status::NotApplicable};
      if (kaonPidMethod == PidMethod::TpcOrTof) {
        pidTrackKa0 = selectorKaon.statusTpcOrTof(candKa0);
        pidTrackKa1 = selectorKaon.statusTpcOrTof(candKa1);
      } else if (kaonPidMethod == PidMethod::TpcAndTof) {
        pidTrackKa0 = selectorKaon.statusTpcAndTof(candKa0);
        pidTrackKa1 = selectorKaon.statusTpcAndTof(candKa1);
      }
      if (HfHelper::selectionBsToJpsiPhiPid(pidTrackKa0, acceptPIDNotApplicable.value) &&
          HfHelper::selectionBsToJpsiPhiPid(pidTrackKa1, acceptPIDNotApplicable.value)) {
        // LOGF(info, "Bs candidate selection failed at PID selection");
        SETBIT(statusBs, SelectionStep::RecoPID);
      } else if (selectionFlagBs >= BIT(SelectionStep::RecoPID) * 2 - 1) {
        return;
      }
    }

    float candidateMlScoreSig = -1;
    if constexpr (WithBsMl) {
      // Bs ML selections
      std::vector<float> inputFeatures = hfMlResponse.getInputFeatures(candidate, candKa0, candKa1);
      if (hfMlResponse.isSelectedMl(inputFeatures, ptCandBs, outputMl)) {
        SETBIT(statusBs, SelectionStep::RecoMl);
      } else if (selectionFlagBs >= BIT(SelectionStep::RecoMl) * 2 - 1) {
        return;
      }
      candidateMlScoreSig = outputMl[1];
    }

    registry.fill(HIST("hMass"), invMassBs, ptCandBs);
    registry.fill(HIST("hMassJpsi"), invMassJpsi, candidate.ptProng0());
    registry.fill(HIST("hMassPhi"), invMassPhi, candidate.ptProng0());
    registry.fill(HIST("hd0K"), candidate.impactParameter1(), candidate.ptProng1());
    if constexpr (DoMc) {
      if (isSignal) {
        registry.fill(HIST("hMassRecSig"), invMassBs, ptCandBs);
        registry.fill(HIST("hMassJpsiRecSig"), invMassJpsi, candidate.ptProng0());
        registry.fill(HIST("hd0KRecSig"), candidate.impactParameter1(), candidate.ptProng1());
      } else if (fillBackground) {
        registry.fill(HIST("hMassRecBg"), invMassBs, ptCandBs);
        registry.fill(HIST("hMassJpsiRecBg"), invMassJpsi, candidate.ptProng0());
        registry.fill(HIST("hd0KRecBg"), candidate.impactParameter1(), candidate.ptProng1());
      }
    }

    float const pseudoRndm = ptJpsi * 1000. - static_cast<int64_t>(ptJpsi * 1000);
    if (ptCandBs >= ptMaxForDownSample || pseudoRndm < downSampleBkgFactor) {
      float ptMother = -1.;
      if constexpr (DoMc) {
        ptMother = candidate.ptMother();
      }

      hfRedCandBsLite(
        // Bs - meson features
        invMassBs,
        ptCandBs,
        candidate.eta(),
        candidate.phi(),
        HfHelper::yBs(candidate),
        candidate.cpa(),
        candidate.cpaXY(),
        candidate.chi2PCA(),
        candidate.decayLength(),
        candidate.decayLengthXY(),
        candidate.decayLengthNormalised(),
        candidate.decayLengthXYNormalised(),
        candidate.ctXY(std::array{o2::constants::physics::MassMuon, o2::constants::physics::MassMuon, o2::constants::physics::MassKPlus, o2::constants::physics::MassKPlus}),
        candidate.impactParameterProduct(),
        candidate.impactParameterProductJpsi(),
        candidate.impactParameterProductPhi(),
        candidate.maxNormalisedDeltaIP(),
        candidateMlScoreSig,
        // J/Psi features
        invMassJpsi,
        ptJpsi,
        invMassPhi,
        candidate.impactParameter0(),
        candidate.impactParameter1(),
        candidate.impactParameter2(),
        candidate.impactParameter3(),
        candJpsi.itsNClsDauPos(),
        candJpsi.tpcNClsCrossedRowsDauPos(),
        candJpsi.itsChi2NClDauPos(),
        candJpsi.tpcChi2NClDauPos(),
        absEta(candJpsi.tglDauPos()),
        candJpsi.itsNClsDauNeg(),
        candJpsi.tpcNClsCrossedRowsDauNeg(),
        candJpsi.itsChi2NClDauNeg(),
        candJpsi.tpcChi2NClDauNeg(),
        absEta(candJpsi.tglDauNeg()),
        // kaon features
        candKa0.pt(),
        candKa0.itsNCls(),
        candKa0.tpcNClsCrossedRows(),
        candKa0.itsChi2NCl(),
        candKa0.tpcChi2NCl(),
        absEta(candKa0.tgl()),
        // candKa.absEtaBach(candKa.tgl()),
        candKa0.tpcNSigmaKa(),
        candKa0.tofNSigmaKa(),
        candKa0.tpcTofNSigmaKa(),
        candKa1.pt(),
        candKa1.itsNCls(),
        candKa1.tpcNClsCrossedRows(),
        candKa1.itsChi2NCl(),
        candKa1.tpcChi2NCl(),
        absEta(candKa1.tgl()),
        candKa1.tpcNSigmaKa(),
        candKa1.tofNSigmaKa(),
        candKa1.tpcTofNSigmaKa(),
        // MC truth
        flagMcMatchRec,
        flagMcDecayChanRec,
        isSignal,
        flagWrongCollision,
        ptMother);
    }
  }

  /// Fill particle histograms (gen MC truth)
  void fillCandMcGen(aod::HfMcGenRedBss::iterator const& particle)
  {
    auto ptParticle = particle.ptTrack();
    auto yParticle = particle.yTrack();
    if (yCandGenMax >= 0. && std::abs(yParticle) > yCandGenMax) {
      return;
    }
    std::array<float, 2> ptProngs = {particle.ptProng0(), particle.ptProng1()};
    std::array<float, 2> etaProngs = {particle.etaProng0(), particle.etaProng1()};
    bool const prongsInAcc = isProngInAcceptance(etaProngs[0], ptProngs[0]) && isProngInAcceptance(etaProngs[1], ptProngs[1]);

    registry.fill(HIST("hPtJpsiGen"), ptProngs[0], ptParticle);
    registry.fill(HIST("hPtKGen"), ptProngs[1], ptParticle);

    // generated Bs with daughters in geometrical acceptance
    if (prongsInAcc) {
      registry.fill(HIST("hYGenWithProngsInAcceptance"), ptParticle, yParticle);
    }
  }

  // Process functions
  void processData(aod::HfRedCandBsToJpsiPhi const& candidates,
                   aod::HfRedJpsis const& candidatesJpsi,
                   aod::HfRedBach0Tracks const& kaon0Tracks,
                   aod::HfRedBach1Tracks const& kaon1Tracks)
  {
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(HfHelper::yBs(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCand<false, false>(candidate, candidatesJpsi, kaon0Tracks, kaon1Tracks);
    } // candidate loop
  } // processData
  PROCESS_SWITCH(HfTaskBsToJpsiPhiReduced, processData, "Process data without ML for Bs", true);

  void processDataWithBsMl(aod::HfRedCandBsToJpsiPhi const& candidates,
                           aod::HfRedJpsis const& candidatesJpsi,
                           aod::HfRedBach0Tracks const& kaon0Tracks,
                           aod::HfRedBach1Tracks const& kaon1Tracks)
  {
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(HfHelper::yBs(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCand<false, true>(candidate, candidatesJpsi, kaon0Tracks, kaon1Tracks);
    } // candidate loop
  } // processDataWithBsMl
  PROCESS_SWITCH(HfTaskBsToJpsiPhiReduced, processDataWithBsMl, "Process data with ML for Bs", false);

  void processMc(soa::Join<aod::HfRedCandBsToJpsiPhi, aod::HfMcRecRedBss> const& candidates,
                 aod::HfMcGenRedBss const& mcParticles,
                 aod::HfRedJpsis const& candidatesJpsi,
                 aod::HfRedBach0Tracks const& kaon0Tracks,
                 aod::HfRedBach1Tracks const& kaon1Tracks)
  {
    // MC rec
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(HfHelper::yBs(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCand<true, false>(candidate, candidatesJpsi, kaon0Tracks, kaon1Tracks);
    } // rec

    // MC gen. level
    for (const auto& particle : mcParticles) {
      fillCandMcGen(particle);
    } // gen
  } // processMc
  PROCESS_SWITCH(HfTaskBsToJpsiPhiReduced, processMc, "Process MC without ML for Bs", false);

  void processMcWithBsMl(soa::Join<aod::HfRedCandBsToJpsiPhi, aod::HfMcRecRedBss> const& candidates,
                         aod::HfMcGenRedBss const& mcParticles,
                         aod::HfRedJpsis const& candidatesJpsi,
                         aod::HfRedBach0Tracks const& kaon0Tracks,
                         aod::HfRedBach1Tracks const& kaon1Tracks)
  {
    // MC rec
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(HfHelper::yBs(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCand<true, true>(candidate, candidatesJpsi, kaon0Tracks, kaon1Tracks);
    } // rec

    // MC gen. level
    for (const auto& particle : mcParticles) {
      fillCandMcGen(particle);
    } // gen
  } // processMcWithBsMl
  PROCESS_SWITCH(HfTaskBsToJpsiPhiReduced, processMcWithBsMl, "Process MC with ML for Bs", false);

}; // struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskBsToJpsiPhiReduced>(cfgc)};
}
