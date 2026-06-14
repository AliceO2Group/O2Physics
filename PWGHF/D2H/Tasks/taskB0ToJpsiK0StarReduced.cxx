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

/// \file taskB0ToJpsiK0StarReduced.cxx
/// \brief B0 → Jpsi K*0 → (µ+ µ-) (K+pi-) analysis task
///
/// \author Fabrizio Chinu <fabrizio.chinu@cern.ch>, Università degli Studi and INFN Torino
/// \author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN

#include "PWGHF/Core/DecayChannels.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/HfMlResponseB0ToJpsiK0StarReduced.h"
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
namespace hf_cand_b0tojpsik0star_lite
{
DECLARE_SOA_COLUMN(PtJpsi, ptJpsi, float);                                             //! Transverse momentum of Jpsi daughter candidate (GeV/c)
DECLARE_SOA_COLUMN(PtBach0, ptBach0, float);                                           //! Transverse momentum of bachelor kaon(<- K*0) (GeV/c)
DECLARE_SOA_COLUMN(PtBach1, ptBach1, float);                                           //! Transverse momentum of bachelor kaon(<- K*0) (GeV/c)
DECLARE_SOA_COLUMN(ItsNClsJpsiDauPos, itsNClsJpsiDauPos, int);                         //! Number of clusters in ITS
DECLARE_SOA_COLUMN(TpcNClsCrossedRowsJpsiDauPos, tpcNClsCrossedRowsJpsiDauPos, int);   //! Number of TPC crossed rows
DECLARE_SOA_COLUMN(ItsChi2NClJpsiDauPos, itsChi2NClJpsiDauPos, float);                 //! ITS chi2 / Number of clusters
DECLARE_SOA_COLUMN(TpcChi2NClJpsiDauPos, tpcChi2NClJpsiDauPos, float);                 //! TPC chi2 / Number of clusters
DECLARE_SOA_COLUMN(AbsEtaJpsiDauPos, absEtaJpsiDauPos, float);                         //! |eta|
DECLARE_SOA_COLUMN(ItsNClsJpsiDauNeg, itsNClsJpsiDauNeg, int);                         //! Number of clusters in ITS
DECLARE_SOA_COLUMN(TpcNClsCrossedRowsJpsiDauNeg, tpcNClsCrossedRowsJpsiDauNeg, int);   //! Number of TPC crossed rows
DECLARE_SOA_COLUMN(ItsChi2NClJpsiDauNeg, itsChi2NClJpsiDauNeg, float);                 //! ITS chi2 / Number of clusters
DECLARE_SOA_COLUMN(TpcChi2NClJpsiDauNeg, tpcChi2NClJpsiDauNeg, float);                 //! TPC chi2 / Number of clusters
DECLARE_SOA_COLUMN(AbsEtaJpsiDauNeg, absEtaJpsiDauNeg, float);                         //! |eta|
DECLARE_SOA_COLUMN(ItsNClsLfTrack0, itsNClsLfTrack0, int);                             //! Number of clusters in ITS
DECLARE_SOA_COLUMN(TpcNClsCrossedRowsLfTrack0, tpcNClsCrossedRowsLfTrack0, int);       //! Number of TPC crossed rows
DECLARE_SOA_COLUMN(ItsChi2NClLfTrack0, itsChi2NClLfTrack0, float);                     //! ITS chi2 / Number of clusters
DECLARE_SOA_COLUMN(TpcChi2NClLfTrack0, tpcChi2NClLfTrack0, float);                     //! TPC chi2 / Number of clusters
DECLARE_SOA_COLUMN(AbsEtaLfTrack0, absEtaLfTrack0, float);                             //! |eta|
DECLARE_SOA_COLUMN(ItsNClsLfTrack1, itsNClsLfTrack1, int);                             //! Number of clusters in ITS
DECLARE_SOA_COLUMN(TpcNClsCrossedRowsLfTrack1, tpcNClsCrossedRowsLfTrack1, int);       //! Number of TPC crossed rows
DECLARE_SOA_COLUMN(ItsChi2NClLfTrack1, itsChi2NClLfTrack1, float);                     //! ITS chi2 / Number of clusters
DECLARE_SOA_COLUMN(TpcChi2NClLfTrack1, tpcChi2NClLfTrack1, float);                     //! TPC chi2 / Number of clusters
DECLARE_SOA_COLUMN(AbsEtaLfTrack1, absEtaLfTrack1, float);                             //! |eta|
DECLARE_SOA_COLUMN(MJpsi, mJpsi, float);                                               //! Invariant mass of Jpsi daughter candidates (GeV/c)
DECLARE_SOA_COLUMN(MK0Star, mK0Star, float);                                           //! Invariant mass of K*0 daughter candidates (GeV/c)
DECLARE_SOA_COLUMN(M, m, float);                                                       //! Invariant mass of candidate (GeV/c2)
DECLARE_SOA_COLUMN(Pt, pt, float);                                                     //! Transverse momentum of candidate (GeV/c)
DECLARE_SOA_COLUMN(PtGen, ptGen, float);                                               //! Transverse momentum of candidate (GeV/c)
DECLARE_SOA_COLUMN(P, p, float);                                                       //! Momentum of candidate (GeV/c)
DECLARE_SOA_COLUMN(Y, y, float);                                                       //! Rapidity of candidate
DECLARE_SOA_COLUMN(Eta, eta, float);                                                   //! Pseudorapidity of candidate
DECLARE_SOA_COLUMN(Phi, phi, float);                                                   //! Azimuth angle of candidate
DECLARE_SOA_COLUMN(E, e, float);                                                       //! Energy of candidate (GeV)
DECLARE_SOA_COLUMN(NSigTpcKaBachelor0, nSigTpcKaBachelor0, float);                     //! TPC Nsigma separation for bachelor 0 with kaon mass hypothesis
DECLARE_SOA_COLUMN(NSigTofKaBachelor0, nSigTofKaBachelor0, float);                     //! TOF Nsigma separation for bachelor 0 with kaon mass hypothesis
DECLARE_SOA_COLUMN(NSigTpcTofKaBachelor0, nSigTpcTofKaBachelor0, float);               //! Combined TPC and TOF Nsigma separation for bachelor 0 with kaon mass hypothesis
DECLARE_SOA_COLUMN(NSigTpcKaBachelor1, nSigTpcKaBachelor1, float);                     //! TPC Nsigma separation for bachelor 1 with kaon mass hypothesis
DECLARE_SOA_COLUMN(NSigTofKaBachelor1, nSigTofKaBachelor1, float);                     //! TOF Nsigma separation for bachelor 1 with kaon mass hypothesis
DECLARE_SOA_COLUMN(NSigTpcTofKaBachelor1, nSigTpcTofKaBachelor1, float);               //! Combined TPC and TOF Nsigma separation for bachelor 1 with kaon mass hypothesis
DECLARE_SOA_COLUMN(NSigTpcPiBachelor0, nSigTpcPiBachelor0, float);                     //! TPC Nsigma separation for bachelor 0 with pion mass hypothesis
DECLARE_SOA_COLUMN(NSigTofPiBachelor0, nSigTofPiBachelor0, float);                     //! TOF Nsigma separation for bachelor 0 with pion mass hypothesis
DECLARE_SOA_COLUMN(NSigTpcTofPiBachelor0, nSigTpcTofPiBachelor0, float);               //! Combined TPC and TOF Nsigma separation for bachelor 0 with pion mass hypothesis
DECLARE_SOA_COLUMN(NSigTpcPiBachelor1, nSigTpcPiBachelor1, float);                     //! TPC Nsigma separation for bachelor 1 with pion mass hypothesis
DECLARE_SOA_COLUMN(NSigTofPiBachelor1, nSigTofPiBachelor1, float);                     //! TOF Nsigma separation for bachelor 1 with pion mass hypothesis
DECLARE_SOA_COLUMN(NSigTpcTofPiBachelor1, nSigTpcTofPiBachelor1, float);               //! Combined TPC and TOF Nsigma separation for bachelor 1 with pion mass hypothesis
DECLARE_SOA_COLUMN(NSigTpcMuJpsiDauPos, nSigTpcMuJpsiDauPos, float);                   //! TPC Nsigma separation for Jpsi DauPos with muon mass hypothesis
DECLARE_SOA_COLUMN(NSigTofMuJpsiDauPos, nSigTofMuJpsiDauPos, float);                   //! TOF Nsigma separation for Jpsi DauPos with muon mass hypothesis
DECLARE_SOA_COLUMN(NSigTpcTofMuJpsiDauPos, nSigTpcTofMuJpsiDauPos, float);             //! Combined TPC and TOF Nsigma separation for Jpsi prong0 with muon mass hypothesis
DECLARE_SOA_COLUMN(NSigTpcMuJpsiDauNeg, nSigTpcMuJpsiDauNeg, float);                   //! TPC Nsigma separation for Jpsi DauNeg with muon mass hypothesis
DECLARE_SOA_COLUMN(NSigTofMuJpsiDauNeg, nSigTofMuJpsiDauNeg, float);                   //! TOF Nsigma separation for Jpsi DauNeg with muon mass hypothesis
DECLARE_SOA_COLUMN(NSigTpcTofMuJpsiDauNeg, nSigTpcTofMuJpsiDauNeg, float);             //! Combined TPC and TOF Nsigma separation for Jpsi prong1 with muon mass hypothesis
DECLARE_SOA_COLUMN(DecayLength, decayLength, float);                                   //! Decay length of candidate (cm)
DECLARE_SOA_COLUMN(DecayLengthXY, decayLengthXY, float);                               //! Transverse decay length of candidate (cm)
DECLARE_SOA_COLUMN(DecayLengthNormalised, decayLengthNormalised, float);               //! Normalised decay length of candidate
DECLARE_SOA_COLUMN(DecayLengthXYNormalised, decayLengthXYNormalised, float);           //! Normalised transverse decay length of candidate
DECLARE_SOA_COLUMN(CtXY, ctXY, float);                                                 //! Pseudo-proper decay length of candidate
DECLARE_SOA_COLUMN(ImpactParameterProduct, impactParameterProduct, float);             //! Impact parameter product of B daughters
DECLARE_SOA_COLUMN(ImpactParameterProductJpsi, impactParameterProductJpsi, float);     //! Impact parameter product of Jpsi daughters
DECLARE_SOA_COLUMN(ImpactParameterProductK0Star, impactParameterProductK0Star, float); //! Impact parameter product of K*0 daughters
DECLARE_SOA_COLUMN(ImpactParameterJpsiDauPos, impactParameterJpsiDauPos, float);       //! Impact parameter of Jpsi daughter candidate
DECLARE_SOA_COLUMN(ImpactParameterJpsiDauNeg, impactParameterJpsiDauNeg, float);       //! Impact parameter of Jpsi daughter candidate
DECLARE_SOA_COLUMN(ImpactParameterLfTrack0, impactParameterLfTrack0, float);           //! Impact parameter of K*0 daughter candidate
DECLARE_SOA_COLUMN(ImpactParameterLfTrack1, impactParameterLfTrack1, float);           //! Impact parameter of K*0 daughter candidate
DECLARE_SOA_COLUMN(Cpa, cpa, float);                                                   //! Cosine pointing angle of candidate
DECLARE_SOA_COLUMN(CpaXY, cpaXY, float);                                               //! Cosine pointing angle of candidate in transverse plane
DECLARE_SOA_COLUMN(CpaJpsi, cpaJpsi, float);                                           //! Cosine pointing angle of Jpsi daughter candidate
DECLARE_SOA_COLUMN(CpaXYJpsi, cpaXYJpsi, float);                                       //! Cosine pointing angle in transverse plane of Jpsi daughter candidate
DECLARE_SOA_COLUMN(MaxNormalisedDeltaIP, maxNormalisedDeltaIP, float);                 //! Maximum normalized difference between measured and expected impact parameter of candidate prongs
DECLARE_SOA_COLUMN(MlScoreSig, mlScoreSig, float);                                     //! ML score for signal class
DECLARE_SOA_COLUMN(FlagWrongCollision, flagWrongCollision, int8_t);                    //! Flag for association with wrong collision
} // namespace hf_cand_b0tojpsik0star_lite

DECLARE_SOA_TABLE(HfRedCandB0Lites, "AOD", "HFREDCANDB0LITE", //! Table with some B0 properties
                  hf_cand_b0tojpsik0star_lite::M,
                  hf_cand_b0tojpsik0star_lite::Pt,
                  hf_cand_b0tojpsik0star_lite::Eta,
                  hf_cand_b0tojpsik0star_lite::Phi,
                  hf_cand_b0tojpsik0star_lite::Y,
                  hf_cand_b0tojpsik0star_lite::Cpa,
                  hf_cand_b0tojpsik0star_lite::CpaXY,
                  hf_cand::Chi2PCA,
                  hf_cand_b0tojpsik0star_lite::DecayLength,
                  hf_cand_b0tojpsik0star_lite::DecayLengthXY,
                  hf_cand_b0tojpsik0star_lite::DecayLengthNormalised,
                  hf_cand_b0tojpsik0star_lite::DecayLengthXYNormalised,
                  hf_cand_b0tojpsik0star_lite::CtXY,
                  hf_cand_b0tojpsik0star_lite::ImpactParameterProduct,
                  hf_cand_b0tojpsik0star_lite::ImpactParameterProductJpsi,
                  hf_cand_b0tojpsik0star_lite::ImpactParameterProductK0Star,
                  hf_cand_b0tojpsik0star_lite::MaxNormalisedDeltaIP,
                  hf_cand_b0tojpsik0star_lite::MlScoreSig,
                  // hf_sel_candidate_bplus::IsSelBsToJpsiPi,
                  //  Jpsi meson features
                  hf_cand_b0tojpsik0star_lite::MJpsi,
                  hf_cand_b0tojpsik0star_lite::PtJpsi,
                  hf_cand_b0tojpsik0star_lite::MK0Star,
                  hf_cand_b0tojpsik0star_lite::ImpactParameterJpsiDauPos,
                  hf_cand_b0tojpsik0star_lite::ImpactParameterJpsiDauNeg,
                  hf_cand_b0tojpsik0star_lite::ImpactParameterLfTrack0,
                  hf_cand_b0tojpsik0star_lite::ImpactParameterLfTrack1,
                  // Jpsi daughter features
                  hf_cand_b0tojpsik0star_lite::ItsNClsJpsiDauPos,
                  hf_cand_b0tojpsik0star_lite::TpcNClsCrossedRowsJpsiDauPos,
                  hf_cand_b0tojpsik0star_lite::ItsChi2NClJpsiDauPos,
                  hf_cand_b0tojpsik0star_lite::TpcChi2NClJpsiDauPos,
                  hf_cand_b0tojpsik0star_lite::AbsEtaJpsiDauPos,
                  hf_cand_b0tojpsik0star_lite::ItsNClsJpsiDauNeg,
                  hf_cand_b0tojpsik0star_lite::TpcNClsCrossedRowsJpsiDauNeg,
                  hf_cand_b0tojpsik0star_lite::ItsChi2NClJpsiDauNeg,
                  hf_cand_b0tojpsik0star_lite::TpcChi2NClJpsiDauNeg,
                  hf_cand_b0tojpsik0star_lite::AbsEtaJpsiDauNeg,
                  // K0* features
                  hf_cand_b0tojpsik0star_lite::PtBach0,
                  hf_cand_b0tojpsik0star_lite::ItsNClsLfTrack0,
                  hf_cand_b0tojpsik0star_lite::TpcNClsCrossedRowsLfTrack0,
                  hf_cand_b0tojpsik0star_lite::ItsChi2NClLfTrack0,
                  hf_cand_b0tojpsik0star_lite::TpcChi2NClLfTrack0,
                  hf_cand_b0tojpsik0star_lite::AbsEtaLfTrack0,
                  hf_cand_b0tojpsik0star_lite::NSigTpcKaBachelor0,
                  hf_cand_b0tojpsik0star_lite::NSigTofKaBachelor0,
                  hf_cand_b0tojpsik0star_lite::NSigTpcTofKaBachelor0,
                  hf_cand_b0tojpsik0star_lite::NSigTpcPiBachelor0,
                  hf_cand_b0tojpsik0star_lite::NSigTofPiBachelor0,
                  hf_cand_b0tojpsik0star_lite::NSigTpcTofPiBachelor0,
                  hf_cand_b0tojpsik0star_lite::PtBach1,
                  hf_cand_b0tojpsik0star_lite::ItsNClsLfTrack1,
                  hf_cand_b0tojpsik0star_lite::TpcNClsCrossedRowsLfTrack1,
                  hf_cand_b0tojpsik0star_lite::ItsChi2NClLfTrack1,
                  hf_cand_b0tojpsik0star_lite::TpcChi2NClLfTrack1,
                  hf_cand_b0tojpsik0star_lite::AbsEtaLfTrack1,
                  hf_cand_b0tojpsik0star_lite::NSigTpcKaBachelor1,
                  hf_cand_b0tojpsik0star_lite::NSigTofKaBachelor1,
                  hf_cand_b0tojpsik0star_lite::NSigTpcTofKaBachelor1,
                  hf_cand_b0tojpsik0star_lite::NSigTpcPiBachelor1,
                  hf_cand_b0tojpsik0star_lite::NSigTofPiBachelor1,
                  hf_cand_b0tojpsik0star_lite::NSigTpcTofPiBachelor1,
                  // MC truth
                  hf_cand_mc_flag::FlagMcMatchRec,
                  hf_cand_mc_flag::FlagMcDecayChanRec,
                  hf_cand_mc_flag::OriginMcRec,
                  hf_cand_b0tojpsik0star_lite::FlagWrongCollision,
                  hf_cand_b0tojpsik0star_lite::PtGen);

} // namespace o2::aod

// string definitions, used for histogram axis labels
const TString stringPt = "#it{p}_{T} (GeV/#it{c})";
const TString stringPtJpsi = "#it{p}_{T}(Jpsi) (GeV/#it{c});";
const TString b0CandTitle = "B^{0} candidates;";
const TString entries = "entries";
const TString b0CandMatch = "B^{0} candidates (matched);";
const TString b0CandUnmatch = "B^{0} candidates (unmatched);";
const TString mcParticleMatched = "MC particles (matched);";

/// B0 analysis task
struct HfTaskB0ToJpsiK0StarReduced {
  Produces<aod::HfRedCandB0Lites> hfRedCandB0Lite;
  // Produces<aod::HfRedB0McCheck> hfRedB0McCheck;

  Configurable<uint8_t> selectionFlagB0{"selectionFlagB0", 1, "Selection Flag for B0"};
  Configurable<double> yCandGenMax{"yCandGenMax", 0.5, "max. gen particle rapidity"};
  Configurable<double> yCandRecoMax{"yCandRecoMax", 0.8, "max. cand. rapidity"};
  Configurable<float> etaTrackMax{"etaTrackMax", 0.8, "max. track pseudo-rapidity"};
  Configurable<float> ptTrackMin{"ptTrackMin", 0.1, "min. track transverse momentum"};
  Configurable<bool> fillBackground{"fillBackground", false, "Flag to enable filling of background histograms/sparses/tree (only MC)"};
  Configurable<float> downSampleBkgFactor{"downSampleBkgFactor", 1., "Fraction of background candidates to keep for ML trainings"};
  Configurable<float> ptMaxForDownSample{"ptMaxForDownSample", 10., "Maximum pt for the application of the downsampling factor"};
  Configurable<bool> useJpsiPdgMass{"useJpsiPdgMass", true, "Whether to use J/Psi PDG mass for B+ candidate mass evaluation or to use the invariant mass of the two prongs"};
  Configurable<bool> useK0StarPdgMass{"useK0StarPdgMass", true, "Whether to use K*0 PDG mass for B+ candidate mass evaluation or to use the invariant mass of the two prongs"};
  // topological cuts
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_b0_to_jpsi_k0star::vecBinsPt}, "pT bin limits"};
  Configurable<LabeledArray<double>> cuts{"cuts", {hf_cuts_b0_to_jpsi_k0star::Cuts[0], hf_cuts_b0_to_jpsi_k0star::NBinsPt, hf_cuts_b0_to_jpsi_k0star::NCutVars, hf_cuts_b0_to_jpsi_k0star::labelsPt, hf_cuts_b0_to_jpsi_k0star::labelsCutVar}, "B0 candidate selection per pT bin"};
  // Enable PID
  Configurable<int> kaonPidMethod{"kaonPidMethod", PidMethod::TpcOrTof, "PID selection method for the bachelor kaon (PidMethod::NoPid: none, PidMethod::TpcOrTof: TPC or TOF, PidMethod::TpcAndTof: TPC and TOF)"};
  Configurable<int> pionPidMethod{"pionPidMethod", PidMethod::TpcOrTof, "PID selection method for the bachelor pion (PidMethod::NoPid: none, PidMethod::TpcOrTof: TPC or TOF, PidMethod::TpcAndTof: TPC and TOF)"};
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
  // B0 ML inference
  Configurable<std::vector<double>> binsPtB0Ml{"binsPtB0Ml", std::vector<double>{hf_cuts_ml::vecBinsPt}, "pT bin limits for ML application"};
  Configurable<std::vector<int>> cutDirB0Ml{"cutDirB0Ml", std::vector<int>{hf_cuts_ml::vecCutDir}, "Whether to reject score values greater or smaller than the threshold"};
  Configurable<LabeledArray<double>> cutsB0Ml{"cutsB0Ml", {hf_cuts_ml::Cuts[0], hf_cuts_ml::NBinsPt, hf_cuts_ml::NCutScores, hf_cuts_ml::labelsPt, hf_cuts_ml::labelsCutScore}, "ML selections per pT bin"};
  Configurable<int> nClassesB0Ml{"nClassesB0Ml", static_cast<int>(hf_cuts_ml::NCutScores), "Number of classes in ML model"};
  Configurable<std::vector<std::string>> namesInputFeatures{"namesInputFeatures", std::vector<std::string>{"feature1", "feature2"}, "Names of ML model input features"};
  // CCDB configuration
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::vector<std::string>> modelPathsCCDB{"modelPathsCCDB", std::vector<std::string>{"path_ccdb/BDT_BS/"}, "Paths of models on CCDB"};
  Configurable<std::vector<std::string>> onnxFileNames{"onnxFileNames", std::vector<std::string>{"ModelHandler_onnx_BSToJpsiPhi.onnx"}, "ONNX file names for each pT bin (if not from CCDB full path)"};
  Configurable<int64_t> timestampCCDB{"timestampCCDB", -1, "timestamp of the ONNX file for ML model used to query in CCDB"};
  Configurable<bool> loadModelsFromCCDB{"loadModelsFromCCDB", false, "Flag to enable or disable the loading of models from CCDB"};

  TrackSelectorKa selectorKaon;
  TrackSelectorPi selectorPion;
  o2::analysis::HfMlResponseB0ToJpsiK0StarReduced<float> hfMlResponse;
  o2::ccdb::CcdbApi ccdbApi;

  using TracksKaon = soa::Join<HfRedTracks, HfRedTracksPid>;
  std::vector<float> outputMl;

  // Filter filterSelectCandidates = (aod::hf_sel_candidate_bplus::isSelBsToJpsiPi >= selectionFlagBs);

  HistogramRegistry registry{"registry"};

  void init(InitContext&)
  {
    std::array<bool, 2> processFuncData{doprocessData, doprocessDataWithB0Ml};
    if ((std::accumulate(processFuncData.begin(), processFuncData.end(), 0)) > 1) {
      LOGP(fatal, "Only one process function for data can be enabled at a time.");
    }
    std::array<bool, 2> processFuncMc{doprocessMc, doprocessMcWithB0Ml};
    if ((std::accumulate(processFuncMc.begin(), processFuncMc.end(), 0)) > 1) {
      LOGP(fatal, "Only one process function for MC can be enabled at a time.");
    }

    if (kaonPidMethod < 0 || kaonPidMethod >= PidMethod::NPidMethods || pionPidMethod < 0 || pionPidMethod >= PidMethod::NPidMethods) {
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

    if (pionPidMethod != PidMethod::NoPid) {
      selectorPion.setRangePtTpc(ptPidTpcMin, ptPidTpcMax);
      selectorPion.setRangeNSigmaTpc(-nSigmaTpcMax, nSigmaTpcMax);
      selectorPion.setRangeNSigmaTpcCondTof(-nSigmaTpcCombinedMax, nSigmaTpcCombinedMax);
      selectorPion.setRangePtTof(ptPidTofMin, ptPidTofMax);
      selectorPion.setRangeNSigmaTof(-nSigmaTofMax, nSigmaTofMax);
      selectorPion.setRangeNSigmaTofCondTpc(-nSigmaTofCombinedMax, nSigmaTofCombinedMax);
    }

    const AxisSpec axisMassB0{150, 4.5, 6.0};
    const AxisSpec axisMassJpsi{600, 2.8f, 3.4f};
    const AxisSpec axisMassK0Star{200, 0.9f, 1.1f};
    const AxisSpec axisPtProng{100, 0., 10.};
    const AxisSpec axisImpactPar{200, -0.05, 0.05};
    const AxisSpec axisPtJpsi{100, 0., 50.};
    const AxisSpec axisRapidity{100, -2., 2.};
    const AxisSpec axisPtB{(std::vector<double>)binsPt, "#it{p}_{T}^{B^{0}} (GeV/#it{c})"};
    const AxisSpec axisPtK0Star{100, 0.f, 10.f};

    registry.add("hMass", b0CandTitle + "inv. mass J/#Psi K*^{0} (GeV/#it{c}^{2});" + stringPt, {HistType::kTH2F, {axisMassB0, axisPtB}});
    registry.add("hMassJpsi", b0CandTitle + "inv. mass #mu^{+}#mu^{#minus} (GeV/#it{c}^{2});" + stringPt, {HistType::kTH2F, {axisMassJpsi, axisPtJpsi}});
    registry.add("hMassK0Star", b0CandTitle + "inv. mass K^{+}#pi^{#minus} (GeV/#it{c}^{2});" + stringPt, {HistType::kTH2F, {axisMassK0Star, axisPtK0Star}});

    // histograms processMC
    if (doprocessMc || doprocessMcWithB0Ml) {
      registry.add("hPtJpsiGen", mcParticleMatched + "J/#Psi #it{p}_{T}^{gen} (GeV/#it{c}); B^{0} " + stringPt, {HistType::kTH2F, {axisPtProng, axisPtB}});
      registry.add("hPtK0StarGen", mcParticleMatched + "K*^{0} #it{p}_{T}^{gen} (GeV/#it{c}); B^{0} " + stringPt, {HistType::kTH2F, {axisPtK0Star, axisPtB}});
      registry.add("hYGenWithProngsInAcceptance", mcParticleMatched + "B^{0} #it{p}_{T}^{gen} (GeV/#it{c}); B^{0} #it{y}", {HistType::kTH2F, {axisPtProng, axisRapidity}});
      registry.add("hMassRecSig", b0CandMatch + "inv. mass J/#Psi K*^{0} (GeV/#it{c}^{2}); B^{0} " + stringPt, {HistType::kTH2F, {axisMassB0, axisPtB}});
      registry.add("hMassJpsiRecSig", b0CandMatch + "inv. mass #mu^{+}#mu^{#minus} (GeV/#it{c}^{2}); J/#Psi " + stringPt, {HistType::kTH2F, {axisMassJpsi, axisPtJpsi}});
      registry.add("hMassK0StarRecSig", b0CandMatch + "inv. mass K^{+}#pi^{#minus} (GeV/#it{c}^{2}); K*^{0} " + stringPt, {HistType::kTH2F, {axisMassK0Star, axisPtK0Star}});
      registry.add("hMassRecBg", b0CandUnmatch + "inv. mass J/#Psi K*^{0} (GeV/#it{c}^{2}); B^{0} " + stringPt, {HistType::kTH2F, {axisMassB0, axisPtB}});
      registry.add("hMassJpsiRecBg", b0CandUnmatch + "inv. mass #mu^{+}#mu^{#minus} (GeV/#it{c}^{2}); J/#Psi " + stringPt, {HistType::kTH2F, {axisMassJpsi, axisPtJpsi}});
      registry.add("hMassK0StarRecBg", b0CandUnmatch + "inv. mass K^{+}#pi^{#minus} (GeV/#it{c}^{2}); K*^{0} " + stringPt, {HistType::kTH2F, {axisMassK0Star, axisPtK0Star}});
    }

    if (doprocessDataWithB0Ml || doprocessMcWithB0Ml) {
      hfMlResponse.configure(binsPtB0Ml, cutsB0Ml, cutDirB0Ml, nClassesB0Ml);
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

  /// Selection of B0 daughter in geometrical acceptance
  /// \param etaProng is the pseudorapidity of B0 prong
  /// \param ptProng is the pT of B0 prong
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
  /// \param withB0Ml is the flag to enable the filling with ML scores for the B0 candidate
  /// \param candidate is the B0 candidate
  /// \param candidatesJpsi is the table with Jpsi candidates
  template <bool DoMc, bool WithB0Ml, typename Cand>
  void fillCand(Cand const& candidate,
                aod::HfRedJpsis const& /*candidatesJpsi*/,
                aod::HfRedBach0Tracks const&,
                aod::HfRedBach1Tracks const&)
  {
    auto ptCandB0 = candidate.pt();
    auto invMassB0KPi = HfHelper::invMassB0ToJpsiK0Star(candidate, useJpsiPdgMass, useK0StarPdgMass, true);
    auto invMassB0PiK = HfHelper::invMassB0ToJpsiK0Star(candidate, useJpsiPdgMass, useK0StarPdgMass, false);
    auto candJpsi = candidate.template jpsi_as<aod::HfRedJpsis>();
    auto candLfDau0 = candidate.template prong0K0Star_as<aod::HfRedBach0Tracks>();
    auto candLfDau1 = candidate.template prong1K0Star_as<aod::HfRedBach1Tracks>();
    auto const pVecMu0 = candidate.pVectorProng0();
    auto const pVecMu1 = candidate.pVectorProng1();
    auto const pVecLfDau0 = candidate.pVectorProng2();
    auto const pVecLfDau1 = candidate.pVectorProng3();
    auto ptJpsi = RecoDecay::pt(pVecMu0, pVecMu1);
    auto ptK0Star = RecoDecay::pt(pVecLfDau0, pVecLfDau1);
    auto invMassJpsi = RecoDecay::m(std::array{pVecMu0, pVecMu1}, std::array{o2::constants::physics::MassMuonPlus, o2::constants::physics::MassMuonMinus});
    auto invMassK0StarKPi = RecoDecay::m(std::array{pVecLfDau0, pVecLfDau1}, std::array{o2::constants::physics::MassKPlus, o2::constants::physics::MassPiPlus});
    auto invMassK0StarPiK = RecoDecay::m(std::array{pVecLfDau0, pVecLfDau1}, std::array{o2::constants::physics::MassPiPlus, o2::constants::physics::MassKPlus});

    int8_t flagMcMatchRec{0}, flagMcDecayChanRec{0}, flagWrongCollision{0};
    bool isSignal = false;
    if constexpr (DoMc) {
      flagMcMatchRec = candidate.flagMcMatchRec();
      flagMcDecayChanRec = candidate.flagMcDecayChanRec();
      flagWrongCollision = candidate.flagWrongCollision();
      isSignal = std::abs(flagMcMatchRec) == o2::hf_decay::hf_cand_beauty::B0ToJpsiPiK &&
                 flagMcDecayChanRec == o2::hf_decay::hf_cand_beauty::B0ToJpsiKstar0;
    }

    bool isSelectedKPi{selectionFlagB0 < BIT(SelectionStep::RecoSkims) * 2 - 1};
    bool isSelectedPiK{selectionFlagB0 < BIT(SelectionStep::RecoSkims) * 2 - 1};

    // topological selection for the two K*0 mass hypotheses
    bool selKPi = HfHelper::selectionB0ToJpsiK0StarTopol(candidate, cuts, binsPt, useJpsiPdgMass, useK0StarPdgMass, true);
    bool selPiK = HfHelper::selectionB0ToJpsiK0StarTopol(candidate, cuts, binsPt, useJpsiPdgMass, useK0StarPdgMass, false);

    if (!selKPi && !selPiK && selectionFlagB0 >= BIT(SelectionStep::RecoTopol) * 2 - 1) {
      return;
    }
    if (selKPi && selectionFlagB0 < BIT(SelectionStep::RecoTopol) * 2 - 1) {
      isSelectedKPi = true;
    }
    if (selPiK && selectionFlagB0 < BIT(SelectionStep::RecoTopol) * 2 - 1) {
      isSelectedPiK = true;
    }

    // track-level PID selection
    // KPi: LfTrack0 as K and LfTrack1 as π ; PiK: LfTrack0 as π and LfTrack1 as K
    if (kaonPidMethod != PidMethod::NoPid || pionPidMethod != PidMethod::NoPid) {
      int pidKa0{TrackSelectorPID::Status::NotApplicable};
      int pidKa1{TrackSelectorPID::Status::NotApplicable};
      int pidPi0{TrackSelectorPID::Status::NotApplicable};
      int pidPi1{TrackSelectorPID::Status::NotApplicable};
      if (kaonPidMethod == PidMethod::TpcOrTof) {
        pidKa0 = selectorKaon.statusTpcOrTof(candLfDau0);
        pidKa1 = selectorKaon.statusTpcOrTof(candLfDau1);
      } else if (kaonPidMethod == PidMethod::TpcAndTof) {
        pidKa0 = selectorKaon.statusTpcAndTof(candLfDau0);
        pidKa1 = selectorKaon.statusTpcAndTof(candLfDau1);
      }
      if (pionPidMethod == PidMethod::TpcOrTof) {
        pidPi0 = selectorPion.statusTpcOrTof(candLfDau0);
        pidPi1 = selectorPion.statusTpcOrTof(candLfDau1);
      } else if (pionPidMethod == PidMethod::TpcAndTof) {
        pidPi0 = selectorPion.statusTpcAndTof(candLfDau0);
        pidPi1 = selectorPion.statusTpcAndTof(candLfDau1);
      }
      bool const pidKPi = HfHelper::selectionB0ToJpsiK0StarPid(pidKa0, acceptPIDNotApplicable.value) &&
                          HfHelper::selectionB0ToJpsiK0StarPid(pidPi1, acceptPIDNotApplicable.value);
      bool const pidPiK = HfHelper::selectionB0ToJpsiK0StarPid(pidPi0, acceptPIDNotApplicable.value) &&
                          HfHelper::selectionB0ToJpsiK0StarPid(pidKa1, acceptPIDNotApplicable.value);
      selKPi = selKPi && pidKPi;
      selPiK = selPiK && pidPiK;
      if (!selKPi && !selPiK && selectionFlagB0 >= BIT(SelectionStep::RecoPID) * 2 - 1) {
        return;
      }
      if (selKPi && selectionFlagB0 < BIT(SelectionStep::RecoPID) * 2 - 1) {
        isSelectedKPi = true;
      }
      if (selPiK && selectionFlagB0 < BIT(SelectionStep::RecoPID) * 2 - 1) {
        isSelectedPiK = true;
      }
    }

    // ML selection, evaluated for each hypothesis that survived topological + PID selection
    float mlScoreSigKPi = -1;
    float mlScoreSigPiK = -1;
    if constexpr (WithB0Ml) {
      if (selKPi) {
        std::vector<float> inputFeaturesKPi = hfMlResponse.getInputFeatures(candidate, candLfDau0, candLfDau1, true);
        outputMl.clear();
        selKPi = hfMlResponse.isSelectedMl(inputFeaturesKPi, ptCandB0, outputMl);
        if (selKPi) {
          mlScoreSigKPi = outputMl[1];
        }
      }
      if (selPiK) {
        std::vector<float> inputFeaturesPiK = hfMlResponse.getInputFeatures(candidate, candLfDau0, candLfDau1, false);
        outputMl.clear();
        selPiK = hfMlResponse.isSelectedMl(inputFeaturesPiK, ptCandB0, outputMl);
        if (selPiK) {
          mlScoreSigPiK = outputMl[1];
        }
      }
      if (!selKPi && !selPiK && selectionFlagB0 >= BIT(SelectionStep::RecoMl) * 2 - 1) {
        return;
      }
    }

    registry.fill(HIST("hMassJpsi"), invMassJpsi, ptJpsi);
    if constexpr (DoMc) {
      if (isSignal) {
        registry.fill(HIST("hMassJpsiRecSig"), invMassJpsi, ptJpsi);
      } else if (fillBackground) {
        registry.fill(HIST("hMassJpsiRecBg"), invMassJpsi, ptJpsi);
      }
    }
    if (selKPi) {
      registry.fill(HIST("hMass"), invMassB0KPi, ptCandB0);
      registry.fill(HIST("hMassK0Star"), invMassK0StarKPi, ptK0Star);
    }
    if (selPiK) {
      registry.fill(HIST("hMass"), invMassB0PiK, ptCandB0);
      registry.fill(HIST("hMassK0Star"), invMassK0StarPiK, ptK0Star);
    }
    if constexpr (DoMc) {
      if (isSignal) {
        if (selKPi) {
          registry.fill(HIST("hMassRecSig"), invMassB0KPi, ptCandB0);
          registry.fill(HIST("hMassK0StarRecSig"), invMassK0StarKPi, ptK0Star);
        }
        if (selPiK) {
          registry.fill(HIST("hMassRecSig"), invMassB0PiK, ptCandB0);
          registry.fill(HIST("hMassK0StarRecSig"), invMassK0StarPiK, ptK0Star);
        }
      } else if (fillBackground) {
        if (selKPi) {
          registry.fill(HIST("hMassRecBg"), invMassB0KPi, ptCandB0);
          registry.fill(HIST("hMassK0StarRecBg"), invMassK0StarKPi, ptK0Star);
        }
        if (selPiK) {
          registry.fill(HIST("hMassRecBg"), invMassB0PiK, ptCandB0);
          registry.fill(HIST("hMassK0StarRecBg"), invMassK0StarPiK, ptK0Star);
        }
      }
    }

    // downsampling for ML trainings
    float const pseudoRndm = ptJpsi * 1000. - static_cast<int64_t>(ptJpsi * 1000);
    if (ptCandB0 < ptMaxForDownSample && pseudoRndm >= downSampleBkgFactor) {
      return;
    }

    float ptMother = -1.;
    if constexpr (DoMc) {
      ptMother = candidate.ptMother();
    }

    auto fillTable = [&](bool isSelKPi) {
      auto ctXY = isSelKPi ? candidate.ctXY(std::array{o2::constants::physics::MassMuon, o2::constants::physics::MassMuon, o2::constants::physics::MassKPlus, o2::constants::physics::MassPiPlus})
                           : candidate.ctXY(std::array{o2::constants::physics::MassMuon, o2::constants::physics::MassMuon, o2::constants::physics::MassPiPlus, o2::constants::physics::MassKPlus});
      auto mlScoreSig = isSelKPi ? mlScoreSigKPi : mlScoreSigPiK;
      auto invMassB0 = isSelKPi ? invMassB0KPi : invMassB0PiK;
      auto invMassK0Star = isSelKPi ? invMassK0StarKPi : invMassK0StarPiK;
      hfRedCandB0Lite(
        // B0 - meson features
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
        ctXY,
        candidate.impactParameterProduct(),
        candidate.impactParameterProductJpsi(),
        candidate.impactParameterProductK0Star(),
        candidate.maxNormalisedDeltaIP(),
        mlScoreSig,
        // J/Psi features
        invMassJpsi,
        ptJpsi,
        invMassK0Star,
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
        // K*0 daughter features
        candLfDau0.pt(),
        candLfDau0.itsNCls(),
        candLfDau0.tpcNClsCrossedRows(),
        candLfDau0.itsChi2NCl(),
        candLfDau0.tpcChi2NCl(),
        absEta(candLfDau0.tgl()),
        candLfDau0.tpcNSigmaKa(),
        candLfDau0.tofNSigmaKa(),
        candLfDau0.tpcTofNSigmaKa(),
        candLfDau0.tpcNSigmaPi(),
        candLfDau0.tofNSigmaPi(),
        candLfDau0.tpcTofNSigmaPi(),
        candLfDau1.pt(),
        candLfDau1.itsNCls(),
        candLfDau1.tpcNClsCrossedRows(),
        candLfDau1.itsChi2NCl(),
        candLfDau1.tpcChi2NCl(),
        absEta(candLfDau1.tgl()),
        candLfDau1.tpcNSigmaKa(),
        candLfDau1.tofNSigmaKa(),
        candLfDau1.tpcTofNSigmaKa(),
        candLfDau1.tpcNSigmaPi(),
        candLfDau1.tofNSigmaPi(),
        candLfDau1.tpcTofNSigmaPi(),
        // MC truth
        flagMcMatchRec,
        flagMcDecayChanRec,
        isSignal,
        flagWrongCollision,
        ptMother);
    };

    if (selKPi) {
      fillTable(true);
    }
    if (selPiK) {
      fillTable(false);
    }
  }

  /// Fill particle histograms (gen MC truth)
  void fillCandMcGen(aod::HfMcGenRedB0s::iterator const& particle)
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
    registry.fill(HIST("hPtK0StarGen"), ptProngs[1], ptParticle);

    // generated B0 with daughters in geometrical acceptance
    if (prongsInAcc) {
      registry.fill(HIST("hYGenWithProngsInAcceptance"), ptParticle, yParticle);
    }
  }

  // Process functions
  void processData(aod::HfRedCandB0ToJpsiK0Star const& candidates,
                   aod::HfRedJpsis const& candidatesJpsi,
                   aod::HfRedBach0Tracks const& kaon0Tracks,
                   aod::HfRedBach1Tracks const& kaon1Tracks)
  {
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(HfHelper::yB0(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCand<false, false>(candidate, candidatesJpsi, kaon0Tracks, kaon1Tracks);
    } // candidate loop
  } // processData
  PROCESS_SWITCH(HfTaskB0ToJpsiK0StarReduced, processData, "Process data without ML for B0", true);

  void processDataWithB0Ml(aod::HfRedCandB0ToJpsiK0Star const& candidates,
                           aod::HfRedJpsis const& candidatesJpsi,
                           aod::HfRedBach0Tracks const& kaon0Tracks,
                           aod::HfRedBach1Tracks const& kaon1Tracks)
  {
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(HfHelper::yB0(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCand<false, true>(candidate, candidatesJpsi, kaon0Tracks, kaon1Tracks);
    } // candidate loop
  } // processDataWithB0Ml
  PROCESS_SWITCH(HfTaskB0ToJpsiK0StarReduced, processDataWithB0Ml, "Process data with ML for B0", false);

  void processMc(soa::Join<aod::HfRedCandB0ToJpsiK0Star, aod::HfMcRecRedB0s> const& candidates,
                 aod::HfMcGenRedB0s const& mcParticles,
                 aod::HfRedJpsis const& candidatesJpsi,
                 aod::HfRedBach0Tracks const& kaon0Tracks,
                 aod::HfRedBach1Tracks const& kaon1Tracks)
  {
    // MC rec
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(HfHelper::yB0(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCand<true, false>(candidate, candidatesJpsi, kaon0Tracks, kaon1Tracks);
    } // rec

    // MC gen. level
    for (const auto& particle : mcParticles) {
      fillCandMcGen(particle);
    } // gen
  } // processMc
  PROCESS_SWITCH(HfTaskB0ToJpsiK0StarReduced, processMc, "Process MC without ML for B0", false);

  void processMcWithB0Ml(soa::Join<aod::HfRedCandB0ToJpsiK0Star, aod::HfMcRecRedB0s> const& candidates,
                         aod::HfMcGenRedB0s const& mcParticles,
                         aod::HfRedJpsis const& candidatesJpsi,
                         aod::HfRedBach0Tracks const& kaon0Tracks,
                         aod::HfRedBach1Tracks const& kaon1Tracks)
  {
    // MC rec
    for (const auto& candidate : candidates) {
      if (yCandRecoMax >= 0. && std::abs(HfHelper::yB0(candidate)) > yCandRecoMax) {
        continue;
      }
      fillCand<true, true>(candidate, candidatesJpsi, kaon0Tracks, kaon1Tracks);
    } // rec

    // MC gen. level
    for (const auto& particle : mcParticles) {
      fillCandMcGen(particle);
    } // gen
  } // processMcWithB0Ml
  PROCESS_SWITCH(HfTaskB0ToJpsiK0StarReduced, processMcWithB0Ml, "Process MC with ML for B0", false);

}; // struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskB0ToJpsiK0StarReduced>(cfgc)};
}
