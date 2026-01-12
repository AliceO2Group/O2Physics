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
//
/// \file qaMatching.cxx
/// \brief Task to compute and evaluate DCA quantities
/// \author Nicolas Bizé <nicolas.bize@cern.ch>, SUBATECH
//
#include "PWGDQ/Core/MuonMatchingMlResponse.h"
#include "PWGDQ/Core/VarManager.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"

#include "Common/DataModel/EventSelection.h"

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "GlobalTracking/MatchGlobalFwd.h"
#include "MFTTracking/Constants.h"

#include <Math/ProbFunc.h>

#include <algorithm>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::aod;

using MyEvents = soa::Join<aod::Collisions, aod::EvSels>;
using MyMuons = soa::Join<aod::FwdTracks, aod::FwdTracksCov>;
using MyMuonsMC = soa::Join<aod::FwdTracks, aod::FwdTracksCov, aod::McFwdTrackLabels>;
using MyMFTs = aod::MFTTracks;
using MyMFTCovariances = aod::MFTTracksCov;
using MyMFTsMC = soa::Join<aod::MFTTracks, aod::McMFTTrackLabels>;

using MyMuon = MyMuons::iterator;
using MyMuonMC = MyMuonsMC::iterator;
using MyMFT = MyMFTs::iterator;
using MyMFTCovariance = MyMFTCovariances::iterator;

using SMatrix55 = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
using SMatrix5 = ROOT::Math::SVector<Double_t, 5>;

static float chi2ToScore(float chi2, int ndf, float chi2max)
{
  double p = -TMath::Log10(ROOT::Math::chisquared_cdf_c(chi2, ndf));
  double pnorm = -TMath::Log10(ROOT::Math::chisquared_cdf_c(chi2max, ndf));
  double result = (1.f / (p / pnorm + 1.f));
  return static_cast<float>(result);
}

struct qaMatching {

  template <class T, int nr, int nc>
  using matrix = std::array<std::array<T, nc>, nr>;

  enum MuonMatchType {
    kMatchTypeTrueLeading = 0,
    kMatchTypeWrongLeading = 1,
    kMatchTypeDecayLeading = 2,
    kMatchTypeFakeLeading = 3,
    kMatchTypeTrueNonLeading = 4,
    kMatchTypeWrongNonLeading = 5,
    kMatchTypeDecayNonLeading = 6,
    kMatchTypeFakeNonLeading = 7,
    kMatchTypeUndefined
  };

  struct MatchingCandidate {
    int64_t collisionId{-1};
    int64_t globalTrackId{-1};
    int64_t muonTrackId{-1};
    int64_t mftTrackId{-1};
    double matchScore{-1};
    double matchChi2{-1};
    int matchRanking{-1};
    double matchScoreProd{-1};
    double matchChi2Prod{-1};
    int matchRankingProd{-1};
    MuonMatchType matchType{kMatchTypeUndefined};
  };

  ////   Variables for selecting muon tracks
  Configurable<float> fPMchLow{"cfgPMchLow", 0.0f, ""};
  Configurable<float> fPtMchLow{"cfgPtMchLow", 0.7f, ""};
  Configurable<float> fEtaMchLow{"cfgEtaMchLow", -4.0f, ""};
  Configurable<float> fEtaMchUp{"cfgEtaMchUp", -2.5f, ""};
  Configurable<float> fRabsLow{"cfgRabsLow", 17.6f, ""};
  Configurable<float> fRabsUp{"cfgRabsUp", 89.5f, ""};
  Configurable<float> fSigmaPdcaUp{"cfgPdcaUp", 6.f, ""};
  Configurable<float> fTrackChi2MchUp{"cfgTrackChi2MchUp", 5.f, ""};
  Configurable<float> fMatchingChi2MchMidUp{"cfgMatchingChi2MchMidUp", 999.f, ""};

  ////   Variables for selecting mft tracks
  Configurable<float> fEtaMftLow{"cfgEtaMftlow", -3.6f, ""};
  Configurable<float> fEtaMftUp{"cfgEtaMftup", -2.5f, ""};
  Configurable<int> fTrackNClustMftLow{"cfgTrackNClustMftLow", 7, ""};
  Configurable<float> fTrackChi2MftUp{"cfgTrackChi2MftUp", 999.f, ""};

  ////   Variables for selecting global tracks
  Configurable<float> fMatchingChi2ScoreMftMchLow{"cfgMatchingChi2ScoreMftMchLow", chi2ToScore(50.f, 5, 50.f), ""};

  ////   Variables for selecting tagged muons
  Configurable<int> fMuonTaggingNCrossedMftPlanesLow{"cfgMuonTaggingNCrossedMftPlanesLow", 5, ""};
  Configurable<float> fMuonTaggingTrackChi2MchUp{"cfgMuonTaggingTrackChi2MchUp", 5.f, ""};
  Configurable<float> fMuonTaggingPMchLow{"cfgMuonTaggingPMchLow", 0.0f, ""};
  Configurable<float> fMuonTaggingPtMchLow{"cfgMuonTaggingPtMchLow", 0.7f, ""};
  Configurable<float> fMuonTaggingEtaMchLow{"cfgMuonTaggingEtaMchLow", -3.6f, ""};
  Configurable<float> fMuonTaggingEtaMchUp{"cfgMuonTaggingEtaMchUp", -2.5f, ""};
  Configurable<float> fMuonTaggingRabsLow{"cfgMuonTaggingRabsLow", 17.6f, ""};
  Configurable<float> fMuonTaggingRabsUp{"cfgMuonTaggingRabsUp", 89.5f, ""};
  Configurable<float> fMuonTaggingSigmaPdcaUp{"cfgMuonTaggingPdcaUp", 4.f, ""};
  Configurable<float> fMuonTaggingChi2DiffLow{"cfgMuonTaggingChi2DiffLow", 100.f, ""};

  ///    Variables to event mixing criteria
  Configurable<float> fSaveMixedMatchingParamsRate{"cfgSaveMixedMatchingParamsRate", 0.002f, ""};
  Configurable<int> fEventMaxDeltaNMFT{"cfgEventMaxDeltaNMFT", 1, ""};
  Configurable<float> fEventMaxDeltaVtxZ{"cfgEventMaxDeltaVtxZ", 1.f, ""};
  Configurable<int> fEventMinDeltaBc{"cfgEventMinDeltaBc", 500, ""};

  ////   Variables for ccdb
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};

  // CCDB connection configurables
  struct : ConfigurableGroup {
    Configurable<std::string> fConfigCcdbUrl{"ccdb-url-", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
    Configurable<int64_t> fConfigNoLaterThan{"ccdb-no-later-than-", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
    Configurable<std::string> fConfigGrpPath{"grpPath-", "GLO/GRP/GRP", "Path of the grp file"};
    Configurable<std::string> fConfigGeoPath{"geoPath-", "GLO/Config/GeometryAligned", "Path of the geometry file"};
    Configurable<std::string> fConfigGrpMagPath{"grpmagPath-", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  } fConfigCCDB;

  ///    Variables for histograms configuration
  Configurable<int> fNCandidatesMax{"nCandidatesMax", 5, ""};

  double mBzAtMftCenter{0};

  o2::globaltracking::MatchGlobalFwd mExtrap;

  using MatchingFunc_t = std::function<std::tuple<double, int>(const o2::dataformats::GlobalFwdTrack& mchtrack, const o2::track::TrackParCovFwd& mfttrack)>;
  std::map<std::string, MatchingFunc_t> mMatchingFunctionMap; ///< MFT-MCH Matching function

  // Chi2 matching interface
  static constexpr int sChi2FunctionsNum = 3;
  struct : ConfigurableGroup {
    std::array<Configurable<std::string>, sChi2FunctionsNum> fFunctionLabel{{
      {"cfgChi2FunctionLabel_0", std::string{"ProdAll"}, "Text label identifying this chi2 matching method"},
      {"cfgChi2FunctionLabel_1", std::string{"MatchXYPhiTanlMom"}, "Text label identifying this chi2 matching method"},
      {"cfgChi2FunctionLabel_2", std::string{"MatchXYPhiTanl"}, "Text label identifying this chi2 matching method"},
    }};
    std::array<Configurable<std::string>, sChi2FunctionsNum> fFunctionName{{{"cfgChi2FunctionNames_0", std::string{"prod"}, "Name of the chi2 matching function"},
                                                                            {"cfgChi2FunctionNames_1", std::string{"matchALL"}, "Name of the chi2 matching function"},
                                                                            {"cfgChi2FunctionNames_2", std::string{"matchXYPhiTanl"}, "Name of the chi2 matching function"}}};
    std::array<Configurable<float>, sChi2FunctionsNum> fMatchingScoreCut{{
      {"cfgChi2FunctionMatchingScoreCut_0", 0.f, "Minimum score value for selecting good matches"},
      {"cfgChi2FunctionMatchingScoreCut_1", 0.5f, "Minimum score value for selecting good matches"},
      {"cfgChi2FunctionMatchingScoreCut_2", 0.5f, "Minimum score value for selecting good matches"},
    }};
    std::array<Configurable<float>, sChi2FunctionsNum> fMatchingPlaneZ{{
      {"cfgChi2FunctionMatchingPlaneZ_0", static_cast<float>(o2::mft::constants::mft::LayerZCoordinate()[9]), "Z position of the matching plane"},
      {"cfgChi2FunctionMatchingPlaneZ_1", static_cast<float>(o2::mft::constants::mft::LayerZCoordinate()[9]), "Z position of the matching plane"},
      {"cfgChi2FunctionMatchingPlaneZ_2", static_cast<float>(o2::mft::constants::mft::LayerZCoordinate()[9]), "Z position of the matching plane"},
    }};
    std::array<Configurable<int>, sChi2FunctionsNum> fMatchingExtrapMethod{{
      {"cfgMatchingExtrapMethod_0", static_cast<int>(0), "Method for MCH track extrapolation to maching plane"},
      {"cfgMatchingExtrapMethod_1", static_cast<int>(0), "Method for MCH track extrapolation to maching plane"},
      {"cfgMatchingExtrapMethod_2", static_cast<int>(0), "Method for MCH track extrapolation to maching plane"},
    }};
  } fConfigChi2MatchingOptions;

  // ML interface
  static constexpr int sMLModelsNum = 2;
  struct : ConfigurableGroup {
    std::array<Configurable<std::string>, sMLModelsNum> fModelLabel{{
      {"cfgMLModelLabel_0", std::string{""}, "Text label identifying this group of ML models"},
      {"cfgMLModelLabel_1", std::string{""}, "Text label identifying this group of ML models"},
    }};
    std::array<Configurable<std::vector<std::string>>, sMLModelsNum> fModelPathsCCDB{{{"cfgMLModelPathsCCDB_0", std::vector<std::string>{"Users/m/mcoquet/MLTest"}, "Paths of models on CCDB"},
                                                                                      {"cfgMLModelPathsCCDB_1", std::vector<std::string>{}, "Paths of models on CCDB"}}};
    std::array<Configurable<std::vector<std::string>>, sMLModelsNum> fInputFeatures{{{"cfgMLInputFeatures_0", std::vector<std::string>{"chi2MCHMFT"}, "Names of ML model input features"},
                                                                                     {"cfgMLInputFeatures_1", std::vector<std::string>{}, "Names of ML model input features"}}};
    std::array<Configurable<std::vector<std::string>>, sMLModelsNum> fModelNames{{{"cfgMLModelNames_0", std::vector<std::string>{"model.onnx"}, "ONNX file names for each pT bin (if not from CCDB full path)"},
                                                                                  {"cfgMLModelNames_1", std::vector<std::string>{}, "ONNX file names for each pT bin (if not from CCDB full path)"}}};
    std::array<Configurable<float>, sMLModelsNum> fMatchingScoreCut{{
      {"cfgMLModelMatchingScoreCut_0", 0.f, "Minimum score value for selecting good matches"},
      {"cfgMLModelMatchingScoreCut_1", 0.f, "Minimum score value for selecting good matches"},
    }};
    std::array<Configurable<float>, sMLModelsNum> fMatchingPlaneZ{{
      {"cfgMLModelMatchingPlaneZ_0", static_cast<float>(o2::mft::constants::mft::LayerZCoordinate()[9]), "Z position of the matching plane"},
      {"cfgMLModelMatchingPlaneZ_1", 0.f, "Z position of the matching plane"},
    }};
    std::array<Configurable<int>, sMLModelsNum> fMatchingExtrapMethod{{
      {"cfgMatchingExtrapMethod_0", static_cast<int>(0), "Method for MCH track extrapolation to maching plane"},
      {"cfgMatchingExtrapMethod_1", static_cast<int>(0), "Method for MCH track extrapolation to maching plane"},
    }};
  } fConfigMlOptions;

  std::vector<double> binsPtMl;
  std::array<double, 1> cutValues;
  std::vector<int> cutDirMl;
  std::map<std::string, o2::analysis::MlResponseMFTMuonMatch<float>> matchingMlResponses;
  std::map<std::string, std::string> matchingChi2Functions;
  std::map<std::string, double> matchingPlanesZ;
  std::map<std::string, double> matchingScoreCuts;
  std::map<std::string, int> matchingExtrapMethod;

  int mRunNumber{0}; // needed to detect if the run changed and trigger update of magnetic field

  Service<o2::ccdb::BasicCCDBManager> ccdbManager;
  o2::ccdb::CcdbApi fCCDBApi;

  o2::aod::rctsel::RCTFlagsChecker rctChecker{"CBT_muon_glo", false, false, true};

  // vector of all MFT-MCH(-MID) matching candidates associated to the same MCH(-MID) track,
  // to be sorted in descending order with respect to the matching score
  // the map key is the MCH(-MID) track global index
  using MatchingCandidates = std::map<int64_t, std::vector<MatchingCandidate>>;

  struct CollisionInfo {
    int64_t index{0};
    uint64_t bc{0};
    // z position of the collision
    double zVertex{0};
    // number of MFT tracks associated to the collision
    int mftTracksMultiplicity{0};
    // vector of MFT track indexes
    std::vector<int64_t> mftTracks;
    // vector of MCH(-MID) track indexes
    std::vector<int64_t> mchTracks;
    // matching candidates
    MatchingCandidates matchingCandidates;
    // vector of MFT-MCH track index pairs belonging to the same MC muon particle
    std::vector<std::pair<int64_t, int64_t>> matchablePairs;
    // vector of MCH track indexes that are expected to have an associated MFT track
    std::vector<int64_t> taggedMuons;
  };

  using CollisionInfos = std::map<int64_t, CollisionInfo>;

  std::unordered_map<int64_t, int32_t> mftTrackCovs;

  std::vector<std::pair<int64_t, int64_t>> fMatchablePairs;
  MatchingCandidates fMatchingCandidates;
  std::vector<int64_t> fTaggedMuons;

  using MuonPair = std::pair<std::pair<int64_t, uint64_t>, std::pair<int64_t, uint64_t>>;
  using GlobalMuonPair = std::pair<std::pair<int64_t, std::vector<MatchingCandidate>>, std::pair<int64_t, std::vector<MatchingCandidate>>>;

  HistogramRegistry registry{"registry", {}};
  HistogramRegistry registryMatching{"registryMatching", {}};
  HistogramRegistry registryMatching0{"registryMatching_0", {}};
  HistogramRegistry registryMatching1{"registryMatching_1", {}};
  HistogramRegistry registryMatching2{"registryMatching_2", {}};
  HistogramRegistry registryMatching3{"registryMatching_3", {}};
  HistogramRegistry registryMatching4{"registryMatching_4", {}};
  HistogramRegistry registryMatching5{"registryMatching_5", {}};
  HistogramRegistry registryMatching6{"registryMatching_6", {}};
  HistogramRegistry registryMatching7{"registryMatching_7", {}};
  HistogramRegistry registryMatching8{"registryMatching_8", {}};
  HistogramRegistry registryMatching9{"registryMatching_9", {}};
  std::vector<HistogramRegistry*> registryMatchingVec{{&registryMatching0,
                                                       &registryMatching1,
                                                       &registryMatching2,
                                                       &registryMatching3,
                                                       &registryMatching4,
                                                       &registryMatching5,
                                                       &registryMatching6,
                                                       &registryMatching7,
                                                       &registryMatching8,
                                                       &registryMatching9}};
  HistogramRegistry registryDimuon{"registryDimuon", {}};

  std::unordered_map<std::string, o2::framework::HistPtr> matchingHistos;
  matrix<o2::framework::HistPtr, 4, 4> dimuonHistos;

  struct EfficiencyPlotter {
    o2::framework::HistPtr p_num;
    o2::framework::HistPtr p_den;
    o2::framework::HistPtr pt_num;
    o2::framework::HistPtr pt_den;
    o2::framework::HistPtr phi_num;
    o2::framework::HistPtr phi_den;
    o2::framework::HistPtr eta_num;
    o2::framework::HistPtr eta_den;

    EfficiencyPlotter(std::string path, std::string title,
                      HistogramRegistry& registry)
    {
      AxisSpec pAxis = {100, 0, 100, "p (GeV/c)"};
      AxisSpec pTAxis = {100, 0, 10, "p_{T} (GeV/c)"};
      AxisSpec etaAxis = {100, -4, -2, "#eta"};
      AxisSpec phiAxis = {90, -180, 180, "#phi (degrees)"};

      std::string histName;
      std::string histTitle;

      // momentum dependence
      histName = path + "p_num";
      histTitle = title + " vs. p - num";
      p_num = registry.add(histName.c_str(), histTitle.c_str(), {HistType::kTH1F, {pAxis}});

      histName = path + "p_den";
      histTitle = title + " vs. p - den";
      p_den = registry.add(histName.c_str(), histTitle.c_str(), {HistType::kTH1F, {pAxis}});

      // pT dependence
      histName = path + "pt_num";
      histTitle = title + " vs. p_{T} - num";
      pt_num = registry.add(histName.c_str(), histTitle.c_str(), {HistType::kTH1F, {pTAxis}});

      histName = path + "pt_den";
      histTitle = title + " vs. p_{T} - den";
      pt_den = registry.add(histName.c_str(), histTitle.c_str(), {HistType::kTH1F, {pTAxis}});

      // eta dependence
      histName = path + "eta_num";
      histTitle = title + " vs. #eta - num";
      eta_num = registry.add(histName.c_str(), histTitle.c_str(), {HistType::kTH1F, {etaAxis}});

      histName = path + "eta_den";
      histTitle = title + " vs. #eta - den";
      eta_den = registry.add(histName.c_str(), histTitle.c_str(), {HistType::kTH1F, {etaAxis}});

      // phi dependence
      histName = path + "phi_num";
      histTitle = title + " vs. #phi - num";
      phi_num = registry.add(histName.c_str(), histTitle.c_str(), {HistType::kTH1F, {phiAxis}});

      histName = path + "phi_den";
      histTitle = title + " vs. #phi - den";
      phi_den = registry.add(histName.c_str(), histTitle.c_str(), {HistType::kTH1F, {phiAxis}});
    }

    template <class T>
    void Fill(const T& track, bool passed)
    {
      double phi = track.phi() * 180 / TMath::Pi();
      std::get<std::shared_ptr<TH1>>(p_den)->Fill(track.p());
      std::get<std::shared_ptr<TH1>>(pt_den)->Fill(track.pt());
      std::get<std::shared_ptr<TH1>>(eta_den)->Fill(track.eta());
      std::get<std::shared_ptr<TH1>>(phi_den)->Fill(phi);

      if (passed) {
        std::get<std::shared_ptr<TH1>>(p_num)->Fill(track.p());
        std::get<std::shared_ptr<TH1>>(pt_num)->Fill(track.pt());
        std::get<std::shared_ptr<TH1>>(eta_num)->Fill(track.eta());
        std::get<std::shared_ptr<TH1>>(phi_num)->Fill(phi);
      }
    }
  };

  struct MatchRankingHistos {
    o2::framework::HistPtr hist;
    o2::framework::HistPtr histVsP;
    o2::framework::HistPtr histVsPt;
    o2::framework::HistPtr histVsMcParticleDz;
    o2::framework::HistPtr histVsMftTrackMult;
    o2::framework::HistPtr histVsMftTrackType;
    o2::framework::HistPtr histVsDeltaChi2;
    o2::framework::HistPtr histVsProdRanking;

    MatchRankingHistos(std::string histName, std::string histTitle, HistogramRegistry* registry)
    {
      AxisSpec pAxis = {100, 0, 100, "p (GeV/c)"};
      AxisSpec ptAxis = {100, 0, 10, "p_{T} (GeV/c)"};
      AxisSpec dzAxis = {100, 0, 50, "#Deltaz (cm)"};
      AxisSpec trackMultAxis = {100, 0, 1000, "MFT track mult."};
      AxisSpec trackTypeAxis = {2, 0, 2, "MFT track type"};
      int matchTypeMax = static_cast<int>(kMatchTypeUndefined);
      AxisSpec matchTypeAxis = {matchTypeMax, 0, static_cast<double>(matchTypeMax), "match type"};
      AxisSpec dchi2Axis = {100, 0, 100, "#Delta#chi^{2}"};
      AxisSpec dqAxis = {3, -1.5, 1.5, "MFT #DeltaQ"};
      AxisSpec indexAxis = {6, 0, 6, "ranking index"};
      AxisSpec indexProdAxis = {6, 0, 6, "ranking index (production)"};

      hist = registry->add(histName.c_str(), histTitle.c_str(), {HistType::kTH1F, {indexAxis}});
      histVsP = registry->add((histName + "VsP").c_str(), (histTitle + " vs. p").c_str(), {HistType::kTH2F, {pAxis, indexAxis}});
      histVsPt = registry->add((histName + "VsPt").c_str(), (histTitle + " vs. p_{T}").c_str(), {HistType::kTH2F, {ptAxis, indexAxis}});
      histVsMcParticleDz = registry->add((histName + "VsMcParticleDz").c_str(), (histTitle + " vs. MC particle #Deltaz").c_str(), {HistType::kTH2F, {dzAxis, indexAxis}});
      histVsMftTrackMult = registry->add((histName + "VsMftTrackMult").c_str(), (histTitle + " vs. MFT track multiplicity").c_str(), {HistType::kTH2F, {trackMultAxis, indexAxis}});
      histVsMftTrackType = registry->add((histName + "VsMftTrackType").c_str(), (histTitle + " vs. MFT track type").c_str(), {HistType::kTH2F, {trackTypeAxis, indexAxis}});
      std::get<std::shared_ptr<TH2>>(histVsMftTrackType)->GetXaxis()->SetBinLabel(1, "Kalman");
      std::get<std::shared_ptr<TH2>>(histVsMftTrackType)->GetXaxis()->SetBinLabel(2, "CA");
      histVsDeltaChi2 = registry->add((histName + "VsDeltaChi2").c_str(), (histTitle + " vs. #Delta#chi^{2}").c_str(), {HistType::kTH2F, {dchi2Axis, indexAxis}});
      histVsProdRanking = registry->add((histName + "VsProdRanking").c_str(), (histTitle + " vs. prod ranking").c_str(), {HistType::kTH2F, {indexProdAxis, indexAxis}});
    }
  };

  struct MatchingPlotter {
    std::unique_ptr<MatchRankingHistos> fMatchRanking;
    std::unique_ptr<MatchRankingHistos> fMatchRankingGoodMCH;
    std::unique_ptr<MatchRankingHistos> fMatchRankingPaired;
    std::unique_ptr<MatchRankingHistos> fMatchRankingPairedGoodMCH;

    //-
    o2::framework::HistPtr fMissedMatches;
    o2::framework::HistPtr fMissedMatchesGoodMCH;
    //-
    o2::framework::HistPtr fMatchRankingWrtProd;
    o2::framework::HistPtr fMatchRankingWrtProdVsP;
    o2::framework::HistPtr fMatchRankingWrtProdVsPt;

    //-
    o2::framework::HistPtr fDecayRankingGoodMatches;
    o2::framework::HistPtr fDecayRankingNonLeadingMatches;
    o2::framework::HistPtr fDecayRankingMissedMatches;

    //-
    o2::framework::HistPtr fScoreGapLeadingTrueMatches;
    o2::framework::HistPtr fScoreGapNonLeadingTrueMatches;

    //-
    o2::framework::HistPtr fMatchType;
    o2::framework::HistPtr fMatchTypeVsP;
    o2::framework::HistPtr fMatchTypeVsPt;

    //-
    o2::framework::HistPtr fMatchScoreVsType;
    o2::framework::HistPtr fMatchScoreVsTypeVsP;
    o2::framework::HistPtr fMatchScoreVsTypeVsPt;
    //-
    o2::framework::HistPtr fMatchChi2VsType;
    o2::framework::HistPtr fMatchChi2VsTypeVsP;
    o2::framework::HistPtr fMatchChi2VsTypeVsPt;
    //-
    o2::framework::HistPtr fMatchScoreVsProd;
    o2::framework::HistPtr fMatchChi2VsProd;
    o2::framework::HistPtr fTrueMatchScoreVsProd;
    o2::framework::HistPtr fTrueMatchChi2VsProd;

    //-
    EfficiencyPlotter fMatchingPurityPlotter;
    EfficiencyPlotter fPairingEfficiencyPlotter;
    EfficiencyPlotter fMatchingEfficiencyPlotter;
    EfficiencyPlotter fFakeMatchingEfficiencyPlotter;

    HistogramRegistry* registry;

    MatchingPlotter(std::string path,
                    HistogramRegistry* reg)
      : fMatchingPurityPlotter(path + "matching-purity/", "Matching purity", *reg),
        fPairingEfficiencyPlotter(path + "pairing-efficiency/", "Pairing efficiency", *reg),
        fMatchingEfficiencyPlotter(path + "matching-efficiency/", "Matching efficiency", *reg),
        fFakeMatchingEfficiencyPlotter(path + "fake-matching-efficiency/", "Fake matching efficiency", *reg)
    {
      registry = reg;
      AxisSpec pAxis = {100, 0, 100, "p (GeV/c)"};
      AxisSpec ptAxis = {100, 0, 10, "p_{T} (GeV/c)"};
      AxisSpec dzAxis = {100, 0, 50, "#Deltaz (cm)"};
      AxisSpec indexAxis = {6, 0, 6, "ranking index"};

      std::string histName = path + "matchRanking";
      std::string histTitle = "True match ranking";

      fMatchRanking = std::make_unique<MatchRankingHistos>(path + "matchRanking", "True match ranking", registry);
      fMatchRankingGoodMCH = std::make_unique<MatchRankingHistos>(path + "matchRankingGoodMCH", "True match ranking (good MCH tracks)", registry);
      fMatchRankingPaired = std::make_unique<MatchRankingHistos>(path + "matchRankingPaired", "True match ranking (paired MCH tracks)", registry);
      fMatchRankingPairedGoodMCH = std::make_unique<MatchRankingHistos>(path + "matchRankingPairedGoodMCH", "True match ranking (good paired MCH tracks)", registry);

      //-
      AxisSpec missedMatchAxis = {5, 0, 5, ""};
      histName = path + "missedMatches";
      histTitle = "Missed matches";
      fMissedMatches = registry->add(histName.c_str(), histTitle.c_str(), {HistType::kTH1F, {missedMatchAxis}});
      std::get<std::shared_ptr<TH1>>(fMissedMatches)->GetXaxis()->SetBinLabel(1, "not paired");
      std::get<std::shared_ptr<TH1>>(fMissedMatches)->GetXaxis()->SetBinLabel(2, "fake MCH");
      std::get<std::shared_ptr<TH1>>(fMissedMatches)->GetXaxis()->SetBinLabel(3, "not stored");
      histName = path + "missedMatchesGoodMCH";
      histTitle = "Missed matches - good MCH tracks";
      fMissedMatchesGoodMCH = registry->add(histName.c_str(), histTitle.c_str(), {HistType::kTH1F, {missedMatchAxis}});
      std::get<std::shared_ptr<TH1>>(fMissedMatchesGoodMCH)->GetXaxis()->SetBinLabel(1, "not paired");
      std::get<std::shared_ptr<TH1>>(fMissedMatchesGoodMCH)->GetXaxis()->SetBinLabel(2, "fake MCH");
      std::get<std::shared_ptr<TH1>>(fMissedMatchesGoodMCH)->GetXaxis()->SetBinLabel(3, "not stored");

      AxisSpec decayRankingAxis = {5, 0, 5, "decay ranking"};
      histName = path + "decayRankingGoodMatches";
      histTitle = "Decay ranking - good matches";
      fDecayRankingGoodMatches = registry->add(histName.c_str(), histTitle.c_str(), {HistType::kTH1F, {decayRankingAxis}});
      histName = path + "decayRankingNonLeadingMatches";
      histTitle = "Decay ranking - non-leading matches";
      fDecayRankingNonLeadingMatches = registry->add(histName.c_str(), histTitle.c_str(), {HistType::kTH1F, {decayRankingAxis}});
      histName = path + "decayRankingMissedMatches";
      histTitle = "Decay ranking - missed matches";
      fDecayRankingMissedMatches = registry->add(histName.c_str(), histTitle.c_str(), {HistType::kTH1F, {decayRankingAxis}});

      AxisSpec scoreGapAxis = {100, 0, 1, "match score difference"};
      histName = path + "scoreGapLeadingTrueMatches";
      histTitle = "Score gap between leading and subleading matches - good matches";
      fScoreGapLeadingTrueMatches = registry->add(histName.c_str(), histTitle.c_str(), {HistType::kTH1F, {scoreGapAxis}});
      histName = path + "scoreGapNonLeadingTrueMatches";
      histTitle = "Score gap between leading and subleading matches - non-leading matches";
      fScoreGapNonLeadingTrueMatches = registry->add(histName.c_str(), histTitle.c_str(), {HistType::kTH1F, {scoreGapAxis}});

      //-
      AxisSpec chi2Axis = {100, 0, 100, "matching #chi^{2}/NDF"};
      AxisSpec scoreAxis = {100, 0, 1, "matching score"};
      int matchTypeMax = static_cast<int>(kMatchTypeUndefined);
      AxisSpec matchTypeAxis = {matchTypeMax, 0, static_cast<double>(matchTypeMax), "match type"};
      histName = path + "matchType";
      histTitle = "Match type";
      fMatchType = registry->add(histName.c_str(), histTitle.c_str(), {HistType::kTH1F, {matchTypeAxis}});
      std::get<std::shared_ptr<TH1>>(fMatchType)->GetXaxis()->SetBinLabel(1, "true (leading)");
      std::get<std::shared_ptr<TH1>>(fMatchType)->GetXaxis()->SetBinLabel(2, "wrong (leading)");
      std::get<std::shared_ptr<TH1>>(fMatchType)->GetXaxis()->SetBinLabel(3, "decay (leading)");
      std::get<std::shared_ptr<TH1>>(fMatchType)->GetXaxis()->SetBinLabel(4, "fake (leading)");
      std::get<std::shared_ptr<TH1>>(fMatchType)->GetXaxis()->SetBinLabel(5, "true (non leading)");
      std::get<std::shared_ptr<TH1>>(fMatchType)->GetXaxis()->SetBinLabel(6, "wrong (non leading)");
      std::get<std::shared_ptr<TH1>>(fMatchType)->GetXaxis()->SetBinLabel(7, "decay (non leading)");
      std::get<std::shared_ptr<TH1>>(fMatchType)->GetXaxis()->SetBinLabel(8, "fake (non leading)");
      histName = path + "matchTypeVsP";
      histTitle = "Match type vs. p";
      fMatchTypeVsP = registry->add(histName.c_str(), histTitle.c_str(), {HistType::kTH2F, {pAxis, matchTypeAxis}});
      std::get<std::shared_ptr<TH2>>(fMatchTypeVsP)->GetYaxis()->SetBinLabel(1, "true (leading)");
      std::get<std::shared_ptr<TH2>>(fMatchTypeVsP)->GetYaxis()->SetBinLabel(2, "wrong (leading)");
      std::get<std::shared_ptr<TH2>>(fMatchTypeVsP)->GetYaxis()->SetBinLabel(3, "decay (leading)");
      std::get<std::shared_ptr<TH2>>(fMatchTypeVsP)->GetYaxis()->SetBinLabel(4, "fake (leading)");
      std::get<std::shared_ptr<TH2>>(fMatchTypeVsP)->GetYaxis()->SetBinLabel(5, "true (non leading)");
      std::get<std::shared_ptr<TH2>>(fMatchTypeVsP)->GetYaxis()->SetBinLabel(6, "wrong (non leading)");
      std::get<std::shared_ptr<TH2>>(fMatchTypeVsP)->GetYaxis()->SetBinLabel(7, "decay (non leading)");
      std::get<std::shared_ptr<TH2>>(fMatchTypeVsP)->GetYaxis()->SetBinLabel(8, "fake (non leading)");
      histName = path + "matchTypeVsPt";
      histTitle = "Match type vs. p_{T}";
      fMatchTypeVsPt = registry->add(histName.c_str(), histTitle.c_str(), {HistType::kTH2F, {ptAxis, matchTypeAxis}});
      std::get<std::shared_ptr<TH2>>(fMatchTypeVsPt)->GetYaxis()->SetBinLabel(1, "true (leading)");
      std::get<std::shared_ptr<TH2>>(fMatchTypeVsPt)->GetYaxis()->SetBinLabel(2, "wrong (leading)");
      std::get<std::shared_ptr<TH2>>(fMatchTypeVsPt)->GetYaxis()->SetBinLabel(3, "decay (leading)");
      std::get<std::shared_ptr<TH2>>(fMatchTypeVsPt)->GetYaxis()->SetBinLabel(4, "fake (leading)");
      std::get<std::shared_ptr<TH2>>(fMatchTypeVsPt)->GetYaxis()->SetBinLabel(5, "true (non leading)");
      std::get<std::shared_ptr<TH2>>(fMatchTypeVsPt)->GetYaxis()->SetBinLabel(6, "wrong (non leading)");
      std::get<std::shared_ptr<TH2>>(fMatchTypeVsPt)->GetYaxis()->SetBinLabel(7, "decay (non leading)");
      std::get<std::shared_ptr<TH2>>(fMatchTypeVsPt)->GetYaxis()->SetBinLabel(8, "fake (non leading)");

      histName = path + "matchChi2VsType";
      histTitle = "Match #chi^{2} vs. match type";
      fMatchChi2VsType = registry->add(histName.c_str(), histTitle.c_str(), {HistType::kTH2F, {matchTypeAxis, chi2Axis}});
      std::get<std::shared_ptr<TH2>>(fMatchChi2VsType)->GetXaxis()->SetBinLabel(1, "true (leading)");
      std::get<std::shared_ptr<TH2>>(fMatchChi2VsType)->GetXaxis()->SetBinLabel(2, "wrong (leading)");
      std::get<std::shared_ptr<TH2>>(fMatchChi2VsType)->GetXaxis()->SetBinLabel(3, "decay (leading)");
      std::get<std::shared_ptr<TH2>>(fMatchChi2VsType)->GetXaxis()->SetBinLabel(4, "fake (leading)");
      std::get<std::shared_ptr<TH2>>(fMatchChi2VsType)->GetXaxis()->SetBinLabel(5, "true (non leading)");
      std::get<std::shared_ptr<TH2>>(fMatchChi2VsType)->GetXaxis()->SetBinLabel(6, "wrong (non leading)");
      std::get<std::shared_ptr<TH2>>(fMatchChi2VsType)->GetXaxis()->SetBinLabel(7, "decay (non leading)");
      std::get<std::shared_ptr<TH2>>(fMatchChi2VsType)->GetXaxis()->SetBinLabel(8, "fake (non leading)");
      histName = path + "matchChi2VsTypeVsP";
      histTitle = "Match #chi^{2} vs. match type vs. p";
      fMatchChi2VsTypeVsP = registry->add(histName.c_str(), histTitle.c_str(), {HistType::kTH3F, {pAxis, matchTypeAxis, chi2Axis}});
      std::get<std::shared_ptr<TH3>>(fMatchChi2VsTypeVsP)->GetYaxis()->SetBinLabel(1, "true (leading)");
      std::get<std::shared_ptr<TH3>>(fMatchChi2VsTypeVsP)->GetYaxis()->SetBinLabel(2, "wrong (leading)");
      std::get<std::shared_ptr<TH3>>(fMatchChi2VsTypeVsP)->GetYaxis()->SetBinLabel(3, "decay (leading)");
      std::get<std::shared_ptr<TH3>>(fMatchChi2VsTypeVsP)->GetYaxis()->SetBinLabel(4, "fake (leading)");
      std::get<std::shared_ptr<TH3>>(fMatchChi2VsTypeVsP)->GetYaxis()->SetBinLabel(5, "true (non leading)");
      std::get<std::shared_ptr<TH3>>(fMatchChi2VsTypeVsP)->GetYaxis()->SetBinLabel(6, "wrong (non leading)");
      std::get<std::shared_ptr<TH3>>(fMatchChi2VsTypeVsP)->GetYaxis()->SetBinLabel(7, "decay (non leading)");
      std::get<std::shared_ptr<TH3>>(fMatchChi2VsTypeVsP)->GetYaxis()->SetBinLabel(8, "fake (non leading)");
      histName = path + "matchChi2VsTypeVsPt";
      histTitle = "Match #chi^{2} vs. match type vs. p_{T}";
      fMatchChi2VsTypeVsPt = registry->add(histName.c_str(), histTitle.c_str(), {HistType::kTH3F, {ptAxis, matchTypeAxis, chi2Axis}});
      std::get<std::shared_ptr<TH3>>(fMatchChi2VsTypeVsPt)->GetYaxis()->SetBinLabel(1, "true (leading)");
      std::get<std::shared_ptr<TH3>>(fMatchChi2VsTypeVsPt)->GetYaxis()->SetBinLabel(2, "wrong (leading)");
      std::get<std::shared_ptr<TH3>>(fMatchChi2VsTypeVsPt)->GetYaxis()->SetBinLabel(3, "decay (leading)");
      std::get<std::shared_ptr<TH3>>(fMatchChi2VsTypeVsPt)->GetYaxis()->SetBinLabel(4, "fake (leading)");
      std::get<std::shared_ptr<TH3>>(fMatchChi2VsTypeVsPt)->GetYaxis()->SetBinLabel(5, "true (non leading)");
      std::get<std::shared_ptr<TH3>>(fMatchChi2VsTypeVsPt)->GetYaxis()->SetBinLabel(6, "wrong (non leading)");
      std::get<std::shared_ptr<TH3>>(fMatchChi2VsTypeVsPt)->GetYaxis()->SetBinLabel(7, "decay (non leading)");
      std::get<std::shared_ptr<TH3>>(fMatchChi2VsTypeVsPt)->GetYaxis()->SetBinLabel(8, "fake (non leading)");
      //-
      histName = path + "matchScoreVsType";
      histTitle = "Match score vs. match type";
      fMatchScoreVsType = registry->add(histName.c_str(), histTitle.c_str(), {HistType::kTH2F, {matchTypeAxis, scoreAxis}});
      std::get<std::shared_ptr<TH2>>(fMatchScoreVsType)->GetXaxis()->SetBinLabel(1, "true (leading)");
      std::get<std::shared_ptr<TH2>>(fMatchScoreVsType)->GetXaxis()->SetBinLabel(2, "wrong (leading)");
      std::get<std::shared_ptr<TH2>>(fMatchScoreVsType)->GetXaxis()->SetBinLabel(3, "decay (leading)");
      std::get<std::shared_ptr<TH2>>(fMatchScoreVsType)->GetXaxis()->SetBinLabel(4, "fake (leading)");
      std::get<std::shared_ptr<TH2>>(fMatchScoreVsType)->GetXaxis()->SetBinLabel(5, "true (non leading)");
      std::get<std::shared_ptr<TH2>>(fMatchScoreVsType)->GetXaxis()->SetBinLabel(6, "wrong (non leading)");
      std::get<std::shared_ptr<TH2>>(fMatchScoreVsType)->GetXaxis()->SetBinLabel(7, "decay (non leading)");
      std::get<std::shared_ptr<TH2>>(fMatchScoreVsType)->GetXaxis()->SetBinLabel(8, "fake (non leading)");
      histName = path + "matchScoreVsTypeVsP";
      histTitle = "Match score vs. match type vs. p";
      fMatchScoreVsTypeVsP = registry->add(histName.c_str(), histTitle.c_str(), {HistType::kTH3F, {pAxis, matchTypeAxis, scoreAxis}});
      std::get<std::shared_ptr<TH3>>(fMatchScoreVsTypeVsP)->GetYaxis()->SetBinLabel(1, "true (leading)");
      std::get<std::shared_ptr<TH3>>(fMatchScoreVsTypeVsP)->GetYaxis()->SetBinLabel(2, "wrong (leading)");
      std::get<std::shared_ptr<TH3>>(fMatchScoreVsTypeVsP)->GetYaxis()->SetBinLabel(3, "decay (leading)");
      std::get<std::shared_ptr<TH3>>(fMatchScoreVsTypeVsP)->GetYaxis()->SetBinLabel(4, "fake (leading)");
      std::get<std::shared_ptr<TH3>>(fMatchScoreVsTypeVsP)->GetYaxis()->SetBinLabel(5, "true (non leading)");
      std::get<std::shared_ptr<TH3>>(fMatchScoreVsTypeVsP)->GetYaxis()->SetBinLabel(6, "wrong (non leading)");
      std::get<std::shared_ptr<TH3>>(fMatchScoreVsTypeVsP)->GetYaxis()->SetBinLabel(7, "decay (non leading)");
      std::get<std::shared_ptr<TH3>>(fMatchScoreVsTypeVsP)->GetYaxis()->SetBinLabel(8, "fake (non leading)");
      histName = path + "matchScoreVsTypeVsPt";
      histTitle = "Match score vs. match type vs. p_{T}";
      fMatchScoreVsTypeVsPt = registry->add(histName.c_str(), histTitle.c_str(), {HistType::kTH3F, {ptAxis, matchTypeAxis, scoreAxis}});
      std::get<std::shared_ptr<TH3>>(fMatchScoreVsTypeVsPt)->GetYaxis()->SetBinLabel(1, "true (leading)");
      std::get<std::shared_ptr<TH3>>(fMatchScoreVsTypeVsPt)->GetYaxis()->SetBinLabel(2, "wrong (leading)");
      std::get<std::shared_ptr<TH3>>(fMatchScoreVsTypeVsPt)->GetYaxis()->SetBinLabel(3, "decay (leading)");
      std::get<std::shared_ptr<TH3>>(fMatchScoreVsTypeVsPt)->GetYaxis()->SetBinLabel(4, "fake (leading)");
      std::get<std::shared_ptr<TH3>>(fMatchScoreVsTypeVsPt)->GetYaxis()->SetBinLabel(5, "true (non leading)");
      std::get<std::shared_ptr<TH3>>(fMatchScoreVsTypeVsPt)->GetYaxis()->SetBinLabel(6, "wrong (non leading)");
      std::get<std::shared_ptr<TH3>>(fMatchScoreVsTypeVsPt)->GetYaxis()->SetBinLabel(7, "decay (non leading)");
      std::get<std::shared_ptr<TH3>>(fMatchScoreVsTypeVsPt)->GetYaxis()->SetBinLabel(8, "fake (non leading)");

      AxisSpec prodScoreAxis = {100, 0, 1, "matching score (prod)"};
      histName = path + "matchScoreVsProd";
      histTitle = "Match score vs. production";
      fMatchScoreVsProd = registry->add(histName.c_str(), histTitle.c_str(), {HistType::kTH2F, {prodScoreAxis, scoreAxis}});
      histName = path + "trueMatchScoreVsProd";
      histTitle = "Match score vs. production - true match";
      fTrueMatchScoreVsProd = registry->add(histName.c_str(), histTitle.c_str(), {HistType::kTH2F, {prodScoreAxis, scoreAxis}});

      AxisSpec prodChi2Axis = {100, 0, 100, "matching #chi^{2}/NDF (prod)"};
      histName = path + "matchChi2VsProd";
      histTitle = "Match #chi^{2} vs. production";
      fMatchChi2VsProd = registry->add(histName.c_str(), histTitle.c_str(), {HistType::kTH2F, {prodChi2Axis, chi2Axis}});
      histName = path + "trueMatchChi2VsProd";
      histTitle = "Match #chi^{2} vs. production - true match";
      fTrueMatchChi2VsProd = registry->add(histName.c_str(), histTitle.c_str(), {HistType::kTH2F, {{100, 0, 10, "matching #chi^{2} (prod)"}, {100, 0, 10, "matching #chi^{2}"}}});
    }
  };
  std::unique_ptr<MatchingPlotter> fChi2MatchingPlotter;
  std::map<std::string, std::unique_ptr<HistogramRegistry>> fMatchingHistogramRegistries;
  std::map<std::string, std::unique_ptr<MatchingPlotter>> fMatchingPlotters;
  std::unique_ptr<MatchingPlotter> fTaggedMuonsMatchingPlotter;
  std::unique_ptr<MatchingPlotter> fSelectedMuonsMatchingPlotter;

  CollisionInfos fCollisionInfos;

  template <typename BC>
  void initCCDB(BC const& bc)
  {
    if (mRunNumber == bc.runNumber())
      return;

    mRunNumber = bc.runNumber();
    std::map<std::string, std::string> metadata;
    auto soreor = o2::ccdb::BasicCCDBManager::getRunDuration(fCCDBApi, mRunNumber);
    auto ts = soreor.first;
    auto grpmag = fCCDBApi.retrieveFromTFileAny<o2::parameters::GRPMagField>(grpmagPath, metadata, ts);
    o2::base::Propagator::initFieldFromGRP(grpmag);
    LOGF(info, "Set field for muons");
    VarManager::SetupMuonMagField();
    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      ccdbManager->get<TGeoManager>(geoPath);
    }
    o2::mch::TrackExtrap::setField();
    auto* fieldB = static_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField());
    if (fieldB) {
      double centerMFT[3] = {0, 0, -61.4}; // Field at center of MFT
      mBzAtMftCenter = fieldB->getBz(centerMFT);
      // std::cout << "fieldB: " << (void*)fieldB << std::endl;
    }
  }

  void CreateMatchingHistosMC()
  {
    AxisSpec chi2Axis = {1000, 0, 1000, "chi^{2}"};
    AxisSpec chi2AxisSmall = {200, 0, 100, "chi^{2}"};
    AxisSpec pAxis = {1000, 0, 100, "p (GeV/c)"};
    AxisSpec pTAxis = {100, 0, 10, "p_{T} (GeV/c)"};
    AxisSpec etaAxis = {100, -4, -2, "#eta"};
    AxisSpec phiAxis = {90, -180, 180, "#phi (degrees)"};
    std::string histPath = "matching/MC/";

    AxisSpec trackPositionXAtMFTAxis = {100, -15, 15, "MFT x (cm)"};
    AxisSpec trackPositionYAtMFTAxis = {100, -15, 15, "MFT y (cm)"};
    registry.add((histPath + "pairedMCHTracksAtMFT").c_str(), "Paired MCH tracks position at MFT end", {HistType::kTH2F, {trackPositionXAtMFTAxis, trackPositionYAtMFTAxis}});
    registry.add((histPath + "pairedMFTTracksAtMFT").c_str(), "Paired MFT tracks position at MFT end", {HistType::kTH2F, {trackPositionXAtMFTAxis, trackPositionYAtMFTAxis}});
    registry.add((histPath + "selectedMCHTracksAtMFT").c_str(), "Selected MCH tracks position at MFT end", {HistType::kTH2F, {trackPositionXAtMFTAxis, trackPositionYAtMFTAxis}});
    registry.add((histPath + "selectedMCHTracksAtMFTTrue").c_str(), "Selected MCH tracks position at MFT end - true", {HistType::kTH2F, {trackPositionXAtMFTAxis, trackPositionYAtMFTAxis}});
    registry.add((histPath + "selectedMCHTracksAtMFTFake").c_str(), "Selected MCH tracks position at MFT end - fake", {HistType::kTH2F, {trackPositionXAtMFTAxis, trackPositionYAtMFTAxis}});

    fChi2MatchingPlotter = std::make_unique<MatchingPlotter>(histPath + "Prod/", &registryMatching);
    int registryIndex = 0;
    for (const auto& [label, func] : matchingChi2Functions) {
      fMatchingPlotters[label] = std::make_unique<MatchingPlotter>(histPath + label + "/", registryMatchingVec[registryIndex]);
      registryIndex += 1;
    }
    for (const auto& [label, response] : matchingMlResponses) {
      fMatchingPlotters[label] = std::make_unique<MatchingPlotter>(histPath + label + "/", (registryMatchingVec[registryIndex]));
      registryIndex += 1;
    }

    fTaggedMuonsMatchingPlotter = std::make_unique<MatchingPlotter>(histPath + "Tagged/", &registryMatching);
    fSelectedMuonsMatchingPlotter = std::make_unique<MatchingPlotter>(histPath + "Selected/", &registryMatching);
  }

  void CreateDimuonHistos()
  {
    AxisSpec invMassAxis = {400, 1, 5, "M_{#mu^{+}#mu^{-}} (GeV/c^{2})"};
    AxisSpec invMassCorrelationAxis = {400, 0, 8, "M_{#mu^{+}#mu^{-}} (GeV/c^{2})"};
    AxisSpec invMassAxisFull = {5000, 0, 100, "M_{#mu^{+}#mu^{-}} (GeV/c^{2})"};
    int matchTypeCombMax = (static_cast<int>(kMatchTypeTrueNonLeading) - 1) * 10 + static_cast<int>(kMatchTypeTrueNonLeading) - 1;
    AxisSpec matchTypeAxis = {matchTypeCombMax + 1, 0, static_cast<double>(matchTypeCombMax + 1), "match type"};

    // MCH-MID tracks with MCH acceptance cuts
    registryDimuon.add("dimuon/invariantMass_MuonKine_MuonCuts", "#mu^{+}#mu^{-} invariant mass (muon cuts)", {HistType::kTH1F, {invMassAxis}});
    // MCH-MID tracks with MFT acceptance cuts
    registryDimuon.add("dimuon/invariantMass_MuonKine_GlobalMuonCuts", "#mu^{+}#mu^{-} invariant mass (global muon cuts)", {HistType::kTH1F, {invMassAxis}});
    // MCH-MID tracks with MFT acceptance cuts vs. muon tracks match type
    registryDimuon.add("dimuon/MC/invariantMass_MuonKine_GlobalMuonCuts_vs_match_type", "#mu^{+}#mu^{-} invariant mass vs. match tye (global muon cuts)", {HistType::kTH2F, {invMassAxis, matchTypeAxis}});
    // MCH-MID tracks with MFT acceptance cuts, good matches
    registryDimuon.add("dimuon/invariantMass_MuonKine_GlobalMuonCuts_GoodMatches", "#mu^{+}#mu^{-} invariant mass (global muon cuts, good matches)", {HistType::kTH1F, {invMassAxis}});
    // MCH-MID tracks with MFT acceptance cuts, good matches + paired muons
    registryDimuon.add("dimuon/MC/invariantMass_MuonKine_GlobalMuonCuts_GoodMatches_vs_match_type", "#mu^{+}#mu^{-} invariant mass vs. match tye (global muon cuts, good matches)", {HistType::kTH2F, {invMassAxis, matchTypeAxis}});

    // scaled kinematics (Hiroshima method)
    // MFT-MCH-MID tracks with MFT acceptance cuts
    registryDimuon.add("dimuon/invariantMass_ScaledMftKine_GlobalMuonCuts", "#mu^{+}#mu^{-} invariant mass (global muon cuts, rescaled MFT)", {HistType::kTH1F, {invMassAxis}});
    // MCH-MID tracks with MFT acceptance cuts vs. muon tracks match type
    registryDimuon.add("dimuon/MC/invariantMass_ScaledMftKine_GlobalMuonCuts_vs_match_type", "#mu^{+}#mu^{-} invariant mass vs. match tye (global muon cuts, rescaled MFT)", {HistType::kTH2F, {invMassAxis, matchTypeAxis}});
    // MFT-MCH-MID tracks with MFT acceptance cuts, good matches
    registryDimuon.add("dimuon/invariantMass_ScaledMftKine_GlobalMuonCuts_GoodMatches", "#mu^{+}#mu^{-} invariant mass (global muon cuts, rescaled MFT, good matches)", {HistType::kTH1F, {invMassAxis}});
    // MFT-MCH-MID tracks with MFT acceptance cuts vs. muon tracks match type, good matches
    registryDimuon.add("dimuon/MC/invariantMass_ScaledMftKine_GlobalMuonCuts_GoodMatches_vs_match_type", "#mu^{+}#mu^{-} invariant mass vs. match tye (global muon cuts, rescaled MFT, good matches)", {HistType::kTH2F, {invMassAxis, matchTypeAxis}});
  }

  void InitMatchingFunctions()
  {
    using SMatrix55Std = ROOT::Math::SMatrix<double, 5>;
    using SMatrix55Sym = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;

    using SVector2 = ROOT::Math::SVector<double, 2>;
    using SVector4 = ROOT::Math::SVector<double, 4>;
    using SVector5 = ROOT::Math::SVector<double, 5>;

    using SMatrix44 = ROOT::Math::SMatrix<double, 4>;
    using SMatrix45 = ROOT::Math::SMatrix<double, 4, 5>;
    using SMatrix22 = ROOT::Math::SMatrix<double, 2>;
    using SMatrix25 = ROOT::Math::SMatrix<double, 2, 5>;

    // Define built-in matching functions
    //________________________________________________________________________________
    mMatchingFunctionMap["matchALL"] = [](const o2::dataformats::GlobalFwdTrack& mchTrack, const o2::track::TrackParCovFwd& mftTrack) -> std::tuple<double, int> {
      // Match two tracks evaluating all parameters: X,Y, phi, tanl & q/pt

      SMatrix55Sym H_k, V_k;
      SVector5 m_k(mftTrack.getX(), mftTrack.getY(), mftTrack.getPhi(),
                   mftTrack.getTanl(), mftTrack.getInvQPt()),
        r_k_kminus1;
      SVector5 GlobalMuonTrackParameters = mchTrack.getParameters();
      SMatrix55Sym GlobalMuonTrackCovariances = mchTrack.getCovariances();
      V_k(0, 0) = mftTrack.getCovariances()(0, 0);
      V_k(1, 1) = mftTrack.getCovariances()(1, 1);
      V_k(2, 2) = mftTrack.getCovariances()(2, 2);
      V_k(3, 3) = mftTrack.getCovariances()(3, 3);
      V_k(4, 4) = mftTrack.getCovariances()(4, 4);
      H_k(0, 0) = 1.0;
      H_k(1, 1) = 1.0;
      H_k(2, 2) = 1.0;
      H_k(3, 3) = 1.0;
      H_k(4, 4) = 1.0;

      // Covariance of residuals
      SMatrix55Std invResCov = (V_k + ROOT::Math::Similarity(H_k, GlobalMuonTrackCovariances));
      invResCov.Invert();

      // Update Parameters
      r_k_kminus1 = m_k - H_k * GlobalMuonTrackParameters; // Residuals of prediction

      auto matchChi2Track = ROOT::Math::Similarity(r_k_kminus1, invResCov);

      // return chi2 and NDF
      return {matchChi2Track, 5};
    };

    //________________________________________________________________________________
    mMatchingFunctionMap["matchXYPhiTanl"] = [](const o2::dataformats::GlobalFwdTrack& mchTrack, const o2::track::TrackParCovFwd& mftTrack) -> std::tuple<double, int> {
      // Match two tracks evaluating positions & angles

      SMatrix45 H_k;
      SMatrix44 V_k;
      SVector4 m_k(mftTrack.getX(), mftTrack.getY(), mftTrack.getPhi(),
                   mftTrack.getTanl()),
        r_k_kminus1;
      SVector5 GlobalMuonTrackParameters = mchTrack.getParameters();
      SMatrix55Sym GlobalMuonTrackCovariances = mchTrack.getCovariances();
      V_k(0, 0) = mftTrack.getCovariances()(0, 0);
      V_k(1, 1) = mftTrack.getCovariances()(1, 1);
      V_k(2, 2) = mftTrack.getCovariances()(2, 2);
      V_k(3, 3) = mftTrack.getCovariances()(3, 3);
      H_k(0, 0) = 1.0;
      H_k(1, 1) = 1.0;
      H_k(2, 2) = 1.0;
      H_k(3, 3) = 1.0;

      // Covariance of residuals
      SMatrix44 invResCov = (V_k + ROOT::Math::Similarity(H_k, GlobalMuonTrackCovariances));
      invResCov.Invert();

      // Residuals of prediction
      r_k_kminus1 = m_k - H_k * GlobalMuonTrackParameters;

      auto matchChi2Track = ROOT::Math::Similarity(r_k_kminus1, invResCov);

      // return chi2 and NDF
      return {matchChi2Track, 4};
    };

    //________________________________________________________________________________
    mMatchingFunctionMap["matchXY"] = [](const o2::dataformats::GlobalFwdTrack& mchTrack, const o2::track::TrackParCovFwd& mftTrack) -> std::tuple<double, int> {
      // Calculate Matching Chi2 - X and Y positions

      SMatrix25 H_k;
      SMatrix22 V_k;
      SVector2 m_k(mftTrack.getX(), mftTrack.getY()), r_k_kminus1;
      SVector5 GlobalMuonTrackParameters = mchTrack.getParameters();
      SMatrix55Sym GlobalMuonTrackCovariances = mchTrack.getCovariances();
      V_k(0, 0) = mftTrack.getCovariances()(0, 0);
      V_k(1, 1) = mftTrack.getCovariances()(1, 1);
      H_k(0, 0) = 1.0;
      H_k(1, 1) = 1.0;

      // Covariance of residuals
      SMatrix22 invResCov = (V_k + ROOT::Math::Similarity(H_k, GlobalMuonTrackCovariances));
      invResCov.Invert();

      // Residuals of prediction
      r_k_kminus1 = m_k - H_k * GlobalMuonTrackParameters;
      auto matchChi2Track = ROOT::Math::Similarity(r_k_kminus1, invResCov);

      // return reduced chi2
      return {matchChi2Track, 2};
    };
  }

  void init(o2::framework::InitContext&)
  {
    // Load geometry
    ccdbManager->setURL(ccdburl);
    ccdbManager->setCaching(true);
    ccdbManager->setLocalObjectValidityChecking();
    fCCDBApi.init(ccdburl);
    mRunNumber = 0;

    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      LOGF(info, "Load geometry from CCDB");
      ccdbManager->get<TGeoManager>(geoPath);
    }

    // Matching functions
    InitMatchingFunctions();
    for (size_t funcId = 0; funcId < sChi2FunctionsNum; funcId++) {
      auto label = fConfigChi2MatchingOptions.fFunctionLabel[funcId].value;
      auto funcName = fConfigChi2MatchingOptions.fFunctionName[funcId].value;
      auto scoreMin = fConfigChi2MatchingOptions.fMatchingScoreCut[funcId].value;
      auto matchingPlaneZ = fConfigChi2MatchingOptions.fMatchingPlaneZ[funcId].value;
      auto extrapMethod = fConfigChi2MatchingOptions.fMatchingExtrapMethod[funcId].value;

      if (label == "" || funcName == "")
        break;

      matchingChi2Functions[label] = funcName;

      matchingScoreCuts[label] = scoreMin;
      matchingPlanesZ[label] = matchingPlaneZ;
      matchingExtrapMethod[label] = extrapMethod;
    }

    // Matching ML models
    // TODO : for now we use hard coded values since the current models use 1 pT bin
    binsPtMl = {-1e-6, 1000.0};
    cutValues = {0.0};
    cutDirMl = {cuts_ml::CutNot};
    o2::framework::LabeledArray<double> mycutsMl(cutValues.data(), 1, 1, std::vector<std::string>{"pT bin 0"}, std::vector<std::string>{"score"});

    for (size_t modelId = 0; modelId < sMLModelsNum; modelId++) {
      auto label = fConfigMlOptions.fModelLabel[modelId].value;
      auto modelPaths = fConfigMlOptions.fModelPathsCCDB[modelId].value;
      auto inputFeatures = fConfigMlOptions.fInputFeatures[modelId].value;
      auto modelNames = fConfigMlOptions.fModelNames[modelId].value;
      auto scoreMin = fConfigMlOptions.fMatchingScoreCut[modelId].value;
      auto matchingPlaneZ = fConfigMlOptions.fMatchingPlaneZ[modelId].value;
      auto extrapMethod = fConfigMlOptions.fMatchingExtrapMethod[modelId].value;

      if (label == "" || modelPaths.empty() || inputFeatures.empty() || modelNames.empty())
        break;

      matchingMlResponses[label].configure(binsPtMl, mycutsMl, cutDirMl, 1);
      matchingMlResponses[label].setModelPathsCCDB(modelNames, fCCDBApi, modelPaths, fConfigCCDB.fConfigNoLaterThan.value);
      matchingMlResponses[label].cacheInputFeaturesIndices(inputFeatures);
      matchingMlResponses[label].init();

      matchingScoreCuts[label] = scoreMin;
      matchingPlanesZ[label] = matchingPlaneZ;
      matchingExtrapMethod[label] = extrapMethod;
    }

    int nTrackTypes = static_cast<int>(o2::aod::fwdtrack::ForwardTrackTypeEnum::MCHStandaloneTrack) + 1;
    AxisSpec trackTypeAxis = {static_cast<int>(nTrackTypes), 0.0, static_cast<double>(nTrackTypes), "track type"};
    registry.add("nTracksPerType", "Number of tracks per type", {HistType::kTH1F, {trackTypeAxis}});

    AxisSpec tracksMultiplicityAxis = {10000, 0, 10000, "tracks multiplicity"};
    registry.add("tracksMultiplicityMFT", "MFT tracks multiplicity", {HistType::kTH1F, {tracksMultiplicityAxis}});
    registry.add("tracksMultiplicityMCH", "MCH tracks multiplicity", {HistType::kTH1F, {tracksMultiplicityAxis}});

    CreateMatchingHistosMC();
    CreateDimuonHistos();
  }

  template <class T, class C>
  bool pDCACut(const T& mchTrack, const C& collision, double nSigmaPDCA)
  {
    static const double sigmaPDCA23 = 80.;
    static const double sigmaPDCA310 = 54.;
    static const double relPRes = 0.0004;
    static const double slopeRes = 0.0005;

    double thetaAbs = TMath::ATan(mchTrack.rAtAbsorberEnd() / 505.) * TMath::RadToDeg();

    // propagate muon track to vertex
    auto mchTrackAtVertex = VarManager::PropagateMuon(mchTrack, collision, VarManager::kToVertex);

    // double pUncorr = mchTrack.p();
    double p = mchTrackAtVertex.getP();

    double pDCA = mchTrack.pDca();
    double sigmaPDCA = (thetaAbs < 3) ? sigmaPDCA23 : sigmaPDCA310;
    double nrp = nSigmaPDCA * relPRes * p;
    double pResEffect = sigmaPDCA / (1. - nrp / (1. + nrp));
    double slopeResEffect = 535. * slopeRes * p;
    double sigmaPDCAWithRes = TMath::Sqrt(pResEffect * pResEffect + slopeResEffect * slopeResEffect);
    if (pDCA > nSigmaPDCA * sigmaPDCAWithRes) {
      return false;
    }

    return true;
  }

  template <class T, class C>
  bool IsGoodMuon(const T& mchTrack, const C& collision,
                  double chi2Cut,
                  double pCut,
                  double pTCut,
                  std::array<double, 2> etaCut,
                  std::array<double, 2> rAbsCut,
                  double nSigmaPdcaCut)
  {
    // chi2 cut
    if (mchTrack.chi2() > chi2Cut)
      return false;

    // momentum cut
    if (mchTrack.p() < pCut) {
      return false; // skip low-momentum tracks
    }

    // transverse momentum cut
    if (mchTrack.pt() < pTCut) {
      return false; // skip low-momentum tracks
    }

    // Eta cut
    double eta = mchTrack.eta();
    if ((eta < etaCut[0] || eta > etaCut[1])) {
      return false;
    }

    // RAbs cut
    double rAbs = mchTrack.rAtAbsorberEnd();
    if ((rAbs < rAbsCut[0] || rAbs > rAbsCut[1])) {
      return false;
    }

    // pDCA cut
    if (!pDCACut(mchTrack, collision, nSigmaPdcaCut)) {
      return false;
    }

    return true;
  }

  template <class T, class C>
  bool IsGoodMuon(const T& muonTrack, const C& collision)
  {
    return IsGoodMuon(muonTrack, collision, fTrackChi2MchUp, fPMchLow, fPtMchLow, {fEtaMchLow, fEtaMchUp}, {fRabsLow, fRabsUp}, fSigmaPdcaUp);
  }

  template <class T, class C>
  bool IsGoodGlobalMuon(const T& muonTrack, const C& collision)
  {
    return IsGoodMuon(muonTrack, collision, fTrackChi2MchUp, fPMchLow, fPtMchLow, {fEtaMftLow, fEtaMftUp}, {fRabsLow, fRabsUp}, fSigmaPdcaUp);
  }

  template <class T>
  bool IsGoodMFT(const T& mftTrack,
                 double chi2Cut,
                 int nClustersCut)
  {
    // std::cout << std::format("Checking MFT track") << std::endl;
    // std::cout << std::format("    chi2={}", mftTrack.chi2()) << std::endl;
    // std::cout << std::format("    nClusters={}", mftTrack.nClusters()) << std::endl;
    //  chi2 cut
    if (mftTrack.chi2() > chi2Cut)
      return false;

    // number of clusters cut
    if (mftTrack.nClusters() < nClustersCut)
      return false;

    return true;
  }

  template <class T>
  bool IsGoodMFT(const T& mftTrack)
  {
    return IsGoodMFT(mftTrack, fTrackChi2MftUp, fTrackNClustMftLow);
  }

  template <class TMUON>
  bool IsGoodGlobalMatching(const TMUON& muonTrack,
                            double matchingScore,
                            double matchingScoreCut)
  {
    if (static_cast<int>(muonTrack.trackType()) > 2)
      return false;

    // MFT-MCH matching score cut
    if (matchingScore < matchingScoreCut)
      return false;

    return true;
  }

  template <class TMUON>
  bool IsGoodGlobalMatching(const TMUON& muonTrack, double matchingScore)
  {
    return IsGoodGlobalMatching(muonTrack, matchingScore, fMatchingChi2ScoreMftMchLow);
  }

  template <class TMUON>
  bool IsTrueGlobalMatching(const TMUON& muonTrack, const std::vector<std::pair<int64_t, int64_t>>& matchablePairs)
  {
    if (static_cast<int>(muonTrack.trackType()) > 2)
      return false;

    int64_t mchTrackId = static_cast<int64_t>(muonTrack.matchMCHTrackId());
    int64_t mftTrackId = static_cast<int64_t>(muonTrack.matchMFTTrackId());

    std::pair<int64_t, int64_t> trackIndexes = std::make_pair(mchTrackId, mftTrackId);

    return (std::find(matchablePairs.begin(), matchablePairs.end(), trackIndexes) != matchablePairs.end());
  }

  bool IsMatchableMCH(int64_t mchTrackId, const std::vector<std::pair<int64_t, int64_t>>& matchablePairs)
  {
    for (auto [id1, id2] : matchablePairs) {
      if (mchTrackId == id1)
        return true;
    }
    return false;
  }

  std::optional<std::pair<int64_t, int64_t>> GetMatchablePairForMCH(int64_t mchTrackId, const std::vector<std::pair<int64_t, int64_t>>& matchablePairs)
  {
    for (auto pair : matchablePairs) {
      if (mchTrackId == pair.first)
        return pair;
    }
    return {};
  }

  template <typename T>
  o2::dataformats::GlobalFwdTrack FwdToTrackPar(const T& track)
  {
    double chi2 = track.chi2();
    SMatrix5 tpars(track.x(), track.y(), track.phi(), track.tgl(), track.signed1Pt());
    std::vector<double> v1{0, 0, 0, 0, 0,
                           0, 0, 0, 0, 0,
                           0, 0, 0, 0, 0};
    SMatrix55 tcovs(v1.begin(), v1.end());
    o2::track::TrackParCovFwd trackparCov{track.z(), tpars, tcovs, chi2};
    o2::dataformats::GlobalFwdTrack fwdtrack;
    fwdtrack.setParameters(trackparCov.getParameters());
    fwdtrack.setZ(trackparCov.getZ());
    fwdtrack.setCovariances(trackparCov.getCovariances());
    return fwdtrack;
  }

  template <typename T, typename C>
  o2::dataformats::GlobalFwdTrack FwdToTrackPar(const T& track, const C& cov)
  {
    double chi2 = track.chi2();
    SMatrix5 tpars(track.x(), track.y(), track.phi(), track.tgl(), track.signed1Pt());
    std::vector<double> v1{cov.cXX(), cov.cXY(), cov.cYY(), cov.cPhiX(), cov.cPhiY(),
                           cov.cPhiPhi(), cov.cTglX(), cov.cTglY(), cov.cTglPhi(), cov.cTglTgl(),
                           cov.c1PtX(), cov.c1PtY(), cov.c1PtPhi(), cov.c1PtTgl(), cov.c1Pt21Pt2()};
    SMatrix55 tcovs(v1.begin(), v1.end());
    o2::track::TrackParCovFwd trackparCov{track.z(), tpars, tcovs, chi2};
    o2::dataformats::GlobalFwdTrack fwdtrack;
    fwdtrack.setParameters(trackparCov.getParameters());
    fwdtrack.setZ(trackparCov.getZ());
    fwdtrack.setCovariances(trackparCov.getCovariances());
    return fwdtrack;
  }

  o2::dataformats::GlobalFwdTrack PropagateToZMCH(const o2::dataformats::GlobalFwdTrack& muon, const double z)
  {
    auto mchTrack = mExtrap.FwdtoMCH(muon);

    float absFront = -90.f;
    float absBack = -505.f;

    if (muon.getZ() < absBack && z > absFront) {
      // extrapolation through the absorber in the upstream direction
      o2::mch::TrackExtrap::extrapToVertexWithoutBranson(mchTrack, z);
    } else {
      // all other cases
      o2::mch::TrackExtrap::extrapToZCov(mchTrack, z);
    }

    auto proptrack = mExtrap.MCHtoFwd(mchTrack);
    o2::dataformats::GlobalFwdTrack propmuon;
    propmuon.setParameters(proptrack.getParameters());
    propmuon.setZ(proptrack.getZ());
    propmuon.setCovariances(proptrack.getCovariances());

    return propmuon;
  }

  template <typename T>
  o2::dataformats::GlobalFwdTrack PropagateToZMCH(const T& muon, const double z)
  {
    double chi2 = muon.chi2();
    SMatrix5 tpars(muon.x(), muon.y(), muon.phi(), muon.tgl(), muon.signed1Pt());
    std::vector<double> v1{muon.cXX(), muon.cXY(), muon.cYY(), muon.cPhiX(), muon.cPhiY(),
                           muon.cPhiPhi(), muon.cTglX(), muon.cTglY(), muon.cTglPhi(), muon.cTglTgl(),
                           muon.c1PtX(), muon.c1PtY(), muon.c1PtPhi(), muon.c1PtTgl(), muon.c1Pt21Pt2()};
    SMatrix55 tcovs(v1.begin(), v1.end());
    o2::track::TrackParCovFwd fwdtrack{muon.z(), tpars, tcovs, chi2};
    o2::dataformats::GlobalFwdTrack track;
    track.setParameters(fwdtrack.getParameters());
    track.setZ(fwdtrack.getZ());
    track.setCovariances(fwdtrack.getCovariances());

    return PropagateToZMCH(track, z);
  }

  o2::dataformats::GlobalFwdTrack PropagateToZMFT(const o2::dataformats::GlobalFwdTrack& mftTrack, const double z)
  {
    o2::dataformats::GlobalFwdTrack fwdtrack{mftTrack};
    fwdtrack.propagateToZ(z, mBzAtMftCenter);
    return fwdtrack;
  }

  template <typename TMFT, typename CMFT>
  o2::dataformats::GlobalFwdTrack PropagateToZMFT(const TMFT& mftTrack, const CMFT& mftCov, const double z)
  {
    o2::dataformats::GlobalFwdTrack fwdtrack = FwdToTrackPar(mftTrack, mftCov);
    return PropagateToZMFT(fwdtrack, z);
  }

  // method 0: standard extrapolation
  // method 1: MFT extrapolation using MCH momentum
  // method 2: MCH track extrapolation constrained to the first MFT track point, MFT extrapolation using MCH momentum
  // method 3: MCH track extrapolation constrained to the collision point, MFT extrapolation using MCH momentum
  template <typename TMCH, typename TMFT, typename CMFT, typename C>
  o2::dataformats::GlobalFwdTrack PropagateToMatchingPlaneMCH(const TMCH& mchTrack, const TMFT& mftTrack, const CMFT& mftTrackCov, const C& collision, const double z, int method)
  {
    if (method == 0 || method == 1) {
      // simple extrapolation upstream through the absorber
      return PropagateToZMCH(mchTrack, z);
    }

    if (method == 2) {
      // extrapolation to the first MFT point and then back to the matching plane
      auto mftTrackPar = FwdToTrackPar(mftTrack, mftTrackCov);
      // std::cout << std::format("[PropagateToMatchingPlaneMCH] extrapolating to MFT: x={:0.3f} y={:0.3f} z={:0.3f}", mftTrackPar.getX(), mftTrackPar.getY(), mftTrackPar.getZ()) << std::endl;
      auto mchTrackAtMFT = PropagateToVertexMCH(FwdToTrackPar(mchTrack, mchTrack),
                                                mftTrackPar.getX(), mftTrackPar.getY(), mftTrackPar.getZ(),
                                                mftTrackPar.getSigma2X(), mftTrackPar.getSigma2Y());
      // std::cout << std::format("[PropagateToMatchingPlaneMCH] extrapolating to z={:0.3f}", z) << std::endl;
      return PropagateToZMCH(mchTrackAtMFT, z);
    }

    if (method == 3) {
      // extrapolation to the vertex and then back to the matching plane
      auto mchTrackAtVertex = VarManager::PropagateMuon(mchTrack, collision, VarManager::kToVertex);
      return PropagateToZMCH(mchTrackAtVertex, z);
    }

    if (method == 4) {
      // extrapolation to the MFT DCA and then back to the matching plane
      auto mftTrackDCA = PropagateToZMFT(FwdToTrackPar(mftTrack, mftTrackCov), collision.posZ());
      auto mchTrackAtDCA = PropagateToVertexMCH(FwdToTrackPar(mchTrack, mchTrack),
                                                mftTrackDCA.getX(), mftTrackDCA.getY(), mftTrackDCA.getZ(),
                                                mftTrackDCA.getSigma2X(), mftTrackDCA.getSigma2Y());
      return PropagateToZMCH(mchTrackAtDCA, z);
    }

    return {};
  }

  template <typename TMCH, typename TMFT, typename CMFT, typename C>
  o2::dataformats::GlobalFwdTrack PropagateToMatchingPlaneMFT(const TMCH& mchTrack, const TMFT& mftTrack, const CMFT& mftTrackCov, const C& collision, const double z, int method)
  {
    if (method == 0) {
      // extrapolation with MFT tools
      return PropagateToZMFT(mftTrack, mftTrackCov, z);
    }

    if (method > 0) {
      // extrapolation with MCH tools
      auto mchTrackAtVertex = VarManager::PropagateMuon(mchTrack, collision, VarManager::kToVertex);
      double pMCH = mchTrackAtVertex.getP();
      double px = pMCH * sin(M_PI / 2 - atan(mftTrack.tgl())) * cos(mftTrack.phi());
      double py = pMCH * sin(M_PI / 2 - atan(mftTrack.tgl())) * sin(mftTrack.phi());
      double pt = std::sqrt(std::pow(px, 2) + std::pow(py, 2));
      double sign = mchTrack.sign();

      o2::dataformats::GlobalFwdTrack track = FwdToTrackPar(mftTrack, mftTrackCov);

      // update momentum in track parameters and errors
      auto newCov = track.getCovariances();
      newCov(4, 4) = mchTrackAtVertex.getSigma2InvQPt();
      track.setCovariances(newCov);
      track.setInvQPt(sign / pt);

      auto trackExt = mExtrap.FwdtoMCH(track);
      o2::mch::TrackExtrap::extrapToZCov(trackExt, z);

      o2::dataformats::GlobalFwdTrack propmuon;
      auto proptrack = mExtrap.MCHtoFwd(trackExt);
      propmuon.setParameters(proptrack.getParameters());
      propmuon.setZ(proptrack.getZ());
      propmuon.setCovariances(proptrack.getCovariances());

      return propmuon;
    }

    return {};
  }

  o2::dataformats::GlobalFwdTrack PropagateToVertexMCH(const o2::dataformats::GlobalFwdTrack& muon,
                                                       const double vx, const double vy, const double vz,
                                                       const double covVx, const double covVy)
  {
    auto mchTrack = mExtrap.FwdtoMCH(muon);

    o2::mch::TrackExtrap::extrapToVertex(mchTrack, vx, vy, vz, covVx, covVy);

    auto proptrack = mExtrap.MCHtoFwd(mchTrack);
    o2::dataformats::GlobalFwdTrack propmuon;
    propmuon.setParameters(proptrack.getParameters());
    propmuon.setZ(proptrack.getZ());
    propmuon.setCovariances(proptrack.getCovariances());

    return propmuon;
  }

  template <class TMCH, class C>
  o2::dataformats::GlobalFwdTrack PropagateToVertexMCH(const TMCH& muon,
                                                       const C& collision)
  {
    return PropagateToVertexMCH(FwdToTrackPar(muon, muon),
                                collision.posX(),
                                collision.posY(),
                                collision.posZ(),
                                collision.covXX(),
                                collision.covYY());
  }

  o2::dataformats::GlobalFwdTrack PropagateToVertexMFT(o2::dataformats::GlobalFwdTrack muon,
                                                       const float vx, const float vy, const float vz,
                                                       const float covVx, const float covVy)
  {
    o2::dataformats::GlobalFwdTrack propmuon;
    auto geoMan = o2::base::GeometryManager::meanMaterialBudget(muon.getX(), muon.getY(), muon.getZ(), vx, vy, vz);
    auto x2x0 = static_cast<float>(geoMan.meanX2X0);
    muon.propagateToVtxhelixWithMCS(vz, {vx, vy}, {covVx, covVy}, mBzAtMftCenter, x2x0);
    propmuon.setParameters(muon.getParameters());
    propmuon.setZ(muon.getZ());
    propmuon.setCovariances(muon.getCovariances());

    return propmuon;
  }

  template <class TMFT, class C>
  o2::dataformats::GlobalFwdTrack PropagateToVertexMFT(const TMFT& muon,
                                                       const C& collision)
  {
    return PropagateToVertexMFT(FwdToTrackPar(muon),
                                collision.posX(),
                                collision.posY(),
                                collision.posZ(),
                                collision.covXX(),
                                collision.covYY());
  }

  template <typename TMCH, typename TMFT, class C>
  o2::dataformats::GlobalFwdTrack PropagateToVertexMFT(const TMFT& muon,
                                                       const TMCH& mchTrack,
                                                       const C& collision)
  {
    // extrapolation with MCH tools
    auto mchTrackAtMFT = mExtrap.FwdtoMCH(FwdToTrackPar(mchTrack));
    o2::mch::TrackExtrap::extrapToVertexWithoutBranson(mchTrackAtMFT, muon.z());

    auto muonTrackProp = mExtrap.FwdtoMCH(FwdToTrackPar(muon));

    // update global track momentum from the MCH track
    double pRatio = muonTrackProp.p() / mchTrackAtMFT.p();
    double newInvBendMom = muonTrackProp.getInverseBendingMomentum() * pRatio;
    muonTrackProp.setInverseBendingMomentum(newInvBendMom);
    muonTrackProp.setCharge(mchTrackAtMFT.getCharge());

    o2::mch::TrackExtrap::extrapToVertex(muonTrackProp,
                                         collision.posX(),
                                         collision.posY(),
                                         collision.posZ(),
                                         collision.covXX(),
                                         collision.covYY());

    return mExtrap.MCHtoFwd(muonTrackProp);
  }

  template <class MCP>
  void GetMotherParticles(MCP const& mcParticle, std::vector<std::pair<int64_t, int64_t>>& motherParticlesVec)
  {
    const auto& motherParticles = mcParticle.template mothers_as<aod::McParticles>();
    if (motherParticles.empty()) {
      return;
    }

    const auto& motherParticle = motherParticles[0];
    motherParticlesVec.emplace_back(std::make_pair(static_cast<int64_t>(motherParticle.pdgCode()), static_cast<int64_t>(motherParticle.globalIndex())));
    GetMotherParticles(motherParticle, motherParticlesVec);
  }

  template <class T>
  std::vector<std::pair<int64_t, int64_t>> GetMotherParticles(T const& track)
  {
    std::vector<std::pair<int64_t, int64_t>> result;
    if (!track.has_mcParticle())
      return result;

    const auto& mcParticle = track.mcParticle();
    result.emplace_back(std::make_pair(static_cast<int64_t>(mcParticle.pdgCode()), static_cast<int64_t>(mcParticle.globalIndex())));

    GetMotherParticles(mcParticle, result);

    return result;
  }

  template <class TMCH, class TMFTs>
  int GetDecayRanking(TMCH const& mchTrack, TMFTs const& mftTracks)
  {
    auto mchMotherParticles = GetMotherParticles(mchTrack);

    int decayRanking = 0;
    // search for an MFT track that is associated to one of the MCH mother particles
    for (const auto& mftTrack : mftTracks) {
      // skip tracks that do not have an associated MC particle
      if (!mftTrack.has_mcParticle())
        continue;
      // get the index associated to the MC particle
      auto mftMcParticle = mftTrack.mcParticle();
      int64_t mftMcTrackIndex = mftMcParticle.globalIndex();

      int ranking = 1;
      for (const auto& mother : mchMotherParticles) {
        if (mother.second == mftMcTrackIndex) {
          decayRanking = ranking;
          break;
        }
        ranking += 1;
      }

      if (decayRanking > 0) {
        break;
      }
    }
    return decayRanking;
  }

  template <class TMUON, class TMFT>
  void FillMatchablePairs(CollisionInfo& collisionInfo,
                          TMUON const& muonTracks,
                          TMFT const& mftTracks)
  {
    collisionInfo.matchablePairs.clear();
    for (const auto& muonTrack : muonTracks) {
      // only consider MCH standalone or MCH-MID matches
      if (static_cast<int>(muonTrack.trackType()) <= 2)
        continue;

      // only consider tracks associated to the current collision
      if (!muonTrack.has_collision())
        continue;
      if (muonTrack.collisionId() != collisionInfo.index)
        continue;

      // skip tracks that do not have an associated MC particle
      if (!muonTrack.has_mcParticle())
        continue;
      // get the index associated to the MC particle
      auto muonMcParticle = muonTrack.mcParticle();
      if (std::abs(muonMcParticle.pdgCode()) != 13)
        continue;

      int64_t muonMcTrackIndex = muonMcParticle.globalIndex();

      for (const auto& mftTrack : mftTracks) {
        // skip tracks that do not have an associated MC particle
        if (!mftTrack.has_mcParticle())
          continue;
        // get the index associated to the MC particle
        auto mftMcParticle = mftTrack.mcParticle();
        int64_t mftMcTrackIndex = mftMcParticle.globalIndex();

        if (muonMcTrackIndex == mftMcTrackIndex) {
          collisionInfo.matchablePairs.emplace_back(std::make_pair(static_cast<int64_t>(muonTrack.globalIndex()),
                                                                   static_cast<int64_t>(mftTrack.globalIndex())));
        }
      }
    }
  }

  template <class TMUON>
  int GetTrueMatchIndex(TMUON const& muonTracks,
                        const std::vector<MatchingCandidate>& matchCandidatesVector,
                        const std::vector<std::pair<int64_t, int64_t>>& matchablePairs)
  {
    // find the index of the matching candidate that corresponds to the true match
    // index=1 corresponds to the leading candidate
    // index=0 means no candidate was found that corresponds to the true match
    int trueMatchIndex = 0;
    for (size_t i = 0; i < matchCandidatesVector.size(); i++) {
      auto const& muonTrack = muonTracks.rawIteratorAt(matchCandidatesVector[i].globalTrackId);

      if (IsTrueGlobalMatching(muonTrack, matchablePairs)) {
        trueMatchIndex = i + 1;
        break;
      }
    }
    return trueMatchIndex;
  }

  template <class TMCH, class TMFT>
  bool IsMuon(const TMCH& mchTrack,
              const TMFT& mftTrack)
  {
    // skip tracks that do not have an associated MC particle
    if (!mchTrack.has_mcParticle())
      return false;
    if (!mftTrack.has_mcParticle())
      return false;

    // get the index associated to the MC particles
    auto mchMcParticle = mchTrack.mcParticle();
    auto mftMcParticle = mftTrack.mcParticle();
    if (mchMcParticle.globalIndex() != mftMcParticle.globalIndex())
      return false;

    if (std::abs(mchMcParticle.pdgCode()) != 13)
      return false;

    return true;
  }

  template <class TMUON, class TMUONS, class TMFTS>
  bool IsMuon(const TMUON& muonTrack,
              TMUONS const& /*muonTracks*/,
              TMFTS const& /*mftTracks*/)
  {
    if (static_cast<int>(muonTrack.trackType()) >= 2)
      return false;

    auto const& mchTrack = muonTrack.template matchMCHTrack_as<TMUONS>();
    auto const& mftTrack = muonTrack.template matchMFTTrack_as<TMFTS>();

    return IsMuon(mchTrack, mftTrack);
  }

  template <class TMUON, class TMUONS, class TMFTS>
  MuonMatchType GetMatchType(const TMUON& muonTrack,
                             TMUONS const& muonTracks,
                             TMFTS const& mftTracks,
                             const std::vector<std::pair<int64_t, int64_t>>& matchablePairs,
                             int ranking)
  {
    if (static_cast<int>(muonTrack.trackType()) > 2)
      return kMatchTypeUndefined;

    auto const& mchTrack = muonTrack.template matchMCHTrack_as<TMUONS>();

    bool isPaired = IsMatchableMCH(mchTrack.globalIndex(), matchablePairs);
    bool isMuon = IsMuon(muonTrack, muonTracks, mftTracks);
    int decayRanking = GetDecayRanking(mchTrack, mftTracks);

    MuonMatchType result{kMatchTypeUndefined};
    if (isPaired) {
      if (isMuon) {
        result = (ranking == 1) ? kMatchTypeTrueLeading : kMatchTypeTrueNonLeading;
      } else {
        result = (ranking == 1) ? kMatchTypeWrongLeading : kMatchTypeWrongNonLeading;
      }
    } else if (decayRanking == 2) {
      result = (ranking == 1) ? kMatchTypeDecayLeading : kMatchTypeDecayNonLeading;
    } else {
      result = (ranking == 1) ? kMatchTypeFakeLeading : kMatchTypeFakeNonLeading;
    }

    if (result == kMatchTypeUndefined) {
      std::cout << std::format("[GetMatchType] isPaired={} isMuon={} decayRanking={} result={}",
                               isPaired, isMuon, decayRanking, static_cast<int>(result))
                << std::endl;
    }

    return result;
  }

  // for each MCH standalone track, collect the associated matching candidates
  template <class TMUON, class C>
  void GetSelectedMuons(const CollisionInfo& collisionInfo,
                        C const& collisions,
                        TMUON const& muonTracks,
                        std::vector<int64_t>& selectedMuons)
  {
    selectedMuons.clear();
    for (auto muonTrack : muonTracks) {

      // only consider MCH-MID matches
      if (static_cast<int>(muonTrack.trackType()) != 3) {
        continue;
      }

      // only select MCH-MID tracks from the current collision
      if (!muonTrack.has_collision())
        continue;
      if (static_cast<int64_t>(muonTrack.collisionId()) != collisionInfo.index)
        continue;

      const auto& collision = collisions.rawIteratorAt(muonTrack.collisionId());

      // select MCH tracks with strict quality cuts
      if (!IsGoodMuon(muonTrack, collision,
                      fMuonTaggingTrackChi2MchUp,
                      fMuonTaggingPMchLow,
                      fMuonTaggingPtMchLow,
                      {fMuonTaggingEtaMchLow, fMuonTaggingEtaMchUp},
                      {fMuonTaggingRabsLow, fMuonTaggingRabsUp},
                      fMuonTaggingSigmaPdcaUp)) {
        continue;
      }

      // propagate MCH track to the vertex
      auto mchTrackAtVertex = VarManager::PropagateMuon(muonTrack, collision, VarManager::kToVertex);

      // propagate the track from the vertex to the first MFT plane
      const auto& extrapToMFTfirst = PropagateToZMCH(mchTrackAtVertex, o2::mft::constants::mft::LayerZCoordinate()[0]);
      double rFront = std::sqrt(extrapToMFTfirst.getX() * extrapToMFTfirst.getX() + extrapToMFTfirst.getY() * extrapToMFTfirst.getY());
      double rMinFront = 3.f;
      double rMaxFront = 9.f;
      if (rFront < rMinFront || rFront > rMaxFront)
        continue;

      // propagate the track from the vertex to the last MFT plane
      const auto& extrapToMFTlast = PropagateToZMCH(mchTrackAtVertex, o2::mft::constants::mft::LayerZCoordinate()[9]);
      double rBack = std::sqrt(extrapToMFTlast.getX() * extrapToMFTlast.getX() + extrapToMFTlast.getY() * extrapToMFTlast.getY());
      double rMinBack = 5.f;
      double rMaxBack = 12.f;
      if (rBack < rMinBack || rBack > rMaxBack)
        continue;

      int64_t muonTrackIndex = muonTrack.globalIndex();
      selectedMuons.emplace_back(muonTrackIndex);
    }
  }

  // for each MCH standalone track, collect the associated matching candidates
  template <class TMUON>
  void GetTaggedMuons(const CollisionInfo& collisionInfo,
                      TMUON const& muonTracks,
                      const std::vector<int64_t>& selectedMuons,
                      std::vector<int64_t>& taggedMuons)
  {
    taggedMuons.clear();
    for (auto [mchIndex, globalTracksVector] : collisionInfo.matchingCandidates) {

      // check if the current muon is selected
      if (std::find(selectedMuons.begin(), selectedMuons.end(), mchIndex) == selectedMuons.end())
        continue;

      // if there is only one candidate, mark the muon as select
      if (globalTracksVector.size() == 1) {
        taggedMuons.emplace_back(mchIndex);
        continue;
      }

      auto const& muonTrack0 = muonTracks.rawIteratorAt(globalTracksVector[0].globalTrackId);
      auto const& muonTrack1 = muonTracks.rawIteratorAt(globalTracksVector[1].globalTrackId);

      double chi2diff = muonTrack1.chi2MatchMCHMFT() - muonTrack0.chi2MatchMCHMFT();
      if (chi2diff < fMuonTaggingChi2DiffLow)
        continue;

      taggedMuons.emplace_back(mchIndex);
    }
  }

  void GetMuonPairs(const CollisionInfo& collisionInfo,
                    std::vector<MuonPair>& muonPairs,
                    std::vector<GlobalMuonPair>& globalMuonPairs)
  {
    // outer loop over muon tracks
    for (auto mchIndex1 : collisionInfo.mchTracks) {

      // inner loop over muon tracks
      for (auto mchIndex2 : collisionInfo.mchTracks) {
        // avoid double-counting of muon pairs
        if (mchIndex2 <= mchIndex1)
          continue;

        MuonPair muonPair{{collisionInfo.index, mchIndex1}, {collisionInfo.index, mchIndex2}};
        muonPairs.emplace_back(muonPair);
      }
    }

    // outer loop over global muon tracks
    for (auto& [mchIndex1, matchingCandidates1] : collisionInfo.matchingCandidates) {

      // inner loop over global muon tracks
      for (auto& [mchIndex2, matchingCandidates2] : collisionInfo.matchingCandidates) {
        // avoid double-counting of muon pairs
        if (mchIndex2 <= mchIndex1)
          continue;

        GlobalMuonPair muonPair{{collisionInfo.index, matchingCandidates1}, {collisionInfo.index, matchingCandidates2}};
        globalMuonPairs.emplace_back(muonPair);
      }
    }
  }

  double GetMuMuInvariantMass(const o2::mch::TrackParam& track1, const o2::mch::TrackParam& track2)
  {
    ROOT::Math::PxPyPzMVector muon1{
      track1.px(),
      track1.py(),
      track1.pz(),
      o2::constants::physics::MassMuon};

    ROOT::Math::PxPyPzMVector muon2{
      track2.px(),
      track2.py(),
      track2.pz(),
      o2::constants::physics::MassMuon};

    auto dimuon = muon1 + muon2;

    return dimuon.M();
  }

  double GetMuMuInvariantMass(const o2::dataformats::GlobalFwdTrack& track1, const o2::dataformats::GlobalFwdTrack& track2)
  {
    ROOT::Math::PxPyPzMVector muon1{
      track1.getPx(),
      track1.getPy(),
      track1.getPz(),
      o2::constants::physics::MassMuon};

    ROOT::Math::PxPyPzMVector muon2{
      track2.getPx(),
      track2.getPy(),
      track2.getPz(),
      o2::constants::physics::MassMuon};

    auto dimuon = muon1 + muon2;

    return dimuon.M();
  }

  template <class EVT, class BC, class TMUON, class TMFT>
  void FillCollisions(EVT const& collisions,
                      BC const& bcs,
                      TMUON const& muonTracks,
                      TMFT const& mftTracks,
                      CollisionInfos& collisionInfos)
  {
    collisionInfos.clear();

    std::vector<int64_t> collisionIds;
    for (const auto& collision : collisions) {
      collisionIds.push_back(collision.globalIndex());
    }

    for (size_t cid = 1; cid < collisionIds.size() - 1; cid++) {
      const auto& collision = collisions.rawIteratorAt(collisionIds[cid]);
      int64_t collisionIndex = collision.globalIndex();
      auto bc = bcs.rawIteratorAt(collision.bcId());

      /*// require a minimum BC gap between the current conllision and the previous/next ones
      const auto& collisionPrev = collisions.rawIteratorAt(collisionIds[cid-1]);
      const auto& collisionNext = collisions.rawIteratorAt(collisionIds[cid+1]);
      auto bcPrev = bcs.rawIteratorAt(collisionPrev.bcId());
      auto bcNext = bcs.rawIteratorAt(collisionNext.bcId());
      int64_t deltaPrev = bc.globalBC() - bcPrev.globalBC();
      int64_t deltaNext = bcNext.globalBC() - bc.globalBC();
      if (deltaPrev < 50 || deltaNext < 50) {
        continue;
      }*/

      // fill collision information for global muon tracks (MFT-MCH-MID matches)
      for (auto muonTrack : muonTracks) {
        if (!muonTrack.has_collision())
          continue;

        if (collisionIndex != muonTrack.collisionId()) {
          continue;
        }

        auto& collisionInfo = collisionInfos[collisionIndex];
        collisionInfo.index = collisionIndex;
        collisionInfo.bc = bc.globalBC();
        collisionInfo.zVertex = collision.posZ();

        if (collisionInfo.matchablePairs.empty()) {
          FillMatchablePairs(collisionInfo, muonTracks, mftTracks);
        }

        if (static_cast<int>(muonTrack.trackType()) > 2) {
          // standalone MCH or MCH-MID tracks
          int64_t mchTrackIndex = muonTrack.globalIndex();
          collisionInfo.mchTracks.push_back(mchTrackIndex);
        } else {
          // global muon tracks (MFT-MCH or MFT-MCH-MID)
          int64_t muonTrackIndex = muonTrack.globalIndex();
          double matchChi2 = muonTrack.chi2MatchMCHMFT() / 5.f;
          double matchScore = chi2ToScore(muonTrack.chi2MatchMCHMFT(), 5, 50.f);
          auto const& mchTrack = muonTrack.template matchMCHTrack_as<TMUON>();
          int64_t mchTrackIndex = mchTrack.globalIndex();
          auto const& mftTrack = muonTrack.template matchMFTTrack_as<TMFT>();
          int64_t mftTrackIndex = mftTrack.globalIndex();

          // check if a vector of global muon candidates is already available for the current MCH index
          // if not, initialize a new one and add the current global muon track
          // bool globalMuonTrackFound = false;
          auto matchingCandidateIterator = collisionInfo.matchingCandidates.find(mchTrackIndex);
          if (matchingCandidateIterator != collisionInfo.matchingCandidates.end()) {
            matchingCandidateIterator->second.emplace_back(MatchingCandidate{
              collisionIndex,
              muonTrackIndex,
              mchTrackIndex,
              mftTrackIndex,
              matchScore,
              matchChi2,
              -1,
              matchScore,
              matchChi2,
              -1,
              kMatchTypeUndefined});
          } else {
            collisionInfo.matchingCandidates[mchTrackIndex].emplace_back(MatchingCandidate{
              collisionIndex,
              muonTrackIndex,
              mchTrackIndex,
              mftTrackIndex,
              matchScore,
              matchChi2,
              -1,
              matchScore,
              matchChi2,
              -1,
              kMatchTypeUndefined});
          }
        }
      }

      // fill collision information for MFT standalone tracks
      for (auto mftTrack : mftTracks) {
        if (!mftTrack.has_collision())
          continue;

        if (collisionIndex != mftTrack.collisionId()) {
          continue;
        }

        int64_t mftTrackIndex = mftTrack.globalIndex();

        auto& collisionInfo = collisionInfos[collisionIndex];
        collisionInfo.index = collisionIndex;
        collisionInfo.bc = bc.globalBC();
        collisionInfo.zVertex = collision.posZ();

        collisionInfo.mftTracks.push_back(mftTrackIndex);
      }
    }

    // sort the vectors of matching candidates in ascending order based on the matching score value
    auto compareMatchingScore = [](const MatchingCandidate& track1, const MatchingCandidate& track2) -> bool {
      return (track1.matchScore > track2.matchScore);
    };

    for (auto& [collisionIndex, collisionInfo] : collisionInfos) {
      for (auto& [mchIndex, globalTracksVector] : collisionInfo.matchingCandidates) {
        std::sort(globalTracksVector.begin(), globalTracksVector.end(), compareMatchingScore);

        int ranking = 1;
        for (auto& candidate : globalTracksVector) {
          const auto& muonTrack = muonTracks.rawIteratorAt(candidate.globalTrackId);

          candidate.matchRanking = ranking;
          candidate.matchRankingProd = ranking;
          candidate.matchType = GetMatchType(muonTrack, muonTracks, mftTracks, collisionInfo.matchablePairs, ranking);
          ranking += 1;
        }
      }
    }
  }

  template <class C, class TMUON, class TMFT>
  void FillMatchingPlotsMC(C const& collision,
                           const CollisionInfo& collisionInfo,
                           TMUON const& muonTracks,
                           TMFT const& mftTracks,
                           const MatchingCandidates& matchingCandidates,
                           const MatchingCandidates& matchingCandidatesProd,
                           const std::vector<std::pair<int64_t, int64_t>>& matchablePairs,
                           double matchingScoreCut,
                           MatchingPlotter* plotter,
                           bool verbose = false)
  {
    int mftTrackMult = collisionInfo.mftTracks.size();

    // ====================================
    // Matching candidates hierarchy

    for (const auto& [mchIndex, globalTracksVector] : matchingCandidates) {
      // check if the MCH track belongs to a matchable pair
      bool isPairedMCH = IsMatchableMCH(static_cast<int64_t>(mchIndex), matchablePairs);

      // get the standalone MCH track
      auto const& mchTrack = muonTracks.rawIteratorAt(mchIndex);
      double mchMom = mchTrack.p();
      double mchPt = mchTrack.pt();

      // MCH track quality flag
      bool isGoodMCH = IsGoodGlobalMuon(mchTrack, collision);

      auto matchablePair = GetMatchablePairForMCH(static_cast<int64_t>(mchIndex), matchablePairs);
      bool hasMatchablePair = matchablePair.has_value();
      int decayRanking = 0;
      int mftTrackType = -1;
      float dchi2 = (globalTracksVector.size() >= 2) ? globalTracksVector[1].matchChi2 - globalTracksVector[0].matchChi2 : -1;
      if (hasMatchablePair) {
        auto const& pairedMftTrack = mftTracks.rawIteratorAt(matchablePair.value().second);
        mftTrackType = pairedMftTrack.isCA() ? 1 : 0;
        decayRanking = GetDecayRanking(mchTrack, mftTracks);
      }

      // find the index of the matching candidate that corresponds to the true match
      // index=1 corresponds to the leading candidate
      // index=0 means no candidate was found that corresponds to the true match
      int trueMatchIndex = GetTrueMatchIndex(muonTracks, globalTracksVector, matchablePairs);
      int trueMatchIndexProd = GetTrueMatchIndex(muonTracks, matchingCandidatesProd.at(mchIndex), matchablePairs);

      float mcParticleDz = -1000;
      if (mchTrack.has_mcParticle()) {
        const auto& mchMcParticle = mchTrack.mcParticle();
        mcParticleDz = collision.posZ() - mchMcParticle.vz();
      }

      std::get<std::shared_ptr<TH1>>(plotter->fMatchRanking->hist)->Fill(trueMatchIndex);
      std::get<std::shared_ptr<TH2>>(plotter->fMatchRanking->histVsP)->Fill(mchMom, trueMatchIndex);
      std::get<std::shared_ptr<TH2>>(plotter->fMatchRanking->histVsPt)->Fill(mchPt, trueMatchIndex);
      std::get<std::shared_ptr<TH2>>(plotter->fMatchRanking->histVsMcParticleDz)->Fill(mcParticleDz, trueMatchIndex);
      std::get<std::shared_ptr<TH2>>(plotter->fMatchRanking->histVsMftTrackMult)->Fill(mftTrackMult, trueMatchIndex);
      std::get<std::shared_ptr<TH2>>(plotter->fMatchRanking->histVsMftTrackType)->Fill(mftTrackType, trueMatchIndex);
      std::get<std::shared_ptr<TH2>>(plotter->fMatchRanking->histVsProdRanking)->Fill(trueMatchIndexProd, trueMatchIndex);
      if (dchi2 >= 0)
        std::get<std::shared_ptr<TH2>>(plotter->fMatchRanking->histVsDeltaChi2)->Fill(dchi2, trueMatchIndex);

      if (isGoodMCH) {
        std::get<std::shared_ptr<TH1>>(plotter->fMatchRankingGoodMCH->hist)->Fill(trueMatchIndex);
        std::get<std::shared_ptr<TH2>>(plotter->fMatchRankingGoodMCH->histVsP)->Fill(mchMom, trueMatchIndex);
        std::get<std::shared_ptr<TH2>>(plotter->fMatchRankingGoodMCH->histVsPt)->Fill(mchPt, trueMatchIndex);
        std::get<std::shared_ptr<TH2>>(plotter->fMatchRankingGoodMCH->histVsMcParticleDz)->Fill(mcParticleDz, trueMatchIndex);
        std::get<std::shared_ptr<TH2>>(plotter->fMatchRankingGoodMCH->histVsMftTrackMult)->Fill(mftTrackMult, trueMatchIndex);
        std::get<std::shared_ptr<TH2>>(plotter->fMatchRankingGoodMCH->histVsMftTrackType)->Fill(mftTrackType, trueMatchIndex);
        std::get<std::shared_ptr<TH2>>(plotter->fMatchRankingGoodMCH->histVsProdRanking)->Fill(trueMatchIndexProd, trueMatchIndex);
        if (dchi2 >= 0)
          std::get<std::shared_ptr<TH2>>(plotter->fMatchRankingGoodMCH->histVsDeltaChi2)->Fill(dchi2, trueMatchIndex);
      }

      if (isPairedMCH) {
        std::get<std::shared_ptr<TH1>>(plotter->fMatchRankingPaired->hist)->Fill(trueMatchIndex);
        std::get<std::shared_ptr<TH2>>(plotter->fMatchRankingPaired->histVsP)->Fill(mchMom, trueMatchIndex);
        std::get<std::shared_ptr<TH2>>(plotter->fMatchRankingPaired->histVsPt)->Fill(mchPt, trueMatchIndex);
        std::get<std::shared_ptr<TH2>>(plotter->fMatchRankingPaired->histVsMcParticleDz)->Fill(mcParticleDz, trueMatchIndex);
        std::get<std::shared_ptr<TH2>>(plotter->fMatchRankingPaired->histVsMftTrackMult)->Fill(mftTrackMult, trueMatchIndex);
        std::get<std::shared_ptr<TH2>>(plotter->fMatchRankingPaired->histVsMftTrackType)->Fill(mftTrackType, trueMatchIndex);
        std::get<std::shared_ptr<TH2>>(plotter->fMatchRankingPaired->histVsProdRanking)->Fill(trueMatchIndexProd, trueMatchIndex);
        if (dchi2 >= 0)
          std::get<std::shared_ptr<TH2>>(plotter->fMatchRankingPaired->histVsDeltaChi2)->Fill(dchi2, trueMatchIndex);
      }

      if (isGoodMCH && isPairedMCH) {
        std::get<std::shared_ptr<TH1>>(plotter->fMatchRankingPairedGoodMCH->hist)->Fill(trueMatchIndex);
        std::get<std::shared_ptr<TH2>>(plotter->fMatchRankingPairedGoodMCH->histVsP)->Fill(mchMom, trueMatchIndex);
        std::get<std::shared_ptr<TH2>>(plotter->fMatchRankingPairedGoodMCH->histVsPt)->Fill(mchPt, trueMatchIndex);
        std::get<std::shared_ptr<TH2>>(plotter->fMatchRankingPairedGoodMCH->histVsMcParticleDz)->Fill(mcParticleDz, trueMatchIndex);
        std::get<std::shared_ptr<TH2>>(plotter->fMatchRankingPairedGoodMCH->histVsMftTrackMult)->Fill(mftTrackMult, trueMatchIndex);
        std::get<std::shared_ptr<TH2>>(plotter->fMatchRankingPairedGoodMCH->histVsMftTrackType)->Fill(mftTrackType, trueMatchIndex);
        std::get<std::shared_ptr<TH2>>(plotter->fMatchRankingPairedGoodMCH->histVsProdRanking)->Fill(trueMatchIndexProd, trueMatchIndex);
        if (dchi2 >= 0)
          std::get<std::shared_ptr<TH2>>(plotter->fMatchRankingPairedGoodMCH->histVsDeltaChi2)->Fill(dchi2, trueMatchIndex);
      }

      if (trueMatchIndex == 0) {
        // missed matches
        if (!isPairedMCH) {
          // the MCH track does not have a corresponding MFT track for matching
          if (mchTrack.has_mcParticle()) {
            // the MCH track is not fake
            std::get<std::shared_ptr<TH1>>(plotter->fMissedMatches)->Fill(0);
            if (isGoodMCH) {
              std::get<std::shared_ptr<TH1>>(plotter->fMissedMatchesGoodMCH)->Fill(0);
            }
          } else {
            // the MCH track is fake
            std::get<std::shared_ptr<TH1>>(plotter->fMissedMatches)->Fill(1);
            if (isGoodMCH) {
              std::get<std::shared_ptr<TH1>>(plotter->fMissedMatchesGoodMCH)->Fill(1);
            }
          }
        } else {
          // the correct match is not among the stored candidates
          std::get<std::shared_ptr<TH1>>(plotter->fMissedMatches)->Fill(2);
          if (isGoodMCH) {
            std::get<std::shared_ptr<TH1>>(plotter->fMissedMatchesGoodMCH)->Fill(2);
          }
        }
      }

      if (globalTracksVector.size() > 1 && trueMatchIndex > 0) {
        // we have a matchable pair with at least two candidates, so we can check the score difference
        // between the leading and the sub-leading
        auto leadingScore = globalTracksVector[0].matchScore;
        auto subleadingScore = globalTracksVector[1].matchScore;
        if (trueMatchIndex == 1) {
          std::get<std::shared_ptr<TH1>>(plotter->fScoreGapLeadingTrueMatches)->Fill(leadingScore - subleadingScore);
        } else {
          std::get<std::shared_ptr<TH1>>(plotter->fScoreGapNonLeadingTrueMatches)->Fill(leadingScore - subleadingScore);
        }
      }

      if (trueMatchIndex == 0) {
        // missed matches
        std::get<std::shared_ptr<TH1>>(plotter->fDecayRankingMissedMatches)->Fill(decayRanking);
      } else if (trueMatchIndex == 1) {
        // good matches
        std::get<std::shared_ptr<TH1>>(plotter->fDecayRankingGoodMatches)->Fill(decayRanking);
      } else {
        // non-leading matches
        std::get<std::shared_ptr<TH1>>(plotter->fDecayRankingNonLeadingMatches)->Fill(decayRanking);
      }
    }

    // ====================================
    // Matching properties

    for (auto [mchIndex, globalTracksVector] : matchingCandidates) {
      if (globalTracksVector.size() < 1)
        continue;

      // get the standalone MCH track
      auto const& mchTrack = muonTracks.rawIteratorAt(mchIndex);

      // skip global muon tracks that do not pass the MCH and MFT quality cuts
      if (!IsGoodGlobalMuon(mchTrack, collision))
        continue;

      double mchMom = mchTrack.p();
      double mchPt = mchTrack.pt();

      // matching score analysis
      for (const auto& candidate : globalTracksVector) {
        int matchType = static_cast<int>(candidate.matchType);
        std::get<std::shared_ptr<TH1>>(plotter->fMatchType)->Fill(matchType);
        std::get<std::shared_ptr<TH2>>(plotter->fMatchTypeVsP)->Fill(mchMom, matchType);
        std::get<std::shared_ptr<TH2>>(plotter->fMatchTypeVsPt)->Fill(mchPt, matchType);

        std::get<std::shared_ptr<TH2>>(plotter->fMatchChi2VsType)->Fill(matchType, candidate.matchChi2);
        std::get<std::shared_ptr<TH3>>(plotter->fMatchChi2VsTypeVsP)->Fill(mchMom, matchType, candidate.matchChi2);
        std::get<std::shared_ptr<TH3>>(plotter->fMatchChi2VsTypeVsPt)->Fill(mchPt, matchType, candidate.matchChi2);

        std::get<std::shared_ptr<TH2>>(plotter->fMatchScoreVsType)->Fill(matchType, candidate.matchScore);
        std::get<std::shared_ptr<TH3>>(plotter->fMatchScoreVsTypeVsP)->Fill(mchMom, matchType, candidate.matchScore);
        std::get<std::shared_ptr<TH3>>(plotter->fMatchScoreVsTypeVsPt)->Fill(mchPt, matchType, candidate.matchScore);
      }
    }

    for (auto [mchIndex, globalTracksVector] : matchingCandidates) {
      if (globalTracksVector.size() < 1)
        continue;

      int trueMatchIndex = GetTrueMatchIndex(muonTracks, globalTracksVector, matchablePairs);

      // loop over candidates
      int candidateIndex = 1;
      for (const auto& candidate : globalTracksVector) {
        auto const& muonTrack = muonTracks.rawIteratorAt(candidate.globalTrackId);

        float matchScore = candidate.matchScore;
        float matchChi2 = candidate.matchChi2;

        float matchChi2Prod = muonTrack.chi2MatchMCHMFT() / 5.f;
        float matchScoreProd = chi2ToScore(muonTrack.chi2MatchMCHMFT(), 5, 50.f);

        std::get<std::shared_ptr<TH2>>(plotter->fMatchScoreVsProd)->Fill(matchScoreProd, matchScore);
        std::get<std::shared_ptr<TH2>>(plotter->fMatchChi2VsProd)->Fill(matchChi2Prod, matchChi2);

        if (candidateIndex == trueMatchIndex) {
          std::get<std::shared_ptr<TH2>>(plotter->fTrueMatchScoreVsProd)->Fill(matchScoreProd, matchScore);
          std::get<std::shared_ptr<TH2>>(plotter->fTrueMatchChi2VsProd)->Fill(matchChi2Prod, matchChi2);
        }

        candidateIndex += 1;
      }
    }

    // ====================================
    // Matching purity
    if (verbose)
      std::cout << std::format("  Filling matching purity plots with score cut {}", matchingScoreCut) << std::endl;
    for (auto [mchIndex, globalTracksVector] : matchingCandidates) {
      if (globalTracksVector.size() < 1)
        continue;

      // get the leading matching candidate
      auto const& muonTrack = muonTracks.rawIteratorAt(globalTracksVector[0].globalTrackId);
      double matchingScore = globalTracksVector[0].matchScore;

      // get the standalone MCH and MFT tracks
      auto const& mchTrack = muonTracks.rawIteratorAt(mchIndex);
      auto const& mftTrack = muonTrack.template matchMFTTrack_as<TMFT>();

      // skip global muon tracks that do not pass the MCH and MFT quality cuts
      if (!IsGoodGlobalMuon(mchTrack, collision))
        continue;
      if (!IsGoodMFT(mftTrack))
        continue;

      // skip  candidates that do not pass the matching quality cuts
      if (!IsGoodGlobalMatching(muonTrack, matchingScore, matchingScoreCut))
        continue;

      // check if the matching candidate is a true one
      bool isTrueMatch = IsTrueGlobalMatching(muonTrack, matchablePairs);

      if (verbose)
        std::cout << std::format("    MCH track #{} -> Muon track #{}, isTrueMatch={}", mchIndex, globalTracksVector[0].globalTrackId, isTrueMatch) << std::endl;
      // fill matching purity plots
      plotter->fMatchingPurityPlotter.Fill(mchTrack, isTrueMatch);
    }

    // ====================================
    // Matching efficiencies

    // outer loop on matchable pairs
    for (auto [matchableMchIndex, matchableMftIndex] : matchablePairs) {
      // get the standalone MCH track
      auto const& mchTrack = muonTracks.rawIteratorAt(matchableMchIndex);

      // skip  track pairs that do not pass the MCH and MFT quality cuts
      // we only consider matchable pairs that fulfill the track quality requirements
      if (!IsGoodGlobalMuon(mchTrack, collision))
        continue;

      bool goodMatchFound = false;
      bool isTrueMatch = false;

      // check if we have some matching candidates for the current matchable MCH track
      if (matchingCandidates.count(matchableMchIndex) > 0) {
        const auto& globalTracksVector = matchingCandidates.at(static_cast<int64_t>(matchableMchIndex));
        if (!globalTracksVector.empty()) {
          // get the leading matching candidate
          auto const& muonTrack = muonTracks.rawIteratorAt(globalTracksVector[0].globalTrackId);
          double matchingScore = globalTracksVector[0].matchScore;

          // get the standalone MFT track
          auto const& mftTrack = muonTrack.template matchMFTTrack_as<TMFT>();
          auto mftIndex = mftTrack.globalIndex();

          goodMatchFound = IsGoodGlobalMatching(muonTrack, matchingScore, matchingScoreCut);
          isTrueMatch = (mftIndex == matchableMftIndex);
        }
      }

      // fill matching efficiency plots
      plotter->fPairingEfficiencyPlotter.Fill(mchTrack, goodMatchFound);
      plotter->fMatchingEfficiencyPlotter.Fill(mchTrack, (goodMatchFound && isTrueMatch));
      plotter->fFakeMatchingEfficiencyPlotter.Fill(mchTrack, (goodMatchFound && !isTrueMatch));
    }
  }

  template <class C, class TMUON, class TMFT>
  void FillDimuonPlotsMC(const CollisionInfo& collisionInfo,
                         C const& collisions,
                         TMUON const& muonTracks,
                         TMFT const& /*mftTracks*/)
  {
    std::vector<MuonPair> muonPairs;
    std::vector<GlobalMuonPair> globalMuonPairs;

    GetMuonPairs(collisionInfo, muonPairs, globalMuonPairs);

    for (auto& [muon1, muon2] : muonPairs) {
      auto const& collision = collisions.rawIteratorAt(muon1.first);

      auto mchIndex1 = muon1.second;
      auto mchIndex2 = muon2.second;
      auto const& muonTrack1 = muonTracks.rawIteratorAt(mchIndex1);
      auto const& muonTrack2 = muonTracks.rawIteratorAt(mchIndex2);
      int sign1 = muonTrack1.sign();
      int sign2 = muonTrack2.sign();

      // only consider opposite-sign pairs
      if ((sign1 * sign2) >= 0)
        continue;

      bool goodMuonTracks = (IsGoodMuon(muonTrack1, collision) && IsGoodMuon(muonTrack2, collision));

      if (goodMuonTracks) {
        double mass = GetMuMuInvariantMass(PropagateToVertexMCH(muonTrack1, collision),
                                           PropagateToVertexMCH(muonTrack2, collision));
        registryDimuon.get<TH1>(HIST("dimuon/invariantMass_MuonKine_MuonCuts"))->Fill(mass);
      }
    }

    for (auto& [muon1, muon2] : globalMuonPairs) {
      auto& candidates1 = muon1.second;
      auto& candidates2 = muon2.second;

      auto const& collision = collisions.rawIteratorAt(muon1.first);

      auto const& muonTrack1 = muonTracks.rawIteratorAt(candidates1[0].globalTrackId);
      auto const& muonTrack2 = muonTracks.rawIteratorAt(candidates2[0].globalTrackId);
      auto matchScore1 = candidates1[0].matchScore;
      auto matchScore2 = candidates2[0].matchScore;
      auto const& mchTrack1 = muonTrack1.template matchMCHTrack_as<TMUON>();
      auto const& mchTrack2 = muonTrack2.template matchMCHTrack_as<TMUON>();
      int sign1 = mchTrack1.sign();
      int sign2 = mchTrack2.sign();

      // only consider opposite-sign pairs
      if ((sign1 * sign2) >= 0)
        continue;

      double p1 = mchTrack1.p();
      double p2 = mchTrack2.p();
      int matchType = -1;
      if (p1 >= p2) {
        matchType = candidates1[0].matchType * 10 + candidates2[0].matchType;
      } else {
        matchType = candidates2[0].matchType * 10 + candidates1[0].matchType;
      }

      bool goodGlobalMuonTracks = (IsGoodGlobalMuon(mchTrack1, collision) && IsGoodGlobalMuon(mchTrack2, collision));
      if (!goodGlobalMuonTracks) {
        continue;
      }

      bool goodGlobalMuonMatches = (IsGoodGlobalMatching(muonTrack1, matchScore1) && IsGoodGlobalMatching(muonTrack2, matchScore2));

      double massMCH = GetMuMuInvariantMass(PropagateToVertexMCH(mchTrack1, collision),
                                            PropagateToVertexMCH(mchTrack2, collision));
      double mass = GetMuMuInvariantMass(PropagateToVertexMCH(muonTrack1, collision),
                                         PropagateToVertexMCH(muonTrack2, collision));
      registryDimuon.get<TH1>(HIST("dimuon/invariantMass_MuonKine_GlobalMuonCuts"))->Fill(massMCH);
      registryDimuon.get<TH1>(HIST("dimuon/invariantMass_ScaledMftKine_GlobalMuonCuts"))->Fill(mass);
      registryDimuon.get<TH2>(HIST("dimuon/MC/invariantMass_MuonKine_GlobalMuonCuts_vs_match_type"))->Fill(massMCH, matchType);
      registryDimuon.get<TH2>(HIST("dimuon/MC/invariantMass_ScaledMftKine_GlobalMuonCuts_vs_match_type"))->Fill(mass, matchType);

      if (goodGlobalMuonMatches) {
        registryDimuon.get<TH1>(HIST("dimuon/invariantMass_MuonKine_GlobalMuonCuts_GoodMatches"))->Fill(massMCH);
        registryDimuon.get<TH1>(HIST("dimuon/invariantMass_ScaledMftKine_GlobalMuonCuts_GoodMatches"))->Fill(mass);
        registryDimuon.get<TH2>(HIST("dimuon/MC/invariantMass_MuonKine_GlobalMuonCuts_GoodMatches_vs_match_type"))->Fill(massMCH, matchType);
        registryDimuon.get<TH2>(HIST("dimuon/MC/invariantMass_ScaledMftKine_GlobalMuonCuts_GoodMatches_vs_match_type"))->Fill(mass, matchType);
      }
    }
  }

  template <class C, class TMUON, class TMFT, class CMFT>
  void RunChi2Matching(C const& collisions,
                       TMUON const& muonTracks,
                       TMFT const& /*mftTracks*/,
                       CMFT const& mftCovs,
                       std::string funcName,
                       float matchingPlaneZ,
                       int extrapMethod,
                       const MatchingCandidates& matchingCandidates,
                       MatchingCandidates& newMatchingCandidates)
  {
    newMatchingCandidates.clear();

    if (funcName == "prod") {
      newMatchingCandidates = matchingCandidates;
      return;
    }

    if (mMatchingFunctionMap.count(funcName) < 1)
      return;
    auto matchingFunc = mMatchingFunctionMap.at(funcName);

    for (auto& [mchIndex, globalTracksVector] : matchingCandidates) {
      auto const& mchTrack = muonTracks.rawIteratorAt(mchIndex);

      for (const auto& candidate : globalTracksVector) {
        auto const& muonTrack = muonTracks.rawIteratorAt(candidate.globalTrackId);
        if (!muonTrack.has_collision())
          continue;

        auto collision = collisions.rawIteratorAt(muonTrack.collisionId());

        // get MCH and MFT standalone tracks
        // auto mchTrack = muonTrack.template matchMCHTrack_as<TMUON>();
        auto const& mftTrack = muonTrack.template matchMFTTrack_as<TMFT>();
        if (mftTrackCovs.count(mftTrack.globalIndex()) < 1) {
          // std::cout << std::format("Covariance matrix for MFT track #{} not found", mftTrack.globalIndex()) << std::endl;
          continue;
        }
        auto const& mftTrackCov = mftCovs.rawIteratorAt(mftTrackCovs[mftTrack.globalIndex()]);

        // get tracks parameters in O2 format
        auto mftTrackProp = FwdToTrackPar(mftTrack, mftTrackCov);
        auto mchTrackProp = FwdToTrackPar(mchTrack, mchTrack);

        if (matchingPlaneZ < 0.) {
          mftTrackProp = PropagateToMatchingPlaneMFT(mchTrack, mftTrack, mftTrackCov, collision, matchingPlaneZ, extrapMethod);
          mchTrackProp = PropagateToMatchingPlaneMCH(mchTrack, mftTrack, mftTrackCov, collision, matchingPlaneZ, extrapMethod);
        }

        // run the chi2 matching function
        auto matchResult = matchingFunc(mchTrackProp, mftTrackProp);
        float matchChi2 = std::get<0>(matchResult) / std::get<1>(matchResult);
        float matchScore = chi2ToScore(std::get<0>(matchResult), std::get<1>(matchResult), 10.f * std::get<1>(matchResult));
        float matchChi2Prod = muonTrack.chi2MatchMCHMFT() / 5.f;
        float matchScoreProd = chi2ToScore(muonTrack.chi2MatchMCHMFT(), 5, 50.f);

        // check if a vector of global muon candidates is already available for the current MCH index
        // if not, initialize a new one and add the current global muon track
        auto matchingCandidateIterator = newMatchingCandidates.find(mchIndex);
        if (matchingCandidateIterator != newMatchingCandidates.end()) {
          matchingCandidateIterator->second.emplace_back(MatchingCandidate{
            muonTrack.collisionId(),
            candidate.globalTrackId,
            mchIndex,
            mftTrack.globalIndex(),
            matchScore,
            matchChi2,
            -1,
            matchScoreProd,
            matchChi2Prod,
            -1,
            kMatchTypeUndefined});
        } else {
          newMatchingCandidates[mchIndex].emplace_back(MatchingCandidate{
            muonTrack.collisionId(),
            candidate.globalTrackId,
            mchIndex,
            mftTrack.globalIndex(),
            matchScore,
            matchChi2,
            -1,
            matchScoreProd,
            matchChi2Prod,
            -1,
            kMatchTypeUndefined});
        }
      }
    }

    // sort the vectors of matching candidates in ascending order based on the matching score value
    auto compareMatchingScore = [](const MatchingCandidate& track1, const MatchingCandidate& track2) -> bool {
      return (track1.matchScore > track2.matchScore);
    };

    for (auto& [mchIndex, globalTracksVector] : newMatchingCandidates) {
      std::sort(globalTracksVector.begin(), globalTracksVector.end(), compareMatchingScore);
    }
  }

  template <class C, class TMUON, class TMFT, class CMFT>
  void RunChi2Matching(C const& collisions,
                       TMUON const& muonTracks,
                       TMFT const& mftTracks,
                       CMFT const& mftCovs,
                       std::string label,
                       const MatchingCandidates& matchingCandidates,
                       MatchingCandidates& newMatchingCandidates)
  {
    newMatchingCandidates.clear();

    auto funcIter = matchingChi2Functions.find(label);
    if (funcIter == matchingChi2Functions.end())
      return;

    auto funcName = funcIter->second;

    if (funcName == "prod") {
      newMatchingCandidates = matchingCandidates;
      return;
    }

    if (mMatchingFunctionMap.count(funcName) < 1)
      return;

    // extrapolation parameters
    auto matchingPlaneZ = matchingPlanesZ[label];
    auto extrapMethod = matchingExtrapMethod[label];

    RunChi2Matching(collisions, muonTracks, mftTracks, mftCovs, funcName, matchingPlaneZ, extrapMethod, matchingCandidates, newMatchingCandidates);
  }

  template <class C, class TMUON, class TMFT, class CMFT>
  void RunMLMatching(C const& collisions,
                     TMUON const& muonTracks,
                     TMFT const& /*mftTracks*/,
                     CMFT const& mftCovs,
                     std::string label,
                     const MatchingCandidates& matchingCandidates,
                     MatchingCandidates& newMatchingCandidates)
  {
    newMatchingCandidates.clear();
    auto mlIter = matchingMlResponses.find(label);
    if (mlIter == matchingMlResponses.end())
      return;

    auto& mlResponse = mlIter->second;
    for (auto& [mchIndex, globalTracksVector] : matchingCandidates) {
      auto const& mchTrack = muonTracks.rawIteratorAt(mchIndex);
      for (const auto& candidate : globalTracksVector) {
        auto const& muonTrack = muonTracks.rawIteratorAt(candidate.globalTrackId);
        if (!muonTrack.has_collision())
          continue;

        auto collision = collisions.rawIteratorAt(muonTrack.collisionId());

        // get MFT standalone track
        auto const& mftTrack = muonTrack.template matchMFTTrack_as<TMFT>();
        if (mftTrackCovs.count(mftTrack.globalIndex()) < 1) {
          // std::cout << std::format("Covariance matrix for MFT track #{} not found", mftTrack.globalIndex()) << std::endl;
          continue;
        }
        // std::cout << fmt::format("Getting covariance matrix for MFT track #{} -> {}", mftTrack.globalIndex(), mftTrackCovs[mftTrack.globalIndex()]) << std::endl;
        auto const& mftTrackCov = mftCovs.rawIteratorAt(mftTrackCovs[mftTrack.globalIndex()]);
        // std::cout << fmt::format("Covariance matrix for MFT track #{} retrieved", mftTrack.globalIndex()) << std::endl;

        // get tracks parameters in O2 format
        auto mftTrackProp = FwdToTrackPar(mftTrack, mftTrackCov);
        auto mchTrackProp = FwdToTrackPar(mchTrack, mchTrack);

        // extrapolate to the matching plane
        auto matchingPlaneZ = matchingPlanesZ[label];
        if (matchingPlaneZ < 0.) {
          mftTrackProp = PropagateToZMFT(mftTrackProp, matchingPlaneZ);
          mchTrackProp = PropagateToZMCH(mchTrackProp, matchingPlaneZ);
        }

        // run the ML model
        std::vector<float> output;
        std::vector<float> inputML = mlResponse.getInputFeaturesGlob(muonTrack, mchTrackProp, mftTrackProp, collision);
        mlResponse.isSelectedMl(inputML, 0, output);
        float matchScore = output[0];
        float matchChi2Prod = muonTrack.chi2MatchMCHMFT() / 5.f;
        float matchScoreProd = chi2ToScore(muonTrack.chi2MatchMCHMFT(), 5, 50.f);
        // std::cout << std::format("Matching score: {}, Chi2: {}", matchingScore, muonTrack.chi2MatchMCHMFT()) << std::endl;

        // check if a vector of global muon candidates is already available for the current MCH index
        // if not, initialize a new one and add the current global muon track
        auto matchingCandidateIterator = newMatchingCandidates.find(mchIndex);
        if (matchingCandidateIterator != newMatchingCandidates.end()) {
          matchingCandidateIterator->second.emplace_back(MatchingCandidate{
            muonTrack.collisionId(),
            candidate.globalTrackId,
            mchIndex,
            mftTrack.globalIndex(),
            matchScore,
            -1,
            -1,
            matchScoreProd,
            matchChi2Prod,
            -1,
            kMatchTypeUndefined});
        } else {
          newMatchingCandidates[mchIndex].emplace_back(MatchingCandidate{
            muonTrack.collisionId(),
            candidate.globalTrackId,
            mchIndex,
            mftTrack.globalIndex(),
            matchScore,
            -1,
            -1,
            matchScoreProd,
            matchChi2Prod,
            -1,
            kMatchTypeUndefined});
        }
      }
    }

    // sort the vectors of matching candidates in ascending order based on the matching score value
    auto compareMatchingScore = [](const MatchingCandidate& track1, const MatchingCandidate& track2) -> bool {
      return (track1.matchScore > track2.matchScore);
    };

    for (auto& [mchIndex, globalTracksVector] : newMatchingCandidates) {
      std::sort(globalTracksVector.begin(), globalTracksVector.end(), compareMatchingScore);
    }
  }

  template <class C, class TMUON, class TMFT, class CMFT>
  void ProcessCollisionMC(const CollisionInfo& collisionInfo,
                          C const& collisions,
                          TMUON const& muonTracks,
                          TMFT const& mftTracks,
                          CMFT const& mftCovs)
  {
    auto collision = collisions.rawIteratorAt(collisionInfo.index);

    registry.get<TH1>(HIST("tracksMultiplicityMFT"))->Fill(collisionInfo.mftTracks.size());
    registry.get<TH1>(HIST("tracksMultiplicityMCH"))->Fill(collisionInfo.mchTracks.size());

    // Chi2-based matching analysis
    FillMatchingPlotsMC(collision, collisionInfo, muonTracks, mftTracks, collisionInfo.matchingCandidates, collisionInfo.matchingCandidates, collisionInfo.matchablePairs, fMatchingChi2ScoreMftMchLow, fChi2MatchingPlotter.get(), false);
    for (auto& [label, func] : matchingChi2Functions) {
      MatchingCandidates matchingCandidates;
      RunChi2Matching(collisions, muonTracks, mftTracks, mftCovs, label, collisionInfo.matchingCandidates, matchingCandidates);

      auto* plotter = fMatchingPlotters.at(label).get();
      double matchingScoreCut = matchingScoreCuts.at(label);

      FillMatchingPlotsMC(collision, collisionInfo, muonTracks, mftTracks, matchingCandidates, collisionInfo.matchingCandidates, collisionInfo.matchablePairs, matchingScoreCut, plotter, false);
    }

    // ML-based matching analysis
    for (auto& [label, mlResponse] : matchingMlResponses) {
      MatchingCandidates matchingCandidates;
      RunMLMatching(collisions, muonTracks, mftTracks, mftCovs, label, collisionInfo.matchingCandidates, matchingCandidates);

      auto* plotter = fMatchingPlotters.at(label).get();
      double matchingScoreCut = matchingScoreCuts.at(label);

      FillMatchingPlotsMC(collision, collisionInfo, muonTracks, mftTracks, matchingCandidates, collisionInfo.matchingCandidates, collisionInfo.matchablePairs, matchingScoreCut, plotter);
    }

    // Muons tagging
    for (auto [mchIndex, mftIndex] : collisionInfo.matchablePairs) {
      auto const& mchTrack = muonTracks.rawIteratorAt(mchIndex);
      if (!mchTrack.has_collision())
        continue;
      auto collision = collisions.rawIteratorAt(mchTrack.collisionId());

      auto const& mftTrack = mftTracks.rawIteratorAt(mftIndex);
      if (mftTrackCovs.count(mftTrack.globalIndex()) < 1) {
        continue;
      }
      auto const& mftTrackCov = mftCovs.rawIteratorAt(mftTrackCovs[mftTrack.globalIndex()]);

      auto mchTrackAtVertex = VarManager::PropagateMuon(mchTrack, collision, VarManager::kToVertex);

      // extrapolate to the matching plane
      auto z = o2::mft::constants::mft::LayerZCoordinate()[9];
      auto mchTrackProp = PropagateToZMCH(mchTrackAtVertex, z);
      auto mftTrackProp = PropagateToZMFT(FwdToTrackPar(mftTrack, mftTrackCov), z);

      registry.get<TH2>(HIST("matching/MC/pairedMCHTracksAtMFT"))->Fill(mchTrackProp.getX(), mchTrackProp.getY());
      registry.get<TH2>(HIST("matching/MC/pairedMFTTracksAtMFT"))->Fill(mftTrackProp.getX(), mftTrackProp.getY());
    }

    std::vector<int64_t> selectedMuons;
    GetSelectedMuons(collisionInfo, collisions, muonTracks, selectedMuons);

    MatchingCandidates selectedMatchingCandidates;
    for (auto [mchIndex, globalTracksVector] : collisionInfo.matchingCandidates) {
      if (std::find(selectedMuons.begin(), selectedMuons.end(), mchIndex) != selectedMuons.end()) {
        selectedMatchingCandidates[mchIndex] = globalTracksVector;
      }
    }
    FillMatchingPlotsMC(collision, collisionInfo, muonTracks, mftTracks, selectedMatchingCandidates, collisionInfo.matchingCandidates, collisionInfo.matchablePairs, fMatchingChi2ScoreMftMchLow, fSelectedMuonsMatchingPlotter.get());

    std::vector<int64_t> taggedMuons;
    GetTaggedMuons(collisionInfo, muonTracks, selectedMuons, taggedMuons);

    MatchingCandidates taggedMatchingCandidates;
    for (auto [mchIndex, globalTracksVector] : collisionInfo.matchingCandidates) {
      if (std::find(taggedMuons.begin(), taggedMuons.end(), mchIndex) != taggedMuons.end()) {
        taggedMatchingCandidates[mchIndex] = globalTracksVector;
      }
    }
    FillMatchingPlotsMC(collision, collisionInfo, muonTracks, mftTracks, taggedMatchingCandidates, collisionInfo.matchingCandidates, collisionInfo.matchablePairs, fMatchingChi2ScoreMftMchLow, fTaggedMuonsMatchingPlotter.get());

    // Di-muon analysis
    FillDimuonPlotsMC(collisionInfo, collisions, muonTracks, mftTracks);
  }

  void processQAMC(MyEvents const& collisions,
                   aod::BCsWithTimestamps const& bcs,
                   MyMuonsMC const& muonTracks,
                   MyMFTsMC const& mftTracks,
                   MyMFTCovariances const& mftCovs,
                   aod::McParticles const& /*mcParticles*/)
  {
    auto bc = bcs.begin();
    initCCDB(bc);

    for (auto& muon : muonTracks) {
      registry.get<TH1>(HIST("nTracksPerType"))->Fill(static_cast<int>(muon.trackType()));
    }

    FillCollisions(collisions, bcs, muonTracks, mftTracks, fCollisionInfos);

    mftTrackCovs.clear();
    for (auto& mftTrackCov : mftCovs) {
      mftTrackCovs[mftTrackCov.matchMFTTrackId()] = mftTrackCov.globalIndex();
    }

    for (auto const& [collisionIndex, collisionInfo] : fCollisionInfos) {
      ProcessCollisionMC(collisionInfo, collisions, muonTracks, mftTracks, mftCovs);
    }
  }

  PROCESS_SWITCH(qaMatching, processQAMC, "process qa MC", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<qaMatching>(cfgc)};
};
