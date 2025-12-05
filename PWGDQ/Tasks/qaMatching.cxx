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

#include <algorithm>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::aod;

using MyEvents = soa::Join<aod::Collisions, aod::EvSels>;
using MyMuonsWithCov = soa::Join<aod::FwdTracks, aod::FwdTracksCov>;
using MyMuonsMC = soa::Join<aod::FwdTracks, aod::FwdTracksCov, aod::McFwdTrackLabels>;
using MyMFTs = aod::MFTTracks;
using MyMFTCovariances = aod::MFTTracksCov;
using MyMFTsMC = soa::Join<aod::MFTTracks, aod::McMFTTrackLabels>;

using MyMuon = MyMuonsWithCov::iterator;
using MyMuonMC = MyMuonsMC::iterator;
using MyMFT = MyMFTs::iterator;
using MyMFTCovariance = MyMFTCovariances::iterator;

using SMatrix55 = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
using SMatrix5 = ROOT::Math::SVector<Double_t, 5>;

static float chi2ToScore(float chi2)
{
  return (1.f / (chi2 / 100.f + 1.f));
}

struct qaMatching {
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
  Configurable<float> fMatchingChi2ScoreMftMchLow{"cfgMatchingChi2ScoreMftMchLow", chi2ToScore(50.f), ""};

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

  using MatchingFunc_t = std::function<double(const o2::dataformats::GlobalFwdTrack& mchtrack, const o2::track::TrackParCovFwd& mfttrack)>;
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
      {"cfgChi2FunctionMatchingScoreCut_1", chi2ToScore(50.f), "Minimum score value for selecting good matches"},
      {"cfgChi2FunctionMatchingScoreCut_2", chi2ToScore(50.f), "Minimum score value for selecting good matches"},
    }};
    std::array<Configurable<float>, sChi2FunctionsNum> fMatchingPlaneZ{{
      {"cfgChi2FunctionMatchingPlaneZ_0", static_cast<float>(o2::mft::constants::mft::LayerZCoordinate()[9]), "Z position of the matching plane"},
      {"cfgChi2FunctionMatchingPlaneZ_1", static_cast<float>(o2::mft::constants::mft::LayerZCoordinate()[9]), "Z position of the matching plane"},
      {"cfgChi2FunctionMatchingPlaneZ_2", static_cast<float>(o2::mft::constants::mft::LayerZCoordinate()[9]), "Z position of the matching plane"},
    }};
  } fConfigChi2MatchingOptions;

  // ML interface
  static constexpr int sMLModelsNum = 2;
  struct : ConfigurableGroup {
    std::array<Configurable<std::string>, sMLModelsNum> fModelLabel{{
      {"cfgMLModelLabel_0", std::string{"TestModel"}, "Text label identifying this group of ML models"},
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
  } fConfigMlOptions;

  std::vector<double> binsPtMl;
  std::array<double, 1> cutValues;
  std::vector<int> cutDirMl;
  std::map<std::string, o2::analysis::MlResponseMFTMuonMatch<float>> matchingMlResponses;
  std::map<std::string, std::string> matchingChi2Functions;
  std::map<std::string, double> matchingPlanesZ;
  std::map<std::string, double> matchingScoreCuts;

  int mRunNumber{0}; // needed to detect if the run changed and trigger update of magnetic field

  Service<o2::ccdb::BasicCCDBManager> ccdbManager;
  o2::ccdb::CcdbApi fCCDBApi;

  o2::aod::rctsel::RCTFlagsChecker rctChecker{"CBT_muon_glo", false, false, true};

  // vector of all MFT-MCH(-MID) matching candidates associated to the same MCH(-MID) track,
  // to be sorted in descending order with respect to the matching score
  // the map key is the MCH(-MID) track global index
  // the elements are pairs og global muon track indexes and associated matching scores
  // for matching candidates computed with the chi2 method, the score is defined as 1/(1+chi2)
  using MatchingCandidates = std::map<int64_t, std::vector<std::pair<int64_t, double>>>;

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
    // vector of MFT-MCH track index pairs belonging to the same MC particle
    std::vector<std::pair<int64_t, int64_t>> matchablePairs;
    // vector of MCH track indexes that are expected to have an associated MFT track
    std::vector<int64_t> taggedMuons;
  };

  using CollisionInfos = std::map<int64_t, CollisionInfo>;

  std::unordered_map<int64_t, int32_t> mftTrackCovs;

  std::vector<std::pair<int64_t, int64_t>> fMatchablePairs;
  MatchingCandidates fMatchingCandidates;
  std::vector<int64_t> fTaggedMuons;

  HistogramRegistry registry{"registry", {}};
  HistogramRegistry registryMatching{"registryMatching", {}};
  HistogramRegistry registryAlignment{"registryAlignment", {}};

  std::unordered_map<std::string, o2::framework::HistPtr> matchingHistos;

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

  struct MatchingPlotter {
    o2::framework::HistPtr fTrueMatchRanking;
    o2::framework::HistPtr fTrueMatchRankingVsP;
    o2::framework::HistPtr fTrueMatchRankingVsPt;
    //-
    o2::framework::HistPtr fTrueMatchRankingGoodMCH;
    o2::framework::HistPtr fTrueMatchRankingGoodMCHVsP;
    o2::framework::HistPtr fTrueMatchRankingGoodMCHVsPt;
    //-
    o2::framework::HistPtr fTrueMatchRankingPairedMCH;
    o2::framework::HistPtr fTrueMatchRankingPairedMCHVsP;
    o2::framework::HistPtr fTrueMatchRankingPairedMCHVsPt;
    //-
    o2::framework::HistPtr fTrueMatchRankingGoodPairedMCH;
    o2::framework::HistPtr fTrueMatchRankingGoodPairedMCHVsP;
    o2::framework::HistPtr fTrueMatchRankingGoodPairedMCHVsPt;
    //-
    o2::framework::HistPtr fTrueMatchRankingGoodPairedMCHMFT;
    o2::framework::HistPtr fTrueMatchRankingGoodPairedMCHMFTVsP;
    o2::framework::HistPtr fTrueMatchRankingGoodPairedMCHMFTVsPt;
    //-
    o2::framework::HistPtr fMissedMatches;
    o2::framework::HistPtr fMissedMatchesGoodMCH;
    o2::framework::HistPtr fMissedMatchesGoodMCHMFT;
    //-
    o2::framework::HistPtr fMatchRankingWrtProd;
    o2::framework::HistPtr fMatchRankingWrtProdVsP;
    o2::framework::HistPtr fMatchRankingWrtProdVsPt;
    //-
    o2::framework::HistPtr fTrueMatchScore;
    o2::framework::HistPtr fTrueMatchScoreVsP;
    o2::framework::HistPtr fTrueMatchScoreVsPt;
    o2::framework::HistPtr fFakeMatchScore;
    o2::framework::HistPtr fFakeMatchScoreVsP;
    o2::framework::HistPtr fFakeMatchScoreVsPt;
    EfficiencyPlotter fMatchingPurityPlotter;
    EfficiencyPlotter fPairingEfficiencyPlotter;
    EfficiencyPlotter fMatchingEfficiencyPlotter;
    EfficiencyPlotter fFakeMatchingEfficiencyPlotter;

    MatchingPlotter(std::string path,
                    HistogramRegistry& registry)
      : fMatchingPurityPlotter(path + "matching-purity/", "Matching purity", registry),
        fPairingEfficiencyPlotter(path + "pairing-efficiency/", "Pairing efficiency", registry),
        fMatchingEfficiencyPlotter(path + "matching-efficiency/", "Matching efficiency", registry),
        fFakeMatchingEfficiencyPlotter(path + "fake-matching-efficiency/", "Fake matching efficiency", registry)
    {
      AxisSpec pAxis = {100, 0, 100, "p (GeV/c)"};
      AxisSpec ptAxis = {100, 0, 10, "p_{T} (GeV/c)"};

      AxisSpec indexAxis = {6, 0, 6, "ranking index"};
      std::string histName = path + "trueMatchRanking";
      std::string histTitle = "True match ranking";
      fTrueMatchRanking = registry.add(histName.c_str(), histTitle.c_str(), {HistType::kTH1F, {indexAxis}});
      histName = path + "trueMatchRankingVsP";
      histTitle = "True match ranking vs. p";
      fTrueMatchRankingVsP = registry.add(histName.c_str(), histTitle.c_str(), {HistType::kTH2F, {pAxis, indexAxis}});
      histName = path + "trueMatchRankingVsPt";
      histTitle = "True match ranking vs. p_{T}";
      fTrueMatchRankingVsPt = registry.add(histName.c_str(), histTitle.c_str(), {HistType::kTH2F, {ptAxis, indexAxis}});
      //-
      histName = path + "trueMatchRankingGoodMCH";
      histTitle = "True match ranking - good MCH tracks";
      fTrueMatchRankingGoodMCH = registry.add(histName.c_str(), histTitle.c_str(), {HistType::kTH1F, {indexAxis}});
      histName = path + "trueMatchRankingGoodMCHVsP";
      histTitle = "True match ranking vs. p - good MCH tracks";
      fTrueMatchRankingGoodMCHVsP = registry.add(histName.c_str(), histTitle.c_str(), {HistType::kTH2F, {pAxis, indexAxis}});
      histName = path + "trueMatchRankingGoodMCHVsPt";
      histTitle = "True match ranking vs. p_{T} - good MCH tracks";
      fTrueMatchRankingGoodMCHVsPt = registry.add(histName.c_str(), histTitle.c_str(), {HistType::kTH2F, {ptAxis, indexAxis}});
      //-
      histName = path + "trueMatchRankingPairedMCH";
      histTitle = "True match ranking - paired MCH tracks";
      fTrueMatchRankingPairedMCH = registry.add(histName.c_str(), histTitle.c_str(), {HistType::kTH1F, {indexAxis}});
      histName = path + "trueMatchRankingPairedMCHVsP";
      histTitle = "True match ranking vs. p - paired MCH tracks";
      fTrueMatchRankingPairedMCHVsP = registry.add(histName.c_str(), histTitle.c_str(), {HistType::kTH2F, {pAxis, indexAxis}});
      histName = path + "trueMatchRankingPairedMCHVsPt";
      histTitle = "True match ranking vs. p_{T} - paired MCH tracks";
      fTrueMatchRankingPairedMCHVsPt = registry.add(histName.c_str(), histTitle.c_str(), {HistType::kTH2F, {ptAxis, indexAxis}});
      //-
      histName = path + "trueMatchRankingGoodPairedMCH";
      histTitle = "True match ranking - good paired MCH tracks";
      fTrueMatchRankingGoodPairedMCH = registry.add(histName.c_str(), histTitle.c_str(), {HistType::kTH1F, {indexAxis}});
      histName = path + "trueMatchRankingGoodPairedMCHVsP";
      histTitle = "True match ranking vs. p - good paired MCH tracks";
      fTrueMatchRankingGoodPairedMCHVsP = registry.add(histName.c_str(), histTitle.c_str(), {HistType::kTH2F, {pAxis, indexAxis}});
      histName = path + "trueMatchRankingGoodPairedMCHVsPt";
      histTitle = "True match ranking vs. p_{T} - good paired MCH tracks";
      fTrueMatchRankingGoodPairedMCHVsPt = registry.add(histName.c_str(), histTitle.c_str(), {HistType::kTH2F, {ptAxis, indexAxis}});
      //-
      histName = path + "trueMatchRankingGoodPairedMCHMFT";
      histTitle = "True match ranking - good paired MFT and MCH tracks";
      fTrueMatchRankingGoodPairedMCHMFT = registry.add(histName.c_str(), histTitle.c_str(), {HistType::kTH1F, {indexAxis}});
      histName = path + "trueMatchRankingGoodPairedMCHMFTVsP";
      histTitle = "True match ranking vs. p - good paired MFT and MCH tracks";
      fTrueMatchRankingGoodPairedMCHMFTVsP = registry.add(histName.c_str(), histTitle.c_str(), {HistType::kTH2F, {pAxis, indexAxis}});
      histName = path + "trueMatchRankingGoodPairedMCHMFTVsPt";
      histTitle = "True match ranking vs. p_{T} - good paired MFT and MCH tracks";
      fTrueMatchRankingGoodPairedMCHMFTVsPt = registry.add(histName.c_str(), histTitle.c_str(), {HistType::kTH2F, {ptAxis, indexAxis}});

      AxisSpec missedMatchAxis = {5, 0, 5, ""};
      histName = path + "missedMatches";
      histTitle = "Missed matches";
      fMissedMatches = registry.add(histName.c_str(), histTitle.c_str(), {HistType::kTH1F, {missedMatchAxis}});
      std::get<std::shared_ptr<TH1>>(fMissedMatches)->GetXaxis()->SetBinLabel(1, "not paired");
      std::get<std::shared_ptr<TH1>>(fMissedMatches)->GetXaxis()->SetBinLabel(2, "not matched");
      std::get<std::shared_ptr<TH1>>(fMissedMatches)->GetXaxis()->SetBinLabel(3, "match missing");
      histName = path + "missedMatchesGoodMCH";
      histTitle = "Missed matches - good MCH tracks";
      fMissedMatchesGoodMCH = registry.add(histName.c_str(), histTitle.c_str(), {HistType::kTH1F, {missedMatchAxis}});
      std::get<std::shared_ptr<TH1>>(fMissedMatchesGoodMCH)->GetXaxis()->SetBinLabel(1, "not paired");
      std::get<std::shared_ptr<TH1>>(fMissedMatchesGoodMCH)->GetXaxis()->SetBinLabel(2, "not matched");
      std::get<std::shared_ptr<TH1>>(fMissedMatchesGoodMCH)->GetXaxis()->SetBinLabel(3, "match missing");
      histName = path + "missedMatchesGoodMCHMFT";
      histTitle = "Missed matches - good MFT and MCH tracks";
      fMissedMatchesGoodMCHMFT = registry.add(histName.c_str(), histTitle.c_str(), {HistType::kTH1F, {missedMatchAxis}});
      std::get<std::shared_ptr<TH1>>(fMissedMatchesGoodMCHMFT)->GetXaxis()->SetBinLabel(1, "not paired");
      std::get<std::shared_ptr<TH1>>(fMissedMatchesGoodMCHMFT)->GetXaxis()->SetBinLabel(2, "not matched");
      std::get<std::shared_ptr<TH1>>(fMissedMatchesGoodMCHMFT)->GetXaxis()->SetBinLabel(3, "match missing");

      AxisSpec scoreAxis = {100, 0, 1, "matching score"};
      histName = path + "trueMatchScore";
      histTitle = "True match score";
      fTrueMatchScore = registry.add(histName.c_str(), histTitle.c_str(), {HistType::kTH1F, {scoreAxis}});
      histName = path + "trueMatchScoreVsP";
      histTitle = "True match score vs. p";
      fTrueMatchScoreVsP = registry.add(histName.c_str(), histTitle.c_str(), {HistType::kTH2F, {pAxis, scoreAxis}});
      histName = path + "trueMatchScoreVsPt";
      histTitle = "True match score vs. p_{T}";
      fTrueMatchScoreVsPt = registry.add(histName.c_str(), histTitle.c_str(), {HistType::kTH2F, {ptAxis, scoreAxis}});

      histName = path + "fakeMatchScore";
      histTitle = "Fake match score";
      fFakeMatchScore = registry.add(histName.c_str(), histTitle.c_str(), {HistType::kTH1F, {scoreAxis}});
      histName = path + "fakeMatchScoreVsP";
      histTitle = "Fake match score vs. p";
      fFakeMatchScoreVsP = registry.add(histName.c_str(), histTitle.c_str(), {HistType::kTH2F, {pAxis, scoreAxis}});
      histName = path + "fakeMatchScoreVsPt";
      histTitle = "Fake match score vs. p_{T}";
      fFakeMatchScoreVsPt = registry.add(histName.c_str(), histTitle.c_str(), {HistType::kTH2F, {ptAxis, scoreAxis}});
    }
  };

  std::unique_ptr<MatchingPlotter> fChi2MatchingPlotter;
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

  void createMatchingHistosMC()
  {
    AxisSpec chi2Axis = {1000, 0, 1000, "chi^{2}"};
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

    AxisSpec pairableType = {2, 0, 2, ""};
    auto pairableTypeHist = registry.add((histPath + "pairableType").c_str(), "Pairable MCH tracks type", {HistType::kTH1F, {pairableType}});
    std::get<std::shared_ptr<TH1>>(pairableTypeHist)->GetXaxis()->SetBinLabel(1, "direct");
    std::get<std::shared_ptr<TH1>>(pairableTypeHist)->GetXaxis()->SetBinLabel(2, "decay");

    fChi2MatchingPlotter = std::make_unique<MatchingPlotter>(histPath + "Prod/", registryMatching);
    for (const auto& [label, func] : matchingChi2Functions) {
      fMatchingPlotters[label] = std::make_unique<MatchingPlotter>(histPath + label + "/", registryMatching);
    }
    for (const auto& [label, response] : matchingMlResponses) {
      fMatchingPlotters[label] = std::make_unique<MatchingPlotter>(histPath + label + "/", registryMatching);
    }

    fTaggedMuonsMatchingPlotter = std::make_unique<MatchingPlotter>(histPath + "Tagged/", registryMatching);
    fSelectedMuonsMatchingPlotter = std::make_unique<MatchingPlotter>(histPath + "Selected/", registryMatching);
  }

  void createAlignmentHistos()
  {
    AxisSpec dxAxis = {200, -10, 10, "#Delta x (cm)"};
    AxisSpec dyAxis = {200, -10, 10, "#Delta y (cm)"};
    AxisSpec pAxis = {100, 0, 100, "p (GeV/c)"};
    std::string histPath = "alignment/";

    registryAlignment.add((histPath + "trackDxAtMFTVsP").c_str(), "Track #Delta x vs. p", {HistType::kTH2F, {pAxis, dxAxis}});
    registryAlignment.add((histPath + "trackDyAtMFTVsP").c_str(), "Track #Delta y vs. p", {HistType::kTH2F, {pAxis, dyAxis}});
    registryAlignment.add((histPath + "trackDxAtMFTVsP_alt").c_str(), "Track #Delta x vs. p (alt method)", {HistType::kTH2F, {pAxis, dxAxis}});
    registryAlignment.add((histPath + "trackDyAtMFTVsP_alt").c_str(), "Track #Delta y vs. p (alt method)", {HistType::kTH2F, {pAxis, dyAxis}});

    registryAlignment.add((histPath + "trackDxAtMFTVsP_fake").c_str(), "Track #Delta x vs. p (fake pairs)", {HistType::kTH2F, {pAxis, dxAxis}});
    registryAlignment.add((histPath + "trackDyAtMFTVsP_fake").c_str(), "Track #Delta y vs. p (fake pairs)", {HistType::kTH2F, {pAxis, dyAxis}});
    registryAlignment.add((histPath + "trackDxAtMFTVsP_alt_fake").c_str(), "Track #Delta x vs. p (alt method, fake pairs)", {HistType::kTH2F, {pAxis, dxAxis}});
    registryAlignment.add((histPath + "trackDyAtMFTVsP_alt_fake").c_str(), "Track #Delta y vs. p (alt method, fake pairs)", {HistType::kTH2F, {pAxis, dyAxis}});
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
    mMatchingFunctionMap["matchALL"] = [](const o2::dataformats::GlobalFwdTrack& mchTrack, const o2::track::TrackParCovFwd& mftTrack) -> double {
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

      return matchChi2Track;
    };

    //________________________________________________________________________________
    mMatchingFunctionMap["matchXYPhiTanl"] = [](const o2::dataformats::GlobalFwdTrack& mchTrack, const o2::track::TrackParCovFwd& mftTrack) -> double {

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

    return matchChi2Track; };

    //________________________________________________________________________________
    mMatchingFunctionMap["matchXY"] = [](const o2::dataformats::GlobalFwdTrack& mchTrack, const o2::track::TrackParCovFwd& mftTrack) -> double {

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

    return matchChi2Track; };
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

      if (label == "" || funcName == "")
        break;

      matchingChi2Functions[label] = funcName;

      matchingScoreCuts[label] = scoreMin;
      matchingPlanesZ[label] = matchingPlaneZ;
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

      if (label == "" || modelPaths.empty() || inputFeatures.empty() || modelNames.empty())
        break;

      matchingMlResponses[label].configure(binsPtMl, mycutsMl, cutDirMl, 1);
      matchingMlResponses[label].setModelPathsCCDB(modelNames, fCCDBApi, modelPaths, fConfigCCDB.fConfigNoLaterThan.value);
      matchingMlResponses[label].cacheInputFeaturesIndices(inputFeatures);
      matchingMlResponses[label].init();

      matchingScoreCuts[label] = scoreMin;
      matchingPlanesZ[label] = matchingPlaneZ;
    }

    int nTrackTypes = static_cast<int>(o2::aod::fwdtrack::ForwardTrackTypeEnum::MCHStandaloneTrack) + 1;
    AxisSpec trackTypeAxis = {static_cast<int>(nTrackTypes), 0.0, static_cast<double>(nTrackTypes), "track type"};
    registry.add("nTracksPerType", "Number of tracks per type", {HistType::kTH1F, {trackTypeAxis}});

    AxisSpec tracksMultiplicityAxis = {10000, 0, 10000, "tracks multiplicity"};
    registry.add("tracksMultiplicityMFT", "MFT tracks multiplicity", {HistType::kTH1F, {tracksMultiplicityAxis}});

    createMatchingHistosMC();
    createAlignmentHistos();
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
      // std::cout << "  extrapToVertexWithoutBranson()" << std::endl;
      o2::mch::TrackExtrap::extrapToVertexWithoutBranson(mchTrack, z);
    } else {
      // all other cases
      // std::cout << "  extrapToZ()" << std::endl;
      o2::mch::TrackExtrap::extrapToZ(mchTrack, z);
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
    auto mchTrack = mExtrap.FwdtoMCH(track);

    float absFront = -90.f;
    float absBack = -505.f;

    if (fwdtrack.getZ() < absBack && z > absFront) {
      // extrapolation through the absorber in the upstream direction
      o2::mch::TrackExtrap::extrapToVertexWithoutBranson(mchTrack, z);
    } else {
      // all other cases
      o2::mch::TrackExtrap::extrapToZ(mchTrack, z);
    }

    auto proptrack = mExtrap.MCHtoFwd(mchTrack);
    o2::dataformats::GlobalFwdTrack propmuon;
    propmuon.setParameters(proptrack.getParameters());
    propmuon.setZ(proptrack.getZ());
    propmuon.setCovariances(proptrack.getCovariances());

    return propmuon;
  }

  o2::dataformats::GlobalFwdTrack PropagateToZMFT(const o2::dataformats::GlobalFwdTrack& mftTrack, const double z)
  {
    o2::dataformats::GlobalFwdTrack fwdtrack{mftTrack};
    fwdtrack.propagateToZhelix(z, mBzAtMftCenter);
    return fwdtrack;
  }

  template <typename TMFT, typename CMFT>
  o2::dataformats::GlobalFwdTrack PropagateToZMFT(const TMFT& mftTrack, const CMFT& mftCov, const double z)
  {
    o2::dataformats::GlobalFwdTrack fwdtrack = FwdToTrackPar(mftTrack, mftCov);
    return PropagateToZMFT(fwdtrack, z);
  }

  template <class EVT, class BC, class TMUON, class TMFT>
  void InitCollisions(EVT const& collisions,
                      BC const& bcs,
                      TMUON const& muonTracks,
                      TMFT const& mftTracks,
                      CollisionInfos& collisionInfos)
  {
    collisionInfos.clear();

    // fill collision information for global muon tracks (MFT-MCH-MID matches)
    for (auto muonTrack : muonTracks) {
      if (!muonTrack.has_collision())
        continue;

      auto collision = collisions.rawIteratorAt(muonTrack.collisionId());
      int64_t collisionIndex = collision.globalIndex();

      // std::cout << std::format("Collision indexes: {} / {}", collisionIndex, muonTrack.collisionId()) << std::endl;

      auto bc = bcs.rawIteratorAt(collision.bcId());

      auto& collisionInfo = collisionInfos[collisionIndex];
      collisionInfo.index = collisionIndex;
      collisionInfo.bc = bc.globalBC();
      collisionInfo.zVertex = collision.posZ();

      if (static_cast<int>(muonTrack.trackType()) > 2) {
        // standalone MCH or MCH-MID tracks
        int64_t mchTrackIndex = muonTrack.globalIndex();
        // if (!IsGoodMuon(muonTrack, collision)) continue;
        collisionInfo.mchTracks.push_back(mchTrackIndex);
      } else {
        // global muon tracks (MFT-MCH or MFT-MCH-MID)
        int64_t muonTrackIndex = muonTrack.globalIndex();
        double matchingScore = chi2ToScore(muonTrack.chi2MatchMCHMFT());
        auto const& mchTrack = muonTrack.template matchMCHTrack_as<TMUON>();
        int64_t mchTrackIndex = mchTrack.globalIndex();

        // check if a vector of global muon candidates is already available for the current MCH index
        // if not, initialize a new one and add the current global muon track
        // bool globalMuonTrackFound = false;
        auto matchingCandidateIterator = collisionInfo.matchingCandidates.find(mchTrackIndex);
        if (matchingCandidateIterator != collisionInfo.matchingCandidates.end()) {
          matchingCandidateIterator->second.push_back(std::make_pair(muonTrackIndex, matchingScore));
          // globalMuonTrackFound = true;
        } else {
          collisionInfo.matchingCandidates[mchTrackIndex].push_back(std::make_pair(muonTrackIndex, matchingScore));
        }
      }
    }

    // fill collision information for MFT standalone tracks
    for (auto mftTrack : mftTracks) {
      if (!mftTrack.has_collision())
        continue;

      auto collision = collisions.rawIteratorAt(mftTrack.collisionId());
      int64_t collisionIndex = collision.globalIndex();

      auto bc = bcs.rawIteratorAt(collision.bcId());

      int64_t mftTrackIndex = mftTrack.globalIndex();

      auto& collisionInfo = collisionInfos[collisionIndex];
      collisionInfo.index = collisionIndex;
      collisionInfo.bc = bc.globalBC();
      collisionInfo.zVertex = collision.posZ();

      collisionInfo.mftTracks.push_back(mftTrackIndex);
    }

    // sort the vectors of matching candidates in ascending order based on the matching score value
    auto compareMatchingScore = [](std::pair<int64_t, double> track1, std::pair<int64_t, double> track2) -> bool {
      return (track1.second > track2.second);
    };

    for (auto& [collisionIndex, collisionInfo] : collisionInfos) {
      for (auto& [mchIndex, globalTracksVector] : collisionInfo.matchingCandidates) {
        std::sort(globalTracksVector.begin(), globalTracksVector.end(), compareMatchingScore);
      }
    }
  }

  template <class TMUON, class TMFT>
  void GetMatchablePairs(const CollisionInfo& collisionInfo,
                         TMUON const& muonTracks,
                         TMFT const& mftTracks,
                         std::vector<std::pair<int64_t, int64_t>>& matchablePairs)
  {
    matchablePairs.clear();
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
      int64_t muonMcTrackIndex = muonMcParticle.globalIndex();

      for (const auto& mftTrack : mftTracks) {
        // skip tracks that do not have an associated MC particle
        if (!mftTrack.has_mcParticle())
          continue;
        // get the index associated to the MC particle
        auto mftMcParticle = mftTrack.mcParticle();
        int64_t mftMcTrackIndex = mftMcParticle.globalIndex();

        if (muonMcTrackIndex == mftMcTrackIndex) {
          matchablePairs.emplace_back(std::make_pair(static_cast<int64_t>(muonTrack.globalIndex()),
                                                     static_cast<int64_t>(mftTrack.globalIndex())));
        } else {
          // check if the muon particle is a decay product of the MFT particle
          for (auto& motherParticle : muonMcParticle.template mothers_as<aod::McParticles>()) {
            if (motherParticle.globalIndex() == mftMcTrackIndex) {
              matchablePairs.emplace_back(std::make_pair(static_cast<int64_t>(muonTrack.globalIndex()),
                                                         static_cast<int64_t>(mftTrack.globalIndex())));
              break;
            }
          }
        }
      }
    }
  }

  template <class TMUON>
  int GetTrueMatchIndex(TMUON const& muonTracks,
                        const std::vector<std::pair<int64_t, double>>& matchCandidatesVector,
                        const std::vector<std::pair<int64_t, int64_t>>& matchablePairs)
  {
    // find the index of the matching candidate that corresponds to the true match
    // index=1 corresponds to the leading candidate
    // index=0 means no candidate was found that corresponds to the true match
    int trueMatchIndex = 0;
    for (size_t i = 0; i < matchCandidatesVector.size(); i++) {
      auto const& muonTrack = muonTracks.rawIteratorAt(matchCandidatesVector[i].first);

      if (IsTrueGlobalMatching(muonTrack, matchablePairs)) {
        trueMatchIndex = i + 1;
        break;
      }
    }
    return trueMatchIndex;
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

      auto const& muonTrack0 = muonTracks.rawIteratorAt(globalTracksVector[0].first);
      auto const& muonTrack1 = muonTracks.rawIteratorAt(globalTracksVector[1].first);

      double chi2diff = muonTrack1.chi2MatchMCHMFT() - muonTrack0.chi2MatchMCHMFT();
      if (chi2diff < fMuonTaggingChi2DiffLow)
        continue;

      taggedMuons.emplace_back(mchIndex);
    }
  }

  template <class C, class TMUON, class TMFT>
  void FillMatchingPlotsMC(C const& collision,
                           TMUON const& muonTracks,
                           TMFT const& mftTracks,
                           const MatchingCandidates& matchingCandidates,
                           const std::vector<std::pair<int64_t, int64_t>>& matchablePairs,
                           double matchingScoreCut,
                           MatchingPlotter* plotter)
  {
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
      bool isGoodPair = false;
      if (hasMatchablePair) {
        auto const& pairedMftTrack = mftTracks.rawIteratorAt(matchablePair.value().second);
        isGoodPair = isGoodMCH && IsGoodMFT(pairedMftTrack);
      }

      // std::cout << std::format("Checking matchable MCH track #{}", mchIndex) << std::endl;

      // find the index of the matching candidate that corresponds to the true match
      // index=1 corresponds to the leading candidate
      // index=0 means no candidate was found that corresponds to the true match
      int trueMatchIndex = GetTrueMatchIndex(muonTracks, globalTracksVector, matchablePairs);

      std::get<std::shared_ptr<TH1>>(plotter->fTrueMatchRanking)->Fill(trueMatchIndex);
      std::get<std::shared_ptr<TH2>>(plotter->fTrueMatchRankingVsP)->Fill(mchMom, trueMatchIndex);
      std::get<std::shared_ptr<TH2>>(plotter->fTrueMatchRankingVsPt)->Fill(mchPt, trueMatchIndex);

      if (isGoodMCH) {
        std::get<std::shared_ptr<TH1>>(plotter->fTrueMatchRankingGoodMCH)->Fill(trueMatchIndex);
        std::get<std::shared_ptr<TH2>>(plotter->fTrueMatchRankingGoodMCHVsP)->Fill(mchMom, trueMatchIndex);
        std::get<std::shared_ptr<TH2>>(plotter->fTrueMatchRankingGoodMCHVsPt)->Fill(mchPt, trueMatchIndex);
      }

      if (isPairedMCH) {
        std::get<std::shared_ptr<TH1>>(plotter->fTrueMatchRankingPairedMCH)->Fill(trueMatchIndex);
        std::get<std::shared_ptr<TH2>>(plotter->fTrueMatchRankingPairedMCHVsP)->Fill(mchMom, trueMatchIndex);
        std::get<std::shared_ptr<TH2>>(plotter->fTrueMatchRankingPairedMCHVsPt)->Fill(mchPt, trueMatchIndex);
      }

      if (isGoodMCH && isPairedMCH) {
        std::get<std::shared_ptr<TH1>>(plotter->fTrueMatchRankingGoodPairedMCH)->Fill(trueMatchIndex);
        std::get<std::shared_ptr<TH2>>(plotter->fTrueMatchRankingGoodPairedMCHVsP)->Fill(mchMom, trueMatchIndex);
        std::get<std::shared_ptr<TH2>>(plotter->fTrueMatchRankingGoodPairedMCHVsPt)->Fill(mchPt, trueMatchIndex);
      }

      if (isGoodPair) {
        std::get<std::shared_ptr<TH1>>(plotter->fTrueMatchRankingGoodPairedMCHMFT)->Fill(trueMatchIndex);
        std::get<std::shared_ptr<TH2>>(plotter->fTrueMatchRankingGoodPairedMCHMFTVsP)->Fill(mchMom, trueMatchIndex);
        std::get<std::shared_ptr<TH2>>(plotter->fTrueMatchRankingGoodPairedMCHMFTVsP)->Fill(mchPt, trueMatchIndex);
      }

      if (trueMatchIndex == 0) {
        // missed matches
        if (!isPairedMCH) {
          // the MCH track does not have a corresponding MFT track for matching
          std::get<std::shared_ptr<TH1>>(plotter->fMissedMatches)->Fill(0);
          if (isGoodMCH) {
            std::get<std::shared_ptr<TH1>>(plotter->fMissedMatchesGoodMCH)->Fill(0);
          }
        } else {
          if (globalTracksVector.empty()) {
            // the MCH track was not matched to any MFT track
            std::get<std::shared_ptr<TH1>>(plotter->fMissedMatches)->Fill(1);
            if (isGoodMCH) {
              std::get<std::shared_ptr<TH1>>(plotter->fMissedMatchesGoodMCH)->Fill(1);
            }
            if (isGoodPair) {
              std::get<std::shared_ptr<TH1>>(plotter->fMissedMatchesGoodMCHMFT)->Fill(1);
            }
          } else {
            // the correct match is not among the stored candidates
            std::get<std::shared_ptr<TH1>>(plotter->fMissedMatches)->Fill(2);
            if (isGoodMCH) {
              std::get<std::shared_ptr<TH1>>(plotter->fMissedMatchesGoodMCH)->Fill(2);
            }
            if (isGoodPair) {
              std::get<std::shared_ptr<TH1>>(plotter->fMissedMatchesGoodMCHMFT)->Fill(2);
            }
          }
        }
      }
    }

    // ====================================
    // Matching properties

    for (auto [mchIndex, globalTracksVector] : matchingCandidates) {
      if (globalTracksVector.size() < 1)
        continue;

      // get the standalone MCH track
      auto const& mchTrack = muonTracks.rawIteratorAt(mchIndex);

      // get leading matching candidate
      auto const& muonTrack = muonTracks.rawIteratorAt(globalTracksVector[0].first);

      auto const& mftTrack = muonTrack.template matchMFTTrack_as<TMFT>();

      // skip global muon tracks that do not pass the MCH and MFT quality cuts
      if (!IsGoodGlobalMuon(mchTrack, collision))
        continue;
      if (!IsGoodMFT(mftTrack))
        continue;

      double matchingScore = globalTracksVector[0].second;
      double mchMom = mchTrack.p();
      double mchPt = mchTrack.pt();

      // check the matching quality, but set the minimum matching score to zero
      bool isGoodMatchNoScore = IsGoodGlobalMatching(muonTrack, matchingScore, 0);
      // bool isGoodMatch = IsGoodGlobalMatching(muonTrack, matchingScore, matchingScoreCut);

      bool isTrueMatch = IsTrueGlobalMatching(muonTrack, matchablePairs);

      // matching score analysis
      if (isGoodMatchNoScore) {
        if (isTrueMatch) {
          std::get<std::shared_ptr<TH1>>(plotter->fTrueMatchScore)->Fill(matchingScore);
          std::get<std::shared_ptr<TH2>>(plotter->fTrueMatchScoreVsP)->Fill(mchMom, matchingScore);
          std::get<std::shared_ptr<TH2>>(plotter->fTrueMatchScoreVsPt)->Fill(mchPt, matchingScore);
        } else {
          std::get<std::shared_ptr<TH1>>(plotter->fFakeMatchScore)->Fill(matchingScore);
          std::get<std::shared_ptr<TH2>>(plotter->fFakeMatchScoreVsP)->Fill(mchMom, matchingScore);
          std::get<std::shared_ptr<TH2>>(plotter->fFakeMatchScoreVsPt)->Fill(mchPt, matchingScore);
        }
      }
    }

    // ====================================
    // Matching purity

    for (auto [mchIndex, globalTracksVector] : matchingCandidates) {
      if (globalTracksVector.size() < 1)
        continue;

      // get the leading matching candidate
      auto const& muonTrack = muonTracks.rawIteratorAt(globalTracksVector[0].first);
      double matchingScore = globalTracksVector[0].second;

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

      // fill matching purity plots
      plotter->fMatchingPurityPlotter.Fill(mchTrack, isTrueMatch);
    }

    // ====================================
    // Matching efficiencies

    // outer loop on matchable pairs
    for (auto [matchableMchIndex, matchableMftIndex] : matchablePairs) {
      // get the standalone MCH track
      // std::cout << std::format("Retrieving paired tracks: {} / {}", matchableMchIndex, matchableMftIndex) << std::endl;
      auto const& mchTrack = muonTracks.rawIteratorAt(matchableMchIndex);
      auto const& pairedMftTrack = mftTracks.rawIteratorAt(matchableMftIndex);
      // std::cout << std::format("... done.") << std::endl;

      // skip  track pairs that do not pass the MCH and MFT quality cuts
      // we only consider matchable pairs that fulfill the track quality requirements
      // std::cout << std::format("Checking tracks...") << std::endl;
      if (!IsGoodGlobalMuon(mchTrack, collision))
        continue;
      // std::cout << std::format("... MCH ok.") << std::endl;
      if (!IsGoodMFT(pairedMftTrack))
        continue;
      // std::cout << std::format("... MFT ok.") << std::endl;

      bool goodMatchFound = false;
      bool isTrueMatch = false;

      // check if we have some matching candidates for the current matchable MCH track
      if (matchingCandidates.count(matchableMchIndex) > 0) {
        // std::cout << std::format("Getting matching candidates for MCH track {}", matchableMchIndex) << std::endl;
        const auto& globalTracksVector = matchingCandidates.at(static_cast<int64_t>(matchableMchIndex));
        // std::cout << std::format("Number of matching candidates: {}", globalTracksVector.size()) << std::endl;
        if (!globalTracksVector.empty()) {
          // get the leading matching candidate
          // std::cout << std::format("Getting leading matching candidate: {}", globalTracksVector[0].first) << std::endl;
          auto const& muonTrack = muonTracks.rawIteratorAt(globalTracksVector[0].first);
          double matchingScore = globalTracksVector[0].second;
          // std::cout << std::format("... done.") << std::endl;

          // get the standalone MFT track
          // std::cout << std::format("Getting MFT track") << std::endl;
          auto const& mftTrack = muonTrack.template matchMFTTrack_as<TMFT>();
          auto mftIndex = mftTrack.globalIndex();
          // std::cout << std::format("... done.") << std::endl;

          // a good match must pass the MFT and matching quality cuts
          // the MCH track quality is already checked in the outer loop
          // std::cout << std::format("Checking matched track quality") << std::endl;
          goodMatchFound = IsGoodMFT(mftTrack) && IsGoodGlobalMatching(muonTrack, matchingScore, matchingScoreCut);
          isTrueMatch = (mftIndex == matchableMftIndex);
          // std::cout << std::format("... done: {}", goodMatchFound) << std::endl;
        }
      }

      // fill matching efficiency plots
      plotter->fPairingEfficiencyPlotter.Fill(mchTrack, goodMatchFound);
      plotter->fMatchingEfficiencyPlotter.Fill(mchTrack, (goodMatchFound && isTrueMatch));
      plotter->fFakeMatchingEfficiencyPlotter.Fill(mchTrack, (goodMatchFound && !isTrueMatch));
    }
  }

  template <class C, class TMUON, class TMFT, class CMFT>
  void RunChi2Matching(C const& collisions,
                       TMUON const& muonTracks,
                       TMFT const& /*mftTracks*/,
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

    // std::cout << std::format("Processing {} matches with chi2 function {}", matchingCandidates.size(), funcName) << std::endl;

    auto matchingFunc = mMatchingFunctionMap.at(funcName);
    for (auto& [mchIndex, globalTracksVector] : matchingCandidates) {
      auto const& mchTrack = muonTracks.rawIteratorAt(mchIndex);

      for (auto& [muonIndex, score] : globalTracksVector) {
        auto const& muonTrack = muonTracks.rawIteratorAt(muonIndex);
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
        auto mchTrackPropAlt = FwdToTrackPar(mchTrack, mchTrack);

        // extrapolate to the matching plane
        auto matchingPlaneZ = matchingPlanesZ[label];
        // std::cout << std::format("Extrapolating tracks to z={}", matchingPlaneZ) << std::endl;
        if (matchingPlaneZ < 0.) {
          mftTrackProp = PropagateToZMFT(mftTrackProp, matchingPlaneZ);
          mchTrackProp = PropagateToZMCH(mchTrackProp, matchingPlaneZ);
          auto mchTrackAtVertex = VarManager::PropagateMuon(mchTrack, collision, VarManager::kToVertex);
          mchTrackPropAlt = PropagateToZMCH(mchTrackAtVertex, matchingPlaneZ);
        }

        // run the chi2 matching function
        float matchingChi2 = matchingFunc(mchTrackProp, mftTrackProp);
        float matchingScore = chi2ToScore(matchingChi2);
        // std::cout << std::format("Matching chi2: {}, original chi2: {}", matchingChi2, muonTrack.chi2MatchMCHMFT()) << std::endl;

        // check if a vector of global muon candidates is already available for the current MCH index
        // if not, initialize a new one and add the current global muon track
        auto matchingCandidateIterator = newMatchingCandidates.find(mchIndex);
        if (matchingCandidateIterator != newMatchingCandidates.end()) {
          matchingCandidateIterator->second.push_back(std::make_pair(muonIndex, matchingScore));
        } else {
          newMatchingCandidates[mchIndex].push_back(std::make_pair(muonIndex, matchingScore));
        }
      }
    }
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
      for (auto& [muonIndex, score] : globalTracksVector) {
        auto const& muonTrack = muonTracks.rawIteratorAt(muonIndex);
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
        float matchingScore = output[0];
        // std::cout << std::format("Matching score: {}, Chi2: {}", matchingScore, muonTrack.chi2MatchMCHMFT()) << std::endl;

        // check if a vector of global muon candidates is already available for the current MCH index
        // if not, initialize a new one and add the current global muon track
        auto matchingCandidateIterator = newMatchingCandidates.find(mchIndex);
        if (matchingCandidateIterator != newMatchingCandidates.end()) {
          matchingCandidateIterator->second.push_back(std::make_pair(muonIndex, matchingScore));
        } else {
          newMatchingCandidates[mchIndex].push_back(std::make_pair(muonIndex, matchingScore));
        }
      }
    }
  }

  template <class C, class TMUON, class TMFT, class CMFT>
  void ProcessCollisionMC(const CollisionInfo& collisionInfo,
                          C const& collisions,
                          TMUON const& muonTracks,
                          TMFT const& mftTracks,
                          CMFT const& mftCovs)
  {
    // std::cout << std::endl << std::format("ProcessMatching() collision #{}", collisionInfo.index) << std::endl;
    auto collision = collisions.rawIteratorAt(collisionInfo.index);

    std::vector<std::pair<int64_t, int64_t>> matchablePairs;
    GetMatchablePairs(collisionInfo, muonTracks, mftTracks, matchablePairs);
    for (auto [mchIndex, mftIndex] : matchablePairs) {
      auto const& muonTrack = muonTracks.rawIteratorAt(mchIndex);
      auto muonMcParticle = muonTrack.mcParticle();
      int64_t muonMcTrackIndex = muonMcParticle.globalIndex();
      double mchMom = muonTrack.p();

      // std::cout << std::format("TOTO1 MCH track #{} type={} p={:0.3} - MFT track #{} - collision #{} - has_collision={}",
      //     mchIndex, muonTrack.trackType(), mchMom, mftIndex, muonTrack.collisionId(), muonTrack.has_collision()) << std::endl;

      auto const& mftTrack = mftTracks.rawIteratorAt(mftIndex);
      auto mftMcParticle = mftTrack.mcParticle();
      int64_t mftMcTrackIndex = mftMcParticle.globalIndex();

      if (muonMcTrackIndex == mftMcTrackIndex) {
        registry.get<TH1>(HIST("matching/MC/pairableType"))->Fill(0);
      } else {
        registry.get<TH1>(HIST("matching/MC/pairableType"))->Fill(1);
      }

      if (mftTrackCovs.count(mftTrack.globalIndex()) < 1) {
        // std::cout << std::format("Covariance matrix for MFT track #{} not found", mftTrack.globalIndex()) << std::endl;
        continue;
      }
      auto const& mftTrackCov = mftCovs.rawIteratorAt(mftTrackCovs[mftTrack.globalIndex()]);

      // propagate tracks at the last MFT plane
      auto mftTrackProp = FwdToTrackPar(mftTrack, mftTrackCov);
      auto mchTrackProp = FwdToTrackPar(muonTrack, muonTrack);
      auto mchTrackPropAlt = FwdToTrackPar(muonTrack, muonTrack);

      auto z = o2::mft::constants::mft::LayerZCoordinate()[9];
      // std::cout << std::format("Extrapolating tracks to z={}", matchingPlaneZ) << std::endl;
      mftTrackProp = PropagateToZMFT(mftTrackProp, z);
      mchTrackProp = PropagateToZMCH(mchTrackProp, z);
      auto mchTrackAtVertex = VarManager::PropagateMuon(muonTrack, collision, VarManager::kToVertex);
      mchTrackPropAlt = PropagateToZMCH(mchTrackAtVertex, z);

      float dx = mchTrackProp.getX() - mftTrackProp.getX();
      float dy = mchTrackProp.getY() - mftTrackProp.getY();
      registryAlignment.get<TH2>(HIST("alignment/trackDxAtMFTVsP"))->Fill(mchMom, dx);
      registryAlignment.get<TH2>(HIST("alignment/trackDyAtMFTVsP"))->Fill(mchMom, dy);

      float dxAlt = mchTrackPropAlt.getX() - mftTrackProp.getX();
      float dyAlt = mchTrackPropAlt.getY() - mftTrackProp.getY();
      registryAlignment.get<TH2>(HIST("alignment/trackDxAtMFTVsP_alt"))->Fill(mchMom, dxAlt);
      registryAlignment.get<TH2>(HIST("alignment/trackDyAtMFTVsP_alt"))->Fill(mchMom, dyAlt);

      // std::cout << std::format("MCH track position: x={} y={}", mchTrackProp.getX(), mchTrackProp.getY()) << std::endl;
      // std::cout << std::format("MCH track pos. alt: x={} y={}", mchTrackPropAlt.getX(), mchTrackPropAlt.getY()) << std::endl;
      // std::cout << std::format("MFT track position: x={} y={}", mftTrackProp.getX(), mftTrackProp.getY()) << std::endl;
    }

    // plot MFT-MCH tracks difference at matching plane for wrong pairs
    // outer loop on MCH tracks associated to the current collision
    for (auto const& mchTrackIndex : collisionInfo.mchTracks) {
      auto matchablePair = GetMatchablePairForMCH(mchTrackIndex, matchablePairs);
      if (!matchablePair.has_value())
        continue;

      auto const& mchTrack = muonTracks.rawIteratorAt(mchTrackIndex);
      double mchMom = mchTrack.p();

      // std::cout << std::format("TOTO2 MCH track #{} type={} p={:0.3} - MFT track #{} - collision #{}",
      //     mchTrackIndex, mchTrack.trackType(), mchMom, matchablePair.value().second, mchTrack.collisionId()) << std::endl;

      // extrapolate the MCH track to the last MFT plane
      auto z = o2::mft::constants::mft::LayerZCoordinate()[9];
      auto mchTrackProp = PropagateToZMCH(FwdToTrackPar(mchTrack, mchTrack), z);
      auto mchTrackAtVertex = VarManager::PropagateMuon(mchTrack, collision, VarManager::kToVertex);
      auto mchTrackPropAlt = PropagateToZMCH(mchTrackAtVertex, z);

      // inner loop on MFT tracks associated to the current collision
      for (auto const& mftTrackIndex : collisionInfo.mftTracks) {
        // skip true matches
        if (mftTrackIndex == matchablePair.value().second)
          continue;
        auto const& mftTrack = mftTracks.rawIteratorAt(mftTrackIndex);

        // get the corresponding covariance matrix
        if (mftTrackCovs.count(mftTrackIndex) < 1) {
          continue;
        }
        auto const& mftTrackCov = mftCovs.rawIteratorAt(mftTrackCovs[mftTrackIndex]);

        auto mftTrackProp = PropagateToZMFT(FwdToTrackPar(mftTrack, mftTrackCov), z);

        float dx = mchTrackProp.getX() - mftTrackProp.getX();
        float dy = mchTrackProp.getY() - mftTrackProp.getY();
        registryAlignment.get<TH2>(HIST("alignment/trackDxAtMFTVsP_fake"))->Fill(mchMom, dx);
        registryAlignment.get<TH2>(HIST("alignment/trackDyAtMFTVsP_fake"))->Fill(mchMom, dy);

        float dxAlt = mchTrackPropAlt.getX() - mftTrackProp.getX();
        float dyAlt = mchTrackPropAlt.getY() - mftTrackProp.getY();
        registryAlignment.get<TH2>(HIST("alignment/trackDxAtMFTVsP_alt_fake"))->Fill(mchMom, dxAlt);
        registryAlignment.get<TH2>(HIST("alignment/trackDyAtMFTVsP_alt_fake"))->Fill(mchMom, dyAlt);
      }
    }

    // Chi2-based matching analysis
    FillMatchingPlotsMC(collision, muonTracks, mftTracks, collisionInfo.matchingCandidates, matchablePairs, fMatchingChi2ScoreMftMchLow, fChi2MatchingPlotter.get());
    for (auto& [label, func] : matchingChi2Functions) {
      MatchingCandidates matchingCandidates;
      RunChi2Matching(collisions, muonTracks, mftTracks, mftCovs, label, collisionInfo.matchingCandidates, matchingCandidates);

      auto* plotter = fMatchingPlotters.at(label).get();
      double matchingScoreCut = matchingScoreCuts.at(label);

      // std::cout << std::format("Calling FillMatchingPlotsMC() for label \"{}\" and score cut {}", label, matchingScoreCut) << std::endl;
      FillMatchingPlotsMC(collision, muonTracks, mftTracks, matchingCandidates, matchablePairs, matchingScoreCut, plotter);
    }

    // ML-based matching analysis
    for (auto& [label, mlResponse] : matchingMlResponses) {
      MatchingCandidates matchingCandidates;
      RunMLMatching(collisions, muonTracks, mftTracks, mftCovs, label, collisionInfo.matchingCandidates, matchingCandidates);

      auto* plotter = fMatchingPlotters.at(label).get();
      double matchingScoreCut = matchingScoreCuts.at(label);

      FillMatchingPlotsMC(collision, muonTracks, mftTracks, matchingCandidates, matchablePairs, matchingScoreCut, plotter);
    }

    // Muons tagging
    for (auto [mchIndex, mftIndex] : matchablePairs) {
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
    for (auto mchIndex : selectedMuons) {
      auto const& mchTrack = muonTracks.rawIteratorAt(mchIndex);
      if (!mchTrack.has_collision())
        continue;
      auto collision = collisions.rawIteratorAt(mchTrack.collisionId());

      auto mchTrackAtVertex = VarManager::PropagateMuon(mchTrack, collision, VarManager::kToVertex);

      // extrapolate to the matching plane
      auto z = o2::mft::constants::mft::LayerZCoordinate()[9];
      auto mchTrackProp = PropagateToZMCH(mchTrackAtVertex, z);

      registry.get<TH2>(HIST("matching/MC/selectedMCHTracksAtMFT"))->Fill(mchTrackProp.getX(), mchTrackProp.getY());

      bool isPairedMCH = IsMatchableMCH(static_cast<int64_t>(mchIndex), matchablePairs);
      if (isPairedMCH) {
        registry.get<TH2>(HIST("matching/MC/selectedMCHTracksAtMFTTrue"))->Fill(mchTrackProp.getX(), mchTrackProp.getY());
      } else {
        registry.get<TH2>(HIST("matching/MC/selectedMCHTracksAtMFTFake"))->Fill(mchTrackProp.getX(), mchTrackProp.getY());
      }
    }

    MatchingCandidates selectedMatchingCandidates;
    for (auto [mchIndex, globalTracksVector] : collisionInfo.matchingCandidates) {
      if (std::find(selectedMuons.begin(), selectedMuons.end(), mchIndex) != selectedMuons.end()) {
        selectedMatchingCandidates[mchIndex] = globalTracksVector;
      }
    }
    FillMatchingPlotsMC(collision, muonTracks, mftTracks, selectedMatchingCandidates, matchablePairs, fMatchingChi2ScoreMftMchLow, fSelectedMuonsMatchingPlotter.get());

    std::vector<int64_t> taggedMuons;
    GetTaggedMuons(collisionInfo, muonTracks, selectedMuons, taggedMuons);

    MatchingCandidates taggedMatchingCandidates;
    for (auto [mchIndex, globalTracksVector] : collisionInfo.matchingCandidates) {
      if (std::find(taggedMuons.begin(), taggedMuons.end(), mchIndex) != taggedMuons.end()) {
        taggedMatchingCandidates[mchIndex] = globalTracksVector;
      }
    }
    FillMatchingPlotsMC(collision, muonTracks, mftTracks, taggedMatchingCandidates, matchablePairs, fMatchingChi2ScoreMftMchLow, fTaggedMuonsMatchingPlotter.get());
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

    InitCollisions(collisions, bcs, muonTracks, mftTracks, fCollisionInfos);

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
