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
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "GlobalTracking/MatchGlobalFwd.h"
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "Common/DataModel/EventSelection.h"
#include "MFTTracking/Constants.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "PWGDQ/Core/VarManager.h"
#include "PWGDQ/Core/MuonMatchingMlResponse.h"

#include <string>
#include <unordered_map>
#include <map>
#include <limits>

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
  Configurable<float> fMatchingChi2ScoreMftMchLow{"cfgMatchingChi2ScoreMftMchLow", (1.f / 51.f), ""};

  ////   Variables for selecting tagged muons
  Configurable<int> fMuonTaggingNCrossedMftPlanesLow{"cfgMuonTaggingNCrossedMftPlanesLow", 5, ""};
  Configurable<float> fMuonTaggingTrackChi2MchUp{"cfgMuonTaggingTrackChi2MchUp", 5.f, ""};
  Configurable<float> fMuonTaggingPMchLow{"cfgMuonTaggingPMchLow", 0.0f, ""};
  Configurable<float> fMuonTaggingPtMchLow{"cfgMuonTaggingPtMchLow", 0.7f, ""};
  Configurable<float> fMuonTaggingEtaMchLow{"cfgMuonTaggingEtaMchLow", -4.0f, ""};
  Configurable<float> fMuonTaggingEtaMchUp{"cfgMuonTaggingEtaMchUp", -2.5f, ""};
  Configurable<float> fMuonTaggingRabsLow{"cfgMuonTaggingRabsLow", 17.6f, ""};
  Configurable<float> fMuonTaggingRabsUp{"cfgMuonTaggingRabsUp", 89.5f, ""};
  Configurable<float> fMuonTaggingSigmaPdcaUp{"cfgMuonTaggingPdcaUp", 6.f, ""};

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

  double mBzAtMftCenter{ 0 };

  o2::globaltracking::MatchGlobalFwd mExtrap;

  using MatchingFunc_t = std::function<double(const o2::dataformats::GlobalFwdTrack& mchtrack, const o2::track::TrackParCovFwd& mfttrack)>;
  std::map<std::string, MatchingFunc_t> mMatchingFunctionMap; ///< MFT-MCH Matching function

  // Chi2 matching interface
  static constexpr int sChi2FunctionsNum = 2;
  struct : ConfigurableGroup {
    std::array<Configurable<std::string>, sChi2FunctionsNum> fFunctionLabel{{
      {"cfgChi2FunctionLabel_0", std::string{"MatchAll"}, "Text label identifying this chi2 matching method"},
      {"cfgChi2FunctionLabel_1", std::string{"MatchXYPhiTanl"}, "Text label identifying this chi2 matching method"},
    }};
    std::array<Configurable<std::string>, sChi2FunctionsNum> fFunctionName{{
      {"cfgChi2FunctionNames_0", std::string{"matchALL"}, "Name of the chi2 matching function"},
      {"cfgChi2FunctionNames_1", std::string{"matchXYPhiTanl"}, "Name of the chi2 matching function"}
    }};
    std::array<Configurable<float>, sChi2FunctionsNum> fMatchingScoreCut{{
      {"cfgChi2FunctionMatchingScoreCut_0", 1.f / 51.f, "Minimum score value for selecting good matches"},
      {"cfgChi2FunctionMatchingScoreCut_1", 1.f / 51.f, "Minimum score value for selecting good matches"},
    }};
    std::array<Configurable<float>, sChi2FunctionsNum> fMatchingPlaneZ{{
      {"cfgChi2FunctionMatchingPlaneZ_0", static_cast<float>(o2::mft::constants::mft::LayerZCoordinate()[9]), "Z position of the matching plane"},
      {"cfgChi2FunctionMatchingPlaneZ_1", static_cast<float>(o2::mft::constants::mft::LayerZCoordinate()[9]), "Z position of the matching plane"},
    }};
  } fConfigChi2MatchingOptions;

  // ML interface
  static constexpr int sMLModelsNum = 2;
  struct : ConfigurableGroup {
    std::array<Configurable<std::string>, sMLModelsNum> fModelLabel{{
      {"cfgMLModelLabel_0", std::string{"TestModel"}, "Text label identifying this group of ML models"},
      {"cfgMLModelLabel_1", std::string{""}, "Text label identifying this group of ML models"},
    }};
    std::array<Configurable<std::vector<std::string>>, sMLModelsNum> fModelPathsCCDB{{
      {"cfgMLModelPathsCCDB_0", std::vector<std::string>{"Users/m/mcoquet/MLTest"}, "Paths of models on CCDB"},
      {"cfgMLModelPathsCCDB_1", std::vector<std::string>{}, "Paths of models on CCDB"}
    }};
    std::array<Configurable<std::vector<std::string>>, sMLModelsNum> fInputFeatures{{
      {"cfgMLInputFeatures_0", std::vector<std::string>{"chi2MCHMFT"}, "Names of ML model input features"},
      {"cfgMLInputFeatures_1", std::vector<std::string>{}, "Names of ML model input features"}
    }};
    std::array<Configurable<std::vector<std::string>>, sMLModelsNum> fModelNames{{
      {"cfgMLModelNames_0", std::vector<std::string>{"model.onnx"}, "ONNX file names for each pT bin (if not from CCDB full path)"},
      {"cfgMLModelNames_1", std::vector<std::string>{}, "ONNX file names for each pT bin (if not from CCDB full path)"}
    }};
    std::array<Configurable<float>, sMLModelsNum> fMatchingScoreCut{{
      {"cfgMLModelMatchingScoreCut_0", 1.f / 41.f, "Minimum score value for selecting good matches"},
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

  int mRunNumber{0};                               // needed to detect if the run changed and trigger update of magnetic field

  Service<o2::ccdb::BasicCCDBManager> ccdbManager;
  o2::ccdb::CcdbApi fCCDBApi;

  o2::aod::rctsel::RCTFlagsChecker rctChecker{"CBT_muon_glo", false, false, true};

  // vector of all MFT-MCH(-MID) matching candidates associated to the same MCH(-MID) track,
  // to be sorted in descending order with respect to the matching score
  // the map key is the MCH(-MID) track global index
  // the elements are pairs og global muon track indexes and associated matching scores
  // for matching candidates computed with the chi2 method, the score is defined as 1/(1+chi2)
  using MatchingCandidates = std::map<uint64_t, std::vector<std::pair<uint64_t, double>>>;

  struct CollisionInfo
  {
    uint64_t index{0};
    uint64_t bc{0};
    // z position of the collision
    double zVertex{0};
    // number of MFT tracks associated to the collision
    int mftTracksMultiplicity{0};
    // vector of MFT track indexes
    std::vector<uint64_t> mftTracks;
    // vector of MCH(-MID) track indexes
    std::vector<uint64_t> mchTracks;
    // matching candidates
    MatchingCandidates matchingCandidates;
    // vector of MFT-MCH track index pairs belonging to the same MC particle
    std::vector<std::pair<uint64_t, uint64_t>> matchablePairs;
    // vector of MCH track indexes that are expected to have an associated MFT track
    std::vector<uint64_t> taggedMuons;
  };

  using CollisionInfos = std::map<uint64_t, CollisionInfo>;

  std::unordered_map<int64_t, int32_t> mftTrackCovs;

  std::vector<std::pair<uint64_t, uint64_t>> fMatchablePairs;
  MatchingCandidates fMatchingCandidates;
  std::vector<uint64_t> fTaggedMuons;

  HistogramRegistry registry{"registry", {}};
  HistogramRegistry registryMatching{"registryMatching", {}};

  std::unordered_map<std::string, o2::framework::HistPtr> matchingHistos;

  struct  EfficiencyPlotter
  {
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
                      //std::unordered_map<std::string, o2::framework::HistPtr> histograms)
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
      //histograms[histName] = p_num;

      histName = path + "p_den";
      histTitle = title + " vs. p - den";
      p_den = registry.add(histName.c_str(), histTitle.c_str(), {HistType::kTH1F, {pAxis}});
      //histograms[histName] = p_den;

      // pT dependence
      histName = path + "pt_num";
      histTitle = title + " vs. p_{T} - num";
      pt_num = registry.add(histName.c_str(), histTitle.c_str(), {HistType::kTH1F, {pTAxis}});
      //histograms[histName] = pt_num;

      histName = path + "pt_den";
      histTitle = title + " vs. p_{T} - den";
      pt_den = registry.add(histName.c_str(), histTitle.c_str(), {HistType::kTH1F, {pTAxis}});
      //histograms[histName] = pt_den;

      // eta dependence
      histName = path + "eta_num";
      histTitle = title + " vs. #eta - num";
      eta_num = registry.add(histName.c_str(), histTitle.c_str(), {HistType::kTH1F, {etaAxis}});
      //histograms[histName] = eta_num;

      histName = path + "eta_den";
      histTitle = title + " vs. #eta - den";
      eta_den = registry.add(histName.c_str(), histTitle.c_str(), {HistType::kTH1F, {etaAxis}});
      //histograms[histName] = eta_den;

      // phi dependence
      histName = path + "phi_num";
      histTitle = title + " vs. #phi - num";
      phi_num = registry.add(histName.c_str(), histTitle.c_str(), {HistType::kTH1F, {phiAxis}});
      //histograms[histName] = phi_num;

      histName = path + "phi_den";
      histTitle = title + " vs. #phi - den";
      phi_den = registry.add(histName.c_str(), histTitle.c_str(), {HistType::kTH1F, {phiAxis}});
      //histograms[histName] = phi_den;
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

  struct MatchingPlotter
  {
    o2::framework::HistPtr fTrueMatchRanking;
    o2::framework::HistPtr fTrueMatchRankingVsP;
    o2::framework::HistPtr fTrueMatchRankingVsPt;
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

      AxisSpec indexAxis = {25, 0, 25, "ranking index"};
      std::string histName = path + "trueMatchRanking";
      std::string histTitle = "True match ranking";
      fTrueMatchRanking = registry.add(histName.c_str(), histTitle.c_str(), {HistType::kTH1F, {indexAxis}});
      histName = path + "trueMatchRankingVsP";
      histTitle = "True match ranking vs. p";
      fTrueMatchRankingVsP = registry.add(histName.c_str(), histTitle.c_str(), {HistType::kTH2F, {pAxis, indexAxis}});
      histName = path + "trueMatchRankingVsPt";
      histTitle = "True match ranking vs. p_{T}";
      fTrueMatchRankingVsPt = registry.add(histName.c_str(), histTitle.c_str(), {HistType::kTH2F, {ptAxis, indexAxis}});

      AxisSpec scoreAxis = {1000, 0, 1, "matching score"};
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
      //std::cout << "fieldB: " << (void*)fieldB << std::endl;
    }
  }

  void createMatchingHistosMC()
  {
    AxisSpec indexAxis = {25, 0, 25, "index"};
    AxisSpec chi2Axis = {1000, 0, 1000, "chi^{2}"};
    AxisSpec pAxis = {1000, 0, 100, "p (GeV/c)"};
    AxisSpec pTAxis = {100, 0, 10, "p_{T} (GeV/c)"};
    AxisSpec etaAxis = {100, -4, -2, "#eta"};
    AxisSpec phiAxis = {90, -180, 180, "#phi (degrees)"};
    std::string histPath = "matching/MC/";

    std::string histName;
/*
    histName = "candidateIdTrueMatch";
    matchingHistos[histName] = registryMatching.add((histPath + histName).c_str(), "candidate index - true matches", {HistType::kTH1F, {indexAxis}});

    histName = "chi2_true";
    matchingHistos[histName] = registryMatching.add((histPath + histName).c_str(), "chi^{2} - true matches", {HistType::kTH1F, {chi2Axis}});
    histName = "chi2VsP_true";
    matchingHistos[histName] = registryMatching.add((histPath + histName).c_str(), "chi^{2} vs. MCH momentum - true matches", {HistType::kTH2F, {pAxis, chi2Axis}});

    histName = "chi2_fake";
    matchingHistos[histName] = registryMatching.add((histPath + histName).c_str(), "chi^{2} - fake matches", {HistType::kTH1F, {chi2Axis}});
    histName = "chi2VsP_fake";
    matchingHistos[histName] = registryMatching.add((histPath + histName).c_str(), "chi^{2} vs. MCH momentum - fake matches", {HistType::kTH2F, {pAxis, chi2Axis}});

    fMatchingPurityPlotter = std::make_unique<EfficiencyPlotter>(histPath + "matching-purity/", "Matching purity", registryMatching);
    fPairingEfficiencyPlotter = std::make_unique<EfficiencyPlotter>(histPath + "pairing-efficiency/", "Pairing efficiency", registryMatching);
    fMatchingEfficiencyPlotter = std::make_unique<EfficiencyPlotter>(histPath + "matching-efficiency/", "Matching efficiency", registryMatching);
    fFakeMatchingEfficiencyPlotter = std::make_unique<EfficiencyPlotter>(histPath + "fake-matching-efficiency/", "Fake matching efficiency", registryMatching);
*/
    fChi2MatchingPlotter = std::make_unique<MatchingPlotter>(histPath + "Chi2/", registryMatching);
    for (const auto& [label, func] : matchingChi2Functions) {
      fMatchingPlotters[label] = std::make_unique<MatchingPlotter>(histPath + label + "/", registryMatching);
    }
    for (const auto& [label, response] : matchingMlResponses) {
      fMatchingPlotters[label] = std::make_unique<MatchingPlotter>(histPath + label + "/", registryMatching);
    }
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
    using SMatrix54 = ROOT::Math::SMatrix<double, 5, 4>;
    using SMatrix22 = ROOT::Math::SMatrix<double, 2>;
    using SMatrix25 = ROOT::Math::SMatrix<double, 2, 5>;
    using SMatrix52 = ROOT::Math::SMatrix<double, 5, 2>;

    // Define built-in matching functions
    //________________________________________________________________________________
    mMatchingFunctionMap["matchALL"] = [](const o2::dataformats::GlobalFwdTrack& mchTrack, const o2::track::TrackParCovFwd& mftTrack) -> double {
      // Match two tracks evaluating all parameters: X,Y, phi, tanl & q/pt

      //SMatrix55Sym I = ROOT::Math::SMatrixIdentity(), H_k, V_k;
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

      // Kalman Gain Matrix
      //SMatrix55Std K_k = GlobalMuonTrackCovariances * ROOT::Math::Transpose(H_k) * invResCov;

      // Update Parameters
      r_k_kminus1 = m_k - H_k * GlobalMuonTrackParameters; // Residuals of prediction

      auto matchChi2Track = ROOT::Math::Similarity(r_k_kminus1, invResCov);

      return matchChi2Track;
    };

    //________________________________________________________________________________
    mMatchingFunctionMap["matchXYPhiTanl"] = [](const o2::dataformats::GlobalFwdTrack& mchTrack, const o2::track::TrackParCovFwd& mftTrack) -> double {

    // Match two tracks evaluating positions & angles

    //SMatrix55Sym I = ROOT::Math::SMatrixIdentity();
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

    // Kalman Gain Matrix
    //SMatrix54 K_k = GlobalMuonTrackCovariances * ROOT::Math::Transpose(H_k) * invResCov;

    // Residuals of prediction
    r_k_kminus1 = m_k - H_k * GlobalMuonTrackParameters;

    auto matchChi2Track = ROOT::Math::Similarity(r_k_kminus1, invResCov);

    return matchChi2Track; };

    //________________________________________________________________________________
    mMatchingFunctionMap["matchXY"] = [](const o2::dataformats::GlobalFwdTrack& mchTrack, const o2::track::TrackParCovFwd& mftTrack) -> double {

    // Calculate Matching Chi2 - X and Y positions

    //SMatrix55Sym I = ROOT::Math::SMatrixIdentity();
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

    // Kalman Gain Matrix
    //SMatrix52 K_k = GlobalMuonTrackCovariances * ROOT::Math::Transpose(H_k) * invResCov;

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

      if (label == "" || funcName == "") break;

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

      if (label == "" || modelPaths.empty() || inputFeatures.empty() || modelNames.empty()) break;

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
  }

  template<class T, class C>
  bool pDCACut(const T& mchTrack, const C& collision, double nSigmaPDCA)
  {
    static const double sigmaPDCA23 = 80.;
    static const double sigmaPDCA310 = 54.;
    static const double relPRes = 0.0004;
    static const double slopeRes = 0.0005;

    double thetaAbs = TMath::ATan(mchTrack.rAtAbsorberEnd() / 505.) * TMath::RadToDeg();

    // propagate muon track to vertex
    auto mchTrackAtVertex = VarManager::PropagateMuon(mchTrack, collision, VarManager::kToVertex);

    //double pUncorr = mchTrack.p();
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

  template<class T, class C>
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

  template<class T, class C>
  bool IsGoodMuon(const T& muonTrack, const C& collision)
  {
    return IsGoodMuon(muonTrack, collision, fTrackChi2MchUp, fPMchLow, fPtMchLow, {fEtaMchLow, fEtaMchUp}, {fRabsLow, fRabsUp}, fSigmaPdcaUp);
  }

  template<class T, class C>
  bool IsGoodGlobalMuon(const T& muonTrack, const C& collision)
  {
    return IsGoodMuon(muonTrack, collision, fTrackChi2MchUp, fPMchLow, fPtMchLow, {fEtaMftLow, fEtaMftUp}, {fRabsLow, fRabsUp}, fSigmaPdcaUp);
  }

  template<class T>
  bool IsGoodMFT(const T& mftTrack,
                 double chi2Cut,
                 int nClustersCut)
  {
    // chi2 cut
    if (mftTrack.chi2() > chi2Cut)
      return false;

    // number of clusters cut
    if (mftTrack.nClusters() < nClustersCut)
      return false;

    return true;
  }

  template<class T>
  bool IsGoodMFT(const T& mftTrack)
  {
    return IsGoodMFT(mftTrack, fTrackChi2MftUp, fTrackNClustMftLow);
  }

  template<class TMUON>
  bool IsGoodGlobalMatching(const TMUON& muonTrack,
                            double matchingScore,
                            double matchingScoreCut)
  {
    if (static_cast<int>(muonTrack.trackType()) >= 2)
      return false;

    // MFT-MCH matching score cut
    if (muonTrack.chi2MatchMCHMFT() > matchingScoreCut)
      return false;

    return true;
  }

  template<class TMUON>
  bool IsGoodGlobalMatching(const TMUON& muonTrack, double matchingScore)
  {
    return IsGoodGlobalMatching(muonTrack, matchingScore, fMatchingChi2ScoreMftMchLow);
  }

  template <class TMUON>
  bool IsTrueGlobalMatching(const TMUON& muonTrack, const std::vector<std::pair<uint64_t, uint64_t>>& matchablePairs)
  {
    if (static_cast<int>(muonTrack.trackType()) > 2)
      return false;

    //auto const& mchTrack = muonTrack.template matchMCHTrack_as<TMUON>();
    //auto const& mftTrack = muonTrack.template matchMFTTrack_as<TMFT>();
    //uint64_t mchTrackId = static_cast<uint64_t>(mchTrack.globalIndex());
    uint64_t mchTrackId = static_cast<uint64_t>(muonTrack.matchMCHTrackId());
    uint64_t mftTrackId = static_cast<uint64_t>(muonTrack.matchMFTTrackId());

    std::pair<uint64_t, uint64_t> trackIndexes = std::make_pair(mchTrackId, mftTrackId);

    return (std::find(matchablePairs.begin(), matchablePairs.end(), trackIndexes) != matchablePairs.end());
  }

  bool IsMatchableMCH(uint64_t mchTrackId, const std::vector<std::pair<uint64_t, uint64_t>>& matchablePairs)
  {
    for (auto [id1, id2] : matchablePairs) {
      if (mchTrackId == id1) return true;
    }
    return false;
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

  template<class EVT, class BC, class TMUON, class TMFT>
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
      uint64_t collisionIndex = collision.globalIndex();

      auto bc = bcs.rawIteratorAt(collision.bcId());

      auto& collisionInfo = collisionInfos[collisionIndex];
      collisionInfo.index = collisionIndex;
      collisionInfo.bc = bc.globalBC();
      collisionInfo.zVertex = collision.posZ();

      if (static_cast<int>(muonTrack.trackType()) > 2) {
        // standalone MCH or MCH-MID tracks
        uint64_t mchTrackIndex = muonTrack.globalIndex();
        if (!IsGoodMuon(muonTrack, collision)) continue;
        collisionInfo.mchTracks.push_back(mchTrackIndex);
      } else {
        // global muon tracks (MFT-MCH or MFT-MCH-MID)
        uint64_t muonTrackIndex = muonTrack.globalIndex();
        double matchingScore = 1.f / (muonTrack.chi2MatchMCHMFT() + 1.f);
        auto const& mchTrack = muonTrack.template matchMCHTrack_as<TMUON>();
        auto const& mftTrack = muonTrack.template matchMFTTrack_as<TMFT>();

        uint64_t mchTrackIndex = mchTrack.globalIndex();

        // check if a vector of global muon candidates is already available for the current MCH index
        // if not, initialize a new one and add the current global muon track
        //bool globalMuonTrackFound = false;
        auto matchingCandidateIterator = collisionInfo.matchingCandidates.find(mchTrackIndex);
        if (matchingCandidateIterator != collisionInfo.matchingCandidates.end()) {
          matchingCandidateIterator->second.push_back(std::make_pair(muonTrackIndex, matchingScore));
          //globalMuonTrackFound = true;
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
      uint64_t collisionIndex = collision.globalIndex();

      auto bc = bcs.rawIteratorAt(collision.bcId());

      uint64_t mftTrackIndex = mftTrack.globalIndex();

      auto& collisionInfo = collisionInfos[collisionIndex];
      collisionInfo.bc = bc.globalBC();
      collisionInfo.zVertex = collision.posZ();

      collisionInfo.mftTracks.push_back(mftTrackIndex);
    }

    // sort the vectors of matching candidates in ascending order based on the matching score value
    auto compareMatchingScore = [](std::pair<uint64_t, double> track1, std::pair<uint64_t, double> track2) -> bool {
      return (track1.second > track2.second);
    };

    for (auto& [collisionIndex, collisionInfo] : collisionInfos) {
      for (auto& [mchIndex, globalTracksVector] : collisionInfo.matchingCandidates) {
        std::sort(globalTracksVector.begin(), globalTracksVector.end(), compareMatchingScore);
      }
    }
  }
/*
  template<class EVT, class BC, class TMUON, class TMFT>
  void InitCollisions(EVT const& collisions,
                      BC const& bcs,
                      TMUON const& muonTracks,
                      TMFT const& mftTracks,
                      std::map<uint64_t, CollisionInfo>& collisionInfos)
  {
    InitCollisions(collisions, bcs, muonTracks, collisionInfos);

    // fill collision information for MFT standalone tracks
    for (auto mftTrack : mftTracks) {
      if (!mftTrack.has_collision())
        continue;

      auto collision = collisions.rawIteratorAt(mftTrack.collisionId());
      uint64_t collisionIndex = collision.globalIndex();

      auto bc = bcs.rawIteratorAt(collision.bcId());

      uint64_t mftTrackIndex = mftTrack.globalIndex();

      auto& collisionInfo = collisionInfos[collisionIndex];
      collisionInfo.bc = bc.globalBC();
      collisionInfo.zVertex = collision.posZ();

      collisionInfo.mftTracks.push_back(mftTrackIndex);
    }
  }
*/

  template <class TMUON, class TMFT>
  void GetMatchablePairs(TMUON const& muonTracks,
                         TMFT const& mftTracks,
                         std::vector<std::pair<uint64_t, uint64_t>>& matchablePairs)
  {
    for (const auto& muonTrack : muonTracks) {
      // only consider MCH standalone or MCH-MID matches
      if (static_cast<int>(muonTrack.trackType()) <= 2) continue;
      // skip tracks that do not have an associated MC particle
      if (!muonTrack.has_mcParticle()) continue;
      // get the index associated to the MC particle
      auto muonMcParticle = muonTrack.mcParticle();
      uint64_t muonMcTrackIndex = muonMcParticle.globalIndex();

      for (const auto& mftTrack : mftTracks) {
        // skip tracks that do not have an associated MC particle
        if (!mftTrack.has_mcParticle()) continue;
        // get the index associated to the MC particle
        auto mftMcParticle = mftTrack.mcParticle();
        uint64_t mftMcTrackIndex = mftMcParticle.globalIndex();

        if (muonMcTrackIndex == mftMcTrackIndex) {
          matchablePairs.emplace_back(std::make_pair(static_cast<uint64_t>(muonTrack.globalIndex()),
                                                     static_cast<uint64_t>(mftTrack.globalIndex())));
        }
      }
    }
  }

  // for each MCH standalone track, collect the associated matching candidates
  template<class TMUON, class TMFT, class COMP>
  void GetMatchingCandidates(TMUON const& muonTracks,
                             TMFT const& mftTracks,
                             COMP comparator,
                             MatchingCandidates& matchingCandidates)
  {
    for (auto muonTrack : muonTracks) {

      if (static_cast<int>(muonTrack.trackType()) > 2) {
        // standalone MCH or MCH-MID tracks
        continue;
      }

      // global muon tracks (MFT-MCH or MFT-MCH-MID)
      uint64_t muonTrackIndex = muonTrack.globalIndex();
      double matchingScore = 1.f / (muonTrack.chi2MatchMCHMFT() + 1.f);
      auto const& mchTrack = muonTrack.template matchMCHTrack_as<TMUON>();
      uint64_t mchTrackIndex = mchTrack.globalIndex();

      // check if a vector of global muon candidates is already available for the current MCH index
      // if not, initialize a new one and add the current global muon track
      auto matchingCandidateIterator = matchingCandidates.find(mchTrackIndex);
      if (matchingCandidateIterator != matchingCandidates.end()) {
        matchingCandidateIterator->second.push_back(std::make_pair(muonTrackIndex, matchingScore));
      } else {
        matchingCandidates[mchTrackIndex].push_back(std::make_pair(muonTrackIndex, matchingScore));
      }
    }

    // sort the vectors of matching candidates in ascending order based on the matching score
    for (auto& [mchIndex, candidates] : matchingCandidates) {
      std::sort(candidates.begin(), candidates.end(),
          [](std::pair<uint64_t, double> track1, std::pair<uint64_t, double> track2) -> bool {
                return (track1.second > track2.second);
      });
    }

    /*// sort the vectors of matching candidates in ascending order based on the user-provided comparator function
    for (auto& [mchIndex, candidates] : matchingCandidates) {
      std::sort(candidates.begin(), candidates.end(),
          [&muonTracks, &comparator](uint64_t trackIndex1, uint64_t trackIndex2) -> bool {
                auto const& track1 = muonTracks.rawIteratorAt(trackIndex1);
                auto const& track2 = muonTracks.rawIteratorAt(trackIndex2);

                return comparator(track1, track2);
      });
    }*/
  }

  // for each MCH standalone track, collect the associated matching candidates
  template<class TMUON, class C>
  void GetTaggedMuons(TMUON const& muonTracks,
                      C const& collisions,
                      std::vector<uint64_t>& taggedMuons)
  {
    for (auto muonTrack : muonTracks) {

      // only consider MCH-MID matches
      if (static_cast<int>(muonTrack.trackType()) != 3) {
        continue;
      }

      //// Get collision information if associated
       if (!muonTrack.has_collision()) {
         continue;
       }
       const auto& collision = collisions.rawIteratorAt(muonTrack.collisionId());

      if (!IsGoodMuon(muonTrack, collision,
                      fMuonTaggingTrackChi2MchUp,
                      fMuonTaggingPMchLow,
                      fMuonTaggingPtMchLow,
                      {fMuonTaggingEtaMchLow, fMuonTaggingEtaMchUp},
                      {fMuonTaggingRabsLow, fMuonTaggingRabsUp},
                      fMuonTaggingSigmaPdcaUp)) {
        continue;
      }
      //std::cout << "[TOTO1] " << (int)muonTrack.globalIndex() << "  " << (int)muonTrack.trackType() << "  " << muonTrack.chi2MatchMCHMFT() << std::endl;

       uint64_t muonTrackIndex = muonTrack.globalIndex();
      //auto const& mchTrack = muonTrack.template matchMCHTrack_as<TMUON>();
      //uint64_t mchTrackIndex = mchTrack.globalIndex();
      //std::cout << "[TOTO1] MCH " << (int)mchTrack.globalIndex() << "  " << (int)mchTrack.trackType() << std::endl;

      // propagate MCH track to the vertex to get the updated momentum
      auto const& mchTrackAtVertex = VarManager::PropagateMuon(muonTrack, collision, VarManager::kToVertex);

      // propagate the extrapolated track to each of the MFT planes
      int nCrossedPlanes = 0;
      for (auto zPlaneMFT : o2::mft::constants::mft::LayerZCoordinate()) {
        const auto& extrapToMFT = PropagateToZMCH(mchTrackAtVertex, zPlaneMFT);
        if (true) {
          nCrossedPlanes += 1;
        }
      }

      // check that the track crosses enough MFT planes
      if (nCrossedPlanes < fMuonTaggingNCrossedMftPlanesLow) {
        continue;
      }

      taggedMuons.emplace_back(muonTrackIndex);
    }
  }

  template <class C, class TMUON, class TMFT>
  void FillMatchingPlotsMC(C const& collision,
                           TMUON const& muonTracks,
                           TMFT const& mftTracks,
                           const MatchingCandidates& matchingCandidates,
                           const std::vector<std::pair<uint64_t, uint64_t>>& matchablePairs,
                           double matchingScoreCut,
                           MatchingPlotter* plotter)
  {
    // ====================================
    // Matching candidates hierarchy

    for (auto [mchIndex, globalTracksVector] : matchingCandidates) {
      // only consider MCH tracks belonging to a matchable pair
      if (!IsMatchableMCH(static_cast<uint64_t>(mchIndex), matchablePairs)) continue;

      // get the standalone MCH track
      auto const& mchTrack = muonTracks.rawIteratorAt(mchIndex);

      // skip global muon tracks that do not pass the MCH quality cuts
      if (!IsGoodGlobalMuon(mchTrack, collision)) continue;

      // find the index of the matching candidate that corresponds to the true match
      // index=1 corresponds to the leading candidate
      // index=0 means no candidate was found that corresponds to the true match
      int trueMatchIndex = 0;
      for (int i = 0; i < globalTracksVector.size(); i++) {
        auto const& muonTrack = muonTracks.rawIteratorAt(globalTracksVector[i].first);
        double matchingScore = globalTracksVector[i].second;

        // get the standalone MFT track
        auto const& mftTrack = muonTrack.template matchMFTTrack_as<TMFT>();

        // only consider global muon tracks that pass the MFT and matching quality cuts
        if (!IsGoodMFT(mftTrack)) continue;
        if (!IsGoodGlobalMatching(muonTrack, matchingScore, matchingScoreCut)) continue;

        if (IsTrueGlobalMatching(muonTrack, matchablePairs)) {
          trueMatchIndex = i + 1;
          break;
        }
      }
      std::get<std::shared_ptr<TH1>>(plotter->fTrueMatchRanking)->Fill(trueMatchIndex);
    }

    // ====================================
    // Matching properties

    for (auto [mchIndex, globalTracksVector] : matchingCandidates) {
      if (globalTracksVector.size() < 1) continue;

      auto const& mchTrack = muonTracks.rawIteratorAt(mchIndex);

      // get leading matching candidate
      auto const& muonTrack = muonTracks.rawIteratorAt(globalTracksVector[0].first);

      auto const& mftTrack = muonTrack.template matchMFTTrack_as<TMFT>();

      // skip global muon tracks that do not pass the MCH and MFT quality cuts
      if (!IsGoodGlobalMuon(muonTrack, collision)) continue;
      if (!IsGoodMFT(mftTrack)) continue;

      double matchingScore = globalTracksVector[0].second;
      double mchMom = mchTrack.p();
      double mchPt = mchTrack.pt();

      // check the matching quality, but set the minimum matching score to zero
      bool isGoodMatchNoScore = IsGoodGlobalMatching(muonTrack, matchingScore, 0);
      //bool isGoodMatch = IsGoodGlobalMatching(muonTrack, matchingScore, matchingScoreCut);

      bool isTrueMatch = IsTrueGlobalMatching(muonTrack, matchablePairs);

      // matching score analysis
      if (isGoodMatchNoScore) {
        if (isTrueMatch) {
          std::get<std::shared_ptr<TH1>>(plotter->fTrueMatchScore)->Fill(matchingScore);
        } else {
          std::get<std::shared_ptr<TH1>>(plotter->fFakeMatchScore)->Fill(matchingScore);
        }
      }
    }

    // ====================================
    // Matching purity

    for (auto [mchIndex, globalTracksVector] : matchingCandidates) {
      if (globalTracksVector.size() < 1) continue;

      // get the leading matching candidate
      auto const& muonTrack = muonTracks.rawIteratorAt(globalTracksVector[0].first);
      double matchingScore = globalTracksVector[0].second;

      // get the standalone MCH and MFT tracks
      auto const& mchTrack = muonTracks.rawIteratorAt(mchIndex);
      auto const& mftTrack = muonTrack.template matchMFTTrack_as<TMFT>();

      // skip global muon tracks that do not pass the MCH and MFT quality cuts
      if (!IsGoodGlobalMuon(muonTrack, collision)) continue;
      if (!IsGoodMFT(mftTrack)) continue;

      // skip  candidates that do not pass the matching quality cuts
      if (!IsGoodGlobalMatching(muonTrack, matchingScore, matchingScoreCut)) continue;

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
      auto const& mchTrack = muonTracks.rawIteratorAt(matchableMchIndex);
      auto const& pairedMftTrack = muonTracks.rawIteratorAt(matchableMftIndex);

      // skip  track pairs that do not pass the MCH and MFT quality cuts
      // we only consider matchable pairs that fulfill the track quality requirements
      if (!IsGoodGlobalMuon(mchTrack, collision)) continue;
      if (!IsGoodMFT(pairedMftTrack)) continue;

      bool goodMatchFound = false;
      bool isTrueMatch = false;

      // check if we have some matching candidates for the current matchable MCH track
      if (matchingCandidates.count(matchableMchIndex) > 0) {
        const auto& globalTracksVector = matchingCandidates.at(static_cast<uint64_t>(matchableMchIndex));
        if (!globalTracksVector.empty()) {
          // get the leading matching candidate
          auto const& muonTrack = muonTracks.rawIteratorAt(globalTracksVector[0].first);
          double matchingScore = globalTracksVector[0].second;

          // get the standalone MFT track
          auto const& mftTrack = muonTrack.template matchMFTTrack_as<TMFT>();
          auto mftIndex = mftTrack.globalIndex();

          // a good match must pass the MFT and matching quality cuts
          // the MCH track quality is already checked in the outer loop
          goodMatchFound = IsGoodMFT(mftTrack) && IsGoodGlobalMatching(muonTrack, matchingScore, matchingScoreCut);
          isTrueMatch = (mftIndex == matchableMftIndex);
        }
      }

      // fill matching efficiency plots
      plotter->fPairingEfficiencyPlotter.Fill(mchTrack, goodMatchFound);
      plotter->fMatchingEfficiencyPlotter.Fill(mchTrack, isTrueMatch);
      plotter->fFakeMatchingEfficiencyPlotter.Fill(mchTrack, !isTrueMatch);
    }
  }

  template <class C, class TMUON, class TMFT, class CMFT>
  void RunChi2Matching(C const& collisions,
                       TMUON const& muonTracks,
                       TMFT const& mftTracks,
                       CMFT const& mftCovs,
                       std::string label,
                       MatchingCandidates matchingCandidates)
  {
    matchingCandidates.clear();
    auto funcIter = matchingChi2Functions.find(label);
    if (funcIter == matchingChi2Functions.end()) return;

    auto funcName = funcIter->second;
    if (mMatchingFunctionMap.count(funcName) < 1) return;
    auto matchingFunc = mMatchingFunctionMap.at(funcName);
    for (auto& [mchIndex, globalTracksVector] : fMatchingCandidates) {
      auto const& mchTrack = muonTracks.rawIteratorAt(mchIndex);

      for (auto& [muonIndex, score] : globalTracksVector) {
        auto const& muonTrack = muonTracks.rawIteratorAt(muonIndex);
        if (!muonTrack.has_collision()) continue;

        auto collision = collisions.rawIteratorAt(muonTrack.collisionId());

        // get MCH and MFT standalone tracks
        //auto mchTrack = muonTrack.template matchMCHTrack_as<TMUON>();
        auto const& mftTrack = muonTrack.template matchMFTTrack_as<TMFT>();
        if (mftTrackCovs.count(mftTrack.globalIndex()) < 1) {
          std::cout << std::format("Covariance matrix for MFT track #{} not found", mftTrack.globalIndex()) << std::endl;
          continue;
        }
        auto const& mftTrackCov = mftCovs.rawIteratorAt(mftTrackCovs[mftTrack.globalIndex()]);

        //continue;
        // get tracks parameters in O2 format
        auto mftprop = FwdToTrackPar(mftTrack, mftTrackCov);
        auto muonprop = FwdToTrackPar(mchTrack, mchTrack);

        // extrapolate to the matching plane
        auto matchingPlaneZ = matchingPlanesZ[label];
        if (matchingPlaneZ < 0.) {
          //mftprop = VarManager::PropagateFwd(mftTrack, mftTrackCov, fConfigVariousOptions.fzMatching.value);
          //muonprop = VarManager::PropagateMuon(mchTrack, collision, VarManager::kToMatching);
        }

        // run the chi2 matching function
        float matchingChi2 = matchingFunc(muonprop, mftprop);
        float matchingScore = 1.f / (matchingChi2 + 1.f);
        std::cout << std::format("Matching chi2: {}, original chi2: {}", matchingChi2, muonTrack.chi2MatchMCHMFT()) << std::endl;

        // check if a vector of global muon candidates is already available for the current MCH index
        // if not, initialize a new one and add the current global muon track
        auto matchingCandidateIterator = matchingCandidates.find(mchIndex);
        if (matchingCandidateIterator != matchingCandidates.end()) {
          matchingCandidateIterator->second.push_back(std::make_pair(muonIndex, matchingScore));
        } else {
          matchingCandidates[mchIndex].push_back(std::make_pair(muonIndex, matchingScore));
        }
      }
    }
  }

  template <class C, class TMUON, class TMFT, class CMFT>
  void RunMLMatching(C const& collisions,
                     TMUON const& muonTracks,
                     TMFT const& mftTracks,
                     CMFT const& mftCovs,
                     std::string label,
                     MatchingCandidates matchingCandidates)
  {
    matchingCandidates.clear();
    auto mlIter = matchingMlResponses.find(label);
    if (mlIter == matchingMlResponses.end()) return;

    auto& mlResponse = mlIter->second;
    for (auto& [mchIndex, globalTracksVector] : fMatchingCandidates) {
      auto const& mchTrack = muonTracks.rawIteratorAt(mchIndex);
      for (auto& [muonIndex, score] : globalTracksVector) {
        auto const& muonTrack = muonTracks.rawIteratorAt(muonIndex);
        if (!muonTrack.has_collision()) continue;

        auto collision = collisions.rawIteratorAt(muonTrack.collisionId());

        // get MCH and MFT standalone tracks
        //auto mchTrack = muonTrack.template matchMCHTrack_as<TMUON>();
        auto const& mftTrack = muonTrack.template matchMFTTrack_as<TMFT>();
        if (mftTrackCovs.count(mftTrack.globalIndex()) < 1) {
          std::cout << std::format("Covariance matrix for MFT track #{} not found", mftTrack.globalIndex()) << std::endl;
          continue;
        }
        std::cout << fmt::format("Getting covariance matrix for MFT track #{} -> {}", mftTrack.globalIndex(), mftTrackCovs[mftTrack.globalIndex()]) << std::endl;
        auto const& mftTrackCov = mftCovs.rawIteratorAt(mftTrackCovs[mftTrack.globalIndex()]);
        std::cout << fmt::format("Covariance matrix for MFT track #{} retrieved", mftTrack.globalIndex()) << std::endl;

        //continue;
        // get tracks parameters in O2 format
        auto mftprop = FwdToTrackPar(mftTrack, mftTrackCov);
        auto muonprop = FwdToTrackPar(mchTrack, mchTrack);

        // extrapolate to the matching plane
        auto matchingPlaneZ = matchingPlanesZ[label];
        if (matchingPlaneZ < 0.) {
          //mftprop = VarManager::PropagateFwd(mftTrack, mftTrackCov, fConfigVariousOptions.fzMatching.value);
          //muonprop = VarManager::PropagateMuon(mchTrack, collision, VarManager::kToMatching);
        }

        // run the ML model
        std::vector<float> output;
        std::vector<float> inputML = mlResponse.getInputFeaturesGlob(muonTrack, muonprop, mftprop, collision);
        mlResponse.isSelectedMl(inputML, 0, output);
        float matchingScore = output[0];
        std::cout << std::format("Matching score: {}, Chi2: {}", matchingScore, muonTrack.chi2MatchMCHMFT()) << std::endl;

        // check if a vector of global muon candidates is already available for the current MCH index
        // if not, initialize a new one and add the current global muon track
        auto matchingCandidateIterator = matchingCandidates.find(mchIndex);
        if (matchingCandidateIterator != matchingCandidates.end()) {
          matchingCandidateIterator->second.push_back(std::make_pair(muonIndex, matchingScore));
        } else {
          matchingCandidates[mchIndex].push_back(std::make_pair(muonIndex, matchingScore));
        }
      }
    }
  }

  template <class C, class TMUON, class TMFT, class CMFT>
  void ProcessMatching(const CollisionInfo& collisionInfo,
                       C const& collisions,
                       TMUON const& muonTracks,
                       TMFT const& mftTracks,
                       CMFT const& mftCovs)
  {
    auto collision = collisions.rawIteratorAt(collisionInfo.index);
    FillMatchingPlotsMC(collision, muonTracks, mftTracks, collisionInfo.matchingCandidates, fMatchablePairs, fMatchingChi2ScoreMftMchLow, fChi2MatchingPlotter.get());

    for (auto& [label, func] : matchingChi2Functions) {
      MatchingCandidates matchingCandidates;
      RunChi2Matching(collisions, muonTracks, mftTracks, mftCovs, label, matchingCandidates);

      auto* plotter = fMatchingPlotters.at(label).get();
      double matchingScoreCut = matchingScoreCuts.at(label);

      FillMatchingPlotsMC(collision, muonTracks, mftTracks, matchingCandidates, fMatchablePairs, matchingScoreCut, plotter);
    }

    for (auto& [label, mlResponse] : matchingMlResponses) {
      MatchingCandidates matchingCandidates;
      RunMLMatching(collisions, muonTracks, mftTracks, mftCovs, label, matchingCandidates);

      auto* plotter = fMatchingPlotters.at(label).get();
      double matchingScoreCut = matchingScoreCuts.at(label);

      FillMatchingPlotsMC(collision, muonTracks, mftTracks, matchingCandidates, fMatchablePairs, matchingScoreCut, plotter);
    }

/*
    for (auto& [label, mlResponse] : matchingMlResponses) {
      MatchingCandidates candidates;
      for (auto& [mchIndex, globalTracksVector] : fMatchingCandidates) {
        for (auto& [muonIndex, score] : globalTracksVector) {
          auto const& muonTrack = muonTracks.rawIteratorAt(muonIndex);
        //auto muonID = muon.matchMCHTrackId();
        auto mchTrack = muonTrack.template matchMCHTrack_as<TMuons>();
        auto mftTrack = muonTrack.template matchMFTTrack_as<TMFTTracks>();
        auto const& mftTrackCov = mfCovs.rawIteratorAt(map_mfttrackcovs[mfttrack.globalIndex()]);
        o2::track::TrackParCovFwd mftprop = VarManager::FwdToTrackPar(mfttrack, mfttrackcov);
        o2::dataformats::GlobalFwdTrack muonprop = VarManager::FwdToTrackPar(muontrack, muontrack);
        auto matchingPlaneZ = matchingPlanesZ[label];
        if (matchingPlaneZ < 0.) {
          mftprop = VarManager::PropagateFwd(mfttrack, mfttrackcov, fConfigVariousOptions.fzMatching.value);
          muonprop = VarManager::PropagateMuon(muontrack, collision, VarManager::kToMatching);
        }
        std::vector<float> output;
        std::vector<float> inputML = matchingMlResponse.getInputFeaturesGlob(muon, muonprop, mftprop, collision);
        matchingMlResponse.isSelectedMl(inputML, 0, output);
        float score = output[0];
      }
    }*/
  }

  void processQAMC(MyEvents const& collisions,
      aod::BCsWithTimestamps const& bcs,
      MyMuonsMC const& muonTracks,
      MyMFTsMC const& mftTracks,
      MyMFTCovariances const& mftCovs,
      aod::McParticles const& mcParticles)
  {
    auto bc = bcs.begin();
    initCCDB(bc);

    InitCollisions(collisions, bcs, muonTracks, mftTracks, fCollisionInfos);

    mftTrackCovs.clear();
    for (auto& mftTrackCov : mftCovs) {
      mftTrackCovs[mftTrackCov.matchMFTTrackId()] = mftTrackCov.globalIndex();
    }

    auto comparatorChi2 = [&](const MyMuonMC& track1, const MyMuonMC& track2) -> bool {
      return (track1.chi2MatchMCHMFT() < track2.chi2MatchMCHMFT());
    };

    fMatchablePairs.clear();
    GetMatchablePairs(muonTracks, mftTracks, fMatchablePairs);

    fMatchingCandidates.clear();
    GetMatchingCandidates(muonTracks, mftTracks, comparatorChi2, fMatchingCandidates);

    fTaggedMuons.clear();
    GetTaggedMuons(muonTracks, collisions, fTaggedMuons);

    for (auto const& [collisionIndex, collisionInfo] : fCollisionInfos) {
    ProcessMatching(collisionInfo, collisions, muonTracks, mftTracks, mftCovs);
    }
    //FillMatchingPlotsMC(muonTracks, mftTracks, fMatchingCandidates, fMatchablePairs, fChi2MatchingPlotter.get());
  }

  PROCESS_SWITCH(qaMatching, processQAMC, "process qa MC", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<qaMatching>(cfgc)};
};
