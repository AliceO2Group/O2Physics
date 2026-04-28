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

#include "Common/CCDB/RCTSelectionFlags.h"
#include "Common/DataModel/EventSelection.h"
#include "Tools/ML/MlResponse.h"

#include <CCDB/BasicCCDBManager.h>
#include <CCDB/CcdbApi.h>
#include <CommonConstants/LHCConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DetectorsBase/GeometryManager.h>
#include <DetectorsBase/Propagator.h>
#include <Field/MagneticField.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Array2D.h>
#include <Framework/Configurable.h>
#include <Framework/DataTypes.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>
#include <GlobalTracking/MatchGlobalFwd.h>
#include <MCHTracking/TrackExtrap.h>
#include <MCHTracking/TrackParam.h>
#include <MFTTracking/Constants.h>
#include <ReconstructionDataFormats/GlobalFwdTrack.h>
#include <ReconstructionDataFormats/TrackFwd.h>

#include <Math/MatrixFunctions.h>
#include <Math/MatrixRepresentationsStatic.h>
#include <Math/ProbFuncMathCore.h>
#include <Math/SMatrix.h>
#include <Math/SVector.h>
#include <Math/Vector4D.h> // IWYU pragma: keep (do not replace with Math/Vector4Dfwd.h)
#include <Math/Vector4Dfwd.h>
#include <TGeoGlobalMagField.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TMath.h>

#include <RtypesCore.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <format>
#include <functional>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include <math.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::aod;

namespace qamatching
{
DECLARE_SOA_COLUMN(P, p, float);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(MatchLabel, matchLabel, int8_t);
DECLARE_SOA_COLUMN(TrackId, trackId, int64_t);
DECLARE_SOA_COLUMN(MatchType, matchType, int8_t);
DECLARE_SOA_COLUMN(MatchScore, matchScore, float);
DECLARE_SOA_COLUMN(MatchRanking, matchRanking, int32_t);
DECLARE_SOA_COLUMN(MftMultiplicity, mftMultiplicity, int32_t);
DECLARE_SOA_COLUMN(TrackType, trackType, int8_t);
DECLARE_SOA_COLUMN(MftMatchAttempts, mftMatchAttempts, int32_t);
DECLARE_SOA_COLUMN(XAtVtx, xAtVtx, float);
DECLARE_SOA_COLUMN(YAtVtx, yAtVtx, float);
DECLARE_SOA_COLUMN(ZAtVtx, zAtVtx, float);
DECLARE_SOA_COLUMN(PxAtVtx, pxAtVtx, float);
DECLARE_SOA_COLUMN(PyAtVtx, pyAtVtx, float);
DECLARE_SOA_COLUMN(PzAtVtx, pzAtVtx, float);
DECLARE_SOA_COLUMN(ColX, colX, float);
DECLARE_SOA_COLUMN(ColY, colY, float);
DECLARE_SOA_COLUMN(ColZ, colZ, float);
} // namespace qamatching

namespace o2::aod
{
DECLARE_SOA_TABLE(QaMatchingEvents, "AOD", "QAMEVT",
                  o2::soa::Index<>,
                  qamatching::MftMultiplicity,
                  qamatching::ColX,
                  qamatching::ColY,
                  qamatching::ColZ);
} // namespace o2::aod

namespace qamatching
{
DECLARE_SOA_INDEX_COLUMN_FULL(ReducedEvent, reducedEvent, int32_t, o2::aod::QaMatchingEvents, "");
} // namespace qamatching

namespace o2::aod
{
DECLARE_SOA_TABLE(QaMatchingMCHTrack, "AOD", "QAMCHTRK",
                  qamatching::ReducedEventId,
                  qamatching::TrackId,
                  qamatching::TrackType,
                  qamatching::P,
                  qamatching::Pt,
                  qamatching::Eta,
                  qamatching::Phi,
                  qamatching::MftMatchAttempts,
                  qamatching::XAtVtx,
                  qamatching::YAtVtx,
                  qamatching::ZAtVtx,
                  qamatching::PxAtVtx,
                  qamatching::PyAtVtx,
                  qamatching::PzAtVtx);
DECLARE_SOA_TABLE(QaMatchingCandidates, "AOD", "QAMCAND",
                  qamatching::ReducedEventId,
                  qamatching::MatchLabel,
                  qamatching::TrackId,
                  qamatching::P, qamatching::Pt, qamatching::Eta, qamatching::Phi,
                  qamatching::MatchType, qamatching::MatchScore, qamatching::MatchRanking,
                  qamatching::XAtVtx,
                  qamatching::YAtVtx,
                  qamatching::ZAtVtx,
                  qamatching::PxAtVtx,
                  qamatching::PyAtVtx,
                  qamatching::PzAtVtx);
} // namespace o2::aod

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
using SMatrix5 = ROOT::Math::SVector<double, 5>;

static float chi2ToScore(float chi2, int ndf, float chi2max)
{
  double p = -std::log10(ROOT::Math::chisquared_cdf_c(chi2, ndf));
  double pnorm = -std::log10(ROOT::Math::chisquared_cdf_c(chi2max, ndf));
  double result = (1.f / (p / pnorm + 1.f));
  return static_cast<float>(result);
}

struct QaMatching {

  template <class T, int nr, int nc>
  using Matrix = std::array<std::array<T, nc>, nr>;

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

  static constexpr int GlobalTrackTypeMax = 2;
  static constexpr int MchMidTrackType = 3;
  static constexpr int FirstDecayMotherRank = 2;
  static constexpr int MftTrackTypeStandard = 0;
  static constexpr int MftTrackTypeCA = 1;
  static constexpr int ThetaAbsBoundaryDeg = 3;
  static constexpr double SlopeResolutionZ = 535.;
  static constexpr int MatchingDegreesOfFreedom = 5;
  static constexpr size_t MinCandidatesForDeltaChi2 = 2;
  static constexpr float MatchingScoreChi2Max = 50.f;
  static constexpr float InvalidDeltaChi2 = -1.f;
  static constexpr int ExtrapolationMethodStandard = 0;
  static constexpr int ExtrapolationMethodMftFirstPoint = 2;
  static constexpr int ExtrapolationMethodVertex = 3;
  static constexpr int ExtrapolationMethodMftDca = 4;
  static constexpr int DecayRankingDirect = 2;

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
    int mftMchMatchAttempts{0};
    MuonMatchType matchType{kMatchTypeUndefined};
  };

  ////   Variables for selecting muon tracks
  Configurable<float> cfgPMchLow{"cfgPMchLow", 0.0f, ""};
  Configurable<float> cfgPtMchLow{"cfgPtMchLow", 0.7f, ""};
  Configurable<float> cfgEtaMchLow{"cfgEtaMchLow", -4.0f, ""};
  Configurable<float> cfgEtaMchUp{"cfgEtaMchUp", -2.5f, ""};
  Configurable<float> cfgRabsLow{"cfgRabsLow", 17.6f, ""};
  Configurable<float> cfgRabsUp{"cfgRabsUp", 89.5f, ""};
  Configurable<float> cfgPdcaUp{"cfgPdcaUp", 6.f, ""};
  Configurable<float> cfgTrackChi2MchUp{"cfgTrackChi2MchUp", 5.f, ""};
  Configurable<float> cfgMatchingChi2MchMidUp{"cfgMatchingChi2MchMidUp", 999.f, ""};

  ////   Variables for selecting mft tracks
  Configurable<float> cfgEtaMftLow{"cfgEtaMftLow", -3.6f, ""};
  Configurable<float> cfgEtaMftUp{"cfgEtaMftUp", -2.5f, ""};
  Configurable<int> cfgTrackNClustMftLow{"cfgTrackNClustMftLow", 7, ""};
  Configurable<float> cfgTrackChi2MftUp{"cfgTrackChi2MftUp", 999.f, ""};

  ////   Variables for selecting global tracks
  Configurable<float> cfgMatchingChi2ScoreMftMchLow{"cfgMatchingChi2ScoreMftMchLow", chi2ToScore(50.f, 5, 50.f), ""};

  ////   Variables for selecting tagged muons
  Configurable<int> cfgMuonTaggingNCrossedMftPlanesLow{"cfgMuonTaggingNCrossedMftPlanesLow", 5, ""};
  Configurable<float> cfgMuonTaggingTrackChi2MchUp{"cfgMuonTaggingTrackChi2MchUp", 5.f, ""};
  Configurable<float> cfgMuonTaggingPMchLow{"cfgMuonTaggingPMchLow", 0.0f, ""};
  Configurable<float> cfgMuonTaggingPtMchLow{"cfgMuonTaggingPtMchLow", 0.7f, ""};
  Configurable<float> cfgMuonTaggingEtaMchLow{"cfgMuonTaggingEtaMchLow", -3.6f, ""};
  Configurable<float> cfgMuonTaggingEtaMchUp{"cfgMuonTaggingEtaMchUp", -2.5f, ""};
  Configurable<float> cfgMuonTaggingRabsLow{"cfgMuonTaggingRabsLow", 17.6f, ""};
  Configurable<float> cfgMuonTaggingRabsUp{"cfgMuonTaggingRabsUp", 89.5f, ""};
  Configurable<float> cfgMuonTaggingPdcaUp{"cfgMuonTaggingPdcaUp", 4.f, ""};
  Configurable<float> cfgMuonTaggingChi2DiffLow{"cfgMuonTaggingChi2DiffLow", 100.f, ""};

  ////   Variables for ccdb
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpMagPath{"grpMagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};

  // CCDB connection configurables
  struct : ConfigurableGroup {
    Configurable<std::string> cfgCcdbUrl{"cfgCcdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
    Configurable<int64_t> cfgCcdbNoLaterThan{"cfgCcdbNoLaterThan", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
    Configurable<std::string> cfgGrpPath{"cfgGrpPath", "GLO/GRP/GRP", "Path of the grp file"};
    Configurable<std::string> cfgGeoPath{"cfgGeoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
    Configurable<std::string> cfgGrpmagPath{"cfgGrpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  } configCcdb;

  struct : ConfigurableGroup {
    Configurable<bool> cfgCreatePdgMomHistograms{"cfgCreatePdgMomHistograms", false, "create matching characteristics plots with particle mom PDG codes"};
  } configQas;

  ///    Variables for histograms configuration
  Configurable<int> cfgNCandidatesMax{"cfgNCandidatesMax", 5, "Number of matching candidates stored for each muon track"};
  Configurable<int> cfgMftTrackMultiplicityMax{"cfgMftTrackMultiplicityMax", 1000, "Maximum number of MFT tracks per collision"};

  Configurable<int> cfgQaMatchingAodDebug{"cfgQaMatchingAodDebug", 0, "If >0, print AO2D filling debug (0=off, N=max collisions)"};

  double mBzAtMftCenter{0};

  o2::globaltracking::MatchGlobalFwd mExtrap;

  using MatchingFunc = std::function<std::tuple<double, int>(const o2::dataformats::GlobalFwdTrack& mchtrack, const o2::track::TrackParCovFwd& mfttrack)>;
  std::map<std::string, MatchingFunc> mMatchingFunctionMap; ///< MFT-MCH Matching function

  // Chi2 matching interface
  static constexpr int Chi2FunctionsNum = 5;
  struct : ConfigurableGroup {
    Configurable<std::string> cfgChi2FunctionLabel1{"cfgChi2FunctionLabel1", std::string{"ProdAll"}, "Text label identifying this chi2 matching method"};
    Configurable<std::string> cfgChi2FunctionLabel2{"cfgChi2FunctionLabel2", std::string{"MatchXYPhiTanlMom"}, "Text label identifying this chi2 matching method"};
    Configurable<std::string> cfgChi2FunctionLabel3{"cfgChi2FunctionLabel3", std::string{"MatchXYPhiTanl"}, "Text label identifying this chi2 matching method"};
    Configurable<std::string> cfgChi2FunctionLabel4{"cfgChi2FunctionLabel4", std::string{""}, "Text label identifying this chi2 matching method"};
    Configurable<std::string> cfgChi2FunctionLabel5{"cfgChi2FunctionLabel5", std::string{""}, "Text label identifying this chi2 matching method"};
    std::array<Configurable<std::string>*, Chi2FunctionsNum> functionLabels{
      &cfgChi2FunctionLabel1, &cfgChi2FunctionLabel2, &cfgChi2FunctionLabel3, &cfgChi2FunctionLabel4, &cfgChi2FunctionLabel5};

    Configurable<std::string> cfgChi2FunctionName1{"cfgChi2FunctionName1", std::string{"prod"}, "Name of the chi2 matching function"};
    Configurable<std::string> cfgChi2FunctionName2{"cfgChi2FunctionName2", std::string{"matchALL"}, "Name of the chi2 matching function"};
    Configurable<std::string> cfgChi2FunctionName3{"cfgChi2FunctionName3", std::string{"matchXYPhiTanl"}, "Name of the chi2 matching function"};
    Configurable<std::string> cfgChi2FunctionName4{"cfgChi2FunctionName4", std::string{""}, "Name of the chi2 matching function"};
    Configurable<std::string> cfgChi2FunctionName5{"cfgChi2FunctionName5", std::string{""}, "Name of the chi2 matching function"};
    std::array<Configurable<std::string>*, Chi2FunctionsNum> functionNames{
      &cfgChi2FunctionName1, &cfgChi2FunctionName2, &cfgChi2FunctionName3, &cfgChi2FunctionName4, &cfgChi2FunctionName5};

    Configurable<float> cfgChi2FunctionMatchingScoreCut1{"cfgChi2FunctionMatchingScoreCut1", 0.f, "Minimum score value for selecting good matches"};
    Configurable<float> cfgChi2FunctionMatchingScoreCut2{"cfgChi2FunctionMatchingScoreCut2", 0.5f, "Minimum score value for selecting good matches"};
    Configurable<float> cfgChi2FunctionMatchingScoreCut3{"cfgChi2FunctionMatchingScoreCut3", 0.5f, "Minimum score value for selecting good matches"};
    Configurable<float> cfgChi2FunctionMatchingScoreCut4{"cfgChi2FunctionMatchingScoreCut4", 0.5f, "Minimum score value for selecting good matches"};
    Configurable<float> cfgChi2FunctionMatchingScoreCut5{"cfgChi2FunctionMatchingScoreCut5", 0.5f, "Minimum score value for selecting good matches"};
    std::array<Configurable<float>*, Chi2FunctionsNum> matchingScoreCuts{
      &cfgChi2FunctionMatchingScoreCut1, &cfgChi2FunctionMatchingScoreCut2, &cfgChi2FunctionMatchingScoreCut3, &cfgChi2FunctionMatchingScoreCut4, &cfgChi2FunctionMatchingScoreCut5};

    Configurable<float> cfgChi2FunctionMatchingPlaneZ1{"cfgChi2FunctionMatchingPlaneZ1", static_cast<float>(o2::mft::constants::mft::LayerZCoordinate()[9]), "Z position of the matching plane"};
    Configurable<float> cfgChi2FunctionMatchingPlaneZ2{"cfgChi2FunctionMatchingPlaneZ2", static_cast<float>(o2::mft::constants::mft::LayerZCoordinate()[9]), "Z position of the matching plane"};
    Configurable<float> cfgChi2FunctionMatchingPlaneZ3{"cfgChi2FunctionMatchingPlaneZ3", static_cast<float>(o2::mft::constants::mft::LayerZCoordinate()[9]), "Z position of the matching plane"};
    Configurable<float> cfgChi2FunctionMatchingPlaneZ4{"cfgChi2FunctionMatchingPlaneZ4", static_cast<float>(o2::mft::constants::mft::LayerZCoordinate()[9]), "Z position of the matching plane"};
    Configurable<float> cfgChi2FunctionMatchingPlaneZ5{"cfgChi2FunctionMatchingPlaneZ5", static_cast<float>(o2::mft::constants::mft::LayerZCoordinate()[9]), "Z position of the matching plane"};
    std::array<Configurable<float>*, Chi2FunctionsNum> matchingPlaneZs{
      &cfgChi2FunctionMatchingPlaneZ1, &cfgChi2FunctionMatchingPlaneZ2, &cfgChi2FunctionMatchingPlaneZ3, &cfgChi2FunctionMatchingPlaneZ4, &cfgChi2FunctionMatchingPlaneZ5};

    Configurable<int> cfgChi2MatchingExtrapMethod1{"cfgChi2MatchingExtrapMethod1", 0, "Method for MCH track extrapolation to maching plane"};
    Configurable<int> cfgChi2MatchingExtrapMethod2{"cfgChi2MatchingExtrapMethod2", 0, "Method for MCH track extrapolation to maching plane"};
    Configurable<int> cfgChi2MatchingExtrapMethod3{"cfgChi2MatchingExtrapMethod3", 0, "Method for MCH track extrapolation to maching plane"};
    Configurable<int> cfgChi2MatchingExtrapMethod4{"cfgChi2MatchingExtrapMethod4", 0, "Method for MCH track extrapolation to maching plane"};
    Configurable<int> cfgChi2MatchingExtrapMethod5{"cfgChi2MatchingExtrapMethod5", 0, "Method for MCH track extrapolation to maching plane"};
    std::array<Configurable<int>*, Chi2FunctionsNum> matchingExtrapMethods{
      &cfgChi2MatchingExtrapMethod1, &cfgChi2MatchingExtrapMethod2, &cfgChi2MatchingExtrapMethod3, &cfgChi2MatchingExtrapMethod4, &cfgChi2MatchingExtrapMethod5};
  } configChi2MatchingOptions;

  // ML interface
  static constexpr int MlModelsNum = 5;
  struct : ConfigurableGroup {
    Configurable<std::string> cfgMlModelLabel1{"cfgMlModelLabel1", std::string{""}, "Text label identifying this group of ML models"};
    Configurable<std::string> cfgMlModelLabel2{"cfgMlModelLabel2", std::string{""}, "Text label identifying this group of ML models"};
    Configurable<std::string> cfgMlModelLabel3{"cfgMlModelLabel3", std::string{""}, "Text label identifying this group of ML models"};
    Configurable<std::string> cfgMlModelLabel4{"cfgMlModelLabel4", std::string{""}, "Text label identifying this group of ML models"};
    Configurable<std::string> cfgMlModelLabel5{"cfgMlModelLabel5", std::string{""}, "Text label identifying this group of ML models"};
    std::array<Configurable<std::string>*, MlModelsNum> modelLabels{
      &cfgMlModelLabel1, &cfgMlModelLabel2, &cfgMlModelLabel3, &cfgMlModelLabel4, &cfgMlModelLabel5};

    Configurable<std::string> cfgMlModelPathCcdb1{"cfgMlModelPathCcdb1", "Users/m/mcoquet/MLTest", "Paths of models on CCDB"};
    Configurable<std::string> cfgMlModelPathCcdb2{"cfgMlModelPathCcdb2", std::string{""}, "Paths of models on CCDB"};
    Configurable<std::string> cfgMlModelPathCcdb3{"cfgMlModelPathCcdb3", std::string{""}, "Paths of models on CCDB"};
    Configurable<std::string> cfgMlModelPathCcdb4{"cfgMlModelPathCcdb4", std::string{""}, "Paths of models on CCDB"};
    Configurable<std::string> cfgMlModelPathCcdb5{"cfgMlModelPathCcdb5", std::string{""}, "Paths of models on CCDB"};
    std::array<Configurable<std::string>*, MlModelsNum> modelPathCcds{
      &cfgMlModelPathCcdb1, &cfgMlModelPathCcdb2, &cfgMlModelPathCcdb3, &cfgMlModelPathCcdb4, &cfgMlModelPathCcdb5};

    Configurable<std::string> cfgMlModelName1{"cfgMlModelName1", "model.onnx", "ONNX file names for each pT bin (if not from CCDB full path)"};
    Configurable<std::string> cfgMlModelName2{"cfgMlModelName2", std::string{""}, "ONNX file names for each pT bin (if not from CCDB full path)"};
    Configurable<std::string> cfgMlModelName3{"cfgMlModelName3", std::string{""}, "ONNX file names for each pT bin (if not from CCDB full path)"};
    Configurable<std::string> cfgMlModelName4{"cfgMlModelName4", std::string{""}, "ONNX file names for each pT bin (if not from CCDB full path)"};
    Configurable<std::string> cfgMlModelName5{"cfgMlModelName5", std::string{""}, "ONNX file names for each pT bin (if not from CCDB full path)"};
    std::array<Configurable<std::string>*, MlModelsNum> modelNames{
      &cfgMlModelName1, &cfgMlModelName2, &cfgMlModelName3, &cfgMlModelName4, &cfgMlModelName5};

    Configurable<std::string> cfgMlInputFeatures1{"cfgMlInputFeatures1", "chi2MCHMFT", "Names of ML model input features"};
    Configurable<std::string> cfgMlInputFeatures2{"cfgMlInputFeatures2", std::string{""}, "Names of ML model input features"};
    Configurable<std::string> cfgMlInputFeatures3{"cfgMlInputFeatures3", std::string{""}, "Names of ML model input features"};
    Configurable<std::string> cfgMlInputFeatures4{"cfgMlInputFeatures4", std::string{""}, "Names of ML model input features"};
    Configurable<std::string> cfgMlInputFeatures5{"cfgMlInputFeatures5", std::string{""}, "Names of ML model input features"};
    std::array<Configurable<std::string>*, MlModelsNum> inputFeatures{
      &cfgMlInputFeatures1, &cfgMlInputFeatures2, &cfgMlInputFeatures3, &cfgMlInputFeatures4, &cfgMlInputFeatures5};

    Configurable<float> cfgMlModelMatchingScoreCut1{"cfgMlModelMatchingScoreCut1", 0.f, "Minimum score value for selecting good matches"};
    Configurable<float> cfgMlModelMatchingScoreCut2{"cfgMlModelMatchingScoreCut2", 0.f, "Minimum score value for selecting good matches"};
    Configurable<float> cfgMlModelMatchingScoreCut3{"cfgMlModelMatchingScoreCut3", 0.f, "Minimum score value for selecting good matches"};
    Configurable<float> cfgMlModelMatchingScoreCut4{"cfgMlModelMatchingScoreCut4", 0.f, "Minimum score value for selecting good matches"};
    Configurable<float> cfgMlModelMatchingScoreCut5{"cfgMlModelMatchingScoreCut5", 0.f, "Minimum score value for selecting good matches"};
    std::array<Configurable<float>*, MlModelsNum> matchingScoreCuts{
      &cfgMlModelMatchingScoreCut1, &cfgMlModelMatchingScoreCut2, &cfgMlModelMatchingScoreCut3, &cfgMlModelMatchingScoreCut4, &cfgMlModelMatchingScoreCut5};

    Configurable<float> cfgMlModelMatchingPlaneZ1{"cfgMlModelMatchingPlaneZ1", static_cast<float>(o2::mft::constants::mft::LayerZCoordinate()[9]), "Z position of the matching plane"};
    Configurable<float> cfgMlModelMatchingPlaneZ2{"cfgMlModelMatchingPlaneZ2", 0.f, "Z position of the matching plane"};
    Configurable<float> cfgMlModelMatchingPlaneZ3{"cfgMlModelMatchingPlaneZ3", 0.f, "Z position of the matching plane"};
    Configurable<float> cfgMlModelMatchingPlaneZ4{"cfgMlModelMatchingPlaneZ4", 0.f, "Z position of the matching plane"};
    Configurable<float> cfgMlModelMatchingPlaneZ5{"cfgMlModelMatchingPlaneZ5", 0.f, "Z position of the matching plane"};
    std::array<Configurable<float>*, MlModelsNum> matchingPlaneZs{
      &cfgMlModelMatchingPlaneZ1, &cfgMlModelMatchingPlaneZ2, &cfgMlModelMatchingPlaneZ3, &cfgMlModelMatchingPlaneZ4, &cfgMlModelMatchingPlaneZ5};

    Configurable<int> cfgMlMatchingExtrapMethod1{"cfgMlMatchingExtrapMethod1", 0, "Method for MCH track extrapolation to maching plane"};
    Configurable<int> cfgMlMatchingExtrapMethod2{"cfgMlMatchingExtrapMethod2", 0, "Method for MCH track extrapolation to maching plane"};
    Configurable<int> cfgMlMatchingExtrapMethod3{"cfgMlMatchingExtrapMethod3", 0, "Method for MCH track extrapolation to maching plane"};
    Configurable<int> cfgMlMatchingExtrapMethod4{"cfgMlMatchingExtrapMethod4", 0, "Method for MCH track extrapolation to maching plane"};
    Configurable<int> cfgMlMatchingExtrapMethod5{"cfgMlMatchingExtrapMethod5", 0, "Method for MCH track extrapolation to maching plane"};
    std::array<Configurable<int>*, MlModelsNum> matchingExtrapMethods{
      &cfgMlMatchingExtrapMethod1, &cfgMlMatchingExtrapMethod2, &cfgMlMatchingExtrapMethod3, &cfgMlMatchingExtrapMethod4, &cfgMlMatchingExtrapMethod5};
  } configMlOptions;

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
  Matrix<o2::framework::HistPtr, 4, 4> dimuonHistos;

  Produces<o2::aod::QaMatchingEvents> qaMatchingEvents;
  Produces<o2::aod::QaMatchingMCHTrack> qaMatchingMCHTrack;
  Produces<o2::aod::QaMatchingCandidates> qaMatchingCandidates;

  struct EfficiencyPlotter {
    o2::framework::HistPtr pNum;
    o2::framework::HistPtr pDen;
    o2::framework::HistPtr pPdgNum;
    o2::framework::HistPtr pPdgDen;
    o2::framework::HistPtr ptNum;
    o2::framework::HistPtr ptDen;
    o2::framework::HistPtr ptPdgNum;
    o2::framework::HistPtr ptPdgDen;
    o2::framework::HistPtr phiNum;
    o2::framework::HistPtr phiDen;
    o2::framework::HistPtr phiPdgNum;
    o2::framework::HistPtr phiPdgDen;
    o2::framework::HistPtr etaNum;
    o2::framework::HistPtr etaDen;
    o2::framework::HistPtr etaPdgNum;
    o2::framework::HistPtr etaPdgDen;

    EfficiencyPlotter(std::string path, std::string title,
                      HistogramRegistry& registry, bool createPdgMomHistograms)
    {
      AxisSpec pAxis = {100, 0, 100, "p (GeV/c)"};
      AxisSpec pTAxis = {100, 0, 10, "p_{T} (GeV/c)"};
      AxisSpec etaAxis = {100, -4, -2, "#eta"};
      AxisSpec phiAxis = {90, -180, 180, "#phi (degrees)"};
      AxisSpec motherPDGAxis{1201, -600.5, 600.5, "Direct mother PDG"};

      std::string histName;
      std::string histTitle;

      // momentum dependence
      histName = path + "p_num";
      histTitle = title + " vs. p - num";
      pNum = registry.add(histName.c_str(), histTitle.c_str(), {HistType::kTH1F, {pAxis}});

      histName = path + "p_den";
      histTitle = title + " vs. p - den";
      pDen = registry.add(histName.c_str(), histTitle.c_str(), {HistType::kTH1F, {pAxis}});

      if (createPdgMomHistograms) {
        histName = path + "p_pdg_num";
        histTitle = title + " vs. p vs pdg ID - num";
        pPdgNum = registry.add(histName.c_str(), histTitle.c_str(), {HistType::kTH2F, {pAxis, motherPDGAxis}});

        histName = path + "p_pdg_den";
        histTitle = title + " vs. p vs pdg ID - den";
        pPdgDen = registry.add(histName.c_str(), histTitle.c_str(), {HistType::kTH2F, {pAxis, motherPDGAxis}});
      }

      // pT dependence
      histName = path + "pt_num";
      histTitle = title + " vs. p_{T} - num";
      ptNum = registry.add(histName.c_str(), histTitle.c_str(), {HistType::kTH1F, {pTAxis}});

      histName = path + "pt_den";
      histTitle = title + " vs. p_{T} - den";
      ptDen = registry.add(histName.c_str(), histTitle.c_str(), {HistType::kTH1F, {pTAxis}});

      if (createPdgMomHistograms) {
        histName = path + "pt_pdg_num";
        histTitle = title + " vs. p_{T} vs pdg ID - num";
        ptPdgNum = registry.add(histName.c_str(), histTitle.c_str(), {HistType::kTH2F, {pTAxis, motherPDGAxis}});

        histName = path + "pt_pdg_den";
        histTitle = title + " vs. p_{T} vs pdg ID - den";
        ptPdgDen = registry.add(histName.c_str(), histTitle.c_str(), {HistType::kTH2F, {pTAxis, motherPDGAxis}});
      }

      // eta dependence
      histName = path + "eta_num";
      histTitle = title + " vs. #eta - num";
      etaNum = registry.add(histName.c_str(), histTitle.c_str(), {HistType::kTH1F, {etaAxis}});

      histName = path + "eta_den";
      histTitle = title + " vs. #eta - den";
      etaDen = registry.add(histName.c_str(), histTitle.c_str(), {HistType::kTH1F, {etaAxis}});

      if (createPdgMomHistograms) {
        histName = path + "eta_pdg_num";
        histTitle = title + " vs. #eta vs pdg ID - num";
        etaPdgNum = registry.add(histName.c_str(), histTitle.c_str(), {HistType::kTH2F, {etaAxis, motherPDGAxis}});

        histName = path + "eta_pdg_den";
        histTitle = title + " vs. #eta vs pdg ID - den";
        etaPdgDen = registry.add(histName.c_str(), histTitle.c_str(), {HistType::kTH2F, {etaAxis, motherPDGAxis}});
      }

      // phi dependence
      histName = path + "phi_num";
      histTitle = title + " vs. #phi - num";
      phiNum = registry.add(histName.c_str(), histTitle.c_str(), {HistType::kTH1F, {phiAxis}});

      histName = path + "phi_den";
      histTitle = title + " vs. #phi - den";
      phiDen = registry.add(histName.c_str(), histTitle.c_str(), {HistType::kTH1F, {phiAxis}});

      if (createPdgMomHistograms) {
        histName = path + "phi_pdg_num";
        histTitle = title + " vs. #phi vs pdg ID - num";
        phiPdgNum = registry.add(histName.c_str(), histTitle.c_str(), {HistType::kTH2F, {phiAxis, motherPDGAxis}});

        histName = path + "phi_pdg_den";
        histTitle = title + " vs. #phi vs pdg ID - den";
        phiPdgDen = registry.add(histName.c_str(), histTitle.c_str(), {HistType::kTH2F, {phiAxis, motherPDGAxis}});
      }
    }

    template <class T>
    void fill(const T& track, bool passed)
    {
      constexpr double RadToDeg = 180. / o2::constants::math::PI;
      double phi = track.phi() * RadToDeg;
      std::get<std::shared_ptr<TH1>>(pDen)->Fill(track.p());
      std::get<std::shared_ptr<TH1>>(ptDen)->Fill(track.pt());
      std::get<std::shared_ptr<TH1>>(etaDen)->Fill(track.eta());
      std::get<std::shared_ptr<TH1>>(phiDen)->Fill(phi);

      if (passed) {
        std::get<std::shared_ptr<TH1>>(pNum)->Fill(track.p());
        std::get<std::shared_ptr<TH1>>(ptNum)->Fill(track.pt());
        std::get<std::shared_ptr<TH1>>(etaNum)->Fill(track.eta());
        std::get<std::shared_ptr<TH1>>(phiNum)->Fill(phi);
      }
    }

    // Study the PDG origin of particles and their effect on the purity score
    template <class T>
    void fill(const T& track, int pdgCode, bool passed)
    {
      constexpr double RadToDeg = 180. / o2::constants::math::PI;
      double phi = track.phi() * RadToDeg;
      std::get<std::shared_ptr<TH2>>(pPdgDen)->Fill(track.p(), pdgCode);
      std::get<std::shared_ptr<TH2>>(ptPdgDen)->Fill(track.pt(), pdgCode);
      std::get<std::shared_ptr<TH2>>(etaPdgDen)->Fill(track.eta(), pdgCode);
      std::get<std::shared_ptr<TH2>>(phiPdgDen)->Fill(phi, pdgCode);

      if (passed) {
        std::get<std::shared_ptr<TH2>>(pPdgNum)->Fill(track.p(), pdgCode);
        std::get<std::shared_ptr<TH2>>(ptPdgNum)->Fill(track.pt(), pdgCode);
        std::get<std::shared_ptr<TH2>>(etaPdgNum)->Fill(track.eta(), pdgCode);
        std::get<std::shared_ptr<TH2>>(phiPdgNum)->Fill(phi, pdgCode);
      }
    }
  };

  struct MatchRankingHistos {
    o2::framework::HistPtr hist;
    o2::framework::HistPtr histVsP;
    o2::framework::HistPtr histVsPt;
    o2::framework::HistPtr histVsMcParticleDz;
    o2::framework::HistPtr histVsMftTrackMult;
    o2::framework::HistPtr histVsMatchAttempts;
    o2::framework::HistPtr histVsMftTrackType;
    o2::framework::HistPtr histVsDeltaChi2;
    o2::framework::HistPtr histVsProdRanking;

    MatchRankingHistos(std::string histName, std::string histTitle, HistogramRegistry* registry, int mftMultMax, int numCandidates)
    {
      AxisSpec pAxis = {100, 0, 100, "p (GeV/c)"};
      AxisSpec ptAxis = {100, 0, 10, "p_{T} (GeV/c)"};
      AxisSpec dzAxis = {100, -1, 4, "#Deltaz (cm)"};
      AxisSpec trackMultAxis = {static_cast<int>(mftMultMax) / 10, 0, static_cast<double>(mftMultMax), "MFT track mult."};
      AxisSpec matchAttemptsAxis = {static_cast<int>(mftMultMax) / 10, 0, static_cast<double>(mftMultMax), "match attempts"};
      AxisSpec trackTypeAxis = {2, 0, 2, "MFT track type"};
      int matchTypeMax = static_cast<int>(kMatchTypeUndefined);
      AxisSpec matchTypeAxis = {matchTypeMax, 0, static_cast<double>(matchTypeMax), "match type"};
      AxisSpec dchi2Axis = {100, 0, 100, "#Delta#chi^{2}"};
      AxisSpec dqAxis = {3, -1.5, 1.5, "MFT #DeltaQ"};
      AxisSpec indexAxis = {numCandidates + 1, 0, static_cast<double>(numCandidates + 1), "ranking index"};
      AxisSpec indexProdAxis = {numCandidates + 1, 0, static_cast<double>(numCandidates + 1), "ranking index (production)"};

      hist = registry->add(histName.c_str(), histTitle.c_str(), {HistType::kTH1F, {indexAxis}});
      histVsP = registry->add((histName + "VsP").c_str(), (histTitle + " vs. p").c_str(), {HistType::kTH2F, {pAxis, indexAxis}});
      histVsPt = registry->add((histName + "VsPt").c_str(), (histTitle + " vs. p_{T}").c_str(), {HistType::kTH2F, {ptAxis, indexAxis}});
      histVsMcParticleDz = registry->add((histName + "VsMcParticleDz").c_str(), (histTitle + " vs. MC particle #Deltaz").c_str(), {HistType::kTH2F, {dzAxis, indexAxis}});
      histVsMftTrackMult = registry->add((histName + "VsMftTrackMult").c_str(), (histTitle + " vs. MFT track multiplicity").c_str(), {HistType::kTH2F, {trackMultAxis, indexAxis}});
      histVsMatchAttempts = registry->add((histName + "VsMatchAttempts").c_str(), (histTitle + " vs. MFT track multiplicity").c_str(), {HistType::kTH2F, {matchAttemptsAxis, indexAxis}});
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
                    HistogramRegistry* reg,
                    bool createPdgMomHistograms,
                    int mftMultMax,
                    int numCandidates)
      : fMatchingPurityPlotter(path + "matching-purity/", "Matching purity", *reg, createPdgMomHistograms),
        fPairingEfficiencyPlotter(path + "pairing-efficiency/", "Pairing efficiency", *reg, createPdgMomHistograms),
        fMatchingEfficiencyPlotter(path + "matching-efficiency/", "Matching efficiency", *reg, createPdgMomHistograms),
        fFakeMatchingEfficiencyPlotter(path + "fake-matching-efficiency/", "Fake matching efficiency", *reg, createPdgMomHistograms)
    {
      registry = reg;
      AxisSpec pAxis = {100, 0, 100, "p (GeV/c)"};
      AxisSpec ptAxis = {100, 0, 10, "p_{T} (GeV/c)"};
      AxisSpec dzAxis = {100, 0, 50, "#Deltaz (cm)"};
      AxisSpec indexAxis = {6, 0, 6, "ranking index"};

      std::string histName = path + "matchRanking";
      std::string histTitle = "True match ranking";

      fMatchRanking = std::make_unique<MatchRankingHistos>(path + "matchRanking", "True match ranking", registry, mftMultMax, numCandidates);
      fMatchRankingGoodMCH = std::make_unique<MatchRankingHistos>(path + "matchRankingGoodMCH", "True match ranking (good MCH tracks)", registry, mftMultMax, numCandidates);
      fMatchRankingPaired = std::make_unique<MatchRankingHistos>(path + "matchRankingPaired", "True match ranking (paired MCH tracks)", registry, mftMultMax, numCandidates);
      fMatchRankingPairedGoodMCH = std::make_unique<MatchRankingHistos>(path + "matchRankingPairedGoodMCH", "True match ranking (good paired MCH tracks)", registry, mftMultMax, numCandidates);

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
  void initCcdb(BC const& bc)
  {
    if (mRunNumber == bc.runNumber())
      return;

    mRunNumber = bc.runNumber();
    std::map<std::string, std::string> metadata;
    auto soreor = o2::ccdb::BasicCCDBManager::getRunDuration(fCCDBApi, mRunNumber);
    auto ts = soreor.first;
    auto grpmag = fCCDBApi.retrieveFromTFileAny<o2::parameters::GRPMagField>(grpMagPath, metadata, ts);
    o2::base::Propagator::initFieldFromGRP(grpmag);
    LOGF(info, "Set field for muons");
    VarManager::SetupMuonMagField();
    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      ccdbManager->get<TGeoManager>(geoPath);
    }
    o2::mch::TrackExtrap::setField();
    auto* fieldB = static_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField());
    if (fieldB) {
      double centerMft[3] = {0, 0, -61.4}; // Field at center of MFT
      mBzAtMftCenter = fieldB->getBz(centerMft);
      // std::cout << "fieldB: " << (void*)fieldB << std::endl;
    }
  }

  void createMatchingHistosMc()
  {
    AxisSpec chi2Axis = {1000, 0, 1000, "chi^{2}"};
    AxisSpec chi2AxisSmall = {200, 0, 100, "chi^{2}"};
    AxisSpec pAxis = {1000, 0, 100, "p (GeV/c)"};
    AxisSpec pTAxis = {100, 0, 10, "p_{T} (GeV/c)"};
    AxisSpec etaAxis = {100, -4, -2, "#eta"};
    AxisSpec phiAxis = {90, -180, 180, "#phi (degrees)"};
    std::string histPath = "matching/MC/";

    AxisSpec trackPositionXAtMftAxis = {100, -15, 15, "MFT x (cm)"};
    AxisSpec trackPositionYAtMftAxis = {100, -15, 15, "MFT y (cm)"};
    registry.add((histPath + "pairedMCHTracksAtMFT").c_str(), "Paired MCH tracks position at MFT end", {HistType::kTH2F, {trackPositionXAtMftAxis, trackPositionYAtMftAxis}});
    registry.add((histPath + "pairedMFTTracksAtMFT").c_str(), "Paired MFT tracks position at MFT end", {HistType::kTH2F, {trackPositionXAtMftAxis, trackPositionYAtMftAxis}});
    registry.add((histPath + "selectedMCHTracksAtMFT").c_str(), "Selected MCH tracks position at MFT end", {HistType::kTH2F, {trackPositionXAtMftAxis, trackPositionYAtMftAxis}});
    registry.add((histPath + "selectedMCHTracksAtMFTTrue").c_str(), "Selected MCH tracks position at MFT end - true", {HistType::kTH2F, {trackPositionXAtMftAxis, trackPositionYAtMftAxis}});
    registry.add((histPath + "selectedMCHTracksAtMFTFake").c_str(), "Selected MCH tracks position at MFT end - fake", {HistType::kTH2F, {trackPositionXAtMftAxis, trackPositionYAtMftAxis}});

    fChi2MatchingPlotter = std::make_unique<MatchingPlotter>(histPath + "Prod/", &registryMatching, configQas.cfgCreatePdgMomHistograms, cfgMftTrackMultiplicityMax, cfgNCandidatesMax);
    int registryIndex = 0;
    for (const auto& [label, func] : matchingChi2Functions) {
      fMatchingPlotters[label] = std::make_unique<MatchingPlotter>(histPath + label + "/", registryMatchingVec[registryIndex], configQas.cfgCreatePdgMomHistograms, cfgMftTrackMultiplicityMax, cfgNCandidatesMax);
      registryIndex += 1;
    }
    for (const auto& [label, response] : matchingMlResponses) {
      fMatchingPlotters[label] = std::make_unique<MatchingPlotter>(histPath + label + "/", (registryMatchingVec[registryIndex]), configQas.cfgCreatePdgMomHistograms, cfgMftTrackMultiplicityMax, cfgNCandidatesMax);
      registryIndex += 1;
    }

    fTaggedMuonsMatchingPlotter = std::make_unique<MatchingPlotter>(histPath + "Tagged/", &registryMatching, configQas.cfgCreatePdgMomHistograms, cfgMftTrackMultiplicityMax, cfgNCandidatesMax);
    fSelectedMuonsMatchingPlotter = std::make_unique<MatchingPlotter>(histPath + "Selected/", &registryMatching, configQas.cfgCreatePdgMomHistograms, cfgMftTrackMultiplicityMax, cfgNCandidatesMax);
  }

  void createDimuonHistos()
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

  void initMatchingFunctions()
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

      SMatrix55Sym hK, vK;
      SVector5 mK(mftTrack.getX(), mftTrack.getY(), mftTrack.getPhi(),
                  mftTrack.getTanl(), mftTrack.getInvQPt()),
        rKKminus1;
      SVector5 globalMuonTrackParameters = mchTrack.getParameters();
      SMatrix55Sym globalMuonTrackCovariances = mchTrack.getCovariances();
      vK(0, 0) = mftTrack.getCovariances()(0, 0);
      vK(1, 1) = mftTrack.getCovariances()(1, 1);
      vK(2, 2) = mftTrack.getCovariances()(2, 2);
      vK(3, 3) = mftTrack.getCovariances()(3, 3);
      vK(4, 4) = mftTrack.getCovariances()(4, 4);
      hK(0, 0) = 1.0;
      hK(1, 1) = 1.0;
      hK(2, 2) = 1.0;
      hK(3, 3) = 1.0;
      hK(4, 4) = 1.0;

      // Covariance of residuals
      SMatrix55Std invResCov = (vK + ROOT::Math::Similarity(hK, globalMuonTrackCovariances));
      invResCov.Invert();

      // Update Parameters
      rKKminus1 = mK - hK * globalMuonTrackParameters; // Residuals of prediction

      auto matchChi2Track = ROOT::Math::Similarity(rKKminus1, invResCov);

      // return chi2 and NDF
      return {matchChi2Track, 5};
    };

    //________________________________________________________________________________
    mMatchingFunctionMap["matchXYPhiTanl"] = [](const o2::dataformats::GlobalFwdTrack& mchTrack, const o2::track::TrackParCovFwd& mftTrack) -> std::tuple<double, int> {
      // Match two tracks evaluating positions & angles

      SMatrix45 hK;
      SMatrix44 vK;
      SVector4 mK(mftTrack.getX(), mftTrack.getY(), mftTrack.getPhi(),
                  mftTrack.getTanl()),
        rKKminus1;
      SVector5 globalMuonTrackParameters = mchTrack.getParameters();
      SMatrix55Sym globalMuonTrackCovariances = mchTrack.getCovariances();
      vK(0, 0) = mftTrack.getCovariances()(0, 0);
      vK(1, 1) = mftTrack.getCovariances()(1, 1);
      vK(2, 2) = mftTrack.getCovariances()(2, 2);
      vK(3, 3) = mftTrack.getCovariances()(3, 3);
      hK(0, 0) = 1.0;
      hK(1, 1) = 1.0;
      hK(2, 2) = 1.0;
      hK(3, 3) = 1.0;

      // Covariance of residuals
      SMatrix44 invResCov = (vK + ROOT::Math::Similarity(hK, globalMuonTrackCovariances));
      invResCov.Invert();

      // Residuals of prediction
      rKKminus1 = mK - hK * globalMuonTrackParameters;

      auto matchChi2Track = ROOT::Math::Similarity(rKKminus1, invResCov);

      // return chi2 and NDF
      return {matchChi2Track, 4};
    };

    //________________________________________________________________________________
    mMatchingFunctionMap["matchXY"] = [](const o2::dataformats::GlobalFwdTrack& mchTrack, const o2::track::TrackParCovFwd& mftTrack) -> std::tuple<double, int> {
      // Calculate Matching Chi2 - X and Y positions

      SMatrix25 hK;
      SMatrix22 vK;
      SVector2 mK(mftTrack.getX(), mftTrack.getY()), rKKminus1;
      SVector5 globalMuonTrackParameters = mchTrack.getParameters();
      SMatrix55Sym globalMuonTrackCovariances = mchTrack.getCovariances();
      vK(0, 0) = mftTrack.getCovariances()(0, 0);
      vK(1, 1) = mftTrack.getCovariances()(1, 1);
      hK(0, 0) = 1.0;
      hK(1, 1) = 1.0;

      // Covariance of residuals
      SMatrix22 invResCov = (vK + ROOT::Math::Similarity(hK, globalMuonTrackCovariances));
      invResCov.Invert();

      // Residuals of prediction
      rKKminus1 = mK - hK * globalMuonTrackParameters;
      auto matchChi2Track = ROOT::Math::Similarity(rKKminus1, invResCov);

      // return reduced chi2
      return {matchChi2Track, 2};
    };
  }

  void init(o2::framework::InitContext&)
  {
    // Load geometry
    ccdbManager->setURL(ccdbUrl);
    ccdbManager->setCaching(true);
    ccdbManager->setLocalObjectValidityChecking();
    fCCDBApi.init(ccdbUrl);
    mRunNumber = 0;

    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      LOGF(info, "Load geometry from CCDB");
      ccdbManager->get<TGeoManager>(geoPath);
    }

    // Matching functions
    initMatchingFunctions();
    for (size_t funcId = 0; funcId < Chi2FunctionsNum; funcId++) {
      auto label = configChi2MatchingOptions.functionLabels[funcId]->value;
      auto funcName = configChi2MatchingOptions.functionNames[funcId]->value;
      auto scoreMin = configChi2MatchingOptions.matchingScoreCuts[funcId]->value;
      auto matchingPlaneZ = configChi2MatchingOptions.matchingPlaneZs[funcId]->value;
      auto extrapMethod = configChi2MatchingOptions.matchingExtrapMethods[funcId]->value;

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

    for (size_t modelId = 0; modelId < MlModelsNum; modelId++) {
      auto label = configMlOptions.modelLabels[modelId]->value;
      auto modelPath = configMlOptions.modelPathCcds[modelId]->value;
      auto inputFeatures = configMlOptions.inputFeatures[modelId]->value;
      auto modelName = configMlOptions.modelNames[modelId]->value;
      auto scoreMin = configMlOptions.matchingScoreCuts[modelId]->value;
      auto matchingPlaneZ = configMlOptions.matchingPlaneZs[modelId]->value;
      auto extrapMethod = configMlOptions.matchingExtrapMethods[modelId]->value;

      if (label == "" || modelPath == "" || inputFeatures == "" || modelName == "")
        break;

      matchingMlResponses[label].configure(binsPtMl, mycutsMl, cutDirMl, 1);
      matchingMlResponses[label].setModelPathsCCDB(std::vector<std::string>{modelName}, fCCDBApi, std::vector<std::string>{modelPath}, configCcdb.cfgCcdbNoLaterThan.value);
      matchingMlResponses[label].cacheInputFeaturesIndices(std::vector<std::string>{inputFeatures});
      matchingMlResponses[label].init();

      matchingScoreCuts[label] = scoreMin;
      matchingPlanesZ[label] = matchingPlaneZ;
      matchingExtrapMethod[label] = extrapMethod;
    }

    int nTrackTypes = static_cast<int>(o2::aod::fwdtrack::ForwardTrackTypeEnum::MCHStandaloneTrack) + 1;
    AxisSpec trackTypeAxis = {static_cast<int>(nTrackTypes), 0.0, static_cast<double>(nTrackTypes), "track type"};
    registry.add("nTracksPerType", "Number of tracks per type", {HistType::kTH1F, {trackTypeAxis}});

    AxisSpec tracksMultiplicityAxis = {cfgMftTrackMultiplicityMax, 0, static_cast<double>(cfgMftTrackMultiplicityMax), "tracks multiplicity"};
    registry.add("tracksMultiplicityMFT", "MFT tracks multiplicity", {HistType::kTH1F, {tracksMultiplicityAxis}});
    registry.add("tracksMultiplicityMCH", "MCH tracks multiplicity", {HistType::kTH1F, {tracksMultiplicityAxis}});

    createMatchingHistosMc();
    createDimuonHistos();
  }

  template <class T, class C>
  bool pDcaCut(const T& mchTrack, const C& collision, double nSigmaPDCA)
  {
    static const double sigmaPDCA23 = 80.;
    static const double sigmaPDCA310 = 54.;
    static const double relPRes = 0.0004;
    static const double slopeRes = 0.0005;

    constexpr double AbsorberEndZ = 505.;
    constexpr double RadToDeg = 180. / o2::constants::math::PI;
    double thetaAbs = std::atan(mchTrack.rAtAbsorberEnd() / AbsorberEndZ) * RadToDeg;

    // propagate muon track to vertex
    auto mchTrackAtVertex = VarManager::PropagateMuon(mchTrack, collision, VarManager::kToVertex);

    // double pUncorr = mchTrack.p();
    double p = mchTrackAtVertex.getP();

    double pDCA = mchTrack.pDca();
    double sigmaPDCA = (thetaAbs < ThetaAbsBoundaryDeg) ? sigmaPDCA23 : sigmaPDCA310;
    double nrp = nSigmaPDCA * relPRes * p;
    double pResEffect = sigmaPDCA / (1. - nrp / (1. + nrp));
    double slopeResEffect = SlopeResolutionZ * slopeRes * p;
    double sigmaPDCAWithRes = std::sqrt(pResEffect * pResEffect + slopeResEffect * slopeResEffect);
    if (pDCA > nSigmaPDCA * sigmaPDCAWithRes) {
      return false;
    }

    return true;
  }

  template <class T, class C>
  bool isGoodMuon(const T& mchTrack, const C& collision,
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
    if (!pDcaCut(mchTrack, collision, nSigmaPdcaCut)) {
      return false;
    }

    return true;
  }

  template <class T, class C>
  bool isGoodMuon(const T& muonTrack, const C& collision)
  {
    return isGoodMuon(muonTrack, collision, cfgTrackChi2MchUp, cfgPMchLow, cfgPtMchLow, {cfgEtaMchLow, cfgEtaMchUp}, {cfgRabsLow, cfgRabsUp}, cfgPdcaUp);
  }

  template <class T, class C>
  bool isGoodGlobalMuon(const T& muonTrack, const C& collision)
  {
    return isGoodMuon(muonTrack, collision, cfgTrackChi2MchUp, cfgPMchLow, cfgPtMchLow, {cfgEtaMftLow, cfgEtaMftUp}, {cfgRabsLow, cfgRabsUp}, cfgPdcaUp);
  }

  template <class T>
  bool isGoodMft(const T& mftTrack,
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
  bool isGoodMft(const T& mftTrack)
  {
    return isGoodMft(mftTrack, cfgTrackChi2MftUp, cfgTrackNClustMftLow);
  }

  template <class TMUON>
  bool isGoodGlobalMatching(const TMUON& muonTrack,
                            double matchingScore,
                            double matchingScoreCut)
  {
    if (static_cast<int>(muonTrack.trackType()) > GlobalTrackTypeMax)
      return false;

    // MFT-MCH matching score cut
    if (matchingScore < matchingScoreCut)
      return false;

    return true;
  }

  template <class TMUON>
  bool isGoodGlobalMatching(const TMUON& muonTrack, double matchingScore)
  {
    return isGoodGlobalMatching(muonTrack, matchingScore, cfgMatchingChi2ScoreMftMchLow);
  }

  template <class TMUON>
  bool isTrueGlobalMatching(const TMUON& muonTrack, const std::vector<std::pair<int64_t, int64_t>>& matchablePairs)
  {
    if (static_cast<int>(muonTrack.trackType()) > GlobalTrackTypeMax)
      return false;

    int64_t mchTrackId = static_cast<int64_t>(muonTrack.matchMCHTrackId());
    int64_t mftTrackId = static_cast<int64_t>(muonTrack.matchMFTTrackId());

    std::pair<int64_t, int64_t> trackIndexes = std::make_pair(mchTrackId, mftTrackId);

    return (std::find(matchablePairs.begin(), matchablePairs.end(), trackIndexes) != matchablePairs.end());
  }

  bool isMatchableMch(int64_t mchTrackId, const std::vector<std::pair<int64_t, int64_t>>& matchablePairs)
  {
    for (const auto& [id1, id2] : matchablePairs) {
      if (mchTrackId == id1)
        return true;
    }
    return false;
  }

  std::optional<std::pair<int64_t, int64_t>> getMatchablePairForMch(int64_t mchTrackId, const std::vector<std::pair<int64_t, int64_t>>& matchablePairs)
  {
    for (const auto& pair : matchablePairs) {
      if (mchTrackId == pair.first)
        return pair;
    }
    return {};
  }

  template <typename T>
  o2::dataformats::GlobalFwdTrack fwdToTrackPar(const T& track)
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
  o2::dataformats::GlobalFwdTrack fwdToTrackPar(const T& track, const C& cov)
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

  o2::dataformats::GlobalFwdTrack propagateToZMch(const o2::dataformats::GlobalFwdTrack& muon, const double z)
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
  o2::dataformats::GlobalFwdTrack propagateToZMch(const T& muon, const double z)
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

    return propagateToZMch(track, z);
  }

  o2::dataformats::GlobalFwdTrack propagateToZMft(const o2::dataformats::GlobalFwdTrack& mftTrack, const double z)
  {
    o2::dataformats::GlobalFwdTrack fwdtrack{mftTrack};
    fwdtrack.propagateToZ(z, mBzAtMftCenter);
    return fwdtrack;
  }

  template <typename TMFT, typename CMFT>
  o2::dataformats::GlobalFwdTrack propagateToZMft(const TMFT& mftTrack, const CMFT& mftCov, const double z)
  {
    o2::dataformats::GlobalFwdTrack fwdtrack = fwdToTrackPar(mftTrack, mftCov);
    return propagateToZMft(fwdtrack, z);
  }

  // method 0: standard extrapolation
  // method 1: MFT extrapolation using MCH momentum
  // method 2: MCH track extrapolation constrained to the first MFT track point, MFT extrapolation using MCH momentum
  // method 3: MCH track extrapolation constrained to the collision point, MFT extrapolation using MCH momentum
  template <typename TMCH, typename TMFT, typename CMFT, typename C>
  o2::dataformats::GlobalFwdTrack propagateToMatchingPlaneMch(const TMCH& mchTrack, const TMFT& mftTrack, const CMFT& mftTrackCov, const C& collision, const double z, int method)
  {
    if (method == ExtrapolationMethodStandard || method == 1) {
      // simple extrapolation upstream through the absorber
      return propagateToZMch(mchTrack, z);
    }

    if (method == ExtrapolationMethodMftFirstPoint) {
      // extrapolation to the first MFT point and then back to the matching plane
      auto mftTrackPar = fwdToTrackPar(mftTrack, mftTrackCov);
      // std::cout << std::format("[propagateToMatchingPlaneMch] extrapolating to MFT: x={:0.3f} y={:0.3f} z={:0.3f}", mftTrackPar.getX(), mftTrackPar.getY(), mftTrackPar.getZ()) << std::endl;
      auto mchTrackAtMFT = propagateToVertexMch(fwdToTrackPar(mchTrack, mchTrack),
                                                mftTrackPar.getX(), mftTrackPar.getY(), mftTrackPar.getZ(),
                                                mftTrackPar.getSigma2X(), mftTrackPar.getSigma2Y());
      // std::cout << std::format("[propagateToMatchingPlaneMch] extrapolating to z={:0.3f}", z) << std::endl;
      return propagateToZMch(mchTrackAtMFT, z);
    }

    if (method == ExtrapolationMethodVertex) {
      // extrapolation to the vertex and then back to the matching plane
      auto mchTrackAtVertex = VarManager::PropagateMuon(mchTrack, collision, VarManager::kToVertex);
      return propagateToZMch(mchTrackAtVertex, z);
    }

    if (method == ExtrapolationMethodMftDca) {
      // extrapolation to the MFT DCA and then back to the matching plane
      auto mftTrackDCA = propagateToZMft(fwdToTrackPar(mftTrack, mftTrackCov), collision.posZ());
      auto mchTrackAtDCA = propagateToVertexMch(fwdToTrackPar(mchTrack, mchTrack),
                                                mftTrackDCA.getX(), mftTrackDCA.getY(), mftTrackDCA.getZ(),
                                                mftTrackDCA.getSigma2X(), mftTrackDCA.getSigma2Y());
      return propagateToZMch(mchTrackAtDCA, z);
    }

    return {};
  }

  template <typename TMCH, typename TMFT, typename CMFT, typename C>
  o2::dataformats::GlobalFwdTrack propagateToMatchingPlaneMft(const TMCH& mchTrack, const TMFT& mftTrack, const CMFT& mftTrackCov, const C& collision, const double z, int method)
  {
    if (method == ExtrapolationMethodStandard) {
      // extrapolation with MFT tools
      return propagateToZMft(mftTrack, mftTrackCov, z);
    }

    if (method > 0) {
      // extrapolation with MCH tools
      auto mchTrackAtVertex = VarManager::PropagateMuon(mchTrack, collision, VarManager::kToVertex);
      double pMCH = mchTrackAtVertex.getP();
      double px = pMCH * std::sin(o2::constants::math::PIHalf - std::atan(mftTrack.tgl())) * std::cos(mftTrack.phi());
      double py = pMCH * std::sin(o2::constants::math::PIHalf - std::atan(mftTrack.tgl())) * std::sin(mftTrack.phi());
      double pt = std::hypot(px, py);
      double sign = mchTrack.sign();

      o2::dataformats::GlobalFwdTrack track = fwdToTrackPar(mftTrack, mftTrackCov);

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

  o2::dataformats::GlobalFwdTrack propagateToVertexMch(const o2::dataformats::GlobalFwdTrack& muon,
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
  o2::dataformats::GlobalFwdTrack propagateToVertexMch(const TMCH& muon,
                                                       const C& collision)
  {
    return propagateToVertexMch(fwdToTrackPar(muon, muon),
                                collision.posX(),
                                collision.posY(),
                                collision.posZ(),
                                collision.covXX(),
                                collision.covYY());
  }

  o2::dataformats::GlobalFwdTrack propagateToVertexMft(o2::dataformats::GlobalFwdTrack muon,
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
  o2::dataformats::GlobalFwdTrack propagateToVertexMft(const TMFT& muon,
                                                       const C& collision)
  {
    return propagateToVertexMft(fwdToTrackPar(muon),
                                collision.posX(),
                                collision.posY(),
                                collision.posZ(),
                                collision.covXX(),
                                collision.covYY());
  }

  template <typename TMCH, typename TMFT, class C>
  o2::dataformats::GlobalFwdTrack propagateToVertexMft(const TMFT& muon,
                                                       const TMCH& mchTrack,
                                                       const C& collision)
  {
    // extrapolation with MCH tools
    auto mchTrackAtMFT = mExtrap.FwdtoMCH(fwdToTrackPar(mchTrack));
    o2::mch::TrackExtrap::extrapToVertexWithoutBranson(mchTrackAtMFT, muon.z());

    auto muonTrackProp = mExtrap.FwdtoMCH(fwdToTrackPar(muon));

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
  void getMotherParticles(MCP const& mcParticle, std::vector<std::pair<int64_t, int64_t>>& motherParticlesVec)
  {
    const auto& motherParticles = mcParticle.template mothers_as<aod::McParticles>();
    if (motherParticles.empty()) {
      return;
    }

    const auto& motherParticle = motherParticles[0];
    motherParticlesVec.emplace_back(std::make_pair(static_cast<int64_t>(motherParticle.pdgCode()), static_cast<int64_t>(motherParticle.globalIndex())));
    getMotherParticles(motherParticle, motherParticlesVec);
  }

  template <class T>
  std::vector<std::pair<int64_t, int64_t>> getMotherParticles(T const& track)
  {
    std::vector<std::pair<int64_t, int64_t>> result;
    if (!track.has_mcParticle())
      return result;

    const auto& mcParticle = track.mcParticle();
    result.emplace_back(std::make_pair(static_cast<int64_t>(mcParticle.pdgCode()), static_cast<int64_t>(mcParticle.globalIndex())));

    getMotherParticles(mcParticle, result);

    return result;
  }

  template <class TMCH, class TMFTs>
  int getDecayRanking(TMCH const& mchTrack, TMFTs const& mftTracks)
  {
    auto mchMotherParticles = getMotherParticles(mchTrack);

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
  void fillMatchablePairs(CollisionInfo& collisionInfo,
                          TMUON const& muonTracks,
                          TMFT const& mftTracks)
  {
    collisionInfo.matchablePairs.clear();
    for (const auto& muonTrack : muonTracks) {
      // only consider MCH standalone or MCH-MID matches
      if (static_cast<int>(muonTrack.trackType()) <= GlobalTrackTypeMax)
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
          collisionInfo.matchablePairs.emplace_back(std::make_pair(static_cast<int64_t>(muonTrack.globalIndex()),
                                                                   static_cast<int64_t>(mftTrack.globalIndex())));
        }
      }
    }
  }

  template <class TMUON>
  int getTrueMatchIndex(TMUON const& muonTracks,
                        const std::vector<MatchingCandidate>& matchCandidatesVector,
                        const std::vector<std::pair<int64_t, int64_t>>& matchablePairs)
  {
    // find the index of the matching candidate that corresponds to the true match
    // index=1 corresponds to the leading candidate
    // index=0 means no candidate was found that corresponds to the true match
    int trueMatchIndex = 0;
    for (size_t i = 0; i < matchCandidatesVector.size(); i++) {
      auto const& muonTrack = muonTracks.rawIteratorAt(matchCandidatesVector[i].globalTrackId);

      if (isTrueGlobalMatching(muonTrack, matchablePairs)) {
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

    if (std::abs(mchMcParticle.pdgCode()) != kMuonMinus)
      return false;

    return true;
  }

  template <class TMUON, class TMUONS, class TMFTS>
  bool IsMuon(const TMUON& muonTrack,
              TMUONS const& /*muonTracks*/,
              TMFTS const& /*mftTracks*/)
  {
    static constexpr int maxGlobalFwdTrackType = 2;
    if (static_cast<int>(muonTrack.trackType()) >= maxGlobalFwdTrackType)
      return false;

    auto const& mchTrack = muonTrack.template matchMCHTrack_as<TMUONS>();
    auto const& mftTrack = muonTrack.template matchMFTTrack_as<TMFTS>();

    return IsMuon(mchTrack, mftTrack);
  }

  template <class TMUON, class TMUONS, class TMFTS>
  MuonMatchType getMatchType(const TMUON& muonTrack,
                             TMUONS const& /*muonTracks*/,
                             TMFTS const& mftTracks,
                             const std::vector<std::pair<int64_t, int64_t>>& matchablePairs,
                             int ranking)
  {
    if (static_cast<int>(muonTrack.trackType()) > GlobalTrackTypeMax)
      return kMatchTypeUndefined;

    auto const& mchTrack = muonTrack.template matchMCHTrack_as<TMUONS>();

    bool isPairable = isMatchableMch(mchTrack.globalIndex(), matchablePairs);
    bool isTrueMatch = isTrueGlobalMatching(muonTrack, matchablePairs);
    int decayRanking = getDecayRanking(mchTrack, mftTracks);

    MuonMatchType result{kMatchTypeUndefined};
    if (isPairable) {
      if (isTrueMatch) {
        result = (ranking == 1) ? kMatchTypeTrueLeading : kMatchTypeTrueNonLeading;
      } else {
        result = (ranking == 1) ? kMatchTypeWrongLeading : kMatchTypeWrongNonLeading;
      }
    } else if (decayRanking == DecayRankingDirect) {
      result = (ranking == 1) ? kMatchTypeDecayLeading : kMatchTypeDecayNonLeading;
    } else {
      result = (ranking == 1) ? kMatchTypeFakeLeading : kMatchTypeFakeNonLeading;
    }

    return result;
  }

  // for each MCH standalone track, collect the associated matching candidates
  template <class TMUON, class C>
  void getSelectedMuons(const CollisionInfo& collisionInfo,
                        C const& collisions,
                        TMUON const& muonTracks,
                        std::vector<int64_t>& selectedMuons)
  {
    selectedMuons.clear();
    for (const auto& muonTrack : muonTracks) {

      // only consider MCH-MID matches
      if (static_cast<int>(muonTrack.trackType()) != MchMidTrackType) {
        continue;
      }

      // only select MCH-MID tracks from the current collision
      if (!muonTrack.has_collision())
        continue;
      if (static_cast<int64_t>(muonTrack.collisionId()) != collisionInfo.index)
        continue;

      const auto& collision = collisions.rawIteratorAt(muonTrack.collisionId());

      // select MCH tracks with strict quality cuts
      if (!isGoodMuon(muonTrack, collision,
                      cfgMuonTaggingTrackChi2MchUp,
                      cfgMuonTaggingPMchLow,
                      cfgMuonTaggingPtMchLow,
                      {cfgMuonTaggingEtaMchLow, cfgMuonTaggingEtaMchUp},
                      {cfgMuonTaggingRabsLow, cfgMuonTaggingRabsUp},
                      cfgMuonTaggingPdcaUp)) {
        continue;
      }

      // propagate MCH track to the vertex
      auto mchTrackAtVertex = VarManager::PropagateMuon(muonTrack, collision, VarManager::kToVertex);

      // propagate the track from the vertex to the first MFT plane
      const auto& extrapToMFTfirst = propagateToZMch(mchTrackAtVertex, o2::mft::constants::mft::LayerZCoordinate()[0]);
      double rFront = std::sqrt(extrapToMFTfirst.getX() * extrapToMFTfirst.getX() + extrapToMFTfirst.getY() * extrapToMFTfirst.getY());
      double rMinFront = 3.f;
      double rMaxFront = 9.f;
      if (rFront < rMinFront || rFront > rMaxFront)
        continue;

      // propagate the track from the vertex to the last MFT plane
      const auto& extrapToMFTlast = propagateToZMch(mchTrackAtVertex, o2::mft::constants::mft::LayerZCoordinate()[9]);
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
  void getTaggedMuons(const CollisionInfo& collisionInfo,
                      TMUON const& muonTracks,
                      const std::vector<int64_t>& selectedMuons,
                      std::vector<int64_t>& taggedMuons)
  {
    taggedMuons.clear();
    for (const auto& [mchIndex, globalTracksVector] : collisionInfo.matchingCandidates) {

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
      if (chi2diff < cfgMuonTaggingChi2DiffLow)
        continue;

      taggedMuons.emplace_back(mchIndex);
    }
  }

  void getMuonPairs(const CollisionInfo& collisionInfo,
                    std::vector<MuonPair>& muonPairs,
                    std::vector<GlobalMuonPair>& globalMuonPairs)
  {
    // outer loop over muon tracks
    for (const auto& mchIndex1 : collisionInfo.mchTracks) {

      // inner loop over muon tracks
      for (const auto& mchIndex2 : collisionInfo.mchTracks) {
        // avoid double-counting of muon pairs
        if (mchIndex2 <= mchIndex1)
          continue;

        MuonPair muonPair{{collisionInfo.index, mchIndex1}, {collisionInfo.index, mchIndex2}};
        muonPairs.emplace_back(muonPair);
      }
    }

    // outer loop over global muon tracks
    for (const auto& [mchIndex1, matchingCandidates1] : collisionInfo.matchingCandidates) {

      // inner loop over global muon tracks
      for (const auto& [mchIndex2, matchingCandidates2] : collisionInfo.matchingCandidates) {
        // avoid double-counting of muon pairs
        if (mchIndex2 <= mchIndex1)
          continue;

        GlobalMuonPair muonPair{{collisionInfo.index, matchingCandidates1}, {collisionInfo.index, matchingCandidates2}};
        globalMuonPairs.emplace_back(muonPair);
      }
    }
  }

  double getMuMuInvariantMass(const o2::mch::TrackParam& track1, const o2::mch::TrackParam& track2)
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

  double getMuMuInvariantMass(const o2::dataformats::GlobalFwdTrack& track1, const o2::dataformats::GlobalFwdTrack& track2)
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

  template <class EVT, class BC, class TMUON, class TMFTS>
  int getMftMchMatchAttempts(EVT const& collisions,
                             BC const& bcs,
                             TMUON const& mchTrack,
                             TMFTS const& mftTracks)
  {
    if (!mchTrack.has_collision()) {
      return 0;
    }
    const auto& collMch = collisions.rawIteratorAt(mchTrack.collisionId());
    const auto& bcMch = bcs.rawIteratorAt(collMch.bcId());

    int attempts{0};
    for (const auto& mftTrack : mftTracks) {
      if (!mftTrack.has_collision()) {
        continue;
      }

      const auto& collMft = collisions.rawIteratorAt(mftTrack.collisionId());
      const auto& bcMft = bcs.rawIteratorAt(collMft.bcId());

      int64_t deltaBc = static_cast<int64_t>(bcMft.globalBC()) - static_cast<int64_t>(bcMch.globalBC());
      double deltaBcNS = o2::constants::lhc::LHCBunchSpacingNS * deltaBc;
      double deltaTrackTime = mftTrack.trackTime() - mchTrack.trackTime() + deltaBcNS;
      double trackTimeResTot = mftTrack.trackTimeRes() + mchTrack.trackTimeRes();

      if (std::fabs(deltaTrackTime) > trackTimeResTot) {
        continue;
      }
      attempts += 1;
    }

    return attempts;
  }

  template <class EVT, class BC, class TMUON, class TMFT>
  void fillCollisions(EVT const& collisions,
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

    if (collisionIds.empty())
      return;

    for (size_t cid = 1; cid < collisionIds.size() - 1; cid++) {
      const auto& collision = collisions.rawIteratorAt(collisionIds[cid]);
      int64_t collisionIndex = collision.globalIndex();
      auto bc = bcs.rawIteratorAt(collision.bcId());

      // fill collision information for global muon tracks (MFT-MCH-MID matches)
      for (const auto& muonTrack : muonTracks) {
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
          fillMatchablePairs(collisionInfo, muonTracks, mftTracks);
        }

        if (static_cast<int>(muonTrack.trackType()) > GlobalTrackTypeMax) {
          // standalone MCH or MCH-MID tracks
          int64_t mchTrackIndex = muonTrack.globalIndex();
          collisionInfo.mchTracks.push_back(mchTrackIndex);
        } else {
          // global muon tracks (MFT-MCH or MFT-MCH-MID)
          int64_t muonTrackIndex = muonTrack.globalIndex();
          double matchChi2 = muonTrack.chi2MatchMCHMFT() / MatchingDegreesOfFreedom;
          double matchScore = chi2ToScore(muonTrack.chi2MatchMCHMFT(), MatchingDegreesOfFreedom, MatchingScoreChi2Max);
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
              0,
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
              0,
              kMatchTypeUndefined});
          }
        }
      }

      // fill collision information for MFT standalone tracks
      for (const auto& mftTrack : mftTracks) {
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

    // sort the vectors of matching candidates in ascending order based on the matching chi2 value
    auto compareMatchingChi2 = [](const MatchingCandidate& track1, const MatchingCandidate& track2) -> bool {
      return (track1.matchChi2 < track2.matchChi2);
    };

    for (auto collisionInfoIt = collisionInfos.begin(); collisionInfoIt != collisionInfos.end(); ++collisionInfoIt) {
      auto& collisionInfo = collisionInfoIt->second;
      for (auto matchingCandidatesIt = collisionInfo.matchingCandidates.begin(); matchingCandidatesIt != collisionInfo.matchingCandidates.end(); ++matchingCandidatesIt) {
        auto& mchIndex = matchingCandidatesIt->first;
        auto& globalTracksVector = matchingCandidatesIt->second;
        std::sort(globalTracksVector.begin(), globalTracksVector.end(), compareMatchingChi2);

        const auto& mchTrack = muonTracks.rawIteratorAt(mchIndex);
        auto mftMchMatchAttempts = getMftMchMatchAttempts(collisions, bcs, mchTrack, mftTracks);
        int ranking = 1;
        for (auto candidateIt = globalTracksVector.begin(); candidateIt != globalTracksVector.end(); ++candidateIt) {
          auto& candidate = *candidateIt;
          const auto& muonTrack = muonTracks.rawIteratorAt(candidate.globalTrackId);

          candidate.matchRanking = ranking;
          candidate.matchRankingProd = ranking;
          candidate.matchType = getMatchType(muonTrack, muonTracks, mftTracks, collisionInfo.matchablePairs, ranking);
          candidate.mftMchMatchAttempts = mftMchMatchAttempts;
          ranking += 1;
        }
      }
    }
  }

  template <class C, class TMUON, class TMFT>
  void fillMatchingPlotsMc(C const& collision,
                           const CollisionInfo& collisionInfo,
                           TMUON const& muonTracks,
                           TMFT const& mftTracks,
                           const MatchingCandidates& matchingCandidates,
                           const MatchingCandidates& matchingCandidatesProd,
                           const std::vector<std::pair<int64_t, int64_t>>& matchablePairs,
                           double matchingScoreCut,
                           MatchingPlotter* plotter,
                           bool /*verbose*/ = false)
  {
    int mftTrackMult = collisionInfo.mftTracks.size();

    // ====================================
    // Matching candidates hierarchy

    for (const auto& [mchIndex, globalTracksVector] : matchingCandidates) {
      // check if the MCH track belongs to a matchable pair
      bool isPairedMCH = isMatchableMch(static_cast<int64_t>(mchIndex), matchablePairs);

      // get the standalone MCH track
      auto const& mchTrack = muonTracks.rawIteratorAt(mchIndex);
      double mchMom = mchTrack.p();
      double mchPt = mchTrack.pt();

      // MCH track quality flag
      bool isGoodMCH = isGoodGlobalMuon(mchTrack, collision);

      auto matchablePair = getMatchablePairForMch(static_cast<int64_t>(mchIndex), matchablePairs);
      bool hasMatchablePair = matchablePair.has_value();
      int decayRanking = 0;
      int mftTrackType = -1;
      float dchi2 = (globalTracksVector.size() >= MinCandidatesForDeltaChi2) ? globalTracksVector[1].matchChi2 - globalTracksVector[0].matchChi2 : InvalidDeltaChi2;
      if (hasMatchablePair) {
        auto const& pairedMftTrack = mftTracks.rawIteratorAt(matchablePair.value().second);
        mftTrackType = pairedMftTrack.isCA() ? MftTrackTypeCA : MftTrackTypeStandard;
        decayRanking = getDecayRanking(mchTrack, mftTracks);
      }

      // find the index of the matching candidate that corresponds to the true match
      // index=1 corresponds to the leading candidate
      // index=0 means no candidate was found that corresponds to the true match
      int trueMatchIndex = getTrueMatchIndex(muonTracks, globalTracksVector, matchablePairs);
      int trueMatchIndexProd = getTrueMatchIndex(muonTracks, matchingCandidatesProd.at(mchIndex), matchablePairs);

      float mcParticleDz = -1000;
      if (mchTrack.has_mcParticle()) {
        const auto& mchMcParticle = mchTrack.mcParticle();
        mcParticleDz = collision.posZ() - mchMcParticle.vz();
      }

      int matchAttempts = globalTracksVector[0].mftMchMatchAttempts;

      std::get<std::shared_ptr<TH1>>(plotter->fMatchRanking->hist)->Fill(trueMatchIndex);
      std::get<std::shared_ptr<TH2>>(plotter->fMatchRanking->histVsP)->Fill(mchMom, trueMatchIndex);
      std::get<std::shared_ptr<TH2>>(plotter->fMatchRanking->histVsPt)->Fill(mchPt, trueMatchIndex);
      std::get<std::shared_ptr<TH2>>(plotter->fMatchRanking->histVsMcParticleDz)->Fill(mcParticleDz, trueMatchIndex);
      std::get<std::shared_ptr<TH2>>(plotter->fMatchRanking->histVsMftTrackMult)->Fill(mftTrackMult, trueMatchIndex);
      std::get<std::shared_ptr<TH2>>(plotter->fMatchRanking->histVsMatchAttempts)->Fill(matchAttempts, trueMatchIndex);
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
        std::get<std::shared_ptr<TH2>>(plotter->fMatchRankingGoodMCH->histVsMatchAttempts)->Fill(matchAttempts, trueMatchIndex);
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
        std::get<std::shared_ptr<TH2>>(plotter->fMatchRankingPaired->histVsMatchAttempts)->Fill(matchAttempts, trueMatchIndex);
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
        std::get<std::shared_ptr<TH2>>(plotter->fMatchRankingPairedGoodMCH->histVsMatchAttempts)->Fill(matchAttempts, trueMatchIndex);
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

    for (const auto& [mchIndex, globalTracksVector] : matchingCandidates) {
      if (globalTracksVector.size() < 1)
        continue;

      // get the standalone MCH track
      auto const& mchTrack = muonTracks.rawIteratorAt(mchIndex);

      // skip global muon tracks that do not pass the MCH and MFT quality cuts
      if (!isGoodGlobalMuon(mchTrack, collision))
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

    for (const auto& [mchIndex, globalTracksVector] : matchingCandidates) {
      if (globalTracksVector.size() < 1)
        continue;

      int trueMatchIndex = getTrueMatchIndex(muonTracks, globalTracksVector, matchablePairs);

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
    for (const auto& [mchIndex, globalTracksVector] : matchingCandidates) {
      if (globalTracksVector.size() < 1)
        continue;

      // get the leading matching candidate
      auto const& muonTrack = muonTracks.rawIteratorAt(globalTracksVector[0].globalTrackId);
      double matchingScore = globalTracksVector[0].matchScore;

      // get the standalone MCH and MFT tracks
      auto const& mchTrack = muonTracks.rawIteratorAt(mchIndex);
      auto const& mftTrack = muonTrack.template matchMFTTrack_as<TMFT>();

      // skip global muon tracks that do not pass the MCH and MFT quality cuts
      if (!isGoodGlobalMuon(mchTrack, collision))
        continue;
      if (!isGoodMft(mftTrack))
        continue;

      // skip  candidates that do not pass the matching quality cuts
      if (!isGoodGlobalMatching(muonTrack, matchingScore, matchingScoreCut))
        continue;

      // check if the matching candidate is a true one
      bool isTrueMatch = isTrueGlobalMatching(muonTrack, matchablePairs);

      // ---- MC ancestry ----
      auto motherParticles = getMotherParticles(muonTrack);
      int motherPDG = 0;
      if (motherParticles.size() > 1) {
        motherPDG = motherParticles[1].first;
      }
      // fill matching purity plots
      plotter->fMatchingPurityPlotter.fill(mchTrack, isTrueMatch);
      if (configQas.cfgCreatePdgMomHistograms) {
        plotter->fMatchingPurityPlotter.fill(mchTrack, motherPDG, isTrueMatch);
      }
    }

    // ====================================
    // Matching efficiencies

    // outer loop on matchable pairs
    for (const auto& [matchableMchIndex, matchableMftIndex] : matchablePairs) {
      // get the standalone MCH track
      auto const& mchTrack = muonTracks.rawIteratorAt(matchableMchIndex);

      // skip  track pairs that do not pass the MCH and MFT quality cuts
      // we only consider matchable pairs that fulfill the track quality requirements
      if (!isGoodGlobalMuon(mchTrack, collision))
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

          goodMatchFound = isGoodGlobalMatching(muonTrack, matchingScore, matchingScoreCut);
          isTrueMatch = (mftIndex == matchableMftIndex);
        }
      }

      if (cfgQaMatchingAodDebug > 0 && goodMatchFound && isTrueMatch) {
        LOGF(info,
             "[good&true] mchId=%lld trackType=%d p=%.3f pt=%.3f eta=%.3f phi=%.3f",
             static_cast<int64_t>(mchTrack.globalIndex()),
             static_cast<int>(mchTrack.trackType()),
             mchTrack.p(),
             mchTrack.pt(),
             mchTrack.eta(),
             mchTrack.phi());
      }

      // ---- MC ancestry ----
      auto motherParticles = getMotherParticles(mchTrack);
      int motherPDG = 0;
      if (motherParticles.size() > 1) {
        motherPDG = motherParticles[1].first;
      }

      // fill matching efficiency plots
      plotter->fPairingEfficiencyPlotter.fill(mchTrack, goodMatchFound);
      if (configQas.cfgCreatePdgMomHistograms) {
        plotter->fPairingEfficiencyPlotter.fill(mchTrack, motherPDG, goodMatchFound);
      }
      plotter->fMatchingEfficiencyPlotter.fill(mchTrack, (goodMatchFound && isTrueMatch));
      if (configQas.cfgCreatePdgMomHistograms) {
        plotter->fMatchingEfficiencyPlotter.fill(mchTrack, motherPDG, (goodMatchFound && isTrueMatch));
      }
      plotter->fFakeMatchingEfficiencyPlotter.fill(mchTrack, (goodMatchFound && !isTrueMatch));
      if (configQas.cfgCreatePdgMomHistograms) {
        plotter->fFakeMatchingEfficiencyPlotter.fill(mchTrack, motherPDG, (goodMatchFound && !isTrueMatch));
      }
    }
  }

  template <class C, class TMUON, class TMFT>
  void fillDimuonPlotsMc(const CollisionInfo& collisionInfo,
                         C const& collisions,
                         TMUON const& muonTracks,
                         TMFT const& /*mftTracks*/)
  {
    std::vector<MuonPair> muonPairs;
    std::vector<GlobalMuonPair> globalMuonPairs;

    getMuonPairs(collisionInfo, muonPairs, globalMuonPairs);

    for (const auto& [muon1, muon2] : muonPairs) {
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

      bool goodMuonTracks = (isGoodMuon(muonTrack1, collision) && isGoodMuon(muonTrack2, collision));

      if (goodMuonTracks) {
        double mass = getMuMuInvariantMass(propagateToVertexMch(muonTrack1, collision),
                                           propagateToVertexMch(muonTrack2, collision));
        registryDimuon.get<TH1>(HIST("dimuon/invariantMass_MuonKine_MuonCuts"))->Fill(mass);
      }
    }

    for (const auto& [muon1, muon2] : globalMuonPairs) {
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

      bool goodGlobalMuonTracks = (isGoodGlobalMuon(mchTrack1, collision) && isGoodGlobalMuon(mchTrack2, collision));
      if (!goodGlobalMuonTracks) {
        continue;
      }

      bool goodGlobalMuonMatches = (isGoodGlobalMatching(muonTrack1, matchScore1) && isGoodGlobalMatching(muonTrack2, matchScore2));

      double massMCH = getMuMuInvariantMass(propagateToVertexMch(mchTrack1, collision),
                                            propagateToVertexMch(mchTrack2, collision));
      double mass = getMuMuInvariantMass(propagateToVertexMch(muonTrack1, collision),
                                         propagateToVertexMch(muonTrack2, collision));
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

  template <class C, class BC, class TMUON, class TMFT, class CMFT>
  void runChi2Matching(C const& collisions,
                       BC const& bcs,
                       TMUON const& muonTracks,
                       TMFT const& mftTracks,
                       CMFT const& mftCovs,
                       std::string funcName,
                       float matchingPlaneZ,
                       int extrapMethod,
                       const std::vector<std::pair<int64_t, int64_t>>& matchablePairs,
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

    for (const auto& [mchIndex, globalTracksVector] : matchingCandidates) {
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
        auto mftTrackProp = fwdToTrackPar(mftTrack, mftTrackCov);
        auto mchTrackProp = fwdToTrackPar(mchTrack, mchTrack);

        if (matchingPlaneZ < 0.) {
          mftTrackProp = propagateToMatchingPlaneMft(mchTrack, mftTrack, mftTrackCov, collision, matchingPlaneZ, extrapMethod);
          mchTrackProp = propagateToMatchingPlaneMch(mchTrack, mftTrack, mftTrackCov, collision, matchingPlaneZ, extrapMethod);
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

    // sort the vectors of matching candidates in ascending order based on the matching chi2 value
    auto compareMatchingChi2 = [](const MatchingCandidate& track1, const MatchingCandidate& track2) -> bool {
      return (track1.matchChi2 < track2.matchChi2);
    };

    for (auto matchingCandidatesIt = newMatchingCandidates.begin(); matchingCandidatesIt != newMatchingCandidates.end(); ++matchingCandidatesIt) {
      auto& mchIndex = matchingCandidatesIt->first;
      auto& globalTracksVector = matchingCandidatesIt->second;
      std::sort(globalTracksVector.begin(), globalTracksVector.end(), compareMatchingChi2);

      const auto& mchTrack = muonTracks.rawIteratorAt(mchIndex);
      auto mftMchMatchAttempts = getMftMchMatchAttempts(collisions, bcs, mchTrack, mftTracks);
      int ranking = 1;
      for (auto candidateIt = globalTracksVector.begin(); candidateIt != globalTracksVector.end(); ++candidateIt) {
        auto& candidate = *candidateIt;
        const auto& muonTrack = muonTracks.rawIteratorAt(candidate.globalTrackId);

        candidate.matchRanking = ranking;
        candidate.matchType = getMatchType(muonTrack, muonTracks, mftTracks, matchablePairs, ranking);
        candidate.mftMchMatchAttempts = mftMchMatchAttempts;
        ranking += 1;
      }
    }
  }

  template <class C, class BC, class TMUON, class TMFT, class CMFT>
  void runChi2Matching(C const& collisions,
                       BC const& bcs,
                       TMUON const& muonTracks,
                       TMFT const& mftTracks,
                       CMFT const& mftCovs,
                       std::string label,
                       const std::vector<std::pair<int64_t, int64_t>>& matchablePairs,
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

    runChi2Matching(collisions, bcs, muonTracks, mftTracks, mftCovs, funcName, matchingPlaneZ, extrapMethod, matchablePairs, matchingCandidates, newMatchingCandidates);
  }

  template <class C, class BC, class TMUON, class TMFT, class CMFT>
  void runMlMatching(C const& collisions,
                     BC const& bcs,
                     TMUON const& muonTracks,
                     TMFT const& mftTracks,
                     CMFT const& mftCovs,
                     std::string label,
                     const std::vector<std::pair<int64_t, int64_t>>& matchablePairs,
                     const MatchingCandidates& matchingCandidates,
                     MatchingCandidates& newMatchingCandidates)
  {
    newMatchingCandidates.clear();
    auto mlIter = matchingMlResponses.find(label);
    if (mlIter == matchingMlResponses.end())
      return;

    auto& mlResponse = mlIter->second;
    for (const auto& [mchIndex, globalTracksVector] : matchingCandidates) {
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
        auto mftTrackProp = fwdToTrackPar(mftTrack, mftTrackCov);
        auto mchTrackProp = fwdToTrackPar(mchTrack, mchTrack);

        // extrapolate to the matching plane
        auto matchingPlaneZ = matchingPlanesZ[label];
        if (matchingPlaneZ < 0.) {
          mftTrackProp = propagateToZMft(mftTrackProp, matchingPlaneZ);
          mchTrackProp = propagateToZMch(mchTrackProp, matchingPlaneZ);
        }

        // run the ML model
        std::vector<float> output;
        std::vector<float> inputML = mlResponse.getInputFeaturesGlob(muonTrack, mchTrackProp, mftTrackProp, collision);
        mlResponse.isSelectedMl(inputML, 0, output);
        float matchScore = output[0];
        float matchChi2Prod = muonTrack.chi2MatchMCHMFT() / MatchingDegreesOfFreedom;
        float matchScoreProd = chi2ToScore(muonTrack.chi2MatchMCHMFT(), MatchingDegreesOfFreedom, MatchingScoreChi2Max);
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

    for (auto matchingCandidatesIt = newMatchingCandidates.begin(); matchingCandidatesIt != newMatchingCandidates.end(); ++matchingCandidatesIt) {
      auto& mchIndex = matchingCandidatesIt->first;
      auto& globalTracksVector = matchingCandidatesIt->second;
      std::sort(globalTracksVector.begin(), globalTracksVector.end(), compareMatchingScore);

      const auto& mchTrack = muonTracks.rawIteratorAt(mchIndex);
      auto mftMchMatchAttempts = getMftMchMatchAttempts(collisions, bcs, mchTrack, mftTracks);
      int ranking = 1;
      for (auto candidateIt = globalTracksVector.begin(); candidateIt != globalTracksVector.end(); ++candidateIt) {
        auto& candidate = *candidateIt;
        const auto& muonTrack = muonTracks.rawIteratorAt(candidate.globalTrackId);

        candidate.matchRanking = ranking;
        candidate.matchType = getMatchType(muonTrack, muonTracks, mftTracks, matchablePairs, ranking);
        candidate.mftMchMatchAttempts = mftMchMatchAttempts;
        ranking += 1;
      }
    }
  }

  template <class C, class BC, class TMUON, class TMFT, class CMFT>
  void processCollisionMc(const CollisionInfo& collisionInfo,
                          C const& collisions,
                          BC const& bcs,
                          TMUON const& muonTracks,
                          TMFT const& mftTracks,
                          CMFT const& mftCovs)
  {
    auto collision = collisions.rawIteratorAt(collisionInfo.index);

    registry.get<TH1>(HIST("tracksMultiplicityMFT"))->Fill(collisionInfo.mftTracks.size());
    registry.get<TH1>(HIST("tracksMultiplicityMCH"))->Fill(collisionInfo.mchTracks.size());

    // Chi2-based matching analysis
    fillMatchingPlotsMc(collision, collisionInfo, muonTracks, mftTracks, collisionInfo.matchingCandidates, collisionInfo.matchingCandidates, collisionInfo.matchablePairs, cfgMatchingChi2ScoreMftMchLow, fChi2MatchingPlotter.get(), false);
    for (const auto& [label, func] : matchingChi2Functions) {
      MatchingCandidates matchingCandidates;
      runChi2Matching(collisions, bcs, muonTracks, mftTracks, mftCovs, label, collisionInfo.matchablePairs, collisionInfo.matchingCandidates, matchingCandidates);

      auto* plotter = fMatchingPlotters.at(label).get();
      double matchingScoreCut = matchingScoreCuts.at(label);

      fillMatchingPlotsMc(collision, collisionInfo, muonTracks, mftTracks, matchingCandidates, collisionInfo.matchingCandidates, collisionInfo.matchablePairs, matchingScoreCut, plotter, false);
    }

    // ML-based matching analysis
    for (const auto& [label, mlResponse] : matchingMlResponses) {
      MatchingCandidates matchingCandidates;
      runMlMatching(collisions, bcs, muonTracks, mftTracks, mftCovs, label, collisionInfo.matchablePairs, collisionInfo.matchingCandidates, matchingCandidates);

      auto* plotter = fMatchingPlotters.at(label).get();
      double matchingScoreCut = matchingScoreCuts.at(label);

      fillMatchingPlotsMc(collision, collisionInfo, muonTracks, mftTracks, matchingCandidates, collisionInfo.matchingCandidates, collisionInfo.matchablePairs, matchingScoreCut, plotter);
    }

    // Muons tagging
    for (const auto& [mchIndex, mftIndex] : collisionInfo.matchablePairs) {
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
      auto mchTrackProp = propagateToZMch(mchTrackAtVertex, z);
      auto mftTrackProp = propagateToZMft(fwdToTrackPar(mftTrack, mftTrackCov), z);

      registry.get<TH2>(HIST("matching/MC/pairedMCHTracksAtMFT"))->Fill(mchTrackProp.getX(), mchTrackProp.getY());
      registry.get<TH2>(HIST("matching/MC/pairedMFTTracksAtMFT"))->Fill(mftTrackProp.getX(), mftTrackProp.getY());
    }

    std::vector<int64_t> selectedMuons;
    getSelectedMuons(collisionInfo, collisions, muonTracks, selectedMuons);

    MatchingCandidates selectedMatchingCandidates;
    for (const auto& [mchIndex, globalTracksVector] : collisionInfo.matchingCandidates) {
      if (std::find(selectedMuons.begin(), selectedMuons.end(), mchIndex) != selectedMuons.end()) {
        selectedMatchingCandidates[mchIndex] = globalTracksVector;
      }
    }
    fillMatchingPlotsMc(collision, collisionInfo, muonTracks, mftTracks, selectedMatchingCandidates, collisionInfo.matchingCandidates, collisionInfo.matchablePairs, cfgMatchingChi2ScoreMftMchLow, fSelectedMuonsMatchingPlotter.get());

    std::vector<int64_t> taggedMuons;
    getTaggedMuons(collisionInfo, muonTracks, selectedMuons, taggedMuons);

    MatchingCandidates taggedMatchingCandidates;
    for (const auto& [mchIndex, globalTracksVector] : collisionInfo.matchingCandidates) {
      if (std::find(taggedMuons.begin(), taggedMuons.end(), mchIndex) != taggedMuons.end()) {
        taggedMatchingCandidates[mchIndex] = globalTracksVector;
      }
    }
    fillMatchingPlotsMc(collision, collisionInfo, muonTracks, mftTracks, taggedMatchingCandidates, collisionInfo.matchingCandidates, collisionInfo.matchablePairs, cfgMatchingChi2ScoreMftMchLow, fTaggedMuonsMatchingPlotter.get());

    // Di-muon analysis
    fillDimuonPlotsMc(collisionInfo, collisions, muonTracks, mftTracks);
  }

  template <class TCOLLISION, class TMUON>
  void fillQaMatchingAodTablesForCollision(TCOLLISION const& collision,
                                           TMUON const& muonTracks,
                                           const MatchingCandidates& matchingCandidates,
                                           int8_t matchLabel,
                                           int32_t reducedEventId)
  {
    for (const auto& [mchIndex, candidates] : matchingCandidates) {
      if (candidates.empty()) {
        continue;
      }

      const auto& mchTrack = muonTracks.rawIteratorAt(mchIndex);
      if (!isGoodGlobalMuon(mchTrack, collision)) {
        continue;
      }

      for (const auto& candidate : candidates) {
        const auto& candidateTrack = muonTracks.rawIteratorAt(candidate.globalTrackId);
        auto candidateTrackAtVertex = VarManager::PropagateMuon(candidateTrack, collision, VarManager::kToVertex);
        qaMatchingCandidates(
          reducedEventId,
          matchLabel,
          mchIndex,
          static_cast<float>(candidateTrack.p()),
          static_cast<float>(candidateTrack.pt()),
          static_cast<float>(candidateTrack.eta()),
          static_cast<float>(candidateTrack.phi()),
          static_cast<int8_t>(candidate.matchType),
          static_cast<float>(candidate.matchScore),
          static_cast<int32_t>(candidate.matchRanking),
          static_cast<float>(candidateTrackAtVertex.getX()),
          static_cast<float>(candidateTrackAtVertex.getY()),
          static_cast<float>(candidateTrackAtVertex.getZ()),
          static_cast<float>(candidateTrackAtVertex.getPx()),
          static_cast<float>(candidateTrackAtVertex.getPy()),
          static_cast<float>(candidateTrackAtVertex.getPz()));
      }
    }
  }

  template <class TCOLLISION>
  void fillQaMatchingAodEventForCollision(const CollisionInfo& collisionInfo,
                                          TCOLLISION const& collision,
                                          int32_t reducedEventId,
                                          int& debugCounter)
  {
    int32_t mftMultiplicity = static_cast<int32_t>(collisionInfo.mftTracks.size());
    qaMatchingEvents(
      mftMultiplicity,
      static_cast<float>(collision.posX()),
      static_cast<float>(collision.posY()),
      static_cast<float>(collision.posZ()));

    if (cfgQaMatchingAodDebug > 0 && debugCounter < cfgQaMatchingAodDebug) {
      LOGF(info, "[AO2D] reducedEvent=%", reducedEventId);
      debugCounter += 1;
    }
  }

  template <class TCOLLISIONS, class TCOLLISION, class TMUON, class TMFT, class TBC>
  void fillQaMatchingMchTracksForCollision(const CollisionInfo& collisionInfo,
                                           TCOLLISIONS const& collisions,
                                           TCOLLISION const& collision,
                                           TMUON const& muonTracks,
                                           TMFT const& mftTracks,
                                           TBC const& bcs,
                                           int32_t reducedEventId)
  {
    std::vector<int64_t> mchIds;
    for (const auto& mchIndex : collisionInfo.mchTracks) {
      if (std::find(mchIds.begin(), mchIds.end(), mchIndex) == mchIds.end()) {
        mchIds.emplace_back(mchIndex);
      }
    }
    for (const auto& [mchIndex, candidates] : collisionInfo.matchingCandidates) {
      (void)candidates;
      if (std::find(mchIds.begin(), mchIds.end(), mchIndex) == mchIds.end()) {
        mchIds.emplace_back(mchIndex);
      }
    }

    for (const auto& mchIndex : mchIds) {
      auto const& mchTrack = muonTracks.rawIteratorAt(mchIndex);
      int mftMchMatchAttempts = getMftMchMatchAttempts(collisions, bcs, mchTrack, mftTracks);
      auto mchTrackAtVertex = VarManager::PropagateMuon(mchTrack, collision, VarManager::kToVertex);
      qaMatchingMCHTrack(
        reducedEventId,
        mchIndex,
        static_cast<int8_t>(mchTrack.trackType()),
        static_cast<float>(mchTrack.p()),
        static_cast<float>(mchTrack.pt()),
        static_cast<float>(mchTrack.eta()),
        static_cast<float>(mchTrack.phi()),
        static_cast<int32_t>(mftMchMatchAttempts),
        static_cast<float>(mchTrackAtVertex.getX()),
        static_cast<float>(mchTrackAtVertex.getY()),
        static_cast<float>(mchTrackAtVertex.getZ()),
        static_cast<float>(mchTrackAtVertex.getPx()),
        static_cast<float>(mchTrackAtVertex.getPy()),
        static_cast<float>(mchTrackAtVertex.getPz()));
    }
  }

  void processQAMC(MyEvents const& collisions,
                   aod::BCsWithTimestamps const& bcs,
                   MyMuonsMC const& muonTracks,
                   MyMFTsMC const& mftTracks,
                   MyMFTCovariances const& mftCovs,
                   aod::McParticles const& /*mcParticles*/)
  {
    auto bc = bcs.begin();
    initCcdb(bc);

    for (const auto& muon : muonTracks) {
      registry.get<TH1>(HIST("nTracksPerType"))->Fill(static_cast<int>(muon.trackType()));
    }

    fillCollisions(collisions, bcs, muonTracks, mftTracks, fCollisionInfos);

    mftTrackCovs.clear();
    for (const auto& mftTrackCov : mftCovs) {
      mftTrackCovs[mftTrackCov.matchMFTTrackId()] = mftTrackCov.globalIndex();
    }

    std::unordered_map<int64_t, int32_t> reducedEventIds;
    int32_t reducedEventCounter = 0;
    for (auto const& [collisionIndex, collisionInfo] : fCollisionInfos) {
      reducedEventIds.emplace(collisionInfo.index, reducedEventCounter);
      reducedEventCounter += 1;
    }

    int debugCounter = 0;
    for (auto const& [collisionIndex, collisionInfo] : fCollisionInfos) {
      auto it = reducedEventIds.find(collisionInfo.index);
      if (it == reducedEventIds.end()) {
        continue;
      }
      int32_t reducedEventId = it->second;
      auto collision = collisions.rawIteratorAt(collisionInfo.index);
      fillQaMatchingAodEventForCollision(collisionInfo, collision, reducedEventId, debugCounter);
      fillQaMatchingMchTracksForCollision(collisionInfo, collisions, collision, muonTracks, mftTracks, bcs, reducedEventId);
    }

    struct AodLabel {
      const char* name;
      int8_t id;
    };
    std::array<AodLabel, 3> aodLabels{{{"ProdAll", 0}, {"MatchXYPhiTanl", 1}, {"MatchXYPhiTanlMom", 2}}};
    for (const auto& aodLabel : aodLabels) {
      if (matchingChi2Functions.find(aodLabel.name) == matchingChi2Functions.end()) {
        LOGF(warn, "[AO2D] Chi2 label not found: %s", aodLabel.name);
        continue;
      }
      debugCounter = 0;
      for (auto const& [collisionIndex, collisionInfo] : fCollisionInfos) {
        auto it = reducedEventIds.find(collisionInfo.index);
        if (it == reducedEventIds.end()) {
          continue;
        }
        int32_t reducedEventId = it->second;
        MatchingCandidates matchingCandidates;
        runChi2Matching(collisions, bcs, muonTracks, mftTracks, mftCovs, aodLabel.name, collisionInfo.matchablePairs, collisionInfo.matchingCandidates, matchingCandidates);
        auto collision = collisions.rawIteratorAt(collisionInfo.index);
        fillQaMatchingAodTablesForCollision(collision, muonTracks, matchingCandidates, aodLabel.id, reducedEventId);
      }
    }

    for (auto const& [collisionIndex, collisionInfo] : fCollisionInfos) {
      processCollisionMc(collisionInfo, collisions, bcs, muonTracks, mftTracks, mftCovs);
    }
  }

  PROCESS_SWITCH(QaMatching, processQAMC, "processQAMC", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<QaMatching>(cfgc)};
};
