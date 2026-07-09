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
/// \file global-muon-matcher.cxx // o2-linter: disable=name/file-cpp,name/workflow-file (legacy workflow executable name)
/// \brief Task for analysis MFT-MCH muon matching
/// \author Andrea Ferrero
///
#include "PWGDQ/Core/MuonMatchingMlResponse.h"
#include "PWGDQ/Core/VarManager.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/fwdtrackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FwdTrackReAlignTables.h"
#include "Tools/ML/MlResponse.h"

#include <CCDB/BasicCCDBManager.h>
#include <CCDB/CcdbApi.h>
#include <CommonConstants/LHCConstants.h>
#include <CommonConstants/MathConstants.h>
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
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>
#include <GPU/GPUROOTCartesianFwd.h>
#include <GlobalTracking/MatchGlobalFwd.h>
#include <MCHBase/TrackerParam.h>
#include <MCHGeometryTransformer/Transformations.h>
#include <MCHTracking/Track.h>
#include <MCHTracking/TrackExtrap.h>
#include <MCHTracking/TrackFitter.h>
#include <MCHTracking/TrackParam.h>
#include <MFTTracking/Constants.h>
#include <MathUtils/Cartesian.h>
#include <ReconstructionDataFormats/TrackFwd.h>

#include <Math/MatrixFunctions.h>
#include <Math/SMatrix.h>
#include <Math/SVector.h>
#include <TGeoGlobalMagField.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <functional>
#include <iterator>
#include <map>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::aod;

namespace o2::aod
{
namespace globalmuonmatching
{
DECLARE_SOA_COLUMN(IsTagged, isTagged, bool);                  //! Whether the MCH(-MID) track passes tagging cuts
DECLARE_SOA_COLUMN(MatchRanking, matchRanking, int32_t);       //! Match candidate ranking (-1 for base MCH entries)
DECLARE_SOA_COLUMN(MixedGroupIndex, mixedGroupIndex, int32_t); //! Mixed-event group index (-1 for same-event candidates)
} // namespace globalmuonmatching

DECLARE_SOA_TABLE(GmmCandFwdTrkExtras, "AOD", "GMMCANDEXTRA", //! Extra info joinable to FwdTracksReAlign
                  globalmuonmatching::IsTagged,
                  globalmuonmatching::MatchRanking,
                  globalmuonmatching::MixedGroupIndex);
} // namespace o2::aod

using MyEvents = soa::Join<aod::Collisions, aod::EvSels>;
using MyMuons = soa::Join<aod::FwdTracks, aod::FwdTracksCov>;
using MyMFTs = aod::MFTTracks;
using MyMFTCovariances = aod::MFTTracksCov;

using SMatrix55Sym = o2::track::SMatrix55Sym;
using SMatrix55Std = o2::track::SMatrix55Std;
using SMatrix5 = o2::track::SMatrix5;

constexpr std::array<int, 10> NDetElemCh = {4, 4, 4, 4, 18, 18, 26, 26, 26, 26};
constexpr std::array<int, 11> SNDetElemCh = {0, 4, 8, 12, 16, 34, 52, 78, 104, 130, 156};

// compute minimum difference between azimuthal angles
static float getDeltaPhi(float phi1, float phi2)
{
  return RecoDecay::constrainAngle(phi1 - phi2, -o2::constants::math::PI);
}

struct GlobalMuonMatching {

  static constexpr int GlobalTrackTypeMax = 2;
  static constexpr int MchMidTrackType = 3;
  static constexpr int NMchChambers = 10;
  static constexpr int MchDetElemNumberingBase = 100;
  static constexpr int NMchDetElems = 156;
  static constexpr int MinRemovableTrackClusters = 10;
  static constexpr int ThetaAbsBoundaryDeg = 3;
  static constexpr double SlopeResolutionZ = 535.;
  static constexpr float MatchingPlaneDefaultZ = -77.5;

  struct MatchingCandidate {
    int64_t muonTrackId{-1};
    int64_t mftTrackId{-1};
    double matchScore{-1};
    double matchChi2{-1};
    int matchRanking{-1};
    int32_t mixedGroupIndex{-1};
  };

  struct MchTrackInfo {
    int nMatchAttempts{-1};
    bool isTagged{false};
    // vector of MFT-MCH matching candidates
    std::vector<MatchingCandidate> matchingCandidates;
    // vector of vectors of MFT-MCH matching candidates from mixed events
    std::vector<std::vector<MatchingCandidate>> mixedMatchingCandidates;
  };

  ////   Variables for selecting tagged muons
  struct : ConfigurableGroup {
    Configurable<int> cfgMuonTaggingNCrossedMftPlanesLow{"cfgMuonTaggingNCrossedMftPlanesLow", 5, ""};
    Configurable<float> cfgMuonTaggingTrackChi2MchUp{"cfgMuonTaggingTrackChi2MchUp", 5.f, ""};
    Configurable<float> cfgMuonTaggingPMchLow{"cfgMuonTaggingPMchLow", 0.0f, ""};
    Configurable<float> cfgMuonTaggingPtMchLow{"cfgMuonTaggingPtMchLow", 0.7f, ""};
    Configurable<float> cfgMuonTaggingEtaMchLow{"cfgMuonTaggingEtaMchLow", -3.6f, ""};
    Configurable<float> cfgMuonTaggingEtaMchUp{"cfgMuonTaggingEtaMchUp", -2.5f, ""};
    Configurable<float> cfgMuonTaggingRabsLow{"cfgMuonTaggingRabsLow", 17.6f, ""};
    Configurable<float> cfgMuonTaggingRabsUp{"cfgMuonTaggingRabsUp", 89.5f, ""};
    Configurable<float> cfgMuonTaggingPdcaUp{"cfgMuonTaggingPdcaUp", 4.f, ""};
    Configurable<float> cfgMuonTaggingRadiusAtMftFrontLow{"cfgMuonTaggingRadiusAtMftFrontLow", 3.f, ""};
    Configurable<float> cfgMuonTaggingRadiusAtMftFrontUp{"cfgMuonTaggingRadiusAtMftFrontUp", 9.f, ""};
    Configurable<float> cfgMuonTaggingRadiusAtMftBackLow{"cfgMuonTaggingRadiusAtMftBackLow", 5.f, ""};
    Configurable<float> cfgMuonTaggingRadiusAtMftBackUp{"cfgMuonTaggingRadiusAtMftBackUp", 12.f, ""};
  } configMuonTagging;

  ////   Variables for MCH realignment
  struct : ConfigurableGroup {
    Configurable<bool> cfgEnableMCHRealign{"cfgEnableMCHRealign", true, "Enable re-alignment of MCH clusters and tracks"};
    Configurable<std::string> cfgGeoRefPath{"cfgGeoRefPath", "GLO/Config/GeometryAligned", "Path of the reference geometry file"};
    Configurable<std::string> cfgGeoNewPath{"cfgGeoNewPath", "GLO/Config/GeometryAligned", "Path of the new geometry file"};
    Configurable<int64_t> cfgCcdbNoLaterThanRef{"cfgCcdbNoLaterThanRef", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object of reference basis"};
    Configurable<int64_t> cfgCcdbNoLaterThanNew{"cfgCcdbNoLaterThanNew", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object of new basis"};
    Configurable<double> cfgChamberResolutionX{"cfgChamberResolutionX", 0.04, "Chamber resolution along X configuration for refit"}; // 0.4cm pp, 0.2cm PbPb
    Configurable<double> cfgChamberResolutionY{"cfgChamberResolutionY", 0.04, "Chamber resolution along Y configuration for refit"}; // 0.4cm pp, 0.2cm PbPb
    Configurable<double> cfgSigmaCutImprove{"cfgSigmaCutImprove", 6., "Sigma cut for track improvement"};                            // 6 for pp, 4 for PbPb
  } configMchRealign;

  ////   Variables for MFT alignment corrections
  struct : ConfigurableGroup {
    Configurable<bool> cfgEnableMftAlignmentCorrections{"cfgEnableMftAlignmentCorrections", true, "Enable alignment corrections for the MFT tracks"};
    // slope corrections
    // Configurable<float> cfgMFTAlignmentCorrXSlopeTop{"cfgMFTAlignmentCorrXSlopeTop", (-0.0006696 - 0.0005621) / 2.f, "MFT X slope correction - top half"};
    // Configurable<float> cfgMFTAlignmentCorrXSlopeBottom{"cfgMFTAlignmentCorrXSlopeBottom", (0.00105 + 0.001007) / 2.f, "MFT X slope correction - bottom half"};
    // Configurable<float> cfgMFTAlignmentCorrYSlopeTop{"cfgMFTAlignmentCorrYSlopeTop", (-0.002299 - 0.002442) / 2.f, "MFT Y slope correction - top half"};
    // Configurable<float> cfgMFTAlignmentCorrYSlopeBottom{"cfgMFTAlignmentCorrYSlopeBottom", (-0.0005339 - 0.0006921) / 2.f, "MFT Y slope correction - bottom half"};
    Configurable<float> cfgMFTAlignmentCorrXSlopeTop{"cfgMFTAlignmentCorrXSlopeTop", 0.f, "MFT X slope correction - top half"};
    Configurable<float> cfgMFTAlignmentCorrXSlopeBottom{"cfgMFTAlignmentCorrXSlopeBottom", 0.f, "MFT X slope correction - bottom half"};
    Configurable<float> cfgMFTAlignmentCorrYSlopeTop{"cfgMFTAlignmentCorrYSlopeTop", 0.f, "MFT Y slope correction - top half"};
    Configurable<float> cfgMFTAlignmentCorrYSlopeBottom{"cfgMFTAlignmentCorrYSlopeBottom", 0.f, "MFT Y slope correction - bottom half"};
    // offset corrections
    Configurable<float> cfgMFTAlignmentCorrXOffsetTop{"cfgMFTAlignmentCorrXOffsetTop", 0.f, "MFT X offset correction - top half"};
    Configurable<float> cfgMFTAlignmentCorrXOffsetBottom{"cfgMFTAlignmentCorrXOffsetBottom", 0.f, "MFT X offset correction - bottom half"};
    Configurable<float> cfgMFTAlignmentCorrYOffsetTop{"cfgMFTAlignmentCorrYOffsetTop", 0.f, "MFT Y offset correction - top half"};
    Configurable<float> cfgMFTAlignmentCorrYOffsetBottom{"cfgMFTAlignmentCorrYOffsetBottom", 0.f, "MFT Y offset correction - bottom half"};
  } configMftAlignmentCorrections;

  // Variables for CCDB objects access and retrieval
  struct : ConfigurableGroup {
    Configurable<std::string> cfgCcdbUrl{"cfgCcdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
    Configurable<int64_t> cfgCcdbNoLaterThan{"cfgCcdbNoLaterThan", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
    Configurable<std::string> cfgGrpPath{"cfgGrpPath", "GLO/GRP/GRP", "Path of the grp file"};
    Configurable<std::string> cfgGeoPath{"cfgGeoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
    Configurable<std::string> cfgGrpMagPath{"cfgGrpMagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  } configCcdb;

  // Matching strategy for the *custom* matches (production baseline is always computed).
  // 0 = chi2 (runChi2Matching), 1 = ML (runMlMatching)
  struct : ConfigurableGroup {
    Configurable<int> cfgCustomMatchingStrategy{"cfgCustomMatchingStrategy", 0, "0=chi2, 1=ML for custom matches"};
    Configurable<bool> cfgIncludeGlobalMuonsInFwdTracks{"cfgIncludeGlobalMuonsInFwdTracks", false, "Include MFT-MCH-MID global muons in GMMCANDTRK table"};
    Configurable<int> cfgMaxCandidatesPerMchTrack{"cfgMaxCandidatesPerMchTrack", -1, "Maximum number of match candidates stored per MCH track (-1: no limit)"};
    Configurable<bool> cfgMatchAllTracks{"cfgMatchAllTracks", false, "If true the matching is performed considering all the MFT tracks for which the covariances are available; if false the matching is performed considering only the global forward tracks stored at production"};
  } configMatching;

  struct : ConfigurableGroup {
    Configurable<int> cfgMixingDepth{"cfgMixingDepth", -1, "Maximum number of mixed candidate groups per MCH track (-1: no limit)"};
    Configurable<int64_t> cfgMinDeltaBc{"cfgMinDeltaBc", 3564, "Minimum DeltaBc between mixed collisions"};
    Configurable<float> cfgMaxDeltaPhi{"cfgMaxDeltaPhi", static_cast<float>(o2::constants::math::PI / 10), "Maximum DelptaPhi between mixed MCH tracks (rad)"};
    Configurable<float> cfgMaxDeltaR{"cfgMaxDeltaR", 10.f, "Maximum DeltaR between mixed MCH tracks"};
    Configurable<float> cfgMaxDeltaAttempts{"cfgMaxDeltaAttempts", 0.1f, "Maximum relative difference in match attempts"};
    Configurable<float> cfgMaxDeltaZ{"cfgMaxDeltaZ", 1.f, "Maximum deltaZ between mixed collisions"};
  } configEventMixing;

  double mBzAtMftCenter{0};

  using MatchingFunc = std::function<std::tuple<double, int>(const o2::track::TrackParCovFwd& mchtrack, const o2::track::TrackParCovFwd& mfttrack)>;
  std::map<std::string, MatchingFunc> mMatchingFunctionMap; ///< MFT-MCH Matching function

  // Chi2 matching interface (single configurable method)
  struct : ConfigurableGroup {
    Configurable<std::string> cfgChi2FunctionLabel{"cfgChi2FunctionLabel", std::string{"ProdAll"}, "Text label identifying the chi2 matching method"};
    Configurable<std::string> cfgChi2FunctionName{"cfgChi2FunctionName", std::string{"prod"}, "Name of the chi2 matching function"};
    Configurable<float> cfgChi2FunctionMatchingPlaneZ{"cfgChi2FunctionMatchingPlaneZ", static_cast<float>(o2::mft::constants::mft::LayerZCoordinate()[9]), "Z position of the matching plane"};
  } configChi2MatchingOptions;

  // ML interface (single configurable model)
  struct : ConfigurableGroup {
    Configurable<std::string> cfgMlModelLabel{"cfgMlModelLabel", std::string{""}, "Text label identifying this ML model"};
    Configurable<std::string> cfgMlModelPathCcdb{"cfgMlModelPathCcdb", "Users/m/mcoquet/MLTest", "Path of model on CCDB"};
    Configurable<std::string> cfgMlModelName{"cfgMlModelName", "model.onnx", "ONNX file name (if not from CCDB full path)"};
    Configurable<std::vector<std::string>> cfgMlInputFeatures{"cfgMlInputFeatures", std::vector<std::string>{"chi2MCHMFT"}, "Names of ML model input features"};
    Configurable<float> cfgMlModelMatchingPlaneZ{"cfgMlModelMatchingPlaneZ", static_cast<float>(o2::mft::constants::mft::LayerZCoordinate()[9]), "Z position of the matching plane"};
  } configMlOptions;

  std::vector<double> binsPtMl;
  std::array<double, 1> cutValues{};
  std::vector<int> cutDirMl;
  bool hasActiveChi2Matching{false};
  std::string activeChi2FunctionName;
  double activeChi2MatchingPlaneZ{0.};

  bool hasActiveMlMatching{false};
  o2::analysis::MlResponseMFTMuonMatch<float> activeMlResponse;
  double activeMlMatchingPlaneZ{0.};

  int mRunNumber{0}; // needed to detect if the run changed and trigger update of magnetic field

  Service<o2::ccdb::BasicCCDBManager> ccdbManager{};
  o2::ccdb::CcdbApi fCCDBApi;

  // vector of all MFT-MCH(-MID) matching candidates associated to the same MCH(-MID) track,
  // to be sorted in descending order with respect to the matching score
  // the map key is the MCH(-MID) track global index
  using MatchingCandidates = std::map<int64_t, std::vector<MatchingCandidate>>;

  class TrackParExt : public o2::track::TrackParCovFwd
  {
   public:
    TrackParExt() = default;
    TrackParExt(const TrackParExt& t) = default;
    explicit TrackParExt(o2::track::TrackParCovFwd const& t, int nc = -1, bool r = false)
      : TrackParCovFwd(t), nClusters(nc), removable(r) {}
    ~TrackParExt() = default;

    TrackParExt& operator=(const TrackParCovFwd& tpf)
    {
      o2::track::TrackParCovFwd::operator=(tpf);
      return *this;
    }
    TrackParExt& operator=(const TrackParExt& tpe)
    {
      o2::track::TrackParCovFwd::operator=(tpe);
      nClusters = tpe.getNClusters();
      removable = tpe.isRemovable();
      return *this;
    }

    void setNClusters(int n) { nClusters = n; }
    [[nodiscard]] int getNClusters() const { return nClusters; }

    void setRemovable() { removable = true; }
    [[nodiscard]] bool isRemovable() const { return removable; }

    [[nodiscard]] o2::track::TrackParCovFwd asTrackParCovFwd() const
    {
      return {static_cast<const o2::track::TrackParCovFwd&>(*this)};
    }

   private:
    int nClusters{-1};
    bool removable{false};
  };

  std::unordered_map<int64_t, TrackParExt> mMchTrackPars;
  std::unordered_map<int64_t, TrackParExt> mMftTrackPars;

  std::unordered_map<int64_t, int32_t> mftTrackCovs;

  Produces<o2::aod::StoredFwdTracksReAlign> gmCandidateFwdTracks;
  Produces<o2::aod::StoredFwdTrksCovReAlign> gmCandidateFwdTracksCov;
  Produces<o2::aod::GmmCandFwdTrkExtras> gmCandidateFwdTrackExtras;
  Produces<o2::aod::AmbiguousFwdTrksReAlign> gmAmbiguousFwdTracksReAlign;

  int32_t mGmmCandFwdTrackRowIndex{0};
  std::unordered_map<int64_t, std::array<int32_t, 2>> mAmbBcSliceByFwdTrackId;
  bool mHasLastMchAmbiguousBcSlice{false};
  std::array<int32_t, 2> mLastMchAmbiguousBcSlice{};

  std::unordered_map<int64_t, MchTrackInfo> mMchTrackInfos;
  std::unordered_map<int64_t, std::vector<MatchingCandidate>> mStoredMatchingCandidates;
  std::unordered_map<int64_t, int32_t> mFwdTrackToGmmCandTrkIndex;

  mch::TrackFitter trackFitter; // Track fitter from MCH tracking library
  mch::geo::TransformationCreator transformation;
  std::map<int, math_utils::Transform3D> transformRef; // reference geometry w.r.t track data
  std::map<int, math_utils::Transform3D> transformNew; // new geometry
  double mImproveCutChi2{0.};                          // Chi2 cut for track improvement.
  TGeoManager* geoNew = nullptr;
  TGeoManager* geoRef = nullptr;
  globaltracking::MatchGlobalFwd mMatching;

  Preslice<aod::FwdTrkCls> perMuon = aod::fwdtrkcl::fwdtrackId;

  template <class T>
  o2::mch::TrackParam fwdToMch(const T& fwdtrack)
  {
    // Convert Forward Track parameters and covariances matrix to the
    // MCH track format.

    // Parameter conversion
    const double x2 = fwdtrack.getPhi();
    const double x3 = fwdtrack.getTanl();
    const double x4 = fwdtrack.getInvQPt();

    const auto sinX2 = std::sin(x2);
    const auto cosX2 = std::cos(x2);

    const double alpha1 = cosX2 / x3;
    const double alpha3 = sinX2 / x3;
    const double alpha4 = x4 / std::sqrt(x3 * x3 + sinX2 * sinX2);

    const auto kNorm = std::sqrt(x3 * x3 + sinX2 * sinX2);
    const auto kNorm3 = kNorm * kNorm * kNorm;

    // Covariances matrix conversion
    SMatrix55Std jacobian;
    SMatrix55Sym covariances;

    covariances(0, 0) = fwdtrack.getCovariances()(0, 0);
    covariances(0, 1) = fwdtrack.getCovariances()(0, 1);
    covariances(0, 2) = fwdtrack.getCovariances()(0, 2);
    covariances(0, 3) = fwdtrack.getCovariances()(0, 3);
    covariances(0, 4) = fwdtrack.getCovariances()(0, 4);

    covariances(1, 1) = fwdtrack.getCovariances()(1, 1);
    covariances(1, 2) = fwdtrack.getCovariances()(1, 2);
    covariances(1, 3) = fwdtrack.getCovariances()(1, 3);
    covariances(1, 4) = fwdtrack.getCovariances()(1, 4);

    covariances(2, 2) = fwdtrack.getCovariances()(2, 2);
    covariances(2, 3) = fwdtrack.getCovariances()(2, 3);
    covariances(2, 4) = fwdtrack.getCovariances()(2, 4);

    covariances(3, 3) = fwdtrack.getCovariances()(3, 3);
    covariances(3, 4) = fwdtrack.getCovariances()(3, 4);

    covariances(4, 4) = fwdtrack.getCovariances()(4, 4);

    jacobian(0, 0) = 1;

    jacobian(1, 2) = -sinX2 / x3;
    jacobian(1, 3) = -cosX2 / (x3 * x3);

    jacobian(2, 1) = 1;

    jacobian(3, 2) = cosX2 / x3;
    jacobian(3, 3) = -sinX2 / (x3 * x3);

    jacobian(4, 2) = -x4 * sinX2 * cosX2 / kNorm3;
    jacobian(4, 3) = -x3 * x4 / kNorm3;
    jacobian(4, 4) = 1 / kNorm;
    // jacobian*covariances*jacobian^T
    covariances = ROOT::Math::Similarity(jacobian, covariances);

    std::array<double, 15> cov = {covariances(0, 0), covariances(1, 0), covariances(1, 1), covariances(2, 0), covariances(2, 1), covariances(2, 2), covariances(3, 0), covariances(3, 1), covariances(3, 2), covariances(3, 3), covariances(4, 0), covariances(4, 1), covariances(4, 2), covariances(4, 3), covariances(4, 4)};
    std::array<double, 5> param = {fwdtrack.getX(), alpha1, fwdtrack.getY(), alpha3, alpha4};

    o2::mch::TrackParam convertedTrack(fwdtrack.getZ(), param.data(), cov.data());
    return {convertedTrack};
  }

  o2::track::TrackParCovFwd mchToFwd(const o2::mch::TrackParam& mchParam)
  {
    // Convert a MCH Track parameters and covariances matrix to the
    // Forward track format. Must be called after propagation though the absorber

    o2::track::TrackParCovFwd convertedTrack;

    // Parameter conversion
    const double alpha1 = mchParam.getNonBendingSlope();
    const double alpha3 = mchParam.getBendingSlope();
    const double alpha4 = mchParam.getInverseBendingMomentum();

    const double x2 = std::atan2(-alpha3, -alpha1);
    const double x3 = -1. / std::sqrt(alpha3 * alpha3 + alpha1 * alpha1);
    const double x4 = alpha4 * -x3 * std::sqrt(1 + alpha3 * alpha3);

    const auto kNorm = alpha1 * alpha1 + alpha3 * alpha3;
    const auto kNorm32 = kNorm * std::sqrt(kNorm);
    const auto slopeLen = std::sqrt(alpha3 * alpha3 + 1);

    // Covariances matrix conversion
    SMatrix55Std jacobian;
    SMatrix55Sym covariances;

    covariances(0, 0) = mchParam.getCovariances()(0, 0);
    covariances(0, 1) = mchParam.getCovariances()(0, 1);
    covariances(0, 2) = mchParam.getCovariances()(0, 2);
    covariances(0, 3) = mchParam.getCovariances()(0, 3);
    covariances(0, 4) = mchParam.getCovariances()(0, 4);

    covariances(1, 1) = mchParam.getCovariances()(1, 1);
    covariances(1, 2) = mchParam.getCovariances()(1, 2);
    covariances(1, 3) = mchParam.getCovariances()(1, 3);
    covariances(1, 4) = mchParam.getCovariances()(1, 4);

    covariances(2, 2) = mchParam.getCovariances()(2, 2);
    covariances(2, 3) = mchParam.getCovariances()(2, 3);
    covariances(2, 4) = mchParam.getCovariances()(2, 4);

    covariances(3, 3) = mchParam.getCovariances()(3, 3);
    covariances(3, 4) = mchParam.getCovariances()(3, 4);

    covariances(4, 4) = mchParam.getCovariances()(4, 4);

    jacobian(0, 0) = 1;

    jacobian(1, 2) = 1;

    jacobian(2, 1) = -alpha3 / kNorm;
    jacobian(2, 3) = alpha1 / kNorm;

    jacobian(3, 1) = alpha1 / kNorm32;
    jacobian(3, 3) = alpha3 / kNorm32;

    jacobian(4, 1) = -alpha1 * alpha4 * slopeLen / kNorm32;
    jacobian(4, 3) = alpha3 * alpha4 * (1 / (std::sqrt(kNorm) * slopeLen) - slopeLen / kNorm32);
    jacobian(4, 4) = slopeLen / std::sqrt(kNorm);

    // jacobian*covariances*jacobian^T
    covariances = ROOT::Math::Similarity(jacobian, covariances);

    // Set output
    convertedTrack.setX(mchParam.getNonBendingCoor());
    convertedTrack.setY(mchParam.getBendingCoor());
    convertedTrack.setZ(mchParam.getZ());
    convertedTrack.setPhi(x2);
    convertedTrack.setTanl(x3);
    convertedTrack.setInvQPt(x4);
    convertedTrack.setCharge(mchParam.getCharge());
    convertedTrack.setCovariances(covariances);

    return convertedTrack;
  }

  int getDetElemId(int iDetElemNumber)
  {
    // make sure detector number is valid
    if (iDetElemNumber < SNDetElemCh[0] ||
        iDetElemNumber >= SNDetElemCh[NMchChambers]) {
      LOGF(fatal, "Invalid detector element number: %d", iDetElemNumber);
    }
    /// get det element number from ID
    // get chamber and element number in chamber
    int iCh = 0;
    int iDet = 0;
    for (int i = 1; i <= NMchChambers; i++) {
      if (iDetElemNumber < SNDetElemCh[i]) {
        iCh = i;
        iDet = iDetElemNumber - SNDetElemCh[i - 1];
        break;
      }
    }

    // make sure detector index is valid
    if (iCh <= 0 || iCh > NMchChambers || iDet >= NDetElemCh[iCh - 1]) {
      LOGF(fatal, "Invalid detector element id: %d", MchDetElemNumberingBase * iCh + iDet);
    }

    // add number of detectors up to this chamber
    return MchDetElemNumberingBase * iCh + iDet;
  }

  bool removeTrack(mch::Track& track)
  {
    // Refit track with re-aligned clusters
    bool shouldRemoveTrack = false;
    try {
      trackFitter.fit(track, false);
    } catch (std::exception const& e) {
      shouldRemoveTrack = true;
      return shouldRemoveTrack;
    }

    auto itStartingParam = std::prev(track.rend());

    while (true) {

      try {
        trackFitter.fit(track, true, false, (itStartingParam == track.rbegin()) ? nullptr : &itStartingParam);
      } catch (std::exception const&) {
        shouldRemoveTrack = true;
        break;
      }

      double worstLocalChi2 = -1.0;

      track.tagRemovableClusters(0x1F, false);

      auto itWorstParam = track.end();

      for (auto itParam = track.begin(); itParam != track.end(); ++itParam) {
        if (itParam->getLocalChi2() > worstLocalChi2) {
          worstLocalChi2 = itParam->getLocalChi2();
          itWorstParam = itParam;
        }
      }

      if (worstLocalChi2 < mImproveCutChi2) {
        break;
      }

      if (!itWorstParam->isRemovable()) {
        shouldRemoveTrack = true;
        track.removable();
        break;
      }

      auto itNextParam = track.removeParamAtCluster(itWorstParam);
      auto itNextToNextParam = (itNextParam == track.end()) ? itNextParam : std::next(itNextParam);
      itStartingParam = track.rbegin();

      if (track.getNClusters() < MinRemovableTrackClusters) {
        shouldRemoveTrack = true;
        break;
      }
      while (itNextToNextParam != track.end()) {
        if (itNextToNextParam->getClusterPtr()->getChamberId() != itNextParam->getClusterPtr()->getChamberId()) {
          itStartingParam = std::make_reverse_iterator(++itNextParam);
          break;
        }
        ++itNextToNextParam;
      }
    }

    if (!shouldRemoveTrack) {
      for (auto& param : track) { // o2-linter: disable=const-ref-in-for-loop (object is modified in loop)
        param.setParameters(param.getSmoothParameters());
        param.setCovariances(param.getSmoothCovariances());
      }
    }

    return shouldRemoveTrack;
  }

  template <typename BC>
  void initCcdb(BC const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }

    mRunNumber = bc.runNumber();
    std::map<std::string, std::string> metadata;
    auto soreor = o2::ccdb::BasicCCDBManager::getRunDuration(fCCDBApi, mRunNumber);
    auto ts = soreor.first;
    auto grpmag = fCCDBApi.retrieveFromTFileAny<o2::parameters::GRPMagField>(configCcdb.cfgGrpMagPath, metadata, ts);
    o2::base::Propagator::initFieldFromGRP(grpmag);
    LOGF(info, "Set field for muons");
    VarManager::SetupMuonMagField();
    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      ccdbManager->get<TGeoManager>(configCcdb.cfgGeoPath);
    }
    mch::TrackExtrap::setField();
    mch::TrackExtrap::useExtrapV2();

    // Load geometry information from CCDB/local
    LOGF(info, "Loading reference aligned geometry from CCDB no later than %d", configMchRealign.cfgCcdbNoLaterThanRef.value);
    ccdbManager->setCreatedNotAfter(configMchRealign.cfgCcdbNoLaterThanRef.value); // this timestamp has to be consistent with what has been used in reco
    geoRef = ccdbManager->getForTimeStamp<TGeoManager>(configMchRealign.cfgGeoRefPath, bc.timestamp());
    ccdbManager->clearCache(configMchRealign.cfgGeoRefPath);
    if (geoRef != nullptr) {
      transformation = mch::geo::transformationFromTGeoManager(*geoRef);
    } else {
      LOGF(fatal, "Reference aligned geometry object is not available in CCDB at timestamp=%llu", bc.timestamp());
    }
    for (int i = 0; i < NMchDetElems; i++) {
      int iDEN = getDetElemId(i);
      transformRef[iDEN] = transformation(iDEN);
    }

    LOGF(info, "Loading new aligned geometry from CCDB no later than %d", configMchRealign.cfgCcdbNoLaterThanNew.value);
    ccdbManager->setCreatedNotAfter(configMchRealign.cfgCcdbNoLaterThanNew.value); // make sure this timestamp can be resolved regarding the reference one
    geoNew = ccdbManager->getForTimeStamp<TGeoManager>(configMchRealign.cfgGeoNewPath, bc.timestamp());
    ccdbManager->clearCache(configMchRealign.cfgGeoNewPath);
    if (geoNew != nullptr) {
      transformation = mch::geo::transformationFromTGeoManager(*geoNew);
    } else {
      LOGF(fatal, "New aligned geometry object is not available in CCDB at timestamp=%llu", bc.timestamp());
    }
    for (int i = 0; i < NMchDetElems; i++) {
      int iDEN = getDetElemId(i);
      transformNew[iDEN] = transformation(iDEN);
    }

    // Init magnetic field for MFT track extrapolation
    auto* fieldB = dynamic_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField());
    if (fieldB) {
      std::array<double, 3> centerMft{0, 0, -61.4}; // Field at center of MFT
      mBzAtMftCenter = fieldB->getBz(centerMft.data());
      // std::cout << "fieldB: " << (void*)fieldB << std::endl;
    }
  }

  void initMatchingFunctions()
  {
    using SVector2 = ROOT::Math::SVector<double, 2>;
    using SVector4 = ROOT::Math::SVector<double, 4>;
    using SVector5 = ROOT::Math::SVector<double, 5>;

    using SMatrix44 = ROOT::Math::SMatrix<double, 4>;
    using SMatrix45 = ROOT::Math::SMatrix<double, 4, 5>;
    using SMatrix22 = ROOT::Math::SMatrix<double, 2>;
    using SMatrix25 = ROOT::Math::SMatrix<double, 2, 5>;

    // Define built-in matching functions
    //________________________________________________________________________________
    mMatchingFunctionMap["matchALL"] = [](const o2::track::TrackParCovFwd& mchTrack, const o2::track::TrackParCovFwd& mftTrack) -> std::tuple<double, int> {
      // Match two tracks evaluating all parameters: X,Y, phi, tanl & q/pt

      SMatrix55Sym hK, vK;
      SVector5 mK(mftTrack.getX(), mftTrack.getY(), mftTrack.getPhi(),
                  mftTrack.getTanl(), mftTrack.getInvQPt()),
        rKKminus1;
      const auto& globalMuonTrackParameters = mchTrack.getParameters();
      const auto& globalMuonTrackCovariances = mchTrack.getCovariances();
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
    mMatchingFunctionMap["matchXYPhiTanl"] = [](const o2::track::TrackParCovFwd& mchTrack, const o2::track::TrackParCovFwd& mftTrack) -> std::tuple<double, int> {
      // Match two tracks evaluating positions & angles

      SMatrix45 hK;
      SMatrix44 vK;
      SVector4 mK(mftTrack.getX(), mftTrack.getY(), mftTrack.getPhi(),
                  mftTrack.getTanl()),
        rKKminus1;
      const auto& globalMuonTrackParameters = mchTrack.getParameters();
      const auto& globalMuonTrackCovariances = mchTrack.getCovariances();
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
    mMatchingFunctionMap["matchXY"] = [](const o2::track::TrackParCovFwd& mchTrack, const o2::track::TrackParCovFwd& mftTrack) -> std::tuple<double, int> {
      // Calculate Matching Chi2 - X and Y positions

      SMatrix25 hK;
      SMatrix22 vK;
      SVector2 mK(mftTrack.getX(), mftTrack.getY()), rKKminus1;
      const auto& globalMuonTrackParameters = mchTrack.getParameters();
      const auto& globalMuonTrackCovariances = mchTrack.getCovariances();
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
    ccdbManager->setURL(configCcdb.cfgCcdbUrl);
    ccdbManager->setCaching(true);
    ccdbManager->setLocalObjectValidityChecking();
    fCCDBApi.init(configCcdb.cfgCcdbUrl);
    mRunNumber = 0;

    // Configuration for track fitter
    const auto& trackerParam = mch::TrackerParam::Instance();
    trackFitter.setBendingVertexDispersion(trackerParam.bendingVertexDispersion);
    trackFitter.setChamberResolution(configMchRealign.cfgChamberResolutionX.value, configMchRealign.cfgChamberResolutionY.value);
    trackFitter.smoothTracks(true);
    trackFitter.useChamberResolution();
    mImproveCutChi2 = 2. * configMchRealign.cfgSigmaCutImprove.value * configMchRealign.cfgSigmaCutImprove.value;

    // Reset matching configuration, then populate only what we need.
    hasActiveChi2Matching = false;
    activeChi2FunctionName.clear();
    activeChi2MatchingPlaneZ = 0.;

    hasActiveMlMatching = false;
    activeMlMatchingPlaneZ = 0.;

    if (configMatching.cfgCustomMatchingStrategy.value == 0) {
      // Matching functions (custom chi2)
      initMatchingFunctions();
      auto label = configChi2MatchingOptions.cfgChi2FunctionLabel.value;
      auto funcName = configChi2MatchingOptions.cfgChi2FunctionName.value;
      auto matchingPlaneZ = configChi2MatchingOptions.cfgChi2FunctionMatchingPlaneZ.value;

      if (!label.empty() && !funcName.empty()) {
        hasActiveChi2Matching = true;
        activeChi2FunctionName = funcName;
        activeChi2MatchingPlaneZ = matchingPlaneZ;
      }
    } else {
      // Matching ML models (custom ML)
      // TODO : for now we use hard coded values since the current models use 1 pT bin
      binsPtMl = {-1e-6, 1000.0};
      cutValues = {0.0};
      cutDirMl = {cuts_ml::CutNot};
      LabeledArray<double> mycutsMl(cutValues.data(), 1, 1, std::vector<std::string>{"pT bin 0"}, std::vector<std::string>{"score"});

      auto label = configMlOptions.cfgMlModelLabel.value;
      auto modelPath = configMlOptions.cfgMlModelPathCcdb.value;
      auto inputFeatures = configMlOptions.cfgMlInputFeatures.value;
      auto modelName = configMlOptions.cfgMlModelName.value;
      auto matchingPlaneZ = configMlOptions.cfgMlModelMatchingPlaneZ.value;

      if (!label.empty() && !modelPath.empty() && !inputFeatures.empty() && !modelName.empty()) {
        activeMlResponse.configure(binsPtMl, mycutsMl, cutDirMl, 1);
        activeMlResponse.setModelPathsCCDB(std::vector<std::string>{modelName}, fCCDBApi, std::vector<std::string>{modelPath}, configCcdb.cfgCcdbNoLaterThan.value);
        activeMlResponse.cacheInputFeaturesIndices(inputFeatures);
        activeMlResponse.init();

        hasActiveMlMatching = true;
        activeMlMatchingPlaneZ = matchingPlaneZ;
      }
    }
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
    return pDCA <= nSigmaPDCA * sigmaPDCAWithRes;
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
    if (mchTrack.chi2() > chi2Cut) {
      return false;
    }

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
    return pDcaCut(mchTrack, collision, nSigmaPdcaCut);
  }

  void storeFwdTrackCovariance(const SMatrix55Sym& cov)
  {
    const float sigX = std::sqrt(cov(0, 0));
    const float sigY = std::sqrt(cov(1, 1));
    const float sigPhi = std::sqrt(cov(2, 2));
    const float sigTgl = std::sqrt(cov(3, 3));
    const float sig1Pt = std::sqrt(cov(4, 4));
    const auto rhoXY = static_cast<int8_t>(128.f * cov(0, 1) / (sigX * sigY));
    const auto rhoPhiX = static_cast<int8_t>(128.f * cov(0, 2) / (sigPhi * sigX));
    const auto rhoPhiY = static_cast<int8_t>(128.f * cov(1, 2) / (sigPhi * sigY));
    const auto rhoTglX = static_cast<int8_t>(128.f * cov(0, 3) / (sigTgl * sigX));
    const auto rhoTglY = static_cast<int8_t>(128.f * cov(1, 3) / (sigTgl * sigY));
    const auto rhoTglPhi = static_cast<int8_t>(128.f * cov(2, 3) / (sigTgl * sigPhi));
    const auto rho1PtX = static_cast<int8_t>(128.f * cov(0, 4) / (sig1Pt * sigX));
    const auto rho1PtY = static_cast<int8_t>(128.f * cov(1, 4) / (sig1Pt * sigY));
    const auto rho1PtPhi = static_cast<int8_t>(128.f * cov(2, 4) / (sig1Pt * sigPhi));
    const auto rho1PtTgl = static_cast<int8_t>(128.f * cov(3, 4) / (sig1Pt * sigTgl));
    gmCandidateFwdTracksCov(sigX, sigY, sigPhi, sigTgl, sig1Pt,
                            rhoXY, rhoPhiY, rhoPhiX, rhoTglX, rhoTglY, rhoTglPhi, rho1PtX, rho1PtY, rho1PtPhi, rho1PtTgl);
  }

  bool isMchTrackTagged(int64_t mchTrackIndex) const
  {
    const auto it = mMchTrackInfos.find(mchTrackIndex);
    return it != mMchTrackInfos.end() && it->second.isTagged;
  }

  template <class TMCH>
  void fillBaseGmmCandFwdTrack(TMCH const& track,
                               TrackParExt const& trackPar,
                               int32_t gmmMchTrackId,
                               float chi2MatchMCHMFT,
                               float matchScoreMCHMFT,
                               bool isTagged)
  {
    const auto collisionId = track.collisionId();
    bool hasBcSlice = false;
    std::array<int32_t, 2> bcSlice{};
    if (collisionId < 0) {
      const auto ambIt = mAmbBcSliceByFwdTrackId.find(track.globalIndex());
      if (ambIt != mAmbBcSliceByFwdTrackId.end()) {
        bcSlice = ambIt->second;
        hasBcSlice = true;
      }
    }

    gmCandidateFwdTracks(
      collisionId,
      track.trackType(),
      trackPar.getX(),
      trackPar.getY(),
      trackPar.getZ(),
      trackPar.getPhi(),
      trackPar.getTgl(),
      trackPar.getInvQPt(),
      trackPar.getNClusters(),
      track.pDca(),
      track.rAtAbsorberEnd(),
      trackPar.isRemovable(),
      trackPar.getTrackChi2(),
      track.chi2MatchMCHMID(),
      chi2MatchMCHMFT,
      matchScoreMCHMFT,
      track.matchMFTTrackId(),
      gmmMchTrackId,
      track.mchBitMap(),
      track.midBitMap(),
      track.midBoards(),
      track.trackTime(),
      track.trackTimeRes());

    storeFwdTrackCovariance(trackPar.getCovariances());
    gmCandidateFwdTrackExtras(isTagged, -1, -1);
    if (hasBcSlice) {
      gmAmbiguousFwdTracksReAlign(mGmmCandFwdTrackRowIndex, bcSlice.data());
    }
    mGmmCandFwdTrackRowIndex += 1;

    mHasLastMchAmbiguousBcSlice = hasBcSlice;
    if (hasBcSlice) {
      mLastMchAmbiguousBcSlice = bcSlice;
    }
  }

  template <class TMCH, class TMFT>
  void fillCandidateFwdTrack(TMCH const& mchTrack,
                             TrackParExt const& mchPar,
                             int32_t gmmMchTrackId,
                             TMFT const& mftTrack,
                             TrackParExt const& mftPar,
                             const MatchingCandidate& candidate)
  {
    using o2::aod::fwdtrack::ForwardTrackTypeEnum;
    using o2::aod::fwdtrackutils::propagationPoint;

    constexpr auto CandidateTrackType = static_cast<uint8_t>(ForwardTrackTypeEnum::GlobalForwardTrack);

    auto propmuonAtMft = fwdToMch(mchPar);
    o2::mch::TrackExtrap::extrapToVertex(propmuonAtMft,
                                         mftPar.getX(),
                                         mftPar.getY(),
                                         mftPar.getZ(),
                                         mftPar.getSigma2X(),
                                         mftPar.getSigma2Y());

    const auto globalMuonRefit = o2::aod::fwdtrackutils::refitGlobalMuonCov(mchToFwd(propmuonAtMft), mftPar);

    const auto nClusters = static_cast<int8_t>(std::min(127, mchPar.getNClusters() + mftPar.getNClusters()));

    const auto chi2 = static_cast<float>(mchTrack.chi2());
    const int32_t collisionId = mchTrack.collisionId();
    bool hasBcSlice = false;
    std::array<int32_t, 2> bcSlice{};
    if (collisionId < 0) {
      if (mHasLastMchAmbiguousBcSlice) {
        bcSlice = mLastMchAmbiguousBcSlice;
        hasBcSlice = true;
      } else {
        const auto ambIt = mAmbBcSliceByFwdTrackId.find(mchTrack.globalIndex());
        if (ambIt != mAmbBcSliceByFwdTrackId.end()) {
          bcSlice = ambIt->second;
          hasBcSlice = true;
        }
      }
    }

    bool isRemovable = mchPar.isRemovable();

    gmCandidateFwdTracks(
      collisionId,
      CandidateTrackType,
      globalMuonRefit.getX(),
      globalMuonRefit.getY(),
      globalMuonRefit.getZ(),
      globalMuonRefit.getPhi(),
      globalMuonRefit.getTgl(),
      globalMuonRefit.getInvQPt(),
      nClusters,
      mchTrack.pDca(),
      mchTrack.rAtAbsorberEnd(),
      isRemovable,
      chi2,
      mchTrack.chi2MatchMCHMID(),
      static_cast<float>(candidate.matchChi2),
      static_cast<float>(candidate.matchScore),
      static_cast<int>(mftTrack.globalIndex()),
      gmmMchTrackId,
      mchTrack.mchBitMap(),
      mchTrack.midBitMap(),
      mchTrack.midBoards(),
      mchTrack.trackTime(),
      mchTrack.trackTimeRes());

    storeFwdTrackCovariance(globalMuonRefit.getCovariances());
    gmCandidateFwdTrackExtras(isMchTrackTagged(mchTrack.globalIndex()), candidate.matchRanking, candidate.mixedGroupIndex);
    if (hasBcSlice) {
      gmAmbiguousFwdTracksReAlign(mGmmCandFwdTrackRowIndex, bcSlice.data());
    }
    mGmmCandFwdTrackRowIndex += 1;
  }

  o2::track::TrackParCovFwd propagateToZMch(const o2::track::TrackParCovFwd& muon, const double z)
  {
    auto mchTrack = fwdToMch(muon);

    float absFront = -90.f;
    float absBack = -505.f;

    if (muon.getZ() < absBack && z > absFront) {
      // extrapolation through the absorber in the upstream direction
      o2::mch::TrackExtrap::extrapToVertexWithoutBranson(mchTrack, z);
    } else {
      // all other cases
      o2::mch::TrackExtrap::extrapToZCov(mchTrack, z);
    }

    return mchToFwd(mchTrack);
  }

  o2::track::TrackParCovFwd propagateToZMft(const o2::track::TrackParCovFwd& mftTrack, const double z)
  {
    o2::track::TrackParCovFwd trackExtrap{mftTrack};
    trackExtrap.propagateToZ(z, mBzAtMftCenter);
    return trackExtrap;
  }

  template <class TMCH, class C>
  o2::track::TrackParCovFwd propagateToVertexMch(const TMCH& muon,
                                                 const C& collision)
  {
    auto mchTrack = fwdToMch(fwdtrackutils::getTrackParCovFwd(muon, muon));
    o2::mch::TrackExtrap::extrapToVertex(mchTrack,
                                         collision.posX(),
                                         collision.posY(),
                                         collision.posZ(),
                                         collision.covXX(),
                                         collision.covYY());
    return mchToFwd(mchTrack);
  }

  // tag muons based on the track quality and the track position at the front and back MFT planes
  template <class TMUON, class C>
  void getTaggedMuons(C const& collisions,
                      TMUON const& muonTracks,
                      std::vector<int64_t>& taggedMuons)
  {
    taggedMuons.clear();
    for (const auto& muonTrack : muonTracks) {

      // only consider MCH-MID matches
      if (static_cast<int>(muonTrack.trackType()) != MchMidTrackType) {
        continue;
      }

      // only select MCH-MID tracks associated to a collision
      if (!muonTrack.has_collision()) {
        continue;
      }

      const auto& collision = collisions.rawIteratorAt(muonTrack.collisionId());

      // select MCH tracks with strict quality cuts
      if (!isGoodMuon(muonTrack, collision,
                      configMuonTagging.cfgMuonTaggingTrackChi2MchUp,
                      configMuonTagging.cfgMuonTaggingPMchLow,
                      configMuonTagging.cfgMuonTaggingPtMchLow,
                      {configMuonTagging.cfgMuonTaggingEtaMchLow, configMuonTagging.cfgMuonTaggingEtaMchUp},
                      {configMuonTagging.cfgMuonTaggingRabsLow, configMuonTagging.cfgMuonTaggingRabsUp},
                      configMuonTagging.cfgMuonTaggingPdcaUp)) {
        continue;
      }

      // propagate MCH track to the vertex
      auto mchTrackAtVertex = propagateToVertexMch(muonTrack, collision);

      // propagate the track from the vertex to the first MFT plane
      const auto& extrapToMFTfirst = propagateToZMch(mchTrackAtVertex, o2::mft::constants::mft::LayerZCoordinate()[0]);
      double rFront = std::sqrt(extrapToMFTfirst.getX() * extrapToMFTfirst.getX() + extrapToMFTfirst.getY() * extrapToMFTfirst.getY());
      if (rFront < configMuonTagging.cfgMuonTaggingRadiusAtMftFrontLow.value || rFront > configMuonTagging.cfgMuonTaggingRadiusAtMftFrontUp.value) {
        continue;
      }

      // propagate the track from the vertex to the last MFT plane
      const auto& extrapToMFTlast = propagateToZMch(mchTrackAtVertex, o2::mft::constants::mft::LayerZCoordinate()[9]);
      double rBack = std::sqrt(extrapToMFTlast.getX() * extrapToMFTlast.getX() + extrapToMFTlast.getY() * extrapToMFTlast.getY());
      if (rBack < configMuonTagging.cfgMuonTaggingRadiusAtMftBackLow.value || rBack > configMuonTagging.cfgMuonTaggingRadiusAtMftBackUp.value) {
        continue;
      }

      int64_t muonTrackIndex = muonTrack.globalIndex();
      taggedMuons.emplace_back(muonTrackIndex);
    }
  }

  template <class EVT, class BC, class TMUON, class TMFT>
  bool isMftMchTimeCompatible(EVT const& collisions,
                              BC const& bcs,
                              TMUON const& mchTrack,
                              TMFT const& mftTrack)
  {
    if (!mchTrack.has_collision() || !mftTrack.has_collision()) {
      return false;
    }

    const auto& collMch = collisions.rawIteratorAt(mchTrack.collisionId());
    const auto& bcMch = bcs.rawIteratorAt(collMch.bcId());
    const auto& collMft = collisions.rawIteratorAt(mftTrack.collisionId());
    const auto& bcMft = bcs.rawIteratorAt(collMft.bcId());

    int64_t deltaBc = static_cast<int64_t>(bcMft.globalBC()) - static_cast<int64_t>(bcMch.globalBC());
    double deltaBcNS = o2::constants::lhc::LHCBunchSpacingNS * deltaBc;
    double deltaTrackTime = mftTrack.trackTime() - mchTrack.trackTime() + deltaBcNS;
    double trackTimeResTot = mftTrack.trackTimeRes() + mchTrack.trackTimeRes();

    return std::fabs(deltaTrackTime) <= trackTimeResTot;
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
  void prepareMatchingCandidates(EVT const& collisions,
                                 BC const& bcs,
                                 TMUON const& muonTracks,
                                 TMFT const& mftTracks,
                                 MyMFTCovariances const& mftCovs)
  {
    mMftTrackPars.clear();
    mMchTrackPars.clear();

    LOGF(info, "Filling matching candidate tables");

    for (const auto& muonTrack : muonTracks) {
      if (static_cast<int>(muonTrack.trackType()) <= GlobalTrackTypeMax) {
        continue;
      }
      auto mchTrackIndex = muonTrack.globalIndex();

      // initialize the MCH track parameters, which will be updated by the realignment if enabled
      mMchTrackPars.try_emplace(mchTrackIndex, TrackParExt(fwdtrackutils::getTrackParCovFwd(muonTrack, muonTrack), muonTrack.nClusters()));

      mMchTrackInfos.try_emplace(mchTrackIndex, MchTrackInfo{
                                                  .nMatchAttempts = 0,
                                                  .matchingCandidates = std::vector<MatchingCandidate>(),
                                                  .mixedMatchingCandidates = std::vector<std::vector<MatchingCandidate>>()});
    }

    for (const auto& mftTrack : mftTracks) {
      auto mftTrackIndex = mftTrack.globalIndex();

      // initialize the MFT track parameters, which will be updated by the alignment corrections if enabled
      if (mftTrackCovs.contains(mftTrackIndex) && !mMftTrackPars.contains(mftTrackIndex)) {
        auto const& mftTrackCov = mftCovs.rawIteratorAt(mftTrackCovs[mftTrackIndex]);
        mMftTrackPars.emplace(mftTrackIndex, TrackParExt(fwdtrackutils::getTrackParCovFwd(mftTrack, mftTrackCov), mftTrack.nClusters()));
      }
    }

    // fill matching candidates table
    if (!configMatching.cfgMatchAllTracks.value) {
      // collect global MFT-MCH or MFT-MCH-MID tracks and associate them to the corresponding MCH(-MID) track
      for (const auto& muonTrack : muonTracks) {
        // skip MCH or MCH-MID tracks
        if (static_cast<int>(muonTrack.trackType()) > GlobalTrackTypeMax) {
          continue;
        }

        auto const& mchTrack = muonTrack.template matchMCHTrack_as<TMUON>();
        int64_t mchTrackIndex = mchTrack.globalIndex();
        auto const& mftTrack = muonTrack.template matchMFTTrack_as<TMFT>();
        int64_t mftTrackIndex = mftTrack.globalIndex();

        if (!mftTrackCovs.contains(mftTrackIndex)) {
          continue;
        }

        auto& mchTrackInfo = mMchTrackInfos[mchTrackIndex];
        mchTrackInfo.matchingCandidates.emplace_back(MatchingCandidate{
          .muonTrackId = muonTrack.globalIndex(),
          .mftTrackId = mftTrackIndex,
          .matchScore = muonTrack.matchScoreMCHMFT(),
          .matchChi2 = muonTrack.chi2MatchMCHMFT()});
      }

      // set the number of match attempts for this track
      for (auto& mchTrackInfo : mMchTrackInfos) { // o2-linter: disable=const-ref-in-for-loop (object is modified in loop)
        const auto& mchTrack = muonTracks.rawIteratorAt(mchTrackInfo.first);
        mchTrackInfo.second.nMatchAttempts = getMftMchMatchAttempts(collisions, bcs, mchTrack, mftTracks);
      }
    } else {
      // build matching candidates from all time-compatible MFT-MCH pairs
      for (const auto& muonTrack : muonTracks) {
        if (static_cast<int>(muonTrack.trackType()) <= GlobalTrackTypeMax) {
          continue;
        }
        auto mchTrackIndex = muonTrack.globalIndex();
        for (const auto& mftTrack : mftTracks) {
          if (!isMftMchTimeCompatible(collisions, bcs, muonTrack, mftTrack)) {
            continue;
          }
          if (!mftTrackCovs.contains(mftTrack.globalIndex())) {
            continue;
          }

          auto& mchTrackInfo = mMchTrackInfos[mchTrackIndex];
          mchTrackInfo.matchingCandidates.emplace_back(MatchingCandidate{
            .mftTrackId = mftTrack.globalIndex()});
        }
      }

      // set the number of match attempts for this track
      for (auto& mchTrackInfo : mMchTrackInfos) { // o2-linter: disable=const-ref-in-for-loop (object is modified in loop)
        mchTrackInfo.second.nMatchAttempts = mchTrackInfo.second.matchingCandidates.size();
      }
    }

    std::vector<int64_t> taggedMuons;
    getTaggedMuons(collisions, muonTracks, taggedMuons);
    for (const auto& mchTrackIndex : taggedMuons) {
      auto it = mMchTrackInfos.find(mchTrackIndex);
      if (it != mMchTrackInfos.end()) {
        it->second.isTagged = true;
      }
    }
  }

  template <class EVT, class BC, class TMUON>
  void prepareEventMixingMatchingCandidates(EVT const& collisions,
                                            BC const& bcs,
                                            TMUON const& muonTracks)
  {
    LOGF(info, "Filling mixed matching candidate tables");

    const int mixingDepth = configEventMixing.cfgMixingDepth.value;
    const int64_t minDeltaBc = configEventMixing.cfgMinDeltaBc.value;
    const float maxDeltaPhi = configEventMixing.cfgMaxDeltaPhi.value;
    const float maxDeltaR = configEventMixing.cfgMaxDeltaR.value;
    const float maxDeltaAttemptsRel = configEventMixing.cfgMaxDeltaAttempts.value;
    const float maxDeltaZ = configEventMixing.cfgMaxDeltaZ.value;

    for (auto& [mchIndex1, mchTrackInfo1] : mMchTrackInfos) { // o2-linter: disable=const-ref-in-for-loop (object is modified in loop)
      const auto& mchTrack1 = muonTracks.rawIteratorAt(mchIndex1);

      if (!mchTrack1.has_collision()) {
        continue;
      }

      const auto& collision1 = collisions.rawIteratorAt(mchTrack1.collisionId());
      const auto& bc1 = bcs.rawIteratorAt(collision1.bcId());

      auto phi1 = std::atan2(mchTrack1.y(), mchTrack1.x());
      auto r1 = std::hypot(mchTrack1.x(), mchTrack1.y());

      auto nAttempts1 = mchTrackInfo1.nMatchAttempts;

      auto vz1 = collision1.posZ();

      for (const auto& [mchIndex2, mchTrackInfo2] : mMchTrackInfos) {
        if (mixingDepth >= 0 && static_cast<int>(mchTrackInfo1.mixedMatchingCandidates.size()) >= mixingDepth) {
          break;
        }

        const auto& mchTrack2 = muonTracks.rawIteratorAt(mchIndex2);

        if (!mchTrack2.has_collision()) {
          continue;
        }

        const auto& collision2 = collisions.rawIteratorAt(mchTrack2.collisionId());
        const auto& bc2 = bcs.rawIteratorAt(collision2.bcId());

        auto deltaBc = std::abs(static_cast<int64_t>(bc2.globalBC()) - static_cast<int64_t>(bc1.globalBC()));
        if (deltaBc < minDeltaBc) {
          continue;
        }

        auto phi2 = std::atan2(mchTrack2.y(), mchTrack2.x());
        auto deltaPhi = std::fabs(getDeltaPhi(phi1, phi2));
        if (deltaPhi > maxDeltaPhi) {
          continue;
        }

        auto r2 = std::hypot(mchTrack2.x(), mchTrack2.y());
        auto deltaR = r2 - r1;
        if (deltaR > maxDeltaR) {
          continue;
        }

        auto nAttempts2 = mchTrackInfo2.nMatchAttempts;

        float deltaAttempts = nAttempts2 - nAttempts1;
        float deltaAttemptsRel = (nAttempts1 > 0) ? deltaAttempts / nAttempts1 : 0;
        if (deltaAttemptsRel > maxDeltaAttemptsRel) {
          continue;
        }

        auto vz2 = collision2.posZ();
        auto deltaZ = std::fabs(vz2 - vz1);
        if (deltaZ > maxDeltaZ) {
          continue;
        }

        // add the candidates of MCH track #2 to the list of mixed candidates of track #1
        mchTrackInfo1.mixedMatchingCandidates.push_back(mchTrackInfo2.matchingCandidates);
        // update the muon track index of the mixed candidates to the index of track #1
        for (auto& candidate : mchTrackInfo1.mixedMatchingCandidates.back()) { // o2-linter: disable=const-ref-in-for-loop (object is modified in loop)
          candidate.muonTrackId = mchIndex1;
        }
      }
    }
  }

  template <typename TMFT, typename TMFTCOV>
  o2::track::TrackParCovFwd transformMft(TMFT& mftTrack, TMFTCOV const& mftTrackCov)
  {
    auto track = fwdToMch(fwdtrackutils::getTrackParCovFwd(mftTrack, mftTrackCov));

    double z = track.getZ();
    // double dZ = zMCH - z;
    double x = track.getNonBendingCoor();
    double y = track.getBendingCoor();
    double xSlope = track.getNonBendingSlope();
    double ySlope = track.getBendingSlope();

    double xSlopeCorrection = (y > 0) ? configMftAlignmentCorrections.cfgMFTAlignmentCorrXSlopeTop : configMftAlignmentCorrections.cfgMFTAlignmentCorrXSlopeBottom;
    double xCorrection = xSlopeCorrection * z +
                         ((y > 0) ? configMftAlignmentCorrections.cfgMFTAlignmentCorrXOffsetTop : configMftAlignmentCorrections.cfgMFTAlignmentCorrXOffsetBottom);
    double xNew = x + xCorrection;
    double xSlopeNew = xSlope + xSlopeCorrection;

    track.setNonBendingCoor(xNew);
    track.setNonBendingSlope(xSlopeNew);

    double ySlopeCorrection = (y > 0) ? configMftAlignmentCorrections.cfgMFTAlignmentCorrYSlopeTop : configMftAlignmentCorrections.cfgMFTAlignmentCorrYSlopeBottom;
    double yCorrection = ySlopeCorrection * z +
                         ((y > 0) ? configMftAlignmentCorrections.cfgMFTAlignmentCorrYOffsetTop : configMftAlignmentCorrections.cfgMFTAlignmentCorrYOffsetBottom);
    track.setBendingCoor(y + yCorrection);
    track.setBendingSlope(ySlope + ySlopeCorrection);

    return mchToFwd(track);
  }

  template <typename TMFTs, typename TMFTCOVs>
  void runMftRealignment(TMFTs const& mftTracks, TMFTCOVs const& mftCovs)
  {
    for (const auto& mftTrack : mftTracks) {
      auto mftTrackIndex = mftTrack.globalIndex();
      if (!mftTrackCovs.contains(mftTrackIndex)) {
        continue;
      }

      auto const& mftTrackCov = mftCovs.rawIteratorAt(mftTrackCovs[mftTrackIndex]);
      mMftTrackPars[mftTrackIndex] = transformMft(mftTrack, mftTrackCov);
    }
  }

  template <typename TMuons, typename TMuonCls>
  void runMuonRealignment(TMuons const& muons, TMuonCls const& clusters)
  {
    // Loop over forward tracks
    for (auto const& muon : muons) {
      int mchIndex = muon.globalIndex();
      // skip global forward matches
      if (muon.trackType() <= GlobalTrackTypeMax) {
        continue;
      }

      // continue;

      auto mchTrackParIt = mMchTrackPars.find(mchIndex);
      if (mchTrackParIt == mMchTrackPars.end()) {
        continue;
      }

      auto clustersSliced = clusters.sliceBy(perMuon, muon.globalIndex()); // Slice clusters by muon id
      mch::Track convertedTrack = mch::Track();                            // Temporary variable to store re-aligned clusters

      int clIndex = -1;
      // Get re-aligned clusters associated to current track
      for (auto const& cluster : clustersSliced) {
        clIndex += 1;

        auto* clusterMCH = new mch::Cluster();

        math_utils::Point3D<double> local;
        math_utils::Point3D<double> master;
        master.SetXYZ(cluster.x(), cluster.y(), cluster.z());

        // Transformation from reference geometry frame to new geometry frame
        transformRef[cluster.deId()].MasterToLocal(master, local);
        transformNew[cluster.deId()].LocalToMaster(local, master);

        clusterMCH->x = master.x();
        clusterMCH->y = master.y();
        clusterMCH->z = master.z();

        const uint32_t clUid = mch::Cluster::buildUniqueId(static_cast<int>(cluster.deId() / 100) - 1, cluster.deId(), clIndex);
        clusterMCH->uid = clUid;
        clusterMCH->ex = cluster.isGoodX() ? 0.2 : 10.0;
        clusterMCH->ey = cluster.isGoodY() ? 0.2 : 10.0;

        // Add transformed cluster into temporary variable
        convertedTrack.createParamAtCluster(*clusterMCH);
        // LOGF(debug, "Track %d, cluster DE%d:  x:%g  y:%g  z:%g", muon.globalIndex(), cluster.deId(), cluster.x(), cluster.y(), cluster.z());
        // LOGF(debug, "Track %d, re-aligned cluster DE%d:  x:%g  y:%g  z:%g", muonRealignId, cluster.deId(), clusterMCH->getX(), clusterMCH->getY(), clusterMCH->getZ());
      }

      // Refit the re-aligned track
      int removable = 0;
      if (convertedTrack.getNClusters() != 0) {
        removable = removeTrack(convertedTrack);
      } else {
        LOGF(fatal, "Muon track %d has no associated clusters.", muon.globalIndex());
      }

      // Get the re-aligned track parameter: track param at the first cluster
      mch::TrackParam trackParam = mch::TrackParam(convertedTrack.first());

      // Convert MCH track to FWD track and store new parameters after realignment
      mchTrackParIt->second = mchToFwd(mch::TrackParam(convertedTrack.first()));
      mchTrackParIt->second.setTrackChi2(trackParam.getTrackChi2() / convertedTrack.getNDF());
      mchTrackParIt->second.setNClusters(convertedTrack.getNClusters());
      if (removable) {
        mchTrackParIt->second.setRemovable();
      }
    }
  }

  void runChi2Matching(const std::string& funcName,
                       float matchingPlaneZ,
                       MatchingCandidates& newMatchingCandidates,
                       bool useMixedMatchingCandidates)
  {
    newMatchingCandidates.clear();

    std::string funcNameEffective = funcName;
    float matchingPlaneZEffective = matchingPlaneZ;
    if (funcName == "prod") {
      funcNameEffective = "matchALL";
      matchingPlaneZEffective = MatchingPlaneDefaultZ;
    }

    if (!mMatchingFunctionMap.contains(funcNameEffective)) {
      return;
    }
    auto matchingFunc = mMatchingFunctionMap.at(funcNameEffective);

    for (const auto& [mchIndex, mchTrackInfo] : mMchTrackInfos) {
      // get the tracks parameters, which have been updated by the realignment if enabled
      const auto mchTrackParIt = mMchTrackPars.find(mchIndex);
      if (mchTrackParIt == mMchTrackPars.end()) {
        continue;
      }

      auto processGroup = [&, mchIndex](const std::vector<MatchingCandidate>& candidatesGroup, int32_t mixedGroupIndex) {
        std::vector<MatchingCandidate> groupResults;
        groupResults.reserve(candidatesGroup.size());

        for (const auto& candidate : candidatesGroup) {
          auto mftTrackParIt = mMftTrackPars.find(candidate.mftTrackId);
          if (mftTrackParIt == mMftTrackPars.end()) {
            continue;
          }

          auto mftTrackProp = mftTrackParIt->second.asTrackParCovFwd();
          auto mchTrackProp = mchTrackParIt->second.asTrackParCovFwd();

          if (matchingPlaneZEffective < 0.) {
            mftTrackProp = propagateToZMft(mftTrackProp, matchingPlaneZ);
            mchTrackProp = propagateToZMch(mchTrackProp, matchingPlaneZ);
          }

          auto matchResult = matchingFunc(mchTrackProp, mftTrackProp);
          float matchChi2 = std::get<0>(matchResult);

          groupResults.emplace_back(MatchingCandidate{
            .muonTrackId = candidate.muonTrackId,
            .mftTrackId = candidate.mftTrackId,
            .matchScore = -1,
            .matchChi2 = matchChi2});
        }

        auto compareMatchingChi2 = [](const MatchingCandidate& track1, const MatchingCandidate& track2) -> bool {
          return (track1.matchChi2 < track2.matchChi2);
        };
        std::sort(groupResults.begin(), groupResults.end(), compareMatchingChi2);

        int ranking = 1;
        for (auto& matchedCandidate : groupResults) { // o2-linter: disable=const-ref-in-for-loop (object is modified in loop)
          matchedCandidate.matchRanking = ranking;
          ranking += 1;
        }

        const int maxCandidates = configMatching.cfgMaxCandidatesPerMchTrack.value;

        auto& storedCandidates = newMatchingCandidates[mchIndex];
        size_t nStoredThisGroup = 0;
        for (auto& result : groupResults) { // o2-linter: disable=const-ref-in-for-loop (object is modified in loop)
          if (maxCandidates >= 0 && nStoredThisGroup >= static_cast<size_t>(maxCandidates)) {
            break;
          }
          result.mixedGroupIndex = useMixedMatchingCandidates ? mixedGroupIndex : -1;
          storedCandidates.push_back(result);
          ++nStoredThisGroup;
        }
      };

      if (useMixedMatchingCandidates) {
        int32_t groupIdx = 0;
        for (const auto& candidatesGroup : mchTrackInfo.mixedMatchingCandidates) {
          processGroup(candidatesGroup, groupIdx);
          groupIdx += 1;
        }
      } else {
        processGroup(mchTrackInfo.matchingCandidates, -1);
      }
    }
  }

  template <class C, class TMUON, class TMFT>
  void runMlMatching(C const& collisions,
                     TMUON const& muonTracks,
                     TMFT const& mftTracks,
                     o2::analysis::MlResponseMFTMuonMatch<float>& mlResponse,
                     float matchingPlaneZ,
                     MatchingCandidates& newMatchingCandidates,
                     bool useMixedMatchingCandidates)
  {
    newMatchingCandidates.clear();
    for (const auto& [mchIndex, mchTrackInfo] : mMchTrackInfos) {
      auto const& mchTrack = muonTracks.rawIteratorAt(mchIndex);
      if (!mchTrack.has_collision()) {
        continue;
      }

      auto collision = collisions.rawIteratorAt(mchTrack.collisionId());

      // get the tracks parameters, which have been updated by the realignment if enabled
      auto mchTrackParIt = mMchTrackPars.find(mchIndex);
      if (mchTrackParIt == mMchTrackPars.end()) {
        continue;
      }

      auto processGroup = [&, mchIndex](const std::vector<MatchingCandidate>& candidatesGroup, int32_t mixedGroupIndex) {
        std::vector<MatchingCandidate> groupResults;
        groupResults.reserve(candidatesGroup.size());

        for (const auto& candidate : candidatesGroup) {
          auto const& muonTrack = (candidate.muonTrackId >= 0) ? muonTracks.rawIteratorAt(candidate.muonTrackId) : mchTrack;
          auto const& mftTrack = mftTracks.rawIteratorAt(candidate.mftTrackId);
          auto mftTrackParIt = mMftTrackPars.find(candidate.mftTrackId);
          if (mftTrackParIt == mMftTrackPars.end()) {
            continue;
          }

          auto mftTrackProp = mftTrackParIt->second.asTrackParCovFwd();
          auto mchTrackProp = mchTrackParIt->second.asTrackParCovFwd();

          if (matchingPlaneZ < 0.) {
            mftTrackProp = propagateToZMft(mftTrackProp, matchingPlaneZ);
            mchTrackProp = propagateToZMch(mchTrackProp, matchingPlaneZ);
          }

          std::vector<float> output;
          std::vector<float> inputML = mlResponse.getInputFeatures(muonTrack, mftTrack, mchTrack, mftTrackProp, mchTrackProp, collision);
          mlResponse.isSelectedMl(inputML, 0, output);
          float matchScore = output[0];

          groupResults.emplace_back(MatchingCandidate{
            .muonTrackId = candidate.muonTrackId,
            .mftTrackId = candidate.mftTrackId,
            .matchScore = matchScore,
            .matchChi2 = -1});
        }

        auto compareMatchingScore = [](const MatchingCandidate& track1, const MatchingCandidate& track2) -> bool {
          return (track1.matchScore > track2.matchScore);
        };
        std::sort(groupResults.begin(), groupResults.end(), compareMatchingScore);

        int ranking = 1;
        for (auto& matchedCandidate : groupResults) { // o2-linter: disable=const-ref-in-for-loop (object is modified in loop)
          matchedCandidate.matchRanking = ranking;
          ranking += 1;
        }

        const int maxCandidates = configMatching.cfgMaxCandidatesPerMchTrack.value;

        auto& storedCandidates = newMatchingCandidates[mchIndex];
        size_t nStoredThisGroup = 0;
        for (auto& result : groupResults) { // o2-linter: disable=const-ref-in-for-loop (object is modified in loop)
          if (maxCandidates >= 0 && nStoredThisGroup >= static_cast<size_t>(maxCandidates)) {
            break;
          }
          result.mixedGroupIndex = useMixedMatchingCandidates ? mixedGroupIndex : -1;
          storedCandidates.push_back(result);
          ++nStoredThisGroup;
        }
      };

      if (useMixedMatchingCandidates) {
        int32_t groupIdx = 0;
        for (const auto& candidatesGroup : mchTrackInfo.mixedMatchingCandidates) {
          processGroup(candidatesGroup, groupIdx);
          groupIdx += 1;
        }
      } else {
        processGroup(mchTrackInfo.matchingCandidates, -1);
      }
    }
  }

  template <class C, class TMUON, class TMFT, class CMFT>
  void processMatchingCandidates(C const& collisions,
                                 TMUON const& muonTracks,
                                 TMFT const& mftTracks,
                                 CMFT const& mftCovs,
                                 aod::FwdTrkCls const& clusters,
                                 bool useMixedMatchingCandidates)
  {
    if (configMchRealign.cfgEnableMCHRealign.value) {
      runMuonRealignment(muonTracks, clusters);
    }

    if (configMftAlignmentCorrections.cfgEnableMftAlignmentCorrections) {
      runMftRealignment(mftTracks, mftCovs);
    }

    MatchingCandidates newMatchingCandidates;
    if (configMatching.cfgCustomMatchingStrategy.value == 0) {
      if (hasActiveChi2Matching) {
        runChi2Matching(activeChi2FunctionName, activeChi2MatchingPlaneZ, newMatchingCandidates, useMixedMatchingCandidates);
      }
    } else {
      if (hasActiveMlMatching) {
        runMlMatching(collisions, muonTracks, mftTracks, activeMlResponse, activeMlMatchingPlaneZ, newMatchingCandidates, useMixedMatchingCandidates);
      }
    }

    for (const auto& [mchIndex, candidates] : newMatchingCandidates) {
      if (candidates.empty()) {
        continue;
      }

      mStoredMatchingCandidates[mchIndex] = candidates;
    }
  }

  int32_t countStoredCandidatesForMchTrack(int64_t mchTrackIndex) const
  {
    const auto candidateIterator = mStoredMatchingCandidates.find(mchTrackIndex);
    if (candidateIterator == mStoredMatchingCandidates.end()) {
      return 0;
    }
    return static_cast<int32_t>(candidateIterator->second.size());
  }

  template <class TMUON, class TMFT>
  void fillGmmCandidateFwdTracks(TMUON const& muonTracks,
                                 TMFT const& mftTracks,
                                 aod::AmbiguousFwdTracks const& ambFwdTracks)
  {
    mFwdTrackToGmmCandTrkIndex.clear();
    mGmmCandFwdTrackRowIndex = 0;
    mHasLastMchAmbiguousBcSlice = false;
    mAmbBcSliceByFwdTrackId.clear();
    for (const auto& ambFwdTrack : ambFwdTracks) {
      const auto bcIds = ambFwdTrack.bcIds();
      mAmbBcSliceByFwdTrackId[ambFwdTrack.fwdtrackId()] = {bcIds[0], bcIds[1]};
    }

    // First pass: assign GMMCANDTRK row indices for MCH/MCH-MID base entries so that
    // MCHTrackId can be remapped consistently even when global muons appear first in FwdTracks.
    int32_t nextGmmCandTrkIndex = 0;
    for (const auto& track : muonTracks) {
      const int trackType = static_cast<int>(track.trackType());
      if (trackType > GlobalTrackTypeMax) {
        mFwdTrackToGmmCandTrkIndex[track.globalIndex()] = nextGmmCandTrkIndex;
        nextGmmCandTrkIndex += 1 + countStoredCandidatesForMchTrack(track.globalIndex());
      } else if (configMatching.cfgIncludeGlobalMuonsInFwdTracks.value) {
        nextGmmCandTrkIndex += 1;
      }
    }

    // Second pass: fill GMMCANDTRK/GMMCANDTRKCOV in FwdTracks order.
    for (const auto& track : muonTracks) {
      const int trackType = static_cast<int>(track.trackType());

      if (trackType > GlobalTrackTypeMax) {
        mHasLastMchAmbiguousBcSlice = false;
        const int64_t mchTrackIndex = track.globalIndex();
        const int32_t gmmMchTrackId = mFwdTrackToGmmCandTrkIndex.at(mchTrackIndex);

        const auto candidateIterator = mStoredMatchingCandidates.find(mchTrackIndex);
        auto mchTrackParIt = mMchTrackPars.find(mchTrackIndex);
        fillBaseGmmCandFwdTrack(track, mchTrackParIt->second, gmmMchTrackId, -1.f, -1.f,
                                isMchTrackTagged(mchTrackIndex));

        if (candidateIterator != mStoredMatchingCandidates.end()) {
          for (const auto& candidate : candidateIterator->second) {
            auto mftTrackParIt = mMftTrackPars.find(candidate.mftTrackId);
            if (mftTrackParIt != mMftTrackPars.end()) {
              const auto& mftTrack = mftTracks.rawIteratorAt(candidate.mftTrackId);
              fillCandidateFwdTrack(track, mchTrackParIt->second, gmmMchTrackId, mftTrack, mftTrackParIt->second, candidate);
            }
          }
        }
      }

      if (configMatching.cfgIncludeGlobalMuonsInFwdTracks.value && trackType <= GlobalTrackTypeMax) {
        int32_t gmmMchTrackId = -1;
        const auto mchIterator = mFwdTrackToGmmCandTrkIndex.find(track.matchMCHTrackId());
        if (mchIterator != mFwdTrackToGmmCandTrkIndex.end()) {
          gmmMchTrackId = mchIterator->second;
        }
        TrackParExt parExt(fwdtrackutils::getTrackParCovFwd(track, track));
        fillBaseGmmCandFwdTrack(track,
                                parExt,
                                gmmMchTrackId,
                                track.chi2MatchMCHMFT(),
                                track.matchScoreMCHMFT(),
                                isMchTrackTagged(track.matchMCHTrackId()));
      }
    }
  }

  void runCandidateProcessing(MyEvents const& collisions,
                              aod::BCsWithTimestamps const& bcs,
                              MyMuons const& muonTracks,
                              MyMFTs const& mftTracks,
                              MyMFTCovariances const& mftCovs,
                              aod::FwdTrkCls const& clusters,
                              aod::AmbiguousFwdTracks const& ambFwdTracks,
                              bool useMixedMatchingCandidates)
  {
    auto bc = bcs.begin();
    initCcdb(bc);

    LOGF(info, "Filling MFT cov");
    mftTrackCovs.clear();
    for (const auto& mftTrackCov : mftCovs) {
      mftTrackCovs[mftTrackCov.matchMFTTrackId()] = mftTrackCov.globalIndex();
    }

    mStoredMatchingCandidates.clear();
    mFwdTrackToGmmCandTrkIndex.clear();
    mMchTrackInfos.clear();

    LOGF(info, "Preparing candidates");
    prepareMatchingCandidates(collisions, bcs, muonTracks, mftTracks, mftCovs);
    if (useMixedMatchingCandidates) {
      prepareEventMixingMatchingCandidates(collisions, bcs, muonTracks);
    }

    LOGF(info, "Processing candidates");
    processMatchingCandidates(collisions, muonTracks, mftTracks, mftCovs, clusters, useMixedMatchingCandidates);

    LOGF(info, "Filling tables");
    fillGmmCandidateFwdTracks(muonTracks, mftTracks, ambFwdTracks);
  }

  void processData(MyEvents const& collisions,
                   aod::BCsWithTimestamps const& bcs,
                   MyMuons const& muonTracks,
                   MyMFTs const& mftTracks,
                   MyMFTCovariances const& mftCovs,
                   aod::FwdTrkCls const& clusters,
                   aod::AmbiguousFwdTracks const& ambFwdTracks)
  {
    runCandidateProcessing(collisions, bcs, muonTracks, mftTracks, mftCovs, clusters, ambFwdTracks, false);
  }

  void processMixedData(MyEvents const& collisions,
                        aod::BCsWithTimestamps const& bcs,
                        MyMuons const& muonTracks,
                        MyMFTs const& mftTracks,
                        MyMFTCovariances const& mftCovs,
                        aod::FwdTrkCls const& clusters,
                        aod::AmbiguousFwdTracks const& ambFwdTracks)
  {
    runCandidateProcessing(collisions, bcs, muonTracks, mftTracks, mftCovs, clusters, ambFwdTracks, true);
  }

  PROCESS_SWITCH(GlobalMuonMatching, processData, "processData", true);
  PROCESS_SWITCH(GlobalMuonMatching, processMixedData, "process event-mixed matching candidates", false);
};

// Extends the fwdtracksrealign table with expression columns
struct GlobalMuonMatchingSpawner {
  Spawns<aod::FwdTrksCovReAlign> realignFwdTrksCov;
  Spawns<aod::FwdTracksReAlign> realignFwdTrks;
  void init(InitContext const&) {}
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<GlobalMuonMatching>(cfgc),
    adaptAnalysisTask<GlobalMuonMatchingSpawner>(cfgc)};
};
