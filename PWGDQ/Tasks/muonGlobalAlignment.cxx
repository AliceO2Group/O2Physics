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
/// \file muonGlobalAlignment.cxx
/// \brief Analysis of global alignment between MFT and MCH-MID
/// \author Andrea Ferrero <andrea.ferrero@cern.ch>, CEA-Saclay

#include "PWGDQ/Core/VarManager.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/RCTSelectionFlags.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/fwdtrackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CCDB/CcdbApi.h>
#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <CommonUtils/ConfigurableParam.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DetectorsBase/GeometryManager.h>
#include <DetectorsBase/Propagator.h>
#include <Field/MagneticField.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>
#include <GPU/GPUROOTCartesianFwd.h>
#include <MCHBase/TrackerParam.h>
#include <MCHGeometryTransformer/Transformations.h>
#include <MCHTracking/Track.h>
#include <MCHTracking/TrackExtrap.h>
#include <MCHTracking/TrackFitter.h>
#include <MCHTracking/TrackParam.h>
#include <MFTTracking/Constants.h>
#include <MathUtils/Cartesian.h>
#include <ReconstructionDataFormats/GlobalFwdTrack.h>
#include <ReconstructionDataFormats/TrackFwd.h>

#include <Math/MatrixFunctions.h>
#include <Math/MatrixRepresentationsStatic.h>
#include <Math/SMatrix.h>
#include <Math/SVector.h>
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>
#include <TGeoGlobalMagField.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <THnSparse.h>
#include <TMath.h>

#include <rapidjson/document.h>
#include <rapidjson/error/error.h>

#include <RtypesCore.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <exception>
#include <format>
#include <iterator>
#include <map>
#include <memory>
#include <random>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::mch;
using namespace o2::framework;
using namespace o2::aod;

using namespace o2::aod::rctsel;

// using MyReducedMuons = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra, aod::ReducedMuonsCov>;
// using MyReducedEvents = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended>;
// using MyReducedEventsVtxCov = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsVtxCov>;

using MyCollisions = aod::Collisions;
using MyBCs = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels>;
using MyEvents = soa::Join<aod::Collisions, aod::EvSels>;
using MyMuonsWithCov = soa::Join<aod::FwdTracks, aod::FwdTracksCov>;
// using MyMuonsWithCov = aod::FwdTracks;
using MyMFTs = aod::MFTTracks;
using MyMFTCovariances = aod::MFTTracksCov;

using MyCollision = MyCollisions::iterator;
using MyBC = MyBCs::iterator;
using MyMUON = MyMuonsWithCov::iterator;
using MyMFT = MyMFTs::iterator;
using MyMFTCovariance = MyMFTCovariances::iterator;

using SMatrix55 = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
using SMatrix5 = ROOT::Math::SVector<Double_t, 5>;

using o2::dataformats::GlobalFwdTrack;
using o2::track::TrackParCovFwd;

namespace o2::aod
{
DECLARE_SOA_TABLE(CompactMFTTracks, "AOD", "COMPACTMFT", //! standalone table for studying alignment
                  collision::PosX, collision::PosY, collision::PosZ,
                  fwdtrack::Signed1Pt, fwdtrack::Tgl, fwdtrack::Phi,
                  fwdtrack::FwdDcaX, fwdtrack::FwdDcaY, o2::aod::fwdtrack::CXX, o2::aod::fwdtrack::CYY, o2::aod::fwdtrack::CXY,
                  fwdtrack::NClusters, fwdtrack::Chi2,
                  fwdtrack::X, fwdtrack::Y, fwdtrack::Z);
using CompactMFTTrack = CompactMFTTracks;
} // namespace o2::aod

struct muonGlobalAlignment { // o2-linter: disable=name/workflow-file,name/struct (exception)

  static constexpr int GlobalTrackTypeMax = 2;
  static constexpr int NMchChambers = 10;
  static constexpr int NMchDetElems = 156;
  static constexpr int ThetaAbsBoundaryDeg = 3;
  static constexpr double SlopeResolutionZ = 535.;
  static constexpr double AbsorberBackZ = -505.f;
  static constexpr double BransonPlaneZ = -466.f;

  Produces<aod::CompactMFTTracks> mftTable;
  Configurable<bool> cfgProduceMFTTable{"cfgProduceMFTTable", false, "flag to produce MFTsa table"};

  ////   Variables for selecting MCH and MFT tracks
  Configurable<float> cfgTrackChi2MchUp{"cfgTrackChi2MchUp", 5.f, ""};
  Configurable<float> cfgPtMchLow{"cfgPtMchLow", 0.7f, ""};
  Configurable<float> cfgEtaMchLow{"cfgEtaMchLow", -4.0f, ""};
  Configurable<float> cfgEtaMchUp{"cfgEtaMchUp", -2.5f, ""};
  Configurable<float> cfgEtaMftLow{"cfgEtaMftLow", -3.6f, ""};
  Configurable<float> cfgEtaMftUp{"cfgEtaMftUp", -2.5f, ""};
  Configurable<float> cfgRabsLow{"cfgRabsLow", 17.6f, ""};
  Configurable<float> cfgRabsUp{"cfgRabsUp", 89.5f, ""};
  Configurable<float> fSigmaPdcaUp{"fSigmaPdcaUp", 6.f, ""};

  Configurable<int> cfgTrackNClustMftLow{"cfgTrackNClustMftLow", 7, ""};
  Configurable<float> cfgTrackChi2MftUp{"cfgTrackChi2MftUp", 999.f, ""};

  Configurable<float> cfgMftDcaMatchChi2Up{"cfgMftDcaMatchChi2Up", 50.f, ""};
  Configurable<float> cfgMftMchResidualsMatchChi2Up{"cfgMftMchResidualsMatchChi2Up", 50.f, ""};
  Configurable<float> cfgDimuonMatchChi2Up{"cfgDimuonMatchChi2Up", 50.f, ""};

  Configurable<float> cfgMftMchResidualsPLow{"cfgMftMchResidualsPLow", 30.f, ""};
  Configurable<float> cfgMftMchResidualsPtLow{"cfgMftMchResidualsPtLow", 4.f, ""};

  Configurable<uint32_t> cfgMftTracksMultiplicityMax{"cfgMftTracksMultiplicityMax", 0, "Maximum number of MFT tracks to be processed per event (zero means no limit)"};

  // Magnetic field position bias
  Configurable<float> cfgFieldOriginBiasZ{"cfgFieldOriginBiasZ", 0.0f, "Bias applied to the magnetic field z position"};
  Configurable<float> cfgDipoleZshift{"cfgDipoleZshift", 0.0f, "Correction to the dipole z position"};
  Configurable<float> cfgVertexZshift{"cfgVertexZshift", 0.0f, "Correction to the vertex z position"};

  ////   Variables for MFT alignment corrections
  struct : ConfigurableGroup {
    Configurable<bool> cfgEnableMFTAlignmentCorrections{"cfgEnableMFTAlignmentCorrections", false, ""};
    // slope corrections
    Configurable<float> cfgMFTAlignmentCorrXSlopeTop{"cfgMFTAlignmentCorrXSlopeTop", (-0.0006696 - 0.0005621) / 2.f, "MFT X slope correction - top half"};
    Configurable<float> cfgMFTAlignmentCorrXSlopeBottom{"cfgMFTAlignmentCorrXSlopeBottom", (0.00105 + 0.001007) / 2.f, "MFT X slope correction - bottom half"};
    Configurable<float> cfgMFTAlignmentCorrYSlopeTop{"cfgMFTAlignmentCorrYSlopeTop", (-0.002299 - 0.002442) / 2.f, "MFT Y slope correction - top half"};
    Configurable<float> cfgMFTAlignmentCorrYSlopeBottom{"cfgMFTAlignmentCorrYSlopeBottom", (-0.0005339 - 0.0006921) / 2.f, "MFT Y slope correction - bottom half"};
    // offset corrections
    Configurable<float> cfgMFTAlignmentCorrXOffsetTop{"cfgMFTAlignmentCorrXOffsetTop", 0.f, "MFT X offset correction - top half"};
    Configurable<float> cfgMFTAlignmentCorrXOffsetBottom{"cfgMFTAlignmentCorrXOffsetBottom", 0.f, "MFT X offset correction - bottom half"};
    Configurable<float> cfgMFTAlignmentCorrYOffsetTop{"cfgMFTAlignmentCorrYOffsetTop", 0.f, "MFT Y offset correction - top half"};
    Configurable<float> cfgMFTAlignmentCorrYOffsetBottom{"cfgMFTAlignmentCorrYOffsetBottom", 0.f, "MFT Y offset correction - bottom half"};
  } configMFTAlignmentCorrections;

  ////   Variables for re-alignment setup
  struct : ConfigurableGroup {
    Configurable<bool> cfgEnableMCHRefit{"cfgEnableMCHRefit", false, "Enable re-fitting of MCH tracks"};
    Configurable<bool> cfgEnableMCHRealign{"cfgEnableMCHRealign", false, "Enable re-alignment of MCH clusters and tracks"};
    Configurable<double> cfgChamberResolutionX{"cfgChamberResolutionX", 0.4, "Chamber resolution along X configuration for refit"}; // 0.4cm pp, 0.2cm PbPb
    Configurable<double> cfgChamberResolutionY{"cfgChamberResolutionY", 0.4, "Chamber resolution along Y configuration for refit"}; // 0.4cm pp, 0.2cm PbPb
    Configurable<double> cfgSigmaCutImprove{"cfgSigmaCutImprove", 6., "Sigma cut for track improvement"};
    Configurable<std::string> cfgMCHRealignCorrections{"cfgMCHRealignCorrections", "", "MCH DE positions/angles corrections in JSON format"};
  } configRealign;

  ////   Variables for ccdb
  struct : ConfigurableGroup {
    Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
    Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
    Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
    Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
    // Configurable<std::string> geoPathRealign{"geoPathRealign", "Users/j/jcastill/GeometryAlignedFix10Fix15ShiftCh1BNew2", "Path of the geometry file"};
    Configurable<std::string> geoPathRealign{"geoPathRealign", "Users/j/jcastill/GeometryAlignedLoczzm4pLHC24anap1sR5a", "Path of the geometry file"};
    Configurable<int64_t> cfgCcdbNoLaterThanRef{"cfgCcdbNoLaterThanRef", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object of reference basis"};
    Configurable<int64_t> cfgCcdbNoLaterThanNew{"cfgCcdbNoLaterThanNew", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object of new basis"};
  } configCCDB;

  Configurable<bool> cfgRequireGoodRCT{"cfgRequireGoodRCT", true, "Require good detector flags in Run Condition Table"};

  Configurable<bool> cfgEnableVertexShiftAnalysis{"cfgEnableVertexShiftAnalysis", false, "Enable the analysis of vertex shift"};
  Configurable<bool> cfgEnableMftDcaAnalysis{"cfgEnableMftDcaAnalysis", false, "Enable the analysis of DCA-based MFT alignment"};
  Configurable<bool> cfgEnableMftDcaExtraPlots{"cfgEnableMftDcaExtraPlots", false, "Enable additional plots for the analysis of DCA-based MFT alignment"};
  Configurable<bool> cfgEnableGlobalFwdDcaAnalysis{"cfgEnableGlobalFwdDcaAnalysis", false, "Enable the analysis of DCA-based MFT alignment using global forward tracks"};
  Configurable<bool> cfgEnableMftMchResidualsAnalysis{"cfgEnableMftMchResidualsAnalysis", true, "Enable the analysis of residuals between MFT tracks and MCH clusters"};
  Configurable<bool> cfgEnableMftMchResidualsExtraPlots{"cfgEnableMftMchResidualsExtraPlots", false, "Enable additional plots for the analysis of residuals between MFT tracks and MCH clusters"};
  Configurable<bool> cfgEnableMftMchMatchingAnalysis{"cfgEnableMftMchMatchingAnalysis", false, "Enable the analysis of residuals between MFT and MCH tracks at reference planes"};
  Configurable<bool> cfgEnableDimuonAnalysis{"cfgEnableDimuonAnalysis", false, "Enable the analysis of di-muon pairs"};

  Configurable<double> cfgRefPlaneZMFT{"cfgRefPlaneZMFT", o2::mft::constants::mft::LayerZCoordinate()[0], "Reference plane on MFT side"};
  Configurable<double> cfgRefPlaneZMCH{"cfgRefPlaneZMCH", -526.0, "Reference plane on MCH side"};

  int mRunNumber{0}; // needed to detect if the run changed and trigger update of magnetic field

  Service<o2::ccdb::BasicCCDBManager> ccdbManager{};
  o2::field::MagneticField* fieldB{nullptr};
  o2::ccdb::CcdbApi ccdbApi;

  // Derived version of mch::Track class that handles the associated clusters as internal objects and deletes them in the destructor
  class TrackRealigned : public mch::Track
  {
   public:
    TrackRealigned() = default;
    ~TrackRealigned()
    {
      // delete the clusters associated to this track
      for (const auto& par : *this) {
        if (par.getClusterPtr()) {
          delete par.getClusterPtr();
        }
      }
    }
  };

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
  std::unordered_map<int64_t, o2::track::TrackParCovFwd> mMftTrackPars;
  std::unordered_map<int64_t, TrackParExt> mMchTrackParsNew;
  std::unordered_map<int64_t, o2::track::TrackParCovFwd> mMftTrackParsNew;

  using MuonPair = std::pair<int64_t, int64_t>;
  using GlobalMuonPair = std::pair<std::vector<int64_t>, std::vector<int64_t>>;

  geo::TransformationCreator transformation;
  std::map<int, math_utils::Transform3D> transformRef; // reference geometry w.r.t track data
  std::map<int, math_utils::Transform3D> transformNew; // new geometry
  TGeoManager* geoNew = nullptr;
  TGeoManager* geoRef = nullptr;
  TrackFitter trackFitter;   // Track fitter from MCH tracking library
  double mImproveCutChi2{0}; // Chi2 cut for track improvement.

  struct AlignmentCorrections {
    double x{0};
    double y{0};
    double z{0};
  };
  std::map<int, AlignmentCorrections> mMchAlignmentCorrections;

  Preslice<aod::FwdTrkCls> perMuon = aod::fwdtrkcl::fwdtrackId;

  o2::aod::rctsel::RCTFlagsChecker rctChecker{"CBT_muon_glo", false, false, true};

  double mBzAtMftCenter{0};

  HistogramRegistry registry{"registry", {}};
  std::array<o2::framework::HistPtr, 10> mMftTrackEffNum;
  std::array<o2::framework::HistPtr, 10> mMftTrackEffDen;

  // vector of all MFT-MCH(-MID) matching candidates associated to the same MCH(-MID) track,
  // to be sorted in descending order with respect to the matching quality
  // the map key is the MCH(-MID) track global index
  using MatchingCandidates = std::map<uint64_t, std::vector<uint64_t>>;

  struct CollisionInfo {
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
    std::map<uint64_t, std::vector<uint64_t>> globalMuonTracks;
  };

  template <typename BC>
  void initCCDB(BC const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }

    mRunNumber = bc.runNumber();
    ccdbManager->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    // auto grpmag = ccdbApi.retrieveFromTFileAny<o2::parameters::GRPMagField>(grpmagPath, metadata, ts);
    auto grpmag = ccdbManager->getForTimeStamp<o2::parameters::GRPMagField>(configCCDB.grpmagPath, bc.timestamp());
    if (grpmag != nullptr) {
      base::Propagator::initFieldFromGRP(grpmag);
      TrackExtrap::setField();
      TrackExtrap::useExtrapV2();
      fieldB = dynamic_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField()); // for MFT
      std::array<double, 3> centerMFT{0, 0, -61.4};                                                 // or use middle point between Vtx and MFT?
      mBzAtMftCenter = fieldB->getBz(centerMFT.data());
    } else {
      LOGF(fatal, "GRP object is not available in CCDB at timestamp=%llu", bc.timestamp());
    }

    // Load geometry information from CCDB/local
    LOGF(info, "Loading reference aligned geometry from CCDB no later than %d", configCCDB.cfgCcdbNoLaterThanRef.value);
    ccdbManager->setCreatedNotAfter(configCCDB.cfgCcdbNoLaterThanRef); // this timestamp has to be consistent with what has been used in reco
    geoRef = ccdbManager->getForTimeStamp<TGeoManager>(configCCDB.geoPath, bc.timestamp());
    ccdbManager->clearCache(configCCDB.geoPath);

    if (configRealign.cfgEnableMCHRealign && (cfgEnableMftMchResidualsAnalysis)) {
      if (geoRef != nullptr) {
        transformation = geo::transformationFromTGeoManager(*geoRef);
      } else {
        LOGF(fatal, "Reference aligned geometry object is not available in CCDB at timestamp=%llu", bc.timestamp());
      }
      for (int i = 0; i < NMchDetElems; i++) {
        int iDEN = GetDetElemId(i);
        transformRef[iDEN] = transformation(iDEN);
      }

      LOGF(info, "Loading new aligned geometry from CCDB no later than %d", configCCDB.cfgCcdbNoLaterThanNew.value);
      ccdbManager->setCreatedNotAfter(configCCDB.cfgCcdbNoLaterThanNew); // make sure this timestamp can be resolved regarding the reference one
      geoNew = ccdbManager->getForTimeStamp<TGeoManager>(configCCDB.geoPathRealign, bc.timestamp());
      ccdbManager->clearCache(configCCDB.geoPathRealign);
      if (geoNew != nullptr) {
        transformation = geo::transformationFromTGeoManager(*geoNew);
      } else {
        LOGF(fatal, "New aligned geometry object is not available in CCDB at timestamp=%llu", bc.timestamp());
      }
      for (int i = 0; i < NMchDetElems; i++) {
        int iDEN = GetDetElemId(i);
        transformNew[iDEN] = transformation(iDEN);
      }
    }
  }

  void init(o2::framework::InitContext&)
  {
    // Load geometry
    ccdbManager->setURL(configCCDB.ccdbUrl);
    ccdbManager->setCaching(true);
    ccdbManager->setLocalObjectValidityChecking();
    ccdbApi.init(configCCDB.ccdbUrl);
    mRunNumber = 0;

    // configure magnetic field position bias
    o2::conf::ConfigurableParam::setValue("FieldOriginBias.z", std::to_string(cfgFieldOriginBiasZ.value));

    // Configuration for track fitter
    const auto& trackerParam = TrackerParam::Instance();
    trackFitter.setBendingVertexDispersion(trackerParam.bendingVertexDispersion);
    trackFitter.setChamberResolution(configRealign.cfgChamberResolutionX, configRealign.cfgChamberResolutionY);
    trackFitter.smoothTracks(true);
    trackFitter.useChamberResolution();
    mImproveCutChi2 = 2. * configRealign.cfgSigmaCutImprove * configRealign.cfgSigmaCutImprove;

    // use the Runge-Kutta extrapolation v2
    TrackExtrap::useExtrapV2();

    // Fill table of MCH alignment corrections
    rapidjson::Document document;
    // Check that the json is parsed correctly
    rapidjson::ParseResult jsonOk = document.Parse(configRealign.cfgMCHRealignCorrections.value.c_str());
    if (jsonOk) {
      for (rapidjson::Value::ConstMemberIterator it = document.MemberBegin(); it != document.MemberEnd(); it++) {
        LOG(info) << "DE" << it->name.GetString() << " alignment corrections:";
        LOG(info) << "  x: " << it->value["x"].GetDouble();
        LOG(info) << "  y: " << it->value["y"].GetDouble();
        LOG(info) << "  z: " << it->value["z"].GetDouble();

        mMchAlignmentCorrections[std::stoi(it->name.GetString())] = AlignmentCorrections{
          .x = it->value["x"].GetDouble(),
          .y = it->value["y"].GetDouble(),
          .z = it->value["z"].GetDouble()};
      }
    } else {
      LOG(error) << "JSON parse error: " << rapidjson::GetParseErrorFunc(jsonOk.Code()) << " (" << jsonOk.Offset() << ")";
    }

    float mftLadderWidth = 1.7;
    AxisSpec dcaxMFTAxis = {400, -0.5, 0.5, "DCA_{x} (cm)"};
    AxisSpec dcayMFTAxis = {400, -0.5, 0.5, "DCA_{y} (cm)"};
    AxisSpec dcaxMCHAxis = {400, -10.0, 10.0, "DCA_{x} (cm)"};
    AxisSpec dcayMCHAxis = {400, -10.0, 10.0, "DCA_{y} (cm)"};
    AxisSpec dcazAxis = {20, -10.0, 10.0, "v_{z} (cm)"};
    AxisSpec txAxis = {30 * 4, -mftLadderWidth * 15.f / 2.f, mftLadderWidth * 15.f / 2.f, "track_{x} (cm)"};
    AxisSpec tyAxis = {24 * 4, -12.f, 12.f, "track_{y} (cm)"};
    AxisSpec txFineAxis = {1500, -15.f, 15.f, "track_{x} (cm)"};
    AxisSpec tyFineAxis = {1500, -15.f, 15.f, "track_{y} (cm)"};
    AxisSpec vxAxis = {400, -0.5, 0.5, "vtx_{x} (cm)"};
    AxisSpec vyAxis = {400, -0.5, 0.5, "vtx_{y} (cm)"};
    AxisSpec vzAxis = {1000, -10.0, 10.0, "vtx_{z} (cm)"};
    AxisSpec phiAxis = {36, -180.0, 180.0, "#phi (degrees)"};
    AxisSpec sxAxis = {50, -0.25, 0.25, "x slope"};
    AxisSpec syAxis = {50, -0.25, 0.25, "y slope"};
    AxisSpec momAxis = {500, 0, 100.0, "p (GeV/c)"};
    AxisSpec nMftClustersAxis = {6, 5, 11, "# of clusters"};
    AxisSpec mftTrackTypeAxis = {2, 0, 2, "track type"};
    AxisSpec mftLayerAxis = {10, 0, 10, "layer"};
    AxisSpec trackChargeSignAxis = {2, 0, 0, "sign"};
    AxisSpec layersPatternAxis = {1024, 0, 1024, "layers pattern"};
    AxisSpec zshiftAxis = {21, -5.25, 5.25, "z shift (mm)"};
    AxisSpec chi2Axis = {500, 0, 500, "chi2"};

    registry.add("vertex_y_vs_x", std::format("Vertex y vs. x").c_str(), {HistType::kTH2F, {vxAxis, vyAxis}});
    registry.add("vertex_z", std::format("Vertex z").c_str(), {HistType::kTH1F, {vzAxis}});

    if (cfgEnableVertexShiftAnalysis || cfgEnableMftDcaAnalysis) {
      registry.add("DCA/MFT/nTracksMFT", std::format("Number of MFT tracks per collision").c_str(), {HistType::kTH1F, {{100, 0, 1000, "# of MFT tracks"}}});
      registry.add("DCA/MFT/DCA_y_vs_x", std::format("DCA y vs. x").c_str(), {HistType::kTH2F, {dcaxMFTAxis, dcayMFTAxis}});
    }

    if (cfgEnableVertexShiftAnalysis) {
      registry.add("DCA/MFT/DCA_x_vs_phi_vs_zshift", std::format("DCA(x) vs. #phi vs. z shift").c_str(), {HistType::kTH3F, {zshiftAxis, phiAxis, dcaxMFTAxis}});
      registry.add("DCA/MFT/DCA_y_vs_phi_vs_zshift", std::format("DCA(y) vs. #phi vs. z shift").c_str(), {HistType::kTH3F, {zshiftAxis, phiAxis, dcayMFTAxis}});

      registry.add("DCA/MFT/DCA_x_vs_slopex_vs_zshift", std::format("DCA(x) vs. x slope vs. z shift").c_str(), {HistType::kTH3F, {zshiftAxis, sxAxis, dcaxMFTAxis}});
      registry.add("DCA/MFT/DCA_x_vs_slopey_vs_zshift", std::format("DCA(x) vs. y slope vs. z shift").c_str(), {HistType::kTH3F, {zshiftAxis, syAxis, dcaxMFTAxis}});
      registry.add("DCA/MFT/DCA_y_vs_slopex_vs_zshift", std::format("DCA(y) vs. x slope vs. z shift").c_str(), {HistType::kTH3F, {zshiftAxis, sxAxis, dcayMFTAxis}});
      registry.add("DCA/MFT/DCA_y_vs_slopey_vs_zshift", std::format("DCA(y) vs. y slope vs. z shift").c_str(), {HistType::kTH3F, {zshiftAxis, syAxis, dcayMFTAxis}});
    }

    if (cfgEnableMftDcaAnalysis) {
      registry.add("DCA/MFT/DCA_x", "DCA(x) vs. vz, tx, ty, nclus",
                   HistType::kTHnSparseF, {dcaxMFTAxis, dcazAxis, txAxis, tyAxis, nMftClustersAxis});
      registry.add("DCA/MFT/DCA_y", "DCA(y) vs. vz, tx, ty, nclus",
                   HistType::kTHnSparseF, {dcayMFTAxis, dcazAxis, txAxis, tyAxis, nMftClustersAxis});

      if (cfgEnableMftDcaExtraPlots) {
        registry.add("DCA/MFT/layers", "Layers vs. tx, ty, nclus",
                     HistType::kTHnSparseF, {mftLayerAxis, txAxis, tyAxis, nMftClustersAxis});
        registry.add("DCA/MFT/trackChi2", "Track #chi^{2} vs. tx, ty, nclus, layer",
                     HistType::kTHnSparseF, {chi2Axis, txAxis, tyAxis, nMftClustersAxis});
        registry.add("DCA/MFT/trackMomentum", "Track momentum vs. tx, ty, nclus, layer",
                     HistType::kTHnSparseF, {momAxis, txAxis, tyAxis, nMftClustersAxis});

        const int nMftLayers = 10;
        for (int i = 0; i < nMftLayers; i++) {
          mMftTrackEffNum[i] = registry.add((std::string("DCA/MFT/mftTrackEffNum_") + std::to_string(i)).c_str(),
                                            (std::string("Track efficiency num. - layer ") + std::to_string(i)).c_str(),
                                            HistType::kTH2F, {txFineAxis, tyFineAxis});
          mMftTrackEffDen[i] = registry.add((std::string("DCA/MFT/mftTrackEffDen_") + std::to_string(i)).c_str(),
                                            (std::string("Track efficiency den. - layer ") + std::to_string(i)).c_str(),
                                            HistType::kTH2F, {txFineAxis, tyFineAxis});
        }
      }
    }

    if (cfgEnableGlobalFwdDcaAnalysis) {
      registry.add("DCA/GlobalFwd/DCA_x", "DCA(x) vs. vz, tx, ty, nclus",
                   HistType::kTHnSparseF, {dcaxMFTAxis, dcazAxis, txAxis, tyAxis, nMftClustersAxis});
      registry.add("DCA/GlobalFwd/DCA_y", "DCA(y) vs. vz, tx, ty, nclus",
                   HistType::kTHnSparseF, {dcayMFTAxis, dcazAxis, txAxis, tyAxis, nMftClustersAxis});
    }

    if (cfgEnableMftMchResidualsAnalysis) {
      AxisSpec dxAxis = {400, -10.0, 10.0, "#Delta x (cm)"};
      AxisSpec dyAxis = {400, -10.0, 10.0, "#Delta y (cm)"};

      registry.add("DCA/MCH/DCA_y_vs_x", std::format("DCA y vs. x").c_str(), {HistType::kTH2F, {dcaxMCHAxis, dcayMCHAxis}});
      registry.add("DCA/MCH/DCA_x_vs_sign_vs_quadrant_vs_mom", std::format("DCA(x) vs. p, quadrant, chargeSign").c_str(), {HistType::kTHnSparseF, {{20, 0, 100.0, "p (GeV/c)"}, {4, 0, 4, "quadrant"}, {2, 0, 2, "sign"}, dcaxMCHAxis}});
      registry.add("DCA/MCH/DCA_y_vs_sign_vs_quadrant_vs_mom", std::format("DCA(y) vs. p, quadrant, chargeSign").c_str(), {HistType::kTHnSparseF, {{20, 0, 100.0, "p (GeV/c)"}, {4, 0, 4, "quadrant"}, {2, 0, 2, "sign"}, dcayMCHAxis}});

      registry.add("DCA/MCH/DCA_x_vs_sign_vs_quadrant_vs_mom_corr", std::format("DCA(x) vs. p, quadrant, chargeSign (with corrections)").c_str(), {HistType::kTHnSparseF, {{20, 0, 100.0, "p (GeV/c)"}, {4, 0, 4, "quadrant"}, {2, 0, 2, "sign"}, dcaxMCHAxis}});
      registry.add("DCA/MCH/DCA_y_vs_sign_vs_quadrant_vs_mom_corr", std::format("DCA(y) vs. p, quadrant, chargeSign (with corrections)").c_str(), {HistType::kTHnSparseF, {{20, 0, 100.0, "p (GeV/c)"}, {4, 0, 4, "quadrant"}, {2, 0, 2, "sign"}, dcayMCHAxis}});

      registry.add("residuals/dx_vs_chamber", "Cluster x residual vs. chamber, quadrant, chargeSign",
                   {HistType::kTHnSparseF, {{10, 1, 11, "chamber"}, {4, 0, 4, "quadrant"}, {2, 0, 2, "sign"}, dxAxis}});
      registry.add("residuals/dy_vs_chamber", "Cluster y residual vs. chamber, quadrant, chargeSign",
                   {HistType::kTHnSparseF, {{10, 1, 11, "chamber"}, {4, 0, 4, "quadrant"}, {2, 0, 2, "sign"}, dyAxis}});

      registry.add("residuals/dx_vs_de", "Cluster x residual vs. DE, quadrant, chargeSign, momentum, xslope",
                   {HistType::kTHnSparseF, {dxAxis, {getNumDE(), 0, static_cast<double>(getNumDE()), "DE"}, {4, 0, 4, "quadrant"}, {2, 0, 2, "sign"}, {20, 0, 100.0, "p (GeV/c)"}, sxAxis}});
      registry.add("residuals/dy_vs_de", "Cluster y residual vs. DE, quadrant, chargeSign, momentum, yslope",
                   {HistType::kTHnSparseF, {dyAxis, {getNumDE(), 0, static_cast<double>(getNumDE()), "DE"}, {4, 0, 4, "quadrant"}, {2, 0, 2, "sign"}, {20, 0, 100.0, "p (GeV/c)"}, syAxis}});

      registry.add("residuals/dx_vs_chamber_corr", "Cluster x residual vs. chamber, quadrant, chargeSign (with corrections)",
                   {HistType::kTHnSparseF, {{10, 1, 11, "chamber"}, {4, 0, 4, "quadrant"}, {2, 0, 2, "sign"}, dxAxis}});
      registry.add("residuals/dy_vs_chamber_corr", "Cluster y residual vs. chamber, quadrant, chargeSign (with corrections)",
                   {HistType::kTHnSparseF, {{10, 1, 11, "chamber"}, {4, 0, 4, "quadrant"}, {2, 0, 2, "sign"}, dyAxis}});

      registry.add("residuals/dx_vs_de_corr", "Cluster x residual vs. DE, quadrant, chargeSign, momentum, xslope (with corrections)",
                   {HistType::kTHnSparseF, {dxAxis, {getNumDE(), 0, static_cast<double>(getNumDE()), "DE"}, {4, 0, 4, "quadrant"}, {2, 0, 2, "sign"}, {20, 0, 100.0, "p (GeV/c)"}, sxAxis}});
      registry.add("residuals/dy_vs_de_corr", "Cluster y residual vs. DE, quadrant, chargeSign, momentum, yslope (with corrections)",
                   {HistType::kTHnSparseF, {dyAxis, {getNumDE(), 0, static_cast<double>(getNumDE()), "DE"}, {4, 0, 4, "quadrant"}, {2, 0, 2, "sign"}, {20, 0, 100.0, "p (GeV/c)"}, syAxis}});

      registry.add("residuals/de_alignment_corrections_x", "DE alignment corrections - X coordinate",
                   {HistType::kTH1F, {{getNumDE(), 0, static_cast<double>(getNumDE()), "DE"}}});
      registry.add("residuals/de_alignment_corrections_y", "DE alignment corrections - Y coordinate",
                   {HistType::kTH1F, {{getNumDE(), 0, static_cast<double>(getNumDE()), "DE"}}});

      for (const auto& [deId, corr] : mMchAlignmentCorrections) {
        auto deIndex = getDEindex(deId);
        registry.get<TH1>(HIST("residuals/de_alignment_corrections_x"))->SetBinContent(deIndex + 1, corr.x);
        registry.get<TH1>(HIST("residuals/de_alignment_corrections_x"))->SetBinError(deIndex + 1, 0.1);
        registry.get<TH1>(HIST("residuals/de_alignment_corrections_y"))->SetBinContent(deIndex + 1, corr.y);
        registry.get<TH1>(HIST("residuals/de_alignment_corrections_y"))->SetBinError(deIndex + 1, 0.1);
      }

      if (cfgEnableMftMchResidualsExtraPlots) {
        registry.add("DCA/MCH/DCA_x_vs_sign_vs_quadrant_vs_vz", std::format("DCA(x) vs. vz, quadrant, chargeSign").c_str(), {HistType::kTHnSparseF, {dcazAxis, {4, 0, 4, "quadrant"}, {2, 0, 2, "sign"}, dcaxMCHAxis}});
        registry.add("DCA/MCH/DCA_y_vs_sign_vs_quadrant_vs_vz", std::format("DCA(y) vs. vz, quadrant, chargeSign").c_str(), {HistType::kTHnSparseF, {dcazAxis, {4, 0, 4, "quadrant"}, {2, 0, 2, "sign"}, dcayMCHAxis}});

        registry.add("residuals/dphi_at_mft", "Track #Delta#phi at MFT",
                     {HistType::kTHnSparseF, {{200, -0.2f, 0.2f, "#Delta#phi"}, {80, -10.f, 10.f, "track_x (cm)"}, {80, -10.f, 10.f, "track_y (cm)"}, {2, 0, 2, "sign"}, {20, 0, 100.0, "p (GeV/c)"}}});
      }
    }

    if (cfgEnableMftMchMatchingAnalysis) {
      AxisSpec dxAxis = {200, -10.0, 10.0, "#Deltax (cm)"};
      AxisSpec dyAxis = {200, -10.0, 10.0, "#Deltay (cm)"};
      AxisSpec dsxAxis = {200, -0.1, 0.1, "#Deltaslope(x) (rad)"};
      AxisSpec dsyAxis = {200, -0.1, 0.1, "#Deltaslope(y) (rad)"};
      AxisSpec dphiAxis = {200, -1.0, 1.0, "#Delta#phi (rad)"};

      // MFT plane
      registry.add("matching/dxAtMFT", "Tracks #Deltax at MFT reference plane",
                   {HistType::kTHnSparseF, {dxAxis, {20, -15.0, 15.0, "track x (cm)"}, {20, -15.0, 15.0, "track y (cm)"}, {4, 0, 4, "quadrant"}, {2, 0, 2, "sign"}, {20, 0, 100.0, "p (GeV/c)"}}});
      registry.add("matching/dyAtMFT", "Tracks #Deltay at MFT reference plane",
                   {HistType::kTHnSparseF, {dyAxis, {20, -15.0, 15.0, "track x (cm)"}, {20, -15.0, 15.0, "track y (cm)"}, {4, 0, 4, "quadrant"}, {2, 0, 2, "sign"}, {20, 0, 100.0, "p (GeV/c)"}}});
      registry.add("matching/dsxAtMFT", "Tracks #Deltaslope(x) at MFT reference plane",
                   {HistType::kTHnSparseF, {dsxAxis, {20, -15.0, 15.0, "track x (cm)"}, {20, -15.0, 15.0, "track y (cm)"}, {4, 0, 4, "quadrant"}, {2, 0, 2, "sign"}, {20, 0, 100.0, "p (GeV/c)"}}});
      registry.add("matching/dsyAtMFT", "Tracks #Deltaslope(y) at MFT reference plane",
                   {HistType::kTHnSparseF, {dsyAxis, {20, -15.0, 15.0, "track x (cm)"}, {20, -15.0, 15.0, "track y (cm)"}, {4, 0, 4, "quadrant"}, {2, 0, 2, "sign"}, {20, 0, 100.0, "p (GeV/c)"}}});
      registry.add("matching/dphiAtMFT", "Tracks #Delta#phi at MFT reference plane",
                   {HistType::kTHnSparseF, {dphiAxis, {20, -15.0, 15.0, "track x (cm)"}, {20, -15.0, 15.0, "track y (cm)"}, {4, 0, 4, "quadrant"}, {2, 0, 2, "sign"}, {20, 0, 100.0, "p (GeV/c)"}}});

      // MCH plane
      registry.add("matching/dxAtMCH", "Tracks #Deltax at MCH reference plane",
                   {HistType::kTHnSparseF, {dxAxis, {20, -100.0, 100.0, "track x (cm)"}, {20, -100.0, 100.0, "track y (cm)"}, {4, 0, 4, "quadrant"}, {2, 0, 2, "sign"}, {20, 0, 100.0, "p (GeV/c)"}}});
      registry.add("matching/dyAtMCH", "Tracks #Deltay at MCH reference plane",
                   {HistType::kTHnSparseF, {dyAxis, {20, -100.0, 100.0, "track x (cm)"}, {20, -100.0, 100.0, "track y (cm)"}, {4, 0, 4, "quadrant"}, {2, 0, 2, "sign"}, {20, 0, 100.0, "p (GeV/c)"}}});
      registry.add("matching/dsxAtMCH", "Tracks #Deltaslope(x) at MCH reference plane",
                   {HistType::kTHnSparseF, {dsxAxis, {20, -100.0, 100.0, "track x (cm)"}, {20, -100.0, 100.0, "track y (cm)"}, {4, 0, 4, "quadrant"}, {2, 0, 2, "sign"}, {20, 0, 100.0, "p (GeV/c)"}}});
      registry.add("matching/dsyAtMCH", "Tracks #Deltaslope(y) at MCH reference plane",
                   {HistType::kTHnSparseF, {dsyAxis, {20, -100.0, 100.0, "track x (cm)"}, {20, -100.0, 100.0, "track y (cm)"}, {4, 0, 4, "quadrant"}, {2, 0, 2, "sign"}, {20, 0, 100.0, "p (GeV/c)"}}});
      registry.add("matching/dphiAtMCH", "Tracks #Delta#phi at MCH reference plane",
                   {HistType::kTHnSparseF, {dphiAxis, {20, -100.0, 100.0, "track x (cm)"}, {20, -100.0, 100.0, "track y (cm)"}, {4, 0, 4, "quadrant"}, {2, 0, 2, "sign"}, {20, 0, 100.0, "p (GeV/c)"}}});
    }

    if (cfgEnableDimuonAnalysis) {
      AxisSpec invMassAxis = {1500, 0, 15, "M_{#mu^{+}#mu^{-}} (GeV/c^{2})"};
      AxisSpec pTAxis = {30, 0, 30, "#mu^{+}#mu^{-} p_{T} (GeV/c)"};
      AxisSpec pAxis = {50, 0, 200, "#mu^{+}#mu^{-} p (GeV/c)"};
      AxisSpec muPosQuadrantAxis = {4, 0, 4, "#mu^{+} quadrant"};
      AxisSpec muNegQuadrantAxis = {4, 0, 4, "#mu^{-} quadrant"};
      AxisSpec angleAxis = {100, 0, 0.5, "#mu^{+}#mu^{-} angle (rad)"};
      AxisSpec angleDiffAxis = {500, -0.05, 0.05, "#mu^{+}#mu^{-} angle difference (rad)"};
      AxisSpec dcaxAxis = {400, -10.0, 10.0, "#mu^{+}#mu^{-} DCA_{x} (cm)"};
      AxisSpec dcayAxis = {400, -10.0, 10.0, "#mu^{+}#mu^{-} DCA_{y} (cm)"};
      AxisSpec dcaMftxAxis = {400, -0.5, 0.5, "#mu^{+}#mu^{-} DCA_{x} (cm)"};
      AxisSpec dcaMftyAxis = {400, -0.5, 0.5, "#mu^{+}#mu^{-} DCA_{y} (cm)"};

      // di-muon invariant mass distributions
      registry.add("dimuon/invariantMass_MuonKine_MuonCuts", "#mu^{+}#mu^{-} invariant mass (muon cuts)", {HistType::kTHnSparseF, {invMassAxis, pAxis, pTAxis, muPosQuadrantAxis, muNegQuadrantAxis}});
      registry.add("dimuon/invariantMass_MuonKine_GlobalMuonCuts_GoodMatches", "#mu^{+}#mu^{-} invariant mass (global muon cuts, good matches)", {HistType::kTHnSparseF, {invMassAxis, pAxis, pTAxis, muPosQuadrantAxis, muNegQuadrantAxis}});
      registry.add("dimuon/invariantMass_ScaledMftKine_GlobalMuonCuts_GoodMatches", "#mu^{+}#mu^{-} invariant mass (global muon cuts, rescaled MFT, good matches)", {HistType::kTHnSparseF, {invMassAxis, pAxis, pTAxis, muPosQuadrantAxis, muNegQuadrantAxis}});
      // difference in mu+mu- opening angle between MCH and global muon tracks
      registry.add("dimuon/angle_GlobalMuonCuts_GoodMatches", "#mu^{+}#mu^{-} opening angle difference (global muon cuts, good matches)", {HistType::kTHnSparseF, {angleDiffAxis, angleAxis, pAxis, pTAxis, muPosQuadrantAxis, muNegQuadrantAxis}});
      // mu+mu- DCA
      registry.add("dimuon/dcax_MuonKine_MuonCuts", "#mu^{+}#mu^{-} DCA_{x} (muon cuts)", {HistType::kTHnSparseF, {dcaxAxis, pAxis, pTAxis, muPosQuadrantAxis, muNegQuadrantAxis}});
      registry.add("dimuon/dcay_MuonKine_MuonCuts", "#mu^{+}#mu^{-} DCA_{y} (muon cuts)", {HistType::kTHnSparseF, {dcayAxis, pAxis, pTAxis, muPosQuadrantAxis, muNegQuadrantAxis}});
      registry.add("dimuon/dcax_MuonKine_GlobalMuonCuts_GoodMatches", "#mu^{+}#mu^{-} DCA_{x} (global muon cuts, good matches)", {HistType::kTHnSparseF, {dcaxAxis, pAxis, pTAxis, muPosQuadrantAxis, muNegQuadrantAxis}});
      registry.add("dimuon/dcay_MuonKine_GlobalMuonCuts_GoodMatches", "#mu^{+}#mu^{-} DCA_{y} (global muon cuts, good matches)", {HistType::kTHnSparseF, {dcayAxis, pAxis, pTAxis, muPosQuadrantAxis, muNegQuadrantAxis}});
      registry.add("dimuon/dcax_ScaledMftKine_GlobalMuonCuts_GoodMatches", "#mu^{+}#mu^{-} DCA_{x} (global muon cuts, rescaled MFT, good matches)", {HistType::kTHnSparseF, {dcaMftxAxis, pAxis, pTAxis, muPosQuadrantAxis, muNegQuadrantAxis}});
      registry.add("dimuon/dcay_ScaledMftKine_GlobalMuonCuts_GoodMatches", "#mu^{+}#mu^{-} DCA_{y} (global muon cuts, rescaled MFT, good matches)", {HistType::kTHnSparseF, {dcaMftyAxis, pAxis, pTAxis, muPosQuadrantAxis, muNegQuadrantAxis}});

      // di-muon invariant mass distributions (realigned/refitted tracks)
      registry.add("dimuon/realign/invariantMass_MuonKine_MuonCuts", "#mu^{+}#mu^{-} invariant mass (muon cuts)", {HistType::kTHnSparseF, {invMassAxis, pAxis, pTAxis, muPosQuadrantAxis, muNegQuadrantAxis}});
      registry.add("dimuon/realign/invariantMass_MuonKine_GlobalMuonCuts_GoodMatches", "#mu^{+}#mu^{-} invariant mass (global muon cuts, good matches)", {HistType::kTHnSparseF, {invMassAxis, pAxis, pTAxis, muPosQuadrantAxis, muNegQuadrantAxis}});
      registry.add("dimuon/realign/invariantMass_ScaledMftKine_GlobalMuonCuts_GoodMatches", "#mu^{+}#mu^{-} invariant mass (global muon cuts, rescaled MFT, good matches)", {HistType::kTHnSparseF, {invMassAxis, pAxis, pTAxis, muPosQuadrantAxis, muNegQuadrantAxis}});
      // difference in mu+mu- opening angle between MCH and global muon tracks (realigned/refitted tracks)
      registry.add("dimuon/realign/angle_GlobalMuonCuts_GoodMatches", "#mu^{+}#mu^{-} opening angle difference (global muon cuts, good matches)", {HistType::kTHnSparseF, {angleDiffAxis, angleAxis, pAxis, pTAxis, muPosQuadrantAxis, muNegQuadrantAxis}});
      // mu+mu- DCA
      registry.add("dimuon/realign/dcax_MuonKine_MuonCuts", "#mu^{+}#mu^{-} DCA_{x} (muon cuts)", {HistType::kTHnSparseF, {dcaxAxis, pAxis, pTAxis, muPosQuadrantAxis, muNegQuadrantAxis}});
      registry.add("dimuon/realign/dcay_MuonKine_MuonCuts", "#mu^{+}#mu^{-} DCA_{y} (muon cuts)", {HistType::kTHnSparseF, {dcayAxis, pAxis, pTAxis, muPosQuadrantAxis, muNegQuadrantAxis}});
      registry.add("dimuon/realign/dcax_MuonKine_GlobalMuonCuts_GoodMatches", "#mu^{+}#mu^{-} DCA_{x} (global muon cuts, good matches)", {HistType::kTHnSparseF, {dcaxAxis, pAxis, pTAxis, muPosQuadrantAxis, muNegQuadrantAxis}});
      registry.add("dimuon/realign/dcay_MuonKine_GlobalMuonCuts_GoodMatches", "#mu^{+}#mu^{-} DCA_{y} (global muon cuts, good matches)", {HistType::kTHnSparseF, {dcayAxis, pAxis, pTAxis, muPosQuadrantAxis, muNegQuadrantAxis}});
      registry.add("dimuon/realign/dcax_ScaledMftKine_GlobalMuonCuts_GoodMatches", "#mu^{+}#mu^{-} DCA_{x} (global muon cuts, rescaled MFT, good matches)", {HistType::kTHnSparseF, {dcaMftxAxis, pAxis, pTAxis, muPosQuadrantAxis, muNegQuadrantAxis}});
      registry.add("dimuon/realign/dcay_ScaledMftKine_GlobalMuonCuts_GoodMatches", "#mu^{+}#mu^{-} DCA_{y} (global muon cuts, rescaled MFT, good matches)", {HistType::kTHnSparseF, {dcaMftyAxis, pAxis, pTAxis, muPosQuadrantAxis, muNegQuadrantAxis}});
    }
  }

  int GetDetElemId(int iDetElemNumber)
  {
    constexpr int fgNCh = 10;
    const std::array<int, fgNCh> fgNDetElemCh{4, 4, 4, 4, 18, 18, 26, 26, 26, 26};
    const std::array<int, fgNCh + 1> fgSNDetElemCh{0, 4, 8, 12, 16, 34, 52, 78, 104, 130, 156};

    // make sure detector number is valid
    if (iDetElemNumber < fgSNDetElemCh[0] ||
        iDetElemNumber >= fgSNDetElemCh[10]) {
      LOGF(fatal, "Invalid detector element number: %d", iDetElemNumber);
    }
    /// get det element number from ID
    // get chamber and element number in chamber
    int iCh = 0;
    int iDet = 0;
    for (int i = 1; i <= NMchChambers; i++) {
      if (iDetElemNumber < fgSNDetElemCh[i]) {
        iCh = i;
        iDet = iDetElemNumber - fgSNDetElemCh[i - 1];
        break;
      }
    }

    // make sure detector index is valid
    if (iCh <= 0 || iCh > NMchChambers || iDet >= fgNDetElemCh[iCh - 1]) {
      LOGF(fatal, "Invalid detector element id: %d", 100 * iCh + iDet);
    }

    // add number of detectors up to this chamber
    return 100 * iCh + iDet;
  }

  int getChamberIndex(int deId)
  {
    return (deId / 100) - 1;
  }

  int getNumDEinChamber(int chIndex)
  {
    int nDE = 0;
    switch (chIndex) {
      case 0:
      case 1:
      case 2:
      case 3:
        nDE = 4;
        break;
      case 4:
      case 5:
        nDE = 18;
        break;
      case 6:
      case 7:
      case 8:
      case 9:
        nDE = 26;
        break;
      default:
        break;
    }
    return nDE;
  }

  int getNumDE()
  {
    static int nDE = 0;
    if (nDE <= 0) {
      for (int c = 0; c < NMchChambers; c++) {
        nDE += getNumDEinChamber(c);
      }
    }

    return nDE;
  }

  int getDEindexInChamber(int deId)
  {
    return (deId - 100) % 100;
  }

  int getChamberOffset(int chIndex)
  {
    int offset = 0;
    for (int c = 0; c < chIndex; c++) {
      offset += getNumDEinChamber(c);
    }
    return offset;
  }

  int getDEindex(int deId)
  {
    auto idx = getDEindexInChamber(deId);
    int offset = getChamberOffset(getChamberIndex(deId));

    return idx + offset;
  }

  int GetQuadrant(float phi)
  {
    if (phi >= 0 && phi < o2::constants::math::PIHalf) {
      return 0;
    }
    if (phi >= o2::constants::math::PIHalf && phi <= o2::constants::math::PI) {
      return 1;
    }
    if (phi >= -o2::constants::math::PI && phi < -o2::constants::math::PIHalf) {
      return 2;
    }
    if (phi >= -o2::constants::math::PIHalf && phi < 0) {
      return 3;
    }
    return -1;
  }

  template <class T>
  int GetQuadrant(const T& track)
  {
    return GetQuadrant(static_cast<float>(track.phi()));
  }

  template <class T>
  static o2::mch::TrackParam FwdtoMCH(const T& fwdtrack)
  {
    using SMatrix55Std = ROOT::Math::SMatrix<double, 5>;
    using SMatrix55Sym = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;

    // Convert Forward Track parameters and covariances matrix to the
    // MCH track format.

    // Parameter conversion
    double x2 = fwdtrack.getPhi();
    double x3 = fwdtrack.getTanl();
    double x4 = fwdtrack.getInvQPt();

    auto sinx2 = TMath::Sin(x2);
    auto cosx2 = TMath::Cos(x2);

    double alpha1 = cosx2 / x3;
    double alpha3 = sinx2 / x3;
    double alpha4 = x4 / TMath::Sqrt(x3 * x3 + sinx2 * sinx2);

    auto K = TMath::Sqrt(x3 * x3 + sinx2 * sinx2);
    auto K3 = K * K * K;

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

    jacobian(1, 2) = -sinx2 / x3;
    jacobian(1, 3) = -cosx2 / (x3 * x3);

    jacobian(2, 1) = 1;

    jacobian(3, 2) = cosx2 / x3;
    jacobian(3, 3) = -sinx2 / (x3 * x3);

    jacobian(4, 2) = -x4 * sinx2 * cosx2 / K3;
    jacobian(4, 3) = -x3 * x4 / K3;
    jacobian(4, 4) = 1 / K;
    // jacobian*covariances*jacobian^T
    covariances = ROOT::Math::Similarity(jacobian, covariances);

    const std::array<Double_t, 15> cov{covariances(0, 0), covariances(1, 0), covariances(1, 1), covariances(2, 0), covariances(2, 1), covariances(2, 2), covariances(3, 0), covariances(3, 1), covariances(3, 2), covariances(3, 3), covariances(4, 0), covariances(4, 1), covariances(4, 2), covariances(4, 3), covariances(4, 4)};
    const std::array<Double_t, 5> param{fwdtrack.getX(), alpha1, fwdtrack.getY(), alpha3, alpha4};

    o2::mch::TrackParam convertedTrack(fwdtrack.getZ(), param.data(), cov.data());
    return {convertedTrack};
  }

  static o2::dataformats::GlobalFwdTrack MCHtoFwd(const o2::mch::TrackParam& mchParam)
  {
    using SMatrix55Std = ROOT::Math::SMatrix<double, 5>;
    using SMatrix55Sym = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;

    // Convert a MCH Track parameters and covariances matrix to the
    // Forward track format. Must be called after propagation though the absorber

    o2::dataformats::GlobalFwdTrack convertedTrack;

    // Parameter conversion
    double alpha1 = mchParam.getNonBendingSlope();
    double alpha3 = mchParam.getBendingSlope();
    double alpha4 = mchParam.getInverseBendingMomentum();

    double x2 = TMath::ATan2(-alpha3, -alpha1);
    double x3 = -1. / TMath::Sqrt(alpha3 * alpha3 + alpha1 * alpha1);
    double x4 = alpha4 * -x3 * TMath::Sqrt(1 + alpha3 * alpha3);

    auto K = alpha1 * alpha1 + alpha3 * alpha3;
    auto K32 = K * TMath::Sqrt(K);
    auto L = TMath::Sqrt(alpha3 * alpha3 + 1);

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

    jacobian(2, 1) = -alpha3 / K;
    jacobian(2, 3) = alpha1 / K;

    jacobian(3, 1) = alpha1 / K32;
    jacobian(3, 3) = alpha3 / K32;

    jacobian(4, 1) = -alpha1 * alpha4 * L / K32;
    jacobian(4, 3) = alpha3 * alpha4 * (1 / (TMath::Sqrt(K) * L) - L / K32);
    jacobian(4, 4) = L / TMath::Sqrt(K);

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

  bool RemoveTrack(mch::Track& track)
  {
    // Refit track with re-aligned clusters
    bool removeTrack = false;
    try {
      trackFitter.fit(track, false);
    } catch (std::exception const& e) {
      removeTrack = true;
      return removeTrack;
    }

    auto itStartingParam = std::prev(track.rend());

    while (true) {

      try {
        trackFitter.fit(track, true, false, (itStartingParam == track.rbegin()) ? nullptr : &itStartingParam);
      } catch (std::exception const&) {
        removeTrack = true;
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
        removeTrack = true;
        track.removable();
        break;
      }

      auto itNextParam = track.removeParamAtCluster(itWorstParam);
      auto itNextToNextParam = (itNextParam == track.end()) ? itNextParam : std::next(itNextParam);
      itStartingParam = track.rbegin();

      if (track.getNClusters() < NMchChambers) {
        removeTrack = true;
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

    if (!removeTrack) {
      for (auto& param : track) { // o2-linter: disable=const-ref-in-for-loop (object is modified in loop)
        param.setParameters(param.getSmoothParameters());
        param.setCovariances(param.getSmoothCovariances());
      }
    }

    return removeTrack;
  }

  template <class T>
  bool IsGoodMFT(const T& mftTrack,
                 double chi2Cut,
                 int nClustersCut)
  {
    // chi2 cut
    if (mftTrack.chi2() > chi2Cut) {
      return false;
    }

    // number of clusters cut
    if (mftTrack.nClusters() < nClustersCut) {
      return false;
    }

    return true;
  }

  template <class T>
  bool IsGoodMFT(const T& mftTrack)
  {
    return IsGoodMFT(mftTrack, cfgTrackChi2MftUp, cfgTrackNClustMftLow);
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
    double sigmaPDCA = (thetaAbs < ThetaAbsBoundaryDeg) ? sigmaPDCA23 : sigmaPDCA310;
    double nrp = nSigmaPDCA * relPRes * p;
    double pResEffect = sigmaPDCA / (1. - nrp / (1. + nrp));
    double slopeResEffect = SlopeResolutionZ * slopeRes * p;
    double sigmaPDCAWithRes = TMath::Sqrt(pResEffect * pResEffect + slopeResEffect * slopeResEffect);

    return (pDCA <= nSigmaPDCA * sigmaPDCAWithRes);
  }

  template <class T, class C>
  bool IsGoodMuon(const T& muonTrack, const C& collision,
                  double chi2Cut,
                  double pCut,
                  double pTCut,
                  std::array<double, 2> etaCut,
                  std::array<double, 2> rAbsCut,
                  double nSigmaPdcaCut)
  {
    auto const& mchTrack = (static_cast<int>(muonTrack.trackType()) <= GlobalTrackTypeMax) ? muonTrack.template matchMCHTrack_as<MyMuonsWithCov>() : muonTrack;

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
    if (!pDCACut(mchTrack, collision, nSigmaPdcaCut)) {
      return false;
    }

    return true;
  }

  template <typename TMUON>
  bool isGoodGlobalMatching(const TMUON& muon,
                            double matchChi2Cut)
  {
    if (static_cast<int>(muon.trackType()) > GlobalTrackTypeMax) {
      return false;
    }

    // MFT-MCH match chi2 cut
    if (muon.chi2MatchMCHMFT() > matchChi2Cut) {
      return false;
    }

    return true;
  }

  template <typename T>
  o2::track::TrackParCovFwd TrackToParCovFwd(const T& track)
  {
    double chi2 = track.chi2();
    SMatrix5 tpars(track.x(), track.y(), track.phi(), track.tgl(), track.signed1Pt());
    std::vector<double> v1{0, 0, 0, 0, 0,
                           0, 0, 0, 0, 0,
                           0, 0, 0, 0, 0};
    SMatrix55 tcovs(v1.begin(), v1.end());
    o2::track::TrackParCovFwd trackparCov{track.z(), tpars, tcovs, chi2};
    return trackparCov;
  }

  template <typename T>
  o2::dataformats::GlobalFwdTrack TrackToGlobalFwd(const T& track)
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

  void TransformMFTPar(o2::mch::TrackParam& track)
  {
    // double zCH10 = -1437.6;
    double z = track.getZ();
    // double dZ = zMCH - z;
    double x = track.getNonBendingCoor();
    double y = track.getBendingCoor();
    double xSlope = track.getNonBendingSlope();
    double ySlope = track.getBendingSlope();

    double xSlopeCorrection = (y > 0) ? configMFTAlignmentCorrections.cfgMFTAlignmentCorrXSlopeTop : configMFTAlignmentCorrections.cfgMFTAlignmentCorrXSlopeBottom;
    double xCorrection = xSlopeCorrection * z +
                         ((y > 0) ? configMFTAlignmentCorrections.cfgMFTAlignmentCorrXOffsetTop : configMFTAlignmentCorrections.cfgMFTAlignmentCorrXOffsetBottom);
    track.setNonBendingCoor(x + xCorrection);
    track.setNonBendingSlope(xSlope + xSlopeCorrection);

    double ySlopeCorrection = (y > 0) ? configMFTAlignmentCorrections.cfgMFTAlignmentCorrYSlopeTop : configMFTAlignmentCorrections.cfgMFTAlignmentCorrYSlopeBottom;
    double yCorrection = ySlopeCorrection * z +
                         ((y > 0) ? configMFTAlignmentCorrections.cfgMFTAlignmentCorrYOffsetTop : configMFTAlignmentCorrections.cfgMFTAlignmentCorrYOffsetBottom);
    track.setBendingCoor(y + yCorrection);
    track.setBendingSlope(ySlope + ySlopeCorrection);
  }

  void TransformMFT(o2::track::TrackParCovFwd& fwdtrack)
  {
    auto mchTrack = FwdtoMCH(fwdtrack);

    TransformMFTPar(mchTrack);

    auto transformedTrack = MCHtoFwd(mchTrack);
    fwdtrack.setParameters(transformedTrack.getParameters());
    fwdtrack.setZ(transformedTrack.getZ());
    fwdtrack.setCovariances(transformedTrack.getCovariances());
  }

  template <typename T>
  T UpdateTrackMomentum(const T& track, const double p, int sign)
  {
    double px = p * std::sin(o2::constants::math::PIHalf - std::atan(track.tgl())) * std::cos(track.phi());
    double py = p * std::sin(o2::constants::math::PIHalf - std::atan(track.tgl())) * std::sin(track.phi());
    double pt = std::sqrt(std::pow(px, 2) + std::pow(py, 2));

    SMatrix5 tpars = {track.x(), track.y(), track.phi(), track.tgl(), sign / pt};
    std::vector<double> v1{0, 0, 0, 0, 0,
                           0, 0, 0, 0, 0,
                           0, 0, 0, 0, 0};
    SMatrix55 tcovs(v1.begin(), v1.end());

    T newTrack;
    newTrack.setParameters(tpars);
    newTrack.setZ(track.z());
    newTrack.setCovariances(tcovs);

    return newTrack;
  }

  template <typename T>
  T UpdateTrackMomentum(const T& track, const o2::mch::TrackParam& track4mom)
  {
    double px = track4mom.p() * std::sin(o2::constants::math::PIHalf - std::atan(track.tgl())) * std::cos(track.phi());
    double py = track4mom.p() * std::sin(o2::constants::math::PIHalf - std::atan(track.tgl())) * std::sin(track.phi());
    double pt = std::sqrt(std::pow(px, 2) + std::pow(py, 2));
    double sign = track4mom.getCharge();

    SMatrix5 tpars = {track.x(), track.y(), track.phi(), track.tgl(), sign / pt};
    std::vector<double> v1{0, 0, 0, 0, 0,
                           0, 0, 0, 0, 0,
                           0, 0, 0, 0, 0};
    SMatrix55 tcovs(v1.begin(), v1.end());

    T newTrack;
    newTrack.setParameters(tpars);
    newTrack.setZ(track.z());
    newTrack.setCovariances(tcovs);

    return track;
  }

  void UpdateTrackMomentum(o2::mch::TrackParam& track, const o2::mch::TrackParam& track4mom)
  {
    double pRatio = track.p() / track4mom.p();
    double newInvBendMom = track.getInverseBendingMomentum() * pRatio;
    track.setInverseBendingMomentum(newInvBendMom);
    track.setCharge(track4mom.getCharge());
  }

  void UpdateTrackMomentum(o2::track::TrackParCovFwd& track, const o2::mch::TrackParam& track4mom)
  {
    auto newTrackMCH = FwdtoMCH(track);
    UpdateTrackMomentum(newTrackMCH, track4mom);
    auto newTrack = MCHtoFwd(newTrackMCH);

    track.setParameters(newTrack.getParameters());
    track.setZ(newTrack.getZ());
    track.setCovariances(newTrack.getCovariances());
  }

  o2::dataformats::GlobalFwdTrack PropagateMCHParam(mch::TrackParam mchTrack, const double z)
  {
    float absFront = -90.f;
    float absBack = -505.f;

    if (mchTrack.getZ() < absBack && z > absFront) {
      // extrapolation through the absorber in the upstream direction
      o2::mch::TrackExtrap::extrapToVertexWithoutBranson(mchTrack, z);
    } else if (z < absBack) {
      // extrapolation downstream of the absorber, correct for dipole longitudinal shift if needed
      if (cfgDipoleZshift.value != 0) {
        mchTrack.setZ(mchTrack.getZ() + cfgDipoleZshift.value);
        o2::mch::TrackExtrap::extrapToZCov(mchTrack, z + cfgDipoleZshift.value);
        mchTrack.setZ(mchTrack.getZ() - cfgDipoleZshift.value);
      } else {
        o2::mch::TrackExtrap::extrapToZCov(mchTrack, z);
      }
    } else {
      // all other cases
      o2::mch::TrackExtrap::extrapToZCov(mchTrack, z);
    }

    auto proptrack = MCHtoFwd(mchTrack);
    o2::dataformats::GlobalFwdTrack propmuon;
    propmuon.setParameters(proptrack.getParameters());
    propmuon.setZ(proptrack.getZ());
    propmuon.setCovariances(proptrack.getCovariances());

    return propmuon;
  }

  o2::dataformats::GlobalFwdTrack PropagateMCH(const o2::dataformats::GlobalFwdTrack& muon, const double z)
  {
    auto mchTrack = FwdtoMCH(muon);
    return PropagateMCHParam(mchTrack, z);
  }

  o2::dataformats::GlobalFwdTrack PropagateMCHRealigned(const mch::Track& muon, const double z)
  {
    mch::TrackParam trackParam = mch::TrackParam(muon.first());
    return PropagateMCHParam(trackParam, z);
  }

  template <typename T>
  o2::dataformats::GlobalFwdTrack PropagateMCH(const T& muon, const double z)
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

    return PropagateMCH(track, z);
  }

  o2::dataformats::GlobalFwdTrack PropagateMCHToVertex(const o2::track::TrackParCovFwd& muon,
                                                       const double vx, const double vy, const double vz,
                                                       const double covVx, const double covVy)
  {
    auto mchTrack = FwdtoMCH(muon);

    o2::mch::TrackExtrap::extrapToVertex(mchTrack, vx, vy, vz, covVx, covVy);

    auto proptrack = MCHtoFwd(mchTrack);
    o2::dataformats::GlobalFwdTrack propmuon;
    propmuon.setParameters(proptrack.getParameters());
    propmuon.setZ(proptrack.getZ());
    propmuon.setCovariances(proptrack.getCovariances());

    return propmuon;
  }

  template <class C>
  o2::dataformats::GlobalFwdTrack PropagateMCHToVertex(const o2::track::TrackParCovFwd& muon,
                                                       const C& collision)
  {
    return PropagateMCHToVertex(muon,
                                collision.posX(),
                                collision.posY(),
                                collision.posZ(),
                                std::sqrt(collision.covXX()),
                                std::sqrt(collision.covYY()));
  }

  template <class TMFT>
  o2::dataformats::GlobalFwdTrack PropagateMFT(const TMFT& mftTrack, float z)
  {
    // static double Bz = -10001;
    double chi2 = mftTrack.chi2();
    SMatrix5 tpars = {mftTrack.x(), mftTrack.y(), mftTrack.phi(), mftTrack.tgl(), mftTrack.signed1Pt()};
    std::vector<double> v1{0, 0, 0, 0, 0,
                           0, 0, 0, 0, 0,
                           0, 0, 0, 0, 0};
    SMatrix55 tcovs(v1.begin(), v1.end());
    o2::track::TrackParCovFwd fwdtrack{mftTrack.z(), tpars, tcovs, chi2};
    o2::dataformats::GlobalFwdTrack propmuon;

    // double propVec[3] = {};
    // propVec[0] = collision.posX() - mftTrack.x();
    // propVec[1] = collision.posY() - mftTrack.y();
    // propVec[2] = collision.posZ() - mftTrack.z();

    // double centerZ[3] = {mftTrack.x() + propVec[0] / 2.,
    //                      mftTrack.y() + propVec[1] / 2.,
    //                      mftTrack.z() + propVec[2] / 2.};
    // if (Bz < -10000) {
    //  double centerZ[3] = {0, 0, (-45.f - 77.5f) / 2.f};
    //  o2::field::MagneticField* field = static_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField());
    //  Bz = field->getBz(centerZ);
    //}
    fwdtrack.propagateToZ(z, mBzAtMftCenter);

    propmuon.setParameters(fwdtrack.getParameters());
    propmuon.setZ(fwdtrack.getZ());
    propmuon.setCovariances(fwdtrack.getCovariances());

    return propmuon;
  }

  template <class C>
  o2::dataformats::GlobalFwdTrack PropagateMFTToDCA(o2::track::TrackParCovFwd mftTrackPar,
                                                    const C& collision,
                                                    float zshift)
  {
    // static double Bz = -10001;

    // double propVec[3] = {};
    // propVec[0] = collision.posX() - mftTrack.x();
    // propVec[1] = collision.posY() - mftTrack.y();
    // propVec[2] = collision.posZ() - mftTrack.z();

    // double centerZ[3] = {mftTrack.x() + propVec[0] / 2.,
    //                      mftTrack.y() + propVec[1] / 2.,
    //                      mftTrack.z() + propVec[2] / 2.};
    // if (Bz < -10000) {
    //  double centerZ[3] = {0, 0, -45.f / 2.f};
    //  o2::field::MagneticField* field = static_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField());
    //  Bz = field->getBz(centerZ);
    // }
    mftTrackPar.propagateToZ(collision.posZ() - zshift, mBzAtMftCenter);

    o2::dataformats::GlobalFwdTrack result;
    result.setParameters(mftTrackPar.getParameters());
    result.setZ(mftTrackPar.getZ());
    result.setCovariances(mftTrackPar.getCovariances());

    return result;
  }

  template <class C>
  o2::dataformats::GlobalFwdTrack PropagateMFTToDCA(o2::track::TrackParCovFwd mftTrackPar,
                                                    const o2::track::TrackParCovFwd& mchTrackPar,
                                                    const C& collision,
                                                    float zshift)
  {
    // static double Bz = -10001;

    // extrapolation with MCH tools
    auto mchTrackAtMFT = FwdtoMCH(mchTrackPar);
    o2::mch::TrackExtrap::extrapToVertex(mchTrackAtMFT,
                                         mftTrackPar.getX(),
                                         mftTrackPar.getY(),
                                         mftTrackPar.getZ(),
                                         std::sqrt(mftTrackPar.getSigma2X()),
                                         std::sqrt(mftTrackPar.getSigma2Y()));
    UpdateTrackMomentum(mftTrackPar, mchTrackAtMFT);

    // double propVec[3] = {};
    // propVec[0] = collision.posX() - mftTrack.x();
    // propVec[1] = collision.posY() - mftTrack.y();
    // propVec[2] = collision.posZ() - mftTrack.z();

    // double centerZ[3] = {mftTrack.x() + propVec[0] / 2.,
    //                      mftTrack.y() + propVec[1] / 2.,
    //                      mftTrack.z() + propVec[2] / 2.};
    // if (Bz < -10000) {
    //  double centerZ[3] = {0, 0, -45.f / 2.f};
    //  o2::field::MagneticField* field = static_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField());
    //  Bz = field->getBz(centerZ);
    // }
    mftTrackPar.propagateToZ(collision.posZ() - zshift, mBzAtMftCenter);

    o2::dataformats::GlobalFwdTrack result;
    result.setParameters(mftTrackPar.getParameters());
    result.setZ(mftTrackPar.getZ());
    result.setCovariances(mftTrackPar.getCovariances());

    return result;
  }

  o2::dataformats::GlobalFwdTrack PropagateMFTtoMCH(const o2::track::TrackParCovFwd& mftTrackPar,
                                                    const o2::mch::TrackParam& mchTrackPar,
                                                    const double z)
  {
    // extrapolation with MCH tools
    auto mchTrackAtMFT = mchTrackPar;
    o2::mch::TrackExtrap::extrapToVertex(mchTrackAtMFT,
                                         mftTrackPar.getX(),
                                         mftTrackPar.getY(),
                                         mftTrackPar.getZ(),
                                         std::sqrt(mftTrackPar.getSigma2X()),
                                         std::sqrt(mftTrackPar.getSigma2Y()));

    auto mftTrackProp = FwdtoMCH(mftTrackPar);
    UpdateTrackMomentum(mftTrackProp, mchTrackAtMFT);
    if (z < AbsorberBackZ) {
      o2::mch::TrackExtrap::extrapToZ(mftTrackProp, BransonPlaneZ);
      UpdateTrackMomentum(mftTrackProp, mchTrackPar);
    }

    if (cfgDipoleZshift.value != 0) {
      // extrapolate to the back of the absorber, taking into account the dipole shift,
      // to avoid that the correction bring the track starting point back into the absorber
      if (cfgDipoleZshift.value < 0) {
        o2::mch::TrackExtrap::extrapToZ(mftTrackProp, AbsorberBackZ);
      } else if (cfgDipoleZshift.value > 0) {
        o2::mch::TrackExtrap::extrapToZ(mftTrackProp, AbsorberBackZ - cfgDipoleZshift.value);
      }
      // shift the track starting point
      mftTrackProp.setZ(mftTrackProp.getZ() + cfgDipoleZshift.value);
      // extrapolate to the final z, corrected for the dipole shift
      o2::mch::TrackExtrap::extrapToZ(mftTrackProp, z + cfgDipoleZshift.value);
      // remove the shift from the extrapolated track
      mftTrackProp.setZ(mftTrackProp.getZ() - cfgDipoleZshift.value);
    } else {
      o2::mch::TrackExtrap::extrapToZ(mftTrackProp, z);
    }

    return MCHtoFwd(mftTrackProp);
  }

  template <class C>
  o2::dataformats::GlobalFwdTrack PropagateMFTToVertex(const o2::track::TrackParCovFwd& mftTrackPar,
                                                       const o2::track::TrackParCovFwd& mchTrackPar,
                                                       const C& collision)
  {
    // extrapolation with MCH tools
    auto mchTrackAtMFT = FwdtoMCH(mchTrackPar);
    o2::mch::TrackExtrap::extrapToVertex(mchTrackAtMFT,
                                         mftTrackPar.getX(),
                                         mftTrackPar.getY(),
                                         mftTrackPar.getZ(),
                                         std::sqrt(mftTrackPar.getSigma2X()),
                                         std::sqrt(mftTrackPar.getSigma2Y()));

    auto fwdTrackRefit = fwdtrackutils::refitGlobalMuonCov(MCHtoFwd(mchTrackAtMFT), mftTrackPar);

    // apply vertex shift correction
    fwdTrackRefit.setZ(fwdTrackRefit.getZ() + cfgVertexZshift.value);

    auto fwdTrackProp = fwdtrackutils::propagateTrackParCovFwd(fwdTrackRefit,
                                                               0,
                                                               collision,
                                                               fwdtrackutils::propagationPoint::kToVertex,
                                                               cfgRefPlaneZMFT,
                                                               mBzAtMftCenter);

    return fwdTrackProp;
  }

  void getMuonPairs(const CollisionInfo& collisionInfo,
                    std::vector<MuonPair>& muonPairs)
  {
    // outer loop over muon tracks
    for (const auto& mchIndex1 : collisionInfo.mchTracks) {
      // inner loop over muon tracks
      for (const auto& mchIndex2 : collisionInfo.mchTracks) {
        // avoid double-counting of muon pairs
        if (mchIndex2 <= mchIndex1) {
          continue;
        }

        muonPairs.emplace_back(mchIndex1, mchIndex2);
      }
    }
  }

  ROOT::Math::PxPyPzMVector getMuMu4Momentum(const o2::dataformats::GlobalFwdTrack& track1, const o2::dataformats::GlobalFwdTrack& track2)
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

    return muon1 + muon2;
  }

  double getMuMuAngle(const o2::dataformats::GlobalFwdTrack& track1, const o2::dataformats::GlobalFwdTrack& track2)
  {
    ROOT::Math::XYZVector muon1{
      track1.getPx(),
      track1.getPy(),
      track1.getPz()};

    ROOT::Math::XYZVector muon2{
      track2.getPx(),
      track2.getPy(),
      track2.getPz()};

    return std::acos(muon1.Unit().Dot(muon2.Unit()));
  }

  double getMuMuInvariantMass(const o2::dataformats::GlobalFwdTrack& track1, const o2::dataformats::GlobalFwdTrack& track2)
  {
    return getMuMu4Momentum(track1, track2).M();
  }

  template <class TMUON, class TCLUS>
  bool MchRealignTrack(const TMUON& mchTrack, const TCLUS& clusters, TrackRealigned& convertedTrack, bool applyCorrections)
  {
    auto mchTrackPar = FwdtoMCH(TrackToGlobalFwd(mchTrack));

    // loop over attached clusters
    int clIndex = -1;
    auto clustersSliced = clusters.sliceBy(perMuon, mchTrack.globalIndex()); // Slice clusters by muon id
    for (auto const& cluster : clustersSliced) {
      clIndex += 1;

      int deId = cluster.deId();
      int chamber = deId / 100 - 1;
      if (chamber < 0 || chamber >= NMchChambers) {
        continue;
      }

      math_utils::Point3D<double> local;
      math_utils::Point3D<double> master;

      master.SetXYZ(cluster.x(), cluster.y(), cluster.z());

      if (configRealign.cfgEnableMCHRealign) {
        // Transformation from reference geometry frame to new geometry frame
        transformRef[cluster.deId()].MasterToLocal(master, local);
        transformNew[cluster.deId()].LocalToMaster(local, master);
      }

      // shift the clusters to correct the longitudinal shift of the dipole
      if (cfgDipoleZshift.value != 0) {
        master.SetZ(master.z() + cfgDipoleZshift.value);
      }

      if (applyCorrections) {
        auto correctionsIt = mMchAlignmentCorrections.find(cluster.deId());
        if (correctionsIt != mMchAlignmentCorrections.end()) {
          const auto& corrections = correctionsIt->second;
          master.SetX(master.x() + corrections.x);
          master.SetY(master.y() + corrections.y);
          master.SetZ(master.z() + corrections.z);
        }
      }

      // realigned MCH cluster
      auto clusterMCH = new mch::Cluster();
      clusterMCH->x = master.x();
      clusterMCH->y = master.y();
      clusterMCH->z = master.z();

      uint32_t ClUId = mch::Cluster::buildUniqueId(static_cast<int>(cluster.deId() / 100) - 1, cluster.deId(), clIndex);
      clusterMCH->uid = ClUId;
      clusterMCH->ex = cluster.isGoodX() ? 0.2 : 10.0;
      clusterMCH->ey = cluster.isGoodY() ? 0.2 : 10.0;

      // Add transformed cluster into temporary variable
      convertedTrack.createParamAtCluster(*clusterMCH);
    }

    bool removable{false};
    // Refit the re-aligned track
    if (convertedTrack.getNClusters() != 0) {
      removable = RemoveTrack(convertedTrack);
    } else {
      LOGF(fatal, "Muon track %d has no associated clusters.", mchTrack.globalIndex());
    }

    // subtract the longitudinal shift of the dipole from the track z
    if (cfgDipoleZshift.value != 0) {
      auto& trackParam = *(convertedTrack.begin());
      trackParam.setZ(trackParam.getZ() - cfgDipoleZshift.value);
    }

    return !removable;
  }

  template <class COLL, class BC, class TMUON>
  void InitCollisions(COLL const& collisions,
                      BC const& bcs,
                      TMUON const& muonTracks,
                      aod::FwdTrkCls const& clusters,
                      std::map<uint64_t, CollisionInfo>& collisionInfos)
  {
    mMchTrackPars.clear();
    mMchTrackParsNew.clear();

    // fill collision information for global muon tracks (MFT-MCH-MID matches)
    for (const auto& muonTrack : muonTracks) {
      if (!muonTrack.has_collision()) {
        continue;
      }

      auto collision = collisions.rawIteratorAt(muonTrack.collisionId());

      if (cfgRequireGoodRCT && !rctChecker(collision)) {
        continue;
      }

      uint64_t collisionIndex = collision.globalIndex();

      auto bc = bcs.rawIteratorAt(collision.bcId());

      auto& collisionInfo = collisionInfos[collisionIndex];
      collisionInfo.bc = bc.globalBC();
      collisionInfo.zVertex = collision.posZ();

      if (static_cast<int>(muonTrack.trackType()) > GlobalTrackTypeMax) {
        // standalone MCH or MCH-MID tracks
        uint64_t mchTrackIndex = muonTrack.globalIndex();
        collisionInfo.mchTracks.push_back(mchTrackIndex);

        // initialize the original MCH track parameters
        mMchTrackPars.try_emplace(mchTrackIndex, TrackParExt(fwdtrackutils::getTrackParCovFwd(muonTrack, muonTrack), muonTrack.nClusters()));

        // refit MCH track if requested
        if (configRealign.cfgEnableMCHRefit || configRealign.cfgEnableMCHRealign) {
          TrackRealigned convertedTrack;
          bool convertedTrackOk = MchRealignTrack(muonTrack, clusters, convertedTrack, !mMchAlignmentCorrections.empty());

          // Get the re-aligned track parameters: track param at the first cluster
          mch::TrackParam trackParam = mch::TrackParam(convertedTrack.first());

          auto mchTrackParIt = mMchTrackParsNew.try_emplace(mchTrackIndex, TrackParExt(MCHtoFwd(trackParam), convertedTrack.getNClusters()));
          if (mchTrackParIt.second) {
            // the insertion succeeded
            mchTrackParIt.first->second.setTrackChi2(trackParam.getTrackChi2() / convertedTrack.getNDF());
            if (!convertedTrackOk) {
              mchTrackParIt.first->second.setRemovable();
            }
          }
        } else {
          // initialize the new MCH track parameters with the original ones, without refitting
          mMchTrackParsNew.try_emplace(mchTrackIndex, TrackParExt(fwdtrackutils::getTrackParCovFwd(muonTrack, muonTrack), muonTrack.nClusters()));
        }
      } else {
        // global muon tracks (MFT-MCH or MFT-MCH-MID)
        uint64_t muonTrackIndex = muonTrack.globalIndex();
        auto const& mchTrack = muonTrack.template matchMCHTrack_as<TMUON>();
        uint64_t mchTrackIndex = mchTrack.globalIndex();

        // check if a vector of global muon candidates is already available for the current MCH index
        // if not, initialize a new one and add the current global muon track
        // bool globalMuonTrackFound = false;
        auto matchingCandidateIterator = collisionInfo.globalMuonTracks.find(mchTrackIndex);
        if (matchingCandidateIterator != collisionInfo.globalMuonTracks.end()) {
          matchingCandidateIterator->second.push_back(muonTrackIndex);
          // globalMuonTrackFound = true;
        } else {
          collisionInfo.globalMuonTracks[mchTrackIndex].push_back(muonTrackIndex);
        }
      }
    }

    // sort the vectors of matching candidates in ascending order based on the matching chi2 value
    auto compareChi2 = [&muonTracks](uint64_t trackIndex1, uint64_t trackIndex2) -> bool {
      auto const& track1 = muonTracks.rawIteratorAt(trackIndex1);
      auto const& track2 = muonTracks.rawIteratorAt(trackIndex2);

      return (track1.chi2MatchMCHMFT() < track2.chi2MatchMCHMFT());
    };

    for (auto& [collisionIndex, collisionInfo] : collisionInfos) {                  // o2-linter: disable=const-ref-in-for-loop (object is modified in loop)
      for (auto& [mchIndex, globalTracksVector] : collisionInfo.globalMuonTracks) { // o2-linter: disable=const-ref-in-for-loop (object is modified in loop)
        std::sort(globalTracksVector.begin(), globalTracksVector.end(), compareChi2);
      }
    }
  }

  void InitCollisions(MyEvents const& collisions,
                      MyBCs const& bcs,
                      MyMuonsWithCov const& muonTracks,
                      aod::FwdTrkCls const& clusters,
                      MyMFTs const& mftTracks,
                      std::map<uint64_t, CollisionInfo>& collisionInfos)
  {
    InitCollisions(collisions, bcs, muonTracks, clusters, collisionInfos);

    mMftTrackPars.clear();
    mMftTrackParsNew.clear();

    // fill collision information for MFT standalone tracks
    for (const auto& mftTrack : mftTracks) {
      if (!mftTrack.has_collision()) {
        continue;
      }

      auto collision = collisions.rawIteratorAt(mftTrack.collisionId());
      uint64_t collisionIndex = collision.globalIndex();

      auto bc = bcs.rawIteratorAt(collision.bcId());

      uint64_t mftTrackIndex = mftTrack.globalIndex();

      auto& collisionInfo = collisionInfos[collisionIndex];
      collisionInfo.bc = bc.globalBC();
      collisionInfo.zVertex = collision.posZ();

      collisionInfo.mftTracks.push_back(mftTrackIndex);

      // initialize the original MFT track parameters
      auto mftTrackFwd = TrackToParCovFwd(mftTrack);
      mMftTrackPars.try_emplace(mftTrackIndex, TrackParExt(mftTrackFwd, mftTrack.nClusters()));

      // initialize the corrected MFT track parameters, if requested
      if (configMFTAlignmentCorrections.cfgEnableMFTAlignmentCorrections) {
        TransformMFT(mftTrackFwd);
        mMftTrackParsNew.try_emplace(mftTrackIndex, TrackParExt(mftTrackFwd, mftTrack.nClusters()));
      } else {
        // initialize the new MFT track parameters with the original ones, without corrections
        mMftTrackParsNew.try_emplace(mftTrackIndex, TrackParExt(mftTrackFwd, mftTrack.nClusters()));
      }
    }
  }

  void FillMftPlots(MyEvents const& collisions,
                    MyBCs const& bcs,
                    MyMuonsWithCov const& muonTracks,
                    MyMFTs const& mftTracks,
                    const std::map<uint64_t, CollisionInfo>& collisionInfos)
  {
    // outer loop over collisions
    for (const auto& [collisionIndex, collisionInfo] : collisionInfos) {
      auto const& collision = collisions.rawIteratorAt(collisionIndex);
      const auto& bc = bcs.rawIteratorAt(collision.bcId());

      // remove TF/ROF borders and ambiguous collisions
      if (!bc.selection_bit(o2::aod::evsel::kNoTimeFrameBorder) ||
          !bc.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
        continue;
      }

      registry.get<TH2>(HIST("vertex_y_vs_x"))->Fill(collision.posX(), collision.posY());
      registry.get<TH1>(HIST("vertex_z"))->Fill(collision.posZ());

      if (cfgEnableVertexShiftAnalysis || cfgEnableMftDcaAnalysis) {
        registry.get<TH1>(HIST("DCA/MFT/nTracksMFT"))->Fill(collisionInfo.mftTracks.size());
      }

      if (cfgEnableVertexShiftAnalysis || cfgEnableMftDcaAnalysis) {
        // loop over MFT tracks
        auto mftTrackIds = collisionInfo.mftTracks;
        if (cfgMftTracksMultiplicityMax > 0 && mftTrackIds.size() > cfgMftTracksMultiplicityMax) {
          auto rng = std::default_random_engine{};
          std::shuffle(std::begin(mftTrackIds), std::end(mftTrackIds), rng);
          mftTrackIds.resize(cfgMftTracksMultiplicityMax);
        }

        for (const auto& mftIndex : mftTrackIds) {
          auto const& mftTrack = mftTracks.rawIteratorAt(mftIndex);

          if (mftTrack.isCA()) {
            continue;
          }

          bool isGoodMFT = IsGoodMFT(mftTrack, 999.f, 5);
          if (!isGoodMFT) {
            continue;
          }

          // get the pre-stored MFT track parameters after corrections
          // if MFT corrections are not enabled, the original MFT track parameters are retrieved
          const auto mftTrackParIt = mMftTrackParsNew.find(mftIndex);
          if (mftTrackParIt == mMftTrackParsNew.end()) {
            continue;
          }
          const auto& mftTrackPar = mftTrackParIt->second;

          auto mftTrackAtDCA = PropagateMFTToDCA(mftTrackPar, collision, cfgVertexZshift);
          double dcax = mftTrackAtDCA.getX() - collision.posX();
          double dcay = mftTrackAtDCA.getY() - collision.posY();
          double phi = mftTrack.phi() * o2::constants::math::Rad2Deg;
          int mftNclusters = mftTrack.nClusters();
          double chi2NDF = static_cast<double>(mftNclusters) * 2 - 5;

          const int nMftLayers = 10;
          std::array<bool, 10> firedLayers{false};
          for (int layer = 0; layer < nMftLayers; layer++) {
            if (((mftTrack.mftClusterSizesAndTrackFlags() >> (layer * 6)) & 0x3F) != 0) {
              firedLayers[layer] = true;
            } else {
              firedLayers[layer] = false;
            }
          }

          if (cfgEnableMftDcaAnalysis) {
            if (cfgEnableMftDcaExtraPlots) {
              for (int i = 0; i < nMftLayers; i++) {
                if (firedLayers[i]) {
                  registry.get<THnSparse>(HIST("DCA/MFT/trackChi2"))->Fill(mftTrack.chi2() / chi2NDF, mftTrack.x(), mftTrack.y(), mftNclusters, i);
                }
              }
            }

            if (mftTrack.chi2() <= cfgTrackChi2MftUp) {
              registry.get<TH2>(HIST("DCA/MFT/DCA_y_vs_x"))->Fill(dcax, dcay);
              registry.get<THnSparse>(HIST("DCA/MFT/DCA_x"))->Fill(dcax, collision.posZ(), mftTrack.x(), mftTrack.y(), mftNclusters);
              registry.get<THnSparse>(HIST("DCA/MFT/DCA_y"))->Fill(dcay, collision.posZ(), mftTrack.x(), mftTrack.y(), mftNclusters);

              if (cfgProduceMFTTable) {
                mftTable(collision.posX(), collision.posY(), collision.posZ(),
                         mftTrack.signed1Pt(), mftTrack.tgl(), mftTrack.phi(),
                         dcax, dcay, mftTrackAtDCA.getSigma2X(), mftTrackAtDCA.getSigma2Y(), mftTrackAtDCA.getSigmaXY(),
                         mftNclusters, mftTrack.chi2(),
                         mftTrack.x(), mftTrack.y(), mftTrack.z());
              }

              if (cfgEnableMftDcaExtraPlots) {
                static constexpr int nMftClustersMin = 6;
                if (mftNclusters >= nMftClustersMin) {
                  for (int i = 0; i < nMftLayers; i++) {
                    auto mftTrackAtLayer = PropagateMFT(mftTrack, o2::mft::constants::mft::LayerZCoordinate()[i]);
                    std::get<std::shared_ptr<TH2>>(mMftTrackEffDen[i])->Fill(mftTrackAtLayer.getX(), mftTrackAtLayer.getY());
                    if (firedLayers[i]) {
                      std::get<std::shared_ptr<TH2>>(mMftTrackEffNum[i])->Fill(mftTrackAtLayer.getX(), mftTrackAtLayer.getY());
                      registry.get<THnSparse>(HIST("DCA/MFT/layers"))->Fill(i, mftTrack.x(), mftTrack.y(), mftNclusters);
                    }
                  }
                }
                for (int i = 0; i < nMftLayers; i++) {
                  if (firedLayers[i]) {
                    registry.get<THnSparse>(HIST("DCA/MFT/trackMomentum"))->Fill(mftTrack.p(), mftTrack.x(), mftTrack.y(), mftNclusters, i);
                  }
                }
              }
            }
          }

          if (cfgEnableVertexShiftAnalysis) {
            static constexpr int nMftClustersMin = 6;
            if (mftTrack.chi2() <= cfgTrackChi2MftUp && std::fabs(collision.posZ()) < 1.f && mftNclusters >= nMftClustersMin) {
              static constexpr int nPoints = 21;
              const std::array<float, nPoints> zshift{// in millimeters
                                                      -5.0, -4.5, -4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0,
                                                      0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0};
              for (int zi = 0; zi < nPoints; zi++) {
                auto mftTrackAtDCAshifted = PropagateMFTToDCA(mftTrackPar, collision, zshift[zi] / 10.f);
                double dcaxShifted = mftTrackAtDCAshifted.getX() - collision.posX();
                double dcayShifted = mftTrackAtDCAshifted.getY() - collision.posY();
                registry.get<TH3>(HIST("DCA/MFT/DCA_x_vs_phi_vs_zshift"))->Fill(zshift[zi], phi, dcaxShifted);
                registry.get<TH3>(HIST("DCA/MFT/DCA_y_vs_phi_vs_zshift"))->Fill(zshift[zi], phi, dcayShifted);

                auto mftTrackAtDCAshiftedPar = FwdtoMCH(mftTrackAtDCAshifted);
                auto slopex = mftTrackAtDCAshiftedPar.getNonBendingSlope();
                auto slopey = mftTrackAtDCAshiftedPar.getBendingSlope();
                registry.get<TH3>(HIST("DCA/MFT/DCA_x_vs_slopex_vs_zshift"))->Fill(zshift[zi], slopex, dcaxShifted);
                registry.get<TH3>(HIST("DCA/MFT/DCA_x_vs_slopey_vs_zshift"))->Fill(zshift[zi], slopey, dcaxShifted);
                registry.get<TH3>(HIST("DCA/MFT/DCA_y_vs_slopex_vs_zshift"))->Fill(zshift[zi], slopex, dcayShifted);
                registry.get<TH3>(HIST("DCA/MFT/DCA_y_vs_slopey_vs_zshift"))->Fill(zshift[zi], slopey, dcayShifted);
              }
            }
          }
        }
      }

      if (cfgEnableGlobalFwdDcaAnalysis) {
        // loop over global muon tracks
        for (const auto& [muonIndex, globalTracksVector] : collisionInfo.globalMuonTracks) {
          auto const& muonTrack = muonTracks.rawIteratorAt(globalTracksVector[0]);
          const auto& mchTrack = muonTrack.template matchMCHTrack_as<MyMuonsWithCov>();
          const auto& mftTrack = muonTrack.template matchMFTTrack_as<MyMFTs>();

          auto mchIndex = mchTrack.globalIndex();
          auto mftIndex = mftTrack.globalIndex();

          if (muonTrack.chi2MatchMCHMFT() < cfgMftDcaMatchChi2Up.value) {
            continue;
          }

          if (globalTracksVector.size() > 1) {
            auto const& muonTrack2 = muonTracks.rawIteratorAt(globalTracksVector[1]);
            double dchi2 = muonTrack2.chi2MatchMCHMFT() - muonTrack.chi2MatchMCHMFT();

            if (dchi2 < cfgMftDcaMatchChi2Up.value) {
              continue;
            }
          }

          if (mftTrack.isCA()) {
            continue;
          }

          bool isGoodMFT = IsGoodMFT(mftTrack, cfgTrackChi2MftUp, 5);
          if (!isGoodMFT) {
            continue;
          }

          bool isGoodMuon = IsGoodMuon(mchTrack, collision, cfgTrackChi2MchUp, 0.f, cfgPtMchLow, {cfgEtaMftLow, cfgEtaMftUp}, {cfgRabsLow, cfgRabsUp}, fSigmaPdcaUp);
          if (!isGoodMuon) {
            continue;
          }

          // get the pre-stored MFT track parameters after corrections
          // if MFT corrections are not enabled, the original MFT track parameters are retrieved
          const auto mftTrackParIt = mMftTrackParsNew.find(mftIndex);
          if (mftTrackParIt == mMftTrackParsNew.end()) {
            continue;
          }
          const auto& mftTrackPar = mftTrackParIt->second;

          // get the pre-stored MCH track parameters
          const auto mchTrackParIt = mMchTrackPars.find(mchIndex);
          if (mchTrackParIt == mMchTrackPars.end()) {
            continue;
          }
          const auto& mchTrackPar = mchTrackParIt->second;

          auto mftTrackAtDCA = PropagateMFTToDCA(mftTrackPar, mchTrackPar, collision, cfgVertexZshift);
          double dcax = mftTrackAtDCA.getX() - collision.posX();
          double dcay = mftTrackAtDCA.getY() - collision.posY();
          int mftNclusters = mftTrack.nClusters();

          if (mftTrack.chi2() <= cfgTrackChi2MftUp) {
            registry.get<THnSparse>(HIST("DCA/GlobalFwd/DCA_x"))->Fill(dcax, collision.posZ(), mftTrack.x(), mftTrack.y(), mftNclusters);
            registry.get<THnSparse>(HIST("DCA/GlobalFwd/DCA_y"))->Fill(dcay, collision.posZ(), mftTrack.x(), mftTrack.y(), mftNclusters);
          }
        }
      }
    }
  }

  void FillMchPlots(MyEvents const& collisions,
                    MyBCs const& bcs,
                    MyMuonsWithCov const& muonTracks,
                    aod::FwdTrkCls const& clusters,
                    const std::map<uint64_t, CollisionInfo>& collisionInfos)
  {
    if (!cfgEnableMftMchResidualsAnalysis && !cfgEnableMftMchMatchingAnalysis) {
      return;
    }

    // loop over collisions
    for (const auto& [collisionIndex, collisionInfo] : collisionInfos) {
      auto const& collision = collisions.rawIteratorAt(collisionIndex);
      const auto& bc = bcs.rawIteratorAt(collision.bcId());

      // remove TF/ROF borders and ambiguous collisions
      if (!bc.selection_bit(o2::aod::evsel::kNoTimeFrameBorder) ||
          !bc.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
        continue;
      }

      // loop over global muon tracks
      for (const auto& [muonIndex, globalTracksVector] : collisionInfo.globalMuonTracks) {
        auto const& muonTrack = muonTracks.rawIteratorAt(globalTracksVector[0]);
        const auto& mchTrack = muonTrack.template matchMCHTrack_as<MyMuonsWithCov>();
        const auto& mftTrack = muonTrack.template matchMFTTrack_as<MyMFTs>();
        // int quadrant = GetQuadrant(mchTrack);
        int quadrant = GetQuadrant(mftTrack);
        int posNeg = (mchTrack.sign() >= 0) ? 0 : 1;

        auto mchIndex = mchTrack.globalIndex();
        auto mftIndex = mftTrack.globalIndex();

        bool isGoodMuon = IsGoodMuon(mchTrack, collision, cfgTrackChi2MchUp, cfgMftMchResidualsPLow, cfgMftMchResidualsPtLow, {cfgEtaMftLow, cfgEtaMftUp}, {cfgRabsLow, cfgRabsUp}, fSigmaPdcaUp);
        if (!isGoodMuon) {
          continue;
        }

        bool isGoodMFT = IsGoodMFT(mftTrack, cfgTrackChi2MftUp, cfgTrackNClustMftLow);
        if (!isGoodMFT) {
          continue;
        }

        // get the pre-stored MFT track parameters
        const auto mftTrackParIt = mMftTrackPars.find(mftIndex);
        if (mftTrackParIt == mMftTrackPars.end()) {
          continue;
        }
        const auto& mftTrackPar = mftTrackParIt->second;

        // get the pre-stored MFT track parameters after corrections
        // if MFT corrections are not enabled, the original MFT track parameters are retrieved
        const auto mftTrackParNewIt = mMftTrackParsNew.find(mftIndex);
        if (mftTrackParNewIt == mMftTrackParsNew.end()) {
          continue;
        }
        const auto& mftTrackParNew = mftTrackParNewIt->second;

        // get the pre-stored MCH track parameters
        const auto mchTrackParIt = mMchTrackPars.find(mchIndex);
        if (mchTrackParIt == mMchTrackPars.end()) {
          continue;
        }
        const auto& mchTrackPar = mchTrackParIt->second;

        // get the pre-stored MCH track parameters after refit
        const auto mchTrackParNewIt = mMchTrackParsNew.find(mchIndex);
        if (mchTrackParNewIt == mMchTrackParsNew.end()) {
          continue;
        }
        const auto& mchTrackParNew = mchTrackParNewIt->second;

        double matchChi2 = muonTrack.chi2MatchMCHMFT();

        // Residuals analysis between MFT tracks and MCH clusters
        if (cfgEnableMftMchResidualsAnalysis && (matchChi2 <= cfgMftMchResidualsMatchChi2Up.value)) {
          // loop over attached clusters
          auto clustersSliced = clusters.sliceBy(perMuon, mchTrack.globalIndex()); // Slice clusters by muon id
          for (auto const& cluster : clustersSliced) {
            int deId = cluster.deId();
            int chamber = deId / 100 - 1;
            if (chamber < 0 || chamber >= NMchChambers) {
              continue;
            }
            int deIndex = getDEindex(deId);

            math_utils::Point3D<double> local;
            math_utils::Point3D<double> master;        // original cluster position
            math_utils::Point3D<double> masterRealign; // cluster position after realignment

            master.SetXYZ(cluster.x(), cluster.y(), cluster.z());
            masterRealign.SetXYZ(cluster.x(), cluster.y(), cluster.z());

            // apply realignment to MCH cluster
            if (configRealign.cfgEnableMCHRealign) {
              // Transformation from reference geometry frame to new geometry frame
              transformRef[cluster.deId()].MasterToLocal(master, local);
              transformNew[cluster.deId()].LocalToMaster(local, masterRealign);
            }

            // apply alignment corrections to MCH cluster (if available)
            if (!mMchAlignmentCorrections.empty()) {
              auto correctionsIt = mMchAlignmentCorrections.find(cluster.deId());
              if (correctionsIt != mMchAlignmentCorrections.end()) {
                const auto& corrections = correctionsIt->second;
                masterRealign.SetX(masterRealign.x() + corrections.x);
                masterRealign.SetY(masterRealign.y() + corrections.y);
                masterRealign.SetZ(masterRealign.z() + corrections.z);
              }
            }

            // MFT-MCH residuals from original alignment
            const auto mftTrackAtCluster = PropagateMFTtoMCH(mftTrackPar, FwdtoMCH(mchTrackPar), master.z());
            const auto mftTrackParamAtCluster = FwdtoMCH(mftTrackAtCluster);

            const std::array<double, 2> xPos{master.x(), mftTrackAtCluster.getX()};
            const std::array<double, 2> yPos{master.y(), mftTrackAtCluster.getY()};

            registry.get<THnSparse>(HIST("residuals/dx_vs_chamber"))->Fill(chamber + 1, quadrant, posNeg, xPos[0] - xPos[1]);
            registry.get<THnSparse>(HIST("residuals/dy_vs_chamber"))->Fill(chamber + 1, quadrant, posNeg, yPos[0] - yPos[1]);

            registry.get<THnSparse>(HIST("residuals/dx_vs_de"))->Fill(xPos[0] - xPos[1], deIndex, quadrant, posNeg, mchTrack.p(), mftTrackParamAtCluster.getNonBendingSlope());
            registry.get<THnSparse>(HIST("residuals/dy_vs_de"))->Fill(yPos[0] - yPos[1], deIndex, quadrant, posNeg, mchTrack.p(), mftTrackParamAtCluster.getBendingSlope());

            // MFT-MCH residuals with realigned and/or corrected MCH clusters
            // if the alignment corrections are available and the refitting is successful, the MFT track is extrpolated
            // by taking the momentum from the MCH track refitted with the alignment corrections and the new
            // alignment (if realignment is enabled)
            if (!mchTrackParNew.isRemovable()) {
              const auto mftTrackAtClusterWithCorr = PropagateMFTtoMCH(mftTrackParNew, FwdtoMCH(mchTrackParNew), masterRealign.z());
              const auto mftTrackParamAtClusterWithCorr = FwdtoMCH(mftTrackAtClusterWithCorr);

              const std::array<double, 2> xPosWithCorr{masterRealign.x(), mftTrackAtClusterWithCorr.getX()};
              const std::array<double, 2> yPosWithCorr{masterRealign.y(), mftTrackAtClusterWithCorr.getY()};

              registry.get<THnSparse>(HIST("residuals/dx_vs_chamber_corr"))->Fill(chamber + 1, quadrant, posNeg, xPosWithCorr[0] - xPosWithCorr[1]);
              registry.get<THnSparse>(HIST("residuals/dy_vs_chamber_corr"))->Fill(chamber + 1, quadrant, posNeg, yPosWithCorr[0] - yPosWithCorr[1]);

              registry.get<THnSparse>(HIST("residuals/dx_vs_de_corr"))->Fill(xPosWithCorr[0] - xPosWithCorr[1], deIndex, quadrant, posNeg, mchTrack.p(), mftTrackParamAtClusterWithCorr.getNonBendingSlope());
              registry.get<THnSparse>(HIST("residuals/dy_vs_de_corr"))->Fill(yPosWithCorr[0] - yPosWithCorr[1], deIndex, quadrant, posNeg, mchTrack.p(), mftTrackParamAtClusterWithCorr.getBendingSlope());
            }
          }

          const auto mchTrackAtDCA = PropagateMCHParam(FwdtoMCH(mchTrackPar), collision.posZ());
          const auto dcax = mchTrackAtDCA.getX() - collision.posX();
          const auto dcay = mchTrackAtDCA.getY() - collision.posY();

          registry.get<TH2>(HIST("DCA/MCH/DCA_y_vs_x"))->Fill(dcax, dcay);
          registry.get<THnSparse>(HIST("DCA/MCH/DCA_x_vs_sign_vs_quadrant_vs_mom"))->Fill(mchTrackPar.getP(), quadrant, posNeg, dcax);
          registry.get<THnSparse>(HIST("DCA/MCH/DCA_y_vs_sign_vs_quadrant_vs_mom"))->Fill(mchTrackPar.getP(), quadrant, posNeg, dcay);

          if (cfgEnableMftMchResidualsExtraPlots) {
            registry.get<THnSparse>(HIST("DCA/MCH/DCA_x_vs_sign_vs_quadrant_vs_vz"))->Fill(collision.posZ(), quadrant, posNeg, dcax);
            registry.get<THnSparse>(HIST("DCA/MCH/DCA_y_vs_sign_vs_quadrant_vs_vz"))->Fill(collision.posZ(), quadrant, posNeg, dcay);
            const auto mchTrackAtMFT = PropagateMCHParam(FwdtoMCH(mchTrackPar), mftTrack.z());
            const double deltaPhi = mchTrackAtMFT.getPhi() - mftTrack.phi();
            registry.get<THnSparse>(HIST("residuals/dphi_at_mft"))->Fill(deltaPhi, mftTrack.x(), mftTrack.y(), posNeg, mchTrackAtMFT.getP());
          }

          if (!mchTrackParNew.isRemovable()) {
            const auto mchTrackAtDCAWithCorr = PropagateMCHParam(FwdtoMCH(mchTrackParNew), collision.posZ());
            const auto dcaxWithCorr = mchTrackAtDCAWithCorr.getX() - collision.posX();
            const auto dcayWithCorr = mchTrackAtDCAWithCorr.getY() - collision.posY();

            registry.get<THnSparse>(HIST("DCA/MCH/DCA_x_vs_sign_vs_quadrant_vs_mom_corr"))->Fill(mchTrackParNew.getP(), quadrant, posNeg, dcaxWithCorr);
            registry.get<THnSparse>(HIST("DCA/MCH/DCA_y_vs_sign_vs_quadrant_vs_mom_corr"))->Fill(mchTrackParNew.getP(), quadrant, posNeg, dcayWithCorr);
          }
        }

        // MFT-MCH track residuals analysis
        if (cfgEnableMftMchMatchingAnalysis && !mchTrackParNew.isRemovable() && (matchChi2 <= cfgMftMchResidualsMatchChi2Up.value)) {
          static constexpr int nRefPlanes = 2;
          const std::array<double, nRefPlanes> refPlaneZ{cfgRefPlaneZMFT, cfgRefPlaneZMCH};

          std::array<std::shared_ptr<THnSparse>, 2> dxPlots{registry.get<THnSparse>(HIST("matching/dxAtMFT")), registry.get<THnSparse>(HIST("matching/dxAtMCH"))};
          std::array<std::shared_ptr<THnSparse>, 2> dyPlots{registry.get<THnSparse>(HIST("matching/dyAtMFT")), registry.get<THnSparse>(HIST("matching/dyAtMCH"))};
          std::array<std::shared_ptr<THnSparse>, 2> dsxPlots{registry.get<THnSparse>(HIST("matching/dsxAtMFT")), registry.get<THnSparse>(HIST("matching/dsxAtMCH"))};
          std::array<std::shared_ptr<THnSparse>, 2> dsyPlots{registry.get<THnSparse>(HIST("matching/dsyAtMFT")), registry.get<THnSparse>(HIST("matching/dsyAtMCH"))};
          std::array<std::shared_ptr<THnSparse>, 2> dphiPlots{registry.get<THnSparse>(HIST("matching/dphiAtMFT")), registry.get<THnSparse>(HIST("matching/dphiAtMCH"))};

          for (int iRefPlane = 0; iRefPlane < nRefPlanes; iRefPlane++) {
            const auto mftTrackAtRefPlane = PropagateMFTtoMCH(mftTrackParNew, FwdtoMCH(mchTrackParNew), refPlaneZ[iRefPlane]);
            const auto mchTrackAtRefPlane = PropagateMCHParam(FwdtoMCH(mchTrackParNew), refPlaneZ[iRefPlane]);
            const auto& refTrackAtRefPlane = (iRefPlane == 0) ? mftTrackAtRefPlane : mchTrackAtRefPlane;

            auto dx = mchTrackAtRefPlane.getX() - mftTrackAtRefPlane.getX();
            dxPlots[iRefPlane]->Fill(dx, refTrackAtRefPlane.getX(), refTrackAtRefPlane.getY(), quadrant, posNeg, mchTrackParNew.getP());
            auto dy = mchTrackAtRefPlane.getY() - mftTrackAtRefPlane.getY();
            dyPlots[iRefPlane]->Fill(dy, refTrackAtRefPlane.getX(), refTrackAtRefPlane.getY(), quadrant, posNeg, mchTrackParNew.getP());

            const auto mftParamAtRefPlane = FwdtoMCH(mftTrackAtRefPlane);
            const auto mchParamAtRefPlane = FwdtoMCH(mchTrackAtRefPlane);

            auto dsx = mchParamAtRefPlane.getNonBendingSlope() - mftParamAtRefPlane.getNonBendingSlope();
            dsxPlots[iRefPlane]->Fill(dsx, refTrackAtRefPlane.getX(), refTrackAtRefPlane.getY(), quadrant, posNeg, mchTrackParNew.getP());
            auto dsy = mchParamAtRefPlane.getBendingSlope() - mftParamAtRefPlane.getBendingSlope();
            dsyPlots[iRefPlane]->Fill(dsy, refTrackAtRefPlane.getX(), refTrackAtRefPlane.getY(), quadrant, posNeg, mchTrackParNew.getP());

            auto dphi = RecoDecay::constrainAngle(mchTrackAtRefPlane.getPhi() - mftTrackAtRefPlane.getPhi(), -o2::constants::math::PI);
            dphiPlots[iRefPlane]->Fill(dphi, refTrackAtRefPlane.getX(), refTrackAtRefPlane.getY(), quadrant, posNeg, mchTrackParNew.getP());
          }
        }
      }
    }
  }

  template <typename T1, typename T2, typename T3, typename T4, typename HistConfigType>
  inline void fillDimuonInvmassPlot(const T1& trackPar1,
                                    const T2& trackPar2,
                                    const T3& trackPar1AtVertex,
                                    const T4& trackPar2AtVertex,
                                    HistConfigType histConfig)
  {
    auto mumu4mom = getMuMu4Momentum(trackPar1AtVertex, trackPar2AtVertex);
    double p = mumu4mom.P();
    double pT = mumu4mom.Pt();
    double mass = mumu4mom.M();
    int quadrant1 = GetQuadrant(static_cast<float>(std::atan2(trackPar1.getY(), trackPar1.getX())));
    int quadrant2 = GetQuadrant(static_cast<float>(std::atan2(trackPar2.getY(), trackPar2.getX())));

    registry.get<THnSparse>(histConfig)->Fill(mass, p, pT, quadrant1, quadrant2);
  }

  template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename HistConfigType1, typename HistConfigType2>
  inline void fillDimuonDcaPlots(const T1& trackPar1,
                                 const T2& trackPar2,
                                 const T3& trackPar1AtVertex,
                                 const T4& trackPar2AtVertex,
                                 const T5& trackPar1AtDca,
                                 const T6& trackPar2AtDca,
                                 HistConfigType1 histConfigX,
                                 HistConfigType2 histConfigY)
  {
    auto mumu4mom = getMuMu4Momentum(trackPar1AtVertex, trackPar2AtVertex);
    double p = mumu4mom.P();
    double pT = mumu4mom.Pt();
    double dcax = trackPar1AtDca.getX() - trackPar2AtDca.getX();
    double dcay = trackPar1AtDca.getY() - trackPar2AtDca.getY();
    int quadrant1 = GetQuadrant(static_cast<float>(std::atan2(trackPar1.getY(), trackPar1.getX())));
    int quadrant2 = GetQuadrant(static_cast<float>(std::atan2(trackPar2.getY(), trackPar2.getX())));

    registry.get<THnSparse>(histConfigX)->Fill(dcax, p, pT, quadrant1, quadrant2);
    registry.get<THnSparse>(histConfigY)->Fill(dcay, p, pT, quadrant1, quadrant2);
  }

  template <typename T1, typename T2, typename T3, typename T4, typename HistConfigType>
  inline void fillDimuonAnglePlot(const T1& trackPar1,
                                  const T2& trackPar2,
                                  const T3& trackPar1AtVertex,
                                  const T4& trackPar2AtVertex,
                                  double mchAngle,
                                  double fwdAngle,
                                  HistConfigType histConfig)
  {
    auto mumu4mom = getMuMu4Momentum(trackPar1AtVertex, trackPar2AtVertex);
    double p = mumu4mom.P();
    double pT = mumu4mom.Pt();
    double dAngle = mchAngle - fwdAngle;
    int quadrant1 = GetQuadrant(static_cast<float>(std::atan2(trackPar1.getY(), trackPar1.getX())));
    int quadrant2 = GetQuadrant(static_cast<float>(std::atan2(trackPar2.getY(), trackPar2.getX())));

    registry.get<THnSparse>(histConfig)->Fill(dAngle, fwdAngle, p, pT, quadrant1, quadrant2);
  }

  void FillDimuonPlots(MyEvents const& collisions,
                       MyMuonsWithCov const& muonTracks,
                       const std::map<uint64_t, CollisionInfo>& collisionInfos)
  {
    if (!cfgEnableDimuonAnalysis) {
      return;
    }

    for (const auto& [collisionIndex, collisionInfo] : collisionInfos) {
      auto const& collision = collisions.rawIteratorAt(collisionIndex);

      std::vector<MuonPair> muonPairs;
      getMuonPairs(collisionInfo, muonPairs);

      for (const auto& [mchIndex1, mchIndex2] : muonPairs) {

        auto const& muonTrack1 = muonTracks.rawIteratorAt(mchIndex1);
        auto const& muonTrack2 = muonTracks.rawIteratorAt(mchIndex2);
        int sign1 = muonTrack1.sign();
        int sign2 = muonTrack2.sign();

        // only consider opposite-sign pairs
        if ((sign1 * sign2) >= 0) {
          continue;
        }

        bool isGoodMuon1 = IsGoodMuon(muonTrack1, collision, cfgTrackChi2MchUp, 0.f, cfgPtMchLow, {cfgEtaMchLow, cfgEtaMchUp}, {cfgRabsLow, cfgRabsUp}, fSigmaPdcaUp);
        bool isGoodMuon2 = IsGoodMuon(muonTrack2, collision, cfgTrackChi2MchUp, 0.f, cfgPtMchLow, {cfgEtaMchLow, cfgEtaMchUp}, {cfgRabsLow, cfgRabsUp}, fSigmaPdcaUp);
        bool goodMuonTracks = (isGoodMuon1 && isGoodMuon2);

        if (!goodMuonTracks) {
          continue;
        }

        // get the pre-stored MCH track parameters
        const auto mchTrackParIt1 = mMchTrackPars.find(mchIndex1);
        if (mchTrackParIt1 == mMchTrackPars.end()) {
          continue;
        }
        const auto mchTrackParIt2 = mMchTrackPars.find(mchIndex2);
        if (mchTrackParIt2 == mMchTrackPars.end()) {
          continue;
        }
        const auto& mchTrackPar1 = mchTrackParIt1->second;
        auto mchTrackPar1AtDca = PropagateMCHParam(FwdtoMCH(mchTrackPar1), collision.posZ());
        auto mchTrackPar1AtVertex = PropagateMCHToVertex(mchTrackPar1, collision);
        const auto& mchTrackPar2 = mchTrackParIt2->second;
        auto mchTrackPar2AtDca = PropagateMCHParam(FwdtoMCH(mchTrackPar2), collision.posZ());
        auto mchTrackPar2AtVertex = PropagateMCHToVertex(mchTrackPar2, collision);

        // get the pre-stored MCH track parameters after refit
        const auto mchTrackParNewIt1 = mMchTrackParsNew.find(mchIndex1);
        if (mchTrackParNewIt1 == mMchTrackParsNew.end()) {
          continue;
        }
        const auto mchTrackParNewIt2 = mMchTrackParsNew.find(mchIndex2);
        if (mchTrackParNewIt2 == mMchTrackParsNew.end()) {
          continue;
        }
        const auto& mchTrackParNew1 = mchTrackParNewIt1->second;
        auto mchTrackParNew1AtDca = PropagateMCHParam(FwdtoMCH(mchTrackParNew1), collision.posZ());
        auto mchTrackParNew1AtVertex = PropagateMCHToVertex(mchTrackParNew1, collision);
        const auto& mchTrackParNew2 = mchTrackParNewIt2->second;
        auto mchTrackParNew2AtDca = PropagateMCHParam(FwdtoMCH(mchTrackParNew2), collision.posZ());
        auto mchTrackParNew2AtVertex = PropagateMCHToVertex(mchTrackParNew2, collision);

        fillDimuonInvmassPlot(mchTrackPar1, mchTrackPar2, mchTrackPar1AtVertex, mchTrackPar2AtVertex, HIST("dimuon/invariantMass_MuonKine_MuonCuts"));
        fillDimuonInvmassPlot(mchTrackParNew1, mchTrackParNew2, mchTrackParNew1AtVertex, mchTrackParNew2AtVertex, HIST("dimuon/realign/invariantMass_MuonKine_MuonCuts"));

        fillDimuonDcaPlots(mchTrackPar1, mchTrackPar2,
                           mchTrackPar1AtVertex, mchTrackPar2AtVertex,
                           mchTrackPar1AtDca, mchTrackPar2AtDca,
                           HIST("dimuon/dcax_MuonKine_MuonCuts"), HIST("dimuon/dcay_MuonKine_MuonCuts"));
        fillDimuonDcaPlots(mchTrackParNew1, mchTrackParNew2,
                           mchTrackParNew1AtVertex, mchTrackParNew2AtVertex,
                           mchTrackParNew1AtDca, mchTrackParNew2AtDca,
                           HIST("dimuon/realign/dcax_MuonKine_MuonCuts"), HIST("dimuon/realign/dcay_MuonKine_MuonCuts"));

        double mchAngle = getMuMuAngle(mchTrackPar1AtVertex, mchTrackPar2AtVertex);
        double mchAngleNew = getMuMuAngle(mchTrackParNew1AtVertex, mchTrackParNew2AtVertex);

        try {
          const auto& candidates1 = collisionInfo.globalMuonTracks.at(mchIndex1);
          const auto& candidates2 = collisionInfo.globalMuonTracks.at(mchIndex2);

          auto fwdIndex1 = candidates1[0];
          auto fwdIndex2 = candidates2[0];

          auto const& fwdTrack1 = muonTracks.rawIteratorAt(fwdIndex1);
          auto const& fwdTrack2 = muonTracks.rawIteratorAt(fwdIndex2);

          bool isGoodGlobalMuon1 = IsGoodMuon(muonTrack1, collision, cfgTrackChi2MchUp, 0.f, cfgPtMchLow, {cfgEtaMftLow, cfgEtaMftUp}, {cfgRabsLow, cfgRabsUp}, fSigmaPdcaUp);
          bool isGoodGlobalMuon2 = IsGoodMuon(muonTrack2, collision, cfgTrackChi2MchUp, 0.f, cfgPtMchLow, {cfgEtaMftLow, cfgEtaMftUp}, {cfgRabsLow, cfgRabsUp}, fSigmaPdcaUp);
          bool goodGlobalMuonTracks = (isGoodGlobalMuon1 && isGoodGlobalMuon2);

          bool isGoodMatch1 = isGoodGlobalMatching(fwdTrack1, cfgDimuonMatchChi2Up.value);
          bool isGoodMatch2 = isGoodGlobalMatching(fwdTrack2, cfgDimuonMatchChi2Up.value);
          bool goodGlobalMuonMatches = (isGoodMatch1 && isGoodMatch2);

          if (!goodGlobalMuonTracks || !goodGlobalMuonMatches) {
            continue;
          }

          fillDimuonInvmassPlot(mchTrackPar1, mchTrackPar2, mchTrackPar1AtVertex, mchTrackPar2AtVertex, HIST("dimuon/invariantMass_MuonKine_GlobalMuonCuts_GoodMatches"));
          fillDimuonInvmassPlot(mchTrackParNew1, mchTrackParNew2, mchTrackParNew1AtVertex, mchTrackParNew2AtVertex, HIST("dimuon/realign/invariantMass_MuonKine_GlobalMuonCuts_GoodMatches"));

          fillDimuonDcaPlots(mchTrackPar1, mchTrackPar2,
                             mchTrackPar1AtVertex, mchTrackPar2AtVertex,
                             mchTrackPar1AtDca, mchTrackPar2AtDca,
                             HIST("dimuon/dcax_MuonKine_GlobalMuonCuts_GoodMatches"), HIST("dimuon/dcay_MuonKine_GlobalMuonCuts_GoodMatches"));
          fillDimuonDcaPlots(mchTrackParNew1, mchTrackParNew2,
                             mchTrackParNew1AtVertex, mchTrackParNew2AtVertex,
                             mchTrackParNew1AtDca, mchTrackParNew2AtDca,
                             HIST("dimuon/realign/dcax_MuonKine_GlobalMuonCuts_GoodMatches"), HIST("dimuon/realign/dcay_MuonKine_GlobalMuonCuts_GoodMatches"));

          auto mftIndex1 = fwdTrack1.matchMFTTrackId();
          auto mftIndex2 = fwdTrack2.matchMFTTrackId();

          const auto mftTrackPar1 = mMftTrackPars.at(mftIndex1);
          auto fwdTrackPar1AtDca = PropagateMFTToDCA(mftTrackPar1, mchTrackPar1, collision, cfgVertexZshift);
          auto fwdTrackPar1AtVertex = PropagateMFTToVertex(mftTrackPar1, mchTrackPar1, collision);
          const auto mftTrackPar2 = mMftTrackPars.at(mftIndex2);
          auto fwdTrackPar2AtDca = PropagateMFTToDCA(mftTrackPar2, mchTrackPar2, collision, cfgVertexZshift);
          auto fwdTrackPar2AtVertex = PropagateMFTToVertex(mftTrackPar2, mchTrackPar2, collision);
          const auto mftTrackParNew1 = mMftTrackParsNew.at(mftIndex1);
          auto fwdTrackParNew1AtDca = PropagateMFTToDCA(mftTrackParNew1, mchTrackParNew1, collision, cfgVertexZshift);
          auto fwdTrackParNew1AtVertex = PropagateMFTToVertex(mftTrackParNew1, mchTrackParNew1, collision);
          const auto mftTrackParNew2 = mMftTrackParsNew.at(mftIndex2);
          auto fwdTrackParNew2AtDca = PropagateMFTToDCA(mftTrackParNew2, mchTrackParNew2, collision, cfgVertexZshift);
          auto fwdTrackParNew2AtVertex = PropagateMFTToVertex(mftTrackParNew2, mchTrackParNew2, collision);

          fillDimuonInvmassPlot(mchTrackPar1, mchTrackPar2, fwdTrackPar1AtVertex, fwdTrackPar2AtVertex, HIST("dimuon/invariantMass_ScaledMftKine_GlobalMuonCuts_GoodMatches"));
          fillDimuonInvmassPlot(mchTrackParNew1, mchTrackParNew2, fwdTrackParNew1AtVertex, fwdTrackParNew2AtVertex, HIST("dimuon/realign/invariantMass_ScaledMftKine_GlobalMuonCuts_GoodMatches"));

          fillDimuonDcaPlots(mchTrackPar1, mchTrackPar2,
                             fwdTrackPar1AtVertex, fwdTrackPar2AtVertex,
                             fwdTrackPar1AtDca, fwdTrackPar2AtDca,
                             HIST("dimuon/dcax_ScaledMftKine_GlobalMuonCuts_GoodMatches"), HIST("dimuon/dcay_ScaledMftKine_GlobalMuonCuts_GoodMatches"));
          fillDimuonDcaPlots(mchTrackParNew1, mchTrackParNew2,
                             fwdTrackParNew1AtVertex, fwdTrackParNew2AtVertex,
                             fwdTrackParNew1AtDca, fwdTrackParNew2AtDca,
                             HIST("dimuon/realign/dcax_ScaledMftKine_GlobalMuonCuts_GoodMatches"), HIST("dimuon/realign/dcay_ScaledMftKine_GlobalMuonCuts_GoodMatches"));

          double fwdAngle = getMuMuAngle(fwdTrackPar1AtVertex, fwdTrackPar2AtVertex);
          double fwdAngleNew = getMuMuAngle(fwdTrackParNew1AtVertex, fwdTrackParNew2AtVertex);

          fillDimuonAnglePlot(mchTrackPar1, mchTrackPar2, fwdTrackPar1AtVertex, fwdTrackPar2AtVertex, mchAngle, fwdAngle, HIST("dimuon/angle_GlobalMuonCuts_GoodMatches"));
          fillDimuonAnglePlot(mchTrackParNew1, mchTrackParNew2, fwdTrackParNew1AtVertex, fwdTrackParNew2AtVertex, mchAngleNew, fwdAngleNew, HIST("dimuon/realign/angle_GlobalMuonCuts_GoodMatches"));
        } catch (const std::exception&) {
          continue;
        }
      }
    }
  }

  void processQA(MyEvents const& collisions,
                 MyBCs const& bcs,
                 MyMuonsWithCov const& muonTracks,
                 MyMFTs const& mftTracks,
                 // MyMFTCovariances const& mftCovariances,
                 aod::FwdTrkCls const& clusters)
  {
    auto bc = bcs.begin();
    if (mRunNumber != bc.runNumber()) {
      initCCDB(bc);
      LOGF(info, "Set field for muons");
      VarManager::SetupMuonMagField();
      mRunNumber = bc.runNumber();
    }

    std::map<uint64_t, CollisionInfo> collisionInfos;
    InitCollisions(collisions, bcs, muonTracks, clusters, mftTracks, collisionInfos);

    FillMftPlots(collisions, bcs, muonTracks, mftTracks, collisionInfos);

    FillMchPlots(collisions, bcs, muonTracks, clusters, collisionInfos);

    FillDimuonPlots(collisions, muonTracks, collisionInfos);
  }

  PROCESS_SWITCH(muonGlobalAlignment, processQA, "processQA", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<muonGlobalAlignment>(cfgc)};
};
