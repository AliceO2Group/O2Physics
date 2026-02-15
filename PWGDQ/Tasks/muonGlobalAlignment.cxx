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
/// \file muonDCA.cxx
/// \brief Task to compute and evaluate DCA quantities
/// \author Nicolas Bizé <nicolas.bize@cern.ch>, SUBATECH
//
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
#include <MCHBase/TrackerParam.h>
#include <MCHGeometryTransformer/Transformations.h>
#include <MCHTracking/Track.h>
#include <MCHTracking/TrackExtrap.h>
#include <MCHTracking/TrackFitter.h>
#include <MCHTracking/TrackParam.h>

#include <algorithm>
#include <limits>
#include <map>
#include <random>
#include <string>
#include <unordered_map>
#include <vector>

using namespace o2;
using namespace o2::mch;
using namespace o2::framework;
using namespace o2::aod;

using namespace o2::aod::rctsel;

using MyReducedMuons = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra, aod::ReducedMuonsCov>;
using MyReducedEvents = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended>;
using MyReducedEventsVtxCov = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsVtxCov>;

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

static o2::globaltracking::MatchGlobalFwd sExtrap;

using o2::dataformats::GlobalFwdTrack;
using o2::track::TrackParCovFwd;

struct muonGlobalAlignment {
  ////   Variables for selecting MCH and MFT tracks
  Configurable<float> fTrackChi2MchUp{"cfgTrackChi2MchUp", 5.f, ""};
  Configurable<float> fPtMchLow{"cfgPtMchLow", 0.7f, ""};
  Configurable<float> fEtaMftLow{"cfgEtaMftlow", -3.6f, ""};
  Configurable<float> fEtaMftUp{"cfgEtaMftup", -2.5f, ""};
  Configurable<float> fRabsLow{"cfgRabsLow", 17.6f, ""};
  Configurable<float> fRabsUp{"cfgRabsUp", 89.5f, ""};
  Configurable<float> fSigmaPdcaUp{"cfgPdcaUp", 6.f, ""};

  Configurable<int> fTrackNClustMftLow{"cfgTrackNClustMftLow", 7, ""};
  Configurable<float> fTrackChi2MftUp{"cfgTrackChi2MftUp", 999.f, ""};

  Configurable<float> fMftMchResidualsPLow{"cfgMftMchResidualsPLow", 30.f, ""};
  Configurable<float> fMftMchResidualsPtLow{"cfgMftMchResidualsPtLow", 4.f, ""};

  Configurable<uint32_t> fMftTracksMultiplicityMax{"cfgMftTracksMultiplicityMax", 0, "Maximum number of MFT tracks to be processed per event (zero means no limit)"};

  Configurable<float> fVertexZshift{"cfgVertexZshift", 0.0f, "Correction to the vertex z position"};

  ////   Variables for MFT alignment corrections
  struct : ConfigurableGroup {
    Configurable<bool> fEnableMFTAlignmentCorrections{"cfgEnableMFTAlignmentCorrections", false, ""};
    // slope corrections
    Configurable<float> fMFTAlignmentCorrXSlopeTop{"cfgMFTAlignmentCorrXSlopeTop", (-0.0006696 - 0.0005621) / 2.f, "MFT X slope correction - top half"};
    Configurable<float> fMFTAlignmentCorrXSlopeBottom{"cfgMFTAlignmentCorrXSlopeBottom", (0.00105 + 0.001007) / 2.f, "MFT X slope correction - bottom half"};
    Configurable<float> fMFTAlignmentCorrYSlopeTop{"cfgMFTAlignmentCorrYSlopeTop", (-0.002299 - 0.002442) / 2.f, "MFT Y slope correction - top half"};
    Configurable<float> fMFTAlignmentCorrYSlopeBottom{"cfgMFTAlignmentCorrYSlopeBottom", (-0.0005339 - 0.0006921) / 2.f, "MFT Y slope correction - bottom half"};
    // offset corrections
    Configurable<float> fMFTAlignmentCorrXOffsetTop{"cfgMFTAlignmentCorrXOffsetTop", 0.f, "MFT X offset correction - top half"};
    Configurable<float> fMFTAlignmentCorrXOffsetBottom{"cfgMFTAlignmentCorrXOffsetBottom", 0.f, "MFT X offset correction - bottom half"};
    Configurable<float> fMFTAlignmentCorrYOffsetTop{"cfgMFTAlignmentCorrYOffsetTop", 0.f, "MFT Y offset correction - top half"};
    Configurable<float> fMFTAlignmentCorrYOffsetBottom{"cfgMFTAlignmentCorrYOffsetBottom", 0.f, "MFT Y offset correction - bottom half"};
  } configMFTAlignmentCorrections;

  ////   Variables for re-alignment setup
  struct : ConfigurableGroup {
    Configurable<bool> fEnableMCHRealign{"cfgEnableMCHRealign", true, "Enable re-alignment of MCH clusters and tracks"};
    Configurable<double> fChamberResolutionX{"cfgChamberResolutionX", 0.4, "Chamber resolution along X configuration for refit"}; // 0.4cm pp, 0.2cm PbPb
    Configurable<double> fChamberResolutionY{"cfgChamberResolutionY", 0.4, "Chamber resolution along Y configuration for refit"}; // 0.4cm pp, 0.2cm PbPb
    Configurable<double> fSigmaCutImprove{"cfgSigmaCutImprove", 6., "Sigma cut for track improvement"};
  } configRealign;

  ////   Variables for ccdb
  struct : ConfigurableGroup {
    Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
    Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
    Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
    Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
    // Configurable<std::string> geoPathRealign{"geoPathRealign", "Users/j/jcastill/GeometryAlignedFix10Fix15ShiftCh1BNew2", "Path of the geometry file"};
    Configurable<std::string> geoPathRealign{"geoPathRealign", "Users/j/jcastill/GeometryAlignedLoczzm4pLHC24anap1sR5a", "Path of the geometry file"};
    Configurable<int64_t> nolaterthan{"ccdb-no-later-than-ref", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object of reference basis"};
    Configurable<int64_t> nolaterthanRealign{"ccdb-no-later-than-new", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object of new basis"};
  } configCCDB;

  Configurable<bool> fRequireGoodRCT{"cfgRequireGoodRCT", true, "Require good detector flags in Run Condition Table"};

  Configurable<bool> fEnableVertexShiftAnalysis{"cfgEnableVertexShiftAnalysis", true, "Enable the analysis of vertex shift"};
  Configurable<bool> fEnableMftDcaAnalysis{"cfgEnableMftDcaAnalysis", true, "Enable the analysis of DCA-based MFT alignment"};
  Configurable<bool> fEnableMftMchResidualsAnalysis{"cfgEnableMftMchResidualsAnalysis", true, "Enable the analysis of residuals between MFT tracks and MCH clusters"};

  int mRunNumber{0}; // needed to detect if the run changed and trigger update of magnetic field

  Service<o2::ccdb::BasicCCDBManager> ccdbManager;
  o2::field::MagneticField* fieldB;
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

  geo::TransformationCreator transformation;
  std::map<int, math_utils::Transform3D> transformRef; // reference geometry w.r.t track data
  std::map<int, math_utils::Transform3D> transformNew; // new geometry
  TGeoManager* geoNew = nullptr;
  TGeoManager* geoRef = nullptr;
  TrackFitter trackFitter; // Track fitter from MCH tracking library
  double mImproveCutChi2;  // Chi2 cut for track improvement.

  Preslice<aod::FwdTrkCl> perMuon = aod::fwdtrkcl::fwdtrackId;

  o2::aod::rctsel::RCTFlagsChecker rctChecker{"CBT_muon_glo", false, false, true};

  double mBzAtMftCenter{0};

  HistogramRegistry registry{"registry", {}};

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

  void InitCollisions(MyEvents const& collisions,
                      MyBCs const& bcs,
                      MyMuonsWithCov const& muonTracks,
                      std::map<uint64_t, CollisionInfo>& collisionInfos)
  {
    // fill collision information for global muon tracks (MFT-MCH-MID matches)
    for (auto muonTrack : muonTracks) {
      if (!muonTrack.has_collision())
        continue;

      auto collision = collisions.rawIteratorAt(muonTrack.collisionId());

      if (fRequireGoodRCT && !rctChecker(collision))
        continue;

      uint64_t collisionIndex = collision.globalIndex();

      auto bc = bcs.rawIteratorAt(collision.bcId());

      auto& collisionInfo = collisionInfos[collisionIndex];
      collisionInfo.bc = bc.globalBC();
      collisionInfo.zVertex = collision.posZ();

      if (static_cast<int>(muonTrack.trackType()) > 2) {
        // standalone MCH or MCH-MID tracks
        uint64_t mchTrackIndex = muonTrack.globalIndex();
        collisionInfo.mchTracks.push_back(mchTrackIndex);
      } else {
        // global muon tracks (MFT-MCH or MFT-MCH-MID)
        uint64_t muonTrackIndex = muonTrack.globalIndex();
        auto const& mchTrack = muonTrack.template matchMCHTrack_as<MyMuonsWithCov>();
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

    for (auto& [collisionIndex, collisionInfo] : collisionInfos) {
      for (auto& [mchIndex, globalTracksVector] : collisionInfo.globalMuonTracks) {
        std::sort(globalTracksVector.begin(), globalTracksVector.end(), compareChi2);
      }
    }
  }

  void InitCollisions(MyEvents const& collisions,
                      MyBCs const& bcs,
                      MyMuonsWithCov const& muonTracks,
                      MyMFTs const& mftTracks,
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

  template <typename BC>
  void initCCDB(BC const& bc)
  {
    if (mRunNumber == bc.runNumber())
      return;

    mRunNumber = bc.runNumber();
    ccdbManager->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    // auto grpmag = ccdbApi.retrieveFromTFileAny<o2::parameters::GRPMagField>(grpmagPath, metadata, ts);
    auto grpmag = ccdbManager->getForTimeStamp<o2::parameters::GRPMagField>(configCCDB.grpmagPath, bc.timestamp());
    if (grpmag != nullptr) {
      base::Propagator::initFieldFromGRP(grpmag);
      TrackExtrap::setField();
      TrackExtrap::useExtrapV2();
      fieldB = static_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField()); // for MFT
      double centerMFT[3] = {0, 0, -61.4};                                                         // or use middle point between Vtx and MFT?
      mBzAtMftCenter = fieldB->getBz(centerMFT);
    } else {
      LOGF(fatal, "GRP object is not available in CCDB at timestamp=%llu", bc.timestamp());
    }

    // Load geometry information from CCDB/local
    LOGF(info, "Loading reference aligned geometry from CCDB no later than %d", configCCDB.nolaterthan.value);
    ccdbManager->setCreatedNotAfter(configCCDB.nolaterthan); // this timestamp has to be consistent with what has been used in reco
    geoRef = ccdbManager->getForTimeStamp<TGeoManager>(configCCDB.geoPath, bc.timestamp());
    ccdbManager->clearCache(configCCDB.geoPath);

    if (configRealign.fEnableMCHRealign && fEnableMftMchResidualsAnalysis) {
      if (geoRef != nullptr) {
        transformation = geo::transformationFromTGeoManager(*geoRef);
      } else {
        LOGF(fatal, "Reference aligned geometry object is not available in CCDB at timestamp=%llu", bc.timestamp());
      }
      for (int i = 0; i < 156; i++) {
        int iDEN = GetDetElemId(i);
        transformRef[iDEN] = transformation(iDEN);
      }

      LOGF(info, "Loading new aligned geometry from CCDB no later than %d", configCCDB.nolaterthanRealign.value);
      ccdbManager->setCreatedNotAfter(configCCDB.nolaterthanRealign); // make sure this timestamp can be resolved regarding the reference one
      geoNew = ccdbManager->getForTimeStamp<TGeoManager>(configCCDB.geoPathRealign, bc.timestamp());
      ccdbManager->clearCache(configCCDB.geoPathRealign);
      if (geoNew != nullptr) {
        transformation = geo::transformationFromTGeoManager(*geoNew);
      } else {
        LOGF(fatal, "New aligned geometry object is not available in CCDB at timestamp=%llu", bc.timestamp());
      }
      for (int i = 0; i < 156; i++) {
        int iDEN = GetDetElemId(i);
        transformNew[iDEN] = transformation(iDEN);
      }
    }
  }

  void init(o2::framework::InitContext&)
  {
    // Load geometry
    ccdbManager->setURL(configCCDB.ccdburl);
    ccdbManager->setCaching(true);
    ccdbManager->setLocalObjectValidityChecking();
    ccdbApi.init(configCCDB.ccdburl);
    mRunNumber = 0;

    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      LOGF(info, "Load geometry from CCDB");
      ccdbManager->get<TGeoManager>(configCCDB.geoPath);
    }

    // Configuration for track fitter
    const auto& trackerParam = TrackerParam::Instance();
    trackFitter.setBendingVertexDispersion(trackerParam.bendingVertexDispersion);
    trackFitter.setChamberResolution(configRealign.fChamberResolutionX, configRealign.fChamberResolutionY);
    trackFitter.smoothTracks(true);
    trackFitter.useChamberResolution();
    mImproveCutChi2 = 2. * configRealign.fSigmaCutImprove * configRealign.fSigmaCutImprove;

    float mftLadderWidth = 1.7;
    AxisSpec dcaxMFTAxis = {400, -0.5, 0.5, "DCA_{x} (cm)"};
    AxisSpec dcayMFTAxis = {400, -0.5, 0.5, "DCA_{y} (cm)"};
    AxisSpec dcaxMCHAxis = {400, -10.0, 10.0, "DCA_{x} (cm)"};
    AxisSpec dcayMCHAxis = {400, -10.0, 10.0, "DCA_{y} (cm)"};
    AxisSpec dcazAxis = {20, -10.0, 10.0, "v_{z} (cm)"};
    AxisSpec txAxis = {30, -mftLadderWidth * 15.f / 2.f, mftLadderWidth * 15.f / 2.f, "track_{x} (cm)"};
    AxisSpec tyAxis = {20, -10.f, 10.f, "track_{y} (cm)"};
    AxisSpec vxAxis = {400, -0.5, 0.5, "vtx_{x} (cm)"};
    AxisSpec vyAxis = {400, -0.5, 0.5, "vtx_{y} (cm)"};
    AxisSpec vzAxis = {1000, -10.0, 10.0, "vtx_{z} (cm)"};
    AxisSpec phiAxis = {36, -180.0, 180.0, "#phi (degrees)"};
    AxisSpec nMftClustersAxis = {6, 5, 11, "# of clusters"};
    AxisSpec mftTrackTypeAxis = {2, 0, 2, "track type"};
    AxisSpec trackChargeSignAxis = {2, 0, 0, "sign"};
    AxisSpec layersPatternAxis = {1024, 0, 1024, "layers pattern"};
    AxisSpec zshiftAxis = {21, -5.25, 5.25, "z shift (mm)"};

    registry.add("vertex_y_vs_x", std::format("Vertex y vs. x").c_str(), {HistType::kTH2F, {vxAxis, vyAxis}});
    registry.add("vertex_z", std::format("Vertex z").c_str(), {HistType::kTH1F, {vzAxis}});

    if (fEnableVertexShiftAnalysis || fEnableMftDcaAnalysis) {
      registry.add("DCA/MFT/nTracksMFT", std::format("Number of MFT tracks per collision").c_str(), {HistType::kTH1F, {{100, 0, 1000, "# of MFT tracks"}}});
      registry.add("DCA/MFT/DCA_y_vs_x", std::format("DCA y vs. x").c_str(), {HistType::kTH2F, {dcaxMFTAxis, dcayMFTAxis}});
    }

    if (fEnableVertexShiftAnalysis) {
      registry.add("DCA/MFT/DCA_x_vs_phi_vs_zshift", std::format("DCA(x) vs. #phi vs. z shift").c_str(), {HistType::kTH3F, {zshiftAxis, phiAxis, dcaxMFTAxis}});
      registry.add("DCA/MFT/DCA_y_vs_phi_vs_zshift", std::format("DCA(y) vs. #phi vs. z shift").c_str(), {HistType::kTH3F, {zshiftAxis, phiAxis, dcayMFTAxis}});
    }

    if (fEnableMftDcaAnalysis) {
      registry.add("DCA/MFT/DCA_x", "DCA(x) vs. vz, tx, ty, nclus, trackType",
                   HistType::kTHnSparseF, {dcaxMFTAxis, dcazAxis, txAxis, tyAxis, nMftClustersAxis, mftTrackTypeAxis});
      registry.add("DCA/MFT/DCA_y", "DCA(y) vs. vz, tx, ty, nclus, trackType",
                   HistType::kTHnSparseF, {dcayMFTAxis, dcazAxis, txAxis, tyAxis, nMftClustersAxis, mftTrackTypeAxis});
      registry.add("DCA/MFT/layers", "Layers pattern vs. tx, ty, nclus, trackType",
                   HistType::kTHnSparseF, {layersPatternAxis, txAxis, tyAxis, nMftClustersAxis, mftTrackTypeAxis});
    }

    if (fEnableMftMchResidualsAnalysis) {
      registry.add("DCA/MCH/DCA_y_vs_x", std::format("DCA y vs. x").c_str(), {HistType::kTH2F, {dcaxMCHAxis, dcayMCHAxis}});
      registry.add("DCA/MCH/DCA_x_vs_sign_vs_quadrant_vs_vz", std::format("DCA(x) vs. vz, quadrant, chargeSign").c_str(), {HistType::kTHnSparseF, {dcazAxis, {4, 0, 4, "quadrant"}, {2, 0, 2, "sign"}, dcaxMCHAxis}});
      registry.add("DCA/MCH/DCA_y_vs_sign_vs_quadrant_vs_vz", std::format("DCA(y) vs. vz, quadrant, chargeSign").c_str(), {HistType::kTHnSparseF, {dcazAxis, {4, 0, 4, "quadrant"}, {2, 0, 2, "sign"}, dcayMCHAxis}});

      //--
      AxisSpec dxAxis = {600, -30.0, 30.0, "#Delta x (cm)"};
      AxisSpec dyAxis = {600, -30.0, 30.0, "#Delta y (cm)"};

      registry.add("residuals/dx_vs_chamber", "Cluster x residual vs. chamber, quadrant, chargeSign",
                   {HistType::kTHnSparseF, {{10, 1, 11, "chamber"}, {4, 0, 4, "quadrant"}, {2, 0, 2, "sign"}, dxAxis}});
      registry.add("residuals/dy_vs_chamber", "Cluster y residual vs. chamber, quadrant, chargeSign",
                   {HistType::kTHnSparseF, {{10, 1, 11, "chamber"}, {4, 0, 4, "quadrant"}, {2, 0, 2, "sign"}, dyAxis}});

      registry.add("residuals/dx_vs_de", "Cluster x residual vs. DE, quadrant, chargeSign",
                   {HistType::kTHnSparseF, {{getNumDE(), 0, static_cast<double>(getNumDE()), "DE"}, {4, 0, 4, "quadrant"}, {2, 0, 2, "sign"}, dxAxis}});
      registry.add("residuals/dy_vs_de", "Cluster y residual vs. DE, quadrant, chargeSign",
                   {HistType::kTHnSparseF, {{getNumDE(), 0, static_cast<double>(getNumDE()), "DE"}, {4, 0, 4, "quadrant"}, {2, 0, 2, "sign"}, dyAxis}});
    }
  }

  int GetDetElemId(int iDetElemNumber)
  {
    const int fgNCh = 10;
    const int fgNDetElemCh[fgNCh] = {4, 4, 4, 4, 18, 18, 26, 26, 26, 26};
    const int fgSNDetElemCh[fgNCh + 1] = {0, 4, 8, 12, 16, 34, 52, 78, 104, 130, 156};

    // make sure detector number is valid
    if (!(iDetElemNumber >= fgSNDetElemCh[0] &&
          iDetElemNumber < fgSNDetElemCh[10])) {
      LOGF(fatal, "Invalid detector element number: %d", iDetElemNumber);
    }
    /// get det element number from ID
    // get chamber and element number in chamber
    int iCh = 0;
    int iDet = 0;
    for (int i = 1; i <= 10; i++) {
      if (iDetElemNumber < fgSNDetElemCh[i]) {
        iCh = i;
        iDet = iDetElemNumber - fgSNDetElemCh[i - 1];
        break;
      }
    }

    // make sure detector index is valid
    if (!(iCh > 0 && iCh <= 10 && iDet < fgNDetElemCh[iCh - 1])) {
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
    static int nDE = -1;
    if (nDE < 0) {
      for (int c = 0; c < 10; c++) {
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

  int GetQuadrant(double phi)
  {
    if (phi >= 0 && phi < 90) {
      return 0;
    }
    if (phi >= 90 && phi <= 180) {
      return 1;
    }
    if (phi >= -180 && phi < -90) {
      return 2;
    }
    if (phi >= -90 && phi < 0) {
      return 3;
    }
    return -1;
  }

  template <class T>
  int GetQuadrant(const T& track)
  {
    double phi = track.phi() * 180 / TMath::Pi();
    return GetQuadrant(phi);
  }

  template <class T>
  static o2::mch::TrackParam FwdtoMCH(const T& fwdtrack)
  {
    using SMatrix55Std = ROOT::Math::SMatrix<double, 5>;
    using SMatrix55Sym = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;

    // Convert Forward Track parameters and covariances matrix to the
    // MCH track format.

    // Parameter conversion
    double alpha1, alpha3, alpha4, x2, x3, x4;

    x2 = fwdtrack.getPhi();
    x3 = fwdtrack.getTanl();
    x4 = fwdtrack.getInvQPt();

    auto sinx2 = TMath::Sin(x2);
    auto cosx2 = TMath::Cos(x2);

    alpha1 = cosx2 / x3;
    alpha3 = sinx2 / x3;
    alpha4 = x4 / TMath::Sqrt(x3 * x3 + sinx2 * sinx2);

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

    double cov[] = {covariances(0, 0), covariances(1, 0), covariances(1, 1), covariances(2, 0), covariances(2, 1), covariances(2, 2), covariances(3, 0), covariances(3, 1), covariances(3, 2), covariances(3, 3), covariances(4, 0), covariances(4, 1), covariances(4, 2), covariances(4, 3), covariances(4, 4)};
    double param[] = {fwdtrack.getX(), alpha1, fwdtrack.getY(), alpha3, alpha4};

    o2::mch::TrackParam convertedTrack(fwdtrack.getZ(), param, cov);
    return o2::mch::TrackParam(convertedTrack);
  }

  static o2::dataformats::GlobalFwdTrack MCHtoFwd(const o2::mch::TrackParam& mchParam)
  {
    using SMatrix55Std = ROOT::Math::SMatrix<double, 5>;
    using SMatrix55Sym = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;

    // Convert a MCH Track parameters and covariances matrix to the
    // Forward track format. Must be called after propagation though the absorber

    o2::dataformats::GlobalFwdTrack convertedTrack;

    // Parameter conversion
    double alpha1, alpha3, alpha4, x2, x3, x4;

    alpha1 = mchParam.getNonBendingSlope();
    alpha3 = mchParam.getBendingSlope();
    alpha4 = mchParam.getInverseBendingMomentum();

    x2 = TMath::ATan2(-alpha3, -alpha1);
    x3 = -1. / TMath::Sqrt(alpha3 * alpha3 + alpha1 * alpha1);
    x4 = alpha4 * -x3 * TMath::Sqrt(1 + alpha3 * alpha3);

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

      if (track.getNClusters() < 10) {
        removeTrack = true;
        break;
      } else {
        while (itNextToNextParam != track.end()) {
          if (itNextToNextParam->getClusterPtr()->getChamberId() != itNextParam->getClusterPtr()->getChamberId()) {
            itStartingParam = std::make_reverse_iterator(++itNextParam);
            break;
          }
          ++itNextToNextParam;
        }
      }
    }

    if (!removeTrack) {
      for (auto& param : track) {
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
  bool IsGoodMuon(const T& muonTrack, const C& collision,
                  double chi2Cut,
                  double pCut,
                  double pTCut,
                  std::array<double, 2> etaCut,
                  std::array<double, 2> rAbsCut,
                  double nSigmaPdcaCut)
  {
    auto const& mchTrack = (static_cast<int>(muonTrack.trackType()) <= 2) ? muonTrack.template matchMCHTrack_as<MyMuonsWithCov>() : muonTrack;

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

  void TransformMFT(o2::mch::TrackParam& track)
  {
    // double zCH10 = -1437.6;
    double z = track.getZ();
    // double dZ = zMCH - z;
    double x = track.getNonBendingCoor();
    double y = track.getBendingCoor();
    double xSlope = track.getNonBendingSlope();
    double ySlope = track.getBendingSlope();

    double xSlopeCorrection = (y > 0) ? configMFTAlignmentCorrections.fMFTAlignmentCorrXSlopeTop : configMFTAlignmentCorrections.fMFTAlignmentCorrXSlopeBottom;
    double xCorrection = xSlopeCorrection * z +
                         ((y > 0) ? configMFTAlignmentCorrections.fMFTAlignmentCorrXOffsetTop : configMFTAlignmentCorrections.fMFTAlignmentCorrXOffsetBottom);
    track.setNonBendingCoor(x + xCorrection);
    track.setNonBendingSlope(xSlope + xSlopeCorrection);

    double ySlopeCorrection = (y > 0) ? configMFTAlignmentCorrections.fMFTAlignmentCorrYSlopeTop : configMFTAlignmentCorrections.fMFTAlignmentCorrYSlopeBottom;
    double yCorrection = ySlopeCorrection * z +
                         ((y > 0) ? configMFTAlignmentCorrections.fMFTAlignmentCorrYOffsetTop : configMFTAlignmentCorrections.fMFTAlignmentCorrYOffsetBottom);
    track.setBendingCoor(y + yCorrection);
    track.setBendingSlope(ySlope + ySlopeCorrection);
    /*
    std::cout << std::format("[TOTO] MFT position:    pos={:0.3f},{:0.3f}", x, y) << std::endl;
    std::cout << std::format("[TOTO] MFT corrections: pos={:0.3f},{:0.3f}  slope={:0.12f},{:0.12f}  angle={:0.12f},{:0.12f}",
        xCorrection, yCorrection, xSlopeCorrection, ySlopeCorrection,
        std::atan2(xSlopeCorrection, 1), std::atan2(ySlopeCorrection, 1)) << std::endl;
    */
  }

  void TransformMFT(o2::dataformats::GlobalFwdTrack& track)
  {
    auto mchTrack = FwdtoMCH(track);

    TransformMFT(mchTrack);

    auto transformedTrack = sExtrap.MCHtoFwd(mchTrack);
    track.setParameters(transformedTrack.getParameters());
    track.setZ(transformedTrack.getZ());
    track.setCovariances(transformedTrack.getCovariances());
  }

  void TransformMFT(o2::track::TrackParCovFwd& fwdtrack)
  {
    o2::dataformats::GlobalFwdTrack track;
    track.setParameters(fwdtrack.getParameters());
    track.setZ(fwdtrack.getZ());
    track.setCovariances(fwdtrack.getCovariances());

    auto mchTrack = FwdtoMCH(track);

    TransformMFT(mchTrack);

    auto transformedTrack = sExtrap.MCHtoFwd(mchTrack);
    fwdtrack.setParameters(transformedTrack.getParameters());
    fwdtrack.setZ(transformedTrack.getZ());
    fwdtrack.setCovariances(transformedTrack.getCovariances());
  }

  o2::dataformats::GlobalFwdTrack PropagateMCHParam(mch::TrackParam mchTrack, const double z)
  {
    float absFront = -90.f;
    float absBack = -505.f;

    if (mchTrack.getZ() < absBack && z > absFront) {
      // extrapolation through the absorber in the upstream direction
      o2::mch::TrackExtrap::extrapToVertexWithoutBranson(mchTrack, z);
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

    float absFront = -90.f;
    float absBack = -505.f;

    if (muon.getZ() < absBack && z > absFront) {
      // extrapolation through the absorber in the upstream direction
      o2::mch::TrackExtrap::extrapToVertexWithoutBranson(mchTrack, z);
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

  template <class TMFT>
  o2::dataformats::GlobalFwdTrack PropagateMFT(const TMFT& mftTrack, float z)
  {
    static double Bz = -10001;
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
    if (Bz < -10000) {
      double centerZ[3] = {0, 0, (-45.f - 77.5f) / 2.f};
      o2::field::MagneticField* field = static_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField());
      Bz = field->getBz(centerZ);
    }
    fwdtrack.propagateToZ(z, Bz);

    propmuon.setParameters(fwdtrack.getParameters());
    propmuon.setZ(fwdtrack.getZ());
    propmuon.setCovariances(fwdtrack.getCovariances());

    return propmuon;
  }

  template <class TMFT, class C>
  o2::dataformats::GlobalFwdTrack PropagateMFTToDCA(const TMFT& mftTrack, const C& collision, float zshift = 0)
  {
    static double Bz = -10001;
    double chi2 = mftTrack.chi2();
    double phiCorrDeg = 0;
    double phiCorr = phiCorrDeg * TMath::Pi() / 180.f;
    double tR = std::hypot(mftTrack.x(), mftTrack.y());
    double tphi = std::atan2(mftTrack.y(), mftTrack.x());
    double tx = std::cos(tphi + phiCorr) * tR;
    double ty = std::sin(tphi + phiCorr) * tR;
    SMatrix5 tpars = {tx, ty, mftTrack.phi() + phiCorr, mftTrack.tgl(), mftTrack.signed1Pt()};
    std::vector<double> v1{0, 0, 0, 0, 0,
                           0, 0, 0, 0, 0,
                           0, 0, 0, 0, 0};
    SMatrix55 tcovs(v1.begin(), v1.end());
    o2::track::TrackParCovFwd fwdtrack{mftTrack.z(), tpars, tcovs, chi2};
    if (configMFTAlignmentCorrections.fEnableMFTAlignmentCorrections) {
      TransformMFT(fwdtrack);
    }
    o2::dataformats::GlobalFwdTrack propmuon;

    // double propVec[3] = {};
    // propVec[0] = collision.posX() - mftTrack.x();
    // propVec[1] = collision.posY() - mftTrack.y();
    // propVec[2] = collision.posZ() - mftTrack.z();

    // double centerZ[3] = {mftTrack.x() + propVec[0] / 2.,
    //                      mftTrack.y() + propVec[1] / 2.,
    //                      mftTrack.z() + propVec[2] / 2.};
    if (Bz < -10000) {
      double centerZ[3] = {0, 0, -45.f / 2.f};
      o2::field::MagneticField* field = static_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField());
      Bz = field->getBz(centerZ);
    }
    fwdtrack.propagateToZ(collision.posZ() - zshift, Bz);

    propmuon.setParameters(fwdtrack.getParameters());
    propmuon.setZ(fwdtrack.getZ());
    propmuon.setCovariances(fwdtrack.getCovariances());

    return propmuon;
  }

  template <typename T>
  T UpdateTrackMomentum(const T& track, const double p, int sign)
  {
    double px = p * sin(M_PI / 2 - atan(track.tgl())) * cos(track.phi());
    double py = p * sin(M_PI / 2 - atan(track.tgl())) * sin(track.phi());
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
    double px = track4mom.p() * sin(M_PI / 2 - atan(track.tgl())) * cos(track.phi());
    double py = track4mom.p() * sin(M_PI / 2 - atan(track.tgl())) * sin(track.phi());
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

  template <typename TMCH, typename TMFT>
  o2::dataformats::GlobalFwdTrack PropagateMFTtoMCH(const TMFT& mftTrack, const TMCH& mchTrack, const double z)
  {
    // extrapolation with MCH tools
    auto mchTrackAtMFT = FwdtoMCH(FwdToTrackPar(mchTrack));
    o2::mch::TrackExtrap::extrapToVertexWithoutBranson(mchTrackAtMFT, mftTrack.z());

    auto mftTrackPar = FwdToTrackPar(mftTrack);
    if (configMFTAlignmentCorrections.fEnableMFTAlignmentCorrections) {
      TransformMFT(mftTrackPar);
    }
    auto mftTrackProp = FwdtoMCH(mftTrackPar);
    UpdateTrackMomentum(mftTrackProp, mchTrackAtMFT);
    if (z < -505.f) {
      o2::mch::TrackExtrap::extrapToZ(mftTrackProp, -466.f);
      UpdateTrackMomentum(mftTrackProp, sExtrap.FwdtoMCH(FwdToTrackPar(mchTrack)));
    }
    o2::mch::TrackExtrap::extrapToZ(mftTrackProp, z);

    return MCHtoFwd(mftTrackProp);
  }

  template <typename TMCH, typename TMFT>
  o2::dataformats::GlobalFwdTrack PropagateMFTtoMCH_(const TMFT& mftTrack, const TMCH& mchTrack, const double z)
  {
    // extrapolation with MCH tools
    auto mftTrackProp = sExtrap.FwdtoMCH(FwdToTrackPar(mftTrack));
    UpdateTrackMomentum(mftTrackProp, sExtrap.FwdtoMCH(FwdToTrackPar(mchTrack)));
    o2::mch::TrackExtrap::extrapToZ(mftTrackProp, z);

    return sExtrap.MCHtoFwd(mftTrackProp);
  }

  void FillDCAPlots(MyEvents const& collisions,
                    MyBCs const& bcs,
                    MyMuonsWithCov const& /*muonTracks*/,
                    MyMFTs const& mftTracks,
                    const std::map<uint64_t, CollisionInfo>& collisionInfos)
  {
    // outer loop over collisions
    for (auto& [collisionIndex, collisionInfo] : collisionInfos) {
      auto const& collision = collisions.rawIteratorAt(collisionIndex);
      const auto& bc = bcs.rawIteratorAt(collision.bcId());

      // remove TF/ROF borders and ambiguous collisions
      if (!bc.selection_bit(o2::aod::evsel::kNoTimeFrameBorder) ||
          !bc.selection_bit(o2::aod::evsel::kNoITSROFrameBorder))
        continue;

      registry.get<TH2>(HIST("vertex_y_vs_x"))->Fill(collision.posX(), collision.posY());
      registry.get<TH1>(HIST("vertex_z"))->Fill(collision.posZ());

      if (fEnableVertexShiftAnalysis || fEnableMftDcaAnalysis) {
        registry.get<TH1>(HIST("DCA/MFT/nTracksMFT"))->Fill(collisionInfo.mftTracks.size());
      }

      if (fEnableVertexShiftAnalysis || fEnableMftDcaAnalysis) {
        // loop over MFT tracks
        auto mftTrackIds = collisionInfo.mftTracks;
        if (fMftTracksMultiplicityMax > 0 && mftTrackIds.size() > fMftTracksMultiplicityMax) {
          auto rng = std::default_random_engine{};
          std::shuffle(std::begin(mftTrackIds), std::end(mftTrackIds), rng);
          mftTrackIds.resize(fMftTracksMultiplicityMax);
        }

        for (auto mftIndex : mftTrackIds) {
          auto const& mftTrack = mftTracks.rawIteratorAt(mftIndex);

          bool isGoodMFT = IsGoodMFT(mftTrack, 999.f, 5);
          if (!isGoodMFT)
            continue;

          auto mftTrackAtDCA = PropagateMFTToDCA(mftTrack, collision, fVertexZshift);
          double dcax = mftTrackAtDCA.getX() - collision.posX();
          double dcay = mftTrackAtDCA.getY() - collision.posY();
          double phi = mftTrack.phi() * 180 / TMath::Pi();
          int mftNclusters = mftTrack.nClusters();
          int mftTrackType = mftTrack.isCA() ? 1 : 0;

          const int nMftLayers = 10;
          int layerPattern = 0;
          for (int layer = 0; layer < nMftLayers; layer++) {
            if ((mftTrack.mftClusterSizesAndTrackFlags() >> (layer * 6)) & 0x3F) {
              layerPattern += (1 << layer);
            }
          }

          if (fEnableMftDcaAnalysis) {
            registry.get<TH2>(HIST("DCA/MFT/DCA_y_vs_x"))->Fill(dcax, dcay);
            registry.get<THnSparse>(HIST("DCA/MFT/DCA_x"))->Fill(dcax, collision.posZ(), mftTrack.x(), mftTrack.y(), mftNclusters, mftTrackType);
            registry.get<THnSparse>(HIST("DCA/MFT/DCA_y"))->Fill(dcay, collision.posZ(), mftTrack.x(), mftTrack.y(), mftNclusters, mftTrackType);

            registry.get<THnSparse>(HIST("DCA/MFT/layers"))->Fill(layerPattern, mftTrack.x(), mftTrack.y(), mftNclusters, mftTrackType);
          }

          if (fEnableVertexShiftAnalysis) {
            if (std::fabs(collision.posZ()) < 1.f) {
              float zshift[21] = {// in millimeters
                                  -5.0, -4.5, -4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0,
                                  0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0};
              for (int zi = 0; zi < 21; zi++) {
                auto mftTrackAtDCAshifted = PropagateMFTToDCA(mftTrack, collision, zshift[zi] / 10.f);
                double dcaxShifted = mftTrackAtDCAshifted.getX() - collision.posX();
                double dcayShifted = mftTrackAtDCAshifted.getY() - collision.posY();
                registry.get<TH3>(HIST("DCA/MFT/DCA_x_vs_phi_vs_zshift"))->Fill(zshift[zi], phi, dcaxShifted);
                registry.get<TH3>(HIST("DCA/MFT/DCA_y_vs_phi_vs_zshift"))->Fill(zshift[zi], phi, dcayShifted);
              }
            }
          }
        }
      }
    }
  }

  void FillResidualsPlots(MyEvents const& collisions,
                          MyBCs const& bcs,
                          MyMuonsWithCov const& muonTracks,
                          aod::FwdTrkCls const& clusters,
                          const std::map<uint64_t, CollisionInfo>& collisionInfos)
  {
    if (!fEnableMftMchResidualsAnalysis) {
      return;
    }

    // loop over collisions
    for (auto& [collisionIndex, collisionInfo] : collisionInfos) {
      auto const& collision = collisions.rawIteratorAt(collisionIndex);
      const auto& bc = bcs.rawIteratorAt(collision.bcId());

      // remove TF/ROF borders and ambiguous collisions
      if (!bc.selection_bit(o2::aod::evsel::kNoTimeFrameBorder) ||
          !bc.selection_bit(o2::aod::evsel::kNoITSROFrameBorder))
        continue;

      // loop over global muon tracks
      for (auto& [muonIndex, globalTracksVector] : collisionInfo.globalMuonTracks) {
        auto const& muonTrack = muonTracks.rawIteratorAt(globalTracksVector[0]);
        const auto& mchTrack = muonTrack.template matchMCHTrack_as<MyMuonsWithCov>();
        const auto& mftTrack = muonTrack.template matchMFTTrack_as<MyMFTs>();
        int quadrantMch = GetQuadrant(mchTrack);
        int posNeg = (mchTrack.sign() >= 0) ? 0 : 1;

        bool isGoodMuon = IsGoodMuon(mchTrack, collision, fTrackChi2MchUp, fMftMchResidualsPLow, fMftMchResidualsPtLow, {fEtaMftLow, fEtaMftUp}, {fRabsLow, fRabsUp}, fSigmaPdcaUp);
        if (!isGoodMuon)
          continue;

        bool isGoodMFT = IsGoodMFT(mftTrack, fTrackChi2MftUp, fTrackNClustMftLow);
        if (!isGoodMFT)
          continue;

        TrackRealigned convertedTrack;
        // loop over attached clusters
        int clIndex = -1;
        auto clustersSliced = clusters.sliceBy(perMuon, mchTrack.globalIndex()); // Slice clusters by muon id
        for (auto const& cluster : clustersSliced) {
          clIndex += 1;

          int deId = cluster.deId();
          int chamber = deId / 100 - 1;
          if (chamber < 0 || chamber > 9)
            continue;
          int deIndex = getDEindex(deId);

          math_utils::Point3D<double> local;
          math_utils::Point3D<double> master;

          master.SetXYZ(cluster.x(), cluster.y(), cluster.z());

          if (configRealign.fEnableMCHRealign) {
            // Transformation from reference geometry frame to new geometry frame
            transformRef[cluster.deId()].MasterToLocal(master, local);
            transformNew[cluster.deId()].LocalToMaster(local, master);

            mch::Cluster* clusterMCH = new mch::Cluster();
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

          auto mftTrackAtCluster = PropagateMFTtoMCH(mftTrack, mchTrack, master.z());

          std::array<double, 2> xPos{master.x(), mftTrackAtCluster.getX()};
          std::array<double, 2> yPos{master.y(), mftTrackAtCluster.getY()};

          registry.get<THnSparse>(HIST("residuals/dx_vs_chamber"))->Fill(chamber + 1, quadrantMch, posNeg, xPos[0] - xPos[1]);
          registry.get<THnSparse>(HIST("residuals/dy_vs_chamber"))->Fill(chamber + 1, quadrantMch, posNeg, yPos[0] - yPos[1]);

          registry.get<THnSparse>(HIST("residuals/dx_vs_de"))->Fill(deIndex, quadrantMch, posNeg, xPos[0] - xPos[1]);
          registry.get<THnSparse>(HIST("residuals/dy_vs_de"))->Fill(deIndex, quadrantMch, posNeg, yPos[0] - yPos[1]);
        }

        bool removable{false};
        if (configRealign.fEnableMCHRealign) {
          // Refit the re-aligned track
          if (convertedTrack.getNClusters() != 0) {
            removable = RemoveTrack(convertedTrack);
          } else {
            LOGF(fatal, "Muon track %d has no associated clusters.", mchTrack.globalIndex());
          }
        }

        if (!removable) {
          auto mchTrackAtDCA = configRealign.fEnableMCHRealign ? PropagateMCHRealigned(convertedTrack, collision.posZ()) : PropagateMCH(mchTrack, collision.posZ());
          auto dcax = mchTrackAtDCA.getX() - collision.posX();
          auto dcay = mchTrackAtDCA.getY() - collision.posY();

          registry.get<TH2>(HIST("DCA/MCH/DCA_y_vs_x"))->Fill(dcax, dcay);
          registry.get<THnSparse>(HIST("DCA/MCH/DCA_x_vs_sign_vs_quadrant_vs_vz"))->Fill(collision.posZ(), quadrantMch, posNeg, dcax);
          registry.get<THnSparse>(HIST("DCA/MCH/DCA_y_vs_sign_vs_quadrant_vs_vz"))->Fill(collision.posZ(), quadrantMch, posNeg, dcay);
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
    InitCollisions(collisions, bcs, muonTracks, mftTracks, collisionInfos);

    FillDCAPlots(collisions, bcs, muonTracks, mftTracks, collisionInfos);

    FillResidualsPlots(collisions, bcs, muonTracks, clusters, collisionInfos);
  }

  PROCESS_SWITCH(muonGlobalAlignment, processQA, "process qa", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<muonGlobalAlignment>(cfgc)};
};
