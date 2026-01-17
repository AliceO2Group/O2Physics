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

#include <algorithm>
#include <limits>
#include <map>
#include <random>
#include <string>
#include <unordered_map>
#include <vector>

using namespace o2;
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

//_________________________________________________________________________________________

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

//_________________________________________________________________________________________

struct muonGlobalAlignment {
  ////   Variables for selecting mft tracks
  Configurable<float> fTrackChi2MchUp{"cfgTrackChi2MchUp", 5.f, ""};
  Configurable<float> fPtMchLow{"cfgPtMchLow", 0.7f, ""};
  Configurable<float> fEtaMftLow{"cfgEtaMftlow", -3.6f, ""};
  Configurable<float> fEtaMftUp{"cfgEtaMftup", -2.5f, ""};
  Configurable<float> fRabsLow{"cfgRabsLow", 17.6f, ""};
  Configurable<float> fRabsUp{"cfgRabsUp", 89.5f, ""};
  Configurable<float> fSigmaPdcaUp{"cfgPdcaUp", 6.f, ""};

  Configurable<int> fTrackNClustMftLow{"cfgTrackNClustMftLow", 7, ""};
  Configurable<float> fTrackChi2MftUp{"cfgTrackChi2MftUp", 999.f, ""};

  Configurable<uint32_t> fMftTracksMultiplicityMax{"cfgMftTracksMultiplicityMax", 0, "Maximum number of MFT tracks to be processed per event (zero means no limit)"};

  Configurable<float> fVertexZshift{"cfgVertexZshift", 0.0f, "Correction to the vertex z position"};

  ////   Variables for ccdb
  Configurable<std::string> ccdburl{"cfgCcdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpmagPath{"cfgGrpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> geoPath{"cfgGeoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};

  Configurable<bool> fEnableVertexShiftAnalysis{"cfgEnableVertexShiftAnalysis", true, "Enable the analysis of vertex shift"};
  Configurable<bool> fEnableMftDcaAnalysis{"cfgEnableMftDcaAnalysis", true, "Enable the analysis of DCA-based MFT alignment"};
  Configurable<bool> fEnableMftMchResidualsAnalysis{"cfgEnableMftMchResidualsAnalysis", true, "Enable the analysis of residuals between MFT tracks and MCH clusters"};

  int mRunNumber{0}; // needed to detect if the run changed and trigger update of magnetic field

  Service<o2::ccdb::BasicCCDBManager> ccdbManager;
  o2::field::MagneticField* fieldB;
  o2::ccdb::CcdbApi ccdbApi;

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
    std::map<std::string, std::string> metadata;
    auto soreor = o2::ccdb::BasicCCDBManager::getRunDuration(ccdbApi, mRunNumber);
    auto ts = soreor.first;
    auto grpmag = ccdbApi.retrieveFromTFileAny<o2::parameters::GRPMagField>(grpmagPath, metadata, ts);
    o2::base::Propagator::initFieldFromGRP(grpmag);
    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      ccdbManager->get<TGeoManager>(geoPath);
    }
    o2::mch::TrackExtrap::setField();
    fieldB = static_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField());
    double centerMFT[3] = {0, 0, -61.4}; // Field at center of MFT
    mBzAtMftCenter = fieldB->getBz(centerMFT);
  }

  void init(o2::framework::InitContext&)
  {
    // Load geometry
    ccdbManager->setURL(ccdburl);
    ccdbManager->setCaching(true);
    ccdbManager->setLocalObjectValidityChecking();
    ccdbApi.init(ccdburl);
    mRunNumber = 0;

    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      LOGF(info, "Load geometry from CCDB");
      ccdbManager->get<TGeoManager>(geoPath);
    }

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

  template <class TMFT, class C>
  o2::dataformats::GlobalFwdTrack PropagateMftToDCA(const TMFT& mftTrack, const C& collision, float zshift = 0)
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
    auto mchTrackAtMFT = sExtrap.FwdtoMCH(FwdToTrackPar(mchTrack));
    o2::mch::TrackExtrap::extrapToVertexWithoutBranson(mchTrackAtMFT, mftTrack.z());

    auto mftTrackProp = sExtrap.FwdtoMCH(FwdToTrackPar(mftTrack));
    UpdateTrackMomentum(mftTrackProp, mchTrackAtMFT);
    if (z < -505.f) {
      o2::mch::TrackExtrap::extrapToZ(mftTrackProp, -466.f);
      UpdateTrackMomentum(mftTrackProp, sExtrap.FwdtoMCH(FwdToTrackPar(mchTrack)));
    }
    o2::mch::TrackExtrap::extrapToZ(mftTrackProp, z);

    return sExtrap.MCHtoFwd(mftTrackProp);
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
                    MyMuonsWithCov const& muonTracks,
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

      if (fEnableMftMchResidualsAnalysis) {
        // loop over muon tracks
        for (auto mchIndex : collisionInfo.mchTracks) {
          auto const& mchTrack = muonTracks.rawIteratorAt(mchIndex);
          int quadrant = GetQuadrant(mchTrack);
          bool isGoodMuon = IsGoodMuon(mchTrack, collision, fTrackChi2MchUp, 30.0, 4.0, {fEtaMftLow, fEtaMftUp}, {fRabsLow, fRabsUp}, fSigmaPdcaUp);
          if (!isGoodMuon)
            continue;
          int sign = (mchTrack.sign() > 0) ? 0 : 1;

          auto mchTrackAtDCA = VarManager::PropagateMuon(mchTrack, collision, VarManager::kToDCA);
          double dcax = mchTrackAtDCA.getX() - collision.posX();
          double dcay = mchTrackAtDCA.getY() - collision.posY();

          registry.get<TH2>(HIST("DCA/MCH/DCA_y_vs_x"))->Fill(dcax, dcay);
          registry.get<THnSparse>(HIST("DCA/MCH/DCA_x_vs_sign_vs_quadrant_vs_vz"))->Fill(collision.posZ(), quadrant, sign, dcax);
          registry.get<THnSparse>(HIST("DCA/MCH/DCA_y_vs_sign_vs_quadrant_vs_vz"))->Fill(collision.posZ(), quadrant, sign, dcay);
        }
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

          auto mftTrackAtDCA = PropagateMftToDCA(mftTrack, collision, fVertexZshift);
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
                auto mftTrackAtDCAshifted = PropagateMftToDCA(mftTrack, collision, zshift[zi] / 10.f);
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
        int quadrant = GetQuadrant(mftTrack);
        int posNeg = (mchTrack.sign() >= 0) ? 0 : 1;

        bool isGoodMuon = IsGoodMuon(mchTrack, collision, fTrackChi2MchUp, 20.0, fPtMchLow, {fEtaMftLow, fEtaMftUp}, {fRabsLow, fRabsUp}, fSigmaPdcaUp);
        if (!isGoodMuon)
          continue;

        bool isGoodMFT = IsGoodMFT(mftTrack, fTrackChi2MftUp, fTrackNClustMftLow);
        if (!isGoodMFT)
          continue;

        // loop over MCH tracks
        for (auto mchIndex : collisionInfo.mchTracks) {
          auto const& mchTrack2 = muonTracks.rawIteratorAt(mchIndex);

          // loop over attached clusters
          for (auto const& cluster : clusters) {

            if (cluster.template fwdtrack_as<MyMuonsWithCov>() != mchTrack2) {
              continue;
            }

            int deId = cluster.deId();
            int chamber = deId / 100 - 1;
            if (chamber < 0 || chamber > 9)
              continue;
            int deIndex = getDEindex(deId);

            double xCluster = cluster.x();
            double yCluster = cluster.y();
            double zCluster = cluster.z();

            auto mftTrackAtCluster = PropagateMFTtoMCH(mftTrack, mchTrack, zCluster);

            std::array<double, 2> xPos{xCluster, mftTrackAtCluster.getX()};
            std::array<double, 2> yPos{yCluster, mftTrackAtCluster.getY()};

            registry.get<THnSparse>(HIST("residuals/dx_vs_chamber"))->Fill(chamber + 1, quadrant, posNeg, xPos[0] - xPos[1]);
            registry.get<THnSparse>(HIST("residuals/dy_vs_chamber"))->Fill(chamber + 1, quadrant, posNeg, yPos[0] - yPos[1]);

            registry.get<THnSparse>(HIST("residuals/dx_vs_de"))->Fill(deIndex, quadrant, posNeg, xPos[0] - xPos[1]);
            registry.get<THnSparse>(HIST("residuals/dy_vs_de"))->Fill(deIndex, quadrant, posNeg, yPos[0] - yPos[1]);
          }
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
      // grpmag = ccdbManager->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, bc.timestamp());
      // if (grpmag != nullptr) {
      //   LOGF(info, "Init field from GRP");
      //   o2::base::Propagator::initFieldFromGRP(grpmag);
      // }
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
