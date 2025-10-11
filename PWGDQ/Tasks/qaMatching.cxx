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

#include <string>
#include <unordered_map>
#include <map>
#include <limits>

using namespace o2;
using namespace o2::framework;
using namespace o2::aod;

using MyCollisions = aod::Collisions;
using MyBCs = soa::Join<aod::BCs, aod::Timestamps>;
using MyEvents = soa::Join<aod::Collisions, aod::EvSels>;
using MyMuonsWithCov = soa::Join<aod::FwdTracks, aod::FwdTracksCov>;
using MyMuonsFullMC = soa::Join<aod::FwdTracks, aod::FwdTracksCov, aod::McFwdTrackLabels>;
using MyMFTs = aod::MFTTracks;
using MyMFTCovariances = aod::MFTTracksCov;
using MyMFTsMC = soa::Join<aod::MFTTracks, aod::McMFTTrackLabels>;

using MyCollision = MyCollisions::iterator;
using MyBC = MyBCs::iterator;
using MyMuon = MyMuonsWithCov::iterator;
using MyMuonMC = MyMuonsFullMC::iterator;
using MyMFT = MyMFTs::iterator;
using MyMFTCovariance = MyMFTCovariances::iterator;

using SMatrix55 = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
using SMatrix5 = ROOT::Math::SVector<Double_t, 5>;

std::unordered_map<int, std::vector<int64_t>> map_mfttracks;
std::unordered_map<int, std::vector<int64_t>> map_muontracks;
std::unordered_map<int, bool> map_collisions;
std::unordered_map<int, bool> map_has_mfttracks_collisions;
std::unordered_map<int, bool> map_has_muontracks_collisions;
std::unordered_map<int, float> map_vtxz;
std::unordered_map<int, int> map_nmfttrack;

const int fgNCh = 10;
const int fgNDetElemCh[fgNCh] = {4, 4, 4, 4, 18, 18, 26, 26, 26, 26};
const int fgSNDetElemCh[fgNCh + 1] = {0, 4, 8, 12, 16, 34, 52, 78, 104, 130, 156};

constexpr double muonMass = 0.1056584;
constexpr double muonMass2 = muonMass * muonMass;

// constexpr static uint32_t gkMuonDCAFillMapWithCov = VarManager::ObjTypes::ReducedMuon | VarManager::ObjTypes::ReducedMuonExtra | VarManager::ObjTypes::ReducedMuonCov | VarManager::ObjTypes::MuonDCA;

constexpr static int toVertex = 0; //VarManager::kToVertex;
constexpr static int toDCA = 1; //VarManager::kToDCA;
constexpr static int toRabs = 2; //VarManager::kToRabs;
constexpr static int toMFT = 3;

constexpr double firstMFTPlaneZ = o2::mft::constants::mft::LayerZCoordinate()[0];
constexpr double lastMFTPlaneZ = o2::mft::constants::mft::LayerZCoordinate()[9];
constexpr double firstMCHPlaneZ = -526.16;

static o2::globaltracking::MatchGlobalFwd sExtrap;

using o2::dataformats::GlobalFwdTrack;
using o2::track::TrackParCovFwd;
typedef std::function<double(const GlobalFwdTrack& mchtrack, const TrackParCovFwd& mfttrack)> MatchingFunc_t;
std::map<std::string, MatchingFunc_t> mMatchingFunctionMap; ///< MFT-MCH Matching function

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
  Configurable<float> fMatchingChi2MftMchUp{"cfgMatchingChi2MftMchUp", 50.f, ""};

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

  ///    Variables for histograms configuration
  Configurable<int> fNCandidatesMax{"nCandidatesMax", 5, ""};

  int mRunNumber{0};                               // needed to detect if the run changed and trigger update of magnetic field

  Service<o2::ccdb::BasicCCDBManager> ccdbManager;
  o2::ccdb::CcdbApi ccdbApi;

  o2::aod::rctsel::RCTFlagsChecker rctChecker{"CBT_muon_glo", false, false, true};

  double mBzAtMftCenter{ 0 };

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

    EfficiencyPlotter(std::string path, std::string prefix, std::string title,
                      HistogramRegistry& registry,
                      std::unordered_map<std::string, o2::framework::HistPtr> histograms)
    {
      AxisSpec pAxis = {1000, 0, 100, "p (GeV/c)"};
      AxisSpec pTAxis = {100, 0, 10, "p_{T} (GeV/c)"};
      AxisSpec etaAxis = {100, -4, -2, "#eta"};
      AxisSpec phiAxis = {90, -180, 180, "#phi (degrees)"};

      std::string histName;
      std::string histTitle;

      histName = path + prefix + "VsPt_num";
      histTitle = title + " vs. p_{T}";
      histograms[histName] = registry.add(histName.c_str(), "good matches vs. p_{T}", {HistType::kTH1F, {pTAxis}});

    }
  };

  std::array<double, 3> zRefPlane{
      firstMFTPlaneZ,
      lastMFTPlaneZ,
      firstMCHPlaneZ
  };
  std::vector<std::pair<std::string, double>> referencePlanes{
      {"MFT-begin", 10.0},
      {"MFT-end", 15.0},
      {"MCH-begin", 100.0}
  };
  std::array<std::string, 4> quadrants{"Q0", "Q1", "Q2", "Q3"};

  double mMatchingPlaneZ{ lastMFTPlaneZ };

  // vector of all MFT-MCH(-MID) matching candidates associated to the same MCH(-MID) track,
  // to be sorted in descending order with respect to the matching quality
  // the map key is the MCH(-MID) track global index
  using MatchingCandidates = std::map<uint64_t, std::vector<uint64_t>>;

  struct CollisionInfo
  {
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
    MatchingCandidates globalMuonTracks;
  };

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

      //std::cout << "[TOTO1] " << (int)muonTrack.globalIndex() << "  " << (int)muonTrack.trackType() << "  " << muonTrack.chi2MatchMCHMFT() << std::endl;

      // global muon tracks (MFT-MCH or MFT-MCH-MID)
      uint64_t muonTrackIndex = muonTrack.globalIndex();
      auto const& mchTrack = muonTrack.template matchMCHTrack_as<TMUON>();
      uint64_t mchTrackIndex = mchTrack.globalIndex();
      //std::cout << "[TOTO1] MCH " << (int)mchTrack.globalIndex() << "  " << (int)mchTrack.trackType() << std::endl;

      // check if a vector of global muon candidates is already available for the current MCH index
      // if not, initialize a new one and add the current global muon track
      //bool globalMuonTrackFound = false;
      auto matchingCandidateIterator = matchingCandidates.find(mchTrackIndex);
      if (matchingCandidateIterator != matchingCandidates.end()) {
        matchingCandidateIterator->second.push_back(muonTrackIndex);
        //globalMuonTrackFound = true;
      } else {
        matchingCandidates[mchTrackIndex].push_back(muonTrackIndex);
      }
    }

    // sort the vectors of matching candidates in ascending order based on the user-provided comparator function
    /*auto comparatorChi2 = [&muonTracks](uint64_t trackIndex1, uint64_t trackIndex2) -> bool {
      auto const& track1 = muonTracks.rawIteratorAt(trackIndex1);
      auto const& track2 = muonTracks.rawIteratorAt(trackIndex2);

      return (track1.chi2MatchMCHMFT() < track2.chi2MatchMCHMFT());
    };*/
    for (auto& [mchIndex, candidates] : matchingCandidates) {
      //std::sort(candidates.begin(), candidates.end(), comparatorChi2);
      std::sort(candidates.begin(), candidates.end(),
          [&muonTracks, &comparator](uint64_t trackIndex1, uint64_t trackIndex2) -> bool {
                auto const& track1 = muonTracks.rawIteratorAt(trackIndex1);
                auto const& track2 = muonTracks.rawIteratorAt(trackIndex2);

                return comparator(track1, track2);
      });
    }
  }

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

  template<class EVT, class BC, class T>
  void InitCollisions(EVT const& collisions,
                      BC const& bcs,
                      T const& muonTracks,
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
        //bool globalMuonTrackFound = false;
        auto matchingCandidateIterator = collisionInfo.globalMuonTracks.find(mchTrackIndex);
        if (matchingCandidateIterator != collisionInfo.globalMuonTracks.end()) {
          matchingCandidateIterator->second.push_back(muonTrackIndex);
          //globalMuonTrackFound = true;
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
    AxisSpec chi2Axis = {1000, 0, 1000, "chi^{2}"};
    AxisSpec pAxis = {1000, 0, 100, "p (GeV/c)"};
    AxisSpec pTAxis = {100, 0, 10, "p_{T} (GeV/c)"};
    AxisSpec etaAxis = {100, -4, -2, "#eta"};
    AxisSpec phiAxis = {90, -180, 180, "#phi (degrees)"};
    std::string histPath = "matching/MC/";

    std::string histName;

    histName = "chi2_true";
    matchingHistos[histName] = registryMatching.add((histPath + histName).c_str(), "chi^{2} - true matches", {HistType::kTH1F, {chi2Axis}});
    histName = "chi2VsP_true";
    matchingHistos[histName] = registryMatching.add((histPath + histName).c_str(), "chi^{2} vs. MCH momentum - true matches", {HistType::kTH2F, {pAxis, chi2Axis}});

    histName = "chi2_fake";
    matchingHistos[histName] = registryMatching.add((histPath + histName).c_str(), "chi^{2} - fake matches", {HistType::kTH1F, {chi2Axis}});
    histName = "chi2VsP_fake";
    matchingHistos[histName] = registryMatching.add((histPath + histName).c_str(), "chi^{2} vs. MCH momentum - fake matches", {HistType::kTH2F, {pAxis, chi2Axis}});

    // pairing purity plots
    // reconstructed tracks
    histName = "goodMatchesVsPt_reco";
    matchingHistos[histName] = registryMatching.add((histPath + histName).c_str(), "good matches vs. p_{T}", {HistType::kTH1F, {pTAxis}});
    histName = "goodMatchesVsP_reco";
    matchingHistos[histName] = registryMatching.add((histPath + histName).c_str(), "good matches vs. p", {HistType::kTH1F, {pAxis}});
    histName = "goodMatchesVsEta_reco";
    matchingHistos[histName] = registryMatching.add((histPath + histName).c_str(), "good matches vs. #eta", {HistType::kTH1F, {etaAxis}});
    histName = "goodMatchesVsPhi_reco";
    matchingHistos[histName] = registryMatching.add((histPath + histName).c_str(), "good matches vs. #phi", {HistType::kTH1F, {phiAxis}});
    // MC true tracks
    histName = "goodMatchesVsPt_true";
    matchingHistos[histName] = registryMatching.add((histPath + histName).c_str(), "good matches vs. p_{T} - true", {HistType::kTH1F, {pTAxis}});
    histName = "goodMatchesVsP_true";
    matchingHistos[histName] = registryMatching.add((histPath + histName).c_str(), "good matches vs. p - true", {HistType::kTH1F, {pAxis}});
    histName = "goodMatchesVsEta_true";
    matchingHistos[histName] = registryMatching.add((histPath + histName).c_str(), "good matches vs. #eta - true", {HistType::kTH1F, {etaAxis}});
    histName = "goodMatchesVsPhi_true";
    matchingHistos[histName] = registryMatching.add((histPath + histName).c_str(), "good matches vs. #phi - true", {HistType::kTH1F, {phiAxis}});
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

    int nTrackTypes = static_cast<int>(o2::aod::fwdtrack::ForwardTrackTypeEnum::MCHStandaloneTrack) + 1;
    AxisSpec trackTypeAxis = {static_cast<int>(nTrackTypes), 0.0, static_cast<double>(nTrackTypes), "track type"};
    registry.add("nTracksPerType", "Number of tracks per type", {HistType::kTH1F, {trackTypeAxis}});

    createMatchingHistosMC();

    // Define built-in matching functions
    //________________________________________________________________________________
    mMatchingFunctionMap["matchALL"] = [](const GlobalFwdTrack& mchTrack, const TrackParCovFwd& mftTrack) -> double {
      // Match two tracks evaluating all parameters: X,Y, phi, tanl & q/pt

      SMatrix55Sym I = ROOT::Math::SMatrixIdentity(), H_k, V_k;
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
      SMatrix55Std K_k = GlobalMuonTrackCovariances * ROOT::Math::Transpose(H_k) * invResCov;

      // Update Parameters
      r_k_kminus1 = m_k - H_k * GlobalMuonTrackParameters; // Residuals of prediction

      auto matchChi2Track = ROOT::Math::Similarity(r_k_kminus1, invResCov);

      return matchChi2Track;
    };

    //________________________________________________________________________________
    mMatchingFunctionMap["matchXYPhiTanl"] = [](const GlobalFwdTrack& mchTrack, const TrackParCovFwd& mftTrack) -> double {

    // Match two tracks evaluating positions & angles

    SMatrix55Sym I = ROOT::Math::SMatrixIdentity();
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
    SMatrix54 K_k = GlobalMuonTrackCovariances * ROOT::Math::Transpose(H_k) * invResCov;

    // Residuals of prediction
    r_k_kminus1 = m_k - H_k * GlobalMuonTrackParameters;

    auto matchChi2Track = ROOT::Math::Similarity(r_k_kminus1, invResCov);

    return matchChi2Track; };

    //________________________________________________________________________________
    mMatchingFunctionMap["matchXY"] = [](const GlobalFwdTrack& mchTrack, const TrackParCovFwd& mftTrack) -> double {

    // Calculate Matching Chi2 - X and Y positions

    SMatrix55Sym I = ROOT::Math::SMatrixIdentity();
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
    SMatrix52 K_k = GlobalMuonTrackCovariances * ROOT::Math::Transpose(H_k) * invResCov;

    // Residuals of prediction
    r_k_kminus1 = m_k - H_k * GlobalMuonTrackParameters;
    auto matchChi2Track = ROOT::Math::Similarity(r_k_kminus1, invResCov);

    return matchChi2Track; };

    //________________________________________________________________________________
    mMatchingFunctionMap["matchNeedsName"] = [this](const GlobalFwdTrack& mchTrack, const TrackParCovFwd& mftTrack) -> double {

    //Hiroshima's Matching function needs a physics-based name

    //Matching constants
    Double_t LAbs = 415.;        //Absorber Length[cm]
    Double_t mumass = 0.106;     //mass of muon [GeV/c^2]
    Double_t absorberPathLength; //the length that extrapolated MCHtrack passes through absorber

    if (mMatchingPlaneZ >= -90.0) {
      absorberPathLength = LAbs;
    } else {
      absorberPathLength = 505.0 + mMatchingPlaneZ;
    }

    //defference between MFTtrack and MCHtrack

    auto dx = mftTrack.getX() - mchTrack.getX();
    auto dy = mftTrack.getY() - mchTrack.getY();
    auto dthetax = TMath::ATan(mftTrack.getPx() / TMath::Abs(mftTrack.getPz())) - TMath::ATan(mchTrack.getPx() / TMath::Abs(mchTrack.getPz()));
    auto dthetay = TMath::ATan(mftTrack.getPy() / TMath::Abs(mftTrack.getPz())) - TMath::ATan(mchTrack.getPy() / TMath::Abs(mchTrack.getPz()));

    //Multiple Scattering(=MS)

    auto pMCH = mchTrack.getP();
    auto lorentzbeta = pMCH / TMath::Sqrt(mumass * mumass + pMCH * pMCH);
    auto zMS = copysign(1.0, mchTrack.getCharge());
    auto thetaMS = 13.6 / (1000.0 * pMCH * lorentzbeta * 1.0) * zMS * TMath::Sqrt(60.0 * absorberPathLength / LAbs) * (1.0 + 0.038 * TMath::Log(60.0 * absorberPathLength / LAbs));
    auto xMS = thetaMS * absorberPathLength / TMath::Sqrt(3.0);

    //normalize by theoritical Multiple Coulomb Scattering width to be momentum-independent
    //make the dx and dtheta dimensionless

    auto dxnorm = dx / xMS;
    auto dynorm = dy / xMS;
    auto dthetaxnorm = dthetax / thetaMS;
    auto dthetaynorm = dthetay / thetaMS;

    //rotate distribution

    auto dxrot = dxnorm * TMath::Cos(TMath::Pi() / 4.0) - dthetaxnorm * TMath::Sin(TMath::Pi() / 4.0);
    auto dthetaxrot = dxnorm * TMath::Sin(TMath::Pi() / 4.0) + dthetaxnorm * TMath::Cos(TMath::Pi() / 4.0);
    auto dyrot = dynorm * TMath::Cos(TMath::Pi() / 4.0) - dthetaynorm * TMath::Sin(TMath::Pi() / 4.0);
    auto dthetayrot = dynorm * TMath::Sin(TMath::Pi() / 4.0) + dthetaynorm * TMath::Cos(TMath::Pi() / 4.0);

    //convert ellipse to circle

    auto k = 0.7; //need to optimize!!
    auto dxcircle = dxrot;
    auto dycircle = dyrot;
    auto dthetaxcircle = dthetaxrot / k;
    auto dthetaycircle = dthetayrot / k;

    //score

    auto scoreX = TMath::Sqrt(dxcircle * dxcircle + dthetaxcircle * dthetaxcircle);
    auto scoreY = TMath::Sqrt(dycircle * dycircle + dthetaycircle * dthetaycircle);
    auto score = TMath::Sqrt(scoreX * scoreX + scoreY * scoreY);

    return score; };
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

  template<class T>
  int GetQuadrant(const T& track)
  {
    double phi = track.phi() * 180 / TMath::Pi();
    return GetQuadrant(phi);
  }

  template<class T, class C>
  double GetDCA(const T& track, const C& collision)
  {
    // propagate muon track to DCA
    auto trackAtDCA = VarManager::PropagateMuon(track, collision, toDCA);
    // Calculate DCA quantities (preferable to do it with VarManager)
    double dcax = trackAtDCA.getX() - collision.posX();
    double dcay = trackAtDCA.getY() - collision.posY();
    return std::sqrt(dcax * dcax + dcay * dcay);
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
    auto mchTrackAtVertex = VarManager::PropagateMuon(mchTrack, collision, toVertex);

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
  bool IsGoodMuon(const T& muonTrack, const C& collision,
                  double chi2Cut,
                  double pCut,
                  double pTCut,
                  std::array<double, 2> etaCut,
                  std::array<double, 2> rAbsCut,
                  double nSigmaPdcaCut)
  {
    auto const& mchTrack = (static_cast<int>(muonTrack.trackType()) <= 2) ?
        muonTrack.template matchMCHTrack_as<MyMuonsWithCov>() :
        muonTrack;

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

  template<class TMUON, class TMFT>
  bool IsGoodGlobalMatching(const TMUON& muonTrack,
                            const TMFT& mftTrack,
                            double chi2CutMFT,
                            int nClustersCutMFT,
                            double matchingChi2Cut)
  {
    if (static_cast<int>(muonTrack.trackType()) >= 2)
      return false;

    //auto const& mftTrack = muonTrack.template matchMFTTrack_as<MyMFTsMC>();

    if (!IsGoodMFT(mftTrack,
                   chi2CutMFT,
                   nClustersCutMFT))
      return false;

    // MFT-MCH matching chi2 cut
    if (muonTrack.chi2MatchMCHMFT() > matchingChi2Cut)
      return false;

    return true;
  }

  template<class TMUON, class TMFT>
  bool IsGoodGlobalMatching(const TMUON& muonTrack, const TMFT& mftTrack)
  {
    return IsGoodGlobalMatching(muonTrack, mftTrack, fTrackChi2MftUp, fTrackNClustMftLow, fMatchingChi2MftMchUp);
  }

  template<class C>
  bool IsSameEvent(const C& c1, const C& c2)
  {
    return (c1.bc == c2.bc);
  }

  template<class C>
  bool IsMixedEvent(const C& c1, const C& c2)
  {
    if (IsSameEvent(c1, c2))
      return false;

    uint64_t bcDiff = (c2.bc > c1.bc) ? (c2.bc - c1.bc) : (c2.bc - c1.bc);
    // in the event mixing case, we require a minimum BC gap between the collisions
    if (bcDiff < static_cast<uint64_t>(fEventMinDeltaBc))
      return false;
    // we also require that the collisions have similar Z positions and multiplicity of MFT tracks
    if (std::fabs(c2.zVertex - c1.zVertex) > fEventMaxDeltaVtxZ)
      return false;
    if (std::abs(c2.mftTracksMultiplicity - c1.mftTracksMultiplicity) > fEventMaxDeltaNMFT)
      return false;

    return true;
  }

  template <class TMUON>
  bool IsTrueGlobalMatching(const TMUON& muonTrack, const std::vector<std::pair<uint64_t, uint64_t>>& matchablePairs)
  {
    if (static_cast<int>(muonTrack.trackType()) >= 2)
      return false;

    //auto const& mchTrack = muonTrack.template matchMCHTrack_as<TMUON>();
    //auto const& mftTrack = muonTrack.template matchMFTTrack_as<TMFT>();
    //uint64_t mchTrackId = static_cast<uint64_t>(mchTrack.globalIndex());
    uint64_t mchTrackId = static_cast<uint64_t>(muonTrack.matchMCHTrackId());
    uint64_t mftTrackId = static_cast<uint64_t>(muonTrack.matchMFTTrackId());

    std::pair<uint64_t, uint64_t> trackIndexes = std::make_pair(mchTrackId, mftTrackId);

    return (std::find(matchablePairs.begin(), matchablePairs.end(), trackIndexes) != matchablePairs.end());
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

    //std::cout << std::format("[TOTO] P1=({:0.2f} [{:0.2f},{:0.2f},{:0.2f}])  P2=({:0.2f} [{:0.2f},{:0.2f},{:0.2f}])  M={:0.2f}",
    //    track1.getP(), track1.getPx(), track1.getPy(), track1.getPz(),
    //    track2.getP(), track2.getPx(), track2.getPy(), track2.getPz(),
    //    dimuon.M()) << std::endl;

    return dimuon.M();
  }

  template<class T, class C>
  double GetMuMuInvariantMass(const T& track1, const T& track2, const C& collision)
  {
    // propagate muon tracks to vertex
    auto const& muonTrack1AtVertex = VarManager::PropagateMuon(track1, collision, toVertex);
    auto const& muonTrack2AtVertex = VarManager::PropagateMuon(track2, collision, toVertex);

    return GetMuMuInvariantMass(muonTrack1AtVertex, muonTrack2AtVertex);
  }

  template<class TMFT, class TMCH, class C>
  o2::dataformats::GlobalFwdTrack PropagateMftToVertex(const TMFT& mftTrack, const TMCH& mchTrack, const C& collision)
  {
    // propagate MCH track to the vertex to get the updated momentum
    auto const& mchTrackAtVertex = VarManager::PropagateMuon(mchTrack, collision, VarManager::kToDCA);

    double px = mchTrackAtVertex.getP() * sin(M_PI / 2 - atan(mftTrack.tgl())) * cos(mftTrack.phi());
    double py = mchTrackAtVertex.getP() * sin(M_PI / 2 - atan(mftTrack.tgl())) * sin(mftTrack.phi());
    double pt = std::sqrt(std::pow(px, 2) + std::pow(py, 2));
    double sign = mchTrack.sign();
    double signed1Pt = sign / pt;

    //std::cout << std::format("[TOTO]   P(MCH)=({:0.2f},{:0.2f})  P(MFT)=({:0.2f},{:0.2f})  Pt(scaled)={:0.2f}",
    //    mchTrackAtVertex.getP(), mchTrackAtVertex.getPt(),
    //    mftTrack.p(), mftTrack.signed1Pt(), signed1Pt) << std::endl;

    double chi2 = mftTrack.chi2();
    SMatrix5 tpars = {mftTrack.x(), mftTrack.y(), mftTrack.phi(), mftTrack.tgl(), signed1Pt};
    std::vector<double> v1{0, 0, 0, 0, 0,
      0, 0, 0, 0, 0,
      0, 0, 0, 0, 0};
    SMatrix55 tcovs(v1.begin(), v1.end());
    o2::track::TrackParCovFwd fwdtrack{mftTrack.z(), tpars, tcovs, chi2};
    o2::dataformats::GlobalFwdTrack propmuon;

    double centerMFT[3] = {0, 0, -61.4};
    o2::field::MagneticField* field = static_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField());
    auto Bz = field->getBz(centerMFT); // Get field at centre of MFT
    auto geoMan = o2::base::GeometryManager::meanMaterialBudget(mftTrack.x(), mftTrack.y(), mftTrack.z(), collision.posX(), collision.posY(), collision.posZ());
    auto x2x0 = static_cast<float>(geoMan.meanX2X0);
    fwdtrack.propagateToVtxhelixWithMCS(collision.posZ(), {collision.posX(), collision.posY()}, {collision.covXX(), collision.covYY()}, Bz, x2x0);
    propmuon.setParameters(fwdtrack.getParameters());
    propmuon.setZ(fwdtrack.getZ());
    propmuon.setCovariances(fwdtrack.getCovariances());

    return propmuon;
  }

  template<class TMFT, class TMCH, class C>
  double GetMuMuInvariantMass(const TMFT& mftTrack1, const TMCH& mchTrack1, const TMFT& mftTrack2, const TMCH& mchTrack2, const C& collision)
  {
    auto mftTrack1AtVertex = PropagateMftToVertex(mftTrack1, mchTrack1, collision);
    auto mftTrack2AtVertex = PropagateMftToVertex(mftTrack2, mchTrack2, collision);

    return GetMuMuInvariantMass(mftTrack1AtVertex, mftTrack2AtVertex);
  }

  using MuonPair = std::pair<std::pair<uint64_t, uint64_t>, std::pair<uint64_t, uint64_t>>;
  using GlobalMuonPair = std::pair<std::pair<uint64_t, std::vector<uint64_t>>, std::pair<uint64_t, std::vector<uint64_t>>>;

  void GetMuonPairs(const std::map<uint64_t, CollisionInfo>& collisionInfos,
                    std::vector<MuonPair>& muonPairs,
                    std::vector<GlobalMuonPair>& globalMuonPairs)
  {
    // muon tracks - outer loop over collisions
    for (auto& [collisionIndex1, collisionInfo1] : collisionInfos) {

      // outer loop over muon tracks
      for (auto mchIndex1 : collisionInfo1.mchTracks) {

        // inner loop over collisions
        for (auto& [collisionIndex2, collisionInfo2] : collisionInfos) {
          // avoid double-counting of collisions
          if (collisionIndex2 < collisionIndex1) continue;

          bool sameEvent = (collisionIndex1 == collisionIndex2);
          bool mixedEvent = IsMixedEvent(collisionInfo1, collisionInfo2);

          if (!sameEvent && !mixedEvent)
            continue;

          // inner loop over muon tracks
          for (auto mchIndex2 : collisionInfo2.mchTracks) {
            // avoid double-counting of muon pairs if we are not mixing events
            if (sameEvent && mchIndex2 <= mchIndex1) continue;

            MuonPair muonPair{{collisionIndex1, mchIndex1}, {collisionIndex2, mchIndex2}};
            muonPairs.emplace_back(muonPair);
          }
        }
      }
    }

    // global muon tracks - outer loop over collisions
    for (auto& [collisionIndex1, collisionInfo1] : collisionInfos) {

      // outer loop over global muon tracks
      for (auto& [mchIndex1, globalTracksVector1] : collisionInfo1.globalMuonTracks) {

        // inner loop over collisions
        for (auto& [collisionIndex2, collisionInfo2] : collisionInfos) {
          // avoid double-counting of collisions
          if (collisionIndex2 < collisionIndex1) continue;

          bool sameEvent = (collisionIndex1 == collisionIndex2);
          bool mixedEvent = IsMixedEvent(collisionInfo1, collisionInfo2);

          if (!sameEvent && !mixedEvent)
            continue;

          // outer loop over global muon tracks
          for (auto& [mchIndex2, globalTracksVector2] : collisionInfo2.globalMuonTracks) {
            // avoid double-counting of muon pairs if we are not mixing events
            if (sameEvent && mchIndex2 <= mchIndex1) continue;

            GlobalMuonPair muonPair{{collisionIndex1, globalTracksVector1}, {collisionIndex2, globalTracksVector2}};
            globalMuonPairs.emplace_back(muonPair);
          }
        }
      }
    }
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
    track.setParameters(tpars);
    track.setZ(fwdtrack.getZ());
    track.setCovariances(tcovs);
    auto mchTrack = sExtrap.FwdtoMCH(track);

    o2::mch::TrackExtrap::extrapToZ(mchTrack, z);

    auto proptrack = sExtrap.MCHtoFwd(mchTrack);
    o2::dataformats::GlobalFwdTrack propmuon;
    propmuon.setParameters(proptrack.getParameters());
    propmuon.setZ(proptrack.getZ());
    propmuon.setCovariances(proptrack.getCovariances());

    return propmuon;
  }

  template <typename TMFT>
  o2::dataformats::GlobalFwdTrack PropagateToZMFT(const TMFT& mftTrack, const double pMCH, int signMCH, const double z)
  {
    double px = pMCH * sin(M_PI / 2 - atan(mftTrack.tgl())) * cos(mftTrack.phi());
    double py = pMCH * sin(M_PI / 2 - atan(mftTrack.tgl())) * sin(mftTrack.phi());
    double pt = std::sqrt(std::pow(px, 2) + std::pow(py, 2));
    double sign = signMCH;

    SMatrix5 tpars = {mftTrack.x(), mftTrack.y(), mftTrack.phi(), mftTrack.tgl(), sign / pt};
    std::vector<double> v1{0, 0, 0, 0, 0,
      0, 0, 0, 0, 0,
      0, 0, 0, 0, 0};
    SMatrix55 tcovs(v1.begin(), v1.end());

    o2::dataformats::GlobalFwdTrack track;
    track.setParameters(tpars);
    track.setZ(mftTrack.z());
    track.setCovariances(tcovs);

    auto mchTrackExt = sExtrap.FwdtoMCH(track);

    o2::mch::TrackExtrap::extrapToZ(mchTrackExt, z);

    o2::dataformats::GlobalFwdTrack propmuon;
    auto proptrack = sExtrap.MCHtoFwd(mchTrackExt);
    propmuon.setParameters(proptrack.getParameters());
    propmuon.setZ(proptrack.getZ());
    propmuon.setCovariances(proptrack.getCovariances());

    return propmuon;
  }

  template<class TMFT, class TMCH>
  double DoMatchingStandard(const TMFT& mftTrack, MyMFTCovariance const& mftTrackCov, const TMCH& mchTrack, std::string matchingFunction, bool correctAlignment)
  {
    ///////////////////////////////////////////////////////////////////////////
    // Extrapolat MFT
    ///////////////////////////////////////////////////////////////////////////

    SMatrix5 tmftpars(mftTrack.x(),
        mftTrack.y(),
        mftTrack.phi(),
        mftTrack.tgl(),
        mftTrack.signed1Pt());

    SMatrix55Sym tmftcovs;
    tmftcovs(0, 0) = mftTrackCov.cXX();
    //std::cout << "mftTrackCov.cXX(): " << mftTrackCov.cXX() << std::endl;
    tmftcovs(0, 1) = mftTrackCov.cXY();
    tmftcovs(0, 2) = mftTrackCov.cPhiX();
    tmftcovs(0, 3) = mftTrackCov.cTglX();
    tmftcovs(0, 4) = mftTrackCov.c1PtX();

    tmftcovs(1, 1) = mftTrackCov.cYY();
    tmftcovs(1, 2) = mftTrackCov.cPhiY();
    tmftcovs(1, 3) = mftTrackCov.cTglY();
    tmftcovs(1, 4) = mftTrackCov.c1PtY();

    tmftcovs(2, 2) = mftTrackCov.cPhiPhi();
    tmftcovs(2, 3) = mftTrackCov.cTglPhi();
    tmftcovs(2, 4) = mftTrackCov.c1PtPhi();

    tmftcovs(3, 3) = mftTrackCov.cTglTgl();
    tmftcovs(3, 4) = mftTrackCov.c1PtTgl();

    tmftcovs(4, 4) = mftTrackCov.c1Pt21Pt2();

    o2::track::TrackParCovFwd extrap_mfttrack{mftTrack.z(),
      tmftpars,
      tmftcovs,
      mftTrack.chi2()};

    float zPlane = o2::mft::constants::mft::LayerZCoordinate()[9];

    extrap_mfttrack.propagateToZ(zPlane, mBzAtMftCenter); // z in cm

    o2::dataformats::GlobalFwdTrack mftTrackAtMatchingPlane;
    mftTrackAtMatchingPlane.setParameters(extrap_mfttrack.getParameters());
    mftTrackAtMatchingPlane.setZ(extrap_mfttrack.getZ());
    mftTrackAtMatchingPlane.setCovariances(extrap_mfttrack.getCovariances());

    ///////////////////////////////////////////////////////////////////////////
    // Extrapolate MCH
    ///////////////////////////////////////////////////////////////////////////

    float cov[15] = {
        mchTrack.cXX(), mchTrack.cXY(), mchTrack.cYY(),
        mchTrack.cPhiX(), mchTrack.cPhiY(), mchTrack.cPhiPhi(),
        mchTrack.cTglX(), mchTrack.cTglY(), mchTrack.cTglPhi(),
        mchTrack.cTglTgl(), mchTrack.c1PtX(), mchTrack.c1PtY(),
        mchTrack.c1PtPhi(), mchTrack.c1PtTgl(), mchTrack.c1Pt21Pt2()};

    SMatrix5 tpars(mchTrack.x(),
        mchTrack.y(),
        mchTrack.phi(),
        mchTrack.tgl(),
        mchTrack.signed1Pt());
    SMatrix55 tcovs(cov, cov + 15);
    double chi2 = mchTrack.chi2();
    //std::cout << "mchTrack.cXX(): " << mchTrack.cXX() << std::endl;

    o2::track::TrackParCovFwd parcovmuontrack{mchTrack.z(), tpars, tcovs, chi2};

    o2::dataformats::GlobalFwdTrack gtrack;
    gtrack.setParameters(tpars);
    gtrack.setZ(parcovmuontrack.getZ());
    gtrack.setCovariances(tcovs);

    auto mchtrack = sExtrap.FwdtoMCH(gtrack);
    o2::mch::TrackExtrap::extrapToVertexWithoutBranson(mchtrack,zPlane);

    auto fwdtrack = sExtrap.MCHtoFwd(mchtrack);

    o2::dataformats::GlobalFwdTrack muonTrackAtMatchingPlane;
    muonTrackAtMatchingPlane.setParameters(fwdtrack.getParameters());
    muonTrackAtMatchingPlane.setZ(fwdtrack.getZ());
    muonTrackAtMatchingPlane.setCovariances(fwdtrack.getCovariances());

    ///////////////////////////////////////////////////////////////////////////
    // CALCULATE CHI2
    ///////////////////////////////////////////////////////////////////////////

    //std::cout << "mchTrackAtMatchingPlane.cXX(): " << muonTrackAtMatchingPlane.getCovariances()(0, 0) << std::endl;
    //std::cout << "mftTrackAtMatchingPlane.cXX(): " << mftTrackAtMatchingPlane.getCovariances()(0, 0) << std::endl;
    double matchingChi2 = mMatchingFunctionMap[matchingFunction](muonTrackAtMatchingPlane, mftTrackAtMatchingPlane);
    //std::cout << std::format("Matching standard: z={:0.2f}  MCH=[{:0.2f} {:0.2f} {:0.2f}]  MFT=[{:0.2f} {:0.2f} {:0.2f}]  chi2={:0.3f}",
    //    zPlane,
    //    muonTrackAtMatchingPlane.getX(), muonTrackAtMatchingPlane.getY(), muonTrackAtMatchingPlane.getZ(),
    //    mftTrackAtMatchingPlane.getX(), mftTrackAtMatchingPlane.getY(), mftTrackAtMatchingPlane.getZ(),
    //    matchingChi2) << std::endl;
    //return mMatchingFunctionMap["matchALL"](muonTrackAtMatchingPlane, mftTrackAtMatchingPlane);
    return matchingChi2;
  }

  template<class TMFT, class TMCH>
  double DoMatchingAlt(const TMFT& mftTrack, MyMFTCovariance const& mftTrackCov, const TMCH& mchTrack, std::string matchingFunction, double zMatchingPlane, bool correctAlignment)
  {
    ///////////////////////////////////////////////////////////////////////////
    // Extrapolate MCH
    ///////////////////////////////////////////////////////////////////////////

    float cov[15] = {
        mchTrack.cXX(), mchTrack.cXY(), mchTrack.cYY(),
        mchTrack.cPhiX(), mchTrack.cPhiY(), mchTrack.cPhiPhi(),
        mchTrack.cTglX(), mchTrack.cTglY(), mchTrack.cTglPhi(),
        mchTrack.cTglTgl(), mchTrack.c1PtX(), mchTrack.c1PtY(),
        mchTrack.c1PtPhi(), mchTrack.c1PtTgl(), mchTrack.c1Pt21Pt2()};

    SMatrix5 tpars(mchTrack.x(),
        mchTrack.y(),
        mchTrack.phi(),
        mchTrack.tgl(),
        mchTrack.signed1Pt());
    SMatrix55 tcovs(cov, cov + 15);
    double chi2 = mchTrack.chi2();
    //std::cout << "mchTrack.cXX(): " << mchTrack.cXX() << std::endl;

    //o2::track::TrackParCovFwd parcovmuontrack{mchTrack.z(), tpars, tcovs, chi2};

    o2::dataformats::GlobalFwdTrack mchTrackPars;
    mchTrackPars.setParameters(tpars);
    mchTrackPars.setZ(mchTrack.z());
    mchTrackPars.setCovariances(tcovs);

    auto mchTrackParsAtMFT = sExtrap.FwdtoMCH(mchTrackPars);

    // extrapolate MCH track to first MFT plane
    o2::mch::TrackExtrap::extrapToVertexWithoutBranson(mchTrackParsAtMFT, mftTrack.z());

    // extrapolate MCH track to matching plane
    auto mchTrackParsTemp = mchTrackParsAtMFT;
    o2::mch::TrackExtrap::extrapToZCov(mchTrackParsTemp, zMatchingPlane);
    //o2::mch::TrackExtrap::extrapToVertexWithoutBranson(mchTrackParsTemp,zMatchingPlane);

    auto mchTrackParsTemp2 = sExtrap.MCHtoFwd(mchTrackParsTemp);

    o2::dataformats::GlobalFwdTrack mchTrackAtMatchingPlane;
    mchTrackAtMatchingPlane.setParameters(mchTrackParsTemp2.getParameters());
    mchTrackAtMatchingPlane.setZ(mchTrackParsTemp2.getZ());
    mchTrackAtMatchingPlane.setCovariances(mchTrackParsTemp2.getCovariances());

    ///////////////////////////////////////////////////////////////////////////
    // Extrapolat MFT
    ///////////////////////////////////////////////////////////////////////////
    double pMCH = mchTrackParsAtMFT.p();
    double signMCH = mchTrack.sign();
    double px = pMCH * sin(M_PI / 2 - atan(mftTrack.tgl())) * cos(mftTrack.phi());
    double py = pMCH * sin(M_PI / 2 - atan(mftTrack.tgl())) * sin(mftTrack.phi());
    double pt = std::sqrt(std::pow(px, 2) + std::pow(py, 2));
    double sign = signMCH;

    SMatrix5 tmftpars(mftTrack.x(),
        mftTrack.y(),
        mftTrack.phi(),
        mftTrack.tgl(),
        sign / pt);

    SMatrix55Sym tmftcovs;
    tmftcovs(0, 0) = mftTrackCov.cXX();
    //std::cout << "mftTrackCov.cXX(): " << mftTrackCov.cXX() << std::endl;
    tmftcovs(0, 1) = mftTrackCov.cXY();
    tmftcovs(0, 2) = mftTrackCov.cPhiX();
    tmftcovs(0, 3) = mftTrackCov.cTglX();
    tmftcovs(0, 4) = mftTrackCov.c1PtX();

    tmftcovs(1, 1) = mftTrackCov.cYY();
    tmftcovs(1, 2) = mftTrackCov.cPhiY();
    tmftcovs(1, 3) = mftTrackCov.cTglY();
    tmftcovs(1, 4) = mftTrackCov.c1PtY();

    tmftcovs(2, 2) = mftTrackCov.cPhiPhi();
    tmftcovs(2, 3) = mftTrackCov.cTglPhi();
    tmftcovs(2, 4) = mftTrackCov.c1PtPhi();

    tmftcovs(3, 3) = mftTrackCov.cTglTgl();
    tmftcovs(3, 4) = mftTrackCov.c1PtTgl();

    tmftcovs(4, 4) = mftTrackCov.c1Pt21Pt2();

    o2::dataformats::GlobalFwdTrack mftTrackPars;
    mftTrackPars.setParameters(tmftpars);
    mftTrackPars.setZ(mftTrack.z());
    mftTrackPars.setCovariances(tmftcovs);

    auto mftTrackParsTemp = sExtrap.FwdtoMCH(mftTrackPars);

    o2::mch::TrackExtrap::extrapToZCov(mftTrackParsTemp, zMatchingPlane);

    auto mftTrackParsTemp2 = sExtrap.MCHtoFwd(mftTrackParsTemp);

    o2::dataformats::GlobalFwdTrack mftTrackAtMatchingPlane;
    mftTrackAtMatchingPlane.setParameters(mftTrackParsTemp2.getParameters());
    mftTrackAtMatchingPlane.setZ(mftTrackParsTemp2.getZ());
    mftTrackAtMatchingPlane.setCovariances(mftTrackParsTemp2.getCovariances());

    ///////////////////////////////////////////////////////////////////////////
    // CALCULATE CHI2
    ///////////////////////////////////////////////////////////////////////////

    //std::cout << "mchTrackAtMatchingPlane.cXX(): " << mchTrackAtMatchingPlane.getCovariances()(0, 0) << std::endl;
    //std::cout << "mftTrackAtMatchingPlane.cXX(): " << mftTrackAtMatchingPlane.getCovariances()(0, 0) << std::endl;
    double matchingChi2 = mMatchingFunctionMap[matchingFunction](mchTrackAtMatchingPlane, mftTrackAtMatchingPlane);
    //std::cout << std::format("Matching alt:      z={:0.2f}  MCH=[{:0.2f} {:0.2f} {:0.2f}]  MFT=[{:0.2f} {:0.2f} {:0.2f}]  chi2={:0.3f}",
    //    zMatchingPlane,
    //    mchTrackAtMatchingPlane.getX(), mchTrackAtMatchingPlane.getY(), mchTrackAtMatchingPlane.getZ(),
    //    mftTrackAtMatchingPlane.getX(), mftTrackAtMatchingPlane.getY(), mftTrackAtMatchingPlane.getZ(),
    //    matchingChi2) << std::endl;
    return matchingChi2;
  }

  template <class TMUON, class TMFT>
  void FillMatchingPlotsMC(TMUON const& muonTracks,
      //MyMuonsWithCov const& muonTracks,
      TMFT const& mftTracks)
  {
    auto comparatorChi2 = [&](const MyMuonMC& track1, const MyMuonMC& track2) -> bool {
      return (track1.chi2MatchMCHMFT() < track2.chi2MatchMCHMFT());
    };

    MatchingCandidates matchingCandidates;
    GetMatchingCandidates(muonTracks, mftTracks, comparatorChi2, matchingCandidates);

    std::vector<std::pair<uint64_t, uint64_t>> matchablePairs;
    GetMatchablePairs(muonTracks, mftTracks, matchablePairs);

    for (auto [mchIndex, globalTracksVector] : matchingCandidates) {
      if (globalTracksVector.size() < 1) continue;
      auto const& muonTrack = muonTracks.rawIteratorAt(globalTracksVector[0]);
      if (static_cast<int>(muonTrack.trackType()) >= 2) continue;

      //std::cout << "[TOTO2] " << (int)muonTrack.globalIndex() << "  " << (int)muonTrack.trackType() << "  " << muonTrack.chi2MatchMCHMFT() << std::endl;

      //auto const& mchTrack = muonTracks.rawIteratorAt(mchIndex);
      auto const& mchTrack = muonTrack.template matchMCHTrack_as<TMUON>();
      //std::cout << "[TOTO2] MCH " << (int)mchTrack.globalIndex() << "  " << (int)mchTrack.trackType() << std::endl;

      double chi2Prod = muonTrack.chi2MatchMCHMFT();
      double mchMom = mchTrack.p();
      double mchPt = mchTrack.pt();

      //continue;

      auto const& mftTrack = muonTrack.template matchMFTTrack_as<TMFT>();
      bool isGoodMatch = IsGoodGlobalMatching(muonTrack, mftTrack);
      bool isGoodMatchNoChi2 = IsGoodGlobalMatching(muonTrack, mftTrack, fTrackChi2MftUp, fTrackNClustMftLow, 1.0e6);

      bool isTrueMatch = IsTrueGlobalMatching(muonTrack, matchablePairs);

      // matching chi2 analysis
      if (isGoodMatchNoChi2) {
        if (isTrueMatch) {
          std::get<std::shared_ptr<TH1>>(matchingHistos["chi2_true"])->Fill(chi2Prod);
          std::get<std::shared_ptr<TH2>>(matchingHistos["chi2VsP_true"])->Fill(mchMom, chi2Prod);
        } else {
          std::get<std::shared_ptr<TH1>>(matchingHistos["chi2_fake"])->Fill(chi2Prod);
          std::get<std::shared_ptr<TH2>>(matchingHistos["chi2VsP_fake"])->Fill(mchMom, chi2Prod);
        }
      }

      // pairing purity analysis
      if (isGoodMatch) {
        std::get<std::shared_ptr<TH1>>(matchingHistos["goodMatchesVsPt_reco"])->Fill(mchPt);
        if (isTrueMatch) {
          std::get<std::shared_ptr<TH1>>(matchingHistos["goodMatchesVsPt_true"])->Fill(mchPt);
        }
      }
    }
  }

  void checkCollisions(MyEvents const& collisions,
                      aod::BCsWithTimestamps const& bcs)
  {
    for (auto& coll : collisions) {
      rctChecker(coll);
    }
  }

  PROCESS_SWITCH(qaMatching, checkCollisions, "check collisions", false);

  void processQAMC(MyEvents const& collisions,
      aod::BCsWithTimestamps const& bcs,
      MyMuonsFullMC const& muonTracks,
      //MyMuonsWithCov const& muonTracks,
      MyMFTsMC const& mftTracks,
      //MyMFTCovariances const& mftCovariances,
      //aod::FwdTrkCls const& clusters,
      aod::McParticles const& mcParticles)
  {
    auto bc = bcs.begin();
    initCCDB(bc);

    FillMatchingPlotsMC(muonTracks, mftTracks);

    //std::map<uint64_t, CollisionInfo> collisionInfos;
    //InitCollisions(collisions, bcs, muonTracks, mftTracks, collisionInfos);
  }

  PROCESS_SWITCH(qaMatching, processQAMC, "process qa MC", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<qaMatching>(cfgc)};
};
