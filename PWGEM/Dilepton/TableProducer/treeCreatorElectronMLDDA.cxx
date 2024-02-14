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
// ========================
//
// This code will create data table for inputs to machine learning for electrons.
//    Please write to: daiki.sekihata@cern.ch

#include <vector>
#include <array>
#include <random>
#include "Math/Vector4D.h"
#include "DCAFitter/DCAFitterN.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/trackUtilities.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "CommonConstants/PhysicsConstants.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"
#include "PWGEM/Dilepton/DataModel/lmeeMLTables.h"
#include "PWGEM/PhotonMeson/Utils/PCMUtilities.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::constants::physics;

using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults>;
using MyCollision = MyCollisions::iterator;

using MyTracks = soa::Join<aod::Tracks, aod::TracksCov, aod::TracksExtra, aod::TracksDCA,
                           aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
                           aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta>;
using MyTrack = MyTracks::iterator;

struct TreeCreatorElectronMLDDA {
  enum class EM_V0_Label : int { // Reconstructed V0
    kUndef = -1,
    kGamma = 0,
    kK0S = 1,
    kLambda = 2,
    kAntiLambda = 3,
  };

  SliceCache cache;
  Produces<o2::aod::EMPrimaryTracks> emprimarytracks; // flat table containing collision + track information

  // Basic checks
  HistogramRegistry registry{
    "registry",
    {
      {"hEventCounter", "hEventCounter", {HistType::kTH1F, {{5, 0.5f, 5.5f}}}},
      {"hAP", "Armenteros Podolanski", {HistType::kTH2F, {{200, -1.f, +1.f}, {250, 0, 0.25}}}},
      {"hCosPA", "V0 cosine of pointing angle", {HistType::kTH1F, {{100, 0.99, 1.f}}}},
      {"hPCA", "V0 distance between 2 legs", {HistType::kTH1F, {{200, 0.f, 2.f}}}},
      {"hMassGamma", "V0 mass gamma", {HistType::kTH1F, {{100, 0, 0.1}}}},
      {"hMassK0Short", "V0 mass K0S", {HistType::kTH1F, {{200, 0.4, 0.6}}}},
      {"hMassLambda", "V0 mass Lambda", {HistType::kTH1F, {{100, 1.05, 1.15}}}},
      {"hMassAntiLambda", "V0 mass AntiLambda", {HistType::kTH1F, {{100, 1.05, 1.15}}}},
      {"hMvsPhiV", "mee vs. phiv", {HistType::kTH2F, {{72, 0, M_PI}, {100, 0, 0.1}}}},
    },
  };

  // Configurables

  // CCDB options
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<bool> skipGRPOquery{"skipGRPOquery", true, "skip grpo query"};
  // Operation and minimisation criteria
  Configurable<double> d_bz_input{"d_bz", -999, "bz field, -999 is automatic"};
  Configurable<bool> d_UseAbsDCA{"d_UseAbsDCA", true, "Use Abs DCAs"};
  Configurable<bool> d_UseWeightedPCA{"d_UseWeightedPCA", false, "Vertices use cov matrices"};
  Configurable<int> useMatCorrType{"useMatCorrType", 0, "0: none, 1: TGeo, 2: LUT"};

  // common variables
  Configurable<float> minpt{"minpt", 0.05, "min pt for global tracks"};
  Configurable<float> maxeta{"maxeta", 0.9, "max. eta for global tracks"};
  Configurable<float> maxdcaXY{"maxdcaXY", 0.5, "max. DCA in XY"};
  Configurable<float> maxdcaZ{"maxdcaZ", 0.5, "max. DCA in Z"};

  // TPC-related variables
  Configurable<int> mincrossedrows{"mincrossedrows", 80, "min. crossed rows"};
  Configurable<float> maxchi2tpc{"maxchi2tpc", 5.0, "max. chi2/NclsTPC"};

  // ITS-related variables
  Configurable<int> minitsncls{"minitsncls", 4, "min. number of ITS clusters"};
  Configurable<float> maxchi2its{"maxchi2its", 6.0, "max. chi2/NclsITS"};

  // for V0s
  Configurable<float> minv0cospa{"minv0cospa", 0.997, "minimum V0 CosPA"};
  Configurable<float> maxdcav0dau{"maxdcav0dau", 0.5, "max distance between V0 Daughters"};
  Configurable<float> downscaling_pion{"downscaling_pion", 0.05, "down scaling factor to store pion"};

  // for phiv vs. mee
  Configurable<float> minTPCNsigmaEl{"minTPCNsigmaEl", -4.0, "min. TPC n sigma for electron inclusion"};
  Configurable<float> maxTPCNsigmaEl{"maxTPCNsigmaEl", 4.0, "max. TPC n sigma for electron inclusion"};
  Configurable<float> slope{"slope", 0.0185, "slope for m vs. phiv"};
  Configurable<float> intercept{"intercept", -0.0280, "intercept for m vs. phiv"};

  int mRunNumber;
  float d_bz;
  float maxSnp;  // max sine phi for propagation
  float maxStep; // max step size (cm) for propagation
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::base::MatLayerCylSet* lut = nullptr;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;
  o2::vertexing::DCAFitterN<2> fitter;

  std::mt19937 engine;
  std::uniform_real_distribution<float> dist01;

  void init(InitContext& context)
  {
    mRunNumber = 0;
    d_bz = 0;
    maxSnp = 0.85f;  // could be changed later
    maxStep = 2.00f; // could be changed later

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    if (useMatCorrType == 1) {
      LOGF(info, "TGeo correction requested, loading geometry");
      if (!o2::base::GeometryManager::isGeometryLoaded()) {
        ccdb->get<TGeoManager>(geoPath);
      }
    }
    if (useMatCorrType == 2) {
      LOGF(info, "LUT correction requested, loading LUT");
      lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(lutPath));
    }

    fitter.setPropagateToPCA(true);
    fitter.setMaxR(200.);
    fitter.setMinParamChange(1e-3);
    fitter.setMinRelChi2Change(0.9);
    fitter.setMaxDZIni(1e9);
    fitter.setMaxChi2(1e9);
    fitter.setUseAbsDCA(d_UseAbsDCA);
    fitter.setWeightedFinalPCA(d_UseWeightedPCA);

    if (useMatCorrType == 1) {
      matCorr = o2::base::Propagator::MatCorrType::USEMatCorrTGeo;
    } else if (useMatCorrType == 2) {
      matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
    }
    fitter.setMatCorrType(matCorr);

    std::random_device seed_gen;
    engine = std::mt19937(seed_gen());
    dist01 = std::uniform_real_distribution<float>(0.0f, 1.0f);
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }

    // In case override, don't proceed, please - no CCDB access required
    if (d_bz_input > -990) {
      d_bz = d_bz_input;
      fitter.setBz(d_bz);
      o2::parameters::GRPMagField grpmag;
      if (fabs(d_bz) > 1e-5) {
        grpmag.setL3Current(30000.f / (d_bz / 5.0f));
      }
      o2::base::Propagator::initFieldFromGRP(&grpmag);
      mRunNumber = bc.runNumber();
      return;
    }

    auto run3grp_timestamp = bc.timestamp();
    o2::parameters::GRPObject* grpo = 0x0;
    o2::parameters::GRPMagField* grpmag = 0x0;
    if (!skipGRPOquery)
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(grpPath, run3grp_timestamp);
    if (grpo) {
      o2::base::Propagator::initFieldFromGRP(grpo);
      // Fetch magnetic field from ccdb for current collision
      d_bz = grpo->getNominalL3Field();
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    } else {
      grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, run3grp_timestamp);
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path " << grpmagPath << " of object GRPMagField and " << grpPath << " of object GRPObject for timestamp " << run3grp_timestamp;
      }
      o2::base::Propagator::initFieldFromGRP(grpmag);
      // Fetch magnetic field from ccdb for current collision
      d_bz = std::lround(5.f * grpmag->getL3Current() / 30000.f);
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    }
    mRunNumber = bc.runNumber();
    // Set magnetic field value once known
    fitter.setBz(d_bz);

    if (useMatCorrType == 2) {
      // setMatLUT only after magfield has been initalized (setMatLUT has implicit and problematic init field call if not)
      o2::base::Propagator::Instance()->setMatLUT(lut);
    }
  }

  EM_V0_Label checkV0(const float alpha, const float qt)
  {
    // Gamma cuts
    const float cutAlphaG = 0.95;
    const float cutQTG = 0.01;

    // K0S cuts
    const float cutQTK0S[2] = {0.1075, 0.215};
    const float cutAPK0S[2] = {0.199, 0.8}; // parameters for curved QT cut

    // Lambda & A-Lambda cuts
    const float cutQTL = 0.03;
    const float cutAlphaL[2] = {0.35, 0.7};
    const float cutAlphaAL[2] = {-0.7, -0.35};
    const float cutAPL[3] = {0.107, -0.69, 0.5}; // parameters for curved QT cut

    // Check for Gamma candidates
    if (std::pow(alpha / cutAlphaG, 2) + std::pow(qt / cutQTG, 2) < std::pow(1.f, 2)) {
      return EM_V0_Label::kGamma;
    }

    // Check for K0S candidates
    float q = cutAPK0S[0] * sqrt(abs(1 - alpha * alpha / (cutAPK0S[1] * cutAPK0S[1])));
    if ((qt > cutQTK0S[0]) && (qt < cutQTK0S[1]) && (qt > q)) {
      return EM_V0_Label::kK0S;
    }

    // Check for Lambda candidates
    q = cutAPL[0] * sqrt(abs(1 - ((alpha + cutAPL[1]) * (alpha + cutAPL[1])) / (cutAPL[2] * cutAPL[2])));
    if ((alpha > cutAlphaL[0]) && (alpha < cutAlphaL[1]) && (qt > cutQTL) && (qt < q)) {
      return EM_V0_Label::kLambda;
    }

    // Check for AntiLambda candidates
    q = cutAPL[0] * sqrt(abs(1 - ((alpha - cutAPL[1]) * (alpha - cutAPL[1])) / (cutAPL[2] * cutAPL[2])));
    if ((alpha > cutAlphaAL[0]) && (alpha < cutAlphaAL[1]) && (qt > cutQTL) && (qt < q)) {
      return EM_V0_Label::kAntiLambda;
    }

    return EM_V0_Label::kUndef;
  }

  template <typename TTrack>
  bool IsSelected(TTrack const& track)
  {
    if (!track.hasITS() || !track.hasTPC()) {
      return false;
    }

    if (track.itsNCls() < minitsncls) {
      return false;
    }
    if (track.itsChi2NCl() < -1) { // if tracks are not reconstructed properly, chi2/ITSncls is set to -999;
      return false;
    }

    if (track.tpcNClsCrossedRows() < mincrossedrows) {
      return false;
    }
    return true;
  }

  template <typename TCollision, typename TTrack>
  void fillTrackTable(TCollision const& collision, TTrack const& track, const int pidlabel)
  {
    if (std::find(stored_trackIds.begin(), stored_trackIds.end(), track.globalIndex()) == stored_trackIds.end()) {
      emprimarytracks(collision.posZ(), collision.numContrib(),
                      track.pt(), track.eta(), track.phi(), track.tgl(), track.signed1Pt(),
                      track.dcaXY(), track.dcaZ(), track.cYY(), track.cZZ(), track.cZY(),
                      track.tpcNClsFindable(), track.tpcNClsFound(), track.tpcNClsCrossedRows(),
                      track.tpcChi2NCl(), track.tpcInnerParam(),
                      track.tpcSignal(), track.tpcNSigmaEl(), track.tpcNSigmaMu(), track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr(),
                      track.beta(), track.tofNSigmaEl(), track.tofNSigmaMu(), track.tofNSigmaPi(), track.tofNSigmaKa(), track.tofNSigmaPr(),
                      track.itsClusterSizes(), track.itsChi2NCl(), track.detectorMap(), pidlabel);
      stored_trackIds.emplace_back(track.globalIndex());
    }
  }

  Filter collisionFilter_common = nabs(o2::aod::collision::posZ) < 10.f && o2::aod::collision::numContrib > (uint16_t)0 && o2::aod::evsel::sel8 == true;
  using filteredMyCollisions = soa::Filtered<MyCollisions>;

  //! type of V0. 0: built solely for cascades (does not pass standard V0 cuts), 1: standard 2, 3: photon-like with TPC-only use. Regular analysis should always use type 1 or 3.
  Filter v0Filter = o2::aod::v0::v0Type == (uint8_t)1;
  using filteredMyV0s = soa::Filtered<aod::V0s>;
  Preslice<filteredMyV0s> perCollision_v0 = o2::aod::v0::collisionId;

  // Filter trackFilter = o2::aod::track::pt > minpt && nabs(o2::aod::track::eta) < maxeta && o2::aod::track::tpcChi2NCl < maxchi2tpc && o2::aod::track::itsChi2NCl < maxchi2its && nabs(o2::aod::track::dcaXY) < maxdcaXY && nabs(o2::aod::track::dcaZ) < maxdcaZ;
  // Filter trackFilter2 = ncheckbit(aod::track::detectorMap, (uint8_t)o2::aod::track::ITS) == true && ncheckbit(aod::track::detectorMap, (uint8_t)o2::aod::track::TPC) == true;
  // Filter eidFilter = minTPCNsigmaEl < o2::aod::pidtpc::tpcNSigmaEl && o2::aod::pidtpc::tpcNSigmaEl < maxTPCNsigmaEl && ((0.95f < o2::aod::pidtofbeta::beta && o2::aod::pidtofbeta::beta < 1.05f) || o2::aod::pidtofbeta::beta < 0.f);
  // using filteredMyTracks = soa::Filtered<MyTracks>;

  // Don't apply filter to tracks, because posTrack_as<>, negTrack_as<> is used.
  Preslice<MyTracks> perCollision_track = o2::aod::track::collisionId;
  Partition<MyTracks> posTracks = o2::aod::track::signed1Pt > 0.f && o2::aod::track::pt > minpt&& nabs(o2::aod::track::eta) < maxeta&& o2::aod::track::tpcChi2NCl < maxchi2tpc&& o2::aod::track::itsChi2NCl < maxchi2its&& nabs(o2::aod::track::dcaXY) < maxdcaXY&& nabs(o2::aod::track::dcaZ) < maxdcaZ&& ncheckbit(aod::track::detectorMap, (uint8_t)o2::aod::track::ITS) == true && ncheckbit(aod::track::detectorMap, (uint8_t)o2::aod::track::TPC) == true && minTPCNsigmaEl < o2::aod::pidtpc::tpcNSigmaEl&& o2::aod::pidtpc::tpcNSigmaEl < maxTPCNsigmaEl && ((0.95f < o2::aod::pidtofbeta::beta && o2::aod::pidtofbeta::beta < 1.05f) || o2::aod::pidtofbeta::beta < 0.f);
  Partition<MyTracks> negTracks = o2::aod::track::signed1Pt < 0.f && o2::aod::track::pt > minpt&& nabs(o2::aod::track::eta) < maxeta&& o2::aod::track::tpcChi2NCl < maxchi2tpc&& o2::aod::track::itsChi2NCl < maxchi2its&& nabs(o2::aod::track::dcaXY) < maxdcaXY&& nabs(o2::aod::track::dcaZ) < maxdcaZ&& ncheckbit(aod::track::detectorMap, (uint8_t)o2::aod::track::ITS) == true && ncheckbit(aod::track::detectorMap, (uint8_t)o2::aod::track::TPC) == true && minTPCNsigmaEl < o2::aod::pidtpc::tpcNSigmaEl&& o2::aod::pidtpc::tpcNSigmaEl < maxTPCNsigmaEl && ((0.95f < o2::aod::pidtofbeta::beta && o2::aod::pidtofbeta::beta < 1.05f) || o2::aod::pidtofbeta::beta < 0.f);
  std::vector<uint64_t> stored_trackIds;

  void processPID(filteredMyCollisions const& collisions, aod::BCsWithTimestamps const&, filteredMyV0s const& v0s, MyTracks const& tracks)
  {
    stored_trackIds.reserve(tracks.size());
    for (auto& collision : collisions) {
      registry.fill(HIST("hEventCounter"), 1.0); // all

      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      std::array<float, 3> pVtx = {collision.posX(), collision.posY(), collision.posZ()};
      auto v0s_coll = v0s.sliceBy(perCollision_v0, collision.globalIndex());
      for (auto& v0 : v0s_coll) {
        auto pos = v0.template posTrack_as<MyTracks>();
        auto neg = v0.template negTrack_as<MyTracks>();
        if (!IsSelected(pos) || !IsSelected(neg)) {
          continue;
        }

        // reconstruct V0s
        auto pTrack = getTrackParCov(pos); // positive
        auto nTrack = getTrackParCov(neg); // negative
        std::array<float, 3> svpos = {0.}; // secondary vertex position
        std::array<float, 3> pvec0 = {0.};
        std::array<float, 3> pvec1 = {0.};

        int nCand = fitter.process(pTrack, nTrack);
        if (nCand != 0) {
          fitter.propagateTracksToVertex();
          const auto& vtx = fitter.getPCACandidate();
          for (int i = 0; i < 3; i++) {
            svpos[i] = vtx[i];
          }
          fitter.getTrack(0).getPxPyPzGlo(pvec0); // positive
          fitter.getTrack(1).getPxPyPzGlo(pvec1); // negative
        } else {
          continue;
        }
        std::array<float, 3> pvxyz{pvec0[0] + pvec1[0], pvec0[1] + pvec1[1], pvec0[2] + pvec1[2]};

        float v0dca = std::sqrt(fitter.getChi2AtPCACandidate()); // distance between 2 legs.
        float v0CosinePA = RecoDecay::cpa(pVtx, svpos, pvxyz);
        registry.fill(HIST("hPCA"), v0dca);
        registry.fill(HIST("hCosPA"), v0CosinePA);

        if (v0dca > maxdcav0dau) {
          continue;
        }
        if (v0CosinePA < minv0cospa) {
          continue;
        }

        float alpha = v0_alpha(pvec0[0], pvec0[1], pvec0[2], pvec1[0], pvec1[1], pvec1[2]);
        float qt = v0_qt(pvec0[0], pvec0[1], pvec0[2], pvec1[0], pvec1[1], pvec1[2]);
        registry.fill(HIST("hAP"), alpha, qt);

        float mGamma = RecoDecay::m(std::array{std::array{pvec0[0], pvec0[1], pvec0[2]}, std::array{pvec1[0], pvec1[1], pvec1[2]}}, std::array{o2::constants::physics::MassElectron, o2::constants::physics::MassElectron});
        float mK0Short = RecoDecay::m(std::array{std::array{pvec0[0], pvec0[1], pvec0[2]}, std::array{pvec1[0], pvec1[1], pvec1[2]}}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassPionCharged});
        float mLambda = RecoDecay::m(std::array{std::array{pvec0[0], pvec0[1], pvec0[2]}, std::array{pvec1[0], pvec1[1], pvec1[2]}}, std::array{o2::constants::physics::MassProton, o2::constants::physics::MassPionCharged});
        float mAntiLambda = RecoDecay::m(std::array{std::array{pvec0[0], pvec0[1], pvec0[2]}, std::array{pvec1[0], pvec1[1], pvec1[2]}}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassProton});

        EM_V0_Label v0id = checkV0(alpha, qt);
        if (v0id == EM_V0_Label::kGamma && (abs(pos.tpcNSigmaEl()) < 5.f && abs(neg.tpcNSigmaEl()) < 5.f)) {
          registry.fill(HIST("hMassGamma"), mGamma);
          if (0.f < mGamma && mGamma < 0.02) {
            fillTrackTable(collision, pos, static_cast<int>(o2::aod::pwgem::dilepton::PID_Label::kElectron));
            fillTrackTable(collision, neg, static_cast<int>(o2::aod::pwgem::dilepton::PID_Label::kElectron));
          }
        } else if (v0id == EM_V0_Label::kK0S && (abs(pos.tpcNSigmaPi()) < 5.f && abs(neg.tpcNSigmaPi()) < 5.f)) {
          registry.fill(HIST("hMassK0Short"), mK0Short);
          if (abs(mK0Short - 0.497) < 0.01) {
            if (dist01(engine) < downscaling_pion) {
              fillTrackTable(collision, pos, static_cast<int>(o2::aod::pwgem::dilepton::PID_Label::kPion));
              fillTrackTable(collision, neg, static_cast<int>(o2::aod::pwgem::dilepton::PID_Label::kPion));
            }
          }
        } else if (v0id == EM_V0_Label::kLambda && (abs(pos.tpcNSigmaPr()) < 5.f && abs(neg.tpcNSigmaPi()) < 5.f)) {
          registry.fill(HIST("hMassLambda"), mLambda);
          if (abs(mLambda - 1.115) < 0.005) {
            fillTrackTable(collision, pos, static_cast<int>(o2::aod::pwgem::dilepton::PID_Label::kProton));
            if (dist01(engine) < downscaling_pion) {
              fillTrackTable(collision, neg, static_cast<int>(o2::aod::pwgem::dilepton::PID_Label::kPion));
            }
          }
        } else if (v0id == EM_V0_Label::kAntiLambda && (abs(pos.tpcNSigmaPi()) < 5.f && abs(neg.tpcNSigmaPr()) < 5.f)) {
          registry.fill(HIST("hMassAntiLambda"), mAntiLambda);
          if (abs(mAntiLambda - 1.115) < 0.005) {
            if (dist01(engine) < downscaling_pion) {
              fillTrackTable(collision, pos, static_cast<int>(o2::aod::pwgem::dilepton::PID_Label::kPion));
            }
            fillTrackTable(collision, neg, static_cast<int>(o2::aod::pwgem::dilepton::PID_Label::kProton));
          }
        }
      } // end of v0 loop

      // extra statistics for electrons tagged by photon conversions
      auto posTracks_per_coll = posTracks->sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), cache);
      auto negTracks_per_coll = negTracks->sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), cache);
      for (auto& [pos, ele] : combinations(CombinationsFullIndexPolicy(posTracks_per_coll, negTracks_per_coll))) {
        if (!IsSelected(pos) || !IsSelected(ele)) {
          continue;
        }

        ROOT::Math::PtEtaPhiMVector v1(ele.pt(), ele.eta(), ele.phi(), o2::constants::physics::MassElectron);
        ROOT::Math::PtEtaPhiMVector v2(pos.pt(), pos.eta(), pos.phi(), o2::constants::physics::MassElectron);
        ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
        float phiv = getPhivPair(pos.px(), pos.py(), pos.pz(), ele.px(), ele.py(), ele.pz(), pos.sign(), ele.sign(), d_bz);
        registry.fill(HIST("hMvsPhiV"), phiv, v12.M());

        if (v12.M() < slope * phiv + intercept) { // photon conversion is found.
          fillTrackTable(collision, ele, static_cast<int>(o2::aod::pwgem::dilepton::PID_Label::kElectron));
          fillTrackTable(collision, pos, static_cast<int>(o2::aod::pwgem::dilepton::PID_Label::kElectron));
        }

      } // end of ULS pair loop

    } // end of collision loop
    stored_trackIds.clear();
    stored_trackIds.shrink_to_fit();
  } // end of process
  PROCESS_SWITCH(TreeCreatorElectronMLDDA, processPID, "produce ML input for single track level", false);

  void processDummy(filteredMyCollisions const& collisions) {}
  PROCESS_SWITCH(TreeCreatorElectronMLDDA, processDummy, "process dummy", true);
};

struct MLTrackQC {
  // Basic checks
  HistogramRegistry registry{
    "registry",
    {
      {"hTPCdEdx_Pin_All", "TPC dE/dx vs. pin;p_{in} (GeV/c);TPC dE/dx", {HistType::kTH2F, {{1000, 0.f, 10.f}, {200, 0, 200}}}},
      {"hTOFbeta_Pin_All", "TOF beta vs. pin;p_{in} (GeV/c);TOF #beta", {HistType::kTH2F, {{1000, 0.f, 10.f}, {220, 0.0, 1.1}}}},
      {"hTPCdEdx_Pin_Electron", "TPC dE/dx vs. pin;p_{in} (GeV/c);TPC dE/dx", {HistType::kTH2F, {{1000, 0.f, 10.f}, {200, 0, 200}}}},
      {"hTOFbeta_Pin_Electron", "TOF beta vs. pin;p_{in} (GeV/c);TOF #beta", {HistType::kTH2F, {{1000, 0.f, 10.f}, {220, 0.0, 1.1}}}},
      {"hTPCdEdx_Pin_Pion", "TPC dE/dx vs. pin;p_{in} (GeV/c);TPC dE/dx", {HistType::kTH2F, {{1000, 0.f, 10.f}, {200, 0, 200}}}},
      {"hTOFbeta_Pin_Pion", "TOF beta vs. pin;p_{in} (GeV/c);TOF #beta", {HistType::kTH2F, {{1000, 0.f, 10.f}, {220, 0.0, 1.1}}}},
      {"hTPCdEdx_Pin_Proton", "TPC dE/dx vs. pin;p_{in} (GeV/c);TPC dE/dx", {HistType::kTH2F, {{1000, 0.f, 10.f}, {200, 0, 200}}}},
      {"hTOFbeta_Pin_Proton", "TOF beta vs. pin;p_{in} (GeV/c);TOF #beta", {HistType::kTH2F, {{1000, 0.f, 10.f}, {220, 0.0, 1.1}}}},
    },
  };

  void process(aod::EMPrimaryTracks const& tracks)
  {
    for (auto& track : tracks) {
      registry.fill(HIST("hTPCdEdx_Pin_All"), track.tpcInnerParam(), track.tpcSignal());
      registry.fill(HIST("hTOFbeta_Pin_All"), track.tpcInnerParam(), track.beta());
      if (track.pidlabel() == static_cast<int>(o2::aod::pwgem::dilepton::PID_Label::kElectron)) {
        registry.fill(HIST("hTPCdEdx_Pin_Electron"), track.tpcInnerParam(), track.tpcSignal());
        registry.fill(HIST("hTOFbeta_Pin_Electron"), track.tpcInnerParam(), track.beta());
      } else if (track.pidlabel() == static_cast<int>(o2::aod::pwgem::dilepton::PID_Label::kPion)) {
        registry.fill(HIST("hTPCdEdx_Pin_Pion"), track.tpcInnerParam(), track.tpcSignal());
        registry.fill(HIST("hTOFbeta_Pin_Pion"), track.tpcInnerParam(), track.beta());
      } else if (track.pidlabel() == static_cast<int>(o2::aod::pwgem::dilepton::PID_Label::kProton)) {
        registry.fill(HIST("hTPCdEdx_Pin_Proton"), track.tpcInnerParam(), track.tpcSignal());
        registry.fill(HIST("hTOFbeta_Pin_Proton"), track.tpcInnerParam(), track.beta());
      }
    } // end of track loop
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<TreeCreatorElectronMLDDA>(cfgc, TaskName{"tree-creator-ele-ml-dda"}),
    adaptAnalysisTask<MLTrackQC>(cfgc, TaskName{"ml-track-qc"})};
}
