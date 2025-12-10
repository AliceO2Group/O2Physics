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
/// \brief The task for muon QA
/// \author Andrea Ferrero <andrea.ferrero@cern.ch>
/// \author Paul Veen <paul.veen@cern.ch>
/// \author Chi Zhang <chi.zhang@cern.ch>

#include "Common/DataModel/EventSelection.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DataFormatsMCH/Cluster.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DetectorsBase/GRPGeomHelper.h>
#include <DetectorsBase/GeometryManager.h>
#include <DetectorsBase/Propagator.h>
#include <Field/MagneticField.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/DataTypes.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
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
#include <ReconstructionDataFormats/GlobalFwdTrack.h>
#include <ReconstructionDataFormats/TrackFwd.h>

#include <Math/MatrixRepresentationsStatic.h>
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
#include <exception>
#include <format>
#include <iterator>
#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <math.h> // FIXME: Replace M_PI

using namespace o2;
using namespace o2::aod;
using namespace o2::mch;
using namespace o2::framework;
using namespace o2::framework::expressions;

using MyEvents = soa::Join<aod::Collisions, aod::EvSels>;
using MyMuonsWithCov = soa::Join<aod::FwdTracks, aod::FwdTracksCov>;
using MyMFTs = aod::MFTTracks;

using SMatrix55 = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
using SMatrix5 = ROOT::Math::SVector<Double_t, 5>;

using MuonPair = std::pair<std::pair<uint64_t, uint64_t>, std::pair<uint64_t, uint64_t>>;
using GlobalMuonPair = std::pair<std::pair<uint64_t, std::vector<uint64_t>>, std::pair<uint64_t, std::vector<uint64_t>>>;

const int fgNCh = 10;
const int fgNDetElemCh[fgNCh] = {4, 4, 4, 4, 18, 18, 26, 26, 26, 26};
const int fgSNDetElemCh[fgNCh + 1] = {0, 4, 8, 12, 16, 34, 52, 78, 104, 130, 156};
const float zAtAbsEnd = -505.;

constexpr double firstMFTPlaneZ = o2::mft::constants::mft::LayerZCoordinate()[0];
constexpr double lastMFTPlaneZ = o2::mft::constants::mft::LayerZCoordinate()[9];

std::array<double, 5> zRefPlane{
  firstMFTPlaneZ,
  lastMFTPlaneZ,
  -90.0,
  -300.0,
  //-505.0,
  -520.0};

std::vector<std::pair<std::string, double>> referencePlanes{
  {"MFT-begin", 10.0},
  {"MFT-end", 15.0},
  {"Absorber-begin", 20.0},
  {"Absorber-mid", 75.0},
  //{"Absorber-end", 100.0},
  {"MCH-begin", 100.0}};

enum MuonExtrapolation {
  // Index used to set different options for muon propagation
  kToVtx = 0, // propagtion to vertex by default
  kToDCA,
  kToAbsEnd,
  kToZ
};

struct VarColl {
  int64_t globalIndex = 0;
  float x = 0.f;
  float y = 0.f;
  float z = 0.f;
  float covXX = 0.f;
  float covYY = 0.f;
  int64_t bc = 0;
  int multMFT = 0;
};

struct VarTrack {
  int64_t collisionId = -1;
  int64_t globalIndex = 0;
  int nClusters = 0; // Only MCH
  int sign = 0;
  int64_t bc = 0;
  int trackType = 0;
  float trackTime = 0.f;

  // Basic kinematics
  float x = 0.f;
  float y = 0.f;
  float z = 0.f;
  float eta = 0.f;
  float phi = 0.f;
  float tgl = 0.f;

  float px = 0.f;
  float py = 0.f;
  float pz = 0.f;
  float pT = 0.f;
  float p = 0.f;

  // Propagation related infos
  float dcaX = 0.f;
  float dcaY = 0.f;
  float dcaXY = 0.f;
  float pDca = 0.f;
  float rabs = 0.f;
  float chi2 = 0.f;
  float chi2matching = 0.f;
};

struct VarClusters {
  std::vector<std::vector<float>> posClusters;   // (x,y,z)
  std::vector<std::vector<float>> errorClusters; // (ex,ey)
  std::vector<int> DEIDs;
};

struct muonQa {
  ////   Variables for enabling QA options
  struct : ConfigurableGroup {
    Configurable<bool> fEnableQAMatching{"cfgEnableQAMatching", false, "Enable MCH-MFT matching QA checks"};
    Configurable<bool> fEnableQAResidual{"cfgEnableQAResidual", false, "Enable residual QA checks"};
    Configurable<bool> fEnableQADCA{"cfgEnableQADCA", false, "Enable DCA QA checks"};
    Configurable<bool> fEnableQADimuon{"cfgEnableQADimuon", false, "Enable dimuon QA checks"};
    Configurable<bool> fEnableQADimuonSameSignDCA{"cfgEnableQADimuonSameSignDCA", false, "Enable same sign dimuon DCA QA checks"};
    Configurable<bool> fEnableSingleMuonDiMuonCorrelations{"cfgEnableMuonDiMuonCorrelations", false, "Enable muon-dimuon QA checks"};
  } configQAs;

  ////   Variables for selecting muon tracks
  struct : ConfigurableGroup {
    Configurable<float> fPMchLow{"cfgPMchLow", 0.0f, ""};
    Configurable<float> fPtMchLow{"cfgPtMchLow", 0.7f, ""};
    Configurable<float> fEtaMchLow{"cfgEtaMchLow", -4.0f, ""};
    Configurable<float> fEtaMchUp{"cfgEtaMchUp", -2.5f, ""};
    Configurable<float> fRabsLow{"cfgRabsLow", 17.6f, ""};
    Configurable<float> fRabsUp{"cfgRabsUp", 89.5f, ""};
    Configurable<float> fSigmaPdcaUp{"cfgPdcaUp", 6.f, ""};
    Configurable<float> fTrackChi2MchUp{"cfgTrackChi2MchUp", 5.f, ""};
    Configurable<float> fMatchingChi2MchMidUp{"cfgMatchingChi2MchMidUp", 999.f, ""};
  } configMuons;

  ////   Variables for selecting mft tracks
  struct : ConfigurableGroup {
    Configurable<float> fEtaMftLow{"cfgEtaMftlow", -3.6f, ""};
    Configurable<float> fEtaMftUp{"cfgEtaMftup", -2.5f, ""};
    Configurable<int> fTrackNClustMftLow{"cfgTrackNClustMftLow", 7, ""};
    Configurable<float> fTrackChi2MftUp{"cfgTrackChi2MftUp", 999.f, ""};
  } configMFTs;

  ////   Variables for selecting global tracks
  Configurable<float> fMatchingChi2MftMchUp{"cfgMatchingChi2MftMchUp", 50.f, ""};
  Configurable<float> fMchPUpForGlobalDCA{"cfgMchPUpForGlobalDCA", 20.f, ""};

  ////   Variables for selecting dimuon DCA candidates
  Configurable<float> fDimuonDCAMassLow{"cfgDimuonDCAMassLow", 2.8f, ""};
  Configurable<float> fDimuonDCAMassHigh{"cfgDimuonDCAMassHigh", 3.4f, ""};

  ////   Variables for alignment corrections
  Configurable<bool> fEnableMFTAlignmentCorrections{"cfgEnableMFTAlignmentCorrections", false, ""};

  ////   Variables for re-alignment setup
  struct : ConfigurableGroup {
    Configurable<bool> fDoRealign{"cfgDoRealign", false, "Switch to apply re-alignment"};
    Configurable<double> fChamberResolutionX{"cfgChamberResolutionX", 0.4, "Chamber resolution along X configuration for refit"}; // 0.4cm pp, 0.2cm PbPb
    Configurable<double> fChamberResolutionY{"cfgChamberResolutionY", 0.4, "Chamber resolution along Y configuration for refit"}; // 0.4cm pp, 0.2cm PbPb
    Configurable<double> fSigmaCutImprove{"cfgSigmaCutImprove", 6., "Sigma cut for track improvement"};
  } configRealign;

  ///    Variables to event mixing criteria
  struct : ConfigurableGroup {
    Configurable<int> fEventMaxDeltaNMFT{"cfgEventMaxDeltaNMFT", 1, ""};
    Configurable<float> fEventMaxDeltaVtxZ{"cfgEventMaxDeltaVtxZ", 1.f, ""};
    Configurable<uint64_t> fEventMinDeltaBc{"cfgEventMinDeltaBc", 500, ""};
  } configMixing;

  ////   Variables for ccdb
  struct : ConfigurableGroup {
    Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
    Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
    Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
    Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
    Configurable<std::string> geoPathRealign{"geoPathRealign", "GLO/Config/GeometryAligned", "Path of the geometry file"};
    Configurable<int64_t> nolaterthan{"ccdb-no-later-than-ref", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object of reference basis"};
    Configurable<int64_t> nolaterthanRealign{"ccdb-no-later-than-new", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object of new basis"};
  } configCCDB;

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

  ////    Variables for histograms configuration
  Configurable<int> fNCandidatesMax{"nCandidatesMax", 5, ""};

  parameters::GRPMagField* grpmag = nullptr;
  TrackFitter trackFitter; // Track fitter from MCH tracking library

  globaltracking::MatchGlobalFwd mMatching;
  int fCurrentRun;        // needed to detect if the run changed and trigger update of calibrations etc.
  double mImproveCutChi2; // Chi2 cut for track improvement.
  Service<ccdb::BasicCCDBManager> ccdb;
  o2::field::MagneticField* fieldB = nullptr;
  double Bz; // Bz for MFT

  geo::TransformationCreator transformation;
  std::map<int, math_utils::Transform3D> transformRef; // reference geometry w.r.t track data
  std::map<int, math_utils::Transform3D> transformNew; // new geometry
  TGeoManager* geoNew = nullptr;
  TGeoManager* geoRef = nullptr;

  Preslice<aod::FwdTrkCl> perMuon = aod::fwdtrkcl::fwdtrackId;
  Preslice<MyMuonsWithCov> fwdtracksPerCollision = aod::fwdtrack::collisionId;
  Preslice<MyMFTs> mftPerCollision = aod::fwdtrack::collisionId;

  HistogramRegistry registry{"registry", {}};
  HistogramRegistry registryDCA{"registryDCA", {}};
  HistogramRegistry registryDCAdiMuons{"registryDCAdiMuons", {}};
  HistogramRegistry registryResiduals{"registryResiduals", {}};
  HistogramRegistry registryResidualsMFT{"registryResidualsMFT", {}};
  HistogramRegistry registryResidualsMCH{"registryResidualsMCH", {}};
  HistogramRegistry registryDimuon{"registryDimuon", {}};

  std::array<std::string, 4> quadrants = {"Q0", "Q1", "Q2", "Q3"};

  std::array<std::array<std::array<std::unordered_map<std::string, o2::framework::HistPtr>, 3>, 4>, 2> dcaHistos;
  std::array<std::array<std::array<std::unordered_map<std::string, o2::framework::HistPtr>, 3>, 4>, 2> dcaHistosMixedEvents;
  std::array<std::array<std::array<std::unordered_map<std::string, o2::framework::HistPtr>, 3>, 4>, 2> dcaHistosGlobal;
  std::array<std::array<std::unordered_map<std::string, o2::framework::HistPtr>, 3>, 4> dcaHistosGlobalSubtracted;
  std::array<std::array<std::unordered_map<std::string, o2::framework::HistPtr>, 3>, 4> dcaHistosGlobalSubtractedMCHpCut;

  std::array<std::array<std::unordered_map<std::string, o2::framework::HistPtr>, 4>, 6> trackResidualsHistos;
  std::array<std::array<std::unordered_map<std::string, o2::framework::HistPtr>, 4>, 6> trackResidualsHistosMixedEvents;

  std::array<std::array<std::unordered_map<std::string, o2::framework::HistPtr>, 10>, 4> residualsHistos;
  std::array<std::array<std::unordered_map<std::string, o2::framework::HistPtr>, 10>, 4> residualsHistosMixedEvents;

  std::array<std::array<std::array<std::unordered_map<std::string, o2::framework::HistPtr>, 10>, 2>, 2> residualsHistosPerDE;
  std::array<std::array<std::array<std::unordered_map<std::string, o2::framework::HistPtr>, 10>, 2>, 2> residualsHistosPerDEMixedEvents;

  std::array<std::array<std::array<std::unordered_map<std::string, o2::framework::HistPtr>, 10>, 2>, 2> mchResidualsHistosPerDE;
  std::array<std::array<std::array<std::unordered_map<std::string, o2::framework::HistPtr>, 10>, 2>, 2> mchResidualsHistosPerDEMixedEvents;

  VarTrack fgValuesMCH;
  VarTrack fgValuesMCHpv;
  VarTrack fgValuesMFT;
  VarTrack fgValuesGlobal;
  std::vector<VarTrack> fgValuesCandidates;

  void CreateBasicHistograms()
  {
    // ======================
    // Muons plots
    // ======================

    AxisSpec chi2Axis = {1000, 0, 1000, "chi2"};
    AxisSpec momentumAxis = {1000, 0, 1000, "p (GeV/c)"};
    AxisSpec transverseMomentumAxis = {1000, 0, 100, "p_{T} (GeV/c)"};
    AxisSpec etaAxis = {80, -5, -1, "#eta"};
    AxisSpec rAbsAxis = {100, 0., 100.0, "R_{abs} (cm)"};
    AxisSpec dcaAxis = {400, 0.0, 20.0, "DCA"};
    AxisSpec pdcaAxis = {5000, 0.0, 5000.0, "p #times DCA"};
    AxisSpec phiAxis = {360, -180.0, 180.0, "#phi (degrees)"};

    registry.add("muons/TrackChi2", "MCH track #chi^{2}", {HistType::kTH1F, {chi2Axis}});
    registry.add("muons/TrackP", "MCH track momentum", {HistType::kTH1F, {momentumAxis}});
    registry.add("muons/TrackPt", "MCH track transverse momentum", {HistType::kTH1F, {transverseMomentumAxis}});
    registry.add("muons/TrackEta", "MCH track #eta", {HistType::kTH1F, {etaAxis}});
    registry.add("muons/TrackRabs", "MCH track R_{abs}", {HistType::kTH1F, {rAbsAxis}});
    registry.add("muons/TrackDCA", "MCH track DCA", {HistType::kTH1F, {dcaAxis}});
    registry.add("muons/TrackPDCA", "MCH track p #times DCA", {HistType::kTH1F, {pdcaAxis}});
    registry.add("muons/TrackPhi", "MCH track #phi", {HistType::kTH1F, {phiAxis}});

    // muon origin from MCH quadrants
    registry.add("muons/TrackEtaPos", "MCH #mu^{+} track #eta", {HistType::kTH1F, {etaAxis}});
    registry.add("muons/TrackEtaNeg", "MCH #mu^{-} track #eta", {HistType::kTH1F, {etaAxis}});
    // -- pT and eta
    registry.add("muons/TrackPt_TrackEtaPos", "track pT and #eta", {HistType::kTH2F, {transverseMomentumAxis, etaAxis}});
    registry.add("muons/TrackPt_TrackEtaNeg", "track pT and #eta", {HistType::kTH2F, {transverseMomentumAxis, etaAxis}});
    // top-bottom
    registry.add("muons/TrackEtaPos_T", "MCH #mu^{+} track #eta, top MCH CH1", {HistType::kTH1F, {etaAxis}});
    registry.add("muons/TrackEtaPos_B", "MCH #mu^{+} track #eta, bottom MCH CH1", {HistType::kTH1F, {etaAxis}});
    registry.add("muons/TrackEtaNeg_T", "MCH #mu^{-} track #eta, top MCH CH1", {HistType::kTH1F, {etaAxis}});
    registry.add("muons/TrackEtaNeg_B", "MCH #mu^{-} track #eta, bottom MCH CH1", {HistType::kTH1F, {etaAxis}});
    // -- pT and eta
    registry.add("muons/TrackPt_TrackEtaPos_T", "track p_{T} and #eta, top MCH CH1", {HistType::kTH2F, {transverseMomentumAxis, etaAxis}});
    registry.add("muons/TrackPt_TrackEtaPos_B", "track p_{T} and #eta, bottom MCH CH1", {HistType::kTH2F, {transverseMomentumAxis, etaAxis}});
    registry.add("muons/TrackPt_TrackEtaNeg_T", "track p_{T} and #eta, top MCH CH1", {HistType::kTH2F, {transverseMomentumAxis, etaAxis}});
    registry.add("muons/TrackPt_TrackEtaNeg_B", "track p_{T} and #eta, bottom MCH CH1", {HistType::kTH2F, {transverseMomentumAxis, etaAxis}});
    // left-right
    registry.add("muons/TrackEtaPos_L", "MCH #mu^{+} track #eta, left MCH CH1", {HistType::kTH1F, {etaAxis}});
    registry.add("muons/TrackEtaPos_R", "MCH #mu^{+} track #eta, right MCH CH1", {HistType::kTH1F, {etaAxis}});
    registry.add("muons/TrackEtaNeg_L", "MCH #mu^{-} track #eta, left MCH CH1", {HistType::kTH1F, {etaAxis}});
    registry.add("muons/TrackEtaNeg_R", "MCH #mu^{-} track #eta, right MCH CH1", {HistType::kTH1F, {etaAxis}});
    // -- pT and eta
    registry.add("muons/TrackPt_TrackEtaPos_L", "track p_{T} and #eta, top MCH CH1", {HistType::kTH2F, {transverseMomentumAxis, etaAxis}});
    registry.add("muons/TrackPt_TrackEtaPos_R", "track p_{T} and #eta, bottom MCH CH1", {HistType::kTH2F, {transverseMomentumAxis, etaAxis}});
    registry.add("muons/TrackPt_TrackEtaNeg_L", "track p_{T} and #eta, top MCH CH1", {HistType::kTH2F, {transverseMomentumAxis, etaAxis}});
    registry.add("muons/TrackPt_TrackEtaNeg_R", "track p_{T} and #eta, bottom MCH CH1", {HistType::kTH2F, {transverseMomentumAxis, etaAxis}});

    // ======================
    // Global muons plots
    // ======================
    int nTrackTypes = static_cast<int>(o2::aod::fwdtrack::ForwardTrackTypeEnum::MCHStandaloneTrack) + 1;
    AxisSpec trackTypeAxis = {static_cast<int>(nTrackTypes), 0.0, static_cast<double>(nTrackTypes), "track type"};
    registry.add("global-muons/nTracksPerType", "Number of tracks per type", {HistType::kTH1F, {trackTypeAxis}});

    AxisSpec nCandidatesAxis = {static_cast<int>(fNCandidatesMax), 0.0, static_cast<double>(fNCandidatesMax), "match candidate rank"};
    registry.add("global-muons/NCandidates", "Number of MFT-MCH match candidates", {HistType::kTH1F, {nCandidatesAxis}});
    registry.add("global-muons/MatchChi2", "MFT-MCH match chi2", {HistType::kTH2F, {chi2Axis, nCandidatesAxis}});

    registry.add("global-muons/TrackChi2", "Muon track #chi^{2}", {HistType::kTH1F, {chi2Axis}});
    registry.add("global-muons/TrackP", "Muon track momentum", {HistType::kTH1F, {momentumAxis}});
    registry.add("global-muons/TrackPt", "Muon track transverse momentum", {HistType::kTH1F, {transverseMomentumAxis}});
    registry.add("global-muons/TrackEta", "Muon track #eta", {HistType::kTH1F, {etaAxis}});
    registry.add("global-muons/TrackRabs", "Muon track R_{abs}", {HistType::kTH1F, {rAbsAxis}});
    registry.add("global-muons/TrackDCA", "Muon track DCA", {HistType::kTH1F, {dcaAxis}});
    registry.add("global-muons/TrackPDCA", "Muon track p #times DCA", {HistType::kTH1F, {pdcaAxis}});
    registry.add("global-muons/TrackPhi", "Muon track #phi", {HistType::kTH1F, {phiAxis}});

    // ======================
    // Global muon plots with matching cuts
    // ======================

    if (configQAs.fEnableQAMatching) {
      AxisSpec dbcAxis = {1000, -500, 500, "#Delta_{BC}"};
      registry.add("global-matches/BCdifference", "MCH-MFT BC difference", {HistType::kTH1F, {dbcAxis}});

      AxisSpec nClustersAxis = {20, 0, 20, "# of MFT clusters per track"};

      registry.add("global-matches/MatchChi2", "MFT-MCH match chi2", {HistType::kTH1F, {chi2Axis}});

      registry.add("global-matches/TrackChi2_MFT", "MFT track #chi^{2}", {HistType::kTH1F, {chi2Axis}});
      registry.add("global-matches/TrackNclusters_MFT", "MFT track Nclusters", {HistType::kTH1F, {nClustersAxis}});

      registry.add("global-matches/TrackChi2", "Muon track #chi^{2}", {HistType::kTH1F, {chi2Axis}});
      registry.add("global-matches/TrackP", "Muon track momentum", {HistType::kTH1F, {momentumAxis}});
      registry.add("global-matches/TrackPt", "Muon track transverse momentum", {HistType::kTH1F, {transverseMomentumAxis}});
      registry.add("global-matches/TrackEta", "Muon track #eta", {HistType::kTH1F, {etaAxis}});
      registry.add("global-matches/TrackRabs", "Muon track R_{abs}", {HistType::kTH1F, {rAbsAxis}});
      registry.add("global-matches/TrackDCA", "Muon track DCA", {HistType::kTH1F, {dcaAxis}});
      registry.add("global-matches/TrackPDCA", "Muon track p #times DCA", {HistType::kTH1F, {pdcaAxis}});
      registry.add("global-matches/TrackPhi", "Muon track #phi", {HistType::kTH1F, {phiAxis}});

      registry.add("global-matches/TrackP_glo", "Global muon track momentum", {HistType::kTH1F, {momentumAxis}});
      registry.add("global-matches/TrackPt_glo", "Global muon track transverse momentum", {HistType::kTH1F, {transverseMomentumAxis}});
      registry.add("global-matches/TrackEta_glo", "Global muon track #eta", {HistType::kTH1F, {etaAxis}});
      registry.add("global-matches/TrackDCA_glo", "Global muon track DCA", {HistType::kTH1F, {dcaAxis}});
      registry.add("global-matches/TrackPhi_glo", "Global muon track #phi", {HistType::kTH1F, {phiAxis}});
    }

    AxisSpec momentumCorrelationAxis = {100, 0, 100, "momentum (GeV/c)"};
    AxisSpec momentumDeltaAxis = {100, -1, 1, "#DeltaP (GeV/c)"};
    // Momentum correlations
    registry.add("global-muons/MomentumCorrelation_Global_vs_Muon",
                 "P_{global} vs. P_{MCH}",
                 {HistType::kTH2F, {momentumCorrelationAxis, momentumCorrelationAxis}});
    registry.add("global-muons/MomentumDifference_Global_vs_Muon",
                 "(P_{global} - P_{MCH}) / P_{MCH} vs. P_{MCH}",
                 {HistType::kTH2F, {momentumCorrelationAxis, momentumDeltaAxis}});
    registry.add("global-muons/MomentumCorrelation_subleading_vs_leading",
                 "P_{subleading_match} vs. P_{leading_match}",
                 {HistType::kTH2F, {momentumCorrelationAxis, momentumCorrelationAxis}});
    registry.add("global-muons/MomentumDifference_subleading_vs_leading",
                 "(P_{subleading_match} - P_{leading_match}) / P_{leading_match} vs. P_{leading_match}",
                 {HistType::kTH2F, {momentumCorrelationAxis, momentumDeltaAxis}});

    // AxisSpec etaAxis = {100, -5.0, -2.0, "#eta"};
    AxisSpec etaCorrelationAxis = {80, -5.0, -1.0, "#eta"};
    AxisSpec etaDeltaAxis = {100, -0.2, 0.2, "#Delta#eta"};
    // Eta correlations
    registry.add("global-muons/EtaCorrelation_Global_vs_Muon",
                 "#eta_{global} vs. #eta_{MCH}",
                 {HistType::kTH2F, {etaCorrelationAxis, etaCorrelationAxis}});
    registry.add("global-muons/EtaDifference_Global_vs_Muon",
                 "(#eta_{global} - #eta_{MCH}) / #eta_{MCH} vs. #eta_{MCH}",
                 {HistType::kTH2F, {etaCorrelationAxis, etaDeltaAxis}});
    registry.add("global-muons/EtaCorrelation_subleading_vs_leading",
                 "#eta_{subleading_match} vs. #eta_{leading_match}",
                 {HistType::kTH2F, {etaCorrelationAxis, etaCorrelationAxis}});
    registry.add("global-muons/EtaDifference_subleading_vs_leading",
                 "(#eta_{subleading_match} - #eta_{leading_match}) / #eta_{leading_match} vs. #eta_{leading_match}",
                 {HistType::kTH2F, {etaCorrelationAxis, etaDeltaAxis}});
  }

  void CreateDetailedHistograms()
  {
    AxisSpec dcaxMFTAxis = {400, -0.5, 0.5, "DCA_{x} (cm)"};
    AxisSpec dcayMFTAxis = {400, -0.5, 0.5, "DCA_{y} (cm)"};
    AxisSpec dcaxyMFTAxis = {400, 0.0, 0.75, "DCA_{xy} (cm)"};
    AxisSpec dcaxMCHAxis = {400, -10.0, 10.0, "DCA_{x} (cm)"};
    AxisSpec dcayMCHAxis = {400, -10.0, 10.0, "DCA_{y} (cm)"};
    AxisSpec dcaxyMCHAxis = {400, 0.0, 15.0, "DCA_{xy} (cm)"};
    AxisSpec dcazAxis = {20, -10.0, 10.0, "DCA_{z} (cm)"};
    AxisSpec dxAxis = {600, -30.0, 30.0, "#Delta x (cm)"};
    AxisSpec dyAxis = {600, -30.0, 30.0, "#Delta y (cm)"};
    AxisSpec thetaxAxis = {10, 0.0, 20.0, "#theta_{x} (degrees)"};
    AxisSpec dThetaxAxis = {500, -5.0, 5.0, "#Delta#theta_{x} (degrees)"};
    AxisSpec thetayAxis = {10, 0.0, 20.0, "#theta_{y} (degrees)"};
    AxisSpec dThetayAxis = {500, -5.0, 5.0, "#Delta#theta_{y} (degrees)"};
    AxisSpec phiAxis = {360, -180.0, 180.0, "#phi (degrees)"};
    AxisSpec dPhiAxis = {200, -20.0, 20.0, "#Delta#phi (degrees)"};

    if (configQAs.fEnableQAResidual) {
      for (size_t i = 0; i < referencePlanes.size(); i++) {
        const auto& refPLane = referencePlanes[i];
        AxisSpec xAxis = {10, 0, refPLane.second, "|x| (cm)"};
        AxisSpec yAxis = {10, 0, refPLane.second, "|y| (cm)"};
        for (size_t j = 0; j < quadrants.size(); j++) {
          const auto& quadrant = quadrants[j];
          std::string histPath = std::string("Alignment/same-event/Residuals/ReferencePlanes/") + refPLane.first + "/" + quadrant + "/";
          trackResidualsHistos[i][j]["dx_vs_x"] = registry.add((histPath + "dx_vs_x").c_str(), std::format("#Delta x vs. |x| - {}", quadrant).c_str(), {HistType::kTH2F, {xAxis, dxAxis}});
          trackResidualsHistos[i][j]["dx_vs_y"] = registry.add((histPath + "dx_vs_y").c_str(), std::format("#Delta x vs. |y| - {}", quadrant).c_str(), {HistType::kTH2F, {yAxis, dxAxis}});
          trackResidualsHistos[i][j]["dy_vs_x"] = registry.add((histPath + "dy_vs_x").c_str(), std::format("#Delta y vs. |x| - {}", quadrant).c_str(), {HistType::kTH2F, {xAxis, dyAxis}});
          trackResidualsHistos[i][j]["dy_vs_y"] = registry.add((histPath + "dy_vs_y").c_str(), std::format("#Delta y vs. |y| - {}", quadrant).c_str(), {HistType::kTH2F, {yAxis, dyAxis}});

          trackResidualsHistos[i][j]["dthetax_vs_x"] = registry.add((histPath + "dthetax_vs_x").c_str(), std::format("#Delta #theta_x vs. |x| - {}", quadrant).c_str(), {HistType::kTH2F, {xAxis, dThetaxAxis}});
          trackResidualsHistos[i][j]["dthetax_vs_y"] = registry.add((histPath + "dthetax_vs_y").c_str(), std::format("#Delta #theta_x vs. |y| - {}", quadrant).c_str(), {HistType::kTH2F, {yAxis, dThetaxAxis}});
          trackResidualsHistos[i][j]["dthetax_vs_thetax"] = registry.add((histPath + "dthetax_vs_thetax").c_str(), std::format("#Delta #theta_x vs. |#theta_x| - {}", quadrant).c_str(), {HistType::kTH2F, {thetaxAxis, dThetaxAxis}});

          trackResidualsHistos[i][j]["dthetay_vs_x"] = registry.add((histPath + "dthetay_vs_x").c_str(), std::format("#Delta #theta_y vs. |x| - {}", quadrant).c_str(), {HistType::kTH2F, {xAxis, dThetayAxis}});
          trackResidualsHistos[i][j]["dthetay_vs_y"] = registry.add((histPath + "dthetay_vs_y").c_str(), std::format("#Delta #theta_y vs. |y| - {}", quadrant).c_str(), {HistType::kTH2F, {yAxis, dThetayAxis}});
          trackResidualsHistos[i][j]["dthetay_vs_thetay"] = registry.add((histPath + "dthetay_vs_thetay").c_str(), std::format("#Delta #theta_y vs. |#theta_y| - {}", quadrant).c_str(), {HistType::kTH2F, {thetayAxis, dThetayAxis}});

          // mixed events
          histPath = std::string("Alignment/mixed-event/Residuals/ReferencePlanes/") + refPLane.first + "/" + quadrant + "/";
          trackResidualsHistosMixedEvents[i][j]["dx_vs_x"] = registry.add((histPath + "dx_vs_x").c_str(), std::format("#Delta x vs. |x| - {}", quadrant).c_str(), {HistType::kTH2F, {xAxis, dxAxis}});
          trackResidualsHistosMixedEvents[i][j]["dx_vs_y"] = registry.add((histPath + "dx_vs_y").c_str(), std::format("#Delta x vs. |y| - {}", quadrant).c_str(), {HistType::kTH2F, {yAxis, dxAxis}});
          trackResidualsHistosMixedEvents[i][j]["dy_vs_x"] = registry.add((histPath + "dy_vs_x").c_str(), std::format("#Delta y vs. |x| - {}", quadrant).c_str(), {HistType::kTH2F, {xAxis, dyAxis}});
          trackResidualsHistosMixedEvents[i][j]["dy_vs_y"] = registry.add((histPath + "dy_vs_y").c_str(), std::format("#Delta y vs. |y| - {}", quadrant).c_str(), {HistType::kTH2F, {yAxis, dyAxis}});

          trackResidualsHistosMixedEvents[i][j]["dthetax_vs_x"] = registry.add((histPath + "dthetax_vs_x").c_str(), std::format("#Delta #theta_x vs. |x| - {}", quadrant).c_str(), {HistType::kTH2F, {xAxis, dThetaxAxis}});
          trackResidualsHistosMixedEvents[i][j]["dthetax_vs_y"] = registry.add((histPath + "dthetax_vs_y").c_str(), std::format("#Delta #theta_x vs. |y| - {}", quadrant).c_str(), {HistType::kTH2F, {yAxis, dThetaxAxis}});
          trackResidualsHistosMixedEvents[i][j]["dthetax_vs_thetax"] = registry.add((histPath + "dthetax_vs_thetax").c_str(), std::format("#Delta #theta_x vs. |#theta_x| - {}", quadrant).c_str(), {HistType::kTH2F, {thetaxAxis, dThetaxAxis}});

          trackResidualsHistosMixedEvents[i][j]["dthetay_vs_x"] = registry.add((histPath + "dthetay_vs_x").c_str(), std::format("#Delta #theta_y vs. |x| - {}", quadrant).c_str(), {HistType::kTH2F, {xAxis, dThetayAxis}});
          trackResidualsHistosMixedEvents[i][j]["dthetay_vs_y"] = registry.add((histPath + "dthetay_vs_y").c_str(), std::format("#Delta #theta_y vs. |y| - {}", quadrant).c_str(), {HistType::kTH2F, {yAxis, dThetayAxis}});
          trackResidualsHistosMixedEvents[i][j]["dthetay_vs_thetay"] = registry.add((histPath + "dthetay_vs_thetay").c_str(), std::format("#Delta #theta_y vs. |#theta_y| - {}", quadrant).c_str(), {HistType::kTH2F, {thetayAxis, dThetayAxis}});
        }
      }

      for (size_t j = 0; j < quadrants.size(); j++) {
        const auto& quadrant = quadrants[j];
        AxisSpec xAxis = {20, 0, 200, "|x| (cm)"};
        AxisSpec yAxis = {10, 0, 200, "|y| (cm)"};
        for (int chamber = 0; chamber < 10; chamber++) {
          std::string histPath = std::string("Alignment/same-event/Residuals/MFT/") + quadrant + "/CH" + std::to_string(chamber + 1) + "/";
          // Delta x at cluster
          residualsHistos[j][chamber]["dx_vs_x"] = registryResiduals.add((histPath + "dx_vs_x").c_str(), "Cluster x residual vs. x", {HistType::kTH2F, {xAxis, dxAxis}});
          residualsHistos[j][chamber]["dx_vs_y"] = registryResiduals.add((histPath + "dx_vs_y").c_str(), "Cluster x residual vs. y", {HistType::kTH2F, {yAxis, dxAxis}});
          residualsHistos[j][chamber]["dy_vs_x"] = registryResiduals.add((histPath + "dy_vs_x").c_str(), "Cluster y residual vs. x", {HistType::kTH2F, {xAxis, dyAxis}});
          residualsHistos[j][chamber]["dy_vs_y"] = registryResiduals.add((histPath + "dy_vs_y").c_str(), "Cluster y residual vs. y", {HistType::kTH2F, {yAxis, dyAxis}});

          // mixed events
          histPath = std::string("Alignment/mixed-event/Residuals/MFT/") + quadrant + "/CH" + std::to_string(chamber + 1) + "/";
          // Delta x at cluster
          residualsHistosMixedEvents[j][chamber]["dx_vs_x"] = registryResiduals.add((histPath + "dx_vs_x").c_str(), "Cluster x residual vs. x", {HistType::kTH2F, {xAxis, dxAxis}});
          residualsHistosMixedEvents[j][chamber]["dx_vs_y"] = registryResiduals.add((histPath + "dx_vs_y").c_str(), "Cluster x residual vs. y", {HistType::kTH2F, {yAxis, dxAxis}});
          residualsHistosMixedEvents[j][chamber]["dy_vs_x"] = registryResiduals.add((histPath + "dy_vs_x").c_str(), "Cluster y residual vs. x", {HistType::kTH2F, {xAxis, dyAxis}});
          residualsHistosMixedEvents[j][chamber]["dy_vs_y"] = registryResiduals.add((histPath + "dy_vs_y").c_str(), "Cluster y residual vs. y", {HistType::kTH2F, {yAxis, dyAxis}});
        }
      }

      for (size_t i = 0; i < 2; i++) {
        std::string topBottom = (i == 0) ? "top" : "bottom";
        AxisSpec deAxis = {26, 0, 26, "DE index"};
        AxisSpec phiAxis = {16, -180, 180, "#phi (degrees)"};
        for (size_t j = 0; j < 2; j++) {
          std::string sign = (j == 0) ? "positive" : "negative";
          for (int chamber = 0; chamber < 10; chamber++) {
            std::string histPath = std::string("Alignment/same-event/Residuals/MFT/MFT_") + topBottom + "/" + sign + "/CH" + std::to_string(chamber + 1) + "/";
            // Delta x and y at cluster
            residualsHistosPerDE[i][j][chamber]["dx_vs_de"] = registryResidualsMFT.add((histPath + "dx_vs_de").c_str(), "Cluster x residual vs. DE index", {HistType::kTH2F, {deAxis, dxAxis}});
            residualsHistosPerDE[i][j][chamber]["dy_vs_de"] = registryResidualsMFT.add((histPath + "dy_vs_de").c_str(), "Cluster y residual vs. DE index", {HistType::kTH2F, {deAxis, dxAxis}});

            residualsHistosPerDE[i][j][chamber]["dx_vs_phi"] = registryResidualsMFT.add((histPath + "dx_vs_phi").c_str(), "Cluster x residual vs. cluster #phi", {HistType::kTH2F, {phiAxis, dxAxis}});
            residualsHistosPerDE[i][j][chamber]["dy_vs_phi"] = registryResidualsMFT.add((histPath + "dy_vs_phi").c_str(), "Cluster y residual vs. cluster #phi", {HistType::kTH2F, {phiAxis, dxAxis}});

            // mixed events
            histPath = std::string("Alignment/mixed-event/Residuals/MFT/MFT_") + topBottom + "/" + sign + "/CH" + std::to_string(chamber + 1) + "/";
            // Delta x and y at cluster
            residualsHistosPerDEMixedEvents[i][j][chamber]["dx_vs_de"] = registryResidualsMFT.add((histPath + "dx_vs_de").c_str(), "Cluster x residual vs. DE index", {HistType::kTH2F, {deAxis, dxAxis}});
            residualsHistosPerDEMixedEvents[i][j][chamber]["dy_vs_de"] = registryResidualsMFT.add((histPath + "dy_vs_de").c_str(), "Cluster y residual vs. DE index", {HistType::kTH2F, {deAxis, dxAxis}});

            residualsHistosPerDEMixedEvents[i][j][chamber]["dx_vs_phi"] = registryResidualsMFT.add((histPath + "dx_vs_phi").c_str(), "Cluster x residual vs. cluster #phi", {HistType::kTH2F, {phiAxis, dxAxis}});
            residualsHistosPerDEMixedEvents[i][j][chamber]["dy_vs_phi"] = registryResidualsMFT.add((histPath + "dy_vs_phi").c_str(), "Cluster y residual vs. cluster #phi", {HistType::kTH2F, {phiAxis, dxAxis}});
          }
        }
      }

      for (size_t i = 0; i < 2; i++) {
        std::string topBottom = (i == 0) ? "top" : "bottom";
        AxisSpec deAxis = {26, 0, 26, "DE index"};
        for (size_t j = 0; j < 2; j++) {
          std::string sign = (j == 0) ? "positive" : "negative";
          for (int chamber = 0; chamber < 10; chamber++) {
            std::string histPath = std::string("Alignment/same-event/Residuals/MCH/MCH_") + topBottom + "/" + sign + "/CH" + std::to_string(chamber + 1) + "/";
            // Delta x and y at cluster
            mchResidualsHistosPerDE[i][j][chamber]["dx_vs_de"] = registryResidualsMCH.add((histPath + "dx_vs_de").c_str(), "Cluster x residual vs. DE index", {HistType::kTH2F, {deAxis, dxAxis}});
            mchResidualsHistosPerDE[i][j][chamber]["dy_vs_de"] = registryResidualsMCH.add((histPath + "dy_vs_de").c_str(), "Cluster y residual vs. DE index", {HistType::kTH2F, {deAxis, dxAxis}});

            // mixed events
            histPath = std::string("Alignment/mixed-event/Residuals/MCH/MCH_") + topBottom + "/" + sign + "/CH" + std::to_string(chamber + 1) + "/";
            // Delta x and y at cluster
            mchResidualsHistosPerDEMixedEvents[i][j][chamber]["dx_vs_de"] = registryResidualsMCH.add((histPath + "dx_vs_de").c_str(), "Cluster x residual vs. DE index", {HistType::kTH2F, {deAxis, dxAxis}});
            mchResidualsHistosPerDEMixedEvents[i][j][chamber]["dy_vs_de"] = registryResidualsMCH.add((histPath + "dy_vs_de").c_str(), "Cluster y residual vs. DE index", {HistType::kTH2F, {deAxis, dxAxis}});
          }
        }
      }
    }

    if (configQAs.fEnableQADCA) {
      for (size_t j = 0; j < quadrants.size(); j++) {
        const auto& quadrant = quadrants[j];
        std::string histPath = std::string("Alignment/same-event/DCA/MFT/") + quadrant + "/";
        dcaHistos[0][j][0]["DCA_x"] = registryDCA.add((histPath + "DCA_x").c_str(), std::format("DCA(x) - {}", quadrant).c_str(), {HistType::kTH1F, {dcaxMFTAxis}});
        dcaHistos[0][j][1]["DCA_x"] = registryDCA.add((histPath + "DCA_x_pos").c_str(), std::format("DCA(x) - {} charge > 0", quadrant).c_str(), {HistType::kTH1F, {dcaxMFTAxis}});
        dcaHistos[0][j][2]["DCA_x"] = registryDCA.add((histPath + "DCA_x_neg").c_str(), std::format("DCA(x) - {} charge < 0", quadrant).c_str(), {HistType::kTH1F, {dcaxMFTAxis}});
        dcaHistos[0][j][0]["DCA_y"] = registryDCA.add((histPath + "DCA_y").c_str(), std::format("DCA(y) - {}", quadrant).c_str(), {HistType::kTH1F, {dcayMFTAxis}});
        dcaHistos[0][j][1]["DCA_y"] = registryDCA.add((histPath + "DCA_y_pos").c_str(), std::format("DCA(y) - {} charge > 0", quadrant).c_str(), {HistType::kTH1F, {dcayMFTAxis}});
        dcaHistos[0][j][2]["DCA_y"] = registryDCA.add((histPath + "DCA_y_neg").c_str(), std::format("DCA(y) - {} charge < 0", quadrant).c_str(), {HistType::kTH1F, {dcayMFTAxis}});
        dcaHistos[0][j][0]["DCA_xy"] = registryDCA.add((histPath + "DCA_xy").c_str(), std::format("DCA(xy) - {}", quadrant).c_str(), {HistType::kTH1F, {dcaxyMFTAxis}});
        dcaHistos[0][j][1]["DCA_xy"] = registryDCA.add((histPath + "DCA_xy_pos").c_str(), std::format("DCA(xy) - {} charge > 0", quadrant).c_str(), {HistType::kTH1F, {dcaxyMFTAxis}});
        dcaHistos[0][j][2]["DCA_xy"] = registryDCA.add((histPath + "DCA_xy_neg").c_str(), std::format("DCA(xy) - {} charge < 0", quadrant).c_str(), {HistType::kTH1F, {dcaxyMFTAxis}});
        dcaHistos[0][j][0]["DCA_x_vs_DCA_y"] = registryDCA.add((histPath + "DCA_y_vs_DCA_x").c_str(), std::format("DCA(x) vs. DCA(y) - {}", quadrant).c_str(), {HistType::kTH2F, {dcaxMFTAxis, dcayMFTAxis}});
        dcaHistos[0][j][0]["DCA_x_vs_z"] = registryDCA.add((histPath + "DCA_x_vs_z").c_str(), std::format("DCA(x) vs. z - {}", quadrant).c_str(), {HistType::kTH2F, {dcazAxis, dcaxMFTAxis}});
        dcaHistos[0][j][0]["DCA_y_vs_z"] = registryDCA.add((histPath + "DCA_y_vs_z").c_str(), std::format("DCA(y) vs. z - {}", quadrant).c_str(), {HistType::kTH2F, {dcazAxis, dcayMFTAxis}});
        dcaHistos[0][j][0]["DCA_xy_vs_z"] = registryDCA.add((histPath + "DCA_xy_vs_z").c_str(), std::format("DCA(xy) vs. z - {}", quadrant).c_str(), {HistType::kTH2F, {dcazAxis, dcaxyMFTAxis}});
        dcaHistos[0][j][0]["DCA_x_vs_phi"] = registryDCA.add((histPath + "DCA_x_vs_phi").c_str(), std::format("DCA(x) vs. phi - {}", quadrant).c_str(), {HistType::kTH2F, {phiAxis, dcaxMFTAxis}});
        dcaHistos[0][j][0]["DCA_y_vs_phi"] = registryDCA.add((histPath + "DCA_y_vs_phi").c_str(), std::format("DCA(y) vs. phi - {}", quadrant).c_str(), {HistType::kTH2F, {phiAxis, dcayMFTAxis}});
        dcaHistos[0][j][0]["DCA_xy_vs_phi"] = registryDCA.add((histPath + "DCA_xy_vs_phi").c_str(), std::format("DCA(xy) vs. phi - {}", quadrant).c_str(), {HistType::kTH2F, {phiAxis, dcaxyMFTAxis}});

        histPath = std::string("Alignment/same-event/DCA/MCH/") + quadrant + "/";
        dcaHistos[1][j][0]["DCA_x"] = registryDCA.add((histPath + "DCA_x").c_str(), std::format("DCA(x) - {}", quadrant).c_str(), {HistType::kTH1F, {dcaxMCHAxis}});
        dcaHistos[1][j][1]["DCA_x"] = registryDCA.add((histPath + "DCA_x_pos").c_str(), std::format("DCA(x) - {} charge > 0", quadrant).c_str(), {HistType::kTH1F, {dcaxMCHAxis}});
        dcaHistos[1][j][2]["DCA_x"] = registryDCA.add((histPath + "DCA_x_neg").c_str(), std::format("DCA(x) - {} charge < 0", quadrant).c_str(), {HistType::kTH1F, {dcaxMCHAxis}});
        dcaHistos[1][j][0]["DCA_y"] = registryDCA.add((histPath + "DCA_y").c_str(), std::format("DCA(y) - {}", quadrant).c_str(), {HistType::kTH1F, {dcayMCHAxis}});
        dcaHistos[1][j][1]["DCA_y"] = registryDCA.add((histPath + "DCA_y_pos").c_str(), std::format("DCA(y) - {} charge > 0", quadrant).c_str(), {HistType::kTH1F, {dcayMCHAxis}});
        dcaHistos[1][j][2]["DCA_y"] = registryDCA.add((histPath + "DCA_y_neg").c_str(), std::format("DCA(y) - {} charge < 0", quadrant).c_str(), {HistType::kTH1F, {dcayMCHAxis}});
        dcaHistos[1][j][0]["DCA_xy"] = registryDCA.add((histPath + "DCA_xy").c_str(), std::format("DCA(xy) - {}", quadrant).c_str(), {HistType::kTH1F, {dcaxyMCHAxis}});
        dcaHistos[1][j][1]["DCA_xy"] = registryDCA.add((histPath + "DCA_xy_pos").c_str(), std::format("DCA(xy) - {} charge > 0", quadrant).c_str(), {HistType::kTH1F, {dcaxyMCHAxis}});
        dcaHistos[1][j][2]["DCA_xy"] = registryDCA.add((histPath + "DCA_xy_neg").c_str(), std::format("DCA(xy) - {} charge < 0", quadrant).c_str(), {HistType::kTH1F, {dcaxyMCHAxis}});
        dcaHistos[1][j][0]["DCA_x_vs_DCA_y"] = registryDCA.add((histPath + "DCA_y_vs_DCA_x").c_str(), std::format("DCA(x) vs. DCA(y) - {}", quadrant).c_str(), {HistType::kTH2F, {dcaxMCHAxis, dcayMCHAxis}});
        // dcaHistos[1][j][0]["DCA_x_vs_z"] = registryDCA.add((histPath + "DCA_x_vs_z").c_str(), std::format("DCA(x) vs. z - {}", quadrant).c_str(), {HistType::kTH2F, {dcazAxis, dcaxMCHAxis}});
        // dcaHistos[1][j][0]["DCA_y_vs_z"] = registryDCA.add((histPath + "DCA_y_vs_z").c_str(), std::format("DCA(y) vs. z - {}", quadrant).c_str(), {HistType::kTH2F, {dcazAxis, dcayMCHAxis}});
        dcaHistos[1][j][0]["DCA_x_vs_phi"] = registryDCA.add((histPath + "DCA_x_vs_phi").c_str(), std::format("DCA(x) vs. phi - {}", quadrant).c_str(), {HistType::kTH2F, {phiAxis, dcaxMCHAxis}});
        dcaHistos[1][j][0]["DCA_y_vs_phi"] = registryDCA.add((histPath + "DCA_y_vs_phi").c_str(), std::format("DCA(y) vs. phi - {}", quadrant).c_str(), {HistType::kTH2F, {phiAxis, dcayMCHAxis}});
        dcaHistos[1][j][0]["DCA_xy_vs_phi"] = registryDCA.add((histPath + "DCA_xy_vs_phi").c_str(), std::format("DCA(xy) vs. phi - {}", quadrant).c_str(), {HistType::kTH2F, {phiAxis, dcaxyMCHAxis}});

        histPath = std::string("Alignment/same-event/DCA/GlobalMuons/MFT/") + quadrant + "/";
        dcaHistosGlobal[0][j][0]["DCA_x"] = registryDCA.add((histPath + "DCA_x").c_str(), std::format("DCA(x) - {}", quadrant).c_str(), {HistType::kTH1F, {dcaxMFTAxis}});
        dcaHistosGlobal[0][j][1]["DCA_x"] = registryDCA.add((histPath + "DCA_x_pos").c_str(), std::format("DCA(x) - {} charge > 0", quadrant).c_str(), {HistType::kTH1F, {dcaxMFTAxis}});
        dcaHistosGlobal[0][j][2]["DCA_x"] = registryDCA.add((histPath + "DCA_x_neg").c_str(), std::format("DCA(x) - {} charge < 0", quadrant).c_str(), {HistType::kTH1F, {dcaxMFTAxis}});
        dcaHistosGlobal[0][j][0]["DCA_y"] = registryDCA.add((histPath + "DCA_y").c_str(), std::format("DCA(y) - {}", quadrant).c_str(), {HistType::kTH1F, {dcayMFTAxis}});
        dcaHistosGlobal[0][j][1]["DCA_y"] = registryDCA.add((histPath + "DCA_y_pos").c_str(), std::format("DCA(y) - {} charge > 0", quadrant).c_str(), {HistType::kTH1F, {dcayMFTAxis}});
        dcaHistosGlobal[0][j][2]["DCA_y"] = registryDCA.add((histPath + "DCA_y_neg").c_str(), std::format("DCA(y) - {} charge < 0", quadrant).c_str(), {HistType::kTH1F, {dcayMFTAxis}});
        dcaHistosGlobal[0][j][0]["DCA_xy"] = registryDCA.add((histPath + "DCA_xy").c_str(), std::format("DCA(xy) - {}", quadrant).c_str(), {HistType::kTH1F, {dcaxyMFTAxis}});
        dcaHistosGlobal[0][j][1]["DCA_xy"] = registryDCA.add((histPath + "DCA_xy_pos").c_str(), std::format("DCA(xy) - {} charge > 0", quadrant).c_str(), {HistType::kTH1F, {dcaxyMFTAxis}});
        dcaHistosGlobal[0][j][2]["DCA_xy"] = registryDCA.add((histPath + "DCA_xy_neg").c_str(), std::format("DCA(xy) - {} charge < 0", quadrant).c_str(), {HistType::kTH1F, {dcaxyMFTAxis}});
        dcaHistosGlobal[0][j][0]["DCA_x_vs_DCA_y"] = registryDCA.add((histPath + "DCA_y_vs_DCA_x").c_str(), std::format("DCA(x) vs. DCA(y) - {}", quadrant).c_str(), {HistType::kTH2F, {dcaxMFTAxis, dcayMFTAxis}});
        dcaHistosGlobal[0][j][0]["DCA_x_vs_z"] = registryDCA.add((histPath + "DCA_x_vs_z").c_str(), std::format("DCA(x) vs. z - {}", quadrant).c_str(), {HistType::kTH2F, {dcazAxis, dcaxMFTAxis}});
        dcaHistosGlobal[0][j][0]["DCA_y_vs_z"] = registryDCA.add((histPath + "DCA_y_vs_z").c_str(), std::format("DCA(y) vs. z - {}", quadrant).c_str(), {HistType::kTH2F, {dcazAxis, dcayMFTAxis}});
        dcaHistosGlobal[0][j][0]["DCA_xy_vs_z"] = registryDCA.add((histPath + "DCA_xy_vs_z").c_str(), std::format("DCA(xy) vs. z - {}", quadrant).c_str(), {HistType::kTH2F, {dcazAxis, dcaxyMFTAxis}});
        dcaHistosGlobal[0][j][0]["DCA_x_vs_phi"] = registryDCA.add((histPath + "DCA_x_vs_phi").c_str(), std::format("DCA(x) vs. phi - {}", quadrant).c_str(), {HistType::kTH2F, {phiAxis, dcaxMFTAxis}});
        dcaHistosGlobal[0][j][0]["DCA_y_vs_phi"] = registryDCA.add((histPath + "DCA_y_vs_phi").c_str(), std::format("DCA(y) vs. phi - {}", quadrant).c_str(), {HistType::kTH2F, {phiAxis, dcayMFTAxis}});
        dcaHistosGlobal[0][j][0]["DCA_xy_vs_phi"] = registryDCA.add((histPath + "DCA_xy_vs_phi").c_str(), std::format("DCA(xy) vs. phi - {}", quadrant).c_str(), {HistType::kTH2F, {phiAxis, dcaxyMFTAxis}});

        histPath = std::string("Alignment/same-event/DCA/GlobalMuons/MCH/") + quadrant + "/";
        dcaHistosGlobal[1][j][0]["DCA_x"] = registryDCA.add((histPath + "DCA_x").c_str(), std::format("DCA(x) - {}", quadrant).c_str(), {HistType::kTH1F, {dcaxMCHAxis}});
        dcaHistosGlobal[1][j][1]["DCA_x"] = registryDCA.add((histPath + "DCA_x_pos").c_str(), std::format("DCA(x) - {} charge > 0", quadrant).c_str(), {HistType::kTH1F, {dcaxMCHAxis}});
        dcaHistosGlobal[1][j][2]["DCA_x"] = registryDCA.add((histPath + "DCA_x_neg").c_str(), std::format("DCA(x) - {} charge < 0", quadrant).c_str(), {HistType::kTH1F, {dcaxMCHAxis}});
        dcaHistosGlobal[1][j][0]["DCA_y"] = registryDCA.add((histPath + "DCA_y").c_str(), std::format("DCA(y) - {}", quadrant).c_str(), {HistType::kTH1F, {dcayMCHAxis}});
        dcaHistosGlobal[1][j][1]["DCA_y"] = registryDCA.add((histPath + "DCA_y_pos").c_str(), std::format("DCA(y) - {} charge > 0", quadrant).c_str(), {HistType::kTH1F, {dcayMCHAxis}});
        dcaHistosGlobal[1][j][2]["DCA_y"] = registryDCA.add((histPath + "DCA_y_neg").c_str(), std::format("DCA(y) - {} charge < 0", quadrant).c_str(), {HistType::kTH1F, {dcayMCHAxis}});
        dcaHistosGlobal[1][j][0]["DCA_xy"] = registryDCA.add((histPath + "DCA_xy").c_str(), std::format("DCA(xy) - {}", quadrant).c_str(), {HistType::kTH1F, {dcaxyMCHAxis}});
        dcaHistosGlobal[1][j][1]["DCA_xy"] = registryDCA.add((histPath + "DCA_xy_pos").c_str(), std::format("DCA(xy) - {} charge > 0", quadrant).c_str(), {HistType::kTH1F, {dcaxyMCHAxis}});
        dcaHistosGlobal[1][j][2]["DCA_xy"] = registryDCA.add((histPath + "DCA_xy_neg").c_str(), std::format("DCA(xy) - {} charge < 0", quadrant).c_str(), {HistType::kTH1F, {dcaxyMCHAxis}});
        dcaHistosGlobal[1][j][0]["DCA_x_vs_DCA_y"] = registryDCA.add((histPath + "DCA_y_vs_DCA_x").c_str(), std::format("DCA(x) vs. DCA(y) - {}", quadrant).c_str(), {HistType::kTH2F, {dcaxMCHAxis, dcayMCHAxis}});
        dcaHistosGlobal[1][j][0]["DCA_x_vs_z"] = registryDCA.add((histPath + "DCA_x_vs_z").c_str(), std::format("DCA(x) vs. z - {}", quadrant).c_str(), {HistType::kTH2F, {dcazAxis, dcaxMCHAxis}});
        dcaHistosGlobal[1][j][0]["DCA_y_vs_z"] = registryDCA.add((histPath + "DCA_y_vs_z").c_str(), std::format("DCA(y) vs. z - {}", quadrant).c_str(), {HistType::kTH2F, {dcazAxis, dcayMCHAxis}});
        dcaHistosGlobal[1][j][0]["DCA_x_vs_phi"] = registryDCA.add((histPath + "DCA_x_vs_phi").c_str(), std::format("DCA(x) vs. phi - {}", quadrant).c_str(), {HistType::kTH2F, {phiAxis, dcaxMCHAxis}});
        dcaHistosGlobal[1][j][0]["DCA_y_vs_phi"] = registryDCA.add((histPath + "DCA_y_vs_phi").c_str(), std::format("DCA(y) vs. phi - {}", quadrant).c_str(), {HistType::kTH2F, {phiAxis, dcayMCHAxis}});
        dcaHistosGlobal[1][j][0]["DCA_xy_vs_phi"] = registryDCA.add((histPath + "DCA_xy_vs_phi").c_str(), std::format("DCA(xy) vs. phi - {}", quadrant).c_str(), {HistType::kTH2F, {phiAxis, dcaxyMCHAxis}});

        histPath = std::string("Alignment/same-event/DCA/GlobalMuons/MFT-MCH/") + quadrant + "/";
        dcaHistosGlobalSubtracted[j][0]["DCA_x"] = registryDCA.add((histPath + "DCA_x").c_str(), std::format("scaled MFT DCA(x) - MCH DCA(x) - {}", quadrant).c_str(), {HistType::kTH1F, {dcaxMCHAxis}});
        dcaHistosGlobalSubtracted[j][1]["DCA_x"] = registryDCA.add((histPath + "DCA_x_pos").c_str(), std::format("scaled MFT DCA(x) - MCH DCA(x) - {} charge > 0", quadrant).c_str(), {HistType::kTH1F, {dcaxMCHAxis}});
        dcaHistosGlobalSubtracted[j][2]["DCA_x"] = registryDCA.add((histPath + "DCA_x_neg").c_str(), std::format("scaled MFT DCA(x) - MCH DCA(x) - {} charge < 0", quadrant).c_str(), {HistType::kTH1F, {dcaxMCHAxis}});
        dcaHistosGlobalSubtracted[j][0]["DCA_y"] = registryDCA.add((histPath + "DCA_y").c_str(), std::format("scaled MFT DCA(y) - MCH DCA(y) - {}", quadrant).c_str(), {HistType::kTH1F, {dcayMCHAxis}});
        dcaHistosGlobalSubtracted[j][1]["DCA_y"] = registryDCA.add((histPath + "DCA_y_pos").c_str(), std::format("scaled MFT DCA(y) - MCH DCA(y) - {} charge > 0", quadrant).c_str(), {HistType::kTH1F, {dcayMCHAxis}});
        dcaHistosGlobalSubtracted[j][2]["DCA_y"] = registryDCA.add((histPath + "DCA_y_neg").c_str(), std::format("scaled MFT DCA(y) - MCH DCA(y) - {} charge < 0", quadrant).c_str(), {HistType::kTH1F, {dcayMCHAxis}});
        dcaHistosGlobalSubtracted[j][0]["DCA_xy"] = registryDCA.add((histPath + "DCA_xy").c_str(), std::format("scaled MFT DCA(xy) - MCH DCA(xy) - {}", quadrant).c_str(), {HistType::kTH1F, {dcaxyMCHAxis}});
        dcaHistosGlobalSubtracted[j][1]["DCA_xy"] = registryDCA.add((histPath + "DCA_xy_pos").c_str(), std::format("scaled MFT DCA(xy) - MCH DCA(xy) - {} charge > 0", quadrant).c_str(), {HistType::kTH1F, {dcaxyMCHAxis}});
        dcaHistosGlobalSubtracted[j][2]["DCA_xy"] = registryDCA.add((histPath + "DCA_xy_neg").c_str(), std::format("scaled MFT DCA(xy) - MCH DCA(xy) - {} charge < 0", quadrant).c_str(), {HistType::kTH1F, {dcaxyMCHAxis}});
        dcaHistosGlobalSubtracted[j][0]["DCA_x_vs_DCA_y"] = registryDCA.add((histPath + "DCA_x_vs_DCA_y").c_str(), std::format("scaled MFT DCA(x) - MCH DCA(x) - vs. scaled MFT DCA(y) - MCH DCA(y) - {}", quadrant).c_str(), {HistType::kTH2F, {dcaxMCHAxis, dcayMCHAxis}});
        dcaHistosGlobalSubtracted[j][0]["DCA_x_vs_z"] = registryDCA.add((histPath + "DCA_x_vs_z").c_str(), std::format("scaled MFT DCA(x) - MCH DCA(x) - vs. z - {}", quadrant).c_str(), {HistType::kTH2F, {dcazAxis, dcaxMCHAxis}});
        dcaHistosGlobalSubtracted[j][0]["DCA_y_vs_z"] = registryDCA.add((histPath + "DCA_y_vs_z").c_str(), std::format("scaled MFT DCA(y) - MCH DCA(y) - vs. z - {}", quadrant).c_str(), {HistType::kTH2F, {dcazAxis, dcayMCHAxis}});
        dcaHistosGlobalSubtracted[j][0]["DCA_xy_vs_z"] = registryDCA.add((histPath + "DCA_xy_vs_z").c_str(), std::format("scaled MFT DCA(xy) - MCH DCA(xy) - vs. z - {}", quadrant).c_str(), {HistType::kTH2F, {dcazAxis, dcaxyMCHAxis}});
        dcaHistosGlobalSubtracted[j][0]["DCA_x_vs_phi"] = registryDCA.add((histPath + "DCA_x_vs_phi").c_str(), std::format("scaled MFT DCA(x) - MCH DCA(x) - vs. phi - {}", quadrant).c_str(), {HistType::kTH2F, {phiAxis, dcaxMCHAxis}});
        dcaHistosGlobalSubtracted[j][0]["DCA_y_vs_phi"] = registryDCA.add((histPath + "DCA_y_vs_phi").c_str(), std::format("scaled MFT DCA(y) - MCH DCA(y) - vs. phi - {}", quadrant).c_str(), {HistType::kTH2F, {phiAxis, dcayMCHAxis}});
        dcaHistosGlobalSubtracted[j][0]["DCA_xy_vs_phi"] = registryDCA.add((histPath + "DCA_xy_vs_phi").c_str(), std::format("scaled MFT DCA(xy) - MCH DCA(xy) - vs. phi - {}", quadrant).c_str(), {HistType::kTH2F, {phiAxis, dcaxyMCHAxis}});

        histPath = std::string("Alignment/same-event/DCA/GlobalMuons/MCHpCut/") + quadrant + "/";
        dcaHistosGlobalSubtractedMCHpCut[j][0]["DCA_x"] = registryDCA.add((histPath + "DCA_x").c_str(), std::format("scaled MFT DCA(x) - MCH DCA(x) - {} for MCH p > cut", quadrant).c_str(), {HistType::kTH1F, {dcaxMCHAxis}});
        dcaHistosGlobalSubtractedMCHpCut[j][1]["DCA_x"] = registryDCA.add((histPath + "DCA_x_pos").c_str(), std::format("scaled MFT DCA(x) - MCH DCA(x) - {} charge > 0 for MCH p > cut", quadrant).c_str(), {HistType::kTH1F, {dcaxMCHAxis}});
        dcaHistosGlobalSubtractedMCHpCut[j][2]["DCA_x"] = registryDCA.add((histPath + "DCA_x_neg").c_str(), std::format("scaled MFT DCA(x) - MCH DCA(x) - {} charge < 0 for MCH p > cut", quadrant).c_str(), {HistType::kTH1F, {dcaxMCHAxis}});
        dcaHistosGlobalSubtractedMCHpCut[j][0]["DCA_y"] = registryDCA.add((histPath + "DCA_y").c_str(), std::format("scaled MFT DCA(y) - MCH DCA(y) - {} for MCH p > cut", quadrant).c_str(), {HistType::kTH1F, {dcayMCHAxis}});
        dcaHistosGlobalSubtractedMCHpCut[j][1]["DCA_y"] = registryDCA.add((histPath + "DCA_y_pos").c_str(), std::format("scaled MFT DCA(y) - MCH DCA(y) - {} charge > 0 for MCH p > cut", quadrant).c_str(), {HistType::kTH1F, {dcayMCHAxis}});
        dcaHistosGlobalSubtractedMCHpCut[j][2]["DCA_y"] = registryDCA.add((histPath + "DCA_y_neg").c_str(), std::format("scaled MFT DCA(y) - MCH DCA(y) - {} charge < 0 for MCH p > cut", quadrant).c_str(), {HistType::kTH1F, {dcayMCHAxis}});
        dcaHistosGlobalSubtractedMCHpCut[j][0]["DCA_xy"] = registryDCA.add((histPath + "DCA_xy").c_str(), std::format("scaled MFT DCA(xy) - MCH DCA(xy) - {} for MCH p > cut", quadrant).c_str(), {HistType::kTH1F, {dcaxyMCHAxis}});
        dcaHistosGlobalSubtractedMCHpCut[j][1]["DCA_xy"] = registryDCA.add((histPath + "DCA_xy_pos").c_str(), std::format("scaled MFT DCA(xy) - MCH DCA(xy) - {} charge > 0 for MCH p > cut", quadrant).c_str(), {HistType::kTH1F, {dcaxyMCHAxis}});
        dcaHistosGlobalSubtractedMCHpCut[j][2]["DCA_xy"] = registryDCA.add((histPath + "DCA_xy_neg").c_str(), std::format("scaled MFT DCA(xy) - MCH DCA(xy) - {} charge < 0 for MCH p > cut", quadrant).c_str(), {HistType::kTH1F, {dcaxyMCHAxis}});
        dcaHistosGlobalSubtractedMCHpCut[j][0]["DCA_x_vs_DCA_y"] = registryDCA.add((histPath + "DCA_x_vs_DCA_y").c_str(), std::format("scaled MFT DCA(x) - MCH DCA(x) - vs. scaled MFT DCA(y) - MCH DCA(y) - {} for MCH p > cut", quadrant).c_str(), {HistType::kTH2F, {dcaxMCHAxis, dcayMCHAxis}});
        dcaHistosGlobalSubtractedMCHpCut[j][0]["DCA_x_vs_z"] = registryDCA.add((histPath + "DCA_x_vs_z").c_str(), std::format("scaled MFT DCA(x) - MCH DCA(x) - vs. z - {} for MCH p > cut", quadrant).c_str(), {HistType::kTH2F, {dcazAxis, dcaxMCHAxis}});
        dcaHistosGlobalSubtractedMCHpCut[j][0]["DCA_y_vs_z"] = registryDCA.add((histPath + "DCA_y_vs_z").c_str(), std::format("scaled MFT DCA(y) - MCH DCA(y) - vs. z - {} for MCH p > cut", quadrant).c_str(), {HistType::kTH2F, {dcazAxis, dcayMCHAxis}});
        dcaHistosGlobalSubtractedMCHpCut[j][0]["DCA_xy_vs_z"] = registryDCA.add((histPath + "DCA_xy_vs_z").c_str(), std::format("scaled MFT DCA(xy) - MCH DCA(xy) - vs. z - {} for MCH p > cut", quadrant).c_str(), {HistType::kTH2F, {dcazAxis, dcaxyMCHAxis}});
        dcaHistosGlobalSubtractedMCHpCut[j][0]["DCA_x_vs_phi"] = registryDCA.add((histPath + "DCA_x_vs_phi").c_str(), std::format("scaled MFT DCA(x) - MCH DCA(x) - vs. phi - {} for MCH p > cut", quadrant).c_str(), {HistType::kTH2F, {phiAxis, dcaxMCHAxis}});
        dcaHistosGlobalSubtractedMCHpCut[j][0]["DCA_y_vs_phi"] = registryDCA.add((histPath + "DCA_y_vs_phi").c_str(), std::format("scaled MFT DCA(y) - MCH DCA(y) - vs. phi - {} for MCH p > cut", quadrant).c_str(), {HistType::kTH2F, {phiAxis, dcayMCHAxis}});
        dcaHistosGlobalSubtractedMCHpCut[j][0]["DCA_xy_vs_phi"] = registryDCA.add((histPath + "DCA_xy_vs_phi").c_str(), std::format("scaled MFT DCA(xy) - MCH DCA(xy) - vs. phi - {} for MCH p > cut", quadrant).c_str(), {HistType::kTH2F, {phiAxis, dcaxyMCHAxis}});

        histPath = std::string("Alignment/mixed-event/DCA/MFT/") + quadrant + "/";
        dcaHistosMixedEvents[0][j][0]["DCA_x"] = registryDCA.add((histPath + "DCA_x").c_str(), std::format("DCA(x) - {}", quadrant).c_str(), {HistType::kTH1F, {dcaxMFTAxis}});
        dcaHistosMixedEvents[0][j][1]["DCA_x"] = registryDCA.add((histPath + "DCA_x_pos").c_str(), std::format("DCA(x) - {} charge > 0", quadrant).c_str(), {HistType::kTH1F, {dcaxMFTAxis}});
        dcaHistosMixedEvents[0][j][2]["DCA_x"] = registryDCA.add((histPath + "DCA_x_neg").c_str(), std::format("DCA(x) - {} charge < 0", quadrant).c_str(), {HistType::kTH1F, {dcaxMFTAxis}});
        dcaHistosMixedEvents[0][j][0]["DCA_y"] = registryDCA.add((histPath + "DCA_y").c_str(), std::format("DCA(y) - {}", quadrant).c_str(), {HistType::kTH1F, {dcayMFTAxis}});
        dcaHistosMixedEvents[0][j][1]["DCA_y"] = registryDCA.add((histPath + "DCA_y_pos").c_str(), std::format("DCA(y) - {} charge > 0", quadrant).c_str(), {HistType::kTH1F, {dcayMFTAxis}});
        dcaHistosMixedEvents[0][j][2]["DCA_y"] = registryDCA.add((histPath + "DCA_y_neg").c_str(), std::format("DCA(y) - {} charge < 0", quadrant).c_str(), {HistType::kTH1F, {dcayMFTAxis}});
        dcaHistosMixedEvents[0][j][0]["DCA_x_vs_z"] = registryDCA.add((histPath + "DCA_x_vs_z").c_str(), std::format("DCA(x) vs. z - {}", quadrant).c_str(), {HistType::kTH2F, {dcazAxis, dcaxMFTAxis}});
        dcaHistosMixedEvents[0][j][0]["DCA_y_vs_z"] = registryDCA.add((histPath + "DCA_y_vs_z").c_str(), std::format("DCA(y) vs. z - {}", quadrant).c_str(), {HistType::kTH2F, {dcazAxis, dcayMFTAxis}});

        histPath = std::string("Alignment/mixed-event/DCA/MCH/") + quadrant + "/";
        dcaHistosMixedEvents[1][j][0]["DCA_x"] = registryDCA.add((histPath + "DCA_x").c_str(), std::format("DCA(x) - {}", quadrant).c_str(), {HistType::kTH1F, {dcaxMCHAxis}});
        dcaHistosMixedEvents[1][j][1]["DCA_x"] = registryDCA.add((histPath + "DCA_x_pos").c_str(), std::format("DCA(x) - {} charge > 0", quadrant).c_str(), {HistType::kTH1F, {dcaxMCHAxis}});
        dcaHistosMixedEvents[1][j][2]["DCA_x"] = registryDCA.add((histPath + "DCA_x_neg").c_str(), std::format("DCA(x) - {} charge < 0", quadrant).c_str(), {HistType::kTH1F, {dcaxMCHAxis}});
        dcaHistosMixedEvents[1][j][0]["DCA_y"] = registryDCA.add((histPath + "DCA_y").c_str(), std::format("DCA(y) - {}", quadrant).c_str(), {HistType::kTH1F, {dcayMCHAxis}});
        dcaHistosMixedEvents[1][j][1]["DCA_y"] = registryDCA.add((histPath + "DCA_y_pos").c_str(), std::format("DCA(y) - {} charge > 0", quadrant).c_str(), {HistType::kTH1F, {dcayMCHAxis}});
        dcaHistosMixedEvents[1][j][2]["DCA_y"] = registryDCA.add((histPath + "DCA_y_neg").c_str(), std::format("DCA(y) - {} charge < 0", quadrant).c_str(), {HistType::kTH1F, {dcayMCHAxis}});
      }
    }

    if (configQAs.fEnableQADimuon) {
      // single muons
      AxisSpec transverseMomentumAxis = {100, 0, 30, "p_{T} (GeV/c)"};
      AxisSpec etaAxis = {40, -5, -1, "#eta"};
      AxisSpec rAbsAxis = {10, 0., 100.0, "R_{abs} (cm)"};
      AxisSpec dcaAxis = {400, -10.0, 10.0, "DCA"};
      AxisSpec dcaAxisReduced = {40, -10.0, 10.0, "DCA"};
      AxisSpec phiAxis = {36, -180.0, 180.0, "#phi (degrees)"};
      // dimuons
      AxisSpec invMassAxis = {400, 1, 5, "M_{#mu^{+}#mu^{-}} (GeV/c^{2})"};
      AxisSpec invMassCorrelationAxis = {80, 0, 8, "M_{#mu^{+}#mu^{-}} (GeV/c^{2})"};
      AxisSpec invMassAxisFull = {5000, 0, 100, "M_{#mu^{+}#mu^{-}} (GeV/c^{2})"};
      AxisSpec yPairAxis = {120, 0.0, 6.0, "#y_{pair}"};
      AxisSpec invMassAxis2D = {750, 0, 15, "M_{#mu^{+}#mu^{-}} (GeV/c^{2})"};
      AxisSpec pTAxis2D = {120, 0, 30, "p_{T} (GeV/c)"};
      // Jpsi candidate DCA histograms
      registryDimuon.add("dimuon/same-event/opposite-sign/DCA/pT_MuPosDCAx_minus_MuNegDCAx_MuonKine_MuonCuts", "DCA_x #mu^{+} minus DCA_x #mu^{-} and #mu^{+}#mu^{-} p_{T}", {HistType::kTH2F, {pTAxis2D, dcaAxis}});
      registryDimuon.add("dimuon/same-event/opposite-sign/DCA/pT_MuPosTDCAx_minus_MuNegTDCAx_MuonKine_MuonCuts", "DCA_x #mu^{+} top minus DCA_x #mu^{-} top and #mu^{+}#mu^{-} p_{T}", {HistType::kTH2F, {pTAxis2D, dcaAxis}});
      registryDimuon.add("dimuon/same-event/opposite-sign/DCA/pT_MuPosTDCAx_minus_MuNegBDCAx_MuonKine_MuonCuts", "DCA_x #mu^{+} top minus DCA_x #mu^{-} bottom and #mu^{+}#mu^{-} p_{T}", {HistType::kTH2F, {pTAxis2D, dcaAxis}});
      registryDimuon.add("dimuon/same-event/opposite-sign/DCA/pT_MuPosBDCAx_minus_MuNegBDCAx_MuonKine_MuonCuts", "DCA_x #mu^{+} bottom minus DCA_x #mu^{-} bottom and #mu^{+}#mu^{-} p_{T}", {HistType::kTH2F, {pTAxis2D, dcaAxis}});
      registryDimuon.add("dimuon/same-event/opposite-sign/DCA/pT_MuPosBDCAx_minus_MuNegTDCAx_MuonKine_MuonCuts", "DCA_x #mu^{+} bottom minus DCA_x #mu^{-} top and #mu^{+}#mu^{-} p_{T}", {HistType::kTH2F, {pTAxis2D, dcaAxis}});
      registryDimuon.add("dimuon/same-event/opposite-sign/DCA/pT_MuPosDCAy_minus_MuNegDCAy_MuonKine_MuonCuts", "DCA_y #mu^{+} minus DCA_y #mu^{-} and #mu^{+}#mu^{-} p_{T}", {HistType::kTH2F, {pTAxis2D, dcaAxis}});
      registryDimuon.add("dimuon/same-event/opposite-sign/DCA/pT_MuPosTDCAy_minus_MuNegTDCAy_MuonKine_MuonCuts", "DCA_y #mu^{+} top minus DCA_y #mu^{-} top and #mu^{+}#mu^{-} p_{T}", {HistType::kTH2F, {pTAxis2D, dcaAxis}});
      registryDimuon.add("dimuon/same-event/opposite-sign/DCA/pT_MuPosTDCAy_minus_MuNegBDCAy_MuonKine_MuonCuts", "DCA_y #mu^{+} top minus DCA_y #mu^{-} bottom and #mu^{+}#mu^{-} p_{T}", {HistType::kTH2F, {pTAxis2D, dcaAxis}});
      registryDimuon.add("dimuon/same-event/opposite-sign/DCA/pT_MuPosBDCAy_minus_MuNegBDCAy_MuonKine_MuonCuts", "DCA_y #mu^{+} bottom minus DCA_y #mu^{-} bottom and #mu^{+}#mu^{-} p_{T}", {HistType::kTH2F, {pTAxis2D, dcaAxis}});
      registryDimuon.add("dimuon/same-event/opposite-sign/DCA/pT_MuPosBDCAy_minus_MuNegTDCAy_MuonKine_MuonCuts", "DCA_y #mu^{+} bottom minus DCA_y #mu^{-} top and #mu^{+}#mu^{-} p_{T}", {HistType::kTH2F, {pTAxis2D, dcaAxis}});
      if (configQAs.fEnableQADimuonSameSignDCA) {
        // mu+mu+
        registryDimuon.add("dimuon/same-event/same-sign-PP/DCA/pT_Mu1DCAx_minus_Mu2DCAx_MuonKine_MuonCuts", "DCA_x #mu_{1} minus DCA_x #mu_{2} and #mu^{+}#mu^{+} p_{T}", {HistType::kTH2F, {pTAxis2D, dcaAxis}});
        registryDimuon.add("dimuon/same-event/same-sign-PP/DCA/pT_Mu1TDCAx_minus_Mu2TDCAx_MuonKine_MuonCuts", "DCA_x #mu_{1} top minus DCA_x #mu_{2} top and #mu^{+}#mu^{+} p_{T}", {HistType::kTH2F, {pTAxis2D, dcaAxis}});
        registryDimuon.add("dimuon/same-event/same-sign-PP/DCA/pT_MuTDCAx_minus_MuBDCAx_MuonKine_MuonCuts", "DCA_x #mu top minus DCA_x #mu bottom and #mu^{+}#mu^{+} p_{T}", {HistType::kTH2F, {pTAxis2D, dcaAxis}});
        registryDimuon.add("dimuon/same-event/same-sign-PP/DCA/pT_Mu1BDCAx_minus_Mu2BDCAx_MuonKine_MuonCuts", "DCA_x #mu_{1} bottom minus DCA_x #mu_{2} bottom and #mu^{+}#mu^{+} p_{T}", {HistType::kTH2F, {pTAxis2D, dcaAxis}});
        registryDimuon.add("dimuon/same-event/same-sign-PP/DCA/pT_Mu1DCAy_minus_Mu2DCAy_MuonKine_MuonCuts", "DCA_y #mu_{1} minus DCA_y #mu_{2} and #mu^{+}#mu^{+} p_{T}", {HistType::kTH2F, {pTAxis2D, dcaAxis}});
        registryDimuon.add("dimuon/same-event/same-sign-PP/DCA/pT_Mu1TDCAy_minus_Mu2TDCAy_MuonKine_MuonCuts", "DCA_y #mu_{1} top minus DCA_y #mu_{2} top and #mu^{+}#mu^{+} p_{T}", {HistType::kTH2F, {pTAxis2D, dcaAxis}});
        registryDimuon.add("dimuon/same-event/same-sign-PP/DCA/pT_MuTDCAy_minus_MuBDCAy_MuonKine_MuonCuts", "DCA_y #mu top minus DCA_y #mu bottom and #mu^{+}#mu^{+} p_{T}", {HistType::kTH2F, {pTAxis2D, dcaAxis}});
        registryDimuon.add("dimuon/same-event/same-sign-PP/DCA/pT_Mu1BDCAy_minus_Mu2BDCAy_MuonKine_MuonCuts", "DCA_y #mu_{1} bottom minus DCA_y #mu_{2} bottom and #mu^{+}#mu^{+} p_{T}", {HistType::kTH2F, {pTAxis2D, dcaAxis}});
        // mu-mu-
        registryDimuon.add("dimuon/same-event/same-sign-MM/DCA/pT_Mu1DCAx_minus_Mu2DCAx_MuonKine_MuonCuts", "DCA_x #mu_{1} minus DCA_x #mu_{2} and #mu^{-}#mu^{-} p_{T}", {HistType::kTH2F, {pTAxis2D, dcaAxis}});
        registryDimuon.add("dimuon/same-event/same-sign-MM/DCA/pT_Mu1TDCAx_minus_Mu2TDCAx_MuonKine_MuonCuts", "DCA_x #mu_{1} top minus DCA_x #mu_{2} top and #mu^{-}#mu^{-} p_{T}", {HistType::kTH2F, {pTAxis2D, dcaAxis}});
        registryDimuon.add("dimuon/same-event/same-sign-MM/DCA/pT_MuTDCAx_minus_MuBDCAx_MuonKine_MuonCuts", "DCA_x #mu top minus DCA_x #mu bottom and #mu^{-}#mu^{-} p_{T}", {HistType::kTH2F, {pTAxis2D, dcaAxis}});
        registryDimuon.add("dimuon/same-event/same-sign-MM/DCA/pT_Mu1BDCAx_minus_Mu2BDCAx_MuonKine_MuonCuts", "DCA_x #mu_{1} bottom minus DCA_x #mu_{2} bottom and #mu^{-}#mu^{-} p_{T}", {HistType::kTH2F, {pTAxis2D, dcaAxis}});
        registryDimuon.add("dimuon/same-event/same-sign-MM/DCA/pT_Mu1DCAy_minus_Mu2DCAy_MuonKine_MuonCuts", "DCA_y #mu_{1} minus DCA_y #mu_{-} and #mu^{-}#mu^{-} p_{T}", {HistType::kTH2F, {pTAxis2D, dcaAxis}});
        registryDimuon.add("dimuon/same-event/same-sign-MM/DCA/pT_Mu1TDCAy_minus_Mu2TDCAy_MuonKine_MuonCuts", "DCA_y #mu_{1} top minus DCA_y #mu_{2} top and #mu^{-}#mu^{-} p_{T}", {HistType::kTH2F, {pTAxis2D, dcaAxis}});
        registryDimuon.add("dimuon/same-event/same-sign-MM/DCA/pT_MuTDCAy_minus_MuBDCAy_MuonKine_MuonCuts", "DCA_y #mu top minus DCA_y #mu bottom and #mu^{-}#mu^{-} p_{T}", {HistType::kTH2F, {pTAxis2D, dcaAxis}});
        registryDimuon.add("dimuon/same-event/same-sign-MM/DCA/pT_Mu1BDCAy_minus_Mu2BDCAy_MuonKine_MuonCuts", "DCA_y #mu_{1} bottom minus DCA_y #mu_{2} bottom and #mu^{-}#mu^{-} p_{T}", {HistType::kTH2F, {pTAxis2D, dcaAxis}});
      }
      registryDimuon.add("dimuon/mixed-event/DCA/pT_MuPosDCAx_minus_MuNegDCAx_MuonKine_MuonCuts", "DCA_x #mu^{+} minus DCA_x #mu^{-} and #mu^{+}#mu^{-} p_{T}", {HistType::kTH2F, {pTAxis2D, dcaAxis}});
      registryDimuon.add("dimuon/mixed-event/DCA/pT_MuPosTDCAx_minus_MuNegTDCAx_MuonKine_MuonCuts", "DCA_x #mu^{+} top minus DCA_x #mu^{-} top and #mu^{+}#mu^{-} p_{T}", {HistType::kTH2F, {pTAxis2D, dcaAxis}});
      registryDimuon.add("dimuon/mixed-event/DCA/pT_MuPosTDCAx_minus_MuNegBDCAx_MuonKine_MuonCuts", "DCA_x #mu^{+} top minus DCA_x #mu^{-} bottom and #mu^{+}#mu^{-} p_{T}", {HistType::kTH2F, {pTAxis2D, dcaAxis}});
      registryDimuon.add("dimuon/mixed-event/DCA/pT_MuPosBDCAx_minus_MuNegBDCAx_MuonKine_MuonCuts", "DCA_x #mu^{+} bottom minus DCA_x #mu^{-} bottom and #mu^{+}#mu^{-} p_{T}", {HistType::kTH2F, {pTAxis2D, dcaAxis}});
      registryDimuon.add("dimuon/mixed-event/DCA/pT_MuPosBDCAx_minus_MuNegTDCAx_MuonKine_MuonCuts", "DCA_x #mu^{+} bottom minus DCA_x #mu^{-} top and #mu^{+}#mu^{-} p_{T}", {HistType::kTH2F, {pTAxis2D, dcaAxis}});
      registryDimuon.add("dimuon/mixed-event/DCA/pT_MuPosDCAy_minus_MuNegDCAy_MuonKine_MuonCuts", "DCA_y #mu^{+} minus DCA_y #mu^{-} and #mu^{+}#mu^{-} p_{T}", {HistType::kTH2F, {pTAxis2D, dcaAxis}});
      registryDimuon.add("dimuon/mixed-event/DCA/pT_MuPosTDCAy_minus_MuNegTDCAy_MuonKine_MuonCuts", "DCA_y #mu^{+} top minus DCA_y #mu^{-} top and #mu^{+}#mu^{-} p_{T}", {HistType::kTH2F, {pTAxis2D, dcaAxis}});
      registryDimuon.add("dimuon/mixed-event/DCA/pT_MuPosTDCAy_minus_MuNegBDCAy_MuonKine_MuonCuts", "DCA_y #mu^{+} top minus DCA_y #mu^{-} bottom and #mu^{+}#mu^{-} p_{T}", {HistType::kTH2F, {pTAxis2D, dcaAxis}});
      registryDimuon.add("dimuon/mixed-event/DCA/pT_MuPosBDCAy_minus_MuNegBDCAy_MuonKine_MuonCuts", "DCA_y #mu^{+} bottom minus DCA_y #mu^{-} bottom and #mu^{+}#mu^{-} p_{T}", {HistType::kTH2F, {pTAxis2D, dcaAxis}});
      registryDimuon.add("dimuon/mixed-event/DCA/pT_MuPosBDCAy_minus_MuNegTDCAy_MuonKine_MuonCuts", "DCA_y #mu^{+} bottom minus DCA_y #mu^{-} top and #mu^{+}#mu^{-} p_{T}", {HistType::kTH2F, {pTAxis2D, dcaAxis}});
      if (configQAs.fEnableSingleMuonDiMuonCorrelations) {
        // Single muons - dimuons correlations
        registryDimuon.add("dimuon/same-event/single-muon-dimuon-correlations/invariantMass_pT_MuPosPt_MuonKine_MuonCuts", "#mu^{+}#mu^{-} and #mu^{+} p_{T}", {HistType::kTH3F, {invMassAxis2D, pTAxis2D, pTAxis2D}});
        registryDimuon.add("dimuon/same-event/single-muon-dimuon-correlations/invariantMass_pT_MuNegPt_MuonKine_MuonCuts", "#mu^{+}#mu^{-} and #mu^{-} p_{T}", {HistType::kTH3F, {invMassAxis2D, pTAxis2D, pTAxis2D}});
        //
        registryDimuon.add("dimuon/same-event/single-muon-dimuon-correlations/invariantMass_pT_MuPosEta_MuonKine_MuonCuts", "#mu^{+}#mu^{-} and #mu^{+} #eta", {HistType::kTH3F, {invMassAxis2D, pTAxis2D, etaAxis}});
        registryDimuon.add("dimuon/same-event/single-muon-dimuon-correlations/invariantMass_pT_MuNegEta_MuonKine_MuonCuts", "#mu^{+}#mu^{-} and #mu^{-} #eta", {HistType::kTH3F, {invMassAxis2D, pTAxis2D, etaAxis}});
        //
        registryDimuon.add("dimuon/same-event/single-muon-dimuon-correlations/invariantMass_pT_MuPosRabs_MuonKine_MuonCuts", "#mu^{+}#mu^{-} and #mu^{+} R_{abs}", {HistType::kTH3F, {invMassAxis2D, pTAxis2D, rAbsAxis}});
        registryDimuon.add("dimuon/same-event/single-muon-dimuon-correlations/invariantMass_pT_MuNegRabs_MuonKine_MuonCuts", "#mu^{+}#mu^{-} and #mu^{-} R_{abs}", {HistType::kTH3F, {invMassAxis2D, pTAxis2D, rAbsAxis}});
        //
        registryDimuon.add("dimuon/same-event/single-muon-dimuon-correlations/invariantMass_pT_MuPosDca_MuonKine_MuonCuts", "#mu^{+}#mu^{-} and #mu^{+} DCA", {HistType::kTH3F, {invMassAxis2D, pTAxis2D, dcaAxisReduced}});
        registryDimuon.add("dimuon/same-event/single-muon-dimuon-correlations/invariantMass_pT_MuNegDca_MuonKine_MuonCuts", "#mu^{+}#mu^{-} and #mu^{-} DCA", {HistType::kTH3F, {invMassAxis2D, pTAxis2D, dcaAxisReduced}});
        //
        registryDimuon.add("dimuon/same-event/single-muon-dimuon-correlations/invariantMass_pT_MuPosPhi_MuonKine_MuonCuts", "#mu^{+}#mu^{-} and #mu^{+} #phi", {HistType::kTH3F, {invMassAxis2D, pTAxis2D, phiAxis}});
        registryDimuon.add("dimuon/same-event/single-muon-dimuon-correlations/invariantMass_pT_MuNegPhi_MuonKine_MuonCuts", "#mu^{+}#mu^{-} and #mu^{-} #phi", {HistType::kTH3F, {invMassAxis2D, pTAxis2D, phiAxis}});
      }
      // MCH-MID tracks with MCH acceptance cuts
      registryDimuon.add("dimuon/same-event/invariantMass_MuonKine_MuonCuts", "#mu^{+}#mu^{-} invariant mass", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/same-event/invariantMassFull_MuonKine_MuonCuts", "#mu^{+}#mu^{-} invariant mass", {HistType::kTH1F, {invMassAxisFull}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_MuonKine_MuonCuts", "#mu^{+}#mu^{-} invariant mass", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/mixed-event/invariantMassFull_MuonKine_MuonCuts", "#mu^{+}#mu^{-} invariant mass", {HistType::kTH1F, {invMassAxisFull}});
      // -- Mass and pT
      registryDimuon.add("dimuon/same-event/invariantMass_pT_MuonKine_MuonCuts", "#mu^{+}#mu^{-} invariant mass", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_pT_MuonKine_MuonCuts", "#mu^{+}#mu^{-} invariant mass", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      // MCH-MID tracks with MCH acceptance cuts and combinations from the top and bottom halfs of MCH
      registryDimuon.add("dimuon/same-event/invariantMass_MuonKine_MuonCuts_TT", "#mu^{+}#mu^{-} invariant mass, top-top", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/same-event/invariantMassFull_MuonKine_MuonCuts_TT", "#mu^{+}#mu^{-} invariant mass, top-top", {HistType::kTH1F, {invMassAxisFull}});
      registryDimuon.add("dimuon/same-event/invariantMass_MuonKine_MuonCuts_TB", "#mu^{+}#mu^{-} invariant mass, top-bottom or bottom-top", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/same-event/invariantMassFull_MuonKine_MuonCuts_TB", "#mu^{+}#mu^{-} invariant mass, top-bottom or bottom-top", {HistType::kTH1F, {invMassAxisFull}});
      registryDimuon.add("dimuon/same-event/invariantMass_MuonKine_MuonCuts_TPBN", "#mu^{+}#mu^{-} invariant mass, #mu^{+} top and #mu^{-} bottom", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/same-event/invariantMass_MuonKine_MuonCuts_TNBP", "#mu^{+}#mu^{-} invariant mass, #mu^{-} top and #mu^{+} bottom", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/same-event/invariantMassFull_MuonKine_MuonCuts_TPBN", "#mu^{+}#mu^{-} invariant mass, #mu^{+} top and #mu^{-} bottom", {HistType::kTH1F, {invMassAxisFull}});
      registryDimuon.add("dimuon/same-event/invariantMassFull_MuonKine_MuonCuts_TNBP", "#mu^{+}#mu^{-} invariant mass, #mu^{-} top and #mu^{+} bottom", {HistType::kTH1F, {invMassAxisFull}});
      registryDimuon.add("dimuon/same-event/invariantMass_MuonKine_MuonCuts_BB", "#mu^{+}#mu^{-} invariant mass, bottom-bottom", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/same-event/invariantMassFull_MuonKine_MuonCuts_BB", "#mu^{+}#mu^{-} invariant mass, bottom-bottom", {HistType::kTH1F, {invMassAxisFull}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_MuonKine_MuonCuts_TT", "#mu^{+}#mu^{-} invariant mass, top-top", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/mixed-event/invariantMassFull_MuonKine_MuonCuts_TT", "#mu^{+}#mu^{-} invariant mass, top-top", {HistType::kTH1F, {invMassAxisFull}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_MuonKine_MuonCuts_TB", "#mu^{+}#mu^{-} invariant mass, top-bottom or bottom-top", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/mixed-event/invariantMassFull_MuonKine_MuonCuts_TB", "#mu^{+}#mu^{-} invariant mass, top-bottom or bottom-top", {HistType::kTH1F, {invMassAxisFull}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_MuonKine_MuonCuts_TPBN", "#mu^{+}#mu^{-} invariant mass, #mu^{+} top and #mu^{-} bottom", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_MuonKine_MuonCuts_TNBP", "#mu^{+}#mu^{-} invariant mass, #mu^{-} top and #mu^{+} bottom", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/mixed-event/invariantMassFull_MuonKine_MuonCuts_TPBN", "#mu^{+}#mu^{-} invariant mass, #mu^{+} top and #mu^{-} bottom", {HistType::kTH1F, {invMassAxisFull}});
      registryDimuon.add("dimuon/mixed-event/invariantMassFull_MuonKine_MuonCuts_TNBP", "#mu^{+}#mu^{-} invariant mass, #mu^{-} top and #mu^{+} bottom", {HistType::kTH1F, {invMassAxisFull}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_MuonKine_MuonCuts_BB", "#mu^{+}#mu^{-} invariant mass, bottom-bottom", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/mixed-event/invariantMassFull_MuonKine_MuonCuts_BB", "#mu^{+}#mu^{-} invariant mass, bottom-bottom", {HistType::kTH1F, {invMassAxisFull}});
      // -- Mass and pT
      registryDimuon.add("dimuon/same-event/invariantMass_pT_MuonKine_MuonCuts_TT", "#mu^{+}#mu^{-} invariant mass, top-top", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/same-event/invariantMass_pT_MuonKine_MuonCuts_TB", "#mu^{+}#mu^{-} invariant mass, top-bottom or bottom-top", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/same-event/invariantMass_pT_MuonKine_MuonCuts_TPBN", "#mu^{+}#mu^{-} invariant mass, #mu^{+} top and #mu^{-} bottom", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/same-event/invariantMass_pT_MuonKine_MuonCuts_TNBP", "#mu^{+}#mu^{-} invariant mass, #mu^{-} top and #mu^{+} bottom", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/same-event/invariantMass_pT_MuonKine_MuonCuts_BB", "#mu^{+}#mu^{-} invariant mass, bottom-bottom", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_pT_MuonKine_MuonCuts_TT", "#mu^{+}#mu^{-} invariant mass, top-top", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_pT_MuonKine_MuonCuts_TB", "#mu^{+}#mu^{-} invariant mass, top-bottom or bottom-top", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_pT_MuonKine_MuonCuts_TPBN", "#mu^{+}#mu^{-} invariant mass, #mu^{+} top and #mu^{-} bottom", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_pT_MuonKine_MuonCuts_TNBP", "#mu^{+}#mu^{-} invariant mass, #mu^{-} top and #mu^{+} bottom", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_pT_MuonKine_MuonCuts_BB", "#mu^{+}#mu^{-} invariant mass, bottom-bottom", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      // MCH-MID tracks with MFT acceptance cuts and combinations from the top and bottom halfs of MCH
      registryDimuon.add("dimuon/same-event/invariantMass_MuonKine_GlobalMuonCuts", "#mu^{+}#mu^{-} invariant mass", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/same-event/invariantMassFull_MuonKine_GlobalMuonCuts", "#mu^{+}#mu^{-} invariant mass", {HistType::kTH1F, {invMassAxisFull}});
      registryDimuon.add("dimuon/same-event/invariantMass_MuonKine_GlobalMuonCuts_TT", "#mu^{+}#mu^{-} invariant mass, top-top", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/same-event/invariantMassFull_MuonKine_GlobalMuonCuts_TT", "#mu^{+}#mu^{-} invariant mass, top-top", {HistType::kTH1F, {invMassAxisFull}});
      registryDimuon.add("dimuon/same-event/invariantMass_MuonKine_GlobalMuonCuts_TB", "#mu^{+}#mu^{-} invariant mass, top-bottom or bottom-top", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/same-event/invariantMassFull_MuonKine_GlobalMuonCuts_TB", "#mu^{+}#mu^{-} invariant mass, top-bottom or bottom-top", {HistType::kTH1F, {invMassAxisFull}});
      registryDimuon.add("dimuon/same-event/invariantMass_MuonKine_GlobalMuonCuts_TPBN", "#mu^{+}#mu^{-} invariant mass, #mu^{+} top and #mu^{-} bottom", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/same-event/invariantMass_MuonKine_GlobalMuonCuts_TNBP", "#mu^{+}#mu^{-} invariant mass, #mu^{-} top and #mu^{+} bottom", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/same-event/invariantMassFull_MuonKine_GlobalMuonCuts_TPBN", "#mu^{+}#mu^{-} invariant mass, #mu^{+} top and #mu^{-} bottom", {HistType::kTH1F, {invMassAxisFull}});
      registryDimuon.add("dimuon/same-event/invariantMassFull_MuonKine_GlobalMuonCuts_TNBP", "#mu^{+}#mu^{-} invariant mass, #mu^{-} top and #mu^{+} bottom", {HistType::kTH1F, {invMassAxisFull}});
      registryDimuon.add("dimuon/same-event/invariantMass_MuonKine_GlobalMuonCuts_BB", "#mu^{+}#mu^{-} invariant mass, bottom-bottom", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/same-event/invariantMassFull_MuonKine_GlobalMuonCuts_BB", "#mu^{+}#mu^{-} invariant mass, bottom-bottom", {HistType::kTH1F, {invMassAxisFull}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_MuonKine_GlobalMuonCuts", "#mu^{+}#mu^{-} invariant mass", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/mixed-event/invariantMassFull_MuonKine_GlobalMuonCuts", "#mu^{+}#mu^{-} invariant mass", {HistType::kTH1F, {invMassAxisFull}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_MuonKine_GlobalMuonCuts_TT", "#mu^{+}#mu^{-} invariant mass, top-top", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/mixed-event/invariantMassFull_MuonKine_GlobalMuonCuts_TT", "#mu^{+}#mu^{-} invariant mass, top-top", {HistType::kTH1F, {invMassAxisFull}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_MuonKine_GlobalMuonCuts_TB", "#mu^{+}#mu^{-} invariant mass, top-bottom or bottom-top", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/mixed-event/invariantMassFull_MuonKine_GlobalMuonCuts_TB", "#mu^{+}#mu^{-} invariant mass, top-bottom or bottom-top", {HistType::kTH1F, {invMassAxisFull}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_MuonKine_GlobalMuonCuts_TPBN", "#mu^{+}#mu^{-} invariant mass, #mu^{+} top and #mu^{-} bottom", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_MuonKine_GlobalMuonCuts_TNBP", "#mu^{+}#mu^{-} invariant mass, #mu^{-} top and #mu^{+} bottom", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/mixed-event/invariantMassFull_MuonKine_GlobalMuonCuts_TPBN", "#mu^{+}#mu^{-} invariant mass, #mu^{+} top and #mu^{-} bottom", {HistType::kTH1F, {invMassAxisFull}});
      registryDimuon.add("dimuon/mixed-event/invariantMassFull_MuonKine_GlobalMuonCuts_TNBP", "#mu^{+}#mu^{-} invariant mass, #mu^{-} top and #mu^{+} bottom", {HistType::kTH1F, {invMassAxisFull}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_MuonKine_GlobalMuonCuts_BB", "#mu^{+}#mu^{-} invariant mass, bottom-bottom", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/mixed-event/invariantMassFull_MuonKine_GlobalMuonCuts_BB", "#mu^{+}#mu^{-} invariant mass, bottom-bottom", {HistType::kTH1F, {invMassAxisFull}});
      // -- Mass and pT
      registryDimuon.add("dimuon/same-event/invariantMass_pT_MuonKine_GlobalMuonCuts", "#mu^{+}#mu^{-} invariant mass", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/same-event/invariantMass_pT_MuonKine_GlobalMuonCuts_TT", "#mu^{+}#mu^{-} invariant mass, top-top", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/same-event/invariantMass_pT_MuonKine_GlobalMuonCuts_TB", "#mu^{+}#mu^{-} invariant mass, top-bottom or bottom-top", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/same-event/invariantMass_pT_MuonKine_GlobalMuonCuts_TPBN", "#mu^{+}#mu^{-} invariant mass, #mu^{+} top and #mu^{-} bottom", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/same-event/invariantMass_pT_MuonKine_GlobalMuonCuts_TNBP", "#mu^{+}#mu^{-} invariant mass, #mu^{-} top and #mu^{+} bottom", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/same-event/invariantMass_pT_MuonKine_GlobalMuonCuts_BB", "#mu^{+}#mu^{-} invariant mass, bottom-bottom", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_pT_MuonKine_GlobalMuonCuts", "#mu^{+}#mu^{-} invariant mass", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_pT_MuonKine_GlobalMuonCuts_TT", "#mu^{+}#mu^{-} invariant mass, top-top", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_pT_MuonKine_GlobalMuonCuts_TB", "#mu^{+}#mu^{-} invariant mass, top-bottom or bottom-top", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_pT_MuonKine_GlobalMuonCuts_TPBN", "#mu^{+}#mu^{-} invariant mass, #mu^{+} top and #mu^{-} bottom", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_pT_MuonKine_GlobalMuonCuts_TNBP", "#mu^{+}#mu^{-} invariant mass, #mu^{-} top and #mu^{+} bottom", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_pT_MuonKine_GlobalMuonCuts_BB", "#mu^{+}#mu^{-} invariant mass, bottom-bottom", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      // MCH-MID tracks with MCH acceptance cuts and combinations from the left and right halfs of MCH
      registryDimuon.add("dimuon/same-event/invariantMass_MuonKine_MuonCuts_LL", "#mu^{+}#mu^{-} invariant mass, left-left", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/same-event/invariantMassFull_MuonKine_MuonCuts_LL", "#mu^{+}#mu^{-} invariant mass, left-left", {HistType::kTH1F, {invMassAxisFull}});
      registryDimuon.add("dimuon/same-event/invariantMass_MuonKine_MuonCuts_LR", "#mu^{+}#mu^{-} invariant mass, left-right or right-left", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/same-event/invariantMassFull_MuonKine_MuonCuts_LR", "#mu^{+}#mu^{-} invariant mass, left-right or right-left", {HistType::kTH1F, {invMassAxisFull}});
      registryDimuon.add("dimuon/same-event/invariantMass_MuonKine_MuonCuts_LPRN", "#mu^{+}#mu^{-} invariant mass, #mu^{+} left and #mu^{-} right", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/same-event/invariantMass_MuonKine_MuonCuts_LNRP", "#mu^{+}#mu^{-} invariant mass, #mu^{-} left and #mu^{+} right", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/same-event/invariantMassFull_MuonKine_MuonCuts_LPRN", "#mu^{+}#mu^{-} invariant mass, #mu^{+} left and #mu^{-} right", {HistType::kTH1F, {invMassAxisFull}});
      registryDimuon.add("dimuon/same-event/invariantMassFull_MuonKine_MuonCuts_LNRP", "#mu^{+}#mu^{-} invariant mass, #mu^{-} left and #mu^{+} right", {HistType::kTH1F, {invMassAxisFull}});
      registryDimuon.add("dimuon/same-event/invariantMass_MuonKine_MuonCuts_RR", "#mu^{+}#mu^{-} invariant mass, right-right", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/same-event/invariantMassFull_MuonKine_MuonCuts_RR", "#mu^{+}#mu^{-} invariant mass, right-right", {HistType::kTH1F, {invMassAxisFull}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_MuonKine_MuonCuts_LL", "#mu^{+}#mu^{-} invariant mass, left-left", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/mixed-event/invariantMassFull_MuonKine_MuonCuts_LL", "#mu^{+}#mu^{-} invariant mass, left-left", {HistType::kTH1F, {invMassAxisFull}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_MuonKine_MuonCuts_LR", "#mu^{+}#mu^{-} invariant mass, left-right or right-left", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/mixed-event/invariantMassFull_MuonKine_MuonCuts_LR", "#mu^{+}#mu^{-} invariant mass, left-right or right-left", {HistType::kTH1F, {invMassAxisFull}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_MuonKine_MuonCuts_LPRN", "#mu^{+}#mu^{-} invariant mass, #mu^{+} left and #mu^{-} right", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_MuonKine_MuonCuts_LNRP", "#mu^{+}#mu^{-} invariant mass, #mu^{-} left and #mu^{+} right", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/mixed-event/invariantMassFull_MuonKine_MuonCuts_LPRN", "#mu^{+}#mu^{-} invariant mass, #mu^{+} left and #mu^{-} right", {HistType::kTH1F, {invMassAxisFull}});
      registryDimuon.add("dimuon/mixed-event/invariantMassFull_MuonKine_MuonCuts_LNRP", "#mu^{+}#mu^{-} invariant mass, #mu^{-} left and #mu^{+} right", {HistType::kTH1F, {invMassAxisFull}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_MuonKine_MuonCuts_RR", "#mu^{+}#mu^{-} invariant mass, right-right", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/mixed-event/invariantMassFull_MuonKine_MuonCuts_RR", "#mu^{+}#mu^{-} invariant mass, right-right", {HistType::kTH1F, {invMassAxisFull}});
      // -- Mass and pT
      registryDimuon.add("dimuon/same-event/invariantMass_pT_MuonKine_MuonCuts_LL", "#mu^{+}#mu^{-} invariant mass, left-left", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/same-event/invariantMass_pT_MuonKine_MuonCuts_LR", "#mu^{+}#mu^{-} invariant mass, left-right or right-left", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/same-event/invariantMass_pT_MuonKine_MuonCuts_LPRN", "#mu^{+}#mu^{-} invariant mass, #mu^{+} left and #mu^{-} right", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/same-event/invariantMass_pT_MuonKine_MuonCuts_LNRP", "#mu^{+}#mu^{-} invariant mass, #mu^{-} left and #mu^{+} right", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/same-event/invariantMass_pT_MuonKine_MuonCuts_RR", "#mu^{+}#mu^{-} invariant mass, right-right", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_pT_MuonKine_MuonCuts_LL", "#mu^{+}#mu^{-} invariant mass, left-left", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_pT_MuonKine_MuonCuts_LR", "#mu^{+}#mu^{-} invariant mass, left-right or left-right", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_pT_MuonKine_MuonCuts_LPRN", "#mu^{+}#mu^{-} invariant mass, #mu^{+} left and #mu^{-} right", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_pT_MuonKine_MuonCuts_LNRP", "#mu^{+}#mu^{-} invariant mass, #mu^{-} left and #mu^{+} right", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_pT_MuonKine_MuonCuts_RR", "#mu^{+}#mu^{-} invariant mass, right-right", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      // MCH-MID tracks with MFT acceptance cuts and combinations from the left and right halfs of MCH
      registryDimuon.add("dimuon/same-event/invariantMass_MuonKine_GlobalMuonCuts_LL", "#mu^{+}#mu^{-} invariant mass, left-left", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/same-event/invariantMassFull_MuonKine_GlobalMuonCuts_LL", "#mu^{+}#mu^{-} invariant mass, left-left", {HistType::kTH1F, {invMassAxisFull}});
      registryDimuon.add("dimuon/same-event/invariantMass_MuonKine_GlobalMuonCuts_LR", "#mu^{+}#mu^{-} invariant mass, left-right or right-left", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/same-event/invariantMassFull_MuonKine_GlobalMuonCuts_LR", "#mu^{+}#mu^{-} invariant mass, left-right or right-left", {HistType::kTH1F, {invMassAxisFull}});
      registryDimuon.add("dimuon/same-event/invariantMass_MuonKine_GlobalMuonCuts_LPRN", "#mu^{+}#mu^{-} invariant mass, #mu^{+} left and #mu^{-} right", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/same-event/invariantMass_MuonKine_GlobalMuonCuts_LNRP", "#mu^{+}#mu^{-} invariant mass, #mu^{-} left and #mu^{+} right", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/same-event/invariantMassFull_MuonKine_GlobalMuonCuts_LPRN", "#mu^{+}#mu^{-} invariant mass, #mu^{+} left and #mu^{-} right", {HistType::kTH1F, {invMassAxisFull}});
      registryDimuon.add("dimuon/same-event/invariantMassFull_MuonKine_GlobalMuonCuts_LNRP", "#mu^{+}#mu^{-} invariant mass, #mu^{-} left and #mu^{+} right", {HistType::kTH1F, {invMassAxisFull}});
      registryDimuon.add("dimuon/same-event/invariantMass_MuonKine_GlobalMuonCuts_RR", "#mu^{+}#mu^{-} invariant mass, right-right", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/same-event/invariantMassFull_MuonKine_GlobalMuonCuts_RR", "#mu^{+}#mu^{-} invariant mass, right-right", {HistType::kTH1F, {invMassAxisFull}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_MuonKine_GlobalMuonCuts_LL", "#mu^{+}#mu^{-} invariant mass, left-left", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/mixed-event/invariantMassFull_MuonKine_GlobalMuonCuts_LL", "#mu^{+}#mu^{-} invariant mass, left-left", {HistType::kTH1F, {invMassAxisFull}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_MuonKine_GlobalMuonCuts_LR", "#mu^{+}#mu^{-} invariant mass, left-right or right-left", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/mixed-event/invariantMassFull_MuonKine_GlobalMuonCuts_LR", "#mu^{+}#mu^{-} invariant mass, left-right or right-left", {HistType::kTH1F, {invMassAxisFull}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_MuonKine_GlobalMuonCuts_LPRN", "#mu^{+}#mu^{-} invariant mass, #mu^{+} left and #mu^{-} right", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_MuonKine_GlobalMuonCuts_LNRP", "#mu^{+}#mu^{-} invariant mass, #mu^{-} left and #mu^{+} right", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/mixed-event/invariantMassFull_MuonKine_GlobalMuonCuts_LPRN", "#mu^{+}#mu^{-} invariant mass, #mu^{+} left and #mu^{-} right", {HistType::kTH1F, {invMassAxisFull}});
      registryDimuon.add("dimuon/mixed-event/invariantMassFull_MuonKine_GlobalMuonCuts_LNRP", "#mu^{+}#mu^{-} invariant mass, #mu^{-} left and #mu^{+} right", {HistType::kTH1F, {invMassAxisFull}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_MuonKine_GlobalMuonCuts_RR", "#mu^{+}#mu^{-} invariant mass, right-right", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/mixed-event/invariantMassFull_MuonKine_GlobalMuonCuts_RR", "#mu^{+}#mu^{-} invariant mass right-right", {HistType::kTH1F, {invMassAxisFull}});
      // -- Mass and pT
      registryDimuon.add("dimuon/same-event/invariantMass_pT_MuonKine_GlobalMuonCuts_LL", "#mu^{+}#mu^{-} invariant mass, left-left", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/same-event/invariantMass_pT_MuonKine_GlobalMuonCuts_LR", "#mu^{+}#mu^{-} invariant mass, left-right or right-left", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/same-event/invariantMass_pT_MuonKine_GlobalMuonCuts_LPRN", "#mu^{+}#mu^{-} invariant mass, #mu^{+} left and #mu^{-} right", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/same-event/invariantMass_pT_MuonKine_GlobalMuonCuts_LNRP", "#mu^{+}#mu^{-} invariant mass, #mu^{-} left and #mu^{+} right", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/same-event/invariantMass_pT_MuonKine_GlobalMuonCuts_RR", "#mu^{+}#mu^{-} invariant mass, right-right", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_pT_MuonKine_GlobalMuonCuts_LL", "#mu^{+}#mu^{-} invariant mass, left-left", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_pT_MuonKine_GlobalMuonCuts_LR", "#mu^{+}#mu^{-} invariant mass, left-right or right-left", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_pT_MuonKine_GlobalMuonCuts_LPRN", "#mu^{+}#mu^{-} invariant mass, #mu^{+} left and #mu^{-} right", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_pT_MuonKine_GlobalMuonCuts_LNRP", "#mu^{+}#mu^{-} invariant mass, #mu^{-} left and #mu^{+} right", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_pT_MuonKine_GlobalMuonCuts_RR", "#mu^{+}#mu^{-} invariant mass, right-right", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      // Good MFT-MCH-MID tracks with MCH parameters and MFT acceptance cuts
      registryDimuon.add("dimuon/same-event/invariantMass_MuonKine_GlobalMatchesCuts", "#mu^{+}#mu^{-} invariant mass", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/same-event/invariantMassFull_MuonKine_GlobalMatchesCuts", "#mu^{+}#mu^{-} invariant mass", {HistType::kTH1F, {invMassAxisFull}});
      registryDimuon.add("dimuon/same-event/invariantMass_MuonKine_GlobalMatchesCuts_TT", "#mu^{+}#mu^{-} invariant mass, top - top", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/same-event/invariantMass_MuonKine_GlobalMatchesCuts_TB", "#mu^{+}#mu^{-} invariant mass, top - bottom or bottom - top", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/same-event/invariantMass_MuonKine_GlobalMatchesCuts_TPBN", "#mu^{+}#mu^{-} invariant mass, #mu^{+} top and #mu^{-} bottom", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/same-event/invariantMass_MuonKine_GlobalMatchesCuts_TNBP", "#mu^{+}#mu^{-} invariant mass, #mu^{-} top and #mu^{+} bottom", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/same-event/invariantMass_MuonKine_GlobalMatchesCuts_BB", "#mu^{+}#mu^{-} invariant mass, bottom - bottom", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_MuonKine_GlobalMatchesCuts", "#mu^{+}#mu^{-} invariant mass", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/mixed-event/invariantMassFull_MuonKine_GlobalMatchesCuts", "#mu^{+}#mu^{-} invariant mass", {HistType::kTH1F, {invMassAxisFull}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_MuonKine_GlobalMatchesCuts_TT", "#mu^{+}#mu^{-} invariant mass, top - top", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_MuonKine_GlobalMatchesCuts_TB", "#mu^{+}#mu^{-} invariant mass, top - bottom or bottom - top", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_MuonKine_GlobalMatchesCuts_TPBN", "#mu^{+}#mu^{-} invariant mass, #mu^{+} top and #mu^{-} bottom", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_MuonKine_GlobalMatchesCuts_TNBP", "#mu^{+}#mu^{-} invariant mass, #mu^{-} top and #mu^{+} bottom", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_MuonKine_GlobalMatchesCuts_BB", "#mu^{+}#mu^{-} invariant mass, bottom - bottom", {HistType::kTH1F, {invMassAxis}});
      // -- Mass and pT
      registryDimuon.add("dimuon/same-event/invariantMass_pT_MuonKine_GlobalMatchesCuts", "#mu^{+}#mu^{-} invariant mass", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/same-event/invariantMass_pT_MuonKine_GlobalMatchesCuts_TT", "#mu^{+}#mu^{-} invariant mass, top - top", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/same-event/invariantMass_pT_MuonKine_GlobalMatchesCuts_TB", "#mu^{+}#mu^{-} invariant mass, top - bottom or bottom - top", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/same-event/invariantMass_pT_MuonKine_GlobalMatchesCuts_TPBN", "#mu^{+}#mu^{-} invariant mass, #mu^{+} top and #mu^{-} bottom", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/same-event/invariantMass_pT_MuonKine_GlobalMatchesCuts_TNBP", "#mu^{+}#mu^{-} invariant mass, #mu^{-} top and #mu^{+} bottom", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/same-event/invariantMass_pT_MuonKine_GlobalMatchesCuts_BB", "#mu^{+}#mu^{-} invariant mass, bottom - bottom", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_pT_MuonKine_GlobalMatchesCuts", "#mu^{+}#mu^{-} invariant mass", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_pT_MuonKine_GlobalMatchesCuts_TT", "#mu^{+}#mu^{-} invariant mass, top - top", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_pT_MuonKine_GlobalMatchesCuts_TB", "#mu^{+}#mu^{-} invariant mass, top - bottom or bottom - top", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_pT_MuonKine_GlobalMatchesCuts_TPBN", "#mu^{+}#mu^{-} invariant mass, #mu^{+} top and #mu^{-} bottom", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_pT_MuonKine_GlobalMatchesCuts_TNBP", "#mu^{+}#mu^{-} invariant mass, #mu^{-} top and #mu^{+} bottom", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_pT_MuonKine_GlobalMatchesCuts_BB", "#mu^{+}#mu^{-} invariant mass, bottom - bottom", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      // Good MFT-MCH-MID tracks with global parameters MFT acceptance cuts
      registryDimuon.add("dimuon/same-event/invariantMass_GlobalMuonKine_GlobalMatchesCuts", "#mu^{+}#mu^{-} invariant mass", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/same-event/invariantMass_GlobalMuonKine_GlobalMatchesCuts_TT", "#mu^{+}#mu^{-} invariant mass, top - top", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/same-event/invariantMass_GlobalMuonKine_GlobalMatchesCuts_TB", "#mu^{+}#mu^{-} invariant mass, top - bottom or bottom-top", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/same-event/invariantMass_GlobalMuonKine_GlobalMatchesCuts_TPBN", "#mu^{+}#mu^{-} invariant mass, #mu^{+} top and #mu^{-} bottom", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/same-event/invariantMass_GlobalMuonKine_GlobalMatchesCuts_TNBP", "#mu^{+}#mu^{-} invariant mass, #mu^{-} top and #mu^{+} bottom", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/same-event/invariantMass_GlobalMuonKine_GlobalMatchesCuts_BB", "#mu^{+}#mu^{-} invariant mass, bottom - bottom", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/same-event/invariantMassFull_GlobalMuonKine_GlobalMatchesCuts", "#mu^{+}#mu^{-} invariant mass", {HistType::kTH1F, {invMassAxisFull}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_GlobalMuonKine_GlobalMatchesCuts", "#mu^{+}#mu^{-} invariant mass", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_GlobalMuonKine_GlobalMatchesCuts_TT", "#mu^{+}#mu^{-} invariant mass, top - top", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_GlobalMuonKine_GlobalMatchesCuts_TB", "#mu^{+}#mu^{-} invariant mass, top - bottom or bottom - top", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_GlobalMuonKine_GlobalMatchesCuts_TPBN", "#mu^{+}#mu^{-} invariant mass, #mu^{+} top and #mu^{-} bottom", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_GlobalMuonKine_GlobalMatchesCuts_TNBP", "#mu^{+}#mu^{-} invariant mass, #mu^{-} top and #mu^{+} bottom", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_GlobalMuonKine_GlobalMatchesCuts_BB", "#mu^{+}#mu^{-} invariant mass, bottom - bottom", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/mixed-event/invariantMassFull_GlobalMuonKine_GlobalMatchesCuts", "#mu^{+}#mu^{-} invariant mass", {HistType::kTH1F, {invMassAxisFull}});
      // -- Mass and pT
      registryDimuon.add("dimuon/same-event/invariantMass_pT_GlobalMuonKine_GlobalMatchesCuts", "#mu^{+}#mu^{-} invariant mass", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/same-event/invariantMass_pT_GlobalMuonKine_GlobalMatchesCuts_TT", "#mu^{+}#mu^{-} invariant mass, top - top", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/same-event/invariantMass_pT_GlobalMuonKine_GlobalMatchesCuts_TB", "#mu^{+}#mu^{-} invariant mass, top - bottom or bottom - top", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/same-event/invariantMass_pT_GlobalMuonKine_GlobalMatchesCuts_TPBN", "#mu^{+}#mu^{-} invariant mass, #mu^{+} top and #mu^{-} bottom", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/same-event/invariantMass_pT_GlobalMuonKine_GlobalMatchesCuts_TNBP", "#mu^{+}#mu^{-} invariant mass, #mu^{-} top and #mu^{+} bottom", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/same-event/invariantMass_pT_GlobalMuonKine_GlobalMatchesCuts_BB", "#mu^{+}#mu^{-} invariant mass, bottom - bottom", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_pT_GlobalMuonKine_GlobalMatchesCuts", "#mu^{+}#mu^{-} invariant mass", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_pT_GlobalMuonKine_GlobalMatchesCuts_TT", "#mu^{+}#mu^{-} invariant mass, top - top", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_pT_GlobalMuonKine_GlobalMatchesCuts_TB", "#mu^{+}#mu^{-} invariant mass, top - bottom or bottom - top", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_pT_GlobalMuonKine_GlobalMatchesCuts_TPBN", "#mu^{+}#mu^{-} invariant mass, #mu^{+} top and #mu^{-} bottom", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_pT_GlobalMuonKine_GlobalMatchesCuts_TNBP", "#mu^{+}#mu^{-} invariant mass, #mu^{-} top and #mu^{+} bottom", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_pT_GlobalMuonKine_GlobalMatchesCuts_BB", "#mu^{+}#mu^{-} invariant mass, bottom - bottom", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      // Good MFT-MCH-MID tracks with re-scaled MFT kinematics and MFT acceptance cuts
      registryDimuon.add("dimuon/same-event/invariantMass_ScaledMftKine_GlobalMatchesCuts", "M_{#mu^{+}#mu^{-}} - rescaled MFT momentum", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/same-event/invariantMassFull_ScaledMftKine_GlobalMatchesCuts", "M_{#mu^{+}#mu^{-}} - rescaled MFT momentum", {HistType::kTH1F, {invMassAxisFull}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_ScaledMftKine_GlobalMatchesCuts", "M_{#mu^{+}#mu^{-}} - rescaled MFT momentum", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/mixed-event/invariantMassFull_ScaledMftKine_GlobalMatchesCuts", "M_{#mu^{+}#mu^{-}} - rescaled MFT momentum", {HistType::kTH1F, {invMassAxisFull}});
      // combinations of tracks from top and bottom halfs of MFT
      registryDimuon.add("dimuon/same-event/invariantMass_ScaledMftKine_GlobalMatchesCuts_TT", "M_{#mu^{+}#mu^{-}} - rescaled MFT momentum, top - top", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/same-event/invariantMass_ScaledMftKine_GlobalMatchesCuts_TB", "M_{#mu^{+}#mu^{-}} - rescaled MFT momentum, top - bottom or bottom - top", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/same-event/invariantMass_ScaledMftKine_GlobalMatchesCuts_TPBN", "#mu^{+}#mu^{-} - rescaled MFT momentum, #mu^{+} top and #mu^{-} bottom", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/same-event/invariantMass_ScaledMftKine_GlobalMatchesCuts_TNBP", "#mu^{+}#mu^{-} - rescaled MFT momentum, #mu^{-} top and #mu^{+} bottom", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/same-event/invariantMass_ScaledMftKine_GlobalMatchesCuts_BB", "M_{#mu^{+}#mu^{-}} - rescaled MFT momentum, bottom-bottom", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_ScaledMftKine_GlobalMatchesCuts_TT", "M_{#mu^{+}#mu^{-}} - rescaled MFT momentum, top - top", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_ScaledMftKine_GlobalMatchesCuts_TB", "M_{#mu^{+}#mu^{-}} - rescaled MFT momentum, top - bottom or bottom - top", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_ScaledMftKine_GlobalMatchesCuts_TPBN", "#mu^{+}#mu^{-} - rescaled MFT momentum, #mu^{+} top and #mu^{-} bottom", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_ScaledMftKine_GlobalMatchesCuts_TNBP", "#mu^{+}#mu^{-} - rescaled MFT momentum, #mu^{-} top and #mu^{+} bottom", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_ScaledMftKine_GlobalMatchesCuts_BB", "M_{#mu^{+}#mu^{-}} - rescaled MFT momentum, bottom - bottom", {HistType::kTH1F, {invMassAxis}});
      // -- Mass and pT
      registryDimuon.add("dimuon/same-event/invariantMass_pT_ScaledMftKine_GlobalMatchesCuts", "M_{#mu^{+}#mu^{-}} - rescaled MFT momentum", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/same-event/invariantMass_pT_ScaledMftKine_GlobalMatchesCuts_TT", "M_{#mu^{+}#mu^{-}} - rescaled MFT momentum, top - top", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/same-event/invariantMass_pT_ScaledMftKine_GlobalMatchesCuts_TB", "M_{#mu^{+}#mu^{-}} - rescaled MFT momentum, top - bottom or bottom - top", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/same-event/invariantMass_pT_ScaledMftKine_GlobalMatchesCuts_TPBN", "#mu^{+}#mu^{-} - rescaled MFT momentum, #mu^{+} top and #mu^{-} bottom", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/same-event/invariantMass_pT_ScaledMftKine_GlobalMatchesCuts_TNBP", "#mu^{+}#mu^{-} - rescaled MFT momentum, #mu^{-} top and #mu^{+} bottom", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/same-event/invariantMass_pT_ScaledMftKine_GlobalMatchesCuts_BB", "M_{#mu^{+}#mu^{-}} - rescaled MFT momentum, bottom - bottom", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_pT_ScaledMftKine_GlobalMatchesCuts", "M_{#mu^{+}#mu^{-}} - rescaled MFT momentum", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_pT_ScaledMftKine_GlobalMatchesCuts_TT", "M_{#mu^{+}#mu^{-}} - rescaled MFT momentum, top - top", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_pT_ScaledMftKine_GlobalMatchesCuts_TB", "M_{#mu^{+}#mu^{-}} - rescaled MFT momentum, top - bottom or bottom - top", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_pT_ScaledMftKine_GlobalMatchesCuts_TPBN", "#mu^{+}#mu^{-} - rescaled MFT momentum, #mu^{+} top and #mu^{-} bottom", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_pT_ScaledMftKine_GlobalMatchesCuts_TNBP", "#mu^{+}#mu^{-} - rescaled MFT momentum, #mu^{-} top and #mu^{+} bottom", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      registryDimuon.add("dimuon/mixed-event/invariantMass_pT_ScaledMftKine_GlobalMatchesCuts_BB", "M_{#mu^{+}#mu^{-}} - rescaled MFT momentum, bottom - bottom", {HistType::kTH2F, {invMassAxis2D, pTAxis2D}});
      // combinations with sub-leading matches
      registryDimuon.add("dimuon/same-event/invariantMass_GlobalMuonKine_GlobalMatchesCuts_leading_subleading", "#mu^{+}#mu^{-} invariant mass", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/same-event/invariantMassFull_GlobalMuonKine_GlobalMatchesCuts_leading_subleading", "#mu^{+}#mu^{-} invariant mass", {HistType::kTH1F, {invMassAxisFull}});
      registryDimuon.add("dimuon/same-event/invariantMass_GlobalMuonKine_GlobalMatchesCuts_subleading_leading", "#mu^{+}#mu^{-} invariant mass", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/same-event/invariantMassFull_GlobalMuonKine_GlobalMatchesCuts_subleading_leading", "#mu^{+}#mu^{-} invariant mass", {HistType::kTH1F, {invMassAxisFull}});
      registryDimuon.add("dimuon/same-event/invariantMass_GlobalMuonKine_GlobalMatchesCuts_subleading_subleading", "#mu^{+}#mu^{-} invariant mass", {HistType::kTH1F, {invMassAxis}});
      registryDimuon.add("dimuon/same-event/invariantMassFull_GlobalMuonKine_GlobalMatchesCuts_subleading_subleading", "#mu^{+}#mu^{-} invariant mass", {HistType::kTH1F, {invMassAxisFull}});
      // invariant mass correlations
      registryDimuon.add("dimuon/same-event/invariantMass_MuonKine_vs_GlobalMuonKine", "M_{#mu^{+}#mu^{-}} - muon tracks vs. global tracks", {HistType::kTH2F, {invMassCorrelationAxis, invMassCorrelationAxis}});
      registryDimuon.add("dimuon/same-event/invariantMass_ScaledMftKine_vs_GlobalMuonKine", "M_{#mu^{+}#mu^{-}} - rescaled MFT tracks vs. global tracks", {HistType::kTH2F, {invMassCorrelationAxis, invMassCorrelationAxis}});
      registryDimuon.add("dimuon/same-event/invariantMass_GlobalMuonKine_subleading_vs_leading", "M_{#mu^{+}#mu^{-}} - subleading vs. leading matches", {HistType::kTH2F, {invMassCorrelationAxis, invMassCorrelationAxis}});

      // pseudorapidity (only for MCH acceptance cuts)
      // MCH-MID tracks with MCH acceptance cuts
      registryDimuon.add("dimuon/same-event/rapPair_MuonKine_MuonCuts", "#eta of dimuon pair", {HistType::kTH1F, {yPairAxis}});
      // -- Mass and eta
      registryDimuon.add("dimuon/same-event/invariantMass_rapPair_MuonKine_MuonCuts", "#mu^{+}#mu^{-} invariant mass", {HistType::kTH2F, {invMassAxis2D, yPairAxis}});
      // -- pT and eta
      registryDimuon.add("dimuon/same-event/pT_rapPair_MuonKine_MuonCuts", "#mu^{+}#mu^{-} p_{T} and #eta", {HistType::kTH2F, {pTAxis2D, yPairAxis}});
      // MCH-MID tracks with MCH acceptance cuts and combinations from the top and bottom halfs of MCH
      registryDimuon.add("dimuon/same-event/rapPair_MuonKine_MuonCuts_TT", "#eta of dimuon pair, top-top", {HistType::kTH1F, {yPairAxis}});
      registryDimuon.add("dimuon/same-event/rapPair_MuonKine_MuonCuts_TB", "#eta of dimuon pair, top-bottom or bottom-top", {HistType::kTH1F, {yPairAxis}});
      registryDimuon.add("dimuon/same-event/rapPair_MuonKine_MuonCuts_BB", "#eta of dimuon pair, bottom-bottom", {HistType::kTH1F, {yPairAxis}});
      // -- Mass and eta
      registryDimuon.add("dimuon/same-event/invariantMass_rapPair_MuonKine_MuonCuts_TT", "#mu^{+}#mu^{-} invariant mass, top-top", {HistType::kTH2F, {invMassAxis2D, yPairAxis}});
      registryDimuon.add("dimuon/same-event/invariantMass_rapPair_MuonKine_MuonCuts_TB", "#mu^{+}#mu^{-} invariant mass, top-bottom or bottom-top", {HistType::kTH2F, {invMassAxis2D, yPairAxis}});
      registryDimuon.add("dimuon/same-event/invariantMass_rapPair_MuonKine_MuonCuts_BB", "#mu^{+}#mu^{-} invariant mass, bottom-bottom", {HistType::kTH2F, {invMassAxis2D, yPairAxis}});
      // -- pT and eta
      registryDimuon.add("dimuon/same-event/pT_rapPair_MuonKine_MuonCuts_TT", "#mu^{+}#mu^{-} p_{T} and #eta, top-top", {HistType::kTH2F, {pTAxis2D, yPairAxis}});
      registryDimuon.add("dimuon/same-event/pT_rapPair_MuonKine_MuonCuts_TB", "#mu^{+}#mu^{-} p_{T} and #eta, top-bottom or bottom-top", {HistType::kTH2F, {pTAxis2D, yPairAxis}});
      registryDimuon.add("dimuon/same-event/pT_rapPair_MuonKine_MuonCuts_BB", "#mu^{+}#mu^{-} p_{T} and #eta, bottom-bottom", {HistType::kTH2F, {pTAxis2D, yPairAxis}});
      // MCH-MID tracks with MCH acceptance cuts and combinations from the left and right halfs of MCH
      registryDimuon.add("dimuon/same-event/rapPair_MuonKine_MuonCuts_LL", "#eta of dimuon pair, left-left", {HistType::kTH1F, {yPairAxis}});
      registryDimuon.add("dimuon/same-event/rapPair_MuonKine_MuonCuts_LR", "#eta of dimuon pair, left-right or right-left", {HistType::kTH1F, {yPairAxis}});
      registryDimuon.add("dimuon/same-event/rapPair_MuonKine_MuonCuts_RR", "#eta of dimuon pair, right-right", {HistType::kTH1F, {yPairAxis}});
      // -- Mass and eta
      registryDimuon.add("dimuon/same-event/invariantMass_rapPair_MuonKine_MuonCuts_LL", "#mu^{+}#mu^{-} invariant mass, left-left", {HistType::kTH2F, {invMassAxis2D, yPairAxis}});
      registryDimuon.add("dimuon/same-event/invariantMass_rapPair_MuonKine_MuonCuts_LR", "#mu^{+}#mu^{-} invariant mass, left-right or right-left", {HistType::kTH2F, {invMassAxis2D, yPairAxis}});
      registryDimuon.add("dimuon/same-event/invariantMass_rapPair_MuonKine_MuonCuts_RR", "#mu^{+}#mu^{-} invariant mass, right-right", {HistType::kTH2F, {invMassAxis2D, yPairAxis}});
      // -- pT and eta
      // registryDimuon.add("dimuon/same-event/pT_rapPair_MuonKine_MuonCuts_LL", "#mu^{+}#mu^{-} p_{T} and #eta, left-left", {HistType::kTH2F, {pTAxis2D, yPairAxis}});
      // registryDimuon.add("dimuon/same-event/pT_rapPair_MuonKine_MuonCuts_LR", "#mu^{+}#mu^{-} p_{T} and #eta, left-right or right-left", {HistType::kTH2F, {pTAxis2D, yPairAxis}});
      // registryDimuon.add("dimuon/same-event/pT_rapPair_MuonKine_MuonCuts_RR", "#mu^{+}#mu^{-} p_{T} and #eta, right-right", {HistType::kTH2F, {pTAxis2D, yPairAxis}});
    }
  }

  void doTransformMFT(o2::mch::TrackParam& track)
  {
    double zCH10 = -1437.6;
    double z = track.getZ();
    // double dZ = zMCH - z;
    double x = track.getNonBendingCoor();
    double y = track.getBendingCoor();
    double xSlope = track.getNonBendingSlope();
    double ySlope = track.getBendingSlope();

    double xShiftMCH = (y > 0) ? 0.8541 : -1.5599;
    double xCorrection = xShiftMCH * z / zCH10;
    track.setNonBendingCoor(x + xCorrection);
    double xSlopeCorrection = xShiftMCH / zCH10;
    track.setNonBendingSlope(xSlope + xSlopeCorrection);

    double yShiftMCH = (y > 0) ? 3.0311 : 0.7588;
    double yCorrection = yShiftMCH * z / zCH10;
    track.setBendingCoor(y + yCorrection);
    double ySlopeCorrection = yShiftMCH / zCH10;
    track.setBendingSlope(ySlope + ySlopeCorrection);
  }

  template <bool GlobalFwdFillMap, typename TTrack>
  void TransformMFT(TTrack& track)
  {
    if constexpr (static_cast<bool>(GlobalFwdFillMap)) {
      auto mchTrack = mMatching.FwdtoMCH(track);
      doTransformMFT(mchTrack);

      auto transformedTrack = mMatching.MCHtoFwd(mchTrack);
      track.setParameters(transformedTrack.getParameters());
      track.setZ(transformedTrack.getZ());
      track.setCovariances(transformedTrack.getCovariances());
    } else {
      o2::dataformats::GlobalFwdTrack fwdtrack;
      fwdtrack.setParameters(track.getParameters());
      fwdtrack.setZ(track.getZ());
      fwdtrack.setCovariances(track.getCovariances());
      auto mchTrack = mMatching.FwdtoMCH(fwdtrack);
      doTransformMFT(mchTrack);

      auto transformedTrack = mMatching.MCHtoFwd(mchTrack);
      track.setParameters(transformedTrack.getParameters());
      track.setZ(transformedTrack.getZ());
      track.setCovariances(transformedTrack.getCovariances());
    }
  }

  int GetDetElemId(int iDetElemNumber)
  {
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

  int GetQuadrantPhi(double phi)
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

  template <typename TTrack>
  int GetQuadrantTrack(TTrack const& track)
  {
    double phi = static_cast<float>(track.phi()) * 180 / TMath::Pi();
    return GetQuadrantPhi(phi);
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

  template <typename Var>
  bool pDCACut(Var const& fgValues, Var const& fgValuesPV, double nSigmaPDCA)
  {
    static const double sigmaPDCA23 = 80.;
    static const double sigmaPDCA310 = 54.;
    static const double relPRes = 0.0004;
    static const double slopeRes = 0.0005;

    double thetaAbs = TMath::ATan(fgValues.rabs / 505.) * TMath::RadToDeg();
    double p = fgValuesPV.p;

    double pDCA = fgValues.pDca;
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

  template <typename Var>
  bool IsMixedEvent(Var const& fgValues1, Var const& fgValues2)
  {
    if (fgValues1.bc == fgValues2.bc) {
      return false;
    }

    uint64_t bcDiff = (fgValues2.bc > fgValues1.bc) ? (fgValues2.bc - fgValues1.bc) : (fgValues1.bc - fgValues2.bc);
    // in the event mixing case, we require a minimum BC gap between the collisions
    if (bcDiff < configMixing.fEventMinDeltaBc)
      return false;

    // we also require that the collisions have similar Z positions and multiplicity of MFT tracks
    if (std::fabs(fgValues2.z - fgValues1.z) > configMixing.fEventMaxDeltaVtxZ) {
      return false;
    }

    if (std::abs(fgValues2.multMFT - fgValues1.multMFT) > configMixing.fEventMaxDeltaNMFT) {
      return false;
    }

    return true;
  }

  template <typename Var>
  bool IsGoodMuon(Var const& fgValues, Var const& fgValuesPV, float fTrackChi2MchUp, float fPMchLow, float fPtMchLow, float fEtaMchLow, float fEtaMchUp, float fRabsLow, float fRabsUp, float fSigmaPdcaUp)
  {
    // chi2 cut
    if (fgValues.chi2 > fTrackChi2MchUp)
      return false;

    // momentum cut
    if (fgValues.p < fPMchLow) {
      return false; // skip low-momentum tracks
    }

    // transverse momentum cut
    if (fgValues.pT < fPtMchLow) {
      return false; // skip low-momentum tracks
    }

    // Eta cut
    if ((fgValues.eta < fEtaMchLow || fgValues.eta > fEtaMchUp)) {
      return false;
    }

    // RAbs cut
    if ((fgValues.rabs < fRabsLow || fgValues.rabs > fRabsUp)) {
      return false;
    }

    // pDCA cut
    if (!pDCACut(fgValues, fgValuesPV, fSigmaPdcaUp)) {
      return false;
    }

    return true;
  }

  template <typename Var>
  bool IsGoodMuon(Var const& fgValues, Var const& fgValuesPV)
  {
    return IsGoodMuon(fgValues, fgValuesPV, configMuons.fTrackChi2MchUp, configMuons.fPMchLow, configMuons.fPtMchLow, configMuons.fEtaMchLow, configMuons.fEtaMchUp, configMuons.fRabsLow, configMuons.fRabsUp, configMuons.fSigmaPdcaUp);
  }

  template <typename Var>
  bool IsGoodGlobalMuon(Var const& fgValues, Var const& fgValuesPV)
  {
    return IsGoodMuon(fgValues, fgValuesPV, configMuons.fTrackChi2MchUp, configMuons.fPMchLow, configMuons.fPtMchLow, configMFTs.fEtaMftLow, configMFTs.fEtaMftUp, configMuons.fRabsLow, configMuons.fRabsUp, configMuons.fSigmaPdcaUp);
  }

  template <typename Var>
  bool IsGoodMFT(Var const& fgValues, float fTrackChi2MftUp, int fTrackNClustMftLow)
  {
    // chi2 cut
    if (fgValues.chi2 > fTrackChi2MftUp) {
      return false;
    }

    // number of clusters cut
    if (fgValues.nClusters < fTrackNClustMftLow) {
      return false;
    }

    return true;
  }

  template <typename Var>
  bool IsGoodGlobalMatching(Var const& fgValues, float fTrackChi2MftUp, int fTrackNClustMftLow, float fMatchingChi2MftMchUp)
  {
    if (!IsGoodMFT(fgValues, fTrackChi2MftUp, fTrackNClustMftLow)) {
      return false;
    }

    if (fgValues.chi2matching > fMatchingChi2MftMchUp) {
      return false;
    }

    return true;
  }

  template <typename Var>
  bool IsGoodGlobalMatching(Var const& fgValues)
  {
    return IsGoodGlobalMatching(fgValues, configMFTs.fTrackChi2MftUp, configMFTs.fTrackNClustMftLow, fMatchingChi2MftMchUp);
  }

  template <typename VarT>
  double GetMuMuInvariantMass(VarT const& track1, VarT const& track2)
  {
    ROOT::Math::PxPyPzMVector muon1{
      track1.px,
      track1.py,
      track1.pz,
      o2::constants::physics::MassMuon};

    ROOT::Math::PxPyPzMVector muon2{
      track2.px,
      track2.py,
      track2.pz,
      o2::constants::physics::MassMuon};

    auto dimuon = muon1 + muon2;

    return dimuon.M();
  }

  template <typename VarT>
  double GetMuMuPt(VarT const& track1, VarT const& track2)
  {
    ROOT::Math::PxPyPzMVector muon1{
      track1.px,
      track1.py,
      track1.pz,
      o2::constants::physics::MassMuon};

    ROOT::Math::PxPyPzMVector muon2{
      track2.px,
      track2.py,
      track2.pz,
      o2::constants::physics::MassMuon};

    auto dimuon = muon1 + muon2;

    return dimuon.Pt();
  }

  template <typename VarT>
  double GetMuMuRap(VarT const& track1, VarT const& track2)
  {
    ROOT::Math::PxPyPzMVector muon1{
      track1.px,
      track1.py,
      track1.pz,
      o2::constants::physics::MassMuon};

    ROOT::Math::PxPyPzMVector muon2{
      track2.px,
      track2.py,
      track2.pz,
      o2::constants::physics::MassMuon};

    auto dimuon = muon1 + muon2;

    return dimuon.Y();
  }

  template <typename TMuons, typename TCandidates>
  void GetMuonPairs(TMuons const& muons, TCandidates const& matchingCandidates, const std::map<uint64_t, VarColl>& collisionInfos,
                    std::vector<MuonPair>& muonPairs,
                    std::vector<GlobalMuonPair>& globalMuonPairs)
  {
    // muon tracks - outer loop over collisions
    for (const auto& [collisionIndex1, collisionInfo1] : collisionInfos) {

      // outer loop over muon tracks
      auto muonCollision1 = muons.sliceBy(fwdtracksPerCollision, collisionInfo1.globalIndex);
      for (const auto& muon1 : muonCollision1) {

        if (muon1.trackType() <= 2) {
          continue;
        }
        auto mchIndex1 = muon1.globalIndex();

        // inner loop over collisions
        for (const auto& [collisionIndex2, collisionInfo2] : collisionInfos) {
          // avoid double-counting of collisions
          if (collisionIndex2 < collisionIndex1)
            continue;

          bool sameEvent = (collisionIndex1 == collisionIndex2);
          bool mixedEvent = IsMixedEvent(collisionInfo1, collisionInfo2);

          if (!sameEvent && !mixedEvent)
            continue;

          // inner loop over muon tracks
          auto muonCollision2 = muons.sliceBy(fwdtracksPerCollision, collisionInfo2.globalIndex);
          for (const auto& muon2 : muonCollision2) {
            if (muon2.trackType() <= 2) {
              continue;
            }
            auto mchIndex2 = muon2.globalIndex();

            // avoid double-counting of muon pairs if we are not mixing events
            if (sameEvent && mchIndex2 <= mchIndex1)
              continue;

            MuonPair muonPair{{collisionIndex1, mchIndex1}, {collisionIndex2, mchIndex2}};
            muonPairs.emplace_back(muonPair);
          }
        }
      }
    }

    // global muon tracks - outer loop over collisions
    for (const auto& [collisionIndex1, collisionInfo1] : collisionInfos) {

      // outer loop over global muon tracks
      auto muonCollision1 = muons.sliceBy(fwdtracksPerCollision, collisionInfo1.globalIndex);
      for (const auto& muon1 : muonCollision1) {

        if (muon1.trackType() <= 2) {
          continue;
        }
        auto mchIndex1 = muon1.globalIndex();
        auto matchingCandidateIt1 = matchingCandidates.find(mchIndex1);
        if (matchingCandidateIt1 == matchingCandidates.end()) {
          continue;
        }

        // inner loop over collisions
        for (const auto& [collisionIndex2, collisionInfo2] : collisionInfos) {
          // avoid double-counting of collisions
          if (collisionIndex2 < collisionIndex1)
            continue;

          bool sameEvent = (collisionIndex1 == collisionIndex2);
          bool mixedEvent = IsMixedEvent(collisionInfo1, collisionInfo2);

          if (!sameEvent && !mixedEvent)
            continue;

          // outer loop over global muon tracks
          auto muonCollision2 = muons.sliceBy(fwdtracksPerCollision, collisionInfo2.globalIndex);
          for (const auto& muon2 : muonCollision2) {

            if (muon2.trackType() <= 2) {
              continue;
            }
            auto mchIndex2 = muon2.globalIndex();
            auto matchingCandidateIt2 = matchingCandidates.find(mchIndex2);
            if (matchingCandidateIt2 == matchingCandidates.end()) {
              continue;
            }

            // avoid double-counting of muon pairs if we are not mixing events
            if (sameEvent && mchIndex2 <= mchIndex1)
              continue;

            GlobalMuonPair muonPair{{collisionIndex1, matchingCandidateIt1->second}, {collisionIndex2, matchingCandidateIt2->second}};
            globalMuonPairs.emplace_back(muonPair);
          }
        }
      }
    }
  }

  template <typename TEvent, typename Var>
  void FillCollision(TEvent const& collision, Var& fgValues)
  {
    fgValues.globalIndex = collision.globalIndex();
    fgValues.x = collision.posX();
    fgValues.y = collision.posY();
    fgValues.z = collision.posZ();
    fgValues.covXX = collision.covXX();
    fgValues.covYY = collision.covYY();
  }

  template <typename Var>
  void FillTrack(mch::Track const& muon, Var& fgValues)
  {
    mch::TrackParam trackParam = mch::TrackParam(muon.first());
    auto proptrack = mMatching.MCHtoFwd(trackParam);

    fgValues.pT = proptrack.getPt();
    fgValues.x = proptrack.getX();
    fgValues.y = proptrack.getY();
    fgValues.z = proptrack.getZ();
    fgValues.eta = proptrack.getEta();
    fgValues.tgl = proptrack.getTgl();
    fgValues.phi = proptrack.getPhi();

    fgValues.p = proptrack.getP();
    fgValues.px = proptrack.getPx();
    fgValues.py = proptrack.getPy();
    fgValues.pz = proptrack.getPz();

    fgValues.chi2 = trackParam.getTrackChi2();
    fgValues.nClusters = muon.getNClusters();
    fgValues.sign = trackParam.getCharge();
  }

  template <bool MuonFillMap, typename TTrack, typename Var>
  void FillTrack(TTrack const& muon, Var& fgValues)
  {
    fgValues.collisionId = muon.collisionId();
    fgValues.globalIndex = muon.globalIndex();
    fgValues.trackTime = muon.trackTime();

    fgValues.pT = muon.pt();
    fgValues.x = muon.x();
    fgValues.y = muon.y();
    fgValues.z = muon.z();
    fgValues.eta = muon.eta();
    fgValues.tgl = muon.tgl();
    fgValues.phi = muon.phi();

    fgValues.p = muon.p();
    fgValues.px = muon.px();
    fgValues.py = muon.py();
    fgValues.pz = muon.pz();

    fgValues.chi2 = muon.chi2();
    fgValues.nClusters = muon.nClusters();
    fgValues.sign = muon.sign();

    if constexpr (static_cast<bool>(MuonFillMap)) {
      // Direct info from AO2D without re-propagation
      fgValues.pDca = muon.pDca();
      fgValues.rabs = muon.rAtAbsorberEnd();
      fgValues.trackType = muon.trackType();
    }
  }

  template <typename TMCHTrack, typename TFwdCls, typename Var>
  bool FillClusters(TMCHTrack const& muon, TFwdCls const& mchcls, Var& fgValues, TrackRealigned& convertedTrack)
  {
    int removable = 0;
    auto clustersSliced = mchcls.sliceBy(perMuon, muon.globalIndex()); // Slice clusters by muon id
    std::vector<std::vector<float>> posClusters;

    int clIndex = -1;
    // Get re-aligned clusters associated to current track
    for (auto const& cluster : clustersSliced) {
      clIndex += 1;

      math_utils::Point3D<double> local;
      math_utils::Point3D<double> master;

      mch::Cluster* clusterMCH = new mch::Cluster();
      master.SetXYZ(cluster.x(), cluster.y(), cluster.z());

      if (configRealign.fDoRealign) {
        // Transformation from reference geometry frame to new geometry frame
        transformRef[cluster.deId()].MasterToLocal(master, local);
        transformNew[cluster.deId()].LocalToMaster(local, master);
      }

      clusterMCH->x = master.x();
      clusterMCH->y = master.y();
      clusterMCH->z = master.z();

      uint32_t ClUId = mch::Cluster::buildUniqueId(static_cast<int>(cluster.deId() / 100) - 1, cluster.deId(), clIndex);
      clusterMCH->uid = ClUId;
      clusterMCH->ex = cluster.isGoodX() ? 0.2 : 10.0;
      clusterMCH->ey = cluster.isGoodY() ? 0.2 : 10.0;

      // Fill temporary values
      std::vector<float> posCls = {clusterMCH->x, clusterMCH->y, clusterMCH->z};
      std::vector<float> eCls = {clusterMCH->ex, clusterMCH->ey};
      posClusters.emplace_back(posCls);
      fgValues.errorClusters.emplace_back(eCls);
      fgValues.DEIDs.emplace_back(cluster.deId());

      // Add transformed cluster into temporary variable
      convertedTrack.createParamAtCluster(*clusterMCH);
    }

    if (configRealign.fDoRealign) {
      // Refit the re-aligned track
      if (convertedTrack.getNClusters() != 0) {
        removable = RemoveTrack(convertedTrack);
      } else {
        LOGF(fatal, "Muon track %d has no associated clusters.", muon.globalIndex());
      }

      for (auto it = convertedTrack.begin(); it != convertedTrack.end(); it++) {
        std::vector<float> pos = {static_cast<float>(it->getNonBendingCoor()), static_cast<float>(it->getBendingCoor()), static_cast<float>(it->getZ())};
        fgValues.posClusters.emplace_back(pos);
      }

    } else {
      fgValues.posClusters = posClusters;
    }

    return !removable;
  }

  template <typename TMuon, typename TMCH, typename TMap>
  void FillMatchingCandidates(TMuon const& muon, TMCH const& mchtrack, TMap& matchingCandidates)
  {
    uint64_t muonId = muon.globalIndex();
    uint64_t mchId = mchtrack.globalIndex();

    //// Save matching candidates index pairs
    auto matchingCandidateIt = matchingCandidates.find(mchId);
    if (matchingCandidateIt != matchingCandidates.end()) {
      matchingCandidateIt->second.push_back(muonId);
    } else {
      matchingCandidates[mchId].push_back(muonId);
    }
  }

  template <typename VarC, typename VarT>
  void FillPropagation(mch::Track const& muon, VarC const& collision, VarT& fgValues, int endPoint = kToVtx, int endZ = 0)
  {
    o2::dataformats::GlobalFwdTrack propmuon;
    mch::TrackParam trackParam = mch::TrackParam(muon.first());
    fgValues.chi2 = trackParam.getTrackChi2();
    fgValues.nClusters = muon.getNClusters();
    fgValues.sign = trackParam.getCharge();

    if (endPoint == kToVtx) {
      o2::mch::TrackExtrap::extrapToVertex(trackParam, collision.x, collision.y, collision.z, collision.covXX, collision.covYY);
    }
    if (endPoint == kToDCA) {
      o2::mch::TrackExtrap::extrapToVertexWithoutBranson(trackParam, collision.z);
    }
    if (endPoint == kToAbsEnd) {
      o2::mch::TrackExtrap::extrapToZ(trackParam, zAtAbsEnd);
    }
    if (endPoint == kToZ) {
      o2::mch::TrackExtrap::extrapToZ(trackParam, endZ);
    }

    auto proptrack = mMatching.MCHtoFwd(trackParam);
    propmuon.setParameters(proptrack.getParameters());
    propmuon.setZ(proptrack.getZ());
    propmuon.setCovariances(proptrack.getCovariances());

    //// Fill propagation informations
    if (endPoint == kToVtx || endPoint == kToZ) {
      fgValues.pT = propmuon.getPt();
      fgValues.x = propmuon.getX();
      fgValues.y = propmuon.getY();
      fgValues.z = propmuon.getZ();
      fgValues.eta = propmuon.getEta();
      fgValues.tgl = propmuon.getTgl();
      fgValues.phi = propmuon.getPhi();

      fgValues.p = propmuon.getP();
      fgValues.px = propmuon.getP() * std::sin(M_PI / 2 - std::atan(propmuon.getTgl())) * std::cos(propmuon.getPhi());
      fgValues.py = propmuon.getP() * std::sin(M_PI / 2 - std::atan(propmuon.getTgl())) * std::sin(propmuon.getPhi());
      fgValues.pz = propmuon.getP() * std::cos(M_PI / 2 - std::atan(propmuon.getTgl()));
    }

    if (endPoint == kToDCA) {
      fgValues.dcaX = (propmuon.getX() - collision.x);
      fgValues.dcaY = (propmuon.getY() - collision.y);
      float dcaXY = std::sqrt(fgValues.dcaX * fgValues.dcaX + fgValues.dcaY * fgValues.dcaY);
      fgValues.dcaXY = dcaXY;

      mch::TrackParam trackParam = mch::TrackParam(muon.first());
      float p = trackParam.p();
      fgValues.pDca = p * dcaXY;
    }

    if (endPoint == kToAbsEnd) {
      double xAbs = propmuon.getX();
      double yAbs = propmuon.getY();
      fgValues.rabs = std::sqrt(xAbs * xAbs + yAbs * yAbs);
    }
  }

  template <bool MuonFillMap, bool Scaled = false, typename TTrack, typename VarC, typename VarT>
  void FillPropagation(TTrack const& muon, VarC const& collision, VarT const& fgValuesMCH, VarT& fgValues, int endPoint = kToVtx, int endZ = 0)
  {
    o2::dataformats::GlobalFwdTrack propmuon;
    double chi2 = muon.chi2();
    fgValues.chi2 = chi2;
    fgValues.nClusters = muon.nClusters();
    fgValues.sign = muon.sign();

    if constexpr (static_cast<bool>(MuonFillMap)) {
      o2::dataformats::GlobalFwdTrack track;

      SMatrix5 tpars(muon.x(), muon.y(), muon.phi(), muon.tgl(), muon.signed1Pt());
      std::vector<double> v1{muon.cXX(), muon.cXY(), muon.cYY(), muon.cPhiX(), muon.cPhiY(),
                             muon.cPhiPhi(), muon.cTglX(), muon.cTglY(), muon.cTglPhi(), muon.cTglTgl(),
                             muon.c1PtX(), muon.c1PtY(), muon.c1PtPhi(), muon.c1PtTgl(), muon.c1Pt21Pt2()};
      SMatrix55 tcovs(v1.begin(), v1.end());
      o2::track::TrackParCovFwd fwdtrack{muon.z(), tpars, tcovs, chi2};

      track.setParameters(tpars);
      track.setZ(fwdtrack.getZ());
      track.setCovariances(tcovs);
      auto mchTrack = mMatching.FwdtoMCH(track);

      if (endPoint == kToVtx) {
        o2::mch::TrackExtrap::extrapToVertex(mchTrack, collision.x, collision.y, collision.z, collision.covXX, collision.covYY);
      }
      if (endPoint == kToDCA) {
        o2::mch::TrackExtrap::extrapToVertexWithoutBranson(mchTrack, collision.z);
      }
      if (endPoint == kToAbsEnd) {
        o2::mch::TrackExtrap::extrapToZ(mchTrack, zAtAbsEnd);
      }
      if (endPoint == kToZ) {
        o2::mch::TrackExtrap::extrapToZ(mchTrack, endZ);
      }

      auto proptrack = mMatching.MCHtoFwd(mchTrack);
      propmuon.setParameters(proptrack.getParameters());
      propmuon.setZ(proptrack.getZ());
      propmuon.setCovariances(proptrack.getCovariances());

    } else {

      o2::dataformats::GlobalFwdTrack track;

      if constexpr (static_cast<bool>(Scaled)) {
        double pMCH = fgValuesMCH.p;
        int sign = fgValuesMCH.sign;

        double px = pMCH * std::sin(M_PI / 2 - std::atan(muon.tgl())) * std::cos(muon.phi());
        double py = pMCH * std::sin(M_PI / 2 - std::atan(muon.tgl())) * std::sin(muon.phi());
        // double pz = pMCH * cos(M_PI / 2 - atan(mft.tgl()));
        double pt = std::sqrt(std::pow(px, 2) + std::pow(py, 2));

        double chi2 = muon.chi2();
        double signed1Pt = endPoint == kToDCA ? muon.signed1Pt() : sign / pt;
        SMatrix5 tpars(muon.x(), muon.y(), muon.phi(), muon.tgl(), signed1Pt);
        std::vector<double> v1{0, 0, 0, 0, 0,
                               0, 0, 0, 0, 0,
                               0, 0, 0, 0, 0};
        SMatrix55 tcovs(v1.begin(), v1.end());
        o2::track::TrackParCovFwd fwdtrack{muon.z(), tpars, tcovs, chi2};
        track.setParameters(tpars);
        track.setZ(fwdtrack.getZ());
        track.setCovariances(tcovs);
      } else {
        SMatrix5 tpars(muon.x(), muon.y(), muon.phi(), muon.tgl(), muon.signed1Pt());
        std::vector<double> v1{muon.cXX(), muon.cXY(), muon.cYY(), muon.cPhiX(), muon.cPhiY(),
                               muon.cPhiPhi(), muon.cTglX(), muon.cTglY(), muon.cTglPhi(), muon.cTglTgl(),
                               muon.c1PtX(), muon.c1PtY(), muon.c1PtPhi(), muon.c1PtTgl(), muon.c1Pt21Pt2()};
        SMatrix55 tcovs(v1.begin(), v1.end());
        o2::track::TrackParCovFwd fwdtrack{muon.z(), tpars, tcovs, chi2};
        track.setParameters(tpars);
        track.setZ(fwdtrack.getZ());
        track.setCovariances(tcovs);
      }

      if (endPoint == kToVtx) {
        if (fEnableMFTAlignmentCorrections) {
          TransformMFT<0>(track);
        }
        auto geoMan = o2::base::GeometryManager::meanMaterialBudget(muon.x(), muon.y(), muon.z(), collision.x, collision.y, collision.z);
        auto x2x0 = static_cast<float>(geoMan.meanX2X0);
        track.propagateToVtxhelixWithMCS(collision.z, {collision.x, collision.y}, {collision.covXX, collision.covYY}, Bz, x2x0);
      }
      if (endPoint == kToDCA) {
        if (fEnableMFTAlignmentCorrections) {
          TransformMFT<0>(track);
        }
        track.propagateToZ(collision.z, Bz);
      }
      if (endPoint == kToZ) {
        auto mchTrackExt = mMatching.FwdtoMCH(track);
        if (fEnableMFTAlignmentCorrections) {
          doTransformMFT(mchTrackExt);
        }
        o2::mch::TrackExtrap::extrapToZ(mchTrackExt, endZ);
        track = mMatching.MCHtoFwd(mchTrackExt);
      }

      propmuon.setParameters(track.getParameters());
      propmuon.setZ(track.getZ());
      propmuon.setCovariances(track.getCovariances());
      if (endPoint == kToDCA) {
        fgValues.dcaX = (propmuon.getX() - collision.x);
        fgValues.dcaY = (propmuon.getY() - collision.y);
        float dcaXY = std::sqrt(fgValues.dcaX * fgValues.dcaX + fgValues.dcaY * fgValues.dcaY);
        fgValues.dcaXY = dcaXY;
      }
    }

    //// Fill propagation informations
    if (endPoint == kToVtx || endPoint == kToZ) {
      fgValues.pT = propmuon.getPt();
      fgValues.x = propmuon.getX();
      fgValues.y = propmuon.getY();
      fgValues.z = propmuon.getZ();
      fgValues.eta = propmuon.getEta();
      fgValues.tgl = propmuon.getTgl();
      fgValues.phi = propmuon.getPhi();

      fgValues.p = propmuon.getP();
      fgValues.px = propmuon.getP() * std::sin(M_PI / 2 - std::atan(propmuon.getTgl())) * std::cos(propmuon.getPhi());
      fgValues.py = propmuon.getP() * std::sin(M_PI / 2 - std::atan(propmuon.getTgl())) * std::sin(propmuon.getPhi());
      fgValues.pz = propmuon.getP() * std::cos(M_PI / 2 - std::atan(propmuon.getTgl()));
    }

    if (endPoint == kToDCA) {
      fgValues.dcaX = (propmuon.getX() - collision.x);
      fgValues.dcaY = (propmuon.getY() - collision.y);
      float dcaXY = std::sqrt(fgValues.dcaX * fgValues.dcaX + fgValues.dcaY * fgValues.dcaY);
      fgValues.dcaXY = dcaXY;

      if constexpr (static_cast<bool>(MuonFillMap)) {
        float p = muon.p();
        fgValues.pDca = p * dcaXY;
      }
    }

    if (endPoint == kToAbsEnd) {
      double xAbs = propmuon.getX();
      double yAbs = propmuon.getY();
      fgValues.rabs = std::sqrt(xAbs * xAbs + yAbs * yAbs);
    }
  }

  template <typename TTrack, typename Var>
  void FillMatching(TTrack const& track, Var& fgValuesMCH, Var& fgValuesMFT)
  {
    fgValuesMCH.chi2matching = track.chi2MatchMCHMID();
    fgValuesMFT.chi2matching = track.chi2MatchMCHMFT();
  }

  template <bool MuonFillMap, bool GlobalMuonFillMap, bool GlobalMatchingFillMap, typename Var, typename VarVector>
  void FillMuonHistograms(Var const& fgValuesMCH, Var const& fgValuesMCHpv, Var const& fgValuesMFT, Var const& fgValuesGlobal, VarVector const& fgValuesCandidates)
  {
    if constexpr (static_cast<bool>(MuonFillMap)) {
      // Muon histograms
      if (IsGoodMuon(fgValuesMCH, fgValuesMCHpv, 1.E10, configMuons.fPMchLow, configMuons.fPtMchLow, configMuons.fEtaMchLow, configMuons.fEtaMchUp, configMuons.fRabsLow, configMuons.fRabsUp, configMuons.fSigmaPdcaUp)) {
        registry.get<TH1>(HIST("muons/TrackChi2"))->Fill(fgValuesMCH.chi2);
      }
      if (IsGoodMuon(fgValuesMCH, fgValuesMCHpv, configMuons.fTrackChi2MchUp, 0., configMuons.fPtMchLow, configMuons.fEtaMchLow, configMuons.fEtaMchUp, configMuons.fRabsLow, configMuons.fRabsUp, configMuons.fSigmaPdcaUp)) {
        registry.get<TH1>(HIST("muons/TrackP"))->Fill(fgValuesMCH.p);
      }
      if (IsGoodMuon(fgValuesMCH, fgValuesMCHpv, configMuons.fTrackChi2MchUp, configMuons.fPMchLow, 0., configMuons.fEtaMchLow, configMuons.fEtaMchUp, configMuons.fRabsLow, configMuons.fRabsUp, configMuons.fSigmaPdcaUp)) {
        registry.get<TH1>(HIST("muons/TrackPt"))->Fill(fgValuesMCH.pT);
      }
      if (IsGoodMuon(fgValuesMCH, fgValuesMCHpv, configMuons.fTrackChi2MchUp, configMuons.fPMchLow, configMuons.fPtMchLow, -1.E10, 1.E10, configMuons.fRabsLow, configMuons.fRabsUp, configMuons.fSigmaPdcaUp)) {
        registry.get<TH1>(HIST("muons/TrackEta"))->Fill(fgValuesMCH.eta);
      }
      if (IsGoodMuon(fgValuesMCH, fgValuesMCHpv, configMuons.fTrackChi2MchUp, configMuons.fPMchLow, configMuons.fPtMchLow, configMuons.fEtaMchLow, configMuons.fEtaMchUp, 0., 1.E10, configMuons.fSigmaPdcaUp)) {
        registry.get<TH1>(HIST("muons/TrackRabs"))->Fill(fgValuesMCH.rabs);
      }
      if (IsGoodMuon(fgValuesMCH, fgValuesMCHpv, configMuons.fTrackChi2MchUp, configMuons.fPMchLow, configMuons.fPtMchLow, configMuons.fEtaMchLow, configMuons.fEtaMchUp, configMuons.fRabsLow, configMuons.fRabsUp, 1.E10)) {
        registry.get<TH1>(HIST("muons/TrackPDCA"))->Fill(fgValuesMCH.pDca);
      }
      if (IsGoodMuon(fgValuesMCH, fgValuesMCHpv, configMuons.fTrackChi2MchUp, configMuons.fPMchLow, configMuons.fPtMchLow, configMuons.fEtaMchLow, configMuons.fEtaMchUp, configMuons.fRabsLow, configMuons.fRabsUp, configMuons.fSigmaPdcaUp)) {
        registry.get<TH1>(HIST("muons/TrackPhi"))->Fill(fgValuesMCH.phi * 180.0 / TMath::Pi());
        registry.get<TH1>(HIST("muons/TrackDCA"))->Fill(std::sqrt(fgValuesMCHpv.dcaX * fgValuesMCHpv.dcaX + fgValuesMCHpv.dcaY * fgValuesMCHpv.dcaY));
      }

      // muon origin for MCH top-bottom and left-right parts
      if (IsGoodMuon(fgValuesMCH, fgValuesMCHpv, configMuons.fTrackChi2MchUp, configMuons.fPMchLow, configMuons.fPtMchLow, -1.E10, 1.E10, configMuons.fRabsLow, configMuons.fRabsUp, configMuons.fSigmaPdcaUp)) {
        int Quadrant = GetQuadrantPhi(fgValuesMCH.phi * 180.0 / TMath::Pi());
        int TopBottom = (Quadrant == 0 || Quadrant == 1) ? 0 : 1;
        int LeftRight = (Quadrant == 0 || Quadrant == 3) ? 0 : 1;
        int PosNeg = fgValuesMCH.sign > 0 ? 0 : 1;
        float eta = fgValuesMCH.eta;
        float pT = fgValuesMCH.pT;

        // same-event case
        if (PosNeg == 0) {
          registry.get<TH1>(HIST("muons/TrackEtaPos"))->Fill(eta);
          registry.get<TH2>(HIST("muons/TrackPt_TrackEtaPos"))->Fill(pT, eta);

          if (TopBottom == 0) {
            registry.get<TH1>(HIST("muons/TrackEtaPos_T"))->Fill(eta);
            registry.get<TH2>(HIST("muons/TrackPt_TrackEtaPos_T"))->Fill(pT, eta);
          } else if (TopBottom == 1) {
            registry.get<TH1>(HIST("muons/TrackEtaPos_B"))->Fill(eta);
            registry.get<TH2>(HIST("muons/TrackPt_TrackEtaPos_B"))->Fill(pT, eta);
          }

          if (LeftRight == 0) {
            registry.get<TH1>(HIST("muons/TrackEtaPos_L"))->Fill(eta);
            registry.get<TH2>(HIST("muons/TrackPt_TrackEtaPos_L"))->Fill(pT, eta);
          } else if (LeftRight == 1) {
            registry.get<TH1>(HIST("muons/TrackEtaPos_R"))->Fill(eta);
            registry.get<TH2>(HIST("muons/TrackPt_TrackEtaPos_R"))->Fill(pT, eta);
          }
        } else if (PosNeg == 1) {
          registry.get<TH1>(HIST("muons/TrackEtaNeg"))->Fill(eta);
          registry.get<TH2>(HIST("muons/TrackPt_TrackEtaNeg"))->Fill(pT, eta);
          if (TopBottom == 0) {
            registry.get<TH1>(HIST("muons/TrackEtaNeg_T"))->Fill(eta);
            registry.get<TH2>(HIST("muons/TrackPt_TrackEtaNeg_T"))->Fill(pT, eta);
          } else if (TopBottom == 1) {
            registry.get<TH1>(HIST("muons/TrackEtaNeg_B"))->Fill(eta);
            registry.get<TH2>(HIST("muons/TrackPt_TrackEtaNeg_B"))->Fill(pT, eta);
          }

          if (LeftRight == 0) {
            registry.get<TH1>(HIST("muons/TrackEtaNeg_L"))->Fill(eta);
            registry.get<TH2>(HIST("muons/TrackPt_TrackEtaNeg_L"))->Fill(pT, eta);
          } else if (LeftRight == 1) {
            registry.get<TH1>(HIST("muons/TrackEtaNeg_R"))->Fill(eta);
            registry.get<TH2>(HIST("muons/TrackPt_TrackEtaNeg_R"))->Fill(pT, eta);
          }
        }
      }
    }

    if constexpr (static_cast<bool>(GlobalMuonFillMap)) {
      // Global muon histograms
      if (IsGoodMuon(fgValuesMCH, fgValuesMCHpv, 1.E10, configMuons.fPMchLow, configMuons.fPtMchLow, configMFTs.fEtaMftLow, configMFTs.fEtaMftUp, configMuons.fRabsLow, configMuons.fRabsUp, configMuons.fSigmaPdcaUp)) {
        registry.get<TH1>(HIST("global-muons/TrackChi2"))->Fill(fgValuesMCH.chi2);
      }
      if (IsGoodMuon(fgValuesMCH, fgValuesMCHpv, configMuons.fTrackChi2MchUp, 0., configMuons.fPtMchLow, configMFTs.fEtaMftLow, configMFTs.fEtaMftUp, configMuons.fRabsLow, configMuons.fRabsUp, configMuons.fSigmaPdcaUp)) {
        registry.get<TH1>(HIST("global-muons/TrackP"))->Fill(fgValuesMCH.p);

        // Momentum correlations
        registry.get<TH2>(HIST("global-muons/MomentumCorrelation_Global_vs_Muon"))->Fill(fgValuesMCH.p, fgValuesGlobal.p);
        if (fgValuesMCH.p != 0) {
          registry.get<TH2>(HIST("global-muons/MomentumDifference_Global_vs_Muon"))->Fill(fgValuesMCH.p, (fgValuesGlobal.p - fgValuesMCH.p) / fgValuesMCH.p);
        }

        if (fgValuesCandidates.size() >= 2) {
          registry.get<TH2>(HIST("global-muons/MomentumCorrelation_subleading_vs_leading"))->Fill(fgValuesCandidates[0].p, fgValuesCandidates[1].p);
          if (fgValuesCandidates[0].p != 0) {
            registry.get<TH2>(HIST("global-muons/MomentumDifference_subleading_vs_leading"))->Fill(fgValuesCandidates[0].p, (fgValuesCandidates[1].p - fgValuesCandidates[0].p) / fgValuesCandidates[0].p);
          }
        }
      }
      if (IsGoodMuon(fgValuesMCH, fgValuesMCHpv, configMuons.fTrackChi2MchUp, configMuons.fPMchLow, 0., configMFTs.fEtaMftLow, configMFTs.fEtaMftUp, configMuons.fRabsLow, configMuons.fRabsUp, configMuons.fSigmaPdcaUp)) {
        registry.get<TH1>(HIST("global-muons/TrackPt"))->Fill(fgValuesMCH.pT);
      }
      if (IsGoodMuon(fgValuesMCH, fgValuesMCHpv, configMuons.fTrackChi2MchUp, configMuons.fPMchLow, configMuons.fPtMchLow, -1.E10, 1.E10, configMuons.fRabsLow, configMuons.fRabsUp, configMuons.fSigmaPdcaUp)) {
        registry.get<TH1>(HIST("global-muons/TrackEta"))->Fill(fgValuesMCH.eta);

        // Eta correlations
        registry.get<TH2>(HIST("global-muons/EtaCorrelation_Global_vs_Muon"))->Fill(fgValuesMCH.eta, fgValuesGlobal.eta);
        if (fgValuesMCH.eta != 0) {
          registry.get<TH2>(HIST("global-muons/EtaDifference_Global_vs_Muon"))->Fill(fgValuesMCH.eta, (fgValuesGlobal.eta - fgValuesMCH.eta) / fgValuesMCH.eta);
        }

        if (fgValuesCandidates.size() >= 2) {
          registry.get<TH2>(HIST("global-muons/EtaCorrelation_subleading_vs_leading"))->Fill(fgValuesCandidates[0].eta, fgValuesCandidates[1].eta);
          if (fgValuesCandidates[0].eta != 0) {
            registry.get<TH2>(HIST("global-muons/EtaDifference_subleading_vs_leading"))->Fill(fgValuesCandidates[0].eta, (fgValuesCandidates[1].eta - fgValuesCandidates[0].eta) / fgValuesCandidates[0].eta);
          }
        }
      }
      if (IsGoodMuon(fgValuesMCH, fgValuesMCHpv, configMuons.fTrackChi2MchUp, configMuons.fPMchLow, configMuons.fPtMchLow, configMFTs.fEtaMftLow, configMFTs.fEtaMftUp, 0., 1.E10, configMuons.fSigmaPdcaUp)) {
        registry.get<TH1>(HIST("global-muons/TrackRabs"))->Fill(fgValuesMCH.rabs);
      }
      if (IsGoodMuon(fgValuesMCH, fgValuesMCHpv, configMuons.fTrackChi2MchUp, configMuons.fPMchLow, configMuons.fPtMchLow, configMFTs.fEtaMftLow, configMFTs.fEtaMftUp, configMuons.fRabsLow, configMuons.fRabsUp, 1.E10)) {
        registry.get<TH1>(HIST("global-muons/TrackPDCA"))->Fill(fgValuesMCH.pDca);
      }
      if (IsGoodMuon(fgValuesMCH, fgValuesMCHpv, configMuons.fTrackChi2MchUp, configMuons.fPMchLow, configMuons.fPtMchLow, configMFTs.fEtaMftLow, configMFTs.fEtaMftUp, configMuons.fRabsLow, configMuons.fRabsUp, configMuons.fSigmaPdcaUp)) {
        registry.get<TH1>(HIST("global-muons/TrackPhi"))->Fill(fgValuesMCH.phi * 180.0 / TMath::Pi());
        registry.get<TH1>(HIST("global-muons/TrackDCA"))->Fill(std::sqrt(fgValuesMCHpv.dcaX * fgValuesMCHpv.dcaX + fgValuesMCHpv.dcaY * fgValuesMCHpv.dcaY));
      }
    }

    if constexpr (static_cast<bool>(GlobalMatchingFillMap)) {
      // Global muon histograms
      if (IsGoodMuon(fgValuesMCH, fgValuesMCHpv, 1.E10, configMuons.fPMchLow, configMuons.fPtMchLow, configMFTs.fEtaMftLow, configMFTs.fEtaMftUp, configMuons.fRabsLow, configMuons.fRabsUp, configMuons.fSigmaPdcaUp)) {
        if (IsGoodGlobalMatching(fgValuesMFT, configMFTs.fTrackChi2MftUp, configMFTs.fTrackNClustMftLow, fMatchingChi2MftMchUp)) {
          registry.get<TH1>(HIST("global-matches/TrackChi2"))->Fill(fgValuesMCH.chi2);
        }
      }
      if (IsGoodMuon(fgValuesMCH, fgValuesMCHpv, configMuons.fTrackChi2MchUp, 0., configMuons.fPtMchLow, configMFTs.fEtaMftLow, configMFTs.fEtaMftUp, configMuons.fRabsLow, configMuons.fRabsUp, configMuons.fSigmaPdcaUp)) {
        if (IsGoodGlobalMatching(fgValuesMFT, configMFTs.fTrackChi2MftUp, configMFTs.fTrackNClustMftLow, fMatchingChi2MftMchUp)) {
          registry.get<TH1>(HIST("global-matches/TrackP"))->Fill(fgValuesMCH.p);
        }
      }
      if (IsGoodMuon(fgValuesMCH, fgValuesMCHpv, configMuons.fTrackChi2MchUp, configMuons.fPMchLow, 0., configMFTs.fEtaMftLow, configMFTs.fEtaMftUp, configMuons.fRabsLow, configMuons.fRabsUp, configMuons.fSigmaPdcaUp)) {
        if (IsGoodGlobalMatching(fgValuesMFT, configMFTs.fTrackChi2MftUp, configMFTs.fTrackNClustMftLow, fMatchingChi2MftMchUp)) {
          registry.get<TH1>(HIST("global-matches/TrackPt"))->Fill(fgValuesMCH.pT);
        }
      }
      if (IsGoodMuon(fgValuesMCH, fgValuesMCHpv, configMuons.fTrackChi2MchUp, configMuons.fPMchLow, configMuons.fPtMchLow, -1.E10, 1.E10, configMuons.fRabsLow, configMuons.fRabsUp, configMuons.fSigmaPdcaUp)) {
        if (IsGoodGlobalMatching(fgValuesMFT, configMFTs.fTrackChi2MftUp, configMFTs.fTrackNClustMftLow, fMatchingChi2MftMchUp)) {
          registry.get<TH1>(HIST("global-matches/TrackEta"))->Fill(fgValuesMCH.eta);
        }
      }
      if (IsGoodMuon(fgValuesMCH, fgValuesMCHpv, configMuons.fTrackChi2MchUp, configMuons.fPMchLow, configMuons.fPtMchLow, configMFTs.fEtaMftLow, configMFTs.fEtaMftUp, 0., 1.E10, configMuons.fSigmaPdcaUp)) {
        if (IsGoodGlobalMatching(fgValuesMFT, configMFTs.fTrackChi2MftUp, configMFTs.fTrackNClustMftLow, fMatchingChi2MftMchUp)) {
          registry.get<TH1>(HIST("global-matches/TrackRabs"))->Fill(fgValuesMCH.rabs);
        }
      }
      if (IsGoodMuon(fgValuesMCH, fgValuesMCHpv, configMuons.fTrackChi2MchUp, configMuons.fPMchLow, configMuons.fPtMchLow, configMFTs.fEtaMftLow, configMFTs.fEtaMftUp, configMuons.fRabsLow, configMuons.fRabsUp, 1.E10)) {
        if (IsGoodGlobalMatching(fgValuesMFT, configMFTs.fTrackChi2MftUp, configMFTs.fTrackNClustMftLow, fMatchingChi2MftMchUp)) {
          registry.get<TH1>(HIST("global-matches/TrackPDCA"))->Fill(fgValuesMCH.pDca);
        }
      }
      if (IsGoodMuon(fgValuesMCH, fgValuesMCHpv, configMuons.fTrackChi2MchUp, configMuons.fPMchLow, configMuons.fPtMchLow, configMFTs.fEtaMftLow, configMFTs.fEtaMftUp, configMuons.fRabsLow, configMuons.fRabsUp, configMuons.fSigmaPdcaUp)) {
        if (IsGoodGlobalMatching(fgValuesMFT, configMFTs.fTrackChi2MftUp, configMFTs.fTrackNClustMftLow, fMatchingChi2MftMchUp)) {
          registry.get<TH1>(HIST("global-matches/TrackPhi"))->Fill(fgValuesMCH.phi * 180.0 / TMath::Pi());
          registry.get<TH1>(HIST("global-matches/TrackDCA"))->Fill(std::sqrt(fgValuesMCHpv.dcaX * fgValuesMCHpv.dcaX + fgValuesMCHpv.dcaY * fgValuesMCHpv.dcaY));

          registry.get<TH1>(HIST("global-matches/TrackP_glo"))->Fill(fgValuesGlobal.p);
          registry.get<TH1>(HIST("global-matches/TrackPt_glo"))->Fill(fgValuesGlobal.pT);
          registry.get<TH1>(HIST("global-matches/TrackEta_glo"))->Fill(fgValuesGlobal.eta);
          registry.get<TH1>(HIST("global-matches/TrackPhi_glo"))->Fill(fgValuesGlobal.phi * 180.0 / TMath::Pi());
          registry.get<TH1>(HIST("global-matches/TrackDCA_glo"))->Fill(std::sqrt(fgValuesGlobal.dcaX * fgValuesGlobal.dcaX + fgValuesGlobal.dcaY * fgValuesGlobal.dcaY));
        }
        if (IsGoodGlobalMatching(fgValuesMFT, 1.E10, configMFTs.fTrackNClustMftLow, fMatchingChi2MftMchUp)) {
          registry.get<TH1>(HIST("global-matches/TrackChi2_MFT"))->Fill(fgValuesMFT.chi2);
        }
        if (IsGoodGlobalMatching(fgValuesMFT, configMFTs.fTrackChi2MftUp, 0, fMatchingChi2MftMchUp)) {
          registry.get<TH1>(HIST("global-matches/TrackNclusters_MFT"))->Fill(fgValuesMFT.nClusters);
        }
        if (IsGoodGlobalMatching(fgValuesMFT, configMFTs.fTrackChi2MftUp, configMFTs.fTrackNClustMftLow, 1.E10)) {
          registry.get<TH1>(HIST("global-matches/MatchChi2"))->Fill(fgValuesMFT.chi2matching);
        }
      }
    }
  }

  template <typename VarVector>
  void FillTrackResidualHistograms(VarVector const& fgVectorsMCH, VarVector const& fgVectorsMFT, int quadrant, bool same, bool mixed)
  {
    std::vector<std::array<double, 2>> xPos;
    std::vector<std::array<double, 2>> yPos;
    std::vector<std::array<double, 2>> thetax;
    std::vector<std::array<double, 2>> thetay;
    for (int zi = 0; zi < static_cast<int>(zRefPlane.size()); zi++) {
      xPos.emplace_back(std::array<double, 2>{fgVectorsMCH[zi].x, fgVectorsMFT[zi].x});
      yPos.emplace_back(std::array<double, 2>{fgVectorsMCH[zi].y, fgVectorsMFT[zi].y});
      thetax.emplace_back(std::array<double, 2>{
        std::atan2(fgVectorsMCH[zi].px, -1.0 * fgVectorsMCH[zi].pz) * 180 / TMath::Pi(),
        std::atan2(fgVectorsMFT[zi].px, -1.0 * fgVectorsMFT[zi].pz) * 180 / TMath::Pi()});
      thetay.emplace_back(std::array<double, 2>{
        std::atan2(fgVectorsMCH[zi].py, -1.0 * fgVectorsMCH[zi].pz) * 180 / TMath::Pi(),
        std::atan2(fgVectorsMFT[zi].py, -1.0 * fgVectorsMFT[zi].pz) * 180 / TMath::Pi()});
    }

    for (int i = 0; i < static_cast<int>(zRefPlane.size()); i++) {
      if (same) {
        std::get<std::shared_ptr<TH2>>(trackResidualsHistos[i][quadrant]["dx_vs_x"])->Fill(std::fabs(xPos[i][1]), xPos[i][0] - xPos[i][1]);
        std::get<std::shared_ptr<TH2>>(trackResidualsHistos[i][quadrant]["dx_vs_y"])->Fill(std::fabs(yPos[i][1]), xPos[i][0] - xPos[i][1]);
        std::get<std::shared_ptr<TH2>>(trackResidualsHistos[i][quadrant]["dy_vs_x"])->Fill(std::fabs(xPos[i][1]), yPos[i][0] - yPos[i][1]);
        std::get<std::shared_ptr<TH2>>(trackResidualsHistos[i][quadrant]["dy_vs_y"])->Fill(std::fabs(yPos[i][1]), yPos[i][0] - yPos[i][1]);

        std::get<std::shared_ptr<TH2>>(trackResidualsHistos[i][quadrant]["dthetax_vs_x"])->Fill(std::fabs(xPos[i][1]), thetax[i][0] - thetax[i][1]);
        std::get<std::shared_ptr<TH2>>(trackResidualsHistos[i][quadrant]["dthetax_vs_y"])->Fill(std::fabs(yPos[i][1]), thetax[i][0] - thetax[i][1]);
        std::get<std::shared_ptr<TH2>>(trackResidualsHistos[i][quadrant]["dthetax_vs_thetax"])->Fill(std::fabs(thetax[i][1]), thetax[i][0] - thetax[i][1]);
        std::get<std::shared_ptr<TH2>>(trackResidualsHistos[i][quadrant]["dthetay_vs_x"])->Fill(std::fabs(xPos[i][1]), thetay[i][0] - thetay[i][1]);
        std::get<std::shared_ptr<TH2>>(trackResidualsHistos[i][quadrant]["dthetay_vs_y"])->Fill(std::fabs(yPos[i][1]), thetay[i][0] - thetay[i][1]);
        std::get<std::shared_ptr<TH2>>(trackResidualsHistos[i][quadrant]["dthetay_vs_thetay"])->Fill(std::fabs(thetay[i][1]), thetay[i][0] - thetay[i][1]);
      }
      if (mixed) {
        std::get<std::shared_ptr<TH2>>(trackResidualsHistosMixedEvents[i][quadrant]["dx_vs_x"])->Fill(std::fabs(xPos[i][1]), xPos[i][0] - xPos[i][1]);
        std::get<std::shared_ptr<TH2>>(trackResidualsHistosMixedEvents[i][quadrant]["dx_vs_y"])->Fill(std::fabs(yPos[i][1]), xPos[i][0] - xPos[i][1]);
        std::get<std::shared_ptr<TH2>>(trackResidualsHistosMixedEvents[i][quadrant]["dy_vs_x"])->Fill(std::fabs(xPos[i][1]), yPos[i][0] - yPos[i][1]);
        std::get<std::shared_ptr<TH2>>(trackResidualsHistosMixedEvents[i][quadrant]["dy_vs_y"])->Fill(std::fabs(yPos[i][1]), yPos[i][0] - yPos[i][1]);

        std::get<std::shared_ptr<TH2>>(trackResidualsHistosMixedEvents[i][quadrant]["dthetax_vs_x"])->Fill(std::fabs(xPos[i][1]), thetax[i][0] - thetax[i][1]);
        std::get<std::shared_ptr<TH2>>(trackResidualsHistosMixedEvents[i][quadrant]["dthetax_vs_y"])->Fill(std::fabs(yPos[i][1]), thetax[i][0] - thetax[i][1]);
        std::get<std::shared_ptr<TH2>>(trackResidualsHistosMixedEvents[i][quadrant]["dthetax_vs_thetax"])->Fill(std::fabs(thetax[i][1]), thetax[i][0] - thetax[i][1]);
        std::get<std::shared_ptr<TH2>>(trackResidualsHistosMixedEvents[i][quadrant]["dthetay_vs_x"])->Fill(std::fabs(xPos[i][1]), thetay[i][0] - thetay[i][1]);
        std::get<std::shared_ptr<TH2>>(trackResidualsHistosMixedEvents[i][quadrant]["dthetay_vs_y"])->Fill(std::fabs(yPos[i][1]), thetay[i][0] - thetay[i][1]);
        std::get<std::shared_ptr<TH2>>(trackResidualsHistosMixedEvents[i][quadrant]["dthetay_vs_thetay"])->Fill(std::fabs(thetay[i][1]), thetay[i][0] - thetay[i][1]);
      }
    }
  }

  template <bool MuonFillMap, bool GlobalMuonFillMap, bool GlobalMatchingFillMap, typename VarT, typename VarC>
  void FillDCAHistograms(VarT const& fgValuesMCH, VarT const& fgValuesMCHpv, VarT const& fgValuesMFT, VarC const& fgValuesColl, int sign, int quadrant, bool same, bool mixed)
  {
    if constexpr (static_cast<bool>(MuonFillMap)) {
      if (same) {
        std::get<std::shared_ptr<TH1>>(dcaHistos[1][quadrant][0]["DCA_x"])->Fill(fgValuesMCH.dcaX);
        std::get<std::shared_ptr<TH1>>(dcaHistos[1][quadrant][0]["DCA_y"])->Fill(fgValuesMCH.dcaY);
        std::get<std::shared_ptr<TH1>>(dcaHistos[1][quadrant][0]["DCA_xy"])->Fill(fgValuesMCH.dcaXY);
        std::get<std::shared_ptr<TH1>>(dcaHistos[1][quadrant][sign]["DCA_x"])->Fill(fgValuesMCH.dcaX);
        std::get<std::shared_ptr<TH1>>(dcaHistos[1][quadrant][sign]["DCA_y"])->Fill(fgValuesMCH.dcaY);
        std::get<std::shared_ptr<TH1>>(dcaHistos[1][quadrant][sign]["DCA_xy"])->Fill(fgValuesMCH.dcaXY);
        std::get<std::shared_ptr<TH2>>(dcaHistos[1][quadrant][0]["DCA_y_vs_DCA_x"])->Fill(fgValuesMCH.dcaX, fgValuesMCH.dcaY);
        // TODO: doesn't work, could check, but not very important
        // std::get<std::shared_ptr<TH2>>(dcaHistos[1][quadrant][0]["DCA_x_vs_z"])->Fill(fgValuesColl.z, fgValuesMCH.dcaX);
        // std::get<std::shared_ptr<TH2>>(dcaHistos[1][quadrant][0]["DCA_y_vs_z"])->Fill(fgValuesColl.z, fgValuesMCH.dcaY);
        std::get<std::shared_ptr<TH2>>(dcaHistos[1][quadrant][0]["DCA_x_vs_phi"])->Fill(fgValuesMCH.phi, fgValuesMCH.dcaX);
        std::get<std::shared_ptr<TH2>>(dcaHistos[1][quadrant][0]["DCA_y_vs_phi"])->Fill(fgValuesMCH.phi, fgValuesMCH.dcaY);
        std::get<std::shared_ptr<TH2>>(dcaHistos[1][quadrant][0]["DCA_xy_vs_phi"])->Fill(fgValuesMCH.phi, fgValuesMCH.dcaXY);
      }

      if (mixed) {
        std::get<std::shared_ptr<TH1>>(dcaHistosMixedEvents[1][quadrant][0]["DCA_x"])->Fill(fgValuesMCH.dcaX);
        std::get<std::shared_ptr<TH1>>(dcaHistosMixedEvents[1][quadrant][0]["DCA_y"])->Fill(fgValuesMCH.dcaY);
        std::get<std::shared_ptr<TH1>>(dcaHistosMixedEvents[1][quadrant][sign]["DCA_x"])->Fill(fgValuesMCH.dcaX);
        std::get<std::shared_ptr<TH1>>(dcaHistosMixedEvents[1][quadrant][sign]["DCA_y"])->Fill(fgValuesMCH.dcaY);
        // TODO: doesn't work, could check, but not very important
        // std::get<std::shared_ptr<TH2>>(dcaHistosMixedEvents[1][quadrant][0]["DCA_x_vs_z"])->Fill(fgValuesColl.z, fgValuesMCH.dcaX);
        // std::get<std::shared_ptr<TH2>>(dcaHistosMixedEvents[1][quadrant][0]["DCA_y_vs_z"])->Fill(fgValuesColl.z, fgValuesMCH.dcaY);
      }
    }

    if constexpr (static_cast<bool>(GlobalMuonFillMap)) {
      if (same) {
        std::get<std::shared_ptr<TH1>>(dcaHistos[0][quadrant][0]["DCA_x"])->Fill(fgValuesMFT.dcaX);
        std::get<std::shared_ptr<TH1>>(dcaHistos[0][quadrant][0]["DCA_y"])->Fill(fgValuesMFT.dcaY);
        std::get<std::shared_ptr<TH1>>(dcaHistos[0][quadrant][0]["DCA_xy"])->Fill(fgValuesMFT.dcaXY);
        std::get<std::shared_ptr<TH1>>(dcaHistos[0][quadrant][sign]["DCA_x"])->Fill(fgValuesMFT.dcaX);
        std::get<std::shared_ptr<TH1>>(dcaHistos[0][quadrant][sign]["DCA_y"])->Fill(fgValuesMFT.dcaY);
        std::get<std::shared_ptr<TH1>>(dcaHistos[0][quadrant][sign]["DCA_xy"])->Fill(fgValuesMFT.dcaXY);
        std::get<std::shared_ptr<TH2>>(dcaHistos[0][quadrant][0]["DCA_y_vs_DCA_x"])->Fill(fgValuesMFT.dcaX, fgValuesMFT.dcaY);
        std::get<std::shared_ptr<TH2>>(dcaHistos[0][quadrant][0]["DCA_x_vs_z"])->Fill(fgValuesColl.z, fgValuesMFT.dcaX);
        std::get<std::shared_ptr<TH2>>(dcaHistos[0][quadrant][0]["DCA_y_vs_z"])->Fill(fgValuesColl.z, fgValuesMFT.dcaY);
        std::get<std::shared_ptr<TH2>>(dcaHistos[0][quadrant][0]["DCA_xy_vs_z"])->Fill(fgValuesColl.z, fgValuesMFT.dcaXY);
        std::get<std::shared_ptr<TH2>>(dcaHistos[0][quadrant][0]["DCA_x_vs_phi"])->Fill(fgValuesMFT.phi, fgValuesMFT.dcaX);
        std::get<std::shared_ptr<TH2>>(dcaHistos[0][quadrant][0]["DCA_y_vs_phi"])->Fill(fgValuesMFT.phi, fgValuesMFT.dcaY);
        std::get<std::shared_ptr<TH2>>(dcaHistos[0][quadrant][0]["DCA_xy_vs_phi"])->Fill(fgValuesMFT.phi, fgValuesMFT.dcaXY);
      }

      if (mixed) {
        std::get<std::shared_ptr<TH1>>(dcaHistosMixedEvents[0][quadrant][0]["DCA_x"])->Fill(fgValuesMFT.dcaX);
        std::get<std::shared_ptr<TH1>>(dcaHistosMixedEvents[0][quadrant][0]["DCA_y"])->Fill(fgValuesMFT.dcaY);
        std::get<std::shared_ptr<TH1>>(dcaHistosMixedEvents[0][quadrant][sign]["DCA_x"])->Fill(fgValuesMFT.dcaX);
        std::get<std::shared_ptr<TH1>>(dcaHistosMixedEvents[0][quadrant][sign]["DCA_y"])->Fill(fgValuesMFT.dcaY);
        std::get<std::shared_ptr<TH2>>(dcaHistosMixedEvents[0][quadrant][0]["DCA_x_vs_z"])->Fill(fgValuesColl.z, fgValuesMFT.dcaX);
        std::get<std::shared_ptr<TH2>>(dcaHistosMixedEvents[0][quadrant][0]["DCA_y_vs_z"])->Fill(fgValuesColl.z, fgValuesMFT.dcaY);
      }
    }

    if constexpr (static_cast<bool>(GlobalMatchingFillMap)) {
      if (same) {
        // In prinicple the MCH DCAs are the same for global and single
        std::get<std::shared_ptr<TH1>>(dcaHistosGlobal[1][quadrant][0]["DCA_x"])->Fill(fgValuesMCH.dcaX);
        std::get<std::shared_ptr<TH1>>(dcaHistosGlobal[1][quadrant][0]["DCA_y"])->Fill(fgValuesMCH.dcaY);
        std::get<std::shared_ptr<TH1>>(dcaHistosGlobal[1][quadrant][0]["DCA_xy"])->Fill(fgValuesMCH.dcaXY);
        std::get<std::shared_ptr<TH1>>(dcaHistosGlobal[1][quadrant][sign]["DCA_x"])->Fill(fgValuesMCH.dcaX);
        std::get<std::shared_ptr<TH1>>(dcaHistosGlobal[1][quadrant][sign]["DCA_y"])->Fill(fgValuesMCH.dcaY);
        std::get<std::shared_ptr<TH1>>(dcaHistosGlobal[1][quadrant][sign]["DCA_xy"])->Fill(fgValuesMCH.dcaXY);
        std::get<std::shared_ptr<TH2>>(dcaHistosGlobal[1][quadrant][0]["DCA_y_vs_DCA_x"])->Fill(fgValuesMCH.dcaX, fgValuesMCH.dcaY);
        std::get<std::shared_ptr<TH2>>(dcaHistosGlobal[1][quadrant][0]["DCA_x_vs_z"])->Fill(fgValuesColl.z, fgValuesMCH.dcaX);
        std::get<std::shared_ptr<TH2>>(dcaHistosGlobal[1][quadrant][0]["DCA_y_vs_z"])->Fill(fgValuesColl.z, fgValuesMCH.dcaY);
        std::get<std::shared_ptr<TH2>>(dcaHistosGlobal[1][quadrant][0]["DCA_xy_vs_z"])->Fill(fgValuesColl.z, fgValuesMCH.dcaXY);
        std::get<std::shared_ptr<TH2>>(dcaHistosGlobal[1][quadrant][0]["DCA_x_vs_phi"])->Fill(fgValuesMCH.phi, fgValuesMCH.dcaX);
        std::get<std::shared_ptr<TH2>>(dcaHistosGlobal[1][quadrant][0]["DCA_y_vs_phi"])->Fill(fgValuesMCH.phi, fgValuesMCH.dcaY);
        std::get<std::shared_ptr<TH2>>(dcaHistosGlobal[1][quadrant][0]["DCA_xy_vs_phi"])->Fill(fgValuesMCH.phi, fgValuesMCH.dcaXY);

        std::get<std::shared_ptr<TH1>>(dcaHistosGlobal[0][quadrant][0]["DCA_x"])->Fill(fgValuesMFT.dcaX);
        std::get<std::shared_ptr<TH1>>(dcaHistosGlobal[0][quadrant][0]["DCA_y"])->Fill(fgValuesMFT.dcaY);
        std::get<std::shared_ptr<TH1>>(dcaHistosGlobal[0][quadrant][0]["DCA_xy"])->Fill(fgValuesMFT.dcaXY);
        std::get<std::shared_ptr<TH1>>(dcaHistosGlobal[0][quadrant][sign]["DCA_x"])->Fill(fgValuesMFT.dcaX);
        std::get<std::shared_ptr<TH1>>(dcaHistosGlobal[0][quadrant][sign]["DCA_y"])->Fill(fgValuesMFT.dcaY);
        std::get<std::shared_ptr<TH1>>(dcaHistosGlobal[0][quadrant][sign]["DCA_xy"])->Fill(fgValuesMFT.dcaXY);
        std::get<std::shared_ptr<TH2>>(dcaHistosGlobal[0][quadrant][0]["DCA_y_vs_DCA_x"])->Fill(fgValuesMFT.dcaX, fgValuesMFT.dcaY);
        std::get<std::shared_ptr<TH2>>(dcaHistosGlobal[0][quadrant][0]["DCA_x_vs_z"])->Fill(fgValuesColl.z, fgValuesMFT.dcaX);
        std::get<std::shared_ptr<TH2>>(dcaHistosGlobal[0][quadrant][0]["DCA_y_vs_z"])->Fill(fgValuesColl.z, fgValuesMFT.dcaY);
        std::get<std::shared_ptr<TH2>>(dcaHistosGlobal[0][quadrant][0]["DCA_xy_vs_z"])->Fill(fgValuesColl.z, fgValuesMFT.dcaXY);
        std::get<std::shared_ptr<TH2>>(dcaHistosGlobal[0][quadrant][0]["DCA_x_vs_phi"])->Fill(fgValuesMFT.phi, fgValuesMFT.dcaX);
        std::get<std::shared_ptr<TH2>>(dcaHistosGlobal[0][quadrant][0]["DCA_y_vs_phi"])->Fill(fgValuesMFT.phi, fgValuesMFT.dcaY);
        std::get<std::shared_ptr<TH2>>(dcaHistosGlobal[0][quadrant][0]["DCA_xy_vs_phi"])->Fill(fgValuesMFT.phi, fgValuesMFT.dcaXY);

        // Calculate the difference (scaled) DCA_MFT - DCA_MCH
        std::get<std::shared_ptr<TH1>>(dcaHistosGlobalSubtracted[quadrant][0]["DCA_x"])->Fill(fgValuesMFT.dcaX - fgValuesMCH.dcaX);
        std::get<std::shared_ptr<TH1>>(dcaHistosGlobalSubtracted[quadrant][0]["DCA_y"])->Fill(fgValuesMFT.dcaY - fgValuesMCH.dcaY);
        std::get<std::shared_ptr<TH1>>(dcaHistosGlobalSubtracted[quadrant][0]["DCA_xy"])->Fill(fgValuesMFT.dcaXY - fgValuesMCH.dcaXY);
        std::get<std::shared_ptr<TH1>>(dcaHistosGlobalSubtracted[quadrant][sign]["DCA_x"])->Fill(fgValuesMFT.dcaX - fgValuesMCH.dcaX);
        std::get<std::shared_ptr<TH1>>(dcaHistosGlobalSubtracted[quadrant][sign]["DCA_y"])->Fill(fgValuesMFT.dcaY - fgValuesMCH.dcaY);
        std::get<std::shared_ptr<TH1>>(dcaHistosGlobalSubtracted[quadrant][sign]["DCA_xy"])->Fill(fgValuesMFT.dcaXY - fgValuesMCH.dcaXY);
        std::get<std::shared_ptr<TH2>>(dcaHistosGlobalSubtracted[quadrant][0]["DCA_y_vs_DCA_x"])->Fill(fgValuesMFT.dcaX, fgValuesMFT.dcaY);
        std::get<std::shared_ptr<TH2>>(dcaHistosGlobalSubtracted[quadrant][0]["DCA_x_vs_z"])->Fill(fgValuesColl.z, fgValuesMFT.dcaX - fgValuesMCH.dcaX);
        std::get<std::shared_ptr<TH2>>(dcaHistosGlobalSubtracted[quadrant][0]["DCA_y_vs_z"])->Fill(fgValuesColl.z, fgValuesMFT.dcaY - fgValuesMCH.dcaY);
        std::get<std::shared_ptr<TH2>>(dcaHistosGlobalSubtracted[quadrant][0]["DCA_xy_vs_z"])->Fill(fgValuesColl.z, fgValuesMFT.dcaXY - fgValuesMCH.dcaXY);
        // Which phi to use for the subtracted histos (from MCH or from MFT? - I guess MFT has better resolution)
        std::get<std::shared_ptr<TH2>>(dcaHistosGlobalSubtracted[quadrant][0]["DCA_x_vs_phi"])->Fill(fgValuesMFT.phi, fgValuesMFT.dcaX - fgValuesMCH.dcaX);
        std::get<std::shared_ptr<TH2>>(dcaHistosGlobalSubtracted[quadrant][0]["DCA_y_vs_phi"])->Fill(fgValuesMFT.phi, fgValuesMFT.dcaY - fgValuesMCH.dcaY);
        std::get<std::shared_ptr<TH2>>(dcaHistosGlobalSubtracted[quadrant][0]["DCA_xy_vs_phi"])->Fill(fgValuesMFT.phi, fgValuesMFT.dcaXY - fgValuesMCH.dcaXY);

        // Fill only when MCH momentum is larger than cut
        if (fgValuesMCHpv.p > fMchPUpForGlobalDCA) {
          std::get<std::shared_ptr<TH1>>(dcaHistosGlobalSubtractedMCHpCut[quadrant][0]["DCA_x"])->Fill(fgValuesMFT.dcaX - fgValuesMCH.dcaX);
          std::get<std::shared_ptr<TH1>>(dcaHistosGlobalSubtractedMCHpCut[quadrant][0]["DCA_y"])->Fill(fgValuesMFT.dcaY - fgValuesMCH.dcaY);
          std::get<std::shared_ptr<TH1>>(dcaHistosGlobalSubtractedMCHpCut[quadrant][0]["DCA_xy"])->Fill(fgValuesMFT.dcaXY - fgValuesMCH.dcaXY);
          std::get<std::shared_ptr<TH1>>(dcaHistosGlobalSubtractedMCHpCut[quadrant][sign]["DCA_x"])->Fill(fgValuesMFT.dcaX - fgValuesMCH.dcaX);
          std::get<std::shared_ptr<TH1>>(dcaHistosGlobalSubtractedMCHpCut[quadrant][sign]["DCA_y"])->Fill(fgValuesMFT.dcaY - fgValuesMCH.dcaY);
          std::get<std::shared_ptr<TH1>>(dcaHistosGlobalSubtractedMCHpCut[quadrant][sign]["DCA_xy"])->Fill(fgValuesMFT.dcaXY - fgValuesMCH.dcaXY);
          std::get<std::shared_ptr<TH2>>(dcaHistosGlobalSubtractedMCHpCut[quadrant][0]["DCA_y_vs_DCA_x"])->Fill(fgValuesMFT.dcaX, fgValuesMFT.dcaY);
          std::get<std::shared_ptr<TH2>>(dcaHistosGlobalSubtractedMCHpCut[quadrant][0]["DCA_x_vs_z"])->Fill(fgValuesColl.z, fgValuesMFT.dcaX - fgValuesMCH.dcaX);
          std::get<std::shared_ptr<TH2>>(dcaHistosGlobalSubtractedMCHpCut[quadrant][0]["DCA_y_vs_z"])->Fill(fgValuesColl.z, fgValuesMFT.dcaY - fgValuesMCH.dcaY);
          std::get<std::shared_ptr<TH2>>(dcaHistosGlobalSubtractedMCHpCut[quadrant][0]["DCA_xy_vs_z"])->Fill(fgValuesColl.z, fgValuesMFT.dcaXY - fgValuesMCH.dcaXY);
          // Which phi to use for the subtracted histos (from MCH or from MFT? - I guess MFT has better resolution)
          std::get<std::shared_ptr<TH2>>(dcaHistosGlobalSubtractedMCHpCut[quadrant][0]["DCA_x_vs_phi"])->Fill(fgValuesMFT.phi, fgValuesMFT.dcaX - fgValuesMCH.dcaX);
          std::get<std::shared_ptr<TH2>>(dcaHistosGlobalSubtractedMCHpCut[quadrant][0]["DCA_y_vs_phi"])->Fill(fgValuesMFT.phi, fgValuesMFT.dcaY - fgValuesMCH.dcaY);
          std::get<std::shared_ptr<TH2>>(dcaHistosGlobalSubtractedMCHpCut[quadrant][0]["DCA_xy_vs_phi"])->Fill(fgValuesMFT.phi, fgValuesMFT.dcaXY - fgValuesMCH.dcaXY);
        }
      }
    }
  }

  template <bool MuonFillMap, bool MFTFillMap, typename Var>
  void FillResidualHistograms(Var const& fgValuesProp, Var const& fgValuesMCH, Var const& fgValuesMCHpv, Var const& fgValuesMFT, float xCls, float yCls, int topBottom, int posNeg, int quadrant, int chamber, int deIndex, bool same, bool mixed)
  {
    std::array<double, 2> xPos{xCls, fgValuesProp.x};
    std::array<double, 2> yPos{yCls, fgValuesProp.y};
    double phiClus = std::atan2(yCls, xCls) * 180 / TMath::Pi();

    if constexpr (static_cast<bool>(MuonFillMap)) {
      if (IsGoodMuon(fgValuesMCH, fgValuesMCHpv, configMuons.fTrackChi2MchUp, 20., configMuons.fPtMchLow, configMFTs.fEtaMftLow, configMFTs.fEtaMftUp, configMuons.fRabsLow, configMuons.fRabsUp, configMuons.fSigmaPdcaUp)) {
        if (same) {
          std::get<std::shared_ptr<TH2>>(mchResidualsHistosPerDE[topBottom][posNeg][chamber]["dx_vs_de"])->Fill(deIndex, xPos[0] - xPos[1]);
          std::get<std::shared_ptr<TH2>>(mchResidualsHistosPerDE[topBottom][posNeg][chamber]["dy_vs_de"])->Fill(deIndex, yPos[0] - yPos[1]);
        }
        if (mixed) {
          std::get<std::shared_ptr<TH2>>(mchResidualsHistosPerDEMixedEvents[topBottom][posNeg][chamber]["dx_vs_de"])->Fill(deIndex, xPos[0] - xPos[1]);
          std::get<std::shared_ptr<TH2>>(mchResidualsHistosPerDEMixedEvents[topBottom][posNeg][chamber]["dy_vs_de"])->Fill(deIndex, yPos[0] - yPos[1]);
        }
      }
    }

    if constexpr (static_cast<bool>(MFTFillMap)) {
      if (IsGoodMuon(fgValuesMCH, fgValuesMCHpv, configMuons.fTrackChi2MchUp, 20., configMuons.fPtMchLow, configMFTs.fEtaMftLow, configMFTs.fEtaMftUp, configMuons.fRabsLow, configMuons.fRabsUp, configMuons.fSigmaPdcaUp)) {
        if (IsGoodMFT(fgValuesMFT, configMFTs.fTrackChi2MftUp, configMFTs.fTrackNClustMftLow)) {
          if (same) {
            std::get<std::shared_ptr<TH2>>(residualsHistos[quadrant][chamber]["dx_vs_x"])->Fill(std::fabs(xPos[1]), xPos[0] - xPos[1]);
            std::get<std::shared_ptr<TH2>>(residualsHistos[quadrant][chamber]["dx_vs_y"])->Fill(std::fabs(yPos[1]), xPos[0] - xPos[1]);
            std::get<std::shared_ptr<TH2>>(residualsHistos[quadrant][chamber]["dy_vs_x"])->Fill(std::fabs(xPos[1]), yPos[0] - yPos[1]);
            std::get<std::shared_ptr<TH2>>(residualsHistos[quadrant][chamber]["dy_vs_y"])->Fill(std::fabs(yPos[1]), yPos[0] - yPos[1]);

            // residuals vs. DE index
            std::get<std::shared_ptr<TH2>>(residualsHistosPerDE[topBottom][posNeg][chamber]["dx_vs_de"])->Fill(deIndex, xPos[0] - xPos[1]);
            std::get<std::shared_ptr<TH2>>(residualsHistosPerDE[topBottom][posNeg][chamber]["dy_vs_de"])->Fill(deIndex, yPos[0] - yPos[1]);

            // residuals vs. cluster phi
            std::get<std::shared_ptr<TH2>>(residualsHistosPerDE[topBottom][posNeg][chamber]["dx_vs_phi"])->Fill(phiClus, xPos[0] - xPos[1]);
            std::get<std::shared_ptr<TH2>>(residualsHistosPerDE[topBottom][posNeg][chamber]["dy_vs_phi"])->Fill(phiClus, yPos[0] - yPos[1]);
          }
          if (mixed) {
            std::get<std::shared_ptr<TH2>>(residualsHistosMixedEvents[quadrant][chamber]["dx_vs_x"])->Fill(std::fabs(xPos[1]), xPos[0] - xPos[1]);
            std::get<std::shared_ptr<TH2>>(residualsHistosMixedEvents[quadrant][chamber]["dx_vs_y"])->Fill(std::fabs(yPos[1]), xPos[0] - xPos[1]);
            std::get<std::shared_ptr<TH2>>(residualsHistosMixedEvents[quadrant][chamber]["dy_vs_x"])->Fill(std::fabs(xPos[1]), yPos[0] - yPos[1]);
            std::get<std::shared_ptr<TH2>>(residualsHistosMixedEvents[quadrant][chamber]["dy_vs_y"])->Fill(std::fabs(yPos[1]), yPos[0] - yPos[1]);

            // residuals vs. DE index
            std::get<std::shared_ptr<TH2>>(residualsHistosPerDEMixedEvents[topBottom][posNeg][chamber]["dx_vs_de"])->Fill(deIndex, xPos[0] - xPos[1]);
            std::get<std::shared_ptr<TH2>>(residualsHistosPerDEMixedEvents[topBottom][posNeg][chamber]["dy_vs_de"])->Fill(deIndex, yPos[0] - yPos[1]);

            // residuals vs. cluster phi
            std::get<std::shared_ptr<TH2>>(residualsHistosPerDEMixedEvents[topBottom][posNeg][chamber]["dx_vs_phi"])->Fill(phiClus, xPos[0] - xPos[1]);
            std::get<std::shared_ptr<TH2>>(residualsHistosPerDEMixedEvents[topBottom][posNeg][chamber]["dy_vs_phi"])->Fill(phiClus, yPos[0] - yPos[1]);
          }
        }
      }
    }
  }

  template <typename Var>
  void resetVar(Var& fgValues)
  {
    fgValues = {};
  }

  void initCCDB(aod::BCsWithTimestamps const& bcs)
  {
    // Update CCDB informations
    if (bcs.size() > 0 && fCurrentRun != bcs.begin().runNumber()) {
      // Load magnetic field information from CCDB/local
      ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
      grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(configCCDB.grpmagPath, bcs.begin().timestamp());
      if (grpmag != nullptr) {
        base::Propagator::initFieldFromGRP(grpmag);
        TrackExtrap::setField();
        TrackExtrap::useExtrapV2();
        fieldB = static_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField()); // for MFT
        double centerMFT[3] = {0, 0, -61.4};                                                         // or use middle point between Vtx and MFT?
        Bz = fieldB->getBz(centerMFT);                                                               // Get field at centre of MFT
      } else {
        LOGF(fatal, "GRP object is not available in CCDB at timestamp=%llu", bcs.begin().timestamp());
      }

      // Load geometry information from CCDB/local
      LOGF(info, "Loading reference aligned geometry from CCDB no later than %d", configCCDB.nolaterthan.value);
      ccdb->setCreatedNotAfter(configCCDB.nolaterthan); // this timestamp has to be consistent with what has been used in reco
      geoRef = ccdb->getForTimeStamp<TGeoManager>(configCCDB.geoPath, bcs.begin().timestamp());
      ccdb->clearCache(configCCDB.geoPath);
      if (geoRef != nullptr) {
        transformation = geo::transformationFromTGeoManager(*geoRef);
      } else {
        LOGF(fatal, "Reference aligned geometry object is not available in CCDB at timestamp=%llu", bcs.begin().timestamp());
      }
      for (int i = 0; i < 156; i++) {
        int iDEN = GetDetElemId(i);
        transformRef[iDEN] = transformation(iDEN);
      }

      if (configRealign.fDoRealign) {
        LOGF(info, "Loading new aligned geometry from CCDB no later than %d", configCCDB.nolaterthanRealign.value);
        ccdb->setCreatedNotAfter(configCCDB.nolaterthanRealign); // make sure this timestamp can be resolved regarding the reference one
        geoNew = ccdb->getForTimeStamp<TGeoManager>(configCCDB.geoPathRealign, bcs.begin().timestamp());
        ccdb->clearCache(configCCDB.geoPathRealign);
        if (geoNew != nullptr) {
          transformation = geo::transformationFromTGeoManager(*geoNew);
        } else {
          LOGF(fatal, "New aligned geometry object is not available in CCDB at timestamp=%llu", bcs.begin().timestamp());
        }
        for (int i = 0; i < 156; i++) {
          int iDEN = GetDetElemId(i);
          transformNew[iDEN] = transformation(iDEN);
        }
      }

      fCurrentRun = bcs.begin().runNumber();
    }
  }

  void init(InitContext const&)
  {
    fCurrentRun = 0;

    // Configuration for CCDB server
    ccdb->setURL(configCCDB.ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    // Configuration for track fitter
    const auto& trackerParam = TrackerParam::Instance();
    trackFitter.setBendingVertexDispersion(trackerParam.bendingVertexDispersion);
    trackFitter.setChamberResolution(configRealign.fChamberResolutionX, configRealign.fChamberResolutionY);
    trackFitter.smoothTracks(true);
    trackFitter.useChamberResolution();
    mImproveCutChi2 = 2. * configRealign.fSigmaCutImprove * configRealign.fSigmaCutImprove;

    CreateBasicHistograms();
    CreateDetailedHistograms();
  }

  template <bool MuonFillMap, bool GlobalMuonFillMap, bool GlobalMatchingFillMap,
            typename TEventMap, typename TMFTTracks, typename TMFTTrack, typename TTrack,
            typename TMCHTrack, typename VarC, typename VarT>
  void runDCA(TEventMap const& collisions, TMFTTracks const& mfts, TMFTTrack const& mfttrack, TTrack const& muon, TMCHTrack const& mchtrack, mch::Track const& mchrealigned, VarC& fgValuesColl, VarT& fgValuesMCH, VarT& fgValuesMCHpv, VarT& fgValuesMFT)
  {
    if constexpr (static_cast<bool>(MuonFillMap)) {

      // track selection
      if (!IsGoodMuon(fgValuesMCH, fgValuesMCHpv, configMuons.fTrackChi2MchUp, 30., 4., configMFTs.fEtaMftLow, configMFTs.fEtaMftUp, configMuons.fRabsLow, configMuons.fRabsUp, configMuons.fSigmaPdcaUp)) {
        return;
      }

      // Loop over collisions
      for (auto& [collisionId, fgValuesColltmp] : collisions) {

        bool sameEvent = (fgValuesColltmp.bc == fgValuesColl.bc);
        bool mixedEvent = IsMixedEvent(fgValuesColltmp, fgValuesColl);
        if (!sameEvent && !mixedEvent) {
          continue;
        }

        // Fill propagation of MCH track to DCA
        if (configRealign.fDoRealign) {
          FillPropagation(mchrealigned, fgValuesColltmp, fgValuesMCH, kToDCA);
        } else {
          FillPropagation<1>(muon, fgValuesColltmp, VarTrack{}, fgValuesMCH, kToDCA);
        }

        double phi = fgValuesMCH.phi * 180 / TMath::Pi();
        int quadrant = GetQuadrantPhi(phi);
        int sign = (fgValuesMCH.sign > 0) ? 1 : 2;

        // Fill DCA QA histograms
        FillDCAHistograms<1, 0, 0>(fgValuesMCH, fgValuesMCHpv, VarTrack{}, fgValuesColltmp, sign, quadrant, sameEvent, mixedEvent);
      }
    }

    if constexpr (static_cast<bool>(GlobalMuonFillMap)) {

      auto mftsThisCollision = mfts.sliceBy(mftPerCollision, fgValuesColl.globalIndex);
      for (auto const& mft : mftsThisCollision) {

        // Fill MFT track
        VarTrack fgValuesMFTtmp;
        FillTrack<0>(mft, fgValuesMFTtmp);

        if (fgValuesMFT.trackTime != fgValuesMFTtmp.trackTime) {
          continue; // if not time compatible
        }

        if (!IsGoodMFT(fgValuesMFTtmp, configMFTs.fTrackChi2MftUp, configMFTs.fTrackNClustMftLow)) {
          continue;
        }

        int quadrant = GetQuadrantTrack(mft);
        if (quadrant < 0) {
          continue;
        }

        int sign = (fgValuesMFTtmp.sign > 0) ? 1 : 2;

        for (auto& [collisionId, fgValuesColltmp] : collisions) {

          bool sameEvent = (fgValuesColltmp.bc == fgValuesColl.bc);
          bool mixedEvent = IsMixedEvent(fgValuesColltmp, fgValuesColl);

          if (!sameEvent && !mixedEvent) {
            continue;
          }

          // Propagate MFT track to DCA
          FillPropagation<0, 1>(mft, fgValuesColltmp, VarTrack{}, fgValuesMFTtmp, kToDCA);

          // Fill DCA QA histograms
          FillDCAHistograms<0, 1, 0>(VarTrack{}, VarTrack{}, fgValuesMFTtmp, fgValuesColltmp, sign, quadrant, sameEvent, mixedEvent);
        }
        resetVar(fgValuesMFTtmp);
      }
    }

    if constexpr (static_cast<bool>(GlobalMatchingFillMap)) {

      // track selection
      if (!IsGoodMuon(fgValuesMCH, fgValuesMCHpv, configMuons.fTrackChi2MchUp, 30., 4., configMFTs.fEtaMftLow, configMFTs.fEtaMftUp, configMuons.fRabsLow, configMuons.fRabsUp, configMuons.fSigmaPdcaUp)) {
        return;
      }
      if (!IsGoodMFT(fgValuesMFT, configMFTs.fTrackChi2MftUp, configMFTs.fTrackNClustMftLow)) {
        return;
      }

      bool goodGlobalMuonTracks = IsGoodGlobalMuon(fgValuesMCH, fgValuesMCHpv);
      bool goodGlobalMuonMatches = IsGoodGlobalMatching(fgValuesMFT);

      if (!goodGlobalMuonTracks || !goodGlobalMuonMatches) {
        return;
      }

      // do collision loop and (again) propagation?
      // propagation in principle already done in runMuonQA?
      for (auto& [collisionId, fgValuesCollGlo] : collisions) {

        bool sameEvent = (fgValuesCollGlo.bc == fgValuesColl.bc);
        bool mixedEvent = IsMixedEvent(fgValuesCollGlo, fgValuesColl);

        if (!sameEvent && !mixedEvent) {
          continue;
        }

        if (mixedEvent) {
          // TODO: possible to implement mixed events in global DCAs, but not foreseen to be useful in the (near) future
          continue;
        }

        // Fill propagation of MCH track to DCA
        if (configRealign.fDoRealign) {
          FillPropagation(mchrealigned, fgValuesCollGlo, fgValuesMCH, kToDCA);
        } else {
          FillPropagation<1>(mchtrack, fgValuesCollGlo, VarTrack{}, fgValuesMCH, kToDCA);
        }

        // Propagate MFT track to DCA
        FillPropagation<0, 1>(mfttrack, fgValuesCollGlo, VarTrack{}, fgValuesMFT, kToDCA);

        double phi = fgValuesMCH.phi * 180 / TMath::Pi();
        int quadrant = GetQuadrantPhi(phi);
        int sign = (fgValuesMCH.sign > 0) ? 1 : 2;

        // Fill DCA QA histograms
        FillDCAHistograms<0, 0, 1>(fgValuesMCH, fgValuesMCHpv, fgValuesMFT, fgValuesCollGlo, sign, quadrant, true, false);
      }
    }
  }

  template <typename TEventMap, typename TMuons, typename TMFTTracks, typename TMCHTrack, typename TMFTTrack, typename TMuonCls, typename VarC, typename VarT>
  void runResidual(TEventMap const& collisions, TMuons const& muons, TMFTTracks const& mfts, TMuonCls const& clusters, TMCHTrack const& mchtrack, mch::Track const& mchrealigned, TMFTTrack const& mfttrack, VarC const& fgValuesColl, VarT const& fgValuesMCH, VarT const& fgValuesMCHpv, VarT const& fgValuesMFT)
  {
    if (!IsGoodMuon(fgValuesMCH, fgValuesMCHpv, configMuons.fTrackChi2MchUp, 20., configMuons.fPtMchLow, configMFTs.fEtaMftLow, configMFTs.fEtaMftUp, configMuons.fRabsLow, configMuons.fRabsUp, configMuons.fSigmaPdcaUp)) {
      return;
    }

    double phi = fgValuesMCH.phi * 180 / TMath::Pi();
    int quadrant = GetQuadrantPhi(phi);

    //// MCH-MFT track residuals
    if (mfttrack.has_collision()) {
      auto& fgValuesCollMatched = collisions.at(mfttrack.collisionId());

      // Do extrapolation for muons to all reference planes
      std::vector<VarTrack> mchTrackExtrap;
      for (const double z : zRefPlane) {
        VarTrack fgValues;
        if (configRealign.fDoRealign) {
          FillPropagation(mchrealigned, VarColl{}, fgValues, kToZ, z);
        } else {
          FillPropagation<1>(mchtrack, VarColl{}, VarTrack{}, fgValues, kToZ, z);
        }
        mchTrackExtrap.emplace_back(fgValues);
      }

      // Loop over MFT tracks
      for (auto const& mft : mfts) {

        if (!mft.has_collision()) {
          continue;
        }

        // Fill MFT track
        VarTrack fgValuesMFTtmp;
        FillTrack<0>(mft, fgValuesMFTtmp);

        // Track selection
        if (!IsGoodMFT(fgValuesMFTtmp, configMFTs.fTrackChi2MftUp, configMFTs.fTrackNClustMftLow)) {
          continue;
        }

        auto& fgValuesCollMFT = collisions.at(mft.collisionId());

        bool sameEvent = (fgValuesCollMFT.bc == fgValuesCollMatched.bc);
        bool mixedEvent = IsMixedEvent(fgValuesCollMFT, fgValuesColl);
        if (!sameEvent && !mixedEvent) {
          continue;
        }

        // Do extrapolation for MFTs to all reference planes
        std::vector<VarTrack> mftTrackExtrap;
        for (const double z : zRefPlane) {
          VarTrack fgValues;
          FillPropagation<0, 1>(mft, fgValuesCollMFT, mchTrackExtrap[1], fgValues, kToZ, z);
          mftTrackExtrap.emplace_back(fgValues);
        }

        //// Fill QA histograms for alignment checks
        FillTrackResidualHistograms(mchTrackExtrap, mftTrackExtrap, quadrant, sameEvent, mixedEvent);

        resetVar(fgValuesMFTtmp);
      }
    }

    //// Track-Cluster residuals
    for (auto const& muon : muons) {
      if (static_cast<int>(muon.trackType()) <= 2) {
        continue;
      }
      if (!muon.has_collision()) {
        continue;
      }
      auto& fgValuesCollMCH = collisions.at(muon.collisionId());

      bool sameEvent = (fgValuesCollMCH.bc == fgValuesColl.bc);
      bool mixedEvent = IsMixedEvent(fgValuesCollMCH, fgValuesColl);
      if (!sameEvent && !mixedEvent) {
        continue;
      }

      //// Fill MCH clusters: do re-alignment if asked
      TrackRealigned mchrealignedTmp;
      VarClusters fgValuesClsTmp;
      if (!FillClusters(muon, clusters, fgValuesClsTmp, mchrealignedTmp)) {
        continue; // Refit is not valid
      }

      // Loop over attached clusters
      for (int iCls = 0; iCls < static_cast<int>(fgValuesClsTmp.posClusters.size()); iCls++) {

        double phiCls = std::atan2(fgValuesClsTmp.posClusters[iCls][1], fgValuesClsTmp.posClusters[iCls][0]) * 180 / TMath::Pi();
        int quadrantCls = GetQuadrantPhi(phiCls);
        int DEId = fgValuesClsTmp.DEIDs[iCls];
        int chamber = DEId / 100 - 1;
        int deIndex = DEId % 100;

        //// MCH residuals
        //// Propagate MCH track to given cluster
        VarTrack fgValuesMCHprop;
        if (configRealign.fDoRealign) {
          FillPropagation(mchrealigned, VarColl{}, fgValuesMCHprop, kToZ, fgValuesClsTmp.posClusters[iCls][2]);
        } else {
          FillPropagation<1>(mchtrack, VarColl{}, VarTrack{}, fgValuesMCHprop, kToZ, fgValuesClsTmp.posClusters[iCls][2]);
        }

        //// Fill residual QA histograms
        int topBottom = (fgValuesMCH.y >= 0) ? 0 : 1;
        int posNeg = (fgValuesMCH.sign >= 0) ? 0 : 1;
        FillResidualHistograms<1, 0>(fgValuesMCHprop, fgValuesMCH, fgValuesMCHpv, fgValuesMFT, fgValuesClsTmp.posClusters[iCls][0], fgValuesClsTmp.posClusters[iCls][1], topBottom, posNeg, quadrantCls, chamber, deIndex, sameEvent, mixedEvent);
        resetVar(fgValuesMCHprop);

        //// MFT residuals
        if (IsGoodMFT(fgValuesMFT, configMFTs.fTrackChi2MftUp, configMFTs.fTrackNClustMftLow)) {
          //// Propagate MFT track to given cluster
          VarTrack fgValuesMFTprop;
          FillPropagation<0, 1>(mfttrack, VarColl{}, fgValuesMCH, fgValuesMFTprop, kToZ, fgValuesClsTmp.posClusters[iCls][2]);

          //// Fill residual QA histograms for MFT
          topBottom = (mfttrack.y() >= 0) ? 0 : 1;
          posNeg = (fgValuesMCH.sign >= 0) ? 0 : 1;
          FillResidualHistograms<0, 1>(fgValuesMFTprop, fgValuesMCH, fgValuesMCHpv, fgValuesMFT, fgValuesClsTmp.posClusters[iCls][0], fgValuesClsTmp.posClusters[iCls][1], topBottom, posNeg, quadrantCls, chamber, deIndex, sameEvent, mixedEvent);
          resetVar(fgValuesMFTprop);
        }
      }
    }
  }

  template <typename TEvents, typename TBcs, typename TFwdTracks, typename TMFTTracks, typename TMap>
  void runEventSelection(TEvents const& collisions, TBcs const& bcs, TFwdTracks const& muons, TMFTTracks const& mfts, TMap& collisionSel)
  {
    for (auto const& collision : collisions) {

      uint64_t collisionIndex = collision.globalIndex();
      auto muonsThisCollision = muons.sliceBy(fwdtracksPerCollision, collisionIndex);
      auto mftsThisCollision = mfts.sliceBy(mftPerCollision, collisionIndex);

      if (muonsThisCollision.size() < 1 && mftsThisCollision.size() < 1) {
        continue;
      }

      auto& fgValuesColl = collisionSel[collisionIndex];
      FillCollision(collision, fgValuesColl);
      fgValuesColl.bc = bcs.rawIteratorAt(collision.bcId()).globalBC();
      fgValuesColl.multMFT = mftsThisCollision.size();
    }
  }

  template <typename TEventMap, typename TCandidateMap, typename TFwdTracks, typename TMFTTracks, typename TMuonCls>
  void runMuonQA(TEventMap const& collisions, TCandidateMap& matchingCandidates, TFwdTracks const& muons, TMFTTracks const& mfts, TMuonCls const& clusters)
  {
    //// First loop over all muon tracks
    for (auto const& muon : muons) {

      //// Get collision information if associated
      VarColl fgValuesColl;
      if (muon.has_collision()) {
        fgValuesColl = collisions.at(muon.collisionId());
      } else {
        continue;
      }

      if (static_cast<int>(muon.trackType()) <= 2) { // MFT-MCH-MID(0) or MFT-MCH(2)

        registry.get<TH1>(HIST("global-muons/nTracksPerType"))->Fill(static_cast<int>(muon.trackType()));

        // auto mfttrack = muon.template matchMFTTrack_as<TMFTTracks>(); // unused parameter?
        auto mchtrack = muon.template matchMCHTrack_as<TFwdTracks>();

        // Fill global matching candidates: global muons per MCH track
        FillMatchingCandidates(muon, mchtrack, matchingCandidates);

      } else { // MCH-MID(3) or MCH(4)

        // Fill MCH tracks
        FillTrack<1>(muon, fgValuesMCH);

        // Propagate MCH to PV
        FillPropagation<1>(muon, fgValuesColl, fgValuesMCH, fgValuesMCHpv); // copied in a separate variable

        //// Fill MCH clusters: re-align clusters if required
        TrackRealigned mchrealigned;
        VarClusters fgValuesCls;
        if (!FillClusters(muon, clusters, fgValuesCls, mchrealigned)) {
          continue; // if refit was not passed
        }

        //// Update MCH tracks kinematics if using realigned muons
        if (configRealign.fDoRealign) {

          // Update track info
          FillTrack(mchrealigned, fgValuesMCH);

          // Update propagate of MCH to PV
          FillPropagation(mchrealigned, fgValuesColl, fgValuesMCHpv);

          // Update pDCA and Rabs values
          FillPropagation(mchrealigned, fgValuesColl, fgValuesMCH, kToAbsEnd);
          FillPropagation(mchrealigned, fgValuesColl, fgValuesMCH, kToDCA);
        }

        //// Fill muon QA histograms
        FillMuonHistograms<1, 0, 0>(fgValuesMCH, fgValuesMCHpv, fgValuesMFT, VarTrack{}, nullptr);

        //// Fill muon DCA QA checks
        if (configQAs.fEnableQADCA) {
          // mchrealigned is a dummy variable in the first argument (but not the second!)
          auto const& dummyMFT = *mfts.begin();
          runDCA<1, 0, 0>(collisions, mfts, dummyMFT, muon, mch::Track(), mchrealigned, fgValuesColl, fgValuesMCH, fgValuesMCHpv, fgValuesMFT);
        }
      }

      resetVar(fgValuesMFT);
      resetVar(fgValuesMCH);
      resetVar(fgValuesMCHpv);
    }

    //// Second loop over global muon tracks
    for (auto& [mchIndex, globalMuonsVector] : matchingCandidates) {

      //// sort matching candidates in ascending order based on the matching chi2
      auto compareChi2 = [&muons](uint64_t trackIndex1, uint64_t trackIndex2) -> bool {
        auto const& track1 = muons.rawIteratorAt(trackIndex1);
        auto const& track2 = muons.rawIteratorAt(trackIndex2);

        return (track1.chi2MatchMCHMFT() < track2.chi2MatchMCHMFT());
      };
      std::sort(globalMuonsVector.begin(), globalMuonsVector.end(), compareChi2);

      //// Get tracks
      auto muontrack = muons.rawIteratorAt(globalMuonsVector[0]);
      auto mchtrack = muontrack.template matchMCHTrack_as<MyMuonsWithCov>();
      auto mfttrack = muontrack.template matchMFTTrack_as<MyMFTs>();

      //// Fill matching chi2
      FillMatching(muontrack, fgValuesMCH, fgValuesMFT);

      //// Fill global informations
      registry.get<TH1>(HIST("global-muons/NCandidates"))->Fill(static_cast<int>(globalMuonsVector.size()));
      for (size_t candidateIndex = 0; candidateIndex < globalMuonsVector.size(); candidateIndex++) {
        auto const& muon = muons.rawIteratorAt(globalMuonsVector[candidateIndex]);
        registry.get<TH2>(HIST("global-muons/MatchChi2"))->Fill(muon.chi2MatchMCHMFT(), candidateIndex);
      }

      //// Fill collision information if avalaible
      auto& fgValuesCollGlo = collisions.at(muontrack.collisionId());
      VarColl fgValuesCollMCH; // in principal should be the same as global muon
      if (mchtrack.has_collision()) {
        fgValuesCollMCH = collisions.at(mchtrack.collisionId());
      }

      //// Fill MCH and MFT tracks: basic info copied from input tables
      FillTrack<1>(mchtrack, fgValuesMCH);
      FillTrack<0>(mfttrack, fgValuesMFT);
      FillTrack<1>(muontrack, fgValuesGlobal);

      //// Propagate MCH to PV
      FillPropagation<1>(mchtrack, fgValuesCollMCH, VarTrack{}, fgValuesMCHpv); // saved in separate variable fgValuesMCHpv

      //// Propagate MFT to PV?
      FillPropagation<0, 1>(mfttrack, fgValuesCollGlo, VarTrack{}, fgValuesMFT, kToDCA);

      //// Fill MCH clusters: re-align clusters if required
      TrackRealigned mchrealigned;
      VarClusters fgValuesCls;
      if (!FillClusters(mchtrack, clusters, fgValuesCls, mchrealigned)) {
        continue; // if refit was not passed
      }

      //// Update MCH tracks kinematics if using realigned muons
      if (configRealign.fDoRealign) {
        // Update track info
        FillTrack(mchrealigned, fgValuesMCH);

        // Update propagation of MCH to PV
        FillPropagation(mchrealigned, fgValuesCollMCH, fgValuesMCHpv);

        // Update pDCA and Rabs values
        FillPropagation(mchrealigned, fgValuesCollMCH, fgValuesMCH, kToAbsEnd);
        FillPropagation(mchrealigned, fgValuesCollMCH, fgValuesMCH, kToDCA);
      }

      //// Fill global muon candidates info
      for (int i = 0; i < static_cast<int>(globalMuonsVector.size()); i++) {
        VarTrack fgValuesTmp;
        auto muonCandidate = muons.rawIteratorAt(globalMuonsVector[i]);
        FillTrack<0>(muonCandidate, fgValuesTmp);
        fgValuesCandidates.emplace_back(fgValuesTmp);
      }

      //// Fill global muons QA : fill global matching QA if required
      if (configQAs.fEnableQAMatching) {

        // Propagate global muon tracks to DCA: treat it as MFT using p from MCH?
        // Use here fgValuesMCHpv or fgValuesMCH?
        FillPropagation<0, 1>(muontrack, fgValuesCollGlo, fgValuesMCHpv, fgValuesGlobal, kToDCA);

        // Fill bc difference of matched MCH and MFT
        if (muontrack.has_collision() && mfttrack.has_collision()) {
          fgValuesMCH.bc = collisions.at(muontrack.collisionId()).bc;
          fgValuesMFT.bc = collisions.at(mfttrack.collisionId()).bc;
          int64_t dbc = fgValuesMCH.bc - fgValuesMFT.bc;
          registry.get<TH1>(HIST("global-matches/BCdifference"))->Fill(dbc);
        }

        // Fill QA histograms including global matching
        FillMuonHistograms<0, 1, 1>(fgValuesMCH, fgValuesMCHpv, fgValuesMFT, fgValuesGlobal, fgValuesCandidates);

      } else {
        // Fill QA histograms
        FillMuonHistograms<0, 1, 0>(fgValuesMCH, fgValuesMCHpv, fgValuesMFT, fgValuesGlobal, fgValuesCandidates);
      }

      //// Fill residual QA checks if required
      if (configQAs.fEnableQAResidual) {
        runResidual(collisions, muons, mfts, clusters, mchtrack, mchrealigned, mfttrack, fgValuesCollGlo, fgValuesMCH, fgValuesMCHpv, fgValuesMFT);
      }

      //// Fill MFT and global muon DCA QA checks if required
      if (configQAs.fEnableQADCA) {
        auto const& dummyMuon = *muons.begin();
        // mchrealigned is a dummy variable in the first call
        runDCA<0, 1, 0>(collisions, mfts, mfttrack, dummyMuon, mchrealigned, mchrealigned, fgValuesCollGlo, fgValuesMCH, fgValuesMCHpv, fgValuesMFT);
        // Now fill global DCAs and compare (scaled) MFT and MCH
        runDCA<0, 0, 1>(collisions, mfts, mfttrack, dummyMuon, mchtrack, mchrealigned, fgValuesCollGlo, fgValuesMCH, fgValuesMCHpv, fgValuesMFT);
      }

      fgValuesCandidates.clear();
      resetVar(fgValuesMFT);
      resetVar(fgValuesMCH);
      resetVar(fgValuesMCHpv);
      resetVar(fgValuesGlobal);
    }
  }

  template <typename TEventMap, typename TCandidateMap, typename TFwdTracks, typename TMuonCls>
  void runDimuonQA(TEventMap const& collisions, TCandidateMap const& matchingCandidates, TFwdTracks const& muonTracks, TMuonCls const& clusters)
  {
    std::vector<MuonPair> muonPairs;
    std::vector<GlobalMuonPair> globalMuonPairs;

    GetMuonPairs(muonTracks, matchingCandidates, collisions, muonPairs, globalMuonPairs);

    for (const auto& [muon1, muon2] : muonPairs) {
      auto collisionIndex1 = muon1.first;
      auto const& collision1 = collisions.at(collisionIndex1);
      auto collisionIndex2 = muon2.first;
      auto const& collision2 = collisions.at(collisionIndex2);

      auto mchIndex1 = muon1.second;
      auto mchIndex2 = muon2.second;
      auto const& muonTrack1 = muonTracks.rawIteratorAt(mchIndex1);
      auto const& muonTrack2 = muonTracks.rawIteratorAt(mchIndex2);

      VarTrack fgValuesMuon1, fgValuesMuonPV1;
      VarTrack fgValuesMuon2, fgValuesMuonPV2;
      TrackRealigned mchrealigned1, mchrealigned2;
      VarClusters fgValuesCls1, fgValuesCls2;
      if (!FillClusters(muonTrack1, clusters, fgValuesCls1, mchrealigned1) || !FillClusters(muonTrack2, clusters, fgValuesCls2, mchrealigned2)) {
        continue; // Refit is not valid
      }

      if (configRealign.fDoRealign) {

        FillTrack(mchrealigned1, fgValuesMuon1);
        FillTrack(mchrealigned2, fgValuesMuon2);

        // Propagate MCH to PV
        FillPropagation(mchrealigned1, collision1, fgValuesMuonPV1);
        FillPropagation(mchrealigned2, collision2, fgValuesMuonPV2);

        // Recalculate pDCA and Rabs values
        FillPropagation(mchrealigned1, collision1, fgValuesMuon1, kToAbsEnd);
        FillPropagation(mchrealigned1, collision1, fgValuesMuon1, kToDCA);
        FillPropagation(mchrealigned2, collision2, fgValuesMuon2, kToAbsEnd);
        FillPropagation(mchrealigned2, collision2, fgValuesMuon2, kToDCA);
      } else {
        FillTrack<1>(muonTrack1, fgValuesMuon1);
        FillTrack<1>(muonTrack2, fgValuesMuon2);

        // Propagate MCH to PV
        FillPropagation<1>(muonTrack1, collision1, fgValuesMuon1, fgValuesMuonPV1);
        FillPropagation<1>(muonTrack2, collision2, fgValuesMuon2, fgValuesMuonPV2);
        // Calculate DCA
        FillPropagation<1>(muonTrack1, collision1, fgValuesMuon1, fgValuesMuonPV1, kToDCA);
        FillPropagation<1>(muonTrack2, collision2, fgValuesMuon2, fgValuesMuonPV2, kToDCA);
      }

      int sign1 = muonTrack1.sign();
      int sign2 = muonTrack2.sign();

      // only consider opposite-sign pairs
      if ((sign1 * sign2) >= 0)
        continue;

      const auto& muonPos = fgValuesMuon1.sign > 0 ? fgValuesMuon1 : fgValuesMuon2;
      const auto& muonNeg = fgValuesMuon1.sign < 0 ? fgValuesMuon1 : fgValuesMuon2;
      // for DCA
      const auto& muonPosPV = fgValuesMuon1.sign > 0 ? fgValuesMuonPV1 : fgValuesMuonPV2;
      const auto& muonNegPV = fgValuesMuon1.sign < 0 ? fgValuesMuonPV1 : fgValuesMuonPV2;
      // μ⁺ variables
      double muPosPt = muonPos.pT;
      double muPosEta = muonPos.eta;
      double muPosRabs = muonPos.rabs;
      double muPosPhi = muonPos.phi * 180.0 / TMath::Pi();
      double muPosDca = std::sqrt(muonPos.dcaX * muonPos.dcaX + muonPos.dcaY * muonPos.dcaY);
      // μ⁻ variables
      double muNegPt = muonNeg.pT;
      double muNegEta = muonNeg.eta;
      double muNegRabs = muonNeg.rabs;
      double muNegPhi = muonNeg.phi * 180.0 / TMath::Pi();
      double muNegDca = std::sqrt(muonNeg.dcaX * muonNeg.dcaX + muonNeg.dcaY * muonNeg.dcaY);

      int Quadrant1 = GetQuadrantPhi(muonTrack1.phi() * 180.0 / TMath::Pi());
      int Quadrant2 = GetQuadrantPhi(muonTrack2.phi() * 180.0 / TMath::Pi());
      int TopBottom1 = (Quadrant1 == 0 || Quadrant1 == 1) ? 0 : 1;
      int TopBottom2 = (Quadrant2 == 0 || Quadrant2 == 1) ? 0 : 1;
      int LeftRight1 = (Quadrant1 == 0 || Quadrant1 == 3) ? 0 : 1;
      int LeftRight2 = (Quadrant2 == 0 || Quadrant2 == 3) ? 0 : 1;

      bool goodMuonTracks = (IsGoodMuon(fgValuesMuon1, fgValuesMuonPV1) && IsGoodMuon(fgValuesMuon2, fgValuesMuonPV2));
      bool goodGlobalMuonTracks = (IsGoodGlobalMuon(fgValuesMuon1, fgValuesMuonPV1) && IsGoodGlobalMuon(fgValuesMuon2, fgValuesMuonPV2));

      bool sameEvent = (collisionIndex1 == collisionIndex2);

      // dimuon variables
      double mass = GetMuMuInvariantMass(fgValuesMuonPV1, fgValuesMuonPV2);
      double pT = GetMuMuPt(fgValuesMuonPV1, fgValuesMuonPV2);
      double yPair = GetMuMuRap(fgValuesMuonPV1, fgValuesMuonPV2);
      double dcaXPair;
      double dcaYPair;
      if (configRealign.fDoRealign) {
        dcaXPair = muonPos.dcaX - muonNeg.dcaX;
        dcaYPair = muonPos.dcaY - muonNeg.dcaY;
      } else {
        dcaXPair = muonPosPV.dcaX - muonNegPV.dcaX;
        dcaYPair = muonPosPV.dcaY - muonNegPV.dcaY;
      }
      if (goodMuonTracks) {
        if (sameEvent) {
          // same-event case
          if (configQAs.fEnableSingleMuonDiMuonCorrelations) {
            // single muons
            registryDimuon.get<TH3>(HIST("dimuon/same-event/single-muon-dimuon-correlations/invariantMass_pT_MuPosPt_MuonKine_MuonCuts"))->Fill(mass, pT, muPosPt);
            registryDimuon.get<TH3>(HIST("dimuon/same-event/single-muon-dimuon-correlations/invariantMass_pT_MuNegPt_MuonKine_MuonCuts"))->Fill(mass, pT, muNegPt);
            //
            registryDimuon.get<TH3>(HIST("dimuon/same-event/single-muon-dimuon-correlations/invariantMass_pT_MuPosEta_MuonKine_MuonCuts"))->Fill(mass, pT, muPosEta);
            registryDimuon.get<TH3>(HIST("dimuon/same-event/single-muon-dimuon-correlations/invariantMass_pT_MuNegEta_MuonKine_MuonCuts"))->Fill(mass, pT, muNegEta);
            //
            registryDimuon.get<TH3>(HIST("dimuon/same-event/single-muon-dimuon-correlations/invariantMass_pT_MuPosRabs_MuonKine_MuonCuts"))->Fill(mass, pT, muPosRabs);
            registryDimuon.get<TH3>(HIST("dimuon/same-event/single-muon-dimuon-correlations/invariantMass_pT_MuNegRabs_MuonKine_MuonCuts"))->Fill(mass, pT, muNegRabs);
            //
            registryDimuon.get<TH3>(HIST("dimuon/same-event/single-muon-dimuon-correlations/invariantMass_pT_MuPosDca_MuonKine_MuonCuts"))->Fill(mass, pT, muPosDca);
            registryDimuon.get<TH3>(HIST("dimuon/same-event/single-muon-dimuon-correlations/invariantMass_pT_MuNegDca_MuonKine_MuonCuts"))->Fill(mass, pT, muNegDca);
            //
            registryDimuon.get<TH3>(HIST("dimuon/same-event/single-muon-dimuon-correlations/invariantMass_pT_MuPosPhi_MuonKine_MuonCuts"))->Fill(mass, pT, muPosPhi);
            registryDimuon.get<TH3>(HIST("dimuon/same-event/single-muon-dimuon-correlations/invariantMass_pT_MuNegPhi_MuonKine_MuonCuts"))->Fill(mass, pT, muNegPhi);
          }
          registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMass_MuonKine_MuonCuts"))->Fill(mass);
          registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMassFull_MuonKine_MuonCuts"))->Fill(mass);
          registryDimuon.get<TH2>(HIST("dimuon/same-event/invariantMass_pT_MuonKine_MuonCuts"))->Fill(mass, pT);
          registryDimuon.get<TH1>(HIST("dimuon/same-event/rapPair_MuonKine_MuonCuts"))->Fill(yPair);
          registryDimuon.get<TH2>(HIST("dimuon/same-event/invariantMass_rapPair_MuonKine_MuonCuts"))->Fill(mass, yPair);
          registryDimuon.get<TH2>(HIST("dimuon/same-event/pT_rapPair_MuonKine_MuonCuts"))->Fill(pT, yPair);

          // dimuon DCA
          if (mass >= fDimuonDCAMassLow && mass <= fDimuonDCAMassHigh) {
            registryDimuon.get<TH2>(HIST("dimuon/same-event/opposite-sign/DCA/pT_MuPosDCAx_minus_MuNegDCAx_MuonKine_MuonCuts"))->Fill(pT, dcaXPair);
            registryDimuon.get<TH2>(HIST("dimuon/same-event/opposite-sign/DCA/pT_MuPosDCAy_minus_MuNegDCAy_MuonKine_MuonCuts"))->Fill(pT, dcaYPair);
          }

          // dimuon top-bottom and left-right separation
          if (TopBottom1 == 0 && TopBottom2 == 0) {
            registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMass_MuonKine_MuonCuts_TT"))->Fill(mass);
            registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMassFull_MuonKine_MuonCuts_TT"))->Fill(mass);
            registryDimuon.get<TH2>(HIST("dimuon/same-event/invariantMass_pT_MuonKine_MuonCuts_TT"))->Fill(mass, pT);
            registryDimuon.get<TH1>(HIST("dimuon/same-event/rapPair_MuonKine_MuonCuts_TT"))->Fill(yPair);
            registryDimuon.get<TH2>(HIST("dimuon/same-event/invariantMass_rapPair_MuonKine_MuonCuts_TT"))->Fill(mass, yPair);
            registryDimuon.get<TH2>(HIST("dimuon/same-event/pT_rapPair_MuonKine_MuonCuts_TT"))->Fill(pT, yPair);
            if (mass >= fDimuonDCAMassLow && mass <= fDimuonDCAMassHigh) {
              registryDimuon.get<TH2>(HIST("dimuon/same-event/opposite-sign/DCA/pT_MuPosTDCAx_minus_MuNegTDCAx_MuonKine_MuonCuts"))->Fill(pT, dcaXPair);
              registryDimuon.get<TH2>(HIST("dimuon/same-event/opposite-sign/DCA/pT_MuPosTDCAy_minus_MuNegTDCAy_MuonKine_MuonCuts"))->Fill(pT, dcaYPair);
            }
          } else if ((TopBottom1 == 0 && TopBottom2 == 1) || (TopBottom1 == 1 && TopBottom2 == 0)) {
            registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMass_MuonKine_MuonCuts_TB"))->Fill(mass);
            registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMassFull_MuonKine_MuonCuts_TB"))->Fill(mass);
            registryDimuon.get<TH2>(HIST("dimuon/same-event/invariantMass_pT_MuonKine_MuonCuts_TB"))->Fill(mass, pT);
            registryDimuon.get<TH1>(HIST("dimuon/same-event/rapPair_MuonKine_MuonCuts_TB"))->Fill(yPair);
            registryDimuon.get<TH2>(HIST("dimuon/same-event/invariantMass_rapPair_MuonKine_MuonCuts_TB"))->Fill(mass, yPair);
            registryDimuon.get<TH2>(HIST("dimuon/same-event/pT_rapPair_MuonKine_MuonCuts_TB"))->Fill(pT, yPair);
            if (TopBottom1 == 0 && TopBottom2 == 1) {
              if (sign1 > 0) {
                registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMass_MuonKine_MuonCuts_TPBN"))->Fill(mass);
                registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMassFull_MuonKine_MuonCuts_TPBN"))->Fill(mass);
                registryDimuon.get<TH2>(HIST("dimuon/same-event/invariantMass_pT_MuonKine_MuonCuts_TPBN"))->Fill(mass, pT);
                if (mass >= fDimuonDCAMassLow && mass <= fDimuonDCAMassHigh) {
                  registryDimuon.get<TH2>(HIST("dimuon/same-event/opposite-sign/DCA/pT_MuPosTDCAx_minus_MuNegBDCAx_MuonKine_MuonCuts"))->Fill(pT, dcaXPair);
                  registryDimuon.get<TH2>(HIST("dimuon/same-event/opposite-sign/DCA/pT_MuPosTDCAy_minus_MuNegBDCAy_MuonKine_MuonCuts"))->Fill(pT, dcaYPair);
                }
              } else {
                registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMass_MuonKine_MuonCuts_TNBP"))->Fill(mass);
                registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMassFull_MuonKine_MuonCuts_TNBP"))->Fill(mass);
                registryDimuon.get<TH2>(HIST("dimuon/same-event/invariantMass_pT_MuonKine_MuonCuts_TNBP"))->Fill(mass, pT);
                if (mass >= fDimuonDCAMassLow && mass <= fDimuonDCAMassHigh) {
                  registryDimuon.get<TH2>(HIST("dimuon/same-event/opposite-sign/DCA/pT_MuPosBDCAx_minus_MuNegTDCAx_MuonKine_MuonCuts"))->Fill(pT, dcaXPair);
                  registryDimuon.get<TH2>(HIST("dimuon/same-event/opposite-sign/DCA/pT_MuPosBDCAy_minus_MuNegTDCAy_MuonKine_MuonCuts"))->Fill(pT, dcaYPair);
                }
              }
            } else if (TopBottom1 == 1 && TopBottom2 == 0) {
              if (sign2 > 0) {
                registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMass_MuonKine_MuonCuts_TPBN"))->Fill(mass);
                registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMassFull_MuonKine_MuonCuts_TPBN"))->Fill(mass);
                registryDimuon.get<TH2>(HIST("dimuon/same-event/invariantMass_pT_MuonKine_MuonCuts_TPBN"))->Fill(mass, pT);
                if (mass >= fDimuonDCAMassLow && mass <= fDimuonDCAMassHigh) {
                  registryDimuon.get<TH2>(HIST("dimuon/same-event/opposite-sign/DCA/pT_MuPosTDCAx_minus_MuNegBDCAx_MuonKine_MuonCuts"))->Fill(pT, dcaXPair);
                  registryDimuon.get<TH2>(HIST("dimuon/same-event/opposite-sign/DCA/pT_MuPosTDCAy_minus_MuNegBDCAy_MuonKine_MuonCuts"))->Fill(pT, dcaYPair);
                }
              } else {
                registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMass_MuonKine_MuonCuts_TNBP"))->Fill(mass);
                registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMassFull_MuonKine_MuonCuts_TNBP"))->Fill(mass);
                registryDimuon.get<TH2>(HIST("dimuon/same-event/invariantMass_pT_MuonKine_MuonCuts_TNBP"))->Fill(mass, pT);
                if (mass >= fDimuonDCAMassLow && mass <= fDimuonDCAMassHigh) {
                  registryDimuon.get<TH2>(HIST("dimuon/same-event/opposite-sign/DCA/pT_MuPosBDCAx_minus_MuNegTDCAx_MuonKine_MuonCuts"))->Fill(pT, dcaXPair);
                  registryDimuon.get<TH2>(HIST("dimuon/same-event/opposite-sign/DCA/pT_MuPosBDCAy_minus_MuNegTDCAy_MuonKine_MuonCuts"))->Fill(pT, dcaYPair);
                }
              }
            }
          } else if (TopBottom1 == 1 && TopBottom2 == 1) {
            registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMass_MuonKine_MuonCuts_BB"))->Fill(mass);
            registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMassFull_MuonKine_MuonCuts_BB"))->Fill(mass);
            registryDimuon.get<TH2>(HIST("dimuon/same-event/invariantMass_pT_MuonKine_MuonCuts_BB"))->Fill(mass, pT);
            registryDimuon.get<TH1>(HIST("dimuon/same-event/rapPair_MuonKine_MuonCuts_BB"))->Fill(yPair);
            registryDimuon.get<TH2>(HIST("dimuon/same-event/invariantMass_rapPair_MuonKine_MuonCuts_BB"))->Fill(mass, yPair);
            registryDimuon.get<TH2>(HIST("dimuon/same-event/pT_rapPair_MuonKine_MuonCuts_BB"))->Fill(pT, yPair);
            if (mass >= fDimuonDCAMassLow && mass <= fDimuonDCAMassHigh) {
              registryDimuon.get<TH2>(HIST("dimuon/same-event/opposite-sign/DCA/pT_MuPosBDCAx_minus_MuNegBDCAx_MuonKine_MuonCuts"))->Fill(pT, dcaXPair);
              registryDimuon.get<TH2>(HIST("dimuon/same-event/opposite-sign/DCA/pT_MuPosBDCAy_minus_MuNegBDCAy_MuonKine_MuonCuts"))->Fill(pT, dcaYPair);
            }
          }

          if (LeftRight1 == 0 && LeftRight2 == 0) {
            registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMass_MuonKine_MuonCuts_LL"))->Fill(mass);
            registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMassFull_MuonKine_MuonCuts_LL"))->Fill(mass);
            registryDimuon.get<TH2>(HIST("dimuon/same-event/invariantMass_pT_MuonKine_MuonCuts_LL"))->Fill(mass, pT);
            registryDimuon.get<TH1>(HIST("dimuon/same-event/rapPair_MuonKine_MuonCuts_LL"))->Fill(yPair);
            registryDimuon.get<TH2>(HIST("dimuon/same-event/invariantMass_rapPair_MuonKine_MuonCuts_LL"))->Fill(mass, yPair);
          } else if ((LeftRight1 == 0 && LeftRight2 == 1) || (LeftRight1 == 1 && LeftRight2 == 0)) {
            registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMass_MuonKine_MuonCuts_LR"))->Fill(mass);
            registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMassFull_MuonKine_MuonCuts_LR"))->Fill(mass);
            registryDimuon.get<TH2>(HIST("dimuon/same-event/invariantMass_pT_MuonKine_MuonCuts_LR"))->Fill(mass, pT);
            registryDimuon.get<TH1>(HIST("dimuon/same-event/rapPair_MuonKine_MuonCuts_LR"))->Fill(yPair);
            registryDimuon.get<TH2>(HIST("dimuon/same-event/invariantMass_rapPair_MuonKine_MuonCuts_LR"))->Fill(mass, yPair);
            if (LeftRight1 == 0 && LeftRight2 == 1) {
              if (sign1 > 0) {
                registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMass_MuonKine_MuonCuts_LPRN"))->Fill(mass);
                registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMassFull_MuonKine_MuonCuts_LPRN"))->Fill(mass);
                registryDimuon.get<TH2>(HIST("dimuon/same-event/invariantMass_pT_MuonKine_MuonCuts_LPRN"))->Fill(mass, pT);
              } else {
                registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMass_MuonKine_MuonCuts_LNRP"))->Fill(mass);
                registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMassFull_MuonKine_MuonCuts_LNRP"))->Fill(mass);
                registryDimuon.get<TH2>(HIST("dimuon/same-event/invariantMass_pT_MuonKine_MuonCuts_LNRP"))->Fill(mass, pT);
              }
            } else if (LeftRight1 == 1 && LeftRight2 == 0) {
              if (sign2 > 0) {
                registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMass_MuonKine_MuonCuts_LPRN"))->Fill(mass);
                registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMassFull_MuonKine_MuonCuts_LPRN"))->Fill(mass);
                registryDimuon.get<TH2>(HIST("dimuon/same-event/invariantMass_pT_MuonKine_MuonCuts_LPRN"))->Fill(mass, pT);
              } else {
                registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMass_MuonKine_MuonCuts_LNRP"))->Fill(mass);
                registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMassFull_MuonKine_MuonCuts_LNRP"))->Fill(mass);
                registryDimuon.get<TH2>(HIST("dimuon/same-event/invariantMass_pT_MuonKine_MuonCuts_LNRP"))->Fill(mass, pT);
              }
            }
          } else if (LeftRight1 == 1 && LeftRight2 == 1) {
            registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMass_MuonKine_MuonCuts_RR"))->Fill(mass);
            registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMassFull_MuonKine_MuonCuts_RR"))->Fill(mass);
            registryDimuon.get<TH2>(HIST("dimuon/same-event/invariantMass_pT_MuonKine_MuonCuts_RR"))->Fill(mass, pT);
            registryDimuon.get<TH1>(HIST("dimuon/same-event/rapPair_MuonKine_MuonCuts_RR"))->Fill(yPair);
            registryDimuon.get<TH2>(HIST("dimuon/same-event/invariantMass_rapPair_MuonKine_MuonCuts_RR"))->Fill(mass, yPair);
          }
        } else {
          // event-mixing case
          registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMass_MuonKine_MuonCuts"))->Fill(mass);
          registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMassFull_MuonKine_MuonCuts"))->Fill(mass);
          registryDimuon.get<TH2>(HIST("dimuon/mixed-event/invariantMass_pT_MuonKine_MuonCuts"))->Fill(mass, pT);

          // dimuon DCA
          if (mass >= fDimuonDCAMassLow && mass <= fDimuonDCAMassHigh) {
            registryDimuon.get<TH2>(HIST("dimuon/mixed-event/DCA/pT_MuPosDCAx_minus_MuNegDCAx_MuonKine_MuonCuts"))->Fill(pT, dcaXPair);
            registryDimuon.get<TH2>(HIST("dimuon/mixed-event/DCA/pT_MuPosDCAy_minus_MuNegDCAy_MuonKine_MuonCuts"))->Fill(pT, dcaYPair);
          }

          // dimuon top-bottom and left-right separation
          if (TopBottom1 == 0 && TopBottom2 == 0) {
            registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMass_MuonKine_MuonCuts_TT"))->Fill(mass);
            registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMassFull_MuonKine_MuonCuts_TT"))->Fill(mass);
            registryDimuon.get<TH2>(HIST("dimuon/mixed-event/invariantMass_pT_MuonKine_MuonCuts_TT"))->Fill(mass, pT);
            if (mass >= fDimuonDCAMassLow && mass <= fDimuonDCAMassHigh) {
              registryDimuon.get<TH2>(HIST("dimuon/mixed-event/DCA/pT_MuPosTDCAx_minus_MuNegTDCAx_MuonKine_MuonCuts"))->Fill(pT, dcaXPair);
              registryDimuon.get<TH2>(HIST("dimuon/mixed-event/DCA/pT_MuPosTDCAy_minus_MuNegTDCAy_MuonKine_MuonCuts"))->Fill(pT, dcaYPair);
            }
          } else if ((TopBottom1 == 0 && TopBottom2 == 1) || (TopBottom1 == 1 && TopBottom2 == 0)) {
            registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMass_MuonKine_MuonCuts_TB"))->Fill(mass);
            registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMassFull_MuonKine_MuonCuts_TB"))->Fill(mass);
            registryDimuon.get<TH2>(HIST("dimuon/mixed-event/invariantMass_pT_MuonKine_MuonCuts_TB"))->Fill(mass, pT);
            if (TopBottom1 == 0 && TopBottom2 == 1) {
              if (sign1 > 0) {
                registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMass_MuonKine_MuonCuts_TPBN"))->Fill(mass);
                registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMassFull_MuonKine_MuonCuts_TPBN"))->Fill(mass);
                registryDimuon.get<TH2>(HIST("dimuon/mixed-event/invariantMass_pT_MuonKine_MuonCuts_TPBN"))->Fill(mass, pT);
                if (mass >= fDimuonDCAMassLow && mass <= fDimuonDCAMassHigh) {
                  registryDimuon.get<TH2>(HIST("dimuon/mixed-event/DCA/pT_MuPosTDCAx_minus_MuNegBDCAx_MuonKine_MuonCuts"))->Fill(pT, dcaXPair);
                  registryDimuon.get<TH2>(HIST("dimuon/mixed-event/DCA/pT_MuPosTDCAy_minus_MuNegBDCAy_MuonKine_MuonCuts"))->Fill(pT, dcaYPair);
                }
              } else {
                registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMass_MuonKine_MuonCuts_TNBP"))->Fill(mass);
                registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMassFull_MuonKine_MuonCuts_TNBP"))->Fill(mass);
                registryDimuon.get<TH2>(HIST("dimuon/mixed-event/invariantMass_pT_MuonKine_MuonCuts_TNBP"))->Fill(mass, pT);
                if (mass >= fDimuonDCAMassLow && mass <= fDimuonDCAMassHigh) {
                  registryDimuon.get<TH2>(HIST("dimuon/mixed-event/DCA/pT_MuPosBDCAx_minus_MuNegTDCAx_MuonKine_MuonCuts"))->Fill(pT, dcaXPair);
                  registryDimuon.get<TH2>(HIST("dimuon/mixed-event/DCA/pT_MuPosBDCAy_minus_MuNegTDCAy_MuonKine_MuonCuts"))->Fill(pT, dcaYPair);
                }
              }
            } else if (TopBottom1 == 1 && TopBottom2 == 0) {
              if (sign2 > 0) {
                registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMass_MuonKine_MuonCuts_TPBN"))->Fill(mass);
                registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMassFull_MuonKine_MuonCuts_TPBN"))->Fill(mass);
                registryDimuon.get<TH2>(HIST("dimuon/mixed-event/invariantMass_pT_MuonKine_MuonCuts_TPBN"))->Fill(mass, pT);
                if (mass >= fDimuonDCAMassLow && mass <= fDimuonDCAMassHigh) {
                  registryDimuon.get<TH2>(HIST("dimuon/mixed-event/DCA/pT_MuPosTDCAx_minus_MuNegBDCAx_MuonKine_MuonCuts"))->Fill(pT, dcaXPair);
                  registryDimuon.get<TH2>(HIST("dimuon/mixed-event/DCA/pT_MuPosTDCAy_minus_MuNegBDCAy_MuonKine_MuonCuts"))->Fill(pT, dcaYPair);
                }
              } else {
                registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMass_MuonKine_MuonCuts_TNBP"))->Fill(mass);
                registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMassFull_MuonKine_MuonCuts_TNBP"))->Fill(mass);
                registryDimuon.get<TH2>(HIST("dimuon/mixed-event/invariantMass_pT_MuonKine_MuonCuts_TNBP"))->Fill(mass, pT);
                if (mass >= fDimuonDCAMassLow && mass <= fDimuonDCAMassHigh) {
                  registryDimuon.get<TH2>(HIST("dimuon/mixed-event/DCA/pT_MuPosBDCAx_minus_MuNegTDCAx_MuonKine_MuonCuts"))->Fill(pT, dcaXPair);
                  registryDimuon.get<TH2>(HIST("dimuon/mixed-event/DCA/pT_MuPosBDCAy_minus_MuNegTDCAy_MuonKine_MuonCuts"))->Fill(pT, dcaYPair);
                }
              }
            }
          } else if (TopBottom1 == 1 && TopBottom2 == 1) {
            registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMass_MuonKine_MuonCuts_BB"))->Fill(mass);
            registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMassFull_MuonKine_MuonCuts_BB"))->Fill(mass);
            registryDimuon.get<TH2>(HIST("dimuon/mixed-event/invariantMass_pT_MuonKine_MuonCuts_BB"))->Fill(mass, pT);
            if (mass >= fDimuonDCAMassLow && mass <= fDimuonDCAMassHigh) {
              registryDimuon.get<TH2>(HIST("dimuon/mixed-event/DCA/pT_MuPosBDCAx_minus_MuNegBDCAx_MuonKine_MuonCuts"))->Fill(pT, dcaXPair);
              registryDimuon.get<TH2>(HIST("dimuon/mixed-event/DCA/pT_MuPosBDCAy_minus_MuNegBDCAy_MuonKine_MuonCuts"))->Fill(pT, dcaYPair);
            }
          }

          if (LeftRight1 == 0 && LeftRight2 == 0) {
            registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMass_MuonKine_MuonCuts_LL"))->Fill(mass);
            registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMassFull_MuonKine_MuonCuts_LL"))->Fill(mass);
            registryDimuon.get<TH2>(HIST("dimuon/mixed-event/invariantMass_pT_MuonKine_MuonCuts_LL"))->Fill(mass, pT);
          } else if ((LeftRight1 == 0 && LeftRight2 == 1) || (LeftRight1 == 1 && LeftRight2 == 0)) {
            registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMass_MuonKine_MuonCuts_LR"))->Fill(mass);
            registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMassFull_MuonKine_MuonCuts_LR"))->Fill(mass);
            registryDimuon.get<TH2>(HIST("dimuon/mixed-event/invariantMass_pT_MuonKine_MuonCuts_LR"))->Fill(mass, pT);
            if (LeftRight1 == 0 && LeftRight2 == 1) {
              if (sign1 > 0) {
                registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMass_MuonKine_MuonCuts_LPRN"))->Fill(mass);
                registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMassFull_MuonKine_MuonCuts_LPRN"))->Fill(mass);
                registryDimuon.get<TH2>(HIST("dimuon/mixed-event/invariantMass_pT_MuonKine_MuonCuts_LPRN"))->Fill(mass, pT);
              } else {
                registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMass_MuonKine_MuonCuts_LNRP"))->Fill(mass);
                registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMassFull_MuonKine_MuonCuts_LNRP"))->Fill(mass);
                registryDimuon.get<TH2>(HIST("dimuon/mixed-event/invariantMass_pT_MuonKine_MuonCuts_LNRP"))->Fill(mass, pT);
              }
            } else if (LeftRight1 == 1 && LeftRight2 == 0) {
              if (sign2 > 0) {
                registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMass_MuonKine_MuonCuts_LPRN"))->Fill(mass);
                registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMassFull_MuonKine_MuonCuts_LPRN"))->Fill(mass);
                registryDimuon.get<TH2>(HIST("dimuon/mixed-event/invariantMass_pT_MuonKine_MuonCuts_LPRN"))->Fill(mass, pT);
              } else {
                registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMass_MuonKine_MuonCuts_LNRP"))->Fill(mass);
                registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMassFull_MuonKine_MuonCuts_LNRP"))->Fill(mass);
                registryDimuon.get<TH2>(HIST("dimuon/mixed-event/invariantMass_pT_MuonKine_MuonCuts_LNRP"))->Fill(mass, pT);
              }
            }
          } else if (LeftRight1 == 1 && LeftRight2 == 1) {
            registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMass_MuonKine_MuonCuts_RR"))->Fill(mass);
            registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMassFull_MuonKine_MuonCuts_RR"))->Fill(mass);
            registryDimuon.get<TH2>(HIST("dimuon/mixed-event/invariantMass_pT_MuonKine_MuonCuts_RR"))->Fill(mass, pT);
          }
        }
      }

      if (goodGlobalMuonTracks) {
        if (sameEvent) {
          // same-event case
          registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMass_MuonKine_GlobalMuonCuts"))->Fill(mass);
          registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMassFull_MuonKine_GlobalMuonCuts"))->Fill(mass);
          registryDimuon.get<TH2>(HIST("dimuon/same-event/invariantMass_pT_MuonKine_GlobalMuonCuts"))->Fill(mass, pT);

          if (TopBottom1 == 0 && TopBottom2 == 0) {
            registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMass_MuonKine_GlobalMuonCuts_TT"))->Fill(mass);
            registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMassFull_MuonKine_GlobalMuonCuts_TT"))->Fill(mass);
            registryDimuon.get<TH2>(HIST("dimuon/same-event/invariantMass_pT_MuonKine_GlobalMuonCuts_TT"))->Fill(mass, pT);
          } else if ((TopBottom1 == 0 && TopBottom2 == 1) || (TopBottom1 == 1 && TopBottom2 == 0)) {
            registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMass_MuonKine_GlobalMuonCuts_TB"))->Fill(mass);
            registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMassFull_MuonKine_GlobalMuonCuts_TB"))->Fill(mass);
            registryDimuon.get<TH2>(HIST("dimuon/same-event/invariantMass_pT_MuonKine_GlobalMuonCuts_TB"))->Fill(mass, pT);
            if (TopBottom1 == 0 && TopBottom2 == 1) {
              if (sign1 > 0) {
                registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMass_MuonKine_GlobalMuonCuts_TPBN"))->Fill(mass);
                registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMassFull_MuonKine_GlobalMuonCuts_TPBN"))->Fill(mass);
                registryDimuon.get<TH2>(HIST("dimuon/same-event/invariantMass_pT_MuonKine_GlobalMuonCuts_TPBN"))->Fill(mass, pT);
              } else {
                registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMass_MuonKine_GlobalMuonCuts_TNBP"))->Fill(mass);
                registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMassFull_MuonKine_GlobalMuonCuts_TNBP"))->Fill(mass);
                registryDimuon.get<TH2>(HIST("dimuon/same-event/invariantMass_pT_MuonKine_GlobalMuonCuts_TNBP"))->Fill(mass, pT);
              }
            } else if (TopBottom1 == 1 && TopBottom2 == 0) {
              if (sign2 > 0) {
                registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMass_MuonKine_GlobalMuonCuts_TPBN"))->Fill(mass);
                registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMassFull_MuonKine_GlobalMuonCuts_TPBN"))->Fill(mass);
                registryDimuon.get<TH2>(HIST("dimuon/same-event/invariantMass_pT_MuonKine_GlobalMuonCuts_TPBN"))->Fill(mass, pT);
              } else {
                registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMass_MuonKine_GlobalMuonCuts_TNBP"))->Fill(mass);
                registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMassFull_MuonKine_GlobalMuonCuts_TNBP"))->Fill(mass);
                registryDimuon.get<TH2>(HIST("dimuon/same-event/invariantMass_pT_MuonKine_GlobalMuonCuts_TNBP"))->Fill(mass, pT);
              }
            }
          } else if (TopBottom1 == 1 && TopBottom2 == 1) {
            registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMass_MuonKine_GlobalMuonCuts_BB"))->Fill(mass);
            registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMassFull_MuonKine_GlobalMuonCuts_BB"))->Fill(mass);
            registryDimuon.get<TH2>(HIST("dimuon/same-event/invariantMass_pT_MuonKine_GlobalMuonCuts_BB"))->Fill(mass, pT);
          }

          if (LeftRight1 == 0 && LeftRight2 == 0) {
            registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMass_MuonKine_GlobalMuonCuts_LL"))->Fill(mass);
            registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMassFull_MuonKine_GlobalMuonCuts_LL"))->Fill(mass);
            registryDimuon.get<TH2>(HIST("dimuon/same-event/invariantMass_pT_MuonKine_GlobalMuonCuts_LL"))->Fill(mass, pT);
          } else if ((LeftRight1 == 0 && LeftRight2 == 1) || (LeftRight1 == 1 && LeftRight2 == 0)) {
            registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMass_MuonKine_GlobalMuonCuts_LR"))->Fill(mass);
            registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMassFull_MuonKine_GlobalMuonCuts_LR"))->Fill(mass);
            registryDimuon.get<TH2>(HIST("dimuon/same-event/invariantMass_pT_MuonKine_GlobalMuonCuts_LR"))->Fill(mass, pT);
            if (TopBottom1 == 0 && TopBottom2 == 1) {
              if (sign1 > 0) {
                registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMass_MuonKine_GlobalMuonCuts_LPRN"))->Fill(mass);
                registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMassFull_MuonKine_GlobalMuonCuts_LPRN"))->Fill(mass);
                registryDimuon.get<TH2>(HIST("dimuon/same-event/invariantMass_pT_MuonKine_GlobalMuonCuts_LPRN"))->Fill(mass, pT);
              } else {
                registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMass_MuonKine_GlobalMuonCuts_LNRP"))->Fill(mass);
                registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMassFull_MuonKine_GlobalMuonCuts_LNRP"))->Fill(mass);
                registryDimuon.get<TH2>(HIST("dimuon/same-event/invariantMass_pT_MuonKine_GlobalMuonCuts_LNRP"))->Fill(mass, pT);
              }
            } else if (TopBottom1 == 1 && TopBottom2 == 0) {
              if (sign2 > 0) {
                registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMass_MuonKine_GlobalMuonCuts_LPRN"))->Fill(mass);
                registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMassFull_MuonKine_GlobalMuonCuts_LPRN"))->Fill(mass);
                registryDimuon.get<TH2>(HIST("dimuon/same-event/invariantMass_pT_MuonKine_GlobalMuonCuts_LPRN"))->Fill(mass, pT);
              } else {
                registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMass_MuonKine_GlobalMuonCuts_LNRP"))->Fill(mass);
                registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMassFull_MuonKine_GlobalMuonCuts_LNRP"))->Fill(mass);
                registryDimuon.get<TH2>(HIST("dimuon/same-event/invariantMass_pT_MuonKine_GlobalMuonCuts_LNRP"))->Fill(mass, pT);
              }
            }
          } else if (LeftRight1 == 1 && LeftRight2 == 1) {
            registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMass_MuonKine_GlobalMuonCuts_RR"))->Fill(mass);
            registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMassFull_MuonKine_GlobalMuonCuts_RR"))->Fill(mass);
            registryDimuon.get<TH2>(HIST("dimuon/same-event/invariantMass_pT_MuonKine_GlobalMuonCuts_RR"))->Fill(mass, pT);
          }
        } else {
          // event-mixing case
          registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMass_MuonKine_GlobalMuonCuts"))->Fill(mass);
          registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMassFull_MuonKine_GlobalMuonCuts"))->Fill(mass);
          registryDimuon.get<TH2>(HIST("dimuon/mixed-event/invariantMass_pT_MuonKine_GlobalMuonCuts"))->Fill(mass, pT);

          if (TopBottom1 == 0 && TopBottom2 == 0) {
            registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMass_MuonKine_GlobalMuonCuts_TT"))->Fill(mass);
            registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMassFull_MuonKine_GlobalMuonCuts_TT"))->Fill(mass);
            registryDimuon.get<TH2>(HIST("dimuon/mixed-event/invariantMass_pT_MuonKine_GlobalMuonCuts_TT"))->Fill(mass, pT);
          } else if ((TopBottom1 == 0 && TopBottom2 == 1) || (TopBottom1 == 1 && TopBottom2 == 0)) {
            registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMass_MuonKine_GlobalMuonCuts_TB"))->Fill(mass);
            registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMassFull_MuonKine_GlobalMuonCuts_TB"))->Fill(mass);
            registryDimuon.get<TH2>(HIST("dimuon/mixed-event/invariantMass_pT_MuonKine_GlobalMuonCuts_TB"))->Fill(mass, pT);
            if (TopBottom1 == 0 && TopBottom2 == 1) {
              if (sign1 > 0) {
                registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMass_MuonKine_GlobalMuonCuts_TPBN"))->Fill(mass);
                registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMassFull_MuonKine_GlobalMuonCuts_TPBN"))->Fill(mass);
                registryDimuon.get<TH2>(HIST("dimuon/mixed-event/invariantMass_pT_MuonKine_GlobalMuonCuts_TPBN"))->Fill(mass, pT);
              } else {
                registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMass_MuonKine_GlobalMuonCuts_TNBP"))->Fill(mass);
                registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMassFull_MuonKine_GlobalMuonCuts_TNBP"))->Fill(mass);
                registryDimuon.get<TH2>(HIST("dimuon/mixed-event/invariantMass_pT_MuonKine_GlobalMuonCuts_TNBP"))->Fill(mass, pT);
              }
            } else if (TopBottom1 == 1 && TopBottom2 == 0) {
              if (sign2 > 0) {
                registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMass_MuonKine_GlobalMuonCuts_TPBN"))->Fill(mass);
                registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMassFull_MuonKine_GlobalMuonCuts_TPBN"))->Fill(mass);
                registryDimuon.get<TH2>(HIST("dimuon/mixed-event/invariantMass_pT_MuonKine_GlobalMuonCuts_TPBN"))->Fill(mass, pT);
              } else {
                registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMass_MuonKine_GlobalMuonCuts_TNBP"))->Fill(mass);
                registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMassFull_MuonKine_GlobalMuonCuts_TNBP"))->Fill(mass);
                registryDimuon.get<TH2>(HIST("dimuon/mixed-event/invariantMass_pT_MuonKine_GlobalMuonCuts_TNBP"))->Fill(mass, pT);
              }
            }
          } else if (TopBottom1 == 1 && TopBottom2 == 1) {
            registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMass_MuonKine_GlobalMuonCuts_BB"))->Fill(mass);
            registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMassFull_MuonKine_GlobalMuonCuts_BB"))->Fill(mass);
            registryDimuon.get<TH2>(HIST("dimuon/mixed-event/invariantMass_pT_MuonKine_GlobalMuonCuts_BB"))->Fill(mass, pT);
          }

          if (LeftRight1 == 0 && LeftRight2 == 0) {
            registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMass_MuonKine_GlobalMuonCuts_LL"))->Fill(mass);
            registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMassFull_MuonKine_GlobalMuonCuts_LL"))->Fill(mass);
            registryDimuon.get<TH2>(HIST("dimuon/mixed-event/invariantMass_pT_MuonKine_GlobalMuonCuts_LL"))->Fill(mass, pT);
          } else if ((LeftRight1 == 0 && LeftRight2 == 1) || (LeftRight1 == 1 && LeftRight2 == 0)) {
            registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMass_MuonKine_GlobalMuonCuts_LR"))->Fill(mass);
            registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMassFull_MuonKine_GlobalMuonCuts_LR"))->Fill(mass);
            registryDimuon.get<TH2>(HIST("dimuon/mixed-event/invariantMass_pT_MuonKine_GlobalMuonCuts_LR"))->Fill(mass, pT);
            if (LeftRight1 == 0 && LeftRight2 == 1) {
              if (sign1 > 0) {
                registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMass_MuonKine_GlobalMuonCuts_LPRN"))->Fill(mass);
                registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMassFull_MuonKine_GlobalMuonCuts_LPRN"))->Fill(mass);
                registryDimuon.get<TH2>(HIST("dimuon/mixed-event/invariantMass_pT_MuonKine_GlobalMuonCuts_LPRN"))->Fill(mass, pT);
              } else {
                registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMass_MuonKine_GlobalMuonCuts_LNRP"))->Fill(mass);
                registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMassFull_MuonKine_GlobalMuonCuts_LNRP"))->Fill(mass);
                registryDimuon.get<TH2>(HIST("dimuon/mixed-event/invariantMass_pT_MuonKine_GlobalMuonCuts_LNRP"))->Fill(mass, pT);
              }
            } else if (LeftRight1 == 1 && LeftRight2 == 0) {
              if (sign2 > 0) {
                registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMass_MuonKine_GlobalMuonCuts_LPRN"))->Fill(mass);
                registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMassFull_MuonKine_GlobalMuonCuts_LPRN"))->Fill(mass);
                registryDimuon.get<TH2>(HIST("dimuon/mixed-event/invariantMass_pT_MuonKine_GlobalMuonCuts_LPRN"))->Fill(mass, pT);
              } else {
                registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMass_MuonKine_GlobalMuonCuts_LNRP"))->Fill(mass);
                registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMassFull_MuonKine_GlobalMuonCuts_LNRP"))->Fill(mass);
                registryDimuon.get<TH2>(HIST("dimuon/mixed-event/invariantMass_pT_MuonKine_GlobalMuonCuts_LNRP"))->Fill(mass, pT);
              }
            }
          } else if (LeftRight1 == 1 && LeftRight2 == 1) {
            registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMass_MuonKine_GlobalMuonCuts_RR"))->Fill(mass);
            registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMassFull_MuonKine_GlobalMuonCuts_RR"))->Fill(mass);
            registryDimuon.get<TH2>(HIST("dimuon/mixed-event/invariantMass_pT_MuonKine_GlobalMuonCuts_RR"))->Fill(mass, pT);
          }
        }
      }
    }

    if (configQAs.fEnableQADimuonSameSignDCA) {
      for (const auto& [muon1, muon2] : muonPairs) {
        auto collisionIndex1 = muon1.first;
        auto const& collision1 = collisions.at(collisionIndex1);
        auto collisionIndex2 = muon2.first;
        auto const& collision2 = collisions.at(collisionIndex2);

        auto mchIndex1 = muon1.second;
        auto mchIndex2 = muon2.second;
        auto const& muonTrack1 = muonTracks.rawIteratorAt(mchIndex1);
        auto const& muonTrack2 = muonTracks.rawIteratorAt(mchIndex2);

        VarTrack fgValuesMuon1, fgValuesMuonPV1;
        VarTrack fgValuesMuon2, fgValuesMuonPV2;
        TrackRealigned mchrealigned1, mchrealigned2;
        VarClusters fgValuesCls1, fgValuesCls2;
        if (!FillClusters(muonTrack1, clusters, fgValuesCls1, mchrealigned1) || !FillClusters(muonTrack2, clusters, fgValuesCls2, mchrealigned2)) {
          continue; // Refit is not valid
        }

        if (configRealign.fDoRealign) {

          FillTrack(mchrealigned1, fgValuesMuon1);
          FillTrack(mchrealigned2, fgValuesMuon2);

          // Propagate MCH to PV
          FillPropagation(mchrealigned1, collision1, fgValuesMuonPV1);
          FillPropagation(mchrealigned2, collision2, fgValuesMuonPV2);

          // Recalculate pDCA and Rabs values
          FillPropagation(mchrealigned1, collision1, fgValuesMuon1, kToAbsEnd);
          FillPropagation(mchrealigned1, collision1, fgValuesMuon1, kToDCA);
          FillPropagation(mchrealigned2, collision2, fgValuesMuon2, kToAbsEnd);
          FillPropagation(mchrealigned2, collision2, fgValuesMuon2, kToDCA);
        } else {
          FillTrack<1>(muonTrack1, fgValuesMuon1);
          FillTrack<1>(muonTrack2, fgValuesMuon2);

          // Propagate MCH to PV
          FillPropagation<1>(muonTrack1, collision1, fgValuesMuon1, fgValuesMuonPV1);
          FillPropagation<1>(muonTrack2, collision2, fgValuesMuon2, fgValuesMuonPV2);
          // Calculate DCA
          FillPropagation<1>(muonTrack1, collision1, fgValuesMuon1, fgValuesMuonPV1, kToDCA);
          FillPropagation<1>(muonTrack2, collision2, fgValuesMuon2, fgValuesMuonPV2, kToDCA);
        }

        int sign1 = muonTrack1.sign();
        int sign2 = muonTrack2.sign();

        // only consider same-sign pairs
        if ((sign1 * sign2) <= 0)
          continue;

        bool isPP = false;
        bool isMM = false;
        if (sign1 > 0 && sign2 > 0)
          isPP = true;
        else
          isMM = true;

        int Quadrant1 = GetQuadrantPhi(muonTrack1.phi() * 180.0 / TMath::Pi());
        int Quadrant2 = GetQuadrantPhi(muonTrack2.phi() * 180.0 / TMath::Pi());
        int TopBottom1 = (Quadrant1 == 0 || Quadrant1 == 1) ? 0 : 1;
        int TopBottom2 = (Quadrant2 == 0 || Quadrant2 == 1) ? 0 : 1;

        bool goodMuonTracks = (IsGoodMuon(fgValuesMuon1, fgValuesMuonPV1) && IsGoodMuon(fgValuesMuon2, fgValuesMuonPV2));

        // dimuon variables
        double mass = GetMuMuInvariantMass(fgValuesMuonPV1, fgValuesMuonPV2);
        double pT = GetMuMuPt(fgValuesMuonPV1, fgValuesMuonPV2);
        double dcaXPair = 0.0;
        double dcaYPair = 0.0;
        if (TopBottom1 != TopBottom2) {             // only mixed pairs
          if (TopBottom1 == 0 && TopBottom2 == 1) { // muon1 = top, muon2 = bottom
            if (configRealign.fDoRealign) {
              dcaXPair = fgValuesMuon1.dcaX - fgValuesMuon2.dcaX;
              dcaYPair = fgValuesMuon1.dcaY - fgValuesMuon2.dcaY;
            } else {
              dcaXPair = fgValuesMuonPV1.dcaX - fgValuesMuonPV2.dcaX;
              dcaYPair = fgValuesMuonPV1.dcaY - fgValuesMuonPV2.dcaY;
            }
          } else if (TopBottom1 == 1 && TopBottom2 == 0) { // muon2 = top, muon1 = bottom
            if (configRealign.fDoRealign) {
              dcaXPair = fgValuesMuon2.dcaX - fgValuesMuon1.dcaX;
              dcaYPair = fgValuesMuon2.dcaY - fgValuesMuon1.dcaY;
            } else {
              dcaXPair = fgValuesMuonPV2.dcaX - fgValuesMuonPV1.dcaX;
              dcaYPair = fgValuesMuonPV2.dcaY - fgValuesMuonPV1.dcaY;
            }
          } else if (configRealign.fDoRealign) { // no redefinition necessary if both on same half
            dcaXPair = fgValuesMuon1.dcaX - fgValuesMuon2.dcaX;
            dcaYPair = fgValuesMuon1.dcaY - fgValuesMuon2.dcaY;
          } else {
            dcaXPair = fgValuesMuonPV1.dcaX - fgValuesMuonPV2.dcaX;
            dcaYPair = fgValuesMuonPV1.dcaY - fgValuesMuonPV2.dcaY;
          }
        }
        if (mass < fDimuonDCAMassLow || mass > fDimuonDCAMassHigh)
          continue;

        if (goodMuonTracks) {
          // dimuon DCA
          if (isPP) {
            registryDimuon.get<TH2>(HIST("dimuon/same-event/same-sign-PP/DCA/pT_Mu1DCAx_minus_Mu2DCAx_MuonKine_MuonCuts"))->Fill(pT, dcaXPair);
            registryDimuon.get<TH2>(HIST("dimuon/same-event/same-sign-PP/DCA/pT_Mu1DCAy_minus_Mu2DCAy_MuonKine_MuonCuts"))->Fill(pT, dcaYPair);
          } else if (isMM) {
            registryDimuon.get<TH2>(HIST("dimuon/same-event/same-sign-MM/DCA/pT_Mu1DCAx_minus_Mu2DCAx_MuonKine_MuonCuts"))->Fill(pT, dcaXPair);
            registryDimuon.get<TH2>(HIST("dimuon/same-event/same-sign-MM/DCA/pT_Mu1DCAy_minus_Mu2DCAy_MuonKine_MuonCuts"))->Fill(pT, dcaYPair);
          }
          // dimuon top-bottom separation
          // TODO: left-right ?
          if (TopBottom1 == 0 && TopBottom2 == 0) {
            if (isPP) {
              registryDimuon.get<TH2>(HIST("dimuon/same-event/same-sign-PP/DCA/pT_Mu1TDCAx_minus_Mu2TDCAx_MuonKine_MuonCuts"))->Fill(pT, dcaXPair);
              registryDimuon.get<TH2>(HIST("dimuon/same-event/same-sign-PP/DCA/pT_Mu1TDCAy_minus_Mu2TDCAy_MuonKine_MuonCuts"))->Fill(pT, dcaYPair);
            } else if (isMM) {
              registryDimuon.get<TH2>(HIST("dimuon/same-event/same-sign-MM/DCA/pT_Mu1TDCAx_minus_Mu2TDCAx_MuonKine_MuonCuts"))->Fill(pT, dcaXPair);
              registryDimuon.get<TH2>(HIST("dimuon/same-event/same-sign-MM/DCA/pT_Mu1TDCAy_minus_Mu2TDCAy_MuonKine_MuonCuts"))->Fill(pT, dcaYPair);
            }
          } else if ((TopBottom1 == 0 && TopBottom2 == 1) || (TopBottom1 == 1 && TopBottom2 == 0)) {
            if (isPP) {
              registryDimuon.get<TH2>(HIST("dimuon/same-event/same-sign-PP/DCA/pT_MuTDCAx_minus_MuBDCAx_MuonKine_MuonCuts"))->Fill(pT, dcaXPair);
              registryDimuon.get<TH2>(HIST("dimuon/same-event/same-sign-PP/DCA/pT_MuTDCAy_minus_MuBDCAy_MuonKine_MuonCuts"))->Fill(pT, dcaYPair);
            } else if (isMM) {
              registryDimuon.get<TH2>(HIST("dimuon/same-event/same-sign-MM/DCA/pT_MuTDCAx_minus_MuBDCAx_MuonKine_MuonCuts"))->Fill(pT, dcaXPair);
              registryDimuon.get<TH2>(HIST("dimuon/same-event/same-sign-MM/DCA/pT_MuTDCAy_minus_MuBDCAy_MuonKine_MuonCuts"))->Fill(pT, dcaYPair);
            }
          } else if (TopBottom1 == 1 && TopBottom2 == 1) {
            if (isPP) {
              registryDimuon.get<TH2>(HIST("dimuon/same-event/same-sign-PP/DCA/pT_Mu1BDCAx_minus_Mu2BDCAx_MuonKine_MuonCuts"))->Fill(pT, dcaXPair);
              registryDimuon.get<TH2>(HIST("dimuon/same-event/same-sign-PP/DCA/pT_Mu1BDCAy_minus_Mu2BDCAy_MuonKine_MuonCuts"))->Fill(pT, dcaYPair);
            } else if (isMM) {
              registryDimuon.get<TH2>(HIST("dimuon/same-event/same-sign-MM/DCA/pT_Mu1BDCAx_minus_Mu2BDCAx_MuonKine_MuonCuts"))->Fill(pT, dcaXPair);
              registryDimuon.get<TH2>(HIST("dimuon/same-event/same-sign-MM/DCA/pT_Mu1BDCAy_minus_Mu2BDCAy_MuonKine_MuonCuts"))->Fill(pT, dcaYPair);
            }
          }
        }
      }
    }

    for (const auto& [muon1, muon2] : globalMuonPairs) {
      auto collisionIndex1 = muon1.first;
      auto collisionIndex2 = muon2.first;
      auto& globalTracksVector1 = muon1.second;
      auto& globalTracksVector2 = muon2.second;

      auto const& collision1 = collisions.at(collisionIndex1);
      auto const& collision2 = collisions.at(collisionIndex2);

      auto const& muonTrack1 = muonTracks.rawIteratorAt(globalTracksVector1[0]);
      auto const& muonTrack2 = muonTracks.rawIteratorAt(globalTracksVector2[0]);
      auto const& mftTrack1 = muonTrack1.template matchMFTTrack_as<MyMFTs>();
      auto const& mftTrack2 = muonTrack2.template matchMFTTrack_as<MyMFTs>();
      auto const& mchTrack1 = muonTrack1.template matchMCHTrack_as<MyMuonsWithCov>();
      auto const& mchTrack2 = muonTrack2.template matchMCHTrack_as<MyMuonsWithCov>();

      VarTrack fgValuesMuon1, fgValuesMuonPV1, fgValuesMCH1, fgValuesMCHpv1, fgValuesMFT1, fgValuesMFTpv1;
      VarTrack fgValuesMuon2, fgValuesMuonPV2, fgValuesMCH2, fgValuesMCHpv2, fgValuesMFT2, fgValuesMFTpv2;

      // Fill MCH and MFT tracks
      FillTrack<1>(mchTrack1, fgValuesMCH1);
      FillTrack<0>(mftTrack1, fgValuesMFT1);
      FillTrack<1>(muonTrack1, fgValuesMuon1);
      FillTrack<1>(mchTrack2, fgValuesMCH2);
      FillTrack<0>(mftTrack2, fgValuesMFT2);
      FillTrack<1>(muonTrack2, fgValuesMuon2);

      //// Fill matching chi2
      FillMatching(muonTrack1, fgValuesMCH1, fgValuesMFT1);
      FillMatching(muonTrack2, fgValuesMCH2, fgValuesMFT2);

      TrackRealigned mchrealigned1, mchrealigned2;
      VarClusters fgValuesCls1, fgValuesCls2;
      if (!FillClusters(mchTrack1, clusters, fgValuesCls1, mchrealigned1) || !FillClusters(mchTrack2, clusters, fgValuesCls2, mchrealigned2)) {
        continue; // Refit is not valid
      }

      if (configRealign.fDoRealign) {

        FillTrack(mchrealigned1, fgValuesMCH1);
        FillTrack(mchrealigned2, fgValuesMCH2);

        // Propagate MCH to PV
        FillPropagation(mchrealigned1, collision1, fgValuesMCHpv1);
        FillPropagation(mchrealigned2, collision2, fgValuesMCHpv2);

        // Recalculate pDCA and Rabs values
        FillPropagation(mchrealigned1, collision1, fgValuesMCH1, kToAbsEnd);
        FillPropagation(mchrealigned1, collision1, fgValuesMCH1, kToDCA);
        FillPropagation(mchrealigned2, collision2, fgValuesMCH2, kToAbsEnd);
        FillPropagation(mchrealigned2, collision2, fgValuesMCH2, kToDCA);

      } else {
        FillTrack<1>(mchTrack1, fgValuesMCH1);
        FillTrack<1>(mchTrack2, fgValuesMCH2);

        // Propagate MCH to PV
        FillPropagation<1>(mchTrack1, collision1, fgValuesMCH1, fgValuesMCHpv1);
        FillPropagation<1>(mchTrack2, collision2, fgValuesMCH2, fgValuesMCHpv2);
      }

      // Propagate global muon tracks to PV
      FillPropagation<0>(muonTrack1, collision1, fgValuesMCH1, fgValuesMuonPV1);
      FillPropagation<0>(muonTrack2, collision2, fgValuesMCH2, fgValuesMuonPV2);

      // Propagate MFT tracks to PV
      FillPropagation<0, 1>(mftTrack1, collision1, fgValuesMCHpv1, fgValuesMFTpv1);
      FillPropagation<0, 1>(mftTrack2, collision2, fgValuesMCHpv2, fgValuesMFTpv2);

      int sign1 = mchTrack1.sign();
      int sign2 = mchTrack2.sign();

      // only consider opposite-sign pairs
      if ((sign1 * sign2) >= 0)
        continue;

      // OLD definition using MFT halves
      // indexes indicating whether the positive and negative tracks come from the top or bottom halves of MFT
      // int posTopBottom = (sign1 > 0) ? ((muonTrack1.y() >= 0) ? 0 : 1) : ((muonTrack2.y() >= 0) ? 0 : 1);
      // int negTopBottom = (sign1 < 0) ? ((muonTrack1.y() >= 0) ? 0 : 1) : ((muonTrack2.y() >= 0) ? 0 : 1);

      // NEW definition using MCH tracks, as is done for MUON (MCH-MID) tracks above
      int Quadrant1 = GetQuadrantPhi(muonTrack1.phi() * 180.0 / TMath::Pi());
      int Quadrant2 = GetQuadrantPhi(muonTrack2.phi() * 180.0 / TMath::Pi());
      int TopBottom1 = (Quadrant1 == 0 || Quadrant1 == 1) ? 0 : 1;
      int TopBottom2 = (Quadrant2 == 0 || Quadrant2 == 1) ? 0 : 1;
      // int LeftRight1 = (Quadrant1 == 0 || Quadrant1 == 3) ? 0 : 1;
      // int LeftRight2 = (Quadrant2 == 0 || Quadrant2 == 3) ? 0 : 1;

      bool goodGlobalMuonTracks = (IsGoodGlobalMuon(fgValuesMCH1, fgValuesMCHpv1) && IsGoodGlobalMuon(fgValuesMCH2, fgValuesMCHpv2));
      bool goodGlobalMuonMatches = (IsGoodGlobalMatching(fgValuesMFT1) && IsGoodGlobalMatching(fgValuesMFT2));

      bool sameEvent = (collisionIndex1 == collisionIndex2);

      if (goodGlobalMuonTracks && goodGlobalMuonMatches) {

        double massMCH = GetMuMuInvariantMass(fgValuesMCHpv1, fgValuesMCHpv2);
        double pTmch = GetMuMuPt(fgValuesMCHpv1, fgValuesMCHpv2);
        double mass = GetMuMuInvariantMass(fgValuesMuonPV1, fgValuesMuonPV2);
        // double pT = GetMuMuPt(fgValuesMuonPV1, fgValuesMuonPV2);
        double massScaled = GetMuMuInvariantMass(fgValuesMFTpv1, fgValuesMFTpv2);
        // double pTscaled = GetMuMuPt(fgValuesMFTpv1, fgValuesMFTpv2);

        if (sameEvent) {
          // same-event case
          registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMass_MuonKine_GlobalMatchesCuts"))->Fill(massMCH);
          registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMassFull_MuonKine_GlobalMatchesCuts"))->Fill(massMCH);
          registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMass_GlobalMuonKine_GlobalMatchesCuts"))->Fill(mass);
          registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMassFull_GlobalMuonKine_GlobalMatchesCuts"))->Fill(mass);
          registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMass_ScaledMftKine_GlobalMatchesCuts"))->Fill(massScaled);
          registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMassFull_ScaledMftKine_GlobalMatchesCuts"))->Fill(massScaled);
          registryDimuon.get<TH2>(HIST("dimuon/same-event/invariantMass_pT_MuonKine_GlobalMatchesCuts"))->Fill(massMCH, pTmch);
          registryDimuon.get<TH2>(HIST("dimuon/same-event/invariantMass_pT_GlobalMuonKine_GlobalMatchesCuts"))->Fill(mass, pTmch);
          registryDimuon.get<TH2>(HIST("dimuon/same-event/invariantMass_pT_ScaledMftKine_GlobalMatchesCuts"))->Fill(massScaled, pTmch);

          if (TopBottom1 == 0 && TopBottom2 == 0) {
            registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMass_MuonKine_GlobalMatchesCuts_TT"))->Fill(massMCH);
            registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMass_GlobalMuonKine_GlobalMatchesCuts_TT"))->Fill(mass);
            registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMass_ScaledMftKine_GlobalMatchesCuts_TT"))->Fill(massScaled);
            registryDimuon.get<TH2>(HIST("dimuon/same-event/invariantMass_pT_MuonKine_GlobalMatchesCuts_TT"))->Fill(massMCH, pTmch);
            registryDimuon.get<TH2>(HIST("dimuon/same-event/invariantMass_pT_GlobalMuonKine_GlobalMatchesCuts_TT"))->Fill(mass, pTmch);
            registryDimuon.get<TH2>(HIST("dimuon/same-event/invariantMass_pT_ScaledMftKine_GlobalMatchesCuts_TT"))->Fill(massScaled, pTmch);
          } else if ((TopBottom1 == 0 && TopBottom2 == 1) || (TopBottom1 == 1 && TopBottom2 == 0)) {
            registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMass_MuonKine_GlobalMatchesCuts_TB"))->Fill(massMCH);
            registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMass_GlobalMuonKine_GlobalMatchesCuts_TB"))->Fill(mass);
            registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMass_ScaledMftKine_GlobalMatchesCuts_TB"))->Fill(massScaled);
            registryDimuon.get<TH2>(HIST("dimuon/same-event/invariantMass_pT_MuonKine_GlobalMatchesCuts_TB"))->Fill(massMCH, pTmch);
            registryDimuon.get<TH2>(HIST("dimuon/same-event/invariantMass_pT_GlobalMuonKine_GlobalMatchesCuts_TB"))->Fill(mass, pTmch);
            registryDimuon.get<TH2>(HIST("dimuon/same-event/invariantMass_pT_ScaledMftKine_GlobalMatchesCuts_TB"))->Fill(massScaled, pTmch);
            if (TopBottom1 == 0 && TopBottom2 == 1) {
              if (sign1 > 0) {
                registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMass_MuonKine_GlobalMatchesCuts_TPBN"))->Fill(massMCH);
                registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMass_GlobalMuonKine_GlobalMatchesCuts_TPBN"))->Fill(mass);
                registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMass_ScaledMftKine_GlobalMatchesCuts_TPBN"))->Fill(massScaled);
                registryDimuon.get<TH2>(HIST("dimuon/same-event/invariantMass_pT_MuonKine_GlobalMatchesCuts_TPBN"))->Fill(massMCH, pTmch);
                registryDimuon.get<TH2>(HIST("dimuon/same-event/invariantMass_pT_GlobalMuonKine_GlobalMatchesCuts_TPBN"))->Fill(mass, pTmch);
                registryDimuon.get<TH2>(HIST("dimuon/same-event/invariantMass_pT_ScaledMftKine_GlobalMatchesCuts_TPBN"))->Fill(massScaled, pTmch);
              } else {
                registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMass_MuonKine_GlobalMatchesCuts_TNBP"))->Fill(massMCH);
                registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMass_GlobalMuonKine_GlobalMatchesCuts_TNBP"))->Fill(mass);
                registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMass_ScaledMftKine_GlobalMatchesCuts_TNBP"))->Fill(massScaled);
                registryDimuon.get<TH2>(HIST("dimuon/same-event/invariantMass_pT_MuonKine_GlobalMatchesCuts_TNBP"))->Fill(massMCH, pTmch);
                registryDimuon.get<TH2>(HIST("dimuon/same-event/invariantMass_pT_GlobalMuonKine_GlobalMatchesCuts_TNBP"))->Fill(mass, pTmch);
                registryDimuon.get<TH2>(HIST("dimuon/same-event/invariantMass_pT_ScaledMftKine_GlobalMatchesCuts_TNBP"))->Fill(massScaled, pTmch);
              }
            } else if (TopBottom1 == 0 && TopBottom2 == 1) {
              if (sign2 > 0) {
                registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMass_MuonKine_GlobalMatchesCuts_TPBN"))->Fill(massMCH);
                registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMass_GlobalMuonKine_GlobalMatchesCuts_TPBN"))->Fill(mass);
                registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMass_ScaledMftKine_GlobalMatchesCuts_TPBN"))->Fill(massScaled);
                registryDimuon.get<TH2>(HIST("dimuon/same-event/invariantMass_pT_MuonKine_GlobalMatchesCuts_TPBN"))->Fill(massMCH, pTmch);
                registryDimuon.get<TH2>(HIST("dimuon/same-event/invariantMass_pT_GlobalMuonKine_GlobalMatchesCuts_TPBN"))->Fill(mass, pTmch);
                registryDimuon.get<TH2>(HIST("dimuon/same-event/invariantMass_pT_ScaledMftKine_GlobalMatchesCuts_TPBN"))->Fill(massScaled, pTmch);
              } else {
                registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMass_MuonKine_GlobalMatchesCuts_TNBP"))->Fill(massMCH);
                registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMass_GlobalMuonKine_GlobalMatchesCuts_TNBP"))->Fill(mass);
                registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMass_ScaledMftKine_GlobalMatchesCuts_TNBP"))->Fill(massScaled);
                registryDimuon.get<TH2>(HIST("dimuon/same-event/invariantMass_pT_MuonKine_GlobalMatchesCuts_TNBP"))->Fill(massMCH, pTmch);
                registryDimuon.get<TH2>(HIST("dimuon/same-event/invariantMass_pT_GlobalMuonKine_GlobalMatchesCuts_TNBP"))->Fill(mass, pTmch);
                registryDimuon.get<TH2>(HIST("dimuon/same-event/invariantMass_pT_ScaledMftKine_GlobalMatchesCuts_TNBP"))->Fill(massScaled, pTmch);
              }
            }
          } else if (TopBottom1 == 1 && TopBottom2 == 1) {
            registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMass_MuonKine_GlobalMatchesCuts_BB"))->Fill(massMCH);
            registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMass_GlobalMuonKine_GlobalMatchesCuts_BB"))->Fill(mass);
            registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMass_ScaledMftKine_GlobalMatchesCuts_BB"))->Fill(massScaled);
            registryDimuon.get<TH2>(HIST("dimuon/same-event/invariantMass_pT_MuonKine_GlobalMatchesCuts_BB"))->Fill(massMCH, pTmch);
            registryDimuon.get<TH2>(HIST("dimuon/same-event/invariantMass_pT_GlobalMuonKine_GlobalMatchesCuts_BB"))->Fill(mass, pTmch);
            registryDimuon.get<TH2>(HIST("dimuon/same-event/invariantMass_pT_ScaledMftKine_GlobalMatchesCuts_BB"))->Fill(massScaled, pTmch);
          }

          // mass correlation
          registryDimuon.get<TH2>(HIST("dimuon/same-event/invariantMass_MuonKine_vs_GlobalMuonKine"))->Fill(mass, massMCH);
          registryDimuon.get<TH2>(HIST("dimuon/same-event/invariantMass_ScaledMftKine_vs_GlobalMuonKine"))->Fill(mass, massScaled);
        } else {
          // event-mixing case
          registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMass_MuonKine_GlobalMatchesCuts"))->Fill(massMCH);
          registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMassFull_MuonKine_GlobalMatchesCuts"))->Fill(massMCH);
          registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMass_GlobalMuonKine_GlobalMatchesCuts"))->Fill(mass);
          registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMassFull_GlobalMuonKine_GlobalMatchesCuts"))->Fill(mass);
          registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMass_ScaledMftKine_GlobalMatchesCuts"))->Fill(massScaled);
          registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMassFull_ScaledMftKine_GlobalMatchesCuts"))->Fill(massScaled);
          registryDimuon.get<TH2>(HIST("dimuon/mixed-event/invariantMass_pT_MuonKine_GlobalMatchesCuts"))->Fill(massMCH, pTmch);
          registryDimuon.get<TH2>(HIST("dimuon/mixed-event/invariantMass_pT_GlobalMuonKine_GlobalMatchesCuts"))->Fill(mass, pTmch);
          registryDimuon.get<TH2>(HIST("dimuon/mixed-event/invariantMass_pT_ScaledMftKine_GlobalMatchesCuts"))->Fill(massScaled, pTmch);

          if (TopBottom1 == 0 && TopBottom2 == 0) {
            registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMass_MuonKine_GlobalMatchesCuts_TT"))->Fill(massMCH);
            registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMass_GlobalMuonKine_GlobalMatchesCuts_TT"))->Fill(mass);
            registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMass_ScaledMftKine_GlobalMatchesCuts_TT"))->Fill(massScaled);
            registryDimuon.get<TH2>(HIST("dimuon/mixed-event/invariantMass_pT_MuonKine_GlobalMatchesCuts_TT"))->Fill(massMCH, pTmch);
            registryDimuon.get<TH2>(HIST("dimuon/mixed-event/invariantMass_pT_GlobalMuonKine_GlobalMatchesCuts_TT"))->Fill(mass, pTmch);
            registryDimuon.get<TH2>(HIST("dimuon/mixed-event/invariantMass_pT_ScaledMftKine_GlobalMatchesCuts_TT"))->Fill(massScaled, pTmch);
          } else if ((TopBottom1 == 0 && TopBottom2 == 1) || (TopBottom1 == 1 && TopBottom2 == 0)) {
            registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMass_MuonKine_GlobalMatchesCuts_TB"))->Fill(massMCH);
            registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMass_GlobalMuonKine_GlobalMatchesCuts_TB"))->Fill(mass);
            registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMass_ScaledMftKine_GlobalMatchesCuts_TB"))->Fill(massScaled);
            registryDimuon.get<TH2>(HIST("dimuon/mixed-event/invariantMass_pT_MuonKine_GlobalMatchesCuts_TB"))->Fill(massMCH, pTmch);
            registryDimuon.get<TH2>(HIST("dimuon/mixed-event/invariantMass_pT_GlobalMuonKine_GlobalMatchesCuts_TB"))->Fill(mass, pTmch);
            registryDimuon.get<TH2>(HIST("dimuon/mixed-event/invariantMass_pT_ScaledMftKine_GlobalMatchesCuts_TB"))->Fill(massScaled, pTmch);
            if (TopBottom1 == 0 && TopBottom2 == 1) {
              if (sign1 > 0) {
                registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMass_MuonKine_GlobalMatchesCuts_TPBN"))->Fill(massMCH);
                registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMass_GlobalMuonKine_GlobalMatchesCuts_TPBN"))->Fill(mass);
                registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMass_ScaledMftKine_GlobalMatchesCuts_TPBN"))->Fill(massScaled);
                registryDimuon.get<TH2>(HIST("dimuon/mixed-event/invariantMass_pT_MuonKine_GlobalMatchesCuts_TPBN"))->Fill(massMCH, pTmch);
                registryDimuon.get<TH2>(HIST("dimuon/mixed-event/invariantMass_pT_GlobalMuonKine_GlobalMatchesCuts_TPBN"))->Fill(mass, pTmch);
                registryDimuon.get<TH2>(HIST("dimuon/mixed-event/invariantMass_pT_ScaledMftKine_GlobalMatchesCuts_TPBN"))->Fill(massScaled, pTmch);
              } else {
                registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMass_MuonKine_GlobalMatchesCuts_TNBP"))->Fill(massMCH);
                registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMass_GlobalMuonKine_GlobalMatchesCuts_TNBP"))->Fill(mass);
                registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMass_ScaledMftKine_GlobalMatchesCuts_TNBP"))->Fill(massScaled);
                registryDimuon.get<TH2>(HIST("dimuon/mixed-event/invariantMass_pT_MuonKine_GlobalMatchesCuts_TNBP"))->Fill(massMCH, pTmch);
                registryDimuon.get<TH2>(HIST("dimuon/mixed-event/invariantMass_pT_GlobalMuonKine_GlobalMatchesCuts_TNBP"))->Fill(mass, pTmch);
                registryDimuon.get<TH2>(HIST("dimuon/mixed-event/invariantMass_pT_ScaledMftKine_GlobalMatchesCuts_TNBP"))->Fill(massScaled, pTmch);
              }
            } else if (TopBottom1 == 0 && TopBottom2 == 1) {
              if (sign2 > 0) {
                registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMass_MuonKine_GlobalMatchesCuts_TPBN"))->Fill(massMCH);
                registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMass_GlobalMuonKine_GlobalMatchesCuts_TPBN"))->Fill(mass);
                registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMass_ScaledMftKine_GlobalMatchesCuts_TPBN"))->Fill(massScaled);
                registryDimuon.get<TH2>(HIST("dimuon/mixed-event/invariantMass_pT_MuonKine_GlobalMatchesCuts_TPBN"))->Fill(massMCH, pTmch);
                registryDimuon.get<TH2>(HIST("dimuon/mixed-event/invariantMass_pT_GlobalMuonKine_GlobalMatchesCuts_TPBN"))->Fill(mass, pTmch);
                registryDimuon.get<TH2>(HIST("dimuon/mixed-event/invariantMass_pT_ScaledMftKine_GlobalMatchesCuts_TPBN"))->Fill(massScaled, pTmch);
              } else {
                registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMass_MuonKine_GlobalMatchesCuts_TNBP"))->Fill(massMCH);
                registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMass_GlobalMuonKine_GlobalMatchesCuts_TNBP"))->Fill(mass);
                registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMass_ScaledMftKine_GlobalMatchesCuts_TNBP"))->Fill(massScaled);
                registryDimuon.get<TH2>(HIST("dimuon/mixed-event/invariantMass_pT_MuonKine_GlobalMatchesCuts_TNBP"))->Fill(massMCH, pTmch);
                registryDimuon.get<TH2>(HIST("dimuon/mixed-event/invariantMass_pT_GlobalMuonKine_GlobalMatchesCuts_TNBP"))->Fill(mass, pTmch);
                registryDimuon.get<TH2>(HIST("dimuon/mixed-event/invariantMass_pT_ScaledMftKine_GlobalMatchesCuts_TNBP"))->Fill(massScaled, pTmch);
              }
            }
          } else if (TopBottom1 == 1 && TopBottom2 == 1) {
            registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMass_MuonKine_GlobalMatchesCuts_BB"))->Fill(massMCH);
            registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMass_GlobalMuonKine_GlobalMatchesCuts_BB"))->Fill(mass);
            registryDimuon.get<TH1>(HIST("dimuon/mixed-event/invariantMass_ScaledMftKine_GlobalMatchesCuts_BB"))->Fill(massScaled);
            registryDimuon.get<TH2>(HIST("dimuon/mixed-event/invariantMass_pT_MuonKine_GlobalMatchesCuts_BB"))->Fill(massMCH, pTmch);
            registryDimuon.get<TH2>(HIST("dimuon/mixed-event/invariantMass_pT_GlobalMuonKine_GlobalMatchesCuts_BB"))->Fill(mass, pTmch);
            registryDimuon.get<TH2>(HIST("dimuon/mixed-event/invariantMass_pT_ScaledMftKine_GlobalMatchesCuts_BB"))->Fill(massScaled, pTmch);
          }
        }
      }

      // plots for sub-leading matches are only filled in the same-event case
      if (sameEvent) {
        if (globalTracksVector1.size() > 1) {
          VarTrack fgValuesMuonb1, fgValuesMuonbpv1, fgValuesMFTb1;
          auto const& muonTrack1b = muonTracks.rawIteratorAt(globalTracksVector1[1]);
          FillTrack<1>(muonTrack1b, fgValuesMuonb1);
          FillPropagation<0>(muonTrack1b, collision1, fgValuesMCH1, fgValuesMuonbpv1);

          goodGlobalMuonTracks = (IsGoodGlobalMuon(fgValuesMCH1, fgValuesMCHpv1) && IsGoodGlobalMuon(fgValuesMCH2, fgValuesMCHpv2));
          goodGlobalMuonMatches = (IsGoodGlobalMatching(fgValuesMFTb1) && IsGoodGlobalMatching(fgValuesMFT2));
          double mass = GetMuMuInvariantMass(fgValuesMuonbpv1, fgValuesMuonPV2);
          if (goodGlobalMuonTracks && goodGlobalMuonMatches) {
            registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMass_GlobalMuonKine_GlobalMatchesCuts_subleading_leading"))->Fill(mass);
            registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMassFull_GlobalMuonKine_GlobalMatchesCuts_subleading_leading"))->Fill(mass);
          }
        }

        if (globalTracksVector2.size() > 1) {
          VarTrack fgValuesMuonb2, fgValuesMuonbpv2, fgValuesMFTb2;
          auto const& muonTrack2b = muonTracks.rawIteratorAt(globalTracksVector2[1]);
          FillTrack<1>(muonTrack2b, fgValuesMuonb2);
          FillPropagation<0>(muonTrack2b, collision2, fgValuesMCH2, fgValuesMuonbpv2);

          goodGlobalMuonTracks = (IsGoodGlobalMuon(fgValuesMCH1, fgValuesMCHpv1) && IsGoodGlobalMuon(fgValuesMCH2, fgValuesMCHpv2));
          goodGlobalMuonMatches = (IsGoodGlobalMatching(fgValuesMFTb2) && IsGoodGlobalMatching(fgValuesMFT1));
          double mass = GetMuMuInvariantMass(fgValuesMuonbpv2, fgValuesMuonPV1);
          if (goodGlobalMuonTracks && goodGlobalMuonMatches) {
            registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMass_GlobalMuonKine_GlobalMatchesCuts_leading_subleading"))->Fill(mass);
            registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMassFull_GlobalMuonKine_GlobalMatchesCuts_leading_subleading"))->Fill(mass);
          }
        }

        if (globalTracksVector1.size() > 1 && globalTracksVector2.size() > 1) {
          VarTrack fgValuesMuonb1, fgValuesMuonbpv1, fgValuesMFTb1;
          VarTrack fgValuesMuonb2, fgValuesMuonbpv2, fgValuesMFTb2;
          auto const& muonTrack1b = muonTracks.rawIteratorAt(globalTracksVector1[1]);
          auto const& muonTrack2b = muonTracks.rawIteratorAt(globalTracksVector2[1]);

          FillTrack<1>(muonTrack1b, fgValuesMuonb1);
          FillPropagation<0>(muonTrack1b, collision1, fgValuesMCH1, fgValuesMuonbpv1);

          FillTrack<1>(muonTrack2b, fgValuesMuonb2);
          FillPropagation<0>(muonTrack2b, collision2, fgValuesMCH2, fgValuesMuonbpv2);

          goodGlobalMuonTracks = (IsGoodGlobalMuon(fgValuesMCH1, fgValuesMCHpv1) && IsGoodGlobalMuon(fgValuesMCH2, fgValuesMCHpv2));
          goodGlobalMuonMatches = (IsGoodGlobalMatching(fgValuesMFTb1) && IsGoodGlobalMatching(fgValuesMFTb2));
          double mass = GetMuMuInvariantMass(fgValuesMuonbpv1, fgValuesMuonbpv2);
          double massLeading = GetMuMuInvariantMass(fgValuesMuonPV1, fgValuesMuonPV2);
          if (goodGlobalMuonTracks && goodGlobalMuonMatches) {
            registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMass_GlobalMuonKine_GlobalMatchesCuts_subleading_subleading"))->Fill(mass);
            registryDimuon.get<TH1>(HIST("dimuon/same-event/invariantMassFull_GlobalMuonKine_GlobalMatchesCuts_subleading_subleading"))->Fill(mass);

            // mass correlation
            registryDimuon.get<TH2>(HIST("dimuon/same-event/invariantMass_GlobalMuonKine_subleading_vs_leading"))->Fill(massLeading, mass);
          }
        }
      }
    }
  }

  void processMuonQa(MyEvents const& collisions, aod::BCsWithTimestamps const& bcs, MyMuonsWithCov const& muontracks, MyMFTs const& mfttracks, aod::FwdTrkCls const& muonclusters)
  {
    std::map<uint64_t, VarColl> collisionSel;
    std::map<uint64_t, std::vector<uint64_t>> matchingCandidates;

    initCCDB(bcs);

    runEventSelection(collisions, bcs, muontracks, mfttracks, collisionSel);

    runMuonQA(collisionSel, matchingCandidates, muontracks, mfttracks, muonclusters);

    if (configQAs.fEnableQADimuon) {
      runDimuonQA(collisionSel, matchingCandidates, muontracks, muonclusters);
    }
  }
  PROCESS_SWITCH(muonQa, processMuonQa, "Process to run muon QA", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<muonQa>(cfgc)};
}
