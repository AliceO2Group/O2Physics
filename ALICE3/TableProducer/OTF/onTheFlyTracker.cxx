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

/// \file onTheFlyTracker.cxx
///
/// \brief LUT-based on-the-fly analysis task-level tracking
///
/// This task allows for the calculation of aod::collisions and aod::Tracks in a synthetic manner,
/// smearing MC particles with very configurable settings. This will allow for the usage of
/// custom LUTs (obtained through separate studies) and the subsequent estimate of the performance
/// of a future detector even in very statistics-hungry analyses.
///
/// \author David Dobrigkeit Chinellato, UNICAMP
/// \author Nicolò Jacazio <nicolo.jacazio@cern.ch>, UniBo
/// \author Roberto Preghenella preghenella@bo.infn.it
///

#include "GeometryContainer.h"

#include "ALICE3/Core/DetLayer.h"
#include "ALICE3/Core/FastTracker.h"
#include "ALICE3/Core/FlatTrackSmearer.h"
#include "ALICE3/Core/TrackUtilities.h"
#include "ALICE3/DataModel/OTFCollision.h"
#include "ALICE3/DataModel/OTFMCParticle.h"
#include "ALICE3/DataModel/OTFStrangeness.h"
#include "ALICE3/DataModel/collisionAlice3.h"
#include "ALICE3/DataModel/tracksAlice3.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <CommonDataFormat/InteractionRecord.h>
#include <CommonDataFormat/TimeStamp.h>
#include <CommonUtils/ConfigurableParam.h>
#include <DCAFitter/DCAFitterN.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DetectorsBase/MatLayerCylSet.h>
#include <DetectorsBase/Propagator.h>
#include <DetectorsVertexing/PVertexer.h>
#include <DetectorsVertexing/PVertexerHelpers.h>
#include <DetectorsVertexing/PVertexerParams.h>
#include <Field/MagneticField.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/DataTypes.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>
#include <MathUtils/Primitive2D.h>
#include <MathUtils/Utils.h>
#include <ReconstructionDataFormats/DCA.h>
#include <ReconstructionDataFormats/GlobalTrackID.h>
#include <ReconstructionDataFormats/PID.h>
#include <ReconstructionDataFormats/PrimaryVertex.h>
#include <ReconstructionDataFormats/Track.h>
#include <ReconstructionDataFormats/TrackParametrization.h>
#include <SimulationDataFormat/InteractionSampler.h>
#include <SimulationDataFormat/MCCompLabel.h>
#include <SimulationDataFormat/MCEventLabel.h>

#include <TGenPhaseSpace.h>
#include <TGeoGlobalMagField.h>
#include <TH1.h>
#include <TH2.h>
#include <TLorentzVector.h>
#include <TMCProcess.h>
#include <TMath.h>
#include <TPDGCode.h>
#include <TRandom3.h>
#include <TString.h>

#include <sys/types.h>

#include <RtypesCore.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <map>
#include <memory>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <math.h>

using namespace o2;
using namespace o2::framework;
using std::array;
#define getHist(type, name) \
  std::get<std::shared_ptr<type>>(histPointers[name])
#define fillHist(type, name, ...)                                           \
  if (histPointers.find(name) != histPointers.end())                        \
    std::get<std::shared_ptr<type>>(histPointers[name])->Fill(__VA_ARGS__); \
  else                                                                      \
    LOG(fatal) << "Histogram " << name << " not found!";
#define insertHist(name, ...) histPointers[name] = histos.add((name).c_str(), __VA_ARGS__);

enum TrackType {
  kNone = 0,
  kRecoPrimary,
  kGhostPrimary,
  kRecoV0Daug,
  kGhostV0Daug,
  kRecoCasc,
  kGenCasc,
  kRecoCascDaug,
  kGenCascDaug,
  kNTrackTypes,
};

struct OnTheFlyTracker {
  Produces<aod::Collisions> tableCollisions;
  Produces<aod::McCollisionLabels> tableMcCollisionLabels;
  Produces<aod::StoredTracks> tableStoredTracks;
  Produces<aod::TracksExtension> tableTracksExtension;
  Produces<aod::StoredTracksCov> tableStoredTracksCov;
  Produces<aod::TracksCovExtension> tableTracksCovExtension;
  Produces<aod::McTrackLabels> tableMcTrackLabels;
  Produces<aod::McTrackWithDauLabels> tableMcTrackWithDauLabels;
  Produces<aod::TracksDCA> tableTracksDCA;
  Produces<aod::TracksDCACov> tableTracksDCACov;
  Produces<aod::CollisionsAlice3> tableCollisionsAlice3;
  Produces<aod::TracksAlice3> tableTracksAlice3;
  Produces<aod::TracksExtraA3> tableTracksExtraA3;
  Produces<aod::UpgradeCascades> tableUpgradeCascades;
  Produces<aod::OTFLUTConfigId> tableOTFLUTConfigId;
  Produces<aod::UpgradeV0s> tableUpgradeV0s;

  // optionally produced, empty (to be tuned later)
  Produces<aod::StoredTracksExtra_002> tableStoredTracksExtra; // base table, extend later
  Produces<aod::TrackSelection> tableTrackSelection;
  Produces<aod::TrackSelectionExtension> tableTrackSelectionExtension;

  Configurable<int> seed{"seed", 0, "TGenPhaseSpace seed"};
  Configurable<float> maxEta{"maxEta", 1.5, "maximum eta to consider viable"};
  Configurable<float> multEtaRange{"multEtaRange", 0.8, "eta range to compute the multiplicity"};
  Configurable<float> minPt{"minPt", 0.1, "minimum pt to consider viable"};
  Configurable<bool> enablePrimarySmearing{"enablePrimarySmearing", false, "Enable smearing of primary particles"};
  Configurable<bool> enableSecondarySmearing{"enableSecondarySmearing", false, "Enable smearing of weak decay daughters"};
  Configurable<bool> enableNucleiSmearing{"enableNucleiSmearing", false, "Enable smearing of nuclei"};
  Configurable<bool> enablePrimaryVertexing{"enablePrimaryVertexing", true, "Enable primary vertexing"};
  Configurable<std::string> primaryVertexOption{"primaryVertexOption", "pvertexer.maxChi2TZDebris=10;pvertexer.acceptableScale2=9;pvertexer.minScale2=2;pvertexer.timeMarginVertexTime=1.3;;pvertexer.maxChi2TZDebris=40;pvertexer.maxChi2Mean=12;pvertexer.maxMultRatDebris=1.;pvertexer.addTimeSigma2Debris=1e-2;pvertexer.meanVertexExtraErrSelection=0.03;", "Option for the primary vertexer"};
  Configurable<bool> interpolateLutEfficiencyVsNch{"interpolateLutEfficiencyVsNch", true, "interpolate LUT efficiency as f(Nch)"};

  Configurable<bool> populateTracksDCA{"populateTracksDCA", true, "populate TracksDCA table"};
  Configurable<bool> populateTracksDCACov{"populateTracksDCACov", false, "populate TracksDCACov table"};
  Configurable<bool> populateTracksExtra{"populateTracksExtra", false, "populate TrackExtra table (legacy)"};
  Configurable<bool> populateTrackSelection{"populateTrackSelection", false, "populate TrackSelection table (legacy)"};

  Configurable<bool> processUnreconstructedTracks{"processUnreconstructedTracks", false, "process (smear) unreco-ed tracks"};
  Configurable<bool> doExtraQA{"doExtraQA", false, "do extra 2D QA plots"};
  Configurable<bool> extraQAwithoutDecayDaughters{"extraQAwithoutDecayDaughters", false, "remove decay daughters from qa plots (yes/no)"};

  struct : ConfigurableGroup {
    ConfigurableAxis axisMomentum{"axisMomentum", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "#it{p} (GeV/#it{c})"};
    ConfigurableAxis axisNVertices{"axisNVertices", {20, -0.5, 19.5}, "N_{vertices}"};
    ConfigurableAxis axisDeltaVtxCoord{"axisDeltaVtxCoord", {100, -5., 5}, "Delta Vtx coord (cm)"};
    ConfigurableAxis axisDeltaMultPVRecoGen{"axisDeltaMultPVRecoGen", {51, -25, 25}, "Delta Multiplicity_{PV} (cm)"};
    ConfigurableAxis axisVtxMult{"axisVtxMult", {101, 0, 100}, "Vertex Multiplicity"};
    ConfigurableAxis axisMultiplicity{"axisMultiplicity", {100, -0.5, 99.5}, "N_{contributors}"};
    ConfigurableAxis axisVertexZ{"axisVertexZ", {40, -20, 20}, "vertex Z (cm)"};
    ConfigurableAxis axisDCA{"axisDCA", {400, -200, 200}, "DCA (#mum)"};
    ConfigurableAxis axisX{"axisX", {250, -50, 200}, "track X (cm)"};
    ConfigurableAxis axisDecayRadius{"axisDecayRadius", {55, 0.01, 100}, "decay radius"};
    ConfigurableAxis axisK0Mass{"axisK0Mass", {200, 0.4f, 0.6f}, ""};
    ConfigurableAxis axisLambdaMass{"axisLambdaMass", {200, 1.101f, 1.131f}, ""};
    ConfigurableAxis axisXiMass{"axisXiMass", {200, 1.22f, 1.42f}, ""};

    ConfigurableAxis axisPtRes{"axisPtRes", {200, -0.4f, 0.4f}, "#Delta p_{T} / Reco p_{T}"};
    ConfigurableAxis axisDeltaPt{"axisDeltaPt", {200, -1.0f, +1.0f}, "#Delta p_{T}"};
    ConfigurableAxis axisDeltaEta{"axisDeltaEta", {200, -0.5f, +0.5f}, "#Delta #eta"};

    ConfigurableAxis axisRadius{"axisRadius", {2500, 0.0f, +250.0f}, "R (cm)"};
    ConfigurableAxis axisZ{"axisZ", {100, -250.0f, +250.0f}, "Z (cm)"};
  } axes;

  // for topo var QA
  struct : ConfigurableGroup {
    std::string prefix = "fastTrackerSettings"; // JSON group name
    Configurable<int> minSiliconHits{"minSiliconHits", 6, "minimum number of silicon hits to accept track"};
    Configurable<int> minSiliconHitsForKinkReco{"minSiliconHitsForKinkReco", 4, "minimum number of silicon hits to accept track"};
    Configurable<int> minSiliconHitsIfTPCUsed{"minSiliconHitsIfTPCUsed", 2, "minimum number of silicon hits to accept track in case TPC info is present"};
    Configurable<int> minTPCClusters{"minTPCClusters", 70, "minimum number of TPC hits necessary to consider minSiliconHitsIfTPCUsed"};
    Configurable<bool> applyZacceptance{"applyZacceptance", false, "apply z limits to detector layers or not"};
    Configurable<bool> applyMSCorrection{"applyMSCorrection", true, "apply ms corrections for secondaries or not"};
    Configurable<bool> applyElossCorrection{"applyElossCorrection", true, "apply eloss corrections for secondaries or not"};
  } fastTrackerSettings; // allows for gap between peak and bg in case someone wants to

  struct : ConfigurableGroup {
    std::string prefix = "fastPrimaryTrackerSettings";
    Configurable<bool> fastTrackPrimaries{"fastTrackPrimaries", false, "Use fasttracker for primary tracks. Enable with care"};
    Configurable<int> minSiliconHits{"minSiliconHits", 4, "minimum number of silicon hits to accept track"};
    Configurable<bool> applyZacceptance{"applyZacceptance", false, "apply z limits to detector layers or not"};
    Configurable<bool> applyMSCorrection{"applyMSCorrection", true, "apply ms corrections for secondaries or not"};
    Configurable<bool> applyElossCorrection{"applyElossCorrection", true, "apply eloss corrections for secondaries or not"};
  } fastPrimaryTrackerSettings;

  struct : ConfigurableGroup {
    std::string prefix = "cascadeDecaySettings"; // Cascade decay settings
    Configurable<bool> decayXi{"decayXi", false, "Manually decay Xi and fill tables with daughters"};
    Configurable<bool> findXi{"findXi", false, "if decayXi on, find Xi and fill Tracks table also with Xi"};
    Configurable<bool> trackXi{"trackXi", false, "if findXi on, attempt to track Xi"};
    Configurable<bool> doXiQA{"doXiQA", false, "QA plots for when treating Xi"};
    Configurable<int> minStraTrackHits{"minStraTrackHits", 1, "if trackXi, set min strangeness tracking hits"};
    Configurable<int> doKinkReco{"doKinkReco", 0, "Flag for kink reco setting: 0 - disabled, 1 - complementary, 2 - only"};
  } cascadeDecaySettings;

  struct : ConfigurableGroup {
    std::string prefix = "v0DecaySettings"; // Cascade decay settings
    Configurable<bool> decayV0{"decayV0", false, "Manually decay V0 and fill tables with daughters"};
    Configurable<bool> findV0{"findV0", false, "if decayV0 on, find V0 and fill Tracks table also with Xi"};
    Configurable<bool> doV0QA{"doV0QA", false, "QA plots for when treating V0"};
  } v0DecaySettings;

  struct : ConfigurableGroup {
    std::string prefix = "cfgFitter";
    Configurable<bool> propagateToPCA{"propagateToPCA", false, "create tracks version propagated to PCA"};
    Configurable<double> maxR{"maxR", 200., "reject PCA's above this radius"};
    Configurable<double> minParamChange{"minParamChange", 1.e-3, "stop iterations if largest change of any X is smaller than this"};
    Configurable<double> minRelChi2Change{"minRelChi2Change", 0.9, "stop iterations is chi2/chi2old > this"};
    Configurable<double> maxDZIni{"maxDZIni", 1e9, "reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
    Configurable<double> maxDXYIni{"maxDXYIni", 4, "reject (if>0) PCA candidate if tracks DXY exceeds threshold"};
    Configurable<double> maxVtxChi2{"maxVtxChi2", 1e9, "reject (if>0) vtx. chi2 above this value"};
    Configurable<bool> useAbsDCA{"useAbsDCA", true, "Minimise abs. distance rather than chi2"};
    Configurable<bool> useWeightedFinalPCA{"useWeightedFinalPCA", false, "Recalculate vertex position using track covariances, effective only if useAbsDCA is true"};
  } cfgFitter;

  using PVertex = o2::dataformats::PrimaryVertex;

  // for secondary vertex finding
  o2::vertexing::DCAFitterN<2> fitter;

  // FastTracker machinery
  std::vector<std::unique_ptr<o2::fastsim::FastTracker>> fastTracker;

  // V0 names for filling histograms
  static constexpr int NtypesV0 = 3;
  static constexpr std::string_view NameV0s[NtypesV0] = {"K0", "Lambda", "AntiLambda"};

  // Class to hold the track information for the O2 vertexing
  class TrackAlice3 : public o2::track::TrackParCov
  {
    using TimeEst = o2::dataformats::TimeStampWithError<float, float>;

   public:
    TrackAlice3() = default;
    ~TrackAlice3() = default;
    TrackAlice3(const TrackAlice3& src) = default;
    TrackAlice3(const o2::track::TrackParCov& src, const int64_t label,
                const float time = 0,
                const float timeError = 1,
                bool decayDauInput = false,
                bool weakDecayDauInput = false,
                int isUsedInCascadingInput = 0,
                int nSiliconHitsInput = 0,
                int nTPCHitsInput = 0,
                TrackType trackTypeInput = TrackType::kRecoPrimary) : o2::track::TrackParCov(src),
                                                                      mcLabel{label},
                                                                      timeEst{time, timeError},
                                                                      isDecayDau(decayDauInput),
                                                                      isWeakDecayDau(weakDecayDauInput),
                                                                      isUsedInCascading(isUsedInCascadingInput),
                                                                      nSiliconHits(nSiliconHitsInput),
                                                                      nTPCHits(nTPCHitsInput),
                                                                      trackType(trackTypeInput) {}
    const TimeEst& getTimeMUS() const { return timeEst; }
    int64_t mcLabel;
    TimeEst timeEst; ///< time estimate in ns
    bool isDecayDau;
    bool isWeakDecayDau;
    int isUsedInCascading; // 0: not at all, 1: is a cascade, 2: is a bachelor, 3: is a pion, 4: is a proton
    int nSiliconHits;
    int nTPCHits;
    TrackType trackType;
  };

  // Helper struct to pass cascade information
  struct cascadecandidate {
    int cascadeTrackId; // track index in the Tracks table
    int positiveId;     // track index in the Tracks table
    int negativeId;     // track index in the Tracks table
    int bachelorId;     // track index in the Tracks table

    float pt;
    float eta;
    float dcaV0dau;
    float dcacascdau;
    float v0radius;
    float cascradius;
    float cascradiusMC;

    // for tracking
    int findableClusters;
    int foundClusters;

    float mLambda;
    float mXi;
  };
  cascadecandidate thisCascade;

  // Helper struct to pass V0 information
  struct v0candidate {
    int positiveId;   // track index in the Tracks table
    int negativeId;   // track index in the Tracks table
    int mcParticleId; // mc particle index

    float pt;

    float dcaV0dau;
    float v0radius;

    float mLambda;
    float mAntiLambda;
    float mK0;
  };
  v0candidate thisV0;
  // Constants
  static constexpr int kv0Prongs = 2;
  static constexpr std::array<int, 3> v0PDGs = {kK0Short,
                                                kLambda0,
                                                kLambda0Bar};

  static constexpr std::array<int, 5> longLivedHandledPDGs = {kElectron,
                                                              kMuonMinus,
                                                              kPiPlus,
                                                              kKPlus,
                                                              kProton};

  static constexpr std::array<int, 4> nucleiPDGs = {o2::constants::physics::kDeuteron,
                                                    o2::constants::physics::kTriton,
                                                    o2::constants::physics::kHelium3,
                                                    o2::constants::physics::kAlpha};

  // Primary vertexing QA
  std::pair<float, float> vertexReconstructionEfficiencyCounters = {0, 0}; // {nVerticesWithMoreThan2Contributors, nVerticesReconstructed}

  // necessary for particle charges
  Service<o2::framework::O2DatabasePDG> pdgDB;

  // for handling basic QA histograms if requested
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  std::map<std::string, HistPtr> histPointers;

  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;

  // Track smearer array, one per geometry
  std::vector<std::unique_ptr<o2::delphes::TrackSmearer>> mSmearer;

  // For processing and vertexing
  std::vector<TrackAlice3> recoPrimaries;
  std::vector<TrackAlice3> ghostPrimaries;
  std::vector<TrackAlice3> recoV0Daugs;
  std::vector<TrackAlice3> recoCascDaugs;
  std::vector<TrackAlice3> tracksAlice3;
  std::vector<TrackAlice3> ghostTracksAlice3;
  std::vector<o2::InteractionRecord> bcData;
  o2::steer::InteractionSampler irSampler;
  o2::vertexing::PVertexer vertexer;
  std::vector<cascadecandidate> cascadesAlice3;
  std::vector<v0candidate> v0sAlice3;

  // For TGenPhaseSpace seed
  TRandom3 rand;
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  // Configuration defined at init time
  o2::fastsim::GeometryContainer mGeoContainer;
  float mMagneticField = 0.0f;
  // Time resolution constants
  const float timeResolutionNs = 100.f; // ns
  const float nsToMus = 1e-3f;
  const float timeResolutionUs = timeResolutionNs * nsToMus; // us

  o2::dataformats::DCA dcaInfo;
  o2::dataformats::VertexBase vtx;

  void init(o2::framework::InitContext& initContext)
  {
    LOG(info) << "Initializing OnTheFlyTracker task";
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setTimestamp(-1);
    mGeoContainer.setCcdbManager(ccdb.operator->());
    mGeoContainer.init(initContext);

    const int nGeometries = mGeoContainer.getNumberOfConfigurations();
    mMagneticField = mGeoContainer.getFloatValue(0, "global", "magneticfield");
    for (int icfg = 0; icfg < nGeometries; ++icfg) {
      const std::string histPath = "Configuration_" + std::to_string(icfg) + "/";
      mSmearer.emplace_back(std::make_unique<o2::delphes::TrackSmearer>());
      mSmearer[icfg]->setCcdbManager(ccdb.operator->());
      std::map<std::string, std::string> globalConfiguration = mGeoContainer.getConfiguration(icfg, "global");
      if (enablePrimarySmearing) {
        // load LUTs for primaries
        for (const auto& entry : globalConfiguration) {
          int pdg = 0;
          if (entry.first.find("lut") != 0) {
            continue;
          }
          if (entry.first.find("lutEl") != std::string::npos) {
            pdg = kElectron;
          } else if (entry.first.find("lutMu") != std::string::npos) {
            pdg = kMuonMinus;
          } else if (entry.first.find("lutPi") != std::string::npos) {
            pdg = kPiPlus;
          } else if (entry.first.find("lutKa") != std::string::npos) {
            pdg = kKPlus;
          } else if (entry.first.find("lutPr") != std::string::npos) {
            pdg = kProton;
          } else if (entry.first.find("lutDe") != std::string::npos) {
            pdg = o2::constants::physics::kDeuteron;
          } else if (entry.first.find("lutTr") != std::string::npos) {
            pdg = o2::constants::physics::kTriton;
          } else if (entry.first.find("lutHe3") != std::string::npos) {
            pdg = o2::constants::physics::kHelium3;
          } else if (entry.first.find("lutAl") != std::string::npos) {
            pdg = o2::constants::physics::kAlpha;
          }

          std::string filename = entry.second;
          if (pdg == 0) {
            LOG(fatal) << "Unknown LUT entry " << entry.first << " for global configuration";
          }
          LOG(info) << "Loading LUT for pdg " << pdg << " for config " << icfg << " from provided file '" << filename << "'";
          if (filename.empty()) {
            LOG(warning) << "No LUT file passed for pdg " << pdg << ", skipping.";
          }
          // strip from leading/trailing spaces
          filename.erase(0, filename.find_first_not_of(" "));
          filename.erase(filename.find_last_not_of(" ") + 1);
          if (filename.empty()) {
            LOG(warning) << "No LUT file passed for pdg " << pdg << ", skipping.";
          }
          bool success = mSmearer[icfg]->loadTable(pdg, filename.c_str());
          if (!success) {
            LOG(fatal) << "Having issue with loading the LUT " << pdg << " " << filename;
          }
        }

        // interpolate efficiencies if requested to do so
        mSmearer[icfg]->interpolateEfficiency(interpolateLutEfficiencyVsNch.value);

        // smear un-reco'ed tracks if asked to do so
        mSmearer[icfg]->skipUnreconstructed(!processUnreconstructedTracks.value);

        insertHist(histPath + "hPtGenerated", "hPtGenerated;#it{p}_{T} (GeV/c);Counts", {kTH1D, {{axes.axisMomentum}}});
        insertHist(histPath + "hPhiGenerated", "hPhiGenerated;#phi (rad);Counts", {kTH1D, {{100, 0.0f, 2 * M_PI, "#phi (rad)"}}});

        insertHist(histPath + "hPtGeneratedEl", "hPtGeneratedEl;Gen #it{p}_{T} (GeV/c);Counts", {kTH1D, {{axes.axisMomentum}}});
        insertHist(histPath + "hPtGeneratedPi", "hPtGeneratedPi;Gen #it{p}_{T} (GeV/c);Counts", {kTH1D, {{axes.axisMomentum}}});
        insertHist(histPath + "hPtGeneratedKa", "hPtGeneratedKa;Gen #it{p}_{T} (GeV/c);Counts", {kTH1D, {{axes.axisMomentum}}});
        insertHist(histPath + "hPtGeneratedPr", "hPtGeneratedPr;Gen #it{p}_{T} (GeV/c);Counts", {kTH1D, {{axes.axisMomentum}}});
        insertHist(histPath + "hPtReconstructed", "hPtReconstructed;Reco #it{p}_{T} (GeV/c);Counts", {kTH1D, {{axes.axisMomentum}}});
        insertHist(histPath + "hPtReconstructedEl", "hPtReconstructedEl;Reco #it{p}_{T} (GeV/c);Counts", {kTH1D, {{axes.axisMomentum}}});
        insertHist(histPath + "hPtReconstructedPi", "hPtReconstructedPi;Reco #it{p}_{T} (GeV/c);Counts", {kTH1D, {{axes.axisMomentum}}});
        insertHist(histPath + "hPtReconstructedKa", "hPtReconstructedKa;Reco #it{p}_{T} (GeV/c);Counts", {kTH1D, {{axes.axisMomentum}}});
        insertHist(histPath + "hPtReconstructedPr", "hPtReconstructedPr;Reco #it{p}_{T} (GeV/c);Counts", {kTH1D, {{axes.axisMomentum}}});
      }
      // Collision QA
      insertHist(histPath + "hPVz", "hPVz;Primary Vertex Z (cm);Counts", {kTH1D, {{axes.axisVertexZ}}});
      insertHist(histPath + "hLUTMultiplicity", "hLUTMultiplicity;dN/d#eta;Counts", {kTH1D, {{axes.axisMultiplicity}}});
      insertHist(histPath + "hSimMultiplicity", "hSimMultiplicity;Gen. multiplicity;Counts", {kTH1D, {{axes.axisMultiplicity}}});
      insertHist(histPath + "hRecoMultiplicity", "hRecoMultiplicity;Reco Multiplicity;Counts", {kTH1D, {{axes.axisMultiplicity}}});

      if (enablePrimaryVertexing) {
        insertHist(histPath + "hDeltaXPVRecoGen", "hDeltaXPVRecoGen;Delta X (reco - gen), cm", {kTH2D, {{axes.axisDeltaVtxCoord, axes.axisMultiplicity}}});
        insertHist(histPath + "hDeltaYPVRecoGen", "hDeltaYPVRecoGen;Delta Y (reco - gen), cm", {kTH2D, {{axes.axisDeltaVtxCoord, axes.axisMultiplicity}}});
        insertHist(histPath + "hDeltaZPVRecoGen", "hDeltaZPVRecoGen;Delta Z (reco - gen), cm", {kTH2D, {{axes.axisDeltaVtxCoord, axes.axisMultiplicity}}});
        insertHist(histPath + "hDeltaMultPVRecoGen", "hDeltaMultPVRecoGen;Delta Multiplicity (reco - gen)", {kTH1D, {{axes.axisDeltaMultPVRecoGen}}});
        insertHist(histPath + "hVtxMultGen", "hVtxMultGen;Generated Vertex Multiplicity", {kTH1D, {{axes.axisVtxMult}}});
        insertHist(histPath + "hVtxMultReco", "hVtxMultReco;Reconstructed Vertex Multiplicity", {kTH1D, {{axes.axisVtxMult}}});
        insertHist(histPath + "hVtxTrials", "hVtxTrials;Vertex Reconstruction Trials", {kTH1D, {{2, -0.5, 1.5}}});
        // Set the bin labels
        getHist(TH1, histPath + "hVtxTrials")->GetXaxis()->SetBinLabel(1, "Tried");
        getHist(TH1, histPath + "hVtxTrials")->GetXaxis()->SetBinLabel(2, "Succeeded");
      }

      if (enableSecondarySmearing) {
        fastTracker.emplace_back(std::make_unique<o2::fastsim::FastTracker>());
        fastTracker[icfg]->SetMagneticField(mMagneticField);
        fastTracker[icfg]->SetApplyZacceptance(fastTrackerSettings.applyZacceptance);
        fastTracker[icfg]->SetApplyMSCorrection(fastTrackerSettings.applyMSCorrection);
        fastTracker[icfg]->SetApplyElossCorrection(fastTrackerSettings.applyElossCorrection);
        fastTracker[icfg]->AddGenericDetector(mGeoContainer.getEntry(icfg), ccdb.operator->());
        fastTracker[icfg]->Print(); // print fastTracker settings

        if (cascadeDecaySettings.doXiQA) {
          insertHist(histPath + "hXiBuilding", "hXiBuilding", {kTH1F, {{10, -0.5f, 9.5f}}});
          getHist(TH1, histPath + "hXiBuilding")->GetXaxis()->SetBinLabel(1, "Generated");
          getHist(TH1, histPath + "hXiBuilding")->GetXaxis()->SetBinLabel(2, "Secondary smearing prong 0");
          getHist(TH1, histPath + "hXiBuilding")->GetXaxis()->SetBinLabel(3, "Secondary smearing prong 1");
          getHist(TH1, histPath + "hXiBuilding")->GetXaxis()->SetBinLabel(4, "Secondary smearing prong 2");
          getHist(TH1, histPath + "hXiBuilding")->GetXaxis()->SetBinLabel(5, "Not Nan");
          getHist(TH1, histPath + "hXiBuilding")->GetXaxis()->SetBinLabel(6, "Start Reco");
          getHist(TH1, histPath + "hXiBuilding")->GetXaxis()->SetBinLabel(7, "V0 fitter ok");
          getHist(TH1, histPath + "hXiBuilding")->GetXaxis()->SetBinLabel(8, "Kink fitter ok");

          insertHist(histPath + "hGenXi", "hGenXi;Decay Radius;#it{p}_{T}", {kTH2F, {axes.axisDecayRadius, axes.axisMomentum}});
          insertHist(histPath + "hRecoXi", "hRecoXi;Decay Radius;#it{p}_{T}", {kTH2F, {axes.axisDecayRadius, axes.axisMomentum}});

          insertHist(histPath + "hGenPiFromXi", "hGenPiFromXi;Decay Radius;#it{p}_{T}", {kTH2F, {axes.axisDecayRadius, axes.axisMomentum}});
          insertHist(histPath + "hGenPiFromLa", "hGenPiFromLa;Decay Radius;#it{p}_{T}", {kTH2F, {axes.axisDecayRadius, axes.axisMomentum}});
          insertHist(histPath + "hGenPrFromLa", "hGenPrFromLa;Decay Radius;#it{p}_{T}", {kTH2F, {axes.axisDecayRadius, axes.axisMomentum}});
          insertHist(histPath + "hRecoPiFromXi", "hRecoPiFromXi;Decay Radius;#it{p}_{T}", {kTH2F, {axes.axisDecayRadius, axes.axisMomentum}});
          insertHist(histPath + "hRecoPiFromLa", "hRecoPiFromLa;Decay Radius;#it{p}_{T}", {kTH2F, {axes.axisDecayRadius, axes.axisMomentum}});
          insertHist(histPath + "hRecoPrFromLa", "hRecoPrFromLa;Decay Radius;#it{p}_{T}", {kTH2F, {axes.axisDecayRadius, axes.axisMomentum}});

          // basic mass histograms to see if we're in business
          insertHist(histPath + "hMassLambda", "hMassLambda;Mass;Counts", {kTH1F, {{axes.axisLambdaMass}}});
          insertHist(histPath + "hMassXi", "hMassXi;Mass;Counts", {kTH1F, {{axes.axisXiMass}}});
          insertHist(histPath + "h2dMassXi", "h2dMassXi;Mass;Counts", {kTH2F, {{axes.axisXiMass}, {axes.axisMomentum}}});

          // OTF strangeness tracking QA
          insertHist(histPath + "hFoundVsFindable", "hFoundVsFindable;Found;Findable", {kTH2F, {{10, -0.5f, 9.5f}, {10, -0.5f, 9.5f}}});

          insertHist(histPath + "h2dDCAxyCascade", "h2dDCAxyCascade;#it{p}_{T};DCA_{xy}", {kTH2F, {axes.axisMomentum, axes.axisDCA}});
          insertHist(histPath + "h2dDCAxyCascadeBachelor", "h2dDCAxyCascadeBachelor;#it{p}_{T};DCA_{xy}", {kTH2F, {axes.axisMomentum, axes.axisDCA}});
          insertHist(histPath + "h2dDCAxyCascadeNegative", "h2dDCAxyCascadeNegative;#it{p}_{T};DCA_{xy}", {kTH2F, {axes.axisMomentum, axes.axisDCA}});
          insertHist(histPath + "h2dDCAxyCascadePositive", "h2dDCAxyCascadePositive;#it{p}_{T};DCA_{xy}", {kTH2F, {axes.axisMomentum, axes.axisDCA}});

          insertHist(histPath + "h2dDCAzCascade", "h2dDCAzCascade;#it{p}_{T};DCA_{z}", {kTH2F, {axes.axisMomentum, axes.axisDCA}});
          insertHist(histPath + "h2dDCAzCascadeBachelor", "h2dDCAzCascadeBachelor;#it{p}_{T};DCA_{z}", {kTH2F, {axes.axisMomentum, axes.axisDCA}});
          insertHist(histPath + "h2dDCAzCascadeNegative", "h2dDCAzCascadeNegative;#it{p}_{T};DCA_{z}", {kTH2F, {axes.axisMomentum, axes.axisDCA}});
          insertHist(histPath + "h2dDCAzCascadePositive", "h2dDCAzCascadePositive;#it{p}_{T};DCA_{z}", {kTH2F, {axes.axisMomentum, axes.axisDCA}});

          insertHist(histPath + "h2dDeltaPtVsPt", "h2dDeltaPtVsPt;Gen p_{T};#Delta p_{T}", {kTH2F, {axes.axisMomentum, axes.axisDeltaPt}});
          insertHist(histPath + "h2dDeltaEtaVsPt", "h2dDeltaEtaVsPt;Gen p_{T};#Delta #eta", {kTH2F, {axes.axisMomentum, axes.axisDeltaEta}});

          insertHist(histPath + "nSiliconHitsCascadeProngs", "nSiliconHitsCascadeProngs", {kTH1F, {{40, -0.5f, 39.5f}}});
          insertHist(histPath + "nTPCHitsCascadeProngs", "nTPCHitsCascadeProngs", {kTH1F, {{10, -0.5f, 9.5f}}});
          insertHist(histPath + "hFastTrackerHits", "hFastTrackerHits", {kTH2F, {axes.axisZ, axes.axisRadius}});
          insertHist(histPath + "hFastTrackerQA", "hFastTrackerQA", {kTH1F, {{8, -0.5f, 7.5f}}});
          getHist(TH1, histPath + "hFastTrackerQA")->GetXaxis()->SetBinLabel(1, "Negative eigenvalue");
          getHist(TH1, histPath + "hFastTrackerQA")->GetXaxis()->SetBinLabel(2, "Failed sanity check");
          getHist(TH1, histPath + "hFastTrackerQA")->GetXaxis()->SetBinLabel(3, "intercept original radius");
          getHist(TH1, histPath + "hFastTrackerQA")->GetXaxis()->SetBinLabel(4, "propagate to original radius");
          getHist(TH1, histPath + "hFastTrackerQA")->GetXaxis()->SetBinLabel(5, "problematic layer");
          getHist(TH1, histPath + "hFastTrackerQA")->GetXaxis()->SetBinLabel(6, "multiple scattering");
          getHist(TH1, histPath + "hFastTrackerQA")->GetXaxis()->SetBinLabel(7, "energy loss");
          getHist(TH1, histPath + "hFastTrackerQA")->GetXaxis()->SetBinLabel(8, "efficiency");
        }
      }

      if (doExtraQA) {
        insertHist(histPath + "h2dPtRes", "h2dPtRes;Gen p_{T};#Delta p_{T} / Reco p_{T}", {kTH2D, {{axes.axisMomentum, axes.axisPtRes}}});
        insertHist(histPath + "h2dPtResAbs", "h2dPtResAbs;Gen p_{T};#Delta p_{T}", {kTH2D, {{axes.axisMomentum, axes.axisPtRes}}});
        insertHist(histPath + "h2dDCAxy", "h2dDCAxy;p_{T};DCA_{xy}", {kTH2D, {{axes.axisMomentum, axes.axisDCA}}});
        insertHist(histPath + "h2dDCAz", "h2dDCAz;p_{T};DCA_{z}", {kTH2D, {{axes.axisMomentum, axes.axisDCA}}});
      }

    } // end config loop

    // Basic QA
    auto hNaN = histos.add<TH2>("hNaNBookkeeping", "hNaNBookkeeping", kTH2F, {{10, -0.5f, 9.5f}, {10, -0.5f, 9.5f}});

    hNaN->GetXaxis()->SetBinLabel(1, "Primary");
    hNaN->GetXaxis()->SetBinLabel(2, "Bachelor");
    hNaN->GetXaxis()->SetBinLabel(3, "Pi from La");
    hNaN->GetXaxis()->SetBinLabel(4, "Pr from La");

    hNaN->GetYaxis()->SetBinLabel(1, "Smear NaN");
    hNaN->GetYaxis()->SetBinLabel(2, "Smear OK");

    auto hCovMatOK = histos.add<TH1>("hCovMatOK", "hCovMatOK", kTH1D, {{2, -0.5f, 1.5f}});
    hCovMatOK->GetXaxis()->SetBinLabel(1, "Not OK");
    hCovMatOK->GetXaxis()->SetBinLabel(2, "OK");

    auto hFitterStatusCode = histos.add<TH1>("hFitterStatusCode", "hFitterStatusCode", kTH1D, {{15, -0.5, 14.5}});
    hFitterStatusCode->GetXaxis()->SetBinLabel(1, "None"); // no status set (should not be possible!)

    /* Good Conditions */
    hFitterStatusCode->GetXaxis()->SetBinLabel(2, "Converged"); // fit converged
    hFitterStatusCode->GetXaxis()->SetBinLabel(3, "MaxIter");   // max iterations reached before fit convergence

    /* Error Conditions */
    hFitterStatusCode->GetXaxis()->SetBinLabel(4, "NoCrossing");       // no reasaonable crossing was found
    hFitterStatusCode->GetXaxis()->SetBinLabel(5, "RejRadius");        // radius of crossing was not acceptable
    hFitterStatusCode->GetXaxis()->SetBinLabel(6, "RejTrackX");        // one candidate track x was below the mimimum required radius
    hFitterStatusCode->GetXaxis()->SetBinLabel(7, "RejTrackRoughZ");   // rejected by rough cut on tracks Z difference
    hFitterStatusCode->GetXaxis()->SetBinLabel(8, "RejChi2Max");       // rejected by maximum chi2 cut
    hFitterStatusCode->GetXaxis()->SetBinLabel(9, "FailProp");         // propagation of at least prong to PCA failed
    hFitterStatusCode->GetXaxis()->SetBinLabel(10, "FailInvCov");      // inversion of cov.-matrix failed
    hFitterStatusCode->GetXaxis()->SetBinLabel(11, "FailInvWeight");   // inversion of Ti weight matrix failed
    hFitterStatusCode->GetXaxis()->SetBinLabel(12, "FailInv2ndDeriv"); // inversion of 2nd derivatives failed
    hFitterStatusCode->GetXaxis()->SetBinLabel(13, "FailCorrTracks");  // correction of tracks to updated x failed
    hFitterStatusCode->GetXaxis()->SetBinLabel(14, "FailCloserAlt");   // alternative PCA is closer
    hFitterStatusCode->GetXaxis()->SetBinLabel(15, "NStatusesDefined");

    if (doExtraQA) {
      histos.add("h2dVerticesVsContributors", "h2dVerticesVsContributors;Multiplicity;N vertices", kTH2F, {axes.axisMultiplicity, axes.axisNVertices});
      histos.add("h1dVerticesNotReco", "h1dVerticesNotReco;Multiplicity;Vertices Not Reco", kTH1F, {axes.axisMultiplicity});
      histos.add("hRecoVsSimMultiplicity", "hRecoVsSimMultiplicity;Reco Mult.;Sim Mult.", kTH2F, {axes.axisMultiplicity, axes.axisMultiplicity});

      histos.add("hSimTrackX", "hSimTrackX", kTH1F, {axes.axisX});
      histos.add("hRecoTrackX", "hRecoTrackX", kTH1F, {axes.axisX});
      histos.add("hTrackXatDCA", "hTrackXatDCA", kTH1F, {axes.axisX});
    }

    if (v0DecaySettings.doV0QA) {
      for (int icfg = 0; icfg < nGeometries; icfg++) {
        std::string v0histPath = "V0Building_Configuration_" + std::to_string(icfg) + "/";
        insertHist(v0histPath + "hV0Building", "hV0Building", kTH1F, {{10, -0.5f, 9.5f}});
        insertHist(v0histPath + "hFastTrackerHits", "hV0Building", kTH2F, {{axes.axisZ, axes.axisRadius}});
        auto h = histos.add<TH1>(v0histPath + "hFastTrackerQA", "hFastTrackerQA", kTH1D, {{8, -0.5f, 7.5f}});
        h->GetXaxis()->SetBinLabel(1, "Negative eigenvalue");
        h->GetXaxis()->SetBinLabel(2, "Failed sanity check");
        h->GetXaxis()->SetBinLabel(3, "intercept original radius");
        h->GetXaxis()->SetBinLabel(4, "propagate to original radius");
        h->GetXaxis()->SetBinLabel(5, "problematic layer");
        h->GetXaxis()->SetBinLabel(6, "multiple scattering");
        h->GetXaxis()->SetBinLabel(7, "energy loss");
        h->GetXaxis()->SetBinLabel(8, "efficiency");
        h->GetXaxis()->SetBinLabel(8, "no layers hit");
        histPointers.insert({v0histPath + "hFastTrackerQA", h});
        // K0s
        insertHist(v0histPath + "K0/hGen", "hGen", kTH2F, {axes.axisDecayRadius, axes.axisMomentum});
        insertHist(v0histPath + "K0/hReco", "hReco", kTH2F, {axes.axisDecayRadius, axes.axisMomentum});
        insertHist(v0histPath + "K0/hGenNegDaughterFromV0", "hGenNegDaughterFromV0", kTH2F, {axes.axisDecayRadius, axes.axisMomentum});
        insertHist(v0histPath + "K0/hGenPosDaughterFromV0", "hGenPosDaughterFromV0", kTH2F, {axes.axisDecayRadius, axes.axisMomentum});
        insertHist(v0histPath + "K0/hRecoNegDaughterFromV0", "hRecoNegDaughterFromV0", kTH2F, {axes.axisDecayRadius, axes.axisMomentum});
        insertHist(v0histPath + "K0/hRecoPosDaughterFromV0", "hRecoPosDaughterFromV0", kTH2F, {axes.axisDecayRadius, axes.axisMomentum});
        insertHist(v0histPath + "K0/hMass", "hMass", kTH2F, {axes.axisK0Mass, axes.axisMomentum});
        // Lambda
        insertHist(v0histPath + "Lambda/hGen", "hGen", kTH2F, {axes.axisDecayRadius, axes.axisMomentum});
        insertHist(v0histPath + "Lambda/hReco", "hReco", kTH2F, {axes.axisDecayRadius, axes.axisMomentum});
        insertHist(v0histPath + "Lambda/hGenNegDaughterFromV0", "hGenNegDaughterFromV0", kTH2F, {axes.axisDecayRadius, axes.axisMomentum});
        insertHist(v0histPath + "Lambda/hGenPosDaughterFromV0", "hGenPosDaughterFromV0", kTH2F, {axes.axisDecayRadius, axes.axisMomentum});
        insertHist(v0histPath + "Lambda/hRecoNegDaughterFromV0", "hRecoNegDaughterFromV0", kTH2F, {axes.axisDecayRadius, axes.axisMomentum});
        insertHist(v0histPath + "Lambda/hRecoPosDaughterFromV0", "hRecoPosDaughterFromV0", kTH2F, {axes.axisDecayRadius, axes.axisMomentum});
        insertHist(v0histPath + "Lambda/hMass", "hMass", kTH2F, {axes.axisLambdaMass, axes.axisMomentum});
        // AntiLambda
        insertHist(v0histPath + "AntiLambda/hGen", "hGen", kTH2F, {axes.axisDecayRadius, axes.axisMomentum});
        insertHist(v0histPath + "AntiLambda/hReco", "hReco", kTH2F, {axes.axisDecayRadius, axes.axisMomentum});
        insertHist(v0histPath + "AntiLambda/hGenNegDaughterFromV0", "hGenNegDaughterFromV0", kTH2F, {axes.axisDecayRadius, axes.axisMomentum});
        insertHist(v0histPath + "AntiLambda/hGenPosDaughterFromV0", "hGenPosDaughterFromV0", kTH2F, {axes.axisDecayRadius, axes.axisMomentum});
        insertHist(v0histPath + "AntiLambda/hRecoNegDaughterFromV0", "hRecoNegDaughterFromV0", kTH2F, {axes.axisDecayRadius, axes.axisMomentum});
        insertHist(v0histPath + "AntiLambda/hRecoPosDaughterFromV0", "hRecoPosDaughterFromV0", kTH2F, {axes.axisDecayRadius, axes.axisMomentum});
        insertHist(v0histPath + "AntiLambda/hMass", "hMass", kTH2F, {axes.axisLambdaMass, axes.axisMomentum});
      }
    }

    LOG(info) << "Initializing magnetic field to value: " << mMagneticField << " kG";
    o2::parameters::GRPMagField grpmag;
    grpmag.setFieldUniformity(true);
    grpmag.setL3Current(30000.f * (mMagneticField / 5.0f));
    auto field = grpmag.getNominalL3Field();
    o2::base::Propagator::initFieldFromGRP(&grpmag);

    auto fieldInstance = static_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField());
    if (!fieldInstance) {
      LOGF(fatal, "Failed to set up magnetic field! Stopping now!");
    }

    // Cross-check
    LOGF(info, "Check field at (0, 0, 0): %.1f kG, nominal: %.1f", static_cast<float>(fieldInstance->GetBz(0, 0, 0)), static_cast<float>(field));
    LOGF(info, "Initializing empty material cylinder LUT - could be better in the future");
    o2::base::MatLayerCylSet* lut = new o2::base::MatLayerCylSet();
    lut->addLayer(200, 200.1, 2, 1.0f, 100.0f);
    LOGF(info, "MatLayerCylSet::optimizePhiSlices()");
    lut->optimizePhiSlices();
    LOGF(info, "Setting lut now...");
    o2::base::Propagator::Instance()->setMatLUT(lut);

    irSampler.setInteractionRate(100);
    irSampler.setFirstIR(o2::InteractionRecord(0, 0));
    irSampler.init();

    vertexer.setValidateWithIR(kFALSE);
    vertexer.setTrackSources(o2::dataformats::GlobalTrackID::ITS);
    vertexer.setBunchFilling(irSampler.getBunchFilling());
    vertexer.setBz(mMagneticField);
    o2::conf::ConfigurableParam::updateFromString("pvertexer.doBCValidation=false;" + primaryVertexOption.value);
    vertexer.init();

    o2::vertexing::PVertexerParams::Instance().printKeyValues();

    // initialize O2 2-prong fitter
    fitter.setPropagateToPCA(cfgFitter.propagateToPCA);
    fitter.setMaxR(cfgFitter.maxR);
    fitter.setMinParamChange(cfgFitter.minParamChange);
    fitter.setMinRelChi2Change(cfgFitter.minRelChi2Change);
    fitter.setMaxDZIni(cfgFitter.maxDZIni);
    fitter.setMaxDXYIni(cfgFitter.maxDXYIni);
    fitter.setMaxChi2(cfgFitter.maxVtxChi2);
    fitter.setUseAbsDCA(cfgFitter.useAbsDCA);
    fitter.setWeightedFinalPCA(cfgFitter.useWeightedFinalPCA);
    fitter.setBz(mMagneticField);
    fitter.setMatCorrType(o2::base::Propagator::MatCorrType::USEMatCorrNONE); // such a light detector here

    // Set seed for TGenPhaseSpace
    rand.SetSeed(seed);
    gRandom->SetSeed(seed);
  }

  /// Function to get the internal PID for a given pdgCode
  /// \param pdgCode pdg code for a common particle (particle handled by the tracker)
  int pdgCodeToPID(int pdgCode) const
  {
    if (std::abs(pdgCode) == PDG_t::kElectron) {
      return o2::track::PID::Electron;
    } else if (std::abs(pdgCode) == PDG_t::kMuonMinus) {
      return o2::track::PID::Muon;
    } else if (std::abs(pdgCode) == PDG_t::kPiPlus) {
      return o2::track::PID::Pion;
    } else if (std::abs(pdgCode) == PDG_t::kKPlus) {
      return o2::track::PID::Kaon;
    } else if (std::abs(pdgCode) == PDG_t::kProton) {
      return o2::track::PID::Proton;
    } else if (std::abs(pdgCode) == PDG_t::kLambda0) {
      return o2::track::PID::Lambda;
    } else if (std::abs(pdgCode) == PDG_t::kXiMinus) {
      return o2::track::PID::XiMinus;
    } else if (std::abs(pdgCode) == PDG_t::kOmegaMinus) {
      return o2::track::PID::OmegaMinus;
    } else {
      return o2::track::PID::Pion; // Default trackParCov assumption
    }
  }

  /// Function to decay the xi
  /// \param particle the particle to decay
  /// \param track track of particle to decay
  /// \param decayDaughters the address of resulting daughters
  /// \param xiDecayVertex the address of the xi decay vertex
  /// \param laDecayVertex the address of the la decay vertex
  template <typename McParticleType>
  void decayCascade(McParticleType particle, o2::track::TrackParCov track, std::vector<TLorentzVector>& decayDaughters, std::vector<double>& xiDecayVertex, std::vector<double>& laDecayVertex)
  {
    const double uXi = rand.Uniform(0, 1);
    const double ctauXi = 4.91; // cm
    const double betaGammaXi = particle.p() / o2::constants::physics::MassXiMinus;
    const double rxyzXi = (-betaGammaXi * ctauXi * std::log(1 - uXi));

    float sna, csa;
    o2::math_utils::CircleXYf_t circleXi;
    track.getCircleParams(mMagneticField, circleXi, sna, csa);
    const double rxyXi = rxyzXi / std::sqrt(1. + track.getTgl() * track.getTgl());
    const double theta = rxyXi / circleXi.rC;
    const double newX = ((particle.vx() - circleXi.xC) * std::cos(theta) - (particle.vy() - circleXi.yC) * std::sin(theta)) + circleXi.xC;
    const double newY = ((particle.vy() - circleXi.yC) * std::cos(theta) + (particle.vx() - circleXi.xC) * std::sin(theta)) + circleXi.yC;
    const double newPx = particle.px() * std::cos(theta) - particle.py() * std::sin(theta);
    const double newPy = particle.py() * std::cos(theta) + particle.px() * std::sin(theta);
    const double newE = std::sqrt(o2::constants::physics::MassXiMinus * o2::constants::physics::MassXiMinus + newPx * newPx + newPy * newPy + particle.pz() * particle.pz());

    xiDecayVertex.push_back(newX);
    xiDecayVertex.push_back(newY);
    xiDecayVertex.push_back(particle.vz() + rxyzXi * (particle.pz() / particle.p()));

    std::vector<double> xiDaughters = {o2::constants::physics::MassLambda, o2::constants::physics::MassPionCharged};
    TLorentzVector xi(newPx, newPy, particle.pz(), newE);
    TGenPhaseSpace xiDecay;
    xiDecay.SetDecay(xi, 2, xiDaughters.data());
    xiDecay.Generate();
    decayDaughters.push_back(*xiDecay.GetDecay(1));
    TLorentzVector la = *xiDecay.GetDecay(0);

    const double uLa = rand.Uniform(0, 1);
    const double ctauLa = 7.845; // cm
    const double betaGammaLa = la.P() / o2::constants::physics::MassLambda;
    const double rxyzLa = (-betaGammaLa * ctauLa * std::log(1 - uLa));
    laDecayVertex.push_back(xiDecayVertex[0] + rxyzLa * (xiDecay.GetDecay(0)->Px() / xiDecay.GetDecay(0)->P()));
    laDecayVertex.push_back(xiDecayVertex[1] + rxyzLa * (xiDecay.GetDecay(0)->Py() / xiDecay.GetDecay(0)->P()));
    laDecayVertex.push_back(xiDecayVertex[2] + rxyzLa * (xiDecay.GetDecay(0)->Pz() / xiDecay.GetDecay(0)->P()));

    std::vector<double> laDaughters = {o2::constants::physics::MassPionCharged, o2::constants::physics::MassProton};
    TGenPhaseSpace laDecay;
    laDecay.SetDecay(la, 2, laDaughters.data());
    laDecay.Generate();
    decayDaughters.push_back(*laDecay.GetDecay(0));
    decayDaughters.push_back(*laDecay.GetDecay(1));
  }

  /// Function to decay the V0
  /// \param particle the particle to decay
  /// \param decayDaughters the address of resulting daughters
  /// \param v0DecayVertex the address of the la decay vertex
  template <typename McParticleType>
  void decayV0Particle(McParticleType particle, std::vector<TLorentzVector>& decayDaughters, std::vector<double>& v0DecayVertex, int pdgCode)
  {
    double u = rand.Uniform(0, 1);
    double v0Mass = -1.;
    double negDauMass = -1.;
    double posDauMass = -1.;
    double ctau = -1.;

    switch (pdgCode) {
      case kK0Short:
      case -kK0Short:
        v0Mass = o2::constants::physics::MassK0Short;
        negDauMass = o2::constants::physics::MassPionCharged;
        posDauMass = o2::constants::physics::MassPionCharged;
        ctau = 2.68;
        break;
      case kLambda0:
        v0Mass = o2::constants::physics::MassLambda;
        negDauMass = o2::constants::physics::MassPionCharged;
        posDauMass = o2::constants::physics::MassProton;
        ctau = 7.845;
        break;
      case kLambda0Bar:
        v0Mass = o2::constants::physics::MassLambda;
        negDauMass = o2::constants::physics::MassProton;
        posDauMass = o2::constants::physics::MassPionCharged;
        ctau = 7.845;
        break;
      default:
        LOG(fatal) << "Trying to decay unsupported V0 with PDG " << pdgCode;
    }

    const double v0BetaGamma = particle.p() / v0Mass;
    const double v0rxyz = (-v0BetaGamma * ctau * std::log(1 - u));
    TLorentzVector v0(particle.px(), particle.py(), particle.pz(), particle.e());

    v0DecayVertex.push_back(particle.vx() + v0rxyz * (particle.px() / particle.p()));
    v0DecayVertex.push_back(particle.vy() + v0rxyz * (particle.py() / particle.p()));
    v0DecayVertex.push_back(particle.vz() + v0rxyz * (particle.pz() / particle.p()));
    std::vector<double> v0Daughters = {negDauMass, posDauMass};

    TGenPhaseSpace v0Decay;
    v0Decay.SetDecay(v0, 2, v0Daughters.data());
    v0Decay.Generate();
    decayDaughters.push_back(*v0Decay.GetDecay(0));
    decayDaughters.push_back(*v0Decay.GetDecay(1));
  }

  /// Function to compute dN/deta for a given set of MC particles
  /// \param dNdEta the address of the variable to fill with the computed dN/deta value
  /// \param mcParticles the set of MC particles to compute dN/deta from
  /// \param histPath the path to the histogram where the computed dN/deta value will be stored for QA purposes
  template <typename McParticleType>
  void computeDNDEta(float& dNdEta, McParticleType const& mcParticles, const std::string histPath)
  {
    for (const auto& mcParticle : mcParticles) {
      if (std::abs(mcParticle.eta()) > multEtaRange) {
        continue;
      }

      if (!mcParticle.isPhysicalPrimary()) {
        continue;
      }

      const auto pdg = std::abs(mcParticle.pdgCode());
      const bool longLivedToBeHandled = std::find(longLivedHandledPDGs.begin(), longLivedHandledPDGs.end(), pdg) != longLivedHandledPDGs.end();
      const bool nucleiToBeHandled = std::find(nucleiPDGs.begin(), nucleiPDGs.end(), pdg) != nucleiPDGs.end();
      const bool pdgsToBeHandled = longLivedToBeHandled || (enableNucleiSmearing && nucleiToBeHandled);
      if (!pdgsToBeHandled) {
        continue;
      }

      const auto& pdgInfo = pdgDB->GetParticle(mcParticle.pdgCode());
      if (!pdgInfo) {
        LOG(warning) << "PDG code " << mcParticle.pdgCode() << " not found in the database";
        continue;
      }
      if (pdgInfo->Charge() == 0) {
        continue;
      }
      dNdEta += 1.f;
    }
    LOG(debug) << "Computed dNch/deta before normalization: " << dNdEta;

    dNdEta /= (multEtaRange * 2.0f);
    getHist(TH1, histPath + "hLUTMultiplicity")->Fill(dNdEta);
  }

  /// Function to study the cascade decay and fill the relevant histograms and output track vector
  /// \param trackTableOffset the offset to apply to the global track indices when filling the output track vector
  /// \param mcParticle the cascade MC particle to study
  /// \param tracksCascadeProngs the address of the vector of output tracks to fill with the cascade daughters
  /// \param primaryVertex the primary vertex of the event, needed for the track propagation in the cascade decay
  /// \param icfg the index of the current configuration, needed for histogram filling
  /// \param dNdEta the dN/deta of the event, needed for the track time smearing
  /// \param eventCollisionTimeNS the collision time of the event in ns, needed for the track time smearing
  template <typename McParticleType>
  void studyCascade(const int trackTableOffset,
                    const McParticleType& mcParticle,
                    std::vector<TrackAlice3>& tracksCascadeProngs,
                    o2::vertexing::PVertex& primaryVertex,
                    int icfg,
                    int dNdEta,
                    float eventCollisionTimeNS)
  {
    o2::track::TrackParCov trackParCov;
    o2::upgrade::convertMCParticleToO2Track(mcParticle, trackParCov, pdgDB);
    const std::string histPath = "Configuration_" + std::to_string(icfg) + "/";

    std::vector<TLorentzVector> cascadeDecayProducts;
    std::vector<double> xiDecayVertex, laDecayVertex;
    static constexpr int kCascProngs = 3;
    std::array<o2::track::TrackParCov, kCascProngs> xiDaughterTrackParCovsPerfect;
    std::array<o2::track::TrackParCov, kCascProngs> xiDaughterTrackParCovsTracked;
    std::array<bool, kCascProngs> isReco;
    std::array<int, kCascProngs> nHitsCascadeProngs;        // total
    std::array<int, kCascProngs> nSiliconHitsCascadeProngs; // silicon type
    std::array<int, kCascProngs> nTPCHitsCascadeProngs;     // TPC type

    o2::track::TrackParCov xiTrackParCov;
    o2::upgrade::convertMCParticleToO2Track(mcParticle, xiTrackParCov, pdgDB);
    decayCascade(mcParticle, xiTrackParCov, cascadeDecayProducts, xiDecayVertex, laDecayVertex);
    if (cascadeDecayProducts.size() != 3) {
      LOG(fatal) << "Xi decay did not produce 3 daughters as expected!";
    }
    double xiDecayRadius2D = std::hypot(xiDecayVertex[0], xiDecayVertex[1]);
    double laDecayRadius2D = std::hypot(laDecayVertex[0], laDecayVertex[1]);

    if (cascadeDecaySettings.doXiQA) {
      getHist(TH2, histPath + "hGenXi")->Fill(xiDecayRadius2D, mcParticle.pt());
      getHist(TH2, histPath + "hGenPiFromXi")->Fill(xiDecayRadius2D, cascadeDecayProducts[0].Pt());
      getHist(TH2, histPath + "hGenPiFromLa")->Fill(laDecayRadius2D, cascadeDecayProducts[1].Pt());
      getHist(TH2, histPath + "hGenPrFromLa")->Fill(laDecayRadius2D, cascadeDecayProducts[2].Pt());
    }

    if (cascadeDecaySettings.doXiQA) {
      getHist(TH1, histPath + "hXiBuilding")->Fill(0.0f);
    }

    o2::upgrade::convertTLorentzVectorToO2Track(PDG_t::kPiMinus, cascadeDecayProducts[0], xiDecayVertex, xiDaughterTrackParCovsPerfect[0], pdgDB);
    xiDaughterTrackParCovsPerfect[0].setPID(pdgCodeToPID(PDG_t::kPiMinus));
    o2::upgrade::convertTLorentzVectorToO2Track(PDG_t::kPiMinus, cascadeDecayProducts[1], laDecayVertex, xiDaughterTrackParCovsPerfect[1], pdgDB);
    xiDaughterTrackParCovsPerfect[1].setPID(pdgCodeToPID(PDG_t::kPiMinus));
    o2::upgrade::convertTLorentzVectorToO2Track(PDG_t::kProton, cascadeDecayProducts[2], laDecayVertex, xiDaughterTrackParCovsPerfect[2], pdgDB);
    xiDaughterTrackParCovsPerfect[2].setPID(pdgCodeToPID(PDG_t::kProton));
    o2::track::TrackParCov perfectCascadeTrack;
    o2::upgrade::convertMCParticleToO2Track(mcParticle, perfectCascadeTrack, pdgDB);
    perfectCascadeTrack.setPID(pdgCodeToPID(PDG_t::kXiMinus)); // FIXME: not OK for omegas

    thisCascade.bachelorId = trackTableOffset + 1;
    thisCascade.negativeId = trackTableOffset + 2; // Lambda daughters
    thisCascade.positiveId = trackTableOffset + 3; // Lambda daughters
    thisCascade.cascadeTrackId = trackTableOffset + 4;

    float trackTime = (eventCollisionTimeNS + gRandom->Gaus(0., timeResolutionNs)) * nsToMus;
    tracksCascadeProngs.push_back(TrackAlice3{xiDaughterTrackParCovsPerfect[0], mcParticle.globalIndex(), trackTime, timeResolutionUs, true, true, 1, -1, -1, TrackType::kGenCascDaug});
    trackTime = (eventCollisionTimeNS + gRandom->Gaus(0., timeResolutionNs)) * nsToMus;
    tracksCascadeProngs.push_back(TrackAlice3{xiDaughterTrackParCovsPerfect[1], mcParticle.globalIndex(), trackTime, timeResolutionUs, true, true, 2, -1, -1, TrackType::kGenCascDaug});
    trackTime = (eventCollisionTimeNS + gRandom->Gaus(0., timeResolutionNs)) * nsToMus;
    tracksCascadeProngs.push_back(TrackAlice3{xiDaughterTrackParCovsPerfect[2], mcParticle.globalIndex(), trackTime, timeResolutionUs, true, true, 3, -1, -1, TrackType::kGenCascDaug});
    trackTime = (eventCollisionTimeNS + gRandom->Gaus(0., timeResolutionNs)) * nsToMus;
    tracksCascadeProngs.push_back(TrackAlice3{perfectCascadeTrack, mcParticle.globalIndex(), trackTime, timeResolutionUs, true, true, 3, -1, -1, TrackType::kGenCascDaug});

    for (int i = 0; i < kCascProngs; i++) {
      isReco[i] = false;
      nHitsCascadeProngs[i] = 0;
      nSiliconHitsCascadeProngs[i] = 0;
      nTPCHitsCascadeProngs[i] = 0;
      if (enableSecondarySmearing) {
        nHitsCascadeProngs[i] = fastTracker[icfg]->FastTrack(xiDaughterTrackParCovsPerfect[i], xiDaughterTrackParCovsTracked[i], dNdEta);
        nSiliconHitsCascadeProngs[i] = fastTracker[icfg]->GetNSiliconPoints();
        nTPCHitsCascadeProngs[i] = fastTracker[icfg]->GetNGasPoints();

        if (nHitsCascadeProngs[i] < 0 && cascadeDecaySettings.doXiQA) { // QA
          getHist(TH1, histPath + "hFastTrackerQA")->Fill(o2::math_utils::abs(nHitsCascadeProngs[i]));
        }

        getHist(TH1, histPath + "nSiliconHitsCascadeProngs")->Fill(nSiliconHitsCascadeProngs[i]);
        getHist(TH1, histPath + "nTPCHitsCascadeProngs")->Fill(nTPCHitsCascadeProngs[i]);
        if (nSiliconHitsCascadeProngs[i] >= fastTrackerSettings.minSiliconHits ||
            (nSiliconHitsCascadeProngs[i] >= fastTrackerSettings.minSiliconHitsIfTPCUsed &&
             nTPCHitsCascadeProngs[i] >= fastTrackerSettings.minTPCClusters)) {
          isReco[i] = true;
        } else {
          continue; // extra sure
        }
        if (cascadeDecaySettings.doXiQA) {
          getHist(TH1, histPath + "hXiBuilding")->Fill(static_cast<float>(i + 1));
        }
        isReco[i] = true;

        for (uint32_t ih = 0; ih < fastTracker[icfg]->GetNHits() && cascadeDecaySettings.doXiQA; ih++) {
          getHist(TH2, histPath + "hFastTrackerHits")->Fill(fastTracker[icfg]->GetHitZ(ih), std::hypot(fastTracker[icfg]->GetHitX(ih), fastTracker[icfg]->GetHitY(ih)));
        }
      } else {
        isReco[i] = true;
        xiDaughterTrackParCovsTracked[i] = xiDaughterTrackParCovsPerfect[i];
        if (cascadeDecaySettings.doXiQA) {
          getHist(TH1, histPath + "hXiBuilding")->Fill(static_cast<float>(i + 1));
        }
      }

      if (TMath::IsNaN(xiDaughterTrackParCovsTracked[i].getZ())) {
        isReco[i] = false;
        continue;
      } else {
        getHist(TH1, histPath + "hXiBuilding")->Fill(4.0f);
        histos.fill(HIST("hNaNBookkeeping"), i + 1, 1.0f);
      }
      trackTime = (eventCollisionTimeNS + gRandom->Gaus(0., timeResolutionNs)) * nsToMus;
      // TODO: add flag for whether it's a ghost track or not, currently assuming all are reconstructed tracks if they pass the fast tracker requirements
      TrackType trackType = isReco[i] ? TrackType::kRecoCascDaug : TrackType::kGenCascDaug;
      tracksCascadeProngs[i] = TrackAlice3{xiDaughterTrackParCovsTracked[i], mcParticle.globalIndex(), trackTime, timeResolutionUs, true, true, i + 2, nSiliconHitsCascadeProngs[i], nTPCHitsCascadeProngs[i], trackType};
    }

    bool tryKinkReco = false;
    if (!isReco[1] || !isReco[2]) {
      tryKinkReco = true; // Lambda outside acceptance, set flag for kink reco to be used if mode 1
    }

    bool reconstructedCascade = false;
    if (isReco[0] && isReco[1] && isReco[2]) {
      reconstructedCascade = true;
    }

    bool fillCascadeTable{false};

    // +-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+
    // combine particles into actual Xi candidate
    // cascade building starts here
    if (cascadeDecaySettings.findXi && reconstructedCascade && cascadeDecaySettings.doKinkReco != 2) {
      if (cascadeDecaySettings.doXiQA) {
        getHist(TH1, histPath + "hXiBuilding")->Fill(3.0f);
      }

      // use DCA fitters
      int nCand = 0;
      bool dcaFitterV0Status = true;
      try {
        nCand = fitter.process(xiDaughterTrackParCovsTracked[1], xiDaughterTrackParCovsTracked[2]);
      } catch (...) {
        LOG(error) << "Exception caught in DCA fitter for V0!";
        dcaFitterV0Status = false;
      }
      if (nCand == 0) {
        // LOG(info) << "No V0 candidate found in DCA fitter!";
        dcaFitterV0Status = false;
      }

      fitter.propagateTracksToVertex();
      if (!fitter.isPropagateTracksToVertexDone()) {
        dcaFitterV0Status = false;
      }

      // V0 found successfully
      if (dcaFitterV0Status) {
        if (cascadeDecaySettings.doXiQA) {
          getHist(TH1, histPath + "hXiBuilding")->Fill(4.0f);
        }

        std::array<float, 3> pos;
        std::array<float, 3> posCascade;
        std::array<float, 3> posP;
        std::array<float, 3> negP;
        std::array<float, 3> bachP;

        o2::track::TrackParCov pTrackAtPCA = fitter.getTrack(1); // proton (positive)
        o2::track::TrackParCov nTrackAtPCA = fitter.getTrack(0); // pion (negative)
        pTrackAtPCA.getPxPyPzGlo(posP);
        nTrackAtPCA.getPxPyPzGlo(negP);

        // get decay vertex coordinates
        const auto& vtx = fitter.getPCACandidate();
        for (int i = 0; i < 3; i++) {
          pos[i] = vtx[i];
        }

        // calculate basic V0 properties here
        // DCA to PV taken care of in daughter tracks already, not necessary
        thisCascade.dcaV0dau = TMath::Sqrt(fitter.getChi2AtPCACandidate());
        thisCascade.v0radius = std::hypot(pos[0], pos[1]);
        thisCascade.mLambda = RecoDecay::m(std::array{std::array{posP[0], posP[1], posP[2]},
                                                      std::array{negP[0], negP[1], negP[2]}},
                                           std::array{o2::constants::physics::MassProton,
                                                      o2::constants::physics::MassPionCharged});

        // go for cascade: create V0 (pseudo)track from reconstructed V0
        std::array<float, 21> covV = {0.};
        constexpr int MomInd[6] = {9, 13, 14, 18, 19, 20}; // cov matrix elements for momentum component
        for (int i = 0; i < 6; i++) {
          covV[MomInd[i]] = 1e-6;
          covV[i] = 1e-6;
        }

        o2::track::TrackParCov v0Track = o2::track::TrackParCov(
          {pos[0], pos[1], pos[2]},
          {posP[0] + negP[0], posP[1] + negP[1], posP[2] + negP[2]},
          covV, 0, true);
        v0Track.setAbsCharge(0);
        v0Track.setPID(pdgCodeToPID(PDG_t::kLambda0));

        // dca fitter step
        nCand = 0;
        bool dcaFitterCascadeStatus = true;
        try {
          nCand = fitter.process(v0Track, xiDaughterTrackParCovsTracked[0]);
        } catch (...) {
          // LOG(error) << "Exception caught in DCA fitter cascade!";
          dcaFitterCascadeStatus = false;
        }
        if (nCand == 0) {
          dcaFitterCascadeStatus = false;
        }

        fitter.propagateTracksToVertex();
        if (!fitter.isPropagateTracksToVertexDone()) {
          dcaFitterCascadeStatus = false;
        }

        const u_int8_t fitterStatusCode = fitter.getFitStatus();
        histos.fill(HIST("hFitterStatusCode"), fitterStatusCode);
        // Cascade found successfully
        if (dcaFitterCascadeStatus) {
          if (cascadeDecaySettings.doXiQA) {
            getHist(TH1, histPath + "hXiBuilding")->Fill(6.0f);
          }

          o2::track::TrackParCov bachelorTrackAtPCA = fitter.getTrack(1);
          const auto& vtxCascade = fitter.getPCACandidate();
          for (int i = 0; i < 3; i++) {
            posCascade[i] = vtxCascade[i];
          }

          // basic properties of the cascade
          thisCascade.dcacascdau = TMath::Sqrt(fitter.getChi2AtPCACandidate());
          thisCascade.cascradius = std::hypot(posCascade[0], posCascade[1]);
          bachelorTrackAtPCA.getPxPyPzGlo(bachP);

          thisCascade.mXi = RecoDecay::m(std::array{std::array{bachP[0], bachP[1], bachP[2]},
                                                    std::array{posP[0] + negP[0], posP[1] + negP[1], posP[2] + negP[2]}},
                                         std::array{o2::constants::physics::MassPionCharged,
                                                    o2::constants::physics::MassLambda});

          // initialize cascade track
          o2::track::TrackParCov cascadeTrack = fitter.createParentTrackParCov();
          cascadeTrack.setAbsCharge(-1);                      // may require more adjustments
          cascadeTrack.setPID(pdgCodeToPID(PDG_t::kXiMinus)); // FIXME: not OK for omegas

          thisCascade.cascradiusMC = xiDecayRadius2D;
          thisCascade.pt = cascadeTrack.getPt();
          thisCascade.eta = cascadeTrack.getEta();
          thisCascade.findableClusters = 0;
          thisCascade.foundClusters = 0;

          if (cascadeDecaySettings.trackXi) {
            // optionally, add the points in the layers before the decay of the Xi
            // will back-track the perfect MC cascade to relevant layers, find hit, smear and add to smeared cascade
            for (int i = fastTracker[icfg]->GetLayers().size() - 1; i >= 0; --i) {
              o2::fastsim::DetLayer layer = fastTracker[icfg]->GetLayer(i);
              if (layer.isInert()) {
                continue; // Not an active tracking layer
              }

              if (thisCascade.cascradiusMC < layer.getRadius()) {
                continue; // Cascade did not reach this layer
              }

              // cascade decayed after the corresponding radius
              thisCascade.findableClusters++; // add to findable

              // find perfect intercept XYZ
              float targetX = 1e+3;
              xiTrackParCov.getXatLabR(layer.getRadius(), targetX, mMagneticField);
              if (targetX > 999) {
                continue; // failed to find intercept
              }

              if (!xiTrackParCov.propagateTo(targetX, mMagneticField)) {
                continue; // failed to propagate
              }

              // get potential cluster position
              std::array<float, 3> posClusterCandidate;
              xiTrackParCov.getXYZGlo(posClusterCandidate);
              float r{std::hypot(posClusterCandidate[0], posClusterCandidate[1])};
              float phi{std::atan2(-posClusterCandidate[1], -posClusterCandidate[0]) + o2::constants::math::PI};

              if (layer.getResolutionRPhi() > 1e-8 && layer.getResolutionZ() > 1e-8) { // catch zero (though should not really happen...)
                phi = gRandom->Gaus(phi, std::asin(layer.getResolutionRPhi() / r));
                posClusterCandidate[0] = r * std::cos(phi);
                posClusterCandidate[1] = r * std::sin(phi);
                posClusterCandidate[2] = gRandom->Gaus(posClusterCandidate[2], layer.getResolutionZ());
              }

              if (std::isnan(phi)) {
                continue; // Catch when getXatLabR misses layer[i]
              }

              // towards adding cluster: move to track alpha
              double alpha = cascadeTrack.getAlpha();
              double xyz1[3]{
                TMath::Cos(alpha) * posClusterCandidate[0] + TMath::Sin(alpha) * posClusterCandidate[1],
                -TMath::Sin(alpha) * posClusterCandidate[0] + TMath::Cos(alpha) * posClusterCandidate[1],
                posClusterCandidate[2]};

              if (!(cascadeTrack.propagateTo(xyz1[0], mMagneticField))) {
                continue;
              }

              const o2::track::TrackParametrization<float>::dim2_t hitpoint = {static_cast<float>(xyz1[1]), static_cast<float>(xyz1[2])};
              const o2::track::TrackParametrization<float>::dim3_t hitpointcov = {layer.getResolutionRPhi() * layer.getResolutionRPhi(), 0.f, layer.getResolutionZ() * layer.getResolutionZ()};
              if (layer.isInDeadPhiRegion(phi)) {
                continue; // No hit for strangeness tracking update
              }

              cascadeTrack.update(hitpoint, hitpointcov);
              thisCascade.foundClusters++; // add to findable
            }

            if (thisCascade.foundClusters < cascadeDecaySettings.minStraTrackHits) {
              return; // We didn't find enough hits for strangeness tracking
            }
          }
          tracksCascadeProngs[kCascProngs] = TrackAlice3{cascadeTrack, mcParticle.globalIndex(), trackTime, timeResolutionUs, false, false, 1, thisCascade.foundClusters, TrackType::kGenCascDaug};
          fillCascadeTable = true;
        }
      }
    } // end cascade building

    if (isReco[0] && ((cascadeDecaySettings.doKinkReco == 1 && tryKinkReco) || cascadeDecaySettings.doKinkReco == 2)) { // mode 1 or 2
      o2::track::TrackParCov trackedCascade;
      const o2::track::TrackParCov& trackedBach = xiDaughterTrackParCovsTracked[0];
      const int nCascHits = fastTracker[icfg]->FastTrack(perfectCascadeTrack, trackedCascade, dNdEta, xiDecayRadius2D);
      reconstructedCascade = (fastTrackerSettings.minSiliconHitsForKinkReco < nCascHits) ? true : false;
      if (reconstructedCascade) {
        std::array<float, 3> pCasc;
        std::array<float, 3> pBach;
        std::array<float, 3> pV0;
        trackedCascade.getPxPyPzGlo(pCasc);
        trackedBach.getPxPyPzGlo(pBach);
        for (size_t i = 0; i < pCasc.size(); ++i) {
          pV0[i] = pCasc[i] - pBach[i];
        }

        // Reset indices
        if (!isReco[1]) {
          thisCascade.negativeId = -1;
        }
        if (!isReco[2]) {
          thisCascade.positiveId = -1;
        }

        int nCand = 0;
        bool kinkFitterOK = true;
        try {
          nCand = fitter.process(trackedCascade, trackedBach);
        } catch (...) {
          kinkFitterOK = false;
        }

        if (nCand == 0) {
          kinkFitterOK = false;
        }

        fitter.propagateTracksToVertex();
        if (!fitter.isPropagateTracksToVertexDone()) {
          kinkFitterOK = false;
        }

        const u_int8_t fitterStatusCode = fitter.getFitStatus();
        histos.fill(HIST("hFitterStatusCode"), fitterStatusCode);
        if (kinkFitterOK) {
          if (cascadeDecaySettings.doXiQA) {
            getHist(TH1, histPath + "hXiBuilding")->Fill(7.0f);
          }

          o2::track::TrackParCov newCascadeTrack = fitter.getTrack(0); // (cascade)
          std::array<float, 3> kinkVtx = {-999, -999, -999};
          kinkVtx = fitter.getPCACandidatePos();
          thisCascade.dcaV0dau = -1.f; // unknown
          thisCascade.v0radius = -1.f; // unknown
          thisCascade.dcacascdau = std::sqrt(fitter.getChi2AtPCACandidate());
          thisCascade.cascradius = std::hypot(kinkVtx[0], kinkVtx[1]);
          thisCascade.cascradiusMC = xiDecayRadius2D;
          thisCascade.mLambda = o2::constants::physics::MassLambda;
          thisCascade.findableClusters = nCascHits;
          thisCascade.foundClusters = nCascHits;
          thisCascade.pt = newCascadeTrack.getPt();
          thisCascade.eta = newCascadeTrack.getEta();
          thisCascade.mXi = RecoDecay::m(std::array{std::array{pBach[0], pBach[1], pBach[2]}, std::array{pV0[0], pV0[1], pV0[2]}},
                                         std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassLambda});
          newCascadeTrack.setPID(pdgCodeToPID(PDG_t::kXiMinus)); // FIXME: not OK for omegas
          float trackTime = (eventCollisionTimeNS + gRandom->Gaus(0., timeResolutionNs)) * nsToMus;
          if (reconstructedCascade) {
            tracksCascadeProngs[kCascProngs + 1] = TrackAlice3{newCascadeTrack, mcParticle.globalIndex(), trackTime, timeResolutionUs, false, false, 1, thisCascade.foundClusters, TrackType::kRecoCascDaug};
          }
          fillCascadeTable = true;
        } // end fitter OK
      } // end cascade found
    } // end cascade kink building

    // +-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+
    if (cascadeDecaySettings.doXiQA) {
      double dcaXY{-1.}, dcaZ{-1.};
      if (reconstructedCascade) {
        getHist(TH2, histPath + "hRecoXi")->Fill(xiDecayRadius2D, mcParticle.pt());
        getHist(TH1, histPath + "hMassLambda")->Fill(thisCascade.mLambda);
        getHist(TH1, histPath + "hMassXi")->Fill(thisCascade.mXi);
        getHist(TH2, histPath + "h2dMassXi")->Fill(thisCascade.mXi, thisCascade.pt);
        getHist(TH2, histPath + "h2dDeltaPtVsPt")->Fill(thisCascade.pt, (mcParticle.pt() - thisCascade.pt) / thisCascade.pt);
        getHist(TH2, histPath + "h2dDeltaEtaVsPt")->Fill(thisCascade.pt, mcParticle.eta() - thisCascade.eta);
        getHist(TH2, histPath + "hFoundVsFindable")->Fill(thisCascade.findableClusters, thisCascade.foundClusters);

        o2::track::TrackParCov trackParametrization(xiTrackParCov);
        trackParametrization.propagateToDCA(primaryVertex, mMagneticField, &dcaInfo);
        getHist(TH2, histPath + "h2dDCAxyCascade")->Fill(trackParametrization.getPt(), dcaXY * 1e+4); // in microns, please
        getHist(TH2, histPath + "h2dDCAzCascade")->Fill(trackParametrization.getPt(), dcaZ * 1e+4);   // in microns, please
      }
      if (isReco[0]) {
        getHist(TH2, histPath + "hRecoPiFromXi")->Fill(xiDecayRadius2D, cascadeDecayProducts[0].Pt());
        o2::track::TrackParCov trackParametrizationCascProng0(xiTrackParCov);
        if (populateTracksDCA && xiTrackParCov.propagateToDCA(primaryVertex, mMagneticField, &dcaInfo)) { // FIXME: this is not the right trackParametrization, need to propagate the bachelor track
          dcaXY = dcaInfo.getY();
          dcaZ = dcaInfo.getZ();
          getHist(TH2, histPath + "h2dDCAxyCascadeBachelor")->Fill(trackParametrizationCascProng0.getPt(), dcaXY * 1e+4); // in microns, please
          getHist(TH2, histPath + "h2dDCAzCascadeBachelor")->Fill(trackParametrizationCascProng0.getPt(), dcaZ * 1e+4);   // in microns, please
        }
      }
      if (isReco[1]) {
        getHist(TH2, histPath + "hRecoPiFromLa")->Fill(laDecayRadius2D, cascadeDecayProducts[1].Pt());
        o2::track::TrackParCov trackParametrizationCascProng1(xiTrackParCov);
        if (populateTracksDCA && xiTrackParCov.propagateToDCA(primaryVertex, mMagneticField, &dcaInfo)) { // FIXME: this is not the right trackParametrization, need to propagate the negative pion track
          dcaXY = dcaInfo.getY();
          dcaZ = dcaInfo.getZ();
          getHist(TH2, histPath + "h2dDCAxyCascadeNegative")->Fill(trackParametrizationCascProng1.getPt(), dcaXY * 1e+4); // in microns, please
          getHist(TH2, histPath + "h2dDCAzCascadeNegative")->Fill(trackParametrizationCascProng1.getPt(), dcaZ * 1e+4);   // in microns, please
        }
      }
      if (isReco[2]) {
        getHist(TH2, histPath + "hRecoPrFromLa")->Fill(laDecayRadius2D, cascadeDecayProducts[2].Pt());
        o2::track::TrackParCov trackParametrizationCascProng2(xiTrackParCov);
        if (populateTracksDCA && xiTrackParCov.propagateToDCA(primaryVertex, mMagneticField, &dcaInfo)) { // FIXME: this is not the right trackParametrization, need to propagate the positive proton track
          dcaXY = dcaInfo.getY();
          dcaZ = dcaInfo.getZ();
          getHist(TH2, histPath + "h2dDCAxyCascadePositive")->Fill(trackParametrizationCascProng2.getPt(), dcaXY * 1e+4); // in microns, please
          getHist(TH2, histPath + "h2dDCAzCascadePositive")->Fill(trackParametrizationCascProng2.getPt(), dcaZ * 1e+4);   // in microns, please
        }
      }
    }

    if (!fillCascadeTable) {
      tracksCascadeProngs.clear(); // clear the tracks added for this cascade since we won't be filling the table
      return;
    }

    // populate Cascades
    tableUpgradeCascades(tableCollisions.lastIndex(),
                         thisCascade.cascadeTrackId,
                         thisCascade.positiveId,
                         thisCascade.negativeId,
                         thisCascade.bachelorId,
                         thisCascade.dcaV0dau,
                         thisCascade.dcacascdau,
                         thisCascade.v0radius,
                         thisCascade.cascradius,
                         thisCascade.cascradiusMC,
                         thisCascade.mLambda,
                         thisCascade.mXi,
                         thisCascade.findableClusters,
                         thisCascade.foundClusters);
  }

  /// Function to study V0s and fill the relevant histograms
  /// \param trackTableOffset offset to be added to the indices of the daughter tracks when filling the V0 table
  /// \param mcParticle the MC particle corresponding to the V0 to be studied
  /// \param tracksV0Daugs vector of tracks to which the V0 daughters will be added
  /// \param icfg index of the current configuration, used for histogram filling
  /// \param dNdEta multiplicity of the event, used for fast tracker smearing
  /// \param eventCollisionTimeNS collision time of the event in nanoseconds, used for track time smearing
  template <typename McParticleType>
  void studyV0(const int trackTableOffset,
               const McParticleType& mcParticle,
               std::vector<TrackAlice3>& tracksV0Daugs,
               int icfg,
               int dNdEta,
               float eventCollisionTimeNS)
  {
    o2::track::TrackParCov trackParCov;
    o2::upgrade::convertMCParticleToO2Track(mcParticle, trackParCov, pdgDB);
    const std::string histPath = "Configuration_" + std::to_string(icfg) + "/";

    std::vector<TLorentzVector> v0DecayProducts;
    std::vector<double> laDecayVertex, v0DecayVertex;
    decayV0Particle(mcParticle, v0DecayProducts, v0DecayVertex, mcParticle.pdgCode());
    if (v0DecayProducts.size() != 2) {
      LOG(fatal) << "V0 decay did not produce 2 daughters as expected!";
    }
    double v0DecayRadius2D = std::hypot(v0DecayVertex[0], v0DecayVertex[1]);

    if (v0DecaySettings.doV0QA) {
      for (size_t indexV0 = 0; indexV0 < v0PDGs.size(); indexV0++) {
        if (mcParticle.pdgCode() != v0PDGs[indexV0]) {
          continue;
        }
        for (int indexDetector = 0; indexDetector < mGeoContainer.getNumberOfConfigurations(); indexDetector++) {
          std::string path = Form("V0Building_Configuration_%i/%s/", indexDetector, NameV0s[indexV0].data());
          fillHist(TH2, path + "hGen", v0DecayRadius2D, mcParticle.pt());
          fillHist(TH2, path + "hGenNegDaughterFromV0", v0DecayRadius2D, v0DecayProducts[0].Pt());
          fillHist(TH2, path + "hGenPosDaughterFromV0", v0DecayRadius2D, v0DecayProducts[1].Pt());
        }
      }
    }

    // V0 handling
    std::vector<o2::track::TrackParCov> v0DaughterTrackParCovsPerfect(2);
    std::vector<o2::track::TrackParCov> v0DaughterTrackParCovsTracked(2);
    std::vector<bool> isV0Reco(kv0Prongs);
    std::vector<int> nV0Hits(kv0Prongs);        // total
    std::vector<int> nV0SiliconHits(kv0Prongs); // silicon type
    std::vector<int> nV0TPCHits(kv0Prongs);     // TPC type
    if (v0DecaySettings.doV0QA) {
      fillHist(TH1, Form("V0Building_Configuration_%i/hV0Building", icfg), 0.0f);
    }
    switch (mcParticle.pdgCode()) {
      case kK0Short:
        o2::upgrade::convertTLorentzVectorToO2Track(kPiMinus, v0DecayProducts[0], v0DecayVertex, v0DaughterTrackParCovsPerfect[0], pdgDB);
        o2::upgrade::convertTLorentzVectorToO2Track(kPiPlus, v0DecayProducts[1], v0DecayVertex, v0DaughterTrackParCovsPerfect[1], pdgDB);
        break;
      case kLambda0:
        o2::upgrade::convertTLorentzVectorToO2Track(kPiMinus, v0DecayProducts[0], v0DecayVertex, v0DaughterTrackParCovsPerfect[0], pdgDB);
        o2::upgrade::convertTLorentzVectorToO2Track(kProton, v0DecayProducts[1], v0DecayVertex, v0DaughterTrackParCovsPerfect[1], pdgDB);
        break;
      case kLambda0Bar:
        o2::upgrade::convertTLorentzVectorToO2Track(kProtonBar, v0DecayProducts[0], v0DecayVertex, v0DaughterTrackParCovsPerfect[0], pdgDB);
        o2::upgrade::convertTLorentzVectorToO2Track(kPiPlus, v0DecayProducts[1], v0DecayVertex, v0DaughterTrackParCovsPerfect[1], pdgDB);
        break;
      default:
        LOG(fatal) << "Unhandled V0 PDG code " << mcParticle.pdgCode();
    }

    // Store not reconstructed daughters, will update them in case reconstruction is successful
    float trackTime = (eventCollisionTimeNS + gRandom->Gaus(0., timeResolutionNs)) * nsToMus;
    tracksV0Daugs.push_back(TrackAlice3{v0DaughterTrackParCovsPerfect[0], mcParticle.globalIndex(), 0.f, timeResolutionUs, true, true, 1, TrackType::kGhostV0Daug});
    trackTime = (eventCollisionTimeNS + gRandom->Gaus(0., timeResolutionNs)) * nsToMus;
    tracksV0Daugs.push_back(TrackAlice3{v0DaughterTrackParCovsPerfect[1], mcParticle.globalIndex(), 0.f, timeResolutionUs, true, true, 1, TrackType::kGhostV0Daug});

    bool fillV0Table{false};
    for (int i = 0; i < kv0Prongs; i++) {
      isV0Reco[i] = false;
      nV0Hits[i] = 0;
      nV0SiliconHits[i] = 0;
      nV0TPCHits[i] = 0;
      if (enableSecondarySmearing) {
        nV0Hits[i] = fastTracker[icfg]->FastTrack(v0DaughterTrackParCovsPerfect[i], v0DaughterTrackParCovsTracked[i], dNdEta);
        nV0SiliconHits[i] = fastTracker[icfg]->GetNSiliconPoints();
        nV0TPCHits[i] = fastTracker[icfg]->GetNGasPoints();

        if (nV0SiliconHits[i] >= fastTrackerSettings.minSiliconHits ||
            (nV0SiliconHits[i] >= fastTrackerSettings.minSiliconHitsIfTPCUsed &&
             nV0TPCHits[i] >= fastTrackerSettings.minTPCClusters)) {
          isV0Reco[i] = true;
        } else {
          continue; // extra sure
        }

        if (v0DecaySettings.doV0QA) { // QA
          if (nV0Hits[i] < 0) {
            fillHist(TH1, Form("V0Building_Configuration_%i/hFastTrackerQA", icfg), o2::math_utils::abs(nV0Hits[i]));
          }
          for (uint32_t ih = 0; ih < fastTracker[icfg]->GetNHits(); ih++) {
            fillHist(TH2, Form("V0Building_Configuration_%i/hFastTrackerHits", icfg), fastTracker[icfg]->GetHitZ(ih), std::hypot(fastTracker[icfg]->GetHitX(ih), fastTracker[icfg]->GetHitY(ih)));
          }
        }
      } else {
        isV0Reco[i] = true;
        v0DaughterTrackParCovsTracked[i] = v0DaughterTrackParCovsPerfect[i];
      }

      trackTime = (eventCollisionTimeNS + gRandom->Gaus(0., timeResolutionNs)) * nsToMus;
      // TODO: flag to separate ghost and reco tracks
      TrackType trackType = isV0Reco[i] ? TrackType::kRecoV0Daug : TrackType::kGhostV0Daug;
      tracksV0Daugs[i] = TrackAlice3{v0DaughterTrackParCovsTracked[i], mcParticle.globalIndex(), trackTime, timeResolutionUs, true, true, i + 2, trackType};
    }
    if (v0DecaySettings.doV0QA) {
      if (isV0Reco[0] && isV0Reco[1]) {
        fillHist(TH1, Form("V0Building_Configuration_%i/hV0Building", icfg), 1.0f);
        for (size_t indexV0 = 0; indexV0 < v0PDGs.size(); indexV0++) {
          if (mcParticle.pdgCode() == v0PDGs[indexV0]) {
            fillHist(TH2, Form("V0Building_Configuration_%i/%s/hReco", icfg, NameV0s[indexV0].data()), v0DecayRadius2D, mcParticle.pt());
          }
        }
      }
      if (isV0Reco[0]) {
        for (size_t indexV0 = 0; indexV0 < v0PDGs.size(); indexV0++) {
          if (mcParticle.pdgCode() == v0PDGs[indexV0]) {
            fillHist(TH2, Form("V0Building_Configuration_%i/%s/hRecoNegDaughterFromV0", icfg, NameV0s[indexV0].data()), v0DecayRadius2D, v0DecayProducts[0].Pt());
          }
        }
      }
      if (isV0Reco[1]) {
        for (size_t indexV0 = 0; indexV0 < v0PDGs.size(); indexV0++) {
          if (mcParticle.pdgCode() == v0PDGs[indexV0]) {
            fillHist(TH2, Form("V0Building_Configuration_%i/%s/hRecoPosDaughterFromV0", icfg, NameV0s[indexV0].data()), v0DecayRadius2D, v0DecayProducts[1].Pt());
          }
        }
      }
    }

    // +-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+
    // combine particles into actual V0 candidate
    // V0 building starts here
    if (v0DecaySettings.findV0 && isV0Reco[0] && isV0Reco[1]) {
      if (v0DecaySettings.doV0QA) {
        fillHist(TH1, Form("V0Building_Configuration_%i/hV0Building", icfg), 2.0f);
      }

      // assign indices of the daughter particles
      thisV0.mcParticleId = mcParticle.globalIndex();
      thisV0.negativeId = trackTableOffset + 1;
      thisV0.positiveId = trackTableOffset + 2;

      // use DCA fitters
      int nCand = 0;
      bool dcaFitterV0Status = true;
      try {
        nCand = fitter.process(v0DaughterTrackParCovsTracked[0], v0DaughterTrackParCovsTracked[1]);
      } catch (...) {
        // LOG(error) << "Exception caught in DCA fitter process call!";
        dcaFitterV0Status = false;
      }
      if (nCand == 0) {
        dcaFitterV0Status = false;
      }

      fitter.propagateTracksToVertex();
      if (!fitter.isPropagateTracksToVertexDone()) {
        dcaFitterV0Status = false;
      }

      const u_int8_t fitterStatusCode = fitter.getFitStatus();
      histos.fill(HIST("hFitterStatusCode"), fitterStatusCode);
      // V0 found successfully
      if (dcaFitterV0Status) {
        if (v0DecaySettings.doV0QA) {
          fillHist(TH1, Form("V0Building_Configuration_%i/hV0Building", icfg), 3.0f);
        }

        std::array<float, 3> pos;
        std::array<float, 3> posP;
        std::array<float, 3> negP;

        o2::track::TrackParCov pTrackAtPCA = fitter.getTrack(1); //  (positive)
        o2::track::TrackParCov nTrackAtPCA = fitter.getTrack(0); //  (negative)
        pTrackAtPCA.getPxPyPzGlo(posP);
        nTrackAtPCA.getPxPyPzGlo(negP);

        // get decay vertex coordinates
        const auto& vtx = fitter.getPCACandidate();
        for (int i = 0; i < 3; i++) {
          pos[i] = vtx[i];
        }

        // calculate basic V0 properties here
        // DCA to PV taken care of in daughter tracks already, not necessary
        thisV0.dcaV0dau = TMath::Sqrt(fitter.getChi2AtPCACandidate());
        thisV0.v0radius = std::hypot(pos[0], pos[1]);
        thisV0.pt = std::hypot(std::cos(v0DaughterTrackParCovsTracked[0].getPhi()) * v0DaughterTrackParCovsTracked[0].getPt() + std::cos(v0DaughterTrackParCovsTracked[1].getPhi()) * v0DaughterTrackParCovsTracked[1].getPt(),
                               std::sin(v0DaughterTrackParCovsTracked[0].getPhi()) * v0DaughterTrackParCovsTracked[0].getPt() + std::sin(v0DaughterTrackParCovsTracked[1].getPhi()) * v0DaughterTrackParCovsTracked[1].getPt());

        thisV0.mK0 = -1, thisV0.mLambda = -1, thisV0.mAntiLambda = -1; // initialize to unphysical values
        if (std::abs(mcParticle.pdgCode()) == kK0Short) {
          thisV0.mK0 = RecoDecay::m(std::array{std::array{posP[0], posP[1], posP[2]},
                                               std::array{negP[0], negP[1], negP[2]}},
                                    std::array{o2::constants::physics::MassPionCharged,
                                               o2::constants::physics::MassPionCharged});
        }

        if (mcParticle.pdgCode() == kLambda0) {
          thisV0.mLambda = RecoDecay::m(std::array{std::array{posP[0], posP[1], posP[2]},
                                                   std::array{negP[0], negP[1], negP[2]}},
                                        std::array{o2::constants::physics::MassProton,
                                                   o2::constants::physics::MassPionCharged});
        }

        if (mcParticle.pdgCode() == kLambda0Bar) {
          thisV0.mAntiLambda = RecoDecay::m(std::array{std::array{posP[0], posP[1], posP[2]},
                                                       std::array{negP[0], negP[1], negP[2]}},
                                            std::array{o2::constants::physics::MassPionCharged,
                                                       o2::constants::physics::MassProton});
        }

        if (v0DecaySettings.doV0QA) {
          fillHist(TH1, Form("V0Building_Configuration_%i/hV0Building", icfg), 4.0f);
          if (std::abs(mcParticle.pdgCode()) == kK0Short) {
            fillHist(TH2, Form("V0Building_Configuration_%i/K0/hMass", icfg), thisV0.mK0, thisV0.pt);
          }
          if (mcParticle.pdgCode() == kLambda0) {
            fillHist(TH2, Form("V0Building_Configuration_%i/Lambda/hMass", icfg), thisV0.mLambda, thisV0.pt);
          }
          if (mcParticle.pdgCode() == kLambda0Bar) {
            fillHist(TH2, Form("V0Building_Configuration_%i/AntiLambda/hMass", icfg), thisV0.mAntiLambda, thisV0.pt);
          }
        }

        LOG(info) << "Will fill V0 table for MC particle with PDG " << mcParticle.pdgCode() << " and global index " << mcParticle.globalIndex();
        fillV0Table = true;
      }
    }

    if (!fillV0Table) {
      tracksV0Daugs.clear(); // clear the tracks added for this V0 since we won't be filling the table
      return;                // don't fill the table if we didn't find a V0 candidate
    }

    // populate V0s
    tableUpgradeV0s(tableCollisions.lastIndex(), // now we know the collision index -> populate table
                    thisV0.mcParticleId,
                    thisV0.positiveId,
                    thisV0.negativeId,
                    thisV0.dcaV0dau,
                    thisV0.v0radius,
                    thisV0.mLambda,
                    thisV0.mAntiLambda,
                    thisV0.mK0,
                    thisV0.pt);
  }

  /// Function to compute the primary vertex position using the provided tracks and MC collision information
  /// \param mcCollision the MC collision information, used to get the true vertex position for comparison
  /// \param prmTrks the vector of tracks to be used for vertex reconstruction
  /// \param primaryVertex the output variable where the computed primary vertex will be stored
  /// \param icfg index of the current configuration, used for histogram filling
  template <typename McCollisionType, typename TrackType>
  void computeVertex(McCollisionType& mcCollision, const std::vector<TrackType>& prmTrks, o2::vertexing::PVertex& primaryVertex, const int icfg)
  {

    if (!enablePrimaryVertexing) {
      primaryVertex.setXYZ(mcCollision.posX(), mcCollision.posY(), mcCollision.posZ());
      return;
    }
    vertexReconstructionEfficiencyCounters.first += 1;
    const std::string histPath = "Configuration_" + std::to_string(icfg) + "/";
    fillHist(TH1, histPath + "hVtxMultGen", prmTrks.size());
    std::vector<o2::MCCompLabel> lblTracks;
    std::vector<o2::vertexing::PVertex> vertices;
    std::vector<o2::vertexing::GIndex> vertexTrackIDs;
    std::vector<o2::vertexing::V2TRef> v2tRefs;
    std::vector<o2::MCEventLabel> lblVtx;
    lblVtx.emplace_back(mcCollision.globalIndex(), 1);
    std::vector<o2::dataformats::GlobalTrackID> idxVec; // store IDs

    idxVec.reserve(prmTrks.size());
    for (unsigned i = 0; i < prmTrks.size(); i++) {
      lblTracks.emplace_back(prmTrks[i].mcLabel, mcCollision.globalIndex(), 1, false);
      idxVec.emplace_back(i, o2::dataformats::GlobalTrackID::ITS); // let's say ITS
    }

    getHist(TH1, histPath + "hVtxTrials")->Fill(0); // Tried vertexing

    // Calculate vertices
    const int n_vertices = vertexer.process(prmTrks, // track array
                                            idxVec,
                                            gsl::span<o2::InteractionRecord>{bcData},
                                            vertices,
                                            vertexTrackIDs,
                                            v2tRefs,
                                            gsl::span<const o2::MCCompLabel>{lblTracks},
                                            lblVtx);

    if (n_vertices < 1) {
      if (doExtraQA) {
        histos.fill(HIST("h1dVerticesNotReco"), prmTrks.size());
      }
      return; // primary vertex not reconstructed
    }
    vertexReconstructionEfficiencyCounters.second += 1;
    getHist(TH1, histPath + "hVtxTrials")->Fill(1); // Succeeded vertexing

    // Find largest vertex
    int largestVertex = 0;
    for (Int_t iv = 1; iv < n_vertices; iv++) {
      if (vertices[iv].getNContributors() > vertices[largestVertex].getNContributors()) {
        largestVertex = iv;
      }
    }
    primaryVertex = vertices[largestVertex];
    if (doExtraQA) {
      histos.fill(HIST("h2dVerticesVsContributors"), primaryVertex.getNContributors(), n_vertices);
    }
    fillHist(TH1, histPath + "hVtxMultReco", primaryVertex.getNContributors());
    fillHist(TH1, histPath + "hDeltaMultPVRecoGen", static_cast<int>(primaryVertex.getNContributors()) - static_cast<int>(prmTrks.size()));
    fillHist(TH2, histPath + "hDeltaXPVRecoGen", primaryVertex.getX() - mcCollision.posX(), primaryVertex.getNContributors());
    fillHist(TH2, histPath + "hDeltaYPVRecoGen", primaryVertex.getY() - mcCollision.posY(), primaryVertex.getNContributors());
    fillHist(TH2, histPath + "hDeltaZPVRecoGen", primaryVertex.getZ() - mcCollision.posZ(), primaryVertex.getNContributors());
  }

  /// Function to fill track information into the relevant tables and histograms
  /// \param tracks the vector of tracks to be processed
  /// \param primaryVertex the primary vertex position, used for DCA calculation
  /// \param icfg index of the current configuration, used for histogram filling
  void fillTracksInfo(std::vector<TrackAlice3> const& tracks, o2::vertexing::PVertex const& primaryVertex, const int icfg)
  {
    const std::string histPath = "Configuration_" + std::to_string(icfg) + "/";

    for (const auto& trackParCov : tracks) {
      // Fixme: collision index could be changeable
      aod::track::TrackTypeEnum trackType = aod::track::Track;

      if (populateTracksDCA) {
        float dcaXY = 1e+10, dcaZ = 1e+10;
        o2::track::TrackParCov trackParametrization(trackParCov);
        if (trackParametrization.propagateToDCA(primaryVertex, mMagneticField, &dcaInfo)) {
          dcaXY = dcaInfo.getY();
          dcaZ = dcaInfo.getZ();
        }
        if (doExtraQA && (!extraQAwithoutDecayDaughters || (extraQAwithoutDecayDaughters && !trackParCov.isDecayDau))) {
          getHist(TH2, histPath + "h2dDCAxy")->Fill(trackParametrization.getPt(), dcaXY * 1e+4);
          getHist(TH2, histPath + "h2dDCAz")->Fill(trackParametrization.getPt(), dcaZ * 1e+4);
          histos.fill(HIST("hTrackXatDCA"), trackParametrization.getX());
        }
        tableTracksDCA(dcaXY, dcaZ);
        if (populateTracksDCACov) {
          tableTracksDCACov(dcaInfo.getSigmaY2(), dcaInfo.getSigmaZ2());
        }
      }
      tableStoredTracks(tableCollisions.lastIndex(), trackType, trackParCov.getX(), trackParCov.getAlpha(), trackParCov.getY(), trackParCov.getZ(), trackParCov.getSnp(), trackParCov.getTgl(), trackParCov.getQ2Pt());
      tableTracksExtension(trackParCov.getPt(), trackParCov.getP(), trackParCov.getEta(), trackParCov.getPhi());

      // TODO do we keep the rho as 0? Also the sigma's are duplicated information
      tableStoredTracksCov(std::sqrt(trackParCov.getSigmaY2()), std::sqrt(trackParCov.getSigmaZ2()), std::sqrt(trackParCov.getSigmaSnp2()),
                           std::sqrt(trackParCov.getSigmaTgl2()), std::sqrt(trackParCov.getSigma1Pt2()), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
      tableTracksCovExtension(trackParCov.getSigmaY2(), trackParCov.getSigmaZY(), trackParCov.getSigmaZ2(), trackParCov.getSigmaSnpY(),
                              trackParCov.getSigmaSnpZ(), trackParCov.getSigmaSnp2(), trackParCov.getSigmaTglY(), trackParCov.getSigmaTglZ(), trackParCov.getSigmaTglSnp(),
                              trackParCov.getSigmaTgl2(), trackParCov.getSigma1PtY(), trackParCov.getSigma1PtZ(), trackParCov.getSigma1PtSnp(), trackParCov.getSigma1PtTgl(),
                              trackParCov.getSigma1Pt2());
      tableMcTrackLabels(trackParCov.mcLabel, 0);
      tableTracksExtraA3(trackParCov.nSiliconHits, trackParCov.nTPCHits, trackParCov.trackType);

      // populate extra tables if required to do so
      if (populateTracksExtra) {
        tableStoredTracksExtra(0.0f, static_cast<uint32_t>(0), static_cast<uint8_t>(0), static_cast<uint8_t>(0), static_cast<uint8_t>(0),
                               static_cast<int8_t>(0), static_cast<int8_t>(0), static_cast<uint8_t>(0), static_cast<uint8_t>(0),
                               0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
                               0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
      }
      if (populateTrackSelection) {
        tableTrackSelection(static_cast<uint8_t>(0), false, false, false, false, false, false);
        tableTrackSelectionExtension(false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false);
      }
      tableTracksAlice3(true);
    }
  }

  void processWithLUTs(aod::McCollision const& mcCollision, aod::McParticles const& mcParticles, const int icfg)
  {
    const std::string histPath = "Configuration_" + std::to_string(icfg) + "/";

    std::vector<int> genCascades;
    std::vector<int> genV0s;
    bcData.clear();
    recoPrimaries.clear();
    ghostPrimaries.clear();

    // *+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*
    // Study collision and perform vertexing
    float dNdEta{0.f}; // Charged particle multiplicity to use in the efficiency evaluation
    computeDNDEta(dNdEta, mcParticles, histPath);
    auto ir = irSampler.generateCollisionTime();
    const float eventCollisionTimeNS = ir.timeInBCNS;

    uint32_t multiplicityCounter = 0;
    // Now that the multiplicity is known, we can process the particles to smear them
    for (const auto& mcParticle : mcParticles) {

      if (!mcParticle.isPhysicalPrimary()) {
        continue;
      }

      if ((std::fabs(mcParticle.eta()) > maxEta) || (mcParticle.pt() < minPt)) {
        continue;
      }

      const bool isCascadeToDecay = (mcParticle.pdgCode() == kXiMinus) && cascadeDecaySettings.decayXi;
      const bool isV0ToDecay = std::find(v0PDGs.begin(), v0PDGs.end(), mcParticle.pdgCode()) != v0PDGs.end() && v0DecaySettings.decayV0;

      const bool longLivedToBeHandled = std::find(longLivedHandledPDGs.begin(), longLivedHandledPDGs.end(), std::abs(mcParticle.pdgCode())) != longLivedHandledPDGs.end();
      const bool nucleiToBeHandled = std::find(nucleiPDGs.begin(), nucleiPDGs.end(), std::abs(mcParticle.pdgCode())) != nucleiPDGs.end();
      const bool pdgsToBeHandled = longLivedToBeHandled ||
                                   (enableNucleiSmearing && nucleiToBeHandled) ||
                                   (isCascadeToDecay) || (isV0ToDecay);
      if (!pdgsToBeHandled) {
        continue;
      }

      if (isCascadeToDecay) {
        genCascades.push_back(mcParticle.globalIndex());
      }
      if (isV0ToDecay) {
        genV0s.push_back(mcParticle.globalIndex());
      }

      multiplicityCounter++;
      o2::track::TrackParCov trackParCov;
      const bool isDecayDaughter = (mcParticle.getProcess() == TMCProcess::kPDecay);
      if (doExtraQA) {
        histos.fill(HIST("hSimTrackX"), trackParCov.getX());
      }

      bool reconstructed = true;
      int nTrkHits = 0;
      if (enablePrimarySmearing) {
        if (fastPrimaryTrackerSettings.fastTrackPrimaries) {
          o2::track::TrackParCov perfectTrackParCov;
          o2::upgrade::convertMCParticleToO2Track(mcParticle, perfectTrackParCov, pdgDB);
          perfectTrackParCov.setPID(pdgCodeToPID(mcParticle.pdgCode()));
          nTrkHits = fastTracker[icfg]->FastTrack(perfectTrackParCov, trackParCov, dNdEta);
          if (nTrkHits < fastPrimaryTrackerSettings.minSiliconHits) {
            reconstructed = false;
          }
        } else {
          o2::upgrade::convertMCParticleToO2Track(mcParticle, trackParCov, pdgDB);
          reconstructed = mSmearer[icfg]->smearTrack(trackParCov, mcParticle.pdgCode(), dNdEta);
          nTrkHits = fastTrackerSettings.minSiliconHits;
        }
        getHist(TH1, histPath + "hPtGenerated")->Fill(mcParticle.pt());
        getHist(TH1, histPath + "hPhiGenerated")->Fill(mcParticle.phi());
        if (std::abs(mcParticle.pdgCode()) == kElectron)
          getHist(TH1, histPath + "hPtGeneratedEl")->Fill(mcParticle.pt());
        if (std::abs(mcParticle.pdgCode()) == kPiPlus)
          getHist(TH1, histPath + "hPtGeneratedPi")->Fill(mcParticle.pt());
        if (std::abs(mcParticle.pdgCode()) == kKPlus)
          getHist(TH1, histPath + "hPtGeneratedKa")->Fill(mcParticle.pt());
        if (std::abs(mcParticle.pdgCode()) == kProton)
          getHist(TH1, histPath + "hPtGeneratedPr")->Fill(mcParticle.pt());

        if (!reconstructed && !processUnreconstructedTracks) {
          continue;
        }

        getHist(TH1, histPath + "hPtReconstructed")->Fill(trackParCov.getPt());
        if (std::abs(mcParticle.pdgCode()) == kElectron)
          getHist(TH1, histPath + "hPtReconstructedEl")->Fill(trackParCov.getPt());
        if (std::abs(mcParticle.pdgCode()) == kPiPlus)
          getHist(TH1, histPath + "hPtReconstructedPi")->Fill(trackParCov.getPt());
        if (std::abs(mcParticle.pdgCode()) == kKPlus)
          getHist(TH1, histPath + "hPtReconstructedKa")->Fill(trackParCov.getPt());
        if (std::abs(mcParticle.pdgCode()) == kProton)
          getHist(TH1, histPath + "hPtReconstructedPr")->Fill(trackParCov.getPt());
      }
      if (doExtraQA) {
        getHist(TH2, histPath + "h2dPtRes")->Fill(trackParCov.getPt(), (trackParCov.getPt() - mcParticle.pt()) / trackParCov.getPt());
        getHist(TH2, histPath + "h2dPtResAbs")->Fill(trackParCov.getPt(), trackParCov.getPt() - mcParticle.pt());
        histos.fill(HIST("hRecoTrackX"), trackParCov.getX());
      }

      if (TMath::IsNaN(trackParCov.getZ())) {
        histos.fill(HIST("hNaNBookkeeping"), 0.0f, 0.0f);
        continue; // capture smearing mistakes
      }
      histos.fill(HIST("hNaNBookkeeping"), 0.0f, 1.0f); // ok!

      // Time associated to the mcParticle: collision time + smearing
      float trackTime = (eventCollisionTimeNS + gRandom->Gaus(0., timeResolutionNs)) * nsToMus;
      TrackType trackType = reconstructed ? TrackType::kRecoPrimary : TrackType::kGhostPrimary;
      if (reconstructed) {
        recoPrimaries.push_back(TrackAlice3{trackParCov, mcParticle.globalIndex(), trackTime, timeResolutionUs, isDecayDaughter, false, 0, nTrkHits, trackType});
      } else {
        ghostPrimaries.push_back(TrackAlice3{trackParCov, mcParticle.globalIndex(), trackTime, timeResolutionUs, isDecayDaughter, false, 0, nTrkHits, trackType});
      }
    }

    // Compute primary vertex with primary tracks
    o2::vertexing::PVertex primaryVertex;
    if (enablePrimaryVertexing && recoPrimaries.size() <= 2) {
      LOG(info) << "Not enough primary tracks (" << recoPrimaries.size() << ") to compute vertex, skipping vertexing and collision population.";
      return;
    }
    computeVertex(mcCollision, recoPrimaries, primaryVertex, icfg);
    getHist(TH1, histPath + "hPVz")->Fill(primaryVertex.getZ());
    // populate collisions
    tableCollisions(-1, // BC is irrelevant in synthetic MC tests for now, could be adjusted in future
                    primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ(),
                    primaryVertex.getSigmaX2(), primaryVertex.getSigmaXY(), primaryVertex.getSigmaY2(),
                    primaryVertex.getSigmaXZ(), primaryVertex.getSigmaYZ(), primaryVertex.getSigmaZ2(),
                    0, primaryVertex.getChi2(), primaryVertex.getNContributors(),
                    eventCollisionTimeNS, 0.f); // For the moment the event collision time is taken as the "GEANT" time, the computation of the event time is done a posteriori from the tracks in the OTF TOF PID task
    tableMcCollisionLabels(mcCollision.globalIndex(), 0);
    tableCollisionsAlice3(dNdEta);
    tableOTFLUTConfigId(icfg); // Populate OTF LUT configuration ID

    // Study V0s and cascades
    int trackTableOffset = tableStoredTracks.lastIndex();
    std::vector<TrackAlice3> tracksCascadeProngs;
    for (size_t iCasc = 0; iCasc < genCascades.size(); iCasc++) {
      auto genCasc = mcParticles.rawIteratorAt(genCascades[iCasc] - mcParticles.offset());
      studyCascade(trackTableOffset, genCasc, tracksCascadeProngs, primaryVertex, icfg, dNdEta, eventCollisionTimeNS);
      fillTracksInfo(tracksCascadeProngs, primaryVertex, icfg);
      trackTableOffset += tracksCascadeProngs.size(); // each cascade takes 4 tracks in the table (cascade + 3 daughters)
      tracksCascadeProngs.clear();
    }

    std::vector<TrackAlice3> tracksV0Daugs;
    for (size_t iV0 = 0; iV0 < genV0s.size(); iV0++) {
      auto genV0 = mcParticles.rawIteratorAt(genV0s[iV0] - mcParticles.offset());
      studyV0(trackTableOffset, genV0, tracksV0Daugs, icfg, dNdEta, eventCollisionTimeNS);
      fillTracksInfo(tracksV0Daugs, primaryVertex, icfg);
      trackTableOffset += tracksV0Daugs.size(); // each V0 takes 2 tracks in the table (2 daughters)
      tracksV0Daugs.clear();
    }

    // *+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*
    // populate tracks
    fillTracksInfo(recoPrimaries, primaryVertex, icfg);
    fillTracksInfo(ghostPrimaries, primaryVertex, icfg);

    // do bookkeeping of fastTracker tracking
    if (enableSecondarySmearing) {
      histos.fill(HIST("hCovMatOK"), 0.0f, fastTracker[icfg]->GetCovMatNotOK());
      histos.fill(HIST("hCovMatOK"), 1.0f, fastTracker[icfg]->GetCovMatOK());
    }
    if (doExtraQA) {
      histos.fill(HIST("hRecoVsSimMultiplicity"), multiplicityCounter, recoPrimaries.size());
      getHist(TH1, histPath + "hSimMultiplicity")->Fill(multiplicityCounter);
      getHist(TH1, histPath + "hRecoMultiplicity")->Fill(recoPrimaries.size());
    }

    LOG(debug) << " <- Finished processing OTF tracking with LUT configuration ID " << icfg;
  } // end process

  void processOnTheFly(aod::McCollision const& mcCollision, aod::McParticles const& mcParticles)
  {
    for (size_t icfg = 0; icfg < mSmearer.size(); ++icfg) {
      LOG(debug) << "  -> Processing OTF tracking with LUT configuration ID " << icfg;
      processWithLUTs(mcCollision, mcParticles, static_cast<int>(icfg));
    }
  }

  void processConfigurationDev(aod::McCollision const& mcCollision, aod::McPartWithDaus const& mcParticles, const int icfg)
  {
    const std::string histPath = "Configuration_" + std::to_string(icfg) + "/";
    tracksAlice3.clear();
    ghostTracksAlice3.clear();
    bcData.clear();
    cascadesAlice3.clear();
    v0sAlice3.clear();

    // generate collision time
    auto ir = irSampler.generateCollisionTime();
    const float eventCollisionTimeNS = ir.timeInBCNS;

    // First we compute the number of charged particles in the event
    float dNdEta{0.f};
    computeDNDEta(dNdEta, mcParticles, histPath);

    uint32_t multiplicityCounter = 0;
    // Now that the multiplicity is known, we can process the particles to smear them
    for (const auto& mcParticle : mcParticles) {
      const bool longLivedToBeHandled = std::find(longLivedHandledPDGs.begin(), longLivedHandledPDGs.end(), std::abs(mcParticle.pdgCode())) != longLivedHandledPDGs.end();
      const bool nucleiToBeHandled = std::find(nucleiPDGs.begin(), nucleiPDGs.end(), std::abs(mcParticle.pdgCode())) != nucleiPDGs.end();
      const bool pdgsToBeHandled = longLivedToBeHandled || (enableNucleiSmearing && nucleiToBeHandled);
      if (!pdgsToBeHandled) {
        continue;
      }

      if (std::fabs(mcParticle.eta()) > maxEta) {
        continue;
      }

      if (mcParticle.pt() < minPt) {
        continue;
      }

      if (enablePrimarySmearing) {
        getHist(TH1, histPath + "hPtGenerated")->Fill(mcParticle.pt());
        getHist(TH1, histPath + "hPhiGenerated")->Fill(mcParticle.phi());
        switch (std::abs(mcParticle.pdgCode())) {
          case kElectron:
            getHist(TH1, histPath + "hPtGeneratedEl")->Fill(mcParticle.pt());
            break;
          case kPiPlus:
            getHist(TH1, histPath + "hPtGeneratedPi")->Fill(mcParticle.pt());
            break;
          case kKPlus:
            getHist(TH1, histPath + "hPtGeneratedKa")->Fill(mcParticle.pt());
            break;
          case kProton:
            getHist(TH1, histPath + "hPtGeneratedPr")->Fill(mcParticle.pt());
            break;
        }
      }

      multiplicityCounter++;
      o2::track::TrackParCov trackParCov;
      const bool isDecayDaughter = (mcParticle.getProcess() == TMCProcess::kPDecay);
      const float time = (eventCollisionTimeNS + gRandom->Gaus(0., timeResolutionNs)) * nsToMus;

      bool reconstructed = false;
      int nTrkHits = 0;
      if (enablePrimarySmearing && mcParticle.isPrimary()) {
        o2::upgrade::convertMCParticleToO2Track(mcParticle, trackParCov, pdgDB);
        reconstructed = mSmearer[icfg]->smearTrack(trackParCov, mcParticle.pdgCode(), dNdEta);
        nTrkHits = fastTrackerSettings.minSiliconHits;
      } else if (enableSecondarySmearing) {
        o2::track::TrackParCov perfectTrackParCov;
        o2::upgrade::convertMCParticleToO2Track(mcParticle, perfectTrackParCov, pdgDB);
        perfectTrackParCov.setPID(pdgCodeToPID(mcParticle.pdgCode()));
        nTrkHits = fastTracker[icfg]->FastTrack(perfectTrackParCov, trackParCov, dNdEta);
        if (nTrkHits < fastTrackerSettings.minSiliconHits) {
          reconstructed = false;
        } else {
          reconstructed = true;
        }
      }

      if (!reconstructed && processUnreconstructedTracks) {
        continue; // failed to reconstruct track
      }

      if (std::isnan(trackParCov.getZ())) {
        histos.fill(HIST("hNaNBookkeeping"), 0.0f, 0.0f);
        continue; // capture smearing mistakes
      }

      histos.fill(HIST("hNaNBookkeeping"), 0.0f, 1.0f);
      if (enablePrimarySmearing) {
        getHist(TH1, histPath + "hPtReconstructed")->Fill(trackParCov.getPt());
        if (std::abs(mcParticle.pdgCode()) == kElectron)
          getHist(TH1, histPath + "hPtReconstructedEl")->Fill(trackParCov.getPt());
        if (std::abs(mcParticle.pdgCode()) == kPiPlus)
          getHist(TH1, histPath + "hPtReconstructedPi")->Fill(trackParCov.getPt());
        if (std::abs(mcParticle.pdgCode()) == kKPlus)
          getHist(TH1, histPath + "hPtReconstructedKa")->Fill(trackParCov.getPt());
        if (std::abs(mcParticle.pdgCode()) == kProton)
          getHist(TH1, histPath + "hPtReconstructedPr")->Fill(trackParCov.getPt());
      }

      if (reconstructed) {
        tracksAlice3.push_back(TrackAlice3{trackParCov, mcParticle.globalIndex(), time, timeResolutionUs, isDecayDaughter, false, 0, nTrkHits, kRecoPrimary});
      } else {
        ghostTracksAlice3.push_back(TrackAlice3{trackParCov, mcParticle.globalIndex(), time, timeResolutionUs, isDecayDaughter, false, 0, nTrkHits, kGhostPrimary});
      }
    }

    o2::vertexing::PVertex primaryVertex;
    if (enablePrimaryVertexing && recoPrimaries.size() <= 2) {
      LOG(info) << "Not enough primary tracks (" << recoPrimaries.size() << ") to compute vertex, skipping vertexing and collision population.";
      return;
    }
    computeVertex(mcCollision, tracksAlice3, primaryVertex, icfg);

    getHist(TH1, histPath + "hSimMultiplicity")->Fill(multiplicityCounter);
    getHist(TH1, histPath + "hRecoMultiplicity")->Fill(tracksAlice3.size());
    getHist(TH1, histPath + "hPVz")->Fill(primaryVertex.getZ());

    // *+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*
    // populate collisions
    tableCollisions(-1, // BC is irrelevant in synthetic MC tests for now, could be adjusted in future
                    primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ(),
                    primaryVertex.getSigmaX2(), primaryVertex.getSigmaXY(), primaryVertex.getSigmaY2(),
                    primaryVertex.getSigmaXZ(), primaryVertex.getSigmaYZ(), primaryVertex.getSigmaZ2(),
                    0, primaryVertex.getChi2(), primaryVertex.getNContributors(),
                    eventCollisionTimeNS, 0.f); // For the moment the event collision time is taken as the "GEANT" time, the computation of the event time is done a posteriori from the tracks in the OTF TOF PID task
    tableMcCollisionLabels(mcCollision.globalIndex(), 0);
    tableCollisionsAlice3(dNdEta);
    tableOTFLUTConfigId(icfg); // Populate OTF LUT configuration ID

    // *+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*
    // populate tracks
    fillTracksInfo(tracksAlice3, primaryVertex, icfg);
    fillTracksInfo(ghostTracksAlice3, primaryVertex, icfg);
  }

  void processDecayer(aod::McCollision const& mcCollision, aod::McPartWithDaus const& mcParticles)
  {
    for (size_t icfg = 0; icfg < mSmearer.size(); ++icfg) {
      processConfigurationDev(mcCollision, mcParticles, static_cast<int>(icfg));
    }
  }

  PROCESS_SWITCH(OnTheFlyTracker, processOnTheFly, "Enable default workflow", true);
  PROCESS_SWITCH(OnTheFlyTracker, processDecayer, "Enable experimental decayer workflow", false);
};

/// Extends TracksExtra if necessary
struct onTheFlyTrackerInitializer {
  Spawns<aod::TracksExtra> tableStoredTracksExtra;
  void init(InitContext const&) {}
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<OnTheFlyTracker>(cfgc),
    adaptAnalysisTask<onTheFlyTrackerInitializer>(cfgc)};
}
