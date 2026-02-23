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

#include "ALICE3/Core/DelphesO2TrackSmearer.h"
#include "ALICE3/Core/DetLayer.h"
#include "ALICE3/Core/FastTracker.h"
#include "ALICE3/Core/TrackUtilities.h"
#include "ALICE3/DataModel/OTFCollision.h"
#include "ALICE3/DataModel/OTFMCParticle.h"
#include "ALICE3/DataModel/OTFStrangeness.h"
#include "ALICE3/DataModel/collisionAlice3.h"
#include "ALICE3/DataModel/tracksAlice3.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DCAFitter/DCAFitterN.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DetectorsBase/Propagator.h>
#include <DetectorsVertexing/PVertexer.h>
#include <DetectorsVertexing/PVertexerHelpers.h>
#include <DetectorsVertexing/PVertexerParams.h>
#include <Field/MagneticField.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/StaticFor.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/DCA.h>
#include <ReconstructionDataFormats/PID.h>
#include <SimulationDataFormat/InteractionSampler.h>

#include <TGenPhaseSpace.h>
#include <TGeoGlobalMagField.h>
#include <TLorentzVector.h>
#include <TPDGCode.h>
#include <TRandom3.h>

#include <array>
#include <map>
#include <string>
#include <utility>
#include <vector>

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

struct OnTheFlyTracker {
  Produces<aod::Collisions> tableCollisions;
  Produces<aod::McCollisionLabels> tableMcCollisionLabels;
  Produces<aod::StoredTracks> tableStoredTracks;
  Produces<aod::TracksExtension> tableTracksExtension;
  Produces<aod::StoredTracksCov> tableStoredTracksCov;
  Produces<aod::TracksCovExtension> tableTracksCovExtension;
  Produces<aod::McTrackLabels> tableMcTrackLabels;
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
                int nTPCHitsInput = 0) : o2::track::TrackParCov(src),
                                         mcLabel{label},
                                         timeEst{time, timeError},
                                         isDecayDau(decayDauInput),
                                         isWeakDecayDau(weakDecayDauInput),
                                         isUsedInCascading(isUsedInCascadingInput),
                                         nSiliconHits(nSiliconHitsInput),
                                         nTPCHits(nTPCHitsInput) {}
    const TimeEst& getTimeMUS() const { return timeEst; }
    int64_t mcLabel;
    TimeEst timeEst; ///< time estimate in ns
    bool isDecayDau;
    bool isWeakDecayDau;
    int isUsedInCascading; // 0: not at all, 1: is a cascade, 2: is a bachelor, 3: is a pion, 4: is a proton
    int nSiliconHits;
    int nTPCHits;
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

  // necessary for particle charges
  Service<o2::framework::O2DatabasePDG> pdgDB;

  // for handling basic QA histograms if requested
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  std::map<std::string, HistPtr> histPointers;

  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;

  // Track smearer array, one per geometry
  std::vector<std::unique_ptr<o2::delphes::DelphesO2TrackSmearer>> mSmearer;

  // For processing and vertexing
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
  void init(o2::framework::InitContext& initContext)
  {
    LOG(info) << "Initializing OnTheFlyTracker task";
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setTimestamp(-1);
    mGeoContainer.init(initContext);

    const int nGeometries = mGeoContainer.getNumberOfConfigurations();
    mMagneticField = mGeoContainer.getFloatValue(0, "global", "magneticfield");
    for (int icfg = 0; icfg < nGeometries; ++icfg) {
      const std::string histPath = "Configuration_" + std::to_string(icfg) + "/";
      mSmearer.emplace_back(std::make_unique<o2::delphes::DelphesO2TrackSmearer>());
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
        insertHist(histPath + "hDeltaXPVRecoGen", "hDeltaXPVRecoGen;Delta X (reco - gen), cm", {kTH1D, {{axes.axisDeltaVtxCoord}}});
        insertHist(histPath + "hDeltaYPVRecoGen", "hDeltaYPVRecoGen;Delta Y (reco - gen), cm", {kTH1D, {{axes.axisDeltaVtxCoord}}});
        insertHist(histPath + "hDeltaZPVRecoGen", "hDeltaZPVRecoGen;Delta Z (reco - gen), cm", {kTH1D, {{axes.axisDeltaVtxCoord}}});
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

          insertHist(histPath + "hGenXi", "hGenXi", {kTH2F, {axes.axisDecayRadius, axes.axisMomentum}});
          insertHist(histPath + "hRecoXi", "hRecoXi", {kTH2F, {axes.axisDecayRadius, axes.axisMomentum}});

          insertHist(histPath + "hGenPiFromXi", "hGenPiFromXi", {kTH2F, {axes.axisDecayRadius, axes.axisMomentum}});
          insertHist(histPath + "hGenPiFromLa", "hGenPiFromLa", {kTH2F, {axes.axisDecayRadius, axes.axisMomentum}});
          insertHist(histPath + "hGenPrFromLa", "hGenPrFromLa", {kTH2F, {axes.axisDecayRadius, axes.axisMomentum}});
          insertHist(histPath + "hRecoPiFromXi", "hRecoPiFromXi", {kTH2F, {axes.axisDecayRadius, axes.axisMomentum}});
          insertHist(histPath + "hRecoPiFromLa", "hRecoPiFromLa", {kTH2F, {axes.axisDecayRadius, axes.axisMomentum}});
          insertHist(histPath + "hRecoPrFromLa", "hRecoPrFromLa", {kTH2F, {axes.axisDecayRadius, axes.axisMomentum}});

          // basic mass histograms to see if we're in business
          insertHist(histPath + "hMassLambda", "hMassLambda", {kTH1F, {{axes.axisLambdaMass}}});
          insertHist(histPath + "hMassXi", "hMassXi", {kTH1F, {{axes.axisXiMass}}});
          insertHist(histPath + "h2dMassXi", "h2dMassXi", {kTH2F, {{axes.axisXiMass}, {axes.axisMomentum}}});

          // OTF strangeness tracking QA
          insertHist(histPath + "hFoundVsFindable", "hFoundVsFindable", {kTH2F, {{10, -0.5f, 9.5f}, {10, -0.5f, 9.5f}}});

          insertHist(histPath + "h2dDCAxyCascade", "h2dDCAxyCascade", {kTH2F, {axes.axisMomentum, axes.axisDCA}});
          insertHist(histPath + "h2dDCAxyCascadeBachelor", "h2dDCAxyCascadeBachelor", {kTH2F, {axes.axisMomentum, axes.axisDCA}});
          insertHist(histPath + "h2dDCAxyCascadeNegative", "h2dDCAxyCascadeNegative", {kTH2F, {axes.axisMomentum, axes.axisDCA}});
          insertHist(histPath + "h2dDCAxyCascadePositive", "h2dDCAxyCascadePositive", {kTH2F, {axes.axisMomentum, axes.axisDCA}});

          insertHist(histPath + "h2dDCAzCascade", "h2dDCAzCascade", {kTH2F, {axes.axisMomentum, axes.axisDCA}});
          insertHist(histPath + "h2dDCAzCascadeBachelor", "h2dDCAzCascadeBachelor", {kTH2F, {axes.axisMomentum, axes.axisDCA}});
          insertHist(histPath + "h2dDCAzCascadeNegative", "h2dDCAzCascadeNegative", {kTH2F, {axes.axisMomentum, axes.axisDCA}});
          insertHist(histPath + "h2dDCAzCascadePositive", "h2dDCAzCascadePositive", {kTH2F, {axes.axisMomentum, axes.axisDCA}});

          insertHist(histPath + "h2dDeltaPtVsPt", "h2dDeltaPtVsPt;Gen p_{T};#Delta p_{T}", {kTH2F, {axes.axisMomentum, axes.axisDeltaPt}});
          insertHist(histPath + "h2dDeltaEtaVsPt", "h2dDeltaEtaVsPt;Gen p_{T};#Delta #eta", {kTH2F, {axes.axisMomentum, axes.axisDeltaEta}});

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
        histPointers.insert({v0histPath + "hFastTrackerQA", h});
        // K0s
        insertHist(v0histPath + "K0/hGen", "hGen", kTH2F, {axes.axisDecayRadius, axes.axisMomentum});
        insertHist(v0histPath + "K0/hReco", "hReco", kTH2F, {axes.axisDecayRadius, axes.axisMomentum});
        insertHist(v0histPath + "K0/hGenNegDaughterFromV0", "hGenNegDaughterFromV0", kTH2F, {axes.axisDecayRadius, axes.axisMomentum});
        insertHist(v0histPath + "K0/hGenPosDaughterFromV0", "hGenPosDaughterFromV0", kTH2F, {axes.axisDecayRadius, axes.axisMomentum});
        insertHist(v0histPath + "K0/hRecoNegDaughterFromV0", "hRecoNegDaughterFromV0", kTH2F, {axes.axisDecayRadius, axes.axisMomentum});
        insertHist(v0histPath + "K0/hRecoPosDaughterFromV0", "hRecoPosDaughterFromV0", kTH2F, {axes.axisDecayRadius, axes.axisMomentum});
        insertHist(v0histPath + "K0/hMass", "hMass", kTH2F, {axes.axisK0Mass, axes.axisMomentum});
        insertHist(v0histPath + "K0/hFinalMass", "hFinalMass", kTH1F, {axes.axisK0Mass});
        // Lambda
        insertHist(v0histPath + "Lambda/hGen", "hGen", kTH2F, {axes.axisDecayRadius, axes.axisMomentum});
        insertHist(v0histPath + "Lambda/hReco", "hReco", kTH2F, {axes.axisDecayRadius, axes.axisMomentum});
        insertHist(v0histPath + "Lambda/hGenNegDaughterFromV0", "hGenNegDaughterFromV0", kTH2F, {axes.axisDecayRadius, axes.axisMomentum});
        insertHist(v0histPath + "Lambda/hGenPosDaughterFromV0", "hGenPosDaughterFromV0", kTH2F, {axes.axisDecayRadius, axes.axisMomentum});
        insertHist(v0histPath + "Lambda/hRecoNegDaughterFromV0", "hRecoNegDaughterFromV0", kTH2F, {axes.axisDecayRadius, axes.axisMomentum});
        insertHist(v0histPath + "Lambda/hRecoPosDaughterFromV0", "hRecoPosDaughterFromV0", kTH2F, {axes.axisDecayRadius, axes.axisMomentum});
        insertHist(v0histPath + "Lambda/hMass", "hMass", kTH2F, {axes.axisLambdaMass, axes.axisMomentum});
        insertHist(v0histPath + "Lambda/hFinalMass", "hFinalMass", kTH1F, {axes.axisLambdaMass});

        // AntiLambda
        insertHist(v0histPath + "AntiLambda/hGen", "hGen", kTH2F, {axes.axisDecayRadius, axes.axisMomentum});
        insertHist(v0histPath + "AntiLambda/hReco", "hReco", kTH2F, {axes.axisDecayRadius, axes.axisMomentum});
        insertHist(v0histPath + "AntiLambda/hGenNegDaughterFromV0", "hGenNegDaughterFromV0", kTH2F, {axes.axisDecayRadius, axes.axisMomentum});
        insertHist(v0histPath + "AntiLambda/hGenPosDaughterFromV0", "hGenPosDaughterFromV0", kTH2F, {axes.axisDecayRadius, axes.axisMomentum});
        insertHist(v0histPath + "AntiLambda/hRecoNegDaughterFromV0", "hRecoNegDaughterFromV0", kTH2F, {axes.axisDecayRadius, axes.axisMomentum});
        insertHist(v0histPath + "AntiLambda/hRecoPosDaughterFromV0", "hRecoPosDaughterFromV0", kTH2F, {axes.axisDecayRadius, axes.axisMomentum});
        insertHist(v0histPath + "AntiLambda/hMass", "hMass", kTH2F, {axes.axisLambdaMass, axes.axisMomentum});
        insertHist(v0histPath + "AntiLambda/hFinalMass", "hFinalMass", kTH1F, {axes.axisLambdaMass});
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
    fitter.setPropagateToPCA(true);
    fitter.setMaxR(200.);
    fitter.setMinParamChange(1e-3);
    fitter.setMinRelChi2Change(0.9);
    fitter.setMaxDZIni(1e9);
    fitter.setMaxDXYIni(4);
    fitter.setMaxChi2(1e9);
    fitter.setUseAbsDCA(true);
    fitter.setWeightedFinalPCA(false);
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

  float dNdEta = 0.f;                                                      // Charged particle multiplicity to use in the efficiency evaluation
  std::pair<float, float> vertexReconstructionEfficiencyCounters = {0, 0}; // {nVerticesWithMoreThan2Contributors, nVerticesReconstructed}
  void processWithLUTs(aod::McCollision const& mcCollision, aod::McParticles const& mcParticles, const int icfg)
  {
    vertexReconstructionEfficiencyCounters.first += 1;
    const int lastTrackIndex = tableStoredTracksCov.lastIndex() + 1; // bookkeep the last added track
    const std::string histPath = "Configuration_" + std::to_string(icfg) + "/";

    tracksAlice3.clear();
    ghostTracksAlice3.clear();
    bcData.clear();
    cascadesAlice3.clear();
    v0sAlice3.clear();

    o2::dataformats::DCA dcaInfo;
    o2::dataformats::VertexBase vtx;

    // generate collision time
    auto ir = irSampler.generateCollisionTime();
    const float eventCollisionTimeNS = ir.timeInBCNS;

    // First we compute the number of charged particles in the event
    dNdEta = 0.f;
    LOG(debug) << "Processing " << mcParticles.size() << " MC particles to compute dNch/deta";
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
      const bool pdgsToBeHandled = longLivedToBeHandled || (enableNucleiSmearing && nucleiToBeHandled) || (cascadeDecaySettings.decayXi && mcParticle.pdgCode() == kXiMinus);
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
    uint32_t multiplicityCounter = 0;
    getHist(TH1, histPath + "hLUTMultiplicity")->Fill(dNdEta);

    // Now that the multiplicity is known, we can process the particles to smear them
    double xiDecayRadius2D = 0;
    double laDecayRadius2D = 0;
    double v0DecayRadius2D = 0;
    std::vector<TLorentzVector> cascadeDecayProducts;
    std::vector<TLorentzVector> v0DecayProducts;
    std::vector<double> xiDecayVertex, laDecayVertex, v0DecayVertex;
    for (const auto& mcParticle : mcParticles) {
      xiDecayRadius2D = 0;
      laDecayRadius2D = 0;
      v0DecayRadius2D = 0;
      xiDecayVertex.clear();
      laDecayVertex.clear();
      v0DecayVertex.clear();

      cascadeDecayProducts.clear();
      v0DecayProducts.clear();
      const bool isCascade = mcParticle.pdgCode() == kXiMinus;
      if (cascadeDecaySettings.decayXi && isCascade) {
        o2::track::TrackParCov xiTrackParCov;
        o2::upgrade::convertMCParticleToO2Track(mcParticle, xiTrackParCov, pdgDB);
        decayCascade(mcParticle, xiTrackParCov, cascadeDecayProducts, xiDecayVertex, laDecayVertex);
        if (cascadeDecayProducts.size() != 3) {
          LOG(fatal) << "Xi decay did not produce 3 daughters as expected!";
        }
        xiDecayRadius2D = std::hypot(xiDecayVertex[0], xiDecayVertex[1]);
        laDecayRadius2D = std::hypot(laDecayVertex[0], laDecayVertex[1]);
      }
      const bool isV0 = std::find(v0PDGs.begin(), v0PDGs.end(), mcParticle.pdgCode()) != v0PDGs.end();

      if (v0DecaySettings.decayV0 && isV0) {
        decayV0Particle(mcParticle, v0DecayProducts, v0DecayVertex, mcParticle.pdgCode());
        if (v0DecayProducts.size() != 2) {
          LOG(fatal) << "V0 decay did not produce 2 daughters as expected!";
        }
        v0DecayRadius2D = std::hypot(v0DecayVertex[0], v0DecayVertex[1]);
      }

      if (!mcParticle.isPhysicalPrimary()) {
        continue;
      }

      const bool longLivedToBeHandled = std::find(longLivedHandledPDGs.begin(), longLivedHandledPDGs.end(), std::abs(mcParticle.pdgCode())) != longLivedHandledPDGs.end();
      const bool nucleiToBeHandled = std::find(nucleiPDGs.begin(), nucleiPDGs.end(), std::abs(mcParticle.pdgCode())) != nucleiPDGs.end();
      const bool pdgsToBeHandled = longLivedToBeHandled || (enableNucleiSmearing && nucleiToBeHandled) || (cascadeDecaySettings.decayXi && isCascade) || (v0DecaySettings.decayV0 && isV0);
      if (!pdgsToBeHandled) {
        continue;
      }

      if (std::fabs(mcParticle.eta()) > maxEta) {
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
      if (cascadeDecaySettings.doXiQA && isCascade) {
        getHist(TH2, histPath + "hGenXi")->Fill(xiDecayRadius2D, mcParticle.pt());
        getHist(TH2, histPath + "hGenPiFromXi")->Fill(xiDecayRadius2D, cascadeDecayProducts[0].Pt());
        getHist(TH2, histPath + "hGenPiFromLa")->Fill(laDecayRadius2D, cascadeDecayProducts[1].Pt());
        getHist(TH2, histPath + "hGenPrFromLa")->Fill(laDecayRadius2D, cascadeDecayProducts[2].Pt());
      }
      if (v0DecaySettings.doV0QA && isV0) {
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
      if (mcParticle.pt() < minPt) {
        continue;
      }

      o2::track::TrackParCov trackParCov;
      o2::upgrade::convertMCParticleToO2Track(mcParticle, trackParCov, pdgDB);

      const bool isDecayDaughter = (mcParticle.getProcess() == TMCProcess::kPDecay);

      multiplicityCounter++;
      const float timeResolutionNs = 100.f; // ns
      const float nsToMus = 1e-3f;
      const float timeResolutionUs = timeResolutionNs * nsToMus; // us
      const float time = (eventCollisionTimeNS + gRandom->Gaus(0., timeResolutionNs)) * nsToMus;
      static constexpr int kCascProngs = 3;
      std::vector<o2::track::TrackParCov> xiDaughterTrackParCovsPerfect(3);
      std::vector<o2::track::TrackParCov> xiDaughterTrackParCovsTracked(3);
      std::vector<bool> isReco(kCascProngs);
      std::vector<int> nHits(kCascProngs);        // total
      std::vector<int> nSiliconHits(kCascProngs); // silicon type
      std::vector<int> nTPCHits(kCascProngs);     // TPC type

      bool tryKinkReco = false;
      if (cascadeDecaySettings.decayXi && isCascade) {
        bool reconstructedCascade = false;
        if (cascadeDecaySettings.doXiQA) {
          getHist(TH1, histPath + "hXiBuilding")->Fill(0.0f);
        }

        o2::upgrade::convertTLorentzVectorToO2Track(PDG_t::kPiMinus, cascadeDecayProducts[0], xiDecayVertex, xiDaughterTrackParCovsPerfect[0], pdgDB);
        xiDaughterTrackParCovsPerfect[0].setPID(pdgCodeToPID(PDG_t::kPiMinus));
        o2::upgrade::convertTLorentzVectorToO2Track(PDG_t::kPiMinus, cascadeDecayProducts[1], laDecayVertex, xiDaughterTrackParCovsPerfect[1], pdgDB);
        xiDaughterTrackParCovsPerfect[1].setPID(pdgCodeToPID(PDG_t::kPiMinus));
        o2::upgrade::convertTLorentzVectorToO2Track(PDG_t::kProton, cascadeDecayProducts[2], laDecayVertex, xiDaughterTrackParCovsPerfect[2], pdgDB);
        xiDaughterTrackParCovsPerfect[2].setPID(pdgCodeToPID(PDG_t::kProton));

        for (int i = 0; i < kCascProngs; i++) {
          isReco[i] = false;
          nHits[i] = 0;
          nSiliconHits[i] = 0;
          nTPCHits[i] = 0;
          if (enableSecondarySmearing) {
            nHits[i] = fastTracker[icfg]->FastTrack(xiDaughterTrackParCovsPerfect[i], xiDaughterTrackParCovsTracked[i], dNdEta);
            nSiliconHits[i] = fastTracker[icfg]->GetNSiliconPoints();
            nTPCHits[i] = fastTracker[icfg]->GetNGasPoints();

            if (nHits[i] < 0 && cascadeDecaySettings.doXiQA) { // QA
              getHist(TH1, histPath + "hFastTrackerQA")->Fill(o2::math_utils::abs(nHits[i]));
            }

            if (nSiliconHits[i] >= fastTrackerSettings.minSiliconHits || (nSiliconHits[i] >= fastTrackerSettings.minSiliconHitsIfTPCUsed && nTPCHits[i] >= fastTrackerSettings.minTPCClusters)) {
              isReco[i] = true;
            } else {
              continue; // extra sure
            }
            for (uint32_t ih = 0; ih < fastTracker[icfg]->GetNHits() && cascadeDecaySettings.doXiQA; ih++) {
              getHist(TH2, histPath + "hFastTrackerHits")->Fill(fastTracker[icfg]->GetHitZ(ih), std::hypot(fastTracker[icfg]->GetHitX(ih), fastTracker[icfg]->GetHitY(ih)));
            }
          } else {
            isReco[i] = true;
            xiDaughterTrackParCovsTracked[i] = xiDaughterTrackParCovsPerfect[i];
          }

          if (TMath::IsNaN(xiDaughterTrackParCovsTracked[i].getZ())) {
            continue;
          } else {
            histos.fill(HIST("hNaNBookkeeping"), i + 1, 1.0f);
          }
          if (isReco[i]) {
            tracksAlice3.push_back(TrackAlice3{xiDaughterTrackParCovsTracked[i], mcParticle.globalIndex(), time, timeResolutionUs, true, true, i + 2, nSiliconHits[i], nTPCHits[i]});
          } else {
            ghostTracksAlice3.push_back(TrackAlice3{xiDaughterTrackParCovsTracked[i], mcParticle.globalIndex(), time, timeResolutionUs, true, true, i + 2});
          }
        }

        if (!isReco[1] || !isReco[2]) {
          tryKinkReco = true; // Lambda outside acceptance, set flag for kink reco to be used if mode 1
        }

        if (isReco[0] && isReco[1] && isReco[2]) {
          reconstructedCascade = true;
        }

        // +-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+
        // combine particles into actual Xi candidate
        // cascade building starts here
        if (cascadeDecaySettings.findXi && isReco[0] && isReco[1] && isReco[2] && cascadeDecaySettings.doKinkReco != 2) {
          if (cascadeDecaySettings.doXiQA) {
            getHist(TH1, histPath + "hXiBuilding")->Fill(3.0f);
          }
          // assign indices of the particles we've used
          // they should be the last ones to be filled, in order:
          // n-1: proton from lambda
          // n-2: pion from lambda
          // n-3: pion from xi
          thisCascade.positiveId = lastTrackIndex + tracksAlice3.size() - 1;
          thisCascade.negativeId = lastTrackIndex + tracksAlice3.size() - 2;
          thisCascade.bachelorId = lastTrackIndex + tracksAlice3.size() - 3;

          // use DCA fitters
          int nCand = 0;
          bool dcaFitterOK_V0 = true;
          try {
            nCand = fitter.process(xiDaughterTrackParCovsTracked[1], xiDaughterTrackParCovsTracked[2]);
          } catch (...) {
            // LOG(error) << "Exception caught in DCA fitter process call!";
            dcaFitterOK_V0 = false;
          }
          if (nCand == 0) {
            dcaFitterOK_V0 = false;
          }
          // V0 found successfully
          if (dcaFitterOK_V0) {
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
            bool dcaFitterOK_Cascade = true;
            try {
              nCand = fitter.process(v0Track, xiDaughterTrackParCovsTracked[0]);
            } catch (...) {
              // LOG(error) << "Exception caught in DCA fitter process call!";
              dcaFitterOK_Cascade = false;
            }
            if (nCand == 0) {
              dcaFitterOK_Cascade = false;
            }

            // Cascade found successfully
            if (dcaFitterOK_Cascade) {
              if (cascadeDecaySettings.doXiQA) {
                getHist(TH1, histPath + "hXiBuilding")->Fill(5.0f);
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
                  trackParCov.getXatLabR(layer.getRadius(), targetX, mMagneticField);
                  if (targetX > 999)
                    continue; // failed to find intercept

                  if (!trackParCov.propagateTo(targetX, mMagneticField)) {
                    continue; // failed to propagate
                  }

                  // get potential cluster position
                  std::array<float, 3> posClusterCandidate;
                  trackParCov.getXYZGlo(posClusterCandidate);
                  float r{std::hypot(posClusterCandidate[0], posClusterCandidate[1])};
                  float phi{std::atan2(-posClusterCandidate[1], -posClusterCandidate[0]) + o2::constants::math::PI};

                  if (layer.getResolutionRPhi() > 1e-8 && layer.getResolutionZ() > 1e-8) { // catch zero (though should not really happen...)
                    phi = gRandom->Gaus(phi, std::asin(layer.getResolutionRPhi() / r));
                    posClusterCandidate[0] = r * std::cos(phi);
                    posClusterCandidate[1] = r * std::sin(phi);
                    posClusterCandidate[2] = gRandom->Gaus(posClusterCandidate[2], layer.getResolutionZ());
                  }

                  if (std::isnan(phi))
                    continue; // Catch when getXatLabR misses layer[i]

                  // towards adding cluster: move to track alpha
                  double alpha = cascadeTrack.getAlpha();
                  double xyz1[3]{
                    TMath::Cos(alpha) * posClusterCandidate[0] + TMath::Sin(alpha) * posClusterCandidate[1],
                    -TMath::Sin(alpha) * posClusterCandidate[0] + TMath::Cos(alpha) * posClusterCandidate[1],
                    posClusterCandidate[2]};

                  if (!(cascadeTrack.propagateTo(xyz1[0], mMagneticField)))
                    continue;
                  const o2::track::TrackParametrization<float>::dim2_t hitpoint = {static_cast<float>(xyz1[1]), static_cast<float>(xyz1[2])};
                  const o2::track::TrackParametrization<float>::dim3_t hitpointcov = {layer.getResolutionRPhi() * layer.getResolutionRPhi(), 0.f, layer.getResolutionZ() * layer.getResolutionZ()};
                  if (layer.isInDeadPhiRegion(phi)) {
                    continue; // No hit for strangeness tracking update
                  }

                  cascadeTrack.update(hitpoint, hitpointcov);
                  thisCascade.foundClusters++; // add to findable
                }
              }

              if (thisCascade.foundClusters < cascadeDecaySettings.minStraTrackHits) {
                continue; // We didn't find enough hits for strangeness tracking
              }

              // add cascade track
              thisCascade.cascadeTrackId = lastTrackIndex + tracksAlice3.size(); // this is the next index to be filled -> should be it
              tracksAlice3.push_back(TrackAlice3{cascadeTrack, mcParticle.globalIndex(), time, timeResolutionUs, false, false, 1, thisCascade.foundClusters});

              // add this cascade to vector (will fill cursor later with collision ID)
              cascadesAlice3.push_back(thisCascade);
            }
          }
        } // end cascade building

        if (isReco[0] && ((cascadeDecaySettings.doKinkReco == 1 && tryKinkReco) || cascadeDecaySettings.doKinkReco == 2)) { // mode 1 or 2
          o2::track::TrackParCov prefectCascadeTrack, trackedCasc;
          const o2::track::TrackParCov& trackedBach = xiDaughterTrackParCovsTracked[0];
          o2::upgrade::convertMCParticleToO2Track(mcParticle, prefectCascadeTrack, pdgDB);

          // back track is already smeared
          prefectCascadeTrack.setPID(pdgCodeToPID(PDG_t::kXiMinus)); // FIXME: not OK for omegas
          int nCascHits = fastTracker[icfg]->FastTrack(prefectCascadeTrack, trackedCasc, dNdEta);
          reconstructedCascade = (fastTrackerSettings.minSiliconHitsForKinkReco < nCascHits) ? false : true;

          if (reconstructedCascade) {
            std::array<float, 3> pCasc;
            std::array<float, 3> pBach;
            std::array<float, 3> pV0;
            trackedCasc.getPxPyPzGlo(pCasc);
            trackedBach.getPxPyPzGlo(pBach);
            for (size_t i = 0; i < pCasc.size(); ++i) {
              pV0[i] = pCasc[i] - pBach[i];
            }

            if (isReco[1] && !isReco[2]) {
              thisCascade.negativeId = lastTrackIndex + tracksAlice3.size() - 1;
              thisCascade.positiveId = -1;
            } else if (!isReco[1] && isReco[2]) {
              thisCascade.negativeId = -1;
              thisCascade.positiveId = lastTrackIndex + tracksAlice3.size() - 1;
            } else if (isReco[1] && isReco[2]) {
              thisCascade.positiveId = lastTrackIndex + tracksAlice3.size() - 1;
              thisCascade.negativeId = lastTrackIndex + tracksAlice3.size() - 2;
            } else {
              thisCascade.positiveId = -1;
              thisCascade.negativeId = -1;
            }

            int nCand = 0;
            bool kinkFitterOK = true;
            try {
              nCand = fitter.process(trackedCasc, trackedBach);
            } catch (...) {
              kinkFitterOK = false;
            }

            if (nCand == 0) {
              kinkFitterOK = false;
            }

            if (kinkFitterOK) {
              if (cascadeDecaySettings.doXiQA) {
                getHist(TH1, histPath + "hXiBuilding")->Fill(6.0f);
              }
            }

            fitter.propagateTracksToVertex(); // propagate e and K to D vertex
            if (!fitter.isPropagateTracksToVertexDone()) {
              kinkFitterOK = false;
            }

            o2::track::TrackParCov newCascadeTrack = fitter.getTrack(0); // (cascade)
            std::array<float, 3> kinkVtx = {-999, -999, -999};
            kinkVtx = fitter.getPCACandidatePos();

            thisCascade.bachelorId = lastTrackIndex + tracksAlice3.size() - isReco.size();
            thisCascade.cascadeTrackId = lastTrackIndex + tracksAlice3.size(); // this should be ok
            thisCascade.dcaV0dau = -1.f;                                       // unknown
            thisCascade.v0radius = -1.f;                                       // unknown
            thisCascade.dcacascdau = std::sqrt(fitter.getChi2AtPCACandidate());
            thisCascade.cascradius = std::hypot(kinkVtx[0], kinkVtx[1]);
            thisCascade.cascradiusMC = xiDecayRadius2D;
            thisCascade.mLambda = o2::constants::physics::MassLambda;
            thisCascade.findableClusters = nCascHits;
            thisCascade.foundClusters = nCascHits;
            thisCascade.pt = newCascadeTrack.getPt();
            thisCascade.eta = newCascadeTrack.getEta();
            thisCascade.mXi = RecoDecay::m(std::array{std::array{pBach[0], pBach[1], pBach[2]},
                                                      std::array{pV0[0], pV0[1], pV0[2]}},
                                           std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassLambda});

            newCascadeTrack.setPID(pdgCodeToPID(PDG_t::kXiMinus)); // FIXME: not OK for omegas
            tracksAlice3.push_back(TrackAlice3{newCascadeTrack, mcParticle.globalIndex(), time, timeResolutionUs, false, false, 1, thisCascade.foundClusters});

            // add this cascade to vector (will fill cursor later with collision ID)
            cascadesAlice3.push_back(thisCascade);
          }
        } // end cascade kink building

        // +-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+
        if (cascadeDecaySettings.doXiQA) {
          if (reconstructedCascade) {
            getHist(TH2, histPath + "hRecoXi")->Fill(xiDecayRadius2D, mcParticle.pt());
            getHist(TH1, histPath + "hMassLambda")->Fill(thisCascade.mLambda);
            getHist(TH1, histPath + "hMassXi")->Fill(thisCascade.mXi);
            getHist(TH2, histPath + "h2dMassXi")->Fill(thisCascade.mXi, thisCascade.pt);
            getHist(TH2, histPath + "h2dDeltaPtVsPt")->Fill(thisCascade.pt, mcParticle.pt() - thisCascade.pt);
            getHist(TH2, histPath + "h2dDeltaEtaVsPt")->Fill(thisCascade.pt, mcParticle.eta() - thisCascade.eta);
            getHist(TH2, histPath + "hFoundVsFindable")->Fill(thisCascade.findableClusters, thisCascade.foundClusters);
          }
          if (isReco[0]) {
            getHist(TH2, histPath + "hRecoPiFromXi")->Fill(xiDecayRadius2D, cascadeDecayProducts[0].Pt());
          }
          if (isReco[1]) {
            getHist(TH2, histPath + "hRecoPiFromLa")->Fill(laDecayRadius2D, cascadeDecayProducts[1].Pt());
          }
          if (isReco[2]) {
            getHist(TH2, histPath + "hRecoPrFromLa")->Fill(laDecayRadius2D, cascadeDecayProducts[2].Pt());
          }
        }

        continue; // Cascade handling done, should not be considered anymore
      }

      // V0 handling
      std::vector<o2::track::TrackParCov> v0DaughterTrackParCovsPerfect(2);
      std::vector<o2::track::TrackParCov> v0DaughterTrackParCovsTracked(2);
      std::vector<bool> isV0Reco(kv0Prongs);
      std::vector<int> nV0Hits(kv0Prongs);        // total
      std::vector<int> nV0SiliconHits(kv0Prongs); // silicon type
      std::vector<int> nV0TPCHits(kv0Prongs);     // TPC type
      if (v0DecaySettings.decayV0 && isV0) {
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
        for (int i = 0; i < kv0Prongs; i++) {
          isV0Reco[i] = false;
          nV0Hits[i] = 0;
          nV0SiliconHits[i] = 0;
          nV0TPCHits[i] = 0;
          if (enableSecondarySmearing) {
            nV0Hits[i] = fastTracker[icfg]->FastTrack(v0DaughterTrackParCovsPerfect[i], v0DaughterTrackParCovsTracked[i], dNdEta);
            nV0SiliconHits[i] = fastTracker[icfg]->GetNSiliconPoints();
            nV0TPCHits[i] = fastTracker[icfg]->GetNGasPoints();

            if (nV0Hits[i] < 0 && v0DecaySettings.doV0QA) { // QA
              fillHist(TH1, Form("V0Building_Configuration_%i/hFastTrackerQA", icfg), o2::math_utils::abs(nV0Hits[i]));
            }

            if (nV0SiliconHits[i] >= fastTrackerSettings.minSiliconHits || (nV0SiliconHits[i] >= fastTrackerSettings.minSiliconHitsIfTPCUsed && nV0TPCHits[i] >= fastTrackerSettings.minTPCClusters)) {
              isReco[i] = true;
            } else {
              continue; // extra sure
            }
            for (uint32_t ih = 0; ih < fastTracker[icfg]->GetNHits() && v0DecaySettings.doV0QA; ih++) {
              fillHist(TH2, Form("V0Building_Configuration_%i/hFastTrackerHits", icfg), fastTracker[icfg]->GetHitZ(ih), std::hypot(fastTracker[icfg]->GetHitX(ih), fastTracker[icfg]->GetHitY(ih)));
            }
          } else {
            isReco[i] = true;
            v0DaughterTrackParCovsTracked[i] = v0DaughterTrackParCovsPerfect[i];
          }

          // if (TMath::IsNaN(v0DaughterTrackParCovsTracked[i].getZ())) {
          //   continue;
          // } else {
          //   histos.fill(HIST("hNaNBookkeeping"), i + 1, 1.0f);
          // }
          if (isReco[i]) {
            tracksAlice3.push_back(TrackAlice3{v0DaughterTrackParCovsTracked[i], mcParticle.globalIndex(), time, timeResolutionUs, true, true, i + 2, nSiliconHits[i], nTPCHits[i]});
          } else {
            ghostTracksAlice3.push_back(TrackAlice3{v0DaughterTrackParCovsTracked[i], mcParticle.globalIndex(), time, timeResolutionUs, true, true, i + 2});
          }
        }
        if (v0DecaySettings.doV0QA) {
          if (isReco[0] && isReco[1]) {
            fillHist(TH1, Form("V0Building_Configuration_%i/hV0Building", icfg), 1.0f);
            for (size_t indexV0 = 0; indexV0 < v0PDGs.size(); indexV0++) {
              if (mcParticle.pdgCode() == v0PDGs[indexV0]) {
                fillHist(TH2, Form("V0Building_Configuration_%i/%s/hReco", icfg, NameV0s[indexV0].data()), v0DecayRadius2D, mcParticle.pt());
              }
            }
          }
          if (isReco[0]) {
            for (size_t indexV0 = 0; indexV0 < v0PDGs.size(); indexV0++) {
              if (mcParticle.pdgCode() == v0PDGs[indexV0]) {
                fillHist(TH2, Form("V0Building_Configuration_%i/%s/hRecoNegDaughterFromV0", icfg, NameV0s[indexV0].data()), v0DecayRadius2D, v0DecayProducts[0].Pt());
              }
            }
          }
          if (isReco[1]) {
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
        if (v0DecaySettings.findV0 && isReco[0] && isReco[1]) {
          if (v0DecaySettings.doV0QA) {
            fillHist(TH1, Form("V0Building_Configuration_%i/hV0Building", icfg), 2.0f);
          }

          // assign indices of the particles we've used
          // they should be the last ones to be filled, in order:
          // n-1: positive Track from V0
          // n-2: negative Track from V0
          thisV0.positiveId = lastTrackIndex + tracksAlice3.size() - 1;
          thisV0.negativeId = lastTrackIndex + tracksAlice3.size() - 2;
          thisV0.mcParticleId = mcParticle.globalIndex();
          // use DCA fitters
          int nCand = 0;
          bool dcaFitterOK_V0 = true;
          try {
            nCand = fitter.process(v0DaughterTrackParCovsTracked[0], v0DaughterTrackParCovsTracked[1]);
          } catch (...) {
            // LOG(error) << "Exception caught in DCA fitter process call!";
            dcaFitterOK_V0 = false;
          }
          if (nCand == 0) {
            dcaFitterOK_V0 = false;
          }
          // V0 found successfully
          if (dcaFitterOK_V0) {
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
            if (std::abs(mcParticle.pdgCode()) == kK0Short) {
              thisV0.mK0 = RecoDecay::m(std::array{std::array{posP[0], posP[1], posP[2]},
                                                   std::array{negP[0], negP[1], negP[2]}},
                                        std::array{o2::constants::physics::MassPionCharged,
                                                   o2::constants::physics::MassPionCharged});
            } else {
              thisV0.mK0 = -1;
            }

            if (mcParticle.pdgCode() == kLambda0) {
              thisV0.mLambda = RecoDecay::m(std::array{std::array{posP[0], posP[1], posP[2]},
                                                       std::array{negP[0], negP[1], negP[2]}},
                                            std::array{o2::constants::physics::MassProton,
                                                       o2::constants::physics::MassPionCharged});
            } else {
              thisV0.mLambda = -1;
            }

            if (mcParticle.pdgCode() == kLambda0Bar) {
              thisV0.mAntiLambda = RecoDecay::m(std::array{std::array{posP[0], posP[1], posP[2]},
                                                           std::array{negP[0], negP[1], negP[2]}},
                                                std::array{o2::constants::physics::MassPionCharged,
                                                           o2::constants::physics::MassProton});
            } else {
              thisV0.mAntiLambda = -1;
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

            // add this V0 to vector (will fill cursor later with collision ID)
            v0sAlice3.push_back(thisV0);
          }
        }
      }

      if (doExtraQA) {
        histos.fill(HIST("hSimTrackX"), trackParCov.getX());
      }
      if (isV0) {
        continue; // V0 handling done, should not be considered anymore
      }

      bool reconstructed = true;
      if (enablePrimarySmearing && !fastPrimaryTrackerSettings.fastTrackPrimaries) {
        reconstructed = mSmearer[icfg]->smearTrack(trackParCov, mcParticle.pdgCode(), dNdEta);
      } else if (fastPrimaryTrackerSettings.fastTrackPrimaries) {
        o2::track::TrackParCov o2Track;
        o2::upgrade::convertMCParticleToO2Track(mcParticle, o2Track, pdgDB);
        o2Track.setPID(pdgCodeToPID(mcParticle.pdgCode()));
        const int nHits = fastTracker[icfg]->FastTrack(o2Track, trackParCov, dNdEta);
        if (nHits < fastPrimaryTrackerSettings.minSiliconHits) {
          reconstructed = false;
        }
      }

      if (!reconstructed && !processUnreconstructedTracks) {
        continue;
      }
      if (TMath::IsNaN(trackParCov.getZ())) {
        // capture rare smearing mistakes / corrupted tracks
        histos.fill(HIST("hNaNBookkeeping"), 0.0f, 0.0f);
        continue;
      } else {
        histos.fill(HIST("hNaNBookkeeping"), 0.0f, 1.0f); // ok!
      }

      // Base QA (note: reco pT here)
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
      if (doExtraQA) {
        getHist(TH2, histPath + "h2dPtRes")->Fill(trackParCov.getPt(), (trackParCov.getPt() - mcParticle.pt()) / trackParCov.getPt());
        histos.fill(HIST("hRecoTrackX"), trackParCov.getX());
      }

      // populate vector with track if we reco-ed it
      if (reconstructed) {
        tracksAlice3.push_back(TrackAlice3{trackParCov, mcParticle.globalIndex(), time, timeResolutionUs, isDecayDaughter});
      } else {
        ghostTracksAlice3.push_back(TrackAlice3{trackParCov, mcParticle.globalIndex(), time, timeResolutionUs, isDecayDaughter});
      }
    }

    // *+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*
    // Calculate primary vertex with tracks from this collision
    // data preparation
    o2::vertexing::PVertex primaryVertex;
    if (enablePrimaryVertexing) {
      LOG(debug) << "Starting primary vertexing with " << tracksAlice3.size() << " tracks.";
      fillHist(TH1, histPath + "hVtxMultGen", tracksAlice3.size());
      std::vector<o2::MCCompLabel> lblTracks;
      std::vector<o2::vertexing::PVertex> vertices;
      std::vector<o2::vertexing::GIndex> vertexTrackIDs;
      std::vector<o2::vertexing::V2TRef> v2tRefs;
      std::vector<o2::MCEventLabel> lblVtx;
      lblVtx.emplace_back(mcCollision.globalIndex(), 1);
      std::vector<o2::dataformats::GlobalTrackID> idxVec; // store IDs

      idxVec.reserve(tracksAlice3.size());
      for (unsigned i = 0; i < tracksAlice3.size(); i++) {
        lblTracks.emplace_back(tracksAlice3[i].mcLabel, mcCollision.globalIndex(), 1, false);
        idxVec.emplace_back(i, o2::dataformats::GlobalTrackID::ITS); // let's say ITS
      }

      getHist(TH1, histPath + "hVtxTrials")->Fill(0); // Tried vertexing

      // Calculate vertices
      const int n_vertices = vertexer.process(tracksAlice3, // track array
                                              idxVec,
                                              gsl::span<o2::InteractionRecord>{bcData},
                                              vertices,
                                              vertexTrackIDs,
                                              v2tRefs,
                                              gsl::span<const o2::MCCompLabel>{lblTracks},
                                              lblVtx);

      LOG(debug) << "Vertex reconstruction efficiency " << vertexReconstructionEfficiencyCounters.second << "/" << vertexReconstructionEfficiencyCounters.first << "=" << vertexReconstructionEfficiencyCounters.second / vertexReconstructionEfficiencyCounters.first;
      if (n_vertices < 1) {
        LOG(debug) << "Vertexing completed, no vtx found.";
        if (doExtraQA) {
          histos.fill(HIST("h1dVerticesNotReco"), tracksAlice3.size());
        }
        return; // primary vertex not reconstructed
      }
      vertexReconstructionEfficiencyCounters.second += 1;
      getHist(TH1, histPath + "hVtxTrials")->Fill(1); // Succeeded vertexing
      LOG(debug) << "Vertexing completed with " << n_vertices << " vertices found.";

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
      fillHist(TH1, histPath + "hDeltaMultPVRecoGen", static_cast<int>(primaryVertex.getNContributors()) - static_cast<int>(tracksAlice3.size()));
      fillHist(TH1, histPath + "hDeltaXPVRecoGen", primaryVertex.getX() - mcCollision.posX());
      fillHist(TH1, histPath + "hDeltaYPVRecoGen", primaryVertex.getY() - mcCollision.posY());
      fillHist(TH1, histPath + "hDeltaZPVRecoGen", primaryVertex.getZ() - mcCollision.posZ());
    } else {
      primaryVertex.setXYZ(mcCollision.posX(), mcCollision.posY(), mcCollision.posZ());
    }
    // *+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*

    // debug / informational
    getHist(TH1, histPath + "hSimMultiplicity")->Fill(multiplicityCounter);
    getHist(TH1, histPath + "hRecoMultiplicity")->Fill(tracksAlice3.size());
    getHist(TH1, histPath + "hPVz")->Fill(primaryVertex.getZ());

    if (doExtraQA) {
      histos.fill(HIST("hRecoVsSimMultiplicity"), multiplicityCounter, tracksAlice3.size());
    }

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

    // *+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*
    // populate tracks
    LOG(debug) << "Populating " << tracksAlice3.size() << " tracks.";
    for (const auto& trackParCov : tracksAlice3) {
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
        if (cascadeDecaySettings.doXiQA) {
          if (trackParCov.isUsedInCascading == 1) {
            getHist(TH2, histPath + "h2dDCAxyCascade")->Fill(trackParametrization.getPt(), dcaXY * 1e+4); // in microns, please
            getHist(TH2, histPath + "h2dDCAzCascade")->Fill(trackParametrization.getPt(), dcaZ * 1e+4);   // in microns, please
          }
          if (trackParCov.isUsedInCascading == 2) {
            getHist(TH2, histPath + "h2dDCAxyCascadeBachelor")->Fill(trackParametrization.getPt(), dcaXY * 1e+4); // in microns, please
            getHist(TH2, histPath + "h2dDCAzCascadeBachelor")->Fill(trackParametrization.getPt(), dcaZ * 1e+4);   // in microns, please
          }
          if (trackParCov.isUsedInCascading == 3) {
            getHist(TH2, histPath + "h2dDCAxyCascadeNegative")->Fill(trackParametrization.getPt(), dcaXY * 1e+4); // in microns, please
            getHist(TH2, histPath + "h2dDCAzCascadeNegative")->Fill(trackParametrization.getPt(), dcaZ * 1e+4);   // in microns, please
          }
          if (trackParCov.isUsedInCascading == 4) {
            getHist(TH2, histPath + "h2dDCAxyCascadePositive")->Fill(trackParametrization.getPt(), dcaXY * 1e+4); // in microns, please
            getHist(TH2, histPath + "h2dDCAzCascadePositive")->Fill(trackParametrization.getPt(), dcaZ * 1e+4);   // in microns, please
          }
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
      tableTracksExtraA3(trackParCov.nSiliconHits, trackParCov.nTPCHits);

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

    // populate ghost tracks
    LOG(debug) << "Populating " << ghostTracksAlice3.size() << " ghost tracks.";
    for (const auto& trackParCov : ghostTracksAlice3) {
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
          histos.fill(HIST("h2dDCAxy"), trackParametrization.getPt(), dcaXY * 1e+4); // in microns, please
          histos.fill(HIST("h2dDCAz"), trackParametrization.getPt(), dcaZ * 1e+4);   // in microns, please
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
      tableTracksExtraA3(trackParCov.nSiliconHits, trackParCov.nTPCHits);

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
      tableTracksAlice3(false);
    }

    // populate Cascades
    LOG(debug) << "Populating " << cascadesAlice3.size() << " cascades.";
    for (const auto& cascade : cascadesAlice3) {
      tableUpgradeCascades(tableCollisions.lastIndex(), // now we know the collision index -> populate table
                           cascade.cascadeTrackId,
                           cascade.positiveId,
                           cascade.negativeId,
                           cascade.bachelorId,
                           cascade.dcaV0dau,
                           cascade.dcacascdau,
                           cascade.v0radius,
                           cascade.cascradius,
                           cascade.cascradiusMC,
                           cascade.mLambda,
                           cascade.mXi,
                           cascade.findableClusters,
                           cascade.foundClusters);
    }

    // populate V0s
    LOG(debug) << "Populating " << v0sAlice3.size() << " V0s.";
    for (const auto& v0 : v0sAlice3) {
      tableUpgradeV0s(tableCollisions.lastIndex(), // now we know the collision index -> populate table
                      v0.mcParticleId,
                      v0.positiveId,
                      v0.negativeId,
                      v0.dcaV0dau,
                      v0.v0radius,
                      v0.mLambda,
                      v0.mAntiLambda,
                      v0.mK0,
                      v0.pt);
    }

    // do bookkeeping of fastTracker tracking
    if (enableSecondarySmearing) {
      histos.fill(HIST("hCovMatOK"), 0.0f, fastTracker[icfg]->GetCovMatNotOK());
      histos.fill(HIST("hCovMatOK"), 1.0f, fastTracker[icfg]->GetCovMatOK());
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

  void processConfigurationDev(aod::McCollision const& mcCollision, aod::McPartsWithDau const& mcParticles, const int icfg)
  {
    // const int lastTrackIndex = tableStoredTracksCov.lastIndex() + 1; // bookkeep the last added track
    const std::string histPath = "Configuration_" + std::to_string(icfg) + "/";
    tracksAlice3.clear();
    ghostTracksAlice3.clear();
    bcData.clear();
    cascadesAlice3.clear();
    v0sAlice3.clear();

    o2::dataformats::DCA dcaInfo;
    o2::dataformats::VertexBase vtx;

    // generate collision time
    auto ir = irSampler.generateCollisionTime();
    const float eventCollisionTimeNS = ir.timeInBCNS;

    // First we compute the number of charged particles in the event
    dNdEta = 0.f;
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

    dNdEta /= (multEtaRange * 2.0f);
    uint32_t multiplicityCounter = 0;
    getHist(TH1, histPath + "hLUTMultiplicity")->Fill(dNdEta);

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
      const float timeResolutionNs = 100.f; // ns
      const float nsToMus = 1e-3f;
      const float timeResolutionUs = timeResolutionNs * nsToMus; // us
      const float time = (eventCollisionTimeNS + gRandom->Gaus(0., timeResolutionNs)) * nsToMus;

      bool reconstructed = false;
      if (enablePrimarySmearing && mcParticle.isPrimary()) {
        o2::upgrade::convertMCParticleToO2Track(mcParticle, trackParCov, pdgDB);
        reconstructed = mSmearer[icfg]->smearTrack(trackParCov, mcParticle.pdgCode(), dNdEta);
      } else if (enableSecondarySmearing) {
        o2::track::TrackParCov perfectTrackParCov;
        o2::upgrade::convertMCParticleToO2Track(mcParticle, perfectTrackParCov, pdgDB);
        perfectTrackParCov.setPID(pdgCodeToPID(mcParticle.pdgCode()));
        const int nHits = fastTracker[icfg]->FastTrack(perfectTrackParCov, trackParCov, dNdEta);
        if (nHits < fastTrackerSettings.minSiliconHits) {
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
        tracksAlice3.push_back(TrackAlice3{trackParCov, mcParticle.globalIndex(), time, timeResolutionUs, isDecayDaughter});
      } else {
        ghostTracksAlice3.push_back(TrackAlice3{trackParCov, mcParticle.globalIndex(), time, timeResolutionUs, isDecayDaughter});
      }
    }

    o2::vertexing::PVertex primaryVertex;
    if (enablePrimaryVertexing) { // disabled for now
      primaryVertex.setXYZ(mcCollision.posX(), mcCollision.posY(), mcCollision.posZ());
    } else {
      primaryVertex.setXYZ(mcCollision.posX(), mcCollision.posY(), mcCollision.posZ());
    }

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
    for (const auto& trackParCov : tracksAlice3) {
      aod::track::TrackTypeEnum trackType = aod::track::Track;

      if (populateTracksDCA) {
        float dcaXY = 1e+10, dcaZ = 1e+10;
        o2::track::TrackParCov trackParametrization(trackParCov);
        if (trackParametrization.propagateToDCA(primaryVertex, mMagneticField, &dcaInfo)) {
          dcaXY = dcaInfo.getY();
          dcaZ = dcaInfo.getZ();
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
      tableTracksExtraA3(trackParCov.nSiliconHits, trackParCov.nTPCHits);

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

    // *+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*
    // populate ghost tracks
    for (const auto& trackParCov : ghostTracksAlice3) {
      aod::track::TrackTypeEnum trackType = aod::track::Track;

      if (populateTracksDCA) {
        float dcaXY = 1e+10, dcaZ = 1e+10;
        o2::track::TrackParCov trackParametrization(trackParCov);
        if (trackParametrization.propagateToDCA(primaryVertex, mMagneticField, &dcaInfo)) {
          dcaXY = dcaInfo.getY();
          dcaZ = dcaInfo.getZ();
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
      tableTracksExtraA3(trackParCov.nSiliconHits, trackParCov.nTPCHits);

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
      tableTracksAlice3(false);
    }
  }

  void processDecayer(aod::McCollision const& mcCollision, aod::McPartsWithDau const& mcParticles)
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
