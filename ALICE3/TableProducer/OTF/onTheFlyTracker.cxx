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
#include "ALICE3/DataModel/OTFStrangeness.h"
#include "ALICE3/DataModel/OTFTracks.h"
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
#include <Field/MagneticField.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/StaticFor.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/DCA.h>
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
#define getHist(type, name) std::get<std::shared_ptr<type>>(histPointers[name])

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
  Configurable<float> magneticField{"magneticField", 20.0f, "magnetic field in kG"};
  Configurable<float> maxEta{"maxEta", 1.5, "maximum eta to consider viable"};
  Configurable<float> multEtaRange{"multEtaRange", 0.8, "eta range to compute the multiplicity"};
  Configurable<float> minPt{"minPt", 0.1, "minimum pt to consider viable"};
  Configurable<bool> enablePrimarySmearing{"enablePrimarySmearing", false, "Enable smearing of primary particles"};
  Configurable<bool> enableSecondarySmearing{"enableSecondarySmearing", false, "Enable smearing of weak decay daughters"};
  Configurable<bool> enableNucleiSmearing{"enableNucleiSmearing", false, "Enable smearing of nuclei"};
  Configurable<bool> enablePrimaryVertexing{"enablePrimaryVertexing", true, "Enable primary vertexing"};
  Configurable<bool> interpolateLutEfficiencyVsNch{"interpolateLutEfficiencyVsNch", true, "interpolate LUT efficiency as f(Nch)"};

  Configurable<bool> populateTracksDCA{"populateTracksDCA", true, "populate TracksDCA table"};
  Configurable<bool> populateTracksDCACov{"populateTracksDCACov", false, "populate TracksDCACov table"};
  Configurable<bool> populateTracksExtra{"populateTracksExtra", false, "populate TrackExtra table (legacy)"};
  Configurable<bool> populateTrackSelection{"populateTrackSelection", false, "populate TrackSelection table (legacy)"};

  Configurable<bool> processUnreconstructedTracks{"processUnreconstructedTracks", false, "process (smear) unreco-ed tracks"};
  Configurable<bool> doExtraQA{"doExtraQA", false, "do extra 2D QA plots"};
  Configurable<bool> extraQAwithoutDecayDaughters{"extraQAwithoutDecayDaughters", false, "remove decay daughters from qa plots (yes/no)"};

  struct : ConfigurableGroup {
    std::string prefix = "lookUpTables"; // JSON group name
    Configurable<std::vector<std::string>> lutEl{"lutEl", std::vector<std::string>{"lutCovm.el.dat"}, "LUT for electrons (if emtpy no LUT is taken)"};
    Configurable<std::vector<std::string>> lutMu{"lutMu", std::vector<std::string>{"lutCovm.mu.dat"}, "LUT for muons (if emtpy no LUT is taken)"};
    Configurable<std::vector<std::string>> lutPi{"lutPi", std::vector<std::string>{"lutCovm.pi.dat"}, "LUT for pions (if emtpy no LUT is taken)"};
    Configurable<std::vector<std::string>> lutKa{"lutKa", std::vector<std::string>{"lutCovm.ka.dat"}, "LUT for kaons (if emtpy no LUT is taken)"};
    Configurable<std::vector<std::string>> lutPr{"lutPr", std::vector<std::string>{"lutCovm.pr.dat"}, "LUT for protons (if emtpy no LUT is taken)"};
    Configurable<std::vector<std::string>> lutDe{"lutDe", std::vector<std::string>{""}, "LUT for deuterons (if emtpy no LUT is taken)"};
    Configurable<std::vector<std::string>> lutTr{"lutTr", std::vector<std::string>{""}, "LUT for tritons (if emtpy no LUT is taken)"};
    Configurable<std::vector<std::string>> lutHe3{"lutHe3", std::vector<std::string>{""}, "LUT for Helium-3 (if emtpy no LUT is taken)"};
    Configurable<std::vector<std::string>> lutAl{"lutAl", std::vector<std::string>{""}, "LUT for Alphas (if emtpy no LUT is taken)"};
  } lookUpTables;

  struct : ConfigurableGroup {
    ConfigurableAxis axisMomentum{"axisMomentum", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "#it{p} (GeV/#it{c})"};
    ConfigurableAxis axisNVertices{"axisNVertices", {20, -0.5, 19.5}, "N_{vertices}"};
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
    Configurable<int> minSiliconHitsIfTPCUsed{"minSiliconHitsIfTPCUsed", 2, "minimum number of silicon hits to accept track in case TPC info is present"};
    Configurable<int> minTPCClusters{"minTPCClusters", 70, "minimum number of TPC hits necessary to consider minSiliconHitsIfTPCUsed"};
    Configurable<std::vector<std::string>> alice3geo{"alice3geo", std::vector<std::string>{"2"}, "0: ALICE 3 v1, 1: ALICE 3 v4, 2: ALICE 3 Sep 2025, or path to ccdb with a3 geo (ccdb:Users/u/user/)"};
    Configurable<bool> applyZacceptance{"applyZacceptance", false, "apply z limits to detector layers or not"};
    Configurable<bool> applyMSCorrection{"applyMSCorrection", true, "apply ms corrections for secondaries or not"};
    Configurable<bool> applyElossCorrection{"applyElossCorrection", true, "apply eloss corrections for secondaries or not"};
    Configurable<bool> applyEffCorrection{"applyEffCorrection", true, "apply efficiency correction or not"};
    Configurable<int> scaleVD{"scaleVD", 1, "scale x0 and xrho in VD layers"};
    Configurable<std::vector<float>> pixelRes{"pixelRes", {0.00025, 0.00025, 0.001, 0.001}, "RPhiIT, ZIT, RPhiOT, ZOT"};
  } fastTrackerSettings; // allows for gap between peak and bg in case someone wants to

  struct : ConfigurableGroup {
    std::string prefix = "fastPrimaryTrackerSettings";
    Configurable<bool> fastTrackPrimaries{"fastTrackPrimaries", false, "Use fasttracker for primary tracks. Enable with care"};
    Configurable<int> minSiliconHits{"minSiliconHits", 4, "minimum number of silicon hits to accept track"};
    Configurable<std::string> alice3geo{"alice3geo", "2", "0: ALICE 3 v1, 1: ALICE 3 v4, 2: ALICE 3 Sep 2025, or path to ccdb with a3 geo"};
    Configurable<bool> applyZacceptance{"applyZacceptance", false, "apply z limits to detector layers or not"};
    Configurable<bool> applyMSCorrection{"applyMSCorrection", true, "apply ms corrections for secondaries or not"};
    Configurable<bool> applyElossCorrection{"applyElossCorrection", true, "apply eloss corrections for secondaries or not"};
    Configurable<bool> applyEffCorrection{"applyEffCorrection", true, "apply efficiency correction or not"};
  } fastPrimaryTrackerSettings;

  struct : ConfigurableGroup {
    std::string prefix = "cascadeDecaySettings"; // Cascade decay settings
    Configurable<bool> decayXi{"decayXi", false, "Manually decay Xi and fill tables with daughters"};
    Configurable<bool> findXi{"findXi", false, "if decayXi on, find Xi and fill Tracks table also with Xi"};
    Configurable<bool> trackXi{"trackXi", false, "if findXi on, attempt to track Xi"};
    Configurable<bool> doXiQA{"doXiQA", false, "QA plots for when treating Xi"};
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
  // o2::fastsim::FastTracker fastTracker;
  std::vector<std::unique_ptr<o2::fastsim::FastTracker>> fastTracker;
  o2::fastsim::FastTracker fastPrimaryTracker;

  // V0 names for filling histograms
  static constexpr std::string_view kV0names[] = {"K0", "Lambda", "AntiLambda"};

  // Class to hold the track information for the O2 vertexing
  class TrackAlice3 : public o2::track::TrackParCov
  {
    using TimeEst = o2::dataformats::TimeStampWithError<float, float>;

   public:
    TrackAlice3() = default;
    ~TrackAlice3() = default;
    TrackAlice3(const TrackAlice3& src) = default;
    TrackAlice3(const o2::track::TrackParCov& src, const int64_t label,
                const float t = 0,
                const float te = 1,
                bool decayDauInput = false,
                bool weakDecayDauInput = false,
                int isUsedInCascadingInput = 0,
                int nSiliconHitsInput = 0,
                int nTPCHitsInput = 0) : o2::track::TrackParCov(src),
                                         mcLabel{label},
                                         timeEst{t, te},
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
    int positiveId; // track index in the Tracks table
    int negativeId; // track index in the Tracks table

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

  // necessary for particle charges
  Service<o2::framework::O2DatabasePDG> pdgDB;

  // for handling basic QA histograms if requested
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  std::map<std::string, HistPtr> histPointers;

  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;

  // Track smearer
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

  static constexpr int kMaxLUTConfigs = 20;
  void init(o2::framework::InitContext&)
  {
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setTimestamp(-1);

    if (enablePrimarySmearing) {
      auto loadLUT = [&](int icfg, int pdg, const std::vector<std::string>& tables) {
        const bool foundNewCfg = static_cast<size_t>(icfg) < tables.size();
        const std::string& lutFile = foundNewCfg ? tables[icfg] : tables.front();
        bool success = mSmearer[icfg]->loadTable(pdg, lutFile.c_str());
        if (!success && !lutFile.empty()) {
          LOG(fatal) << "Having issue with loading the LUT " << pdg << " " << lutFile;
        }

        return foundNewCfg;
      };

      for (int icfg = 0; icfg < kMaxLUTConfigs; ++icfg) {
        mSmearer.emplace_back(std::make_unique<o2::delphes::DelphesO2TrackSmearer>());
        mSmearer[icfg]->setCcdbManager(ccdb.operator->());

        // check if more configs were provided, fall back to first entry
        bool newLUTLoaded = false;
        newLUTLoaded |= loadLUT(icfg, kElectron, lookUpTables.lutEl.value);
        newLUTLoaded |= loadLUT(icfg, kMuonMinus, lookUpTables.lutMu.value);
        newLUTLoaded |= loadLUT(icfg, kPiPlus, lookUpTables.lutPi.value);
        newLUTLoaded |= loadLUT(icfg, kKPlus, lookUpTables.lutKa.value);
        newLUTLoaded |= loadLUT(icfg, kProton, lookUpTables.lutPr.value);
        newLUTLoaded |= loadLUT(icfg, o2::constants::physics::kDeuteron, lookUpTables.lutDe.value);
        newLUTLoaded |= loadLUT(icfg, o2::constants::physics::kTriton, lookUpTables.lutTr.value);
        newLUTLoaded |= loadLUT(icfg, o2::constants::physics::kHelium3, lookUpTables.lutHe3.value);
        newLUTLoaded |= loadLUT(icfg, o2::constants::physics::kAlpha, lookUpTables.lutAl.value);

        if (!newLUTLoaded) {
          mSmearer.pop_back();
          break;
        }

        // interpolate efficiencies if requested to do so
        mSmearer[icfg]->interpolateEfficiency(static_cast<bool>(interpolateLutEfficiencyVsNch));

        // smear un-reco'ed tracks if asked to do so
        mSmearer[icfg]->skipUnreconstructed(static_cast<bool>(!processUnreconstructedTracks));

        std::string histPath = "Configuration_" + std::to_string(icfg) + "/";
        histPointers.insert({histPath + "hPtGenerated", histos.add((histPath + "hPtGenerated").c_str(), "hPtGenerated", {kTH1D, {{axes.axisMomentum}}})});
        histPointers.insert({histPath + "hPhiGenerated", histos.add((histPath + "hPhiGenerated").c_str(), "hPhiGenerated", {kTH1D, {{100, 0.0f, 2 * M_PI, "#phi (rad)"}}})});

        histPointers.insert({histPath + "hPtGeneratedEl", histos.add((histPath + "hPtGeneratedEl").c_str(), "hPtGeneratedEl", {kTH1D, {{axes.axisMomentum}}})});
        histPointers.insert({histPath + "hPtGeneratedPi", histos.add((histPath + "hPtGeneratedPi").c_str(), "hPtGeneratedPi", {kTH1D, {{axes.axisMomentum}}})});
        histPointers.insert({histPath + "hPtGeneratedKa", histos.add((histPath + "hPtGeneratedKa").c_str(), "hPtGeneratedKa", {kTH1D, {{axes.axisMomentum}}})});
        histPointers.insert({histPath + "hPtGeneratedPr", histos.add((histPath + "hPtGeneratedPr").c_str(), "hPtGeneratedPr", {kTH1D, {{axes.axisMomentum}}})});
        histPointers.insert({histPath + "hPtReconstructed", histos.add((histPath + "hPtReconstructed").c_str(), "hPtReconstructed", {kTH1D, {{axes.axisMomentum}}})});
        histPointers.insert({histPath + "hPtReconstructedEl", histos.add((histPath + "hPtReconstructedEl").c_str(), "hPtReconstructedEl", {kTH1D, {{axes.axisMomentum}}})});
        histPointers.insert({histPath + "hPtReconstructedPi", histos.add((histPath + "hPtReconstructedPi").c_str(), "hPtReconstructedPi", {kTH1D, {{axes.axisMomentum}}})});
        histPointers.insert({histPath + "hPtReconstructedKa", histos.add((histPath + "hPtReconstructedKa").c_str(), "hPtReconstructedKa", {kTH1D, {{axes.axisMomentum}}})});
        histPointers.insert({histPath + "hPtReconstructedPr", histos.add((histPath + "hPtReconstructedPr").c_str(), "hPtReconstructedPr", {kTH1D, {{axes.axisMomentum}}})});

        // Collision QA
        histPointers.insert({histPath + "hPVz", histos.add((histPath + "hPVz").c_str(), "hPVz", {kTH1D, {{axes.axisVertexZ}}})});
        histPointers.insert({histPath + "hLUTMultiplicity", histos.add((histPath + "hLUTMultiplicity").c_str(), "hLUTMultiplicity", {kTH1D, {{axes.axisMultiplicity}}})});
        histPointers.insert({histPath + "hSimMultiplicity", histos.add((histPath + "hSimMultiplicity").c_str(), "hSimMultiplicity", {kTH1D, {{axes.axisMultiplicity}}})});
        histPointers.insert({histPath + "hRecoMultiplicity", histos.add((histPath + "hRecoMultiplicity").c_str(), "hRecoMultiplicity", {kTH1D, {{axes.axisMultiplicity}}})});

        if (enableSecondarySmearing) {
          fastTracker.emplace_back(std::make_unique<o2::fastsim::FastTracker>());
          fastTracker[icfg]->SetMagneticField(magneticField);
          fastTracker[icfg]->SetApplyZacceptance(fastTrackerSettings.applyZacceptance);
          fastTracker[icfg]->SetApplyMSCorrection(fastTrackerSettings.applyMSCorrection);
          fastTracker[icfg]->SetApplyElossCorrection(fastTrackerSettings.applyElossCorrection);

          if (fastTrackerSettings.alice3geo.value[icfg] == "0") {
            fastTracker[icfg]->AddSiliconALICE3v2(fastTrackerSettings.pixelRes);
          } else if (fastTrackerSettings.alice3geo.value[icfg] == "1") {
            fastTracker[icfg]->AddSiliconALICE3v4(fastTrackerSettings.pixelRes);
            fastTracker[icfg]->AddTPC(0.1, 0.1);
          } else if (fastTrackerSettings.alice3geo.value[icfg] == "2") {
            fastTracker[icfg]->AddSiliconALICE3(fastTrackerSettings.scaleVD, fastTrackerSettings.pixelRes);
          } else {
            fastTracker[icfg]->AddGenericDetector(fastTrackerSettings.alice3geo.value[icfg], ccdb.operator->());
          }

          // print fastTracker settings
          fastTracker[icfg]->Print();
          histPointers.insert({histPath + "hMassXi", histos.add((histPath + "hMassXi").c_str(), "hMassXi", {kTH1D, {{axes.axisXiMass}}})});
        }

        if (doExtraQA) {
          histPointers.insert({histPath + "h2dPtRes", histos.add((histPath + "h2dPtRes").c_str(), "h2dPtRes", {kTH2D, {{axes.axisMomentum, axes.axisPtRes}}})});
          histPointers.insert({histPath + "h2dDCAxy", histos.add((histPath + "h2dDCAxy").c_str(), "h2dDCAxy", {kTH2D, {{axes.axisMomentum, axes.axisDCA}}})});
          histPointers.insert({histPath + "h2dDCAz", histos.add((histPath + "h2dDCAz").c_str(), "h2dDCAz", {kTH2D, {{axes.axisMomentum, axes.axisDCA}}})});
        }

      } // end config loop
    }

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
      histos.add("h2dVerticesVsContributors", "h2dVerticesVsContributors", kTH2F, {axes.axisMultiplicity, axes.axisNVertices});
      histos.add("hRecoVsSimMultiplicity", "hRecoVsSimMultiplicity", kTH2F, {axes.axisMultiplicity, axes.axisMultiplicity});

      histos.add("hSimTrackX", "hSimTrackX", kTH1F, {axes.axisX});
      histos.add("hRecoTrackX", "hRecoTrackX", kTH1F, {axes.axisX});
      histos.add("hTrackXatDCA", "hTrackXatDCA", kTH1F, {axes.axisX});
    }

    if (cascadeDecaySettings.doXiQA) {
      histos.add("hXiBuilding", "hXiBuilding", kTH1F, {{10, -0.5f, 9.5f}});

      histos.add("hGenXi", "hGenXi", kTH2F, {axes.axisDecayRadius, axes.axisMomentum});
      histos.add("hRecoXi", "hRecoXi", kTH2F, {axes.axisDecayRadius, axes.axisMomentum});

      histos.add("hGenPiFromXi", "hGenPiFromXi", kTH2F, {axes.axisDecayRadius, axes.axisMomentum});
      histos.add("hGenPiFromLa", "hGenPiFromLa", kTH2F, {axes.axisDecayRadius, axes.axisMomentum});
      histos.add("hGenPrFromLa", "hGenPrFromLa", kTH2F, {axes.axisDecayRadius, axes.axisMomentum});
      histos.add("hRecoPiFromXi", "hRecoPiFromXi", kTH2F, {axes.axisDecayRadius, axes.axisMomentum});
      histos.add("hRecoPiFromLa", "hRecoPiFromLa", kTH2F, {axes.axisDecayRadius, axes.axisMomentum});
      histos.add("hRecoPrFromLa", "hRecoPrFromLa", kTH2F, {axes.axisDecayRadius, axes.axisMomentum});

      // basic mass histograms to see if we're in business
      histos.add("hMassLambda", "hMassLambda", kTH1F, {axes.axisLambdaMass});
      histos.add("hMassXi", "hMassXi", kTH1F, {axes.axisXiMass});

      // OTF strangeness tracking QA
      histos.add("hFoundVsFindable", "hFoundVsFindable", kTH2F, {{10, -0.5f, 9.5f}, {10, -0.5f, 9.5f}});

      histos.add("h2dDCAxyCascade", "h2dDCAxyCascade", kTH2F, {axes.axisMomentum, axes.axisDCA});
      histos.add("h2dDCAxyCascadeBachelor", "h2dDCAxyCascadeBachelor", kTH2F, {axes.axisMomentum, axes.axisDCA});
      histos.add("h2dDCAxyCascadeNegative", "h2dDCAxyCascadeNegative", kTH2F, {axes.axisMomentum, axes.axisDCA});
      histos.add("h2dDCAxyCascadePositive", "h2dDCAxyCascadePositive", kTH2F, {axes.axisMomentum, axes.axisDCA});

      histos.add("h2dDCAzCascade", "h2dDCAzCascade", kTH2F, {axes.axisMomentum, axes.axisDCA});
      histos.add("h2dDCAzCascadeBachelor", "h2dDCAzCascadeBachelor", kTH2F, {axes.axisMomentum, axes.axisDCA});
      histos.add("h2dDCAzCascadeNegative", "h2dDCAzCascadeNegative", kTH2F, {axes.axisMomentum, axes.axisDCA});
      histos.add("h2dDCAzCascadePositive", "h2dDCAzCascadePositive", kTH2F, {axes.axisMomentum, axes.axisDCA});

      histos.add("h2dDeltaPtVsPt", "h2dDeltaPtVsPt", kTH2F, {axes.axisMomentum, axes.axisDeltaPt});
      histos.add("h2dDeltaEtaVsPt", "h2dDeltaEtaVsPt", kTH2F, {axes.axisMomentum, axes.axisDeltaEta});

      histos.add("hFastTrackerHits", "hFastTrackerHits", kTH2F, {axes.axisZ, axes.axisRadius});
      auto hFastTrackerQA = histos.add<TH1>("hFastTrackerQA", "hFastTrackerQA", kTH1D, {{8, -0.5f, 7.5f}});
      hFastTrackerQA->GetXaxis()->SetBinLabel(1, "Negative eigenvalue");
      hFastTrackerQA->GetXaxis()->SetBinLabel(2, "Failed sanity check");
      hFastTrackerQA->GetXaxis()->SetBinLabel(3, "intercept original radius");
      hFastTrackerQA->GetXaxis()->SetBinLabel(4, "propagate to original radius");
      hFastTrackerQA->GetXaxis()->SetBinLabel(5, "problematic layer");
      hFastTrackerQA->GetXaxis()->SetBinLabel(6, "multiple scattering");
      hFastTrackerQA->GetXaxis()->SetBinLabel(7, "energy loss");
      hFastTrackerQA->GetXaxis()->SetBinLabel(8, "efficiency");
    }
    if (v0DecaySettings.doV0QA) {
      histos.add("V0Building/hV0Building", "hV0Building", kTH1F, {{10, -0.5f, 9.5f}});

      histos.add("V0Building/K0/hGen", "hGen", kTH2F, {axes.axisDecayRadius, axes.axisMomentum});
      histos.add("V0Building/K0/hReco", "hReco", kTH2F, {axes.axisDecayRadius, axes.axisMomentum});
      histos.add("V0Building/K0/hGenNegDaughterFromV0", "hGenNegDaughterFromV0", kTH2F, {axes.axisDecayRadius, axes.axisMomentum});
      histos.add("V0Building/K0/hGenPosDaughterFromV0", "hGenPosDaughterFromV0", kTH2F, {axes.axisDecayRadius, axes.axisMomentum});
      histos.add("V0Building/K0/hRecoNegDaughterFromV0", "hRecoNegDaughterFromV0", kTH2F, {axes.axisDecayRadius, axes.axisMomentum});
      histos.add("V0Building/K0/hRecoPosDaughterFromV0", "hRecoPosDaughterFromV0", kTH2F, {axes.axisDecayRadius, axes.axisMomentum});

      // histos.add("V0Building/K0/h2dDCAxyV0Negative", "h2dDCAxyV0Negative", kTH2F, {axes.axisMomentum, axes.axisDCA});
      // histos.add("V0Building/K0/h2dDCAxyV0Positive", "h2dDCAxyV0Positive", kTH2F, {axes.axisMomentum, axes.axisDCA});

      // histos.add("V0Building/K0/h2dDCAzV0Negative", "h2dDCAzV0Negative", kTH2F, {axes.axisMomentum, axes.axisDCA});
      // histos.add("V0Building/K0/h2dDCAzV0Positive", "h2dDCAzV0Positive", kTH2F, {axes.axisMomentum, axes.axisDCA});

      histos.addClone("V0Building/K0/", "V0Building/Lambda/");
      histos.addClone("V0Building/K0/", "V0Building/AntiLambda/");

      histos.add("V0Building/K0/hMass", "hMass", kTH2F, {axes.axisK0Mass, axes.axisMomentum});
      histos.add("V0Building/K0/hFinalMass", "hMass", kTH1F, {axes.axisK0Mass});
      histos.add("V0Building/Lambda/hMass", "hMass", kTH2F, {axes.axisLambdaMass, axes.axisMomentum});
      histos.add("V0Building/AntiLambda/hMass", "hMass", kTH2F, {axes.axisLambdaMass, axes.axisMomentum});

      histos.add("V0Building/hFastTrackerHits", "hFastTrackerHits", kTH2F, {axes.axisZ, axes.axisRadius});
      auto hFastTrackerQA = histos.add<TH1>("V0Building/hFastTrackerQA", "hFastTrackerQA", kTH1D, {{8, -0.5f, 7.5f}});
      hFastTrackerQA->GetXaxis()->SetBinLabel(1, "Negative eigenvalue");
      hFastTrackerQA->GetXaxis()->SetBinLabel(2, "Failed sanity check");
      hFastTrackerQA->GetXaxis()->SetBinLabel(3, "intercept original radius");
      hFastTrackerQA->GetXaxis()->SetBinLabel(4, "propagate to original radius");
      hFastTrackerQA->GetXaxis()->SetBinLabel(5, "problematic layer");
      hFastTrackerQA->GetXaxis()->SetBinLabel(6, "multiple scattering");
      hFastTrackerQA->GetXaxis()->SetBinLabel(7, "energy loss");
      hFastTrackerQA->GetXaxis()->SetBinLabel(8, "efficiency");
    }

    LOGF(info, "Initializing magnetic field to value: %.3f kG", static_cast<float>(magneticField));
    o2::parameters::GRPMagField grpmag;
    grpmag.setFieldUniformity(true);
    grpmag.setL3Current(30000.f * (magneticField / 5.0f));
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
    vertexer.setBunchFilling(irSampler.getBunchFilling());
    vertexer.init();

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

    // Configure FastTracker for primaries
    if (fastPrimaryTrackerSettings.fastTrackPrimaries) {
      fastPrimaryTracker.SetMagneticField(magneticField);
      fastPrimaryTracker.SetApplyZacceptance(fastPrimaryTrackerSettings.applyZacceptance);
      fastPrimaryTracker.SetApplyMSCorrection(fastPrimaryTrackerSettings.applyMSCorrection);
      fastPrimaryTracker.SetApplyElossCorrection(fastPrimaryTrackerSettings.applyElossCorrection);

      if (fastPrimaryTrackerSettings.alice3geo.value == "0") {
        fastPrimaryTracker.AddSiliconALICE3v2({0.00025, 0.00025, 0.001, 0.001});
      } else if (fastPrimaryTrackerSettings.alice3geo.value == "1") {
        fastPrimaryTracker.AddSiliconALICE3v4({0.00025, 0.00025, 0.001, 0.001});
        fastPrimaryTracker.AddTPC(0.1, 0.1);
      } else if (fastPrimaryTrackerSettings.alice3geo.value == "2") {
        fastPrimaryTracker.AddSiliconALICE3(1., {0.00025, 0.00025, 0.001, 0.001});
      } else {
        fastPrimaryTracker.AddGenericDetector(fastPrimaryTrackerSettings.alice3geo, ccdb.operator->());
      }

      // print fastTracker settings
      fastPrimaryTracker.Print();
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
    track.getCircleParams(magneticField, circleXi, sna, csa);
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

    if (std::abs(pdgCode) == kK0Short) {
      v0Mass = o2::constants::physics::MassK0Short;
      negDauMass = o2::constants::physics::MassPionCharged;
      posDauMass = o2::constants::physics::MassPionCharged;
      ctau = 2.68;
    } else if (std::abs(pdgCode) == kLambda0) {
      v0Mass = o2::constants::physics::MassLambda;
      negDauMass = o2::constants::physics::MassPionCharged;
      posDauMass = o2::constants::physics::MassProton;
      ctau = 7.845;
    } else if (std::abs(pdgCode) == kLambda0Bar) {
      v0Mass = o2::constants::physics::MassLambda;
      negDauMass = o2::constants::physics::MassProton;
      posDauMass = o2::constants::physics::MassPionCharged;
      ctau = 7.845;
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

  float dNdEta = 0.f; // Charged particle multiplicity to use in the efficiency evaluation
  void processWithLUTs(aod::McCollision const& mcCollision, aod::McParticles const& mcParticles, int const& icfg)
  {
    int lastTrackIndex = tableStoredTracksCov.lastIndex() + 1; // bookkeep the last added track
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
    const float eventCollisionTime = ir.timeInBCNS;

    constexpr std::array<int, 5> longLivedHandledPDGs = {kElectron,
                                                         kMuonMinus,
                                                         kPiPlus,
                                                         kKPlus,
                                                         kProton};

    constexpr std::array<int, 4> nucleiPDGs = {o2::constants::physics::kDeuteron,
                                               o2::constants::physics::kTriton,
                                               o2::constants::physics::kHelium3,
                                               o2::constants::physics::kAlpha};

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

    dNdEta /= (multEtaRange * 2.0f);
    uint32_t multiplicityCounter = 0;
    getHist(TH1, histPath + "hLUTMultiplicity")->Fill(dNdEta);

    gRandom->SetSeed(seed);

    for (const auto& mcParticle : mcParticles) {
      double xiDecayRadius2D = 0;
      double laDecayRadius2D = 0;
      double v0DecayRadius2D = 0;
      std::vector<TLorentzVector> decayProducts;
      std::vector<TLorentzVector> v0DecayProducts;
      std::vector<double> xiDecayVertex, laDecayVertex, v0DecayVertex;
      if (cascadeDecaySettings.decayXi) {
        if (mcParticle.pdgCode() == kXiMinus) {
          o2::track::TrackParCov xiTrackParCov;
          o2::upgrade::convertMCParticleToO2Track(mcParticle, xiTrackParCov, pdgDB);
          decayCascade(mcParticle, xiTrackParCov, decayProducts, xiDecayVertex, laDecayVertex);
          xiDecayRadius2D = std::hypot(xiDecayVertex[0], xiDecayVertex[1]);
          laDecayRadius2D = std::hypot(laDecayVertex[0], laDecayVertex[1]);
        }
      }
      const auto pdg = std::abs(mcParticle.pdgCode());
      bool isV0 = std::find(v0PDGs.begin(), v0PDGs.end(), pdg) != v0PDGs.end();
      bool isK0 = (std::abs(pdg) == kK0Short);
      bool isLambda = (std::abs(pdg) == kLambda0);
      bool isAntiLambda = (std::abs(pdg) == kLambda0Bar);

      if (v0DecaySettings.decayV0 && isV0) {
        decayV0Particle(mcParticle, v0DecayProducts, v0DecayVertex, pdg);
        v0DecayRadius2D = std::hypot(v0DecayVertex[0], v0DecayVertex[1]);
      }

      if (!mcParticle.isPhysicalPrimary()) {
        continue;
      }

      const bool longLivedToBeHandled = std::find(longLivedHandledPDGs.begin(), longLivedHandledPDGs.end(), pdg) != longLivedHandledPDGs.end();
      const bool nucleiToBeHandled = std::find(nucleiPDGs.begin(), nucleiPDGs.end(), pdg) != nucleiPDGs.end();
      const bool pdgsToBeHandled = longLivedToBeHandled || (enableNucleiSmearing && nucleiToBeHandled) || (cascadeDecaySettings.decayXi && mcParticle.pdgCode() == kXiMinus) || (v0DecaySettings.decayV0 && isV0);
      if (!pdgsToBeHandled) {
        continue;
      }

      if (std::fabs(mcParticle.eta()) > maxEta) {
        continue;
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

      if (cascadeDecaySettings.doXiQA && mcParticle.pdgCode() == kXiMinus) {
        histos.fill(HIST("hGenXi"), xiDecayRadius2D, mcParticle.pt());
        histos.fill(HIST("hGenPiFromXi"), xiDecayRadius2D, decayProducts[0].Pt());
        histos.fill(HIST("hGenPiFromLa"), laDecayRadius2D, decayProducts[1].Pt());
        histos.fill(HIST("hGenPrFromLa"), laDecayRadius2D, decayProducts[2].Pt());
      }
      if (v0DecaySettings.doV0QA && isV0) {
        static_for<0, 2>([&](auto i) {
          constexpr int Index = i.value;
          histos.fill(HIST("V0Building/") + HIST(kV0names[Index]) + HIST("/hGen"), v0DecayRadius2D, mcParticle.pt());
          histos.fill(HIST("V0Building/") + HIST(kV0names[Index]) + HIST("/hGenNegDaughterFromV0"), v0DecayRadius2D, v0DecayProducts[0].Pt());
          histos.fill(HIST("V0Building/") + HIST(kV0names[Index]) + HIST("/hGenPosDaughterFromV0"), v0DecayRadius2D, v0DecayProducts[1].Pt());
        });
      }
      if (mcParticle.pt() < minPt) {
        continue;
      }

      o2::track::TrackParCov trackParCov;
      o2::upgrade::convertMCParticleToO2Track(mcParticle, trackParCov, pdgDB);

      bool isDecayDaughter = false;
      if (mcParticle.getProcess() == TMCProcess::kPDecay)
        isDecayDaughter = true;

      multiplicityCounter++;
      const float t = (eventCollisionTime + gRandom->Gaus(0., 100.)) * 1e-3;
      static constexpr int kCascProngs = 3;
      std::vector<o2::track::TrackParCov> xiDaughterTrackParCovsPerfect(3);
      std::vector<o2::track::TrackParCov> xiDaughterTrackParCovsTracked(3);
      std::vector<bool> isReco(kCascProngs);
      std::vector<int> nHits(kCascProngs);        // total
      std::vector<int> nSiliconHits(kCascProngs); // silicon type
      std::vector<int> nTPCHits(kCascProngs);     // TPC type
      if (cascadeDecaySettings.decayXi && mcParticle.pdgCode() == kXiMinus) {
        if (cascadeDecaySettings.doXiQA) {
          histos.fill(HIST("hXiBuilding"), 0.0f);
        }

        o2::upgrade::convertTLorentzVectorToO2Track(kPiMinus, decayProducts[0], xiDecayVertex, xiDaughterTrackParCovsPerfect[0], pdgDB);
        o2::upgrade::convertTLorentzVectorToO2Track(kPiMinus, decayProducts[1], laDecayVertex, xiDaughterTrackParCovsPerfect[1], pdgDB);
        o2::upgrade::convertTLorentzVectorToO2Track(kProton, decayProducts[2], laDecayVertex, xiDaughterTrackParCovsPerfect[2], pdgDB);

        for (int i = 0; i < kCascProngs; i++) {
          isReco[i] = false;
          nHits[i] = 0;
          nSiliconHits[i] = 0;
          nTPCHits[i] = 0;
          if (enableSecondarySmearing) {
            nHits[i] = fastTracker[icfg]->FastTrack(xiDaughterTrackParCovsPerfect[i], xiDaughterTrackParCovsTracked[i], dNdEta);
            nSiliconHits[i] = fastTracker[icfg]->GetNSiliconPoints();
            nTPCHits[i] = fastTracker[icfg]->GetNGasPoints();

            if (nHits[i] < 0) { // QA
              histos.fill(HIST("hFastTrackerQA"), o2::math_utils::abs(nHits[i]));
            }

            if (nSiliconHits[i] >= fastTrackerSettings.minSiliconHits || (nSiliconHits[i] >= fastTrackerSettings.minSiliconHitsIfTPCUsed && nTPCHits[i] >= fastTrackerSettings.minTPCClusters)) {
              isReco[i] = true;
            } else {
              continue; // extra sure
            }
            for (uint32_t ih = 0; ih < fastTracker[icfg]->GetNHits(); ih++) {
              histos.fill(HIST("hFastTrackerHits"), fastTracker[icfg]->GetHitZ(ih), std::hypot(fastTracker[icfg]->GetHitX(ih), fastTracker[icfg]->GetHitY(ih)));
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
            tracksAlice3.push_back(TrackAlice3{xiDaughterTrackParCovsTracked[i], mcParticle.globalIndex(), t, 100.f * 1e-3, true, true, i + 2, nSiliconHits[i], nTPCHits[i]});
          } else {
            ghostTracksAlice3.push_back(TrackAlice3{xiDaughterTrackParCovsTracked[i], mcParticle.globalIndex(), t, 100.f * 1e-3, true, true, i + 2});
          }
        }

        if (cascadeDecaySettings.doXiQA && mcParticle.pdgCode() == kXiMinus) {
          if (isReco[0] && isReco[1] && isReco[2]) {
            histos.fill(HIST("hXiBuilding"), 2.0f);
            histos.fill(HIST("hRecoXi"), xiDecayRadius2D, mcParticle.pt());
          }
          if (isReco[0])
            histos.fill(HIST("hRecoPiFromXi"), xiDecayRadius2D, decayProducts[0].Pt());
          if (isReco[1])
            histos.fill(HIST("hRecoPiFromLa"), laDecayRadius2D, decayProducts[1].Pt());
          if (isReco[2])
            histos.fill(HIST("hRecoPrFromLa"), laDecayRadius2D, decayProducts[2].Pt());
        }

        // +-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+
        // combine particles into actual Xi candidate
        // cascade building starts here
        if (cascadeDecaySettings.findXi && mcParticle.pdgCode() == kXiMinus && isReco[0] && isReco[1] && isReco[2]) {
          if (cascadeDecaySettings.doXiQA)
            histos.fill(HIST("hXiBuilding"), 3.0f);
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
            if (cascadeDecaySettings.doXiQA)
              histos.fill(HIST("hXiBuilding"), 4.0f);
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
            v0Track.setPID(o2::track::PID::Lambda);

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
              if (cascadeDecaySettings.doXiQA)
                histos.fill(HIST("hXiBuilding"), 5.0f);
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
              cascadeTrack.setAbsCharge(-1);                // may require more adjustments
              cascadeTrack.setPID(o2::track::PID::XiMinus); // FIXME: not OK for omegas

              thisCascade.cascradiusMC = xiDecayRadius2D;
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
                  trackParCov.getXatLabR(layer.getRadius(), targetX, magneticField);
                  if (targetX > 999)
                    continue; // failed to find intercept

                  if (!trackParCov.propagateTo(targetX, magneticField)) {
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

                  if (!(cascadeTrack.propagateTo(xyz1[0], magneticField)))
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

              // add cascade track
              thisCascade.cascadeTrackId = lastTrackIndex + tracksAlice3.size(); // this is the next index to be filled -> should be it
              tracksAlice3.push_back(TrackAlice3{cascadeTrack, mcParticle.globalIndex(), t, 100.f * 1e-3, false, false, 1, thisCascade.foundClusters});

              if (cascadeDecaySettings.doXiQA) {
                histos.fill(HIST("hXiBuilding"), 6.0f);
                histos.fill(HIST("h2dDeltaPtVsPt"), trackParCov.getPt(), cascadeTrack.getPt() - trackParCov.getPt());
                histos.fill(HIST("h2dDeltaEtaVsPt"), trackParCov.getPt(), cascadeTrack.getEta() - trackParCov.getEta());

                histos.fill(HIST("hMassLambda"), thisCascade.mLambda);
                histos.fill(HIST("hMassXi"), thisCascade.mXi);
                histos.fill(HIST("hFoundVsFindable"), thisCascade.findableClusters, thisCascade.foundClusters);
                getHist(TH1, histPath + "hMassXi")->Fill(thisCascade.mXi);
              }

              // add this cascade to vector (will fill cursor later with collision ID)
              cascadesAlice3.push_back(thisCascade);
            }
          }
        } // end cascade building
        // +-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+
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
          histos.fill(HIST("V0Building/hV0Building"), 0.0f);
        }
        if (isK0) {
          o2::upgrade::convertTLorentzVectorToO2Track(kPiMinus, v0DecayProducts[0], v0DecayVertex, v0DaughterTrackParCovsPerfect[0], pdgDB);
          o2::upgrade::convertTLorentzVectorToO2Track(kPiPlus, v0DecayProducts[1], v0DecayVertex, v0DaughterTrackParCovsPerfect[1], pdgDB);
        } else if (isLambda) {
          o2::upgrade::convertTLorentzVectorToO2Track(kPiMinus, v0DecayProducts[0], v0DecayVertex, v0DaughterTrackParCovsPerfect[0], pdgDB);
          o2::upgrade::convertTLorentzVectorToO2Track(kProton, v0DecayProducts[1], v0DecayVertex, v0DaughterTrackParCovsPerfect[1], pdgDB);
        } else if (isAntiLambda) {
          o2::upgrade::convertTLorentzVectorToO2Track(kProtonBar, v0DecayProducts[0], v0DecayVertex, v0DaughterTrackParCovsPerfect[0], pdgDB);
          o2::upgrade::convertTLorentzVectorToO2Track(kPiPlus, v0DecayProducts[1], v0DecayVertex, v0DaughterTrackParCovsPerfect[1], pdgDB);
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

            if (nV0Hits[i] < 0) { // QA
              histos.fill(HIST("V0Building/hFastTrackerQA"), o2::math_utils::abs(nV0Hits[i]));
            }

            if (nV0SiliconHits[i] >= fastTrackerSettings.minSiliconHits || (nV0SiliconHits[i] >= fastTrackerSettings.minSiliconHitsIfTPCUsed && nV0TPCHits[i] >= fastTrackerSettings.minTPCClusters)) {
              isReco[i] = true;
            } else {
              continue; // extra sure
            }
            for (uint32_t ih = 0; ih < fastTracker[icfg]->GetNHits(); ih++) {
              histos.fill(HIST("V0Building/hFastTrackerHits"), fastTracker[icfg]->GetHitZ(ih), std::hypot(fastTracker[icfg]->GetHitX(ih), fastTracker[icfg]->GetHitY(ih)));
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
            tracksAlice3.push_back(TrackAlice3{v0DaughterTrackParCovsTracked[i], mcParticle.globalIndex(), t, 100.f * 1e-3, true, true, i + 2, nSiliconHits[i], nTPCHits[i]});
          } else {
            ghostTracksAlice3.push_back(TrackAlice3{v0DaughterTrackParCovsTracked[i], mcParticle.globalIndex(), t, 100.f * 1e-3, true, true, i + 2});
          }
        }
        if (v0DecaySettings.doV0QA) {
          static_for<0, 2>([&](auto i) {
            constexpr int Index = i.value;
            if (pdg == v0PDGs[Index]) {
              if (isReco[0] && isReco[1]) {
                histos.fill(HIST("V0Building/") + HIST(kV0names[Index]) + HIST("/hReco"), v0DecayRadius2D, mcParticle.pt());
              }
              if (isReco[0])
                histos.fill(HIST("V0Building/") + HIST(kV0names[Index]) + HIST("/hRecoNegDaughterFromV0"), v0DecayRadius2D, v0DecayProducts[0].Pt());
              if (isReco[1])
                histos.fill(HIST("V0Building/") + HIST(kV0names[Index]) + HIST("/hRecoPosDaughterFromV0"), v0DecayRadius2D, v0DecayProducts[1].Pt());
            }
          });
          if (isReco[0] && isReco[1]) {
            histos.fill(HIST("V0Building/hV0Building"), 1.0f);
          }
        }
        // +-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+-~-+
        // combine particles into actual V0 candidate
        // V0 building starts here
        if (v0DecaySettings.findV0 && isReco[0] && isReco[1]) {
          if (v0DecaySettings.doV0QA)
            histos.fill(HIST("V0Building/hV0Building"), 2.0f);
          // assign indices of the particles we've used
          // they should be the last ones to be filled, in order:
          // n-1: positive Track from V0
          // n-2: negative Track from V0
          thisV0.positiveId = lastTrackIndex + tracksAlice3.size() - 1;
          thisV0.negativeId = lastTrackIndex + tracksAlice3.size() - 2;
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
            if (v0DecaySettings.doV0QA)
              histos.fill(HIST("V0Building/hV0Building"), 3.0f);
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
            // thisV0.mLambda = RecoDecay::m(std::array{std::array{posP[0], posP[1], posP[2]},
            //                                                 std::array{negP[0], negP[1], negP[2]}},
            //                                      std::array{o2::constants::physics::MassProton,
            //                                                 o2::constants::physics::MassPionCharged});
            if (isK0) {
              thisV0.mK0 = RecoDecay::m(std::array{std::array{posP[0], posP[1], posP[2]},
                                                   std::array{negP[0], negP[1], negP[2]}},
                                        std::array{o2::constants::physics::MassPionCharged,
                                                   o2::constants::physics::MassPionCharged});
            } else {
              thisV0.mK0 = -1;
            }

            if (isLambda) {
              thisV0.mLambda = RecoDecay::m(std::array{std::array{posP[0], posP[1], posP[2]},
                                                       std::array{negP[0], negP[1], negP[2]}},
                                            std::array{o2::constants::physics::MassPionCharged,
                                                       o2::constants::physics::MassProton});
            } else {
              thisV0.mLambda = -1;
            }

            if (isAntiLambda) {
              thisV0.mAntiLambda = RecoDecay::m(std::array{std::array{posP[0], posP[1], posP[2]},
                                                           std::array{negP[0], negP[1], negP[2]}},
                                                std::array{o2::constants::physics::MassProton,
                                                           o2::constants::physics::MassPionCharged});
            } else {
              thisV0.mAntiLambda = -1;
            }

            if (v0DecaySettings.doV0QA) {
              histos.fill(HIST("V0Building/hV0Building"), 4.0f);

              histos.fill(HIST("V0Building/K0/hMass"), thisV0.mK0, thisV0.pt);
              histos.fill(HIST("V0Building/Lambda/hMass"), thisV0.mLambda, thisV0.pt);
              histos.fill(HIST("V0Building/AntiLambda/hMass"), thisV0.mAntiLambda, thisV0.pt);
            }

            // add this V0 to vector (will fill cursor later with collision ID)
            v0sAlice3.push_back(thisV0);
          }
        }
      }

      if (doExtraQA) {
        histos.fill(HIST("hSimTrackX"), trackParCov.getX());
      }

      bool reconstructed = true;
      if (enablePrimarySmearing && !fastPrimaryTrackerSettings.fastTrackPrimaries) {
        reconstructed = mSmearer[icfg]->smearTrack(trackParCov, mcParticle.pdgCode(), dNdEta);
      } else if (fastPrimaryTrackerSettings.fastTrackPrimaries) {
        o2::track::TrackParCov o2Track;
        o2::upgrade::convertMCParticleToO2Track(mcParticle, o2Track, pdgDB);
        int nHits = fastPrimaryTracker.FastTrack(o2Track, trackParCov, dNdEta);
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
      getHist(TH1, histPath + "hPtReconstructed")->Fill(trackParCov.getPt());
      if (std::abs(mcParticle.pdgCode()) == kElectron)
        getHist(TH1, histPath + "hPtReconstructedEl")->Fill(trackParCov.getPt());
      if (std::abs(mcParticle.pdgCode()) == kPiPlus)
        getHist(TH1, histPath + "hPtReconstructedPi")->Fill(trackParCov.getPt());
      if (std::abs(mcParticle.pdgCode()) == kKPlus)
        getHist(TH1, histPath + "hPtReconstructedKa")->Fill(trackParCov.getPt());
      if (std::abs(mcParticle.pdgCode()) == kProton)
        getHist(TH1, histPath + "hPtReconstructedPr")->Fill(trackParCov.getPt());

      if (doExtraQA) {
        getHist(TH2, histPath + "h2dPtRes")->Fill(trackParCov.getPt(), (trackParCov.getPt() - mcParticle.pt()) / trackParCov.getPt());
        histos.fill(HIST("hRecoTrackX"), trackParCov.getX());
      }

      // populate vector with track if we reco-ed it
      if (reconstructed) {
        tracksAlice3.push_back(TrackAlice3{trackParCov, mcParticle.globalIndex(), t, 100.f * 1e-3, isDecayDaughter});
      } else {
        ghostTracksAlice3.push_back(TrackAlice3{trackParCov, mcParticle.globalIndex(), t, 100.f * 1e-3, isDecayDaughter});
      }
    }

    // *+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*
    // Calculate primary vertex with tracks from this collision
    // data preparation
    o2::vertexing::PVertex primaryVertex;
    if (enablePrimaryVertexing) {
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

      // Calculate vertices
      const int n_vertices = vertexer.process(tracksAlice3, // track array
                                              idxVec,
                                              gsl::span<o2::InteractionRecord>{bcData},
                                              vertices,
                                              vertexTrackIDs,
                                              v2tRefs,
                                              gsl::span<const o2::MCCompLabel>{lblTracks},
                                              lblVtx);

      if (n_vertices < 1) {
        return; // primary vertex not reconstructed
      }

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
                    eventCollisionTime, 0.f); // For the moment the event collision time is taken as the "GEANT" time, the computation of the event time is done a posteriori from the tracks in the OTF TOF PID task
    tableMcCollisionLabels(mcCollision.globalIndex(), 0);
    tableCollisionsAlice3(dNdEta);
    // *+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*

    // *+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*
    // populate tracks
    for (const auto& trackParCov : tracksAlice3) {
      // Fixme: collision index could be changeable
      aod::track::TrackTypeEnum trackType = aod::track::Track;

      if (populateTracksDCA) {
        float dcaXY = 1e+10, dcaZ = 1e+10;
        o2::track::TrackParCov trackParametrization(trackParCov);
        if (trackParametrization.propagateToDCA(primaryVertex, magneticField, &dcaInfo)) {
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
            histos.fill(HIST("h2dDCAxyCascade"), trackParametrization.getPt(), dcaXY * 1e+4); // in microns, please
            histos.fill(HIST("h2dDCAzCascade"), trackParametrization.getPt(), dcaZ * 1e+4);   // in microns, please
          }
          if (trackParCov.isUsedInCascading == 2) {
            histos.fill(HIST("h2dDCAxyCascadeBachelor"), trackParametrization.getPt(), dcaXY * 1e+4); // in microns, please
            histos.fill(HIST("h2dDCAzCascadeBachelor"), trackParametrization.getPt(), dcaZ * 1e+4);   // in microns, please
          }
          if (trackParCov.isUsedInCascading == 3) {
            histos.fill(HIST("h2dDCAxyCascadeNegative"), trackParametrization.getPt(), dcaXY * 1e+4); // in microns, please
            histos.fill(HIST("h2dDCAzCascadeNegative"), trackParametrization.getPt(), dcaZ * 1e+4);   // in microns, please
          }
          if (trackParCov.isUsedInCascading == 4) {
            histos.fill(HIST("h2dDCAxyCascadePositive"), trackParametrization.getPt(), dcaXY * 1e+4); // in microns, please
            histos.fill(HIST("h2dDCAzCascadePositive"), trackParametrization.getPt(), dcaZ * 1e+4);   // in microns, please
          }
        }
        tableTracksDCA(dcaXY, dcaZ);
        if (populateTracksDCACov) {
          tableTracksDCACov(dcaInfo.getSigmaY2(), dcaInfo.getSigmaZ2());
        }
      }
      tableOTFLUTConfigId(icfg);
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
    for (const auto& trackParCov : ghostTracksAlice3) {
      // Fixme: collision index could be changeable
      aod::track::TrackTypeEnum trackType = aod::track::Track;

      if (populateTracksDCA) {
        float dcaXY = 1e+10, dcaZ = 1e+10;
        o2::track::TrackParCov trackParametrization(trackParCov);
        if (trackParametrization.propagateToDCA(primaryVertex, magneticField, &dcaInfo)) {
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
    for (const auto& v0 : v0sAlice3) {
      if (v0.mK0 > 0)
        histos.fill(HIST("V0Building/K0/hFinalMass"), v0.mK0);
      tableUpgradeV0s(tableCollisions.lastIndex(), // now we know the collision index -> populate table
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
  } // end process

  void process(aod::McCollision const& mcCollision, aod::McParticles const& mcParticles)
  {
    for (size_t icfg = 0; icfg < mSmearer.size(); ++icfg) {
      processWithLUTs(mcCollision, mcParticles, static_cast<int>(icfg));
    }
  }
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
