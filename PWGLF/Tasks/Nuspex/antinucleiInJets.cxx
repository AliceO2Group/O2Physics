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
///
/// \file antinucleiInJets.cxx
///
/// \brief task for analysis of antinuclei in jets using Fastjet
/// \author Alberto Caliva (alberto.caliva@cern.ch), Chiara Pinto (chiara.pinto@cern.ch), Francesca Casillo (francesca.casillo@cern.ch)
/// \since February 13, 2025

#include "PWGJE/Core/JetBkgSubUtils.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/Core/JetUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/Zorro.h"
#include "Common/Core/ZorroSummary.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CCDB/CcdbApi.h"
#include "DataFormatsTPC/BetheBlochAleph.h"
#include "Framework/ASoA.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/DataTypes.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/Logger.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/PID.h"
#include "ReconstructionDataFormats/Track.h"

#include "TGrid.h"
#include <Math/GenVector/Boost.h>
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>
#include <TList.h>
#include <TMath.h>
#include <TPDGCode.h>
#include <TRandom3.h>
#include <TVector2.h>
#include <TVector3.h>

#include <fastjet/AreaDefinition.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include <fastjet/GhostedAreaSpec.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/Selector.hh>
#include <fastjet/tools/JetMedianBackgroundEstimator.hh>
#include <fastjet/tools/Subtractor.hh>

#include <chrono>
#include <cmath>
#include <memory>
#include <random>
#include <string>
#include <vector>

using namespace std;
using namespace o2;
using namespace o2::soa;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;
using namespace o2::constants::math;
using std::array;

// Define convenient aliases for commonly used table joins
using SelectedCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms>;
using RecCollisionsMc = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::McCollisionLabels>;
using GenCollisionsMc = aod::McCollisions;
using AntiNucleiTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TrackSelectionExtension, aod::TracksDCA, aod::pidTPCFullPr, aod::pidTPCFullDe, aod::pidTPCFullHe, aod::pidTOFFullPr, aod::pidTOFFullDe, aod::pidTOFFullHe>;
using AntiNucleiTracksMc = soa::Join<AntiNucleiTracks, aod::McTrackLabels>;

using LorentzVector = ROOT::Math::PxPyPzEVector;

// Lightweight particle container for fast kinematic access
struct ReducedParticle {
  double px;
  double py;
  double pz;
  int pdgCode;
  int mcIndex;
  bool used;

  // Pseudorapidity
  double eta() const
  {
    double p = std::sqrt(px * px + py * py + pz * pz);
    if (p == std::abs(pz)) {
      return (pz >= 0) ? 1e10 : -1e10;
    }
    return 0.5 * std::log((p + pz) / (p - pz));
  }

  // Azimuthal Angle
  double phi() const
  {
    double angle = PI + std::atan2(-py, -px);
    return angle;
  }

  // Transverse Momentum
  double pt() const
  {
    return std::sqrt(px * px + py * py);
  }
};

struct AntinucleiInJets {

  // Histogram registries for data, MC, quality control, multiplicity and correlations
  HistogramRegistry registryData{"registryData", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry registryMC{"registryMC", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry registryQC{"registryQC", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry registryMult{"registryMult", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry registryCorr{"registryCorr", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // Random generator for subsample assignment
  TRandom3 mRand;

  // Event selection criteria
  Configurable<bool> rejectITSROFBorder{"rejectITSROFBorder", true, "Reject events near the ITS ROF border"};
  Configurable<bool> rejectTFBorder{"rejectTFBorder", true, "Reject events near the TF border"};
  Configurable<bool> requireVtxITSTPC{"requireVtxITSTPC", true, "Require at least one ITS-TPC matched track"};
  Configurable<bool> rejectSameBunchPileup{"rejectSameBunchPileup", true, "Reject events with same-bunch pileup collisions"};
  Configurable<bool> requireIsGoodZvtxFT0VsPV{"requireIsGoodZvtxFT0VsPV", true, "Require consistent FT0 vs PV z-vertex"};
  Configurable<bool> requireIsVertexTOFmatched{"requireIsVertexTOFmatched", false, "Require at least one vertex track matched to TOF"};

  // Skimmed data flag and list of active triggers for processing
  Configurable<bool> cfgSkimmedProcessing{"cfgSkimmedProcessing", false, "Skimmed dataset processing"};
  Configurable<std::string> triggerList{"triggerList", "fHe", "Trigger list"};

  // Jet selection and event filtering parameters
  Configurable<double> minJetPt{"minJetPt", 10.0, "Minimum pt of the jet after bkg subtraction"};
  Configurable<double> ptLeadingMin{"ptLeadingMin", 5.0, "pt Leading Min"};
  Configurable<double> rJet{"rJet", 0.3, "Jet resolution parameter R"};
  Configurable<double> zVtx{"zVtx", 10.0, "Maximum zVertex"};
  Configurable<bool> applyAreaCut{"applyAreaCut", true, "apply area cut"};
  Configurable<double> maxNormalizedJetArea{"maxNormalizedJetArea", 1.0, "area cut"};
  Configurable<double> deltaEtaEdge{"deltaEtaEdge", 0.05, "eta gap from the edge"};
  Configurable<int> nSyst{"nSyst", 50, "number of systematic variations"};
  Configurable<int> nSubsamples{"nSubsamples", 50, "number of subsamples"};

  // Track quality, kinematic, and PID selection parameters
  Configurable<bool> requirePvContributor{"requirePvContributor", false, "require that the track is a PV contributor"};
  Configurable<bool> applyItsPid{"applyItsPid", false, "apply ITS PID"};
  Configurable<int> minItsNclusters{"minItsNclusters", 5, "minimum number of ITS clusters"};
  Configurable<int> minTpcNcrossedRows{"minTpcNcrossedRows", 100, "minimum number of TPC crossed pad rows"};
  Configurable<double> minChiSquareTpc{"minChiSquareTpc", 0.0, "minimum TPC chi^2/Ncls"};
  Configurable<double> maxChiSquareTpc{"maxChiSquareTpc", 4.0, "maximum TPC chi^2/Ncls"};
  Configurable<double> maxChiSquareIts{"maxChiSquareIts", 36.0, "maximum ITS chi^2/Ncls"};
  Configurable<double> minPt{"minPt", 0.3, "minimum pt of the tracks"};
  Configurable<double> minEta{"minEta", -0.8, "minimum eta"};
  Configurable<double> maxEta{"maxEta", +0.8, "maximum eta"};
  Configurable<double> maxDcaxy{"maxDcaxy", 0.05, "Maximum DCAxy"};
  Configurable<double> maxDcaz{"maxDcaz", 0.05, "Maximum DCAz"};
  Configurable<double> minNsigmaTpc{"minNsigmaTpc", -3.0, "Minimum nsigma TPC"};
  Configurable<double> maxNsigmaTpc{"maxNsigmaTpc", +3.0, "Maximum nsigma TPC"};
  Configurable<double> minNsigmaTof{"minNsigmaTof", -3.0, "Minimum nsigma TOF"};
  Configurable<double> maxNsigmaTof{"maxNsigmaTof", +3.5, "Maximum nsigma TOF"};
  Configurable<double> ptMaxItsPidProt{"ptMaxItsPidProt", 1.0, "maximum pt for ITS PID for protons"};
  Configurable<double> ptMaxItsPidDeut{"ptMaxItsPidDeut", 1.0, "maximum pt for ITS PID for deuterons"};
  Configurable<double> ptMaxItsPidHel{"ptMaxItsPidHel", 1.0, "maximum pt for ITS PID for helium"};
  Configurable<double> nSigmaItsMin{"nSigmaItsMin", -3.0, "nSigmaITS min"};
  Configurable<double> nSigmaItsMax{"nSigmaItsMax", +3.0, "nSigmaITS max"};
  Configurable<bool> setMCDefaultItsParams{"setMCDefaultItsParams", true, "set MC default parameters"};
  Configurable<bool> cfgCompensatePIDinTracking{"cfgCompensatePIDinTracking", false, "If true, divide tpcInnerParam by the electric charge"};
  Configurable<std::array<double, 5>> cfgBetheBlochParams{"cfgBetheBlochParams", {0.6539, 1.591, 0.8225, 2.363, 0.09}, "TPC Bethe-Bloch parameterisation for He3"};

  // Configuration parameters for CCDB access and reweighting input files
  Configurable<bool> applyReweighting{"applyReweighting", true, "enable reweighting for efficiency"};
  Configurable<std::string> urlToCcdb{"urlToCcdb", "http://alice-ccdb.cern.ch", "url of the personal ccdb"};
  Configurable<std::string> pathToFile{"pathToFile", "Users/a/alcaliva/reweightingHistogramsAntipInJet", "path to file"};
  Configurable<std::string> weightsProton{"weightsProton", "", "weightsProton"};
  Configurable<std::string> weightsLambda{"weightsLambda", "", "weightsLambda"};
  Configurable<std::string> weightsSigma{"weightsSigma", "", "weightsSigma"};
  Configurable<std::string> weightsXi{"weightsXi", "", "weightsXi"};
  Configurable<std::string> weightsOmega{"weightsOmega", "", "weightsOmega"};
  Configurable<std::string> weightsJet{"weightsJet", "", "weightsJet"};
  Configurable<std::string> weightsUe{"weightsUe", "", "weightsUe"};

  // Number of events
  Configurable<int> shrinkInterval{"shrinkInterval", 1000, "variable that controls how often shrinking happens"};

  // Coalescence momentum
  Configurable<double> coalescenceMomentum{"coalescenceMomentum", 0.15, "p0 (GeV/c)"};

  // Reweighting histograms
  TH1F* primaryAntiprotons;
  TH1F* primaryAntiLambda;
  TH1F* primaryAntiSigma;
  TH1F* primaryAntiXi;
  TH1F* primaryAntiOmega;
  TH1F* antiprotonsInsideJets;
  TH1F* antiprotonsPerpCone;

  // CCDB manager service for accessing condition data
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  // Direct interface to the CCDB API for manual data access
  o2::ccdb::CcdbApi ccdbApi;

  // Instantiate the main Zorro processing object and define an output to store summary information
  Zorro zorro;
  OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};

  // Utility object for jet background subtraction methods
  JetBkgSubUtils backgroundSub;

  // Initialize ITS PID Response object
  o2::aod::ITSResponse itsResponse;

  // Initialize CCDB access and histogram registry for Zorro processing
  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (cfgSkimmedProcessing) {
      zorro.initCCDB(ccdb.service, bc.runNumber(), bc.timestamp(), triggerList);
      zorro.populateHistRegistry(registryData, bc.runNumber());
    }
  }

  void init(InitContext const&)
  {
    // Set summary object if processing skimmed data
    if (cfgSkimmedProcessing) {
      zorroSummary.setObject(zorro.getZorroSummary());
    }

    // Set default MC parametrization for ITS response
    if (setMCDefaultItsParams) {
      itsResponse.setMCDefaultParameters();
    }

    // Initialize random seed using high-resolution clock to ensure unique sequences across parallel Grid jobs
    auto timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    mRand.SetSeed(timeSeed);

    // Load reweighting histograms from CCDB if antinuclei efficiency processing is enabled
    if (doprocessAntinucleiEfficiency || doprocessJetsMCgen || doprocessJetsMCrec) {
      ccdb->setURL(urlToCcdb.value);
      ccdb->setCaching(true);
      ccdb->setLocalObjectValidityChecking();
      ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
      ccdb->setFatalWhenNull(false);
      getReweightingHistograms(ccdb, TString(pathToFile), TString(weightsProton), TString(weightsLambda), TString(weightsSigma), TString(weightsXi), TString(weightsOmega), TString(weightsJet), TString(weightsUe));
    }

    // Binning
    double min = 0.0;
    double max = 6.0;
    int nbins = 120;

    // Quality control histograms for jet/UE topology and multiplicity
    if (doprocessQC) {
      registryQC.add("deltaEta_deltaPhi_jet", "deltaEta_deltaPhi_jet", HistType::kTH2F, {{200, -0.5, 0.5, "#Delta#eta"}, {200, 0, PIHalf, "#Delta#phi"}});
      registryQC.add("deltaEta_deltaPhi_ue", "deltaEta_deltaPhi_ue", HistType::kTH2F, {{200, -0.5, 0.5, "#Delta#eta"}, {200, 0, PIHalf, "#Delta#phi"}});
      registryQC.add("eta_phi_jet", "eta_phi_jet", HistType::kTH2F, {{200, -0.5, 0.5, "#eta_{jet}"}, {200, 0, TwoPI, "#phi_{jet}"}});
      registryQC.add("eta_phi_ue", "eta_phi_ue", HistType::kTH2F, {{200, -0.5, 0.5, "#eta_{UE}"}, {200, 0, TwoPI, "#phi_{UE}"}});
      registryQC.add("NchJetCone", "NchJetCone", HistType::kTH1F, {{100, 0, 100, "#it{N}_{ch}"}});
      registryQC.add("NchJet", "NchJet", HistType::kTH1F, {{100, 0, 100, "#it{N}_{ch}"}});
      registryQC.add("NchUE", "NchUE", HistType::kTH1F, {{100, 0, 100, "#it{N}_{ch}"}});
      registryQC.add("sumPtJetCone", "sumPtJetCone", HistType::kTH1F, {{500, 0, 50, "#it{p}_{T} (GeV/#it{c})"}});
      registryQC.add("sumPtJet", "sumPtJet", HistType::kTH1F, {{500, 0, 50, "#it{p}_{T} (GeV/#it{c})"}});
      registryQC.add("sumPtUE", "sumPtUE", HistType::kTH1F, {{500, 0, 50, "#it{p}_{T} (GeV/#it{c})"}});
      registryQC.add("nJetsFound", "nJetsFound", HistType::kTH1F, {{50, 0, 50, "#it{n}_{Jet}"}});
      registryQC.add("nJetsInAcceptance", "nJetsInAcceptance", HistType::kTH1F, {{50, 0, 50, "#it{n}_{Jet}"}});
      registryQC.add("nJetsSelectedHighPt", "nJetsSelectedHighPt", HistType::kTH1F, {{50, 0, 50, "#it{n}_{Jet}"}});
      registryQC.add("jetPtDifference", "jetPtDifference", HistType::kTH1F, {{200, -1, 1, "#Deltap_{T}^{jet}"}});
      registryQC.add("ptDistributionJetCone", "ptDistributionJetCone", HistType::kTH1F, {{2000, 0, 200, "#it{p}_{T} (GeV/#it{c})"}});
      registryQC.add("ptDistributionJet", "ptDistributionJet", HistType::kTH1F, {{2000, 0, 200, "#it{p}_{T} (GeV/#it{c})"}});
    }

    // Multiplicity histograms: check charged-particle multiplicity in events with selected jets (Run 3) or ptTrigger > threshold (Run 2-like)
    if (doprocessMultEvents) {
      registryMult.add("multiplicityEvtsPtLeading", "multiplicityEvtsPtLeading", HistType::kTH1F, {{1000, 0, 1000, "#it{N}_{ch}"}});
      registryMult.add("multiplicityEvtsWithJet", "multiplicityEvtsWithJet", HistType::kTH1F, {{1000, 0, 1000, "#it{N}_{ch}"}});
    }

    // Histograms for real data
    if (doprocessData) {

      // Event counters
      registryData.add("number_of_events_data", "number of events in data", HistType::kTH1F, {{20, 0, 20, "counter"}});

      // Configuration
      registryData.add("settingData", "settingData", HistType::kTH2F, {{100, 0.0, 50.0, "min #it{p}^{jet}_{T} (GeV/#it{c})"}, {20, 0.0, 1.0, "#it{R}_{jet}"}});

      // Jet effective area over piR^2
      registryData.add("jetEffectiveAreaOverPiR2", "jet effective area / piR^2", HistType::kTH1F, {{2000, 0, 2, "Area/#piR^{2}"}});

      // angle between track and jet axis
      registryData.add("theta_track_jet", "theta_track_jet", HistType::kTH2F, {{100, 0, 100, "#it{p}^{jet}_{T} (GeV/#it{c})"}, {400, 0, 20.0, "#theta_{track-jet} (deg)"}});

      // Antiprotons
      registryData.add("antiproton_jet_tpc", "antiproton_jet_tpc", HistType::kTH2F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{TPC}"}});
      registryData.add("antiproton_jet_tof", "antiproton_jet_tof", HistType::kTH2F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{TOF}"}});
      registryData.add("antiproton_ue_tpc", "antiproton_ue_tpc", HistType::kTH2F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{TPC}"}});
      registryData.add("antiproton_ue_tof", "antiproton_ue_tof", HistType::kTH2F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{TOF}"}});
      registryData.add("antiproton_dca_jet", "antiproton_dca_jet", HistType::kTH2F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}, {200, -1.0, 1.0, "DCA_{xy} (cm)"}});
      registryData.add("antiproton_dca_ue", "antiproton_dca_ue", HistType::kTH2F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}, {200, -1.0, 1.0, "DCA_{xy} (cm)"}});

      // Antideuterons
      registryData.add("antideuteron_jet_tpc", "antideuteron_jet_tpc", HistType::kTH2F, {{nbins, 2 * min, 2 * max, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{TPC}"}});
      registryData.add("antideuteron_jet_tof", "antideuteron_jet_tof", HistType::kTH2F, {{nbins, 2 * min, 2 * max, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{TOF}"}});
      registryData.add("antideuteron_ue_tpc", "antideuteron_ue_tpc", HistType::kTH2F, {{nbins, 2 * min, 2 * max, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{TPC}"}});
      registryData.add("antideuteron_ue_tof", "antideuteron_ue_tof", HistType::kTH2F, {{nbins, 2 * min, 2 * max, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{TOF}"}});

      // Deuterons
      registryData.add("deuteron_jet_tpc", "deuteron_jet_tpc", HistType::kTH2F, {{nbins, 2 * min, 2 * max, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{TPC}"}});
      registryData.add("deuteron_jet_tof", "deuteron_jet_tof", HistType::kTH2F, {{nbins, 2 * min, 2 * max, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{TOF}"}});
      registryData.add("deuteron_ue_tpc", "deuteron_ue_tpc", HistType::kTH2F, {{nbins, 2 * min, 2 * max, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{TPC}"}});
      registryData.add("deuteron_ue_tof", "deuteron_ue_tof", HistType::kTH2F, {{nbins, 2 * min, 2 * max, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{TOF}"}});

      // Antihelium-3
      registryData.add("antihelium3_jet_tpc", "antihelium3_jet_tpc", HistType::kTH2F, {{nbins, 3 * min, 3 * max, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{TPC}"}});
      registryData.add("antihelium3_ue_tpc", "antihelium3_ue_tpc", HistType::kTH2F, {{nbins, 3 * min, 3 * max, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{TPC}"}});

      // Helium-3
      registryData.add("helium3_jet_tpc", "helium3_jet_tpc", HistType::kTH2F, {{nbins, 3 * min, 3 * max, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{TPC}"}});
      registryData.add("helium3_ue_tpc", "helium3_ue_tpc", HistType::kTH2F, {{nbins, 3 * min, 3 * max, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{TPC}"}});

      // nsigmaITS for antiproton candidates
      registryData.add("antiproton_nsigma_its_data", "antiproton_nsigma_its_data", HistType::kTH2F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{ITS}"}});

      // custom nsigma for He3 - needed for 24 pp data
      registryData.add("tpcsignal_data", "Specific energy loss", HistType::kTH2F, {{600, -6., 6., "#it{p} (GeV/#it{c})"}, {1400, 0, 1400, "d#it{E} / d#it{X} (a. u.)"}});
      registryData.add("antihelium3_jet_tpc_custom", "antihelium3_jet_tpc_custom", HistType::kTH2F, {{nbins, 3 * min, 3 * max, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{TPC} custom"}});
      registryData.add("antihelium3_ue_tpc_custom", "antihelium3_ue_tpc_custom", HistType::kTH2F, {{nbins, 3 * min, 3 * max, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{TPC} custom"}});
    }

    // Generated antiproton spectra in jets and UE from MC truth
    if (doprocessJetsMCgen) {

      // Event counter
      registryMC.add("genEvents", "number of generated events in mc", HistType::kTH1F, {{10, 0, 10, "counter"}});
      registryMC.add("genJets", "number of generated jets", HistType::kTH1F, {{10, 0, 10, "counter"}});

      // Size of particle array
      registryMC.add("sizeParticleArray", "number of particles", HistType::kTH1F, {{1000, 0, 1000, "#it{N}_{part}"}});

      // Generated spectra of antiprotons
      registryMC.add("antiproton_gen_jet", "antiproton_gen_jet", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("antiproton_gen_ue", "antiproton_gen_ue", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("antiproton_gen_full", "antiproton_gen_full", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});

      // Normalization histogram
      registryMC.add("antiproton_deltay_deltaphi_jet", "antiproton_deltay_deltaphi_jet", HistType::kTH2F, {{2000, -1.0, 1.0, "#Delta#it{y}"}, {2000, 0.0, 2.0, "#Delta#phi"}});
      registryMC.add("antiproton_deltay_deltaphi_ue", "antiproton_deltay_deltaphi_ue", HistType::kTH2F, {{2000, -1.0, 1.0, "#Delta#it{y}"}, {2000, 0.0, 2.0, "#Delta#phi"}});

      // 2d kinematic distributions (eta,pt) in jets and UE
      registryMC.add("antiproton_eta_pt_jet", "antiproton_eta_pt_jet", HistType::kTH2F, {{500, 0.0, 5.0, "#it{p}_{T} (GeV/#it{c})"}, {18, -0.9, 0.9, "#eta"}});
      registryMC.add("antiproton_eta_pt_ue", "antiproton_eta_pt_ue", HistType::kTH2F, {{500, 0.0, 5.0, "#it{p}_{T} (GeV/#it{c})"}, {18, -0.9, 0.9, "#eta"}});
    }

    // Reconstructed antiproton spectra in jets and UE (MC-matched) with TPC/TOF PID
    if (doprocessJetsMCrec) {

      // Event counter
      registryMC.add("recEvents", "number of reconstructed events in mc", HistType::kTH1F, {{20, 0, 20, "counter"}});
      registryMC.add("recJets", "number of reconstructed jets", HistType::kTH1F, {{10, 0, 10, "counter"}});

      // Configuration
      registryMC.add("settingMC", "settingMC", HistType::kTH2F, {{100, 0.0, 50.0, "min #it{p}^{jet}_{T} (GeV/#it{c})"}, {20, 0.0, 1.0, "#it{R}_{jet}"}});

      // Reconstructed spectra of antiprotons
      registryMC.add("antiproton_rec_tpc_jet", "antiproton_rec_tpc_jet", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("antiproton_rec_tof_jet", "antiproton_rec_tof_jet", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("antiproton_rec_tpc_ue", "antiproton_rec_tpc_ue", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("antiproton_rec_tof_ue", "antiproton_rec_tof_ue", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("antiproton_rec_tpc_full", "antiproton_rec_tpc_full", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("antiproton_rec_tof_full", "antiproton_rec_tof_full", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});

      // Fraction of primary antiprotons
      registryMC.add("antiproton_prim_jet", "antiproton_prim_jet", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("antiproton_incl_jet", "antiproton_incl_jet", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("antiproton_prim_ue", "antiproton_prim_ue", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("antiproton_incl_ue", "antiproton_incl_ue", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});

      // DCA templates
      registryMC.add("antiproton_prim_dca_jet", "antiproton_prim_dca_jet", HistType::kTH2F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}, {200, -1, 1, "DCA_{xy} (cm)"}});
      registryMC.add("antiproton_prim_dca_ue", "antiproton_prim_dca_ue", HistType::kTH2F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}, {200, -1, 1, "DCA_{xy} (cm)"}});
      registryMC.add("antiproton_all_dca_jet", "antiproton_all_dca_jet", HistType::kTH2F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}, {200, -1, 1, "DCA_{xy} (cm)"}});
      registryMC.add("antiproton_all_dca_ue", "antiproton_all_dca_ue", HistType::kTH2F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}, {200, -1, 1, "DCA_{xy} (cm)"}});

      // nsigmaTOF for antiprotons
      registryMC.add("antiproton_nsigma_tof_jet_mc", "antiproton_nsigma_tof_jet_mc", HistType::kTH2F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{TOF}"}});
    }

    // Efficiency of antinuclei
    if (doprocessAntinucleiEfficiency) {

      // Event counter MC
      registryMC.add("number_of_events_mc_nuclei_efficiency", "number of events in mc", HistType::kTH1F, {{20, 0, 20, "counter"}});

      // Generated spectra of (anti)protons
      registryMC.add("antip_gen_jet", "antip_gen_jet", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("antip_gen_ue", "antip_gen_ue", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});

      // Generated spectra of (anti)deuterons
      registryMC.add("deuteron_gen_jet", "deuteron_gen_jet", HistType::kTH1F, {{nbins, 2 * min, 2 * max, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("deuteron_gen_ue", "deuteron_gen_ue", HistType::kTH1F, {{nbins, 2 * min, 2 * max, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("antideuteron_gen_jet", "antideuteron_gen_jet", HistType::kTH1F, {{nbins, 2 * min, 2 * max, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("antideuteron_gen_ue", "antideuteron_gen_ue", HistType::kTH1F, {{nbins, 2 * min, 2 * max, "#it{p}_{T} (GeV/#it{c})"}});

      // Generated spectra of (anti)helium3
      registryMC.add("helium3_gen_jet", "helium3_gen_jet", HistType::kTH1F, {{nbins, 3 * min, 3 * max, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("helium3_gen_ue", "helium3_gen_ue", HistType::kTH1F, {{nbins, 3 * min, 3 * max, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("antihelium3_gen_jet", "antihelium3_gen_jet", HistType::kTH1F, {{nbins, 3 * min, 3 * max, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("antihelium3_gen_ue", "antihelium3_gen_ue", HistType::kTH1F, {{nbins, 3 * min, 3 * max, "#it{p}_{T} (GeV/#it{c})"}});

      // Reconstructed spectra of antiprotons
      registryMC.add("antip_rec_tpc_jet", "antip_rec_tpc_jet", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("antip_rec_tof_jet", "antip_rec_tof_jet", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("antip_rec_tpc_ue", "antip_rec_tpc_ue", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("antip_rec_tof_ue", "antip_rec_tof_ue", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});

      // Reconstructed spectra of (anti)deuterons
      registryMC.add("deuteron_rec_tpc_jet", "deuteron_rec_tpc_jet", HistType::kTH1F, {{nbins, 2 * min, 2 * max, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("deuteron_rec_tof_jet", "deuteron_rec_tof_jet", HistType::kTH1F, {{nbins, 2 * min, 2 * max, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("deuteron_rec_tpc_ue", "deuteron_rec_tpc_ue", HistType::kTH1F, {{nbins, 2 * min, 2 * max, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("deuteron_rec_tof_ue", "deuteron_rec_tof_ue", HistType::kTH1F, {{nbins, 2 * min, 2 * max, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("antideuteron_rec_tpc_jet", "antideuteron_rec_tpc_jet", HistType::kTH1F, {{nbins, 2 * min, 2 * max, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("antideuteron_rec_tof_jet", "antideuteron_rec_tof_jet", HistType::kTH1F, {{nbins, 2 * min, 2 * max, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("antideuteron_rec_tpc_ue", "antideuteron_rec_tpc_ue", HistType::kTH1F, {{nbins, 2 * min, 2 * max, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("antideuteron_rec_tof_ue", "antideuteron_rec_tof_ue", HistType::kTH1F, {{nbins, 2 * min, 2 * max, "#it{p}_{T} (GeV/#it{c})"}});

      // Reconstructed spectra of (anti)helium3
      registryMC.add("helium3_rec_tpc_jet", "helium3_rec_tpc_jet", HistType::kTH1F, {{nbins, 3 * min, 3 * max, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("helium3_rec_tpc_ue", "helium3_rec_tpc_ue", HistType::kTH1F, {{nbins, 3 * min, 3 * max, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("antihelium3_rec_tpc_jet", "antihelium3_rec_tpc_jet", HistType::kTH1F, {{nbins, 3 * min, 3 * max, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("antihelium3_rec_tpc_ue", "antihelium3_rec_tpc_ue", HistType::kTH1F, {{nbins, 3 * min, 3 * max, "#it{p}_{T} (GeV/#it{c})"}});

      // Generated spectra needed for reweighting
      registryMC.add("protonBar", "protonBar", HistType::kTH1F, {{1000, 0, 10, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("lambdaBar", "lambdaBar", HistType::kTH1F, {{1000, 0, 10, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("xiBar", "xiBar", HistType::kTH1F, {{1000, 0, 10, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("omegaBar", "omegaBar", HistType::kTH1F, {{1000, 0, 10, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("sigmaBar", "sigmaBar", HistType::kTH1F, {{1000, 0, 10, "#it{p}_{T} (GeV/#it{c})"}});

      // nsigmaITS for antiproton candidates
      registryMC.add("antiproton_nsigma_its_mc", "antiproton_nsigma_its_mc", HistType::kTH2F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{ITS}"}});

      // Systematics on the fraction of primary antiprotons
      registryMC.add("antip_prim_pythia", "antip_prim_pythia", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("antip_prim_std", "antip_prim_std", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("antip_prim_up", "antip_prim_up", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("antip_prim_low", "antip_prim_low", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});

      registryMC.add("antip_sec_pythia", "antip_sec_pythia", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("antip_sec_std", "antip_sec_std", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("antip_sec_up", "antip_sec_up", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("antip_sec_low", "antip_sec_low", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});
    }

    // Coalescence
    if (doprocessCoalescence) {
      registryMC.add("genEventsCoalescence", "genEventsCoalescence", HistType::kTH1F, {{20, 0, 20, "counter"}});
      registryMC.add("antideuteron_coal_fullEvent", "antideuteron_coal_fullEvent", HistType::kTH1F, {{nbins, 2 * min, 2 * max, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("antideuteron_coal_jet", "antideuteron_coal_jet", HistType::kTH1F, {{nbins, 2 * min, 2 * max, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("antideuteron_coal_ue", "antideuteron_coal_ue", HistType::kTH1F, {{nbins, 2 * min, 2 * max, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("antiproton_coal_fullEvent", "antiproton_coal_fullEvent", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("antiproton_coal_jet", "antiproton_coal_jet", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("antiproton_coal_ue", "antiproton_coal_ue", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});
    }

    // Systematic uncertainties (Data)
    if (doprocessSystData) {
      registryData.add("number_of_events_data_syst", "event counter", HistType::kTH1F, {{20, 0, 20, "counter"}});

      registryData.add("antiproton_tpc_syst", "antiproton_tpc_syst", HistType::kTH3F, {{50, 0, 50, "systematic uncertainty"}, {nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{TPC}"}});
      registryData.add("antiproton_tof_syst", "antiproton_tof_syst", HistType::kTH3F, {{50, 0, 50, "systematic uncertainty"}, {nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{TOF}"}});
      registryData.add("antideuteron_tpc_syst", "antideuteron_tpc_syst", HistType::kTH3F, {{50, 0, 50, "systematic uncertainty"}, {nbins, 2 * min, 2 * max, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{TPC}"}});
      registryData.add("antideuteron_tof_syst", "antideuteron_tof_syst", HistType::kTH3F, {{50, 0, 50, "systematic uncertainty"}, {nbins, 2 * min, 2 * max, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{TOF}"}});
      registryData.add("antihelium3_tpc_syst", "antihelium3_tpc_syst", HistType::kTH3F, {{50, 0, 50, "systematic uncertainty"}, {nbins, 3 * min, 3 * max, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{TPC}"}});
    }

    // Systematic uncertainties (MC)
    if (doprocessSystEff) {

      // Histograms for generated antiparticles
      registryMC.add("antiproton_gen_syst", "antiproton_gen_syst", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("antideuteron_gen_syst", "antideuteron_gen_syst", HistType::kTH1F, {{nbins, 2 * min, 2 * max, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("antihelium3_gen_syst", "antihelium3_gen_syst", HistType::kTH1F, {{nbins, 3 * min, 3 * max, "#it{p}_{T} (GeV/#it{c})"}});

      // Histograms for reconstructed antiparticles
      registryMC.add("antiproton_rec_tpc_syst", "antiproton_rec_tpc_syst", HistType::kTH2F, {{50, 0, 50, "systematic uncertainty"}, {nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("antiproton_rec_tof_syst", "antiproton_rec_tof_syst", HistType::kTH2F, {{50, 0, 50, "systematic uncertainty"}, {nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("antideuteron_rec_tpc_syst", "antideuteron_rec_tpc_syst", HistType::kTH2F, {{50, 0, 50, "systematic uncertainty"}, {nbins, 2 * min, 2 * max, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("antideuteron_rec_tof_syst", "antideuteron_rec_tof_syst", HistType::kTH2F, {{50, 0, 50, "systematic uncertainty"}, {nbins, 2 * min, 2 * max, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("antihelium3_rec_tpc_syst", "antihelium3_rec_tpc_syst", HistType::kTH2F, {{50, 0, 50, "systematic uncertainty"}, {nbins, 3 * min, 3 * max, "#it{p}_{T} (GeV/#it{c})"}});

      // Histograms for primary antiprotons
      registryMC.add("antiproton_incl_syst", "antiproton_incl_syst", HistType::kTH2F, {{50, 0, 50, "systematic uncertainty"}, {nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("antiproton_prim_syst", "antiproton_prim_syst", HistType::kTH2F, {{50, 0, 50, "systematic uncertainty"}, {nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});
    }

    // Correlation analysis
    if (doprocessCorr) {

      // Axes definitions for multidimensional histogram binning
      const AxisSpec multiplicityAxis{100, 0.0, 100.0, "multiplicity percentile"};
      const AxisSpec ptPerNucleonAxis{5, 0.4, 0.9, "{p}_{T}/A (GeV/#it{c})"};
      const AxisSpec nAntideuteronsAxis{10, 0.0, 10.0, "N_{#bar{d}}"};
      const AxisSpec nAntiprotonsAxis{10, 0.0, 10.0, "N_{#bar{p}}"};
      const AxisSpec nBarD2Axis{100, 0.0, 100.0, "N_{#bar{d}}^{i} #times N_{#bar{d}}^{j}"};
      const AxisSpec nBarP2Axis{100, 0.0, 100.0, "N_{#bar{p}}^{i} #times N_{#bar{p}}^{j}"};
      const AxisSpec nBarDnBarPAxis{100, 0.0, 100.0, "N_{#bar{d}}^{i} #times N_{#bar{p}}^{j}"};
      const AxisSpec subsampleAxis{nSubsamples, 0, static_cast<double>(nSubsamples), "Subsample Index"};

      // Event counter
      registryCorr.add("eventCounter", "number of events", HistType::kTH1F, {{20, 0, 20, "counter"}});
      registryCorr.add("eventCounter_centrality_fullEvent", "Number of events per centrality (Full Event)", HistType::kTH2F, {multiplicityAxis, subsampleAxis});
      // registryCorr.add("eventCounter_centrality_jet", "Number of events per centrality (Jet)", HistType::kTH1F, {multiplicityAxis});
      // registryCorr.add("eventCounter_centrality_ue", "Number of events per centrality (Underlying Event)",  HistType::kTH1F, {multiplicityAxis});

      // Correlation histograms: antiproton vs. antideuteron number vs. event multiplicity
      // registryCorr.add("rho_jet", "rho_jet", HistType::kTH3F, {nAntideuteronsAxis, nAntiprotonsAxis, multiplicityAxis});
      // registryCorr.add("rho_ue", "rho_ue", HistType::kTH3F, {nAntideuteronsAxis, nAntiprotonsAxis, multiplicityAxis});
      registryCorr.add("rho_fullEvent", "rho_fullEvent", HistType::kTHnSparseD, {nAntideuteronsAxis, nAntiprotonsAxis, multiplicityAxis, subsampleAxis});

      // Correlation histograms: net antiproton vs. net antideuteron numbers
      // registryCorr.add("rho_netP_netD_jet", "rho_netP_netD_jet", HistType::kTH2F, {nAntideuteronsAxis, nAntiprotonsAxis});
      // registryCorr.add("rho_netP_netD_ue", "rho_netP_netD_ue", HistType::kTH2F, {nAntideuteronsAxis, nAntiprotonsAxis});
      registryCorr.add("rho_netP_netD_fullEvent", "rho_netP_netD_fullEvent", HistType::kTH3F, {nAntideuteronsAxis, nAntiprotonsAxis, subsampleAxis});

      // Efficiency histograms jet
      // registryCorr.add("q1d_jet", "q1d_jet", HistType::kTH3F, {nAntideuteronsAxis, ptPerNucleonAxis, multiplicityAxis});
      // registryCorr.add("q1p_jet", "q1p_jet", HistType::kTH3F, {nAntiprotonsAxis, ptPerNucleonAxis, multiplicityAxis});
      // registryCorr.add("q1d_square_jet", "q1d_square_jet", HistType::kTHnSparseD, {ptPerNucleonAxis, ptPerNucleonAxis, nBarD2Axis, multiplicityAxis});
      // registryCorr.add("q1p_square_jet", "q1p_square_jet", HistType::kTHnSparseD, {ptPerNucleonAxis, ptPerNucleonAxis, nBarP2Axis, multiplicityAxis});
      // registryCorr.add("q1d_q1p_jet", "q1d_q1p_jet", HistType::kTHnSparseD, {ptPerNucleonAxis, ptPerNucleonAxis, nBarDnBarPAxis, multiplicityAxis});

      // Efficiency histograms UE
      // registryCorr.add("q1d_ue", "q1d_ue", HistType::kTH3F, {nAntideuteronsAxis, ptPerNucleonAxis, multiplicityAxis});
      // registryCorr.add("q1p_ue", "q1p_ue", HistType::kTH3F, {nAntiprotonsAxis, ptPerNucleonAxis, multiplicityAxis});
      // registryCorr.add("q1d_square_ue", "q1d_square_ue", HistType::kTHnSparseD, {ptPerNucleonAxis, ptPerNucleonAxis, nBarD2Axis, multiplicityAxis});
      // registryCorr.add("q1p_square_ue", "q1p_square_ue", HistType::kTHnSparseD, {ptPerNucleonAxis, ptPerNucleonAxis, nBarP2Axis, multiplicityAxis});
      // registryCorr.add("q1d_q1p_ue", "q1d_q1p_ue", HistType::kTHnSparseD, {ptPerNucleonAxis, ptPerNucleonAxis, nBarDnBarPAxis, multiplicityAxis});

      // Efficiency histograms full event
      registryCorr.add("q1d_fullEvent", "q1d_fullEvent", HistType::kTHnSparseD, {nAntideuteronsAxis, ptPerNucleonAxis, multiplicityAxis, subsampleAxis});
      registryCorr.add("q1p_fullEvent", "q1p_fullEvent", HistType::kTHnSparseD, {nAntiprotonsAxis, ptPerNucleonAxis, multiplicityAxis, subsampleAxis});
      registryCorr.add("q1d_square_fullEvent", "q1d_square_fullEvent", HistType::kTHnSparseD, {ptPerNucleonAxis, ptPerNucleonAxis, nBarD2Axis, multiplicityAxis, subsampleAxis});
      registryCorr.add("q1p_square_fullEvent", "q1p_square_fullEvent", HistType::kTHnSparseD, {ptPerNucleonAxis, ptPerNucleonAxis, nBarP2Axis, multiplicityAxis, subsampleAxis});
      registryCorr.add("q1d_q1p_fullEvent", "q1d_q1p_fullEvent", HistType::kTHnSparseD, {ptPerNucleonAxis, ptPerNucleonAxis, nBarDnBarPAxis, multiplicityAxis, subsampleAxis});
    }
      
      // Systematic uncertainties on correlation analysis
      if (doprocessCorrSyst) {
          
          // Axes definitions for multidimensional histogram binning
          const AxisSpec multiplicityAxis{100, 0.0, 100.0, "multiplicity percentile"};
          const AxisSpec ptPerNucleonAxis{5, 0.4, 0.9, "{p}_{T}/A (GeV/#it{c})"};
          const AxisSpec nAntideuteronsAxis{10, 0.0, 10.0, "N_{#bar{d}}"};
          const AxisSpec nAntiprotonsAxis{10, 0.0, 10.0, "N_{#bar{p}}"};
          const AxisSpec nBarD2Axis{100, 0.0, 100.0, "N_{#bar{d}}^{i} #times N_{#bar{d}}^{j}"};
          const AxisSpec nBarP2Axis{100, 0.0, 100.0, "N_{#bar{p}}^{i} #times N_{#bar{p}}^{j}"};
          const AxisSpec nBarDnBarPAxis{100, 0.0, 100.0, "N_{#bar{d}}^{i} #times N_{#bar{p}}^{j}"};
          const AxisSpec systAxis{nSyst, 0, static_cast<double>(nSyst), "Systematic Variation Index"};

          registryCorr.add("eventCounter_syst", "number of events syst", HistType::kTH1F, {{20, 0, 20, "counter"}});
          registryCorr.add("eventCounter_centrality_fullEvent_syst", "Number of events per centrality (Full Event) Syst", HistType::kTH2F, {multiplicityAxis, systAxis});

          // Correlation histograms
          registryCorr.add("rho_fullEvent_syst", "rho_fullEvent_syst", HistType::kTHnSparseD, {nAntideuteronsAxis, nAntiprotonsAxis, multiplicityAxis, systAxis});
          registryCorr.add("rho_netP_netD_fullEvent_syst", "rho_netP_netD_fullEvent_syst", HistType::kTH3F, {nAntideuteronsAxis, nAntiprotonsAxis, systAxis});

          // Efficiency histograms full event
          registryCorr.add("q1d_fullEvent_syst", "q1d_fullEvent_syst", HistType::kTHnSparseD, {nAntideuteronsAxis, ptPerNucleonAxis, multiplicityAxis, systAxis});
          registryCorr.add("q1p_fullEvent_syst", "q1p_fullEvent_syst", HistType::kTHnSparseD, {nAntiprotonsAxis, ptPerNucleonAxis, multiplicityAxis, systAxis});
          registryCorr.add("q1d_square_fullEvent_syst", "q1d_square_fullEvent_syst", HistType::kTHnSparseD, {ptPerNucleonAxis, ptPerNucleonAxis, nBarD2Axis, multiplicityAxis, systAxis});
          registryCorr.add("q1p_square_fullEvent_syst", "q1p_square_fullEvent_syst", HistType::kTHnSparseD, {ptPerNucleonAxis, ptPerNucleonAxis, nBarP2Axis, multiplicityAxis, systAxis});
          registryCorr.add("q1d_q1p_fullEvent_syst", "q1d_q1p_fullEvent_syst", HistType::kTHnSparseD, {ptPerNucleonAxis, ptPerNucleonAxis, nBarDnBarPAxis, multiplicityAxis, systAxis});
        }
  }

  void getReweightingHistograms(o2::framework::Service<o2::ccdb::BasicCCDBManager> const& ccdbObj, TString filepath, TString antip, TString antilambda, TString antisigma, TString antixi, TString antiomega, TString jet, TString ue)
  {
    TList* list = ccdbObj->get<TList>(filepath.Data());
    if (!list) {
      LOGP(error, "Could not retrieve the list from CCDB");
      return;
    }

    // Get reweighting histograms for primary fraction
    primaryAntiprotons = static_cast<TH1F*>(list->FindObject(antip));
    primaryAntiLambda = static_cast<TH1F*>(list->FindObject(antilambda));
    primaryAntiSigma = static_cast<TH1F*>(list->FindObject(antisigma));
    primaryAntiXi = static_cast<TH1F*>(list->FindObject(antixi));
    primaryAntiOmega = static_cast<TH1F*>(list->FindObject(antiomega));

    if (!primaryAntiprotons || !primaryAntiSigma || !primaryAntiLambda || !primaryAntiXi || !primaryAntiOmega) {
      LOGP(error, "Missing one or more reweighting histograms for primary fraction in CCDB list");
    }

    // Get reweighting histograms for efficiency
    antiprotonsInsideJets = static_cast<TH1F*>(list->FindObject(jet));
    antiprotonsPerpCone = static_cast<TH1F*>(list->FindObject(ue));
    if (!antiprotonsInsideJets || !antiprotonsPerpCone) {
      LOGP(error, "Missing one or more reweighting histograms for efficiency in CCDB list");
    }

    LOGP(info, "Successfully loaded reweighting histograms from CCDB path");
  }

  // Get first ancestor
  aod::McParticle getFirstAncestor(aod::McParticle const& particle, aod::McParticles const& mcParticles)
  {
    auto current = particle;

    // If already physical primary, return it
    if (current.isPhysicalPrimary())
      return current;

    while (current.has_mothers()) {
      auto motherId = current.mothersIds()[0];

      // Stop if motherId is invalid
      if (motherId < 0 || motherId >= mcParticles.size()) {
        break;
      }

      // Move up the chain
      current = mcParticles.iteratorAt(motherId);

      if (current.isPhysicalPrimary())
        break;
    }

    return current;
  }

  // Check if particle is a physical primary or a decay product of a heavy-flavor hadron
  bool isPhysicalPrimaryOrFromHF(aod::McParticle const& particle, aod::McParticles const& mcParticles)
  {
    // Keep only pi, K, p, e, mu
    int pdg = std::abs(particle.pdgCode());
    if (!(pdg == PDG_t::kPiPlus || pdg == PDG_t::kKPlus || pdg == PDG_t::kProton || pdg == PDG_t::kElectron || pdg == PDG_t::kMuonMinus))
      return false;

    // Constants for identifying heavy-flavor (charm and bottom) content from PDG codes
    static constexpr int CharmQuark = 4;
    static constexpr int BottomQuark = 5;
    static constexpr int Hundreds = 100;
    static constexpr int Thousands = 1000;

    // Check if particle is from heavy-flavor decay
    bool fromHF = false;
    if (particle.has_mothers()) {
      auto mother = mcParticles.iteratorAt(particle.mothersIds()[0]);
      int motherPdg = std::abs(mother.pdgCode());
      fromHF = (motherPdg / Hundreds == CharmQuark || motherPdg / Hundreds == BottomQuark || motherPdg / Thousands == CharmQuark || motherPdg / Thousands == BottomQuark);
    }

    // Select only physical primary particles or from heavy-flavor
    return (particle.isPhysicalPrimary() || fromHF);
  }

  // Evaluate proton–neutron coalescence for deuteron formation
  template <typename ReducedPart>
  bool passDeuteronCoalescence(const ReducedPart& p, const ReducedPart& n, double p0, TRandom3& mRand)
  {
    // Nucleon masses
    const double mp = o2::constants::physics::MassProton;
    const double mn = o2::constants::physics::MassNeutron;

    // Spin-statistical factor for deuteron formation (S = 1)
    static constexpr double SpinFactor = 3.0 / 4.0;

    // Require proton and neutron to be both matter or both antimatter
    const int signP = p.pdgCode / std::abs(p.pdgCode);
    const int signN = n.pdgCode / std::abs(n.pdgCode);
    if (signP != signN) {
      return false;
    }

    // Build on-shell nucleon four-momenta
    const double ep = std::sqrt(p.px * p.px + p.py * p.py + p.pz * p.pz + mp * mp);
    const double en = std::sqrt(n.px * n.px + n.py * n.py + n.pz * n.pz + mn * mn);

    LorentzVector p4p(p.px, p.py, p.pz, ep);
    LorentzVector p4n(n.px, n.py, n.pz, en);

    // Total four-momentum of the nucleon pair
    const LorentzVector p4tot = p4p + p4n;

    // Boost individual nucleons to the pair center-of-mass frame
    ROOT::Math::Boost boostToCM(p4tot.BoostToCM());
    const LorentzVector p4pcm = boostToCM(p4p);
    // const LorentzVector p4ncm = boostToCM(p4n);

    // Relative momentum magnitude in the CM frame
    const double relativeMomentum = p4pcm.P();

    // Momentum-space coalescence condition
    if (relativeMomentum > p0) {
      return false;
    }

    // Spin-statistical acceptance
    if (mRand.Uniform() > SpinFactor) {
      return false;
    }
    return true;
  }

  // Compute two transverse directions orthogonal to vector p
  void getPerpendicularDirections(const TVector3& p, TVector3& u1, TVector3& u2)
  {
    // Get momentum components
    double px = p.X();
    double py = p.Y();
    double pz = p.Z();

    // Precompute squared terms
    double px2 = px * px;
    double py2 = py * py;
    double pz2 = pz * pz;
    double pz4 = pz2 * pz2;

    // Case 1: vector along z-axis -> undefined perpendiculars
    if (px == 0 && py == 0) {
      u1.SetXYZ(0, 0, 0);
      u2.SetXYZ(0, 0, 0);
      return;
    }

    // Case 2: px = 0 -> avoid division by zero
    if (px == 0 && py != 0) {
      double ux = std::sqrt(py2 - pz4 / py2);
      double uy = -pz2 / py;
      u1.SetXYZ(ux, uy, pz);
      u2.SetXYZ(-ux, uy, pz);
      return;
    }

    // Case 3: py = 0 -> avoid division by zero
    if (py == 0 && px != 0) {
      double ux = -pz2 / px;
      double uy = std::sqrt(px2 - pz4 / px2);
      u1.SetXYZ(ux, uy, pz);
      u2.SetXYZ(ux, -uy, pz);
      return;
    }

    // General case: solve quadratic for perpendicular vectors
    double a = px2 + py2;
    double b = 2.0 * px * pz2;
    double c = pz4 - py2 * py2 - px2 * py2;
    double delta = b * b - 4.0 * a * c;

    // Invalid or degenerate solutions
    if (delta < 0 || a == 0) {
      u1.SetXYZ(0, 0, 0);
      u2.SetXYZ(0, 0, 0);
      return;
    }

    // Solution 1
    double u1x = (-b + std::sqrt(delta)) / (2.0 * a);
    double u1y = (-pz2 - px * u1x) / py;
    u1.SetXYZ(u1x, u1y, pz);

    // Solution 2
    double u2x = (-b - std::sqrt(delta)) / (2.0 * a);
    double u2y = (-pz2 - px * u2x) / py;
    u2.SetXYZ(u2x, u2y, pz);
  }

  // Compute delta phi
  double getDeltaPhi(double a1, double a2)
  {
    double deltaPhi(0);
    double phi1 = TVector2::Phi_0_2pi(a1);
    double phi2 = TVector2::Phi_0_2pi(a2);
    double diff = std::fabs(phi1 - phi2);

    if (diff <= PI)
      deltaPhi = diff;
    if (diff > PI)
      deltaPhi = TwoPI - diff;

    return deltaPhi;
  }

  // Find bin
  int findBin(const std::vector<double>& edges, double value)
  {
    auto it = std::upper_bound(edges.begin(), edges.end(), value);
    int index = static_cast<int>(it - edges.begin()) - 1;
    if (index < 0 || index >= static_cast<int>(edges.size()) - 1) {
      return -1; // value is out of bounds
    }
    return index;
  }

  // ITS hit
  template <typename TrackIts>
  bool hasITSHit(const TrackIts& track, int layer)
  {
    int ibit = layer - 1;
    return (track.itsClusterMap() & (1 << ibit));
  }

  // Single-track selection for particles inside jets
  template <typename JetTrack>
  bool passedTrackSelectionForJetReconstruction(const JetTrack& track)
  {
    static constexpr int MinTpcCr = 70;
    static constexpr double MaxChi2Tpc = 4.0;
    static constexpr double MaxChi2Its = 36.0;
    static constexpr double MinPtTrack = 0.1;
    static constexpr double DcaxyMaxTrackPar0 = 0.0105;
    static constexpr double DcaxyMaxTrackPar1 = 0.035;
    static constexpr double DcaxyMaxTrackPar2 = 1.1;
    static constexpr double DcazMaxTrack = 2.0;

    if (!track.hasITS())
      return false;
    if ((!hasITSHit(track, 1)) && (!hasITSHit(track, 2)) && (!hasITSHit(track, 3)))
      return false;
    if (!track.hasTPC())
      return false;
    if (track.tpcNClsCrossedRows() < MinTpcCr)
      return false;
    if (track.tpcChi2NCl() > MaxChi2Tpc)
      return false;
    if (track.itsChi2NCl() > MaxChi2Its)
      return false;
    if (std::fabs(track.eta()) > maxEta)
      return false;
    if (track.pt() < MinPtTrack)
      return false;
    if (std::fabs(track.dcaXY()) > (DcaxyMaxTrackPar0 + DcaxyMaxTrackPar1 / std::pow(track.pt(), DcaxyMaxTrackPar2)))
      return false;
    if (std::fabs(track.dcaZ()) > DcazMaxTrack)
      return false;
    return true;
  }

  // Single-track selection for antinuclei
  template <typename AntinucleusTrack>
  bool passedTrackSelection(const AntinucleusTrack& track)
  {
    if (requirePvContributor && !(track.isPVContributor()))
      return false;
    if (!track.hasITS())
      return false;
    if ((!hasITSHit(track, 1)) && (!hasITSHit(track, 2)) && (!hasITSHit(track, 3)))
      return false;
    if (track.itsNCls() < minItsNclusters)
      return false;
    if (!track.hasTPC())
      return false;
    if (track.tpcNClsCrossedRows() < minTpcNcrossedRows)
      return false;
    if (track.tpcChi2NCl() < minChiSquareTpc)
      return false;
    if (track.tpcChi2NCl() > maxChiSquareTpc)
      return false;
    if (track.itsChi2NCl() > maxChiSquareIts)
      return false;
    if (track.eta() < minEta || track.eta() > maxEta)
      return false;
    if (track.pt() < minPt)
      return false;

    return true;
  }

  // Single-track selection for antinuclei candidates with systematic variations
  template <typename AntinucleusTrack>
  bool passedTrackSelectionSyst(const AntinucleusTrack& track, int isyst)
  {
    // Define cut settings
    static std::vector<int> minItsNclustersSyst = {
      3, 7, 6, 6, 6, 4, 5, 6, 7, 4,
      4, 3, 6, 3, 7, 5, 4, 6, 5, 7,
      6, 5, 3, 5, 4, 3, 6, 6, 4, 7,
      3, 4, 3, 5, 7, 6, 6, 4, 3, 5,
      4, 7, 3, 6, 4, 5, 6, 3, 7, 5};

    static std::vector<int> minTpcNcrossedRowsSyst = {
      90, 108, 112, 119, 92, 111, 98, 105, 86, 117,
      118, 101, 87, 116, 82, 109, 80, 115, 89, 97,
      107, 120, 104, 94, 100, 93, 103, 84, 102, 85,
      108, 96, 113, 117, 91, 88, 99, 110, 106, 83,
      118, 95, 112, 114, 109, 89, 116, 92, 98, 120};

    static std::vector<double> maxChiSquareTpcSyst = {
      4.28, 4.81, 4.43, 4.02, 3.38, 3.58, 3.11, 4.17, 3.51, 4.53,
      4.90, 3.07, 3.20, 4.86, 4.62, 3.91, 3.98, 4.38, 3.66, 3.84,
      3.03, 3.14, 4.96, 4.07, 4.75, 4.32, 3.31, 3.78, 4.11, 3.23,
      3.87, 3.70, 4.99, 4.48, 4.69, 4.25, 3.93, 3.45, 4.58, 3.35,
      3.18, 3.60, 4.21, 3.75, 4.64, 4.35, 3.26, 3.42, 4.15, 3.09};

    static std::vector<double> maxChiSquareItsSyst = {
      42.84, 48.66, 39.27, 34.09, 43.73, 36.98, 30.23, 49.11, 37.67, 35.10,
      44.55, 46.79, 38.92, 40.66, 47.14, 33.46, 30.88, 41.32, 45.90, 39.68,
      31.42, 32.71, 43.17, 36.04, 49.80, 33.95, 31.89, 38.37, 48.08, 35.87,
      47.61, 44.02, 32.15, 46.21, 34.75, 40.17, 37.14, 30.55, 45.42, 42.30,
      41.79, 33.21, 39.12, 47.98, 36.52, 31.58, 49.44, 38.02, 35.56, 43.49};

    // Track Selection
    if (requirePvContributor && !(track.isPVContributor()))
      return false;
    if (!track.hasITS())
      return false;
    if (track.itsNCls() < minItsNclustersSyst[isyst])
      return false;
    if (!track.hasTPC())
      return false;
    if (track.tpcNClsCrossedRows() < minTpcNcrossedRowsSyst[isyst])
      return false;
    if (track.tpcChi2NCl() > maxChiSquareTpcSyst[isyst])
      return false;
    if (track.itsChi2NCl() > maxChiSquareItsSyst[isyst])
      return false;
    if (track.eta() < minEta || track.eta() > maxEta)
      return false;
    if (track.pt() < minPt)
      return false;

    return true;
  }

  // Selection of high-purity antiproton sample
  template <typename AntiprotonTrack>
  bool isHighPurityAntiproton(const AntiprotonTrack& track)
  {
    // variables
    double nsigmaTPCPr = track.tpcNSigmaPr();
    double nsigmaTOFPr = track.tofNSigmaPr();
    double pt = track.pt();
    static constexpr double PtThreshold = 0.5;
    static constexpr double NsigmaMaxPr = 2.0;

    if (pt < PtThreshold && std::fabs(nsigmaTPCPr) < NsigmaMaxPr)
      return true;
    if (pt >= PtThreshold && std::fabs(nsigmaTPCPr) < NsigmaMaxPr && track.hasTOF() && std::fabs(nsigmaTOFPr) < NsigmaMaxPr)
      return true;
    return false;
  }

  // Selection of (anti)protons
  template <typename ProtonTrack>
  bool isProton(const ProtonTrack& track)
  {
    // Constants
    static constexpr double PtThreshold = 0.6;
    static constexpr double NsigmaMax = 3.0;

    // PID variables and transverse momentum of the track
    const double nsigmaTPC = track.tpcNSigmaPr();
    const double nsigmaTOF = track.tofNSigmaPr();
    const double pt = track.pt();

    // Apply TPC PID cut
    if (std::abs(nsigmaTPC) > NsigmaMax)
      return false;

    // Low-pt: TPC PID is sufficient
    if (pt < PtThreshold)
      return true;

    // High-pt: require valid TOF match and pass TOF PID
    return (track.hasTOF() && std::abs(nsigmaTOF) < NsigmaMax);
  }

  // Selection of (anti)deuterons
  template <typename DeuteronTrack>
  bool isDeuteron(const DeuteronTrack& track)
  {
    // Constants
    static constexpr double PtThreshold = 1.0;
    static constexpr double NsigmaMax = 3.0;

    // PID variables and transverse momentum of the track
    const double nsigmaTPC = track.tpcNSigmaDe();
    const double nsigmaTOF = track.tofNSigmaDe();
    const double pt = track.pt();

    // Apply TPC PID cut
    if (std::abs(nsigmaTPC) > NsigmaMax)
      return false;

    // Low-pt: TPC PID is sufficient
    if (pt < PtThreshold)
      return true;

    // High-pt: require valid TOF match and pass TOF PID
    return (track.hasTOF() && std::abs(nsigmaTOF) < NsigmaMax);
  }
    
    // Selection of (anti)protons with systematic variations
      template <typename ProtonTrack>
      bool isProtonSyst(const ProtonTrack& track, double minSigTPC, double maxSigTPC, double minSigTOF, double maxSigTOF)
      {
        // Constant
        static constexpr double kPtThreshold = 0.6;
          
        // PID variables and transverse momentum of the track
        const double nsigmaTPC = track.tpcNSigmaPr();
        const double nsigmaTOF = track.tofNSigmaPr();
        const double pt = track.pt();

         // Apply TPC PID cut (with systematic variations)
        if (nsigmaTPC < minSigTPC || nsigmaTPC > maxSigTPC)
          return false;

        // Low-pt: TPC PID is sufficient
        if (pt < kPtThreshold)
          return true;

        // High-pt: require valid TOF match and pass TOF PID (with systematic variations)
        return (track.hasTOF() && nsigmaTOF > minSigTOF && nsigmaTOF < maxSigTOF);
      }

      // Selection of (anti)deuterons with systematic variations
      template <typename DeuteronTrack>
      bool isDeuteronSyst(const DeuteronTrack& track, double minSigTPC, double maxSigTPC, double minSigTOF, double maxSigTOF)
      {
        // Constant
        static constexpr double kPtThreshold = 1.0;
          
        // PID variables and transverse momentum of the track
        const double nsigmaTPC = track.tpcNSigmaDe();
        const double nsigmaTOF = track.tofNSigmaDe();
        const double pt = track.pt();

        // Apply TPC PID cut (with systematic variations)
        if (nsigmaTPC < minSigTPC || nsigmaTPC > maxSigTPC)
          return false;

        // Low-pt: TPC PID is sufficient
        if (pt < kPtThreshold)
          return true;

        // High-pt: require valid TOF match and pass TOF PID (with systematic variations)
        return (track.hasTOF() && nsigmaTOF > minSigTOF && nsigmaTOF < maxSigTOF);
      }

  // Process Data
  void processData(SelectedCollisions::iterator const& collision, AntiNucleiTracks const& tracks, aod::BCsWithTimestamps const&)
  {
    // Event counter: before event selection
    registryData.fill(HIST("number_of_events_data"), 0.5);
    registryData.fill(HIST("settingData"), minJetPt.value, rJet.value);

    // Retrieve the bunch crossing information with timestamps from the collision
    auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);

    // If skimmed processing is enabled, apply Zorro trigger selection
    if (cfgSkimmedProcessing && !zorro.isSelected(collision.template bc_as<aod::BCsWithTimestamps>().globalBC())) {
      return;
    }
    registryData.fill(HIST("number_of_events_data"), 1.5);

    // Apply standard event selection
    if (!collision.sel8() || std::fabs(collision.posZ()) > zVtx)
      return;

    // Event counter: after event selection
    registryData.fill(HIST("number_of_events_data"), 2.5);

    // Reject events near the ITS Read-Out Frame border
    if (rejectITSROFBorder && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder))
      return;
    registryData.fill(HIST("number_of_events_data"), 3.5);

    // Reject events at the Time Frame border
    if (rejectTFBorder && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder))
      return;
    registryData.fill(HIST("number_of_events_data"), 4.5);

    // Require at least one ITS-TPC matched track
    if (requireVtxITSTPC && !collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC))
      return;
    registryData.fill(HIST("number_of_events_data"), 5.5);

    // Reject events with same-bunch pileup
    if (rejectSameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup))
      return;
    registryData.fill(HIST("number_of_events_data"), 6.5);

    // Require consistent FT0 vs PV z-vertex
    if (requireIsGoodZvtxFT0VsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV))
      return;
    registryData.fill(HIST("number_of_events_data"), 7.5);

    // Require TOF match for at least one vertex track
    if (requireIsVertexTOFmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTOFmatched))
      return;
    registryData.fill(HIST("number_of_events_data"), 8.5);

    // Loop over reconstructed tracks
    int id(-1);
    std::vector<fastjet::PseudoJet> fjParticles;
    for (auto const& track : tracks) {
      id++;
      if (!passedTrackSelectionForJetReconstruction(track))
        continue;

      // 4-momentum representation of a particle
      fastjet::PseudoJet fourMomentum(track.px(), track.py(), track.pz(), track.energy(MassPionCharged));
      fourMomentum.set_user_index(id);
      fjParticles.emplace_back(fourMomentum);

      // Fill TPC signal vs p*sign for PID calibration
      bool heliumPID = track.pidForTracking() == o2::track::PID::Helium3;
      float correctedTpcInnerParam = (heliumPID && cfgCompensatePIDinTracking) ? track.tpcInnerParam() / 2 : track.tpcInnerParam();
      registryData.fill(HIST("tpcsignal_data"), correctedTpcInnerParam * track.sign(), track.tpcSignal());
    }

    // Reject empty events
    if (fjParticles.empty())
      return;
    registryData.fill(HIST("number_of_events_data"), 9.5);

    // Cluster particles using the anti-kt algorithm
    fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, rJet);
    fastjet::AreaDefinition areaDef(fastjet::active_area, fastjet::GhostedAreaSpec(1.0));
    fastjet::ClusterSequenceArea cs(fjParticles, jetDef, areaDef);
    std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cs.inclusive_jets());
    auto [rhoPerp, rhoMPerp] = jetutilities::estimateRhoPerpCone(fjParticles, jets[0], rJet);

    // Loop over reconstructed jets
    bool isAtLeastOneJetSelected = false;
    for (const auto& jet : jets) {

      // Jet must be fully contained in the acceptance
      if ((std::fabs(jet.eta()) + rJet) > (maxEta - deltaEtaEdge))
        continue;

      // Jet pt must be larger than threshold
      auto jetForSub = jet;
      fastjet::PseudoJet jetMinusBkg = backgroundSub.doRhoAreaSub(jetForSub, rhoPerp, rhoMPerp);
      if (jetMinusBkg.pt() < minJetPt)
        continue;

      // Apply area cut if required
      double normalizedJetArea = jet.area() / (PI * rJet * rJet);
      if (applyAreaCut && normalizedJetArea > maxNormalizedJetArea)
        continue;
      isAtLeastOneJetSelected = true;

      // Perpendicular cones
      double coneRadius = std::sqrt(jet.area() / PI);
      TVector3 jetAxis(jet.px(), jet.py(), jet.pz());
      TVector3 ueAxis1(0, 0, 0), ueAxis2(0, 0, 0);
      getPerpendicularDirections(jetAxis, ueAxis1, ueAxis2);
      if (ueAxis1.Mag() == 0 || ueAxis2.Mag() == 0) {
        continue;
      }

      // Fill histogram with jet effective area / piR^2
      registryData.fill(HIST("jetEffectiveAreaOverPiR2"), jet.area() / (PI * rJet * rJet));

      // Get jet constituents
      std::vector<fastjet::PseudoJet> jetConstituents = jet.constituents();

      // Loop over jet constituents
      for (const auto& particle : jetConstituents) {

        // Get corresponding track and apply track selection criteria
        auto const& track = tracks.iteratorAt(particle.user_index());
        if (!passedTrackSelection(track))
          continue;

        // Define variables
        double nsigmaTPCPr = track.tpcNSigmaPr();
        double nsigmaTOFPr = track.tofNSigmaPr();
        double nsigmaTPCDe = track.tpcNSigmaDe();
        double nsigmaTOFDe = track.tofNSigmaDe();
        double nsigmaTPCHe = track.tpcNSigmaHe();
        double pt = track.pt();
        double dcaxy = track.dcaXY();
        double dcaz = track.dcaZ();

        // Fill DCA distribution for antiprotons
        if (track.sign() < 0 && isHighPurityAntiproton(track) && std::fabs(dcaz) < maxDcaz) {
          registryData.fill(HIST("antiproton_dca_jet"), pt, dcaxy);
        }

        // Apply DCA selections
        if (std::fabs(dcaxy) > maxDcaxy || std::fabs(dcaz) > maxDcaz)
          continue;

        // Fill angular distribution of tracks wrt jet axis
        TVector3 trackDirection(track.px(), track.py(), track.pz());
        double thetaTrackJet = (180.0 / PI) * jetAxis.Angle(trackDirection);
        registryData.fill(HIST("theta_track_jet"), jet.pt(), thetaTrackJet);

        // Particle identification using the ITS cluster size
        bool passedItsPidProt(true), passedItsPidDeut(true), passedItsPidHel(true);
        double nSigmaITSprot = static_cast<double>(itsResponse.nSigmaITS<o2::track::PID::Proton>(track));
        double nSigmaITSdeut = static_cast<double>(itsResponse.nSigmaITS<o2::track::PID::Deuteron>(track));
        double nSigmaITShel3 = static_cast<double>(itsResponse.nSigmaITS<o2::track::PID::Helium3>(track));

        // Fill nsigmaITS for antiproton candidates
        if (isHighPurityAntiproton(track)) {
          registryData.fill(HIST("antiproton_nsigma_its_data"), pt, nSigmaITSprot);
        }

        if (applyItsPid && pt < ptMaxItsPidProt && (nSigmaITSprot < nSigmaItsMin || nSigmaITSprot > nSigmaItsMax)) {
          passedItsPidProt = false;
        }
        if (applyItsPid && pt < ptMaxItsPidDeut && (nSigmaITSdeut < nSigmaItsMin || nSigmaITSdeut > nSigmaItsMax)) {
          passedItsPidDeut = false;
        }
        if (applyItsPid && (2.0 * pt) < ptMaxItsPidHel && (nSigmaITShel3 < nSigmaItsMin || nSigmaITShel3 > nSigmaItsMax)) {
          passedItsPidHel = false;
        }

        // Fill histograms for antimatter
        if (track.sign() < 0) {
          if (passedItsPidProt) {
            registryData.fill(HIST("antiproton_jet_tpc"), pt, nsigmaTPCPr);
            if (nsigmaTPCPr > minNsigmaTpc && nsigmaTPCPr < maxNsigmaTpc && track.hasTOF())
              registryData.fill(HIST("antiproton_jet_tof"), pt, nsigmaTOFPr);
          }
          if (passedItsPidDeut) {
            registryData.fill(HIST("antideuteron_jet_tpc"), pt, nsigmaTPCDe);
            if (nsigmaTPCDe > minNsigmaTpc && nsigmaTPCDe < maxNsigmaTpc && track.hasTOF())
              registryData.fill(HIST("antideuteron_jet_tof"), pt, nsigmaTOFDe);
          }
          if (passedItsPidHel) {
            registryData.fill(HIST("antihelium3_jet_tpc"), 2.0 * pt, nsigmaTPCHe);
            // custom nsigma He3 based on bethe bloch fit of TPC signal
            double tpcSignal = track.tpcSignal();
            double expectedSignalHe3 = tpc::BetheBlochAleph(static_cast<double>(track.tpcInnerParam() * 2. / o2::constants::physics::MassHelium3), cfgBetheBlochParams.value[0], cfgBetheBlochParams.value[1], cfgBetheBlochParams.value[2], cfgBetheBlochParams.value[3], cfgBetheBlochParams.value[4]);
            double nSigmaTPCHe3Custom = ((tpcSignal / expectedSignalHe3) - 1.) / 0.045;
            registryData.fill(HIST("antihelium3_jet_tpc_custom"), 2.0 * pt, nSigmaTPCHe3Custom);
          }
        }

        // Fill histograms for matter
        if (track.sign() > 0) {
          if (passedItsPidDeut) {
            registryData.fill(HIST("deuteron_jet_tpc"), pt, nsigmaTPCDe);
            if (nsigmaTPCDe > minNsigmaTpc && nsigmaTPCDe < maxNsigmaTpc && track.hasTOF())
              registryData.fill(HIST("deuteron_jet_tof"), pt, nsigmaTOFDe);
          }
          if (passedItsPidHel) {
            registryData.fill(HIST("helium3_jet_tpc"), 2.0 * pt, nsigmaTPCHe);
          }
        }
      }

      // Loop over tracks in the underlying event
      for (auto const& track : tracks) {

        // Get corresponding track and apply track selection criteria
        if (!passedTrackSelection(track))
          continue;

        // Calculate the angular distance between the track and underlying event axes in eta-phi space
        double deltaEtaUe1 = track.eta() - ueAxis1.Eta();
        double deltaPhiUe1 = getDeltaPhi(track.phi(), ueAxis1.Phi());
        double deltaRUe1 = std::sqrt(deltaEtaUe1 * deltaEtaUe1 + deltaPhiUe1 * deltaPhiUe1);
        double deltaEtaUe2 = track.eta() - ueAxis2.Eta();
        double deltaPhiUe2 = getDeltaPhi(track.phi(), ueAxis2.Phi());
        double deltaRUe2 = std::sqrt(deltaEtaUe2 * deltaEtaUe2 + deltaPhiUe2 * deltaPhiUe2);

        // Determine the maximum allowed distance from UE axes for particle selection
        double maxConeRadius = coneRadius;
        if (applyAreaCut) {
          maxConeRadius = std::sqrt(maxNormalizedJetArea) * rJet;
        }

        // Reject tracks that lie outside the maxConeRadius from both UE axes
        if (deltaRUe1 > maxConeRadius && deltaRUe2 > maxConeRadius)
          continue;

        // Define variables
        double nsigmaTPCPr = track.tpcNSigmaPr();
        double nsigmaTOFPr = track.tofNSigmaPr();
        double nsigmaTPCDe = track.tpcNSigmaDe();
        double nsigmaTOFDe = track.tofNSigmaDe();
        double nsigmaTPCHe = track.tpcNSigmaHe();
        double pt = track.pt();
        double dcaxy = track.dcaXY();
        double dcaz = track.dcaZ();

        // Fill DCA distribution for antiprotons
        if (track.sign() < 0 && isHighPurityAntiproton(track) && std::fabs(dcaz) < maxDcaz) {
          registryData.fill(HIST("antiproton_dca_ue"), pt, dcaxy);
        }

        // Apply DCA selections
        if (std::fabs(dcaxy) > maxDcaxy || std::fabs(dcaz) > maxDcaz)
          continue;

        // Particle identification using the ITS cluster size
        bool passedItsPidProt(true), passedItsPidDeut(true), passedItsPidHel(true);
        double nSigmaITSprot = static_cast<double>(itsResponse.nSigmaITS<o2::track::PID::Proton>(track));
        double nSigmaITSdeut = static_cast<double>(itsResponse.nSigmaITS<o2::track::PID::Deuteron>(track));
        double nSigmaITShel3 = static_cast<double>(itsResponse.nSigmaITS<o2::track::PID::Helium3>(track));

        if (applyItsPid && pt < ptMaxItsPidProt && (nSigmaITSprot < nSigmaItsMin || nSigmaITSprot > nSigmaItsMax)) {
          passedItsPidProt = false;
        }
        if (applyItsPid && pt < ptMaxItsPidDeut && (nSigmaITSdeut < nSigmaItsMin || nSigmaITSdeut > nSigmaItsMax)) {
          passedItsPidDeut = false;
        }
        if (applyItsPid && (2.0 * pt) < ptMaxItsPidHel && (nSigmaITShel3 < nSigmaItsMin || nSigmaITShel3 > nSigmaItsMax)) {
          passedItsPidHel = false;
        }

        // Fill histograms for antimatter
        if (track.sign() < 0) {
          if (passedItsPidProt) {
            registryData.fill(HIST("antiproton_ue_tpc"), pt, nsigmaTPCPr);
            if (nsigmaTPCPr > minNsigmaTpc && nsigmaTPCPr < maxNsigmaTpc && track.hasTOF())
              registryData.fill(HIST("antiproton_ue_tof"), pt, nsigmaTOFPr);
          }
          if (passedItsPidDeut) {
            registryData.fill(HIST("antideuteron_ue_tpc"), pt, nsigmaTPCDe);
            if (nsigmaTPCDe > minNsigmaTpc && nsigmaTPCDe < maxNsigmaTpc && track.hasTOF())
              registryData.fill(HIST("antideuteron_ue_tof"), pt, nsigmaTOFDe);
          }
          if (passedItsPidHel) {
            registryData.fill(HIST("antihelium3_ue_tpc"), 2.0 * pt, nsigmaTPCHe);
            // custom nsigma He3 based on bethe bloch fit of TPC signal
            double tpcSignal = track.tpcSignal();
            double expectedSignalHe3 = tpc::BetheBlochAleph(static_cast<double>(track.tpcInnerParam() * 2. / o2::constants::physics::MassHelium3), cfgBetheBlochParams.value[0], cfgBetheBlochParams.value[1], cfgBetheBlochParams.value[2], cfgBetheBlochParams.value[3], cfgBetheBlochParams.value[4]);
            double nSigmaTPCHe3Custom = ((tpcSignal / expectedSignalHe3) - 1.) / 0.045;
            registryData.fill(HIST("antihelium3_ue_tpc_custom"), 2.0 * pt, nSigmaTPCHe3Custom);
          }
        }

        //  Fill histograms for matter
        if (track.sign() > 0) {
          if (passedItsPidDeut) {
            registryData.fill(HIST("deuteron_ue_tpc"), pt, nsigmaTPCDe);
            if (nsigmaTPCDe > minNsigmaTpc && nsigmaTPCDe < maxNsigmaTpc && track.hasTOF())
              registryData.fill(HIST("deuteron_ue_tof"), pt, nsigmaTOFDe);
          }
          if (passedItsPidHel) {
            registryData.fill(HIST("helium3_ue_tpc"), 2.0 * pt, nsigmaTPCHe);
          }
        }
      }
    }
    // Event counter: events with at least one jet selected
    if (isAtLeastOneJetSelected) {
      registryData.fill(HIST("number_of_events_data"), 10.5);
    }
  }
  PROCESS_SWITCH(AntinucleiInJets, processData, "Process Data", true);

  // Charged-particle multiplicity in events with selected jets (Run 3) or ptTrigger > threshold (Run 2-like)
  void processMultEvents(SelectedCollisions::iterator const& collision, AntiNucleiTracks const& tracks)
  {
    // Apply event selection
    if (!collision.sel8() || std::fabs(collision.posZ()) > zVtx)
      return;
    if (rejectITSROFBorder && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder))
      return;
    if (rejectTFBorder && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder))
      return;
    if (requireVtxITSTPC && !collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC))
      return;
    if (rejectSameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup))
      return;
    if (requireIsGoodZvtxFT0VsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV))
      return;
    if (requireIsVertexTOFmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTOFmatched))
      return;

    // Initialize variable to store the maximum pt in the event
    double ptMax(0.0);

    // Loop over reconstructed tracks
    int id(-1);
    std::vector<fastjet::PseudoJet> fjParticles;
    for (auto const& track : tracks) {
      id++;
      if (!passedTrackSelectionForJetReconstruction(track))
        continue;
      if (track.pt() > ptMax) {
        ptMax = track.pt();
      }

      // 4-momentum representation of a particle
      fastjet::PseudoJet fourMomentum(track.px(), track.py(), track.pz(), track.energy(MassPionCharged));
      fourMomentum.set_user_index(id);
      fjParticles.emplace_back(fourMomentum);
    }

    // Reject empty events
    if (fjParticles.empty()) {
      return;
    }

    // Fill charged-particle multiplicity for events with leading track having pt>threshold
    if (ptMax > ptLeadingMin) {
      registryMult.fill(HIST("multiplicityEvtsPtLeading"), fjParticles.size());
    }

    // Cluster particles using the anti-kt algorithm
    fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, rJet);
    fastjet::AreaDefinition areaDef(fastjet::active_area, fastjet::GhostedAreaSpec(1.0));
    fastjet::ClusterSequenceArea cs(fjParticles, jetDef, areaDef);
    std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cs.inclusive_jets());
    auto [rhoPerp, rhoMPerp] = jetutilities::estimateRhoPerpCone(fjParticles, jets[0], rJet);

    // Loop over reconstructed jets
    bool isAtLeastOneJetSelected = false;
    for (const auto& jet : jets) {

      // Jet must be fully contained in the acceptance
      if ((std::fabs(jet.eta()) + rJet) > (maxEta - deltaEtaEdge))
        continue;

      // Jet pt must be larger than threshold
      auto jetForSub = jet;
      fastjet::PseudoJet jetMinusBkg = backgroundSub.doRhoAreaSub(jetForSub, rhoPerp, rhoMPerp);
      if (jetMinusBkg.pt() < minJetPt)
        continue;

      // Apply area cut if required
      double normalizedJetArea = jet.area() / (PI * rJet * rJet);
      if (applyAreaCut && normalizedJetArea > maxNormalizedJetArea)
        continue;
      isAtLeastOneJetSelected = true;
    }

    // Fill histogram of charged-particle multiplicity for events containing at least one selected jet
    if (isAtLeastOneJetSelected) {
      registryMult.fill(HIST("multiplicityEvtsWithJet"), fjParticles.size());
    }
  }
  PROCESS_SWITCH(AntinucleiInJets, processMultEvents, "Process Mult Events", false);

  // Process QC
  void processQC(SelectedCollisions::iterator const& collision, AntiNucleiTracks const& tracks)
  {
    // Apply event selection
    if (!collision.sel8() || std::fabs(collision.posZ()) > zVtx)
      return;
    if (rejectITSROFBorder && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder))
      return;
    if (rejectTFBorder && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder))
      return;
    if (requireVtxITSTPC && !collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC))
      return;
    if (rejectSameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup))
      return;
    if (requireIsGoodZvtxFT0VsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV))
      return;
    if (requireIsVertexTOFmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTOFmatched))
      return;

    // Loop over reconstructed tracks
    std::vector<fastjet::PseudoJet> fjParticles;
    for (auto const& track : tracks) {
      if (!passedTrackSelectionForJetReconstruction(track))
        continue;

      // 4-momentum representation of a particle
      fastjet::PseudoJet fourMomentum(track.px(), track.py(), track.pz(), track.energy(MassPionCharged));
      fjParticles.emplace_back(fourMomentum);
    }

    // Reject empty events
    if (fjParticles.empty())
      return;

    // Cluster particles using the anti-kt algorithm
    fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, rJet);
    fastjet::AreaDefinition areaDef(fastjet::active_area, fastjet::GhostedAreaSpec(1.0));
    fastjet::ClusterSequenceArea cs(fjParticles, jetDef, areaDef);
    std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cs.inclusive_jets());
    auto [rhoPerp, rhoMPerp] = jetutilities::estimateRhoPerpCone(fjParticles, jets[0], rJet);

    // Loop over reconstructed jets
    int njetsInAcc(0);
    int njetsHighPt(0);
    for (const auto& jet : jets) {

      // Jet must be fully contained in the acceptance
      if ((std::fabs(jet.eta()) + rJet) > (maxEta - deltaEtaEdge))
        continue;
      njetsInAcc++;
      registryQC.fill(HIST("sumPtJetCone"), jet.pt());
      double ptJetBeforeSub = jet.pt();

      // Jet pt must be larger than threshold
      auto jetForSub = jet;
      fastjet::PseudoJet jetMinusBkg = backgroundSub.doRhoAreaSub(jetForSub, rhoPerp, rhoMPerp);
      double ptJetAfterSub = jetForSub.pt();
      registryQC.fill(HIST("jetPtDifference"), ptJetAfterSub - ptJetBeforeSub);
      registryQC.fill(HIST("ptDistributionJetCone"), ptJetBeforeSub);
      registryQC.fill(HIST("ptDistributionJet"), ptJetAfterSub);

      if (jetMinusBkg.pt() < minJetPt)
        continue;
      njetsHighPt++;
      registryQC.fill(HIST("sumPtJet"), jet.pt());

      // Jet properties and perpendicular cone
      std::vector<fastjet::PseudoJet> jetConstituents = jet.constituents();
      TVector3 jetAxis(jet.px(), jet.py(), jet.pz());
      double coneRadius = std::sqrt(jet.area() / PI);
      TVector3 ueAxis1(0, 0, 0), ueAxis2(0, 0, 0);
      getPerpendicularDirections(jetAxis, ueAxis1, ueAxis2);
      if (ueAxis1.Mag() == 0 || ueAxis2.Mag() == 0) {
        continue;
      }

      registryQC.fill(HIST("NchJetCone"), static_cast<int>(jetConstituents.size()));

      // Loop over jet constituents
      for (const auto& particle : jetConstituents) {

        double deltaEta = particle.eta() - jetAxis.Eta();
        double deltaPhi = getDeltaPhi(particle.phi(), jetAxis.Phi());
        registryQC.fill(HIST("deltaEta_deltaPhi_jet"), deltaEta, deltaPhi);
        registryQC.fill(HIST("eta_phi_jet"), particle.eta(), particle.phi());
      }

      // Loop over particles in perpendicular cones
      double nParticlesPerp(0);
      double ptPerp(0);
      for (auto const& track : tracks) {

        if (!passedTrackSelectionForJetReconstruction(track))
          continue;

        double deltaEtaUe1 = track.eta() - ueAxis1.Eta();
        double deltaPhiUe1 = getDeltaPhi(track.phi(), ueAxis1.Phi());
        double deltaRUe1 = std::sqrt(deltaEtaUe1 * deltaEtaUe1 + deltaPhiUe1 * deltaPhiUe1);
        double deltaEtaUe2 = track.eta() - ueAxis2.Eta();
        double deltaPhiUe2 = getDeltaPhi(track.phi(), ueAxis2.Phi());
        double deltaRUe2 = std::sqrt(deltaEtaUe2 * deltaEtaUe2 + deltaPhiUe2 * deltaPhiUe2);
        if (deltaRUe1 > coneRadius && deltaRUe2 > coneRadius)
          continue;

        ptPerp = ptPerp + track.pt();
        nParticlesPerp++;
        registryQC.fill(HIST("deltaEta_deltaPhi_ue"), deltaEtaUe1, deltaPhiUe1);
        registryQC.fill(HIST("deltaEta_deltaPhi_ue"), deltaEtaUe2, deltaPhiUe2);
        registryQC.fill(HIST("eta_phi_ue"), track.eta(), track.phi());
      }
      registryQC.fill(HIST("NchUE"), 0.5 * nParticlesPerp);
      registryQC.fill(HIST("NchJet"), static_cast<double>(jetConstituents.size()) - 0.5 * nParticlesPerp);
      registryQC.fill(HIST("sumPtUE"), 0.5 * ptPerp);
    }
    registryQC.fill(HIST("nJetsFound"), static_cast<int>(jets.size()));
    registryQC.fill(HIST("nJetsInAcceptance"), njetsInAcc);
    registryQC.fill(HIST("nJetsSelectedHighPt"), njetsHighPt);
  }
  PROCESS_SWITCH(AntinucleiInJets, processQC, "Process QC", false);

  // Define preslices to group MC tracks and MC particles by their associated MC collision
  Preslice<AntiNucleiTracksMc> mcTracksPerMcCollision = o2::aod::track::collisionId;
  Preslice<aod::McParticles> mcParticlesPerMcCollision = o2::aod::mcparticle::mcCollisionId;

  // Antinuclei reconstruction efficiency
  void processAntinucleiEfficiency(GenCollisionsMc const& genCollisions, RecCollisionsMc const& recCollisions, AntiNucleiTracksMc const& mcTracks, aod::McParticles const& mcParticles)
  {
    // Loop over generated collisions
    for (const auto& collision : genCollisions) {

      // Apply event selection: require vertex position to be within the allowed z range
      if (std::fabs(collision.posZ()) > zVtx)
        continue;

      // Get particles in this MC collision
      const auto mcParticlesThisMcColl = mcParticles.sliceBy(mcParticlesPerMcCollision, collision.globalIndex());

      // Loop over all generated Monte Carlo particles for the selected event
      for (const auto& particle : mcParticlesThisMcColl) {

        // Select primary particles
        if (!particle.isPhysicalPrimary())
          continue;

        // Select particles within the specified pseudorapidity interval
        if (particle.eta() < minEta || particle.eta() > maxEta)
          continue;

        // Process different particle species based on PDG code
        switch (particle.pdgCode()) {
          case PDG_t::kProtonBar:
            registryMC.fill(HIST("antip_gen_jet"), particle.pt());
            registryMC.fill(HIST("antip_gen_ue"), particle.pt());
            registryMC.fill(HIST("protonBar"), particle.pt());
            break;
          case o2::constants::physics::Pdg::kDeuteron:
            registryMC.fill(HIST("deuteron_gen_jet"), particle.pt());
            registryMC.fill(HIST("deuteron_gen_ue"), particle.pt());
            break;
          case -o2::constants::physics::Pdg::kDeuteron:
            registryMC.fill(HIST("antideuteron_gen_jet"), particle.pt());
            registryMC.fill(HIST("antideuteron_gen_ue"), particle.pt());
            break;
          case o2::constants::physics::Pdg::kHelium3:
            registryMC.fill(HIST("helium3_gen_jet"), particle.pt());
            registryMC.fill(HIST("helium3_gen_ue"), particle.pt());
            break;
          case -o2::constants::physics::Pdg::kHelium3:
            registryMC.fill(HIST("antihelium3_gen_jet"), particle.pt());
            registryMC.fill(HIST("antihelium3_gen_ue"), particle.pt());
            break;
          // Histograms for re-weighting
          case PDG_t::kLambda0Bar:
            registryMC.fill(HIST("lambdaBar"), particle.pt());
            break;
          case PDG_t::kXiPlusBar:
            registryMC.fill(HIST("xiBar"), particle.pt());
            break;
          case PDG_t::kOmegaPlusBar:
            registryMC.fill(HIST("omegaBar"), particle.pt());
            break;
          case PDG_t::kSigmaBarMinus:
            registryMC.fill(HIST("sigmaBar"), particle.pt());
            break;
        }
      }
    }

    // Loop over all reconstructed collisions
    for (const auto& collision : recCollisions) {

      // Count all generated events before applying any event selection criteria
      registryMC.fill(HIST("number_of_events_mc_nuclei_efficiency"), 0.5);

      // Apply event selection: require sel8 and vertex position within the allowed z range
      if (!collision.sel8() || std::fabs(collision.posZ()) > zVtx)
        continue;

      // Count events that pass the selection criteria
      registryMC.fill(HIST("number_of_events_mc_nuclei_efficiency"), 1.5);

      // Reject events near the ITS Read-Out Frame border
      if (rejectITSROFBorder && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder))
        continue;
      registryMC.fill(HIST("number_of_events_mc_nuclei_efficiency"), 2.5);

      // Reject events at the Time Frame border
      if (rejectTFBorder && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder))
        continue;
      registryMC.fill(HIST("number_of_events_mc_nuclei_efficiency"), 3.5);

      // Require at least one ITS-TPC matched track
      if (requireVtxITSTPC && !collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC))
        continue;
      registryMC.fill(HIST("number_of_events_mc_nuclei_efficiency"), 4.5);

      // Reject events with same-bunch pileup
      if (rejectSameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup))
        continue;
      registryMC.fill(HIST("number_of_events_mc_nuclei_efficiency"), 5.5);

      // Require consistent FT0 vs PV z-vertex
      if (requireIsGoodZvtxFT0VsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV))
        continue;
      registryMC.fill(HIST("number_of_events_mc_nuclei_efficiency"), 6.5);

      // Require TOF match for at least one vertex track
      if (requireIsVertexTOFmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTOFmatched))
        continue;
      registryMC.fill(HIST("number_of_events_mc_nuclei_efficiency"), 7.5);

      // Get tracks in this MC collision
      const auto mcTracksThisMcColl = mcTracks.sliceBy(mcTracksPerMcCollision, collision.globalIndex());

      // Loop over all reconstructed MC tracks
      for (auto const& track : mcTracksThisMcColl) {

        // Apply standard track selection criteria
        if (!passedTrackSelection(track))
          continue;

        // Cut on transverse and longitudinal distance of closest approach
        if (std::fabs(track.dcaXY()) > maxDcaxy || std::fabs(track.dcaZ()) > maxDcaz)
          continue;

        // Skip tracks that are not associated with a true MC particle
        if (!track.has_mcParticle())
          continue;
        const auto particle = track.mcParticle();

        // ****************************************************************************************************************

        // Systematic uncertainty on primary fraction
        if (track.sign() < 0 && particle.pdgCode() == PDG_t::kProtonBar) {

          // Primary antiprotons
          if (particle.isPhysicalPrimary()) {

            // Initialize weights
            double wPrimStd(1.0), wPrimUp(1.0), wPrimLow(1.0);

            // Weight assignment
            if (particle.pt() < primaryAntiprotons->GetXaxis()->GetXmax()) {
              int ipt = primaryAntiprotons->FindBin(particle.pt());
              wPrimStd = primaryAntiprotons->GetBinContent(ipt);
              wPrimUp = wPrimStd + primaryAntiprotons->GetBinError(ipt);
              wPrimLow = wPrimStd - primaryAntiprotons->GetBinError(ipt);
            }

            // Fill histograms
            registryMC.fill(HIST("antip_prim_pythia"), track.pt());
            registryMC.fill(HIST("antip_prim_std"), track.pt(), wPrimStd);
            registryMC.fill(HIST("antip_prim_up"), track.pt(), wPrimUp);
            registryMC.fill(HIST("antip_prim_low"), track.pt(), wPrimLow);
          }

          // Secondary antiprotons from material
          if (!particle.isPhysicalPrimary() && !particle.has_mothers()) {

            // Fill histograms
            registryMC.fill(HIST("antip_sec_pythia"), track.pt());
            registryMC.fill(HIST("antip_sec_std"), track.pt());
            registryMC.fill(HIST("antip_sec_up"), track.pt());
            registryMC.fill(HIST("antip_sec_low"), track.pt());
          }

          // Secondary antiprotons from weak decays
          if (!particle.isPhysicalPrimary() && particle.has_mothers()) {

            // Get first ancestor
            auto ancestor = getFirstAncestor(particle, mcParticles);
            double wSecStd(1.0), wSecUp(1.0), wSecLow(1.0);

            // Antiprotons from antiSigma
            if (ancestor.pdgCode() == PDG_t::kSigmaBarMinus && ancestor.pt() < primaryAntiSigma->GetXaxis()->GetXmax()) {
              int ipt = primaryAntiSigma->FindBin(ancestor.pt());
              wSecStd = primaryAntiSigma->GetBinContent(ipt);
              wSecUp = wSecStd + primaryAntiSigma->GetBinError(ipt);
              wSecLow = wSecStd - primaryAntiSigma->GetBinError(ipt);
            }

            // Antiprotons from antiLambda0
            if (ancestor.pdgCode() == PDG_t::kLambda0Bar && ancestor.pt() < primaryAntiLambda->GetXaxis()->GetXmax()) {
              int ipt = primaryAntiLambda->FindBin(ancestor.pt());
              wSecStd = primaryAntiLambda->GetBinContent(ipt);
              wSecUp = wSecStd + primaryAntiLambda->GetBinError(ipt);
              wSecLow = wSecStd - primaryAntiLambda->GetBinError(ipt);
            }

            // Antiprotons from antiXi
            if (ancestor.pdgCode() == PDG_t::kXiPlusBar && ancestor.pt() < primaryAntiXi->GetXaxis()->GetXmax()) {
              int ipt = primaryAntiXi->FindBin(ancestor.pt());
              wSecStd = primaryAntiXi->GetBinContent(ipt);
              wSecUp = wSecStd + primaryAntiXi->GetBinError(ipt);
              wSecLow = wSecStd - primaryAntiXi->GetBinError(ipt);
            }

            // Antiprotons from antiXi0
            if (ancestor.pdgCode() == -o2::constants::physics::Pdg::kXi0 && ancestor.pt() < primaryAntiXi->GetXaxis()->GetXmax()) {
              int ipt = primaryAntiXi->FindBin(ancestor.pt());
              wSecStd = primaryAntiXi->GetBinContent(ipt);
              wSecUp = wSecStd + primaryAntiXi->GetBinError(ipt);
              wSecLow = wSecStd - primaryAntiXi->GetBinError(ipt);
            }

            // Antiprotons from antiOmega
            if (ancestor.pdgCode() == PDG_t::kOmegaPlusBar && ancestor.pt() < primaryAntiOmega->GetXaxis()->GetXmax()) {
              int ipt = primaryAntiOmega->FindBin(ancestor.pt());
              wSecStd = primaryAntiOmega->GetBinContent(ipt);
              wSecUp = wSecStd + primaryAntiOmega->GetBinError(ipt);
              wSecLow = wSecStd - primaryAntiOmega->GetBinError(ipt);
            }

            // Fill histograms
            registryMC.fill(HIST("antip_sec_pythia"), track.pt());
            registryMC.fill(HIST("antip_sec_std"), track.pt(), wSecStd);
            registryMC.fill(HIST("antip_sec_up"), track.pt(), wSecUp);
            registryMC.fill(HIST("antip_sec_low"), track.pt(), wSecLow);
          }
        }

        // ****************************************************************************************************************

        // Select only physical primary particles
        if (!particle.isPhysicalPrimary())
          continue;

        // Retrieve PID responses from TPC and TOF detectors for proton, deuteron, helium-3
        double nsigmaTPCPr = track.tpcNSigmaPr();
        double nsigmaTOFPr = track.tofNSigmaPr();
        double nsigmaTPCDe = track.tpcNSigmaDe();
        double nsigmaTOFDe = track.tofNSigmaDe();
        double nsigmaTPCHe = track.tpcNSigmaHe();
        double pt = track.pt();

        // particle identification using the ITS cluster size
        bool passedItsPidProt(true), passedItsPidDeut(true), passedItsPidHel(true);
        double nSigmaITSprot = static_cast<double>(itsResponse.nSigmaITS<o2::track::PID::Proton>(track));
        double nSigmaITSdeut = static_cast<double>(itsResponse.nSigmaITS<o2::track::PID::Deuteron>(track));
        double nSigmaITShel3 = static_cast<double>(itsResponse.nSigmaITS<o2::track::PID::Helium3>(track));

        if (applyItsPid && pt < ptMaxItsPidProt && (nSigmaITSprot < nSigmaItsMin || nSigmaITSprot > nSigmaItsMax)) {
          passedItsPidProt = false;
        }
        if (applyItsPid && pt < ptMaxItsPidDeut && (nSigmaITSdeut < nSigmaItsMin || nSigmaITSdeut > nSigmaItsMax)) {
          passedItsPidDeut = false;
        }
        if (applyItsPid && (2.0 * pt) < ptMaxItsPidHel && (nSigmaITShel3 < nSigmaItsMin || nSigmaITShel3 > nSigmaItsMax)) {
          passedItsPidHel = false;
        }

        // Fill nsigmaITS for antiproton candidates
        if (track.sign() < 0 && particle.pdgCode() == PDG_t::kProtonBar && isHighPurityAntiproton(track)) {
          registryMC.fill(HIST("antiproton_nsigma_its_mc"), pt, nSigmaITSprot);
        }

        // Fill histograms of antiprotons
        if (track.sign() < 0 && particle.pdgCode() == PDG_t::kProtonBar && passedItsPidProt) {
          if (nsigmaTPCPr > minNsigmaTpc && nsigmaTPCPr < maxNsigmaTpc) {
            registryMC.fill(HIST("antip_rec_tpc_jet"), track.pt());
            registryMC.fill(HIST("antip_rec_tpc_ue"), track.pt());

            if (track.hasTOF() && nsigmaTOFPr > minNsigmaTof && nsigmaTOFPr < maxNsigmaTof) {
              registryMC.fill(HIST("antip_rec_tof_jet"), track.pt());
              registryMC.fill(HIST("antip_rec_tof_ue"), track.pt());
            }
          }
        }

        // Fill histograms of antideuterons
        if (track.sign() < 0 && particle.pdgCode() == -o2::constants::physics::Pdg::kDeuteron && passedItsPidDeut) {
          if (nsigmaTPCDe > minNsigmaTpc && nsigmaTPCDe < maxNsigmaTpc) {
            registryMC.fill(HIST("antideuteron_rec_tpc_jet"), track.pt());
            registryMC.fill(HIST("antideuteron_rec_tpc_ue"), track.pt());

            if (track.hasTOF() && nsigmaTOFDe > minNsigmaTof && nsigmaTOFDe < maxNsigmaTof) {
              registryMC.fill(HIST("antideuteron_rec_tof_jet"), track.pt());
              registryMC.fill(HIST("antideuteron_rec_tof_ue"), track.pt());
            }
          }
        }

        // Fill histograms of deuterons
        if (track.sign() > 0 && particle.pdgCode() == o2::constants::physics::Pdg::kDeuteron && passedItsPidDeut) {
          if (nsigmaTPCDe > minNsigmaTpc && nsigmaTPCDe < maxNsigmaTpc) {
            registryMC.fill(HIST("deuteron_rec_tpc_jet"), track.pt());
            registryMC.fill(HIST("deuteron_rec_tpc_ue"), track.pt());
            if (track.hasTOF() && nsigmaTOFDe > minNsigmaTof && nsigmaTOFDe < maxNsigmaTof) {
              registryMC.fill(HIST("deuteron_rec_tof_jet"), track.pt());
              registryMC.fill(HIST("deuteron_rec_tof_ue"), track.pt());
            }
          }
        }

        // Fill histograms of antihelium3
        if (track.sign() < 0 && particle.pdgCode() == -o2::constants::physics::Pdg::kHelium3 && passedItsPidHel) {
          if (nsigmaTPCHe > minNsigmaTpc && nsigmaTPCHe < maxNsigmaTpc) {
            registryMC.fill(HIST("antihelium3_rec_tpc_jet"), 2.0 * track.pt());
            registryMC.fill(HIST("antihelium3_rec_tpc_ue"), 2.0 * track.pt());
          }
        }

        // Fill histograms of helium3
        if (track.sign() > 0 && particle.pdgCode() == o2::constants::physics::Pdg::kHelium3 && passedItsPidHel) {
          if (nsigmaTPCHe > minNsigmaTpc && nsigmaTPCHe < maxNsigmaTpc) {
            registryMC.fill(HIST("helium3_rec_tpc_jet"), 2.0 * track.pt());
            registryMC.fill(HIST("helium3_rec_tpc_ue"), 2.0 * track.pt());
          }
        }
      }
    }
  }
  PROCESS_SWITCH(AntinucleiInJets, processAntinucleiEfficiency, "process antinuclei efficiency", false);

  // Generated events
  void processJetsMCgen(GenCollisionsMc const& collisions, aod::McParticles const& mcParticles)
  {
    // Define per-event particle containers
    std::vector<fastjet::PseudoJet> fjParticles;
    std::vector<TVector3> protonMomentum;

    // Event counter
    int eventCounter = 0;

    // Jet and area definitions
    fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, rJet);
    fastjet::AreaDefinition areaDef(fastjet::active_area, fastjet::GhostedAreaSpec(1.0));

    // Loop over all simulated collisions
    for (const auto& collision : collisions) {

      // Increment event counter
      eventCounter++;

      // Clear containers at the start of the event loop
      fjParticles.clear();
      protonMomentum.clear();

      // Event counter: before event selection
      registryMC.fill(HIST("genEvents"), 0.5);

      // Apply event selection: require vertex position to be within the allowed z range
      if (std::fabs(collision.posZ()) > zVtx)
        continue;

      // Event counter: after event selection
      registryMC.fill(HIST("genEvents"), 1.5);

      // Get particles in this MC collision
      const auto mcParticlesThisMcColl = mcParticles.sliceBy(mcParticlesPerMcCollision, collision.globalIndex());

      // Loop over MC particles
      for (const auto& particle : mcParticlesThisMcColl) {

        // Select physical primary particles or HF decay products
        if (!isPhysicalPrimaryOrFromHF(particle, mcParticles))
          continue;

        // Select particles within acceptance
        static constexpr double MinPtParticle = 0.1;
        if (particle.eta() < minEta || particle.eta() > maxEta || particle.pt() < MinPtParticle)
          continue;

        // Store 3-momentum vectors of antiprotons for further analysis
        if (particle.pdgCode() == PDG_t::kProtonBar) {
          TVector3 pVec(particle.px(), particle.py(), particle.pz());
          protonMomentum.emplace_back(pVec);
          registryMC.fill(HIST("antiproton_gen_full"), particle.pt());
        }

        // 4-momentum representation of a particle
        double energy = std::sqrt(particle.p() * particle.p() + MassPionCharged * MassPionCharged);
        fastjet::PseudoJet fourMomentum(particle.px(), particle.py(), particle.pz(), energy);
        fourMomentum.set_user_index(particle.pdgCode());
        fjParticles.emplace_back(fourMomentum);
      }

      // Reject empty events
      if (fjParticles.empty())
        continue;
      registryMC.fill(HIST("genEvents"), 2.5);

      // Size of particle array
      registryMC.fill(HIST("sizeParticleArray"), fjParticles.size());

      // Cluster MC particles into jets using anti-kt algorithm
      fastjet::ClusterSequenceArea cs(fjParticles, jetDef, areaDef);
      std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cs.inclusive_jets());
      auto [rhoPerp, rhoMPerp] = jetutilities::estimateRhoPerpCone(fjParticles, jets[0], rJet);

      // Loop over clustered jets
      bool isAtLeastOneJetSelected = false;
      for (const auto& jet : jets) {

        // Jet must be fully contained in the acceptance
        if ((std::fabs(jet.eta()) + rJet) > (maxEta - deltaEtaEdge))
          continue;

        // Jet pt must be larger than threshold
        auto jetForSub = jet;
        fastjet::PseudoJet jetMinusBkg = backgroundSub.doRhoAreaSub(jetForSub, rhoPerp, rhoMPerp);
        if (jetMinusBkg.pt() < minJetPt)
          continue;

        // Apply area cut if required
        double normalizedJetArea = jet.area() / (PI * rJet * rJet);
        if (applyAreaCut && normalizedJetArea > maxNormalizedJetArea)
          continue;
        isAtLeastOneJetSelected = true;

        // Generated jets
        registryMC.fill(HIST("genJets"), 0.5);

        // Analyze jet constituents
        std::vector<fastjet::PseudoJet> jetConstituents = jet.constituents();
        for (const auto& particle : jetConstituents) {
          if (particle.user_index() != PDG_t::kProtonBar)
            continue;

          if (particle.eta() < minEta || particle.eta() > maxEta)
            continue;

          // Fill normalization histogram
          registryMC.fill(HIST("antiproton_deltay_deltaphi_jet"), particle.eta() - jet.eta(), getDeltaPhi(particle.phi(), jet.phi()));

          // Calculate weight
          double weightJet(1.0);
          if (applyReweighting && particle.pt() < antiprotonsInsideJets->GetXaxis()->GetXmax()) {
            int ipt = antiprotonsInsideJets->FindBin(particle.pt());
            weightJet = antiprotonsInsideJets->GetBinContent(ipt);
          }

          // Fill histogram for generated antiprotons
          registryMC.fill(HIST("antiproton_gen_jet"), particle.pt(), weightJet);

          // Fill 2d (pt,eta) distribution of antiprotons
          registryMC.fill(HIST("antiproton_eta_pt_jet"), particle.pt(), particle.eta(), weightJet);
        }

        // Set up two perpendicular cone axes for underlying event estimation
        TVector3 jetAxis(jet.px(), jet.py(), jet.pz());
        double coneRadius = std::sqrt(jet.area() / PI);
        TVector3 ueAxis1(0, 0, 0), ueAxis2(0, 0, 0);
        getPerpendicularDirections(jetAxis, ueAxis1, ueAxis2);
        if (ueAxis1.Mag() == 0 || ueAxis2.Mag() == 0) {
          continue;
        }

        // Loop over MC particles to analyze underlying event region
        for (const auto& protonVec : protonMomentum) {

          // Compute distance of particle from both perpendicular cone axes
          double deltaEtaUe1 = protonVec.Eta() - ueAxis1.Eta();
          double deltaPhiUe1 = getDeltaPhi(protonVec.Phi(), ueAxis1.Phi());
          double deltaRUe1 = std::sqrt(deltaEtaUe1 * deltaEtaUe1 + deltaPhiUe1 * deltaPhiUe1);
          double deltaEtaUe2 = protonVec.Eta() - ueAxis2.Eta();
          double deltaPhiUe2 = getDeltaPhi(protonVec.Phi(), ueAxis2.Phi());
          double deltaRUe2 = std::sqrt(deltaEtaUe2 * deltaEtaUe2 + deltaPhiUe2 * deltaPhiUe2);

          // Determine the maximum allowed distance from UE axes for particle selection
          double maxConeRadius = coneRadius;
          if (applyAreaCut) {
            maxConeRadius = std::sqrt(maxNormalizedJetArea) * rJet;
          }

          // Reject tracks that lie outside the maxConeRadius from both UE axes
          if (deltaRUe1 > maxConeRadius && deltaRUe2 > maxConeRadius)
            continue;

          // Fill normalization histogram
          registryMC.fill(HIST("antiproton_deltay_deltaphi_ue"), protonVec.Eta() - ueAxis1.Eta(), getDeltaPhi(protonVec.Phi(), ueAxis1.Phi()));
          registryMC.fill(HIST("antiproton_deltay_deltaphi_ue"), protonVec.Eta() - ueAxis2.Eta(), getDeltaPhi(protonVec.Phi(), ueAxis2.Phi()));

          // Calculate weight
          double weightUe(1.0);
          if (applyReweighting && protonVec.Pt() < antiprotonsPerpCone->GetXaxis()->GetXmax()) {
            int ipt = antiprotonsPerpCone->FindBin(protonVec.Pt());
            weightUe = antiprotonsPerpCone->GetBinContent(ipt);
          }

          // Fill histogram for antiprotons in the UE
          registryMC.fill(HIST("antiproton_gen_ue"), protonVec.Pt(), weightUe);

          // Fill 2d (pt,eta) distribution of antiprotons
          registryMC.fill(HIST("antiproton_eta_pt_ue"), protonVec.Pt(), protonVec.Eta(), weightUe);
        }
      }
      if (isAtLeastOneJetSelected) {
        registryMC.fill(HIST("genEvents"), 3.5);
      }

      // Shrink large vectors
      if (eventCounter % shrinkInterval == 0) {
        std::vector<fastjet::PseudoJet>().swap(fjParticles);
        std::vector<TVector3>().swap(protonMomentum);
      }
    }
  }
  PROCESS_SWITCH(AntinucleiInJets, processJetsMCgen, "process jets mc gen", false);

  // Reconstructed events
  void processJetsMCrec(RecCollisionsMc const& collisions, AntiNucleiTracksMc const& mcTracks, McParticles const&)
  {
    // Define per-event containers
    std::vector<fastjet::PseudoJet> fjParticles;
    std::vector<int> antiprotonTrackIndex;

    // Jet and area definitions
    fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, rJet);
    fastjet::AreaDefinition areaDef(fastjet::active_area, fastjet::GhostedAreaSpec(1.0));

    // Event counter
    int eventCounter = 0;

    // Loop over all reconstructed collisions
    for (const auto& collision : collisions) {

      // Configuration
      registryMC.fill(HIST("settingMC"), minJetPt.value, rJet.value);

      // Increment event counter
      eventCounter++;

      // Clear containers at the start of the event loop
      fjParticles.clear();
      antiprotonTrackIndex.clear();

      // Event counter: before event selection
      registryMC.fill(HIST("recEvents"), 0.5);

      // Apply event selection: require sel8 and vertex position to be within the allowed z range
      if (!collision.sel8() || std::fabs(collision.posZ()) > zVtx)
        continue;

      // Event counter: after event selection
      registryMC.fill(HIST("recEvents"), 1.5);

      // Reject events near the ITS Read-Out Frame border
      if (rejectITSROFBorder && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder))
        continue;
      registryMC.fill(HIST("recEvents"), 2.5);

      // Reject events at the Time Frame border
      if (rejectTFBorder && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder))
        continue;
      registryMC.fill(HIST("recEvents"), 3.5);

      // Require at least one ITS-TPC matched track
      if (requireVtxITSTPC && !collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC))
        continue;
      registryMC.fill(HIST("recEvents"), 4.5);

      // Reject events with same-bunch pileup
      if (rejectSameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup))
        continue;
      registryMC.fill(HIST("recEvents"), 5.5);

      // Require consistent FT0 vs PV z-vertex
      if (requireIsGoodZvtxFT0VsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV))
        continue;
      registryMC.fill(HIST("recEvents"), 6.5);

      // Require TOF match for at least one vertex track
      if (requireIsVertexTOFmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTOFmatched))
        continue;
      registryMC.fill(HIST("recEvents"), 7.5);

      // Get tracks in this MC collision
      const auto mcTracksThisMcColl = mcTracks.sliceBy(mcTracksPerMcCollision, collision.globalIndex());

      // Loop over reconstructed tracks
      int id(-1);
      for (auto const& track : mcTracksThisMcColl) {
        id++;

        // Get corresponding MC particle
        if (!track.has_mcParticle())
          continue;
        const auto mcparticle = track.mcParticle();

        // Store track index for antiproton tracks
        if (passedTrackSelection(track) && track.sign() < 0 && mcparticle.pdgCode() == PDG_t::kProtonBar) {
          antiprotonTrackIndex.emplace_back(id);

          double nsigmaTPCPr = track.tpcNSigmaPr();
          double nsigmaTOFPr = track.tofNSigmaPr();
          double pt = track.pt();
          double dcaxy = track.dcaXY();
          double dcaz = track.dcaZ();

          if (mcparticle.isPhysicalPrimary() && std::fabs(dcaxy) < maxDcaxy && std::fabs(dcaz) < maxDcaz && nsigmaTPCPr > minNsigmaTpc && nsigmaTPCPr < maxNsigmaTpc) {
            registryMC.fill(HIST("antiproton_rec_tpc_full"), pt);
            if (track.hasTOF() && nsigmaTOFPr > minNsigmaTof && nsigmaTOFPr < maxNsigmaTof) {
              registryMC.fill(HIST("antiproton_rec_tof_full"), pt);
            }
          }
        }

        // Apply track selection for jet reconstruction
        if (!passedTrackSelectionForJetReconstruction(track))
          continue;

        // 4-momentum representation of a particle
        fastjet::PseudoJet fourMomentum(track.px(), track.py(), track.pz(), track.energy(MassPionCharged));
        fourMomentum.set_user_index(id);
        fjParticles.emplace_back(fourMomentum);
      }

      // Reject empty events
      if (fjParticles.empty())
        continue;
      registryMC.fill(HIST("recEvents"), 8.5);

      // Cluster particles using the anti-kt algorithm
      fastjet::ClusterSequenceArea cs(fjParticles, jetDef, areaDef);
      std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cs.inclusive_jets());
      auto [rhoPerp, rhoMPerp] = jetutilities::estimateRhoPerpCone(fjParticles, jets[0], rJet);

      // Loop over reconstructed jets
      bool isAtLeastOneJetSelected = false;
      for (const auto& jet : jets) {

        // Jet must be fully contained in the acceptance
        if ((std::fabs(jet.eta()) + rJet) > (maxEta - deltaEtaEdge))
          continue;

        // Jet pt must be larger than threshold
        auto jetForSub = jet;
        fastjet::PseudoJet jetMinusBkg = backgroundSub.doRhoAreaSub(jetForSub, rhoPerp, rhoMPerp);
        if (jetMinusBkg.pt() < minJetPt)
          continue;

        // Apply area cut if required
        double normalizedJetArea = jet.area() / (PI * rJet * rJet);
        if (applyAreaCut && normalizedJetArea > maxNormalizedJetArea)
          continue;
        isAtLeastOneJetSelected = true;

        // Reconstructed jets
        registryMC.fill(HIST("recJets"), 0.5);

        // Set up two perpendicular cone axes for underlying event estimation
        double coneRadius = std::sqrt(jet.area() / PI);
        TVector3 jetAxis(jet.px(), jet.py(), jet.pz());
        TVector3 ueAxis1(0, 0, 0), ueAxis2(0, 0, 0);
        getPerpendicularDirections(jetAxis, ueAxis1, ueAxis2);
        if (ueAxis1.Mag() == 0 || ueAxis2.Mag() == 0) {
          continue;
        }

        // Get jet constituents
        std::vector<fastjet::PseudoJet> jetConstituents = jet.constituents();

        // Loop over jet constituents
        for (const auto& particle : jetConstituents) {

          // Get corresponding track and apply track selection criteria
          auto const& track = mcTracksThisMcColl.iteratorAt(particle.user_index());
          if (!passedTrackSelection(track))
            continue;

          // Antimatter selection
          if (track.sign() > 0)
            continue;

          // Get corresponding MC particle
          if (!track.has_mcParticle())
            continue;
          const auto mcparticle = track.mcParticle();

          // Define variables
          double nsigmaTPCPr = track.tpcNSigmaPr();
          double nsigmaTOFPr = track.tofNSigmaPr();
          double pt = track.pt();
          double dcaxy = track.dcaXY();
          double dcaz = track.dcaZ();

          // Fill nsigmaTOF template
          if (track.hasTOF() && std::fabs(dcaxy) < maxDcaxy && std::fabs(dcaz) < maxDcaz && mcparticle.isPhysicalPrimary() && nsigmaTPCPr > minNsigmaTpc && nsigmaTPCPr < maxNsigmaTpc && mcparticle.pdgCode() != PDG_t::kProtonBar) {
            registryMC.fill(HIST("antiproton_nsigma_tof_jet_mc"), pt, nsigmaTOFPr);
          }

          // Antiproton selection based on the PDG
          if (mcparticle.pdgCode() != PDG_t::kProtonBar)
            continue;

          // Fill DCA templates
          if (std::fabs(dcaz) < maxDcaz) {
            if (mcparticle.isPhysicalPrimary()) {
              registryMC.fill(HIST("antiproton_prim_dca_jet"), pt, dcaxy);
            } else {
              registryMC.fill(HIST("antiproton_all_dca_jet"), pt, dcaxy);
            }
          }

          // Apply DCA selections
          if (std::fabs(dcaxy) > maxDcaxy || std::fabs(dcaz) > maxDcaz)
            continue;

          // nsigmaITS for antiprotons
          double nSigmaITSprot = static_cast<double>(itsResponse.nSigmaITS<o2::track::PID::Proton>(track));

          // Particle identification using the ITS cluster size
          bool passedItsPidProt(true);
          if (applyItsPid && pt < ptMaxItsPidProt && (nSigmaITSprot < nSigmaItsMin || nSigmaITSprot > nSigmaItsMax)) {
            passedItsPidProt = false;
          }

          // Fill inclusive antiproton spectrum
          registryMC.fill(HIST("antiproton_incl_jet"), pt);

          // Select physical primary antiprotons
          if (!mcparticle.isPhysicalPrimary())
            continue;

          // Fill antiproton spectrum for physical primaries
          registryMC.fill(HIST("antiproton_prim_jet"), pt);

          // Calculate weight
          double weightJet(1.0);
          if (applyReweighting && mcparticle.pt() < antiprotonsInsideJets->GetXaxis()->GetXmax()) {
            int ipt = antiprotonsInsideJets->FindBin(mcparticle.pt());
            weightJet = antiprotonsInsideJets->GetBinContent(ipt);
          }

          // Fill histograms (TPC and TOF) only for selected candidates
          if (passedItsPidProt && nsigmaTPCPr > minNsigmaTpc && nsigmaTPCPr < maxNsigmaTpc) {
            registryMC.fill(HIST("antiproton_rec_tpc_jet"), pt, weightJet);
            if (track.hasTOF() && nsigmaTOFPr > minNsigmaTof && nsigmaTOFPr < maxNsigmaTof) {
              registryMC.fill(HIST("antiproton_rec_tof_jet"), pt, weightJet);
            }
          }
        }

        // Loop over tracks in the underlying event
        for (auto const& index : antiprotonTrackIndex) {

          // retrieve track associated to index
          auto const& track = mcTracksThisMcColl.iteratorAt(index);

          // Get corresponding MC particle
          if (!track.has_mcParticle())
            continue;
          const auto mcparticle = track.mcParticle();

          // Define variables
          double nsigmaTPCPr = track.tpcNSigmaPr();
          double nsigmaTOFPr = track.tofNSigmaPr();
          double pt = track.pt();
          double dcaxy = track.dcaXY();
          double dcaz = track.dcaZ();

          // Fill DCA templates
          if (std::fabs(dcaz) < maxDcaz) {
            if (mcparticle.isPhysicalPrimary()) {
              registryMC.fill(HIST("antiproton_prim_dca_ue"), pt, dcaxy);
            } else {
              registryMC.fill(HIST("antiproton_all_dca_ue"), pt, dcaxy);
            }
          }

          // Apply DCA selection
          if (std::fabs(dcaxy) > maxDcaxy || std::fabs(dcaz) > maxDcaz)
            continue;

          // Calculate the angular distance between the track and underlying event axes in eta-phi space
          double deltaEtaUe1 = track.eta() - ueAxis1.Eta();
          double deltaPhiUe1 = getDeltaPhi(track.phi(), ueAxis1.Phi());
          double deltaRUe1 = std::sqrt(deltaEtaUe1 * deltaEtaUe1 + deltaPhiUe1 * deltaPhiUe1);
          double deltaEtaUe2 = track.eta() - ueAxis2.Eta();
          double deltaPhiUe2 = getDeltaPhi(track.phi(), ueAxis2.Phi());
          double deltaRUe2 = std::sqrt(deltaEtaUe2 * deltaEtaUe2 + deltaPhiUe2 * deltaPhiUe2);

          // Determine the maximum allowed distance from UE axes for particle selection
          double maxConeRadius = coneRadius;
          if (applyAreaCut) {
            maxConeRadius = std::sqrt(maxNormalizedJetArea) * rJet;
          }

          // Reject tracks that lie outside the maxConeRadius from both UE axes
          if (deltaRUe1 > maxConeRadius && deltaRUe2 > maxConeRadius)
            continue;

          // Particle identification using the ITS cluster size
          bool passedItsPidProt(true);
          double nSigmaITSprot = static_cast<double>(itsResponse.nSigmaITS<o2::track::PID::Proton>(track));
          if (applyItsPid && pt < ptMaxItsPidProt && (nSigmaITSprot < nSigmaItsMin || nSigmaITSprot > nSigmaItsMax)) {
            passedItsPidProt = false;
          }

          // Fill inclusive antiproton spectrum
          registryMC.fill(HIST("antiproton_incl_ue"), pt);

          // Select physical primary antiprotons
          if (!mcparticle.isPhysicalPrimary())
            continue;

          // Fill antiproton spectrum for physical primaries
          registryMC.fill(HIST("antiproton_prim_ue"), pt);

          // Calculate weight
          double weightUe(1.0);
          if (applyReweighting && mcparticle.pt() < antiprotonsPerpCone->GetXaxis()->GetXmax()) {
            int ipt = antiprotonsPerpCone->FindBin(mcparticle.pt());
            weightUe = antiprotonsPerpCone->GetBinContent(ipt);
          }

          // Fill histograms (TPC and TOF) only for selected candidates
          if (passedItsPidProt && nsigmaTPCPr > minNsigmaTpc && nsigmaTPCPr < maxNsigmaTpc) {
            registryMC.fill(HIST("antiproton_rec_tpc_ue"), pt, weightUe);
            if (track.hasTOF() && nsigmaTOFPr > minNsigmaTof && nsigmaTOFPr < maxNsigmaTof) {
              registryMC.fill(HIST("antiproton_rec_tof_ue"), pt, weightUe);
            }
          }
        }
      }
      if (isAtLeastOneJetSelected) {
        registryMC.fill(HIST("recEvents"), 9.5);
      }

      // Shrink large vectors
      if (eventCounter % shrinkInterval == 0) {
        std::vector<fastjet::PseudoJet>().swap(fjParticles);
        std::vector<int>().swap(antiprotonTrackIndex);
      }
    }
  }
  PROCESS_SWITCH(AntinucleiInJets, processJetsMCrec, "process jets MC rec", false);

  // Process real data with systematic variations of analysis parameters
  void processSystData(SelectedCollisions::iterator const& collision, AntiNucleiTracks const& tracks)
  {
    // Event counter: before event selection
    registryData.fill(HIST("number_of_events_data_syst"), 0.5);

    // Apply standard event selection
    if (!collision.sel8() || std::fabs(collision.posZ()) > zVtx)
      return;

    // Event counter: after event selection
    registryData.fill(HIST("number_of_events_data_syst"), 1.5);

    // Reject events near the ITS Read-Out Frame border
    if (rejectITSROFBorder && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder))
      return;
    registryData.fill(HIST("number_of_events_data_syst"), 2.5);

    // Reject events at the Time Frame border
    if (rejectTFBorder && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder))
      return;
    registryData.fill(HIST("number_of_events_data_syst"), 3.5);

    // Require at least one ITS-TPC matched track
    if (requireVtxITSTPC && !collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC))
      return;
    registryData.fill(HIST("number_of_events_data_syst"), 4.5);

    // Reject events with same-bunch pileup
    if (rejectSameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup))
      return;
    registryData.fill(HIST("number_of_events_data_syst"), 5.5);

    // Require consistent FT0 vs PV z-vertex
    if (requireIsGoodZvtxFT0VsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV))
      return;
    registryData.fill(HIST("number_of_events_data_syst"), 6.5);

    // Require TOF match for at least one vertex track
    if (requireIsVertexTOFmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTOFmatched))
      return;
    registryData.fill(HIST("number_of_events_data_syst"), 7.5);

    // Cut settings
    static std::vector<double> maxDcaxySyst = {
      0.071, 0.060, 0.066, 0.031, 0.052, 0.078, 0.045, 0.064, 0.036, 0.074,
      0.079, 0.043, 0.067, 0.059, 0.032, 0.070, 0.048, 0.077, 0.062, 0.034,
      0.057, 0.055, 0.073, 0.038, 0.050, 0.075, 0.041, 0.061, 0.033, 0.069,
      0.035, 0.044, 0.076, 0.049, 0.037, 0.054, 0.072, 0.046, 0.058, 0.040,
      0.068, 0.042, 0.056, 0.039, 0.047, 0.065, 0.051, 0.053, 0.063, 0.030};

    static std::vector<double> maxDcazSyst = {
      0.064, 0.047, 0.032, 0.076, 0.039, 0.058, 0.043, 0.069, 0.050, 0.035,
      0.074, 0.061, 0.045, 0.033, 0.068, 0.055, 0.037, 0.071, 0.042, 0.053,
      0.077, 0.038, 0.065, 0.049, 0.036, 0.059, 0.044, 0.067, 0.041, 0.034,
      0.073, 0.052, 0.040, 0.063, 0.046, 0.031, 0.070, 0.054, 0.037, 0.062,
      0.048, 0.035, 0.075, 0.051, 0.039, 0.066, 0.043, 0.060, 0.032, 0.056};

    static std::vector<double> nSigmaItsMinSyst = {
      -2.9, -2.8, -3.0, -3.4, -2.7, -3.3, -3.0, -3.1, -3.4, -3.1,
      -3.0, -2.8, -3.2, -2.6, -2.7, -3.4, -2.9, -3.0, -3.0, -2.7,
      -2.9, -3.3, -3.0, -3.1, -3.2, -3.0, -2.9, -2.7, -3.3, -3.0,
      -2.8, -3.3, -2.7, -3.3, -2.8, -3.4, -2.8, -3.4, -2.9, -3.1,
      -3.2, -2.6, -3.1, -2.9, -3.1, -2.8, -2.9, -3.3, -3.0, -2.8};

    static std::vector<double> nSigmaItsMaxSyst = {
      2.9, 2.8, 3.0, 3.4, 2.7, 3.3, 3.0, 3.1, 3.4, 3.1,
      3.0, 2.8, 3.2, 2.6, 2.7, 3.4, 2.9, 3.0, 3.0, 2.7,
      2.9, 3.3, 3.0, 3.1, 3.2, 3.0, 2.9, 2.7, 3.3, 3.0,
      2.8, 3.3, 2.7, 3.3, 2.8, 3.4, 2.8, 3.4, 2.9, 3.1,
      3.2, 2.6, 3.1, 2.9, 3.1, 2.8, 2.9, 3.3, 3.0, 2.8};

    static std::vector<double> minNsigmaTpcSyst = {
      -3.2, -2.9, -3.1, -2.9, -3.5, -2.6, -3.3, -3.0, -3.5, -2.7,
      -3.0, -2.6, -3.3, -3.4, -2.8, -3.1, -2.6, -3.2, -3.1, -2.8,
      -3.4, -2.7, -3.4, -2.9, -3.0, -2.5, -3.3, -2.8, -3.1, -2.7,
      -3.4, -2.8, -3.3, -2.6, -3.1, -2.5, -3.4, -3.0, -3.2, -2.6,
      -3.4, -2.8, -3.1, -2.6, -3.3, -2.7, -3.2, -2.7, -3.4, -2.9};

    static std::vector<double> maxNsigmaTpcSyst = {
      3.2, 2.9, 3.1, 2.9, 3.5, 2.6, 3.3, 3.0, 3.5, 2.7,
      3.0, 2.6, 3.3, 3.4, 2.8, 3.1, 2.6, 3.2, 3.1, 2.8,
      3.4, 2.7, 3.4, 2.9, 3.0, 2.5, 3.3, 2.8, 3.1, 2.7,
      3.4, 2.8, 3.3, 2.6, 3.1, 2.5, 3.4, 3.0, 3.2, 2.6,
      3.4, 2.8, 3.1, 2.6, 3.3, 2.7, 3.2, 2.7, 3.4, 2.9};

    // Loop over reconstructed tracks
    for (auto const& track : tracks) {

      // Select only antimatter
      if (track.sign() > 0)
        continue;

      // Loop over different cut settings
      for (int isyst = 0; isyst < nSyst; isyst++) {

        // Apply track selection
        if (!passedTrackSelectionSyst(track, isyst))
          continue;

        // Define variables
        double nsigmaTPCPr = track.tpcNSigmaPr();
        double nsigmaTOFPr = track.tofNSigmaPr();
        double nsigmaTPCDe = track.tpcNSigmaDe();
        double nsigmaTOFDe = track.tofNSigmaDe();
        double nsigmaTPCHe = track.tpcNSigmaHe();
        double pt = track.pt();
        double dcaxy = track.dcaXY();
        double dcaz = track.dcaZ();

        // Apply DCA selections
        if (std::fabs(dcaxy) > maxDcaxySyst[isyst] || std::fabs(dcaz) > maxDcazSyst[isyst])
          continue;

        // Particle identification using the ITS cluster size (vary also PID ITS)
        bool passedItsPidProt(true), passedItsPidDeut(true), passedItsPidHel(true);
        double nSigmaITSprot = static_cast<double>(itsResponse.nSigmaITS<o2::track::PID::Proton>(track));
        double nSigmaITSdeut = static_cast<double>(itsResponse.nSigmaITS<o2::track::PID::Deuteron>(track));
        double nSigmaITShel3 = static_cast<double>(itsResponse.nSigmaITS<o2::track::PID::Helium3>(track));

        if (applyItsPid && pt < ptMaxItsPidProt && (nSigmaITSprot < nSigmaItsMinSyst[isyst] || nSigmaITSprot > nSigmaItsMaxSyst[isyst])) {
          passedItsPidProt = false;
        }
        if (applyItsPid && pt < ptMaxItsPidDeut && (nSigmaITSdeut < nSigmaItsMinSyst[isyst] || nSigmaITSdeut > nSigmaItsMaxSyst[isyst])) {
          passedItsPidDeut = false;
        }
        if (applyItsPid && (2.0 * pt) < ptMaxItsPidHel && (nSigmaITShel3 < nSigmaItsMinSyst[isyst] || nSigmaITShel3 > nSigmaItsMaxSyst[isyst])) {
          passedItsPidHel = false;
        }

        // Fill histograms
        if (passedItsPidProt) {
          registryData.fill(HIST("antiproton_tpc_syst"), isyst, pt, nsigmaTPCPr);
          if (nsigmaTPCPr > minNsigmaTpcSyst[isyst] && nsigmaTPCPr < maxNsigmaTpcSyst[isyst] && track.hasTOF())
            registryData.fill(HIST("antiproton_tof_syst"), isyst, pt, nsigmaTOFPr);
        }
        if (passedItsPidDeut) {
          registryData.fill(HIST("antideuteron_tpc_syst"), isyst, pt, nsigmaTPCDe);
          if (nsigmaTPCDe > minNsigmaTpcSyst[isyst] && nsigmaTPCDe < maxNsigmaTpcSyst[isyst] && track.hasTOF())
            registryData.fill(HIST("antideuteron_tof_syst"), isyst, pt, nsigmaTOFDe);
        }
        if (passedItsPidHel) {
          registryData.fill(HIST("antihelium3_tpc_syst"), isyst, 2.0 * pt, nsigmaTPCHe);
        }
      }
    }
  }
  PROCESS_SWITCH(AntinucleiInJets, processSystData, "Process syst data", false);

  // Process MC with systematic variations of analysis parameters
  void processSystEff(GenCollisionsMc const& genCollisions, RecCollisionsMc const& recCollisions, AntiNucleiTracksMc const& mcTracks, aod::McParticles const& mcParticles)
  {
    // Cut settings
    static std::vector<double> maxDcaxySyst = {
      0.071, 0.060, 0.066, 0.031, 0.052, 0.078, 0.045, 0.064, 0.036, 0.074,
      0.079, 0.043, 0.067, 0.059, 0.032, 0.070, 0.048, 0.077, 0.062, 0.034,
      0.057, 0.055, 0.073, 0.038, 0.050, 0.075, 0.041, 0.061, 0.033, 0.069,
      0.035, 0.044, 0.076, 0.049, 0.037, 0.054, 0.072, 0.046, 0.058, 0.040,
      0.068, 0.042, 0.056, 0.039, 0.047, 0.065, 0.051, 0.053, 0.063, 0.030};

    static std::vector<double> maxDcazSyst = {
      0.064, 0.047, 0.032, 0.076, 0.039, 0.058, 0.043, 0.069, 0.050, 0.035,
      0.074, 0.061, 0.045, 0.033, 0.068, 0.055, 0.037, 0.071, 0.042, 0.053,
      0.077, 0.038, 0.065, 0.049, 0.036, 0.059, 0.044, 0.067, 0.041, 0.034,
      0.073, 0.052, 0.040, 0.063, 0.046, 0.031, 0.070, 0.054, 0.037, 0.062,
      0.048, 0.035, 0.075, 0.051, 0.039, 0.066, 0.043, 0.060, 0.032, 0.056};

    static std::vector<double> nSigmaItsMinSyst = {
      -2.9, -2.8, -3.0, -3.4, -2.7, -3.3, -3.0, -3.1, -3.4, -3.1,
      -3.0, -2.8, -3.2, -2.6, -2.7, -3.4, -2.9, -3.0, -3.0, -2.7,
      -2.9, -3.3, -3.0, -3.1, -3.2, -3.0, -2.9, -2.7, -3.3, -3.0,
      -2.8, -3.3, -2.7, -3.3, -2.8, -3.4, -2.8, -3.4, -2.9, -3.1,
      -3.2, -2.6, -3.1, -2.9, -3.1, -2.8, -2.9, -3.3, -3.0, -2.8};

    static std::vector<double> nSigmaItsMaxSyst = {
      2.9, 2.8, 3.0, 3.4, 2.7, 3.3, 3.0, 3.1, 3.4, 3.1,
      3.0, 2.8, 3.2, 2.6, 2.7, 3.4, 2.9, 3.0, 3.0, 2.7,
      2.9, 3.3, 3.0, 3.1, 3.2, 3.0, 2.9, 2.7, 3.3, 3.0,
      2.8, 3.3, 2.7, 3.3, 2.8, 3.4, 2.8, 3.4, 2.9, 3.1,
      3.2, 2.6, 3.1, 2.9, 3.1, 2.8, 2.9, 3.3, 3.0, 2.8};

    static std::vector<double> minNsigmaTpcSyst = {
      -3.2, -2.9, -3.1, -2.9, -3.5, -2.6, -3.3, -3.0, -3.5, -2.7,
      -3.0, -2.6, -3.3, -3.4, -2.8, -3.1, -2.6, -3.2, -3.1, -2.8,
      -3.4, -2.7, -3.4, -2.9, -3.0, -2.5, -3.3, -2.8, -3.1, -2.7,
      -3.4, -2.8, -3.3, -2.6, -3.1, -2.5, -3.4, -3.0, -3.2, -2.6,
      -3.4, -2.8, -3.1, -2.6, -3.3, -2.7, -3.2, -2.7, -3.4, -2.9};

    static std::vector<double> maxNsigmaTpcSyst = {
      3.2, 2.9, 3.1, 2.9, 3.5, 2.6, 3.3, 3.0, 3.5, 2.7,
      3.0, 2.6, 3.3, 3.4, 2.8, 3.1, 2.6, 3.2, 3.1, 2.8,
      3.4, 2.7, 3.4, 2.9, 3.0, 2.5, 3.3, 2.8, 3.1, 2.7,
      3.4, 2.8, 3.3, 2.6, 3.1, 2.5, 3.4, 3.0, 3.2, 2.6,
      3.4, 2.8, 3.1, 2.6, 3.3, 2.7, 3.2, 2.7, 3.4, 2.9};

    static std::vector<double> minNsigmaTofSyst = {
      -3.2, -2.9, -3.1, -2.9, -3.5, -2.6, -3.3, -3.0, -3.5, -2.7,
      -3.0, -2.6, -3.3, -3.4, -2.8, -3.1, -2.6, -3.2, -3.1, -2.8,
      -3.4, -2.7, -3.4, -2.9, -3.0, -2.5, -3.3, -2.8, -3.1, -2.7,
      -3.4, -2.8, -3.3, -2.6, -3.1, -2.5, -3.4, -3.0, -3.2, -2.6,
      -3.4, -2.8, -3.1, -2.6, -3.3, -2.7, -3.2, -2.7, -3.4, -2.9};

    static std::vector<double> maxNsigmaTofSyst = {
      3.9, 3.6, 3.8, 3.2, 3.2, 3.5, 3.1, 3.8, 3.5, 3.4,
      3.9, 3.8, 3.7, 3.0, 3.6, 3.1, 3.7, 3.4, 4.0, 3.0,
      3.7, 3.3, 3.9, 3.1, 3.3, 3.5, 3.6, 3.2, 3.5, 3.3,
      3.9, 3.0, 3.4, 3.2, 3.1, 3.9, 3.6, 3.1, 3.2, 4.0,
      3.1, 3.7, 3.6, 3.1, 3.3, 3.5, 3.3, 3.4, 3.1, 3.8};

    // Loop over generated collisions
    for (const auto& collision : genCollisions) {

      // Apply event selection: require vertex position to be within the allowed z range
      if (std::fabs(collision.posZ()) > zVtx)
        continue;

      // Get particles in this MC collision
      const auto mcParticlesThisMcColl = mcParticles.sliceBy(mcParticlesPerMcCollision, collision.globalIndex());

      // Loop over all generated Monte Carlo particles for the selected event
      for (const auto& particle : mcParticlesThisMcColl) {

        // Select primary particles
        if (!particle.isPhysicalPrimary())
          continue;

        // Select particles within the specified pseudorapidity interval
        if (particle.eta() < minEta || particle.eta() > maxEta)
          continue;

        // Process different particle species based on PDG code
        switch (particle.pdgCode()) {
          case PDG_t::kProtonBar:
            registryMC.fill(HIST("antiproton_gen_syst"), particle.pt());
            break;
          case -o2::constants::physics::Pdg::kDeuteron:
            registryMC.fill(HIST("antideuteron_gen_syst"), particle.pt());
            break;
          case -o2::constants::physics::Pdg::kHelium3:
            registryMC.fill(HIST("antihelium3_gen_syst"), particle.pt());
            break;
        }
      }
    }

    // Loop over reconstructed collisions
    for (const auto& collision : recCollisions) {

      // Apply standard event selection
      if (!collision.sel8() || std::fabs(collision.posZ()) > zVtx)
        continue;

      // Reject events near the ITS Read-Out Frame border
      if (rejectITSROFBorder && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder))
        continue;

      // Reject events at the Time Frame border
      if (rejectTFBorder && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder))
        continue;

      // Require at least one ITS-TPC matched track
      if (requireVtxITSTPC && !collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC))
        continue;

      // Reject events with same-bunch pileup
      if (rejectSameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup))
        continue;

      // Require consistent FT0 vs PV z-vertex
      if (requireIsGoodZvtxFT0VsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV))
        continue;

      // Require TOF match for at least one vertex track
      if (requireIsVertexTOFmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTOFmatched))
        continue;

      // Get tracks in this MC collision
      const auto mcTracksThisMcColl = mcTracks.sliceBy(mcTracksPerMcCollision, collision.globalIndex());

      // Loop over reconstructed tracks
      for (auto const& track : mcTracksThisMcColl) {

        // Select only antimatter
        if (track.sign() > 0)
          continue;

        // Get corresponding MC particle
        if (!track.has_mcParticle())
          continue;
        const auto mcparticle = track.mcParticle();

        // Loop over different cut settings
        for (int isyst = 0; isyst < nSyst; isyst++) {

          // Apply track selection
          if (!passedTrackSelectionSyst(track, isyst))
            continue;

          // Define variables
          double nsigmaTPCPr = track.tpcNSigmaPr();
          double nsigmaTOFPr = track.tofNSigmaPr();
          double nsigmaTPCDe = track.tpcNSigmaDe();
          double nsigmaTOFDe = track.tofNSigmaDe();
          double nsigmaTPCHe = track.tpcNSigmaHe();
          double pt = track.pt();
          double dcaxy = track.dcaXY();
          double dcaz = track.dcaZ();

          // Apply DCA selections
          if (std::fabs(dcaxy) > maxDcaxySyst[isyst] || std::fabs(dcaz) > maxDcazSyst[isyst])
            continue;

          // Fill inclusive antiproton spectrum
          registryMC.fill(HIST("antiproton_incl_syst"), isyst, pt);

          // Select physical primary antiprotons
          if (!mcparticle.isPhysicalPrimary())
            continue;

          // Fill antiproton spectrum for physical primaries
          registryMC.fill(HIST("antiproton_prim_syst"), isyst, pt);

          // Particle identification using the ITS cluster size (vary also PID ITS)
          bool passedItsPidProt(true), passedItsPidDeut(true), passedItsPidHel(true);
          double nSigmaITSprot = static_cast<double>(itsResponse.nSigmaITS<o2::track::PID::Proton>(track));
          double nSigmaITSdeut = static_cast<double>(itsResponse.nSigmaITS<o2::track::PID::Deuteron>(track));
          double nSigmaITShel3 = static_cast<double>(itsResponse.nSigmaITS<o2::track::PID::Helium3>(track));

          if (applyItsPid && pt < ptMaxItsPidProt && (nSigmaITSprot < nSigmaItsMinSyst[isyst] || nSigmaITSprot > nSigmaItsMaxSyst[isyst])) {
            passedItsPidProt = false;
          }
          if (applyItsPid && pt < ptMaxItsPidDeut && (nSigmaITSdeut < nSigmaItsMinSyst[isyst] || nSigmaITSdeut > nSigmaItsMaxSyst[isyst])) {
            passedItsPidDeut = false;
          }
          if (applyItsPid && (2.0 * pt) < ptMaxItsPidHel && (nSigmaITShel3 < nSigmaItsMinSyst[isyst] || nSigmaITShel3 > nSigmaItsMaxSyst[isyst])) {
            passedItsPidHel = false;
          }

          // Fill histograms for antiprotons
          if (passedItsPidProt && mcparticle.pdgCode() == PDG_t::kProtonBar && nsigmaTPCPr > minNsigmaTpcSyst[isyst] && nsigmaTPCPr < maxNsigmaTpcSyst[isyst]) {
            registryMC.fill(HIST("antiproton_rec_tpc_syst"), isyst, pt);
            if (track.hasTOF() && nsigmaTOFPr > minNsigmaTofSyst[isyst] && nsigmaTOFPr < maxNsigmaTofSyst[isyst])
              registryMC.fill(HIST("antiproton_rec_tof_syst"), isyst, pt);
          }
          // Fill histograms for antideuterons
          if (passedItsPidDeut && mcparticle.pdgCode() == -o2::constants::physics::Pdg::kDeuteron && nsigmaTPCDe > minNsigmaTpcSyst[isyst] && nsigmaTPCDe < maxNsigmaTpcSyst[isyst]) {
            registryMC.fill(HIST("antideuteron_rec_tpc_syst"), isyst, pt);
            if (track.hasTOF() && nsigmaTOFDe > minNsigmaTofSyst[isyst] && nsigmaTOFDe < maxNsigmaTofSyst[isyst])
              registryMC.fill(HIST("antideuteron_rec_tof_syst"), isyst, pt);
          }
          // Fill histograms for antihelium3
          if (passedItsPidHel && mcparticle.pdgCode() == -o2::constants::physics::Pdg::kHelium3 && nsigmaTPCHe > minNsigmaTpcSyst[isyst] && nsigmaTPCHe < maxNsigmaTpcSyst[isyst]) {
            registryMC.fill(HIST("antihelium3_rec_tpc_syst"), isyst, 2.0 * pt);
          }
        }
      }
    }
  }
  PROCESS_SWITCH(AntinucleiInJets, processSystEff, "process syst mc", false);

  // Process correlation
  void processCorr(SelectedCollisions::iterator const& collision, AntiNucleiTracks const& tracks)
  {
    // Event counter: before event selection
    registryCorr.fill(HIST("eventCounter"), 0.5);

    // Apply standard event selection
    if (!collision.sel8() || std::fabs(collision.posZ()) > zVtx)
      return;

    // Event counter: after event selection
    registryCorr.fill(HIST("eventCounter"), 1.5);

    // Reject events near the ITS Read-Out Frame border
    if (rejectITSROFBorder && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder))
      return;
    registryCorr.fill(HIST("eventCounter"), 2.5);

    // Reject events at the Time Frame border
    if (rejectTFBorder && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder))
      return;
    registryCorr.fill(HIST("eventCounter"), 3.5);

    // Require at least one ITS-TPC matched track
    if (requireVtxITSTPC && !collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC))
      return;
    registryCorr.fill(HIST("eventCounter"), 4.5);

    // Reject events with same-bunch pileup
    if (rejectSameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup))
      return;
    registryCorr.fill(HIST("eventCounter"), 5.5);

    // Require consistent FT0 vs PV z-vertex
    if (requireIsGoodZvtxFT0VsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV))
      return;
    registryCorr.fill(HIST("eventCounter"), 6.5);

    // Require TOF match for at least one vertex track
    if (requireIsVertexTOFmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTOFmatched))
      return;
    registryCorr.fill(HIST("eventCounter"), 7.5);

    // Assign event to a random subsample (0-19)
    double sampleId = mRand.Integer(nSubsamples) + 0.5;

    // Multiplicity percentile
    const float multiplicity = collision.centFT0M();

    // Fill event counter vs centrality (full Event region)
    registryCorr.fill(HIST("eventCounter_centrality_fullEvent"), multiplicity, sampleId);

    // pt/A bins
    std::vector<double> ptOverAbins = {0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
    const int nBins = ptOverAbins.size() - 1;

    // Particle counters
    std::vector<int> nAntiprotonFullEvent(nBins, 0);
    std::vector<int> nAntideuteronFullEvent(nBins, 0);
    int nTotProtonFullEvent(0);
    int nTotDeuteronFullEvent(0);
    int nTotAntiprotonFullEvent(0);
    int nTotAntideuteronFullEvent(0);

    // Loop over reconstructed tracks
    for (auto const& track : tracks) {

      // Apply track selection
      if (!passedTrackSelection(track))
        continue;

      // Apply DCA selections
      if (std::fabs(track.dcaXY()) > maxDcaxy || std::fabs(track.dcaZ()) > maxDcaz)
        continue;

      // Particle identification using the ITS cluster size
      bool passedItsPidProt(true), passedItsPidDeut(true);
      double nSigmaITSprot = static_cast<double>(itsResponse.nSigmaITS<o2::track::PID::Proton>(track));
      double nSigmaITSdeut = static_cast<double>(itsResponse.nSigmaITS<o2::track::PID::Deuteron>(track));

      if (applyItsPid && track.pt() < ptMaxItsPidProt && (nSigmaITSprot < nSigmaItsMin || nSigmaITSprot > nSigmaItsMax)) {
        passedItsPidProt = false;
      }
      if (applyItsPid && track.pt() < ptMaxItsPidDeut && (nSigmaITSdeut < nSigmaItsMin || nSigmaITSdeut > nSigmaItsMax)) {
        passedItsPidDeut = false;
      }

      // Kinematic range selection
      if (isProton(track) && passedItsPidProt) {
        if (track.pt() < ptOverAbins[0] || track.pt() >= ptOverAbins[nBins]) {
          continue;
        }
      } else if (isDeuteron(track) && passedItsPidDeut) {
        double ptPerNucleon = 0.5 * track.pt();
        if (ptPerNucleon < ptOverAbins[0] || ptPerNucleon >= ptOverAbins[nBins]) {
          continue;
        }
      } else {
        continue;
      }

      // (Anti)protons
      if (isProton(track) && passedItsPidProt) {
        if (track.sign() > 0) {
          nTotProtonFullEvent++;
        } else if (track.sign() < 0) {
          nTotAntiprotonFullEvent++;
          int ibin = findBin(ptOverAbins, track.pt());
          nAntiprotonFullEvent[ibin]++;
        }
      }

      // (Anti)deuterons
      if (isDeuteron(track) && passedItsPidDeut) {
        const double ptPerNucleon = 0.5 * track.pt();

        if (track.sign() > 0) {
          nTotDeuteronFullEvent++;
        } else if (track.sign() < 0) {
          nTotAntideuteronFullEvent++;
          int ibin = findBin(ptOverAbins, ptPerNucleon);
          nAntideuteronFullEvent[ibin]++;
        }
      }
    }

    // Fill correlation histograms
    int netProtonFullEvent = nTotProtonFullEvent - nTotAntiprotonFullEvent;
    int netDeuteronFullEvent = nTotDeuteronFullEvent - nTotAntideuteronFullEvent;

    registryCorr.fill(HIST("rho_fullEvent"), nTotAntideuteronFullEvent, nTotAntiprotonFullEvent, multiplicity, sampleId);
    registryCorr.fill(HIST("rho_netP_netD_fullEvent"), netDeuteronFullEvent, netProtonFullEvent, sampleId);

    // Fill efficiency histograms
    for (int i = 0; i < nBins; i++) {
      double ptAcenteri = 0.5 * (ptOverAbins[i] + ptOverAbins[i + 1]);

      registryCorr.fill(HIST("q1d_fullEvent"), nAntideuteronFullEvent[i], ptAcenteri, multiplicity, sampleId);
      registryCorr.fill(HIST("q1p_fullEvent"), nAntiprotonFullEvent[i], ptAcenteri, multiplicity, sampleId);
      for (int j = 0; j < nBins; j++) {
        double ptAcenterj = 0.5 * (ptOverAbins[j] + ptOverAbins[j + 1]);
        registryCorr.fill(HIST("q1d_square_fullEvent"), ptAcenteri, ptAcenterj, (nAntideuteronFullEvent[i] * nAntideuteronFullEvent[j]), multiplicity, sampleId);
        registryCorr.fill(HIST("q1p_square_fullEvent"), ptAcenteri, ptAcenterj, (nAntiprotonFullEvent[i] * nAntiprotonFullEvent[j]), multiplicity, sampleId);
        registryCorr.fill(HIST("q1d_q1p_fullEvent"), ptAcenteri, ptAcenterj, (nAntideuteronFullEvent[i] * nAntiprotonFullEvent[j]), multiplicity, sampleId);
      }
    }

    /*
    // Loop over reconstructed tracks (refactoring: this part can be incorporated above)
    int id(-1);
    std::vector<fastjet::PseudoJet> fjParticles;
    for (auto const& track : tracks) {
      id++;
      if (!passedTrackSelectionForJetReconstruction(track))
        continue;

      // 4-momentum representation of a particle
      fastjet::PseudoJet fourMomentum(track.px(), track.py(), track.pz(), track.energy(MassPionCharged));
      fourMomentum.set_user_index(id);
      fjParticles.emplace_back(fourMomentum);
    }

    // Reject empty events
    if (fjParticles.empty())
      return;
    registryCorr.fill(HIST("eventCounter"), 8.5);

    // Cluster particles using the anti-kt algorithm
    fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, rJet);
    fastjet::AreaDefinition areaDef(fastjet::active_area, fastjet::GhostedAreaSpec(1.0));
    fastjet::ClusterSequenceArea cs(fjParticles, jetDef, areaDef);
    std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cs.inclusive_jets());
    auto [rhoPerp, rhoMPerp] = jetutilities::estimateRhoPerpCone(fjParticles, jets[0], rJet);

    // Loop over reconstructed jets
    bool isAtLeastOneJetSelected = false;
    for (const auto& jet : jets) {

      // Jet must be fully contained in the acceptance
      if ((std::fabs(jet.eta()) + rJet) > (maxEta - deltaEtaEdge))
        continue;

      // Jet pt must be larger than threshold
      auto jetForSub = jet;
      fastjet::PseudoJet jetMinusBkg = backgroundSub.doRhoAreaSub(jetForSub, rhoPerp, rhoMPerp);
      if (jetMinusBkg.pt() < minJetPt)
        continue;

      // Apply area cut if required
      double normalizedJetArea = jet.area() / (PI * rJet * rJet);
      if (applyAreaCut && normalizedJetArea > maxNormalizedJetArea)
        continue;
      isAtLeastOneJetSelected = true;

      // Perpendicular cones
      double coneRadius = std::sqrt(jet.area() / PI);
      TVector3 jetAxis(jet.px(), jet.py(), jet.pz());
      TVector3 ueAxis1(0, 0, 0), ueAxis2(0, 0, 0);
      getPerpendicularDirections(jetAxis, ueAxis1, ueAxis2);
      if (ueAxis1.Mag() == 0 || ueAxis2.Mag() == 0) {
        continue;
      }

      // Get jet constituents
      std::vector<fastjet::PseudoJet> jetConstituents = jet.constituents();

      // Particle counters
      std::vector<int> nAntiprotonJet(nBins, 0);
      std::vector<int> nAntideuteronJet(nBins, 0);
      int nTotProtonJet(0);
      int nTotDeuteronJet(0);
      int nTotAntiprotonJet(0);
      int nTotAntideuteronJet(0);

      // Loop over jet constituents
      for (const auto& particle : jetConstituents) {

        // Get corresponding track and apply track selection criteria
        auto const& track = tracks.iteratorAt(particle.user_index());
        if (!passedTrackSelection(track))
          continue;

        // Apply DCA selections
        if (std::fabs(track.dcaXY()) > maxDcaxy || std::fabs(track.dcaZ()) > maxDcaz)
          continue;

        // Particle identification using the ITS cluster size
        bool passedItsPidProt(true), passedItsPidDeut(true);
        double nSigmaITSprot = static_cast<double>(itsResponse.nSigmaITS<o2::track::PID::Proton>(track));
        double nSigmaITSdeut = static_cast<double>(itsResponse.nSigmaITS<o2::track::PID::Deuteron>(track));

        if (applyItsPid && track.pt() < ptMaxItsPidProt && (nSigmaITSprot < nSigmaItsMin || nSigmaITSprot > nSigmaItsMax)) {
          passedItsPidProt = false;
        }
        if (applyItsPid && track.pt() < ptMaxItsPidDeut && (nSigmaITSdeut < nSigmaItsMin || nSigmaITSdeut > nSigmaItsMax)) {
          passedItsPidDeut = false;
        }

        // Kinematic range selection
        if (isProton(track) && passedItsPidProt) {
          if (track.pt() < ptOverAbins[0] || track.pt() >= ptOverAbins[nBins]) {
            continue;
          }
        } else if (isDeuteron(track) && passedItsPidDeut) {
          double ptPerNucleon = 0.5 * track.pt();
          if (ptPerNucleon < ptOverAbins[0] || ptPerNucleon >= ptOverAbins[nBins]) {
            continue;
          }
        } else {
          continue;
        }

        // (Anti)protons
        if (isProton(track) && passedItsPidProt) {
          if (track.sign() > 0) {
            nTotProtonJet++;
          } else if (track.sign() < 0) {
            nTotAntiprotonJet++;
            int ibin = findBin(ptOverAbins, track.pt());
            nAntiprotonJet[ibin]++;
          }
        }

        // (Anti)deuterons
        if (isDeuteron(track) && passedItsPidDeut) {
          const double ptPerNucleon = 0.5 * track.pt();

          if (track.sign() > 0) {
            nTotDeuteronJet++;
          } else if (track.sign() < 0) {
            nTotAntideuteronJet++;
            int ibin = findBin(ptOverAbins, ptPerNucleon);
            nAntideuteronJet[ibin]++;
          }
        }
      } // end of loop over constituents

      // Fill correlation histograms
      int netProtonJet = nTotProtonJet - nTotAntiprotonJet;
      int netDeuteronJet = nTotDeuteronJet - nTotAntideuteronJet;
      registryCorr.fill(HIST("rho_jet"), nTotAntideuteronJet, nTotAntiprotonJet, multiplicity);
      registryCorr.fill(HIST("rho_netP_netD_jet"), netDeuteronJet, netProtonJet);

      // Fill efficiency histograms
      for (int i = 0; i < nBins; i++) {
        double ptAcenteri = 0.5 * (ptOverAbins[i] + ptOverAbins[i + 1]);

        registryCorr.fill(HIST("q1d_jet"), nAntideuteronJet[i], ptAcenteri, multiplicity);
        registryCorr.fill(HIST("q1p_jet"), nAntiprotonJet[i], ptAcenteri, multiplicity);
        for (int j = 0; j < nBins; j++) {
          double ptAcenterj = 0.5 * (ptOverAbins[j] + ptOverAbins[j + 1]);
          registryCorr.fill(HIST("q1d_square_jet"), ptAcenteri, ptAcenterj, nAntideuteronJet[i] * nAntideuteronJet[j], multiplicity);
          registryCorr.fill(HIST("q1p_square_jet"), ptAcenteri, ptAcenterj, nAntiprotonJet[i] * nAntiprotonJet[j], multiplicity);
          registryCorr.fill(HIST("q1d_q1p_jet"), ptAcenteri, ptAcenterj, nAntideuteronJet[i] * nAntiprotonJet[j], multiplicity);
        }
      }

      // Particle counters
      std::vector<int> nAntiprotonUE(nBins, 0);
      std::vector<int> nAntideuteronUE(nBins, 0);
      int nTotProtonUE(0);
      int nTotDeuteronUE(0);
      int nTotAntiprotonUE(0);
      int nTotAntideuteronUE(0);

      // Loop over tracks in the underlying event
      for (auto const& track : tracks) {

        // Get corresponding track and apply track selection criteria
        if (!passedTrackSelection(track))
          continue;

        // Apply DCA selections
        if (std::fabs(track.dcaXY()) > maxDcaxy || std::fabs(track.dcaZ()) > maxDcaz)
          continue;

        // Calculate the angular distance between the track and underlying event axes in eta-phi space
        double deltaEtaUe1 = track.eta() - ueAxis1.Eta();
        double deltaPhiUe1 = getDeltaPhi(track.phi(), ueAxis1.Phi());
        double deltaRUe1 = std::sqrt(deltaEtaUe1 * deltaEtaUe1 + deltaPhiUe1 * deltaPhiUe1);
        double deltaEtaUe2 = track.eta() - ueAxis2.Eta();
        double deltaPhiUe2 = getDeltaPhi(track.phi(), ueAxis2.Phi());
        double deltaRUe2 = std::sqrt(deltaEtaUe2 * deltaEtaUe2 + deltaPhiUe2 * deltaPhiUe2);

        // Determine the maximum allowed distance from UE axes for particle selection
        double maxConeRadius = coneRadius;
        if (applyAreaCut) {
          maxConeRadius = std::sqrt(maxNormalizedJetArea) * rJet;
        }

        // Reject tracks that lie outside the maxConeRadius from both UE axes
        if (deltaRUe1 > maxConeRadius && deltaRUe2 > maxConeRadius)
          continue;

        // Particle identification using the ITS cluster size
        bool passedItsPidProt(true), passedItsPidDeut(true);
        double nSigmaITSprot = static_cast<double>(itsResponse.nSigmaITS<o2::track::PID::Proton>(track));
        double nSigmaITSdeut = static_cast<double>(itsResponse.nSigmaITS<o2::track::PID::Deuteron>(track));

        if (applyItsPid && track.pt() < ptMaxItsPidProt && (nSigmaITSprot < nSigmaItsMin || nSigmaITSprot > nSigmaItsMax)) {
          passedItsPidProt = false;
        }
        if (applyItsPid && track.pt() < ptMaxItsPidDeut && (nSigmaITSdeut < nSigmaItsMin || nSigmaITSdeut > nSigmaItsMax)) {
          passedItsPidDeut = false;
        }

        // Kinematic range selection
        if (isProton(track) && passedItsPidProt) {
          if (track.pt() < ptOverAbins[0] || track.pt() >= ptOverAbins[nBins]) {
            continue;
          }
        } else if (isDeuteron(track) && passedItsPidDeut) {
          double ptPerNucleon = 0.5 * track.pt();
          if (ptPerNucleon < ptOverAbins[0] || ptPerNucleon >= ptOverAbins[nBins]) {
            continue;
          }
        } else {
          continue;
        }

        // (Anti)protons
        if (isProton(track) && passedItsPidProt) {
          if (track.sign() > 0) {
            nTotProtonUE++;
          } else if (track.sign() < 0) {
            nTotAntiprotonUE++;
            int ibin = findBin(ptOverAbins, track.pt());
            nAntiprotonUE[ibin]++;
          }
        }

        // (Anti)deuterons
        if (isDeuteron(track) && passedItsPidDeut) {
          const double ptPerNucleon = 0.5 * track.pt();

          if (track.sign() > 0) {
            nTotDeuteronUE++;
          } else if (track.sign() < 0) {
            nTotAntideuteronUE++;
            int ibin = findBin(ptOverAbins, ptPerNucleon);
            nAntideuteronUE[ibin]++;
          }
        }
      }

      // Fill correlation histograms
      int netProtonUE = nTotProtonUE - nTotAntiprotonUE;
      int netDeuteronUE = nTotDeuteronUE - nTotAntideuteronUE;
      registryCorr.fill(HIST("rho_ue"), nTotAntideuteronUE, nTotAntiprotonUE, multiplicity);
      registryCorr.fill(HIST("rho_netP_netD_ue"), netDeuteronUE, netProtonUE);

      // Fill efficiency histograms
      for (int i = 0; i < nBins; i++) {
        double ptAcenteri = 0.5 * (ptOverAbins[i] + ptOverAbins[i + 1]);

        registryCorr.fill(HIST("q1d_ue"), nAntideuteronUE[i], ptAcenteri, multiplicity);
        registryCorr.fill(HIST("q1p_ue"), nAntiprotonUE[i], ptAcenteri, multiplicity);
        for (int j = 0; j < nBins; j++) {
          double ptAcenterj = 0.5 * (ptOverAbins[j] + ptOverAbins[j + 1]);
          registryCorr.fill(HIST("q1d_square_ue"), ptAcenteri, ptAcenterj, nAntideuteronUE[i] * nAntideuteronUE[j], multiplicity);
          registryCorr.fill(HIST("q1p_square_ue"), ptAcenteri, ptAcenterj, nAntiprotonUE[i] * nAntiprotonUE[j], multiplicity);
          registryCorr.fill(HIST("q1d_q1p_ue"), ptAcenteri, ptAcenterj, nAntideuteronUE[i] * nAntiprotonUE[j], multiplicity);
        }
      }
      // Fill event counter vs centrality (ue region)
      registryCorr.fill(HIST("eventCounter_centrality_ue"), multiplicity);
    }
    // Event counter: events with at least one jet selected
    if (isAtLeastOneJetSelected) {
      registryCorr.fill(HIST("eventCounter"), 9.5);
      registryCorr.fill(HIST("eventCounter_centrality_jet"), multiplicity);
    }
    */
  }
  PROCESS_SWITCH(AntinucleiInJets, processCorr, "Process Correlation analysis", false);
    
    // Process correlation analysis with systematic variations of analysis parameters
    void processCorrSyst(SelectedCollisions::iterator const& collision, AntiNucleiTracks const& tracks)
    {
      // cut settings (from processSystData)
      static std::vector<double> maxDcaxySyst = {
        0.071, 0.060, 0.066, 0.031, 0.052, 0.078, 0.045, 0.064, 0.036, 0.074,
        0.079, 0.043, 0.067, 0.059, 0.032, 0.070, 0.048, 0.077, 0.062, 0.034,
        0.057, 0.055, 0.073, 0.038, 0.050, 0.075, 0.041, 0.061, 0.033, 0.069,
        0.035, 0.044, 0.076, 0.049, 0.037, 0.054, 0.072, 0.046, 0.058, 0.040,
        0.068, 0.042, 0.056, 0.039, 0.047, 0.065, 0.051, 0.053, 0.063, 0.030};

      static std::vector<double> maxDcazSyst = {
        0.064, 0.047, 0.032, 0.076, 0.039, 0.058, 0.043, 0.069, 0.050, 0.035,
        0.074, 0.061, 0.045, 0.033, 0.068, 0.055, 0.037, 0.071, 0.042, 0.053,
        0.077, 0.038, 0.065, 0.049, 0.036, 0.059, 0.044, 0.067, 0.041, 0.034,
        0.073, 0.052, 0.040, 0.063, 0.046, 0.031, 0.070, 0.054, 0.037, 0.062,
        0.048, 0.035, 0.075, 0.051, 0.039, 0.066, 0.043, 0.060, 0.032, 0.056};

      static std::vector<double> nSigmaItsMinSyst = {
        -2.9, -2.8, -3.0, -3.4, -2.7, -3.3, -3.0, -3.1, -3.4, -3.1,
        -3.0, -2.8, -3.2, -2.6, -2.7, -3.4, -2.9, -3.0, -3.0, -2.7,
        -2.9, -3.3, -3.0, -3.1, -3.2, -3.0, -2.9, -2.7, -3.3, -3.0,
        -2.8, -3.3, -2.7, -3.3, -2.8, -3.4, -2.8, -3.4, -2.9, -3.1,
        -3.2, -2.6, -3.1, -2.9, -3.1, -2.8, -2.9, -3.3, -3.0, -2.8};

      static std::vector<double> nSigmaItsMaxSyst = {
        2.9, 2.8, 3.0, 3.4, 2.7, 3.3, 3.0, 3.1, 3.4, 3.1,
        3.0, 2.8, 3.2, 2.6, 2.7, 3.4, 2.9, 3.0, 3.0, 2.7,
        2.9, 3.3, 3.0, 3.1, 3.2, 3.0, 2.9, 2.7, 3.3, 3.0,
        2.8, 3.3, 2.7, 3.3, 2.8, 3.4, 2.8, 3.4, 2.9, 3.1,
        3.2, 2.6, 3.1, 2.9, 3.1, 2.8, 2.9, 3.3, 3.0, 2.8};

      static std::vector<double> minNsigmaTpcSyst = {
        -3.2, -2.9, -3.1, -2.9, -3.5, -2.6, -3.3, -3.0, -3.5, -2.7,
        -3.0, -2.6, -3.3, -3.4, -2.8, -3.1, -2.6, -3.2, -3.1, -2.8,
        -3.4, -2.7, -3.4, -2.9, -3.0, -2.5, -3.3, -2.8, -3.1, -2.7,
        -3.4, -2.8, -3.3, -2.6, -3.1, -2.5, -3.4, -3.0, -3.2, -2.6,
        -3.4, -2.8, -3.1, -2.6, -3.3, -2.7, -3.2, -2.7, -3.4, -2.9};

      static std::vector<double> maxNsigmaTpcSyst = {
        3.2, 2.9, 3.1, 2.9, 3.5, 2.6, 3.3, 3.0, 3.5, 2.7,
        3.0, 2.6, 3.3, 3.4, 2.8, 3.1, 2.6, 3.2, 3.1, 2.8,
        3.4, 2.7, 3.4, 2.9, 3.0, 2.5, 3.3, 2.8, 3.1, 2.7,
        3.4, 2.8, 3.3, 2.6, 3.1, 2.5, 3.4, 3.0, 3.2, 2.6,
        3.4, 2.8, 3.1, 2.6, 3.3, 2.7, 3.2, 2.7, 3.4, 2.9};

      static std::vector<double> minNsigmaTofSyst = {
        -3.2, -2.9, -3.1, -2.9, -3.5, -2.6, -3.3, -3.0, -3.5, -2.7,
        -3.0, -2.6, -3.3, -3.4, -2.8, -3.1, -2.6, -3.2, -3.1, -2.8,
        -3.4, -2.7, -3.4, -2.9, -3.0, -2.5, -3.3, -2.8, -3.1, -2.7,
        -3.4, -2.8, -3.3, -2.6, -3.1, -2.5, -3.4, -3.0, -3.2, -2.6,
        -3.4, -2.8, -3.1, -2.6, -3.3, -2.7, -3.2, -2.7, -3.4, -2.9};

      static std::vector<double> maxNsigmaTofSyst = {
        3.9, 3.6, 3.8, 3.2, 3.2, 3.5, 3.1, 3.8, 3.5, 3.4,
        3.9, 3.8, 3.7, 3.0, 3.6, 3.1, 3.7, 3.4, 4.0, 3.0,
        3.7, 3.3, 3.9, 3.1, 3.3, 3.5, 3.6, 3.2, 3.5, 3.3,
        3.9, 3.0, 3.4, 3.2, 3.1, 3.9, 3.6, 3.1, 3.2, 4.0,
        3.1, 3.7, 3.6, 3.1, 3.3, 3.5, 3.3, 3.4, 3.1, 3.8};

      // Event counter: before event selection
      registryCorr.fill(HIST("eventCounter_syst"), 0.5);

      // Apply standard event selection (Same as processCorr)
      if (!collision.sel8() || std::fabs(collision.posZ()) > zVtx)
        return;

      registryCorr.fill(HIST("eventCounter_syst"), 1.5);

      if (rejectITSROFBorder && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder))
        return;
      registryCorr.fill(HIST("eventCounter_syst"), 2.5);

      if (rejectTFBorder && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder))
        return;
      registryCorr.fill(HIST("eventCounter_syst"), 3.5);

      if (requireVtxITSTPC && !collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC))
        return;
      registryCorr.fill(HIST("eventCounter_syst"), 4.5);

      if (rejectSameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup))
        return;
      registryCorr.fill(HIST("eventCounter_syst"), 5.5);

      if (requireIsGoodZvtxFT0VsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV))
        return;
      registryCorr.fill(HIST("eventCounter_syst"), 6.5);

      if (requireIsVertexTOFmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTOFmatched))
        return;
      registryCorr.fill(HIST("eventCounter_syst"), 7.5);

      const float multiplicity = collision.centFT0M();
        
        
      // Bins setup
      static const std::vector<double> ptOverAbins = {0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
      const int nBins = ptOverAbins.size() - 1;

      // Loop over systematic variations
      for (int isyst = 0; isyst < nSyst; isyst++) {

        // Fill event counter for this systematic
        registryCorr.fill(HIST("eventCounter_centrality_fullEvent_syst"), multiplicity, static_cast<double>(isyst));

        // Particle counters for this specific cut setting
        std::vector<int> nAntiprotonFullEvent(nBins, 0);
        std::vector<int> nAntideuteronFullEvent(nBins, 0);
        int nTotProtonFullEvent(0);
        int nTotDeuteronFullEvent(0);
        int nTotAntiprotonFullEvent(0);
        int nTotAntideuteronFullEvent(0);

        // Loop over tracks
        for (auto const& track : tracks) {

            // Apply track selection (from processSystData)
          if (!passedTrackSelectionSyst(track, isyst))
            continue;

            // Apply DCA selections (from processSystData)
          if (std::fabs(track.dcaXY()) > maxDcaxySyst[isyst] || std::fabs(track.dcaZ()) > maxDcazSyst[isyst])
            continue;

          // Particle identification using the ITS cluster size (vary also PID ITS) (from processSystData)
          bool passedItsPidProt(true), passedItsPidDeut(true);
          double nSigmaITSprot = static_cast<double>(itsResponse.nSigmaITS<o2::track::PID::Proton>(track));
          double nSigmaITSdeut = static_cast<double>(itsResponse.nSigmaITS<o2::track::PID::Deuteron>(track));

          if (applyItsPid && track.pt() < ptMaxItsPidProt && (nSigmaITSprot < nSigmaItsMinSyst[isyst] || nSigmaITSprot > nSigmaItsMaxSyst[isyst])) {
            passedItsPidProt = false;
          }
          if (applyItsPid && track.pt() < ptMaxItsPidDeut && (nSigmaITSdeut < nSigmaItsMinSyst[isyst] || nSigmaITSdeut > nSigmaItsMaxSyst[isyst])) {
            passedItsPidDeut = false;
          }

          // Check Identity (Proton/Deuteron) using TPC/TOF
          bool isProt = isProtonSyst(track, minNsigmaTpcSyst[isyst], maxNsigmaTpcSyst[isyst], minNsigmaTofSyst[isyst], maxNsigmaTofSyst[isyst]);
          bool isDeut = isDeuteronSyst(track, minNsigmaTpcSyst[isyst], maxNsigmaTpcSyst[isyst], minNsigmaTofSyst[isyst], maxNsigmaTofSyst[isyst]);

          // Kinematic range selection & counting
          // (Anti)protons
          if (isProt && passedItsPidProt) {
            if (track.pt() >= ptOverAbins[0] && track.pt() < ptOverAbins[nBins]) {
              if (track.sign() > 0) {
                nTotProtonFullEvent++;
              } else if (track.sign() < 0) {
                nTotAntiprotonFullEvent++;
                int ibin = findBin(ptOverAbins, track.pt());
                nAntiprotonFullEvent[ibin]++;
              }
            }
          }

          // (Anti)deuterons
          if (isDeut && passedItsPidDeut) {
            double ptPerNucleon = 0.5 * track.pt();
            if (ptPerNucleon >= ptOverAbins[0] && ptPerNucleon < ptOverAbins[nBins]) {
              if (track.sign() > 0) {
                nTotDeuteronFullEvent++;
              } else if (track.sign() < 0) {
                nTotAntideuteronFullEvent++;
                int ibin = findBin(ptOverAbins, ptPerNucleon);
                nAntideuteronFullEvent[ibin]++;
              }
            }
          }
        }

        // Fill Correlation Histograms for systematic variations
        int netProtonFullEvent = nTotProtonFullEvent - nTotAntiprotonFullEvent;
        int netDeuteronFullEvent = nTotDeuteronFullEvent - nTotAntideuteronFullEvent;

        registryCorr.fill(HIST("rho_fullEvent_syst"), nTotAntideuteronFullEvent, nTotAntiprotonFullEvent, multiplicity, static_cast<double>(isyst));
        registryCorr.fill(HIST("rho_netP_netD_fullEvent_syst"), netDeuteronFullEvent, netProtonFullEvent, static_cast<double>(isyst));

        // Fill efficiency histograms for systematic variations
        for (int i = 0; i < nBins; i++) {
          double ptAcenteri = 0.5 * (ptOverAbins[i] + ptOverAbins[i + 1]);

          registryCorr.fill(HIST("q1d_fullEvent_syst"), nAntideuteronFullEvent[i], ptAcenteri, multiplicity, static_cast<double>(isyst));
          registryCorr.fill(HIST("q1p_fullEvent_syst"), nAntiprotonFullEvent[i], ptAcenteri, multiplicity, static_cast<double>(isyst));
          for (int j = 0; j < nBins; j++) {
            double ptAcenterj = 0.5 * (ptOverAbins[j] + ptOverAbins[j + 1]);
            registryCorr.fill(HIST("q1d_square_fullEvent_syst"), ptAcenteri, ptAcenterj, (nAntideuteronFullEvent[i] * nAntideuteronFullEvent[j]), multiplicity, static_cast<double>(isyst));
            registryCorr.fill(HIST("q1p_square_fullEvent_syst"), ptAcenteri, ptAcenterj, (nAntiprotonFullEvent[i] * nAntiprotonFullEvent[j]), multiplicity, static_cast<double>(isyst));
            registryCorr.fill(HIST("q1d_q1p_fullEvent_syst"), ptAcenteri, ptAcenterj, (nAntideuteronFullEvent[i] * nAntiprotonFullEvent[j]), multiplicity, static_cast<double>(isyst));
          }
        }

      }
    }
    PROCESS_SWITCH(AntinucleiInJets, processCorrSyst, "Process Correlation systematics", false);
    
    // Process coalescence
    void processCoalescence(GenCollisionsMc const& collisions, aod::McParticles const& mcParticles)
    {
      // Deuteron Mass and minimum pt
      double massDeut = o2::constants::physics::MassDeuteron;
      static constexpr double MinPtParticle = 0.1;

      // Define per-event particle containers
      std::vector<ReducedParticle> chargedParticles;
      std::vector<ReducedParticle> protons;
      std::vector<ReducedParticle> neutrons;
      std::vector<fastjet::PseudoJet> fjParticles;

      // Jet and area definitions
      fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, rJet);
      fastjet::AreaDefinition areaDef(fastjet::active_area, fastjet::GhostedAreaSpec(1.0));

      // Loop over all simulated collisions
      for (const auto& collision : collisions) {

        // Clear containers at the start of the event loop
        chargedParticles.clear();
        protons.clear();
        neutrons.clear();
        fjParticles.clear();

        // Event counter: before event selection
        registryMC.fill(HIST("genEventsCoalescence"), 0.5);

        // Apply event selection: require vertex position to be within the allowed z range
        if (std::fabs(collision.posZ()) > zVtx)
          continue;

        // Event counter: after event selection
        registryMC.fill(HIST("genEventsCoalescence"), 1.5);

        // Get particles in this MC collision
        const auto mcParticlesThisMcColl = mcParticles.sliceBy(mcParticlesPerMcCollision, collision.globalIndex());

        // Loop over MC particles
        for (const auto& particle : mcParticlesThisMcColl) {

          // Monte Carlo index
          int mcId = particle.globalIndex();

          // Store Protons
          if (particle.isPhysicalPrimary() && std::abs(particle.pdgCode()) == PDG_t::kProton) {
            protons.push_back({particle.px(), particle.py(), particle.pz(), particle.pdgCode(), mcId, false});
          }

          // Store Neutrons
          if (particle.isPhysicalPrimary() && std::abs(particle.pdgCode()) == PDG_t::kNeutron) {
            neutrons.push_back({particle.px(), particle.py(), particle.pz(), particle.pdgCode(), mcId, false});
          }

          // Select physical primary particles or HF decay products
          if (!isPhysicalPrimaryOrFromHF(particle, mcParticles))
            continue;

          // Select particles within acceptance
          if (particle.eta() < minEta || particle.eta() > maxEta || particle.pt() < MinPtParticle)
            continue;
          chargedParticles.push_back({particle.px(), particle.py(), particle.pz(), particle.pdgCode(), mcId, false});
        }

        // Reject empty events
        if (chargedParticles.empty())
          continue;
        registryMC.fill(HIST("genEventsCoalescence"), 2.5);

        // Build deuterons
        for (int ip = 0; ip < static_cast<int>(protons.size()); ip++) {
          auto& proton = protons[ip];
          if (proton.used)
            continue;

          for (int in = 0; in < static_cast<int>(neutrons.size()); in++) {
            auto& neutron = neutrons[in];
            if (neutron.used)
              continue;

            if (passDeuteronCoalescence(proton, neutron, coalescenceMomentum, mRand)) {

              int sign = (proton.pdgCode > 0) ? +1 : -1;
              int deuteronPdg = sign * o2::constants::physics::Pdg::kDeuteron;

              double pxDeut = proton.px + neutron.px;
              double pyDeut = proton.py + neutron.py;
              double pzDeut = proton.pz + neutron.pz;
              double energyDeut = std::sqrt(pxDeut * pxDeut + pyDeut * pyDeut + pzDeut * pzDeut + massDeut * massDeut);
              LorentzVector pd(pxDeut, pyDeut, pzDeut, energyDeut);
              if (pd.Eta() < minEta || pd.Eta() > maxEta || pd.Pt() < MinPtParticle)
                continue;

              // Store Deuteron
              chargedParticles.push_back({pxDeut, pyDeut, pzDeut, deuteronPdg, proton.mcIndex, false});

              neutron.used = true;
              proton.used = true;
              break;
            }
          }
        }
          
          for (const auto& part : chargedParticles) {
              if (part.used)
                  continue;
              
              // Fill histograms for Full Event
              if (part.pdgCode == PDG_t::kProtonBar) {
                  registryMC.fill(HIST("antiproton_coal_fullEvent"), part.pt());
              }
              if (part.pdgCode == -o2::constants::physics::Pdg::kDeuteron) {
                  registryMC.fill(HIST("antideuteron_coal_fullEvent"), part.pt());
              }
          }

        // Fill particle array to feed to Fastjet
        for (const auto& part : chargedParticles) {

          if (part.used)
            continue;

          double energy = std::sqrt(part.px * part.px + part.py * part.py + part.pz * part.pz + MassPionCharged * MassPionCharged);
          fastjet::PseudoJet fourMomentum(part.px, part.py, part.pz, energy);
          fourMomentum.set_user_index(part.pdgCode);
          fjParticles.emplace_back(fourMomentum);
        }

        // Reject empty events
        if (fjParticles.empty())
          continue;
        registryMC.fill(HIST("genEventsCoalescence"), 3.5);

        // Cluster MC particles into jets using anti-kt algorithm
        fastjet::ClusterSequenceArea cs(fjParticles, jetDef, areaDef);
        std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cs.inclusive_jets());
        if (jets.empty())
          continue;
        registryMC.fill(HIST("genEventsCoalescence"), 4.5);
        auto [rhoPerp, rhoMPerp] = jetutilities::estimateRhoPerpCone(fjParticles, jets[0], rJet);

        // Loop over clustered jets
        bool isAtLeastOneJetSelected = false;
        for (const auto& jet : jets) {

          // Jet must be fully contained in the acceptance
          if ((std::fabs(jet.eta()) + rJet) > (maxEta - deltaEtaEdge))
            continue;

          // Jet pt must be larger than threshold
          auto jetForSub = jet;
          fastjet::PseudoJet jetMinusBkg = backgroundSub.doRhoAreaSub(jetForSub, rhoPerp, rhoMPerp);
          if (jetMinusBkg.pt() < minJetPt)
            continue;

          // Apply area cut if required
          double normalizedJetArea = jet.area() / (PI * rJet * rJet);
          if (applyAreaCut && normalizedJetArea > maxNormalizedJetArea)
            continue;
          isAtLeastOneJetSelected = true;

          // Analyze jet constituents
          std::vector<fastjet::PseudoJet> jetConstituents = jet.constituents();
          for (const auto& particle : jetConstituents) {
            if (particle.user_index() == PDG_t::kProtonBar)
              registryMC.fill(HIST("antiproton_coal_jet"), particle.pt());
            if (particle.user_index() == -o2::constants::physics::Pdg::kDeuteron)
              registryMC.fill(HIST("antideuteron_coal_jet"), particle.pt());
          }

          // Set up two perpendicular cone axes for underlying event estimation
          TVector3 jetAxis(jet.px(), jet.py(), jet.pz());
          double coneRadius = std::sqrt(jet.area() / PI);
          TVector3 ueAxis1(0, 0, 0), ueAxis2(0, 0, 0);
          getPerpendicularDirections(jetAxis, ueAxis1, ueAxis2);
          if (ueAxis1.Mag() == 0 || ueAxis2.Mag() == 0) {
            continue;
          }

          // Loop over MC particles to analyze underlying event region
          for (const auto& chParticle : chargedParticles) {

            // Skip used particles
            if (chParticle.used)
              continue;

            // Compute distance of particle from both perpendicular cone axes
            double deltaEtaUe1 = chParticle.eta() - ueAxis1.Eta();
            double deltaPhiUe1 = getDeltaPhi(chParticle.phi(), ueAxis1.Phi());
            double deltaRUe1 = std::sqrt(deltaEtaUe1 * deltaEtaUe1 + deltaPhiUe1 * deltaPhiUe1);
            double deltaEtaUe2 = chParticle.eta() - ueAxis2.Eta();
            double deltaPhiUe2 = getDeltaPhi(chParticle.phi(), ueAxis2.Phi());
            double deltaRUe2 = std::sqrt(deltaEtaUe2 * deltaEtaUe2 + deltaPhiUe2 * deltaPhiUe2);

            // Determine the maximum allowed distance from UE axes for particle selection
            double maxConeRadius = coneRadius;
            if (applyAreaCut) {
              maxConeRadius = std::sqrt(maxNormalizedJetArea) * rJet;
            }

            // Reject tracks that lie outside the maxConeRadius from both UE axes
            if (deltaRUe1 > maxConeRadius && deltaRUe2 > maxConeRadius)
              continue;

            // Fill histograms for UE
            if (chParticle.pdgCode == PDG_t::kProtonBar)
              registryMC.fill(HIST("antiproton_coal_ue"), chParticle.pt());
            if (chParticle.pdgCode == -o2::constants::physics::Pdg::kDeuteron)
              registryMC.fill(HIST("antideuteron_coal_ue"), chParticle.pt());
          }
        }
        if (isAtLeastOneJetSelected) {
          registryMC.fill(HIST("genEventsCoalescence"), 5.5);
        }
      }
    }
    PROCESS_SWITCH(AntinucleiInJets, processCoalescence, "process coalescence", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<AntinucleiInJets>(cfgc)};
}
