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
/// \author Alberto Caliva (alberto.caliva@cern.ch), Chiara Pinto (chiara.pinto@cern.ch)
/// \since February 13, 2025

#include "PWGJE/Core/JetBkgSubUtils.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "EventFiltering/Zorro.h"
#include "EventFiltering/ZorroSummary.h"

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
#include <TList.h>
#include <TPDGCode.h>
#include <TRandom.h>
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

struct AntinucleiInJets {

  // Histogram registries for data, MC, quality control, multiplicity and correlations
  HistogramRegistry registryData{"registryData", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry registryMC{"registryMC", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry registryQC{"registryQC", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry registryMult{"registryMult", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry registryCorr{"registryCorr", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

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
  Configurable<bool> applyAreaCut{"applyAreaCut", false, "apply area cut"};
  Configurable<double> maxNormalizedJetArea{"maxNormalizedJetArea", 0.6, "area cut"};
  Configurable<double> deltaEtaEdge{"deltaEtaEdge", 0.05, "eta gap from the edge"};
  Configurable<int> nSyst{"nSyst", 50, "number of systematic variations"};

  // Track quality, kinematic, and PID selection parameters
  Configurable<bool> requirePvContributor{"requirePvContributor", false, "require that the track is a PV contributor"};
  Configurable<bool> applyItsPid{"applyItsPid", true, "apply ITS PID"};
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
  Configurable<bool> setMCDefaultItsParams{"setMCDefaultItsParams", false, "set MC default parameters"};
  Configurable<bool> cfgCompensatePIDinTracking{"cfgCompensatePIDinTracking", false, "If true, divide tpcInnerParam by the electric charge"};
  Configurable<std::array<double, 5>> cfgBetheBlochParams{"cfgBetheBlochParams", {0.6539, 1.591, 0.8225, 2.363, 0.09}, "TPC Bethe-Bloch parameterisation for He3"};

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

      // Jet effective area over piR^2
      registryData.add("jetEffectiveAreaOverPiR2", "jet effective area / piR^2", HistType::kTH1F, {{2000, 0, 2, "Area/#piR^{2}"}});

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

      // Generated spectra of antiprotons
      registryMC.add("antiproton_gen_jet", "antiproton_gen_jet", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("antiproton_gen_ue", "antiproton_gen_ue", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});

      // Normalization histogram
      registryMC.add("antiproton_deltay_deltaphi_jet", "antiproton_deltay_deltaphi_jet", HistType::kTH2F, {{2000, -1.0, 1.0, "#Delta#it{y}"}, {2000, 0.0, 2.0, "#Delta#phi"}});
      registryMC.add("antiproton_deltay_deltaphi_ue", "antiproton_deltay_deltaphi_ue", HistType::kTH2F, {{2000, -1.0, 1.0, "#Delta#it{y}"}, {2000, 0.0, 2.0, "#Delta#phi"}});
    }

    // Reconstructed antiproton spectra in jets and UE (MC-matched) with TPC/TOF PID
    if (doprocessJetsMCrec) {

      // Event counter
      registryMC.add("recEvents", "number of reconstructed events in mc", HistType::kTH1F, {{20, 0, 20, "counter"}});
      registryMC.add("recJets", "number of reconstructed jets", HistType::kTH1F, {{10, 0, 10, "counter"}});

      // Reconstructed spectra of antiprotons
      registryMC.add("antiproton_rec_tpc_jet", "antiproton_rec_tpc_jet", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("antiproton_rec_tof_jet", "antiproton_rec_tof_jet", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("antiproton_rec_tpc_ue", "antiproton_rec_tpc_ue", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("antiproton_rec_tof_ue", "antiproton_rec_tof_ue", HistType::kTH1F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}});

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

      // nsigmaITS for antiproton candidates
      registryMC.add("antiproton_nsigma_its_mc", "antiproton_nsigma_its_mc", HistType::kTH2F, {{nbins, min, max, "#it{p}_{T} (GeV/#it{c})"}, {400, -20.0, 20.0, "n#sigma_{ITS}"}});

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
      registryMC.add("protonBar", "protonBar", HistType::kTH1F, {{5000, 0, 5, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("lambdaBar", "lambdaBar", HistType::kTH1F, {{5000, 0, 5, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("xiBar", "xiBar", HistType::kTH1F, {{5000, 0, 5, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("omegaBar", "omegaBar", HistType::kTH1F, {{5000, 0, 5, "#it{p}_{T} (GeV/#it{c})"}});
      registryMC.add("sigmaBar", "sigmaBar", HistType::kTH1F, {{5000, 0, 5, "#it{p}_{T} (GeV/#it{c})"}});
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

      // Event counter
      registryCorr.add("eventCounter", "number of events", HistType::kTH1F, {{20, 0, 20, "counter"}});

      // Correlation histograms: antiproton vs. antideuteron number vs. event multiplicity
      registryCorr.add("rho_jet", "rho_jet", HistType::kTH3F, {nAntideuteronsAxis, nAntiprotonsAxis, multiplicityAxis});
      registryCorr.add("rho_ue", "rho_ue", HistType::kTH3F, {nAntideuteronsAxis, nAntiprotonsAxis, multiplicityAxis});
      registryCorr.add("rho_fullEvent", "rho_fullEvent", HistType::kTH3F, {nAntideuteronsAxis, nAntiprotonsAxis, multiplicityAxis});

      // Correlation histograms: net antiproton vs. net antideuteron numbers
      registryCorr.add("rho_netP_netD_jet", "rho_netP_netD_jet", HistType::kTH2F, {nAntideuteronsAxis, nAntiprotonsAxis});
      registryCorr.add("rho_netP_netD_ue", "rho_netP_netD_ue", HistType::kTH2F, {nAntideuteronsAxis, nAntiprotonsAxis});
      registryCorr.add("rho_netP_netD_fullEvent", "rho_netP_netD_fullEvent", HistType::kTH2F, {nAntideuteronsAxis, nAntiprotonsAxis});

      // Efficiency histograms jet
      registryCorr.add("q1d_jet", "q1d_jet", HistType::kTH2F, {nAntideuteronsAxis, ptPerNucleonAxis});
      registryCorr.add("q1p_jet", "q1p_jet", HistType::kTH2F, {nAntiprotonsAxis, ptPerNucleonAxis});
      registryCorr.add("q1d_square_jet", "q1d_square_jet", HistType::kTH3F, {ptPerNucleonAxis, ptPerNucleonAxis, nBarD2Axis});
      registryCorr.add("q1p_square_jet", "q1p_square_jet", HistType::kTH3F, {ptPerNucleonAxis, ptPerNucleonAxis, nBarP2Axis});
      registryCorr.add("q1d_q1p_jet", "q1d_q1p_jet", HistType::kTH3F, {ptPerNucleonAxis, ptPerNucleonAxis, nBarDnBarPAxis});

      // Efficiency histograms UE
      registryCorr.add("q1d_ue", "q1d_ue", HistType::kTH2F, {nAntideuteronsAxis, ptPerNucleonAxis});
      registryCorr.add("q1p_ue", "q1p_ue", HistType::kTH2F, {nAntiprotonsAxis, ptPerNucleonAxis});
      registryCorr.add("q1d_square_ue", "q1d_square_ue", HistType::kTH3F, {ptPerNucleonAxis, ptPerNucleonAxis, nBarD2Axis});
      registryCorr.add("q1p_square_ue", "q1p_square_ue", HistType::kTH3F, {ptPerNucleonAxis, ptPerNucleonAxis, nBarP2Axis});
      registryCorr.add("q1d_q1p_ue", "q1d_q1p_ue", HistType::kTH3F, {ptPerNucleonAxis, ptPerNucleonAxis, nBarDnBarPAxis});

      // Efficiency histograms full event
      registryCorr.add("q1d_fullEvent", "q1d_fullEvent", HistType::kTH2F, {nAntideuteronsAxis, ptPerNucleonAxis});
      registryCorr.add("q1p_fullEvent", "q1p_fullEvent", HistType::kTH2F, {nAntiprotonsAxis, ptPerNucleonAxis});
      registryCorr.add("q1d_square_fullEvent", "q1d_square_fullEvent", HistType::kTH3F, {ptPerNucleonAxis, ptPerNucleonAxis, nBarD2Axis});
      registryCorr.add("q1p_square_fullEvent", "q1p_square_fullEvent", HistType::kTH3F, {ptPerNucleonAxis, ptPerNucleonAxis, nBarP2Axis});
      registryCorr.add("q1d_q1p_fullEvent", "q1d_q1p_fullEvent", HistType::kTH3F, {ptPerNucleonAxis, ptPerNucleonAxis, nBarDnBarPAxis});
    }
  }

  // Compute two unit vectors perpendicular to p
  void getPerpendicularAxis(const TVector3& p, TVector3& u, double sign)
  {
    double px = p.X();
    double py = p.Y();
    double pz = p.Z();

    double px2 = px * px;
    double py2 = py * py;
    double pz2 = pz * pz;
    double pz4 = pz2 * pz2;

    // px and py are both zero
    if (px == 0 && py == 0) {
      u.SetXYZ(0, 0, 0);
      return;
    }

    // protection 1
    if (px == 0 && py != 0) {
      double ux = sign * std::sqrt(py2 - pz4 / py2);
      double uy = -pz2 / py;
      u.SetXYZ(ux, uy, pz);
      return;
    }

    // protection 2
    if (py == 0 && px != 0) {
      double ux = -pz2 / px;
      double uy = sign * std::sqrt(px2 - pz4 / px2);
      u.SetXYZ(ux, uy, pz);
      return;
    }

    // General case
    double a = px2 + py2;
    double b = 2.0 * px * pz2;
    double c = pz4 - py2 * py2 - px2 * py2;

    double delta = b * b - 4.0 * a * c;

    if (delta < 0 || a == 0) {
      LOGP(warn, "Invalid input in getPerpendicularAxis: delta = {}, a = {}", delta, a);
      u.SetXYZ(0, 0, 0);
      return;
    }

    double ux = (-b + sign * std::sqrt(delta)) / (2.0 * a);
    double uy = (-pz2 - px * ux) / py;
    u.SetXYZ(ux, uy, pz);
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

    static std::vector<double> minEtaSyst = {
      -0.804, -0.844, -0.751, -0.784, -0.819, -0.823, -0.768, -0.781, -0.845, -0.787,
      -0.758, -0.828, -0.776, -0.842, -0.808, -0.763, -0.849, -0.770, -0.799, -0.754,
      -0.825, -0.847, -0.806, -0.783, -0.796, -0.835, -0.777, -0.752, -0.838, -0.813,
      -0.785, -0.802, -0.795, -0.846, -0.780, -0.829, -0.817, -0.773, -0.765, -0.789,
      -0.800, -0.839, -0.758, -0.820, -0.833, -0.762, -0.792, -0.809, -0.827, -0.751};

    static std::vector<double> maxEtaSyst = {
      0.804, 0.844, 0.751, 0.784, 0.819, 0.823, 0.768, 0.781, 0.845, 0.787,
      0.758, 0.828, 0.776, 0.842, 0.808, 0.763, 0.849, 0.770, 0.799, 0.754,
      0.825, 0.847, 0.806, 0.783, 0.796, 0.835, 0.777, 0.752, 0.838, 0.813,
      0.785, 0.802, 0.795, 0.846, 0.780, 0.829, 0.817, 0.773, 0.765, 0.789,
      0.800, 0.839, 0.758, 0.820, 0.833, 0.762, 0.792, 0.809, 0.827, 0.751};

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
    if (track.eta() < minEtaSyst[isyst] || track.eta() > maxEtaSyst[isyst])
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
    static constexpr double kPtThreshold = 0.6;
    static constexpr double kNsigmaMax = 3.0;

    // PID variables and transverse momentum of the track
    const double nsigmaTPC = track.tpcNSigmaPr();
    const double nsigmaTOF = track.tofNSigmaPr();
    const double pt = track.pt();

    // Apply TPC PID cut
    if (std::abs(nsigmaTPC) > kNsigmaMax)
      return false;

    // Low-pt: TPC PID is sufficient
    if (pt < kPtThreshold)
      return true;

    // High-pt: require valid TOF match and pass TOF PID
    return (track.hasTOF() && std::abs(nsigmaTOF) < kNsigmaMax);
  }

  // Selection of (anti)deuterons
  template <typename DeuteronTrack>
  bool isDeuteron(const DeuteronTrack& track)
  {
    // Constants
    static constexpr double kPtThreshold = 1.0;
    static constexpr double kNsigmaMax = 3.0;

    // PID variables and transverse momentum of the track
    const double nsigmaTPC = track.tpcNSigmaDe();
    const double nsigmaTOF = track.tofNSigmaDe();
    const double pt = track.pt();

    // Apply TPC PID cut
    if (std::abs(nsigmaTPC) > kNsigmaMax)
      return false;

    // Low-pt: TPC PID is sufficient
    if (pt < kPtThreshold)
      return true;

    // High-pt: require valid TOF match and pass TOF PID
    return (track.hasTOF() && std::abs(nsigmaTOF) < kNsigmaMax);
  }

  // Process Data
  void processData(SelectedCollisions::iterator const& collision, AntiNucleiTracks const& tracks, aod::BCsWithTimestamps const&)
  {
    // Event counter: before event selection
    registryData.fill(HIST("number_of_events_data"), 0.5);

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
    auto [rhoPerp, rhoMPerp] = backgroundSub.estimateRhoPerpCone(fjParticles, jets);

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
      TVector3 ueAxis1(0, 0, 0);
      TVector3 ueAxis2(0, 0, 0);
      getPerpendicularAxis(jetAxis, ueAxis1, +1);
      getPerpendicularAxis(jetAxis, ueAxis2, -1);

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
            double sigmaHe3 = expectedSignalHe3 * cfgBetheBlochParams.value[4];
            double nSigmaTPCHe3Custom = (tpcSignal - expectedSignalHe3) / sigmaHe3;
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
            double sigmaHe3 = expectedSignalHe3 * cfgBetheBlochParams.value[4];
            double nSigmaTPCHe3Custom = (tpcSignal - expectedSignalHe3) / sigmaHe3;
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
    auto [rhoPerp, rhoMPerp] = backgroundSub.estimateRhoPerpCone(fjParticles, jets);

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
    auto [rhoPerp, rhoMPerp] = backgroundSub.estimateRhoPerpCone(fjParticles, jets);

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
      TVector3 ueAxis1(0, 0, 0);
      TVector3 ueAxis2(0, 0, 0);
      getPerpendicularAxis(jetAxis, ueAxis1, +1);
      getPerpendicularAxis(jetAxis, ueAxis2, -1);

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

  // Antinuclei reconstruction efficiency
  void processAntinucleiEfficiency(RecCollisionsMc const& collisions, AntiNucleiTracksMc const& mcTracks, aod::McParticles const& mcParticles)
  {
    // Loop over all simulated collision events
    for (const auto& collision : collisions) {

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

      // Loop over all generated Monte Carlo particles for the selected event
      for (const auto& particle : mcParticles) {

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

      // Loop over all reconstructed MC tracks
      for (auto const& track : mcTracks) {

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
    // Loop over all simulated collisions
    for (const auto& collision : collisions) {

      // Event counter: before event selection
      registryMC.fill(HIST("genEvents"), 0.5);

      // Apply event selection: require vertex position to be within the allowed z range
      if (std::fabs(collision.posZ()) > zVtx)
        continue;

      // Event counter: after event selection
      registryMC.fill(HIST("genEvents"), 1.5);

      // Loop over all MC particles
      std::vector<fastjet::PseudoJet> fjParticles;
      std::vector<TVector3> protonMomentum;
      for (const auto& particle : mcParticles) {

        // Select physical primaries within acceptance
        if (!particle.isPhysicalPrimary())
          continue;
        static constexpr double MinPtParticle = 0.1;
        if (particle.eta() < minEta || particle.eta() > maxEta || particle.pt() < MinPtParticle)
          continue;

        // Store 3-momentum vectors of antiprotons for further analysis
        if (particle.pdgCode() == PDG_t::kProtonBar) {
          TVector3 pVec(particle.px(), particle.py(), particle.pz());
          protonMomentum.push_back(pVec);
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

      // Cluster MC particles into jets using anti-kt algorithm
      fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, rJet);
      fastjet::AreaDefinition areaDef(fastjet::active_area, fastjet::GhostedAreaSpec(1.0));
      fastjet::ClusterSequenceArea cs(fjParticles, jetDef, areaDef);
      std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cs.inclusive_jets());
      auto [rhoPerp, rhoMPerp] = backgroundSub.estimateRhoPerpCone(fjParticles, jets);

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

          // Fill histogram for generated antiprotons
          registryMC.fill(HIST("antiproton_gen_jet"), particle.pt());
        }

        // Set up two perpendicular cone axes for underlying event estimation
        TVector3 jetAxis(jet.px(), jet.py(), jet.pz());
        double coneRadius = std::sqrt(jet.area() / PI);
        TVector3 ueAxis1(0, 0, 0), ueAxis2(0, 0, 0);
        getPerpendicularAxis(jetAxis, ueAxis1, +1);
        getPerpendicularAxis(jetAxis, ueAxis2, -1);

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
          if (deltaRUe1 < maxConeRadius) {
            registryMC.fill(HIST("antiproton_deltay_deltaphi_ue"), protonVec.Eta() - ueAxis1.Eta(), getDeltaPhi(protonVec.Phi(), ueAxis1.Phi()));
          }
          if (deltaRUe2 < maxConeRadius) {
            registryMC.fill(HIST("antiproton_deltay_deltaphi_ue"), protonVec.Eta() - ueAxis2.Eta(), getDeltaPhi(protonVec.Phi(), ueAxis2.Phi()));
          }

          // Fill histogram for antiprotons in the UE
          registryMC.fill(HIST("antiproton_gen_ue"), protonVec.Pt());
        }
      }
      if (isAtLeastOneJetSelected) {
        registryMC.fill(HIST("genEvents"), 3.5);
      }
    }
  }
  PROCESS_SWITCH(AntinucleiInJets, processJetsMCgen, "process jets mc gen", false);

  // Reconstructed events
  void processJetsMCrec(RecCollisionsMc const& collisions, AntiNucleiTracksMc const& mcTracks, McParticles const&)
  {
    // Loop over all reconstructed collisions
    for (const auto& collision : collisions) {

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

      // Loop over reconstructed tracks
      int id(-1);
      std::vector<fastjet::PseudoJet> fjParticles;
      std::vector<int> antiprotonTrackIndex;
      for (auto const& track : mcTracks) {
        id++;

        // Get corresponding MC particle
        if (!track.has_mcParticle())
          continue;
        const auto mcparticle = track.mcParticle();

        // Store track index for antiproton tracks
        if (passedTrackSelection(track) && track.sign() < 0 && mcparticle.pdgCode() == PDG_t::kProtonBar) {
          antiprotonTrackIndex.emplace_back(id);
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
      fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, rJet);
      fastjet::AreaDefinition areaDef(fastjet::active_area, fastjet::GhostedAreaSpec(1.0));
      fastjet::ClusterSequenceArea cs(fjParticles, jetDef, areaDef);
      std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cs.inclusive_jets());
      auto [rhoPerp, rhoMPerp] = backgroundSub.estimateRhoPerpCone(fjParticles, jets);

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
        getPerpendicularAxis(jetAxis, ueAxis1, +1);
        getPerpendicularAxis(jetAxis, ueAxis2, -1);

        // Get jet constituents
        std::vector<fastjet::PseudoJet> jetConstituents = jet.constituents();

        // Loop over jet constituents
        for (const auto& particle : jetConstituents) {

          // Get corresponding track and apply track selection criteria
          auto const& track = mcTracks.iteratorAt(particle.user_index());
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

          // Fill nsigmaITS for antiproton candidates
          if (isHighPurityAntiproton(track)) {
            registryMC.fill(HIST("antiproton_nsigma_its_mc"), pt, nSigmaITSprot);
          }

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

          // Fill histograms (TPC and TOF) only for selected candidates
          if (passedItsPidProt && nsigmaTPCPr > minNsigmaTpc && nsigmaTPCPr < maxNsigmaTpc) {
            registryMC.fill(HIST("antiproton_rec_tpc_jet"), pt);
            if (track.hasTOF() && nsigmaTOFPr > minNsigmaTof && nsigmaTOFPr < maxNsigmaTof) {
              registryMC.fill(HIST("antiproton_rec_tof_jet"), pt);
            }
          }
        }

        // Loop over tracks in the underlying event
        for (auto const& index : antiprotonTrackIndex) {

          // retrieve track associated to index
          auto const& track = mcTracks.iteratorAt(index);

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

          // Fill histograms (TPC and TOF) only for selected candidates
          if (passedItsPidProt && nsigmaTPCPr > minNsigmaTpc && nsigmaTPCPr < maxNsigmaTpc) {
            registryMC.fill(HIST("antiproton_rec_tpc_ue"), pt);
            if (track.hasTOF() && nsigmaTOFPr > minNsigmaTof && nsigmaTOFPr < maxNsigmaTof) {
              registryMC.fill(HIST("antiproton_rec_tof_ue"), pt);
            }
          }
        }
      }
      if (isAtLeastOneJetSelected) {
        registryMC.fill(HIST("recEvents"), 9.5);
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

      // Loop over all generated Monte Carlo particles for the selected event
      for (const auto& particle : mcParticles) {

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

      // Loop over reconstructed tracks
      for (auto const& track : mcTracks) {

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

    // Multiplicity percentile
    const float multiplicity = collision.centFT0M();

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
    registryCorr.fill(HIST("rho_fullEvent"), nTotAntideuteronFullEvent, nTotAntiprotonFullEvent, multiplicity);
    registryCorr.fill(HIST("rho_netP_netD_fullEvent"), netDeuteronFullEvent, netProtonFullEvent);

    // Fill efficiency histograms
    for (int i = 0; i < nBins; i++) {
      double ptAcenteri = 0.5 * (ptOverAbins[i] + ptOverAbins[i + 1]);

      registryCorr.fill(HIST("q1d_fullEvent"), nAntideuteronFullEvent[i], ptAcenteri);
      registryCorr.fill(HIST("q1p_fullEvent"), nAntiprotonFullEvent[i], ptAcenteri);
      for (int j = 0; j < nBins; j++) {
        double ptAcenterj = 0.5 * (ptOverAbins[j] + ptOverAbins[j + 1]);
        registryCorr.fill(HIST("q1d_square_fullEvent"), ptAcenteri, ptAcenterj, nAntideuteronFullEvent[i] * nAntideuteronFullEvent[j]);
        registryCorr.fill(HIST("q1p_square_fullEvent"), ptAcenteri, ptAcenterj, nAntiprotonFullEvent[i] * nAntiprotonFullEvent[j]);
        registryCorr.fill(HIST("q1d_q1p_fullEvent"), ptAcenteri, ptAcenterj, nAntideuteronFullEvent[i] * nAntiprotonFullEvent[j]);
      }
    }

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
    auto [rhoPerp, rhoMPerp] = backgroundSub.estimateRhoPerpCone(fjParticles, jets);

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
      TVector3 ueAxis1(0, 0, 0);
      TVector3 ueAxis2(0, 0, 0);
      getPerpendicularAxis(jetAxis, ueAxis1, +1);
      getPerpendicularAxis(jetAxis, ueAxis2, -1);

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

        registryCorr.fill(HIST("q1d_jet"), nAntideuteronJet[i], ptAcenteri);
        registryCorr.fill(HIST("q1p_jet"), nAntiprotonJet[i], ptAcenteri);
        for (int j = 0; j < nBins; j++) {
          double ptAcenterj = 0.5 * (ptOverAbins[j] + ptOverAbins[j + 1]);
          registryCorr.fill(HIST("q1d_square_jet"), ptAcenteri, ptAcenterj, nAntideuteronJet[i] * nAntideuteronJet[j]);
          registryCorr.fill(HIST("q1p_square_jet"), ptAcenteri, ptAcenterj, nAntiprotonJet[i] * nAntiprotonJet[j]);
          registryCorr.fill(HIST("q1d_q1p_jet"), ptAcenteri, ptAcenterj, nAntideuteronJet[i] * nAntiprotonJet[j]);
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

        registryCorr.fill(HIST("q1d_ue"), nAntideuteronUE[i], ptAcenteri);
        registryCorr.fill(HIST("q1p_ue"), nAntiprotonUE[i], ptAcenteri);
        for (int j = 0; j < nBins; j++) {
          double ptAcenterj = 0.5 * (ptOverAbins[j] + ptOverAbins[j + 1]);
          registryCorr.fill(HIST("q1d_square_ue"), ptAcenteri, ptAcenterj, nAntideuteronUE[i] * nAntideuteronUE[j]);
          registryCorr.fill(HIST("q1p_square_ue"), ptAcenteri, ptAcenterj, nAntiprotonUE[i] * nAntiprotonUE[j]);
          registryCorr.fill(HIST("q1d_q1p_ue"), ptAcenteri, ptAcenterj, nAntideuteronUE[i] * nAntiprotonUE[j]);
        }
      }
    }
    // Event counter: events with at least one jet selected
    if (isAtLeastOneJetSelected) {
      registryCorr.fill(HIST("eventCounter"), 9.5);
    }
  }
  PROCESS_SWITCH(AntinucleiInJets, processCorr, "Process Correlation analysis", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<AntinucleiInJets>(cfgc)};
}
