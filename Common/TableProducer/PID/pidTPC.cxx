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
/// \file   pidTPC.cxx
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \author Christian Sonnabend christian.sonnabend@cern.ch
/// \author Annalena Kalteyer annalena.sophie.kalteyer@cern.ch
/// \author Jeremy Wilkinson jeremy.wilkinson@cern.ch
/// \brief  Task to produce PID tables for TPC split for each particle.
///         Only the tables for the mass hypotheses requested are filled, and only for the requested table size ("Full" or "Tiny"). The others are sent empty.
///
#include <utility>
#include <map>
#include <memory>
#include <string>
#include <vector>
// ROOT includes
#include "TFile.h"
#include "TRandom.h"
#include "TSystem.h"

// O2 includes
#include "CCDB/BasicCCDBManager.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "CCDB/CcdbApi.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/Core/PID/TPCPIDResponse.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/EventSelection.h"
#include "TableHelper.h"
#include "Tools/ML/model.h"
#include "pidTPCBase.h"
#include "MetadataHelper.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::pid;
using namespace o2::pid::tpc;
using namespace o2::framework::expressions;
using namespace o2::track;
using namespace o2::ml;

MetadataHelper metadataInfo; // Metadata helper

void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  std::vector<ConfigParamSpec> options{{"add-qa", VariantType::Int, 0, {"Legacy. No effect."}}};
  std::swap(workflowOptions, options);
}

/// Task to produce the response table
struct tpcPid {
  using Trks = soa::Join<aod::Tracks, aod::TracksExtra>;
  using Coll = soa::Join<aod::Collisions, aod::PIDMults, aod::EvSels>;

  using TrksMC = soa::Join<aod::Tracks, aod::TracksExtra, aod::McTrackLabels>;
  using CollMC = soa::Join<aod::Collisions, aod::PIDMults, aod::McCollisionLabels, aod::EvSels>;

  // Tables to produce
  Produces<o2::aod::pidTPCFullEl> tablePIDFullEl;
  Produces<o2::aod::pidTPCFullMu> tablePIDFullMu;
  Produces<o2::aod::pidTPCFullPi> tablePIDFullPi;
  Produces<o2::aod::pidTPCFullKa> tablePIDFullKa;
  Produces<o2::aod::pidTPCFullPr> tablePIDFullPr;
  Produces<o2::aod::pidTPCFullDe> tablePIDFullDe;
  Produces<o2::aod::pidTPCFullTr> tablePIDFullTr;
  Produces<o2::aod::pidTPCFullHe> tablePIDFullHe;
  Produces<o2::aod::pidTPCFullAl> tablePIDFullAl;

  Produces<o2::aod::pidTPCEl> tablePIDTinyEl;
  Produces<o2::aod::pidTPCMu> tablePIDTinyMu;
  Produces<o2::aod::pidTPCPi> tablePIDTinyPi;
  Produces<o2::aod::pidTPCKa> tablePIDTinyKa;
  Produces<o2::aod::pidTPCPr> tablePIDTinyPr;
  Produces<o2::aod::pidTPCDe> tablePIDTinyDe;
  Produces<o2::aod::pidTPCTr> tablePIDTinyTr;
  Produces<o2::aod::pidTPCHe> tablePIDTinyHe;
  Produces<o2::aod::pidTPCAl> tablePIDTinyAl;
  Produces<o2::aod::mcTPCTuneOnData> tableTuneOnData;

  // TPC PID Response
  o2::pid::tpc::Response* response;

  // Network correction for TPC PID response
  OnnxModel network;
  o2::ccdb::CcdbApi ccdbApi;
  std::map<std::string, std::string> metadata;
  std::map<std::string, std::string> nullmetadata;
  std::map<std::string, std::string> headers;
  std::vector<int> speciesNetworkFlags = std::vector<int>(9);
  std::string networkVersion;

  // Input parameters
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Configurable<std::string> paramfile{"param-file", "", "Path to the parametrization object, if empty the parametrization is not taken from file"};
  Configurable<std::string> url{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> ccdbPath{"ccdbPath", "Analysis/PID/TPC/Response", "Path of the TPC parametrization on the CCDB"};
  Configurable<std::string> recoPass{"recoPass", "", "Reconstruction pass name for CCDB query (automatically takes latest object for timestamp if blank)"};
  Configurable<int64_t> ccdbTimestamp{"ccdb-timestamp", 0, "timestamp of the object used to query in CCDB the detector response. Exceptions: -1 gets the latest object, 0 gets the run dependent timestamp"};
  // Parameters for loading network from a file / downloading the file
  Configurable<bool> useNetworkCorrection{"useNetworkCorrection", 0, "(bool) Wether or not to use the network correction for the TPC dE/dx signal"};
  Configurable<bool> autofetchNetworks{"autofetchNetworks", 1, "(bool) Automatically fetches networks from CCDB for the correct run number"};
  Configurable<bool> skipTPCOnly{"skipTPCOnly", false, "Flag to skip TPC only tracks (faster but affects the analyses that use TPC only tracks)"};
  Configurable<std::string> networkPathLocally{"networkPathLocally", "network.onnx", "(std::string) Path to the local .onnx file. If autofetching is enabled, then this is where the files will be downloaded"};
  Configurable<std::string> networkPathCCDB{"networkPathCCDB", "Analysis/PID/TPC/ML", "Path on CCDB"};
  Configurable<bool> enableNetworkOptimizations{"enableNetworkOptimizations", 1, "(bool) If the neural network correction is used, this enables GraphOptimizationLevel::ORT_ENABLE_EXTENDED in the ONNX session"};
  Configurable<int> networkSetNumThreads{"networkSetNumThreads", 0, "Especially important for running on a SLURM cluster. Sets the number of threads used for execution."};
  // Configuration flags to include and exclude particle hypotheses
  Configurable<int> pidFullEl{"pid-full-el", -1, {"Produce PID information for the Electron mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  Configurable<int> pidFullMu{"pid-full-mu", -1, {"Produce PID information for the Muon mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  Configurable<int> pidFullPi{"pid-full-pi", -1, {"Produce PID information for the Pion mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  Configurable<int> pidFullKa{"pid-full-ka", -1, {"Produce PID information for the Kaon mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  Configurable<int> pidFullPr{"pid-full-pr", -1, {"Produce PID information for the Proton mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  Configurable<int> pidFullDe{"pid-full-de", -1, {"Produce PID information for the Deuterons mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  Configurable<int> pidFullTr{"pid-full-tr", -1, {"Produce PID information for the Triton mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  Configurable<int> pidFullHe{"pid-full-he", -1, {"Produce PID information for the Helium3 mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  Configurable<int> pidFullAl{"pid-full-al", -1, {"Produce PID information for the Alpha mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  Configurable<int> pidTinyEl{"pid-tiny-el", -1, {"Produce PID information for the Electron mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  Configurable<int> pidTinyMu{"pid-tiny-mu", -1, {"Produce PID information for the Muon mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  Configurable<int> pidTinyPi{"pid-tiny-pi", -1, {"Produce PID information for the Pion mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  Configurable<int> pidTinyKa{"pid-tiny-ka", -1, {"Produce PID information for the Kaon mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  Configurable<int> pidTinyPr{"pid-tiny-pr", -1, {"Produce PID information for the Proton mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  Configurable<int> pidTinyDe{"pid-tiny-de", -1, {"Produce PID information for the Deuterons mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  Configurable<int> pidTinyTr{"pid-tiny-tr", -1, {"Produce PID information for the Triton mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  Configurable<int> pidTinyHe{"pid-tiny-he", -1, {"Produce PID information for the Helium3 mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  Configurable<int> pidTinyAl{"pid-tiny-al", -1, {"Produce PID information for the Alpha mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  Configurable<int> enableTuneOnDataTable{"enableTuneOnDataTable", -1, {"Produce tuned dE/dx signal table for MC to be used as raw signal in other tasks (default -1, 'only if needed'"}};
  Configurable<int> useNetworkEl{"useNetworkEl", 1, {"Switch for applying neural network on the electron mass hypothesis (if network enabled) (set to 0 to disable)"}};
  Configurable<int> useNetworkMu{"useNetworkMu", 1, {"Switch for applying neural network on the muon mass hypothesis (if network enabled) (set to 0 to disable)"}};
  Configurable<int> useNetworkPi{"useNetworkPi", 1, {"Switch for applying neural network on the pion mass hypothesis (if network enabled) (set to 0 to disable)"}};
  Configurable<int> useNetworkKa{"useNetworkKa", 1, {"Switch for applying neural network on the kaon mass hypothesis (if network enabled) (set to 0 to disable)"}};
  Configurable<int> useNetworkPr{"useNetworkPr", 1, {"Switch for applying neural network on the proton mass hypothesis (if network enabled) (set to 0 to disable)"}};
  Configurable<int> useNetworkDe{"useNetworkDe", 1, {"Switch for applying neural network on the deuteron mass hypothesis (if network enabled) (set to 0 to disable)"}};
  Configurable<int> useNetworkTr{"useNetworkTr", 1, {"Switch for applying neural network on the triton mass hypothesis (if network enabled) (set to 0 to disable)"}};
  Configurable<int> useNetworkHe{"useNetworkHe", 1, {"Switch for applying neural network on the helium3 mass hypothesis (if network enabled) (set to 0 to disable)"}};
  Configurable<int> useNetworkAl{"useNetworkAl", 1, {"Switch for applying neural network on the alpha mass hypothesis (if network enabled) (set to 0 to disable)"}};
  Configurable<float> networkBetaGammaCutoff{"networkBetaGammaCutoff", 0.45, {"Lower value of beta-gamma to override the NN application"}};

  // Parametrization configuration
  bool useCCDBParam = false;

  void init(o2::framework::InitContext& initContext)
  {
    // Protection for process flags
    if ((doprocessStandard && doprocessMcTuneOnData) || (!doprocessStandard && !doprocessMcTuneOnData)) {
      LOG(fatal) << "pid-tpc must have only one of the options 'processStandard' OR 'processMcTuneOnData' enabled. Please check your configuration.";
    }
    response = new o2::pid::tpc::Response();
    // Checking the tables are requested in the workflow and enabling them
    auto enableFlag = [&](const std::string particle, Configurable<int>& flag) {
      enableFlagIfTableRequired(initContext, "pidTPC" + particle, flag);
    };
    enableFlag("FullEl", pidFullEl);
    enableFlag("FullMu", pidFullMu);
    enableFlag("FullPi", pidFullPi);
    enableFlag("FullKa", pidFullKa);
    enableFlag("FullPr", pidFullPr);
    enableFlag("FullDe", pidFullDe);
    enableFlag("FullTr", pidFullTr);
    enableFlag("FullHe", pidFullHe);
    enableFlag("FullAl", pidFullAl);

    enableFlag("El", pidTinyEl);
    enableFlag("Mu", pidTinyMu);
    enableFlag("Pi", pidTinyPi);
    enableFlag("Ka", pidTinyKa);
    enableFlag("Pr", pidTinyPr);
    enableFlag("De", pidTinyDe);
    enableFlag("Tr", pidTinyTr);
    enableFlag("He", pidTinyHe);
    enableFlag("Al", pidTinyAl);

    if (doprocessMcTuneOnData) {
      enableFlagIfTableRequired(initContext, "mcTPCTuneOnData", enableTuneOnDataTable);
    }

    speciesNetworkFlags[0] = useNetworkEl;
    speciesNetworkFlags[1] = useNetworkMu;
    speciesNetworkFlags[2] = useNetworkPi;
    speciesNetworkFlags[3] = useNetworkKa;
    speciesNetworkFlags[4] = useNetworkPr;
    speciesNetworkFlags[5] = useNetworkDe;
    speciesNetworkFlags[6] = useNetworkTr;
    speciesNetworkFlags[7] = useNetworkHe;
    speciesNetworkFlags[8] = useNetworkAl;

    // Initialise metadata object for CCDB calls from AO2D metadata
    if (recoPass.value == "") {
      if (metadataInfo.isFullyDefined()) {
        metadata["RecoPassName"] = metadataInfo.get("RecoPassName");
        LOGP(info, "Automatically setting reco pass for TPC Response to {} from AO2D", metadata["RecoPassName"]);
      }
    } else {
      LOGP(info, "Setting reco pass for TPC response to user-defined name {}", recoPass.value);
      metadata["RecoPassName"] = recoPass.value;
    }

    /// TPC PID Response
    const TString fname = paramfile.value;
    if (fname != "") { // Loading the parametrization from file
      LOGP(info, "Loading TPC response from file {}", fname.Data());
      try {
        std::unique_ptr<TFile> f(TFile::Open(fname, "READ"));
        f->GetObject("Response", response);
      } catch (...) {
        LOGF(fatal, "Loading the TPC PID Response from file {} failed!", fname.Data());
      }
      response->PrintAll();
    } else {
      useCCDBParam = true;
      const std::string path = ccdbPath.value;
      const auto time = ccdbTimestamp.value;
      ccdb->setURL(url.value);
      ccdb->setFatalWhenNull(false); // manual fallback in case ccdb entry empty
      ccdb->setCaching(true);
      ccdb->setLocalObjectValidityChecking();
      ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
      ccdbApi.init(url);
      if (time != 0) {
        LOGP(info, "Initialising TPC PID response for fixed timestamp {} and reco pass {}:", time, recoPass.value);
        ccdb->setTimestamp(time);
        response = ccdb->getSpecific<o2::pid::tpc::Response>(path, time, metadata);
        headers = ccdbApi.retrieveHeaders(path, metadata, time);
        if (!response) {
          LOGF(warning, "Unable to find TPC parametrisation for specified pass name - falling back to latest object");
          response = ccdb->getForTimeStamp<o2::pid::tpc::Response>(path, time);
          headers = ccdbApi.retrieveHeaders(path, metadata, time);
          networkVersion = headers["NN-Version"];
          if (!response) {
            LOGF(fatal, "Unable to find any TPC object corresponding to timestamp {}!", time);
          }
        }
        LOG(info) << "Successfully retrieved TPC PID object from CCDB for timestamp " << time << ", period " << headers["LPMProductionTag"] << ", recoPass " << headers["RecoPassName"];
        metadata["RecoPassName"] = headers["RecoPassName"]; // Force pass number for NN request to match retrieved BB
        response->PrintAll();
      }
    }

    /// Neural network init for TPC PID

    if (!useNetworkCorrection) {
      return;
    } else {
      /// CCDB and auto-fetching

      if (!autofetchNetworks) {
        if (ccdbTimestamp > 0) {
          /// Fetching network for specific timestamp
          LOG(info) << "Fetching network for timestamp: " << ccdbTimestamp.value;
          bool retrieveSuccess = ccdbApi.retrieveBlob(networkPathCCDB.value, ".", metadata, ccdbTimestamp.value, false, networkPathLocally.value);
          headers = ccdbApi.retrieveHeaders(networkPathCCDB.value, metadata, ccdbTimestamp.value);
          networkVersion = headers["NN-Version"];
          if (retrieveSuccess) {
            network.initModel(networkPathLocally.value, enableNetworkOptimizations.value, networkSetNumThreads.value, strtoul(headers["Valid-From"].c_str(), NULL, 0), strtoul(headers["Valid-Until"].c_str(), NULL, 0));
            std::vector<float> dummyInput(network.getNumInputNodes(), 1.);
            network.evalModel(dummyInput); /// Init the model evaluations
            LOGP(info, "Retrieved NN corrections for production tag {}, pass number {}, and NN-Version {}", headers["LPMProductionTag"], headers["RecoPassName"], headers["NN-Version"]);
          } else {
            LOG(fatal) << "No valid NN object found matching retrieved Bethe-Bloch parametrisation for pass " << metadata["RecoPassName"] << ". Please ensure that the requested pass has dedicated NN corrections available";
          }
        } else {
          /// Taking the network from local file
          if (networkPathLocally.value == "") {
            LOG(fatal) << "Local path must be set (flag networkPathLocally)! Aborting...";
          }
          LOG(info) << "Using local file [" << networkPathLocally.value << "] for the TPC PID response correction.";
          network.initModel(networkPathLocally.value, enableNetworkOptimizations.value, networkSetNumThreads.value);
          std::vector<float> dummyInput(network.getNumInputNodes(), 1.);
          network.evalModel(dummyInput); // This is an initialisation and might reduce the overhead of the model
        }
      } else {
        return;
      }
    }
  }

  Partition<Trks> notTPCStandaloneTracks = (aod::track::tpcNClsFindable > static_cast<uint8_t>(0)) && ((aod::track::itsClusterSizes > static_cast<uint32_t>(0)) || (aod::track::trdPattern > static_cast<uint8_t>(0)) || (aod::track::tofExpMom > 0.f && aod::track::tofChi2 > 0.f)); // To count number of tracks for use in NN array
  Partition<Trks> tracksWithTPC = (aod::track::tpcNClsFindable > (uint8_t)0);

  template <typename C, typename T, typename B>
  std::vector<float> createNetworkPrediction(C const& collisions, T const& tracks, B const& bcs, const size_t size)
  {

    std::vector<float> network_prediction;

    auto start_network_total = std::chrono::high_resolution_clock::now();
    if (autofetchNetworks) {
      const auto& bc = bcs.begin();
      // Initialise correct TPC response object before NN setup (for NCl normalisation)
      if (useCCDBParam && ccdbTimestamp.value == 0 && !ccdb->isCachedObjectValid(ccdbPath.value, bc.timestamp())) { // Updating parametrisation only if the initial timestamp is 0
        if (recoPass.value == "") {
          LOGP(info, "Retrieving latest TPC response object for timestamp {}:", bc.timestamp());
        } else {
          LOGP(info, "Retrieving TPC Response for timestamp {} and recoPass {}:", bc.timestamp(), recoPass.value);
        }
        response = ccdb->getSpecific<o2::pid::tpc::Response>(ccdbPath.value, bc.timestamp(), metadata);
        headers = ccdbApi.retrieveHeaders(ccdbPath.value, metadata, bc.timestamp());
        networkVersion = headers["NN-Version"];
        if (!response) {
          LOGP(warning, "!! Could not find a valid TPC response object for specific pass name {}! Falling back to latest uploaded object.", metadata["RecoPassName"]);
          headers = ccdbApi.retrieveHeaders(ccdbPath.value, nullmetadata, bc.timestamp());
          response = ccdb->getForTimeStamp<o2::pid::tpc::Response>(ccdbPath.value, bc.timestamp());
          if (!response) {
            LOGP(fatal, "Could not find ANY TPC response object for the timestamp {}!", bc.timestamp());
          }
        }
        LOG(info) << "Successfully retrieved TPC PID object from CCDB for timestamp " << bc.timestamp() << ", period " << headers["LPMProductionTag"] << ", recoPass " << headers["RecoPassName"];
        metadata["RecoPassName"] = headers["RecoPassName"]; // Force pass number for NN request to match retrieved BB
        response->PrintAll();
      }

      if (bc.timestamp() < network.getValidityFrom() || bc.timestamp() > network.getValidityUntil()) { // fetches network only if the runnumbers change
        LOG(info) << "Fetching network for timestamp: " << bc.timestamp();
        bool retrieveSuccess = ccdbApi.retrieveBlob(networkPathCCDB.value, ".", metadata, bc.timestamp(), false, networkPathLocally.value);
        headers = ccdbApi.retrieveHeaders(networkPathCCDB.value, metadata, bc.timestamp());
        networkVersion = headers["NN-Version"];
        if (retrieveSuccess) {
          network.initModel(networkPathLocally.value, enableNetworkOptimizations.value, networkSetNumThreads.value, strtoul(headers["Valid-From"].c_str(), NULL, 0), strtoul(headers["Valid-Until"].c_str(), NULL, 0));
          std::vector<float> dummyInput(network.getNumInputNodes(), 1.);
          network.evalModel(dummyInput);
          LOGP(info, "Retrieved NN corrections for production tag {}, pass number {}, NN-Version number{}", headers["LPMProductionTag"], headers["RecoPassName"], headers["NN-Version"]);
        } else {
          LOG(fatal) << "No valid NN object found matching retrieved Bethe-Bloch parametrisation for pass " << metadata["RecoPassName"] << ". Please ensure that the requested pass has dedicated NN corrections available";
        }
      }
    }

    // Defining some network parameters
    int input_dimensions = network.getNumInputNodes();
    int output_dimensions = network.getNumOutputNodes();
    const uint64_t track_prop_size = input_dimensions * size;
    const uint64_t prediction_size = output_dimensions * size;

    network_prediction = std::vector<float>(prediction_size * 9); // For each mass hypotheses
    const float nNclNormalization = response->GetNClNormalization();
    float duration_network = 0;

    std::vector<float> track_properties(track_prop_size);
    uint64_t counter_track_props = 0;
    int loop_counter = 0;

    // Filling a std::vector<float> to be evaluated by the network
    // Evaluation on single tracks brings huge overhead: Thus evaluation is done on one large vector
    for (int i = 0; i < 9; i++) { // Loop over particle number for which network correction is used
      for (auto const& trk : tracks) {
        if (!trk.hasTPC()) {
          continue;
        }
        if (skipTPCOnly) {
          if (!trk.hasITS() && !trk.hasTRD() && !trk.hasTOF()) {
            continue;
          }
        }
        track_properties[counter_track_props] = trk.tpcInnerParam();
        track_properties[counter_track_props + 1] = trk.tgl();
        track_properties[counter_track_props + 2] = trk.signed1Pt();
        track_properties[counter_track_props + 3] = o2::track::pid_constants::sMasses[i];
        track_properties[counter_track_props + 4] = trk.has_collision() ? collisions.iteratorAt(trk.collisionId()).multTPC() / 11000. : 1.;
        track_properties[counter_track_props + 5] = std::sqrt(nNclNormalization / trk.tpcNClsFound());
        if (input_dimensions == 7 && networkVersion == "2") {
          track_properties[counter_track_props + 6] = trk.has_collision() ? collisions.iteratorAt(trk.collisionId()).ft0cOccupancyInTimeRange() / 60000. : 1.;
        }
        counter_track_props += input_dimensions;
      }

      auto start_network_eval = std::chrono::high_resolution_clock::now();
      float* output_network = network.evalModel(track_properties);
      auto stop_network_eval = std::chrono::high_resolution_clock::now();
      duration_network += std::chrono::duration<float, std::ratio<1, 1000000000>>(stop_network_eval - start_network_eval).count();
      for (uint64_t i = 0; i < prediction_size; i += output_dimensions) {
        for (int j = 0; j < output_dimensions; j++) {
          network_prediction[i + j + prediction_size * loop_counter] = output_network[i + j];
        }
      }

      counter_track_props = 0;
      loop_counter += 1;
    }
    track_properties.clear();

    auto stop_network_total = std::chrono::high_resolution_clock::now();
    LOG(debug) << "Neural Network for the TPC PID response correction: Time per track (eval ONNX): " << duration_network / (size * 9) << "ns ; Total time (eval ONNX): " << duration_network / 1000000000 << " s";
    LOG(debug) << "Neural Network for the TPC PID response correction: Time per track (eval + overhead): " << std::chrono::duration<float, std::ratio<1, 1000000000>>(stop_network_total - start_network_total).count() / (size * 9) << "ns ; Total time (eval + overhead): " << std::chrono::duration<float, std::ratio<1, 1000000000>>(stop_network_total - start_network_total).count() / 1000000000 << " s";

    return network_prediction;
  }

  template <typename C, typename T, typename NSF, typename NST>
  void makePidTables(const int flagFull, NSF& tableFull, const int flagTiny, NST& tableTiny, const o2::track::PID::ID pid, const float tpcSignal, const T& trk, const C& collisions, const std::vector<float>& network_prediction, const int& count_tracks, const int& tracksForNet_size)
  {
    if (flagFull != 1 && flagTiny != 1) {
      return;
    }
    if (!trk.hasTPC() || tpcSignal < 0.f) {
      if (flagFull)
        tableFull(-999.f, -999.f);
      if (flagTiny)
        tableTiny(aod::pidtpc_tiny::binning::underflowBin);
      return;
    }
    if (skipTPCOnly) {
      if (!trk.hasITS() && !trk.hasTRD() && !trk.hasTOF()) {
        if (flagFull)
          tableFull(-999.f, -999.f);
        if (flagTiny)
          tableTiny(aod::pidtpc_tiny::binning::underflowBin);
        return;
      }
    }
    auto expSignal = response->GetExpectedSignal(trk, pid);
    auto expSigma = trk.has_collision() ? response->GetExpectedSigma(collisions.iteratorAt(trk.collisionId()), trk, pid) : 0.07 * expSignal; // use default sigma value of 7% if no collision information to estimate resolution
    if (expSignal < 0. || expSigma < 0.) {                                                                                                   // skip if expected signal invalid
      if (flagFull)
        tableFull(-999.f, -999.f);
      if (flagTiny)
        tableTiny(aod::pidtpc_tiny::binning::underflowBin);
      return;
    }

    float nSigma = -999.f;
    float bg = trk.tpcInnerParam() / o2::track::pid_constants::sMasses[pid]; // estimated beta-gamma for network cutoff
    if (useNetworkCorrection && speciesNetworkFlags[pid] && trk.has_collision() && bg > networkBetaGammaCutoff) {

      // Here comes the application of the network. The output--dimensions of the network determine the application: 1: mean, 2: sigma, 3: sigma asymmetric
      // For now only the option 2: sigma will be used. The other options are kept if there would be demand later on
      if (network.getNumOutputNodes() == 1) { // Expected mean correction; no sigma correction
        nSigma = (tpcSignal - network_prediction[count_tracks + tracksForNet_size * pid] * expSignal) / expSigma;
      } else if (network.getNumOutputNodes() == 2) { // Symmetric sigma correction
        expSigma = (network_prediction[2 * (count_tracks + tracksForNet_size * pid) + 1] - network_prediction[2 * (count_tracks + tracksForNet_size * pid)]) * expSignal;
        nSigma = (tpcSignal / expSignal - network_prediction[2 * (count_tracks + tracksForNet_size * pid)]) / (network_prediction[2 * (count_tracks + tracksForNet_size * pid) + 1] - network_prediction[2 * (count_tracks + tracksForNet_size * pid)]);
      } else if (network.getNumOutputNodes() == 3) { // Asymmetric sigma corection
        if (tpcSignal / expSignal >= network_prediction[3 * (count_tracks + tracksForNet_size * pid)]) {
          expSigma = (network_prediction[3 * (count_tracks + tracksForNet_size * pid) + 1] - network_prediction[3 * (count_tracks + tracksForNet_size * pid)]) * expSignal;
          nSigma = (tpcSignal / expSignal - network_prediction[3 * (count_tracks + tracksForNet_size * pid)]) / (network_prediction[3 * (count_tracks + tracksForNet_size * pid) + 1] - network_prediction[3 * (count_tracks + tracksForNet_size * pid)]);
        } else {
          expSigma = (network_prediction[3 * (count_tracks + tracksForNet_size * pid)] - network_prediction[3 * (count_tracks + tracksForNet_size * pid) + 2]) * expSignal;
          nSigma = (tpcSignal / expSignal - network_prediction[3 * (count_tracks + tracksForNet_size * pid)]) / (network_prediction[3 * (count_tracks + tracksForNet_size * pid)] - network_prediction[3 * (count_tracks + tracksForNet_size * pid) + 2]);
        }
      } else {
        LOGF(fatal, "Network output-dimensions incompatible!");
      }
    } else {
      nSigma = response->GetNumberOfSigmaMCTuned(collisions.iteratorAt(trk.collisionId()), trk, pid, tpcSignal);
    }
    if (flagFull)
      tableFull(expSigma, nSigma);
    if (flagTiny)
      aod::pidutils::packInTable<aod::pidtpc_tiny::binning>(nSigma, tableTiny);
  };

  void processStandard(Coll const& collisions, Trks const& tracks, aod::BCsWithTimestamps const& bcs)
  {

    const uint64_t outTable_size = tracks.size();

    auto reserveTable = [&outTable_size](const Configurable<int>& flag, auto& table) {
      if (flag.value != 1) {
        return;
      }
      table.reserve(outTable_size);
    };

    // Prepare memory for enabled tables
    reserveTable(pidFullEl, tablePIDFullEl);
    reserveTable(pidFullMu, tablePIDFullMu);
    reserveTable(pidFullPi, tablePIDFullPi);
    reserveTable(pidFullKa, tablePIDFullKa);
    reserveTable(pidFullPr, tablePIDFullPr);
    reserveTable(pidFullDe, tablePIDFullDe);
    reserveTable(pidFullTr, tablePIDFullTr);
    reserveTable(pidFullHe, tablePIDFullHe);
    reserveTable(pidFullAl, tablePIDFullAl);

    reserveTable(pidTinyEl, tablePIDTinyEl);
    reserveTable(pidTinyMu, tablePIDTinyMu);
    reserveTable(pidTinyPi, tablePIDTinyPi);
    reserveTable(pidTinyKa, tablePIDTinyKa);
    reserveTable(pidTinyPr, tablePIDTinyPr);
    reserveTable(pidTinyDe, tablePIDTinyDe);
    reserveTable(pidTinyTr, tablePIDTinyTr);
    reserveTable(pidTinyHe, tablePIDTinyHe);
    reserveTable(pidTinyAl, tablePIDTinyAl);

    const uint64_t tracksForNet_size = (skipTPCOnly) ? notTPCStandaloneTracks.size() : tracksWithTPC.size();
    std::vector<float> network_prediction;

    if (useNetworkCorrection) {
      network_prediction = createNetworkPrediction(collisions, tracks, bcs, tracksForNet_size);
    }

    uint64_t count_tracks = 0;

    for (auto const& trk : tracks) {
      // Loop on Tracks

      const auto& bc = trk.has_collision() ? collisions.iteratorAt(trk.collisionId()).bc_as<aod::BCsWithTimestamps>() : bcs.begin();
      if (useCCDBParam && ccdbTimestamp.value == 0 && !ccdb->isCachedObjectValid(ccdbPath.value, bc.timestamp())) { // Updating parametrisation only if the initial timestamp is 0
        if (recoPass.value == "") {
          LOGP(info, "Retrieving latest TPC response object for timestamp {}:", bc.timestamp());
        } else {
          LOGP(info, "Retrieving TPC Response for timestamp {} and recoPass {}:", bc.timestamp(), recoPass.value);
        }
        response = ccdb->getSpecific<o2::pid::tpc::Response>(ccdbPath.value, bc.timestamp(), metadata);
        headers = ccdbApi.retrieveHeaders(ccdbPath.value, metadata, bc.timestamp());
        if (!response) {
          LOGP(warning, "!! Could not find a valid TPC response object for specific pass name {}! Falling back to latest uploaded object.", metadata["RecoPassName"]);
          response = ccdb->getForTimeStamp<o2::pid::tpc::Response>(ccdbPath.value, bc.timestamp());
          headers = ccdbApi.retrieveHeaders(ccdbPath.value, nullmetadata, bc.timestamp());
          if (!response) {
            LOGP(fatal, "Could not find ANY TPC response object for the timestamp {}!", bc.timestamp());
          }
        }
        LOG(info) << "Successfully retrieved TPC PID object from CCDB for timestamp " << bc.timestamp() << ", period " << headers["LPMProductionTag"] << ", recoPass " << headers["RecoPassName"];
        response->PrintAll();
      }

      auto makePidTablesDefault = [&trk, &collisions, &network_prediction, &count_tracks, &tracksForNet_size, this](const int flagFull, auto& tableFull, const int flagTiny, auto& tableTiny, const o2::track::PID::ID pid) {
        makePidTables(flagFull, tableFull, flagTiny, tableTiny, pid, trk.tpcSignal(), trk, collisions, network_prediction, count_tracks, tracksForNet_size);
      };

      makePidTablesDefault(pidFullEl, tablePIDFullEl, pidTinyEl, tablePIDTinyEl, o2::track::PID::Electron);
      makePidTablesDefault(pidFullMu, tablePIDFullMu, pidTinyMu, tablePIDTinyMu, o2::track::PID::Muon);
      makePidTablesDefault(pidFullPi, tablePIDFullPi, pidTinyPi, tablePIDTinyPi, o2::track::PID::Pion);
      makePidTablesDefault(pidFullKa, tablePIDFullKa, pidTinyKa, tablePIDTinyKa, o2::track::PID::Kaon);
      makePidTablesDefault(pidFullPr, tablePIDFullPr, pidTinyPr, tablePIDTinyPr, o2::track::PID::Proton);
      makePidTablesDefault(pidFullDe, tablePIDFullDe, pidTinyDe, tablePIDTinyDe, o2::track::PID::Deuteron);
      makePidTablesDefault(pidFullTr, tablePIDFullTr, pidTinyTr, tablePIDTinyTr, o2::track::PID::Triton);
      makePidTablesDefault(pidFullHe, tablePIDFullHe, pidTinyHe, tablePIDTinyHe, o2::track::PID::Helium3);
      makePidTablesDefault(pidFullAl, tablePIDFullAl, pidTinyAl, tablePIDTinyAl, o2::track::PID::Alpha);

      if (trk.hasTPC() && (!skipTPCOnly || trk.hasITS() || trk.hasTRD() || trk.hasTOF())) {
        count_tracks++; // Increment network track counter only if track has TPC, and (not skipping TPConly) or (is not TPConly)
      }
    }
  }

  PROCESS_SWITCH(tpcPid, processStandard, "Creating PID tables without MC TuneOnData", true);

  Partition<TrksMC> mcnotTPCStandaloneTracks = (aod::track::tpcNClsFindable > static_cast<uint8_t>(0)) && ((aod::track::itsClusterSizes > static_cast<uint32_t>(0)) || (aod::track::trdPattern > static_cast<uint8_t>(0)) || (aod::track::tofExpMom > 0.f && aod::track::tofChi2 > 0.f)); // To count number of tracks for use in NN array
  Partition<TrksMC> mctracksWithTPC = (aod::track::tpcNClsFindable > (uint8_t)0);

  void processMcTuneOnData(CollMC const& collisionsMc, TrksMC const& tracksMc, aod::BCsWithTimestamps const& bcs, aod::McParticles const&)
  {
    gRandom->SetSeed(0); // Ensure unique seed from UUID for each process call
    const uint64_t outTable_size = tracksMc.size();

    auto reserveTable = [&outTable_size](const Configurable<int>& flag, auto& table) {
      if (flag.value != 1) {
        return;
      }
      table.reserve(outTable_size);
    };

    // Prepare memory for enabled tables
    reserveTable(pidFullEl, tablePIDFullEl);
    reserveTable(pidFullMu, tablePIDFullMu);
    reserveTable(pidFullPi, tablePIDFullPi);
    reserveTable(pidFullKa, tablePIDFullKa);
    reserveTable(pidFullPr, tablePIDFullPr);
    reserveTable(pidFullDe, tablePIDFullDe);
    reserveTable(pidFullTr, tablePIDFullTr);
    reserveTable(pidFullHe, tablePIDFullHe);
    reserveTable(pidFullAl, tablePIDFullAl);

    reserveTable(pidTinyEl, tablePIDTinyEl);
    reserveTable(pidTinyMu, tablePIDTinyMu);
    reserveTable(pidTinyPi, tablePIDTinyPi);
    reserveTable(pidTinyKa, tablePIDTinyKa);
    reserveTable(pidTinyPr, tablePIDTinyPr);
    reserveTable(pidTinyDe, tablePIDTinyDe);
    reserveTable(pidTinyTr, tablePIDTinyTr);
    reserveTable(pidTinyHe, tablePIDTinyHe);
    reserveTable(pidTinyAl, tablePIDTinyAl);

    reserveTable(enableTuneOnDataTable, tableTuneOnData); // Only produce the table of tuned dE/dx if the signal is requested by another task

    const uint64_t tracksForNet_size = (skipTPCOnly) ? mcnotTPCStandaloneTracks.size() : mctracksWithTPC.size();
    std::vector<float> network_prediction;

    if (useNetworkCorrection) {
      network_prediction = createNetworkPrediction(collisionsMc, tracksMc, bcs, tracksForNet_size);
    }

    uint64_t count_tracks = 0;

    for (auto const& trk : tracksMc) {
      // Loop on Tracks
      const auto& bc = trk.has_collision() ? collisionsMc.iteratorAt(trk.collisionId()).bc_as<aod::BCsWithTimestamps>() : bcs.begin();
      if (useCCDBParam && ccdbTimestamp.value == 0 && !ccdb->isCachedObjectValid(ccdbPath.value, bc.timestamp())) { // Updating parametrisation only if the initial timestamp is 0
        if (recoPass.value == "") {
          LOGP(info, "Retrieving latest TPC response object for timestamp {}:", bc.timestamp());
        } else {
          LOGP(info, "Retrieving TPC Response for timestamp {} and recoPass {}:", bc.timestamp(), recoPass.value);
        }
        response = ccdb->getSpecific<o2::pid::tpc::Response>(ccdbPath.value, bc.timestamp(), metadata);
        if (!response) {
          LOGP(warning, "!! Could not find a valid TPC response object for specific pass name {}! Falling back to latest uploaded object.", metadata["RecoPassName"]);
          response = ccdb->getForTimeStamp<o2::pid::tpc::Response>(ccdbPath.value, bc.timestamp());
          if (!response) {
            LOGP(fatal, "Could not find ANY TPC response object for the timestamp {}!", bc.timestamp());
          }
        }
        response->PrintAll();
      }

      // Perform TuneOnData sampling for MC dE/dx
      float mcTunedTPCSignal = 0.;
      if (!trk.hasTPC()) {
        mcTunedTPCSignal = -999.f;
      } else {
        if (skipTPCOnly) {
          if (!trk.hasITS() && !trk.hasTRD() && !trk.hasTOF()) {
            mcTunedTPCSignal = -999.f;
          }
        }
        int pid = getPIDIndex(trk.mcParticle().pdgCode());

        auto expSignal = response->GetExpectedSignal(trk, pid);
        auto expSigma = response->GetExpectedSigma(collisionsMc.iteratorAt(trk.collisionId()), trk, pid);
        if (expSignal < 0. || expSigma < 0.) { // if expectation invalid then give undefined signal
          mcTunedTPCSignal = -999.f;
        }
        float bg = trk.tpcInnerParam() / o2::track::pid_constants::sMasses[pid]; // estimated beta-gamma for network cutoff

        if (useNetworkCorrection && speciesNetworkFlags[pid] && trk.has_collision() && bg > networkBetaGammaCutoff) {
          auto mean = network_prediction[2 * (count_tracks + tracksForNet_size * pid)] * expSignal; // Absolute mean, i.e. the mean dE/dx value of the data in that slice, not the mean of the NSigma distribution
          auto sigma = (network_prediction[2 * (count_tracks + tracksForNet_size * pid) + 1] - network_prediction[2 * (count_tracks + tracksForNet_size * pid)]) * expSignal;
          if (mean < 0.f || sigma < 0.f) {
            mcTunedTPCSignal = -999.f;
          } else {
            mcTunedTPCSignal = gRandom->Gaus(mean, sigma);
          }
        } else {
          mcTunedTPCSignal = gRandom->Gaus(expSignal, expSigma);
        }
      }
      if (enableTuneOnDataTable)
        tableTuneOnData(mcTunedTPCSignal);

      // Check and fill enabled nsigma tables

      auto makePidTablesMCTune = [&trk, &collisionsMc, &network_prediction, &count_tracks, &tracksForNet_size, &mcTunedTPCSignal, this](const int flagFull, auto& tableFull, const int flagTiny, auto& tableTiny, const o2::track::PID::ID pid) {
        makePidTables(flagFull, tableFull, flagTiny, tableTiny, pid, mcTunedTPCSignal, trk, collisionsMc, network_prediction, count_tracks, tracksForNet_size);
      };

      makePidTablesMCTune(pidFullEl, tablePIDFullEl, pidTinyEl, tablePIDTinyEl, o2::track::PID::Electron);
      makePidTablesMCTune(pidFullMu, tablePIDFullMu, pidTinyMu, tablePIDTinyMu, o2::track::PID::Muon);
      makePidTablesMCTune(pidFullPi, tablePIDFullPi, pidTinyPi, tablePIDTinyPi, o2::track::PID::Pion);
      makePidTablesMCTune(pidFullKa, tablePIDFullKa, pidTinyKa, tablePIDTinyKa, o2::track::PID::Kaon);
      makePidTablesMCTune(pidFullPr, tablePIDFullPr, pidTinyPr, tablePIDTinyPr, o2::track::PID::Proton);
      makePidTablesMCTune(pidFullDe, tablePIDFullDe, pidTinyDe, tablePIDTinyDe, o2::track::PID::Deuteron);
      makePidTablesMCTune(pidFullTr, tablePIDFullTr, pidTinyTr, tablePIDTinyTr, o2::track::PID::Triton);
      makePidTablesMCTune(pidFullHe, tablePIDFullHe, pidTinyHe, tablePIDTinyHe, o2::track::PID::Helium3);
      makePidTablesMCTune(pidFullAl, tablePIDFullAl, pidTinyAl, tablePIDTinyAl, o2::track::PID::Alpha);

      if (trk.hasTPC() && (!skipTPCOnly || trk.hasITS() || trk.hasTRD() || trk.hasTOF())) {
        count_tracks++; // Increment network track counter only if track has TPC, and (not skipping TPConly) or (is not TPConly)
      }
    }
  }

  PROCESS_SWITCH(tpcPid, processMcTuneOnData, "Creating PID tables with MC TuneOnData", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  metadataInfo.initMetadata(cfgc); // Parse AO2D metadata
  return WorkflowSpec{adaptAnalysisTask<tpcPid>(cfgc)};
}
