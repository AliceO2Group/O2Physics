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

/// \file pidTPCModule.h
/// \brief  Task to produce PID tables for TPC split for each particle.
///         Only the tables for the mass hypotheses requested are filled, and only for the requested table size
///         ("Full" or "Tiny"). The others are sent empty.
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \author Christian Sonnabend christian.sonnabend@cern.ch
/// \author Annalena Kalteyer annalena.sophie.kalteyer@cern.ch
/// \author Jeremy Wilkinson jeremy.wilkinson@cern.ch

#ifndef COMMON_TOOLS_PIDTPCMODULE_H_
#define COMMON_TOOLS_PIDTPCMODULE_H_

#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>
// ROOT includes
#include "TFile.h"
#include "TRandom.h"
#include "TSystem.h"

// O2 includes
#include "MetadataHelper.h"
#include "TableHelper.h"
#include "pidTPCBase.h"

#include "Common/Core/PID/TPCPIDResponse.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Tools/ML/model.h"

#include "CCDB/BasicCCDBManager.h"
#include "CCDB/CcdbApi.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

namespace o2::aod
{
namespace pid
{
struct pidTPCProducts : o2::framework::ProducesGroup {
  // Intermediate tables (provide only if requested)
  o2::framework::Produces<aod::DEdxsCorrected> dEdxCorrected;
  o2::framework::Produces<aod::PIDMults> mult;

  // Tables produced by TPC component 
  o2::framework::Produces<o2::aod::pidTPCFullEl> tablePIDFullEl;
  o2::framework::Produces<o2::aod::pidTPCFullMu> tablePIDFullMu;
  o2::framework::Produces<o2::aod::pidTPCFullPi> tablePIDFullPi;
  o2::framework::Produces<o2::aod::pidTPCFullKa> tablePIDFullKa;
  o2::framework::Produces<o2::aod::pidTPCFullPr> tablePIDFullPr;
  o2::framework::Produces<o2::aod::pidTPCFullDe> tablePIDFullDe;
  o2::framework::Produces<o2::aod::pidTPCFullTr> tablePIDFullTr;
  o2::framework::Produces<o2::aod::pidTPCFullHe> tablePIDFullHe;
  o2::framework::Produces<o2::aod::pidTPCFullAl> tablePIDFullAl;

  o2::framework::Produces<o2::aod::pidTPCEl> tablePIDTinyEl;
  o2::framework::Produces<o2::aod::pidTPCMu> tablePIDTinyMu;
  o2::framework::Produces<o2::aod::pidTPCPi> tablePIDTinyPi;
  o2::framework::Produces<o2::aod::pidTPCKa> tablePIDTinyKa;
  o2::framework::Produces<o2::aod::pidTPCPr> tablePIDTinyPr;
  o2::framework::Produces<o2::aod::pidTPCDe> tablePIDTinyDe;
  o2::framework::Produces<o2::aod::pidTPCTr> tablePIDTinyTr;
  o2::framework::Produces<o2::aod::pidTPCHe> tablePIDTinyHe;
  o2::framework::Produces<o2::aod::pidTPCAl> tablePIDTinyAl;
  o2::framework::Produces<o2::aod::mcTPCTuneOnData> tableTuneOnData;
};

struct pidTPCConfigurables : o2::framework::ConfigurableGroup {
  std::string prefix = "pidTPC";
  o2::framework::Configurable<std::string> paramfile{"param-file", "", "Path to the parametrization object, if empty the parametrization is not taken from file"};
  o2::framework::Configurable<std::string> ccdbPath{"ccdbPath", "Analysis/PID/TPC/Response", "Path of the TPC parametrization on the CCDB"};
  o2::framework::Configurable<std::string> recoPass{"recoPass", "", "Reconstruction pass name for CCDB query (automatically takes latest object for timestamp if blank)"};
  o2::framework::Configurable<int64_t> ccdbTimestamp{"ccdb-timestamp", 0, "timestamp of the object used to query in CCDB the detector response. Exceptions: -1 gets the latest object, 0 gets the run dependent timestamp"};
  // Parameters for loading network from a file / downloading the file
  o2::framework::Configurable<bool> useNetworkCorrection{"useNetworkCorrection", 0, "(bool) Wether or not to use the network correction for the TPC dE/dx signal"};
  o2::framework::Configurable<bool> autofetchNetworks{"autofetchNetworks", 1, "(bool) Automatically fetches networks from CCDB for the correct run number"};
  o2::framework::Configurable<bool> skipTPCOnly{"skipTPCOnly", false, "Flag to skip TPC only tracks (faster but affects the analyses that use TPC only tracks)"};
  o2::framework::Configurable<std::string> networkPathLocally{"networkPathLocally", "network.onnx", "(std::string) Path to the local .onnx file. If autofetching is enabled, then this is where the files will be downloaded"};
  o2::framework::Configurable<std::string> networkPathCCDB{"networkPathCCDB", "Analysis/PID/TPC/ML", "Path on CCDB"};
  o2::framework::Configurable<bool> enableNetworkOptimizations{"enableNetworkOptimizations", 1, "(bool) If the neural network correction is used, this enables GraphOptimizationLevel::ORT_ENABLE_EXTENDED in the ONNX session"};
  o2::framework::Configurable<int> networkSetNumThreads{"networkSetNumThreads", 0, "Especially important for running on a SLURM cluster. Sets the number of threads used for execution."};
  // Configuration flags to include and exclude particle hypotheses
  o2::framework::Configurable<int> savedEdxsCorrected{"savedEdxsCorrected", -1, {"Save table with corrected dE/dx calculated on the spot. 0: off, 1: on, -1: auto"}};
  o2::framework::Configurable<bool> useCorrecteddEdx{"useCorrecteddEdx", false, "(bool) If true, use corrected dEdx value in Nsigma calculation instead of the one in the AO2D"};

  o2::framework::Configurable<int> pidFullEl{"pid-full-el", -1, {"Produce PID information for the Electron mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  o2::framework::Configurable<int> pidFullMu{"pid-full-mu", -1, {"Produce PID information for the Muon mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  o2::framework::Configurable<int> pidFullPi{"pid-full-pi", -1, {"Produce PID information for the Pion mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  o2::framework::Configurable<int> pidFullKa{"pid-full-ka", -1, {"Produce PID information for the Kaon mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  o2::framework::Configurable<int> pidFullPr{"pid-full-pr", -1, {"Produce PID information for the Proton mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  o2::framework::Configurable<int> pidFullDe{"pid-full-de", -1, {"Produce PID information for the Deuterons mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  o2::framework::Configurable<int> pidFullTr{"pid-full-tr", -1, {"Produce PID information for the Triton mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  o2::framework::Configurable<int> pidFullHe{"pid-full-he", -1, {"Produce PID information for the Helium3 mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  o2::framework::Configurable<int> pidFullAl{"pid-full-al", -1, {"Produce PID information for the Alpha mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  o2::framework::Configurable<int> pidTinyEl{"pid-tiny-el", -1, {"Produce PID information for the Electron mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  o2::framework::Configurable<int> pidTinyMu{"pid-tiny-mu", -1, {"Produce PID information for the Muon mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  o2::framework::Configurable<int> pidTinyPi{"pid-tiny-pi", -1, {"Produce PID information for the Pion mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  o2::framework::Configurable<int> pidTinyKa{"pid-tiny-ka", -1, {"Produce PID information for the Kaon mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  o2::framework::Configurable<int> pidTinyPr{"pid-tiny-pr", -1, {"Produce PID information for the Proton mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  o2::framework::Configurable<int> pidTinyDe{"pid-tiny-de", -1, {"Produce PID information for the Deuterons mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  o2::framework::Configurable<int> pidTinyTr{"pid-tiny-tr", -1, {"Produce PID information for the Triton mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  o2::framework::Configurable<int> pidTinyHe{"pid-tiny-he", -1, {"Produce PID information for the Helium3 mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  o2::framework::Configurable<int> pidTinyAl{"pid-tiny-al", -1, {"Produce PID information for the Alpha mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  o2::framework::Configurable<int> enableTuneOnDataTable{"enableTuneOnDataTable", -1, {"Produce tuned dE/dx signal table for MC to be used as raw signal in other tasks (default -1, 'only if needed'"}};
  o2::framework::Configurable<int> useNetworkEl{"useNetworkEl", 1, {"Switch for applying neural network on the electron mass hypothesis (if network enabled) (set to 0 to disable)"}};
  o2::framework::Configurable<int> useNetworkMu{"useNetworkMu", 1, {"Switch for applying neural network on the muon mass hypothesis (if network enabled) (set to 0 to disable)"}};
  o2::framework::Configurable<int> useNetworkPi{"useNetworkPi", 1, {"Switch for applying neural network on the pion mass hypothesis (if network enabled) (set to 0 to disable)"}};
  o2::framework::Configurable<int> useNetworkKa{"useNetworkKa", 1, {"Switch for applying neural network on the kaon mass hypothesis (if network enabled) (set to 0 to disable)"}};
  o2::framework::Configurable<int> useNetworkPr{"useNetworkPr", 1, {"Switch for applying neural network on the proton mass hypothesis (if network enabled) (set to 0 to disable)"}};
  o2::framework::Configurable<int> useNetworkDe{"useNetworkDe", 1, {"Switch for applying neural network on the deuteron mass hypothesis (if network enabled) (set to 0 to disable)"}};
  o2::framework::Configurable<int> useNetworkTr{"useNetworkTr", 1, {"Switch for applying neural network on the triton mass hypothesis (if network enabled) (set to 0 to disable)"}};
  o2::framework::Configurable<int> useNetworkHe{"useNetworkHe", 1, {"Switch for applying neural network on the helium3 mass hypothesis (if network enabled) (set to 0 to disable)"}};
  o2::framework::Configurable<int> useNetworkAl{"useNetworkAl", 1, {"Switch for applying neural network on the alpha mass hypothesis (if network enabled) (set to 0 to disable)"}};
  o2::framework::Configurable<float> networkBetaGammaCutoff{"networkBetaGammaCutoff", 0.45, {"Lower value of beta-gamma to override the NN application"}};
};

// helper getter - FIXME should be separate
int getPIDIndex(const int pdgCode) // Get O2 PID index corresponding to MC PDG code
{
  switch (abs(pdgCode)) {
    case 11:
      return o2::track::PID::Electron;
    case 13:
      return o2::track::PID::Muon;
    case 211:
      return o2::track::PID::Pion;
    case 321:
      return o2::track::PID::Kaon;
    case 2212:
      return o2::track::PID::Proton;
    case 1000010020:
      return o2::track::PID::Deuteron;
    case 1000010030:
      return o2::track::PID::Triton;
    case 1000020030:
      return o2::track::PID::Helium3;
    case 1000020040:
      return o2::track::PID::Alpha;
    default: // treat as pion if not any of the above
      return o2::track::PID::Pion;
  }
}

class pidTPCModule
{
 public:
  pidTPCModule()
  {
    // constructor
  }
  o2::aod::pid::pidTPCConfigurables pidTPCopts;

  // TPC PID Response
  o2::pid::tpc::Response* response;

  // Network correction for TPC PID response
  ml::OnnxModel network;
  std::map<std::string, std::string> metadata;
  std::map<std::string, std::string> nullmetadata;
  std::map<std::string, std::string> headers;
  std::vector<int> speciesNetworkFlags = std::vector<int>(9);
  std::string networkVersion;

  // Parametrization configuration
  bool useCCDBParam = false;

  //__________________________________________________
  template <typename TCCDB, typename TCCDBApi, typename TContext, typename TpidTPCOpts, typename TMetadataInfo>
  void init(TCCDB& ccdb, TCCDBApi& ccdbApi, TContext& context, TpidTPCOpts const& external_pidtpcopts, TMetadataInfo const& metadataInfo)
  {
    // read in configurations from the task where it's used
    pidTPCopts = external_pidtpcopts;

    // initialize PID response
    response = new o2::pid::tpc::Response();

    enableFlagIfTableRequired(context, "DEdxsCorrected", pidTPCopts.savedEdxsCorrected);

    // Checking the tables are requested in the workflow and enabling them
    auto enableFlag = [&](const std::string particle, o2::framework::Configurable<int>& flag) {
      enableFlagIfTableRequired(context, "pidTPC" + particle, flag);
    };
    enableFlag("FullEl", pidTPCopts.pidFullEl);
    enableFlag("FullMu", pidTPCopts.pidFullMu);
    enableFlag("FullPi", pidTPCopts.pidFullPi);
    enableFlag("FullKa", pidTPCopts.pidFullKa);
    enableFlag("FullPr", pidTPCopts.pidFullPr);
    enableFlag("FullDe", pidTPCopts.pidFullDe);
    enableFlag("FullTr", pidTPCopts.pidFullTr);
    enableFlag("FullHe", pidTPCopts.pidFullHe);
    enableFlag("FullAl", pidTPCopts.pidFullAl);

    enableFlag("El", pidTPCopts.pidTinyEl);
    enableFlag("Mu", pidTPCopts.pidTinyMu);
    enableFlag("Pi", pidTPCopts.pidTinyPi);
    enableFlag("Ka", pidTPCopts.pidTinyKa);
    enableFlag("Pr", pidTPCopts.pidTinyPr);
    enableFlag("De", pidTPCopts.pidTinyDe);
    enableFlag("Tr", pidTPCopts.pidTinyTr);
    enableFlag("He", pidTPCopts.pidTinyHe);
    enableFlag("Al", pidTPCopts.pidTinyAl);

    if (metadataInfo.isMC()) {
      enableFlagIfTableRequired(context, "mcTPCTuneOnData", pidTPCopts.enableTuneOnDataTable);
    }

    speciesNetworkFlags[0] = pidTPCopts.useNetworkEl;
    speciesNetworkFlags[1] = pidTPCopts.useNetworkMu;
    speciesNetworkFlags[2] = pidTPCopts.useNetworkPi;
    speciesNetworkFlags[3] = pidTPCopts.useNetworkKa;
    speciesNetworkFlags[4] = pidTPCopts.useNetworkPr;
    speciesNetworkFlags[5] = pidTPCopts.useNetworkDe;
    speciesNetworkFlags[6] = pidTPCopts.useNetworkTr;
    speciesNetworkFlags[7] = pidTPCopts.useNetworkHe;
    speciesNetworkFlags[8] = pidTPCopts.useNetworkAl;

    // Initialise metadata object for CCDB calls from AO2D metadata
    if (pidTPCopts.recoPass.value == "") {
      if (metadataInfo.isFullyDefined()) {
        metadata["RecoPassName"] = metadataInfo.get("RecoPassName");
        LOGP(info, "Automatically setting reco pass for TPC Response to {} from AO2D", metadata["RecoPassName"]);
      }
    } else {
      LOGP(info, "Setting reco pass for TPC response to user-defined name {}", pidTPCopts.recoPass.value);
      metadata["RecoPassName"] = pidTPCopts.recoPass.value;
    }

    /// TPC PID Response
    const TString fname = pidTPCopts.paramfile.value;
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
      const std::string path = pidTPCopts.ccdbPath.value;
      const auto time = pidTPCopts.ccdbTimestamp.value;
      if (time != 0) {
        LOGP(info, "Initialising TPC PID response for fixed timestamp {} and reco pass {}:", time, pidTPCopts.recoPass.value);
        ccdb->setTimestamp(time);
        response = ccdb->template getSpecific<o2::pid::tpc::Response>(path, time, metadata);
        headers = ccdbApi.retrieveHeaders(path, metadata, time);
        if (!response) {
          LOGF(warning, "Unable to find TPC parametrisation for specified pass name - falling back to latest object");
          response = ccdb->template getForTimeStamp<o2::pid::tpc::Response>(path, time);
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
    if (!pidTPCopts.useNetworkCorrection) {
      return;
    } else {
      /// CCDB and auto-fetching
      if (!pidTPCopts.autofetchNetworks) {
        if (pidTPCopts.ccdbTimestamp > 0) {
          /// Fetching network for specific timestamp
          LOG(info) << "Fetching network for timestamp: " << pidTPCopts.ccdbTimestamp.value;
          bool retrieveSuccess = ccdbApi.retrieveBlob(pidTPCopts.networkPathCCDB.value, ".", metadata, pidTPCopts.ccdbTimestamp.value, false, pidTPCopts.networkPathLocally.value);
          headers = ccdbApi.retrieveHeaders(pidTPCopts.networkPathCCDB.value, metadata, pidTPCopts.ccdbTimestamp.value);
          networkVersion = headers["NN-Version"];
          if (retrieveSuccess) {
            network.initModel(pidTPCopts.networkPathLocally.value, pidTPCopts.enableNetworkOptimizations.value, pidTPCopts.networkSetNumThreads.value, strtoul(headers["Valid-From"].c_str(), NULL, 0), strtoul(headers["Valid-Until"].c_str(), NULL, 0));
            std::vector<float> dummyInput(network.getNumInputNodes(), 1.);
            network.evalModel(dummyInput); /// Init the model evaluations
            LOGP(info, "Retrieved NN corrections for production tag {}, pass number {}, and NN-Version {}", headers["LPMProductionTag"], headers["RecoPassName"], headers["NN-Version"]);
          } else {
            LOG(fatal) << "No valid NN object found matching retrieved Bethe-Bloch parametrisation for pass " << metadata["RecoPassName"] << ". Please ensure that the requested pass has dedicated NN corrections available";
          }
        } else {
          /// Taking the network from local file
          if (pidTPCopts.networkPathLocally.value == "") {
            LOG(fatal) << "Local path must be set (flag networkPathLocally)! Aborting...";
          }
          LOG(info) << "Using local file [" << pidTPCopts.networkPathLocally.value << "] for the TPC PID response correction.";
          network.initModel(pidTPCopts.networkPathLocally.value, pidTPCopts.enableNetworkOptimizations.value, pidTPCopts.networkSetNumThreads.value);
          std::vector<float> dummyInput(network.getNumInputNodes(), 1.);
          network.evalModel(dummyInput); // This is an initialisation and might reduce the overhead of the model
        }
      } else {
        return;
      }
    }
  } // end init 

  //__________________________________________________
  template <typename TCCDB, typename TCCDBApi, typename C, typename M, typename T, typename B>
  std::vector<float> createNetworkPrediction(TCCDB& ccdb, TCCDBApi& ccdbApi, C const& collisions, M const& mults, T const& tracks, B const& bcs, const size_t size)
  {

    std::vector<float> network_prediction;

    auto start_network_total = std::chrono::high_resolution_clock::now();
    if (pidTPCopts.autofetchNetworks) {
      const auto& bc = bcs.begin();
      // Initialise correct TPC response object before NN setup (for NCl normalisation)
      if (useCCDBParam && pidTPCopts.ccdbTimestamp.value == 0 && !ccdb->isCachedObjectValid(pidTPCopts.ccdbPath.value, bc.timestamp())) { // Updating parametrisation only if the initial timestamp is 0
        if (pidTPCopts.recoPass.value == "") {
          LOGP(info, "Retrieving latest TPC response object for timestamp {}:", bc.timestamp());
        } else {
          LOGP(info, "Retrieving TPC Response for timestamp {} and recoPass {}:", bc.timestamp(), pidTPCopts.recoPass.value);
        }
        response = ccdb->template getSpecific<o2::pid::tpc::Response>(pidTPCopts.ccdbPath.value, bc.timestamp(), metadata);
        headers = ccdbApi.retrieveHeaders(pidTPCopts.ccdbPath.value, metadata, bc.timestamp());
        networkVersion = headers["NN-Version"];
        if (!response) {
          LOGP(warning, "!! Could not find a valid TPC response object for specific pass name {}! Falling back to latest uploaded object.", metadata["RecoPassName"]);
          headers = ccdbApi.retrieveHeaders(pidTPCopts.ccdbPath.value, nullmetadata, bc.timestamp());
          response = ccdb->template getForTimeStamp<o2::pid::tpc::Response>(pidTPCopts.ccdbPath.value, bc.timestamp());
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
        bool retrieveSuccess = ccdbApi.retrieveBlob(pidTPCopts.networkPathCCDB.value, ".", metadata, bc.timestamp(), false, pidTPCopts.networkPathLocally.value);
        headers = ccdbApi.retrieveHeaders(pidTPCopts.networkPathCCDB.value, metadata, bc.timestamp());
        networkVersion = headers["NN-Version"];
        if (retrieveSuccess) {
          network.initModel(pidTPCopts.networkPathLocally.value, pidTPCopts.enableNetworkOptimizations.value, pidTPCopts.networkSetNumThreads.value, strtoul(headers["Valid-From"].c_str(), NULL, 0), strtoul(headers["Valid-Until"].c_str(), NULL, 0));
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
        if (pidTPCopts.skipTPCOnly) {
          if (!trk.hasITS() && !trk.hasTRD() && !trk.hasTOF()) {
            continue;
          }
        }
        track_properties[counter_track_props] = trk.tpcInnerParam();
        track_properties[counter_track_props + 1] = trk.tgl();
        track_properties[counter_track_props + 2] = trk.signed1Pt();
        track_properties[counter_track_props + 3] = o2::track::pid_constants::sMasses[i];
        track_properties[counter_track_props + 4] = trk.has_collision() ? mults[trk.collisionId()] / 11000. : 1.;
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

  //__________________________________________________
  template <typename C, typename T, typename NSF, typename NST>
  void makePidTables(const int flagFull, NSF& tableFull, const int flagTiny, NST& tableTiny, const o2::track::PID::ID pid, const float tpcSignal, const T& trk, const C& collisions, const long multTPC, const std::vector<float>& network_prediction, const int& count_tracks, const int& tracksForNet_size)
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
    if (pidTPCopts.skipTPCOnly) {
      if (!trk.hasITS() && !trk.hasTRD() && !trk.hasTOF()) {
        if (flagFull)
          tableFull(-999.f, -999.f);
        if (flagTiny)
          tableTiny(aod::pidtpc_tiny::binning::underflowBin);
        return;
      }
    }
    auto expSignal = response->GetExpectedSignal(trk, pid);
    auto expSigma = trk.has_collision() ? response->GetExpectedSigma(collisions.iteratorAt(trk.collisionId()), multTPC, trk, pid) : 0.07 * expSignal; // use default sigma value of 7% if no collision information to estimate resolution
    if (expSignal < 0. || expSigma < 0.) {                                                                                                   // skip if expected signal invalid
      if (flagFull)
        tableFull(-999.f, -999.f);
      if (flagTiny)
        tableTiny(aod::pidtpc_tiny::binning::underflowBin);
      return;
    }

    float nSigma = -999.f;
    float bg = trk.tpcInnerParam() / o2::track::pid_constants::sMasses[pid]; // estimated beta-gamma for network cutoff
    if (pidTPCopts.useNetworkCorrection && speciesNetworkFlags[pid] && trk.has_collision() && bg > pidTPCopts.networkBetaGammaCutoff) {

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
      nSigma = response->GetNumberOfSigmaMCTuned(collisions.iteratorAt(trk.collisionId()), multTPC, trk, pid, tpcSignal);
    }
    if (flagFull)
      tableFull(expSigma, nSigma);
    if (flagTiny)
      aod::pidtpc_tiny::binning::packInTable(nSigma, tableTiny);
  };

  //__________________________________________________
  template <typename TCCDB, typename TCCDBApi, typename TBCs, typename TCollisions, typename TTracks, typename TProducts>
  void process(TCCDB& ccdb, TCCDBApi& ccdbApi, TBCs const& bcs, TCollisions const& cols, TTracks const& tracks, TProducts& products)
  {
    // preparatory step: we need the multiplicities for each collision
    std::vector<int> pidmults; 
    long totalTPCtracks = 0;
    long totalTPCnotStandalone = 0;
    pidmults.resize(cols.size(), 0);

    // faster counting
    for (const auto& track : tracks) {
      if(track.hasTPC()){ 
        pidmults[track.collisionId()]++;
        totalTPCtracks++;
        if(track.hasITS()||track.hasTOF()||track.hasTRD()){
          totalTPCnotStandalone++;
        }
      }
    }

    const uint64_t outTable_size = tracks.size();

    auto reserveTable = [&outTable_size](const o2::framework::Configurable<int>& flag, auto& table) {
      if (flag.value != 1) {
        return;
      }
      table.reserve(outTable_size);
    };

    // Prepare memory for enabled tables
    reserveTable(pidTPCopts.pidFullEl, products.tablePIDFullEl);
    reserveTable(pidTPCopts.pidFullMu, products.tablePIDFullMu);
    reserveTable(pidTPCopts.pidFullPi, products.tablePIDFullPi);
    reserveTable(pidTPCopts.pidFullKa, products.tablePIDFullKa);
    reserveTable(pidTPCopts.pidFullPr, products.tablePIDFullPr);
    reserveTable(pidTPCopts.pidFullDe, products.tablePIDFullDe);
    reserveTable(pidTPCopts.pidFullTr, products.tablePIDFullTr);
    reserveTable(pidTPCopts.pidFullHe, products.tablePIDFullHe);
    reserveTable(pidTPCopts.pidFullAl, products.tablePIDFullAl);

    reserveTable(pidTPCopts.pidTinyEl, products.tablePIDTinyEl);
    reserveTable(pidTPCopts.pidTinyMu, products.tablePIDTinyMu);
    reserveTable(pidTPCopts.pidTinyPi, products.tablePIDTinyPi);
    reserveTable(pidTPCopts.pidTinyKa, products.tablePIDTinyKa);
    reserveTable(pidTPCopts.pidTinyPr, products.tablePIDTinyPr);
    reserveTable(pidTPCopts.pidTinyDe, products.tablePIDTinyDe);
    reserveTable(pidTPCopts.pidTinyTr, products.tablePIDTinyTr);
    reserveTable(pidTPCopts.pidTinyHe, products.tablePIDTinyHe);
    reserveTable(pidTPCopts.pidTinyAl, products.tablePIDTinyAl);

    const uint64_t tracksForNet_size = (pidTPCopts.skipTPCOnly) ? totalTPCnotStandalone : totalTPCtracks;
    std::vector<float> network_prediction;

    if (pidTPCopts.useNetworkCorrection) {
      network_prediction = createNetworkPrediction(ccdb, ccdbApi, cols, pidmults, tracks, bcs, tracksForNet_size);
    }

    uint64_t count_tracks = 0;

    for (auto const& trk : tracks) {
      // Loop on Tracks

      const auto& bc = trk.has_collision() ? cols.iteratorAt(trk.collisionId()).template bc_as<aod::BCsWithTimestamps>() : bcs.begin();
      if (useCCDBParam && pidTPCopts.ccdbTimestamp.value == 0 && !ccdb->isCachedObjectValid(pidTPCopts.ccdbPath.value, bc.timestamp())) { // Updating parametrisation only if the initial timestamp is 0
        if (pidTPCopts.recoPass.value == "") {
          LOGP(info, "Retrieving latest TPC response object for timestamp {}:", bc.timestamp());
        } else {
          LOGP(info, "Retrieving TPC Response for timestamp {} and recoPass {}:", bc.timestamp(), pidTPCopts.recoPass.value);
        }
        response = ccdb->template getSpecific<o2::pid::tpc::Response>(pidTPCopts.ccdbPath.value, bc.timestamp(), metadata);
        headers = ccdbApi.retrieveHeaders(pidTPCopts.ccdbPath.value, metadata, bc.timestamp());
        if (!response) {
          LOGP(warning, "!! Could not find a valid TPC response object for specific pass name {}! Falling back to latest uploaded object.", metadata["RecoPassName"]);
          response = ccdb->template getForTimeStamp<o2::pid::tpc::Response>(pidTPCopts.ccdbPath.value, bc.timestamp());
          headers = ccdbApi.retrieveHeaders(pidTPCopts.ccdbPath.value, nullmetadata, bc.timestamp());
          if (!response) {
            LOGP(fatal, "Could not find ANY TPC response object for the timestamp {}!", bc.timestamp());
          }
        }
        LOG(info) << "Successfully retrieved TPC PID object from CCDB for timestamp " << bc.timestamp() << ", period " << headers["LPMProductionTag"] << ", recoPass " << headers["RecoPassName"];
        response->PrintAll();
      }

      auto makePidTablesDefault = [&trk, &cols, &pidmults, &network_prediction, &count_tracks, &tracksForNet_size, this](const int flagFull, auto& tableFull, const int flagTiny, auto& tableTiny, const o2::track::PID::ID pid) {
        makePidTables(flagFull, tableFull, flagTiny, tableTiny, pid, trk.tpcSignal(), trk, cols, pidmults[trk.collisionId()], network_prediction, count_tracks, tracksForNet_size);
      };

      makePidTablesDefault(pidTPCopts.pidFullEl, products.tablePIDFullEl, pidTPCopts.pidTinyEl, products.tablePIDTinyEl, o2::track::PID::Electron);
      makePidTablesDefault(pidTPCopts.pidFullMu, products.tablePIDFullMu, pidTPCopts.pidTinyMu, products.tablePIDTinyMu, o2::track::PID::Muon);
      makePidTablesDefault(pidTPCopts.pidFullPi, products.tablePIDFullPi, pidTPCopts.pidTinyPi, products.tablePIDTinyPi, o2::track::PID::Pion);
      makePidTablesDefault(pidTPCopts.pidFullKa, products.tablePIDFullKa, pidTPCopts.pidTinyKa, products.tablePIDTinyKa, o2::track::PID::Kaon);
      makePidTablesDefault(pidTPCopts.pidFullPr, products.tablePIDFullPr, pidTPCopts.pidTinyPr, products.tablePIDTinyPr, o2::track::PID::Proton);
      makePidTablesDefault(pidTPCopts.pidFullDe, products.tablePIDFullDe, pidTPCopts.pidTinyDe, products.tablePIDTinyDe, o2::track::PID::Deuteron);
      makePidTablesDefault(pidTPCopts.pidFullTr, products.tablePIDFullTr, pidTPCopts.pidTinyTr, products.tablePIDTinyTr, o2::track::PID::Triton);
      makePidTablesDefault(pidTPCopts.pidFullHe, products.tablePIDFullHe, pidTPCopts.pidTinyHe, products.tablePIDTinyHe, o2::track::PID::Helium3);
      makePidTablesDefault(pidTPCopts.pidFullAl, products.tablePIDFullAl, pidTPCopts.pidTinyAl, products.tablePIDTinyAl, o2::track::PID::Alpha);

      if (trk.hasTPC() && (!pidTPCopts.skipTPCOnly || trk.hasITS() || trk.hasTRD() || trk.hasTOF())) {
        count_tracks++; // Increment network track counter only if track has TPC, and (not skipping TPConly) or (is not TPConly)
      }
    }
  } // end process function
};

} // namespace pid
} // namespace o2::aod

#endif // COMMON_TOOLS_PIDTPCMODULE_H_
