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
/// \file   pidTPCFull.cxx
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \author Christian Sonnabend christian.sonnabend@cern.ch
/// \author Annalena Kalteyer annalena.sophie.kalteyer@cern.ch
/// \brief  Task to produce PID tables for TPC split for each particle.
///         Only the tables for the mass hypotheses requested are filled, the others are sent empty.
///         QA histograms for the TPC PID can be produced by adding `--add-qa 1` to the workflow
///

// ROOT includes
#include "TFile.h"
#include "TSystem.h"

// O2 includes
#include <CCDB/BasicCCDBManager.h>
#include "Framework/AnalysisTask.h"
#include "ReconstructionDataFormats/Track.h"
#include "CCDB/CcdbApi.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/Core/PID/TPCPIDResponse.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/DataModel/Multiplicity.h"
#include "TableHelper.h"
#include "Tools/ML/model.h"
#include "pidTPCBase.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::pid;
using namespace o2::pid::tpc;
using namespace o2::framework::expressions;
using namespace o2::track;
using namespace o2::ml;

void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  std::vector<ConfigParamSpec> options{{"add-qa", VariantType::Int, 0, {"Legacy. No effect."}}};
  std::swap(workflowOptions, options);
}

#include "Framework/runDataProcessing.h"

/// Task to produce the response table
struct tpcPidFull {
  using Trks = soa::Join<aod::Tracks, aod::TracksExtra>;
  using Coll = soa::Join<aod::Collisions, aod::PIDMults>;

  // Tables to produce
  Produces<o2::aod::pidTPCFullEl> tablePIDEl;
  Produces<o2::aod::pidTPCFullMu> tablePIDMu;
  Produces<o2::aod::pidTPCFullPi> tablePIDPi;
  Produces<o2::aod::pidTPCFullKa> tablePIDKa;
  Produces<o2::aod::pidTPCFullPr> tablePIDPr;
  Produces<o2::aod::pidTPCFullDe> tablePIDDe;
  Produces<o2::aod::pidTPCFullTr> tablePIDTr;
  Produces<o2::aod::pidTPCFullHe> tablePIDHe;
  Produces<o2::aod::pidTPCFullAl> tablePIDAl;

  // TPC PID Response
  o2::pid::tpc::Response* response;

  // Network correction for TPC PID response
  OnnxModel network;
  o2::ccdb::CcdbApi ccdbApi;
  std::map<std::string, std::string> metadata;
  std::map<std::string, std::string> headers;

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
  Configurable<int> pidEl{"pid-el", -1, {"Produce PID information for the Electron mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  Configurable<int> pidMu{"pid-mu", -1, {"Produce PID information for the Muon mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  Configurable<int> pidPi{"pid-pi", -1, {"Produce PID information for the Pion mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  Configurable<int> pidKa{"pid-ka", -1, {"Produce PID information for the Kaon mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  Configurable<int> pidPr{"pid-pr", -1, {"Produce PID information for the Proton mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  Configurable<int> pidDe{"pid-de", -1, {"Produce PID information for the Deuterons mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  Configurable<int> pidTr{"pid-tr", -1, {"Produce PID information for the Triton mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  Configurable<int> pidHe{"pid-he", -1, {"Produce PID information for the Helium3 mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};
  Configurable<int> pidAl{"pid-al", -1, {"Produce PID information for the Alpha mass hypothesis, overrides the automatic setup: the corresponding table can be set off (0) or on (1)"}};

  // Paramatrization configuration
  bool useCCDBParam = false;

  void init(o2::framework::InitContext& initContext)
  {
    response = new o2::pid::tpc::Response();
    // Checking the tables are requested in the workflow and enabling them
    auto enableFlag = [&](const std::string particle, Configurable<int>& flag) {
      enableFlagIfTableRequired(initContext, "pidTPCFull" + particle, flag);
    };
    enableFlag("El", pidEl);
    enableFlag("Mu", pidMu);
    enableFlag("Pi", pidPi);
    enableFlag("Ka", pidKa);
    enableFlag("Pr", pidPr);
    enableFlag("De", pidDe);
    enableFlag("Tr", pidTr);
    enableFlag("He", pidHe);
    enableFlag("Al", pidAl);

    // Initialise metadata object for CCDB calls
    if (recoPass.value == "") {
      LOGP(info, "Reco pass not specified; CCDB will take latest available object");
    } else {
      LOGP(info, "CCDB object will be requested for reconstruction pass {}", recoPass.value);
      metadata["RecoPassName"] = recoPass.value;
    }

    /// TPC PID Response
    const TString fname = paramfile.value;
    if (fname != "") { // Loading the parametrization from file
      LOGP(info, "Loading TPC response from file {}", fname);
      try {
        std::unique_ptr<TFile> f(TFile::Open(fname, "READ"));
        f->GetObject("Response", response);
      } catch (...) {
        LOGF(fatal, "Loading the TPC PID Response from file {} failed!", fname);
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
      if (time != 0) {
        LOGP(info, "Initialising TPC PID response for fixed timestamp {} and reco pass {}:", time, recoPass.value);
        ccdb->setTimestamp(time);
        response = ccdb->getSpecific<o2::pid::tpc::Response>(path, time, metadata);
        if (!response) {
          LOGF(warning, "Unable to find TPC parametrisation for specified pass name - falling back to latest object");
          response = ccdb->getForTimeStamp<o2::pid::tpc::Response>(path, time);
          if (!response) {
            LOGF(fatal, "Unable to find any TPC object corresponding to timestamp {}!", time);
          }
        }
        response->PrintAll();
      }
    }

    /// Neural network init for TPC PID

    if (!useNetworkCorrection) {
      return;
    } else {

      /// CCDB and auto-fetching
      ccdbApi.init(url);
      if (!autofetchNetworks) {
        if (ccdbTimestamp > 0) {
          /// Fetching network for specific timestamp
          LOG(info) << "Fetching network for timestamp: " << ccdbTimestamp.value;
          bool retrieveSuccess = ccdbApi.retrieveBlob(networkPathCCDB.value, ".", metadata, ccdbTimestamp.value, false, networkPathLocally.value);
          headers = ccdbApi.retrieveHeaders(networkPathCCDB.value, metadata, ccdbTimestamp.value);
          if (retrieveSuccess) {
            network.initModel(networkPathLocally.value, enableNetworkOptimizations.value, networkSetNumThreads.value, strtoul(headers["Valid-From"].c_str(), NULL, 0), strtoul(headers["Valid-Until"].c_str(), NULL, 0));
            std::vector<float> dummyInput(network.getNumInputNodes(), 1.);
            network.evalModel(dummyInput); /// Init the model evaluations
          } else {
            LOG(fatal) << "Error encountered while fetching/loading the network from CCDB! Maybe the network doesn't exist yet for this runnumber/timestamp?";
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

  Partition<Trks> notTPCStandaloneTracks = (aod::track::tpcNClsFindable > (uint8_t)0) && ((aod::track::itsClusterSizes > (uint32_t)0) || (aod::track::trdPattern > (uint8_t)0) || (aod::track::tofExpMom > 0.f && aod::track::tofChi2 > 0.f)); // To count number of tracks for use in NN array
  Partition<Trks> tracksWithTPC = (aod::track::tpcNClsFindable > (uint8_t)0);

  void process(Coll const& collisions, Trks const& tracks,
               aod::BCsWithTimestamps const&)
  {

    const uint64_t outTable_size = tracks.size();
    auto reserveTable = [&outTable_size](const Configurable<int>& flag, auto& table) {
      if (flag.value != 1) {
        return;
      }
      table.reserve(outTable_size);
    };
    // Prepare memory for enabled tables
    reserveTable(pidEl, tablePIDEl);
    reserveTable(pidMu, tablePIDMu);
    reserveTable(pidPi, tablePIDPi);
    reserveTable(pidKa, tablePIDKa);
    reserveTable(pidPr, tablePIDPr);
    reserveTable(pidDe, tablePIDDe);
    reserveTable(pidTr, tablePIDTr);
    reserveTable(pidHe, tablePIDHe);
    reserveTable(pidAl, tablePIDAl);

    std::vector<float> network_prediction;
    const uint64_t tracksForNet_size = (skipTPCOnly) ? notTPCStandaloneTracks.size() : tracksWithTPC.size();

    if (useNetworkCorrection) {

      auto start_network_total = std::chrono::high_resolution_clock::now();
      if (autofetchNetworks) {
        auto bc = collisions.iteratorAt(0).bc_as<aod::BCsWithTimestamps>();
        // Initialise correct TPC response object before NN setup (for NCl normalisation)
        if (useCCDBParam && ccdbTimestamp.value == 0 && !ccdb->isCachedObjectValid(ccdbPath.value, bc.timestamp())) { // Updating parametrisation only if the initial timestamp is 0
          if (recoPass.value == "") {
            LOGP(info, "Retrieving latest TPC response object for timestamp {}:", bc.timestamp());
          } else {
            LOGP(info, "Retrieving TPC Response for timestamp {} and recoPass {}:", bc.timestamp(), recoPass.value);
          }
          response = ccdb->getSpecific<o2::pid::tpc::Response>(ccdbPath.value, bc.timestamp(), metadata);
          if (!response) {
            LOGP(warning, "!! Could not find a valid TPC response object for specific pass name {}! Falling back to latest uploaded object.", recoPass.value);
            response = ccdb->getForTimeStamp<o2::pid::tpc::Response>(ccdbPath.value, bc.timestamp());
            if (!response) {
              LOGP(fatal, "Could not find ANY TPC response object for the timestamp {}!", bc.timestamp());
            }
          }
          response->PrintAll();
        }

        if (bc.timestamp() < network.getValidityFrom() || bc.timestamp() > network.getValidityUntil()) { // fetches network only if the runnumbers change
          LOG(info) << "Fetching network for timestamp: " << bc.timestamp();
          bool retrieveSuccess = ccdbApi.retrieveBlob(networkPathCCDB.value, ".", metadata, bc.timestamp(), false, networkPathLocally.value);
          headers = ccdbApi.retrieveHeaders(networkPathCCDB.value, metadata, bc.timestamp());
          if (retrieveSuccess) {
            network.initModel(networkPathLocally.value, enableNetworkOptimizations.value, networkSetNumThreads.value, strtoul(headers["Valid-From"].c_str(), NULL, 0), strtoul(headers["Valid-Until"].c_str(), NULL, 0));
            std::vector<float> dummyInput(network.getNumInputNodes(), 1.);
            network.evalModel(dummyInput);
          } else {
            LOG(fatal) << "Error encountered while fetching/loading the network from CCDB! Maybe the network doesn't exist yet for this runnumber/timestamp?";
          }
        }
      }

      // Defining some network parameters
      int input_dimensions = network.getNumInputNodes();
      int output_dimensions = network.getNumOutputNodes();
      const uint64_t track_prop_size = input_dimensions * tracksForNet_size;
      const uint64_t prediction_size = output_dimensions * tracksForNet_size;

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
          track_properties[counter_track_props + 4] = collisions.iteratorAt(trk.collisionId()).multTPC() / 11000.;
          track_properties[counter_track_props + 5] = std::sqrt(nNclNormalization / trk.tpcNClsFound());
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
      LOG(debug) << "Neural Network for the TPC PID response correction: Time per track (eval ONNX): " << duration_network / (tracksForNet_size * 9) << "ns ; Total time (eval ONNX): " << duration_network / 1000000000 << " s";
      LOG(debug) << "Neural Network for the TPC PID response correction: Time per track (eval + overhead): " << std::chrono::duration<float, std::ratio<1, 1000000000>>(stop_network_total - start_network_total).count() / (tracksForNet_size * 9) << "ns ; Total time (eval + overhead): " << std::chrono::duration<float, std::ratio<1, 1000000000>>(stop_network_total - start_network_total).count() / 1000000000 << " s";
    }

    uint64_t count_tracks = 0;

    for (auto const& trk : tracks) {
      // Loop on Tracks
      if (trk.has_collision()) {
        const auto& bc = collisions.iteratorAt(trk.collisionId()).bc_as<aod::BCsWithTimestamps>();
        if (useCCDBParam && ccdbTimestamp.value == 0 && !ccdb->isCachedObjectValid(ccdbPath.value, bc.timestamp())) { // Updating parametrisation only if the initial timestamp is 0
          if (recoPass.value == "") {
            LOGP(info, "Retrieving latest TPC response object for timestamp {}:", bc.timestamp());
          } else {
            LOGP(info, "Retrieving TPC Response for timestamp {} and recoPass {}:", bc.timestamp(), recoPass.value);
          }
          response = ccdb->getSpecific<o2::pid::tpc::Response>(ccdbPath.value, bc.timestamp(), metadata);
          if (!response) {
            LOGP(warning, "!! Could not find a valid TPC response object for specific pass name {}! Falling back to latest uploaded object.", recoPass.value);
            response = ccdb->getForTimeStamp<o2::pid::tpc::Response>(ccdbPath.value, bc.timestamp());
            if (!response) {
              LOGP(fatal, "Could not find ANY TPC response object for the timestamp {}!", bc.timestamp());
            }
          }
          response->PrintAll();
        }
      }
      // Check and fill enabled tables
      auto makeTable = [&trk, &collisions, &network_prediction, &count_tracks, &tracksForNet_size, this](const Configurable<int>& flag, auto& table, const o2::track::PID::ID pid) {
        if (flag.value != 1) {
          return;
        }
        if (!trk.hasTPC()) {
          table(-999.f, -999.f);
          return;
        }
        if (skipTPCOnly) {
          if (!trk.hasITS() && !trk.hasTRD() && !trk.hasTOF()) {
            table(-999.f, -999.f);
            return;
          }
        }
        auto expSignal = response->GetExpectedSignal(trk, pid);
        auto expSigma = response->GetExpectedSigma(collisions.iteratorAt(trk.collisionId()), trk, pid);
        if (expSignal < 0. || expSigma < 0.) { // skip if expected signal invalid
          table(-999.f, -999.f);
          return;
        }

        if (useNetworkCorrection) {

          // Here comes the application of the network. The output--dimensions of the network dtermine the application: 1: mean, 2: sigma, 3: sigma asymmetric
          // For now only the option 2: sigma will be used. The other options are kept if there would be demand later on
          if (network.getNumOutputNodes() == 1) {
            table(expSigma,
                  (trk.tpcSignal() - network_prediction[count_tracks + tracksForNet_size * pid] * expSignal) / expSigma);
          } else if (network.getNumOutputNodes() == 2) {
            table((network_prediction[2 * (count_tracks + tracksForNet_size * pid) + 1] - network_prediction[2 * (count_tracks + tracksForNet_size * pid)]) * expSignal,
                  (trk.tpcSignal() / expSignal - network_prediction[2 * (count_tracks + tracksForNet_size * pid)]) / (network_prediction[2 * (count_tracks + tracksForNet_size * pid) + 1] - network_prediction[2 * (count_tracks + tracksForNet_size * pid)]));
          } else if (network.getNumOutputNodes() == 3) {
            if (trk.tpcSignal() / expSignal >= network_prediction[3 * (count_tracks + tracksForNet_size * pid)]) {
              table((network_prediction[3 * (count_tracks + tracksForNet_size * pid) + 1] - network_prediction[3 * (count_tracks + tracksForNet_size * pid)]) * expSignal,
                    (trk.tpcSignal() / expSignal - network_prediction[3 * (count_tracks + tracksForNet_size * pid)]) / (network_prediction[3 * (count_tracks + tracksForNet_size * pid) + 1] - network_prediction[3 * (count_tracks + tracksForNet_size * pid)]));
            } else {
              table((network_prediction[3 * (count_tracks + tracksForNet_size * pid)] - network_prediction[3 * (count_tracks + tracksForNet_size * pid) + 2]) * expSignal,
                    (trk.tpcSignal() / expSignal - network_prediction[3 * (count_tracks + tracksForNet_size * pid)]) / (network_prediction[3 * (count_tracks + tracksForNet_size * pid)] - network_prediction[3 * (count_tracks + tracksForNet_size * pid) + 2]));
            }
          } else {
            LOGF(fatal, "Network output-dimensions incompatible!");
          }
        } else {
          table(expSigma,
                response->GetNumberOfSigma(collisions.iteratorAt(trk.collisionId()), trk, pid));
        }
      };

      makeTable(pidEl, tablePIDEl, o2::track::PID::Electron);
      makeTable(pidMu, tablePIDMu, o2::track::PID::Muon);
      makeTable(pidPi, tablePIDPi, o2::track::PID::Pion);
      makeTable(pidKa, tablePIDKa, o2::track::PID::Kaon);
      makeTable(pidPr, tablePIDPr, o2::track::PID::Proton);
      makeTable(pidDe, tablePIDDe, o2::track::PID::Deuteron);
      makeTable(pidTr, tablePIDTr, o2::track::PID::Triton);
      makeTable(pidHe, tablePIDHe, o2::track::PID::Helium3);
      makeTable(pidAl, tablePIDAl, o2::track::PID::Alpha);

      if (trk.hasTPC() && (!skipTPCOnly || trk.hasITS() || trk.hasTRD() || trk.hasTOF())) {
        count_tracks++; // Increment network track counter only if track has TPC, and (not skipping TPConly) or (is not TPConly)
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<tpcPidFull>(cfgc)}; }
