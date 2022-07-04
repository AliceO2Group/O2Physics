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

// O2 includes
#include "Framework/AnalysisTask.h"
#include "ReconstructionDataFormats/Track.h"
#include <CCDB/BasicCCDBManager.h>
#include "Common/DataModel/PIDResponse.h"
#include "Common/Core/PID/TPCPIDResponse.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/DataModel/Multiplicity.h"
#include "TableHelper.h"
#include "Common/TableProducer/PID/pidTPCML.h"
#include "DPG/Tasks/qaPIDTPC.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::pid;
using namespace o2::pid::tpc;
using namespace o2::framework::expressions;
using namespace o2::track;

void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  std::vector<ConfigParamSpec> options{{"add-qa", VariantType::Int, 0, {"Produce TPC PID QA histograms"}}};
  std::swap(workflowOptions, options);
}

#include "Framework/runDataProcessing.h"

/// Task to produce the response table
struct tpcPidFull {
  using Trks = soa::Join<aod::Tracks, aod::TracksExtra>;
  using Coll = soa::Join<aod::Collisions, aod::Mults>;

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
  o2::pid::tpc::Response response;
  o2::pid::tpc::Response* responseptr = nullptr;
  // Network correction for TPC PID response
  Network network;
  std::string currentRunNumber = "-1";

  // Input parameters
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Configurable<std::string> paramfile{"param-file", "", "Path to the parametrization object, if emtpy the parametrization is not taken from file"};
  Configurable<std::string> url{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> ccdbPath{"ccdbPath", "Analysis/PID/TPC/Response", "Path of the TPC parametrization on the CCDB"};
  Configurable<long> ccdbTimestamp{"ccdb-timestamp", 0, "timestamp of the object used to query in CCDB the detector response. Exceptions: -1 gets the latest object, 0 gets the run dependent timestamp"};
  // Parameters for loading network from a file / downloading the file
  Configurable<bool> useNetworkCorrection{"useNetworkCorrection", 0, "(bool) Wether or not to use the network correction for the TPC dE/dx signal"};
  Configurable<bool> downloadNetworkFromAlien{"downloadNetworkFromAlien", 1, "(bool) Download network from AliEn. If true and networkPathAlien==false, then the networks are fetched automatically using the runnumber of the collisions"};
  Configurable<std::string> networkPathLocally{"networkPathLocally", "network.onnx", "(std::string) Path to local .onnx file containing the network, if downloadNetworkFromAlien is True, the network will be downloaded to this directory"};
  Configurable<std::string> networkPathAlien{"networkPathAlien", "", "(std::string) Set this variable if you want to use a network from a specific path on AliEn for the full evaluation, generally: alien:///alice/cern.ch/user/c/csonnabe/TPC/NN_PID/network_[runnumber].onnx. downloadNetworkFromAlien must be set to true."};
  Configurable<std::string> networkSetAlienDir{"networkSetAlienDir", "alien:///alice/cern.ch/user/c/csonnabe/TPC/NN_PID", "(std::string) Set the AliEn directory from which the network-files are fetched. Default: alien:///alice/cern.ch/user/c/csonnabe/TPC/NN_PID"};
  Configurable<bool> enableNetworkOptimizations{"enableNetworkOptimizations", 1, "(bool) If the neural network correction is used, this enables GraphOptimizationLevel::ORT_ENABLE_EXTENDED in the ONNX session"};
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

    const TString fname = paramfile.value;
    if (fname != "") { // Loading the parametrization from file
      LOGP(info, "Loading TPC response from file {}", fname);
      try {
        std::unique_ptr<TFile> f(TFile::Open(fname, "READ"));
        f->GetObject("Response", responseptr);
        response.SetParameters(responseptr);
      } catch (...) {
        LOGF(fatal, "Loading the TPC PID Response from file {} failed!", fname);
      };
    } else {
      useCCDBParam = true;
      const std::string path = ccdbPath.value;
      const auto time = ccdbTimestamp.value;
      ccdb->setURL(url.value);

      ccdb->setCaching(true);
      ccdb->setLocalObjectValidityChecking();
      ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
      if (time != 0) {
        LOGP(info, "Initialising TPC PID response for fixed timestamp {}:", time);
        ccdb->setTimestamp(time);
        response.SetParameters(ccdb->getForTimeStamp<o2::pid::tpc::Response>(path, time));
      } else
        LOGP(info, "Initialising default TPC PID response:");
      response.PrintAll();
    }

    if (!useNetworkCorrection) {
      return;
    } else {
      if (downloadNetworkFromAlien.value && networkPathAlien.value == "") {
        return;
      } else {
        Network temp_net(networkPathLocally.value,
                         downloadNetworkFromAlien.value,
                         networkPathAlien.value,
                         enableNetworkOptimizations.value);
        network = temp_net;

        network.evalNetwork(std::vector<float>(network.getInputDimensions(), 1.)); // This is an initialisation and might reduce the overhead of the model
      }
    }
  }

  void process(Coll const& collisions, Trks const& tracks,
               aod::BCsWithTimestamps const&)
  {
    auto reserveTable = [&tracks](const Configurable<int>& flag, auto& table) {
      if (flag.value != 1) {
        return;
      }
      table.reserve(tracks.size());
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

    if (useNetworkCorrection) {

      auto start_overhead = std::chrono::high_resolution_clock::now();

      if (downloadNetworkFromAlien.value && networkPathAlien.value == "") {

        auto bc = collisions.iteratorAt(0).bc_as<aod::BCsWithTimestamps>();

        if (currentRunNumber.c_str() != std::to_string(bc.runNumber())) { // fetches network only if the runnumbers change
          currentRunNumber = std::to_string(bc.runNumber());
          LOG(info) << "Fetching network for runnumber: " << currentRunNumber;
          Network temp_net(networkPathLocally.value,
                           downloadNetworkFromAlien.value,
                           networkSetAlienDir.value + "/network_" + currentRunNumber + ".onnx",
                           enableNetworkOptimizations.value);
          network = temp_net;

          network.evalNetwork(std::vector<float>(network.getInputDimensions(), 1.)); // This is an initialisation and might reduce the overhead of the model
        }
      }

      // Defining some network parameters
      const int inputDimensions = network.getInputDimensions();
      const int outputDimensions = network.getOutputDimensions();
      const unsigned long track_prop_size = tracks.size() * 9;

      network_prediction = std::vector<float>(outputDimensions * track_prop_size);

      std::vector<float> track_properties(inputDimensions * track_prop_size);
      unsigned long counter_track_props = 0;

      // Filling a std::vector<float> to be evaluated by the network
      // Evaluation on single tracks brings huge overhead: Thus evaluation is done on one large vector
      for (int i = 0; i < 9; i++) { // Loop over particle number for which network correction is used
        for (auto const& trk : tracks) {
          track_properties[counter_track_props] = trk.tpcInnerParam();
          track_properties[counter_track_props + 1] = trk.tgl();
          track_properties[counter_track_props + 2] = trk.signed1Pt();
          track_properties[counter_track_props + 3] = o2::track::pid_constants::sMasses[i];
          track_properties[counter_track_props + 4] = collisions.iteratorAt(trk.collisionId()).multTPC() / 11000.;
          track_properties[counter_track_props + 5] = std::sqrt(159. / trk.tpcNClsFound());
          counter_track_props += inputDimensions;
        }
      }

      auto stop_overhead = std::chrono::high_resolution_clock::now();
      float duration_overhead = std::chrono::duration<float, std::ratio<1, 1000000000>>(stop_overhead - start_overhead).count();
      float time_per_track_overhead = duration_overhead / track_prop_size; // There are n (typically n=7) variables in each track which are being extracted in track_properties. Each network evaluation takes time_per_track_overhead/9 nano-seconds
      LOG(info) << "Neural Network for the TPC PID response: Time per track (overhead): " << time_per_track_overhead << "ns ; Overhead total (all tracks): " << duration_overhead / 1000000000 << "s";

      auto start_network = std::chrono::high_resolution_clock::now();
      float* output_network = network.evalNetwork(track_properties);
      track_properties.clear();
      for (unsigned long i = 0; i < (outputDimensions * track_prop_size); i += outputDimensions) {
        for (int j = 0; j < outputDimensions; j++) {
          network_prediction[i + j] = output_network[i + j];
        }
      }
      auto stop_network = std::chrono::high_resolution_clock::now();
      float duration_network = std::chrono::duration<float, std::ratio<1, 1000000000>>(stop_network - start_network).count();
      float time_per_track_net = duration_network / track_prop_size;
      LOG(info) << "Neural Network for the TPC PID response: Time per track (eval): " << time_per_track_net << "ns ; Evaluation total (all tracks): " << duration_network / 1000000000 << "s"; // The time per track but with 9 particle mass hypotheses: So actual time per track is (time_per_track_net / 9)
    }

    int lastCollisionId = -1; // Last collision ID analysed
    unsigned long count_tracks = 0;
    const int tracks_size = tracks.size();

    for (auto const& trk : tracks) {
      // Loop on Tracks
      if (useCCDBParam && ccdbTimestamp.value == 0 && trk.has_collision() && trk.collisionId() != lastCollisionId) { // Updating parametrization only if the initial timestamp is 0
        lastCollisionId = trk.collisionId();
        const auto& bc = collisions.iteratorAt(trk.collisionId()).bc_as<aod::BCsWithTimestamps>();
        response.SetParameters(ccdb->getForTimeStamp<o2::pid::tpc::Response>(ccdbPath.value, bc.timestamp()));
      }
      // Check and fill enabled tables
      auto makeTable = [&trk, &collisions, &network_prediction, &count_tracks, &tracks_size, this](const Configurable<int>& flag, auto& table, const o2::track::PID::ID pid) {
        if (flag.value != 1) {
          return;
        }

        if (useNetworkCorrection) {

          // Here comes the application of the network. The output--dimensions of the network dtermine the application: 1: mean, 2: sigma, 3: sigma asymmetric
          // For now only the option 2: sigma will be used. The other options are kept if there would be demand later on
          if (network.getOutputDimensions() == 1) {
            table(response.GetExpectedSigma(collisions.iteratorAt(trk.collisionId()), trk, pid),
                  (trk.tpcSignal() - network_prediction[count_tracks + tracks_size * pid] * response.GetExpectedSignal(trk, pid)) / response.GetExpectedSigma(collisions.iteratorAt(trk.collisionId()), trk, pid));
          } else if (network.getOutputDimensions() == 2) {
            table((network_prediction[2 * (count_tracks + tracks_size * pid) + 1] - network_prediction[2 * (count_tracks + tracks_size * pid)]) * response.GetExpectedSignal(trk, pid),
                  (trk.tpcSignal() / response.GetExpectedSignal(trk, pid) - network_prediction[2 * (count_tracks + tracks_size * pid)]) / (network_prediction[2 * (count_tracks + tracks_size * pid) + 1] - network_prediction[2 * (count_tracks + tracks_size * pid)]));
          } else if (network.getOutputDimensions() == 3) {
            if (trk.tpcSignal() / response.GetExpectedSignal(trk, pid) >= network_prediction[3 * (count_tracks + tracks_size * pid)]) {
              table((network_prediction[3 * (count_tracks + tracks_size * pid) + 1] - network_prediction[3 * (count_tracks + tracks_size * pid)]) * response.GetExpectedSignal(trk, pid),
                    (trk.tpcSignal() / response.GetExpectedSignal(trk, pid) - network_prediction[3 * (count_tracks + tracks_size * pid)]) / (network_prediction[3 * (count_tracks + tracks_size * pid) + 1] - network_prediction[3 * (count_tracks + tracks_size * pid)]));
            } else {
              table((network_prediction[3 * (count_tracks + tracks_size * pid)] - network_prediction[3 * (count_tracks + tracks_size * pid) + 2]) * response.GetExpectedSignal(trk, pid),
                    (trk.tpcSignal() / response.GetExpectedSignal(trk, pid) - network_prediction[3 * (count_tracks + tracks_size * pid)]) / (network_prediction[3 * (count_tracks + tracks_size * pid)] - network_prediction[3 * (count_tracks + tracks_size * pid) + 2]));
            }
          } else {
            LOGF(fatal, "Network output-dimensions incompatible!");
          }
        } else {
          table(response.GetExpectedSigma(collisions.iteratorAt(trk.collisionId()), trk, pid),
                response.GetNumberOfSigma(collisions.iteratorAt(trk.collisionId()), trk, pid));
        }
      };

      // const o2::pid::tpc::Response& response;
      makeTable(pidEl, tablePIDEl, o2::track::PID::Electron);
      makeTable(pidMu, tablePIDMu, o2::track::PID::Muon);
      makeTable(pidPi, tablePIDPi, o2::track::PID::Pion);
      makeTable(pidKa, tablePIDKa, o2::track::PID::Kaon);
      makeTable(pidPr, tablePIDPr, o2::track::PID::Proton);
      makeTable(pidDe, tablePIDDe, o2::track::PID::Deuteron);
      makeTable(pidTr, tablePIDTr, o2::track::PID::Triton);
      makeTable(pidHe, tablePIDHe, o2::track::PID::Helium3);
      makeTable(pidAl, tablePIDAl, o2::track::PID::Alpha);

      count_tracks++;
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  auto workflow = WorkflowSpec{adaptAnalysisTask<tpcPidFull>(cfgc)};
  if (cfgc.options().get<int>("add-qa")) {
    workflow.push_back(adaptAnalysisTask<tpcPidQa>(cfgc));
  }
  return workflow;
}