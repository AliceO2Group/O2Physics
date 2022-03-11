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
/// \file   handleParamTPCBetheBloch.cxx
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \since  2020-06-22
/// \brief  A simple tool to produce Bethe Bloch parametrization objects for the TPC PID Response
///

#include "Common/Core/PID/BetheBloch.h"
#include "Common/Core/PID/TPCReso.h"
#include "handleParamBase.h"

using namespace o2::pid::tpc;

bool initOptionsAndParse(bpo::options_description& options, int argc, char* argv[])
{
  options.add_options()(
    "url,u", bpo::value<std::string>()->default_value("http://alice-ccdb.cern.ch"), "URL of the CCDB database e.g. http://ccdb-test.cern.ch:8080 or http://alice-ccdb.cern.ch")(
    "ccdb-path,c", bpo::value<std::string>()->default_value("Analysis/PID/TPC"), "CCDB path for storage/retrieval")(
    "rct-path", bpo::value<std::string>()->default_value("RCT/RunInformation"), "path to the ccdb RCT objects for the SOR/EOR timestamps")(
    "start,s", bpo::value<long>()->default_value(0), "Start timestamp of object validity. If 0 and runnumber != 0 it will be set to the run SOR")(
    "stop,S", bpo::value<long>()->default_value(0), "Stop timestamp of object validity. If 0 and runnumber != 0 it will be set to the run EOR")(
    "timestamp,T", bpo::value<long>()->default_value(-1), "Timestamp of the object to retrieve, used in alternative to the run number")(
    "runnumber,R", bpo::value<unsigned int>()->default_value(0), "Timestamp of the object to retrieve, used in alternative to the timestamp (if 0 using the timestamp)")(
    "delete-previous,delete_previous,d", bpo::value<int>()->default_value(0), "Flag to delete previous versions of converter objects in the CCDB before uploading the new one so as to avoid proliferation on CCDB")(
    "save-to-file,file,f,o", bpo::value<std::string>()->default_value(""), "Option to save parametrization to file instead of uploading to ccdb")(
    "read-from-file,i", bpo::value<std::string>()->default_value(""), "Option to get parametrization from a file")(
    "exp-name,n", bpo::value<std::string>()->default_value("BetheBloch"), "Name of the parametrization object")(
    "reso-name,n", bpo::value<std::string>()->default_value("TPCReso"), "Name of the parametrization object")(
    "mode,m", bpo::value<unsigned int>()->default_value(1), "Working mode: 0 push 1 pull and test")(
    "p0", bpo::value<float>()->default_value(0.0320981), "Parameter 0 of the TPC expected value")(
    "p1", bpo::value<float>()->default_value(19.9768), "Parameter 1 of the TPC expected value")(
    "p2", bpo::value<float>()->default_value(2.52666e-16), "Parameter 2 of the TPC expected value")(
    "p3", bpo::value<float>()->default_value(2.72123), "Parameter 3 of the TPC expected value")(
    "p4", bpo::value<float>()->default_value(6.08092), "Parameter 4 of the TPC expected value")(
    "p5", bpo::value<float>()->default_value(50.f), "Parameter 5 of the TPC expected value")(
    "p6", bpo::value<float>()->default_value(2.3), "Parameter 6 of the TPC expected value")(
    "dryrun,D", bpo::value<int>()->default_value(0), "Perform a dryrun check before uploading")(
    "verbose,v", bpo::value<int>()->default_value(0), "Verbose level 0, 1")(
    "help,h", "Produce help message.");
  try {
    bpo::store(parse_command_line(argc, argv, options), arguments);

    // help
    if (arguments.count("help")) {
      LOG(info) << options;
      return false;
    }

    bpo::notify(arguments);
  } catch (const bpo::error& e) {
    LOG(error) << e.what() << "\n";
    LOG(error) << "Error parsing command line arguments; Available options:";
    LOG(error) << options;
    return false;
  }
  return true;
}

int main(int argc, char* argv[])
{
  bpo::options_description options("Allowed options");
  if (!initOptionsAndParse(options, argc, argv)) {
    return 1;
  }

  // Fetch options
  const auto mode = arguments["mode"].as<unsigned int>();
  const auto runnumber = arguments["runnumber"].as<unsigned int>();
  auto timestamp = arguments["timestamp"].as<long>();
  const auto path = arguments["ccdb-path"].as<std::string>();
  auto start = arguments["start"].as<long>();
  auto stop = arguments["stop"].as<long>();

  // Init CCDB
  const std::string url = arguments["url"].as<std::string>();
  api.init(url);
  if (!api.isHostReachable()) {
    LOG(warning) << "CCDB host " << url << " is not reacheable, cannot go forward";
    return 1;
  }

  // Init timestamps
  setupTimestamps(timestamp, start, stop);

  BetheBloch* bb = nullptr;
  TPCReso* reso = nullptr;
  const std::string exp_name = arguments["exp-name"].as<std::string>();
  const std::string reso_name = arguments["reso-name"].as<std::string>();
  if (mode == 0) { // Push mode
    LOG(info) << "Handling TPC parametrization in create mode";
    const std::string input_file_name = arguments["read-from-file"].as<std::string>();
    if (!input_file_name.empty()) {
      TFile f(input_file_name.data(), "READ");
      if (!f.IsOpen()) {
        LOG(warning) << "Input file " << input_file_name << " is not reacheable, cannot get param from file";
      }
      f.GetObject(exp_name.c_str(), bb);
      f.GetObject(reso_name.c_str(), reso);
      f.Close();
    }
    if (!bb) {
      bb = new BetheBloch();
      const std::vector<float> bbparams = {arguments["p0"].as<float>(), arguments["p1"].as<float>(), arguments["p2"].as<float>(), arguments["p3"].as<float>(), arguments["p4"].as<float>(), arguments["p5"].as<float>(), arguments["p6"].as<float>()};
      bb->SetParameters(bbparams);
    }
    bb->Print();
    if (!reso) {
      reso = new TPCReso();
      const std::vector<float> resoparams = {0.07, 0.0};
      reso->SetParameters(resoparams);
    }
    reso->Print();
    const std::string fname = arguments["save-to-file"].as<std::string>();
    if (!fname.empty()) { // Saving it to file
      LOG(info) << "Saving parametrization to file " << fname;
      TFile f(fname.data(), "RECREATE");
      bb->Write();
      bb->GetParameters().Write();
      reso->Write();
      reso->GetParameters().Write();
      f.ls();
      f.Close();
    } else { // Saving it to CCDB
      LOG(info) << "Saving parametrization to CCDB " << path << " with validity " << start << " "
                << stop;

      std::map<std::string, std::string> metadata;
      if (runnumber != 0) {
        metadata["runnumber"] = Form("%i", runnumber);
      }
      // Storing parametrization objects
      storeOnCCDB(path + "/" + exp_name, metadata, start, stop, bb);
      storeOnCCDB(path + "/" + reso_name, metadata, start, stop, reso);
      // Storing parametrization parameters
      o2::pid::Parameters* params;
      bb->GetParameters(params);
      storeOnCCDB(path + "/Parameters/" + exp_name, metadata, start, stop, params);
      reso->GetParameters(params);
      storeOnCCDB(path + "/Parameters/" + reso_name, metadata, start, stop, params);
    }
  } else { // Pull and test mode
    LOG(info) << "Handling TPC parametrization in test mode for timestamp "
              << timestamp << " -> " << timeStampToHReadble(timestamp);
    const float x[2] = {1, 1};
    bb = retrieveFromCCDB<BetheBloch>(path + "/" + exp_name, timestamp);
    bb->Print();
    LOG(info) << "BetheBloch " << bb->operator()(x);

    reso = retrieveFromCCDB<TPCReso>(path + "/" + reso_name, timestamp);
    reso->Print();
    LOG(info) << "TPCReso " << reso->operator()(x);
  }

  return 0;
}
