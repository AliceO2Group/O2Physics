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
/// \file handleParamTPCResponse.cxx
/// \author Jeremy Wilkinson
/// \brief exec for writing and reading Response object

#include "CCDB/CcdbApi.h"
#include <boost/program_options.hpp>
#include <FairLogger.h>
#include <array>
#include "TFile.h"
#include "Common/Core/PID/PIDResponse.h"
#include "Common/Core/PID/TPCPIDResponse.h"
using namespace o2::pid::tpc;
namespace bpo = boost::program_options;

bool initOptionsAndParse(bpo::options_description& options, int argc, char* argv[], bpo::variables_map& vm)
{
  options.add_options()(
    "url,u", bpo::value<std::string>()->default_value("http://alice-ccdb.cern.ch"), "URL of the CCDB database")(
    "ccdb-path,c", bpo::value<std::string>()->default_value("Analysis/PID/TPC/Response"), "CCDB path for object")(
    "start,s", bpo::value<long>()->default_value(0), "Start timestamp for calibration validity")(
    "stop,S", bpo::value<long>()->default_value(4108971600000), "End timestamp for calibration validity")(
    "delete-previous", bpo::value<int>()->default_value(0), "delete previous entry from CCDB (1 = true)")(
    "save-to-file,file,f,o", bpo::value<std::string>()->default_value(""), "Option to save parametrization to file instead of uploading to ccdb")(
    "read-from-file,i", bpo::value<std::string>()->default_value(""), "Option to get parametrization from a file")(
    "objname,n", bpo::value<std::string>()->default_value("Response"), "Object name to be stored in file")(
    "bb0", bpo::value<float>()->default_value(0.0320981f), "Bethe-Bloch parameter 0")(
    "bb1", bpo::value<float>()->default_value(19.9768f), "Bethe-Bloch parameter 1")(
    "bb2", bpo::value<float>()->default_value(2.52666e-16f), "Bethe-Bloch parameter 2")(
    "bb3", bpo::value<float>()->default_value(2.72123f), "Bethe-Bloch parameter 3")(
    "bb4", bpo::value<float>()->default_value(6.08092f), "Bethe-Bloch parameter 4")(
    "sig0", bpo::value<float>()->default_value(0.07f), "Sigma parameter 0")(
    "sig1", bpo::value<float>()->default_value(0.f), "Sigma parameter 1")(
    "paramMIP", bpo::value<float>()->default_value(50.f), "MIP parameter value")(
    "paramChargeFactor", bpo::value<float>()->default_value(2.3f), "Charge factor value")(
    "mode", bpo::value<string>()->default_value(""), "Running mode ('read' from file, 'write' to file, 'pull' from CCDB, 'push' to CCDB)")(
    "help,h", "Print this help.");
  try {
    bpo::store(parse_command_line(argc, argv, options), vm);

    // help
    if (vm.count("help")) {
      LOG(info) << options;
      return false;
    }
    bpo::notify(vm);
  } catch (const bpo::error& e) {
    LOG(error) << e.what() << "\n";
    LOG(error) << "Error parsing command line arguments; Available options:";
    LOG(error) << options;
    return false;
  }
  return true;

} //initOptionsAndParse

int main(int argc, char* argv[])
{
  bpo::options_description options("Allowed options");
  bpo::variables_map vm;
  if (!initOptionsAndParse(options, argc, argv, vm)) {
    return 1;
  }

  Response* tpc = nullptr;

  const std::string urlCCDB = vm["url"].as<std::string>();
  const std::string pathCCDB = vm["ccdb-path"].as<std::string>();
  const long startTime = vm["start"].as<long>();
  const long endTime = vm["stop"].as<long>();
  const int optDelete = vm["delete-previous"].as<int>();

  const std::string outFilename = vm["save-to-file"].as<std::string>();
  const std::string inFilename = vm["read-from-file"].as<std::string>();
  const std::string objname = vm["objname"].as<std::string>();

  const float bb0 = vm["bb0"].as<float>();
  const float bb1 = vm["bb1"].as<float>();
  const float bb2 = vm["bb2"].as<float>();
  const float bb3 = vm["bb3"].as<float>();
  const float bb4 = vm["bb4"].as<float>();
  const float sig0 = vm["sig0"].as<float>();
  const float sig1 = vm["sig1"].as<float>();
  const float mipval = vm["paramMIP"].as<float>();
  const float chargefacval = vm["paramChargeFactor"].as<float>();
  const std::string optMode = vm["mode"].as<std::string>();
  if (optMode.empty()) {
    LOG(error) << "--mode must be specified (read, write, pull, push)";
    return 1;
  }
  // create parameter arrays from commandline options
  std::array<float, 5> BBparams = {bb0, bb1, bb2, bb3, bb4};
  std::array<float, 2> sigparams = {sig0, sig1};
  
  // initialise CCDB API
  std::map<std::string, std::string> metadata;
  std::map<std::string, std::string>* headers;
  o2::ccdb::CcdbApi api;
  if (optMode.compare("push") == 0 || optMode.compare("pull") == 0) {  // Initialise CCDB if in push/pull mode
    api.init(urlCCDB);
    if (!api.isHostReachable()) {
      LOG(warning) << "CCDB mode (push/pull) enabled but host " << urlCCDB << " is unreachable.";
      return 1;
      }
    }
  
  if (optMode.compare("read") == 0 ) { // Read existing object from local file
    if (inFilename.empty()) {
      LOG(error) << "read mode defined with no input file, please set --read-from-file";
      return 1; 
    }
    
    TFile fin(inFilename.data(), "READ");
    if (!fin.IsOpen()) {
      LOG(error) << "Input file " << inFilename << " could not be read";
      return 1;
    }

    fin.GetObject(objname.c_str(), tpc);
    if (!tpc) {
      LOG(error) << "Object with name " << objname << " could not be found in file " << inFilename;
      return 1;
    }
    LOG(info) << "Reading existing TPCPIDResponse object " << objname << " from " << inFilename << ":" ;
    tpc->PrintAll();
    return 0;
  }

  else if (optMode.compare("write") == 0 || optMode.compare("push") == 0 )  // Create new object to write to local file or push to CCDB
  {
    LOG(info) << "Creating new TPCPIDResponse object with defined parameters:";
    
    tpc = new Response();
    tpc->SetBetheBlochParams(BBparams);
    tpc->SetResolutionParams(sigparams);
    tpc->SetMIP(mipval);
    tpc->SetChargeFactor(chargefacval);
    tpc->PrintAll();

    if (optMode.compare("write") == 0) {
      if (outFilename.empty()) {
        LOG(error) << "'write' mode specified, but no output filename. Quitting";
        return 1;
      }

      LOG(info) << "Writing to output file " << outFilename;
      TFile fout(outFilename.data(), "RECREATE");
      if (!fout.IsOpen()) {
        LOG(error) << "Output file " << inFilename << " could not be written to";
        return 1;
      }
      fout.cd();
      //   tpc->Print();
      fout.WriteObject(tpc, objname.c_str());
      fout.Close();
      LOG(info) << "File successfully written";
      return 0;
    }

    if (optMode.compare("push") == 0) {
      LOG(info) << "Attempting to push object to CCDB";
      
      if (optDelete) {
        api.truncate(pathCCDB);
      }
      api.storeAsTFileAny(tpc, pathCCDB + "/" + objname, metadata, startTime, endTime);
    }

  }

  if (optMode.compare("pull") == 0) { // pull existing from CCDB; write out to file if requested
      LOG(info) << "Attempting to pull object from CCDB";

      tpc = api.retrieveFromTFileAny<Response>(pathCCDB + "/" + objname, metadata, -1, headers);
      tpc->PrintAll();

      if (!outFilename.empty()) {
        LOG(info) << "Writing pulled object to local file";
        TFile fout(outFilename.data(), "RECREATE");
        if (!fout.IsOpen()) {
          LOG(error) << "Output file " << inFilename << " could not be written to";
          return 1;
        }
        fout.cd();
        //   tpc->Print();
        fout.WriteObject(tpc, objname.c_str());
        fout.Close();
        LOG(info) << "File successfully written";
      }
      return 0;
  }

  else {
    LOG(error) << "Invalid mode specified! (must be 'read', 'write', 'pull' or 'push')";
    return 1;
  }
} //main
