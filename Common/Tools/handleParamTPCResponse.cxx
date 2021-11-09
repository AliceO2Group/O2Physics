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
/// \brief exec for writing and reading TPCPIDResponse object

#include "CCDB/CcdbApi.h"
#include <boost/program_options.hpp>
#include <FairLogger.h>
#include "TFile.h"
#include "Common/Core/PID/PIDResponse.h"
#include "Common/Core/PID/TPCPIDResponse.h"

using namespace o2::pid::tpc;
namespace bpo = boost::program_options;

bool initOptionsAndParse(bpo::options_description& options, int argc, char* argv[], bpo::variables_map& vm)
{
  options.add_options()(
    "save-to-file,file,f,o", bpo::value<std::string>()->default_value(""), "Option to save parametrization to file instead of uploading to ccdb")(
    "read-from-file,i", bpo::value<std::string>()->default_value(""), "Option to get parametrization from a file")(
    "objname,n", bpo::value<std::string>()->default_value("TPCRespCustom"), "Object name to be stored in file")(
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

} // initOptionsAndParse

int main(int argc, char* argv[])
{
  bpo::options_description options("Allowed options");
  bpo::variables_map vm;
  if (!initOptionsAndParse(options, argc, argv, vm)) {
    return 1;
  }

  TPCPIDResponse* tpc = nullptr;
  const std::string outFilename = vm["save-to-file"].as<std::string>();
  const std::string inFilename = vm["read-from-file"].as<std::string>();
  const std::string objname = vm["objname"].as<std::string>();
  if (!outFilename.empty() && !inFilename.empty()) {
    LOG(error) << "Cannot both read and write at the same time!";
    return 1;
  }
  if (outFilename.empty() && inFilename.empty()) {
    LOG(error) << "neither input nor output defined!";
    return 1;
  }

  if (!outFilename.empty()) { // Write new object to file
    TFile fout(outFilename.data(), "RECREATE");
    tpc = new TPCPIDResponse();
    fout.cd();
    //   tpc->Print();
    tpc->Write();
    fout.Close();
  }

  if (!inFilename.empty()) { // Read and execute object from file
    TFile fin(inFilename.data(), "READ");
    if (!fin.IsOpen()) {
      LOG(warning) << "Input file " << inFilename << " could not be read";
    }

    fin.GetObject(objname.c_str(), tpc);
    //   tpc->Print();
  }

} // main