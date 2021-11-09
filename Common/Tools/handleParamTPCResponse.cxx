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
#include <array>
#include "TFile.h"
#include "Common/Core/PID/PIDResponse.h"
#include "Common/Core/PID/TPCPIDResponse.h"
// std::array<float,5>( {0.0320981, 19.9768, 2.52666e-16, 2.72123, 6.08092})
using namespace o2::pid::tpc;
namespace bpo = boost::program_options;

bool initOptionsAndParse(bpo::options_description& options, int argc, char* argv[], bpo::variables_map& vm)
{
  options.add_options()(
    "save-to-file,file,f,o", bpo::value<std::string>()->default_value(""), "Option to save parametrization to file instead of uploading to ccdb")(
    "read-from-file,i", bpo::value<std::string>()->default_value(""), "Option to get parametrization from a file")(
    "objname,n", bpo::value<std::string>()->default_value("TPCRespCustom"), "Object name to be stored in file")(
    "bb0", bpo::value<float>()->default_value(0.0320981f), "Bethe-Bloch parameter 0")(
    "bb1", bpo::value<float>()->default_value(19.9768f), "Bethe-Bloch parameter 1")(
    "bb2", bpo::value<float>()->default_value(2.52666e-16f), "Bethe-Bloch parameter 2")(
    "bb3", bpo::value<float>()->default_value(2.72123f), "Bethe-Bloch parameter 3")
    "bb4", bpo::value<float>()->default_value(6.08092f), "Bethe-Bloch parameter 4")(
    "sig0", bpo::value<float>()->default_value(0.07f), "Sigma parameter 0")(
    "sig1", bpo::value<float>()->default_value(0.f), "Sigma parameter 1")(
    "paramMIP", bpo::value<float>()->default_value(50.f), "MIP parameter value")(
    "paramChargeFactor", bpo::value<float>()->default_value(2.3f), "Charge factor value")(
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

  TPCPIDResponse* tpc = nullptr;
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

  std::array<float, 5> BBparams = {bb0, bb1, bb2, bb3, bb4};
  std::array<float, 2> sigparams = {sig0, sig1};

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
    tpc->SetBetheBlochParams(BBparams);
    tpc->SetResolutionParams(sigparams);
    tpc->SetMIP(mipval);
    tpc->SetChargeFactor(chargefacval);
    fout.cd();
    //   tpc->Print();
    tpc->Write();
    fout.Close();
  }

  if (!inFilename.empty()) { // Read and execute object from file
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
    tpc->PrintAll();
    //   tpc->Print();
  }

  return 0;
} //main
