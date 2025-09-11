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
/// \brief exec for writing and reading TPC PID Response object

#include "Common/Core/PID/TPCPIDResponse.h"
#include "Common/Tools/PID/handleParamBase.h"

#include <Algorithm/RangeTokenizer.h>
#include <Framework/Logger.h>

#include <TFile.h>
#include <TString.h>

#include <boost/algorithm/string/predicate.hpp>
#include <boost/program_options.hpp>

#include <array>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using namespace o2::pid::tpc;

bool initOptionsAndParse(bpo::options_description& options, int argc, char* argv[])
{
  options.add_options()(
    "ccdb-path,c", bpo::value<std::string>()->default_value("Analysis/PID/TPC"), "CCDB path for storage/retrieval")(
    "objname,n", bpo::value<std::string>()->default_value("Response"), "Object name to be stored in file")(
    "inobjname", bpo::value<std::string>()->default_value("Response"), "Object name to be read from file in 'push' mode")(
    "recopass", bpo::value<std::string>()->default_value(""), "Reconstruction pass name")(
    "period", bpo::value<std::string>()->default_value(""), "Period name")(
    "jiraticket", bpo::value<std::string>()->default_value(""), "JIRA ticket")(
    "comment", bpo::value<std::string>()->default_value(""), "Comment for metadata")(
    "bb0", bpo::value<float>()->default_value(0.03209809958934784f), "Bethe-Bloch parameter 0")(
    "bb1", bpo::value<float>()->default_value(19.9768009185791f), "Bethe-Bloch parameter 1")(
    "bb2", bpo::value<float>()->default_value(2.5266601063857674e-16f), "Bethe-Bloch parameter 2")(
    "bb3", bpo::value<float>()->default_value(2.7212300300598145f), "Bethe-Bloch parameter 3")(
    "bb4", bpo::value<float>()->default_value(6.080920219421387f), "Bethe-Bloch parameter 4")(
    "s0", bpo::value<float>()->default_value(0.07), "resolution parameter 0")(
    "s1", bpo::value<float>()->default_value(0.), "resolution parameter 1")(
    "reso-param-path", bpo::value<std::string>()->default_value(""), "Path to resolution parameter file")(
    "sigmaGlobal", bpo::value<std::string>()->default_value("5.43799e-7,0.053044,0.667584,0.0142667,0.00235175,1.22482,2.3501e-7,0.031585"), "Sigma parameters global")(
    "paramMIP", bpo::value<float>()->default_value(50.f), "MIP parameter value")(
    "paramChargeFactor", bpo::value<float>()->default_value(2.299999952316284f), "Charge factor value")(
    "paramMultNormalization", bpo::value<float>()->default_value(11000.), "Multiplicity Normalization")(
    "paramnClNormalization", bpo::value<float>()->default_value(152.), "Maximum nClusters for normalisation (159 for run 2, 152 for run 3)")(
    "useDefaultParam", bpo::value<bool>()->default_value(true), "Use default sigma parametrisation")(
    "mode", bpo::value<std::string>()->default_value(""), "Running mode ('read' from file, 'write' to file, 'pull' from CCDB, 'push' to CCDB)")(
    "help,h", "Produce help message.");
  setStandardOpt(options);
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

} // initOptionsAndParse

int main(int argc, char* argv[])
{
  bpo::options_description options("Allowed options");
  if (!initOptionsAndParse(options, argc, argv)) { //, arguments
    return 1;
  }

  Response* tpc = nullptr;

  const std::string urlCCDB = arguments["url"].as<std::string>();
  const auto pathCCDB = arguments["ccdb-path"].as<std::string>();

  const std::string outFilename = arguments["save-to-file"].as<std::string>();
  const std::string inFilename = arguments["read-from-file"].as<std::string>();
  const std::string objname = arguments["objname"].as<std::string>();
  const std::string inobjname = arguments["inobjname"].as<std::string>();
  const std::string recopass = arguments["recopass"].as<std::string>();
  const std::string periodname = arguments["period"].as<std::string>();
  const std::string jiraticket = arguments["jiraticket"].as<std::string>();
  const std::string comment = arguments["comment"].as<std::string>();

  const float bb0 = arguments["bb0"].as<float>();
  const float bb1 = arguments["bb1"].as<float>();
  const float bb2 = arguments["bb2"].as<float>();
  const float bb3 = arguments["bb3"].as<float>();
  const float bb4 = arguments["bb4"].as<float>();
  const float s0 = arguments["s0"].as<float>();
  const float s1 = arguments["s1"].as<float>();
  const std::string pathResoParam = arguments["reso-param-path"].as<std::string>();
  const std::string sigmaGlobal = arguments["sigmaGlobal"].as<std::string>();
  const float mipval = arguments["paramMIP"].as<float>();
  const float chargefacval = arguments["paramChargeFactor"].as<float>();
  const float multNormval = arguments["paramMultNormalization"].as<float>();
  const float nClNormval = arguments["paramnClNormalization"].as<float>();
  const bool useDefaultParam = arguments["useDefaultParam"].as<bool>();
  const std::string optMode = arguments["mode"].as<std::string>();
  if (optMode.empty()) {
    LOG(error) << "--mode must be specified (read, write, pull, push)";
    return 1;
  }
  // create parameter arrays from commandline options
  std::array<float, 5> BBparams = {bb0, bb1, bb2, bb3, bb4};
  std::array<float, 2> sparams = {s0, s1};
  std::vector<double> sigparamsGlobal;

  if (pathResoParam != "") { // Read resolution parameters from file
    std::ifstream infile(pathResoParam);
    std::string paramstring;
    std::vector<double> sigmaParams;
    std::getline(infile, paramstring);
    std::istringstream ss(paramstring);
    sigmaParams.push_back({});
    double param = 0.;
    while (ss >> param) {
      sigmaParams.push_back(param);
    }
    sigparamsGlobal = sigmaParams;
  } else {
    sigparamsGlobal = o2::RangeTokenizer::tokenize<double>(sigmaGlobal);
  }

  if (optMode.compare("push") == 0 || optMode.compare("pull") == 0) { // Initialise CCDB if in push/pull mode
    initCCDBApi();
    setupTimestamps(ccdbTimestamp, validityStart, validityStop);
  }

  if (optMode.compare("read") == 0) { // Read existing object from local file
    if (inFilename.empty()) {
      LOG(error) << "Read mode defined with no input file, please set --read-from-file";
      return 1;
    }

    TFile fin(inFilename.data(), "READ");
    if (!fin.IsOpen()) {
      LOG(error) << "Input file " << inFilename << " could not be read";
      return 1;
    }

    tpc = fin.Get<Response>(objname.c_str());
    if (!tpc) {
      LOG(error) << "Object with name " << objname << " could not be found in file " << inFilename;
      return 1;
    }
    LOG(info) << "Reading existing TPCPIDResponse object " << objname << " from " << inFilename << ":";
    tpc->PrintAll();
    return 0;

  } else if (optMode.compare("write") == 0 || optMode.compare("push") == 0) { // Create new object to write to local file or push to CCDB

    if (!inFilename.empty()) { // Read from existing file to push to CCDB
      LOG(info) << "Reading from existing file to write to CCDB:";
      TFile fin(inFilename.data(), "READ");
      if (!fin.IsOpen()) {
        LOG(error) << "Input file " << inFilename << " could not be read";
        return 1;
      }

      tpc = fin.Get<Response>(inobjname.c_str());
      if (!tpc) {
        LOG(error) << "Object with name " << objname << " could not be found in file " << inFilename;
        return 1;
      }

      tpc->PrintAll();
    } else { // Create new object if file not specified
      LOG(info) << "Creating new TPCPIDResponse object with defined parameters:";

      tpc = new Response();
      tpc->SetBetheBlochParams(BBparams);
      tpc->SetResolutionParamsDefault(sparams);
      tpc->SetResolutionParams(sigparamsGlobal);
      tpc->SetMIP(mipval);
      tpc->SetChargeFactor(chargefacval);
      tpc->SetMultiplicityNormalization(multNormval);
      tpc->SetNClNormalization(nClNormval);
      tpc->SetUseDefaultResolutionParam(useDefaultParam);
      tpc->PrintAll();
    }

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
      fout.WriteObject(tpc, objname.c_str());
      fout.Close();
      LOG(info) << "File successfully written";
      return 0;
    }

    if (optMode.compare("push") == 0) {
      LOG(info) << "Attempting to push object to CCDB";
      std::map<std::string, std::string> metadata;

      // Sanity check: Request confirmation if any of recoPass, period name, comment, jira not specified
      int missingPass = 0, missingPeriod = 0, missingComment = 0, missingJira = 0;
      if (recopass.empty()) {
        missingPass = 1;
      } else {
        metadata["RecoPassName"] = recopass;
      }

      if (periodname.empty()) {
        missingPeriod = 1;
      } else {
        metadata["LPMProductionTag"] = periodname;
      }

      if (comment.empty()) {
        missingComment = 1;
      } else {
        metadata["Comment"] = comment;
      }
      if (jiraticket.empty()) {
        missingJira = 1;
      } else {
        metadata["JIRA"] = jiraticket;
      }

      if (missingPass || missingPeriod || missingComment || missingJira) {
        LOG(info) << "WARNING: Attempting to push an object with missing metadata elements:";
        if (missingPass)
          LOG(info) << "\t- Pass name";
        if (missingPeriod)
          LOG(info) << "\t- Period name";
        if (missingComment)
          LOG(info) << "\t- Comment";
        if (missingJira)
          LOG(info) << "\t- JIRA ticket";

        // Request interactive confirmation to upload
        LOG(info) << "Continue with object upload anyway? (Y/n)";
        std::string confirm;
        std::cin >> confirm;
        if (boost::iequals(confirm.substr(0, 1), "y")) {
          LOG(info) << "Continuing with object upload";
        } else {
          LOG(info) << "Aborting upload";
          return 1;
        }
      }

      // Fill metadata map for start/end run number
      if (minRunNumber != 0) {
        metadata["min-runnumber"] = Form("%i", minRunNumber);
        metadata["max-runnumber"] = Form("%i", maxRunNumber);
      }

      storeOnCCDB(pathCCDB + "/" + objname, metadata, validityStart, validityStop, tpc);
    }
  } else if (optMode.compare("pull") == 0) { // pull existing from CCDB; write out to file if requested
    LOG(info) << "Attempting to pull object from CCDB (" << urlCCDB << "): " << pathCCDB << "/" << objname << (recopass.empty() ? "" : "/recoPass=" + recopass);
    std::map<std::string, std::string> metadata;
    if (!recopass.empty()) {
      metadata["RecoPassName"] = recopass;
    }

    tpc = retrieveFromCCDB<Response>(pathCCDB + "/" + objname, ccdbTimestamp, metadata);
    if (!tpc) {
      LOG(error) << "No CCDB object found with requested parameters";
      return 1;
    }
    tpc->PrintAll();

    if (!outFilename.empty()) {
      LOG(info) << "Writing pulled object to local file";
      TFile fout(outFilename.data(), "RECREATE");
      if (!fout.IsOpen()) {
        LOG(error) << "Output file " << inFilename << " could not be written to";
        return 1;
      }
      fout.cd();
      fout.WriteObject(tpc, objname.c_str());
      fout.Close();
      LOG(info) << "File successfully written";
    }
    return 0;
  } else {
    LOG(error) << "Invalid mode specified! (must be 'read', 'write', 'pull' or 'push')";
    return 1;
  }
} // main
