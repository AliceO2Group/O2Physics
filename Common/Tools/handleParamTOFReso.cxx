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
/// \file   handleParamTOFReso.cxx
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \since  2020-06-22
/// \brief  A simple tool to produce resolution parametrization objects for the TOF PID Response
///

#include "TGraph.h"
#include "TSystem.h"
#include "TCanvas.h"

#include "PID/ParamBase.h"
#include "PID/DetectorResponse.h"
#include "PID/PIDTOF.h"
#include "PID/TOFReso.h"
#include "handleParamBase.h"

#include <chrono>
using namespace std::chrono;
using namespace o2::pid::tof;
using namespace o2::pid;

// Utility class
struct DebugTrack { // Track that mimics the O2 data structure
  float mp = 0.1f;
  float p() const { return mp; }
  bool isEvTimeDefined() const { return true; }
  bool hasTOF() const { return true; }
  float tofEvTime() const { return 0.f; }
  float tofEvTimeErr() const { return 20.f; }
  float length() const { return 400.f; }
  float tofSignal() const { return length() * 1.01 / kCSPEED; }
  int trackType() const { return 0; };
  float tofExpMom() const
  {
    constexpr float mass = o2::constants::physics::MassPionCharged; // default pid = pion
    const float expBeta = (length() / (tofSignal() * kCSPEED));
    return mass * expBeta / std::sqrt(1.f - expBeta * expBeta);
  }

} debugTrack;

bool initOptionsAndParse(bpo::options_description& options, int argc, char* argv[])
{
  options.add_options()(
    "ccdb-path,c", bpo::value<std::string>()->default_value("Analysis/PID/TOF"), "CCDB path for storage/retrieval")(
    "reso-name,n", bpo::value<std::string>()->default_value("TOFResoParams"), "Name of the parametrization object")(
    "mode,m", bpo::value<unsigned int>()->default_value(1), "Working mode: 0 push, 1 pull and test, 2 create and performance")(
    "p0", bpo::value<float>()->default_value(0.008f), "Parameter 0 of the TOF resolution")(
    "p1", bpo::value<float>()->default_value(0.008f), "Parameter 1 of the TOF resolution")(
    "p2", bpo::value<float>()->default_value(0.002f), "Parameter 2 of the TOF resolution")(
    "p3", bpo::value<float>()->default_value(40.0f), "Parameter 3 of the TOF resolution")(
    "p4", bpo::value<float>()->default_value(60.0f), "Parameter 4 of the TOF resolution: average TOF resolution");
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
}

int main(int argc, char* argv[])
{
  bpo::options_description options("Allowed options");
  if (!initOptionsAndParse(options, argc, argv)) {
    return 1;
  }

  // Fetch options
  const auto mode = arguments["mode"].as<unsigned int>();
  const auto ccdbPath = arguments["ccdb-path"].as<std::string>() + "/" + arguments["reso-name"].as<std::string>();

  // Init CCDB
  initCCDBApi();

  // Init timestamps
  setupTimestamps(ccdbTimestamp, validityStart, validityStop);

  TOFReso* resoOld = nullptr;
  TOFResoParams* reso = nullptr;
  if (mode == 0) { // Push mode
    LOG(info) << "Handling TOF parametrization in create mode";
    const std::string input_file_name = arguments["read-from-file"].as<std::string>();
    reso = new TOFResoParams();
    if (!input_file_name.empty()) { // Load parameters from input file
      LOG(info) << "Loading parameters from file";
      reso->LoadParamFromFile(input_file_name.c_str(), arguments["reso-name"].as<std::string>().c_str());
    } else { // Create new object
      LOG(info) << "Loading parameters from command line";
      reso->SetParameters(std::array<float, 5>{arguments["p0"].as<float>(),
                                               arguments["p1"].as<float>(),
                                               arguments["p2"].as<float>(),
                                               arguments["p3"].as<float>(),
                                               arguments["p4"].as<float>()});
    }
    reso->Print();
    const std::string fname = arguments["save-to-file"].as<std::string>();
    if (!fname.empty()) { // Saving it to file
      LOG(info) << "Saving parametrization to file " << fname;
      TFile f(fname.data(), "RECREATE");
      reso->Write();
      f.ls();
      f.Close();
    } else { // Saving it to CCDB
      LOG(info) << "Saving parametrization to CCDB " << ccdbPath << " with validity " << validityStart << " "
                << validityStop;

      std::map<std::string, std::string> metadata;
      if (minRunNumber != 0) {
        metadata["min-runnumber"] = Form("%i", minRunNumber);
        metadata["max-runnumber"] = Form("%i", maxRunNumber);
      }
      reso->AddToMetadata(metadata);
      // Storing parametrization parameters
      storeOnCCDB(ccdbPath, metadata, validityStart, validityStop, reso);
    }
  } else if (mode == 1) { // Pull and test mode
    LOG(info) << "Handling TOF parametrization in test mode for timestamp "
              << ccdbTimestamp << " -> " << timeStampToHReadble(ccdbTimestamp);
    reso = retrieveFromCCDB<TOFResoParams>(ccdbPath, ccdbTimestamp);
    reso->Print();
    using RespImp = ExpTimes<DebugTrack, 2>;
    LOG(info) << "TOF expected resolution at p=" << debugTrack.p() << " GeV/c and mass " << RespImp::mMassZ << ":" << RespImp::GetExpectedSigma(*reso, debugTrack);
  } else { // Create and test + performance
    LOG(info) << "Creating TOF parametrization and testing";
    reso = new TOFResoParams();
    reso->Print();
    DetectorResponse response;
    if (!resoOld) {
      resoOld = new TOFReso();
      const std::vector<float> resoparams = {arguments["p0"].as<float>(), arguments["p1"].as<float>(), arguments["p2"].as<float>(), arguments["p3"].as<float>(), arguments["p4"].as<float>()};
      resoOld->SetParameters(resoparams);
    }
    response.LoadParam(DetectorResponse::kSigma, resoOld);
    // Draw it
    using RespImp = ExpTimes<DebugTrack, 2>;
    //
    std::map<std::string, TGraph*> graphs;
    graphs["ExpSigma"] = new TGraph();
    graphs["NSigma"] = new TGraph();
    graphs["durationExpSigma"] = new TGraph();
    graphs["durationNSigma"] = new TGraph();
    //
    graphs["ExpSigmaOld"] = new TGraph();
    graphs["NSigmaOld"] = new TGraph();
    graphs["durationExpSigmaOld"] = new TGraph();
    graphs["durationNSigmaOld"] = new TGraph();
    const int nsamp = 1000;
    for (int i = 0; i < nsamp; i++) {
      debugTrack.mp += 0.01f;
      //
      auto start = high_resolution_clock::now();
      RespImp::GetExpectedSigma(*reso, debugTrack);
      auto stop = high_resolution_clock::now();
      auto duration = duration_cast<nanoseconds>(stop - start).count();
      graphs["ExpSigma"]->SetPoint(i, debugTrack.p(), RespImp::GetExpectedSigma(*reso, debugTrack));
      graphs["durationExpSigma"]->SetPoint(i + 1, i, duration);
      //
      start = high_resolution_clock::now();
      RespImp::GetSeparation(*reso, debugTrack);
      stop = high_resolution_clock::now();
      duration = duration_cast<nanoseconds>(stop - start).count();
      graphs["NSigma"]->SetPoint(i, debugTrack.p(), RespImp::GetSeparation(*reso, debugTrack));
      graphs["durationNSigma"]->SetPoint(i + 1, i, duration);
      //
      start = high_resolution_clock::now();
      RespImp::GetExpectedSigma(response, debugTrack);
      stop = high_resolution_clock::now();
      duration = duration_cast<nanoseconds>(stop - start).count();
      graphs["ExpSigmaOld"]->SetPoint(i, debugTrack.p(), RespImp::GetExpectedSigma(*reso, debugTrack));
      graphs["durationExpSigmaOld"]->SetPoint(i + 1, i, duration);
      //
      start = high_resolution_clock::now();
      RespImp::GetSeparation(response, debugTrack);
      stop = high_resolution_clock::now();
      duration = duration_cast<nanoseconds>(stop - start).count();
      graphs["NSigmaOld"]->SetPoint(i, debugTrack.p(), RespImp::GetSeparation(*reso, debugTrack));
      graphs["durationNSigmaOld"]->SetPoint(i + 1, i, duration);
    }
    TFile fdebug("/tmp/tofParamDebug.root", "UPDATE");
    TString dn = Form("%i", fdebug.GetListOfKeys()->GetEntries());
    LOG(info) << "Saving performance graphs to " << fdebug.GetName() << " iteration " << dn;
    fdebug.mkdir(dn);
    fdebug.cd(dn);
    for (const auto& i : graphs) {
      i.second->SetName(i.first.c_str());
      i.second->SetTitle(i.first.c_str());
      i.second->Write();
    }
    fdebug.Close();
    LOG(info) << "TOF expected resolution at p=" << debugTrack.p() << " GeV/c and mass " << RespImp::mMassZ << ": " << RespImp::GetExpectedSigma(*reso, debugTrack) << " ps";
  }

  return 0;
}
