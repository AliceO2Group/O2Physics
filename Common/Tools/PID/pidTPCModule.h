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

#ifndef COMMON_TOOLS_PID_PIDTPCMODULE_H_
#define COMMON_TOOLS_PID_PIDTPCMODULE_H_

#include "Common/CCDB/ctpRateFetcher.h"
#include "Common/Core/CollisionTypeHelper.h"
#include "Common/Core/PID/TPCPIDResponse.h"
#include "Common/Core/TableHelper.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/TableProducer/PID/pidTPCBase.h" // IWYU pragma: keep
#include "Tools/ML/model.h"

#include <DataFormatsParameters/GRPLHCIFData.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/Configurable.h>
#include <Framework/DeviceSpec.h>
#include <Framework/RunningWorkflowInfo.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/PID.h>

#include <TFile.h>
#include <TMatrixD.h> // IWYU pragma: keep (do not replace with TMatrixDfwd.h)
#include <TMatrixDfwd.h>
#include <TRandom.h>

#include <chrono>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <map>
#include <memory>
#include <ratio>
#include <string>
#include <vector>

#include <math.h>

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
  o2::framework::Configurable<int> skipTPCOnly{"skipTPCOnly", -1, "Flag to skip TPC only tracks (faster but affects the analyses that use TPC only tracks). 0: do not skip, 1: skip, -1: check if needed by specific tasks"};
  o2::framework::Configurable<std::vector<std::string>> devicesRequiringTPCOnlyPID{"devicesRequiringTPCOnlyPID", std::vector<std::string>{"photon-conversion-builder"}, "List of device names of tasks requiring TPC-only tracks to have TPC PID calculated"};
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
  o2::framework::Configurable<std::string> cfgPathGrpLhcIf{"ccdb-path-grplhcif", "GLO/Config/GRPLHCIF", "Path on the CCDB for the GRPLHCIF object"};
};

// helper getter - FIXME should be separate
int getPIDIndex(const int pdgCode) // Get O2 PID index corresponding to MC PDG code
{
  switch (std::abs(pdgCode)) {
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

typedef struct Str_dEdx_correction {
  TMatrixD fMatrix;
  bool warning = true;

  // void init(std::vector<double>& params)
  void init()
  {
    double elements[32] = {0.99091, -0.015053, 0.0018912, -0.012305,
                           0.081387, 0.003205, -0.0087404, -0.0028608,
                           0.013066, 0.017012, -0.0018469, -0.0052177,
                           -0.0035655, 0.0017846, 0.0019127, -0.00012964,
                           0.0049428, 0.0055592, -0.0010618, -0.0016134,
                           -0.0059098, 0.0013335, 0.00052133, 3.1119e-05,
                           -0.004882, 0.00077317, -0.0013827, 0.003249,
                           -0.00063689, 0.0016218, -0.00045215, -1.5815e-05};
    fMatrix.ResizeTo(4, 8);
    fMatrix.SetMatrixArray(elements);
  }

  float fReal_fTPCSignalN(std::vector<float> vec1, std::vector<float> vec2)
  {
    float result = 0.f;
    // push 1.
    vec1.insert(vec1.begin(), 1.0);
    vec2.insert(vec2.begin(), 1.0);
    for (int i = 0; i < fMatrix.GetNrows(); i++) {
      for (int j = 0; j < fMatrix.GetNcols(); j++) {
        double param = fMatrix(i, j);
        double value1 = i >= static_cast<int>(vec1.size()) ? 0 : vec1[i];
        double value2 = j >= static_cast<int>(vec2.size()) ? 0 : vec2[j];
        result += param * value1 * value2;
      }
    }
    return result;
  }
} Str_dEdx_correction;

class pidTPCModule
{
 public:
  o2::aod::pid::pidTPCConfigurables pidTPCopts;

  // TPC PID Response
  o2::pid::tpc::Response* response{nullptr};

  // Network correction for TPC PID response
  ml::OnnxModel network;
  std::map<std::string, std::string> metadata;
  std::map<std::string, std::string> headers;
  std::vector<int> speciesNetworkFlags = std::vector<int>(9);
  std::string networkVersion;

  // To get automatically the proper Hadronic Rate
  std::string irSource = "";
  o2::common::core::CollisionSystemType::collType collsys = o2::common::core::CollisionSystemType::kCollSysUndef;

  // Parametrization configuration
  bool useCCDBParam = false;

  // for dEdx correction
  ctpRateFetcher mRateFetcher;
  Str_dEdx_correction str_dedx_correction;

  //__________________________________________________
  template <typename TCCDB, typename TContext, typename TpidTPCOpts, typename TMetadataInfo>
  void init(TCCDB& ccdb, TContext& context, TpidTPCOpts const& external_pidtpcopts, TMetadataInfo const& metadataInfo)
  {
    // read in configurations from the task where it's used
    pidTPCopts = external_pidtpcopts;

    if (pidTPCopts.useCorrecteddEdx.value) {
      LOGF(warning, "***************************************************");
      LOGF(warning, " WARNING: YOU HAVE SWITCHED ON 'corrected dEdx!");
      LOGF(warning, " This mode is still in development and it is meant");
      LOGF(warning, " ONLY FOR EXPERTS at this time. Please switch ");
      LOGF(warning, " this option off UNLESS you are absolutely SURE");
      LOGF(warning, " of what you're doing! You've been warned!");
      LOGF(warning, "***************************************************");
    }

    if (pidTPCopts.skipTPCOnly.value == -1) {
      LOGF(info, "***************************************************");
      LOGF(info, " the skipTPConly flag has a value of -1! ");
      LOGF(info, " ---> autodetecting TPC-only track necessity now ");
      LOGF(info, "***************************************************");
      // print list of devices that are being checked for
      for (std::size_t devIdx{0}; devIdx < pidTPCopts.devicesRequiringTPCOnlyPID->size(); devIdx++) {
        LOGF(info, "Will search for #%i device requiring TPC PID for TPC only: %s", devIdx, pidTPCopts.devicesRequiringTPCOnlyPID->at(devIdx));
      }
      LOGF(info, "***************************************************");

      // assume that TPC tracks are not needed, but check if tasks
      // requiring them are present in the chain
      pidTPCopts.skipTPCOnly.value = 1;

      // loop over devices in this execution
      auto const& workflows = context.services().template get<o2::framework::RunningWorkflowInfo const>();
      for (o2::framework::DeviceSpec const& device : workflows.devices) {
        // Look for propagation service
        if (device.name.compare("propagation-service") == 0) {
          LOGF(info, " ---> propagation service detected, checking if photons enabled...");
          for (auto const& option : device.options) {
            // check for photon generation enabled or not
            if (option.name.compare("v0BuilderOpts.generatePhotonCandidates") == 0) {
              if (option.defaultValue.get<bool>()) {
                LOGF(info, " ---> propagation service: photons enabled, will calculate TPC PID for TPC only tracks.");
                pidTPCopts.skipTPCOnly.value = 0;
              } else {
                LOGF(info, " ---> propagation service: photons disabled, TPC PID not required for TPC-only tracks");
              }
            }
          }
        }

        // Check 2: specific tasks that require TPC PID based on configurable
        for (std::size_t devIdx{0}; devIdx < pidTPCopts.devicesRequiringTPCOnlyPID->size(); devIdx++) {
          if (device.name.compare(pidTPCopts.devicesRequiringTPCOnlyPID->at(devIdx)) == 0) {
            LOGF(info, " ---> %s detected! ", pidTPCopts.devicesRequiringTPCOnlyPID->at(devIdx));
            LOGF(info, " ---> enabling TPC only track TPC PID calculations now.");
            pidTPCopts.skipTPCOnly.value = 0;
          }
        }
      }

      if (pidTPCopts.skipTPCOnly.value == 1) {
        LOGF(info, "***************************************************");
        LOGF(info, "No need for TPC only information detected. Will not generate Nsigma for TPC only tracks");
        LOGF(info, "If this is unexpected behaviour and a necessity was not identified, please add the");
        LOGF(info, "corresponding task to the configurable 'devicesRequiringTPCOnlyPID' of this task");
        LOGF(info, "To do that, please get in touch with core service wagon maintainers and ask:");
        LOGF(info, "It is always best to use core service wagons instead of private copies");
      }
      LOGF(info, "***************************************************");
    }

    // initialize PID response
    response = new o2::pid::tpc::Response();

    o2::common::core::enableFlagIfTableRequired(context, "DEdxsCorrected", pidTPCopts.savedEdxsCorrected);

    // Checking the tables are requested in the workflow and enabling them
    auto enableFlag = [&](const std::string particle, o2::framework::Configurable<int>& flag) {
      o2::common::core::enableFlagIfTableRequired(context, "pidTPC" + particle, flag);
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
      o2::common::core::enableFlagIfTableRequired(context, "mcTPCTuneOnData", pidTPCopts.enableTuneOnDataTable);
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
        response = ccdb->template getSpecific<o2::pid::tpc::Response>(path, time, metadata, &headers);
        if (!response) {
          LOGF(warning, "Unable to find TPC parametrisation for specified pass name - falling back to latest object");
          response = ccdb->template getForTimeStamp<o2::pid::tpc::Response>(path, time, &headers);
          if (!response) {
            LOGF(fatal, "Unable to find any TPC object corresponding to timestamp {}!", time);
          }
        }
        networkVersion = headers["NN-Version"];
        LOG(info) << "Successfully retrieved TPC PID object from CCDB for timestamp " << time << ", period " << headers["LPMProductionTag"] << ", recoPass " << headers["RecoPassName"];
        metadata["RecoPassName"] = headers["RecoPassName"]; // Force pass number for NN request to match retrieved BB
        o2::parameters::GRPLHCIFData* grpo = ccdb->template getForTimeStamp<o2::parameters::GRPLHCIFData>(pidTPCopts.cfgPathGrpLhcIf.value, time);
        if (grpo) {
          LOG(info) << " collision type::" << CollisionSystemType::getCollisionTypeFromGrp(grpo);
          collsys = CollisionSystemType::getCollisionTypeFromGrp(grpo);
          if (collsys == CollisionSystemType::kCollSyspp) {
            irSource = std::string("T0VTX");
          } else {
            irSource = std::string("ZNC hadronic");
          }
        } else {
          LOGF(info, "No grpo object found. irSource will remain undefined.");
        }
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
          bool retrieveSuccess = ccdb->getCCDBAccessor().retrieveBlob(pidTPCopts.networkPathCCDB.value, ".", metadata, pidTPCopts.ccdbTimestamp.value, false, pidTPCopts.networkPathLocally.value, "", "", &headers);
          networkVersion = headers["NN-Version"];
          if (retrieveSuccess) {
            network.initModel(pidTPCopts.networkPathLocally.value, pidTPCopts.enableNetworkOptimizations.value, pidTPCopts.networkSetNumThreads.value, strtoul(headers["Valid-From"].c_str(), NULL, 0), strtoul(headers["Valid-Until"].c_str(), NULL, 0));
            std::vector<float> dummyInput(network.getNumInputNodes(), 1.);
            network.evalModel(dummyInput); /// Init the model evaluations
            setupColumnInputNetwork();
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
          setupColumnInputNetwork();
        }
      } else {
        return;
      }
    }

    if (pidTPCopts.useCorrecteddEdx.value && networkVersion != "5") {
      LOGF(fatal, "Using corrected dE/dx with a network version other than 5 will not work. Crashing now.");
    }
  } // end init

  //__________________________________________________
  void setupColumnInputNetwork()
  {
    using PI = o2::ml::OnnxModel::PreprocInput;
    using PF = o2::ml::OnnxModel::PreprocFeature;
    const int nFeat = network.getNumInputNodes(); // # network features (6..9), original model

    // Raw graph inputs (this order defines the tensor feeding order in
    // createNetworkPrediction). All per-track columns are wrapped zero-copy from
    // the Arrow buffers; nclNorm/hrDivisor/hrFallback are per-DF runtime scalars.
    std::vector<PI> in;
    in.push_back({"tpcInnerParam", PI::Type::TrackFloat});
    in.push_back({"tgl", PI::Type::TrackFloat});
    in.push_back({"signed1Pt", PI::Type::TrackFloat});
    in.push_back({"mass", PI::Type::ScalarFloat});
    in.push_back({"collisionId", PI::Type::TrackInt32});
    in.push_back({"multArray", PI::Type::CollisionFloat});
    in.push_back({"nclNorm", PI::Type::ScalarFloat});
    in.push_back({"nclsFindable", PI::Type::TrackUint8});
    in.push_back({"nclsFMF", PI::Type::TrackInt8});
    if (nFeat >= 7) {
      in.push_back({"occArray", PI::Type::CollisionFloat});
    }
    if (nFeat >= 8) {
      in.push_back({"hrArray", PI::Type::CollisionFloat});
      in.push_back({"hrDivisor", PI::Type::ScalarFloat});
      in.push_back({"hrFallback", PI::Type::ScalarFloat});
    }
    if (nFeat >= 9) {
      in.push_back({"phi", PI::Type::TrackFloat});
    }
    in.push_back({"validMask", PI::Type::TrackBool});

    // Per-feature preprocessing (exactly nFeat entries, in the training order).
    std::vector<PF> feat;
    {
      PF f;
      f.op = PF::Op::Passthrough;
      f.a = "tpcInnerParam";
      feat.push_back(f);
    }
    {
      PF f;
      f.op = PF::Op::Passthrough;
      f.a = "tgl";
      feat.push_back(f);
    }
    {
      PF f;
      f.op = PF::Op::Passthrough;
      f.a = "signed1Pt";
      feat.push_back(f);
    }
    {
      PF f;
      f.op = PF::Op::BroadcastScalar;
      f.a = "mass";
      f.shapeRef = "collisionId";
      feat.push_back(f);
    }
    {
      PF f;
      f.op = PF::Op::GatherNormWhere;
      f.a = "multArray";
      f.b = "collisionId";
      f.c = {11000.f, 1.f, 0.f};
      feat.push_back(f);
    }
    {
      PF f;
      f.op = PF::Op::NClsSqrtRecip;
      f.a = "nclsFindable";
      f.b = "nclsFMF";
      f.scaleInput = "nclNorm";
      feat.push_back(f);
    }
    if (nFeat >= 7) {
      PF f;
      f.op = PF::Op::GatherNormWhere;
      f.a = "occArray";
      f.b = "collisionId";
      f.c = {60000.f, 1.f, 0.f};
      feat.push_back(f);
    }
    if (nFeat >= 8) {
      PF f;
      f.op = PF::Op::GatherNormWhere;
      f.a = "hrArray";
      f.b = "collisionId";
      f.scaleInput = "hrDivisor";
      f.fallbackInput = "hrFallback";
      feat.push_back(f);
    }
    if (nFeat >= 9) {
      PF f;
      f.op = PF::Op::Mod2;
      f.a = "phi";
      f.c = {2.f * static_cast<float>(M_PI), 2.f * static_cast<float>(M_PI), static_cast<float>(M_PI) / 9.0f};
      feat.push_back(f);
    }

    network.setupColumnInputs(in, feat, "validMask");
  }

  //__________________________________________________
  template <typename TCCDB, typename M, typename T, typename B>
  std::vector<float> createNetworkPrediction(TCCDB& ccdb, soa::Join<aod::Collisions, aod::EvSels> const& collisions, M const& mults, T const& tracks, B const& bcs, const size_t size)
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
        response = ccdb->template getSpecific<o2::pid::tpc::Response>(pidTPCopts.ccdbPath.value, bc.timestamp(), metadata, &headers);
        networkVersion = headers["NN-Version"];
        if (!response) {
          LOGP(warning, "!! Could not find a valid TPC response object for specific pass name {}! Falling back to latest uploaded object.", metadata["RecoPassName"]);
          response = ccdb->template getForTimeStamp<o2::pid::tpc::Response>(pidTPCopts.ccdbPath.value, bc.timestamp(), &headers);
          if (!response) {
            LOGP(fatal, "Could not find ANY TPC response object for the timestamp {}!", bc.timestamp());
          }
        }
        LOG(info) << "Successfully retrieved TPC PID object from CCDB for timestamp " << bc.timestamp() << ", period " << headers["LPMProductionTag"] << ", recoPass " << headers["RecoPassName"];
        metadata["RecoPassName"] = headers["RecoPassName"]; // Force pass number for NN request to match retrieved BB
        o2::parameters::GRPLHCIFData* grpo = ccdb->template getForTimeStamp<o2::parameters::GRPLHCIFData>(pidTPCopts.cfgPathGrpLhcIf.value, bc.timestamp());
        if (grpo) {
          LOG(info) << "Collision type::" << CollisionSystemType::getCollisionTypeFromGrp(grpo);
          collsys = CollisionSystemType::getCollisionTypeFromGrp(grpo);
          if (collsys == CollisionSystemType::kCollSyspp) {
            irSource = std::string("T0VTX");
          } else {
            irSource = std::string("ZNC hadronic");
          }
        } else {
          LOGF(info, "No grpo object found. irSource will remain undefined.");
        }
        response->PrintAll();
      }

      if (bc.timestamp() < network.getValidityFrom() || bc.timestamp() > network.getValidityUntil()) { // fetches network only if the runnumbers change
        LOG(info) << "Fetching network for timestamp: " << bc.timestamp();
        bool retrieveSuccess = ccdb->getCCDBAccessor().retrieveBlob(pidTPCopts.networkPathCCDB.value, ".", metadata, bc.timestamp(), false, pidTPCopts.networkPathLocally.value, "", "", &headers);
        networkVersion = headers["NN-Version"];
        if (retrieveSuccess) {
          network.initModel(pidTPCopts.networkPathLocally.value, pidTPCopts.enableNetworkOptimizations.value, pidTPCopts.networkSetNumThreads.value, strtoul(headers["Valid-From"].c_str(), NULL, 0), strtoul(headers["Valid-Until"].c_str(), NULL, 0));
          std::vector<float> dummyInput(network.getNumInputNodes(), 1.);
          network.evalModel(dummyInput);
          setupColumnInputNetwork();
          LOGP(info, "Retrieved NN corrections for production tag {}, pass number {}, NN-Version number {}", headers["LPMProductionTag"], headers["RecoPassName"], headers["NN-Version"]);
        } else {
          LOG(fatal) << "No valid NN object found matching retrieved Bethe-Bloch parametrisation for pass " << metadata["RecoPassName"] << ". Please ensure that the requested pass has dedicated NN corrections available";
        }
      }
    }

    // Defining some network parameters
    const int nFeat = network.getNumFeatures();
    int output_dimensions = network.getNumOutputNodes();
    const uint64_t prediction_size = output_dimensions * size;

    network_prediction = std::vector<float>(prediction_size * 9); // For each mass hypotheses
    const float nNclNormalization = response->GetNClNormalization();
    float duration_network = 0;

    // To load the Hadronic rate once for each collision
    float hadronicRateBegin = 0.;
    std::vector<float> hadronicRateForCollision(collisions.size(), 0.0f);
    size_t i = 0;
    for (const auto& collision : collisions) {
      const auto& bc = collision.template bc_as<B>();
      if (irSource.compare("") != 0) {
        hadronicRateForCollision[i] = mRateFetcher.fetch(ccdb.service, bc.timestamp(), bc.runNumber(), irSource) * 1.e-3;
      } else {
        hadronicRateForCollision[i] = 0.0f;
      }
      i++;
    }
    auto bc = bcs.begin();
    if (irSource.compare("") != 0) {
      hadronicRateBegin = mRateFetcher.fetch(ccdb.service, bc.timestamp(), bc.runNumber(), irSource) * 1.e-3; // kHz
    } else {
      hadronicRateBegin = 0.0f;
    }

    static constexpr int NParticleTypes = 9;
    constexpr int ExpectedInputDimensionsNNV2 = 7;
    constexpr int ExpectedInputDimensionsNNV3 = 8;
    constexpr int ExpectedInputDimensionsNNV4 = 9;

    const float hadronicRateDivisor = (collsys == CollisionSystemType::kCollSyspp) ? 1500.f : 50.f;

    // Per-collision arrays (O(nColl)); gathered per track inside the model via the
    // collisionId column, then normalised in-graph.
    const int64_t nColl = static_cast<int64_t>(collisions.size());
    std::vector<float> multArray(nColl);
    std::vector<float> occArray(nFeat >= ExpectedInputDimensionsNNV2 ? nColl : 0);
    {
      int64_t c = 0;
      for (const auto& col : collisions) {
        multArray[c] = static_cast<float>(mults[c]);
        if (nFeat >= ExpectedInputDimensionsNNV2) {
          occArray[c] = col.ft0cOccupancyInTimeRange();
        }
        ++c;
      }
    }

    // Raw per-track Arrow column buffers (zero-copy; one chunk per DataFrame).
    auto arrowTable = tracks.asArrowTable();
    auto chunk0 = [&](const char* name) -> std::shared_ptr<arrow::Array> {
      const int idx = arrowTable->schema()->GetFieldIndex(name);
      if (idx < 0) {
        LOG(fatal) << "createNetworkPrediction: column '" << name << "' not found in tracks table";
      }
      auto col = arrowTable->column(idx);
      if (col->num_chunks() != 1) {
        LOG(fatal) << "createNetworkPrediction: column '" << name << "' has " << col->num_chunks()
                   << " chunks; a single chunk per DataFrame is required for zero-copy input";
      }
      return col->chunk(0);
    };
    const int64_t nTrk = static_cast<int64_t>(tracks.size());
    const float* pTpcInner = std::static_pointer_cast<arrow::FloatArray>(chunk0("fTPCInnerParam"))->raw_values();
    const float* pTgl = std::static_pointer_cast<arrow::FloatArray>(chunk0("fTgl"))->raw_values();
    const float* pSigned1Pt = std::static_pointer_cast<arrow::FloatArray>(chunk0("fSigned1Pt"))->raw_values();
    const int32_t* pCollId = std::static_pointer_cast<arrow::Int32Array>(chunk0("fIndexCollisions"))->raw_values();
    const uint8_t* pFindable = std::static_pointer_cast<arrow::UInt8Array>(chunk0("fTPCNClsFindable"))->raw_values();
    const int8_t* pFMF = std::static_pointer_cast<arrow::Int8Array>(chunk0("fTPCNClsFindableMinusFound"))->raw_values();
    const float* pPhi = (nFeat >= ExpectedInputDimensionsNNV4)
                          ? std::static_pointer_cast<arrow::FloatArray>(chunk0("fPhi"))->raw_values()
                          : nullptr;

    // Single boolean mask of the tracks the network runs on; the model Compress'es
    // to exactly these rows so the output is compact and the consumer's
    // count_tracks indexing is unchanged. Condition matches process()'s counter.
    std::vector<uint8_t> validMask(nTrk);
    {
      int64_t t = 0;
      for (auto const& trk : tracks) {
        bool valid = trk.hasTPC();
        if (valid && pidTPCopts.skipTPCOnly && !trk.hasITS() && !trk.hasTRD() && !trk.hasTOF()) {
          valid = false;
        }
        validMask[t++] = valid ? 1 : 0;
      }
    }

    auto memInfo = Ort::MemoryInfo::CreateCpu(OrtAllocatorType::OrtArenaAllocator, OrtMemType::OrtMemTypeDefault);
    const int64_t one = 1;
    float massVal = 0.f;
    float nclNormVal = nNclNormalization;
    float hrDivisorVal = hadronicRateDivisor;
    float hrFallbackVal = hadronicRateBegin / hadronicRateDivisor;

    // Evaluate once per mass hypothesis; only the mass scalar input changes.
    for (int j = 0; j < NParticleTypes; j++) {
      massVal = o2::track::pid_constants::sMasses[j];

      std::vector<Ort::Value> inputTensors;
      inputTensors.emplace_back(Ort::Value::CreateTensor<float>(memInfo, const_cast<float*>(pTpcInner), nTrk, &nTrk, 1));
      inputTensors.emplace_back(Ort::Value::CreateTensor<float>(memInfo, const_cast<float*>(pTgl), nTrk, &nTrk, 1));
      inputTensors.emplace_back(Ort::Value::CreateTensor<float>(memInfo, const_cast<float*>(pSigned1Pt), nTrk, &nTrk, 1));
      inputTensors.emplace_back(Ort::Value::CreateTensor<float>(memInfo, &massVal, 1, &one, 1));
      inputTensors.emplace_back(Ort::Value::CreateTensor<int32_t>(memInfo, const_cast<int32_t*>(pCollId), nTrk, &nTrk, 1));
      inputTensors.emplace_back(Ort::Value::CreateTensor<float>(memInfo, multArray.data(), nColl, &nColl, 1));
      inputTensors.emplace_back(Ort::Value::CreateTensor<float>(memInfo, &nclNormVal, 1, &one, 1));
      inputTensors.emplace_back(Ort::Value::CreateTensor<uint8_t>(memInfo, const_cast<uint8_t*>(pFindable), nTrk, &nTrk, 1));
      inputTensors.emplace_back(Ort::Value::CreateTensor<int8_t>(memInfo, const_cast<int8_t*>(pFMF), nTrk, &nTrk, 1));
      if (nFeat >= ExpectedInputDimensionsNNV2) {
        inputTensors.emplace_back(Ort::Value::CreateTensor<float>(memInfo, occArray.data(), nColl, &nColl, 1));
      }
      if (nFeat >= ExpectedInputDimensionsNNV3) {
        inputTensors.emplace_back(Ort::Value::CreateTensor<float>(memInfo, hadronicRateForCollision.data(), nColl, &nColl, 1));
        inputTensors.emplace_back(Ort::Value::CreateTensor<float>(memInfo, &hrDivisorVal, 1, &one, 1));
        inputTensors.emplace_back(Ort::Value::CreateTensor<float>(memInfo, &hrFallbackVal, 1, &one, 1));
      }
      if (nFeat >= ExpectedInputDimensionsNNV4) {
        inputTensors.emplace_back(Ort::Value::CreateTensor<float>(memInfo, const_cast<float*>(pPhi), nTrk, &nTrk, 1));
      }
      inputTensors.emplace_back(Ort::Value::CreateTensor<bool>(memInfo, reinterpret_cast<bool*>(validMask.data()), nTrk, &nTrk, 1));

      auto start_network_eval = std::chrono::high_resolution_clock::now();
      float* output_network = network.evalModel<float>(inputTensors);
      auto stop_network_eval = std::chrono::high_resolution_clock::now();
      duration_network += std::chrono::duration<float, std::ratio<1, 1000000000>>(stop_network_eval - start_network_eval).count();
      for (uint64_t k = 0; k < prediction_size; k += output_dimensions) {
        for (int l = 0; l < output_dimensions; l++) {
          network_prediction[k + l + prediction_size * j] = output_network[k + l];
        }
      }
    }
    auto stop_network_total = std::chrono::high_resolution_clock::now();
    LOG(debug) << "Neural Network for the TPC PID response correction: Time per track (eval ONNX): " << duration_network / (size * 9) << "ns ; Total time (eval ONNX): " << duration_network / 1000000000 << " s";
    LOG(debug) << "Neural Network for the TPC PID response correction: Time per track (eval + overhead): " << std::chrono::duration<float, std::ratio<1, 1000000000>>(stop_network_total - start_network_total).count() / (size * 9) << "ns ; Total time (eval + overhead): " << std::chrono::duration<float, std::ratio<1, 1000000000>>(stop_network_total - start_network_total).count() / 1000000000 << " s";

    return network_prediction;
  }

  //__________________________________________________
  template <typename T, typename NSF, typename NST>
  void makePidTables(const int flagFull, NSF& tableFull, const int flagTiny, NST& tableTiny, const o2::track::PID::ID pid, const float tpcSignal, const T& trk, const int64_t multTPC, const std::vector<float>& network_prediction, const int& count_tracks, const int& tracksForNet_size)
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
    auto expSigma = trk.has_collision() ? response->GetExpectedSigmaAtMultiplicity(multTPC, trk, pid) : 0.07 * expSignal; // use default sigma value of 7% if no collision information to estimate resolution
    if (expSignal < 0. || expSigma < 0.) {                                                                                // skip if expected signal invalid
      if (flagFull)
        tableFull(-999.f, -999.f);
      if (flagTiny)
        tableTiny(aod::pidtpc_tiny::binning::underflowBin);
      return;
    }

    float nSigma = -999.f;
    float bg = trk.tpcInnerParam() / o2::track::pid_constants::sMasses[pid]; // estimated beta-gamma for network cutoff
    constexpr int NumOutputNodesSymmetricSigma = 2;
    constexpr int NumOutputNodesAsymmetricSigma = 3;
    if (pidTPCopts.useNetworkCorrection && speciesNetworkFlags[pid] && trk.has_collision() && bg > pidTPCopts.networkBetaGammaCutoff) {

      // Here comes the application of the network. The output--dimensions of the network determine the application: 1: mean, 2: sigma, 3: sigma asymmetric
      // For now only the option 2: sigma will be used. The other options are kept if there would be demand later on
      if (network.getNumOutputNodes() == 1) { // Expected mean correction; no sigma correction
        nSigma = (tpcSignal - network_prediction[count_tracks + tracksForNet_size * pid] * expSignal) / expSigma;
      } else if (network.getNumOutputNodes() == NumOutputNodesSymmetricSigma) { // Symmetric sigma correction
        expSigma = (network_prediction[NumOutputNodesSymmetricSigma * (count_tracks + tracksForNet_size * pid) + 1] - network_prediction[NumOutputNodesSymmetricSigma * (count_tracks + tracksForNet_size * pid)]) * expSignal;
        nSigma = (tpcSignal / expSignal - network_prediction[NumOutputNodesSymmetricSigma * (count_tracks + tracksForNet_size * pid)]) / (network_prediction[NumOutputNodesSymmetricSigma * (count_tracks + tracksForNet_size * pid) + 1] - network_prediction[NumOutputNodesSymmetricSigma * (count_tracks + tracksForNet_size * pid)]);
      } else if (network.getNumOutputNodes() == NumOutputNodesAsymmetricSigma) { // Asymmetric sigma corection
        if (tpcSignal / expSignal >= network_prediction[NumOutputNodesAsymmetricSigma * (count_tracks + tracksForNet_size * pid)]) {
          expSigma = (network_prediction[NumOutputNodesAsymmetricSigma * (count_tracks + tracksForNet_size * pid) + 1] - network_prediction[NumOutputNodesAsymmetricSigma * (count_tracks + tracksForNet_size * pid)]) * expSignal;
          nSigma = (tpcSignal / expSignal - network_prediction[NumOutputNodesAsymmetricSigma * (count_tracks + tracksForNet_size * pid)]) / (network_prediction[NumOutputNodesAsymmetricSigma * (count_tracks + tracksForNet_size * pid) + 1] - network_prediction[NumOutputNodesAsymmetricSigma * (count_tracks + tracksForNet_size * pid)]);
        } else {
          expSigma = (network_prediction[NumOutputNodesAsymmetricSigma * (count_tracks + tracksForNet_size * pid)] - network_prediction[NumOutputNodesAsymmetricSigma * (count_tracks + tracksForNet_size * pid) + 2]) * expSignal;
          nSigma = (tpcSignal / expSignal - network_prediction[NumOutputNodesAsymmetricSigma * (count_tracks + tracksForNet_size * pid)]) / (network_prediction[NumOutputNodesAsymmetricSigma * (count_tracks + tracksForNet_size * pid)] - network_prediction[NumOutputNodesAsymmetricSigma * (count_tracks + tracksForNet_size * pid) + 2]);
        }
      } else {
        LOGF(fatal, "Network output dimensions incompatible!");
      }
    } else {
      nSigma = response->GetNumberOfSigmaMCTunedAtMultiplicity(multTPC, trk, pid, tpcSignal);
    }
    if (flagFull)
      tableFull(expSigma, nSigma);
    if (flagTiny)
      aod::pidtpc_tiny::binning::packInTable(nSigma, tableTiny);
  };

  //__________________________________________________
  template <typename TCCDB, typename TBCs, typename TTracks, typename TTracksQA, typename TProducts>
  void process(TCCDB& ccdb, TBCs const& bcs, soa::Join<aod::Collisions, aod::EvSels> const& cols, TTracks const& tracks, TTracksQA const& tracksQA, TProducts& products)
  {
    if (tracks.size() == 0) {
      return; // empty protection
    }
    auto trackiterator = tracks.begin();
    if constexpr (requires { trackiterator.mcParticleId(); }) {
      gRandom->SetSeed(0); // Ensure unique seed from UUID for each process call
    }

    // preparatory step: we need the multiplicities for each collision
    std::vector<int64_t> pidmults;
    int64_t totalTPCtracks = 0;
    int64_t totalTPCnotStandalone = 0;
    pidmults.resize(cols.size(), 0);

    // faster counting
    for (const auto& track : tracks) {
      if (track.hasTPC()) {
        if (track.collisionId() > -1) {
          pidmults[track.collisionId()]++;
        }
        totalTPCtracks++;
        if (track.hasITS() || track.hasTOF() || track.hasTRD()) {
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
      network_prediction = createNetworkPrediction(ccdb, cols, pidmults, tracks, bcs, tracksForNet_size);
    }

    uint64_t count_tracks = 0;

    //_______________________________________
    // process tracksQA in case present
    std::vector<int64_t> indexTrack2TrackQA(outTable_size, -1);
    if constexpr (soa::is_table<TTracksQA>) {
      for (const auto& trackQA : tracksQA) {
        indexTrack2TrackQA[trackQA.trackId()] = trackQA.globalIndex();
      }
    }
    //_______________________________________

    // Fill Hadronic rate per collision in case CorrectedDEdx is requested
    std::vector<float> hadronicRateForCollision(cols.size(), 0.0f);
    float hadronicRateBegin = 0.0f;
    if (pidTPCopts.useCorrecteddEdx) {
      size_t i = 0;
      for (const auto& collision : cols) {
        const auto& bc = collision.template bc_as<aod::BCsWithTimestamps>();
        if (irSource.compare("") != 0) {
          hadronicRateForCollision[i] = mRateFetcher.fetch(ccdb.service, bc.timestamp(), bc.runNumber(), irSource) * 1.e-3;
        } else {
          hadronicRateForCollision[i] = 0.0f;
        }
        i++;
      }
      auto bc = bcs.begin();
      if (irSource.compare("") != 0) {
        hadronicRateBegin = mRateFetcher.fetch(ccdb.service, bc.timestamp(), bc.runNumber(), irSource) * 1.e-3; // kHz
      } else {
        hadronicRateBegin = 0.0f;
      }
    }

    for (auto const& trk : tracks) {
      // get the TPC signal to be used in the PID
      float tpcSignalToEvaluatePID = trk.tpcSignal();

      int64_t multTPC = 0;
      if (trk.has_collision()) {
        multTPC = pidmults[trk.collisionId()];
      }

      // if corrected dE/dx is requested, correct it here on the spot and use that
      if (pidTPCopts.useCorrecteddEdx) {

        //_________________________________________________________
        // bypass TPC signal in case TracksQA information present
        if constexpr (soa::is_table<TTracksQA>) {
          tpcSignalToEvaluatePID = -999.f;
          if (indexTrack2TrackQA[trk.globalIndex()] != -1) {
            auto trackQA = tracksQA.rawIteratorAt(indexTrack2TrackQA[trk.globalIndex()]);
            tpcSignalToEvaluatePID = trackQA.tpcdEdxNorm();
          }
        }
        //_________________________________________________________

        double hadronicRate;
        int occupancy;
        if (trk.has_collision()) {
          auto collision = cols.iteratorAt(trk.collisionId());
          hadronicRate = hadronicRateForCollision[trk.collisionId()];
          occupancy = collision.trackOccupancyInTimeRange();
        } else {
          hadronicRate = hadronicRateBegin;
          occupancy = 0;
        }

        constexpr float kExpectedTPCSignalMIP = 50.0f;
        constexpr float kMaxAllowedRatio = 1.05f;
        constexpr float kMinAllowedRatio = 0.05f;
        constexpr float kMaxAllowedOcc = 12.0f;

        float fTPCSignal = tpcSignalToEvaluatePID;
        float fNormMultTPC = multTPC / 11000.;

        float fTrackOccN = occupancy / 1000.;
        float fOccTPCN = fNormMultTPC * 10; //(fNormMultTPC*10).clip(0,12)
        if (fOccTPCN > kMaxAllowedOcc)
          fOccTPCN = kMaxAllowedOcc;
        else if (fOccTPCN < 0)
          fOccTPCN = 0;

        float fTrackOccMeanN = hadronicRate / 5;
        float side = trk.tgl() > 0 ? 1 : 0;
        float a1pt = std::abs(trk.signed1Pt());
        float a1pt2 = a1pt * a1pt;
        float atgl = std::abs(trk.tgl());
        float mbb0R = kExpectedTPCSignalMIP / fTPCSignal;
        if (mbb0R > kMaxAllowedRatio)
          mbb0R = kMaxAllowedRatio;
        else if (mbb0R < kMinAllowedRatio)
          mbb0R = kMinAllowedRatio;
        // float mbb0R =  max(0.05,  min(50 / fTPCSignal, 1.05));
        float a1ptmbb0R = a1pt * mbb0R;
        float atglmbb0R = atgl * mbb0R;

        std::vector<float> vec_occu = {fTrackOccN, fOccTPCN, fTrackOccMeanN};
        std::vector<float> vec_track = {mbb0R, a1pt, atgl, atglmbb0R, a1ptmbb0R, side, a1pt2};

        float fTPCSignalN_CR0 = str_dedx_correction.fReal_fTPCSignalN(vec_occu, vec_track);

        float mbb0R1 = kExpectedTPCSignalMIP / (fTPCSignal / fTPCSignalN_CR0);
        if (mbb0R1 > kMaxAllowedRatio)
          mbb0R1 = kMaxAllowedRatio;
        else if (mbb0R1 < kMinAllowedRatio)
          mbb0R1 = kMinAllowedRatio;

        std::vector<float> vec_track1 = {mbb0R1, a1pt, atgl, atgl * mbb0R1, a1pt * mbb0R1, side, a1pt2};
        float fTPCSignalN_CR1 = str_dedx_correction.fReal_fTPCSignalN(vec_occu, vec_track1);

        // change the signal used for PID
        tpcSignalToEvaluatePID = fTPCSignal / fTPCSignalN_CR1;

        if (pidTPCopts.savedEdxsCorrected) {
          // populated cursor if requested or autodetected
          products.dEdxCorrected(tpcSignalToEvaluatePID);
        }
      }

      const auto& bc = trk.has_collision() ? cols.rawIteratorAt(trk.collisionId()).template bc_as<aod::BCsWithTimestamps>() : bcs.begin();
      if (useCCDBParam && pidTPCopts.ccdbTimestamp.value == 0 && !ccdb->isCachedObjectValid(pidTPCopts.ccdbPath.value, bc.timestamp())) { // Updating parametrisation only if the initial timestamp is 0
        if (pidTPCopts.recoPass.value == "") {
          LOGP(info, "Retrieving latest TPC response object for timestamp {}:", bc.timestamp());
        } else {
          LOGP(info, "Retrieving TPC Response for timestamp {} and recoPass {}:", bc.timestamp(), pidTPCopts.recoPass.value);
        }
        response = ccdb->template getSpecific<o2::pid::tpc::Response>(pidTPCopts.ccdbPath.value, bc.timestamp(), metadata, &headers);
        if (!response) {
          LOGP(warning, "!! Could not find a valid TPC response object for specific pass name {}! Falling back to latest uploaded object.", metadata["RecoPassName"]);
          response = ccdb->template getForTimeStamp<o2::pid::tpc::Response>(pidTPCopts.ccdbPath.value, bc.timestamp(), &headers);
          if (!response) {
            LOGP(fatal, "Could not find ANY TPC response object for the timestamp {}!", bc.timestamp());
          }
        }
        LOG(info) << "Successfully retrieved TPC PID object from CCDB for timestamp " << bc.timestamp() << ", period " << headers["LPMProductionTag"] << ", recoPass " << headers["RecoPassName"];
        o2::parameters::GRPLHCIFData* grpo = ccdb->template getForTimeStamp<o2::parameters::GRPLHCIFData>(pidTPCopts.cfgPathGrpLhcIf.value, bc.timestamp());
        if (grpo) {
          LOG(info) << "Collisions type::" << CollisionSystemType::getCollisionTypeFromGrp(grpo);
          collsys = CollisionSystemType::getCollisionTypeFromGrp(grpo);
          if (collsys == CollisionSystemType::kCollSyspp) {
            irSource = std::string("T0VTX");
          } else {
            irSource = std::string("ZNC hadronic");
          }
        } else {
          LOGF(info, "No grpo object found. irSource will remain undefined.");
        }
        response->PrintAll();
      }

      // if this is a MC process function, go for MC tune on data processing
      if constexpr (requires { trk.mcParticleId(); }) {
        // Perform TuneOnData sampling for MC dE/dx
        if (!trk.has_mcParticle()) {
          products.tableTuneOnData(-999.f);
          tpcSignalToEvaluatePID = -999.f; // pass this for further eval
        } else {
          float mcTunedTPCSignal = 0.;
          if (!trk.hasTPC()) {
            mcTunedTPCSignal = -999.f;
          } else {
            if (pidTPCopts.skipTPCOnly) {
              if (!trk.hasITS() && !trk.hasTRD() && !trk.hasTOF()) {
                mcTunedTPCSignal = -999.f;
              }
            }
            int pid = getPIDIndex(trk.mcParticle().pdgCode());

            auto expSignal = response->GetExpectedSignal(trk, pid);
            auto expSigma = response->GetExpectedSigmaAtMultiplicity(multTPC, trk, pid);
            if (expSignal < 0. || expSigma < 0.) { // if expectation invalid then give undefined signal
              mcTunedTPCSignal = -999.f;
            }
            float bg = trk.tpcInnerParam() / o2::track::pid_constants::sMasses[pid]; // estimated beta-gamma for network cutoff

            if (pidTPCopts.useNetworkCorrection && speciesNetworkFlags[pid] && trk.has_collision() && bg > pidTPCopts.networkBetaGammaCutoff) {
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
          tpcSignalToEvaluatePID = mcTunedTPCSignal; // pass this for further eval
          if (pidTPCopts.enableTuneOnDataTable)
            products.tableTuneOnData(mcTunedTPCSignal);
        }
      }

      auto makePidTablesDefault = [&trk, &tpcSignalToEvaluatePID, &multTPC, &network_prediction, &count_tracks, &tracksForNet_size, this](const int flagFull, auto& tableFull, const int flagTiny, auto& tableTiny, const o2::track::PID::ID pid) {
        this->makePidTables(flagFull, tableFull, flagTiny, tableTiny, pid, tpcSignalToEvaluatePID, trk, multTPC, network_prediction, count_tracks, tracksForNet_size);
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

#endif // COMMON_TOOLS_PID_PIDTPCMODULE_H_
