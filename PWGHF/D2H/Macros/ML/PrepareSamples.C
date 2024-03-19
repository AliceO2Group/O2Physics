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

/// \file PrepareSamples.C
/// \brief Macro to prepare data samples for ML training
/// \author Alexandre Bigot <alexandre.bigot@cern.ch>, Strasbourg University

// if .h file not found, please include your local rapidjson/document.h and rapidjson/filereadstream.h here
#include <rapidjson/document.h>
#include <rapidjson/filereadstream.h>
#include <typeinfo>

#include "ROOT/RDataFrame.hxx"

#include <TChain.h>
#include <TFile.h>
#include <TString.h>
#include <TTree.h>

// define constants
static constexpr uint8_t idBkg = 0;
static constexpr uint8_t idPrompt = 1;
static constexpr uint8_t idNonPrompt = 2;
static constexpr std::array<uint8_t, 3> ids = {idBkg, idPrompt, idNonPrompt};

static constexpr float nSigTolerance = 0.1;
static constexpr float nSigDummy = -999. + nSigTolerance;


///Function to retrieve information from JSON in a vector<T>
/// \param jsonEntry is the entry loaded from JSON
/// \return vec a vector of type T
template <typename T, typename TypeJsonEntry>
const std::vector<T> getVectorFromJson(const TypeJsonEntry& jsonEntry)
{
  std::vector<T> vec;
  const auto& jsonArray = jsonEntry.GetArray();

  if constexpr (std::is_same<T, std::string>::value) {
    for (const auto& entry : jsonArray) {
      vec.emplace_back(entry.GetString());
    }
  } else if constexpr (std::is_same<T, TString>::value) {
    for (const auto& entry : jsonArray) {
      vec.emplace_back(TString(entry.GetString()));
    }
  } else if constexpr (std::is_same<T, uint8_t>::value) {
    for (const auto& entry : jsonArray) {
      vec.emplace_back((uint8_t)entry.GetInt());
    }
  } else if constexpr (std::is_same<T, int>::value) {
    for (const auto& entry : jsonArray) {
      vec.emplace_back((int)entry.GetInt());
    }
  } else if constexpr (std::is_same<T, float>::value) {
    for (const auto& entry : jsonArray) {
      vec.emplace_back(entry.GetFloat());
    }
  } else if constexpr (std::is_same<T, double>::value) {
    for (const auto& entry : jsonArray) {
      vec.emplace_back(entry.GetDouble());
    }
  } else {
    std::cerr << "ERROR: type not available!" << std::endl;
  }

  return vec;
}

/// Function to remove track if it has no TOF nor TPC information and to compute combined nSigma
/// \param nSigTpc nSigma of TPC
/// \param nSigTof nSigma of TOF
/// \return combined nSigma
float combinePid(float nSigTpc, float nSigTof)
{
  bool hasTpc = true ? nSigTpc > nSigDummy : false;
  bool hasTof = true ? nSigTof > nSigDummy : false;

  if (hasTpc && hasTof) { // TPC and TOF
    return std::sqrt(0.5 * (nSigTpc*nSigTpc + nSigTof*nSigTof));
  } else if (hasTpc) { // TPC only
    return std::abs(nSigTpc);
  } else if (hasTof) { // TOF only
    return std::abs(nSigTpc);
  }

  return -1.f; // no TPC nor TOF
}

/// Function to get the path of the input TTrees
/// \param nameInputFiles a vector with input files names
/// \param nameTree the (common) name of the TTree
/// \return pathInputTrees a vector with paths to the input TTrees
template <typename T1, typename T2>
std::vector<TString> getPathInputTrees(const T1& nameInputFiles, const T2& nameTree)
{
  const std::string nameDf = "DF";
  const std::string delimiter = ";";
  std::vector<TString> pathInputTrees;
  for (const auto& nameFile : nameInputFiles) {
    std::unique_ptr<TFile> myFile( TFile::Open(nameFile.data()) );
    std::vector<std::string> nameDirs;
    for (const auto& key : *myFile->GetListOfKeys()) {
      std::string name = key->GetName();
      // fill vector only if DF in the dir name (no parent files dir)
      if (name.find(nameDf) != std::string::npos ) {
        nameDirs.emplace_back(name);
      }
    }
    myFile->Close();

    if (nameDirs.size() > 2) {
      std::cerr << "ERROR: too many DFs in <name of root file>! Run o2-aod-merger with --max-size <great number> to get only one DF." << std::endl;
      return std::vector<TString>{""};
    } else if (nameDirs.size() == 2) {
      // Check if the two DFs are be the same (o2-aod-merger artefact)
      auto name = nameDirs[0];
      std::string name0 = nameDirs[0].substr(0, nameDirs[0].find(delimiter));
      std::string name1 = nameDirs[1].substr(0, nameDirs[1].find(delimiter));

      if (name0 != name1) {
        std::cerr << "ERROR: too many DFs in " << nameFile << " ! Run o2-aod-merger with --max-size <great number> --output AO2D_merged.root to get only one DF. Run o2-aod-merger --help for more information" << std::endl;
        return std::vector<TString>{""};
      }
    }

    pathInputTrees.emplace_back(TString(nameFile + "/" + nameDirs[0] + "/" + nameTree)); // TODO add enforceTrailingSlash method
  }

  return pathInputTrees;
}

/// Main function
/// \param nameCfgFile name of the JSON confifuration file
void PrepareSamples(TString nameCfgFile="./config_preparation_DplusToPiKPi.json")
{
  // load configuration file
  FILE* configFile = fopen(nameCfgFile.Data(), "r");
  if (!configFile) {
    std::cerr << "ERROR: Missing configuration json file: " << nameCfgFile << std::endl;
    return;
  }

  rapidjson::Document config;
  char readBuffer[65536];
  rapidjson::FileReadStream is(configFile, readBuffer, sizeof(readBuffer));
  config.ParseStream(is);
  fclose(configFile);

  // import configurables
  std::string channel = config["channel"].GetString();

  const auto labels = getVectorFromJson<std::string>(config["labels"]);
  const auto nameInputFiles = getVectorFromJson<std::string>(config["prepare_samples"]["input"]["files"]);
  const std::string nameTree = config["prepare_samples"]["input"]["tree_name"].GetString();
  const auto preSelections = getVectorFromJson<std::string>(config["preselections"]);

  const bool doPidCombination = config["pid"]["combine_vars"].GetBool();
  const std::vector<TString> massHypos = getVectorFromJson<TString>(config["pid"]["mass_hypos"]);
  const uint8_t nProngs = (uint8_t)config["pid"]["n_prongs"].GetInt();

  std::vector<std::string> colsToRemove;
  if (!config["cols_to_remove"].IsNull()) {
    colsToRemove = getVectorFromJson<std::string>(config["cols_to_remove"]);
  }

  const float downscaleBkg = config["downscale_bkg"].GetFloat();
  const std::string invMassSideBands = config["filt_bkg_mass"].GetString();
  const uint8_t seedSplit = (uint8_t)config["seed_split"].GetInt();

  // output
  const bool force = config["output"]["force"].GetBool();
  const auto nameOutputDirs = getVectorFromJson<std::string>(config["output"]["dirs"]);
  const TString nameOutputTree = config["output"]["tree_name"].GetString();

  // configure access to the input TTrees
  auto pathInputTrees = getPathInputTrees(nameInputFiles, nameTree);
  // safety
  if (pathInputTrees.empty()) {
    return;
  }

  uint8_t counter_outdir = 0;
  for (const auto& pathInputTree : pathInputTrees) {
    TChain chain;
    chain.Add(pathInputTree);

    // define dataframe from the input TTrees
    ROOT::EnableImplicitMT(32); // tell ROOT to go parallel
    ROOT::RDataFrame dataFrame(chain);

    std::vector<std::string> colsToKeep = dataFrame.GetColumnNames();
    if (colsToRemove.size() == 0) {
      std::remove(colsToKeep.begin(), colsToKeep.end(), "fOriginMcRec"); // just remove cand_type if not explicitly kept
      colsToKeep.pop_back();
    }

    // apply preselections, if enabled
    auto dfPreselDummy = dataFrame.Filter("fM > 0"); // trick to get a filtered dataframe type
    std::array<decltype(dfPreselDummy), 1> dataFrames{dfPreselDummy};
    if (preSelections.size() != 0) {
      for (const auto& preSelection : preSelections) {
        dataFrames[0] = dataFrame.Filter(preSelection);
      }
    }

    // combine PID variables, if enabled
    if (doPidCombination) {

      const TString fNSigTpc("fNSigTpc");
      const TString fNSigTof("fNSigTof");
      const TString fNSigTpcTof("fNSigTpcTof");

      std::vector<TString> suffixPid;
      for (const auto& massHypo : massHypos) {
        for (uint8_t iProng{0}; iProng < nProngs; ++iProng) {
          suffixPid.emplace_back(massHypo + TString::Format("%d", iProng));
        }
      }

      for (const auto& suffix : suffixPid) {
        // define column names
        TString nameColTpc = fNSigTpc + suffix;
        TString nameColTof = fNSigTpc + suffix;
        TString nameColTpcTof = fNSigTpcTof + suffix;
        TString filterCombinedNSig = nameColTpcTof + ">= 0"; // apply TPC or TOF logic
        // compute combined nSigma
        auto df = dataFrames[0].Define(nameColTpcTof, combinePid, {nameColTpc.Data(), nameColTof.Data()});
        dataFrames[0] = df.Filter(filterCombinedNSig.Data());

        colsToKeep.emplace_back(nameColTpcTof);
      }
    }

    // define total preselected dataframe
    auto dfTot = dataFrames[0];

    // divide dataframe into classes and save them in flagged .root files
    for (const auto& id : ids) {
      if (id == idBkg) {
        dfTot.Filter(TString::Format("fOriginMcRec == %d", id).Data()).Filter(invMassSideBands)
          .Snapshot(nameOutputTree, TString::Format("%s/%s_%s.root", nameOutputDirs[counter_outdir].data(), labels[id].data(), channel.data()), colsToKeep);
      } else {
        dfTot.Filter(TString::Format("fOriginMcRec == %d", id).Data())
          .Snapshot(nameOutputTree, TString::Format("%s/%s_%s.root", nameOutputDirs[counter_outdir].data(), labels[id].data(), channel.data()), colsToKeep);
      }
    }
    counter_outdir++;
  }
}
