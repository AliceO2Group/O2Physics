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

/// \file   testMetadataHelper.C
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \brief  Test the MetadataHelper functionality

#include "Common/Core/MetadataHelper.h"

#include <CCDB/CcdbApi.h>
#include <Framework/ConfigContext.h>
#include <Framework/ConfigParamRegistry.h>
#include <Framework/ConfigParamStore.h>
#include <Framework/ServiceRegistry.h>
#include <Framework/ServiceRegistryRef.h>

#include <TFile.h>
#include <TMap.h>
#include <TObjString.h>
#include <TSystem.h>

#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

// Taken from O2/Framework/AnalysisSupport/src/Plugin.cxx
auto readMetadata(std::unique_ptr<TFile>& currentFile) -> std::vector<o2::framework::ConfigParamSpec>
{
  // Get the metadata, if any
  auto m = (TMap*)currentFile->Get("metaData");
  if (!m) {
    return {};
  }
  std::vector<o2::framework::ConfigParamSpec> results;
  auto it = m->MakeIterator();

  // Serialise metadata into a ; separated string with : separating key and value
  bool first = true;
  while (auto obj = it->Next()) {
    if (first) {
      LOGP(info, "Metadata for file \"{}\":", currentFile->GetName());
      first = false;
    }
    auto objString = (TObjString*)m->GetValue(obj);
    std::string key = "aod-metadata-" + std::string(obj->GetName());
    LOGP(info, "- {}: {} goes into key {}", obj->GetName(), objString->String().Data(), key);
    char const* value = strdup(objString->String());
    results.push_back(o2::framework::ConfigParamSpec{key, o2::framework::VariantType::String, value, {"Metadata in AOD"}});
  }
  return results;
}

// Create a file with all the versions of the O2 software with alienv q
void createO2VersionFile()
{
  // Can do this only if on lxplus
  std::string host = gSystem->HostName() ? gSystem->HostName() : "";
  if (host.find("lxplus") == std::string::npos) {
    LOG(warn) << "Not on lxplus (" << host << "); skipping creation of /tmp/o2version.txt";
    return;
  }
  // If file exists, do nothing
  std::ifstream infile("/tmp/o2version.txt");
  if (infile.is_open()) {
    return;
  }
  gSystem->Exec("alienv q | grep VO_ALICE@O2:: > /tmp/o2version.txt");
}

std::map<std::string, bool> buildMapForCommitHash(const std::string& hash)
{
  // Change directory to /tmp
  std::map<std::string, bool> results;
  std::ifstream infileO2Versions("/tmp/o2version.txt");
  std::string lineOfO2Version;
  const std::string fileContainingCommit = "/tmp/branches_" + hash + ".txt";
  std::ifstream infileO2VersionsWithHash(fileContainingCommit);
  if (!infileO2VersionsWithHash.is_open()) {
    gSystem->cd("/tmp/");
    gSystem->Exec("git clone git@github.com:AliceO2Group/AliceO2.git");
    gSystem->cd("AliceO2");
    std::string cmd = Form("git branch -r --contains %s > %s 2>&1", hash.c_str(), fileContainingCommit.c_str());
    LOG(info) << "Executing command " << cmd;
    gSystem->Exec(cmd.c_str());
  }
  std::string lineOfO2VersionsWithHash;
  while (std::getline(infileO2Versions, lineOfO2Version)) {
    // Extract the tag
    int stripSize = 4;
    std::string tag = lineOfO2Version.substr(lineOfO2Version.find("O2::") + stripSize);
    // Strip a trailing "-1" (some alienv entries append this)
    stripSize = 2;
    if (tag.size() >= stripSize && tag.compare(tag.size() - stripSize, stripSize, "-1") == 0) {
      tag.resize(tag.size() - stripSize);
    }
    LOG(debug) << "Checking tag '" << lineOfO2Version << "' tag (" << tag << ")";
    bool found = false;
    infileO2VersionsWithHash.open(fileContainingCommit);
    while (std::getline(infileO2VersionsWithHash, lineOfO2VersionsWithHash)) {
      // LOG(info) << "Comparing " << lineOfO2Version << " with " << lineOfO2VersionsWithHash;
      if (lineOfO2VersionsWithHash.find(tag) != std::string::npos) {
        LOG(info) << "Tag " << tag << " contains hash " << hash;
        found = true;
        break;
      }
    }
    infileO2VersionsWithHash.close();
    results[tag] = found;
  }
  return results;
}

void populateCCDBWithCommitAvailability(std::map<string, bool> hasHashMap,
                                        const std::string commitHash const std::string ccdbUrl = "http://ccdb-test.cern.ch:8080/")
{
  // First, init the CCDB manager to test if the ccdb is already populated
  o2::ccdb::CcdbApi api;
  api.init(ccdbUrl);
  if (!api.isHostReachable()) {
    LOG(fatal) << "CCDB host " << ccdbUrl << " is not reacheable, cannot go forward";
  }
  for (const auto& entry : hasHashMap) {
    if (!entry.second) { // Version of the code does not have the hash
      continue;
    }
    LOG(info) << "Populating CCDB with information that commit hash " << commitHash << " is contained in software tag " << entry.first;
    std::map<std::string, std::string> metadata;
    metadata["O2Version"] = entry.first;
    const std::string ccdbPath = "O2Version/CommitHash/" + commitHash;
    auto headers = api.retrieveHeaders(ccdbPath, metadata, -1);
    if (headers.size() != 0) {
      LOG(info) << "Entry in CCDB already present for commit hash " << commitHash << ", skipping creation";
      continue;
    }
    LOG(info) << "No entry in CCDB for commit hash " << commitHash << ", creating it";
    std::string s = "available";
    api.storeAsTFileAny<std::string>(&s, ccdbPath, metadata);
  }
}

void testMetadataHelper(std::string aod = "/tmp/AO2D.root")
{
  createO2VersionFile();
  const std::string commitHash = "63bc2e3893851ef0f849bb4c98c65eae1ba21e47";
  const std::map<std::string, bool> hasHashMap = buildMapForCommitHash(commitHash);
  populateCCDBWithCommitAvailability(hasHashMap, commitHash);

  TFile* file = TFile::Open(aod.c_str());
  if (!file || file->IsZombie()) {
    LOG(fatal) << "Could not open file " << aod;
  }
  std::unique_ptr<TFile> currentFile{file};
  std::vector<o2::framework::ConfigParamSpec> specs = readMetadata(currentFile);

  std::vector<std::unique_ptr<o2::framework::ParamRetriever>> retrievers;
  auto paramStore = std::make_unique<o2::framework::ConfigParamStore>(specs, std::move(retrievers));
  paramStore->preload();
  paramStore->activate();
  o2::framework::ConfigParamRegistry paramRegistry(std::move(paramStore));
  o2::framework::ServiceRegistry serviceRegistry;
  o2::framework::ServiceRegistryRef services(serviceRegistry);
  o2::framework::ConfigContext aodCfg(paramRegistry, services, 0, nullptr);
  LOG(info) << "Loaded " << aodCfg.options().specs().size() << " configuration entries from file " << aod;
  aodCfg.options().get<std::string>("aod-metadata-DataType");
  o2::common::core::MetadataHelper metadataInfo;
  metadataInfo.initMetadata(aodCfg);
  metadataInfo.set("O2Version", "epn-20250715"); // Override the O2 version to a known one
  metadataInfo.print();
  LOG(info) << "Metadata label: " << metadataInfo.makeMetadataLabel();

  // Check if the hash is in the software tag
  const std::string v = metadataInfo.getO2Version();
  if (hasHashMap.find(v) == hasHashMap.end()) {
    LOG(fatal) << "Software tag " << v << " not found in available O2 versions";
  }
  if (hasHashMap.at(v)) {
    LOG(info) << "Hash " << commitHash << " is contained in software tag " << v;
  } else {
    LOG(warn) << "Hash " << commitHash << " is NOT contained in software tag " << v;
  }
  if (metadataInfo.isCommitInSoftwareTag(commitHash)) {
    LOG(info) << "MetadataHelper confirms that hash " << commitHash << " is contained in software tag " << v;
  } else {
    LOG(warn) << "MetadataHelper confirms that hash " << commitHash << " is NOT contained in software tag " << v;
  }
}
