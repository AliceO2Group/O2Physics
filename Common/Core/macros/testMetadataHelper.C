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

#include <Framework/ConfigContext.h>
#include <Framework/ConfigParamRegistry.h>
#include <Framework/ConfigParamStore.h>
#include <Framework/ServiceRegistry.h>
#include <Framework/ServiceRegistryRef.h>

#include <TFile.h>
#include <TMap.h>
#include <TObjString.h>

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

void testMetadataHelper(std::string aod = "/tmp/AO2D.root")
{

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
  metadataInfo.print();
  LOG(info) << "Metadata label: " << metadataInfo.makeMetadataLabel();
}
