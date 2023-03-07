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

#include "CCDB/CcdbApi.h"
#include "TObjArray.h"
#include "TriggerAliases.h"
#include "TTree.h"
#include "TString.h"
#include <fstream>
#include <map>
#include <string>
#include <vector>
using std::map;
using std::string;

void createDefaultAliases(map<int, TString>& mAliases)
{
  mAliases[kTVXinTRD] = "minbias_TVX";
  mAliases[kTVXinEMC] = "minbias_TVX_L0";
}

void upload_trigger_aliases_run3()
{
  map<int, TString> mAliases;
  createDefaultAliases(mAliases);

  TObjArray* classNames[kNaliases];
  for (auto& al : mAliases) {
    classNames[al.first] = (al.second).Tokenize(",");
  }

  o2::ccdb::CcdbApi ccdb;
  ccdb.init("http://alice-ccdb.cern.ch");
  // ccdb.init("http://ccdb-test.cern.ch:8080");
  // ccdb.truncate("EventSelection/TriggerAliases");

  map<string, string> metadata, metadataRCT, header;

  // read list of runs from text file
  std::ifstream f("runs_run3.txt");
  std::vector<int> runs;
  int r = 0;
  while (f >> r) {
    runs.push_back(r);
  }

  if (1) {
    ULong64_t sor = 1543767116001;
    ULong64_t eor = 1669611662530;
    TriggerAliases* aliases = new TriggerAliases();
    metadata["runNumber"] = "default";
    ccdb.storeAsTFileAny(aliases, "EventSelection/TriggerAliases", metadata, sor, eor);
  }

  for (auto& run : runs) {
    LOGP(info, "run = {}", run);
    if (run < 519903)
      continue; // no CTP info
    if (run == 527349)
      continue; // no CTP info
    if (run == 527963)
      continue; // no CTP info
    if (run == 528537)
      continue; // no CTP info
    if (run == 528543)
      continue; // no CTP info
    // read SOR and EOR timestamps from RCT CCDB
    header = ccdb.retrieveHeaders(Form("RCT/Info/RunInformation/%i", run), metadataRCT, -1);
    ULong64_t sor = atol(header["SOR"].c_str());
    ULong64_t eor = atol(header["EOR"].c_str());
    ULong64_t ts = sor;
    // read CTP config
    metadata["runNumber"] = Form("%d", run);
    auto ctpcfg = ccdb.retrieveFromTFileAny<o2::ctp::CTPConfiguration>("CTP/Config/Config", metadata, ts);
    if (!ctpcfg)
      continue;

    if (run == 529414) { // adding tolerance to sor for this run
      sor = 1668809980000;
    }

    std::vector<o2::ctp::CTPClass> classes = ctpcfg->getCTPClasses();
    // ctpcfg->printConfigString();
    // create trigger aliases
    TriggerAliases* aliases = new TriggerAliases();
    for (auto& al : mAliases) {
      int aliasId = al.first;
      for (const auto& cl : classes) {
        int classId = cl.getIndex();
        LOGP(debug, "class index = {}, name = {}, cluster = {}", classId, cl.name, cl.cluster->name);
        if (cl.name == al.second) {
          if (aliasId == kTVXinTRD && cl.cluster->name != "trd") { // workaround for configs with ambiguous class names
            continue;
          }
          if (aliasId == kTVXinEMC && cl.cluster->name != "emc") { // workaround for configs with ambiguous class names
            continue;
          }
          aliases->AddClassIdToAlias(aliasId, classId);
        }
      }
    }
    aliases->Print();
    ccdb.storeAsTFileAny(aliases, "EventSelection/TriggerAliases", metadata, sor, eor + 10000); // adding tolerance of 10s to eor
  }
}
