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
  mAliases[kINT7] = "CINT7-B-NOPF-CENT,CV0L7-B-NOPF-CENT,CINT7-B-NOPF-CENTNOTRD,CINT7ZAC-B-NOPF-CENTNOPMD";
  mAliases[kEMC7] = "CEMC7-B-NOPF-CENTNOPMD,CDMC7-B-NOPF-CENTNOPMD";
  mAliases[kINT7inMUON] = "CINT7-B-NOPF-MUFAST";
  mAliases[kMuonSingleLowPt7] = "CMSL7-B-NOPF-MUFAST";
  mAliases[kMuonUnlikeLowPt7] = "CMUL7-B-NOPF-MUFAST";
  mAliases[kMuonLikeLowPt7] = "CMLL7-B-NOPF-MUFAST";
  mAliases[kMuonSingleHighPt7] = "CMSH7-B-NOPF-MUFAST";
  mAliases[kCUP8] = "CCUP8-B-NOPF-CENTNOTRD";
  mAliases[kCUP9] = "CCUP9-B-NOPF-CENTNOTRD";
  mAliases[kMUP10] = "CMUP10-B-NOPF-MUFAST";
  mAliases[kMUP11] = "CMUP11-B-NOPF-MUFAST";
}

void upload_trigger_aliases()
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
  std::ifstream f("runs.txt");
  std::vector<int> runs;
  int r = 0;
  while (f >> r) {
    runs.push_back(r);
  }

  for (auto& run : runs) {
    LOGP(info, "run = {}", run);
    // read list of trigger classes
    map<string, int>* pClassNameToIndexMap = ccdb.retrieveFromTFileAny<map<string, int>>("CTP/ClassNameToIndexMap", metadata, run);
    // fill map of alias-to-class-indices
    TriggerAliases* aliases = new TriggerAliases();
    for (auto& al : mAliases) {
      int aliasId = al.first;
      LOGP(debug, " aliadId = {}", aliasId);
      for (const auto& className : *(classNames[aliasId])) {
        string sname = className->GetName();
        int classId = (*pClassNameToIndexMap)[sname] - 1;
        if (classId < 0) { // class doesn't exist in this run
          continue;
        }
        LOGP(debug, " className = {}", sname.data());
        aliases->AddClassIdToAlias(aliasId, classId);
      }
    }

    // read SOR and EOR timestamps from RCT CCDB
    header = ccdb.retrieveHeaders(Form("RCT/Info/RunInformation/%i", run), metadataRCT, -1);
    ULong64_t sor = atol(header["SOR"].c_str());
    ULong64_t eor = atol(header["EOR"].c_str());
    // add safety margins to avoid edge effects due to SOR/EOR time differences in DCS and CTP
    sor -= 60000;
    eor += 300000;
    ccdb.storeAsTFileAny(aliases, "EventSelection/TriggerAliases", metadata, sor, eor);
  }
}
