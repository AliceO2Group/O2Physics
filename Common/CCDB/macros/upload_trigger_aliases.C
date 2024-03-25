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
#include "CCDB/BasicCCDBManager.h"
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
  mAliases[kINT7] = "CINT7-B-NOPF-CENT,CINT7-B-NOPF-FAST,CV0L7-B-NOPF-CENT,CINT7-B-NOPF-CENTNOTRD,CINT7ZAC-B-NOPF-CENTNOPMD,CINT7-B-NOPF-ALLNOTRD,CINT7-I-NOPF-ALLNOTRD,CINT7-S-NOPF-ALLNOTRD,CMBAC-B-NOPF-ALL,CMBACS2-B-NOPF-ALLNOTRD,CMBACS2-B-NOPF-ALL";
  mAliases[kEMC7] = "CEMC7-B-NOPF-CENTNOPMD,CEMC7-B-NOPF-CENT,CEMC7-B-NOPF-CENTNOTRD,CEMC7-B-NOPF-ALLNOTRD,CEMC7-S-NOPF-ALLNOTRD";
  mAliases[kINT7inMUON] = "CINT7-B-NOPF-MUFAST";
  mAliases[kMuonSingleLowPt7] = "CMSL7-B-NOPF-MUFAST,CMSL7-B-NOPF-MUON,CMSL7-S-NOPF-MUON,CMSL7-SC-NOPF-MUON,CPBI1MSL-B-NOPF-MUON,CMUS7-B-NOPF-MUON";
  mAliases[kMuonUnlikeLowPt7] = "CMUL7-B-NOPF-MUFAST,CMUL7-B-NOPF-MUON,CMUL7-S-NOPF-MUON,CMUL7-S-NOPF-ALLNOTRD,CPBI1MUL-B-NOPF-MUON,CMUU7-B-NOPF-MUON,CMUU7-B-NOPF-ALLNOTRD";
  mAliases[kMuonLikeLowPt7] = "CMLL7-B-NOPF-MUFAST,CMLL7-B-NOPF-MUON,CMLL7-S-NOPF-MUON,CPBI1MLL-B-NOPF-MUON";
  mAliases[kMuonSingleHighPt7] = "CMSH7-B-NOPF-MUFAST,CMSH7-B-NOPF-MUON,CMSH7-S-NOPF-MUON,CMSH7-S-NOPF-ALLNOTRD,CPBI1MSH-B-NOPF-MUON,CMUSH7-B-NOPF-MUON";
  mAliases[kCUP8] = "CCUP8-B-NOPF-CENTNOTRD";
  mAliases[kCUP9] = "CCUP9-B-NOPF-CENTNOTRD";
  mAliases[kMUP10] = "CMUP10-B-NOPF-MUFAST";
  mAliases[kMUP11] = "CMUP11-B-NOPF-MUFAST";
  mAliases[kINT1] = "CINT1B-ABCE-NOPF-ALL,CINT1-B-NOPF-ALLNOTRD,CINT1-B-NOPF-ALLNOTRD";
  mAliases[kUnbiased] = "CBEAMB-ABCE-NOPF-ALL,CBEAMB-B-NOPF-ALLNOTRD,CTRUE-B-NOPF-ALLNOTRD,CTRUE-S-NOPF-ALLNOTRD,CTRUE-B-NOPF-CENTNOTRD";
  mAliases[kDMC7] = "CDMC7-B-NOPF-CENTNOPMD,CDMC7-B-NOPF-CENT,CDMC7-B-NOPF-CENTNOTRD,CDMC7-B-NOPF-ALLNOTRD";
  mAliases[kEG1] = "CINT7EG1-B-NOPF-CENTNOPMD,CEMC7EG1-B-NOPF-CENT,CEMC7EG1-B-NOPF-CENTNOTRD,CEMC7EG1-B-NOPF-CENTNOPMD,CEMC7EG1-B-NOPF-ALLNOTRD,CEMC7EGA-B-NOPF-CENTNOTRD,CEMC7EGA-S-NOPF-CENTNOTRD,CEMC7EGA-S-NOPF-ALLNOTRD,CPBI2EGA-B-NOPF-CENTNOTRD";
  mAliases[kEJ1] = "CINT7EJ1-B-NOPF-CENTNOTRD,CEMC7EJ1-B-NOPF-CENT,CEMC7EJ1-B-NOPF-CENTNOTRD,CEMC7EJ1-B-NOPF-CENTNOPMD,CEMC7EJ1-B-NOPF-ALLNOTRD,CEMC7EJE-B-NOPF-CENTNOTRD,CEMC7EJE-S-NOPF-CENTNOTRD,CEMC7EJE-S-NOPF-ALLNOTRD,CPBI2EJE-B-NOPF-CENTNOTRD";
  mAliases[kEG2] = "CINT7EG2-B-NOPF-CENTNOPMD,CEMC7EG2-B-NOPF-CENT,CEMC7EG2-B-NOPF-CENTNOTRD,CEMC7EG2-B-NOPF-CENTNOPMD,CEMC7EG2-B-NOPF-ALLNOTRD,CEMC7EG2PER-B-NOPF-CENTNOPMD";
  mAliases[kEJ2] = "CINT7EJ2-B-NOPF-CENTNOPMD,CEMC7EJ2-B-NOPF-CENT,CEMC7EJ2-B-NOPF-CENTNOTRD,CEMC7EJ2-B-NOPF-CENTNOPMD,CEMC7EJ2-B-NOPF-ALLNOTRD";
  mAliases[kDG1] = "CINT7DG1-B-NOPF-CENTNOPMD,CDMC7DG1-B-NOPF-CENT,CDMC7DG1-B-NOPF-CENTNOTRD,CDMC7DG1-B-NOPF-CENTNOPMD";
  mAliases[kDJ1] = "CINT7DJ1-B-NOPF-CENTNOPMD,CDMC7DJ1-B-NOPF-CENT,CDMC7DJ1-B-NOPF-CENTNOTRD,CDMC7DJ1-B-NOPF-CENTNOPMD";
  mAliases[kDG2] = "CINT7DG2-B-NOPF-CENTNOPMD,CDMC7DG2-B-NOPF-CENT,CDMC7DG2-B-NOPF-CENTNOTRD,CDMC7DG2-B-NOPF-CENTNOPMD,CDMC7DG2PER-B-NOPF-CENTNOPMD";
  mAliases[kDJ2] = "CINT7DJ2-B-NOPF-CENTNOPMD,CDMC7DJ2-B-NOPF-CENT,CDMC7DJ2-B-NOPF-CENTNOTRD,CDMC7DJ2-B-NOPF-CENTNOPMD";
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

  map<string, string> metadata;

  // read list of runs from text file
  std::ifstream f("runs_run1.txt");
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

    // read SOR and EOR timestamps from RCT CCDB (via utility function)
    auto soreor = o2::ccdb::BasicCCDBManager::getRunDuration(ccdb, run);
    ULong64_t sor = soreor.first;
    ULong64_t eor = soreor.second;

    // add safety margins to avoid edge effects due to SOR/EOR time differences in DCS and CTP
    sor -= 60000;
    eor += 300000;
    ccdb.storeAsTFileAny(aliases, "EventSelection/TriggerAliases", metadata, sor, eor);
  }
}
