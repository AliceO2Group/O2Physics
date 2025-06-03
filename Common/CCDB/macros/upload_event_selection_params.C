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
#include "TString.h"
#include <map>
#include <string>
#include "EventSelectionParams.h"
using std::map;
using std::string;

void upload_event_selection_params()
{
  o2::ccdb::CcdbApi ccdb;

  // ccdb.init("http://ccdb-test.cern.ch:8080");
  // ccdb.truncate("EventSelection/EventSelectionParams");
  ccdb.init("https://alice-ccdb.cern.ch");

  const int nPeriodsMax = 100;
  EventSelectionParams* par[nPeriodsMax];
  string period[nPeriodsMax];
  int runFirst[nPeriodsMax];
  int runLast[nPeriodsMax];
  bool isNew[nPeriodsMax] = {0};

  int n = 0;
  period[n] = "pp2010";
  par[n] = new EventSelectionParams(0, 1);
  runFirst[n] = 114786;
  runLast[n] = 136377;

  n++;
  period[n] = "PbPb_10h";
  par[n] = new EventSelectionParams(3, 1);
  runFirst[n] = 137135;
  runLast[n] = 139517;

  n++;
  period[n] = "pp2011";
  par[n] = new EventSelectionParams(0, 1);
  runFirst[n] = 144871;
  runLast[n] = 159582;

  n++;
  period[n] = "PbPb_11h";
  par[n] = new EventSelectionParams(3, 1);
  runFirst[n] = 167915;
  runLast[n] = 170593;

  n++;
  period[n] = "pp2012";
  par[n] = new EventSelectionParams(0, 2);
  runFirst[n] = 176661;
  runLast[n] = 193752;

  n++;
  period[n] = "pPb2013";
  par[n] = new EventSelectionParams(1, 1);
  runFirst[n] = 195344;
  runLast[n] = 196310;

  n++;
  period[n] = "Pbp2013";
  par[n] = new EventSelectionParams(2, 1);
  runFirst[n] = 196433;
  runLast[n] = 197388;

  n++;
  period[n] = "defaultRun2";
  par[n] = new EventSelectionParams(0);
  runFirst[n] = 209122;
  runLast[n] = 297624;

  n++;
  period[n] = "lhc15f_isolated_bunches";
  par[n] = new EventSelectionParams(0);
  par[n]->DisableOutOfBunchPileupCuts();
  runFirst[n] = 225000;
  runLast[n] = 225719;

  n++;
  period[n] = "lhc15f_isolated_bunches2";
  par[n] = new EventSelectionParams(0);
  par[n]->DisableOutOfBunchPileupCuts();
  runFirst[n] = 226062;
  runLast[n] = 226500;

  n++;
  period[n] = "lhc15f_50ns_trains1";
  par[n] = new EventSelectionParams(0);
  runFirst[n] = 225753;
  runLast[n] = 225768;
  par[n]->SetOnVsOfParams(-372.579114, 9.415265, -6.65857, 0.546801);

  n++;
  period[n] = "lhc15f_50ns_trains2_missing_V0C3";
  par[n] = new EventSelectionParams(0);
  par[n]->SetOnVsOfParams(-372.579114, 9.415265, -6.65857, 0.546801);
  par[n]->fV0CasymA = 0;
  par[n]->fV0CasymB = 0;
  runFirst[n] = 226530;
  runLast[n] = 226606;

  n++;
  period[n] = "lhc15h";
  par[n] = new EventSelectionParams(0);
  par[n]->SetOnVsOfParams(-245.12, 6.86754, -6.65857, 0.546801);
  runFirst[n] = 232914;
  runLast[n] = 233859;

  n++;
  period[n] = "lhc15h_isolated_bunches";
  par[n] = new EventSelectionParams(0);
  par[n]->DisableOutOfBunchPileupCuts();
  runFirst[n] = 233912;
  runLast[n] = 234050;

  n++;
  period[n] = "lhc15i";
  par[n] = new EventSelectionParams(0);
  par[n]->SetOnVsOfParams(-223.155660, 7.117266, -6.218793, 0.543201);
  runFirst[n] = 235196;
  runLast[n] = 236866;

  n++;
  period[n] = "lhc15j";
  par[n] = new EventSelectionParams(0);
  par[n]->SetOnVsOfParams(-222.631866, 7.431432, -6.610850, 0.587165);
  runFirst[n] = 236892;
  runLast[n] = 238621;

  n++;
  period[n] = "lhc15l";
  par[n] = new EventSelectionParams(0);
  par[n]->SetOnVsOfParams(-198.639921, 7.454714, -5.018572, 0.585245);
  runFirst[n] = 239188;
  runLast[n] = 241544;

  n++;
  period[n] = "lhc15n";
  par[n] = new EventSelectionParams(0);
  par[n]->SetOnVsOfParams(-336.279729, 10.694535, -4.144493, 0.851104);
  runFirst[n] = 244340;
  runLast[n] = 244628;

  n++;
  period[n] = "lhc15o1";
  par[n] = new EventSelectionParams(3);
  runFirst[n] = 244824;
  runLast[n] = 245705;

  n++;
  period[n] = "lhc15o2";
  par[n] = new EventSelectionParams(3);
  runFirst[n] = 245829;
  runLast[n] = 246994;

  n++;
  period[n] = "lhc15o_common_zna_tdc";
  par[n] = new EventSelectionParams(3);
  par[n]->fZNSumMean = -123.1;
  par[n]->fZNDifMean = 123.1;
  runFirst[n] = 245729;
  runLast[n] = 245793;

  n++;
  period[n] = "lhc16do";
  par[n] = new EventSelectionParams(0);
  par[n]->SetOnVsOfParams(-65.42, 7.43, -5.62, 0.85);
  runFirst[n] = 252235;
  runLast[n] = 264035;

  n++;
  period[n] = "lhc16r_pPb";
  par[n] = new EventSelectionParams(1);
  runFirst[n] = 265304;
  runLast[n] = 266318;

  n++;
  period[n] = "lhc16s_Pbp";
  par[n] = new EventSelectionParams(2);
  runFirst[n] = 266405;
  runLast[n] = 267131;

  n++;
  period[n] = "lhc16t_pPb";
  par[n] = new EventSelectionParams(1);
  runFirst[n] = 267161;
  runLast[n] = 267166;

  n++;
  period[n] = "pp2017";
  par[n] = new EventSelectionParams(0);
  par[n]->SetOnVsOfParams(-70, 5.2, -3.0, 0.76);
  runFirst[n] = 270531;
  runLast[n] = 280140;

  n++;
  period[n] = "lhc17n_XeXe";
  par[n] = new EventSelectionParams(3);
  runFirst[n] = 280234;
  runLast[n] = 280235;

  n++;
  period[n] = "pp2017pqr";
  par[n] = new EventSelectionParams(0);
  par[n]->SetOnVsOfParams(-70, 5.2, -3.0, 0.76);
  runFirst[n] = 282008;
  runLast[n] = 282704;

  n++;
  period[n] = "lhc18b_isolated_bunches";
  par[n] = new EventSelectionParams(0);
  par[n]->DisableOutOfBunchPileupCuts();
  runFirst[n] = 284706;
  runLast[n] = 285015;

  n++;
  period[n] = "lhc18b_trains";
  par[n] = new EventSelectionParams(0);
  par[n]->SetOnVsOfParams(-65., 4.3, -5.62, 0.85);
  runFirst[n] = 285064;
  runLast[n] = 285203;

  n++;
  period[n] = "lhc18qr_PbPb";
  par[n] = new EventSelectionParams(3);
  runFirst[n] = 295581;
  runLast[n] = 297624;

  n++;
  map<string, string> metadata;
  for (int i = 0; i < n; i++) {
    printf("%s ", period[i].c_str());
    if (!isNew[i]) {
      printf(" .... is not new, skipping\n");
      continue;
    }
    auto sor = o2::ccdb::BasicCCDBManager::getRunDuration(ccdb, runFirst[i]).first;
    auto eor = o2::ccdb::BasicCCDBManager::getRunDuration(ccdb, runLast[i]).second;

    printf("sor=%llu eor=%llu\n", sor, eor);
    metadata["period"] = period[i];
    metadata["run_first"] = Form("%d", runFirst[i]);
    metadata["run_last"] = Form("%d", runLast[i]);
    ccdb.storeAsTFileAny(par[i], "EventSelection/EventSelectionParams", metadata, sor, eor);
  }

  if (1) { // Default Run 3 object
    ULong64_t sorRun3 = 1543767116001;
    ULong64_t eorRun3 = 1893445200000;
    metadata["period"] = "Default Run 3";
    n++;
    par[n] = new EventSelectionParams();
    par[n]->fV0ABBlower = -3.0;  // ns
    par[n]->fV0ABBupper = +2.0;  // ns
    par[n]->fV0ABGlower = 2.0;   // ns
    par[n]->fV0ABGupper = 5.0;   // ns
    par[n]->fFDABBlower = -3.0;  // ns
    par[n]->fFDABBupper = +3.0;  // ns
    par[n]->fFDABGlower = 10.0;  // ns
    par[n]->fFDABGupper = 13.0;  // ns
    par[n]->fFDCBBlower = -3.0;  // ns
    par[n]->fFDCBBupper = +3.0;  // ns
    par[n]->fFDCBGlower = -10.0; // ns
    par[n]->fFDCBGupper = -3.0;  // ns
    par[n]->fT0ABBlower = -1.0;  // ns
    par[n]->fT0ABBupper = +1.0;  // ns
    par[n]->fT0CBBlower = -1.0;  // ns
    par[n]->fT0CBBupper = +1.0;  // ns

    ccdb.storeAsTFileAny(par[n], "EventSelection/EventSelectionParams", metadata, sorRun3, eorRun3);
  }
}
