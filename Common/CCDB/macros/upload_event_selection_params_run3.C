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
#include "TString.h"
#include <map>
#include <string>
#include "EventSelectionParams.h"
using std::map;
using std::string;

void fillMapOfCustomOrbitShifts(std::map<int, int>& mapOrbitShift);

void upload_event_selection_params_run3()
{
  o2::ccdb::CcdbApi ccdb;
  ccdb.init("https://alice-ccdb.cern.ch");

  std::map<int, int> mapOrbitShift;
  fillMapOfCustomOrbitShifts(mapOrbitShift);
  map<string, string> metadata, metadataRCT, header;

  EventSelectionParams* params = new EventSelectionParams();
  params->fV0ABBlower = -3.0;  // ns
  params->fV0ABBupper = +2.0;  // ns
  params->fV0ABGlower = 2.0;   // ns
  params->fV0ABGupper = 5.0;   // ns
  params->fFDABBlower = -3.0;  // ns
  params->fFDABBupper = +3.0;  // ns
  params->fFDABGlower = 10.0;  // ns
  params->fFDABGupper = 13.0;  // ns
  params->fFDCBBlower = -3.0;  // ns
  params->fFDCBBupper = +3.0;  // ns
  params->fFDCBGlower = -10.0; // ns
  params->fFDCBGupper = -3.0;  // ns
  params->fT0ABBlower = -1.0;  // ns
  params->fT0ABBupper = +1.0;  // ns
  params->fT0CBBlower = -1.0;  // ns
  params->fT0CBBupper = +1.0;  // ns

  // default object
  bool uploadDefaultParameters = 1;
  if (uploadDefaultParameters) {
    LOGP(info, "Uploading default object");
    ULong64_t sorRun3 = 1543767116001;
    ULong64_t eorRun3 = 1893445200000;
    metadata["runNumber"] = "Default Run 3";
    ccdb.storeAsTFileAny(params, "EventSelection/EventSelectionParams", metadata, sorRun3, eorRun3);
  }

  // fill param objects for runs with custom orbit shifts
  for (auto m : mapOrbitShift) {
    int run = m.first;
    int orbitShift = m.second;
    LOGP(info, "Uploading custom object for run = {}: orbitShift = {}", run, orbitShift);
    params->fTimeFrameOrbitShift = orbitShift;
    header = ccdb.retrieveHeaders(Form("RCT/Info/RunInformation/%i", run), metadataRCT, -1);
    ULong64_t sor = atol(header["SOR"].c_str());
    ULong64_t eor = atol(header["EOR"].c_str());
    metadata["runNumber"] = Form("%d", run);
    ccdb.storeAsTFileAny(params, "EventSelection/EventSelectionParams", metadata, sor, eor + 10000); // adding tolerance of 10s to eor
  }
}

// map of custom orbit shifts of the first orbit in TF wrt (SOR-OrbitReset)%TFDuration (in orbits)
void fillMapOfCustomOrbitShifts(std::map<int, int>& mapOrbitShift)
{
  mapOrbitShift[517619] = 109;
  mapOrbitShift[517620] = 109;
  mapOrbitShift[517623] = 109;
  mapOrbitShift[517677] = 127;
  mapOrbitShift[517678] = 127;
  mapOrbitShift[517679] = 127;
  mapOrbitShift[517685] = 127;
  mapOrbitShift[517690] = 127;
  mapOrbitShift[517693] = 127;
  mapOrbitShift[517737] = 127;
  mapOrbitShift[517748] = 127;
  mapOrbitShift[517751] = 127;
  mapOrbitShift[517753] = 127;
  mapOrbitShift[517758] = 127;
  mapOrbitShift[517767] = 127;
  mapOrbitShift[518541] = 40;
  mapOrbitShift[518543] = 92;
  mapOrbitShift[518546] = 124;
  mapOrbitShift[518547] = 47;
  mapOrbitShift[519041] = 59;
  mapOrbitShift[519043] = 59;
  mapOrbitShift[519045] = 59;
  mapOrbitShift[519497] = 86;
  mapOrbitShift[519498] = 86;
  mapOrbitShift[519499] = 86;
  mapOrbitShift[519502] = 86;
  mapOrbitShift[519503] = 86;
  mapOrbitShift[519504] = 86;
  mapOrbitShift[519506] = 86;
  mapOrbitShift[519507] = 86;
  mapOrbitShift[519903] = 62;
  mapOrbitShift[519904] = 62;
  mapOrbitShift[519905] = 62;
  mapOrbitShift[519906] = 62;
  mapOrbitShift[520259] = 76;
  mapOrbitShift[520294] = 76;
  mapOrbitShift[520471] = 46;
  mapOrbitShift[520472] = 46;
  mapOrbitShift[520473] = 46;
  mapOrbitShift[523142] = 127;
  mapOrbitShift[523148] = 127;
  mapOrbitShift[523182] = 127;
  mapOrbitShift[523186] = 127;
  mapOrbitShift[523298] = 28;
  mapOrbitShift[523306] = 28;
  mapOrbitShift[523308] = 28;
  mapOrbitShift[523309] = 28;
  mapOrbitShift[523397] = 110;
  mapOrbitShift[523399] = 110;
  mapOrbitShift[523401] = 110;
  mapOrbitShift[523441] = 117;
  mapOrbitShift[523541] = 103;
  mapOrbitShift[523559] = 103;
  mapOrbitShift[523669] = 39;
  mapOrbitShift[523671] = 39;
  mapOrbitShift[523677] = 39;
  mapOrbitShift[523728] = 113;
  mapOrbitShift[523731] = 113;
  mapOrbitShift[523779] = 41;
  mapOrbitShift[523783] = 41;
  mapOrbitShift[523786] = 41;
  mapOrbitShift[523788] = 41;
  mapOrbitShift[523789] = 41;
  mapOrbitShift[523792] = 41;
  mapOrbitShift[523797] = 41;
  mapOrbitShift[523821] = 36;
  mapOrbitShift[523897] = 38;
}