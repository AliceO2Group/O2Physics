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
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/LFStrangenessPIDTables.h"

using namespace o2;
using namespace o2::framework;

// Converts V0 version 001 to 002
struct strangeDerivedConverter {
  Produces<aod::StraRawCents_001> straRawCents_001;
  Produces<aod::StraRawCents_003> straRawCents_003;
  Produces<aod::V0TOFs> v0tofs;
  Produces<aod::CascTOFs> casctofs;

  void processStraRawCents000to001(aod::StraRawCents_000 const& straRawCents_000)
  {
    for (auto& values : straRawCents_000) {
      straRawCents_001(values.multFT0A(), values.multFT0C(), values.multFV0A(), values.multNTracksPVeta1(), 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
    }
  }
  void processStraRawCents002to003(aod::StraRawCents_002 const& straRawCents_002)
  {
    for (auto& values : straRawCents_002) {
      straRawCents_003(values.multFT0A(),
                       values.multFT0C(),
                       values.multFT0A(),
                       values.multNTracksPVeta1(),
                       0, 0,
                       values.multNTracksITSTPC(),
                       values.multAllTracksTPCOnly(),
                       values.multAllTracksITSTPC(),
                       values.multZNA(),
                       values.multZNC(),
                       values.multZEM1(),
                       values.multZEM2(),
                       values.multZPA(),
                       values.multZPC());
    }
  }

  void processGenerateV0TOFs(soa::Join<aod::V0Cores, aod::V0Extras>::iterator const& v0, soa::Join<aod::DauTrackExtras, aod::DauTrackTOFPIDs> const&)
  {
    // de-reference interlink and populate V0-joinable table (legacy storage)
    auto pTra = v0.posTrackExtra_as<soa::Join<aod::DauTrackExtras, aod::DauTrackTOFPIDs>>();
    auto nTra = v0.negTrackExtra_as<soa::Join<aod::DauTrackExtras, aod::DauTrackTOFPIDs>>();
    v0tofs(pTra.length(), nTra.length(),
           pTra.tofSignal(), nTra.tofSignal(),
           pTra.tofEvTime(), nTra.tofEvTime());
  }

  void processGenerateCascTOFs(soa::Join<aod::CascCores, aod::CascExtras>::iterator const& casc, soa::Join<aod::DauTrackExtras, aod::DauTrackTOFPIDs> const&)
  {
    // de-reference interlink and populate V0-joinable table (legacy storage)
    auto pTra = casc.posTrackExtra_as<soa::Join<aod::DauTrackExtras, aod::DauTrackTOFPIDs>>();
    auto nTra = casc.negTrackExtra_as<soa::Join<aod::DauTrackExtras, aod::DauTrackTOFPIDs>>();
    auto bTra = casc.bachTrackExtra_as<soa::Join<aod::DauTrackExtras, aod::DauTrackTOFPIDs>>();
    casctofs(pTra.length(), nTra.length(), bTra.length(),
             pTra.tofSignal(), nTra.tofSignal(), bTra.tofSignal(),
             pTra.tofEvTime(), nTra.tofEvTime(), bTra.tofEvTime());
  }

  PROCESS_SWITCH(strangeDerivedConverter, processStraRawCents000to001, "from StraRawCents 000 to 001", false);
  PROCESS_SWITCH(strangeDerivedConverter, processStraRawCents002to003, "from StraRawCents 002 to 003", false);
  PROCESS_SWITCH(strangeDerivedConverter, processGenerateV0TOFs, "from DauTrackTOFPIDs to V0TOFs", false);
  PROCESS_SWITCH(strangeDerivedConverter, processGenerateCascTOFs, "from DauTrackTOFPIDs to CascTOFs", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<strangeDerivedConverter>(cfgc)};
}
