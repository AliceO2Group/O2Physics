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
//
/// \file upcTestRctTables.cxx
/// \brief Tests Rct Tables in UD tabl;es
///
/// \author Roman Lavicka <roman.lavicka@cern.ch>, Austrian Academy of Sciences & SMI
/// \since  27.06.2025
//

// O2 headers
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"

// O2Physics headers
#include "PWGUD/Core/SGSelector.h"
#include "PWGUD/Core/UPCTauCentralBarrelHelperRL.h"
#include "PWGUD/DataModel/UDTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct UpcTestRctTables {

  // Global varialbes
  SGSelector sgSelector;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // declare configurables
  Configurable<bool> verboseInfo{"verboseInfo", false, {"Print general info to terminal; default it false."}};

  using FullSGUDCollisions = soa::Join<aod::UDCollisions, aod::UDCollisionsSels, aod::UDCollisionSelExtras, aod::SGCollisions>;
  using FullSGUDCollision = FullSGUDCollisions::iterator;

  // init
  void init(InitContext&)
  {
    if (verboseInfo)
      printMediumMessage("INIT METHOD");

    histos.add("OutputTable/hRCTflags", ";RCTflag (-);Number of passed collision (-)", HistType::kTH1D, {{5, -0.5, 4.5}});

  } // end init

  void processDataSG(FullSGUDCollision const& reconstructedCollision)
  {

    histos.get<TH1>(HIST("OutputTable/hRCTflags"))->Fill(0);

    if (sgSelector.isCBTOk(reconstructedCollision))
      histos.get<TH1>(HIST("OutputTable/hRCTflags"))->Fill(1);

    if (sgSelector.isCBTZdcOk(reconstructedCollision))
      histos.get<TH1>(HIST("OutputTable/hRCTflags"))->Fill(2);

    if (sgSelector.isCBTHadronOk(reconstructedCollision))
      histos.get<TH1>(HIST("OutputTable/hRCTflags"))->Fill(3);

    if (sgSelector.isCBTHadronZdcOk(reconstructedCollision))
      histos.get<TH1>(HIST("OutputTable/hRCTflags"))->Fill(4);

  } // end processDataSG

  PROCESS_SWITCH(UpcTestRctTables, processDataSG, "Iterate UD tables with measured data created by SG-Candidate-Producer.", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<UpcTestRctTables>(cfgc)};
}
