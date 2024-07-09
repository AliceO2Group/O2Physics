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

// temporary task to produce HF and DQ tables needed when making D0 jet derived data - should become obsolete when tables are able to be prouduced based on a configurable
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/runDataProcessing.h"
#include "PWGJE/DataModel/JetReducedData.h"
#include "PWGHF/DataModel/DerivedTables.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct JetDerivedDataProducerDummyD0Task {

  Produces<aod::Hf3PCollBases> lcCollisionsTable;
  Produces<aod::Hf3PMcRCollIds> lcCollisionsMatchingTable;
  Produces<aod::Hf3PBases> lcsTable;
  Produces<aod::Hf3PPars> lcParsTable;
  Produces<aod::Hf3PParEs> lcParExtrasTable;
  Produces<aod::Hf3PSels> lcSelsTable;
  Produces<aod::Hf3PMls> lcMlsTable;
  Produces<aod::Hf3PMcs> lcMcsTable;
  Produces<aod::Hf3PMcCollBases> lcMcCollisionsTable;
  Produces<aod::Hf3PPBases> lcParticlesTable;

  Produces<aod::ReducedEvents> dielectronCollisionsTable;
  Produces<aod::Dielectrons> dielectronTable;

  void init(InitContext const&)
  {
  }

  void processDummy(aod::JDummys const&)
  {
  }
  PROCESS_SWITCH(JetDerivedDataProducerDummyD0Task, processDummy, "leaves all tables empty", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<JetDerivedDataProducerDummyD0Task>(cfgc, TaskName{"jet-deriveddata-producer-dummy-d0"})};
}
