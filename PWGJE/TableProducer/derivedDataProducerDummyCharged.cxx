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

// temporary task to produce HF, DQ and EMCal tables needed when making inclusive derived data - should become obsolete when tables are able to be prouduced based on a configurable
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

struct JetDerivedDataProducerDummyChargedTask {

  Produces<aod::HfD0CollBases> d0CollisionsTable;
  Produces<aod::HfD0McRCollIds> d0CollisionsMatchingTable;
  Produces<aod::HfD0Bases> d0sTable;
  Produces<aod::HfD0Pars> d0ParsTable;
  Produces<aod::HfD0ParEs> d0ParExtrasTable;
  Produces<aod::HfD0Sels> d0SelsTable;
  Produces<aod::HfD0Mls> d0MlsTable;
  Produces<aod::HfD0Mcs> d0McsTable;
  Produces<aod::HfD0McCollBases> d0McCollisionsTable;
  Produces<aod::HfD0PBases> d0ParticlesTable;

  Produces<aod::HfLcCollBases> lcCollisionsTable;
  Produces<aod::HfLcMcRCollIds> lcCollisionsMatchingTable;
  Produces<aod::HfLcBases> lcsTable;
  Produces<aod::HfLcPars> lcParsTable;
  Produces<aod::HfLcParEs> lcParExtrasTable;
  Produces<aod::HfLcSels> lcSelsTable;
  Produces<aod::HfLcMls> lcMlsTable;
  Produces<aod::HfLcMcs> lcMcsTable;
  Produces<aod::HfLcMcCollBases> lcMcCollisionsTable;
  Produces<aod::HfLcPBases> lcParticlesTable;

  Produces<aod::HfBplusCollBases> bplusCollisionsTable;
  Produces<aod::HfBplusMcRCollIds> bplusCollisionsMatchingTable;
  Produces<aod::HfBplusBases> bplussTable;
  Produces<aod::HfBplusPars> bplusParsTable;
  Produces<aod::HfBplusParEs> bplusParExtrasTable;
  Produces<aod::HfBplusParD0s> bplusParD0sTable;
  Produces<aod::HfBplusSels> bplusSelsTable;
  Produces<aod::HfBplusMls> bplusMlsTable;
  Produces<aod::HfBplusMlD0s> bplusMlD0sTable;
  Produces<aod::HfBplusMcs> bplusMcsTable;
  Produces<aod::HfBplusMcCollBases> bplusMcCollisionsTable;
  Produces<aod::HfBplusPBases> bplusParticlesTable;

  Produces<aod::ReducedEvents> dielectronCollisionsTable;
  Produces<aod::Dielectrons> dielectronTable;

  Produces<aod::JClustersCorrectedEnergies> jClustersCorrectedEnergiesTable;

  void init(InitContext const&)
  {
  }

  void processDummy(aod::JDummys const&)
  {
  }
  PROCESS_SWITCH(JetDerivedDataProducerDummyChargedTask, processDummy, "leaves all tables empty", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<JetDerivedDataProducerDummyChargedTask>(cfgc, TaskName{"jet-deriveddata-producer-dummy-charged"})};
}
