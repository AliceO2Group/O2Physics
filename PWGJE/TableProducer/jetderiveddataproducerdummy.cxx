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

// temporary task to produce HF tables needed when making inclusive derived data - should become obsolete when tables are able to be prouduced based on a configurable
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/runDataProcessing.h"
#include "PWGJE/DataModel/JetReducedData.h"
#include "PWGHF/DataModel/DerivedTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct JetDerivedDataProducerDummyTask {

  Produces<aod::HfD0CollBases> d0CollisionsTable;
  Produces<aod::HfD0Bases> d0sTable;
  Produces<aod::HfD0Pars> d0ParsTable;
  Produces<aod::HfD0ParEs> d0ParExtrasTable;
  Produces<aod::HfD0Sels> d0SelsTable;
  Produces<aod::HfD0Mls> d0MlsTable;
  Produces<aod::HfD0Mcs> d0McsTable;
  Produces<aod::HfD0PBases> d0ParticlesTable;

  Produces<aod::Hf3PCollBases> LcCollisionsTable;
  Produces<aod::Hf3PBases> LcsTable;
  Produces<aod::Hf3PPars> LcParsTable;
  Produces<aod::Hf3PParEs> LcParExtrasTable;
  Produces<aod::Hf3PSels> LcSelsTable;
  Produces<aod::Hf3PMls> LcMlsTable;
  Produces<aod::Hf3PMcs> LcMcsTable;
  Produces<aod::Hf3PPBases> LcParticlesTable;

  void init(InitContext const&)
  {
  }

  void processDummy(aod::JDummys const&)
  {
  }
  PROCESS_SWITCH(JetDerivedDataProducerDummyTask, processDummy, "leaves all tables empty", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<JetDerivedDataProducerDummyTask>(cfgc, TaskName{"jet-deriveddata-producer-dummy"})};
}
