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
// ========================
//
// This code produces a table with an index between mc collision and rec. collision.
//    Please write to: daiki.sekihata@cern.ch

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"

#include "PWGEM/Dilepton/DataModel/dileptonTables.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

struct associateMCcollision {
  Produces<aod::MostProbableEMEventIdsInMC> mpemeventIds;

  HistogramRegistry fRegistry{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  void init(InitContext&)
  {
    addhistograms();
  }

  ~associateMCcollision() {}

  void addhistograms()
  {
    fRegistry.add("hReccollsPerMCcoll", "Rec. colls per MC coll;Rec. colls per MC coll;Number of MC collisions", kTH1F, {{21, -0.5, 20.5}}, false);
  }

  using MyCollisions = soa::Join<aod::EMEvents, aod::EMMCEventLabels>;
  using MyCollision = MyCollisions::iterator;
  PresliceUnsorted<MyCollisions> recColperMcCollision = aod::emmceventlabel::emmceventId;

  void processNcontrib(aod::EMMCEvents const& mccollisions, MyCollisions const& collisions)
  {
    for (auto& mccollision : mccollisions) {
      auto rec_colls_per_mccoll = collisions.sliceBy(recColperMcCollision, mccollision.globalIndex());
      fRegistry.fill(HIST("hReccollsPerMCcoll"), rec_colls_per_mccoll.size());
      uint32_t maxNumContrib = 0;
      int rec_col_globalIndex = -999;
      for (auto& rec_col : rec_colls_per_mccoll) {
        if (rec_col.numContrib() > maxNumContrib) {
          rec_col_globalIndex = rec_col.globalIndex();
          maxNumContrib = rec_col.numContrib(); // assign mc collision to collision where the number of contibutors is the lagest. LF/MM recommendation
        }
      }
      // LOGF(info, "rec_col_globalIndex = %d", rec_col_globalIndex);
      mpemeventIds(rec_col_globalIndex);
    } // end of mc collision
  } // end of process
  PROCESS_SWITCH(associateMCcollision, processNcontrib, "produce most probable emeventId based on Ncontrib to PV", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<associateMCcollision>(cfgc, TaskName{"associate-mccollision-to-collision"})};
}
