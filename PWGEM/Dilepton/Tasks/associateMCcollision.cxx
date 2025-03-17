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
    if (doprocessNcontrib && doprocessNcontrib_Derived) {
      LOGF(fatal, "Please select only 1 process function.");
    }
    addhistograms();
  }

  ~associateMCcollision() {}

  void addhistograms()
  {
    fRegistry.add("hReccollsPerMCcoll", "Rec. colls per MC coll;Rec. colls per MC coll;Number of MC collisions", kTH1D, {{21, -0.5, 20.5}}, false);
  }

  template <typename TMCCollisions, typename TCollisions, typename TPreslice>
  void runMC(TMCCollisions const& mcCollisions, TCollisions const& collisions, TPreslice const& perMCCollision)
  {

    for (auto& mcCollision : mcCollisions) {
      auto rec_colls_per_mccoll = collisions.sliceBy(perMCCollision, mcCollision.globalIndex());
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

  } // end of runMC

  using MyCollisions = soa::Join<aod::Collisions, aod::McCollisionLabels>;
  using MyCollision = MyCollisions::iterator;
  PresliceUnsorted<aod::McCollisionLabels> recColperMcCollision = aod::mccollisionlabel::mcCollisionId;

  using MyEMCollisions = soa::Join<aod::EMEvents, aod::EMMCEventLabels>;
  using MyEMCollision = MyEMCollisions::iterator;
  PresliceUnsorted<aod::EMMCEventLabels> recColperMcCollision_derived = aod::emmceventlabel::emmceventId;

  void processNcontrib_Derived(aod::EMMCEvents const& mcCollisions, MyEMCollisions const& collisions)
  {
    runMC(mcCollisions, collisions, recColperMcCollision_derived);
  }
  PROCESS_SWITCH(associateMCcollision, processNcontrib_Derived, "produce most probable emeventId based on Ncontrib to PV for derived AOD", true);

  void processNcontrib(aod::McCollisions const& mcCollisions, MyCollisions const& collisions)
  {
    runMC(mcCollisions, collisions, recColperMcCollision);
  }
  PROCESS_SWITCH(associateMCcollision, processNcontrib, "produce most probable emeventId based on Ncontrib to PV for original AOD", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<associateMCcollision>(cfgc, TaskName{"associate-mccollision-to-collision"})};
}
