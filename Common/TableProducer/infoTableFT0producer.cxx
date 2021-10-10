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
#include "Framework/ConfigParamSpec.h"

using namespace o2;
using namespace o2::framework;

// custom configurable for switching between run2 and run3 selection types
/*void customize(std::vector<ConfigParamSpec>& workflowOptions)
{
  workflowOptions.push_back(ConfigParamSpec{"ccdb-path-ft0", o2::framework::VariantType::String, "http://o2-ccdb.internal/", {"URL of the CCDB database"}});
  }*/

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/EventSelection.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/CCDB/TriggerAliases.h"
#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/LHCConstants.h"
#include "Common/DataModel/T0infoTable.h"
#include "iostream"
#include <string>

using namespace o2::aod;
struct infoTableFT0producer {
  Produces<o2::aod::T0info> table;
  using BCsWithMatchings = soa::Join<aod::BCs, aod::Run3MatchedToBCSparse>;
  Service<o2::ccdb::BasicCCDBManager> ccdb; /// Object manager in CCDB
  o2::ccdb::CcdbApi ccdb_api;               /// API to access CCDB
  Configurable<std::string> url{"ccdb-url", "http://alice-ccdb.cern.ch", "URL of the CCDB database"};
  using CollisionEvSel = soa::Join<aod::Collisions, aod::EvSels>::iterator;
  /*  void init(o2::framework::InitContext&)
  {
    LOGF(info, "Initializing CCDB");
    //   ccdb.setURL(url.value); // Setting URL of CCDB manager from configuration
    //    ccdb_api.init(url.value);
    //    if (!ccdb_api.isHostReachable()) {
    //      LOGF(fatal, "CCDB host %s is n
    ot reacheable, cannot go forward", url.value.data());
    //    }
  }
  */
  void process(BCsWithMatchings const& bcs, soa::Join<aod::Collisions, aod::EvSels> const& collisions, aod::FT0s const& ft0s)
  {
    //   auto caliboffsets = ccdb.get<o2::ft0::GlobalOffsetsCalibrationObject>("FT0/Calibration/GlobalOffsets");
    LOG(INFO) << "infoTableFT0producer";
    float c = 29.9792458; // cm/ns
    for (auto& collision : collisions) {
      int64_t foundFT0 = collision.foundFT0();
      float vertexPV = collision.posZ();
      float vertex_corr = vertexPV / c;
      float t0A, t0C, t0AC, vertexT0;
      t0A = t0C = t0AC = vertexT0 = 32000;
      if (foundFT0 != -1) {
        auto ft0 = ft0s.iteratorAt(foundFT0);
        t0A = ft0.timeA();
        t0C = ft0.timeC();
        int triggersignals = ft0.triggerMask();
        bool ora = (triggersignals & (1 << 0)) != 0;
        bool orc = (triggersignals & (1 << 1)) != 0;
        bool vertexTrigger = (triggersignals & (1 << 2)) != 0;
        LOGF(info, "triggers OrA %i OrC %i vertexTrigger %i", ora, orc, vertexTrigger);
        LOGF(info, " T0A = %f, T0C %f, vertex_corr %f, triggersignals %i", t0A, t0C, vertex_corr, triggersignals);
        if (vertexTrigger) {
          vertexT0 = c * (t0C - t0A) / 2;
          t0AC = (t0A + t0C) / 2;
        }
        if (ora) {
          t0A += vertex_corr;
        }
        if (orc) {
          t0C -= vertex_corr;
        }
      }
      LOGF(info, " T0 collision time = %f, T0A = %f, T0C = %f, T0 vertex = %f PV = %f", t0A, t0C, t0AC, vertexT0, vertexPV);
      table(t0A, t0C, t0AC, vertexT0);
    }
  }
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  LOG(INFO) << "WorkflowSpec defineDataProcessin";
  return WorkflowSpec{adaptAnalysisTask<infoTableFT0producer>(cfgc, TaskName{"ft0-table"})};
}
