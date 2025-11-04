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
// This code is for bc counter.
//    Please write to: daiki.sekihata@cern.ch

#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/Qvectors.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "MathUtils/Utils.h"

#include "TString.h"

#include <algorithm>
#include <optional>
#include <string>
#include <type_traits>
#include <vector>

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

using MyBCs = soa::Join<aod::BCsWithTimestamps, aod::BcSels>;
using MyCollisions = soa::Join<aod::Collisions, aod::EvSels>;
using MyMCCollisions = soa::Join<MyCollisions, aod::McCollisionLabels>;

struct bcCounter {
  // Configurables
  // Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  // Service<o2::ccdb::BasicCCDBManager> ccdb;

  HistogramRegistry fRegistry{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  void init(InitContext&)
  {
    // ccdb->setURL(ccdburl);
    // ccdb->setCaching(true);
    // ccdb->setLocalObjectValidityChecking();
    // ccdb->setFatalWhenNull(false);
    addhistograms();
  }

  ~bcCounter() {}

  void addhistograms()
  {
    // event info

    const int nbin_ev = 3;
    auto hBCCounter = fRegistry.add<TH1>("Data/hBCCounter", "bc counter;;Number of bcs", kTH1D, {{nbin_ev, 0.5, nbin_ev + 0.5}}, false);
    hBCCounter->GetXaxis()->SetBinLabel(1, "all");
    hBCCounter->GetXaxis()->SetBinLabel(2, "FT0AND");
    hBCCounter->GetXaxis()->SetBinLabel(3, "FT0AND && vertex found");
    fRegistry.add("Data/hNcollsPerBC", "Number of rec. collisions per BC", kTH1D, {{21, -0.5, 20.5}}, false);

    fRegistry.addClone("Data/", "MC/");
  }

  SliceCache cache;
  PresliceUnsorted<MyCollisions> preslice_collisions_per_bc = o2::aod::evsel::foundBCId;
  // std::unordered_map<uint64_t, int> map_ncolls_per_bc;

  void processData(MyBCs const& bcs, MyCollisions const& collisions)
  {
    // first count the number of collisions per bc
    for (const auto& bc : bcs) {
      auto collisions_per_bc = collisions.sliceBy(preslice_collisions_per_bc, bc.globalIndex());
      // map_ncolls_per_bc[bc.globalIndex()] = collisions_per_bc.size();
      fRegistry.fill(HIST("Data/hNcollsPerBC"), collisions_per_bc.size());
      // LOGF(info, "bc-loop | bc.globalIndex() = %d , collisions_per_bc.size() = %d", bc.globalIndex(), collisions_per_bc.size());

      fRegistry.fill(HIST("Data/hBCCounter"), 1.0);
      if (bc.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
        fRegistry.fill(HIST("Data/hBCCounter"), 2.0);

        if (collisions_per_bc.size() > 0) { // at least 1 reconstructed vertex exists.
          fRegistry.fill(HIST("Data/hBCCounter"), 3.0);
        }
      }
    } // end of bc loop

    // for (const auto& collision : collisions) {
    //   auto bc = collision.template foundBC_as<MyBCs>();
    //   // LOGF(info, "collision-loop | bc.globalIndex() = %d, ncolls_per_bc = %d", bc.globalIndex(), map_ncolls_per_bc[bc.globalIndex()]);
    // } // end of collision loop

    // map_ncolls_per_bc.clear();
  }
  PROCESS_SWITCH(bcCounter, processData, "process Data", true);

  // void processMC(MyBCs const& bcs, MyMCCollisions const& collisions, aod::McCollisions const& mccollisions)
  // {

  //   // first count the number of collisions per bc
  //   for (const auto& bc : bcs) {
  //     auto collisions_per_bc = collisions.sliceBy(preslice_collisions_per_bc, bc.globalIndex());
  //     // map_ncolls_per_bc[bc.globalIndex()] = collisions_per_bc.size();
  //     fRegistry.fill(HIST("hNcollsPerBC"), collisions_per_bc.size());
  //     // LOGF(info, "bc-loop | bc.globalIndex() = %d , collisions_per_bc.size() = %d", bc.globalIndex(), collisions_per_bc.size());

  //     fRegistry.fill(HIST("hBCCounter"), 1.0);
  //     if (bc.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
  //       fRegistry.fill(HIST("hBCCounter"), 2.0);

  //       if (collisions_per_bc.size() > 0) { // at least 1 reconstructed vertex exists.
  //         fRegistry.fill(HIST("hBCCounter"), 3.0);
  //       }
  //     }
  //   } // end of bc loop
  // }
  // PROCESS_SWITCH(bcCounter, processMC, "process MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<bcCounter>(cfgc, TaskName{"bc-counter"})};
}
