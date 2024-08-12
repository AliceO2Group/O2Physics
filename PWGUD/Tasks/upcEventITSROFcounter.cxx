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
#include "ITSMFTBase/DPLAlpideParam.h"
#include "CCDB/BasicCCDBManager.h"

#include "Common/DataModel/EventSelection.h"

#include <TH1D.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using BCsWithRun3Matchings = soa::Join<aod::BCs, aod::Timestamps, aod::Run3MatchedToBCSparse>;
using CCs = soa::Join<aod::Collisions, aod::EvSels>;

struct UpcEventITSROFcounter {
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<int> nTracksForUPCevent{"nTracksForUPCevent", 16, {"Maximum of tracks defining a UPC collision"}};

  void init(InitContext&)
  {

    histos.add("Events/hCountCollisionsExactMatching", ";;Number of collision (-)", HistType::kTH1D, {{11, -0.5, 10.5}});
    histos.add("Events/hCountUPCcollisionsExactMatching", ";;Number of UPC (mult < 17) collision (-)", HistType::kTH1D, {{11, -0.5, 10.5}});
    histos.add("Events/hCountCollisionsInROFborderMatching", ";;Number of collision (-)", HistType::kTH1D, {{11, -0.5, 10.5}});
    histos.add("Events/hCountUPCcollisionsInROFborderMatching", ";;Number of UPC (mult < 17) collision (-)", HistType::kTH1D, {{11, -0.5, 10.5}});

  } // end init

  void process(BCsWithRun3Matchings const& bcs, CCs const& collisions)
  {
    int nAllColls = 0;
    int nUPCcolls = 0;
    uint16_t previousBCinITSROF = 0;
    std::vector<std::pair<uint16_t, uint16_t>> vecITSROFborders;
    bool isFirst = true;
    uint16_t firstBCglobalIndex = 0;
    uint16_t previousBCglobalIndex = 0;

    // extract ITS time frame parameters
    int64_t ts = bcs.iteratorAt(0).timestamp();
    auto alppar = ccdb->getForTimeStamp<o2::itsmft::DPLAlpideParam<0>>("ITS/Config/AlpideParam", ts);

    for (auto bc : bcs) {
      uint64_t globalBC = bc.globalBC();
      uint64_t globalIndex = bc.globalIndex();
      if (isFirst) {
        firstBCglobalIndex = globalIndex;
        isFirst = false;
      }
      uint16_t bcInITSROF = (globalBC + o2::constants::lhc::LHCMaxBunches - alppar->roFrameBiasInBC) % alppar->roFrameLengthInBC;

      if (bcInITSROF - previousBCinITSROF < 0) {
        histos.get<TH1>(HIST("Events/hCountCollisionsExactMatching"))->Fill(nAllColls);
        histos.get<TH1>(HIST("Events/hCountUPCcollisionsExactMatching"))->Fill(nUPCcolls);
        nAllColls = 0;
        nUPCcolls = 0;
        vecITSROFborders.push_back(std::make_pair(firstBCglobalIndex, previousBCglobalIndex));
        firstBCglobalIndex = globalIndex;
      }
      previousBCinITSROF = bcInITSROF;
      previousBCglobalIndex = globalIndex;
      // next is based on exact matching of bc and collision
      for (auto& collision : collisions) {
        if (collision.has_foundBC()) {
          if (collision.foundBCId() == bc.globalIndex()) {
            nAllColls++;
            if (collision.numContrib() < nTracksForUPCevent + 1) {
              nUPCcolls++;
            }
          }
        } else if (collision.bcId() == bc.globalIndex()) {
          nAllColls++;
          if (collision.numContrib() < nTracksForUPCevent + 1) {
            nUPCcolls++;
          }
        }
      } // end loop over collisions
    }   // end loop over bcs

    int arrAllColls[1000] = {0};
    int arrUPCcolls[1000] = {0};

    // next is based on matching of collision bc within ITSROF range in bcs
    for (auto& itsrofBorder : vecITSROFborders) {
      int nAllCollsInROF = 0;
      int nUpcCollsInROF = 0;
      for (auto& collision : collisions) {
        if ((itsrofBorder.first < collision.bcId()) && (collision.bcId() < itsrofBorder.second)) {
          nAllCollsInROF++;
          if (collision.numContrib() < nTracksForUPCevent + 1) {
            nUpcCollsInROF++;
          }
        }
      } // end loop over collisions
      arrAllColls[nAllCollsInROF]++;
      arrUPCcolls[nUpcCollsInROF]++;
    } // end loop over ITSROFs

    for (int ncol = 0; ncol < 12; ncol++) {
      histos.get<TH1>(HIST("Events/hCountCollisionsInROFborderMatching"))->Fill(ncol, arrAllColls[ncol]);
      histos.get<TH1>(HIST("Events/hCountUPCcollisionsInROFborderMatching"))->Fill(ncol, arrUPCcolls[ncol]);
    }
  }

  PROCESS_SWITCH(UpcEventITSROFcounter, process, "Counts number of collisions within ITSROF", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<UpcEventITSROFcounter>(cfgc, TaskName{"upc-event-itsrof-counter"})};
}
