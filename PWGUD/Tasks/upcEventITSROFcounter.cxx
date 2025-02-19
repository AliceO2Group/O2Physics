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

/// \file upcEventITSROFcounter.cxx
/// \brief Personal task to analyze tau events from UPC collisions
///
/// \author Roman Lavicka <roman.lavicka@cern.ch>, Austrian Academy of Sciences & SMI
/// \since  09.08.2024

#include <utility>
#include <vector>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "ITSMFTBase/DPLAlpideParam.h"
#include "CCDB/BasicCCDBManager.h"
#include "ReconstructionDataFormats/Vertex.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/DataModel/EventSelection.h"

#include "PWGUD/DataModel/UDTables.h"
#include "PWGUD/Core/SGSelector.h"

#include <TH1D.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::dataformats;

using BCsWithRun3Matchings = soa::Join<aod::BCs, aod::Timestamps, aod::Run3MatchedToBCSparse>;
using CCs = soa::Join<aod::Collisions, aod::EvSels>;
using FullSGUDCollision = soa::Join<aod::UDCollisions, aod::UDCollisionsSels, aod::SGCollisions, aod::UDZdcsReduced>::iterator;

struct UpcEventITSROFcounter {
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  SGSelector sgSelector;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<int> nTracksForUPCevent{"nTracksForUPCevent", 16, {"Maximum of tracks defining a UPC collision"}};

  Configurable<bool> useTrueGap{"useTrueGap", true, {"Calculate gapSide for a given FV0/FT0/ZDC thresholds"}};
  Configurable<float> cutMyGapSideFV0{"cutMyGapSideFV0", -1, "FV0A threshold for SG selector"};
  Configurable<float> cutMyGapSideFT0A{"cutMyGapSideFT0A", 150., "FT0A threshold for SG selector"};
  Configurable<float> cutMyGapSideFT0C{"cutMyGapSideFT0C", 50., "FT0C threshold for SG selector"};
  Configurable<float> cutMyGapSideZDC{"cutMyGapSideZDC", 10., "ZDC threshold for SG selector"};
  ConfigurableAxis axisRunNumbers{"axisRunNumbers", {1400, 544000.5, 545400.5}, "Range of run numbers"};

  void init(InitContext&)
  {

    histos.add("Events/hCountCollisionsExactMatching", ";;Number of collision (-)", HistType::kTH1D, {{11, -0.5, 10.5}});
    histos.add("Events/hCountUPCcollisionsExactMatching", ";;Number of UPC (mult < 17) collision (-)", HistType::kTH1D, {{11, -0.5, 10.5}});
    histos.add("Events/hCountCollisionsInROFborderMatching", ";;Number of collision (-)", HistType::kTH1D, {{11, -0.5, 10.5}});
    histos.add("Events/hCountUPCcollisionsInROFborderMatching", ";;Number of UPC (mult < 17) collision (-)", HistType::kTH1D, {{11, -0.5, 10.5}});

    histos.add("Events/hPVcontribsVsCollisionsPerITSROFstd", "Collisions reconstructed with standard mode;Number of vertex contributors (-); Number of collisions in one ITSROF (-)", HistType::kTH2D, {{101, -0.5, 100.5}, {11, -0.5, 10.5}});
    histos.add("Events/hPVcontribsVsCollisionsPerITSROFupc", "Collisions reconstructed with upc mode;Number of vertex contributors (-); Number of collisions in one ITSROF (-)", HistType::kTH2D, {{101, -0.5, 100.5}, {11, -0.5, 10.5}});

    histos.add("Runs/hStdModeCollDG", ";Run number;Number of events (-)", HistType::kTH1D, {axisRunNumbers});
    histos.add("Runs/hUpcModeCollDG", ";Run number;Number of events (-)", HistType::kTH1D, {axisRunNumbers});
    histos.add("Runs/hStdModeCollSG1", ";Run number;Number of events (-)", HistType::kTH1D, {axisRunNumbers});
    histos.add("Runs/hUpcModeCollSG1", ";Run number;Number of events (-)", HistType::kTH1D, {axisRunNumbers});
    histos.add("Runs/hStdModeCollSG0", ";Run number;Number of events (-)", HistType::kTH1D, {axisRunNumbers});
    histos.add("Runs/hUpcModeCollSG0", ";Run number;Number of events (-)", HistType::kTH1D, {axisRunNumbers});
    histos.add("Runs/hStdModeCollNG", ";Run number;Number of events (-)", HistType::kTH1D, {axisRunNumbers});
    histos.add("Runs/hUpcModeCollNG", ";Run number;Number of events (-)", HistType::kTH1D, {axisRunNumbers});

  } // end init

  void processCounterPerITSROF(BCsWithRun3Matchings const& bcs, CCs const& collisions)
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

    for (const auto& bc : bcs) {
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
      for (const auto& collision : collisions) {
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
    } // end loop over bcs

    int arrAllColls[1000] = {0};
    int arrUPCcolls[1000] = {0};

    // next is based on matching of collision bc within ITSROF range in bcs
    for (const auto& itsrofBorder : vecITSROFborders) {
      int nAllCollsInROF = 0;
      int nUpcCollsInROF = 0;
      for (const auto& collision : collisions) {
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

    // TEST vertex contributors per reconstruction flag (std vs upc)
    // matching of collision bc within ITSROF range in bcs
    for (const auto& itsrofBorder : vecITSROFborders) {
      std::vector<int> vecNumContribsStd;
      std::vector<int> vecNumContribsUpc;
      for (const auto& collision : collisions) {
        if ((itsrofBorder.first < collision.bcId()) && (collision.bcId() < itsrofBorder.second)) {
          if (collision.flags() & dataformats::Vertex<o2::dataformats::TimeStamp<int>>::Flags::UPCMode) {
            vecNumContribsUpc.push_back(collision.numContrib());
          } else {
            vecNumContribsStd.push_back(collision.numContrib());
          }
        }
      } // end loop over collisions

      for (const auto& numContribs : vecNumContribsStd) {
        histos.get<TH2>(HIST("Events/hPVcontribsVsCollisionsPerITSROFstd"))->Fill(numContribs, vecNumContribsStd.size());
      }
      for (const auto& numContribs : vecNumContribsUpc) {
        histos.get<TH2>(HIST("Events/hPVcontribsVsCollisionsPerITSROFupc"))->Fill(numContribs, vecNumContribsUpc.size());
      }

    } // end loop over ITSROFs
  }

  void processCounterPerRun(FullSGUDCollision const& coll)
  {

    int gapSide = coll.gapSide();
    int trueGapSide = sgSelector.trueGap(coll, cutMyGapSideFV0, cutMyGapSideFT0A, cutMyGapSideFT0C, cutMyGapSideZDC);
    if (useTrueGap) {
      gapSide = trueGapSide;
    }

    if (coll.flags() == 0) {
      if (gapSide == 2) {
        histos.get<TH1>(HIST("Runs/hStdModeCollDG"))->Fill(coll.runNumber());
      } else if (gapSide == 1) {
        histos.get<TH1>(HIST("Runs/hStdModeCollSG1"))->Fill(coll.runNumber());
      } else if (gapSide == 0) {
        histos.get<TH1>(HIST("Runs/hStdModeCollSG0"))->Fill(coll.runNumber());
      } else {
        histos.get<TH1>(HIST("Runs/hStdModeCollNG"))->Fill(coll.runNumber());
      }
    } else {
      if (gapSide == 2) {
        histos.get<TH1>(HIST("Runs/hUpcModeCollDG"))->Fill(coll.runNumber());
      } else if (gapSide == 1) {
        histos.get<TH1>(HIST("Runs/hUpcModeCollSG1"))->Fill(coll.runNumber());
      } else if (gapSide == 0) {
        histos.get<TH1>(HIST("Runs/hUpcModeCollSG0"))->Fill(coll.runNumber());
      } else {
        histos.get<TH1>(HIST("Runs/hUpcModeCollNG"))->Fill(coll.runNumber());
      }
    }
  }

  PROCESS_SWITCH(UpcEventITSROFcounter, processCounterPerITSROF, "Counts number of collisions per ITSROF", false);
  PROCESS_SWITCH(UpcEventITSROFcounter, processCounterPerRun, "Counts number of whatever per RUN", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<UpcEventITSROFcounter>(cfgc)};
}
