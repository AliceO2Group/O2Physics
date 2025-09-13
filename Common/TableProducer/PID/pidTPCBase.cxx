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

///
/// \file   pidTPCBase.cxx
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \brief  Base to build tasks for TPC PID tasks.
///

#include "pidTPCBase.h"

#include "Common/CCDB/ctpRateFetcher.h"
#include "Common/Core/TableHelper.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"

#include <CCDB/BasicCCDBManager.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/InitContext.h>
#include <Framework/SliceCache.h>
#include <Framework/runDataProcessing.h>

#include <cstdint>
#include <cstdlib>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::track;

struct PidMultiplicity {
  SliceCache cache;
  Produces<aod::PIDMults> mult;

  bool enableTable = false;
  void init(InitContext& initContext)
  {
    LOG(info) << "Initializing PID Mult Task";
    // Checking that the table is requested in the workflow and enabling it
    enableTable = isTableRequiredInWorkflow(initContext, "PIDMults");
    if (enableTable) {
      LOG(info) << "Table TPC PID Multiplicity enabled!";
    }
    if (doprocessStandard == true && doprocessIU == true) {
      LOG(fatal) << "Both processStandard and processIU are enabled, pick one!";
    }
  }

  using TrksIU = soa::Join<aod::TracksIU, aod::TracksExtra>;
  Partition<TrksIU> tracksWithTPCIU = (aod::track::tpcNClsFindable > (uint8_t)0);
  void processIU(aod::Collision const& collision, TrksIU const&)
  {
    if (!enableTable) {
      return;
    }
    auto tracksGrouped = tracksWithTPCIU->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    mult(tracksGrouped.size());
  }
  PROCESS_SWITCH(PidMultiplicity, processIU, "Process with IU tracks, faster but works on Run3 only", false);

  using Trks = soa::Join<aod::Tracks, aod::TracksExtra>;
  Partition<Trks> tracksWithTPC = (aod::track::tpcNClsFindable > (uint8_t)0);
  void processStandard(aod::Collision const& collision, Trks const&)
  {
    if (!enableTable) {
      return;
    }
    auto tracksGrouped = tracksWithTPC->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    mult(tracksGrouped.size());
  }
  PROCESS_SWITCH(PidMultiplicity, processStandard, "Process with tracks, needs propagated tracks", true);
};

struct DeDxCorrection {
  Produces<aod::DEdxsCorrected> dEdxCorrected;
  using BCsRun3 = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;
  using ColEvSels = soa::Join<aod::Collisions, aod::EvSels, aod::Mults>;
  using FullTracksIU = soa::Join<aod::TracksIU, aod::TracksExtra>;

  uint64_t minGlobalBC = 0;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  ctpRateFetcher mRateFetcher;

  Str_dEdx_correction str_dedx_correction;

  // void init(InitContext& initContext)
  void init(o2::framework::InitContext&)
  {
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    str_dedx_correction.init();
  }

  void processRun3(
    ColEvSels const& cols,
    FullTracksIU const& tracks,
    aod::BCsWithTimestamps const& bcs)
  {
    const uint64_t outTable_size = tracks.size();
    dEdxCorrected.reserve(outTable_size);

    for (auto const& trk : tracks) {
      double hadronicRate;
      int multTPC;
      int occupancy;
      if (trk.has_collision()) {
        auto collision = cols.iteratorAt(trk.collisionId());
        auto bc = collision.bc_as<aod::BCsWithTimestamps>();
        const int runnumber = bc.runNumber();
        hadronicRate = mRateFetcher.fetch(ccdb.service, bc.timestamp(), runnumber, "ZNC hadronic") * 1.e-3; // kHz
        multTPC = collision.multTPC();
        occupancy = collision.trackOccupancyInTimeRange();
      } else {
        auto bc = bcs.begin();
        const int runnumber = bc.runNumber();
        hadronicRate = mRateFetcher.fetch(ccdb.service, bc.timestamp(), runnumber, "ZNC hadronic") * 1.e-3; // kHz
        multTPC = 0;
        occupancy = 0;
      }

      float fTPCSignal = trk.tpcSignal();
      float fNormMultTPC = multTPC / 11000.;

      float fTrackOccN = occupancy / 1000.;
      float fOccTPCN = fNormMultTPC * 10; //(fNormMultTPC*10).clip(0,12)
      if (fOccTPCN > 12)
        fOccTPCN = 12;
      else if (fOccTPCN < 0)
        fOccTPCN = 0;

      float fTrackOccMeanN = hadronicRate / 5;
      float side = trk.tgl() > 0 ? 1 : 0;
      float a1pt = std::abs(trk.signed1Pt());
      float a1pt2 = a1pt * a1pt;
      float atgl = std::abs(trk.tgl());
      float mbb0R = 50 / fTPCSignal;
      if (mbb0R > 1.05)
        mbb0R = 1.05;
      else if (mbb0R < 0.05)
        mbb0R = 0.05;
      // float mbb0R =  max(0.05,  min(50 / fTPCSignal, 1.05));
      float a1ptmbb0R = a1pt * mbb0R;
      float atglmbb0R = atgl * mbb0R;

      std::vector<float> vec_occu = {fTrackOccN, fOccTPCN, fTrackOccMeanN};
      std::vector<float> vec_track = {mbb0R, a1pt, atgl, atglmbb0R, a1ptmbb0R, side, a1pt2};

      float fTPCSignalN_CR0 = str_dedx_correction.fReal_fTPCSignalN(vec_occu, vec_track);

      float mbb0R1 = 50 / (fTPCSignal / fTPCSignalN_CR0);
      if (mbb0R1 > 1.05)
        mbb0R1 = 1.05;
      else if (mbb0R1 < 0.05)
        mbb0R1 = 0.05;

      std::vector<float> vec_track1 = {mbb0R1, a1pt, atgl, atgl * mbb0R1, a1pt * mbb0R1, side, a1pt2};
      float fTPCSignalN_CR1 = str_dedx_correction.fReal_fTPCSignalN(vec_occu, vec_track1);

      float corrected_dEdx = fTPCSignal / fTPCSignalN_CR1;
      dEdxCorrected(corrected_dEdx);
    }
  }
  PROCESS_SWITCH(DeDxCorrection, processRun3, "dEdx correction process", false);

  void processDummy(aod::Collisions const&) {}

  PROCESS_SWITCH(DeDxCorrection, processDummy, "Do nothing", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<PidMultiplicity>(cfgc),
                      adaptAnalysisTask<DeDxCorrection>(cfgc)};
}
