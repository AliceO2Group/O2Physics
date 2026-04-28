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

#include "tableHMPID.h"

#include "Common/Core/PID/PIDTOF.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/TableProducer/PID/pidTOFBase.h"

#include "CCDB/BasicCCDBManager.h"
#include "CCDB/CcdbApi.h"
#include <Framework/ASoA.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/RunningWorkflowInfo.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/DCA.h>
#include <ReconstructionDataFormats/PID.h>
#include <ReconstructionDataFormats/Track.h>
#include <ReconstructionDataFormats/TrackParametrization.h>

#include <TTree.h>

#include <string>
#include <unordered_set>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

struct HmpidTableProducer {
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  const AxisSpec axisEvtCounter{1, 0, +1, ""};

  // CCDB configurable
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  struct : ConfigurableGroup {
    Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "URL of the CCDB repository"};
  } ccdbConfig;

  Produces<aod::HmpidAnalysis> hmpidAnalysis;

  // configurable for quality requirements
  Configurable<bool> requireITS{"requireITS", true, "Require ITS track"};
  Configurable<bool> requireTPC{"requireTPC", true, "Require TPC track"};
  Configurable<bool> requireTOF{"requireTOF", true, "Require TOF track"};

  using CollisionCandidates = o2::soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFV0As>;

  using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
                                    aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullDe,
                                    aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFFullDe>;

  void init(o2::framework::InitContext&)
  {
    // Configure CCDB
    ccdb->setURL(ccdbConfig.ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    histos.add("eventCounter", "eventCounter", kTH1F, {axisEvtCounter});
    histos.add("goodEventCounter", "goodEventCounter", kTH1F, {axisEvtCounter});
    histos.add("eventsHmpid", "eventsWithHmpid", kTH1F, {axisEvtCounter});
  }

  // function to manage ccdb
  int mCCDBRunNumber = 0;
  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mCCDBRunNumber == bc.runNumber()) {
      return;
    }
    mCCDBRunNumber = bc.runNumber();
  }

  void processEvent(CollisionCandidates::iterator const& col,
                    aod::BCsWithTimestamps const&)
  {
    histos.fill(HIST("eventCounter"), 0.5);
    if (col.sel8()) {
      histos.fill(HIST("goodEventCounter"), 0.5);
    }

    // initialize CCDB for current BC
    initCCDB(col.bc_as<aod::BCsWithTimestamps>());
  }
  PROCESS_SWITCH(HmpidTableProducer, processEvent, "Process event level - collisions", true);

  void processHmpid(
    aod::HMPIDs const& hmpids,
    TrackCandidates const&,
    CollisionCandidates const&,
    aod::BCsWithTimestamps const&)
  {
    // --- Static set to track unique collisions with HMPID tracks ---
    static std::unordered_set<uint32_t> collisionsWithHmpid;

    for (auto const& t : hmpids) {

      // Access the global track associated to the HMPID track
      const auto& globalTrack = t.track_as<TrackCandidates>();

      if (!globalTrack.has_collision())
        continue;

      // Access the associated collision
      const auto& col = globalTrack.collision_as<CollisionCandidates>();
      initCCDB(col.bc_as<aod::BCsWithTimestamps>());
      uint32_t collId = col.globalIndex();

      // --- Track quality selection ---
      if ((requireITS && !globalTrack.hasITS()) ||
          (requireTPC && !globalTrack.hasTPC()) ||
          (requireTOF && !globalTrack.hasTOF())) {
        continue;
      }

      // Count collisions with at least one valid HMPID track
      if (collisionsWithHmpid.insert(collId).second) {
        histos.fill(HIST("eventsHmpid"), 0.5);
      }

      float centrality = col.centFV0A();

      float hmpidPhotsCharge2[o2::aod::kDimPhotonsCharge];

      for (int i = 0; i < o2::aod::kDimPhotonsCharge; i++) {
        hmpidPhotsCharge2[i] = t.hmpidPhotsCharge()[i];
      }

      /////FILL HMPID CUSTOM TABLE
      hmpidAnalysis(t.hmpidSignal(), t.hmpidMom(),
                    globalTrack.p(), t.hmpidXTrack(), t.hmpidYTrack(), t.hmpidXMip(),
                    t.hmpidYMip(), t.hmpidNPhotons(), t.hmpidQMip(), (t.hmpidClusSize() % 1000000) / 1000,
                    t.hmpidClusSize() / 1000000, hmpidPhotsCharge2, globalTrack.eta(), globalTrack.phi(),
                    globalTrack.px(), globalTrack.py(), globalTrack.pz(), globalTrack.itsNCls(),
                    globalTrack.tpcNClsFound(), globalTrack.tpcNClsCrossedRows(), globalTrack.tpcChi2NCl(), globalTrack.itsChi2NCl(),
                    globalTrack.dcaXY(), globalTrack.dcaZ(), globalTrack.tpcNSigmaPi(), globalTrack.tofNSigmaPi(),
                    globalTrack.tpcNSigmaKa(), globalTrack.tofNSigmaKa(), globalTrack.tpcNSigmaPr(), globalTrack.tofNSigmaPr(),
                    globalTrack.tpcNSigmaDe(), globalTrack.tofNSigmaDe(), centrality);
    } // end loop on hmpid table entries
  }

  PROCESS_SWITCH(HmpidTableProducer, processHmpid, "Process hmpid entries - tracks", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg) { return WorkflowSpec{adaptAnalysisTask<HmpidTableProducer>(cfg)}; }
