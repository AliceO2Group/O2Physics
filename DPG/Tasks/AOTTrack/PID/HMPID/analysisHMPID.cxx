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

// O2 includes
#include "Common/Core/PID/PIDTOF.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/TableProducer/PID/pidTOFBase.h"

#include "Framework/ASoA.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/PID.h"
#include "ReconstructionDataFormats/Track.h"
#include "ReconstructionDataFormats/TrackParametrization.h"
// #include "ReconstructionDataFormats/MatchInfoHMP.h"

#include "tableHMPID.h"

#include "CCDB/BasicCCDBManager.h"
#include "CCDB/CcdbApi.h"

#include <TTree.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

struct pidHmpidAnalysis {

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  const AxisSpec axisEvtCounter{1, 0, +1, ""};

  // CCDB configurable
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  struct : ConfigurableGroup {
    Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  } ccdbConfig;

  Produces<aod::HMPID_analysis> HMPID_analysis;

  // using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>;

  using CollisionCandidates = o2::soa::Join<o2::aod::Collisions, o2::aod::EvSels>;

  using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
                                    aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullDe,
                                    aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFFullDe>;

  void init(o2::framework::InitContext&)
  {
    // Configure CCDB
    ccdb->setURL(ccdbConfig.ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    histos.add("eventCounter", "eventCounter", kTH1F, {axisEvtCounter});
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

  void process(aod::Collision const& event,
               const aod::HMPIDs& hmpids,
               TrackCandidates const&,
               aod::BCsWithTimestamps const&)
  {
    // counter for all events processed (move into tracks loop if good event counter is needed)
    histos.fill(HIST("eventCounter"), 0.5);

    initCCDB(event.bc_as<aod::BCsWithTimestamps>());

    // loop sulle tracce hmpid
    for (const auto& t : hmpids) {
      // global tracks associated to hmpid tracks
      const auto& global_track = t.track_as<TrackCandidates>();
      if (!global_track.isGlobalTrack())
        continue;
      if (!global_track.hasITS() || !global_track.hasTPC() || !global_track.hasTOF())
        continue;

      float hmpidPhotsCharge2[10];

      for (int i = 0; i < 10; i++) {
        hmpidPhotsCharge2[i] = t.hmpidPhotsCharge()[i];
      }

      /////FILL TABLE
      HMPID_analysis(t.hmpidSignal(), global_track.phi(), global_track.eta(), t.hmpidMom(),
                     global_track.p(), t.hmpidXTrack(), t.hmpidYTrack(), t.hmpidXMip(),
                     t.hmpidYMip(), t.hmpidNPhotons(), t.hmpidQMip(), (t.hmpidClusSize() % 1000000) / 1000, t.hmpidClusSize() / 1000000,
                     hmpidPhotsCharge2, global_track.eta(), global_track.phi(), global_track.px(), global_track.py(), global_track.pz(),
                     global_track.itsNCls(), global_track.tpcNClsFound(), global_track.tpcNClsCrossedRows(), global_track.tpcChi2NCl(), global_track.itsChi2NCl(),
                     global_track.dcaXY(), global_track.dcaZ(), global_track.tpcNSigmaPi(), global_track.tofNSigmaPi(), global_track.tpcNSigmaKa(), global_track.tofNSigmaKa(),
                     global_track.tpcNSigmaPr(), global_track.tofNSigmaPr(), global_track.tpcNSigmaDe(), global_track.tofNSigmaDe());
    } // end hmpid tracks loop
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg) { return WorkflowSpec{adaptAnalysisTask<pidHmpidAnalysis>(cfg)}; }
