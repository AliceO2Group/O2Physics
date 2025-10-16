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
#include "tableHMPIDpp.h"

#include "Common/Core/PID/PIDTOF.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
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

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

struct HmpidTableProducerPp {

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  const AxisSpec axisEvtCounter{1, 0, +1, ""};

  // CCDB configurable
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  struct : ConfigurableGroup {
    Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "URL of the CCDB repository"};
  } ccdbConfig;

  Produces<aod::HmpidAnalysisPp> hmpidAnalysis;

  // using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>;

  using CollisionCandidates = o2::soa::Join<o2::aod::Collisions, o2::aod::EvSels>;

  using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
                                    aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullDe,
                                    aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFFullDe>;

  void init(o2::framework::InitContext&)
  {
    // Configure CCDB
    ccdb->setURL(ccdbConfig.ccdbUrl);
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

    // loop on hmpid tracks
    for (const auto& t : hmpids) {
      // global tracks associated to hmpid tracks
      const auto& globalTrack = t.track_as<TrackCandidates>();
      if (!globalTrack.isGlobalTrack())
        continue;
      if (!globalTrack.hasITS() || !globalTrack.hasTPC() || !globalTrack.hasTOF())
        continue;

      float hmpidPhotsCharge2[o2::aod::kDimPhotonsCharge];

      for (int i = 0; i < o2::aod::kDimPhotonsCharge; i++) {
        hmpidPhotsCharge2[i] = t.hmpidPhotsCharge()[i];
      }

      /////FILL TABLE
      hmpidAnalysis(
        t.hmpidSignal(), globalTrack.phi(), globalTrack.eta(), t.hmpidMom(),
        globalTrack.p(), t.hmpidXTrack(), t.hmpidYTrack(), t.hmpidXMip(),
        t.hmpidYMip(), t.hmpidNPhotons(), t.hmpidQMip(), (t.hmpidClusSize() % 1000000) / 1000,
        t.hmpidClusSize() / 1000000, hmpidPhotsCharge2, globalTrack.eta(), globalTrack.phi(),
        globalTrack.px(), globalTrack.py(), globalTrack.pz(), globalTrack.itsNCls(),
        globalTrack.tpcNClsFound(), globalTrack.tpcNClsCrossedRows(), globalTrack.tpcChi2NCl(), globalTrack.itsChi2NCl(),
        globalTrack.dcaXY(), globalTrack.dcaZ(), globalTrack.tpcNSigmaPi(), globalTrack.tofNSigmaPi(),
        globalTrack.tpcNSigmaKa(), globalTrack.tofNSigmaKa(), globalTrack.tpcNSigmaPr(), globalTrack.tofNSigmaPr(),
        globalTrack.tpcNSigmaDe(), globalTrack.tofNSigmaDe());
    } // end hmpid tracks loop
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg) { return WorkflowSpec{adaptAnalysisTask<HmpidTableProducerPp>(cfg)}; }
