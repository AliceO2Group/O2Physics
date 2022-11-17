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

#include <CCDB/BasicCCDBManager.h>
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "PWGCF/TwoParticleCorrelations/DataModel/TwoParticleCorrelationsSkimmed.h"
#include "PWGCF/TwoParticleCorrelations/Core/FilterAndAnalysisFramework.h"
#include "Framework/runDataProcessing.h"
#include "DataFormatsParameters/GRPObject.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::soa;
using namespace o2::framework::expressions;
using namespace o2::analysis;

namespace o2::analysis::cfskim
{
#define LOGTRACKTRACKS debug
#ifdef INCORPORATEBAYESIANPID
using pidTables = soa::Join<aod::pidTPCEl, aod::pidTPCMu, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr,
                            aod::pidTOFEl, aod::pidTOFMu, aod::pidTOFPi, aod::pidTOFKa, aod::pidTOFPr,
                            aod::pidBayesEl, aod::pidBayesMu, aod::pidBayesPi, aod::pidBayesKa, aod::pidBayesPr>;
#else
using pidTables = soa::Join<aod::pidTPCEl, aod::pidTPCMu, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr,
                            aod::pidTOFEl, aod::pidTOFMu, aod::pidTOFPi, aod::pidTOFKa, aod::pidTOFPr>;
#endif

PWGCF::FilterAndAnalysisFramework* fFilterFramework = nullptr;
} // namespace o2::analysis::cfskim

using namespace cfskim;

struct TwoParticleCorrelationsTrackSkimming {
  /* skimmed data tables */
  Produces<aod::CFTrackMask> trackmask;
  Produces<aod::CFTrackPIDs> skimmtrackpid;
  Produces<aod::CFMCPartMask> particlemask;

#include "PWGCF/TwoParticleCorrelations/TableProducer/Productions/skimmingconf_20221115.h"

  void init(InitContext const&)
  {
    using namespace cfskim;

    LOGF(info, "DptDptSkimTask::init()");

    /* collision filtering configuration */
    PWGCF::EventSelectionConfigurable eventsel(eventfilter.bfield, eventfilter.centmultsel, {}, eventfilter.zvtxsel, eventfilter.pileuprej);
    /* track filtering configuration */
    PWGCF::TrackSelectionConfigurable trksel(trackfilter.ttype, trackfilter.nclstpc, trackfilter.nxrtpc, trackfilter.nclsits, trackfilter.chi2clustpc,
                                             trackfilter.chi2clusits, trackfilter.xrofctpc, trackfilter.dcaxy, trackfilter.dcaz, trackfilter.ptrange, trackfilter.etarange);
#ifdef INCORPORATEBAYESIANPID
    PWGCF::PIDSelectionConfigurable pidsel(pidfilter.pidtpcfilter.tpcel, pidfilter.pidtpcfilter.tpcmu, pidfilter.pidtpcfilter.tpcpi, pidfilter.pidtpcfilter.tpcka, pidfilter.pidtpcfilter.tpcpr,
                                           pidfilter.pidtoffilter.tpcel, pidfilter.pidtoffilter.tpcmu, pidfilter.pidtoffilter.tpcpi, pidfilter.pidtoffilter.tpcka, pidfilter.pidtoffilter.tpcpr,
                                           pidfilter.pidbayesfilter.bayel, pidfilter.pidbayesfilter.baymu, pidfilter.pidbayesfilter.baypi, pidfilter.pidbayesfilter.bayka, pidfilter.pidbayesfilter.baypr);
#else
    PWGCF::PIDSelectionConfigurable pidsel(pidfilter.pidtpcfilter.tpcel, pidfilter.pidtpcfilter.tpcmu, pidfilter.pidtpcfilter.tpcpi, pidfilter.pidtpcfilter.tpcka, pidfilter.pidtpcfilter.tpcpr,
                                           pidfilter.pidtoffilter.tpcel, pidfilter.pidtoffilter.tpcmu, pidfilter.pidtoffilter.tpcpi, pidfilter.pidtoffilter.tpcka, pidfilter.pidtoffilter.tpcpr);
#endif

    fFilterFramework = new PWGCF::FilterAndAnalysisFramework(filterccdb.ccdburl.value, filterccdb.ccdbpath.value, filterccdb.filterdate.value);
    fFilterFramework->SetConfiguration(eventsel, trksel, pidsel, PWGCF::SelectionFilterAndAnalysis::kFilter);
    fFilterFramework->Init();

    /* TODO: upload the cuts signatures to the CCDB */
    LOGF(info, "Collision skimming signature: %s", fFilterFramework->getEventFilterCutStringSignature().Data());
    LOGF(info, "Track skimming signature: %s", fFilterFramework->getTrackFilterCutStringSignature().Data());
    LOGF(info, "PID skimming signature: %s", fFilterFramework->getPIDFilterCutStringSignature().Data());
  }

  void processRun2(soa::Join<aod::FullTracks, aod::TracksDCA, pidTables> const& tracks, soa::Join<aod::Collisions, aod::CFCollMasks> const&)
  {
    trackmask.reserve(tracks.size());
    skimmtrackpid.reserve(tracks.size());

    int nfilteredtracks = 0;
    for (auto const& track : tracks) {
      if (!track.has_collision()) {
        /* track not assigned to any collision */
        trackmask(0UL);
        skimmtrackpid(0UL);
      } else {
        auto trkmask = fFilterFramework->FilterTrack(track);
        auto pidmask = fFilterFramework->FilterTrackPID(track);
        trackmask(trkmask);
        skimmtrackpid(pidmask);
        if (trkmask != 0UL) {
          nfilteredtracks++;
        }
      }
    }
    LOGF(info, "Filtered %d tracks out of %d", nfilteredtracks, tracks.size());
  }
  PROCESS_SWITCH(TwoParticleCorrelationsTrackSkimming, processRun2, "Process on Run 1 or Run 2 data, i.e. do not store derived data ", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<TwoParticleCorrelationsTrackSkimming>(cfgc)};
  return workflow;
}
