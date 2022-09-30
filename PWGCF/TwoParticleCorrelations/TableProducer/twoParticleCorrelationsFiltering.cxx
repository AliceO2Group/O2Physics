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

#include "Framework/AnalysisTask.h"
#include "PWGCF/TwoParticleCorrelations/Core/EventSelectionFilterAndAnalysis.h"
#include "PWGCF/TwoParticleCorrelations/Core/TrackSelectionFilterAndAnalysis.h"
#include "PWGCF/TwoParticleCorrelations/DataModel/TwoParticleCorrelationsSkimmed.h"
#include "PWGCF/TwoParticleCorrelations/DataModel/TwoParticleCorrelationsFiltered.h"
#include "Framework/runDataProcessing.h"
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TParameter.h>
#include <TList.h>
#include <TDirectory.h>
#include <TFolder.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TProfile3D.h>

#include <cmath>

using namespace o2;
using namespace o2::framework;
using namespace o2::soa;
using namespace o2::framework::expressions;
using namespace o2::analysis;

namespace o2::analysis::twopfilter
{
#define TWOPFILTERLOGCOLLISIONS info
#define TWOPFILTERLOGTRACKS info

PWGCF::EventSelectionFilterAndAnalysis* fCollisionFilter = nullptr;
PWGCF::TrackSelectionFilterAndAnalysis* fTrackFilter = nullptr;
PWGCF::PIDSelectionFilterAndAnalysis* fPIDFilter = nullptr;
} // namespace o2::analysis::twopfilter

using namespace o2::aod::cfskim;

struct TwoParticleCorrelationsFilter {
  Produces<aod::TwoPAcceptedCollisions> acceptedcollisions;
  Produces<aod::TwoPFilteredTracks> accepteddtracks;
  Produces<aod::TwoPAcceptedGenCollisions> acceptedgencollisions;
  Produces<aod::TwoPFilteredParticles> acceptedgentracks;

#include "skimmingconf.h"

  int nReportedTracks;
  //  HistogramRegistry historeg;

  uint64_t collisionmask = 0UL;
  uint64_t collisionmask_opt = 0UL;
  uint64_t collisionmask_forced = 0UL;
  uint64_t trackmask = 0UL;
  uint64_t trackmask_opt = 0UL;
  uint64_t trackmask_forced = 0UL;
  uint64_t pidmask = 0UL;
  uint64_t pidmask_opt = 0UL;
  uint64_t pidmask_forced = 0UL;

  void init(InitContext const&)
  {
    using namespace twopfilter;

    LOGF(info, "TwoParticleCorrelationsFilter::init()");

    /* collision filtering configuration */
    PWGCF::EventSelectionConfigurable eventsel(eventfilter.centmultsel, {}, eventfilter.zvtxsel, {});
    fCollisionFilter = new PWGCF::EventSelectionFilterAndAnalysis(eventsel, PWGCF::SelectionFilterAndAnalysis::kAnalysis);

    /* track filtering configuration */
    PWGCF::TrackSelectionConfigurable trksel(trackfilter.ttype, trackfilter.nclstpc, trackfilter.nxrtpc, trackfilter.nclsits, trackfilter.chi2clustpc,
                                             trackfilter.chi2clusits, trackfilter.xrofctpc, trackfilter.dcaxy, trackfilter.dcaz, trackfilter.ptrange, trackfilter.etarange);
    fTrackFilter = new PWGCF::TrackSelectionFilterAndAnalysis(trksel, PWGCF::SelectionFilterAndAnalysis::kAnalysis);
    PWGCF::PIDSelectionConfigurable pidsel(pidfilter.pidtpcfilter.tpcel, pidfilter.pidtpcfilter.tpcmu, pidfilter.pidtpcfilter.tpcpi, pidfilter.pidtpcfilter.tpcka, pidfilter.pidtpcfilter.tpcpr,
                                           pidfilter.pidtoffilter.tpcel, pidfilter.pidtoffilter.tpcmu, pidfilter.pidtoffilter.tpcpi, pidfilter.pidtoffilter.tpcka, pidfilter.pidtoffilter.tpcpr,
                                           pidfilter.pidbayesfilter.bayel, pidfilter.pidbayesfilter.baymu, pidfilter.pidbayesfilter.baypi, pidfilter.pidbayesfilter.bayka, pidfilter.pidbayesfilter.baypr);
    fPIDFilter = new PWGCF::PIDSelectionFilterAndAnalysis(pidsel, PWGCF::SelectionFilterAndAnalysis::kFilter);

    nReportedTracks = 0;
    collisionmask = fCollisionFilter->getMask();
    collisionmask_opt = fCollisionFilter->getOptMask();
    collisionmask_forced = fCollisionFilter->getForcedMask();
    trackmask = fTrackFilter->getMask();
    trackmask_opt = fTrackFilter->getOptMask();
    trackmask_forced = fTrackFilter->getForcedMask();
    pidmask = fPIDFilter->getMask();
    pidmask_opt = fPIDFilter->getOptMask();
    pidmask_forced = fPIDFilter->getForcedMask();
    LOGF(info, "TwoParticleCorrelationsFilter::init(), collision selection masks 0x%08lx, 0x%08lx, and 0x%08lx ", collisionmask, collisionmask_opt, collisionmask_forced);
    LOGF(info, "TwoParticleCorrelationsFilter::init(), track selection masks 0x%08lx, 0x%08lx, and 0x%08lx ", trackmask, trackmask_opt, trackmask_forced);
    LOGF(info, "TwoParticleCorrelationsFilter::init(), PID selection masks 0x%08lx, 0x%08lx, and 0x%08lx ", pidmask, pidmask_opt, pidmask_forced);
    if (collisionmask == uint64_t(0) or trackmask == uint64_t(0)) {
      LOGF(fatal, "TwoParticleCorrelationsFilter::init() null masks, selecting everything!!!");
    }
    /* TODO: check the cuts signatures against the CCDB contents */
    LOGF(info, "Collision skimming signature: %s", fCollisionFilter->getCutStringSignature().Data());
    LOGF(info, "Track skimming signature: %s", fTrackFilter->getCutStringSignature().Data());
    LOGF(info, "PID skimming signature: %s", fPIDFilter->getCutStringSignature().Data());
  }

  Filter onlyacceptedcolls = ((aod::cfskim::selflags & static_cast<uint64_t>(collisionmask_forced)) == static_cast<uint64_t>(collisionmask_forced));
  Filter onlyacceptedtracks = ((aod::cfskim::trackflags & static_cast<uint64_t>(trackmask_forced)) == static_cast<uint64_t>(trackmask_forced));

  void processRun2(soa::Filtered<aod::CFCollisions>::iterator const& collision, soa::Filtered<aod::CFTracks> const& tracks)
  {
    LOGF(TWOPFILTERLOGCOLLISIONS, "Received collision with mask 0x%016lx and %ld tracks", collision.selflags(), tracks.size());

    /* for some reason we cannot apply this condition in the filter, it does not work */
    if ((collision.selflags() & collisionmask_opt) != 0UL) {
      acceptedcollisions(collision.posZ(), collision.centmult()[0]);
      int nAcceptedTracks = 0;
      for (const auto& track : tracks) {
        /* for some reason we cannot apply this condition in the filter, it does not work */
        if ((track.trackflags() & trackmask_opt) != 0UL) {
          accepteddtracks(acceptedcollisions.lastIndex(), 0, track.pt(), track.eta(), track.phi());
          nAcceptedTracks++;
        }
      }
      LOGF(TWOPFILTERLOGCOLLISIONS, ">> Accepted collision with mask 0x%08lx and %d accepted tracks", collision.selflags(), nAcceptedTracks);
    }
  }
  PROCESS_SWITCH(TwoParticleCorrelationsFilter, processRun2, "Process Run 2 two particle correlations filtering", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<TwoParticleCorrelationsFilter>(cfgc)};
  return workflow;
}
