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

/// \file twoParticleCorrelationsFiltering.cxx
/// \brief Filters collisions and tracks
/// \author victor.gonzalez.sebastian@gmail.com

#include "PWGCF/TwoParticleCorrelations/Core/EventSelectionFilterAndAnalysis.h"
#include "PWGCF/TwoParticleCorrelations/Core/FilterAndAnalysisFramework.h"
#include "PWGCF/TwoParticleCorrelations/Core/PIDSelectionFilterAndAnalysis.h"
#include "PWGCF/TwoParticleCorrelations/Core/SelectionFilterAndAnalysis.h"
#include "PWGCF/TwoParticleCorrelations/Core/TrackSelectionFilterAndAnalysis.h"
#include "PWGCF/TwoParticleCorrelations/DataModel/TwoParticleCorrelationsFiltered.h"
#include "PWGCF/TwoParticleCorrelations/DataModel/TwoParticleCorrelationsSkimmed.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <cmath>
#include <cstdint>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::soa;
using namespace o2::framework::expressions;
using namespace o2::analysis;

namespace o2::analysis::twopfilter
{
#define TWOPFILTERLOGCOLLISIONS debug
#define TWOPFILTERLOGTRACKS debug

uint64_t collisionmask = 0UL;
std::vector<uint64_t> collisionmaskOptional;
uint64_t collisionmaskForced = 0UL;
uint64_t trackmask = 0UL;
std::vector<uint64_t> trackmaskOptional;
uint64_t trackmaskForced = 0UL;
uint64_t pidmask = 0UL;
std::vector<uint64_t> pidmaskOptional;
uint64_t pidmaskForced = 0UL;

pwgcf::FilterAndAnalysisFramework* fFilterFramework = nullptr;

int fMultiplicityIndex = -1; //! the index to the multiplicity values array
} // namespace o2::analysis::twopfilter

using namespace o2::aod::cfskim;

struct TwoParticleCorrelationsFiltering {
  Produces<aod::TwoPAcceptedCollisions> acceptedcollisions;
  Produces<aod::TwoPFilteredTracks> accepteddtracks;
  Produces<aod::TwoPAcceptedGenCollisions> acceptedgencollisions;
  Produces<aod::TwoPFilteredParticles> acceptedgentracks;

#include "PWGCF/TwoParticleCorrelations/TableProducer/Productions/skimmingconf_20221115.h" // NOLINT

  int nReportedTracks;
  //  HistogramRegistry historeg;

  void init(InitContext const&)
  {
    using namespace twopfilter;

    LOGF(info, "TwoParticleCorrelationsFiltering::init()");

    /* collision filtering configuration */
    pwgcf::EventSelectionConfigurable eventsel(eventfilter.bfield, eventfilter.centmultsel, {}, eventfilter.zvtxsel, eventfilter.pileuprej);
    /* track filtering configuration */
    pwgcf::TrackSelectionConfigurable trksel(trackfilter.ttype, trackfilter.nclstpc, trackfilter.nxrtpc, trackfilter.nclsits, trackfilter.chi2clustpc,
                                             trackfilter.chi2clusits, trackfilter.xrofctpc, trackfilter.dcaxy, trackfilter.dcaz, trackfilter.ptrange, trackfilter.etarange);
#ifdef INCORPORATEBAYESIANPID
    pwgcf::PIDSelectionConfigurable pidsel(pidfilter.pidtpcfilter.tpcel, pidfilter.pidtpcfilter.tpcmu, pidfilter.pidtpcfilter.tpcpi, pidfilter.pidtpcfilter.tpcka, pidfilter.pidtpcfilter.tpcpr,
                                           pidfilter.pidtoffilter.tpcel, pidfilter.pidtoffilter.tpcmu, pidfilter.pidtoffilter.tpcpi, pidfilter.pidtoffilter.tpcka, pidfilter.pidtoffilter.tpcpr,
                                           pidfilter.pidbayesfilter.bayel, pidfilter.pidbayesfilter.baymu, pidfilter.pidbayesfilter.baypi, pidfilter.pidbayesfilter.bayka, pidfilter.pidbayesfilter.baypr);
#else
    pwgcf::PIDSelectionConfigurable pidsel(pidfilter.pidtpcfilter.tpcel, pidfilter.pidtpcfilter.tpcmu, pidfilter.pidtpcfilter.tpcpi, pidfilter.pidtpcfilter.tpcka, pidfilter.pidtpcfilter.tpcpr,
                                           pidfilter.pidtoffilter.tpcel, pidfilter.pidtoffilter.tpcmu, pidfilter.pidtoffilter.tpcpi, pidfilter.pidtoffilter.tpcka, pidfilter.pidtoffilter.tpcpr);
#endif
    fFilterFramework = new pwgcf::FilterAndAnalysisFramework(filterccdb.ccdburl.value, filterccdb.ccdbpath.value, filterccdb.filterdate.value);
    fFilterFramework->setConfiguration(eventsel, trksel, pidsel, pwgcf::SelectionFilterAndAnalysis::KAnalysis);
    fFilterFramework->init();
    nReportedTracks = 0;
    collisionmask = fFilterFramework->getCollisionMask();
    collisionmaskOptional = fFilterFramework->getCollisionOptMask();
    collisionmaskForced = fFilterFramework->getCollisionForcedMask();
    fMultiplicityIndex = fFilterFramework->getCollisionMultiplicityIndex();
    trackmask = fFilterFramework->getTrackMask();
    trackmaskOptional = fFilterFramework->getTrackOptMask();
    trackmaskForced = fFilterFramework->getTrackForcedMask();
    pidmask = fFilterFramework->getPIDMask();
    pidmaskOptional = fFilterFramework->getPIDOptMask();
    pidmaskForced = fFilterFramework->getPIDForcedMask();
    LOGF(info, "TwoParticleCorrelationsFiltering::init(), collision selection masks 0x%016lx, %s, and 0x%016lx and multiplicity index %d", collisionmask, fFilterFramework->printCollisionOptionalMasks().Data(), collisionmaskForced, fMultiplicityIndex);
    LOGF(info, "TwoParticleCorrelationsFiltering::init(), track selection masks 0x%016lx, %s, and 0x%016lx ", trackmask, fFilterFramework->printTrackOptionalMasks().Data(), trackmaskForced);
    LOGF(info, "TwoParticleCorrelationsFiltering::init(), PID selection masks 0x%016lx, %s, and 0x%016lx ", pidmask, fFilterFramework->printPIDOptionalMasks().Data(), pidmaskForced);
    if (collisionmask == static_cast<uint64_t>(0) || trackmask == static_cast<uint64_t>(0)) {
      LOGF(fatal, "TwoParticleCorrelationsFiltering::init() null masks, selecting everything!!!");
    }
    /* TODO: check the cuts signatures against the CCDB contents */
    LOGF(info, "Collision skimming signature: %s", fFilterFramework->getEventFilterCutStringSignature().Data());
    LOGF(info, "Track skimming signature: %s", fFilterFramework->getTrackFilterCutStringSignature().Data());
    LOGF(info, "PID skimming signature: %s", fFilterFramework->getPIDFilterCutStringSignature().Data());
  }

  void processRun2(soa::Join<aod::Collisions, aod::CFCollMasks>::iterator const& collision,
                   soa::Join<aod::FullTracks, aod::CFTrackMasks> const& tracks)
  {
    using namespace twopfilter;
    LOGF(TWOPFILTERLOGCOLLISIONS, "Received collision with mask 0x%016lx and %ld tracks", collision.selFlags(), tracks.size());
    auto passOptions = [](const auto& options, auto mask) {
      bool all = true;
      for (const auto& option : options) {
        all = all && ((option & mask) != 0UL);
      }
      return all;
    };

    if ((collision.selFlags() & collisionmaskForced) == collisionmaskForced && passOptions(collisionmaskOptional, collision.selFlags())) {
      LOGF(TWOPFILTERLOGCOLLISIONS, ">> Accepted collision with mask 0x%016lx and %ld unfiltered tracks", collision.selFlags(), tracks.size());
      acceptedcollisions(collision.centMult()[fMultiplicityIndex], uint8_t(true));
      int nAcceptedTracks = 0;
      for (const auto& track : tracks) {
        if ((track.trackFlags() & trackmaskForced) == trackmaskForced && passOptions(trackmaskOptional, track.trackFlags())) {
          accepteddtracks(0); // TODO: the kind of accepted track
          nAcceptedTracks++;
        } else {
          accepteddtracks(-1);
        }
      }
      LOGF(TWOPFILTERLOGCOLLISIONS, ">> Accepted collision with mask 0x%016lx and %d accepted tracks", collision.selFlags(), nAcceptedTracks);
    } else {
      acceptedcollisions(collision.centMult()[fMultiplicityIndex], uint8_t(false));
    }
  }
  PROCESS_SWITCH(TwoParticleCorrelationsFiltering, processRun2, "Process Run 2, i.e. over NOT stored derived data, two particle correlations filtering", true);

  void processRun3(aod::CFCollision const& collision, aod::CFTracks const& tracks)
  {
    using namespace twopfilter;
    LOGF(TWOPFILTERLOGCOLLISIONS, "Received collision with mask 0x%016lx and %ld tracks", collision.selFlags(), tracks.size());
    auto passOptions = [](const auto& options, auto mask) {
      bool all = true;
      for (const auto& option : options) {
        all = all && ((option & mask) != 0UL);
      }
      return all;
    };

    if ((collision.selFlags() & collisionmaskForced) == collisionmaskForced && passOptions(collisionmaskOptional, collision.selFlags())) {
      LOGF(TWOPFILTERLOGCOLLISIONS, ">> Accepted collision with mask 0x%016lx and %ld unfiltered tracks", collision.selFlags(), tracks.size());
      acceptedcollisions(collision.centMult()[fMultiplicityIndex], uint8_t(true));
      int nAcceptedTracks = 0;
      for (const auto& track : tracks) {
        if ((track.trackFlags() & trackmaskForced) == trackmaskForced && passOptions(trackmaskOptional, track.trackFlags())) {
          accepteddtracks(0); // TODO: the kind of accepted track
          nAcceptedTracks++;
        } else {
          accepteddtracks(-1);
        }
      }
      LOGF(TWOPFILTERLOGCOLLISIONS, ">> Accepted collision with mask 0x%016lx and %d accepted tracks", collision.selFlags(), nAcceptedTracks);
    } else {
      acceptedcollisions(collision.centMult()[fMultiplicityIndex], uint8_t(false));
    }
  }
  PROCESS_SWITCH(TwoParticleCorrelationsFiltering, processRun3, "Process Run 3, i.e. over stored derived data, two particle correlations filtering", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<TwoParticleCorrelationsFiltering>(cfgc)};
  return workflow;
}
