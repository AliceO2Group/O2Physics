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

#include "PWGCF/TwoParticleCorrelations/Core/FilterAndAnalysisFramework.h"
#include "PWGCF/TwoParticleCorrelations/DataModel/TwoParticleCorrelationsSkimmed.h"

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "DataFormatsParameters/GRPObject.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include <CCDB/BasicCCDBManager.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::soa;
using namespace o2::framework::expressions;
using namespace o2::analysis;

#include <string>

namespace o2::analysis::cfskim
{
#define LOGTRACKCOLLISIONS debug
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

void setEventCutsLabels(std::shared_ptr<TH1> h)
{
  using namespace aod::run2;

  /* labels taken from O2/Framework DataTypes.h */
  const char* evcutslabel[kTRDHEE + 1] = {
    "kINELgtZERO",
    "kPileupInMultBins",
    "kConsistencySPDandTrackVertices",
    "kTrackletsVsClusters",
    "kNonZeroNContribs",
    "kIncompleteDAQ",
    "kPileUpMV",
    "kTPCPileUp",
    "kTimeRangeCut",
    "kEMCALEDCut",
    "kAliEventCutsAccepted",
    "kIsPileupFromSPD",
    "kIsV0PFPileup",
    "kIsTPCHVdip",
    "kIsTPCLaserWarmUp",
    "kTRDHCO",
    "kTRDHJT",
    "kTRDHSE",
    "kTRDHQU",
    "kTRDHEE"};

  for (int bit = kINELgtZERO; bit <= kTRDHEE; ++bit) {
    h->GetXaxis()->SetBinLabel(bit + 1, evcutslabel[bit]);
  }
}

void reportEventCuts(std::shared_ptr<TH1> h, uint32_t eventcuts)
{
  using namespace aod::run2;
  auto entries = h->GetEntries();
  for (int bit = kINELgtZERO; bit <= kTRDHEE; ++bit) {
    if (TESTBIT(eventcuts, bit)) {
      h->Fill(bit + 0.5);
    }
  }
  h->SetEntries(entries + 1);
}

struct TwoParticleCorrelationsCollisionSkimming {
  /* skimmed data tables */
  Produces<aod::CFCollMask> collisionmask;
  Produces<aod::CFMCCollMasks> gencollisionmask;

  Service<o2::ccdb::BasicCCDBManager> ccdb;

#include "PWGCF/TwoParticleCorrelations/TableProducer/Productions/skimmingconf_20221115.h" // NOLINT

  int nReportedTracks;
  int runNumber = 0;
  int bfield = 0;
  HistogramRegistry historeg;

  int getMagneticField(std::string ccdbpath, uint64_t timestamp)
  {
    // TODO done only once (and not per run). Will be replaced by CCDBConfigurable
    static o2::parameters::GRPObject* grpo = nullptr;
    if (grpo == nullptr) {
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(ccdbpath, timestamp);
      if (grpo == nullptr) {
        LOGF(fatal, "GRP object not found for timestamp %llu", timestamp);
        return 0;
      }
      LOGF(info, "Retrieved GRP for timestamp %llu with magnetic field of %d kG", timestamp, grpo->getNominalL3Field());
    }
    return grpo->getNominalL3Field();
  }

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
    nReportedTracks = 0;

    /* TODO: upload the cuts signatures to the CCDB */
    LOGF(info, "Collision skimming signature: %s", fFilterFramework->getEventFilterCutStringSignature().Data());
    LOGF(info, "Track skimming signature: %s", fFilterFramework->getTrackFilterCutStringSignature().Data());
    LOGF(info, "PID skimming signature: %s", fFilterFramework->getPIDFilterCutStringSignature().Data());

    historeg.add("EventCuts", "EventCuts", {HistType::kTH1F, {{aod::run2::kTRDHEE + 1, 0, aod::run2::kTRDHEE + 1}}});
    setEventCutsLabels(historeg.get<TH1>(HIST("EventCuts")));

    /* initialize access to the CCDB */
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
  }

  template <typename Coll, typename AssociatedTracks, typename BcInfo>
  uint64_t filterRun2Collision(Coll const& collision, AssociatedTracks const& tracks, BcInfo const& bcinfo)
  {
    using namespace aod::run2;
    using namespace aod::collision;

    uint32_t eventcuts = bcinfo.eventCuts();
    bool accepted = true;

    // NOT CONFIGURABLE EVENT SELECTION
    /* incomplete data acquisition */
    if (!TESTBIT(eventcuts, kIncompleteDAQ)) {
      accepted = false;
    }
    /* pile-up */
    /* TODO: check if this is also valid for Run 1 data */
    if (!TESTBIT(eventcuts, kPileupInMultBins) ||
        !TESTBIT(eventcuts, kTrackletsVsClusters) ||
        !TESTBIT(eventcuts, kPileUpMV) ||
        !TESTBIT(eventcuts, kTimeRangeCut) ||
        !TESTBIT(eventcuts, kTPCPileUp) ||
        TESTBIT(eventcuts, kIsPileupFromSPD) ||
        TESTBIT(eventcuts, kIsV0PFPileup)) {
      accepted = false;
    }
    /* TPC issues*/
    if (TESTBIT(eventcuts, kIsTPCHVdip) ||
        TESTBIT(eventcuts, kIsTPCLaserWarmUp)) {
      accepted = false;
    }
    /* vertex */
    if (!TESTBIT(eventcuts, kNonZeroNContribs) ||
        (((collision.flags() & Run2VertexerZ) == Run2VertexerZ) && collision.covZZ() < 0.25)) {
      accepted = false;
    }
    reportEventCuts(historeg.get<TH1>(HIST("EventCuts")), eventcuts);

    // CONFIGURABLE EVENT SELECTION
    /* update magnetic field if needed */
    if (bcinfo.runNumber() != runNumber) {
      bfield = getMagneticField("GLO/GRP/GRP", bcinfo.timestamp());
      runNumber = bcinfo.runNumber();
    }
    if (accepted) {
      return fFilterFramework->FilterCollision(collision, tracks, bfield);
    } else {
      return 0UL;
    }
  }

  void processRun2(soa::Join<aod::Collisions, aod::CentRun2V0Ms, aod::CentRun2CL0s, aod::CentRun2CL1s, aod::Mults>::iterator const& collision,
                   soa::Join<aod::BCs, aod::Timestamps, aod::Run2BCInfos> const&, soa::Join<aod::FullTracks, aod::TracksDCA, pidTables> const& tracks)
  {
    /* for the time being this will apply only to Run 1+2 converted data */
    LOGF(LOGTRACKCOLLISIONS, "Got a new collision with zvtx %.2f and V0M %.2f, CL0 %.2f, CL1 %.2f", collision.posZ(), collision.centRun2V0M(), collision.centRun2CL0(), collision.centRun2CL1());

    auto bc = collision.bc_as<soa::Join<aod::BCs, aod::Timestamps, aod::Run2BCInfos>>();
    auto colmask = filterRun2Collision(collision, tracks, bc);
    LOGF(LOGTRACKCOLLISIONS, "Got mask 0x%16lx", colmask);

    collisionmask(colmask, fFilterFramework->GetCollisionMultiplicities());
  }
  PROCESS_SWITCH(TwoParticleCorrelationsCollisionSkimming, processRun2, "Process on Run 1 or Run 2 data, i.e. do not store derived data ", true);
};

struct TwoParticleCorrelationsTrackSkimming {
  /* skimmed data tables */
  Produces<aod::CFTrackMask> trackmask;
  Produces<aod::CFTrackPIDs> skimmtrackpid;
  Produces<aod::CFMCPartMask> particlemask;

#include "PWGCF/TwoParticleCorrelations/TableProducer/Productions/skimmingconf_20221115.h" // NOLINT

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
    adaptAnalysisTask<TwoParticleCorrelationsCollisionSkimming>(cfgc),
    adaptAnalysisTask<TwoParticleCorrelationsTrackSkimming>(cfgc)};
  return workflow;
}
