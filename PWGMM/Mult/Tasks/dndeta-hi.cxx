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

#include <cmath>
#include <iostream>
#include <chrono>
#include "bestCollisionTable.h"
#include "CCDB/BasicCCDBManager.h"
#include "Common/CCDB/EventSelectionParams.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "CommonConstants/MathConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/Configurable.h"
#include "Framework/RuntimeError.h"
#include "Framework/runDataProcessing.h"
#include "Index.h"
#include "ReconstructionDataFormats/GlobalTrackID.h"
#include "TDatabasePDG.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::track;
using BCsRun3 = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;

enum {
  kECbegin = 0,
  kDATA = 1,
  kINEL,
  kECend
};
enum {
  kTrigbegin = 0,
  kMBAND = 1,
  kTrigend
};

AxisSpec ZAxis = {60, -30, 30, "zaxis"};
AxisSpec DeltaZAxis = {61, -6.1, 6.1};
AxisSpec DCAAxis = {601, -3.01, 3.01};
AxisSpec EtaAxis = {80, -4.0, 4.0, "etaaxis"};
AxisSpec MultAxis = {1001, -0.5, 1000.5};
AxisSpec PhiAxis = {629, 0, 2 * M_PI};
AxisSpec PtAxis = {2401, -0.005, 24.005};
AxisSpec EvtClassAxis = {kECend - 1, kECbegin + 0.5, kECend - 0.5, "eventclass"};
AxisSpec TrigClassAxis = {kTrigend - 1, kTrigbegin + 0.5, kTrigend - 0.5, "triggclass"};
std::vector<double> centBinning = {0., 20, 60., 90., 100};
AxisSpec CentAxis = {centBinning, "centrality"};

static constexpr TrackSelectionFlags::flagtype trackSelectionITS =
  TrackSelectionFlags::kITSNCls | TrackSelectionFlags::kITSChi2NDF |
  TrackSelectionFlags::kITSHits;

static constexpr TrackSelectionFlags::flagtype trackSelectionTPC =
  TrackSelectionFlags::kTPCNCls |
  TrackSelectionFlags::kTPCCrossedRowsOverNCls |
  TrackSelectionFlags::kTPCChi2NDF;

static constexpr TrackSelectionFlags::flagtype trackSelectionDCA =
  TrackSelectionFlags::kDCAz | TrackSelectionFlags::kDCAxy;

using LabeledTracks = soa::Join<aod::Tracks, aod::McTrackLabels>;

struct MultiplicityCounter {
  Service<TDatabasePDG> pdg;

  Configurable<float> estimatorEta{"estimatorEta", 1.0, "eta range for INEL>0 sample definition"};
  Configurable<bool> useEvSel{"useEvSel", true, "use event selection"};
  Configurable<bool> isMC{"isMC", false, "check if MC"};
  Service<ccdb::BasicCCDBManager> ccdb;
  Configurable<std::string> path{"ccdb-path", "Users/s/sherrman/My/Object", "base path to the ccdb object"};
  Configurable<std::string> url{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<int64_t> nolaterthan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};

  HistogramRegistry registry{
    "registry",
    {
      {"Events/Selection", ";status;events", {HistType::kTH1F, {{7, 0.5, 7.5}}}},                                                                          //
      {"hrecdndeta", "evntclass; triggerclass; centrality, zvtex, eta", {HistType::kTHnSparseD, {EvtClassAxis, TrigClassAxis, CentAxis, ZAxis, EtaAxis}}}, //
      {"hgendndeta", "evntclass; centrality, zvtex, eta", {HistType::kTHnSparseD, {EvtClassAxis, CentAxis, ZAxis, EtaAxis}}},                              //
      {"hreczvtx", "evntclass; triggerclass; centrality, zvtex", {HistType::kTHnSparseD, {EvtClassAxis, TrigClassAxis, CentAxis, ZAxis}}},                 //
      {"hgenzvtx", "evntclass; centrality, zvtex", {HistType::kTHnSparseD, {EvtClassAxis, CentAxis, ZAxis}}},                                              //
      {"PhiEta", "; #varphi; #eta; tracks", {HistType::kTH2F, {PhiAxis, EtaAxis}}},                                                                        //
      {"DCAXY", " ; DCA_{XY} (cm)", {HistType::kTH1F, {DCAAxis}}},                                                                                         //
      {"DCAZ", " ; DCA_{Z} (cm)", {HistType::kTH1F, {DCAAxis}}},                                                                                           //
      {"Multiplicity", " ; FV0A (#); FT0A (#); FT0C (#) ", {HistType::kTHnSparseD, {MultAxis, MultAxis, MultAxis}}},                                       //
      {"FV0A", " ; FV0A (%)", {HistType::kTH1F, {{500, 0, 1e3}}}},                                                                                         //
      {"FT0A", " ; FT0A (%)", {HistType::kTH1F, {{500, 0, 1e3}}}},                                                                                         //
      {"FT0C", " ; FT0C (%)", {HistType::kTH1F, {{500, 0, 1e3}}}}                                                                                          //
    }};

  std::vector<int> usedTracksIds;

  void init(InitContext&)
  {
    ccdb->setURL(url.value);
    // Enabling object caching, otherwise each call goes to the CCDB server
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    // Not later than now, will be replaced by the value of the train creation
    // This avoids that users can replace objects **while** a train is running
    ccdb->setCreatedNotAfter(nolaterthan.value);
    LOGF(info, "Getting object %s", path.value.data());
  }

  using FullBCs = soa::Join<aod::BCsWithTimestamps, aod::BcSels>;
  void processEventStat(
    FullBCs const& bcs,
    soa::Join<aod::Collisions, aod::EvSels> const& collisions)
  {
    std::vector<typename std::decay_t<decltype(collisions)>::iterator> cols;
    for (auto& bc : bcs) {
      if (!useEvSel || (bc.selection()[evsel::kIsBBT0A] &
                        bc.selection()[evsel::kIsBBT0C]) != 0) {
        registry.fill(HIST("Events/Selection"), 5.);
        cols.clear();
        for (auto& collision : collisions) {
          if (collision.has_foundBC()) {
            if (collision.foundBCId() == bc.globalIndex()) {
              cols.emplace_back(collision);
            }
          } else if (collision.bcId() == bc.globalIndex()) {
            cols.emplace_back(collision);
          }
        }
        LOGP(debug, "BC {} has {} collisions", bc.globalBC(), cols.size());
        if (!cols.empty()) {
          registry.fill(HIST("Events/Selection"), 6.);
          if (cols.size() > 1) {
            registry.fill(HIST("Events/Selection"), 7.);
          }
        }
      }
    }
  }

  PROCESS_SWITCH(MultiplicityCounter, processEventStat, "Collect event sample stats", false);

  expressions::Filter trackSelectionProper = ((aod::track::trackCutFlag & trackSelectionITS) == trackSelectionITS) &&
                                             ifnode((aod::track::detectorMap & (uint8_t)o2::aod::track::TPC) == (uint8_t)o2::aod::track::TPC,
                                                    (aod::track::trackCutFlag & trackSelectionTPC) == trackSelectionTPC,
                                                    true) &&
                                             ((aod::track::trackCutFlag & trackSelectionDCA) == trackSelectionDCA);

  using ExTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA>;
  using FiTracks = soa::Filtered<ExTracks>;

  expressions::Filter atrackFilter = (aod::track::bestCollisionId >= 0) &&
                                     (nabs(aod::track::bestDCAZ) <= 2.f) &&
                                     (nabs(aod::track::bestDCAXY) <= ((0.0105f + 0.0350f / npow(aod::track::pts, 1.1f))));

  std::vector<Double_t> tracketas;
  void processCounting(
    soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision,
    BCsRun3 const& bcs,
    aod::FT0s const& ft0s,
    aod::FV0As const& fv0as,
    FiTracks const& tracks,
    soa::SmallGroups<aod::ReassignedTracksCore> const& atracks) // soa::Join<aod::AmbiguousTracks, aod::BestCollisions>
  {

    const auto& foundBC = collision.foundBC_as<BCsRun3>();

    float multT0A = 0;
    float multT0C = 0;
    if (foundBC.has_ft0()) {
      for (auto amplitude : foundBC.ft0().amplitudeA()) {
        multT0A += amplitude;
      }
      for (auto amplitude : foundBC.ft0().amplitudeC()) {
        multT0C += amplitude;
      }
    } else {
      multT0A = multT0C = -999;
    }
    float multV0A = 0;
    if (foundBC.has_fv0a()) {
      for (auto amplitude : foundBC.fv0a().amplitude()) {
        multV0A += amplitude;
      }
    } else {
      multV0A = -999;
    }

    registry.fill(HIST("FT0A"), multT0A);
    registry.fill(HIST("FT0C"), multT0C);
    registry.fill(HIST("FV0A"), multV0A);
    registry.fill(HIST("Multiplicity"), multV0A, multT0A, multT0C);

    registry.fill(HIST("Events/Selection"), 1.);
    if (!useEvSel || collision.sel8()) {
      registry.fill(HIST("Events/Selection"), 2.);
      auto z = collision.posZ();

      registry.fill(HIST("hreczvtx"), Double_t(kDATA), Double_t(kMBAND), 50., z);
      usedTracksIds.clear();

      tracketas.clear();
      for (auto& track : atracks) {
        auto otrack = track.track_as<FiTracks>();
        // tracketas.push_back(track.etas());
        tracketas.push_back(otrack.eta());
      }
      for (auto& track : tracks) {
        if (std::find(usedTracksIds.begin(), usedTracksIds.end(), track.globalIndex()) != usedTracksIds.end()) {
          continue;
        }
        if (z > -10 && z < 10) {
          registry.fill(HIST("PhiEta"), track.phi(), track.eta());
          registry.fill(HIST("DCAXY"), track.dcaXY());
          registry.fill(HIST("DCAZ"), track.dcaZ());
        }
        tracketas.push_back(track.eta());
      }

      for (auto eta : tracketas) {
        registry.fill(HIST("hrecdndeta"), Double_t(kDATA), Double_t(kMBAND), 50., z, eta);
        // registry.fill(HIST("hrecdndeta"), Double_t(kINEL), Double_t(kMBAND), 50., z, 1);
      }
    }
  }
  PROCESS_SWITCH(MultiplicityCounter, processCounting, "Count tracks", false);

  using Particles = soa::Filtered<aod::McParticles>;
  using LabeledTracksEx = soa::Join<LabeledTracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA>;
  using Particle = Particles::iterator;
  using ParticlesI = soa::Join<aod::McParticles, aod::ParticlesToTracks>;
  expressions::Filter primaries = (aod::mcparticle::flags & (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary) == (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary;
  Partition<Particles> mcSample = nabs(aod::mcparticle::eta) < estimatorEta;
  Partition<ParticlesI> primariesI = ((aod::mcparticle::flags & (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary) == (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary) &&
                                     (nabs(aod::mcparticle::eta) < estimatorEta);

  Preslice<FiTracks> perCol = aod::track::collisionId;

  Partition<soa::Filtered<LabeledTracksEx>> lsample = nabs(aod::track::eta) < estimatorEta;
  void processMCCounting(
    soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels> const& collisions,
    aod::McCollisions const&, Particles const& mcParticles,
    soa::Filtered<LabeledTracksEx> const&,
    soa::SmallGroups<aod::ReassignedTracksCore> const& atracks)
  {
    for (auto& collision : collisions) {
      if (useEvSel && !collision.sel8()) {
        continue;
      }
      if (!collision.has_mcCollision()) {
        continue;
      }
      auto z = collision.posZ();
      registry.fill(HIST("hreczvtx"), Double_t(kINEL), Double_t(kMBAND), 50., z);
      auto mcCollision = collision.mcCollision();
      auto particles = mcSample->sliceByCached(aod::mcparticle::mcCollisionId, mcCollision.globalIndex());
      auto tracks = lsample->sliceByCached(aod::track::collisionId, collision.globalIndex());
      tracks.bindExternalIndices(&mcParticles);

      usedTracksIds.clear();
      for (auto& track : atracks) {
        auto ttrack = track.track_as<soa::Filtered<LabeledTracksEx>>();
        usedTracksIds.emplace_back(ttrack.globalIndex());
        if (ttrack.has_mcParticle()) {
          registry.fill(HIST("hrecdndeta"), Double_t(kINEL), Double_t(kMBAND), 50., z, ttrack.mcParticle_as<Particles>().eta());
        } else {
          // when secondary
        }
      }
      for (auto& track : tracks) {
        if (std::find(usedTracksIds.begin(), usedTracksIds.end(), track.globalIndex()) != usedTracksIds.end()) {
          continue;
        }
        if (track.has_mcParticle()) {
          registry.fill(HIST("hrecdndeta"), Double_t(kINEL), Double_t(kMBAND), 50., z, track.mcParticle_as<Particles>().eta());
        } else {
          // when secondary
        }
      }
    }
  }
  PROCESS_SWITCH(MultiplicityCounter, processMCCounting, "MC Count tracks", false);

  void processGen(
    aod::McCollisions::iterator const& mcCollision,
    o2::soa::SmallGroups<soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>> const& collisions,
    Particles const& particles, FiTracks const& tracks)
  {
    auto perCollisionMCSample = mcSample->sliceByCached(aod::mcparticle::mcCollisionId, mcCollision.globalIndex());
    auto genz = mcCollision.posZ();
    registry.fill(HIST("hgenzvtx"), Double_t(kINEL), 50., genz);
    for (auto& particle : perCollisionMCSample) {
      auto p = pdg->GetParticle(particle.pdgCode());
      if (p != nullptr) {
        if (std::abs(p->Charge()) >= 3) {
          registry.fill(HIST("hgendndeta"), Double_t(kINEL), 50., genz, particle.eta());
        }
      }
    }
  }
  PROCESS_SWITCH(MultiplicityCounter, processGen, "Process generator-level info", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<MultiplicityCounter>(cfgc)};
}
