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

#include <Math/Vector4D.h>
#include <cmath>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TPDGCode.h>
#include <TDatabasePDG.h>

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/Configurable.h"
#include "Framework/RuntimeError.h"
#include "Framework/runDataProcessing.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "CCDB/BasicCCDBManager.h"

#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/CCDB/EventSelectionParams.h"
#include "CommonConstants/MathConstants.h"
#include "bestCollisionTable.h"
#include "ReconstructionDataFormats/GlobalTrackID.h"
#include "ReconstructionDataFormats/Track.h"
#include "Index.h"

#include "PWGLF/DataModel/LFStrangenessTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::track;

using BCsRun3 = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;
using MyCollisions = soa::Join<aod::Collisions, aod::EvSels>;
using MyCollisionsCent = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs>;
using FullBCs = soa::Join<aod::BCsWithTimestamps, aod::BcSels>;
using ExTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA>;
using FiTracks = soa::Filtered<ExTracks>;
using Particles = soa::Filtered<aod::McParticles>;
using Particle = Particles::iterator;
using LabeledTracks = soa::Join<aod::Tracks, aod::McTrackLabels>;
using LabeledTracksEx = soa::Join<LabeledTracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA>;
using MyTracks = soa::Join<aod::Tracks, aod::McTrackLabels, aod::TracksExtra, aod::pidTPCPr>;

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
enum {
  kSpeciesbegin = 0,
  kK0short = 1,
  kLambda,
  kAntilambda,
  kSpeciesend
};

AxisSpec ZAxis = {60, -30, 30, "", "zaxis"};
AxisSpec DeltaZAxis = {61, -6.1, 6.1, "", "deltaz axis"};
AxisSpec DCAAxis = {601, -3.01, 3.01, "", "DCA axis"};
AxisSpec EtaAxis = {80, -4.0, 4.0, "", "eta axis"};
AxisSpec PhiAxis = {629, 0, 2 * M_PI, "", "phi axis"};
AxisSpec PtAxis = {2401, -0.005, 24.005, "", "P_{T} axis"};
AxisSpec EvtClassAxis = {kECend - 1, kECbegin + 0.5, kECend - 0.5, "", "event class"};
AxisSpec TrigClassAxis = {kTrigend - 1, kTrigbegin + 0.5, kTrigend - 0.5, "", "trigger class"};
std::vector<double> centBinning = {0, 10., 20., 30., 40., 50., 60., 70., 80., 100};
AxisSpec CentAxis = {centBinning, "", "centrality"};
AxisSpec SpeciesAxis = {kSpeciesend - 1, kSpeciesbegin + 0.5, kSpeciesend - 0.5, "", "species class"};
AxisSpec massAxis = {600, 0.0f, 3.0f, "", "Inv. Mass (GeV/c^{2})"};

static constexpr TrackSelectionFlags::flagtype trackSelectionITS =
  TrackSelectionFlags::kITSNCls | TrackSelectionFlags::kITSChi2NDF |
  TrackSelectionFlags::kITSHits;

static constexpr TrackSelectionFlags::flagtype trackSelectionTPC =
  TrackSelectionFlags::kTPCNCls |
  TrackSelectionFlags::kTPCCrossedRowsOverNCls |
  TrackSelectionFlags::kTPCChi2NDF;

static constexpr TrackSelectionFlags::flagtype trackSelectionDCA =
  TrackSelectionFlags::kDCAz | TrackSelectionFlags::kDCAxy;

struct MultiplicityCounter {
  Service<O2DatabasePDG> pdg;

  Configurable<float> estimatorEta{"estimatorEta", 1.0, "eta range for INEL>0 sample definition"};
  Configurable<bool> useEvSel{"useEvSel", true, "use event selection"};
  Configurable<bool> isMC{"isMC", false, "check if MC"};

  ConfigurableAxis multBinning{"multBinning", {301, -0.5, 300.5}, ""};
  AxisSpec MultAxis = {multBinning, "N"};
  Configurable<bool> CentOn{"CentOn", false, "centrality switch ; watch your ccdb timestamps"}; // on debuging..

  HistogramRegistry registry{
    "registry",
    {{"Events/Selection", ";status;events", {HistType::kTH1F, {{7, 0.5, 7.5}}}},                                                                          //
     {"hrecdndeta", "evntclass; triggerclass; centrality, zvtex, eta", {HistType::kTHnSparseD, {EvtClassAxis, TrigClassAxis, CentAxis, ZAxis, EtaAxis}}}, //
     {"hgendndeta", "evntclass; centrality, zvtex, eta", {HistType::kTHnSparseD, {EvtClassAxis, CentAxis, ZAxis, EtaAxis}}},                              //
     {"hreczvtx", "evntclass; triggerclass; centrality, zvtex", {HistType::kTHnSparseD, {EvtClassAxis, TrigClassAxis, CentAxis, ZAxis}}},                 //
     {"hgenzvtx", "evntclass; centrality, zvtex", {HistType::kTHnSparseD, {EvtClassAxis, CentAxis, ZAxis}}},                                              //
     {"PhiEta", "; #varphi; #eta; tracks", {HistType::kTHnSparseD, {EvtClassAxis, PhiAxis, EtaAxis}}},                                                    //
     {"DCAXY", " ; DCA_{XY} (cm)", {HistType::kTHnSparseD, {EvtClassAxis, DCAAxis}}},                                                                     //
     {"DCAZ", " ; DCA_{Z} (cm)", {HistType::kTHnSparseD, {EvtClassAxis, DCAAxis}}},                                                                       //
     {"Multiplicity", " ; FV0A (#); FT0A (#); FT0C (#) ", {HistType::kTHnSparseD, {MultAxis, MultAxis, MultAxis}}},                                       //
     {"Centrality", " ; centrality_FT0C (%) ", {HistType::kTH1F, {CentAxis}}},                                                                            //
     {"Centrality_MBAND", " ; centrality_MBAND_FT0C (%) ", {HistType::kTH1F, {CentAxis}}},                                                                //
     {"hrecpt", " eventclass; centrality; pt_gen; pt_rec ", {HistType::kTHnSparseD, {EvtClassAxis, CentAxis, PtAxis, PtAxis}}},                           //
     {"hgenpt", " eventclass; centrality; pt;  ", {HistType::kTHnSparseD, {EvtClassAxis, CentAxis, PtAxis}}},                                             //
     {"hV0Mass", "species ; evntclass; K0shortMass; LambdaMass; AntiLambdaMass", {HistType::kTHnSparseD, {SpeciesAxis, EvtClassAxis, massAxis, massAxis, massAxis}}}}};

  std::vector<int> usedTracksIds;

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
  expressions::Filter atrackFilter = (aod::track::bestCollisionId >= 0) &&
                                     (nabs(aod::track::bestDCAZ) <= 2.f) &&
                                     (nabs(aod::track::bestDCAXY) <= ((0.0105f + 0.0350f / npow(aod::track::pts, 1.1f))));
  std::vector<Double_t> tracketas;

  void processCounting(
    MyCollisionsCent::iterator const& collision,
    BCsRun3 const& bcs,
    aod::FT0s const& ft0s,
    aod::FV0As const& fv0as,
    FiTracks const& tracks,
    aod::V0Datas const& fullV0s,
    soa::SmallGroups<aod::ReassignedTracksCore> const& atracks)
  {
    const auto& foundBC = collision.foundBC_as<BCsRun3>();
    float multT0A = 0;
    float multT0C = 0;
    float multV0A = 0;
    auto cent = collision.centFT0C();

    registry.fill(HIST("Centrality"), cent);
    registry.fill(HIST("Events/Selection"), 1.);
    auto z = collision.posZ();

    if (!useEvSel || collision.sel8()) {
      if (std::abs(z) < 10) {
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

        if (foundBC.has_fv0a()) {
          for (auto amplitude : foundBC.fv0a().amplitude()) {
            multV0A += amplitude;
          }
        } else {
          multV0A = -999;
        }
        registry.fill(HIST("Multiplicity"), multV0A, multT0A, multT0C);
        registry.fill(HIST("Centrality_MBAND"), cent);

        registry.fill(HIST("Events/Selection"), 2.);

        registry.fill(HIST("hreczvtx"), Double_t(kDATA), Double_t(kMBAND), cent, z);
        usedTracksIds.clear();

        tracketas.clear();
        for (auto& track : atracks) {
          auto otrack = track.track_as<FiTracks>();
          tracketas.push_back(otrack.eta());
          registry.fill(HIST("PhiEta"), Double_t(kDATA), otrack.phi(), otrack.eta());
          registry.fill(HIST("DCAXY"), Double_t(kDATA), otrack.dcaXY());
          registry.fill(HIST("DCAZ"), Double_t(kDATA), otrack.dcaZ());
          registry.fill(HIST("hrecpt"), Double_t(kDATA), cent, -1, otrack.pt());
        }
        for (auto& track : tracks) {
          if (std::find(usedTracksIds.begin(), usedTracksIds.end(), track.globalIndex()) != usedTracksIds.end()) {
            continue;
          }
          registry.fill(HIST("PhiEta"), Double_t(kDATA), track.phi(), track.eta());
          registry.fill(HIST("DCAXY"), Double_t(kDATA), track.dcaXY());
          registry.fill(HIST("DCAZ"), Double_t(kDATA), track.dcaZ());
          registry.fill(HIST("hrecpt"), Double_t(kDATA), cent, -1, track.pt());
          tracketas.push_back(track.eta());
        }

        for (auto eta : tracketas) {
          registry.fill(HIST("hrecdndeta"), Double_t(kDATA), Double_t(kMBAND), cent, z, eta);
        }
        for (auto& v0 : fullV0s) {
          registry.fill(HIST("hV0Mass"), Double_t(kDATA), Double_t(kK0short), v0.mK0Short(), -1, -1);
          registry.fill(HIST("hV0Mass"), Double_t(kDATA), Double_t(kLambda), -1, v0.mLambda(), -1);
          registry.fill(HIST("hV0Mass"), Double_t(kDATA), Double_t(kAntilambda), -1, -1, v0.mAntiLambda());
        }
      }
    }
  }
  PROCESS_SWITCH(MultiplicityCounter, processCounting, "Count tracks", false);

  void processCountingPP(
    MyCollisions::iterator const& collision,
    BCsRun3 const& bcs,
    aod::FT0s const& ft0s,
    aod::FV0As const& fv0as,
    FiTracks const& tracks,
    aod::V0Datas const& fullV0s,
    soa::SmallGroups<aod::ReassignedTracksCore> const& atracks)
  {
    const auto& foundBC = collision.foundBC_as<BCsRun3>();
    float multT0A = 0;
    float multT0C = 0;
    float multV0A = 0;
    auto cent = 50.;

    registry.fill(HIST("Centrality"), cent);
    registry.fill(HIST("Events/Selection"), 1.);
    auto z = collision.posZ();

    if (!useEvSel || collision.sel8()) {
      if (std::abs(z) < 10) {
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

        if (foundBC.has_fv0a()) {
          for (auto amplitude : foundBC.fv0a().amplitude()) {
            multV0A += amplitude;
          }
        } else {
          multV0A = -999;
        }
        registry.fill(HIST("Multiplicity"), multV0A, multT0A, multT0C);
        registry.fill(HIST("Centrality_MBAND"), cent);

        registry.fill(HIST("Events/Selection"), 2.);

        registry.fill(HIST("hreczvtx"), Double_t(kDATA), Double_t(kMBAND), cent, z);
        usedTracksIds.clear();

        tracketas.clear();
        for (auto& track : atracks) {
          auto otrack = track.track_as<FiTracks>();
          tracketas.push_back(otrack.eta());
          registry.fill(HIST("PhiEta"), Double_t(kDATA), otrack.phi(), otrack.eta());
          registry.fill(HIST("DCAXY"), Double_t(kDATA), otrack.dcaXY());
          registry.fill(HIST("DCAZ"), Double_t(kDATA), otrack.dcaZ());
          registry.fill(HIST("hrecpt"), Double_t(kDATA), cent, -1, otrack.pt());
        }
        for (auto& track : tracks) {
          if (std::find(usedTracksIds.begin(), usedTracksIds.end(), track.globalIndex()) != usedTracksIds.end()) {
            continue;
          }
          registry.fill(HIST("PhiEta"), Double_t(kDATA), track.phi(), track.eta());
          registry.fill(HIST("DCAXY"), Double_t(kDATA), track.dcaXY());
          registry.fill(HIST("DCAZ"), Double_t(kDATA), track.dcaZ());
          registry.fill(HIST("hrecpt"), Double_t(kDATA), cent, -1, track.pt());
          tracketas.push_back(track.eta());
        }

        for (auto eta : tracketas) {
          registry.fill(HIST("hrecdndeta"), Double_t(kDATA), Double_t(kMBAND), cent, z, eta);
        }
        for (auto& v0 : fullV0s) {
          registry.fill(HIST("hV0Mass"), Double_t(kDATA), Double_t(kK0short), v0.mK0Short(), -1, -1);
          registry.fill(HIST("hV0Mass"), Double_t(kDATA), Double_t(kLambda), -1, v0.mLambda(), -1);
          registry.fill(HIST("hV0Mass"), Double_t(kDATA), Double_t(kAntilambda), -1, -1, v0.mAntiLambda());
        }
      }
    }
  }
  PROCESS_SWITCH(MultiplicityCounter, processCountingPP, "Count tracks", false);

  expressions::Filter primaries = (aod::mcparticle::flags & (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary) == (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary;
  Partition<Particles> mcSample = nabs(aod::mcparticle::eta) < estimatorEta;
  Preslice<FiTracks> perCol = aod::track::collisionId;
  Partition<soa::Filtered<LabeledTracksEx>> lsample = nabs(aod::track::eta) < estimatorEta;

  void processMCCounting(
    soa::Join<MyCollisions, aod::McCollisionLabels> const& collisions,
    aod::McCollisions const&,
    aod::V0Datas const& fullV0s,
    Particles const& mcParticles,
    soa::Filtered<LabeledTracksEx> const&,
    soa::SmallGroups<aod::ReassignedTracksCore> const& atracks)
  {
    auto cent = 50.;

    for (auto& collision : collisions) {

      auto z = collision.posZ();
      if (useEvSel && !collision.sel8()) {
        if (std::abs(z) > 10) {
          continue;
        }
      }
      if (!collision.has_mcCollision()) {
        continue;
      }
      for (auto& v0 : fullV0s) {
        registry.fill(HIST("hV0Mass"), Double_t(kINEL), Double_t(kK0short), v0.mK0Short(), -1, -1);
        registry.fill(HIST("hV0Mass"), Double_t(kINEL), Double_t(kLambda), -1, v0.mLambda(), -1);
        registry.fill(HIST("hV0Mass"), Double_t(kINEL), Double_t(kAntilambda), -1, -1, v0.mAntiLambda());
      }

      registry.fill(HIST("hreczvtx"), Double_t(kINEL), Double_t(kMBAND), cent, z);
      auto mcCollision = collision.mcCollision();
      auto particles = mcSample->sliceByCached(aod::mcparticle::mcCollisionId, mcCollision.globalIndex());
      auto tracks = lsample->sliceByCached(aod::track::collisionId, collision.globalIndex());
      tracks.bindExternalIndices(&mcParticles);

      usedTracksIds.clear();
      for (auto& track : atracks) {
        auto ttrack = track.track_as<soa::Filtered<LabeledTracksEx>>();
        usedTracksIds.emplace_back(ttrack.globalIndex());
        if (ttrack.has_mcParticle()) {
          registry.fill(HIST("hrecdndeta"), Double_t(kINEL), Double_t(kMBAND), cent, z, ttrack.mcParticle_as<Particles>().eta());
          registry.fill(HIST("hrecpt"), Double_t(kINEL), cent, ttrack.mcParticle_as<Particles>().pt(), ttrack.pt());

          registry.fill(HIST("PhiEta"), Double_t(kINEL), ttrack.phi(), ttrack.eta());
          registry.fill(HIST("DCAXY"), Double_t(kINEL), ttrack.dcaXY());
          registry.fill(HIST("DCAZ"), Double_t(kINEL), ttrack.dcaZ());
        } else {
          // when secondary
        }
      }
      for (auto& track : tracks) {
        if (std::find(usedTracksIds.begin(), usedTracksIds.end(), track.globalIndex()) != usedTracksIds.end()) {
          continue;
        }
        if (track.has_mcParticle()) {
          registry.fill(HIST("hrecdndeta"), Double_t(kINEL), Double_t(kMBAND), cent, z, track.mcParticle_as<Particles>().eta());
          registry.fill(HIST("hrecpt"), Double_t(kINEL), cent, track.mcParticle_as<Particles>().pt(), track.pt());

          registry.fill(HIST("PhiEta"), Double_t(kINEL), track.phi(), track.eta());
          registry.fill(HIST("DCAXY"), Double_t(kINEL), track.dcaXY());
          registry.fill(HIST("DCAZ"), Double_t(kINEL), track.dcaZ());
        } else {
          // when secondary
        }
      }
      for (auto& particle : particles) {
        auto p = pdg->GetParticle(particle.pdgCode());
        if (p != nullptr) {
          if (std::abs(p->Charge()) >= 3) {
            registry.fill(HIST("hgenpt"), Double_t(kINEL), cent, particle.pt());
          }
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
    auto cent = 50.;
    auto perCollisionMCSample = mcSample->sliceByCached(aod::mcparticle::mcCollisionId, mcCollision.globalIndex());
    auto genz = mcCollision.posZ();
    registry.fill(HIST("hgenzvtx"), Double_t(kINEL), cent, genz);
    for (auto& particle : perCollisionMCSample) {
      auto p = pdg->GetParticle(particle.pdgCode());
      if (p != nullptr) {
        if (std::abs(p->Charge()) >= 3) {
          registry.fill(HIST("hgendndeta"), Double_t(kINEL), cent, genz, particle.eta());
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
