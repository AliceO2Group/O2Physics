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
#include "Framework/Configurable.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/RuntimeError.h"
#include "Framework/runDataProcessing.h"

#include "ReconstructionDataFormats/GlobalTrackID.h"
#include "Common/CCDB/EventSelectionParams.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "CommonConstants/MathConstants.h"
#include "TDatabasePDG.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct PseudorapidityDensity {
  Service<TDatabasePDG> pdg;

  Configurable<float> estimatorEta{"estimatorEta", 1.0, "eta range for INEL>0 sample definition"};

  Configurable<bool> useEvSel{"useEvSel", true, "use event selection"};
  Configurable<bool> useDCAZ{"useDCAZ", true, "use DCAZ cut"};
  Configurable<bool> useDCAXY{"useDCAXY", false, "use DCAXY cut"};
  Configurable<bool> usePtDCAXY{"usePtDCAXY", false, "use pt-dependent DCAXY"};
  Configurable<float> maxDCAXY{"maxDCAXY", 2.4, "max allowed transverse DCA"};
  Configurable<float> maxDCAZ{"maxDCAZ", 3.2, "max allowed longitudal DCA"};

  HistogramRegistry registry{
    "registry",
    {
      {"EventsNtrkZvtx", "; N_{trk}; Z_{vtx}; events", {HistType::kTH2F, {{301, -0.5, 300.5}, {201, -20.1, 20.1}}}}, //
      {"TracksEtaZvtx", "; #eta; Z_{vtx}; tracks", {HistType::kTH2F, {{21, -2.1, 2.1}, {201, -20.1, 20.1}}}},        //
      {"TracksEtaZvtx_gt0", "; #eta; Z_{vtx}; tracks", {HistType::kTH2F, {{21, -2.1, 2.1}, {201, -20.1, 20.1}}}},    //
      {"TracksPhiEta", "; #varphi; #eta; tracks", {HistType::kTH2F, {{600, 0, 2 * M_PI}, {21, -2.1, 2.1}}}},         //
      {"EventSelection", ";status;events", {HistType::kTH1F, {{7, 0.5, 7.5}}}}                                       //
    }                                                                                                                //
  };

  void init(InitContext&)
  {
    auto hstat = registry.get<TH1>(HIST("EventSelection"));
    auto* x = hstat->GetXaxis();
    x->SetBinLabel(1, "All");
    x->SetBinLabel(2, "Selected");
    x->SetBinLabel(3, "Selected INEL>0");
    x->SetBinLabel(4, "Rejected");
    x->SetBinLabel(5, "Good BCs");
    x->SetBinLabel(6, "BCs with collisions");
    x->SetBinLabel(7, "BCs with pile-up/splitting");

    if (doprocessGen) {
      registry.add({"EventsNtrkZvtxGen", "; N_{trk}; Z_{vtx}; events", {HistType::kTH2F, {{301, -0.5, 300.5}, {201, -20.1, 20.1}}}});
      registry.add({"EventsNtrkZvtxGen_t", "; N_{part}; Z_{vtx}; events", {HistType::kTH2F, {{301, -0.5, 300.5}, {201, -20.1, 20.1}}}});
      registry.add({"TracksEtaZvtxGen", "; #eta; Z_{vtx}; tracks", {HistType::kTH2F, {{21, -2.1, 2.1}, {201, -20.1, 20.1}}}});
      registry.add({"TracksEtaZvtxGen_t", "; #eta; Z_{vtx}; tracks", {HistType::kTH2F, {{21, -2.1, 2.1}, {201, -20.1, 20.1}}}});
      registry.add({"TracksEtaZvtxGen_gt0", "; #eta; Z_{vtx}; tracks", {HistType::kTH2F, {{21, -2.1, 2.1}, {201, -20.1, 20.1}}}});
      registry.add({"TracksEtaZvtxGen_gt0t", "; #eta; Z_{vtx}; tracks", {HistType::kTH2F, {{21, -2.1, 2.1}, {201, -20.1, 20.1}}}});

      registry.add({"TracksPhiEtaGen", "; #varphi; #eta; tracks", {HistType::kTH2F, {{600, 0, 2 * M_PI}, {21, -2.1, 2.1}}}});
      registry.add({"EventEfficiency", "; status; events", {HistType::kTH1F, {{5, 0.5, 5.5}}}});
      registry.add({"NotFoundEventZvtx", " ; Z_{vtx}", {HistType::kTH1F, {{201, -20.1, 20.1}}}});

      auto heff = registry.get<TH1>(HIST("EventEfficiency"));
      x = heff->GetXaxis();
      x->SetBinLabel(1, "Generated");
      x->SetBinLabel(2, "Generated INEL>0");
      x->SetBinLabel(3, "Reconstructed");
      x->SetBinLabel(4, "Selected");
      x->SetBinLabel(5, "Selected INEL>0");
    }
  }

  using FullBCs = soa::Join<aod::BCsWithTimestamps, aod::BcSels>;
  void processTagging(FullBCs const& bcs, soa::Join<aod::Collisions, aod::EvSels> const& collisions)
  {
    std::vector<typename std::decay_t<decltype(collisions)>::iterator> cols;
    for (auto& bc : bcs) {
      if (!useEvSel || (useEvSel && ((bc.selection()[evsel::kIsBBT0A] & bc.selection()[evsel::kIsBBT0C]) != 0))) {
        registry.fill(HIST("EventSelection"), 5.);
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
          registry.fill(HIST("EventSelection"), 6.);
          if (cols.size() > 1) {
            registry.fill(HIST("EventSelection"), 7.);
          }
        }
      }
    }
  }

  PROCESS_SWITCH(PseudorapidityDensity, processTagging, "Collect event sample stats", false);

  expressions::Filter ITStracks = (aod::track::detectorMap & (uint8_t)o2::aod::track::ITS) != (uint8_t)0;
  expressions::Filter DCAFilterZ = ifnode(useDCAZ.node(), nabs(aod::track::dcaZ) <= maxDCAZ, framework::expressions::LiteralNode{true});
  expressions::Filter DCAFilterXY = ifnode(useDCAXY.node(), ifnode(usePtDCAXY.node(), nabs(aod::track::dcaXY) <= 0.0105f + 0.0350f / npow(aod::track::pt, 1.1f), nabs(aod::track::dcaXY) <= maxDCAXY), framework::expressions::LiteralNode{true});

  using Trks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksExtended>>;
  Partition<Trks> sample = nabs(aod::track::eta) < estimatorEta;

  void process(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, Trks const& tracks)
  {
    registry.fill(HIST("EventSelection"), 1.);
    if (!useEvSel || (useEvSel && collision.sel8())) {
      registry.fill(HIST("EventSelection"), 2.);
      auto z = collision.posZ();
      auto perCollisionSample = sample->sliceByCached(aod::track::collisionId, collision.globalIndex());
      if (perCollisionSample.size() > 0) {
        registry.fill(HIST("EventSelection"), 3.);
      }
      registry.fill(HIST("EventsNtrkZvtx"), perCollisionSample.size(), z);
      for (auto& track : tracks) {
        registry.fill(HIST("TracksEtaZvtx"), track.eta(), z);
        registry.fill(HIST("TracksPhiEta"), track.phi(), track.eta());
        if (perCollisionSample.size() > 0) {
          registry.fill(HIST("TracksEtaZvtx_gt0"), track.eta(), z);
        }
      }
    } else {
      registry.fill(HIST("EventSelection"), 4.);
    }
  }

  using Particles = aod::McParticles;
  expressions::Filter primaries = (aod::mcparticle::flags & (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary) == (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary;
  Partition<soa::Filtered<Particles>> mcSample = nabs(aod::mcparticle::eta) < estimatorEta;

  void processGen(aod::McCollisions::iterator const& mcCollision, o2::soa::SmallGroups<soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>> const& collisions, soa::Filtered<Particles> const& particles, Trks const& /*tracks*/)
  {
    auto perCollisionMCSample = mcSample->sliceByCached(aod::mcparticle::mcCollisionId, mcCollision.globalIndex());
    auto nCharged = 0;
    for (auto& particle : perCollisionMCSample) {
      auto charge = 0;
      auto p = pdg->GetParticle(particle.pdgCode());
      if (p != nullptr) {
        charge = (int)p->Charge();
      }
      if (charge != 0) {
        nCharged++;
      }
    }
    registry.fill(HIST("EventsNtrkZvtxGen_t"), nCharged, mcCollision.posZ());
    registry.fill(HIST("EventEfficiency"), 1.);

    if (nCharged > 0) {
      registry.fill(HIST("EventEfficiency"), 2.);
    }
    bool atLeastOne = false;
    bool atLeastOne_gt0 = false;
    LOGP(debug, "MC col {} has {} reco cols", mcCollision.globalIndex(), collisions.size());
    for (auto& collision : collisions) {
      registry.fill(HIST("EventEfficiency"), 3.);
      if (!useEvSel || (useEvSel && collision.sel8())) {
        atLeastOne = true;
        auto perCollisionSample = sample->sliceByCached(aod::track::collisionId, collision.globalIndex());
        registry.fill(HIST("EventEfficiency"), 4.);
        if (perCollisionSample.size() > 0) {
          atLeastOne_gt0 = true;
          registry.fill(HIST("EventEfficiency"), 5.);
        }
        registry.fill(HIST("EventsNtrkZvtxGen"), perCollisionSample.size(), collision.posZ());
      }
    }
    if (collisions.size() == 0) {
      registry.fill(HIST("NotFoundEventZvtx"), mcCollision.posZ());
    }
    for (auto& particle : particles) {
      auto p = pdg->GetParticle(particle.pdgCode());
      auto charge = 0;
      if (p != nullptr) {
        charge = (int)p->Charge();
      }
      if (charge != 0) {
        registry.fill(HIST("TracksEtaZvtxGen_t"), particle.eta(), mcCollision.posZ());
        if (perCollisionMCSample.size() > 0) {
          registry.fill(HIST("TracksEtaZvtxGen_gt0t"), particle.eta(), mcCollision.posZ());
        }
        if (atLeastOne) {
          registry.fill(HIST("TracksEtaZvtxGen"), particle.eta(), mcCollision.posZ());
          if (atLeastOne_gt0) {
            registry.fill(HIST("TracksEtaZvtxGen_gt0"), particle.eta(), mcCollision.posZ());
          }
        }
        registry.fill(HIST("TracksPhiEtaGen"), particle.phi(), particle.eta());
      }
    }
  }

  PROCESS_SWITCH(PseudorapidityDensity, processGen, "Process generator-level info", false);
};

WorkflowSpec
  defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<PseudorapidityDensity>(cfgc)};
}
