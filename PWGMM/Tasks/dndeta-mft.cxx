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
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "CommonConstants/MathConstants.h"
#include "TDatabasePDG.h"
#include "MathUtils/Utils.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct PseudorapidityDensityMFT {
  Service<TDatabasePDG> pdg;

  Configurable<bool> useEvSel{"useEvSel", true, "use event selection"};

  HistogramRegistry registry{
    "registry",
    {
      {"EventsNtrkZvtx", "; N_{trk}; Z_{vtx}; events", {HistType::kTH2F, {{301, -0.5, 300.5}, {201, -20.1, 20.1}}}}, //
      {"TracksEtaZvtx", "; #eta; Z_{vtx}; tracks", {HistType::kTH2F, {{18, -4.6, -1.}, {201, -20.1, 20.1}}}},        //
      {"TracksPhiEta", "; #varphi; #eta; tracks", {HistType::kTH2F, {{600, 0, 2 * M_PI}, {18, -4.6, -1.}}}},         //
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
      registry.add({"EventsNtrkZvtxGen_t", "; N_{trk}; Z_{vtx}; events", {HistType::kTH2F, {{301, -0.5, 300.5}, {201, -20.1, 20.1}}}});
      registry.add({"TracksEtaZvtxGen", "; #eta; Z_{vtx}; tracks", {HistType::kTH2F, {{18, -4.6, -1.}, {201, -20.1, 20.1}}}});
      registry.add({"TracksEtaZvtxGen_t", "; #eta; Z_{vtx}; tracks", {HistType::kTH2F, {{18, -4.6, -1.}, {201, -20.1, 20.1}}}});
      registry.add({"TracksPhiEtaGen", "; #varphi; #eta; tracks", {HistType::kTH2F, {{600, 0, 2 * M_PI}, {18, -4.6, -1.}}}});
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

  PROCESS_SWITCH(PseudorapidityDensityMFT, processTagging, "Collect event sample stats", true);

  void process(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, aod::MFTTracks const& tracks)
  {
    registry.fill(HIST("EventSelection"), 1.);
    if (!useEvSel || (useEvSel && collision.sel8())) {
      registry.fill(HIST("EventSelection"), 2.);
      auto z = collision.posZ();

      auto groupedTracks = tracks.sliceBy(o2::aod::fwdtrack::collisionId, collision.globalIndex());

      registry.fill(HIST("EventsNtrkZvtx"), groupedTracks.size(), z);
      for (auto& track : tracks) {
        registry.fill(HIST("TracksEtaZvtx"), track.eta(), z);
        float phi = track.phi();
        o2::math_utils::bringTo02Pi(phi);
        registry.fill(HIST("TracksPhiEta"), phi, track.eta());
      }
    } else {
      registry.fill(HIST("EventSelection"), 4.);
    }
  }

  using Particles = aod::McParticles;
  expressions::Filter primaries = (aod::mcparticle::flags & (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary) == (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary;

  void processGen(aod::McCollisions::iterator const& mcCollision, o2::soa::SmallGroups<soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>> const& collisions, soa::Filtered<Particles> const& particles, aod::MFTTracks const& tracks)
  {
    registry.fill(HIST("EventEfficiency"), 1.);

    bool atLeastOne = false;

    LOGP(debug, "MC col {} has {} reco cols", mcCollision.globalIndex(), collisions.size());
    for (auto& collision : collisions) {
      registry.fill(HIST("EventEfficiency"), 3.);
      if (!useEvSel || (useEvSel && collision.sel8())) {
        atLeastOne = true;
        auto groupedTracks = tracks.sliceBy(o2::aod::fwdtrack::collisionId, collision.globalIndex());

        registry.fill(HIST("EventEfficiency"), 4.);

        registry.fill(HIST("EventsNtrkZvtxGen"), groupedTracks.size(), collision.posZ());
      }
    }
    if (collisions.size() == 0) {
      registry.fill(HIST("NotFoundEventZvtx"), mcCollision.posZ());
    }
    auto nCharged = 0;
    for (auto& particle : particles) {
      auto p = pdg->GetParticle(particle.pdgCode());
      auto charge = 0;
      if (p != nullptr) {
        charge = (int)p->Charge();
      }
      if (std::abs(charge) < 3.) {
        continue;
      }

      registry.fill(HIST("TracksEtaZvtxGen_t"), particle.eta(), mcCollision.posZ());
      if (atLeastOne) {
        registry.fill(HIST("TracksEtaZvtxGen"), particle.eta(), mcCollision.posZ());
      }
      registry.fill(HIST("TracksPhiEtaGen"), particle.phi(), particle.eta());

      nCharged++;
    }
    registry.fill(HIST("EventsNtrkZvtxGen_t"), nCharged, mcCollision.posZ());
  }

  PROCESS_SWITCH(PseudorapidityDensityMFT, processGen, "Process generator-level info", false);
};

WorkflowSpec
  defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<PseudorapidityDensityMFT>(cfgc)};
}
