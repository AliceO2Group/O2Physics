// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
#ifndef DNDETA_H
#define DNDETA_H
#include "Framework/Configurable.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"

#include "ReconstructionDataFormats/GlobalTrackID.h"

#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/Core/MC.h"

#include "TDatabasePDG.h"

namespace o2::pwgmm::multiplicity
{
namespace
{
template <typename C, uint8_t TRACKTYPE>
inline auto collisionSelector(C const& collision)
{
  if constexpr (TRACKTYPE == aod::track::TrackTypeEnum::Run2Tracklet) {
    return collision.sel7();
  } else if constexpr (TRACKTYPE == o2::dataformats::GlobalTrackID::ITS) {
    return collision.sel8();
  } else {
    return false;
  }
}
} // namespace

using namespace o2::framework;
template <uint8_t TRACKTYPE>
struct PseudorapidityDensity {
  o2::framework::Service<TDatabasePDG> pdg;

  o2::framework::Configurable<float> etaMax{"etaMax", 2.0, "max eta value"};
  o2::framework::Configurable<float> etaMin{"etaMin", -2.0, "min eta value"};
  o2::framework::Configurable<float> vtxZMax{"vtxZMax", 15, "max z vertex"};
  o2::framework::Configurable<float> vtxZMin{"vtxZMin", -15, "min z vertex"};

  o2::framework::HistogramRegistry registry{
    "registry",
    {
      {"EventsNtrkZvtx", "; N_{trk}; Z_{vtx}; events", {o2::framework::HistType::kTH2F, {{301, -0.5, 300.5}, {201, -20.1, 20.1}}}},    //
      {"TracksEtaZvtx", "; #eta; Z_{vtx}; tracks", {o2::framework::HistType::kTH2F, {{21, -2.1, 2.1}, {201, -20.1, 20.1}}}},           //
      {"TracksPhiEta", "; #varphi; #eta; tracks", {o2::framework::HistType::kTH2F, {{600, 0, 2 * M_PI}, {21, -2.1, 2.1}}}},            //
      {"EventSelection", ";status;events", {o2::framework::HistType::kTH1F, {{3, 0.5, 3.5}}}},                                         //
      {"EventsNtrkZvtxGen", "; N_{trk}; Z_{vtx}; events", {o2::framework::HistType::kTH2F, {{301, -0.5, 300.5}, {201, -20.1, 20.1}}}}, //
      {"TracksEtaZvtxGen", "; #eta; Z_{vtx}; tracks", {o2::framework::HistType::kTH2F, {{21, -2.1, 2.1}, {201, -20.1, 20.1}}}},        //
      {"TracksPhiEtaGen", "; #varphi; #eta; tracks", {o2::framework::HistType::kTH2F, {{600, 0, 2 * M_PI}, {21, -2.1, 2.1}}}},         //
      {"EventEfficiency", "; status; events", {o2::framework::HistType::kTH1F, {{3, 0.5, 3.5}}}}                                       //
    }                                                                                                                                  //
  };

  void init(o2::framework::InitContext&)
  {
    auto hstat = registry.get<TH1>(HIST("EventSelection"));
    auto* x = hstat->GetXaxis();
    x->SetBinLabel(1, "All");
    x->SetBinLabel(2, "Selected");
    x->SetBinLabel(3, "Rejected");

    auto heff = registry.get<TH1>(HIST("EventEfficiency"));
    x = heff->GetXaxis();
    x->SetBinLabel(1, "Generated");
    x->SetBinLabel(2, "Reconstructed");
    x->SetBinLabel(3, "Selected");
  }

  template <typename C>
  bool select(C const& collision)
  {
    return collisionSelector<C, TRACKTYPE>(collision);
  }

  o2::framework::expressions::Filter etaFilter = (aod::track::eta < etaMax) && (aod::track::eta > etaMin);
  o2::framework::expressions::Filter trackTypeFilter = (aod::track::trackType == TRACKTYPE);
  o2::framework::expressions::Filter posZFilter = (aod::collision::posZ < vtxZMax) && (aod::collision::posZ > vtxZMin);
  o2::framework::expressions::Filter posZFilterMC = (aod::mccollision::posZ < vtxZMax) && (aod::mccollision::posZ > vtxZMin);

  void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& collision, soa::Filtered<aod::Tracks> const& tracks)
  {
    registry.fill(HIST("EventSelection"), 1.);
    if (select(collision)) {
      registry.fill(HIST("EventSelection"), 2.);
      auto z = collision.posZ();
      registry.fill(HIST("EventsNtrkZvtx"), tracks.size(), z);
      for (auto& track : tracks) {
        registry.fill(HIST("TracksEtaZvtx"), track.eta(), z);
        registry.fill(HIST("TracksPhiEta"), track.phi(), track.eta());
      }
    } else {
      registry.fill(HIST("EventSelection"), 3.);
    }
  }

  using Particles = aod::McParticles;

  void processGen(soa::Filtered<aod::McCollisions>::iterator const& collision, Particles const& particles)
  {
    registry.fill(HIST("EventEfficiency"), 1.);
    for (auto& particle : particles) {
      auto p = pdg->GetParticle(particle.pdgCode());
      int charge = 0;
      if (p == nullptr) {
        LOGF(WARN, "[%d] Unknown particle with PDG code %d", particle.globalIndex(), particle.pdgCode());
      } else {
        charge = p->Charge();
      }
      if (charge != 0 && MC::isPhysicalPrimary(particle) && (particle.eta() < etaMax) && (particle.eta() > etaMin)) {
        registry.fill(HIST("TracksEtaZvtxGen"), particle.eta(), collision.posZ());
        registry.fill(HIST("TracksPhiEtaGen"), particle.phi(), particle.eta());
      }
    }
  }

  PROCESS_SWITCH(PseudorapidityDensity<TRACKTYPE>, processGen, "Process generator-level info", true);

  void processMatching(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>>::iterator const& collision, soa::Filtered<aod::Tracks> const& tracks, aod::McCollisions const&)
  {
    registry.fill(HIST("EventEfficiency"), 2.);
    if (select(collision)) {
      registry.fill(HIST("EventEfficiency"), 3.);
      auto z = collision.mcCollision().posZ();
      registry.fill(HIST("EventsNtrkZvtxGen"), tracks.size(), z);
    }
  }

  PROCESS_SWITCH(PseudorapidityDensity<TRACKTYPE>, processMatching, "Process generator-level info matched to reco", true);
};
} // namespace o2::pwgmm::multiplicity

#endif // DNDETA_H
