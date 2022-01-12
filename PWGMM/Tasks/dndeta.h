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
#ifndef DNDETA_H
#define DNDETA_H
#include "Framework/Configurable.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/RuntimeError.h"

#include "ReconstructionDataFormats/GlobalTrackID.h"

#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/MC.h"

#include "CommonConstants/MathConstants.h"

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
    throw framework::runtime_error("Unsupported track selection!");
  }
}

template <typename B, uint8_t TRACKTYPE>
inline auto BCSelector(B const& bc)
{
  if constexpr (TRACKTYPE == aod::track::TrackTypeEnum::Run2Tracklet) {
    return true;
  } else if constexpr (TRACKTYPE == o2::dataformats::GlobalTrackID::ITS) {
    return bc.selection()[aod::EventSelectionFlags::kIsBBT0A] & bc.selection()[aod::EventSelectionFlags::kIsBBT0C];
  } else {
    throw framework::runtime_error("Unsupported track selection!");
  }
}
} // namespace

using namespace o2::framework;

template <uint8_t TRACKTYPE>
struct PseudorapidityDensity {
  Service<TDatabasePDG> pdg;

  Configurable<float> estimatorEta{"estimatorEta", 1.0, "eta range for INEL>0 sample definition"};

  Configurable<bool> useEvSel{"useEvSel", true, "use event selection"};
  Configurable<bool> useDCA{"useDCA", false, "use DCA cuts"};
  Configurable<float> maxDCAXY{"maxDCAXY", 2.4, "max allowed transverse DCA"};
  Configurable<float> maxDCAZ{"maxDCAZ", 3.2, "max allowed longitudal DCA"};

  ConfigurableAxis percentileBinning{"pBins",
                                     {VARIABLE_WIDTH, 0., 0.01, 0.1, 0.5, 1, 5, 10, 15, 20, 30, 40, 50, 70, 100},
                                     "Centrality/multiplicity percentile binning"};

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

    if (doprocessBinned) {
      registry.add({"EventsNtrkZvtxBin", "; N_{trk}; Z_{vtx}; Percentile", {HistType::kTH3F, {{301, -0.5, 300.5}, {201, -20.1, 20.1}, percentileBinning}}});
      registry.add({"TracksEtaZvtxBin", "; #eta; Z_{vtx}; Percentile", {HistType::kTH3F, {{21, -2.1, 2.1}, {201, -20.1, 20.1}, percentileBinning}}});
      registry.add({"TracksPhiEtaBin", "; #varphi; #eta; Percentile", {HistType::kTH3F, {{600, 0, 2 * M_PI}, {21, -2.1, 2.1}, percentileBinning}}});
      registry.add({"EventSelectionBin", "; status; Percentile; events", {HistType::kTH2F, {{3, 0.5, 3.5}, percentileBinning}}});

      hstat = registry.get<TH2>(HIST("EventSelectionBin"));
      x = hstat->GetXaxis();
      x->SetBinLabel(1, "All");
      x->SetBinLabel(2, "Selected");
      x->SetBinLabel(3, "Rejected");
    }
  }

  template <typename C>
  bool select(C const& collision)
  {
    if (useEvSel) {
      return collisionSelector<C, TRACKTYPE>(collision);
    }
    return true;
  }

  template <typename B>
  bool selectBC(B const& bc)
  {
    return BCSelector<B, TRACKTYPE>(bc);
  }

  using FullBCs = soa::Join<aod::BCsWithTimestamps, aod::BcSels>;
  void processTagging(FullBCs const& bcs, soa::Join<aod::Collisions, aod::EvSels> const& collisions)
  {
    std::vector<typename std::decay_t<decltype(collisions)>::iterator> cols;
    for (auto& bc : bcs) {
      if (selectBC(bc)) {
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

  PROCESS_SWITCH(PseudorapidityDensity<TRACKTYPE>, processTagging, "Collect event sample stats", false);
  expressions::Filter trackTypeFilter = (aod::track::trackType == TRACKTYPE);
  expressions::Filter DCAFilter = ifnode(useDCA.node(), nabs(aod::track::dcaXY) <= maxDCAXY && nabs(aod::track::dcaZ) <= maxDCAZ, framework::expressions::LiteralNode{true});

  using Trks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtended>>;
  Partition<Trks> sample = nabs(aod::track::eta) < estimatorEta;

  void process(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, Trks const& tracks)
  {
    registry.fill(HIST("EventSelection"), 1.);
    if (select(collision)) {
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

  void processBinned(soa::Join<aod::Collisions, aod::EvSels, aod::CentV0Ms>::iterator const& collision, soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtended>> const& tracks)
  {
    auto p = collision.centV0M();
    registry.fill(HIST("EventSelectionBin"), 1., p);
    if (select(collision)) {
      registry.fill(HIST("EventSelectionBin"), 2., p);
      auto z = collision.posZ();
      auto localSample = sample->sliceByCached(aod::track::collisionId, collision.globalIndex());
      registry.fill(HIST("EventsNtrkZvtxBin"), localSample.size(), z, p);
      for (auto& track : tracks) {
        registry.fill(HIST("TracksEtaZvtxBin"), track.eta(), z, p);
        registry.fill(HIST("TracksPhiEtaBin"), track.phi(), track.eta(), p);
      }
    } else {
      registry.fill(HIST("EventSelectionBin"), 3., p);
    }
  }

  PROCESS_SWITCH(PseudorapidityDensity<TRACKTYPE>, processBinned, "Process centrality/mult. percentile binned", false);

  using Particles = aod::McParticles;
  expressions::Filter primaries = (aod::mcparticle::flags & (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary) == (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary;
  Partition<soa::Filtered<Particles>> mcSample = nabs(aod::mcparticle::eta) < estimatorEta;

  void processGen(aod::McCollisions::iterator const& mcCollision, o2::soa::SmallGroups<soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>> const& collisions, soa::Filtered<Particles> const& particles, soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtended>> const& /*tracks*/)
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
      if (select(collision)) {
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

  PROCESS_SWITCH(PseudorapidityDensity<TRACKTYPE>, processGen, "Process generator-level info", false);
};
} // namespace o2::pwgmm::multiplicity

#endif // DNDETA_H
