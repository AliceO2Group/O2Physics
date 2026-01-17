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

/// \file pccQa.cxx
/// \brief Task producing DCA distributions with and without particle-composition correction.
/// \author Mario Krüger <mario.kruger@cern.ch>

#include "PWGLF/DataModel/particleCompositionCorrectionTable.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <Framework/AnalysisTask.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/runDataProcessing.h>

#include <TMCProcess.h>

#include <vector>

using namespace o2;
using namespace o2::framework;
using aod::track::TrackSelectionFlags;

struct PccQa {
  HistogramRegistry histos;
  Service<o2::framework::O2DatabasePDG> pdg;

  static constexpr float MaxVtxZ = 10.f;

  void init(InitContext const&);

  template <bool IS_MC, typename C, typename T>
  void processMeas(const C& collision, const T& tracks);

  using CollisionTableData = soa::Join<aod::Collisions, aod::EvSels>;
  using TrackTableData = soa::Join<aod::FullTracks, aod::TracksDCA, aod::TrackSelection>;
  void processData(CollisionTableData::iterator const& collision, TrackTableData const& tracks);
  PROCESS_SWITCH(PccQa, processData, "process data", false);

  using CollisionTableMCTrue = aod::McCollisions;
  using CollisionTableMC = soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions, aod::EvSels>>;
  using TrackTableMC = soa::Join<aod::FullTracks, aod::TracksDCA, aod::TrackSelection, aod::McTrackLabels>;
  using ParticleTableMC = soa::Join<aod::McParticles, aod::ParticleCompositionCorrection>;
  Preslice<TrackTableMC> perCollision = aod::track::collisionId;
  void processMC(CollisionTableMCTrue::iterator const& mcCollision, TrackTableMC const& tracks, CollisionTableMC const& collisions, ParticleTableMC const&);
  PROCESS_SWITCH(PccQa, processMC, "process mc", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<PccQa>(cfgc)};
}

void PccQa::init(InitContext const&)
{
  histos.add("eventCounter", "", kTH1D, {{1, 0.5, 1.5}});
  const AxisSpec dcaAxis = {1000, -1., 1., "#it{DCA}_{xy}", "dca"};
  std::vector<double> ptBinEdges = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0};
  const AxisSpec ptAxis{ptBinEdges, "#it{p}_{T} (GeV/#it{c})", "pt"};

  histos.add("DCAxyVsPt", "", kTH2D, {ptAxis, dcaAxis});

  if (doprocessMC) {
    histos.add("DCAxyVsPt_weighted", "", kTH2D, {ptAxis, dcaAxis});
    histos.add("prim/DCAxyVsPt", "", kTH2D, {ptAxis, dcaAxis});
    histos.add("prim/DCAxyVsPt_weighted", "", kTH2D, {ptAxis, dcaAxis});
    histos.add("sec/DCAxyVsPt", "", kTH2D, {ptAxis, dcaAxis});
    histos.add("sec/DCAxyVsPt_weighted", "", kTH2D, {ptAxis, dcaAxis});
    histos.add("sec/dec/DCAxyVsPt", "", kTH2D, {ptAxis, dcaAxis});
    histos.add("sec/dec/DCAxyVsPt_weighted", "", kTH2D, {ptAxis, dcaAxis});
    histos.add("sec/mat/DCAxyVsPt", "", kTH2D, {ptAxis, dcaAxis});
    histos.add("sec/mat/DCAxyVsPt_weighted", "", kTH2D, {ptAxis, dcaAxis});
  }
}

void PccQa::processData(CollisionTableData::iterator const& collision, TrackTableData const& tracks)
{
  processMeas<false>(collision, tracks);
}
void PccQa::processMC(CollisionTableMCTrue::iterator const&, TrackTableMC const& tracks, CollisionTableMC const& collisions, ParticleTableMC const&)
{
  for (const auto& collision : collisions) {
    auto curTracks = tracks.sliceBy(perCollision, collision.globalIndex());
    processMeas<true>(collision, curTracks);
    break;
  }
}

template <bool IS_MC, typename C, typename T>
void PccQa::processMeas(const C& collision, const T& tracks)
{
  if ((std::abs(collision.posZ()) > MaxVtxZ) || !collision.sel8()) {
    return;
  }
  histos.fill(HIST("eventCounter"), 1);

  for (const auto& track : tracks) {
    if (!TrackSelectionFlags::checkFlag(track.trackCutFlag(), TrackSelectionFlags::kGlobalTrackWoDCA)) {
      continue;
    }
    histos.fill(HIST("DCAxyVsPt"), track.pt(), track.dcaXY());

    if constexpr (IS_MC) {
      if (!track.has_mcParticle()) {
        continue;
      }
      const auto& particle = track.template mcParticle_as<ParticleTableMC>();

      histos.fill(HIST("DCAxyVsPt_weighted"), track.pt(), track.dcaXY(), particle.pccWeight());

      if (particle.isPhysicalPrimary()) {
        histos.fill(HIST("prim/DCAxyVsPt"), track.pt(), track.dcaXY());
        histos.fill(HIST("prim/DCAxyVsPt_weighted"), track.pt(), track.dcaXY(), particle.pccWeight());
      } else {
        histos.fill(HIST("sec/DCAxyVsPt"), track.pt(), track.dcaXY());
        histos.fill(HIST("sec/DCAxyVsPt_weighted"), track.pt(), track.dcaXY(), particle.pccWeight());
        if (particle.getProcess() == TMCProcess::kPDecay) {
          histos.fill(HIST("sec/dec/DCAxyVsPt"), track.pt(), track.dcaXY());
          histos.fill(HIST("sec/dec/DCAxyVsPt_weighted"), track.pt(), track.dcaXY(), particle.pccWeight());
        } else if (particle.getProcess() == TMCProcess::kPHInhelastic || particle.getProcess() == TMCProcess::kPHadronic || particle.getProcess() == TMCProcess::kPHElastic) {
          histos.fill(HIST("sec/mat/DCAxyVsPt"), track.pt(), track.dcaXY());
          histos.fill(HIST("sec/mat/DCAxyVsPt_weighted"), track.pt(), track.dcaXY(), particle.pccWeight());
        }
      }
    }
  }
}
