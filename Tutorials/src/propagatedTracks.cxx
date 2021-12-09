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

//
// Task to use tables of track parameters propagated to the primary vertex
//

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/TrackPropagation.h"
#include "Common/Core/trackUtilities.h"
#include "ReconstructionDataFormats/DCA.h"
#include "DetectorsBase/Propagator.h"
#include "CommonUtils/NameConf.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct PropagatedTracksExtended {

  HistogramRegistry registry{
    "registryTracks",
    {{"hcollx", "collision x position (m);entries", {HistType::kTH1F, {{200, -0.5, 0.5}}}},
     {"hpt", "track #it{p}_{T} propagated (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 20.}}}},
     {"hphi", "track #phi propagated;entries", {HistType::kTH1F, {{140, -1., 8.}}}},
     {"heta", "track #eta propagated;entries", {HistType::kTH1F, {{200, -2., 2.}}}},
     {"hdcaXY", "propagated track dcaXY;entries", {HistType::kTH1F, {{160, -2., 2.}}}}}};

  void process(aod::Collision const& collision, soa::Join<aod::TracksPropagated, aod::TracksExtended> const& tracks)
  {
    registry.fill(HIST("hcollx"), collision.posX());
    for (auto& track : tracks) {
      registry.fill(HIST("hpt"), track.pt());
      registry.fill(HIST("hphi"), track.phi());
      registry.fill(HIST("heta"), track.eta());
      registry.fill(HIST("hdcaXY"), track.dcaXY());
    }
  }
};

struct PropagatedTracksExtra {

  HistogramRegistry registry{
    "registryTracks",
    {{"hpt", "track #it{p}_{T} propagated (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 20.}}}},
     {"hphi", "track #phi propagated;entries", {HistType::kTH1F, {{140, -1., 8.}}}},
     {"heta", "track #eta propagated;entries", {HistType::kTH1F, {{200, -2., 2.}}}}}};

  void process(aod::Collision const& collision, soa::Join<aod::TracksPropagated, aod::TracksExtra> const& tracks)
  {
    for (auto& track : tracks) {
      registry.fill(HIST("hpt"), track.pt());
      registry.fill(HIST("hphi"), track.phi());
      registry.fill(HIST("heta"), track.eta());
    }
  }
};

struct PropagatedTracksPar {

  HistogramRegistry registry{
    "registryTracks",
    {{"hX", "propagated track X position;entries", {HistType::kTH1F, {{240, -12., 12.}}}},
     {"hY", "propagated track Y position", {HistType::kTH1F, {{240, -12., 12.}}}},
     {"hZ", "propagated track Z position", {HistType::kTH1F, {{240, -12., 12.}}}}}};

  void process(aod::TrackParPropagated const& track)
  {
    registry.fill(HIST("hX"), track.x());
    registry.fill(HIST("hY"), track.y());
    registry.fill(HIST("hZ"), track.z());
  }
};

//****************************************************************************************
/**
 * Workflow definition.
 */
//****************************************************************************************
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<PropagatedTracksExtended>(cfgc),
    adaptAnalysisTask<PropagatedTracksExtra>(cfgc),
    adaptAnalysisTask<PropagatedTracksPar>(cfgc),
  };
}