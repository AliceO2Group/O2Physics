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
// Task to use tables of track parameters propagated to the collision vertex
// Needs o2-analysis-track-propagation
//

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/trackUtilities.h"
#include "ReconstructionDataFormats/DCA.h"
#include "DetectorsBase/Propagator.h"
#include "CommonUtils/NameConf.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::constants::math;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct PropagatedTracksQa {

  HistogramRegistry registry{
    "registryTracks",
    {
      {"hcollx", ";collision x position (m);entries", {HistType::kTH1F, {{200, -0.5, 0.5}}}},
      {"hpt", ";track sign*#it{p}_{T} propagated (GeV/#it{c});", {HistType::kTH1F, {{400, -20., 20.}}}},
      {"hphi", ";track #varphi propagated;entries", {HistType::kTH1F, {{140, 0., TwoPI}}}},
      {"heta", ";track #eta propagated;entries", {HistType::kTH1F, {{200, -2., 2.}}}},
      {"hptIU", ";track sign*#it{p}_{T} unpropagated (GeV/#it{c});entries", {HistType::kTH1F, {{400, -20., 20.}}}},
      {"hphiIU", ";track #varphi unpropagated;entries", {HistType::kTH1F, {{140, 0., TwoPI}}}},
      {"hetaIU", ";track #eta unpropagated;entries", {HistType::kTH1F, {{200, -2., 2.}}}},
      {"hdcaXY", ";propagated track dcaXY;entries", {HistType::kTH1F, {{160, -2., 2.}}}},
      {"hx2D", "x unpropagated vs. x propagated", {HistType::kTH2F, {{100, -10., 10.}, {100, -10., 10.}}}},
      {"hphi2D", "#varphi unpropagated vs. #varphi propagated", {HistType::kTH2F, {{140, 0., TwoPI}, {140, 0., TwoPI}}}},
      {"heta2D", "#eta unpropagated vs. #eta propagated", {HistType::kTH2F, {{200, -2., 2.}, {200, -2., 2.}}}},
      {"hpt2D", "sign*#it{p}_{T} unpropagated vs. sign*#it{p}_{T} propagated", {HistType::kTH2F, {{400, -20., 20.}, {400, -20., 20.}}}},
      {"hdeltaphivspT", "#Delta#varphi vs. #it{p}_{T}", {HistType::kTH2F, {{200, 0, 20}, {200, -1., 1}}}},
    }};

  void process(aod::Collision const& collision, aod::TracksIU const& tracksIU, soa::Join<aod::Tracks, aod::TracksDCA> const& tracks)
  {
    registry.fill(HIST("hcollx"), collision.posX());
    for (auto& track : tracksIU) {
      registry.fill(HIST("hptIU"), track.pt() * track.sign());
      registry.fill(HIST("hphiIU"), track.phi());
      registry.fill(HIST("hetaIU"), track.eta());
    }
    for (auto& track : tracks) {
      registry.fill(HIST("hpt"), track.pt() * track.sign());
      registry.fill(HIST("hphi"), track.phi());
      registry.fill(HIST("heta"), track.eta());
      registry.fill(HIST("hdcaXY"), track.dcaXY());
    }
    for (int i = 0; i < tracks.size(); i++) {
      registry.fill(HIST("hx2D"), tracks.iteratorAt(i).x(), tracksIU.iteratorAt(i).x());
      registry.fill(HIST("hphi2D"), tracks.iteratorAt(i).phi(), tracksIU.iteratorAt(i).phi());
      registry.fill(HIST("heta2D"), tracks.iteratorAt(i).eta(), tracksIU.iteratorAt(i).eta());
      registry.fill(HIST("hpt2D"), tracks.iteratorAt(i).pt() * tracks.iteratorAt(i).sign(), tracksIU.iteratorAt(i).pt() * tracksIU.iteratorAt(i).sign());
      registry.fill(HIST("hdeltaphivspT"), tracks.iteratorAt(i).pt(), tracks.iteratorAt(i).phi() - tracksIU.iteratorAt(i).phi());
    }
  }
};

struct PropagatedTracksExtra {

  HistogramRegistry registry{
    "registryTracks",
    {
      {"hpt", ";#it{p}_{T} propagated (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 20.}}}},
      {"hcrossedrows", ";track crossed rows;entries", {HistType::kTH1F, {{160, -0.5, 159.5}}}},
    }};

  void process(aod::Collision const& collision, soa::Join<aod::Tracks, aod::TracksExtra> const& tracks)
  {
    for (auto& track : tracks) {
      registry.fill(HIST("hpt"), track.pt());
      registry.fill(HIST("hcrossedrows"), track.tpcNClsCrossedRows());
    }
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
    adaptAnalysisTask<PropagatedTracksQa>(cfgc),
    adaptAnalysisTask<PropagatedTracksExtra>(cfgc),
  };
}