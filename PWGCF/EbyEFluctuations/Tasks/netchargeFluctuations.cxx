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

// author Nida Malik (nida.malik@cern.ch)
// Department of Physics, Aligarh Muslim University, India
// to study the net charge fluctuations by observable, #nu_dyn

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace std;

namespace o2::aod
{
using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms, aod::Mults>;
using MyCollision = MyCollisions::iterator;

using MyTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>;
using MyTrack = MyTracks::iterator;
} // namespace o2::aod

struct NetchargeFluctuations {
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  void init(o2::framework::InitContext&)
  {
    AxisSpec vtxZAxis = {100, -20, 20, "Z (cm)"};
    AxisSpec dcaAxis = {1000, -100, 100, "DCA_{xy} (cm)"};
    AxisSpec dcazAxis = {1000, -100, 100, "DCA_{z} (cm)"};
    AxisSpec ptAxis = {40, 0.0, 4.0, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec etaAxis = {30, -1.5, 1.5, "#eta"};
    AxisSpec centAxis = {100, 0., 100., "centrality"};
    AxisSpec multAxis = {2000, 0., 2000., "multiplicity"};

    histos.add("hVertexZ_bef", "", kTH1F, {vtxZAxis});
    histos.add("hVertexZ_aft", "", kTH1F, {vtxZAxis});
    histos.add("hVertexZ_aft_sel", "", kTH1D, {vtxZAxis});
    histos.add("hDCAxy_bef", "", kTH1D, {dcaAxis});
    histos.add("hDCAxy_aft", "", kTH1D, {dcaAxis});
    histos.add("hDCAz_bef", "", kTH1D, {dcazAxis});
    histos.add("hDCAz_aft", "", kTH1D, {dcazAxis});
    histos.add("hCentrality", "", kTH1D, {centAxis});
    histos.add("hMultiplicity", "", kTH1D, {multAxis});
    histos.add("hEta", "", kTH1F, {etaAxis});
    histos.add("hPt", "", kTH1F, {ptAxis});
  }

  void process(aod::MyCollision const& coll, aod::MyTracks const& inputTracks)
  {
    histos.fill(HIST("hVertexZ_bef"), coll.posZ());

    if (std::fabs(coll.posZ()) > 10.f) {
      return;
    }
    histos.fill(HIST("hVertexZ_aft"), coll.posZ());

    if (!coll.sel7()) {
      return;
    }
    histos.fill(HIST("hVertexZ_aft_sel"), coll.posZ());
    histos.fill(HIST("hCentrality"), coll.centRun2V0M());
    histos.fill(HIST("hMultiplicity"), coll.multFV0M());

    for (auto const& track : inputTracks) {
      if (std::fabs(track.eta()) > 0.8)
        continue;
      if (!(track.pt() > 0.2 && track.pt() < 2.))
        continue;

      histos.fill(HIST("hDCAxy_bef"), track.dcaXY());
      histos.fill(HIST("hDCAz_bef"), track.dcaZ());

      if (!track.isGlobalTrack())
        continue;

      histos.fill(HIST("hDCAxy_aft"), track.dcaXY());
      histos.fill(HIST("hDCAz_aft"), track.dcaZ());
      histos.fill(HIST("hPt"), track.pt());
      histos.fill(HIST("hEta"), track.eta());
    }
  }
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<NetchargeFluctuations>(cfgc)};
  return workflow;
}
