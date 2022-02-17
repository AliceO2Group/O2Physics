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
///
/// \brief
/// \author Paul Buehler, paul.buehler@oeaw.ac.at
/// \since  01.10.2021

#include "PWGUD/Tasks/diffMCHelpers.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct DiffQA {

  // define histograms
  HistogramRegistry registry{
    "registry",
    {
      {"timeResolution", "#timeResolution", {HistType::kTH1F, {{200, 0., 1.E3}}}},
      {"tResvsrTOFTracks", "#tResvsrTOFTracks", {HistType::kTH2F, {{200, 0., 1.E3}, {101, -0.01, 1.01}}}},
      {"numberBCs", "#numberBCs", {HistType::kTH1F, {{101, -0.5, 100.5}}}},
      {"numberTracks", "#numberTracks", {HistType::kTH1F, {{3001, -0.5, 3000.5}}}},
      {"numberVtxTracks", "#numberVtxTracks", {HistType::kTH1F, {{151, -0.5, 150.5}}}},
      {"numberGlobalTracks", "#numberGlobalTracks", {HistType::kTH1F, {{101, -0.5, 100.5}}}},
      {"numberFWDTracks", "#numberFWDTracks", {HistType::kTH1F, {{21, -0.5, 20.5}}}},
      {"VtxvsGlobalTracks", "#VtxvsGlobalTracks", {HistType::kTH2F, {{101, -0.5, 100.5}, {101, -0.5, 100.5}}}},
      {"CDEresolution", "#CDEresolution", {HistType::kTH2F, {{2, 0., 2.}, {2, 0., 2.}}}},
      {"CDEtResvsrTOFTracks", "#CDEtResvsrTOFTracks", {HistType::kTH2F, {{200, 0., 1.E3}, {101, -0.01, 1.01}}}},
      {"CDEnumberBCs", "#CDEnumberBCs", {HistType::kTH1F, {{101, -0.5, 100.5}}}},
      {"CDEnumberTracks", "#CDEnumberTracks", {HistType::kTH1F, {{3001, -0.5, 3000.5}}}},
      {"CDEnumberVtxTracks", "#CDEnumberVtxTracks", {HistType::kTH1F, {{151, -0.5, 150.5}}}},
      {"CDEnumberGlobalTracks", "#CDEnumberGlobalTracks", {HistType::kTH1F, {{101, -0.5, 100.5}}}},
      {"CDEnumberFWDTracks", "#CDEnumberFWDTracks", {HistType::kTH1F, {{21, -0.5, 20.5}}}},
      {"CDEVtxvsGlobalTracks", "#CDEVtxvsGlobalTracks", {HistType::kTH2F, {{101, -0.5, 100.5}, {101, -0.5, 100.5}}}},
      {"CDEtimeResolution", "#CDEtimeResolution", {HistType::kTH1F, {{200, 0., 1.E3}}}},
    }};

  using CCs = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;
  using CC = CCs::iterator;
  using BCs = soa::Join<aod::BCs, aod::Run3MatchedToBCSparse>;
  using TCs = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection>;
  using FWs = aod::FwdTracks;

  void process(CC const& collision, BCs const& bct0s,
               TCs& tracks, FWs& fwdtracks, aod::FT0s& ft0s, aod::FV0As& fv0as, aod::FDDs& fdds,
               aod::McCollisions& McCols, aod::McParticles& McParts)
  {

    // is this a central diffractive events?
    auto MCCol = collision.mcCollision();
    auto MCPartSlice = McParts.sliceBy(aod::mcparticle::mcCollisionId, MCCol.globalIndex());
    auto isCDE = isMcCDE(MCPartSlice);

    // obtain slice of compatible BCs
    auto bcSlice = getMCCompatibleBCs(collision, 4, bct0s);
    LOGF(debug, "<DiffQA> Number of compatible BCs: %i", bcSlice.size());

    // global tracks
    Partition<TCs> goodTracks = aod::track::isGlobalTrack > uint8_t(0);
    goodTracks.bindTable(tracks);

    // number of gobal tracks with TOF hit
    float rtrwTOF = 0.;
    for (auto& track : goodTracks) {
      if (track.hasTOF()) {
        rtrwTOF += 1.;
      }
    }
    rtrwTOF /= goodTracks.size();

    // update histograms
    if (isCDE) {
      registry.get<TH1>(HIST("CDEtimeResolution"))->Fill(collision.collisionTimeRes());
      registry.get<TH2>(HIST("CDEtResvsrTOFTracks"))->Fill(collision.collisionTimeRes(), rtrwTOF);
      registry.get<TH1>(HIST("CDEnumberBCs"))->Fill(bcSlice.size());
      registry.get<TH1>(HIST("CDEnumberTracks"))->Fill(tracks.size());
      registry.get<TH1>(HIST("CDEnumberVtxTracks"))->Fill(collision.numContrib());
      registry.get<TH1>(HIST("CDEnumberGlobalTracks"))->Fill(goodTracks.size());
      registry.get<TH1>(HIST("CDEnumberFWDTracks"))->Fill(fwdtracks.size());
      registry.get<TH2>(HIST("CDEVtxvsGlobalTracks"))->Fill(collision.numContrib(), goodTracks.size());
    } else {
      registry.get<TH1>(HIST("timeResolution"))->Fill(collision.collisionTimeRes());
      registry.get<TH2>(HIST("tResvsrTOFTracks"))->Fill(collision.collisionTimeRes(), rtrwTOF);
      registry.get<TH1>(HIST("numberBCs"))->Fill(bcSlice.size());
      registry.get<TH1>(HIST("numberTracks"))->Fill(tracks.size());
      registry.get<TH1>(HIST("numberVtxTracks"))->Fill(collision.numContrib());
      registry.get<TH1>(HIST("numberGlobalTracks"))->Fill(goodTracks.size());
      registry.get<TH1>(HIST("numberFWDTracks"))->Fill(fwdtracks.size());
      registry.get<TH2>(HIST("VtxvsGlobalTracks"))->Fill(collision.numContrib(), goodTracks.size());
    }

    // is it a DG candidate?
    // DG = no FIT signal in compatible BCs
    //    & number of forward tracks = 0
    //    & number of vertex tracks <= ntr
    //    & number of global tracks <= ntr
    bool isDGcandidate = true;
    int ntr = 2;

    // no FIT signal in compatible BCs
    for (auto& bc : bcSlice) {
      if (bc.has_ft0() || bc.has_fv0a() || bc.has_fdd()) {
        isDGcandidate = false;
        continue;
      }
    }

    // number of forward tracks = 0
    isDGcandidate &= (fwdtracks.size() == 0);

    // number of vertex tracks <= n
    isDGcandidate &= (collision.numContrib() <= ntr);

    // number of global tracks <= ntr
    isDGcandidate &= (goodTracks.size() <= ntr);

    if (isCDE) {
      if (isDGcandidate) {
        registry.get<TH2>(HIST("CDEresolution"))->Fill(0.5, 0.5);
      } else {
        registry.get<TH2>(HIST("CDEresolution"))->Fill(0.5, 1.5);
      }
    } else {
      if (isDGcandidate) {
        registry.get<TH2>(HIST("CDEresolution"))->Fill(1.5, 0.5);
      } else {
        registry.get<TH2>(HIST("CDEresolution"))->Fill(1.5, 1.5);
      }
    }
  };
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<DiffQA>(cfgc, TaskName{"diffqa"}),
  };
}
