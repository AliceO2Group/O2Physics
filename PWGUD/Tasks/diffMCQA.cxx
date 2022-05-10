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

#include "TLorentzVector.h"
#include "CommonConstants/PhysicsConstants.h"
#include "PWGUD/Tasks/diffMCHelpers.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct DiffQA {

  Configurable<int> ntrMin{"minNumberTracks", 2, "Minimum number of tracks"};
  Configurable<int> ntrMax{"maxNumberTracks", 2, "Maximum number of tracks"};
  Configurable<int> nchMin{"minNetCharge", 0, "Minimum net charge"};
  Configurable<int> nchMax{"maxNetCharge", 0, "Maximum net charge"};
  Configurable<float> ptMin{"minpt", 0., "Minimum pt of particles"};
  Configurable<float> ptMax{"maxpt", 3., "Maximum pt of particles"};
  Configurable<float> etaMin{"minEta", -1.5, "Minimum eta of particles"};
  Configurable<float> etaMax{"maxEta", 1.5, "Maximum eta of particles"};
  Configurable<float> massMin{"minMass", 0., "Minimum invariant mass"};
  Configurable<float> massMax{"maxMass", 3., "Maximum invariant mass"};
  Configurable<int> pidHypo{"pidHypoyhesis", 211, "PID of decay particles"};

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
      {"TCH", "#CDETCH", {HistType::kTH1F, {{21, -10.5, 10.5}}}},
      {"IVM", "#IVM", {HistType::kTH2F, {{150, 0., 3.0}, {7, -0.5, 6.5}}}},
      {"pt", "#pt", {HistType::kTH2F, {{150, 0., 3.0}, {7, -0.5, 6.5}}}},
      {"eta", "#eta", {HistType::kTH2F, {{150, -1.5, 1.5}, {7, -0.5, 6.5}}}},

      {"CDEtimeResolution", "#CDEtimeResolution", {HistType::kTH1F, {{200, 0., 1.E3}}}},
      {"CDEtResvsrTOFTracks", "#CDEtResvsrTOFTracks", {HistType::kTH2F, {{200, 0., 1.E3}, {101, -0.01, 1.01}}}},
      {"CDEnumberBCs", "#CDEnumberBCs", {HistType::kTH1F, {{101, -0.5, 100.5}}}},
      {"CDEnumberTracks", "#CDEnumberTracks", {HistType::kTH1F, {{3001, -0.5, 3000.5}}}},
      {"CDEnumberVtxTracks", "#CDEnumberVtxTracks", {HistType::kTH1F, {{151, -0.5, 150.5}}}},
      {"CDEnumberGlobalTracks", "#CDEnumberGlobalTracks", {HistType::kTH1F, {{101, -0.5, 100.5}}}},
      {"CDEnumberFWDTracks", "#CDEnumberFWDTracks", {HistType::kTH1F, {{21, -0.5, 20.5}}}},
      {"CDEVtxvsGlobalTracks", "#CDEVtxvsGlobalTracks", {HistType::kTH2F, {{101, -0.5, 100.5}, {101, -0.5, 100.5}}}},
      {"CDETCH", "#CDETCH", {HistType::kTH1F, {{21, -10.5, 10.5}}}},
      {"CDEIVM", "#CDEIVM", {HistType::kTH2F, {{150, 0., 3.0}, {7, -0.5, 6.5}}}},
      {"CDEpt", "#CDEpt", {HistType::kTH2F, {{150, 0., 3.0}, {7, -0.5, 6.5}}}},
      {"CDEeta", "#CDEeta", {HistType::kTH2F, {{150, -1.5, 1.5}, {7, -0.5, 6.5}}}},

      {"Stat", "#Stat", {HistType::kTH1F, {{12, -0.5, 11.5}}}},
      {"Efficiency", "#Efficiency", {HistType::kTH2F, {{3, -0.5, 2.5}, {2, -0.5, 1.5}}}},
    }};

  using CCs = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;
  using CC = CCs::iterator;
  using BCs = soa::Join<aod::BCs, aod::Run3MatchedToBCSparse>;
  using TCs = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::McTrackLabels>;
  using FWs = aod::FwdTracks;

  void process(CC const& collision, BCs const& bct0s,
               TCs& tracks, FWs& fwdtracks, aod::FT0s& ft0s, aod::FV0As& fv0as, aod::FDDs& fdds,
               aod::McCollisions& McCols, aod::McParticles const& McParts)
  {
    LOGF(debug, "<DiffQA> Start");

    // is this a central diffractive event?
    auto MCCol = collision.mcCollision();
    auto MCPartSlice = McParts.sliceBy(aod::mcparticle::mcCollisionId, MCCol.globalIndex());
    auto isPythiaDiff = isPythiaCDE(MCPartSlice);
    auto isGraniittiDiff = isGraniittiCDE(MCPartSlice);

    // obtain slice of compatible BCs
    auto bcSlice = getMCCompatibleBCs(collision, 4, bct0s);
    LOGF(debug, "<DiffQA> Number of compatible BCs: %i", bcSlice.size());

    // global tracks
    Partition<TCs> goodTracks = requireGlobalTrackInFilter();
    goodTracks.bindTable(tracks);
    LOGF(debug, "<DiffQA> Number of good tracks: %i", goodTracks.size());

    // PV tracks
    // Partition<TCs> vtxTracks = aod::track::isPVContributor == true;
    // vtxTracks.bindTable(tracks);
    // LOGF(debug, "<DiffQA> Number of vtx tracks: %i", vtxTracks.size());

    // number of gobal tracks with TOF hit
    float rgtrwTOF = 0.;
    if (goodTracks.size() > 0) {
      for (auto& track : goodTracks) {
        if (track.hasTOF()) {
          rgtrwTOF += 1.;
        }
      }
      rgtrwTOF /= goodTracks.size();
    }
    LOGF(debug, "<DiffQA> Good tracks with TOF: %f [1]", rgtrwTOF);

    // update histograms
    if (isGraniittiDiff) {
      registry.get<TH1>(HIST("CDEtimeResolution"))->Fill(collision.collisionTimeRes());
      registry.get<TH2>(HIST("CDEtResvsrTOFTracks"))->Fill(collision.collisionTimeRes(), rgtrwTOF);
      registry.get<TH1>(HIST("CDEnumberBCs"))->Fill(bcSlice.size());
      registry.get<TH1>(HIST("CDEnumberTracks"))->Fill(tracks.size());
      registry.get<TH1>(HIST("CDEnumberVtxTracks"))->Fill(collision.numContrib());
      registry.get<TH1>(HIST("CDEnumberGlobalTracks"))->Fill(goodTracks.size());
      registry.get<TH1>(HIST("CDEnumberFWDTracks"))->Fill(fwdtracks.size());
      registry.get<TH2>(HIST("CDEVtxvsGlobalTracks"))->Fill(collision.numContrib(), goodTracks.size());
    } else {
      registry.get<TH1>(HIST("timeResolution"))->Fill(collision.collisionTimeRes());
      registry.get<TH2>(HIST("tResvsrTOFTracks"))->Fill(collision.collisionTimeRes(), rgtrwTOF);
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
    //    & ntrMin <= number of vertex tracks <= ntrMax
    //    & no global track which is not a vertex track
    bool isDGcandidate = true;

    // no FIT signal in compatible BCs
    for (auto& bc : bcSlice) {
      if (bc.has_ft0() || bc.has_fv0a() || bc.has_fdd()) {
        isDGcandidate = false;
        continue;
      }
    }
    registry.get<TH1>(HIST("Stat"))->Fill(0., isDGcandidate * 1.);

    // number of forward tracks = 0
    isDGcandidate &= (fwdtracks.size() == 0);
    registry.get<TH1>(HIST("Stat"))->Fill(1., isDGcandidate * 1.);

    // no global tracks which are no vtx tracks
    for (auto& track : tracks) {
      if (track.isGlobalTrack() && !track.isPVContributor()) {
        isDGcandidate = false;
        continue;
      }
    }
    registry.get<TH1>(HIST("Stat"))->Fill(2., isDGcandidate * 1.);
    LOGF(debug, "<DiffQA> isDGcandidate: %i / %i", isPythiaDiff || isGraniittiDiff, isDGcandidate);

    // number of vertex tracks <= n
    for (int ii = 2; ii <= 10; ii++) {
      registry.get<TH1>(HIST("Stat"))->Fill(ii + 1, isDGcandidate * (collision.numContrib() <= (12 - ii)) * 1.);
    }
    isDGcandidate &= (collision.numContrib() >= ntrMin);
    isDGcandidate &= (collision.numContrib() <= ntrMax);

    // invariant mass
    if (isDGcandidate) {

      // which particle hypothesis?
      auto mass2Use = constants::physics::MassPionCharged;
      if (pidHypo == 321) {
        mass2Use = constants::physics::MassKaonCharged;
      }

      auto netCharge = 0;
      auto lvtmp = TLorentzVector();
      auto ivm = TLorentzVector();
      for (auto& track : tracks) {
        if (track.isPVContributor()) {
          lvtmp.SetXYZM(track.px(), track.py(), track.pz(), mass2Use);
          if (lvtmp.Perp() < ptMin || lvtmp.Perp() > ptMax) {
            isDGcandidate = false;
            break;
          }
          if (lvtmp.Eta() < etaMin || lvtmp.Eta() > etaMax) {
            isDGcandidate = false;
            break;
          }
          netCharge += track.sign();
          ivm += lvtmp;

          // if (track.isPVContributor() && track.mcParticleId()>0) {
          // auto mcpart = track.mcParticle();
          // auto plv = TLorentzVector(mcpart.px(),mcpart.py(),mcpart.pz(),mcpart.e());
          // registry.get<TH1>(HIST("CDEpt"))->Fill(plv.Perp());
          // ivm += plv;
        }
      }
      isDGcandidate &= (netCharge >= nchMin);
      isDGcandidate &= (netCharge <= nchMax);
      isDGcandidate &= (ivm.M() >= massMin);
      isDGcandidate &= (ivm.M() <= massMax);

      if (isDGcandidate) {
        LOGF(debug, "<DiffQA> Invariant mass: %f", ivm.M());

        // invariant mass
        if (!isGraniittiDiff) {
          registry.get<TH1>(HIST("TCH"))->Fill(netCharge);
          registry.get<TH2>(HIST("IVM"))->Fill(ivm.M(), collision.numContrib());
        } else {
          registry.get<TH1>(HIST("CDETCH"))->Fill(netCharge);
          registry.get<TH2>(HIST("CDEIVM"))->Fill(ivm.M(), collision.numContrib());
        }

        // track pt
        for (auto& track : tracks) {
          if (track.isPVContributor()) {
            lvtmp.SetXYZM(track.px(), track.py(), track.pz(), mass2Use);
            if (!isGraniittiDiff) {
              registry.get<TH2>(HIST("pt"))->Fill(lvtmp.Perp(), collision.numContrib());
              registry.get<TH2>(HIST("eta"))->Fill(lvtmp.Eta(), collision.numContrib());
            } else {
              registry.get<TH2>(HIST("CDEpt"))->Fill(lvtmp.Perp(), collision.numContrib());
              registry.get<TH2>(HIST("CDEeta"))->Fill(lvtmp.Eta(), collision.numContrib());
            }
          }
        }
      }
    }

    // update Efficiency
    registry.get<TH2>(HIST("Efficiency"))->Fill(0. + isPythiaDiff * 1. + isGraniittiDiff * 2., 0. + isDGcandidate * 1.);

    LOGF(debug, "<DiffQA> End");
  };
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<DiffQA>(cfgc, TaskName{"diffqa"}),
  };
}
