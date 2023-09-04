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

#include "Framework/AnalysisTask.h"
#include "PWGUD/DataModel/McPIDTable.h"
#include "PWGUD/Core/UDHelpers.h"

using namespace o2;
using namespace o2::framework;

void customize(std::vector<ConfigParamSpec>& workflowOptions)
{
  workflowOptions.push_back(ConfigParamSpec{"runCase", VariantType::Int, 0, {"runCase: 0 - histos,  1 - mcTruth, 2 - mc only, else - tree"}});
}

#include "Framework/runDataProcessing.h"

using namespace o2::framework::expressions;

// Loop over collisions
// find for each collision the number of compatible BCs
struct CompatibleBCs {
  using CCs = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;
  using CC = CCs::iterator;

  void process(CC const& collision, aod::BCs const& bcs)
  {
    auto bcSlice = udhelpers::MCcompatibleBCs(collision, 4, bcs);
    LOGF(debug, "  Number of possible BCs: %i", bcSlice.size());
    for (auto& bc : bcSlice) {
      LOGF(debug, "    This collision may belong to BC %lld", bc.globalBC());
    }
  }
};

// Fill histograms with collision and compatible BCs related information
// runCase = 0
struct collisionsInfo {
  int cnt = 0;
  HistogramRegistry registry{
    "registry",
    {
      {"timeResolution", "#timeResolution", {HistType::kTH1F, {{200, 0., 1.E3}}}},
      {"numberBCs", "#numberBCs", {HistType::kTH1F, {{101, -0.5, 100.5}}}},
      {"DGCandidate", "#DGCandidate", {HistType::kTH1F, {{2, -0.5, 1.5}}}},
      {"numberTracks", "#numberTracks", {HistType::kTH1F, {{5001, -0.5, 5000.5}}}},
      {"numberVtxTracks", "#numberVtxTracks", {HistType::kTH1F, {{101, -0.5, 100.5}}}},
      {"numberGlobalTracks", "#numberGlobalTracks", {HistType::kTH1F, {{101, -0.5, 100.5}}}},
      {"netCharge", "#netCharge", {HistType::kTH1F, {{3, -1.5, 1.5}}}},
      //      {"numberMFTTracks", "#numberMFTTracks", {HistType::kTH1F, {{101, -0.5, 100.5}}}},
      //      {"etaMFTAll", "#etaMFTAll", {HistType::kTH1F, {{100, -5.0, 5.0}}}},
      //      {"etaMFTDG", "#etaMFTDG", {HistType::kTH1F, {{100, -5.0, 5.0}}}},
      {"numberFWDTracks", "#numberFWDTracks", {HistType::kTH1F, {{101, -0.5, 100.5}}}},
      {"etaFWDAll", "#etaFWDAll", {HistType::kTH1F, {{100, -5.0, 5.0}}}},
      {"etaFWDDG", "#etaFWDDG", {HistType::kTH1F, {{100, -5.0, 5.0}}}},
      //      {"globalVsMFTAll", "#globalVsMFTAll", {HistType::kTH2F, {{21, -0.5, 20.5}, {21, -0.5, 20.5}}}},
      //      {"globalVsMFTDG", "#globalVsMFTDG", {HistType::kTH2F, {{21, -0.5, 20.5}, {21, -0.5, 20.5}}}},
      {"VtxvsTracks", "#VtxvsTracks", {HistType::kTH2F, {{101, -0.5, 100.5}, {5001, -0.5, 5000.5}}}},
      {"VtxvsGlobalTracks", "#VtxvsGlobalTracks", {HistType::kTH2F, {{101, -0.5, 100.5}, {101, -0.5, 100.5}}}},
    }};

  void init(o2::framework::InitContext&)
  {
    registry.get<TH1>(HIST("timeResolution"))->GetXaxis()->SetTitle("Time resolution [ns]");
    registry.get<TH1>(HIST("numberBCs"))->GetXaxis()->SetTitle("Number of compatible BCs");
    registry.get<TH1>(HIST("numberTracks"))->GetXaxis()->SetTitle("Number of tracks");
    registry.get<TH1>(HIST("numberVtxTracks"))->GetXaxis()->SetTitle("Number of Vtx tracks");
    registry.get<TH1>(HIST("numberGlobalTracks"))->GetXaxis()->SetTitle("Number of global tracks");
    registry.get<TH1>(HIST("netCharge"))->GetXaxis()->SetTitle("Sign of net charge");
    //    registry.get<TH1>(HIST("numberMFTTracks"))->GetXaxis()->SetTitle("Number of MFT tracks");
    //    registry.get<TH1>(HIST("etaMFTAll"))->GetXaxis()->SetTitle("Pseudo rapidity");
    //    registry.get<TH1>(HIST("etaMFTDG"))->GetXaxis()->SetTitle("Pseudo rapidity");
    registry.get<TH1>(HIST("numberFWDTracks"))->GetXaxis()->SetTitle("Number of FWD tracks");
    registry.get<TH1>(HIST("etaFWDAll"))->GetXaxis()->SetTitle("Pseudo rapidity");
    registry.get<TH1>(HIST("etaFWDDG"))->GetXaxis()->SetTitle("Pseudo rapidity");
    //    registry.get<TH2>(HIST("globalVsMFTAll"))->GetXaxis()->SetTitle("Number of global tracks");
    //    registry.get<TH2>(HIST("globalVsMFTAll"))->GetYaxis()->SetTitle("Number of MFT tracks");
    //    registry.get<TH2>(HIST("globalVsMFTDG"))->GetXaxis()->SetTitle("Number of global tracks");
    //    registry.get<TH2>(HIST("globalVsMFTDG"))->GetYaxis()->SetTitle("Number of MFT tracks");
    registry.get<TH2>(HIST("VtxvsTracks"))->GetXaxis()->SetTitle("Number of Vtx tracks");
    registry.get<TH2>(HIST("VtxvsTracks"))->GetYaxis()->SetTitle("Number of tracks");
    registry.get<TH2>(HIST("VtxvsGlobalTracks"))->GetXaxis()->SetTitle("Number of Vtx tracks");
    registry.get<TH2>(HIST("VtxvsGlobalTracks"))->GetYaxis()->SetTitle("Number of global tracks");
  }

  using CCs = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;
  using CC = CCs::iterator;
  using BCs = soa::Join<aod::BCs, aod::BcSels, aod::Run3MatchedToBCSparse>;
  using TCs = soa::Join<aod::Tracks, aod::TrackSelection>;
  using MFs = aod::MFTTracks;
  using FWs = aod::FwdTracks;

  void process(CC const& collision, BCs const& bct0s,
               TCs& tracks, /* MFs& mfttracks,*/ FWs& fwdtracks, aod::FT0s& ft0s, aod::FV0As& fv0as, aod::FDDs& fdds,
               aod::McCollisions& McCols, aod::McParticles& McParts)
  {

    // obtain slice of compatible BCs
    auto bcSlice = udhelpers::MCcompatibleBCs(collision, 4, bct0s);
    LOGF(info, "  Number of compatible BCs: %i", bcSlice.size());
    registry.get<TH1>(HIST("numberBCs"))->Fill(bcSlice.size());

    // check that there are no FIT signals in any of the compatible BCs
    std::vector<float> const lims{0.0, 0.0, 0.0, 0.0, 0.0}; // amplitude thresholds: FV0A, FT0A, FT0C, FDDA, FDDC
    auto isDGcandidate = true;
    for (auto& bc : bcSlice) {
      if (!udhelpers::cleanFIT(bc, 4., lims)) {
        isDGcandidate = false;
        break;
      }
    }

    // global tracks
    Partition<TCs> goodTracks = requireGlobalTrackInFilter();
    goodTracks.bindTable(tracks);
    int cntGlobal = goodTracks.size();

    // count tracks
    int cntAll = 0;
    int netCharge = 0;
    for (auto track : tracks) {
      cntAll++;
      netCharge += track.sign();
    }
    // netCharge /= abs(netCharge);

    LOGF(info, "Global tracks: %i", cntGlobal);

    if (isDGcandidate) {
      LOGF(info, "  This is a DG candidate with %d tracks and %d net charge.", tracks.size(), netCharge);
      registry.get<TH1>(HIST("DGCandidate"))->Fill(1.);
    } else {
      registry.get<TH1>(HIST("DGCandidate"))->Fill(0.);
    }

    // update histograms with track information
    LOGF(debug, "Number of tracks: Vertex %d, total %d, global %d", collision.numContrib(), cntAll, cntGlobal);
    registry.get<TH1>(HIST("numberTracks"))->Fill(cntAll);
    registry.get<TH1>(HIST("numberVtxTracks"))->Fill(collision.numContrib());
    registry.get<TH1>(HIST("numberGlobalTracks"))->Fill(cntGlobal);
    registry.get<TH1>(HIST("netCharge"))->Fill(netCharge);
    //    registry.get<TH1>(HIST("numberMFTTracks"))->Fill(mfttracks.size());
    registry.get<TH1>(HIST("numberFWDTracks"))->Fill(fwdtracks.size());
    //    registry.get<TH2>(HIST("globalVsMFTAll"))->Fill(cntGlobal, mfttracks.size());
    registry.get<TH2>(HIST("VtxvsTracks"))->Fill(collision.numContrib(), cntAll);
    registry.get<TH2>(HIST("VtxvsGlobalTracks"))->Fill(collision.numContrib(), cntGlobal);

    // loop over MFT tracks
    //    LOGF(debug, "MFT tracks: %i", mfttracks.size());
    //    for (auto mfttrack : mfttracks) {
    //      registry.get<TH1>(HIST("etaMFTAll"))->Fill(mfttrack.eta());
    //    }

    // loop over FWD tracks
    LOGF(info, "FWD tracks: %i", fwdtracks.size());
    for (auto fwdtrack : fwdtracks) {
      registry.get<TH1>(HIST("etaFWDAll"))->Fill(fwdtrack.eta());
    }

    // update timeResolution
    registry.get<TH1>(HIST("timeResolution"))->Fill(collision.collisionTimeRes());

    if (isDGcandidate) {
      // loop over MFT tracks
      //      for (auto mfttrack : mfttracks) {
      //        registry.get<TH1>(HIST("etaMFTDG"))->Fill(mfttrack.eta());
      //      }
      //      registry.get<TH2>(HIST("globalVsMFTDG"))->Fill(cntGlobal, mfttracks.size());

      // loop over FWD tracks
      for (auto fwdtrack : fwdtracks) {
        registry.get<TH1>(HIST("etaFWDDG"))->Fill(fwdtrack.eta());
      }
    }
    cnt++;
    LOGF(info, "#Collisions: %d", cnt);
  }
};

// Loop over BCs
// check aliases, selection, and FIT signals per BC
// runCase = 0
struct BCInfo {
  int cnt = 0;
  HistogramRegistry registry{
    "registry",
    {{"numberCollisions", "#numberCollisions", {HistType::kTH1F, {{11, -0.5, 10.5}}}},
     {"numberCollisionsGT", "#numberCollisionsGT", {HistType::kTH1F, {{11, -0.5, 10.5}}}},
     {"Aliases", "#Aliases", {HistType::kTH1F, {{kNaliases, 0., kNaliases}}}},
     {"Selection", "#Selection", {HistType::kTH1F, {{aod::evsel::kNsel, 0., aod::evsel::kNsel}}}},
     {"DetectorSignals", "#DetectorSignals", {HistType::kTH1F, {{6, 0., 6}}}}}};

  void init(o2::framework::InitContext&)
  {
    registry.get<TH1>(HIST("numberCollisions"))->GetXaxis()->SetTitle("#Collisions per BC");
    registry.get<TH1>(HIST("numberCollisionsGT"))->GetXaxis()->SetTitle("#Collisions with good time per BC");
  }

  using BBs = soa::Join<aod::BCs, aod::BcSels, aod::Run3MatchedToBCSparse>;
  using BB = BBs::iterator;

  void process(BB const& bc, aod::Collisions const& cols)
  {
    LOGF(debug, "BC: %d number of collisions: %d", bc.globalBC(), cols.size());
    registry.get<TH1>(HIST("numberCollisions"))->Fill(cols.size());

    // count collisions with good time resoluton
    auto nColGT = 0;
    for (auto col : cols) {
      if (col.collisionTimeRes() <= 20.) {
        nColGT++;
      }
    }
    registry.get<TH1>(HIST("numberCollisionsGT"))->Fill(nColGT);

    // update Aliases
    for (auto ii = 0; ii < kNaliases; ii++) {
      registry.get<TH1>(HIST("Aliases"))->Fill(ii, bc.alias_bit(ii));
    }

    // update Selection
    for (auto ii = 0; ii < aod::evsel::kNsel; ii++) {
      registry.get<TH1>(HIST("Selection"))->Fill(ii, bc.selection_bit(ii));
    }

    // FIT detector signals
    if (!bc.has_foundFT0()) {
      registry.get<TH1>(HIST("DetectorSignals"))->Fill(0., 1.);
    }
    if (!bc.has_foundFV0()) {
      registry.get<TH1>(HIST("DetectorSignals"))->Fill(1., 1.);
    }
    if (!bc.has_foundFDD()) {
      registry.get<TH1>(HIST("DetectorSignals"))->Fill(2., 1.);
    }
    if (!bc.has_zdc()) {
      registry.get<TH1>(HIST("DetectorSignals"))->Fill(3., 1.);
    }
    auto noFIT = !bc.has_foundFV0() && !bc.has_foundFT0() && !bc.has_foundFDD();
    if (noFIT) {
      registry.get<TH1>(HIST("DetectorSignals"))->Fill(4., 1.);
    }
    if (noFIT && !bc.has_zdc()) {
      registry.get<TH1>(HIST("DetectorSignals"))->Fill(5., 1.);
    }
    cnt++;
    LOGF(debug, "#BCs: %d", cnt);
  }
};

// Loop over tracks
// Make histograms with track type and time resolution
// runCase = 0
struct TrackTypes {
  HistogramRegistry registry{
    "registry",
    {
      {"nTracks", "#nTracks", {HistType::kTH2F, {{6, -0.5, 5.5}, {2, 0., 2.}}}},
      {"timeRes", "#timeRes", {HistType::kTH2F, {{6, -0.5, 5.5}, {2, 0., 2.}}}},
      {"FwdType", "#FwdType", {HistType::kTH2F, {{7, -0.5, 6.5}, {1, -0.5, 0.5}}}},
    }};

  using TCs = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection>;
  using MTs = aod::MFTTracks;
  using FTs = aod::FwdTracks;

  void process(TCs const& tracks, /*MTs const& mfttracks,*/ FTs const& fwdtracks)
  {
    for (auto track : tracks) {
      LOGF(debug, "isGlobal %i Detector map %i %i %i %i time resolution %f", track.isGlobalTrack(),
           track.hasITS(), track.hasTPC(), track.hasTRD(), track.hasTOF(), track.trackTimeRes());

      float isGlobal = track.isGlobalTrack() ? 1. : 0.;

      registry.get<TH2>(HIST("nTracks"))->Fill(0., isGlobal, 1.);
      registry.get<TH2>(HIST("timeRes"))->Fill(0., isGlobal, track.trackTimeRes());

      // has associated collision
      if (track.collisionId() >= 0) {
        registry.get<TH2>(HIST("nTracks"))->Fill(1., isGlobal, 1.);
        registry.get<TH2>(HIST("timeRes"))->Fill(1., isGlobal, track.trackTimeRes());
      }

      // has ITS hit
      if (track.hasITS()) {
        registry.get<TH2>(HIST("nTracks"))->Fill(2., isGlobal, 1.);
        registry.get<TH2>(HIST("timeRes"))->Fill(2., isGlobal, track.trackTimeRes());
      }

      // has TPC hit
      if (track.hasTPC()) {
        registry.get<TH2>(HIST("nTracks"))->Fill(3., isGlobal, 1.);
        registry.get<TH2>(HIST("timeRes"))->Fill(3., isGlobal, track.trackTimeRes());
      }

      // has TRD hit
      if (track.hasTRD()) {
        registry.get<TH2>(HIST("nTracks"))->Fill(4., isGlobal, 1.);
        registry.get<TH2>(HIST("timeRes"))->Fill(4., isGlobal, track.trackTimeRes());
      }

      // has TOF hit
      if (track.hasTOF()) {
        registry.get<TH2>(HIST("nTracks"))->Fill(5., isGlobal, 1.);
        registry.get<TH2>(HIST("timeRes"))->Fill(5., isGlobal, track.trackTimeRes());
      }
    }

    // ForwardTrackTypeEnum has 5 values
    auto nTypes = 5;
    for (auto fwdtrack : fwdtracks) {
      registry.get<TH2>(HIST("FwdType"))->Fill(0., 0., 1.);
      if (fwdtrack.collisionId() >= 0) {
        registry.get<TH2>(HIST("FwdType"))->Fill(1., 0., 1.);
      }
      for (auto ii = 0; ii < nTypes; ii++) {
        if ((fwdtrack.trackType() & (1 << ii)) > 0) {
          registry.get<TH2>(HIST("FwdType"))->Fill((Double_t)(ii + 2), 0., 1.);
        }
      }
    }
  }
};

// MCTruth tracks
// runCase = 1
struct MCTracks {

  using CCs = soa::Join<aod::Collisions, aod::McCollisionLabels>;
  using CC = CCs::iterator;

  Preslice<aod::McParticles> perMcCollision = aod::mcparticle::mcCollisionId;

  void process(CCs const& collisions, aod::McCollisions& McCols, aod::McParticles& McParts)
  {

    for (auto collision : collisions) {

      // get McCollision which belongs to collision
      auto MCCol = collision.mcCollision();

      LOGF(info, "Collision %i MC collision %i Type %i",
           collision.globalIndex(), MCCol.globalIndex(), MCCol.generatorsID());

      // get MCParticles which belong to MCCol
      auto MCPartSlice = McParts.sliceBy(perMcCollision, MCCol.globalIndex());
      LOGF(info, "  Number of McParticles %i", MCPartSlice.size());

      // loop over particles
      float etot = 0.0;
      float px = 0.0, py = 0.0, pz = 0.0;
      bool hasDiff = false;
      int prongs = 0;

      for (auto mcpart : MCPartSlice) {
        LOGF(info, " MCPart: %i %i %i %i %i - %i", mcpart.mcCollisionId(), mcpart.isPhysicalPrimary(), mcpart.getProcess(), mcpart.getGenStatusCode(), mcpart.globalIndex(), mcpart.pdgCode());
        if (mcpart.pdgCode() == 9900110) {
          LOGF(info, "  rho_diff0 energy: %f", mcpart.e());
          hasDiff = true;
        }

        // retrieve mothers
        auto mothers = mcpart.mothers_as<aod::McParticles>();

        LOGF(info, "Mothers %i", mothers.size());
        // for (auto mother : mothers) {
        //   LOGF(info, "  %i %i %i", mother.globalIndex(), mother.pdgCode(), mother.getGenStatusCode());
        // }

        if (hasDiff && mothers.size() > 1) {
          auto mom1 = mothers[0];
          auto mom2 = mothers[1];
          if (mcpart.isPhysicalPrimary() &&
              (mcpart.getGenStatusCode() == 1 || mcpart.getGenStatusCode() == 2) &&
              mom1.globalIndex() != mom2.globalIndex() &&
              mom2.globalIndex() > 0) {

            prongs++;
            etot += mcpart.e();
            px += mcpart.px();
            py += mcpart.py();
            pz += mcpart.pz();
          }
        }
      }
      if (hasDiff) {
        auto mass = TMath::Sqrt(etot * etot - (px * px + py * py + pz * pz));
        LOGF(info, "  mass of X: %f, prongs: %i", mass, prongs);
      }
    }
  }
};

// MC only
// runCase = 2
struct MConly {

  void process(aod::McCollision const& mcCollision, aod::McParticles const& mcParticles)
  {
    LOGF(info, "Number of MC particles %i", mcParticles.size());
  }
};

// TPC nSigma
// runCase > 2
struct TPCnSigma {
  Produces<aod::UDnSigmas> nSigmas;

  using TCs = soa::Join<aod::Tracks, aod::TrackSelection, aod::McTrackLabels>;
  using TCwPIDs = soa::Join<TCs, aod::pidTPCEl, aod::pidTPCMu, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr>;

  void process(TCwPIDs& tracks, aod::McParticles const& mcParticles)
  {
    for (auto track : tracks) {
      if (track.isGlobalTrack()) {
        nSigmas(track.mcParticle_as<aod::McParticles>().pdgCode(), track.mcParticle_as<aod::McParticles>().pt(),
                track.tpcNSigmaEl(), track.tpcNSigmaMu(), track.tpcNSigmaPi(),
                track.tpcNSigmaKa(), track.tpcNSigmaPr());
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  if (cfgc.options().get<int>("runCase") == 0) {
    return WorkflowSpec{
      // adaptAnalysisTask<CompatibleBCs>(cfgc),
      adaptAnalysisTask<collisionsInfo>(cfgc, TaskName{"collisioninformation"}),
      adaptAnalysisTask<BCInfo>(cfgc, TaskName{"bcinformation"}),
      adaptAnalysisTask<TrackTypes>(cfgc, TaskName{"tracktypes"}),
    };
  } else if (cfgc.options().get<int>("runCase") == 1) {
    return WorkflowSpec{
      adaptAnalysisTask<MCTracks>(cfgc, TaskName{"mctracks"}),
    };
  } else if (cfgc.options().get<int>("runCase") == 2) {
    return WorkflowSpec{
      adaptAnalysisTask<MConly>(cfgc, TaskName{"mconly"}),
    };
  } else {
    return WorkflowSpec{
      adaptAnalysisTask<TPCnSigma>(cfgc, TaskName{"tpcnsigma"}),
    };
  }
}
