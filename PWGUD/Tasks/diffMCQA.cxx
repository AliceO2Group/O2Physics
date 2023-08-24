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
/// \brief A QA task for DG events
/// \author Paul Buehler, paul.buehler@oeaw.ac.at
/// \since  01.10.2021

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "TLorentzVector.h"
#include "ReconstructionDataFormats/BCRange.h"
#include "CommonConstants/PhysicsConstants.h"
#include "PWGUD/Core/UDHelpers.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct DiffMCQA {
  SliceCache cache;
  Preslice<aod::Zdcs> perBCzdc = aod::zdc::bcId;
  Preslice<aod::Calos> perBCcalo = aod::calo::bcId;

  float maxdEdxTPC;
  float maxdEdxTOF;

  // get a DGCutparHolder
  DGCutparHolder diffCuts = DGCutparHolder();
  Configurable<DGCutparHolder> DGCuts{"DGCuts", {}, "DG event cuts"};
  Configurable<bool> withAmbTrackAnalysis{"ambiguousTracks", false, "with ambiguous tracks analysis"};
  Configurable<bool> withAmbFwdTrackAnalysis{"ambiguousFwdTracks", false, "with ambiguous forward tracks analysis"};
  Configurable<bool> doCleanFITBC{"doCleanFITBC", false, "Require cleanFIT in compatible BCs"};

  // structures to hold information about the possible BCs the ambiguous tracks/FwdTracks belong to
  o2::dataformats::bcRanges abcrs = o2::dataformats::bcRanges("ambiguous_tracks");
  o2::dataformats::bcRanges afbcrs = o2::dataformats::bcRanges("ambiguous_fwdtracks");

  // 3 different versions of histograms:
  //  MBRDiff: Pythia MBR
  //  Diff   : GRANIITTI
  //  nonDiff: Rest (Pythia MB)
  //
  // initialize HistogramRegistry
  HistogramRegistry registry{
    "registry",
    {}};

  // define abbreviations
  using CCs = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;
  using CC = CCs::iterator;
  using BCs = soa::Join<aod::BCs, aod::BcSels, aod::Run3MatchedToBCSparse>;
  using TCs = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TrackSelectionExtension, aod::TOFSignal>;
  using FWs = soa::Join<aod::FwdTracks, aod::McFwdTrackLabels>;
  using ATs = aod::AmbiguousTracks;
  using AFTs = aod::AmbiguousFwdTracks;

  void init(InitContext& context)
  {
    maxdEdxTPC = 0.;
    maxdEdxTOF = 0.;
    diffCuts = (DGCutparHolder)DGCuts;

    // add histograms for the different process functions
    if (context.mOptions.get<bool>("processMCTruth")) {
      registry.add("MCTruth/Stat", "Simulated event type; event type; Entries", {HistType::kTH1F, {{4, 0.5, 4.5}}});
      registry.add("MCTruth/recCols", "Number of reconstructed collisions; Number of reconstructed collisions; Entries", {HistType::kTH1F, {{31, -0.5, 30.5}}});
    }

    if (context.mOptions.get<bool>("processMain")) {
      registry.add("all/mcCols", "Has MC truth collision; Has MC truth collision; Entries", {HistType::kTH1F, {{2, -0.5, 1.5}}});

      // non diffractive events
      registry.add("nonDiff/Stat", "Cut statistics; Selection criterion; Collisions", {HistType::kTH1F, {{20, -0.5, 19.5}}});
      registry.add("nonDiff/cleanFIT", "Statistics of collisions with empty FIT; Multiple of collision time resolution; FIT status; Collisions", {HistType::kTH2F, {{10, -0.5, 9.5}, {2, -0.5, 1.5}}});
      registry.add("nonDiff/Tracks", "Number of tracks; Number of tracks; Collisions", {HistType::kTH1F, {{50, 0.5, 50.5}}});
      registry.add("nonDiff/PVTracks", "Number of PV tracks; Number of PV tracks; Collisions", {HistType::kTH1F, {{50, 0.5, 50.5}}});
      registry.add("nonDiff/globalTracks", "Number of global tracks; Number of global tracks; Collisions", {HistType::kTH1F, {{50, 0.5, 50.5}}});
      registry.add("nonDiff/notPVTracks", "Not PV tracks; Track status bit; Not PV tracks", {HistType::kTH1F, {{17, -0.5, 16.5}}});
      registry.add("nonDiff/tResvsrTOFTracks", "#tResvsrTOFTracks", {HistType::kTH2F, {{1000, 0., 1.E3}, {101, -0.01, 1.01}}});
      registry.add("nonDiff/PVposxy", "Vertex position in x and y direction; V_x; V_y; Collisions", {HistType::kTH2F, {{200, -2., 2.}, {200, -2., 2.}}});
      registry.add("nonDiff/PVposz", "Vertex position in z direction; V_z; Collisions", {HistType::kTH1F, {{1000, -100., 100.}}});
      registry.add("nonDiff/etapt", "eta versus pT of all tracks; eta of track; p_T of track [MeV/c^2]; Tracks", {HistType::kTH2F, {{80, -2., 2.}, {100, 0., 5.}}});
      registry.add("nonDiff/dEdxTPC", "TPC signal versus signed track momentum; Signed track momentum [GeV/c]; TPC signal; Tracks", {HistType::kTH2F, {{100, 0., 5.0}, {3000, 0., 30000.}}});
      registry.add("nonDiff/dEdxTOF", "TOF signal versus signed track momentum; Signed track momentum [GeV/c]; TOF signal; Tracks", {HistType::kTH2F, {{100, 0., 5.0}, {1000, 0., 500000.}}});
      registry.add("nonDiff/PVposxyDG", "DG: Vertex position in x and y direction; V_x; V_y; Collisions", {HistType::kTH2F, {{200, -2., 2.}, {200, -2., 2.}}});
      registry.add("nonDiff/PVposzDG", "DG: Vertex position in z direction; V_z; Collisions", {HistType::kTH1F, {{1000, -100., 100.}}});
      registry.add("nonDiff/etaptDG", "DG: eta versus pT of all tracks; eta of track; p_T of track [MeV/c^2]; Tracks", {HistType::kTH2F, {{80, -2., 2.}, {100, 0., 5.}}});
      registry.add("nonDiff/dEdxTPCDG", "DG: TPC signal versus signed track momentum; Signed track momentum [GeV/c]; TPC signal; Tracks", {HistType::kTH2F, {{100, 0., 5.0}, {3000, 0., 30000.}}});
      registry.add("nonDiff/dEdxTOFDG", "DG: TOF signal versus signed track momentum; Signed track momentum [GeV/c]; TOF signal; Tracks", {HistType::kTH2F, {{100, 0., 5.0}, {1000, 0., 500000.}}});
      registry.add("nonDiff/netChargeDG", "DG: Net charge; Net charge; DG collisions", {HistType::kTH1F, {{21, -10.5, 10.5}}});
      registry.add("nonDiff/IVMptSysDG", "DG: Invariant mass versus p_{T, system}; Invarian mass [GeV/c^2]; p_{T, system} [GeV/c]; DG collisions", {HistType::kTH2F, {{100, 0., 2.5}, {350, 0., 3.5}}});
      registry.add("nonDiff/IVMptTrkDG", "DG: Invariant mass versus p_{T, tracks}; Invarian mass [GeV/c^2]; p_{T, tracks} [GeV/c]; DG collisions", {HistType::kTH2F, {{100, 0., 2.5}, {350, 0., 3.5}}});
      // PYTHIA8 diffractive events
      registry.add("MBRDiff/Stat", "Cut statistics; Selection criterion; Collisions", {HistType::kTH1F, {{20, -0.5, 19.5}}});
      registry.add("MBRDiff/cleanFIT", "Statistics of collisions with empty FIT; Multiple of collision time resolution; FIT status; Collisions", {HistType::kTH2F, {{10, -0.5, 9.5}, {2, -0.5, 1.5}}});
      registry.add("MBRDiff/Tracks", "Number of tracks; Number of tracks; Collisions", {HistType::kTH1F, {{50, 0.5, 50.5}}});
      registry.add("MBRDiff/PVTracks", "Number of PV tracks; Number of PV tracks; Collisions", {HistType::kTH1F, {{50, 0.5, 50.5}}});
      registry.add("MBRDiff/globalTracks", "Number of global tracks; Number of global tracks; Collisions", {HistType::kTH1F, {{50, 0.5, 50.5}}});
      registry.add("MBRDiff/notPVTracks", "Not PV tracks; Track status bit; Not PV tracks", {HistType::kTH1F, {{17, -0.5, 16.5}}});
      registry.add("MBRDiff/tResvsrTOFTracks", "#tResvsrTOFTracks", {HistType::kTH2F, {{1000, 0., 1.E3}, {101, -0.01, 1.01}}});
      registry.add("MBRDiff/PVposxy", "Vertex position in x and y direction; V_x; V_y; Collisions", {HistType::kTH2F, {{200, -2., 2.}, {200, -2., 2.}}});
      registry.add("MBRDiff/PVposz", "Vertex position in z direction; V_z; Collisions", {HistType::kTH1F, {{1000, -100., 100.}}});
      registry.add("MBRDiff/etapt", "eta versus pT of all tracks; eta of track; p_T of track [MeV/c^2]; Tracks", {HistType::kTH2F, {{80, -2., 2.}, {100, 0., 5.}}});
      registry.add("MBRDiff/dEdxTPC", "TPC signal versus signed track momentum; Signed track momentum [GeV/c]; TPC signal; Tracks", {HistType::kTH2F, {{100, 0., 5.0}, {3000, 0., 30000.}}});
      registry.add("MBRDiff/dEdxTOF", "TOF signal versus signed track momentum; Signed track momentum [GeV/c]; TOF signal; Tracks", {HistType::kTH2F, {{100, 0., 5.0}, {1000, 0., 500000.}}});
      registry.add("MBRDiff/PVposxyDG", "DG: Vertex position in x and y direction; V_x; V_y; Collisions", {HistType::kTH2F, {{200, -2., 2.}, {200, -2., 2.}}});
      registry.add("MBRDiff/PVposzDG", "DG: Vertex position in z direction; V_z; Collisions", {HistType::kTH1F, {{1000, -100., 100.}}});
      registry.add("MBRDiff/etaptDG", "DG: eta versus pT of all tracks; eta of track; p_T of track [MeV/c^2]; Tracks", {HistType::kTH2F, {{80, -2., 2.}, {100, 0., 5.}}});
      registry.add("MBRDiff/dEdxTPCDG", "DG: TPC signal versus signed track momentum; Signed track momentum [GeV/c]; TPC signal; Tracks", {HistType::kTH2F, {{100, 0., 5.0}, {3000, 0., 30000.}}});
      registry.add("MBRDiff/dEdxTOFDG", "DG: TOF signal versus signed track momentum; Signed track momentum [GeV/c]; TOF signal; Tracks", {HistType::kTH2F, {{100, 0., 5.0}, {1000, 0., 500000.}}});
      registry.add("MBRDiff/netChargeDG", "DG: Net charge; Net charge; DG collisions", {HistType::kTH1F, {{21, -10.5, 10.5}}});
      registry.add("MBRDiff/IVMptSysDG", "DG: Invariant mass versus p_{T, system}; Invarian mass [GeV/c^2]; p_{T, system} [GeV/c]; DG collisions", {HistType::kTH2F, {{100, 0., 2.5}, {350, 0., 3.5}}});
      registry.add("MBRDiff/IVMptTrkDG", "DG: Invariant mass versus p_{T, tracks}; Invarian mass [GeV/c^2]; p_{T, tracks} [GeV/c]; DG collisions", {HistType::kTH2F, {{100, 0., 2.5}, {350, 0., 3.5}}});
      // GRANIITTI diffractive events
      registry.add("Diff/Stat", "Cut statistics; Selection criterion; Collisions", {HistType::kTH1F, {{20, -0.5, 19.5}}});
      registry.add("Diff/cleanFIT", "Statistics of collisions with empty FIT; Multiple of collision time resolution; FIT status; Collisions", {HistType::kTH2F, {{10, -0.5, 9.5}, {2, -0.5, 1.5}}});
      registry.add("Diff/Tracks", "Number of tracks; Number of tracks; Collisions", {HistType::kTH1F, {{50, 0.5, 50.5}}});
      registry.add("Diff/PVTracks", "Number of PV tracks; Number of PV tracks; Collisions", {HistType::kTH1F, {{50, 0.5, 50.5}}});
      registry.add("Diff/globalTracks", "Number of global tracks; Number of global tracks; Collisions", {HistType::kTH1F, {{50, 0.5, 50.5}}});
      registry.add("Diff/notPVTracks", "Not PV tracks; Track status bit; Not PV tracks", {HistType::kTH1F, {{17, -0.5, 16.5}}});
      registry.add("Diff/tResvsrTOFTracks", "#tResvsrTOFTracks", {HistType::kTH2F, {{1000, 0., 1.E3}, {101, -0.01, 1.01}}});
      registry.add("Diff/PVposxy", "Vertex position in x and y direction; V_x; V_y; Collisions", {HistType::kTH2F, {{200, -2., 2.}, {200, -2., 2.}}});
      registry.add("Diff/PVposz", "Vertex position in z direction; V_z; Collisions", {HistType::kTH1F, {{1000, -100., 100.}}});
      registry.add("Diff/etapt", "eta versus pT of all tracks; eta of track; p_T of track [MeV/c^2]; Tracks", {HistType::kTH2F, {{80, -2., 2.}, {100, 0., 5.}}});
      registry.add("Diff/dEdxTPC", "TPC signal versus signed track momentum; Signed track momentum [GeV/c]; TPC signal; Tracks", {HistType::kTH2F, {{100, 0., 5.0}, {3000, 0., 30000.}}});
      registry.add("Diff/dEdxTOF", "TOF signal versus signed track momentum; Signed track momentum [GeV/c]; TOF signal; Tracks", {HistType::kTH2F, {{100, 0., 5.0}, {1000, 0., 500000.}}});
      registry.add("Diff/PVposxyDG", "DG: Vertex position in x and y direction; V_x; V_y; Collisions", {HistType::kTH2F, {{200, -2., 2.}, {200, -2., 2.}}});
      registry.add("Diff/PVposzDG", "DG: Vertex position in z direction; V_z; Collisions", {HistType::kTH1F, {{1000, -100., 100.}}});
      registry.add("Diff/etaptDG", "DG: eta versus pT of all tracks; eta of track; p_T of track [MeV/c^2]; Tracks", {HistType::kTH2F, {{80, -2., 2.}, {100, 0., 5.}}});
      registry.add("Diff/dEdxTPCDG", "DG: TPC signal versus signed track momentum; Signed track momentum [GeV/c]; TPC signal; Tracks", {HistType::kTH2F, {{100, 0., 5.0}, {3000, 0., 30000.}}});
      registry.add("Diff/dEdxTOFDG", "DG: TOF signal versus signed track momentum; Signed track momentum [GeV/c]; TOF signal; Tracks", {HistType::kTH2F, {{100, 0., 5.0}, {1000, 0., 500000.}}});
      registry.add("Diff/netChargeDG", "DG: Net charge; Net charge; DG collisions", {HistType::kTH1F, {{21, -10.5, 10.5}}});
      registry.add("Diff/IVMptSysDG", "DG: Invariant mass versus p_{T, system}; Invarian mass [GeV/c^2]; p_{T, system} [GeV/c]; DG collisions", {HistType::kTH2F, {{100, 0., 2.5}, {350, 0., 3.5}}});
      registry.add("Diff/IVMptTrkDG", "DG: Invariant mass versus p_{T, tracks}; Invarian mass [GeV/c^2]; p_{T, tracks} [GeV/c]; DG collisions", {HistType::kTH2F, {{100, 0., 2.5}, {350, 0., 3.5}}});
    }

    if (context.mOptions.get<bool>("processFV0")) {
      registry.add("FV0/FV0Aamp", "Amplitudes of FV0A channels; FV0A amplitude; FV0A channel; Entries", {HistType::kTH2F, {{48, -0.5, 47.5}, {1000, 0., 1000.}}});
    }

    if (context.mOptions.get<bool>("processFT0")) {
      registry.add("FT0/FT0Aamp", "Amplitudes of FT0A channels; FT0A amplitude; FT0A channel; Entries", {HistType::kTH2F, {{96, -0.5, 95.5}, {100, 0., 200.}}});
      registry.add("FT0/FT0Camp", "Amplitudes of FT0C channels; FT0C amplitude; FT0C channel; Entries", {HistType::kTH2F, {{112, -0.5, 111.5}, {100, 0., 200.}}});
    }

    if (context.mOptions.get<bool>("processFDD")) {
      registry.add("FDD/FDDAamp", "Amplitudes of FDDA channels; FDDA amplitude; FDDA channel; Entries", {HistType::kTH2F, {{8, -0.5, 7.5}, {100, 0., 100.}}});
      registry.add("FDD/FDDCamp", "Amplitudes of FDDC channels; FDDC amplitude; FDDC channel; Entries", {HistType::kTH2F, {{8, -0.5, 7.5}, {100, 0., 100.}}});
    }

    if (context.mOptions.get<bool>("processZDC")) {
      registry.add("ZdcEnergies", "Various ZDC energies; Energy type; ZDC energy; Entries", {HistType::kTH2F, {{22, -0.5, 21.5}, {100, 0., 1000.}}});
    }

    if (context.mOptions.get<bool>("processCalo")) {
      registry.add("CaloCell", "Calorimeter cell; Calorimeter cell number; Entries", {HistType::kTH1I, {{18000, -0.5, 17999.5}}});
      registry.add("CaloAmplitude", "Calorimeter amplitudes; Calorimeter amplitude; Entries", {HistType::kTH1F, {{100, 0, 10.}}});
    }
  }

  void run(ProcessingContext& pc)
  {
    // get ambiguous tracks table
    auto t1 = pc.inputs().get<TableConsumer>("BCs")->asArrowTable();
    auto t2 = pc.inputs().get<TableConsumer>("BcSels")->asArrowTable();
    auto t3 = pc.inputs().get<TableConsumer>("Run3MatchedToBCSparse")->asArrowTable();
    auto bcs = BCs({t1, t2, t3});

    if (withAmbTrackAnalysis) {
      auto t4 = pc.inputs().get<TableConsumer>("AmbiguousTracks")->asArrowTable();
      auto ambtracks = ATs({t4});
      ambtracks.bindExternalIndices(&bcs);

      // make sorted list of BC ranges which are associated with an ambiguous track.
      // This is used to efficiently check whether a given BC is contained in one of these ranges
      abcrs.reset();
      LOGF(info, "<DiffMCQA> size of ambiguous tracks table %i", ambtracks.size());
      for (auto const& ambtrack : ambtracks) {
        auto bcfirst = ambtrack.bc().rawIteratorAt(0);
        auto bclast = ambtrack.bc().rawIteratorAt(ambtrack.bc().size() - 1);
        abcrs.add(bcfirst.globalIndex(), bclast.globalIndex());
      }
      abcrs.merge();
    }

    if (withAmbFwdTrackAnalysis) {
      // get ambiguous FwdTracks table
      auto t5 = pc.inputs().get<TableConsumer>("AmbiguousFwdTracks")->asArrowTable();
      auto ambfwdtracks = AFTs({t5});
      ambfwdtracks.bindExternalIndices(&bcs);

      // make sorted list of BC ranges which are associated with an ambiguous FwdTrack.
      afbcrs.reset();
      LOGF(info, "<DiffMCQA> size of ambiguous fwd tracks table %i", ambfwdtracks.size());
      for (auto const& ambfwdtrack : ambfwdtracks) {
        auto bcfirst = ambfwdtrack.bc().rawIteratorAt(0);
        auto bclast = ambfwdtrack.bc().rawIteratorAt(ambfwdtrack.bc().size() - 1);
        afbcrs.add(bcfirst.globalIndex(), bclast.globalIndex());
      }
      afbcrs.merge();
    }
    LOGF(info, "<DiffMCQA> Size of abcrs %i and afbcrs %i", abcrs.size(), afbcrs.size());
  }

  Preslice<aod::McParticles> partPerMcCollision = aod::mcparticle::mcCollisionId;

  // ...............................................................................................................
  void processMCTruth(aod::McCollision const& mccollision, soa::SmallGroups<CCs> const& collisions, aod::McParticles const& McParts)
  {
    LOGF(info, "Number of collisions %d", collisions.size());
    LOGF(info, "Number of McParts %d", McParts.size());

    // is this a central diffractive event?
    // by default it is assumed to be a MB event
    bool isPythiaDiff = false;
    bool isGraniittiDiff = false;
    registry.get<TH1>(HIST("MCTruth/Stat"))->Fill(1., 1.);
    isPythiaDiff = udhelpers::isPythiaCDE(McParts);
    isGraniittiDiff = udhelpers::isGraniittiCDE(McParts);
    if (isPythiaDiff) {
      registry.get<TH1>(HIST("MCTruth/Stat"))->Fill(3., 1.);
    } else if (isGraniittiDiff) {
      registry.get<TH1>(HIST("MCTruth/Stat"))->Fill(4., 1.);
    } else {
      registry.get<TH1>(HIST("MCTruth/Stat"))->Fill(2., 1.);
    }

    // number of reconstructed collision
    registry.get<TH1>(HIST("MCTruth/recCols"))->Fill(collisions.size(), 1.);
  }
  PROCESS_SWITCH(DiffMCQA, processMCTruth, "Process MC truth", true);

  // ...............................................................................................................
  void processMain(CC const& collision, BCs const& bct0s,
                   TCs const& tracks, FWs const& fwdtracks, ATs const& ambtracks, AFTs const& ambfwdtracks,
                   aod::FT0s const& ft0s, aod::FV0As const& fv0as, aod::FDDs const& fdds,
                   aod::Zdcs& zdcs, aod::Calos& calos,
                   aod::V0s const& v0s, aod::Cascades const& cascades,
                   aod::McCollisions const& McCols, aod::McParticles const& McParts)
  {
    bool isDGcandidate = true;

    // is this a central diffractive event?
    // by default it is assumed to be a MB event
    bool isPythiaDiff = false;
    bool isGraniittiDiff = false;
    if (collision.has_mcCollision()) {
      registry.get<TH1>(HIST("all/mcCols"))->Fill(1., 1.);

      auto MCCol = collision.mcCollision();
      auto MCPartSlice = McParts.sliceBy(partPerMcCollision, MCCol.globalIndex());
      isPythiaDiff = udhelpers::isPythiaCDE(MCPartSlice);
      isGraniittiDiff = udhelpers::isGraniittiCDE(MCPartSlice);
    } else {
      registry.get<TH1>(HIST("all/mcCols"))->Fill(0., 1.);
    }

    // global tracks
    Partition<TCs> goodTracks = requireGlobalTrackInFilter();
    goodTracks.bindTable(tracks);

    // update collision histograms
    if (isPythiaDiff) {
      registry.get<TH1>(HIST("MBRDiff/Stat"))->Fill(0., 1.);
      registry.get<TH2>(HIST("MBRDiff/PVposxy"))->Fill(collision.posX(), collision.posY());
      registry.get<TH1>(HIST("MBRDiff/PVposz"))->Fill(collision.posZ());
      registry.get<TH1>(HIST("MBRDiff/Tracks"))->Fill(tracks.size());
      registry.get<TH1>(HIST("MBRDiff/PVTracks"))->Fill(collision.numContrib());
      registry.get<TH1>(HIST("MBRDiff/globalTracks"))->Fill(goodTracks.size());
    } else if (isGraniittiDiff) {
      registry.get<TH1>(HIST("Diff/Stat"))->Fill(0., 1.);
      registry.get<TH2>(HIST("Diff/PVposxy"))->Fill(collision.posX(), collision.posY());
      registry.get<TH1>(HIST("Diff/PVposz"))->Fill(collision.posZ());
      registry.get<TH1>(HIST("Diff/Tracks"))->Fill(tracks.size());
      registry.get<TH1>(HIST("Diff/PVTracks"))->Fill(collision.numContrib());
      registry.get<TH1>(HIST("Diff/globalTracks"))->Fill(goodTracks.size());
    } else {
      registry.get<TH1>(HIST("nonDiff/Stat"))->Fill(0., 1.);
      registry.get<TH2>(HIST("nonDiff/PVposxy"))->Fill(collision.posX(), collision.posY());
      registry.get<TH1>(HIST("nonDiff/PVposz"))->Fill(collision.posZ());
      registry.get<TH1>(HIST("nonDiff/Tracks"))->Fill(tracks.size());
      registry.get<TH1>(HIST("nonDiff/PVTracks"))->Fill(collision.numContrib());
      registry.get<TH1>(HIST("nonDiff/globalTracks"))->Fill(goodTracks.size());
    }

    // test influence of BCrange width
    for (int minBcs = 0; minBcs < 10; minBcs++) {
      auto bcSlice = udhelpers::MCcompatibleBCs(collision, 0, bct0s, minBcs);
      isDGcandidate = true;
      for (auto& bc : bcSlice) {
        isDGcandidate &= udhelpers::cleanFIT(bc, diffCuts.maxFITtime(), diffCuts.FITAmpLimits());
      }
      if (isPythiaDiff) {
        registry.get<TH2>(HIST("MBRDiff/cleanFIT"))->Fill(minBcs, isDGcandidate * 1.);
      } else if (isGraniittiDiff) {
        registry.get<TH2>(HIST("Diff/cleanFIT"))->Fill(minBcs, isDGcandidate * 1.);
      } else {
        registry.get<TH2>(HIST("nonDiff/cleanFIT"))->Fill(minBcs, isDGcandidate * 1.);
      }
    }

    // loop over all tracks
    float rgtrwTOF = 0.;
    for (auto const& track : tracks) {
      // update eta vs pt and dEdx histograms histogram
      if (isPythiaDiff) {
        registry.get<TH2>(HIST("MBRDiff/etapt"))->Fill(track.eta(), track.pt());
        registry.get<TH2>(HIST("MBRDiff/dEdxTPC"))->Fill(track.p(), track.tpcSignal());
      } else if (isGraniittiDiff) {
        registry.get<TH2>(HIST("Diff/etapt"))->Fill(track.eta(), track.pt());
        registry.get<TH2>(HIST("Diff/dEdxTPC"))->Fill(track.p(), track.tpcSignal());
      } else {
        registry.get<TH2>(HIST("nonDiff/etapt"))->Fill(track.eta(), track.pt());
        registry.get<TH2>(HIST("nonDiff/dEdxTPC"))->Fill(track.p(), track.tpcSignal());
      }
      if (track.tpcSignal() > maxdEdxTPC) {
        maxdEdxTPC = track.tpcSignal();
        LOGF(info, "<DiffMCQA> New maxdEdx TPC %f", maxdEdxTPC);
      }

      // TOF hit?
      if (track.hasTOF()) {
        if (isPythiaDiff) {
          registry.get<TH2>(HIST("MBRDiff/dEdxTOF"))->Fill(track.p(), track.tofSignal());
        } else if (isGraniittiDiff) {
          registry.get<TH2>(HIST("Diff/dEdxTOF"))->Fill(track.p(), track.tofSignal());
        } else {
          registry.get<TH2>(HIST("nonDiff/dEdxTOF"))->Fill(track.p(), track.tofSignal());
        }
        if (track.tofSignal() > maxdEdxTOF) {
          maxdEdxTOF = track.tofSignal();
          LOGF(info, "<DiffMCQA> New maxdEdx TOF %f", maxdEdxTOF);
        }

        // vertex track with TOF hit?
        if (track.isPVContributor()) {
          rgtrwTOF += 1.;
        }
      }
    }
    if (collision.numContrib() > 0) {
      rgtrwTOF /= collision.numContrib();
    }
    LOGF(debug, "<DiffMCQA> Vertex tracks with TOF: %f [1]", rgtrwTOF);
    if (isPythiaDiff) {
      registry.get<TH2>(HIST("MBRDiff/tResvsrTOFTracks"))->Fill(collision.collisionTimeRes(), rgtrwTOF);
    } else if (isGraniittiDiff) {
      registry.get<TH2>(HIST("Diff/tResvsrTOFTracks"))->Fill(collision.collisionTimeRes(), rgtrwTOF);
    } else {
      registry.get<TH2>(HIST("nonDiff/tResvsrTOFTracks"))->Fill(collision.collisionTimeRes(), rgtrwTOF);
    }

    // is it a DG candidate?
    // DG = no FIT signal in compatible BCs
    //    & no ZDC signal in compatible BCs
    //    & no Calo signal in compatible BCs
    //    & no V0s
    //    & no Cascades
    //    & number of forward tracks = 0
    //    & no global track which is not a vertex track
    //    & ntrMin <= number of vertex tracks <= ntrMax
    isDGcandidate = true;

    // get BCrange to test for FIT signals
    auto bcSlice = udhelpers::MCcompatibleBCs(collision, diffCuts.NDtcoll(), bct0s, diffCuts.minNBCs());

    // no FIT signal in bcSlice / collision
    if (doCleanFITBC) {
      for (auto const& bc : bcSlice) {
        if (!udhelpers::cleanFIT(bc, diffCuts.maxFITtime(), diffCuts.FITAmpLimits())) {
          isDGcandidate = false;
          break;
        }
      }
    } else {
      if (!udhelpers::cleanFITCollision(collision, diffCuts.maxFITtime(), diffCuts.FITAmpLimits())) {
        isDGcandidate = false;
      }
    }
    if (isPythiaDiff) {
      registry.get<TH1>(HIST("MBRDiff/Stat"))->Fill(1., isDGcandidate * 1.);
    } else if (isGraniittiDiff) {
      registry.get<TH1>(HIST("Diff/Stat"))->Fill(1., isDGcandidate * 1.);
    } else {
      registry.get<TH1>(HIST("nonDiff/Stat"))->Fill(1., isDGcandidate * 1.);
    }

    // no Zdc signal in bcSlice
    std::vector<float> lims(10, 0.);
    for (auto const& bc : bcSlice) {
      if (!udhelpers::cleanZDC(bc, zdcs, lims, cache)) {
        isDGcandidate = false;
        break;
      }
    }
    if (isPythiaDiff) {
      registry.get<TH1>(HIST("MBRDiff/Stat"))->Fill(2., isDGcandidate * 1.);
    } else if (isGraniittiDiff) {
      registry.get<TH1>(HIST("Diff/Stat"))->Fill(2., isDGcandidate * 1.);
    } else {
      registry.get<TH1>(HIST("nonDiff/Stat"))->Fill(2., isDGcandidate * 1.);
    }

    // no Calo signal in bcSlice
    for (auto const& bc : bcSlice) {
      if (!udhelpers::cleanCalo(bc, calos, lims, cache)) {
        isDGcandidate = false;
        break;
      }
    }
    if (isPythiaDiff) {
      registry.get<TH1>(HIST("MBRDiff/Stat"))->Fill(3., isDGcandidate * 1.);
    } else if (isGraniittiDiff) {
      registry.get<TH1>(HIST("Diff/Stat"))->Fill(3., isDGcandidate * 1.);
    } else {
      registry.get<TH1>(HIST("nonDiff/Stat"))->Fill(3., isDGcandidate * 1.);
    }

    // no V0s
    isDGcandidate &= (v0s.size() == 0);
    if (isPythiaDiff) {
      registry.get<TH1>(HIST("MBRDiff/Stat"))->Fill(4., isDGcandidate * 1.);
    } else if (isGraniittiDiff) {
      registry.get<TH1>(HIST("Diff/Stat"))->Fill(4., isDGcandidate * 1.);
    } else {
      registry.get<TH1>(HIST("nonDiff/Stat"))->Fill(4., isDGcandidate * 1.);
    }

    // no Cascades
    isDGcandidate &= (cascades.size() == 0);
    if (isPythiaDiff) {
      registry.get<TH1>(HIST("MBRDiff/Stat"))->Fill(5., isDGcandidate * 1.);
    } else if (isGraniittiDiff) {
      registry.get<TH1>(HIST("Diff/Stat"))->Fill(5., isDGcandidate * 1.);
    } else {
      registry.get<TH1>(HIST("nonDiff/Stat"))->Fill(5., isDGcandidate * 1.);
    }

    // number of forward tracks = 0
    isDGcandidate &= (fwdtracks.size() == 0);
    if (isPythiaDiff) {
      registry.get<TH1>(HIST("MBRDiff/Stat"))->Fill(6., isDGcandidate * 1.);
    } else if (isGraniittiDiff) {
      registry.get<TH1>(HIST("Diff/Stat"))->Fill(6., isDGcandidate * 1.);
    } else {
      registry.get<TH1>(HIST("nonDiff/Stat"))->Fill(6., isDGcandidate * 1.);
    }

    // no global tracks which are no vtx tracks
    bool globalAndVtx = isDGcandidate;
    bool vtxAndGlobal = isDGcandidate;
    for (auto const& track : tracks) {
      if (track.isGlobalTrack() && !track.isPVContributor()) {
        globalAndVtx = false;
      }
      if (track.isPVContributor() && !track.isGlobalTrack()) {
        vtxAndGlobal = false;
      }
    }
    if (isPythiaDiff) {
      registry.get<TH1>(HIST("MBRDiff/Stat"))->Fill(7., globalAndVtx * 1.);
      registry.get<TH1>(HIST("MBRDiff/Stat"))->Fill(8., vtxAndGlobal * 1.);
    } else if (isGraniittiDiff) {
      registry.get<TH1>(HIST("Diff/Stat"))->Fill(7., globalAndVtx * 1.);
      registry.get<TH1>(HIST("Diff/Stat"))->Fill(8., vtxAndGlobal * 1.);
    } else {
      registry.get<TH1>(HIST("nonDiff/Stat"))->Fill(7., globalAndVtx * 1.);
      registry.get<TH1>(HIST("nonDiff/Stat"))->Fill(8., vtxAndGlobal * 1.);
    }
    isDGcandidate &= globalAndVtx;
    if (diffCuts.globalTracksOnly()) {
      isDGcandidate &= vtxAndGlobal;
    }

    // check a given bc for possible ambiguous Tracks
    auto noAmbTracks = isDGcandidate;
    for (auto const& bc : bcSlice) {
      if (abcrs.isInRange(bc.globalIndex())) {
        noAmbTracks = false;
        break;
      }
    }
    if (isPythiaDiff) {
      registry.get<TH1>(HIST("MBRDiff/Stat"))->Fill(9., noAmbTracks * 1.);
    } else if (isGraniittiDiff) {
      registry.get<TH1>(HIST("Diff/Stat"))->Fill(9., noAmbTracks * 1.);
    } else {
      registry.get<TH1>(HIST("nonDiff/Stat"))->Fill(9., noAmbTracks * 1.);
    }

    // check a given bc for possible ambiguous FwdTracks
    auto noAmbFwdTracks = isDGcandidate;
    for (auto const& bc : bcSlice) {
      if (afbcrs.isInRange(bc.globalIndex())) {
        noAmbFwdTracks = false;
        break;
      }
    }
    if (isPythiaDiff) {
      registry.get<TH1>(HIST("MBRDiff/Stat"))->Fill(10., noAmbFwdTracks * 1.);
    } else if (isGraniittiDiff) {
      registry.get<TH1>(HIST("Diff/Stat"))->Fill(10., noAmbFwdTracks * 1.);
    } else {
      registry.get<TH1>(HIST("nonDiff/Stat"))->Fill(10., noAmbFwdTracks * 1.);
    }

    // at least one vtx track with TOF hit
    isDGcandidate &= (rgtrwTOF >= diffCuts.minRgtrwTOF());
    if (isPythiaDiff) {
      registry.get<TH1>(HIST("MBRDiff/Stat"))->Fill(11., isDGcandidate * 1.);
    } else if (isGraniittiDiff) {
      registry.get<TH1>(HIST("Diff/Stat"))->Fill(11., isDGcandidate * 1.);
    } else {
      registry.get<TH1>(HIST("nonDiff/Stat"))->Fill(11., isDGcandidate * 1.);
    }

    // number of vertex tracks <= n
    isDGcandidate &= (collision.numContrib() >= diffCuts.minNTracks());
    if (isPythiaDiff) {
      registry.get<TH1>(HIST("MBRDiff/Stat"))->Fill(12., isDGcandidate * 1.);
    } else if (isGraniittiDiff) {
      registry.get<TH1>(HIST("Diff/Stat"))->Fill(12., isDGcandidate * 1.);
    } else {
      registry.get<TH1>(HIST("nonDiff/Stat"))->Fill(12., isDGcandidate * 1.);
    }
    isDGcandidate &= (collision.numContrib() <= diffCuts.maxNTracks());
    if (isPythiaDiff) {
      registry.get<TH1>(HIST("MBRDiff/Stat"))->Fill(13., isDGcandidate * 1.);
    } else if (isGraniittiDiff) {
      registry.get<TH1>(HIST("Diff/Stat"))->Fill(13., isDGcandidate * 1.);
    } else {
      registry.get<TH1>(HIST("nonDiff/Stat"))->Fill(13., isDGcandidate * 1.);
    }

    // net charge and invariant mass
    bool goodetas = true;
    bool goodpts = true;
    bool goodnchs = true;
    auto netCharge = 0;
    auto lvtmp = TLorentzVector();
    auto ivm = TLorentzVector();
    if (isDGcandidate) {

      // which particle hypothesis?
      auto mass2Use = constants::physics::MassPionCharged;
      if (diffCuts.pidHypothesis() == 321) {
        mass2Use = constants::physics::MassKaonCharged;
      }

      // check also pt and eta of tracks
      for (auto& track : tracks) {
        // update histogram notPVTracks
        if (!track.isPVContributor()) {
          if (isPythiaDiff) {
            registry.get<TH1>(HIST("MBRDiff/notPVTracks"))->Fill(0., 1.);
            registry.get<TH1>(HIST("MBRDiff/notPVTracks"))->Fill(1., track.isGlobalTrackSDD() * 1.);
            registry.get<TH1>(HIST("MBRDiff/notPVTracks"))->Fill(2., track.passedTrackType() * 1.);
            registry.get<TH1>(HIST("MBRDiff/notPVTracks"))->Fill(3., track.passedPtRange() * 1.);
            registry.get<TH1>(HIST("MBRDiff/notPVTracks"))->Fill(4., track.passedEtaRange() * 1.);
            registry.get<TH1>(HIST("MBRDiff/notPVTracks"))->Fill(5., track.passedTPCNCls() * 1.);
            registry.get<TH1>(HIST("MBRDiff/notPVTracks"))->Fill(6., track.passedTPCCrossedRows() * 1.);
            registry.get<TH1>(HIST("MBRDiff/notPVTracks"))->Fill(7., track.passedTPCCrossedRowsOverNCls() * 1.);
            registry.get<TH1>(HIST("MBRDiff/notPVTracks"))->Fill(8., track.passedTPCChi2NDF() * 1.);
            registry.get<TH1>(HIST("MBRDiff/notPVTracks"))->Fill(9., track.passedTPCRefit() * 1.);
            registry.get<TH1>(HIST("MBRDiff/notPVTracks"))->Fill(10., track.passedITSNCls() * 1.);
            registry.get<TH1>(HIST("MBRDiff/notPVTracks"))->Fill(11., track.passedITSChi2NDF() * 1.);
            registry.get<TH1>(HIST("MBRDiff/notPVTracks"))->Fill(12., track.passedITSRefit() * 1.);
            registry.get<TH1>(HIST("MBRDiff/notPVTracks"))->Fill(13., track.passedITSHits() * 1.);
            registry.get<TH1>(HIST("MBRDiff/notPVTracks"))->Fill(14., track.passedGoldenChi2() * 1.);
            registry.get<TH1>(HIST("MBRDiff/notPVTracks"))->Fill(15., track.passedDCAxy() * 1.);
            registry.get<TH1>(HIST("MBRDiff/notPVTracks"))->Fill(16., track.passedDCAz() * 1.);
          } else if (isGraniittiDiff) {
            registry.get<TH1>(HIST("Diff/notPVTracks"))->Fill(0., 1.);
            registry.get<TH1>(HIST("Diff/notPVTracks"))->Fill(1., track.isGlobalTrackSDD() * 1.);
            registry.get<TH1>(HIST("Diff/notPVTracks"))->Fill(2., track.passedTrackType() * 1.);
            registry.get<TH1>(HIST("Diff/notPVTracks"))->Fill(3., track.passedPtRange() * 1.);
            registry.get<TH1>(HIST("Diff/notPVTracks"))->Fill(4., track.passedEtaRange() * 1.);
            registry.get<TH1>(HIST("Diff/notPVTracks"))->Fill(5., track.passedTPCNCls() * 1.);
            registry.get<TH1>(HIST("Diff/notPVTracks"))->Fill(6., track.passedTPCCrossedRows() * 1.);
            registry.get<TH1>(HIST("Diff/notPVTracks"))->Fill(7., track.passedTPCCrossedRowsOverNCls() * 1.);
            registry.get<TH1>(HIST("Diff/notPVTracks"))->Fill(8., track.passedTPCChi2NDF() * 1.);
            registry.get<TH1>(HIST("Diff/notPVTracks"))->Fill(9., track.passedTPCRefit() * 1.);
            registry.get<TH1>(HIST("Diff/notPVTracks"))->Fill(10., track.passedITSNCls() * 1.);
            registry.get<TH1>(HIST("Diff/notPVTracks"))->Fill(11., track.passedITSChi2NDF() * 1.);
            registry.get<TH1>(HIST("Diff/notPVTracks"))->Fill(12., track.passedITSRefit() * 1.);
            registry.get<TH1>(HIST("Diff/notPVTracks"))->Fill(13., track.passedITSHits() * 1.);
            registry.get<TH1>(HIST("Diff/notPVTracks"))->Fill(14., track.passedGoldenChi2() * 1.);
            registry.get<TH1>(HIST("Diff/notPVTracks"))->Fill(15., track.passedDCAxy() * 1.);
            registry.get<TH1>(HIST("Diff/notPVTracks"))->Fill(16., track.passedDCAz() * 1.);
          } else {
            registry.get<TH1>(HIST("nonDiff/notPVTracks"))->Fill(0., 1.);
            registry.get<TH1>(HIST("nonDiff/notPVTracks"))->Fill(1., track.isGlobalTrackSDD() * 1.);
            registry.get<TH1>(HIST("nonDiff/notPVTracks"))->Fill(2., track.passedTrackType() * 1.);
            registry.get<TH1>(HIST("nonDiff/notPVTracks"))->Fill(3., track.passedPtRange() * 1.);
            registry.get<TH1>(HIST("nonDiff/notPVTracks"))->Fill(4., track.passedEtaRange() * 1.);
            registry.get<TH1>(HIST("nonDiff/notPVTracks"))->Fill(5., track.passedTPCNCls() * 1.);
            registry.get<TH1>(HIST("nonDiff/notPVTracks"))->Fill(6., track.passedTPCCrossedRows() * 1.);
            registry.get<TH1>(HIST("nonDiff/notPVTracks"))->Fill(7., track.passedTPCCrossedRowsOverNCls() * 1.);
            registry.get<TH1>(HIST("nonDiff/notPVTracks"))->Fill(8., track.passedTPCChi2NDF() * 1.);
            registry.get<TH1>(HIST("nonDiff/notPVTracks"))->Fill(9., track.passedTPCRefit() * 1.);
            registry.get<TH1>(HIST("nonDiff/notPVTracks"))->Fill(10., track.passedITSNCls() * 1.);
            registry.get<TH1>(HIST("nonDiff/notPVTracks"))->Fill(11., track.passedITSChi2NDF() * 1.);
            registry.get<TH1>(HIST("nonDiff/notPVTracks"))->Fill(12., track.passedITSRefit() * 1.);
            registry.get<TH1>(HIST("nonDiff/notPVTracks"))->Fill(13., track.passedITSHits() * 1.);
            registry.get<TH1>(HIST("nonDiff/notPVTracks"))->Fill(14., track.passedGoldenChi2() * 1.);
            registry.get<TH1>(HIST("nonDiff/notPVTracks"))->Fill(15., track.passedDCAxy() * 1.);
            registry.get<TH1>(HIST("nonDiff/notPVTracks"))->Fill(16., track.passedDCAz() * 1.);
          }
          continue;
        }

        lvtmp.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), mass2Use);
        LOGF(debug, "mass %f track pt %f/%f eta %f/%f", mass2Use, track.pt(), lvtmp.Perp(), track.eta(), lvtmp.Eta());
        if (track.pt() <= diffCuts.minPt() || track.pt() >= diffCuts.maxPt()) {
          goodpts = false;
        }
        if (track.eta() <= diffCuts.minEta() || track.eta() >= diffCuts.maxEta()) {
          goodetas = false;
        }
        netCharge += track.sign();
        ivm += lvtmp;
      }
    }
    isDGcandidate &= goodpts;
    if (isPythiaDiff) {
      registry.get<TH1>(HIST("MBRDiff/Stat"))->Fill(14., isDGcandidate * 1.);
    } else if (isGraniittiDiff) {
      registry.get<TH1>(HIST("Diff/Stat"))->Fill(14., isDGcandidate * 1.);
    } else {
      registry.get<TH1>(HIST("nonDiff/Stat"))->Fill(14., isDGcandidate * 1.);
    }
    isDGcandidate &= goodetas;
    if (isPythiaDiff) {
      registry.get<TH1>(HIST("MBRDiff/Stat"))->Fill(15., isDGcandidate * 1.);
    } else if (isGraniittiDiff) {
      registry.get<TH1>(HIST("Diff/Stat"))->Fill(15., isDGcandidate * 1.);
    } else {
      registry.get<TH1>(HIST("nonDiff/Stat"))->Fill(15., isDGcandidate * 1.);
    }
    auto netChargeValues = diffCuts.netCharges();
    if (std::find(netChargeValues.begin(), netChargeValues.end(), netCharge) == netChargeValues.end()) {
      goodnchs = false;
    }
    isDGcandidate &= goodnchs;
    if (isPythiaDiff) {
      registry.get<TH1>(HIST("MBRDiff/Stat"))->Fill(16., isDGcandidate * 1.);
    } else if (isGraniittiDiff) {
      registry.get<TH1>(HIST("Diff/Stat"))->Fill(16., isDGcandidate * 1.);
    } else {
      registry.get<TH1>(HIST("nonDiff/Stat"))->Fill(16., isDGcandidate * 1.);
    }
    if (isPythiaDiff) {
      registry.get<TH1>(HIST("MBRDiff/Stat"))->Fill(17., isDGcandidate * 1.);
    } else if (isGraniittiDiff) {
      registry.get<TH1>(HIST("Diff/Stat"))->Fill(17., isDGcandidate * 1.);
    } else {
      registry.get<TH1>(HIST("nonDiff/Stat"))->Fill(17., isDGcandidate * 1.);
    }
    isDGcandidate &= (ivm.M() >= diffCuts.minIVM());
    if (isPythiaDiff) {
      registry.get<TH1>(HIST("MBRDiff/Stat"))->Fill(18., isDGcandidate * 1.);
    } else if (isGraniittiDiff) {
      registry.get<TH1>(HIST("Diff/Stat"))->Fill(18., isDGcandidate * 1.);
    } else {
      registry.get<TH1>(HIST("nonDiff/Stat"))->Fill(18., isDGcandidate * 1.);
    }
    isDGcandidate &= (ivm.M() <= diffCuts.maxIVM());
    if (isPythiaDiff) {
      registry.get<TH1>(HIST("MBRDiff/Stat"))->Fill(19., isDGcandidate * 1.);
    } else if (isGraniittiDiff) {
      registry.get<TH1>(HIST("Diff/Stat"))->Fill(19., isDGcandidate * 1.);
    } else {
      registry.get<TH1>(HIST("nonDiff/Stat"))->Fill(19., isDGcandidate * 1.);
    }

    // update some DG histograms
    if (isDGcandidate) {
      if (isPythiaDiff) {
        registry.get<TH2>(HIST("MBRDiff/PVposxyDG"))->Fill(collision.posX(), collision.posY());
        registry.get<TH1>(HIST("MBRDiff/PVposzDG"))->Fill(collision.posZ());
        registry.get<TH1>(HIST("MBRDiff/netChargeDG"))->Fill(netCharge);
        registry.get<TH2>(HIST("MBRDiff/IVMptSysDG"))->Fill(ivm.M(), ivm.Perp());
      } else if (isGraniittiDiff) {
        registry.get<TH2>(HIST("Diff/PVposxyDG"))->Fill(collision.posX(), collision.posY());
        registry.get<TH1>(HIST("Diff/PVposzDG"))->Fill(collision.posZ());
        registry.get<TH1>(HIST("Diff/netChargeDG"))->Fill(netCharge);
        registry.get<TH2>(HIST("Diff/IVMptSysDG"))->Fill(ivm.M(), ivm.Perp());
      } else {
        registry.get<TH2>(HIST("nonDiff/PVposxyDG"))->Fill(collision.posX(), collision.posY());
        registry.get<TH1>(HIST("nonDiff/PVposzDG"))->Fill(collision.posZ());
        registry.get<TH1>(HIST("nonDiff/netChargeDG"))->Fill(netCharge);
        registry.get<TH2>(HIST("nonDiff/IVMptSysDG"))->Fill(ivm.M(), ivm.Perp());
      }

      // fill dEdx of DG event tracks
      for (auto& track : tracks) {
        if (track.isPVContributor()) {
          if (isPythiaDiff) {
            registry.get<TH2>(HIST("MBRDiff/etaptDG"))->Fill(track.eta(), track.pt());
            registry.get<TH2>(HIST("MBRDiff/dEdxTPCDG"))->Fill(track.p(), track.tpcSignal());
            registry.get<TH2>(HIST("MBRDiff/IVMptTrkDG"))->Fill(ivm.M(), track.pt());
          } else if (isGraniittiDiff) {
            registry.get<TH2>(HIST("Diff/etaptDG"))->Fill(track.eta(), track.pt());
            registry.get<TH2>(HIST("Diff/dEdxTPCDG"))->Fill(track.p(), track.tpcSignal());
            registry.get<TH2>(HIST("Diff/IVMptTrkDG"))->Fill(ivm.M(), track.pt());
          } else {
            registry.get<TH2>(HIST("nonDiff/etaptDG"))->Fill(track.eta(), track.pt());
            registry.get<TH2>(HIST("nonDiff/dEdxTPCDG"))->Fill(track.p(), track.tpcSignal());
            registry.get<TH2>(HIST("nonDiff/IVMptTrkDG"))->Fill(ivm.M(), track.pt());
          }
          if (track.hasTOF()) {
            if (isPythiaDiff) {
              registry.get<TH2>(HIST("MBRDiff/dEdxTOFDG"))->Fill(track.p(), track.tofSignal());
            } else if (isGraniittiDiff) {
              registry.get<TH2>(HIST("Diff/dEdxTOFDG"))->Fill(track.p(), track.tofSignal());
            } else {
              registry.get<TH2>(HIST("nonDiff/dEdxTOFDG"))->Fill(track.p(), track.tofSignal());
            }
          }
        }
      }
    }
  }
  PROCESS_SWITCH(DiffMCQA, processMain, "Process Main", true);

  void processFV0(aod::FV0A const& fv0)
  {
    // side A
    for (size_t ind = 0; ind < fv0.channel().size(); ind++) {
      registry.get<TH2>(HIST("FV0/FV0Aamp"))->Fill((fv0.channel())[ind], (fv0.amplitude())[ind]);
    }
  }
  PROCESS_SWITCH(DiffMCQA, processFV0, "Process FV0", true);

  void processFT0(aod::FT0 const& ft0)
  {
    // side A
    for (size_t ind = 0; ind < ft0.channelA().size(); ind++) {
      registry.get<TH2>(HIST("FT0/FT0Aamp"))->Fill((ft0.channelA())[ind], (ft0.amplitudeA())[ind]);
    }

    // side C
    for (size_t ind = 0; ind < ft0.channelC().size(); ind++) {
      registry.get<TH2>(HIST("FT0/FT0Camp"))->Fill((ft0.channelC())[ind], (ft0.amplitudeC())[ind]);
    }
  }
  PROCESS_SWITCH(DiffMCQA, processFT0, "Process FT0", true);

  void processFDD(aod::FDD const& fdd)
  {
    // side A
    for (auto ind = 0; ind < 8; ind++) {
      registry.get<TH2>(HIST("FDD/FDDAamp"))->Fill(ind, (fdd.chargeA())[ind]);
    }

    // side C
    for (auto ind = 0; ind < 8; ind++) {
      registry.get<TH2>(HIST("FDD/FDDCamp"))->Fill(ind, (fdd.chargeC())[ind]);
    }
  }
  PROCESS_SWITCH(DiffMCQA, processFDD, "Process FDD", true);

  // energies:
  //  0: energyZEM1
  //  1: energyZEM2
  //  2: energyCommonZNA
  //  3: energyCommonZNC
  //  4: energyCommonZPA
  //  5: energyCommonZPC
  //  6: energySectorZNA[0]
  //  7: energySectorZNA[1]
  //  8: energySectorZNA[2]
  //  9: energySectorZNA[3]
  // 10: energySectorZNC[0]
  // 11: energySectorZNC[1]
  // 12: energySectorZNC[2]
  // 13: energySectorZNC[3]
  // 14: energySectorZPA[0]
  // 15: energySectorZPA[1]
  // 16: energySectorZPA[2]
  // 17: energySectorZPA[3]
  // 18: energySectorZPC[0]
  // 19: energySectorZPC[1]
  // 20: energySectorZPC[2]
  // 21: energySectorZPC[3]
  void processZDC(aod::Zdc const& zdc)
  {
    // Zdc energies
    registry.get<TH2>(HIST("ZdcEnergies"))->Fill(0., zdc.energyZEM1());
    registry.get<TH2>(HIST("ZdcEnergies"))->Fill(1., zdc.energyZEM2());
    registry.get<TH2>(HIST("ZdcEnergies"))->Fill(2., zdc.energyCommonZNA());
    registry.get<TH2>(HIST("ZdcEnergies"))->Fill(3., zdc.energyCommonZNC());
    registry.get<TH2>(HIST("ZdcEnergies"))->Fill(4., zdc.energyCommonZPA());
    registry.get<TH2>(HIST("ZdcEnergies"))->Fill(5., zdc.energyCommonZPC());
    registry.get<TH2>(HIST("ZdcEnergies"))->Fill(6., (zdc.energySectorZNA())[0]);
    registry.get<TH2>(HIST("ZdcEnergies"))->Fill(7., (zdc.energySectorZNA())[1]);
    registry.get<TH2>(HIST("ZdcEnergies"))->Fill(8., (zdc.energySectorZNA())[2]);
    registry.get<TH2>(HIST("ZdcEnergies"))->Fill(9., (zdc.energySectorZNA())[3]);
    registry.get<TH2>(HIST("ZdcEnergies"))->Fill(10., (zdc.energySectorZNC())[0]);
    registry.get<TH2>(HIST("ZdcEnergies"))->Fill(11., (zdc.energySectorZNC())[1]);
    registry.get<TH2>(HIST("ZdcEnergies"))->Fill(12., (zdc.energySectorZNC())[2]);
    registry.get<TH2>(HIST("ZdcEnergies"))->Fill(13., (zdc.energySectorZNC())[2]);
    registry.get<TH2>(HIST("ZdcEnergies"))->Fill(14., (zdc.energySectorZPA())[0]);
    registry.get<TH2>(HIST("ZdcEnergies"))->Fill(15., (zdc.energySectorZPA())[1]);
    registry.get<TH2>(HIST("ZdcEnergies"))->Fill(16., (zdc.energySectorZPA())[2]);
    registry.get<TH2>(HIST("ZdcEnergies"))->Fill(17., (zdc.energySectorZPA())[3]);
    registry.get<TH2>(HIST("ZdcEnergies"))->Fill(18., (zdc.energySectorZPC())[0]);
    registry.get<TH2>(HIST("ZdcEnergies"))->Fill(19., (zdc.energySectorZPC())[1]);
    registry.get<TH2>(HIST("ZdcEnergies"))->Fill(20., (zdc.energySectorZPC())[2]);
    registry.get<TH2>(HIST("ZdcEnergies"))->Fill(21., (zdc.energySectorZPC())[3]);
  }
  PROCESS_SWITCH(DiffMCQA, processZDC, "Process ZDC", true);

  void processCalo(aod::Calo const& calo)
  {
    // cell number
    registry.get<TH1>(HIST("CaloCell"))->Fill(calo.cellNumber());

    // amplitude
    registry.get<TH1>(HIST("CaloAmplitude"))->Fill(calo.amplitude());
  }
  PROCESS_SWITCH(DiffMCQA, processCalo, "Process Calorimeter", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<DiffMCQA>(cfgc, TaskName{"diffmcqa"})};
}
