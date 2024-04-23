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
/// \brief A task for Asynchronus Quality Control for Ultra-perimpheral and Diffraction (AQC-UD)
/// \author Anisa Khatun, anisa.khatun@cern.ch
/// \author Paul Buehler, paul.buehler@oeaw.ac.at
/// \since  04.05.2023

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "ReconstructionDataFormats/BCRange.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Common/DataModel/FT0Corrected.h"
#include "PWGUD/Core/UDHelpers.h"
#include "Framework/StaticFor.h"
#include "TLorentzVector.h"
#include "TMath.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct UDQC {

  SliceCache cache;
  Preslice<aod::Zdcs> perBCzdc = aod::zdc::bcId;
  Preslice<aod::Calos> perBCcalo = aod::calo::bcId;

  static constexpr std::string_view hcFIT1s[5] = {"CleanFIT1/FV0A", "CleanFIT1/FT0A", "CleanFIT1/FT0C", "CleanFIT1/FDDA", "CleanFIT1/FDDC"};
  static constexpr std::string_view hcFIT2s[5] = {"CleanFIT2/FV0A", "CleanFIT2/FT0A", "CleanFIT2/FT0C", "CleanFIT2/FDDA", "CleanFIT2/FDDC"};
  static constexpr std::string_view hcRelBCs[5] = {"CleanFIT2/BCFV0A", "CleanFIT2/BCFT0A", "CleanFIT2/BCFT0C", "CleanFIT2/BCFDDA", "CleanFIT2/BCFDDC"};

  // global variables
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

  // inivinitialize HistogramRegistry
  HistogramRegistry registry{
    "registry",
    {}};

  // define abbreviations
  using CCs = soa::Join<aod::Collisions, aod::EvSels>;
  using CC = CCs::iterator;
  using BCs = soa::Join<aod::BCs, aod::BcSels, aod::Run3MatchedToBCSparse>;
  using TCs = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA, aod::TrackSelectionExtension, aod::pidTPCFullPi, aod::TOFSignal, aod::pidTOFbeta>;
  using FWs = aod::FwdTracks;
  using ATs = aod::AmbiguousTracks;
  using AFTs = aod::AmbiguousFwdTracks;

  void init(InitContext& context)
  {
    // initialize global variables
    maxdEdxTPC = 0.;
    maxdEdxTOF = 0.;
    diffCuts = (DGCutparHolder)DGCuts;

    // add histograms for the different process functions
    if (context.mOptions.get<bool>("processMain")) {

      // collisions
      registry.add("collisions/Stat", "Cut statistics; Selection criterion; Collisions", {HistType::kTH1F, {{20, -0.5, 19.5}}});
      registry.add("collisions/Tracks", "Number of tracks; Number of tracks; Collisions", {HistType::kTH1F, {{300, 0.5, 300.5}}});
      registry.add("collisions/vtxTracks", "Number of vertex tracks; Number of contributors; Collisions", {HistType::kTH1F, {{300, 0.5, 300.5}}});
      registry.add("collisions/globalTracks", "Number of global tracks; Number of global tracks; Collisions", {HistType::kTH1F, {{300, 0.5, 300.5}}});
      registry.add("collisions/posxy", "Vertex position in x and y direction; V_x; V_y; Collisions", {HistType::kTH2F, {{100, -0.5, 0.5}, {100, -0.5, 0.5}}});
      registry.add("collisions/posz", "Vertex position in z direction; V_z; Collisions", {HistType::kTH1F, {{1000, -100., 100.}}});
      registry.add("collisions/netCharge", "All: net charge; Net charge; collisions with PV tracks", {HistType::kTH1F, {{21, -10.5, 10.5}}});
      registry.add("collisions/notPVTracks", "Not PV tracks; Track status bit; Not PV tracks", {HistType::kTH1F, {{17, -0.5, 16.5}}});
      registry.add("collisions/tResvsrTOFTracks", "Number of PV tracks with TOF hit versus collision time resolution; Collision time resolution [ns]; Fraction of PV tracks with TOF hit; Collisions", {HistType::kTH2F, {{1000, 0., 1.E3}, {101, -0.01, 1.01}}});
      registry.add("collisions/tResvsTOFTrkNoPV", "Number of No PV tracks with TOF hit versus collision time resolution; Collision time resolution [ns]; Fraction of No PV tracks with TOF hit; Collisions", {HistType::kTH2F, {{1000, 0., 1.E3}, {101, -0.01, 1.01}}});

      // tracks
      registry.add("tracks/Stat", "Track bits as function of pT; Track pT; Track bit; Tracks", {HistType::kTH2F, {{100, 0., 5.}, {8, 0.5, 8.5}}});
      registry.add("tracks/etapt", "eta versus pT of all tracks; eta of track; p_T of track [MeV/c^2]; Tracks", {HistType::kTH2F, {{80, -2., 2.}, {100, 0., 5.}}});
      registry.add("tracks/etapt2", "eta versus pT of all quality tracks; eta of track; p_T of track [MeV/c^2]; Tracks", {HistType::kTH2F, {{80, -2., 2.}, {100, 0., 5.}}});
      registry.add("tracks/etapt3", "eta versus pT of all global tracks; eta of track; p_T of track [MeV/c^2]; Tracks", {HistType::kTH2F, {{80, -2., 2.}, {100, 0., 5.}}});
      registry.add("tracks/etapt4", "eta versus pT of tracks with ITS Only; eta of track; p_T of track [MeV/c^2]; Tracks", {HistType::kTH2F, {{80, -2., 2.}, {100, 0., 5.}}});
      registry.add("tracks/etapt5", "eta versus pT of all tracks without TOF; eta of track; p_T of track [MeV/c^2]; Tracks", {HistType::kTH2F, {{80, -2., 2.}, {100, 0., 5.}}});
      registry.add("tracks/etapt6", "eta versus pT of tracks with TPC Only; eta of track; p_T of track [MeV/c^2]; Tracks", {HistType::kTH2F, {{80, -2., 2.}, {100, 0., 5.}}});
      registry.add("tracks/etapt7", "eta versus pT of tracks with TRD Only; eta of track; p_T of track [MeV/c^2]; Tracks", {HistType::kTH2F, {{80, -2., 2.}, {100, 0., 5.}}});
      registry.add("tracks/etapt8", "eta versus pT of tracks with TOF Only; eta of track; p_T of track [MeV/c^2]; Tracks", {HistType::kTH2F, {{80, -2., 2.}, {100, 0., 5.}}});

      registry.add("tracks/dEdxTPC", "TPC signal versus signed track momentum; Signed track momentum [GeV/c]; TPC signal; Tracks", {HistType::kTH2F, {{120, -6., 6.}, {1000, 0., 1000.}}});
      registry.add("tracks/dEdxTOF", "TOF signal versus signed track momentum; Signed track momentum [GeV/c]; TOF signal; Tracks", {HistType::kTH2F, {{1000, 0., 10.}, {1000, 0., 10.}}});

      // DG
      registry.add("DG/PVposxy", "DG: Vertex position in x and y direction; V_x [mm]; V_y [mm]; DG collisions", {HistType::kTH2F, {{100, -0.5, 0.5}, {100, -0.5, 0.5}}});
      registry.add("DG/PVposz", "DG: Vertex position in z direction; V_z; DG collisions", {HistType::kTH1F, {{1000, -100., 100.}}});
      registry.add("DG/netCharge", "DG: net charge; Net charge; DG collisions", {HistType::kTH1F, {{21, -10.5, 10.5}}});
      registry.add("DG/TrackStat", "Track bits as function of pT; Track pT; Track bit; Tracks", {HistType::kTH2F, {{100, 0., 5.}, {10, 0.5, 10.5}}});

      registry.add("DG/etapt", "DG: eta versus pT of all tracks; eta of track; p_T of track [GeV/c^2]; Tracks", {HistType::kTH2F, {{80, -2., 2.}, {100, 0., 5.}}});
      registry.add("DG/etapt2", "DG: eta versus pT of all quality tracks; eta of track; p_T of track [MeV/c^2]; Tracks", {HistType::kTH2F, {{80, -2., 2.}, {100, 0., 5.}}});
      registry.add("DG/etapt3", "DG: eta versus pT of all global tracks; eta of track; p_T of track [MeV/c^2]; Tracks", {HistType::kTH2F, {{80, -2., 2.}, {100, 0., 5.}}});
      registry.add("DG/etapt4", "DG: eta versus pT of all TPC clusers; eta of track; p_T of track [MeV/c^2]; Tracks", {HistType::kTH2F, {{80, -2., 2.}, {100, 0., 5.}}});
      registry.add("DG/etapt5", "DG: eta versus pT of frac.TPCSharedClusters; eta of track; p_T of track [MeV/c^2]; Tracks", {HistType::kTH2F, {{80, -2., 2.}, {100, 0., 5.}}});
      registry.add("DG/etapt6", "DG: eta versus pT of all tracks with ITS; eta of track; p_T of track [MeV/c^2]; Tracks", {HistType::kTH2F, {{80, -2., 2.}, {100, 0., 5.}}});
      registry.add("DG/etapt7", "DG: eta versus pT of all tracks with TPC; eta of track; p_T of track [MeV/c^2]; Tracks", {HistType::kTH2F, {{80, -2., 2.}, {100, 0., 5.}}});
      registry.add("DG/etapt8", "DG: eta versus pT of all tracks with TRD; eta of track; p_T of track [MeV/c^2]; Tracks", {HistType::kTH2F, {{80, -2., 2.}, {100, 0., 5.}}});
      registry.add("DG/etapt9", "DG: eta versus pT of all tracks with TOF; eta of track; p_T of track [MeV/c^2]; Tracks", {HistType::kTH2F, {{80, -2., 2.}, {100, 0., 5.}}});

      registry.add("DG/dEdxTPC", "DG: TPC signal versus signed track momentum; Signed track momentum [GeV/c]; TPC signal; Tracks", {HistType::kTH2F, {{120, -6., 6.}, {1000, 0., 1000.}}});
      registry.add("DG/dEdxTOF", "DG: TOF signal versus signed track momentum; Signed track momentum [GeV/c]; TOF signal; Tracks", {HistType::kTH2F, {{1000, 0., 10.}, {1000, 0., 10.}}});
      registry.add("DG/IVMptSys", "DG: Invariant mass versus p_{T, system}; Invarian mass [GeV/c^2]; p_{T, system} [GeV/c]; DG collisions", {HistType::kTH2F, {{100, 0., 5.}, {350, 0., 3.5}}});
      registry.add("DG/IVMptSys2PVtrk", "DG: Invariant mass versus p_{T, system}; Invarian mass [GeV/c^2]; p_{T, system} [GeV/c]; DG collisions 2 PV tracks", {HistType::kTH2F, {{1000, 0., 10.}, {350, 0., 3.5}}});
      registry.add("DG/etaplus", "DG: eta of positive tracks; eta of track", {HistType::kTH2F, {{1000, -5., 5.}, {120, -6., 6., "#phi"}}});
      registry.add("DG/etaminus", "DG: eta of negative tracks; eta of track", {HistType::kTH2F, {{1000, -5., 5.}, {120, -6., 6., "#phi"}}});
      registry.add("DG/hMass", "DG: Invariant mass of pions; Invarian mass [GeV/c^2]", {HistType::kTH1F, {{1000, 0., 10.}}});
      registry.add("DG/trkDCAxy", "DG: Track DCA of XY; DCAxy", {HistType::kTH1F, {{100, -1., 1.}}});
      registry.add("DG/trkDCAz", "DG: Track DCA of XY;  DCAz", {HistType::kTH1F, {{100, -1., 1.}}});
    }
    if (context.mOptions.get<bool>("processFewProng")) {
      registry.add("fpStat", "#fpStat", {HistType::kTH1F, {{2, 0.5, 2.5}}});
      registry.add("allPVC", "#allPVC", {HistType::kTH1F, {{200, 0.5, 200.5}}});
      registry.add("fpPVC", "#fpPVC", {HistType::kTH1F, {{200, 0.5, 200.5}}});
    }
    if (context.mOptions.get<bool>("processCleanFIT1")) {
      registry.add("CleanFIT1/cleanFIT1", "#cleanFIT1", {HistType::kTH2F, {{20, -0.5, 19.5}, {1000, -0.5, 999.5}}});
      registry.add("CleanFIT1/cF1FV0Aamp", "#cF1FV0Aamp", {HistType::kTH2F, {{20, -0.5, 19.5}, {1000, -0.5, 999.5}}});
      registry.add("CleanFIT1/cF1FT0Aamp", "#cF1FT0Aamp", {HistType::kTH2F, {{20, -0.5, 19.5}, {1000, -0.5, 999.5}}});
      registry.add("CleanFIT1/cF1FT0Camp", "#cF1FT0Camp", {HistType::kTH2F, {{20, -0.5, 19.5}, {1000, -0.5, 999.5}}});
      registry.add("CleanFIT1/cF1FDDAamp", "#cF1FDDAamp", {HistType::kTH2F, {{20, -0.5, 19.5}, {1000, -0.5, 999.5}}});
      registry.add("CleanFIT1/cF1FDDCamp", "#cF1FDDCamp", {HistType::kTH2F, {{20, -0.5, 19.5}, {1000, -0.5, 999.5}}});

      for (auto n{0}; n < 5; n++) {
        registry.add(hcFIT1s[n].data(), hcFIT1s[n].data(), {HistType::kTH2F, {{20, -0.5, 19.5}, {2, -0.5, 1.5}}});
      }
    }
    if (context.mOptions.get<bool>("processCleanFIT2")) {
      registry.add("CleanFIT2/cleanFIT2", "#cleanFIT2", {HistType::kTH2F, {{20, -0.5, 19.5}, {1000, -0.5, 999.5}}});
      registry.add("CleanFIT2/cF2FV0Aamp", "#cF2FV0Aamp", {HistType::kTH2F, {{20, -0.5, 19.5}, {1000, -0.5, 999.5}}});
      registry.add("CleanFIT2/cF2FT0Aamp", "#cF2FT0Aamp", {HistType::kTH2F, {{20, -0.5, 19.5}, {1000, -0.5, 999.5}}});
      registry.add("CleanFIT2/cF2FT0Camp", "#cF2FT0Camp", {HistType::kTH2F, {{20, -0.5, 19.5}, {1000, -0.5, 999.5}}});
      registry.add("CleanFIT2/cF2FDDAamp", "#cF2FDDAamp", {HistType::kTH2F, {{20, -0.5, 19.5}, {1000, -0.5, 999.5}}});
      registry.add("CleanFIT2/cF2FDDCamp", "#cF2FDDCamp", {HistType::kTH2F, {{20, -0.5, 19.5}, {1000, -0.5, 999.5}}});

      for (auto n{0}; n < 5; n++) {
        registry.add(hcFIT2s[n].data(), hcFIT2s[n].data(), {HistType::kTH2F, {{20, -0.5, 19.5}, {2, -0.5, 1.5}}});
        registry.add(hcRelBCs[n].data(), hcRelBCs[n].data(), {HistType::kTH1F, {{3564, -0.5, 3563.5}}});
      }
    }
    if (context.mOptions.get<bool>("processFV0")) {
      registry.add("FV0/FV0A", "#FV0A", {HistType::kTH2F, {{48, -0.5, 47.5}, {2000, 0., 2000.}}});
      registry.add("FV0/hV0A", "Time FV0A", {HistType::kTH1F, {{500, -5.0, 5.0}}});
    }
    if (context.mOptions.get<bool>("processFT0")) {
      registry.add("FT0/FT0A", "#FT0A", {HistType::kTH2F, {{96, -0.5, 95.5}, {400, 0., 400.}}});
      registry.add("FT0/FT0C", "#FT0C", {HistType::kTH2F, {{112, -0.5, 111.5}, {400, 0., 400.}}});
      registry.add("FT0/hT0A", "Time FT0 A side", {HistType::kTH1F, {{500, -5.0, 5.0}}});
      registry.add("FT0/hT0C", "Time FT0 C side", {HistType::kTH1F, {{500, -5.0, 5.0}}});
      registry.add("FT0/hT0ACorr", "Corrected Time FT0 A side", {HistType::kTH1F, {{500, -5.0, 5.0}}});
      registry.add("FT0/hT0CCorr", "Corrected Time FT0 C side", {HistType::kTH1F, {{500, -5.0, 5.0}}});
      registry.add("FT0/hT0AC", "Average Time FT0", {HistType::kTH1F, {{500, -5.0, 5.0}}});
    }
    if (context.mOptions.get<bool>("processFDD")) {
      registry.add("FDD/FDDA", "#FDDA", {HistType::kTH2F, {{8, -0.5, 7.5}, {100, 0., 100.}}});
      registry.add("FDD/FDDC", "#FDDC", {HistType::kTH2F, {{8, -0.5, 7.5}, {100, 0., 100.}}});
      registry.add("FDD/hFDDA", " Avg Time FDD A side", {HistType::kTH1F, {{500, -5.0, 5.0}}});
      registry.add("FDD/hFDDC", " Avg Time FDD C side", {HistType::kTH1F, {{500, -5.0, 5.0}}});
    }
    if (context.mOptions.get<bool>("processZDC")) {
      registry.add("ZdcEnergies", "#ZdcEnergies", {HistType::kTH2F, {{22, -0.5, 21.5}, {100, 0., 1000.}}});
    }
  }

  // ...............................................................................................................
  void processMain(CC const& collision, BCs const& bct0s,
                   TCs const& tracks, FWs const& fwdtracks, ATs const& /*ambtracks*/, AFTs const& /*ambfwdtracks*/,
                   aod::FT0s const& /*ft0s*/, aod::FV0As const& /*fv0as*/, aod::FDDs const& /*fdds*/,
                   aod::Zdcs& zdcs, aod::Calos& calos,
                   aod::V0s const& v0s, aod::Cascades const& cascades)
  {
    LOGF(debug, "<UDQC. Collision %d", collision.globalIndex());
    LOGF(debug, "<UDQC> Start %i", abcrs.size());

    bool isDGcandidate = true;
    registry.get<TH1>(HIST("collisions/Stat"))->Fill(0., isDGcandidate * 1.);

    auto netChargeAll = 0;
    for (auto const& trk : tracks) {
      if (trk.isPVContributor()) {
        netChargeAll += trk.sign();
      }
    }

    // update collision histograms
    // vertex position
    registry.get<TH2>(HIST("collisions/posxy"))->Fill(collision.posX(), collision.posY());
    registry.get<TH1>(HIST("collisions/posz"))->Fill(collision.posZ());
    registry.get<TH1>(HIST("collisions/netCharge"))->Fill(netChargeAll);
    // tracks
    registry.get<TH1>(HIST("collisions/Tracks"))->Fill(tracks.size());
    // vertex tracks normally gives PV contributors from collisions
    registry.get<TH1>(HIST("collisions/vtxTracks"))->Fill(collision.numContrib());
    // global tracks
    Partition<TCs> goodTracks = requireGlobalTrackInFilter();
    goodTracks.bindTable(tracks);
    registry.get<TH1>(HIST("collisions/globalTracks"))->Fill(goodTracks.size());

    // loop over all tracks
    float rgtrwTOF = 0.;
    float norgtrwTOF = 0.;
    for (auto const& track : tracks) {
      // update PV track stats
      if (track.isPVContributor()) {
        registry.get<TH2>(HIST("tracks/Stat"))->Fill(track.pt(), 1., 1.);
        registry.get<TH2>(HIST("tracks/Stat"))->Fill(track.pt(), 2., track.isQualityTrack() * 1.);
        registry.get<TH2>(HIST("tracks/Stat"))->Fill(track.pt(), 3., track.isGlobalTrack() * 1.);
        registry.get<TH2>(HIST("tracks/Stat"))->Fill(track.pt(), 4., (track.hasITS() && !track.hasTPC() && !track.hasTRD() && !track.hasTOF()) * 1.); // tracks with only ITS
        registry.get<TH2>(HIST("tracks/Stat"))->Fill(track.pt(), 5., (track.hasITS() && track.hasTPC() && track.hasTRD() && !track.hasTOF()) * 1.);   // tracks without TOF
        registry.get<TH2>(HIST("tracks/Stat"))->Fill(track.pt(), 6., (!track.hasITS() && track.hasTPC() && !track.hasTRD() && !track.hasTOF()) * 1.); // tracks with only TPC
        registry.get<TH2>(HIST("tracks/Stat"))->Fill(track.pt(), 7., (!track.hasITS() && !track.hasTPC() && track.hasTRD() && !track.hasTOF()) * 1.); // tracks with only TRD
        registry.get<TH2>(HIST("tracks/Stat"))->Fill(track.pt(), 8., (!track.hasITS() && !track.hasTPC() && !track.hasTRD() && track.hasTOF()) * 1.); // tracks with only TOF

        // update eta vs pt histograms
        registry.get<TH2>(HIST("tracks/etapt"))->Fill(track.eta(), track.pt(), 1.);
        registry.get<TH2>(HIST("tracks/etapt2"))->Fill(track.eta(), track.pt(), track.isQualityTrack() * 1.);
        registry.get<TH2>(HIST("tracks/etapt3"))->Fill(track.eta(), track.pt(), track.isGlobalTrack() * 1.);
        registry.get<TH2>(HIST("tracks/etapt4"))->Fill(track.eta(), track.pt(), (track.hasITS() && !track.hasTPC() && !track.hasTRD() && !track.hasTOF()) * 1.); // tracks with only ITS
        registry.get<TH2>(HIST("tracks/etapt5"))->Fill(track.eta(), track.pt(), (track.hasITS() && track.hasTPC() && track.hasTRD() && !track.hasTOF()) * 1.);   // tracks without TOF
        registry.get<TH2>(HIST("tracks/etapt6"))->Fill(track.eta(), track.pt(), (!track.hasITS() && track.hasTPC() && !track.hasTRD() && !track.hasTOF()) * 1.); // tracks with only TPC
        registry.get<TH2>(HIST("tracks/etapt7"))->Fill(track.eta(), track.pt(), (!track.hasITS() && !track.hasTPC() && track.hasTRD() && !track.hasTOF()) * 1.); // tracks with only TRD
        registry.get<TH2>(HIST("tracks/etapt8"))->Fill(track.eta(), track.pt(), (!track.hasITS() && !track.hasTPC() && !track.hasTRD() && track.hasTOF()) * 1.); // tracks with only TOF

        // update dEdx histograms
        registry.get<TH2>(HIST("tracks/dEdxTPC"))->Fill(track.tpcInnerParam() / track.sign(), track.tpcSignal());

        if (track.tpcSignal() > maxdEdxTPC) {
          maxdEdxTPC = track.tpcSignal();
          LOGF(debug, "<UDQC> New maxdEdx TPC %f", maxdEdxTPC);
        }

        // TOF hit?
        if (track.hasTOF()) {
          registry.get<TH2>(HIST("tracks/dEdxTOF"))->Fill(track.p() / track.sign(), track.beta());
          if (track.tofSignal() > maxdEdxTOF) {
            maxdEdxTOF = track.tofSignal();
            LOGF(debug, "<UDQC> New maxdEdx TOF %f", maxdEdxTOF);
          }

          // No vertex track with TOF hit?
          if (!track.isPVContributor()) {
            norgtrwTOF += 1.;
          }

          // vertex track with TOF hit?
          if (track.isPVContributor()) {
            rgtrwTOF += 1.;
          }
        }
      }
    } // closing track loop

    // fraction of No PV tracks with TOF hit
    if (collision.numContrib() > 0) {
      norgtrwTOF /= collision.numContrib();
    }

    // fraction of PV tracks with TOF hit
    if (collision.numContrib() > 0) {
      rgtrwTOF /= collision.numContrib();
    }
    LOGF(debug, "<UDQC> PV tracks with TOF: %f [1]", rgtrwTOF);
    registry.get<TH2>(HIST("collisions/tResvsrTOFTracks"))->Fill(collision.collisionTimeRes(), rgtrwTOF);
    registry.get<TH2>(HIST("collisions/tResvsTOFTrkNoPV"))->Fill(collision.collisionTimeRes(), norgtrwTOF);

    // is it a DG candidate?
    // 1. DG = no FIT signal in compatible BCs
    // 2. & no ZDC signal in compatible BCs
    // 3. & no Calo signal in compatible BCs
    // 4. & no V0s
    // 5. & no Cascades
    // 6. & number of forward tracks = 0
    // 7. & no global track which is not a vertex track
    // 8. & no vertex track which is not a global track
    // 9. & ntrMin <= number of vertex tracks <= ntrMax
    isDGcandidate = true;

    // get BCrange to test for FIT signals
    auto bcSlice = udhelpers::compatibleBCs(collision, diffCuts.NDtcoll(), bct0s, diffCuts.minNBCs());

    // 1. no FIT signal in bcSlice / collision
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
    registry.get<TH1>(HIST("collisions/Stat"))->Fill(1., isDGcandidate * 1.);

    // 2. no Zdc signal in bcSlice
    std::vector<float> lims(10, 0.);
    for (auto const& bc : bcSlice) {
      if (!udhelpers::cleanZDC(bc, zdcs, lims, cache)) {
        isDGcandidate = false;
        break;
      }
    }
    registry.get<TH1>(HIST("collisions/Stat"))->Fill(2., isDGcandidate * 1.);

    // 3. no Calo signal in bcSlice
    for (auto const& bc : bcSlice) {
      if (!udhelpers::cleanCalo(bc, calos, lims, cache)) {
        isDGcandidate = false;
        break;
      }
    }
    registry.get<TH1>(HIST("collisions/Stat"))->Fill(3., isDGcandidate * 1.);

    // 4. no V0s
    isDGcandidate &= (v0s.size() == 0);
    registry.get<TH1>(HIST("collisions/Stat"))->Fill(4., isDGcandidate * 1.);

    // 5. no Cascades
    isDGcandidate &= (cascades.size() == 0);
    registry.get<TH1>(HIST("collisions/Stat"))->Fill(5., isDGcandidate * 1.);

    // 6. number of forward tracks = 0
    isDGcandidate &= (fwdtracks.size() == 0);
    registry.get<TH1>(HIST("collisions/Stat"))->Fill(6., isDGcandidate * 1.);

    // 7. Check for global tracks which are no vtx tracks
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
    registry.get<TH1>(HIST("collisions/Stat"))->Fill(7., globalAndVtx * 1.);

    // 8. check a given bc for possible ambiguous Tracks
    auto noAmbTracks = isDGcandidate;
    for (auto const& bc : bcSlice) {
      if (abcrs.isInRange(bc.globalIndex())) {
        noAmbTracks = false;
        break;
      }
    }
    registry.get<TH1>(HIST("collisions/Stat"))->Fill(8., noAmbTracks * 1.); // noAmbTracks

    // 9. check a given bc for possible ambiguous FwdTracks
    auto noAmbFwdTracks = isDGcandidate;
    for (auto const& bc : bcSlice) {
      if (afbcrs.isInRange(bc.globalIndex())) {
        noAmbFwdTracks = false;
        break;
      }
    }
    registry.get<TH1>(HIST("collisions/Stat"))->Fill(9., noAmbFwdTracks * 1.); // noAmbFwdTracks

    // 10. number of vertex tracks <= n
    isDGcandidate &= (collision.numContrib() >= diffCuts.minNTracks());
    registry.get<TH1>(HIST("collisions/Stat"))->Fill(10., isDGcandidate * 1.);
    isDGcandidate &= (collision.numContrib() <= diffCuts.maxNTracks());
    registry.get<TH1>(HIST("collisions/Stat"))->Fill(11., isDGcandidate * 1.);

    // 11. fraction of PV tracks with TOF hit
    isDGcandidate &= (rgtrwTOF >= diffCuts.minRgtrwTOF());
    registry.get<TH1>(HIST("collisions/Stat"))->Fill(12., isDGcandidate * 1.);

    // 12. net charge and invariant mass
    bool goodetas = true;
    bool goodpts = true;
    bool goodnchs = true;
    auto netCharge = 0;
    auto lvtmp = TLorentzVector();
    auto ivm = TLorentzVector();
    if (isDGcandidate) {

      // which particle hypothesis? // Pion hypothesis has been used currently
      auto mass2Use = constants::physics::MassPionCharged;
      if (diffCuts.pidHypothesis() == 321) {
        mass2Use = constants::physics::MassKaonCharged;
      }

      // check also pt and eta of tracks
      for (auto const& track : tracks) {
        // update histogram for rejectedTracks/notPVtracks
        if (!track.isPVContributor()) {
          registry.get<TH1>(HIST("collisions/notPVTracks"))->Fill(0., 1.);
          registry.get<TH1>(HIST("collisions/notPVTracks"))->Fill(1., track.isGlobalTrackSDD() * 1.);
          registry.get<TH1>(HIST("collisions/notPVTracks"))->Fill(2., track.passedTrackType() * 1.);
          registry.get<TH1>(HIST("collisions/notPVTracks"))->Fill(3., track.passedPtRange() * 1.);
          registry.get<TH1>(HIST("collisions/notPVTracks"))->Fill(4., track.passedEtaRange() * 1.);
          registry.get<TH1>(HIST("collisions/notPVTracks"))->Fill(5., track.passedTPCNCls() * 1.);
          registry.get<TH1>(HIST("collisions/notPVTracks"))->Fill(6., track.passedTPCCrossedRows() * 1.);
          registry.get<TH1>(HIST("collisions/notPVTracks"))->Fill(7., track.passedTPCCrossedRowsOverNCls() * 1.);
          registry.get<TH1>(HIST("collisions/notPVTracks"))->Fill(8., track.passedTPCChi2NDF() * 1.);
          registry.get<TH1>(HIST("collisions/notPVTracks"))->Fill(9., track.passedTPCRefit() * 1.);
          registry.get<TH1>(HIST("collisions/notPVTracks"))->Fill(10., track.passedITSNCls() * 1.);
          registry.get<TH1>(HIST("collisions/notPVTracks"))->Fill(11., track.passedITSChi2NDF() * 1.);
          registry.get<TH1>(HIST("collisions/notPVTracks"))->Fill(12., track.passedITSRefit() * 1.);
          registry.get<TH1>(HIST("collisions/notPVTracks"))->Fill(13., track.passedITSHits() * 1.);
          registry.get<TH1>(HIST("collisions/notPVTracks"))->Fill(14., track.passedGoldenChi2() * 1.);
          registry.get<TH1>(HIST("collisions/notPVTracks"))->Fill(15., track.passedDCAxy() * 1.);
          registry.get<TH1>(HIST("collisions/notPVTracks"))->Fill(16., track.passedDCAz() * 1.);
          continue;
        }
        // define Lorentz vector to create invariant mass
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
    registry.get<TH1>(HIST("collisions/Stat"))->Fill(13., isDGcandidate * 1.);
    isDGcandidate &= goodetas;
    registry.get<TH1>(HIST("collisions/Stat"))->Fill(14., isDGcandidate * 1.);

    auto netChargeValues = diffCuts.netCharges();

    if (std::find(netChargeValues.begin(), netChargeValues.end(), netCharge) == netChargeValues.end()) {
      goodnchs = false;
    }
    isDGcandidate &= goodnchs;
    registry.get<TH1>(HIST("collisions/Stat"))->Fill(15., isDGcandidate * 1.);
    registry.get<TH1>(HIST("collisions/Stat"))->Fill(16., isDGcandidate * 1.);
    isDGcandidate &= (ivm.M() >= diffCuts.minIVM());
    registry.get<TH1>(HIST("collisions/Stat"))->Fill(17., isDGcandidate * 1.);
    isDGcandidate &= (ivm.M() <= diffCuts.maxIVM());
    registry.get<TH1>(HIST("collisions/Stat"))->Fill(18., isDGcandidate * 1.);

    // 8. check for vertex tracks which are no global tracks
    isDGcandidate &= globalAndVtx;
    if (diffCuts.globalTracksOnly()) {
      isDGcandidate &= vtxAndGlobal;
    }
    registry.get<TH1>(HIST("collisions/Stat"))->Fill(19., isDGcandidate * 1.);

    // update some DG histograms
    if (isDGcandidate) {
      // vertex position of DG events
      registry.get<TH2>(HIST("DG/PVposxy"))->Fill(collision.posX(), collision.posY());
      registry.get<TH1>(HIST("DG/PVposz"))->Fill(collision.posZ());
      registry.get<TH1>(HIST("DG/netCharge"))->Fill(netCharge);

      // Invariant mass with 2 PV contributors and all contributors
      if (collision.numContrib() == 2) {
        bool ispipiCand = true;
        for (auto const& track : tracks) {
          if (track.isPVContributor()) {
            if (std::abs(track.tpcNSigmaPi()) > diffCuts.maxNSigmaTPC()) {
              ispipiCand = false;
            }
          }
        }
        if (ispipiCand) {
          registry.get<TH1>(HIST("DG/hMass"))->Fill(ivm.M());
          registry.get<TH2>(HIST("DG/IVMptSys2PVtrk"))->Fill(ivm.M(), ivm.Pt());
        }
      } else {
        registry.get<TH2>(HIST("DG/IVMptSys"))->Fill(ivm.M(), ivm.Pt());
      }
      // fill dEdx and kinematics of DG event tracks
      for (auto const& track : tracks) {
        if (track.isPVContributor()) {
          registry.get<TH2>(HIST("DG/TrackStat"))->Fill(track.pt(), 1., 1.);
          registry.get<TH2>(HIST("DG/TrackStat"))->Fill(track.pt(), 2., track.isQualityTrack() * 1.);
          registry.get<TH2>(HIST("DG/TrackStat"))->Fill(track.pt(), 3., track.isGlobalTrack() * 1.);
          registry.get<TH2>(HIST("DG/TrackStat"))->Fill(track.pt(), 4., track.tpcNClsFound() * 1.);
          registry.get<TH2>(HIST("DG/TrackStat"))->Fill(track.pt(), 5., track.tpcFractionSharedCls() * 1.);
          registry.get<TH2>(HIST("DG/TrackStat"))->Fill(track.pt(), 6., (track.hasITS() && !track.hasTPC() && !track.hasTRD() && !track.hasTOF()) * 1.);
          registry.get<TH2>(HIST("DG/TrackStat"))->Fill(track.pt(), 7., (!track.hasITS() && track.hasTPC() && !track.hasTRD() && !track.hasTOF()) * 1.);
          registry.get<TH2>(HIST("DG/TrackStat"))->Fill(track.pt(), 8., (!track.hasITS() && !track.hasTPC() && track.hasTRD() && !track.hasTOF()) * 1.);
          registry.get<TH2>(HIST("DG/TrackStat"))->Fill(track.pt(), 9., (!track.hasITS() && !track.hasTPC() && !track.hasTRD() && track.hasTOF()) * 1.);
          registry.get<TH2>(HIST("DG/TrackStat"))->Fill(track.pt(), 10, track.dcaXY() * 1.);
          registry.get<TH2>(HIST("DG/TrackStat"))->Fill(track.pt(), 11, track.dcaZ() * 1.);

          registry.get<TH2>(HIST("DG/etapt"))->Fill(track.eta(), track.pt(), 1.);
          registry.get<TH2>(HIST("DG/etapt2"))->Fill(track.eta(), track.pt(), track.isQualityTrack() * 1.);
          registry.get<TH2>(HIST("DG/etapt3"))->Fill(track.eta(), track.pt(), track.isGlobalTrack() * 1.);
          registry.get<TH2>(HIST("DG/etapt4"))->Fill(track.eta(), track.pt(), track.tpcNClsFound() * 1.);
          registry.get<TH2>(HIST("DG/etapt5"))->Fill(track.eta(), track.pt(), track.tpcFractionSharedCls() * 1.);
          registry.get<TH2>(HIST("DG/etapt6"))->Fill(track.eta(), track.pt(), (track.hasITS() && !track.hasTPC() && !track.hasTRD() && !track.hasTOF()) * 1.);
          registry.get<TH2>(HIST("DG/etapt7"))->Fill(track.eta(), track.pt(), (!track.hasITS() && track.hasTPC() && !track.hasTRD() && !track.hasTOF()) * 1.);
          registry.get<TH2>(HIST("DG/etapt8"))->Fill(track.eta(), track.pt(), (!track.hasITS() && !track.hasTPC() && track.hasTRD() && !track.hasTOF()) * 1.);
          registry.get<TH2>(HIST("DG/etapt9"))->Fill(track.eta(), track.pt(), (!track.hasITS() && !track.hasTPC() && !track.hasTRD() && track.hasTOF()) * 1.);
          // eta distribution of positive and negative charged tracks
          if (track.sign() > 0) {
            registry.get<TH2>(HIST("DG/etaplus"))->Fill(track.eta(), track.phi());
          }
          if (track.sign() < 0) {
            registry.get<TH2>(HIST("DG/etaminus"))->Fill(track.eta(), track.phi());
          }

          LOGF(debug, "dEdx TPC %f TOF %i %f", track.tpcSignal(), track.hasTOF(), track.hasTOF() ? track.tofSignal() : 0.);
          if (collision.numContrib() == 2) {
            registry.get<TH2>(HIST("DG/dEdxTPC"))->Fill(track.tpcInnerParam() / track.sign(), track.tpcSignal());
            registry.get<TH1>(HIST("DG/trkDCAxy"))->Fill(track.dcaXY());
            registry.get<TH1>(HIST("DG/trkDCAz"))->Fill(track.dcaZ());
          }

          if (track.hasTOF()) {
            registry.get<TH2>(HIST("DG/dEdxTOF"))->Fill(track.p() / track.sign(), track.beta());
          }
        }
      }
    }
  }
  PROCESS_SWITCH(UDQC, processMain, "Process Main", true);

  // ...............................................................................................................
  // Distribution of number of PV contributors for all collisions and those with empty FT0
  void processFewProng(CC const& collision, BCs const& /*bct0s*/,
                       aod::FT0s const& /*ft0s*/, aod::FV0As const& /*fv0as*/, aod::FDDs const& /*fdds*/)
  {
    // count collisions
    registry.get<TH1>(HIST("fpStat"))->Fill(1., 1.);
    registry.get<TH1>(HIST("allPVC"))->Fill(collision.numContrib(), 1.);

    // check FT0 to be empty
    auto bc = collision.foundBC_as<BCs>();
    if (udhelpers::cleanFT0(bc, diffCuts.maxFITtime(), 0., 0.)) {
      // only collisions with empty FT0 arrive here
      registry.get<TH1>(HIST("fpStat"))->Fill(2., 1.);

      // update #PV contributors in collisions with empty FT0
      registry.get<TH1>(HIST("fpPVC"))->Fill(collision.numContrib(), 1.);
    }
  }
  PROCESS_SWITCH(UDQC, processFewProng, "Process FewProng", true);

  // .............................................................................................................................................
  void processCleanFIT1(CC const& collision, BCs const& bct0s,
                        aod::FT0s const& /*ft0s*/, aod::FV0As const& /*fv0as*/, aod::FDDs const& /*fdds*/)
  {
    LOGF(debug, "<CleanFit. Collision %d", collision.globalIndex());

    // test influence of BCrange width using a series of NDtcoll
    float ampFV0A, ampFT0A, ampFT0C, ampFDDA, ampFDDC;
    auto FITlims = std::vector<float>(5, 1000000.);
    bool isDGcandidate = true;
    for (int NDtcoll = 0; NDtcoll < 20; NDtcoll++) {
      auto bcSlice = udhelpers::compatibleBCs(collision, NDtcoll, bct0s, 0);

      // do for diffCuts.FITAmpLimits
      ampFV0A = ampFT0A = ampFT0C = ampFDDA = ampFDDC = 0.;
      isDGcandidate = true;
      for (auto const& bc : bcSlice) {
        isDGcandidate &= udhelpers::cleanFIT(bc, diffCuts.maxFITtime(), diffCuts.FITAmpLimits());

        if (bc.has_foundFV0()) {
          ampFV0A += udhelpers::FV0AmplitudeA(bc.foundFV0());
        }
        if (bc.has_foundFT0()) {
          ampFT0A += udhelpers::FT0AmplitudeA(bc.foundFT0());
          ampFT0C += udhelpers::FT0AmplitudeA(bc.foundFT0());
        }
        if (bc.has_foundFDD()) {
          ampFDDA += udhelpers::FDDAmplitudeA(bc.foundFDD());
          ampFDDC += udhelpers::FDDAmplitudeA(bc.foundFDD());
        }
      }
      registry.get<TH2>(HIST("CleanFIT1/cleanFIT1"))->Fill(NDtcoll, isDGcandidate * 1.);
      if (isDGcandidate) {
        registry.get<TH2>(HIST("CleanFIT1/cF1FV0Aamp"))->Fill(NDtcoll, ampFV0A);
        registry.get<TH2>(HIST("CleanFIT1/cF1FT0Aamp"))->Fill(NDtcoll, ampFT0A);
        registry.get<TH2>(HIST("CleanFIT1/cF1FT0Camp"))->Fill(NDtcoll, ampFT0C);
        registry.get<TH2>(HIST("CleanFIT1/cF1FDDAamp"))->Fill(NDtcoll, ampFDDA);
        registry.get<TH2>(HIST("CleanFIT1/cF1FDDCamp"))->Fill(NDtcoll, ampFDDC);
      }

      // loop over single detectors
      static_for<0, 4>([&](auto n) {
        FITlims[n] = 0.;
        isDGcandidate = true;
        for (auto const& bc : bcSlice) {
          isDGcandidate &= udhelpers::cleanFIT(bc, diffCuts.maxFITtime(), FITlims);
        }
        constexpr int index = n.value;
        registry.fill(HIST(hcFIT1s[index]), NDtcoll, isDGcandidate * 1.);
        FITlims[n] = 1000000.;
      });
    }
  }

  PROCESS_SWITCH(UDQC, processCleanFIT1, "Process CleanFitTest1", true);
  // .............................................................................................................................................

  void processCleanFIT2(CC const& collision, BCs const& bct0s,
                        aod::FT0s const& /*ft0s*/, aod::FV0As const& /*fv0as*/, aod::FDDs const& /*fdds*/)
  {
    LOGF(debug, "<CleanFit. Collision %d", collision.globalIndex());
    uint64_t bcnum = 0;
    if (collision.has_foundBC()) {
      auto collbc = collision.foundBC_as<BCs>();
      bcnum = collbc.globalBC() % o2::constants::lhc::LHCMaxBunches;
    }
    // test influence of BCrange width using a series of nMinBC
    float ampFV0A, ampFT0A, ampFT0C, ampFDDA, ampFDDC;
    auto FITlims = std::vector<float>(5, 1000000.);
    bool isDGcandidate = true;
    for (int nMinBC = 0; nMinBC < 20; nMinBC++) {
      auto bcSlice = udhelpers::compatibleBCs(collision, 0, bct0s, nMinBC);
      ampFV0A = ampFT0A = ampFT0C = ampFDDA = ampFDDC = 0.;
      isDGcandidate = true;
      for (auto const& bc : bcSlice) {
        isDGcandidate &= udhelpers::cleanFIT(bc, diffCuts.maxFITtime(), diffCuts.FITAmpLimits());

        if (bc.has_foundFV0()) {
          ampFV0A += udhelpers::FV0AmplitudeA(bc.foundFV0());
        }
        if (bc.has_foundFT0()) {
          ampFT0A += udhelpers::FT0AmplitudeA(bc.foundFT0());
          ampFT0C += udhelpers::FT0AmplitudeA(bc.foundFT0());
        }
        if (bc.has_foundFDD()) {
          ampFDDA += udhelpers::FDDAmplitudeA(bc.foundFDD());
          ampFDDC += udhelpers::FDDAmplitudeA(bc.foundFDD());
        }
      }
      registry.get<TH2>(HIST("CleanFIT2/cleanFIT2"))->Fill(nMinBC, isDGcandidate * 1.);

      if (isDGcandidate) {
        registry.get<TH2>(HIST("CleanFIT2/cF2FV0Aamp"))->Fill(nMinBC, ampFV0A);
        registry.get<TH2>(HIST("CleanFIT2/cF2FT0Aamp"))->Fill(nMinBC, ampFT0A);
        registry.get<TH2>(HIST("CleanFIT2/cF2FT0Camp"))->Fill(nMinBC, ampFT0C);
        registry.get<TH2>(HIST("CleanFIT2/cF2FDDAamp"))->Fill(nMinBC, ampFDDA);
        registry.get<TH2>(HIST("CleanFIT2/cF2FDDCamp"))->Fill(nMinBC, ampFDDC);
      }

      // loop over single detectors
      static_for<0, 4>([&](auto n) {
        FITlims[n] = 0.;
        isDGcandidate = true;
        for (auto const& bc : bcSlice) {
          isDGcandidate &= udhelpers::cleanFIT(bc, diffCuts.maxFITtime(), FITlims); // DG
        }
        constexpr int index = n.value;
        registry.fill(HIST(hcFIT2s[index]), nMinBC, isDGcandidate * 1.);
        registry.fill(HIST(hcRelBCs[index]), static_cast<float>(bcnum), isDGcandidate * 1.);
        FITlims[n] = 1000000.;
      });
    }
  }

  PROCESS_SWITCH(UDQC, processCleanFIT2, "Process CleanFitTest2", true);

  // ...............................................................................................................
  void processFV0(aod::FV0As const& fv0s, BCs const&)
  {
    LOGF(info, "<FV0Signals> %d", fv0s.size());
    if (fv0s.size() <= 0) {
      return;
    }

    for (auto fv0 : fv0s) {
      registry.get<TH1>(HIST("FV0/hV0A"))->Fill(fv0.time());
      // side A
      for (size_t ind = 0; ind < fv0.channel().size(); ind++) {
        registry.get<TH2>(HIST("FV0/FV0A"))->Fill((fv0.channel())[ind], (fv0.amplitude())[ind]);
      }
    }
  };
  PROCESS_SWITCH(UDQC, processFV0, "Process FV0", true);

  // ...............................................................................................................
  void processFT0(aod::FT0s const& ft0s, aod::FT0sCorrected const& ft0scorr, BCs const&)
  {
    LOGF(debug, "<processFT0> %d", ft0s.size());
    for (auto const& collision : ft0scorr) {

      if (collision.t0ACorrectedValid()) {
        registry.get<TH1>(HIST("FT0/hT0ACorr"))->Fill(collision.t0ACorrected());
      }
      if (collision.t0CCorrectedValid()) {
        registry.get<TH1>(HIST("FT0/hT0CCorr"))->Fill(collision.t0CCorrected());
      }

      if (collision.t0CCorrectedValid() && collision.t0ACorrectedValid()) {
        registry.get<TH1>(HIST("FT0/hT0AC"))->Fill(collision.t0AC());
      }
    }
    for (auto ft0 : ft0s) {
      registry.get<TH1>(HIST("FT0/hT0A"))->Fill(ft0.timeA());
      registry.get<TH1>(HIST("FT0/hT0C"))->Fill(ft0.timeC());

      // side A
      for (size_t ind = 0; ind < ft0.channelA().size(); ind++) {
        registry.get<TH2>(HIST("FT0/FT0A"))->Fill((ft0.channelA())[ind], (ft0.amplitudeA())[ind]);
      }

      // side C
      for (size_t ind = 0; ind < ft0.channelC().size(); ind++) {
        registry.get<TH2>(HIST("FT0/FT0C"))->Fill((ft0.channelC())[ind], (ft0.amplitudeC())[ind]);
      }
    }
  };
  PROCESS_SWITCH(UDQC, processFT0, "Process FT0", true);

  // ...............................................................................................................
  void processFDD(aod::FDDs const& fdds, BCs const&)
  {
    LOGF(debug, "<FDDSignals> %d", fdds.size());

    for (auto fdd : fdds) {

      registry.get<TH1>(HIST("FDD/hFDDA"))->Fill(fdd.timeA());
      registry.get<TH1>(HIST("FDD/hFDDC"))->Fill(fdd.timeC());
      // side A
      for (auto ind = 0; ind < 8; ind++) {
        registry.get<TH2>(HIST("FDD/FDDA"))->Fill(ind, (fdd.chargeA())[ind]);
      }

      // side C
      for (auto ind = 0; ind < 8; ind++) {
        registry.get<TH2>(HIST("FDD/FDDC"))->Fill(ind, (fdd.chargeC())[ind]);
      }
    }
  };
  PROCESS_SWITCH(UDQC, processFDD, "Process FDD", true);

  // ...............................................................................................................
  void processZDC(aod::Zdc const& zdc)
  {
    LOGF(debug, "<ZDCSignals> %d", zdc.size());

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
    registry.get<TH2>(HIST("ZdcEnergies"))->Fill(13., (zdc.energySectorZNC())[3]);
    registry.get<TH2>(HIST("ZdcEnergies"))->Fill(14., (zdc.energySectorZPA())[0]);
    registry.get<TH2>(HIST("ZdcEnergies"))->Fill(15., (zdc.energySectorZPA())[1]);
    registry.get<TH2>(HIST("ZdcEnergies"))->Fill(16., (zdc.energySectorZPA())[2]);
    registry.get<TH2>(HIST("ZdcEnergies"))->Fill(17., (zdc.energySectorZPA())[3]);
    registry.get<TH2>(HIST("ZdcEnergies"))->Fill(18., (zdc.energySectorZPC())[0]);
    registry.get<TH2>(HIST("ZdcEnergies"))->Fill(19., (zdc.energySectorZPC())[1]);
    registry.get<TH2>(HIST("ZdcEnergies"))->Fill(20., (zdc.energySectorZPC())[2]);
    registry.get<TH2>(HIST("ZdcEnergies"))->Fill(21., (zdc.energySectorZPC())[3]);
  };
  PROCESS_SWITCH(UDQC, processZDC, "Process ZDC", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<UDQC>(cfgc, TaskName{"udQC"}),
  };
}
