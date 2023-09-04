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
/// \since  20.05.2022

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "TLorentzVector.h"
#include "CommonConstants/LHCConstants.h"
#include "ReconstructionDataFormats/BCRange.h"
#include "CommonConstants/PhysicsConstants.h"
#include "PWGUD/Core/UDHelpers.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct DiffQA {
  SliceCache cache;
  Preslice<aod::Zdcs> perBCzdc = aod::zdc::bcId;
  Preslice<aod::Calos> perBCcalo = aod::calo::bcId;

  // constants
  static const int ns = 20;    // number of BCs to save (used in processF[V0, T0, DD])
  static const int ncmin = 20; // minimum length of series of empty BCs  (used in processF[V0, T0, DD])

  static constexpr std::string_view hFV0A[ns + 1] = {"FV0/A00", "FV0/A01", "FV0/A02", "FV0/A03", "FV0/A04", "FV0/A05", "FV0/A06", "FV0/A07", "FV0/A08", "FV0/A09", "FV0/A10", "FV0/A11", "FV0/A12", "FV0/A13", "FV0/A14", "FV0/A15", "FV0/A16", "FV0/A17", "FV0/A18", "FV0/A19", "FV0/A20"};
  static constexpr std::string_view hFT0A[ns + 1] = {"FT0/A00", "FT0/A01", "FT0/A02", "FT0/A03", "FT0/A04", "FT0/A05", "FT0/A06", "FT0/A07", "FT0/A08", "FT0/A09", "FT0/A10", "FT0/A11", "FT0/A12", "FT0/A13", "FT0/A14", "FT0/A15", "FT0/A16", "FT0/A17", "FT0/A18", "FT0/A19", "FT0/A20"};
  static constexpr std::string_view hFT0C[ns + 1] = {"FT0/C00", "FT0/C01", "FT0/C02", "FT0/C03", "FT0/C04", "FT0/C05", "FT0/C06", "FT0/C07", "FT0/C08", "FT0/C09", "FT0/C10", "FT0/C11", "FT0/C12", "FT0/C13", "FT0/C14", "FT0/C15", "FT0/C16", "FT0/C17", "FT0/C18", "FT0/C19", "FT0/C20"};
  static constexpr std::string_view hFDDA[ns + 1] = {"FDD/A00", "FDD/A01", "FDD/A02", "FDD/A03", "FDD/A04", "FDD/A05", "FDD/A06", "FDD/A07", "FDD/A08", "FDD/A09", "FDD/A10", "FDD/A11", "FDD/A12", "FDD/A13", "FDD/A14", "FDD/A15", "FDD/A16", "FDD/A17", "FDD/A18", "FDD/A19", "FDD/A20"};
  static constexpr std::string_view hFDDC[ns + 1] = {"FDD/C00", "FDD/C01", "FDD/C02", "FDD/C03", "FDD/C04", "FDD/C05", "FDD/C06", "FDD/C07", "FDD/C08", "FDD/C09", "FDD/C10", "FDD/C11", "FDD/C12", "FDD/C13", "FDD/C14", "FDD/C15", "FDD/C16", "FDD/C17", "FDD/C18", "FDD/C19", "FDD/C20"};

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

  // initialize HistogramRegistry
  HistogramRegistry registry{
    "registry",
    {}};

  // define abbreviations
  using CCs = soa::Join<aod::Collisions, aod::EvSels>;
  using CC = CCs::iterator;
  using BCs = soa::Join<aod::BCs, aod::BcSels, aod::Run3MatchedToBCSparse>;
  using TCs = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TrackSelectionExtension, aod::TOFSignal>;
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
      registry.add("collisions/PVposxy", "Vertex position in x and y direction; V_x; V_y; Collisions", {HistType::kTH2F, {{100, -0.5, 0.5}, {100, -0.5, 0.5}}});
      registry.add("collisions/PVposz", "Vertex position in z direction; V_z; Collisions", {HistType::kTH1F, {{1000, -100., 100.}}});
      registry.add("collisions/Tracks", "Number of tracks; Number of tracks; Collisions", {HistType::kTH1F, {{50, 0.5, 50.5}}});
      registry.add("collisions/PVTracks", "Number of PV tracks; Number of PV tracks; Collisions", {HistType::kTH1F, {{50, 0.5, 50.5}}});
      registry.add("collisions/globalTracks", "Number of global tracks; Number of global tracks; Collisions", {HistType::kTH1F, {{50, 0.5, 50.5}}});
      registry.add("collisions/notPVTracks", "Not PV tracks; Track status bit; Not PV tracks", {HistType::kTH1F, {{17, -0.5, 16.5}}});
      registry.add("collisions/tResvsrTOFTracks", "Number of PV tracks with TOF hit versus collision time resolution; Collision time resolution [ns]; Fraction of PV tracks with TOF hit; Collisions", {HistType::kTH2F, {{1000, 0., 1.E3}, {101, -0.01, 1.01}}});

      // tracks
      registry.add("tracks/Stat", "Track bits as function of pT; Track pT; Track bit; Tracks", {HistType::kTH2F, {{100, 0., 5.}, {8, 0.5, 8.5}}});

      registry.add("tracks/etapt", "eta versus pT of all tracks; eta of track; p_T of track [MeV/c^2]; Tracks", {HistType::kTH2F, {{80, -2., 2.}, {100, 0., 5.}}});
      registry.add("tracks/etapt2", "eta versus pT of all quality tracks; eta of track; p_T of track [MeV/c^2]; Tracks", {HistType::kTH2F, {{80, -2., 2.}, {100, 0., 5.}}});
      registry.add("tracks/etapt3", "eta versus pT of all global tracks; eta of track; p_T of track [MeV/c^2]; Tracks", {HistType::kTH2F, {{80, -2., 2.}, {100, 0., 5.}}});
      registry.add("tracks/etapt4", "eta versus pT of all PV tracks; eta of track; p_T of track [MeV/c^2]; Tracks", {HistType::kTH2F, {{80, -2., 2.}, {100, 0., 5.}}});
      registry.add("tracks/etapt5", "eta versus pT of all tracks with ITS; eta of track; p_T of track [MeV/c^2]; Tracks", {HistType::kTH2F, {{80, -2., 2.}, {100, 0., 5.}}});
      registry.add("tracks/etapt6", "eta versus pT of all tracks with TPC; eta of track; p_T of track [MeV/c^2]; Tracks", {HistType::kTH2F, {{80, -2., 2.}, {100, 0., 5.}}});
      registry.add("tracks/etapt7", "eta versus pT of all tracks with TRD; eta of track; p_T of track [MeV/c^2]; Tracks", {HistType::kTH2F, {{80, -2., 2.}, {100, 0., 5.}}});
      registry.add("tracks/etapt8", "eta versus pT of all tracks with TOF; eta of track; p_T of track [MeV/c^2]; Tracks", {HistType::kTH2F, {{80, -2., 2.}, {100, 0., 5.}}});

      registry.add("tracks/dEdxTPC", "TPC signal versus signed track momentum; Signed track momentum [GeV/c]; TPC signal; Tracks", {HistType::kTH2F, {{120, -6., 6.}, {1000, 0., 1000.}}});
      registry.add("tracks/dEdxTOF", "TOF signal versus signed track momentum; Signed track momentum [GeV/c]; TOF signal; Tracks", {HistType::kTH2F, {{120, -6., 6.}, {1000, 0., 30000.}}});

      // DG
      registry.add("DG/PVposxy", "DG: Vertex position in x and y direction; V_x [mm]; V_y [mm]; DG collisions", {HistType::kTH2F, {{100, -0.5, 0.5}, {100, -0.5, 0.5}}});
      registry.add("DG/PVposz", "DG: Vertex position in z direction; V_z; DG collisions", {HistType::kTH1F, {{1000, -100., 100.}}});
      registry.add("DG/netCharge", "DG: net charge; Net charge; DG collisions", {HistType::kTH1F, {{21, -10.5, 10.5}}});
      registry.add("DG/TrackStat", "Track bits as function of pT; Track pT; Track bit; Tracks", {HistType::kTH2F, {{100, 0., 5.}, {8, 0.5, 8.5}}});

      registry.add("DG/etapt", "DG: eta versus pT of all tracks; eta of track; p_T of track [GeV/c^2]; Tracks", {HistType::kTH2F, {{80, -2., 2.}, {100, 0., 5.}}});
      registry.add("DG/etapt2", "DG: eta versus pT of all quality tracks; eta of track; p_T of track [MeV/c^2]; Tracks", {HistType::kTH2F, {{80, -2., 2.}, {100, 0., 5.}}});
      registry.add("DG/etapt3", "DG: eta versus pT of all global tracks; eta of track; p_T of track [MeV/c^2]; Tracks", {HistType::kTH2F, {{80, -2., 2.}, {100, 0., 5.}}});
      registry.add("DG/etapt4", "DG: eta versus pT of all PV tracks; eta of track; p_T of track [MeV/c^2]; Tracks", {HistType::kTH2F, {{80, -2., 2.}, {100, 0., 5.}}});
      registry.add("DG/etapt5", "DG: eta versus pT of all tracks with ITS; eta of track; p_T of track [MeV/c^2]; Tracks", {HistType::kTH2F, {{80, -2., 2.}, {100, 0., 5.}}});
      registry.add("DG/etapt6", "DG: eta versus pT of all tracks with TPC; eta of track; p_T of track [MeV/c^2]; Tracks", {HistType::kTH2F, {{80, -2., 2.}, {100, 0., 5.}}});
      registry.add("DG/etapt7", "DG: eta versus pT of all tracks with TRD; eta of track; p_T of track [MeV/c^2]; Tracks", {HistType::kTH2F, {{80, -2., 2.}, {100, 0., 5.}}});
      registry.add("DG/etapt8", "DG: eta versus pT of all tracks with TOF; eta of track; p_T of track [MeV/c^2]; Tracks", {HistType::kTH2F, {{80, -2., 2.}, {100, 0., 5.}}});

      registry.add("DG/dEdxTPC", "DG: TPC signal versus signed track momentum; Signed track momentum [GeV/c]; TPC signal; Tracks", {HistType::kTH2F, {{120, -6., 6.}, {1000, 0., 1000.}}});
      registry.add("DG/dEdxTOF", "DG: TOF signal versus signed track momentum; Signed track momentum [GeV/c]; TOF signal; Tracks", {HistType::kTH2F, {{120, -6., 6.}, {1000, 0., 30000.}}});
      registry.add("DG/IVMptSys", "DG: Invariant mass versus p_{T, system}; Invarian mass [GeV/c^2]; p_{T, system} [GeV/c]; DG collisions", {HistType::kTH2F, {{100, 0., 5.}, {350, 0., 3.5}}});
      registry.add("DG/IVMptTrk", "DG: Invariant mass versus p_{T, tracks}; Invarian mass [GeV/c^2]; p_{T, tracks} [GeV/c]; DG collisions", {HistType::kTH2F, {{100, 0., 5.}, {350, 0., 3.5}}});
    }
    if (context.mOptions.get<bool>("processCleanFT0")) {
      registry.add("cleanFT0/Stat", "Statistics of collisions with clean FT0; FT0 status; Collisions", {HistType::kTH1F, {{2, 0.5, 2.5}}});
      registry.add("cleanFT0/PVTracks", "Distribution of number of PV contributors for all collisions and for collisions with clean FT0; Number of PV contributors;", {HistType::kTH2F, {{100, 0.5, 100.5}, {2, 0.5, 2.5}}});
    }
    if (context.mOptions.get<bool>("processCleanFIT1")) {
      registry.add("cleanFIT1/Stat", "Statistics of collisions with empty FT0; Multiple of collision time resolution; FT0 status; Collisions", {HistType::kTH2F, {{20, -0.5, 19.5}, {2, -0.5, 1.5}}});
      registry.add("cleanFIT1/FV0Aamp", "Amplitude of FV0A in collisions with empty FT0; Multiple of collision time resolution; FV0A amplitude; Collisions", {HistType::kTH2F, {{20, -0.5, 19.5}, {1000, -0.5, 999.5}}});
      registry.add("cleanFIT1/FT0Aamp", "Amplitude of FT0A in collisions with empty FT0; Multiple of collision time resolution; FT0A amplitude; Collisions", {HistType::kTH2F, {{20, -0.5, 19.5}, {1000, -0.5, 999.5}}});
      registry.add("cleanFIT1/FT0Camp", "Amplitude of FT0C in collisions with empty FT0; Multiple of collision time resolution; FT0C amplitude; Collisions", {HistType::kTH2F, {{20, -0.5, 19.5}, {1000, -0.5, 999.5}}});
      registry.add("cleanFIT1/FDDAamp", "Amplitude of FDDA in collisions with empty FT0; Multiple of collision time resolution; FDDA amplitude; Collisions", {HistType::kTH2F, {{20, -0.5, 19.5}, {1000, -0.5, 999.5}}});
      registry.add("cleanFIT1/FDDCamp", "Amplitude of FDDC in collisions with empty FT0; Multiple of collision time resolution; FDDC amplitude; Collisions", {HistType::kTH2F, {{20, -0.5, 19.5}, {1000, -0.5, 999.5}}});
    }
    if (context.mOptions.get<bool>("processCleanFIT2")) {
      registry.add("cleanFIT2/Stat", "Statistics of collisions with empty FIT; Number of neighbouring BCs; FIT status; Collisions", {HistType::kTH2F, {{20, -0.5, 19.5}, {2, -0.5, 1.5}}});
      registry.add("cleanFIT2/FV0Aamp", "Amplitude of FV0A in collisions with empty FIT; Number of neighbouring BCs; FV0A amplitude; Collisions", {HistType::kTH2F, {{20, -0.5, 19.5}, {1000, -0.5, 999.5}}});
      registry.add("cleanFIT2/FT0Aamp", "Amplitude of FT0A in collisions with empty FIT; Number of neighbouring BCs; FT0A amplitude; Collisions", {HistType::kTH2F, {{20, -0.5, 19.5}, {1000, -0.5, 999.5}}});
      registry.add("cleanFIT2/FT0Camp", "Amplitude of FT0C in collisions with empty FIT; Number of neighbouring BCs; FT0C amplitude; Collisions", {HistType::kTH2F, {{20, -0.5, 19.5}, {1000, -0.5, 999.5}}});
      registry.add("cleanFIT2/FDDAamp", "Amplitude of FDDA in collisions with empty FIT; Number of neighbouring BCs; FDDA amplitude; Collisions", {HistType::kTH2F, {{20, -0.5, 19.5}, {1000, -0.5, 999.5}}});
      registry.add("cleanFIT2/FDDCamp", "Amplitude of FDDC in collisions with empty FIT; Number of neighbouring BCs; FDDC amplitude; Collisions", {HistType::kTH2F, {{20, -0.5, 19.5}, {1000, -0.5, 999.5}}});
    }
    if (context.mOptions.get<bool>("processFV0")) {
      registry.add("FV0/FV0Aamp", "FV0A amplitudes; FV0A channel; FV0A amplitude; Entries", {HistType::kTH2F, {{48, -0.5, 47.5}, {2000, 0., 2000.}}});
      registry.add("FV0/emptyBCs", "Distribution of number of consecutive BCs with empty FV0A; Number of consecutive BCs with empty FV0A; Entries", {HistType::kTH1F, {{3564, -0.5, 3563.5}}});
    }
    if (context.mOptions.get<bool>("processFT0")) {
      registry.add("FT0/FT0Aamp", "FT0A amplitudes; FT0A channel; FT0A amplitude; Entries", {HistType::kTH2F, {{96, -0.5, 95.5}, {400, 0., 400.}}});
      registry.add("FT0/FT0Camp", "FT0C amplitudes; FT0C channel; FT0C amplitude; Entries", {HistType::kTH2F, {{112, -0.5, 111.5}, {400, 0., 400.}}});
      registry.add("FT0/AP2BC", "P2 BCs with FT0A signal; P2 BC; Entries", {HistType::kTH1F, {{3564, -0.5, 3563.5}}});
      registry.add("FT0/CP2BC", "P2 BCs with FT0C signal; P2 BC; Entries", {HistType::kTH1F, {{3564, -0.5, 3563.5}}});
      registry.add("FT0/emptyBCs", "Distribution of number of consecutive BCs with empty FT0; Number of consecutive BCs with empty FT0; Entries", {HistType::kTH1F, {{3564, -0.5, 3563.5}}});

      // add amplitude histograms
      for (auto n{0}; n <= ns; n++) {
        registry.add(hFT0A[n].data(), hFT0A[n].data(), {HistType::kTH1F, {{1000, 0., 1000.}}});
        registry.add(hFT0C[n].data(), hFT0C[n].data(), {HistType::kTH1F, {{1000, 0., 1000.}}});
      }
    }
    if (context.mOptions.get<bool>("processFDD")) {
      registry.add("FDD/FDDAamp", "FDDA amplitudes; FDDA channel; FDDA amplitude; Entries", {HistType::kTH2F, {{8, -0.5, 7.5}, {100, 0., 100.}}});
      registry.add("FDD/FDDCamp", "FDDA amplitudes; FDDA channel; FDDA amplitude; Entries", {HistType::kTH2F, {{8, -0.5, 7.5}, {100, 0., 100.}}});
      registry.add("FDD/emptyBCs", "Distribution of number of consecutive BCs with empty FDD; Number of consecutive BCs with empty FDD; Entries", {HistType::kTH1F, {{3564, -0.5, 3563.5}}});
    }
    if (context.mOptions.get<bool>("processZDC")) {
      registry.add("ZDC/Energies", "Registered energies in various ZDC channels; ZDC channel; Energy; Entries", {HistType::kTH2F, {{22, -0.5, 21.5}, {100, 0., 1000.}}});
    }
    // if (context.mOptions.get<bool>("processCalo")) {
    //   registry.add("CaloCell", "#CaloCell", {HistType::kTH1I, {{18000, -0.5, 17999.5}}});
    //   registry.add("CaloAmplitude", "#CaloAmplitude", {HistType::kTH1F, {{100, 0, 10.}}});
    // }
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
      LOGF(debug, "<DiffQA> size of ambiguous tracks table %i", ambtracks.size());
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
      LOGF(debug, "<DiffQA> size of ambiguous fwd tracks table %i", ambfwdtracks.size());
      for (auto const& ambfwdtrack : ambfwdtracks) {
        auto bcfirst = ambfwdtrack.bc().rawIteratorAt(0);
        auto bclast = ambfwdtrack.bc().rawIteratorAt(ambfwdtrack.bc().size() - 1);
        afbcrs.add(bcfirst.globalIndex(), bclast.globalIndex());
      }
      afbcrs.merge();
    }
  }

  // ...............................................................................................................
  void processMain(CC const& collision, BCs const& bct0s,
                   TCs const& tracks, FWs const& fwdtracks, ATs const& ambtracks, AFTs const& ambfwdtracks,
                   aod::FT0s const& ft0s, aod::FV0As const& fv0as, aod::FDDs const& fdds,
                   aod::Zdcs& zdcs, aod::Calos& calos,
                   aod::V0s const& v0s, aod::Cascades const& cascades)
  {
    LOGF(debug, "<DiffQA> Collision %d", collision.globalIndex());
    LOGF(debug, "<DiffQA> Start %i", abcrs.size());

    bool isDGcandidate = true;
    registry.get<TH1>(HIST("collisions/Stat"))->Fill(0., isDGcandidate * 1.);

    // update collision histograms
    // vertex position
    registry.get<TH2>(HIST("collisions/PVposxy"))->Fill(collision.posX(), collision.posY());
    registry.get<TH1>(HIST("collisions/PVposz"))->Fill(collision.posZ());
    // tracks
    registry.get<TH1>(HIST("collisions/Tracks"))->Fill(tracks.size());
    // vertex tracks
    registry.get<TH1>(HIST("collisions/PVTracks"))->Fill(collision.numContrib());
    // global tracks
    Partition<TCs> goodTracks = requireGlobalTrackInFilter();
    goodTracks.bindTable(tracks);
    registry.get<TH1>(HIST("collisions/globalTracks"))->Fill(goodTracks.size());

    // loop over all tracks
    float rgtrwTOF = 0.;
    for (auto const& track : tracks) {
      // update track stats
      registry.get<TH2>(HIST("tracks/Stat"))->Fill(track.pt(), 1., 1.);
      registry.get<TH2>(HIST("tracks/Stat"))->Fill(track.pt(), 2., track.isQualityTrack() * 1.);
      registry.get<TH2>(HIST("tracks/Stat"))->Fill(track.pt(), 3., track.isGlobalTrack() * 1.);
      registry.get<TH2>(HIST("tracks/Stat"))->Fill(track.pt(), 4., track.isPVContributor() * 1.);
      registry.get<TH2>(HIST("tracks/Stat"))->Fill(track.pt(), 5., track.hasITS() * 1.);
      registry.get<TH2>(HIST("tracks/Stat"))->Fill(track.pt(), 6., track.hasTPC() * 1.);
      registry.get<TH2>(HIST("tracks/Stat"))->Fill(track.pt(), 7., track.hasTRD() * 1.);
      registry.get<TH2>(HIST("tracks/Stat"))->Fill(track.pt(), 8., track.hasTOF() * 1.);

      // update eta vs pt histograms
      registry.get<TH2>(HIST("tracks/etapt"))->Fill(track.eta(), track.pt(), 1.);
      registry.get<TH2>(HIST("tracks/etapt2"))->Fill(track.eta(), track.pt(), track.isQualityTrack() * 1.);
      registry.get<TH2>(HIST("tracks/etapt3"))->Fill(track.eta(), track.pt(), track.isGlobalTrack() * 1.);
      registry.get<TH2>(HIST("tracks/etapt4"))->Fill(track.eta(), track.pt(), track.isPVContributor() * 1.);
      registry.get<TH2>(HIST("tracks/etapt5"))->Fill(track.eta(), track.pt(), track.hasITS() * 1.);
      registry.get<TH2>(HIST("tracks/etapt6"))->Fill(track.eta(), track.pt(), track.hasTPC() * 1.);
      registry.get<TH2>(HIST("tracks/etapt7"))->Fill(track.eta(), track.pt(), track.hasTRD() * 1.);
      registry.get<TH2>(HIST("tracks/etapt8"))->Fill(track.eta(), track.pt(), track.hasTOF() * 1.);

      // update dEdx histograms
      registry.get<TH2>(HIST("tracks/dEdxTPC"))->Fill(track.tpcInnerParam() * track.sign(), track.tpcSignal());
      if (track.tpcSignal() > maxdEdxTPC) {
        maxdEdxTPC = track.tpcSignal();
        LOGF(debug, "<DiffQA> New maxdEdx TPC %f", maxdEdxTPC);
      }

      // TOF hit?
      if (track.hasTOF()) {
        registry.get<TH2>(HIST("tracks/dEdxTOF"))->Fill(track.p() * track.sign(), track.tofSignal());
        if (track.tofSignal() > maxdEdxTOF) {
          maxdEdxTOF = track.tofSignal();
          LOGF(debug, "<DiffQA> New maxdEdx TOF %f", maxdEdxTOF);
        }

        // vertex track with TOF hit?
        if (track.isPVContributor()) {
          rgtrwTOF += 1.;
        }
      }
    }
    // fraction of PV tracks with TOF hit
    if (collision.numContrib() > 0) {
      rgtrwTOF /= collision.numContrib();
    }
    LOGF(debug, "<DiffQA> PV tracks with TOF: %f [1]", rgtrwTOF);
    registry.get<TH2>(HIST("collisions/tResvsrTOFTracks"))->Fill(collision.collisionTimeRes(), rgtrwTOF);

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
    auto bcSlice = udhelpers::compatibleBCs(collision, diffCuts.NDtcoll(), bct0s, diffCuts.minNBCs());

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
    registry.get<TH1>(HIST("collisions/Stat"))->Fill(1., isDGcandidate * 1.);

    // no Zdc signal in bcSlice
    std::vector<float> lims(10, 0.);
    for (auto const& bc : bcSlice) {
      if (!udhelpers::cleanZDC(bc, zdcs, lims, cache)) {
        isDGcandidate = false;
        break;
      }
    }
    registry.get<TH1>(HIST("collisions/Stat"))->Fill(2., isDGcandidate * 1.);

    // no Calo signal in bcSlice
    for (auto const& bc : bcSlice) {
      if (!udhelpers::cleanCalo(bc, calos, lims, cache)) {
        isDGcandidate = false;
        break;
      }
    }
    registry.get<TH1>(HIST("collisions/Stat"))->Fill(3., isDGcandidate * 1.);

    // no V0s
    isDGcandidate &= (v0s.size() == 0);
    registry.get<TH1>(HIST("collisions/Stat"))->Fill(4., isDGcandidate * 1.);

    // no Cascades
    isDGcandidate &= (cascades.size() == 0);
    registry.get<TH1>(HIST("collisions/Stat"))->Fill(5., isDGcandidate * 1.);

    // number of forward tracks = 0
    isDGcandidate &= (fwdtracks.size() == 0);
    registry.get<TH1>(HIST("collisions/Stat"))->Fill(6., isDGcandidate * 1.);

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
    registry.get<TH1>(HIST("collisions/Stat"))->Fill(7., globalAndVtx * 1.);
    registry.get<TH1>(HIST("collisions/Stat"))->Fill(8., vtxAndGlobal * 1.);
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
    registry.get<TH1>(HIST("collisions/Stat"))->Fill(9., noAmbTracks * 1.);

    // check a given bc for possible ambiguous FwdTracks
    auto noAmbFwdTracks = isDGcandidate;
    for (auto const& bc : bcSlice) {
      if (afbcrs.isInRange(bc.globalIndex())) {
        noAmbFwdTracks = false;
        break;
      }
    }
    registry.get<TH1>(HIST("collisions/Stat"))->Fill(10., noAmbFwdTracks * 1.);

    // fraction of PV tracks with TOF hit
    isDGcandidate &= (rgtrwTOF >= diffCuts.minRgtrwTOF());
    registry.get<TH1>(HIST("collisions/Stat"))->Fill(11., isDGcandidate * 1.);

    // number of vertex tracks <= n
    isDGcandidate &= (collision.numContrib() >= diffCuts.minNTracks());
    registry.get<TH1>(HIST("collisions/Stat"))->Fill(12., isDGcandidate * 1.);
    isDGcandidate &= (collision.numContrib() <= diffCuts.maxNTracks());
    registry.get<TH1>(HIST("collisions/Stat"))->Fill(13., isDGcandidate * 1.);

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
      for (auto const& track : tracks) {
        // update histogram notPVTracks
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
    registry.get<TH1>(HIST("collisions/Stat"))->Fill(14., isDGcandidate * 1.);
    isDGcandidate &= goodetas;
    registry.get<TH1>(HIST("collisions/Stat"))->Fill(15., isDGcandidate * 1.);
    auto netChargeValues = diffCuts.netCharges();
    if (std::find(netChargeValues.begin(), netChargeValues.end(), netCharge) == netChargeValues.end()) {
      goodnchs = false;
    }
    isDGcandidate &= goodnchs;
    registry.get<TH1>(HIST("collisions/Stat"))->Fill(16., isDGcandidate * 1.);
    registry.get<TH1>(HIST("collisions/Stat"))->Fill(17., isDGcandidate * 1.);
    isDGcandidate &= (ivm.M() >= diffCuts.minIVM());
    registry.get<TH1>(HIST("collisions/Stat"))->Fill(18., isDGcandidate * 1.);
    isDGcandidate &= (ivm.M() <= diffCuts.maxIVM());
    registry.get<TH1>(HIST("collisions/Stat"))->Fill(19., isDGcandidate * 1.);

    // update some DG histograms
    if (isDGcandidate) {
      // vertex position of DG events
      registry.get<TH2>(HIST("DG/PVposxy"))->Fill(collision.posX(), collision.posY());
      registry.get<TH1>(HIST("DG/PVposz"))->Fill(collision.posZ());
      registry.get<TH1>(HIST("DG/netCharge"))->Fill(netCharge);
      registry.get<TH2>(HIST("DG/IVMptSys"))->Fill(ivm.M(), ivm.Perp());

      // fill track status, eta vs pt, and dEdx of DG event tracks
      for (auto const& track : tracks) {
        if (track.isPVContributor()) {
          registry.get<TH2>(HIST("DG/TrackStat"))->Fill(track.pt(), 1., 1.);
          registry.get<TH2>(HIST("DG/TrackStat"))->Fill(track.pt(), 2., track.isQualityTrack() * 1.);
          registry.get<TH2>(HIST("DG/TrackStat"))->Fill(track.pt(), 3., track.isGlobalTrack() * 1.);
          registry.get<TH2>(HIST("DG/TrackStat"))->Fill(track.pt(), 4., 1.);
          registry.get<TH2>(HIST("DG/TrackStat"))->Fill(track.pt(), 5., track.hasITS() * 1.);
          registry.get<TH2>(HIST("DG/TrackStat"))->Fill(track.pt(), 6., track.hasTPC() * 1.);
          registry.get<TH2>(HIST("DG/TrackStat"))->Fill(track.pt(), 7., track.hasTRD() * 1.);
          registry.get<TH2>(HIST("DG/TrackStat"))->Fill(track.pt(), 8., track.hasTOF() * 1.);

          registry.get<TH2>(HIST("DG/etapt"))->Fill(track.eta(), track.pt(), 1.);
          registry.get<TH2>(HIST("DG/etapt2"))->Fill(track.eta(), track.pt(), track.isQualityTrack() * 1.);
          registry.get<TH2>(HIST("DG/etapt3"))->Fill(track.eta(), track.pt(), track.isGlobalTrack() * 1.);
          registry.get<TH2>(HIST("DG/etapt4"))->Fill(track.eta(), track.pt(), 1.);
          registry.get<TH2>(HIST("DG/etapt5"))->Fill(track.eta(), track.pt(), track.hasITS() * 1.);
          registry.get<TH2>(HIST("DG/etapt6"))->Fill(track.eta(), track.pt(), track.hasTPC() * 1.);
          registry.get<TH2>(HIST("DG/etapt7"))->Fill(track.eta(), track.pt(), track.hasTRD() * 1.);
          registry.get<TH2>(HIST("DG/etapt8"))->Fill(track.eta(), track.pt(), track.hasTOF() * 1.);

          LOGF(debug, "dEdx TPC %f TOF %i %f", track.tpcSignal(), track.hasTOF(), track.hasTOF() ? track.tofSignal() : 0.);
          registry.get<TH2>(HIST("DG/dEdxTPC"))->Fill(track.tpcInnerParam() * track.sign(), track.tpcSignal());
          if (track.hasTOF()) {
            registry.get<TH2>(HIST("DG/dEdxTOF"))->Fill(track.p() * track.sign(), track.tofSignal());
          }
          registry.get<TH2>(HIST("DG/IVMptTrk"))->Fill(ivm.M(), track.pt());
        }
      }
    }
  }
  PROCESS_SWITCH(DiffQA, processMain, "Process Main", true);

  // ...............................................................................................................
  // Fraction of collisions with empty FIT as function of NDtcoll
  void processCleanFIT1(CC const& collision, BCs const& bct0s,
                        aod::FT0s const& ft0s, aod::FV0As const& fv0as, aod::FDDs const& fdds)
  {
    LOGF(debug, "<CleanFit1. Collision %d", collision.globalIndex());

    // test influence of BCrange width using a series of NDtcoll
    float ampFV0A, ampFT0A, ampFT0C, ampFDDA, ampFDDC;
    bool isDGcandidate;
    for (int NDtcoll = 0; NDtcoll < 20; NDtcoll++) {
      auto bcSlice = udhelpers::compatibleBCs(collision, NDtcoll, bct0s, 0);
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
      registry.get<TH2>(HIST("cleanFIT1/Stat"))->Fill(NDtcoll, isDGcandidate * 1.);
      if (isDGcandidate) {
        registry.get<TH2>(HIST("cleanFIT1/FV0Aamp"))->Fill(NDtcoll, ampFV0A);
        registry.get<TH2>(HIST("cleanFIT1/FT0Aamp"))->Fill(NDtcoll, ampFT0A);
        registry.get<TH2>(HIST("cleanFIT1/FT0Camp"))->Fill(NDtcoll, ampFT0C);
        registry.get<TH2>(HIST("cleanFIT1/FDDAamp"))->Fill(NDtcoll, ampFDDA);
        registry.get<TH2>(HIST("cleanFIT1/FDDCamp"))->Fill(NDtcoll, ampFDDC);
      }
    }
  }
  PROCESS_SWITCH(DiffQA, processCleanFIT1, "Process CleanFit1", true);

  // ...............................................................................................................
  void processCleanFIT2(CC const& collision, BCs const& bct0s,
                        aod::FT0s const& ft0s, aod::FV0As const& fv0as, aod::FDDs const& fdds)
  {
    LOGF(debug, "<CleanFit2. Collision %d", collision.globalIndex());

    // test influence of BCrange width using a series of nMinBC
    float ampFV0A, ampFT0A, ampFT0C, ampFDDA, ampFDDC;
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
      registry.get<TH2>(HIST("cleanFIT2/Stat"))->Fill(nMinBC, isDGcandidate * 1.);
      if (isDGcandidate) {
        registry.get<TH2>(HIST("cleanFIT2/FV0Aamp"))->Fill(nMinBC, ampFV0A);
        registry.get<TH2>(HIST("cleanFIT2/FT0Aamp"))->Fill(nMinBC, ampFT0A);
        registry.get<TH2>(HIST("cleanFIT2/FT0Camp"))->Fill(nMinBC, ampFT0C);
        registry.get<TH2>(HIST("cleanFIT2/FDDAamp"))->Fill(nMinBC, ampFDDA);
        registry.get<TH2>(HIST("cleanFIT2/FDDCamp"))->Fill(nMinBC, ampFDDC);
      }
    }
  }
  PROCESS_SWITCH(DiffQA, processCleanFIT2, "Process CleanFit2", true);

  // ...............................................................................................................
  // Distribution of number of PV contributors for all collisions and those with empty FT0
  void processCleanFT0(CC const& collision, BCs const& bct0s,
                       aod::FT0s const& ft0s, aod::FV0As const& fv0as, aod::FDDs const& fdds)
  {
    // count collisions
    registry.get<TH1>(HIST("cleanFT0/Stat"))->Fill(1., 1.);
    registry.get<TH2>(HIST("cleanFT0/PVTracks"))->Fill(collision.numContrib(), 1., 1.);

    // check FT0 to be empty
    auto bc = collision.foundBC_as<BCs>();
    if (udhelpers::cleanFT0(bc, diffCuts.maxFITtime(), 0., 0.)) {
      // only collisions with empty FT0 arrive here
      registry.get<TH1>(HIST("cleanFT0/Stat"))->Fill(2., 1.);

      // update #PV contributors in collisions with empty FT0
      registry.get<TH2>(HIST("cleanFT0/PVTracks"))->Fill(collision.numContrib(), 2., 1.);
    }
  }
  PROCESS_SWITCH(DiffQA, processCleanFT0, "Process CleanFIT1", true);

  // ...............................................................................................................
  void processFV0(aod::FV0As const& fv0s, BCs const&)
  {
    LOGF(info, "<FV0Signals> %d", fv0s.size());
    if (fv0s.size() <= 0) {
      return;
    }
    int64_t lastBCwFV0 = fv0s.begin().bc_as<BCs>().globalBC();
    auto lastOrbit = lastBCwFV0 / o2::constants::lhc::LHCMaxBunches;

    for (auto fv0 : fv0s) {

      // side A
      for (size_t ind = 0; ind < fv0.channel().size(); ind++) {
        registry.get<TH2>(HIST("FV0/FV0Aamp"))->Fill((fv0.channel())[ind], (fv0.amplitude())[ind]);
      }

      // length of series of empty BCs
      auto bc = fv0.bc_as<BCs>();
      int64_t aBC = bc.globalBC();
      auto aOrbit = aBC / o2::constants::lhc::LHCMaxBunches;
      auto ampA = udhelpers::FV0AmplitudeA(bc.foundFV0());
      if (ampA > 0.) {
        // require both BCs to be in same orbit
        if (aOrbit == lastOrbit)
          registry.get<TH1>(HIST("FV0/emptyBCs"))->Fill((bc.globalBC() - lastBCwFV0));
        lastBCwFV0 = aBC;
        lastOrbit = aOrbit;
      }
    }
  };
  PROCESS_SWITCH(DiffQA, processFV0, "Process FV0", true);

  // ...............................................................................................................
  void processFT0(aod::FT0s const& ft0s, BCs const&)
  {
    LOGF(debug, "<processFT0> %d", ft0s.size());
    int nc = 0;
    int64_t fBC = 0; // first BC with FIT activity
    int64_t aBC = 0; // actually processed BC
    float minAmpA = 0., fAmpA = 0.;
    float minAmpC = 0., fAmpC = 0.;

    int64_t lastBCwFT0 = ft0s.begin().bc_as<BCs>().globalBC();
    auto lastOrbit = lastBCwFT0 / o2::constants::lhc::LHCMaxBunches;
    for (auto ft0 : ft0s) {

      // side A
      for (size_t ind = 0; ind < ft0.channelA().size(); ind++) {
        registry.get<TH2>(HIST("FT0/FT0Aamp"))->Fill((ft0.channelA())[ind], (ft0.amplitudeA())[ind]);
      }

      // side C
      for (size_t ind = 0; ind < ft0.channelC().size(); ind++) {
        registry.get<TH2>(HIST("FT0/FT0Camp"))->Fill((ft0.channelC())[ind], (ft0.amplitudeC())[ind]);
      }

      // sequence of BCs
      auto bc = ft0.bc_as<BCs>();
      aBC = bc.globalBC();
      auto gBC = aBC % o2::constants::lhc::LHCMaxBunches;
      auto ampA = udhelpers::FT0AmplitudeA(bc.foundFT0());
      auto ampC = udhelpers::FT0AmplitudeC(bc.foundFT0());

      // update AP2BC
      if (ampA > 0.) {
        registry.get<TH1>(HIST("FT0/AP2BC"))->Fill(gBC, 1.);
      }

      // update FT0/CP2BC
      if (ampC > 0.) {
        registry.get<TH1>(HIST("FT0/CP2BC"))->Fill(gBC, 1.);
      }

      // update dFT0BCNUM
      // require both BCs to be in same orbit
      auto aOrbit = aBC / o2::constants::lhc::LHCMaxBunches;
      LOGF(debug, "lastOrbit %d aOrbit %d", lastOrbit, aOrbit);
      if (ampA > 0. || ampC > 0.) {
        if (aOrbit == lastOrbit)
          registry.get<TH1>(HIST("FT0/emptyBCs"))->Fill((aBC - lastBCwFT0));
      } else {
        continue;
      }

      // amplitude distributions in BCs following a long series of empty BCs
      nc = aBC - lastBCwFT0 - 1;
      if (nc >= ncmin) {
        fBC = aBC;
        fAmpA = ampA;
        fAmpC = ampC;
      }
      auto dBC = static_cast<int>(aBC - fBC);
      if (dBC >= 0 && dBC <= ns) {
        LOGF(debug, "<processFT0> dBC %d ampA %f ampC %f", dBC, ampA, ampC);
        switch (dBC) {
          case 0:
            if (fAmpA > minAmpA)
              registry.get<TH1>(HIST(hFT0A[0].data()))->Fill(ampA, 1.);
            if (fAmpC > minAmpC)
              registry.get<TH1>(HIST(hFT0C[0].data()))->Fill(ampC, 1.);
            break;
          case 1:
            if (fAmpA > minAmpA)
              registry.get<TH1>(HIST(hFT0A[1].data()))->Fill(ampA, 1.);
            if (fAmpC > minAmpC)
              registry.get<TH1>(HIST(hFT0C[1].data()))->Fill(ampC, 1.);
            break;
          case 2:
            if (fAmpA > minAmpA)
              registry.get<TH1>(HIST(hFT0A[2].data()))->Fill(ampA, 1.);
            if (fAmpC > minAmpC)
              registry.get<TH1>(HIST(hFT0C[2].data()))->Fill(ampC, 1.);
            break;
          case 3:
            if (fAmpA > minAmpA)
              registry.get<TH1>(HIST(hFT0A[3].data()))->Fill(ampA, 1.);
            if (fAmpC > minAmpC)
              registry.get<TH1>(HIST(hFT0C[3].data()))->Fill(ampC, 1.);
            break;
          case 4:
            if (fAmpA > minAmpA)
              registry.get<TH1>(HIST(hFT0A[4].data()))->Fill(ampA, 1.);
            if (fAmpC > minAmpC)
              registry.get<TH1>(HIST(hFT0C[4].data()))->Fill(ampC, 1.);
            break;
          case 5:
            if (fAmpA > minAmpA)
              registry.get<TH1>(HIST(hFT0A[5].data()))->Fill(ampA, 1.);
            if (fAmpC > minAmpC)
              registry.get<TH1>(HIST(hFT0C[5].data()))->Fill(ampC, 1.);
            break;
          case 6:
            if (fAmpA > minAmpA)
              registry.get<TH1>(HIST(hFT0A[6].data()))->Fill(ampA, 1.);
            if (fAmpC > minAmpC)
              registry.get<TH1>(HIST(hFT0C[6].data()))->Fill(ampC, 1.);
            break;
          case 7:
            if (fAmpA > minAmpA)
              registry.get<TH1>(HIST(hFT0A[7].data()))->Fill(ampA, 1.);
            if (fAmpC > minAmpC)
              registry.get<TH1>(HIST(hFT0C[7].data()))->Fill(ampC, 1.);
            break;
          case 8:
            if (fAmpA > minAmpA)
              registry.get<TH1>(HIST(hFT0A[8].data()))->Fill(ampA, 1.);
            if (fAmpC > minAmpC)
              registry.get<TH1>(HIST(hFT0C[8].data()))->Fill(ampC, 1.);
            break;
          case 9:
            if (fAmpA > minAmpA)
              registry.get<TH1>(HIST(hFT0A[9].data()))->Fill(ampA, 1.);
            if (fAmpC > minAmpC)
              registry.get<TH1>(HIST(hFT0C[9].data()))->Fill(ampC, 1.);
            break;
          case 10:
            if (fAmpA > minAmpA)
              registry.get<TH1>(HIST(hFT0A[10].data()))->Fill(ampA, 1.);
            if (fAmpC > minAmpC)
              registry.get<TH1>(HIST(hFT0C[10].data()))->Fill(ampC, 1.);
            break;
          case 11:
            if (fAmpA > minAmpA)
              registry.get<TH1>(HIST(hFT0A[11].data()))->Fill(ampA, 1.);
            if (fAmpC > minAmpC)
              registry.get<TH1>(HIST(hFT0C[11].data()))->Fill(ampC, 1.);
            break;
          case 12:
            if (fAmpA > minAmpA)
              registry.get<TH1>(HIST(hFT0A[12].data()))->Fill(ampA, 1.);
            if (fAmpC > minAmpC)
              registry.get<TH1>(HIST(hFT0C[12].data()))->Fill(ampC, 1.);
            break;
          case 13:
            if (fAmpA > minAmpA)
              registry.get<TH1>(HIST(hFT0A[13].data()))->Fill(ampA, 1.);
            if (fAmpC > minAmpC)
              registry.get<TH1>(HIST(hFT0C[13].data()))->Fill(ampC, 1.);
            break;
          case 14:
            if (fAmpA > minAmpA)
              registry.get<TH1>(HIST(hFT0A[14].data()))->Fill(ampA, 1.);
            if (fAmpC > minAmpC)
              registry.get<TH1>(HIST(hFT0C[14].data()))->Fill(ampC, 1.);
            break;
          case 15:
            if (fAmpA > minAmpA)
              registry.get<TH1>(HIST(hFT0A[15].data()))->Fill(ampA, 1.);
            if (fAmpC > minAmpC)
              registry.get<TH1>(HIST(hFT0C[15].data()))->Fill(ampC, 1.);
            break;
          case 16:
            if (fAmpA > minAmpA)
              registry.get<TH1>(HIST(hFT0A[16].data()))->Fill(ampA, 1.);
            if (fAmpC > minAmpC)
              registry.get<TH1>(HIST(hFT0C[16].data()))->Fill(ampC, 1.);
            break;
          case 17:
            if (fAmpA > minAmpA)
              registry.get<TH1>(HIST(hFT0A[17].data()))->Fill(ampA, 1.);
            if (fAmpC > minAmpC)
              registry.get<TH1>(HIST(hFT0C[17].data()))->Fill(ampC, 1.);
            break;
          case 18:
            if (fAmpA > minAmpA)
              registry.get<TH1>(HIST(hFT0A[18].data()))->Fill(ampA, 1.);
            if (fAmpC > minAmpC)
              registry.get<TH1>(HIST(hFT0C[18].data()))->Fill(ampC, 1.);
            break;
          case 19:
            if (fAmpA > minAmpA)
              registry.get<TH1>(HIST(hFT0A[19].data()))->Fill(ampA, 1.);
            if (fAmpC > minAmpC)
              registry.get<TH1>(HIST(hFT0C[19].data()))->Fill(ampC, 1.);
            break;
          case 20:
            if (fAmpA > minAmpA)
              registry.get<TH1>(HIST(hFT0A[20].data()))->Fill(ampA, 1.);
            if (fAmpC > minAmpC)
              registry.get<TH1>(HIST(hFT0C[20].data()))->Fill(ampC, 1.);
        }
      }
      lastBCwFT0 = aBC;
      lastOrbit = aOrbit;
    }
  };
  PROCESS_SWITCH(DiffQA, processFT0, "Process FT0", true);

  // ...............................................................................................................
  void processFDD(aod::FDDs const& fdds, BCs const&)
  {
    LOGF(debug, "<FDDSignals> %d", fdds.size());

    int64_t lastBCwFDD = fdds.begin().bc_as<BCs>().globalBC();
    auto lastOrbit = lastBCwFDD / o2::constants::lhc::LHCMaxBunches;
    for (auto fdd : fdds) {

      // side A
      for (auto ind = 0; ind < 8; ind++) {
        registry.get<TH2>(HIST("FDD/FDDAamp"))->Fill(ind, (fdd.chargeA())[ind]);
      }

      // side C
      for (auto ind = 0; ind < 8; ind++) {
        registry.get<TH2>(HIST("FDD/FDDCamp"))->Fill(ind, (fdd.chargeC())[ind]);
      }

      // sequence of BCs
      auto bc = fdd.bc_as<BCs>();
      auto aBC = bc.globalBC();
      int64_t aOrbit = aBC / o2::constants::lhc::LHCMaxBunches;
      auto ampA = udhelpers::FDDAmplitudeA(bc.foundFDD());
      auto ampC = udhelpers::FDDAmplitudeC(bc.foundFDD());
      if (ampA > 0. || ampC > 0.) {
        // require both BCs to be in same orbit
        if (aOrbit == lastOrbit)
          registry.get<TH1>(HIST("FDD/emptyBCs"))->Fill((bc.globalBC() - lastBCwFDD));
        lastBCwFDD = aBC;
        lastOrbit = aOrbit;
      }
    }
  };
  PROCESS_SWITCH(DiffQA, processFDD, "Process FDD", true);

  // ...............................................................................................................
  void processZDC(aod::Zdc const& zdc)
  {
    LOGF(debug, "<ZDCSignals> %d", zdc.size());

    // Zdc energies
    registry.get<TH2>(HIST("ZDC/Energies"))->Fill(0., zdc.energyZEM1());
    registry.get<TH2>(HIST("ZDC/Energies"))->Fill(1., zdc.energyZEM2());
    registry.get<TH2>(HIST("ZDC/Energies"))->Fill(2., zdc.energyCommonZNA());
    registry.get<TH2>(HIST("ZDC/Energies"))->Fill(3., zdc.energyCommonZNC());
    registry.get<TH2>(HIST("ZDC/Energies"))->Fill(4., zdc.energyCommonZPA());
    registry.get<TH2>(HIST("ZDC/Energies"))->Fill(5., zdc.energyCommonZPC());
    registry.get<TH2>(HIST("ZDC/Energies"))->Fill(6., (zdc.energySectorZNA())[0]);
    registry.get<TH2>(HIST("ZDC/Energies"))->Fill(7., (zdc.energySectorZNA())[1]);
    registry.get<TH2>(HIST("ZDC/Energies"))->Fill(8., (zdc.energySectorZNA())[2]);
    registry.get<TH2>(HIST("ZDC/Energies"))->Fill(9., (zdc.energySectorZNA())[3]);
    registry.get<TH2>(HIST("ZDC/Energies"))->Fill(10., (zdc.energySectorZNC())[0]);
    registry.get<TH2>(HIST("ZDC/Energies"))->Fill(11., (zdc.energySectorZNC())[1]);
    registry.get<TH2>(HIST("ZDC/Energies"))->Fill(12., (zdc.energySectorZNC())[2]);
    registry.get<TH2>(HIST("ZDC/Energies"))->Fill(13., (zdc.energySectorZNC())[2]);
    registry.get<TH2>(HIST("ZDC/Energies"))->Fill(14., (zdc.energySectorZPA())[0]);
    registry.get<TH2>(HIST("ZDC/Energies"))->Fill(15., (zdc.energySectorZPA())[1]);
    registry.get<TH2>(HIST("ZDC/Energies"))->Fill(16., (zdc.energySectorZPA())[2]);
    registry.get<TH2>(HIST("ZDC/Energies"))->Fill(17., (zdc.energySectorZPA())[3]);
    registry.get<TH2>(HIST("ZDC/Energies"))->Fill(18., (zdc.energySectorZPC())[0]);
    registry.get<TH2>(HIST("ZDC/Energies"))->Fill(19., (zdc.energySectorZPC())[1]);
    registry.get<TH2>(HIST("ZDC/Energies"))->Fill(20., (zdc.energySectorZPC())[2]);
    registry.get<TH2>(HIST("ZDC/Energies"))->Fill(21., (zdc.energySectorZPC())[3]);
  };
  PROCESS_SWITCH(DiffQA, processZDC, "Process ZDC", true);

  // ...............................................................................................................
  void processTest(CCs const& collisions, BCs const& bcs)
  {
    uint64_t bc1, bc2, bc3;
    for (auto col : collisions) {
      bc1 = -1;
      bc2 = -2;
      bc3 = -3;
      if (col.has_foundBC()) {
        auto bc = col.foundBC_as<BCs>();
        bc1 = bc.globalBC();
      }
      if (col.has_bc()) {
        auto bc = col.bc_as<BCs>();
        bc2 = bc.globalBC();
      }
      auto bc = bcs.rawIteratorAt(col.globalIndex());
      bc3 = bc.globalBC();

      if (bc1 != bc2 || bc1 != bc3) {
        LOGF(info, "BC missmatch: %d %d %d", bc1, bc2, bc3);
      }
    }
  };
  PROCESS_SWITCH(DiffQA, processTest, "Process test", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<DiffQA>(cfgc, TaskName{"diffqa"}),
  };
}
