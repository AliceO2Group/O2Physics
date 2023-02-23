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
#include "ReconstructionDataFormats/BCRange.h"
#include "CommonConstants/PhysicsConstants.h"
#include "PWGUD/Core/UDHelpers.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct DiffQA {

  // constants
  static const int nBCpOrbit = 3564;
  static const int ns = 20;    // number of BCs to save (used in processF[V0, T0, DD])
  static const int ncmin = 20; // minimum length of series of empty BCs  (used in processF[V0, T0, DD])

  static constexpr std::string_view hFV0A[ns + 1] = {"fv0A00", "fv0A01", "fv0A02", "fv0A03", "fv0A04", "fv0A05", "fv0A06", "fv0A07", "fv0A08", "fv0A09", "fv0A10", "fv0A11", "fv0A12", "fv0A13", "fv0A14", "fv0A15", "fv0A16", "fv0A17", "fv0A18", "fv0A19", "fv0A20"};
  static constexpr std::string_view hFT0A[ns + 1] = {"ft0A00", "ft0A01", "ft0A02", "ft0A03", "ft0A04", "ft0A05", "ft0A06", "ft0A07", "ft0A08", "ft0A09", "ft0A10", "ft0A11", "ft0A12", "ft0A13", "ft0A14", "ft0A15", "ft0A16", "ft0A17", "ft0A18", "ft0A19", "ft0A20"};
  static constexpr std::string_view hFT0C[ns + 1] = {"ft0C00", "ft0C01", "ft0C02", "ft0C03", "ft0C04", "ft0C05", "ft0C06", "ft0C07", "ft0C08", "ft0C09", "ft0C10", "ft0C11", "ft0C12", "ft0C13", "ft0C14", "ft0C15", "ft0C16", "ft0C17", "ft0C18", "ft0C19", "ft0C20"};
  static constexpr std::string_view hFDDA[ns + 1] = {"fddA00", "fddA01", "fddA02", "fddA03", "fddA04", "fddA05", "fddA06", "fddA07", "fddA08", "fddA09", "fddA10", "fddA11", "fddA12", "fddA13", "fddA14", "fddA15", "fddA16", "fddA17", "fddA18", "fddA19", "fddA20"};
  static constexpr std::string_view hFDDC[ns + 1] = {"fddC00", "fddC01", "fddC02", "fddC03", "fddC04", "fddC05", "fddC06", "fddC07", "fddC08", "fddC09", "fddC10", "fddC11", "fddC12", "fddC13", "fddC14", "fddC15", "fddC16", "fddC17", "fddC18", "fddC19", "fddC20"};

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
      registry.add("Stat", "#Stat", {HistType::kTH1F, {{20, -0.5, 19.5}}});
      registry.add("Tracks", "#Tracks", {HistType::kTH1F, {{50, 0.5, 50.5}}});
      registry.add("vtxTracks", "#vtxTracks", {HistType::kTH1F, {{50, 0.5, 50.5}}});
      registry.add("globalTracks", "#globalTracks", {HistType::kTH1F, {{50, 0.5, 50.5}}});
      registry.add("rejectedTracks", "#rejectedTracks", {HistType::kTH1F, {{17, -0.5, 16.5}}});
      registry.add("tResvsrTOFTracks", "#tResvsrTOFTracks", {HistType::kTH2F, {{1000, 0., 1.E3}, {101, -0.01, 1.01}}});
      registry.add("vtxPosxy", "#vtxPosxy", {HistType::kTH2F, {{200, -2., 2.}, {200, -2., 2.}}});
      registry.add("vtxPosz", "#vtxPosz", {HistType::kTH1F, {{1000, -100., 100.}}});
      registry.add("etapt", "#etapt", {HistType::kTH2F, {{80, -2., 2.}, {100, 0., 5.}}});
      registry.add("dEdxTPC", "#dEdxTPC", {HistType::kTH2F, {{120, -6., 6.}, {1000, 0., 1000.}}});
      registry.add("dEdxTOF", "#dEdxTOF", {HistType::kTH2F, {{100, 0., 5.0}, {1000, 0., 500000.}}});
      registry.add("vtxPosxyDG", "#vtxPosxyDG", {HistType::kTH2F, {{200, -2., 2.}, {200, -2., 2.}}});
      registry.add("vtxPoszDG", "#vtxPoszDG", {HistType::kTH1F, {{1000, -100., 100.}}});
      registry.add("etaptDG", "#etaptDG", {HistType::kTH2F, {{80, -2., 2.}, {100, 0., 5.}}});
      registry.add("dEdxTPCDG", "#dEdxTPCDG", {HistType::kTH2F, {{120, -6., 6.0}, {1000, 0., 1000.}}});
      registry.add("dEdxTOFDG", "#dEdxTOFDG", {HistType::kTH2F, {{100, 0., 5.0}, {1000, 0., 500000.}}});
      registry.add("netChargeDG", "#netChargeDG", {HistType::kTH1F, {{21, -10.5, 10.5}}});
      registry.add("IVMptSysDG", "#IVMptSysDG", {HistType::kTH2F, {{100, 0., 5.}, {350, 0., 3.5}}});
      registry.add("IVMptTrkDG", "#IVMptTrkDG", {HistType::kTH2F, {{100, 0., 5.}, {350, 0., 3.5}}});
    }
    if (context.mOptions.get<bool>("processFewProng")) {
      registry.add("fpStat", "#fpStat", {HistType::kTH1F, {{2, 0.5, 2.5}}});
      registry.add("allPVC", "#allPVC", {HistType::kTH1F, {{100, 0.5, 100.5}}});
      registry.add("fpPVC", "#fpPVC", {HistType::kTH1F, {{100, 0.5, 100.5}}});
    }
    if (context.mOptions.get<bool>("processCleanFIT1")) {
      registry.add("cleanFIT1", "#cleanFIT1", {HistType::kTH2F, {{20, -0.5, 19.5}, {2, -0.5, 1.5}}});
      registry.add("cF1FV0Aamp", "#cF1FV0Aamp", {HistType::kTH2F, {{20, -0.5, 19.5}, {1000, -0.5, 999.5}}});
      registry.add("cF1FT0Aamp", "#cF1FT0Aamp", {HistType::kTH2F, {{20, -0.5, 19.5}, {1000, -0.5, 999.5}}});
      registry.add("cF1FT0Camp", "#cF1FT0Camp", {HistType::kTH2F, {{20, -0.5, 19.5}, {1000, -0.5, 999.5}}});
      registry.add("cF1FDDAamp", "#cF1FDDAamp", {HistType::kTH2F, {{20, -0.5, 19.5}, {1000, -0.5, 999.5}}});
      registry.add("cF1FDDCamp", "#cF1FDDCamp", {HistType::kTH2F, {{20, -0.5, 19.5}, {1000, -0.5, 999.5}}});
    }
    if (context.mOptions.get<bool>("processCleanFIT2")) {
      registry.add("cleanFIT2", "#cleanFIT2", {HistType::kTH2F, {{20, -0.5, 19.5}, {2, -0.5, 1.5}}});
      registry.add("cF2FV0Aamp", "#cF2FV0Aamp", {HistType::kTH2F, {{20, -0.5, 19.5}, {1000, -0.5, 999.5}}});
      registry.add("cF2FT0Aamp", "#cF2FT0Aamp", {HistType::kTH2F, {{20, -0.5, 19.5}, {1000, -0.5, 999.5}}});
      registry.add("cF2FT0Camp", "#cF2FT0Camp", {HistType::kTH2F, {{20, -0.5, 19.5}, {1000, -0.5, 999.5}}});
      registry.add("cF2FDDAamp", "#cF2FDDAamp", {HistType::kTH2F, {{20, -0.5, 19.5}, {1000, -0.5, 999.5}}});
      registry.add("cF2FDDCamp", "#cF2FDDCamp", {HistType::kTH2F, {{20, -0.5, 19.5}, {1000, -0.5, 999.5}}});
    }
    if (context.mOptions.get<bool>("processFV0")) {
      registry.add("FV0A", "#FV0A", {HistType::kTH2F, {{48, -0.5, 47.5}, {2000, 0., 2000.}}});
      registry.add("FV0BCNUM", "#FV0BCNUM", {HistType::kTH1F, {{3564, -0.5, 3563.5}}});
    }
    if (context.mOptions.get<bool>("processFT0")) {
      registry.add("FT0A", "#FT0A", {HistType::kTH2F, {{96, -0.5, 95.5}, {400, 0., 400.}}});
      registry.add("FT0C", "#FT0C", {HistType::kTH2F, {{112, -0.5, 111.5}, {400, 0., 400.}}});
      registry.add("FT0ABCNUM", "#FT0ABCNUM", {HistType::kTH1F, {{3564, -0.5, 3563.5}}});
      registry.add("FT0CBCNUM", "#FT0CBCNUM", {HistType::kTH1F, {{3564, -0.5, 3563.5}}});
      registry.add("dFT0BCNUM", "#dFT0BCNUM", {HistType::kTH1F, {{3564, -0.5, 3563.5}}});

      // add amplitude histograms
      for (auto n{0}; n <= ns; n++) {
        registry.add(hFT0A[n].data(), hFT0A[n].data(), {HistType::kTH1F, {{1000, 0., 1000.}}});
        registry.add(hFT0C[n].data(), hFT0C[n].data(), {HistType::kTH1F, {{1000, 0., 1000.}}});
      }
    }
    if (context.mOptions.get<bool>("processFDD")) {
      registry.add("FDDA", "#FDDA", {HistType::kTH2F, {{8, -0.5, 7.5}, {100, 0., 100.}}});
      registry.add("FDDC", "#FDDC", {HistType::kTH2F, {{8, -0.5, 7.5}, {100, 0., 100.}}});
      registry.add("FDDBCNUM", "#FDDBCNUM", {HistType::kTH1F, {{3564, -0.5, 3563.5}}});
    }
    if (context.mOptions.get<bool>("processZDC")) {
      registry.add("ZdcEnergies", "#ZdcEnergies", {HistType::kTH2F, {{22, -0.5, 21.5}, {100, 0., 1000.}}});
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
    LOGF(debug, "<DiffQA. Collision %d", collision.globalIndex());
    LOGF(debug, "<DiffQA> Start %i", abcrs.size());

    bool isDGcandidate = true;
    registry.get<TH1>(HIST("Stat"))->Fill(0., isDGcandidate * 1.);

    // update collision histograms
    // vertex position
    registry.get<TH2>(HIST("vtxPosxy"))->Fill(collision.posX(), collision.posY());
    registry.get<TH1>(HIST("vtxPosz"))->Fill(collision.posZ());
    // tracks
    registry.get<TH1>(HIST("Tracks"))->Fill(tracks.size());
    // vertex tracks
    registry.get<TH1>(HIST("vtxTracks"))->Fill(collision.numContrib());
    // global tracks
    Partition<TCs> goodTracks = requireGlobalTrackInFilter();
    goodTracks.bindTable(tracks);
    registry.get<TH1>(HIST("globalTracks"))->Fill(goodTracks.size());

    // number of vertex tracks with TOF hit
    float rgtrwTOF = 0.;
    for (auto const& track : tracks) {
      // update eta vs pt histogram
      registry.get<TH2>(HIST("etapt"))->Fill(track.eta(), track.pt());
      // update dEdx histograms
      registry.get<TH2>(HIST("dEdxTPC"))->Fill(track.p() * track.sign(), track.tpcSignal());
      if (track.tpcSignal() > maxdEdxTPC) {
        maxdEdxTPC = track.tpcSignal();
        LOGF(debug, "<DiffQA> New maxdEdx TPC %f", maxdEdxTPC);
      }

      // TOF hit?
      if (track.hasTOF()) {
        registry.get<TH2>(HIST("dEdxTOF"))->Fill(track.pt(), track.tofSignal());
        if (track.tofSignal() > maxdEdxTOF) {
          maxdEdxTOF = track.tofSignal();
          LOGF(debug, "<DiffQA> New maxdEdx tOF %f", maxdEdxTOF);
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
    LOGF(debug, "<DiffQA> Vertex tracks with TOF: %f [1]", rgtrwTOF);
    registry.get<TH2>(HIST("tResvsrTOFTracks"))->Fill(collision.collisionTimeRes(), rgtrwTOF);

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
        if (!udhelpers::cleanFIT(bc, diffCuts.FITAmpLimits())) {
          isDGcandidate = false;
          break;
        }
      }
    } else {
      if (!udhelpers::cleanFITCollision(collision, diffCuts.FITAmpLimits())) {
        isDGcandidate = false;
      }
    }
    registry.get<TH1>(HIST("Stat"))->Fill(1., isDGcandidate * 1.);

    // no Zdc signal in bcSlice
    std::vector<float> lims(10, 0.);
    for (auto const& bc : bcSlice) {
      if (!udhelpers::cleanZDC(bc, zdcs, lims)) {
        isDGcandidate = false;
        break;
      }
    }
    registry.get<TH1>(HIST("Stat"))->Fill(2., isDGcandidate * 1.);

    // no Calo signal in bcSlice
    for (auto const& bc : bcSlice) {
      if (!udhelpers::cleanCalo(bc, calos, lims)) {
        isDGcandidate = false;
        break;
      }
    }
    registry.get<TH1>(HIST("Stat"))->Fill(3., isDGcandidate * 1.);

    // no V0s
    isDGcandidate &= (v0s.size() == 0);
    registry.get<TH1>(HIST("Stat"))->Fill(4., isDGcandidate * 1.);

    // no Cascades
    isDGcandidate &= (cascades.size() == 0);
    registry.get<TH1>(HIST("Stat"))->Fill(5., isDGcandidate * 1.);

    // number of forward tracks = 0
    isDGcandidate &= (fwdtracks.size() == 0);
    registry.get<TH1>(HIST("Stat"))->Fill(6., isDGcandidate * 1.);

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
    registry.get<TH1>(HIST("Stat"))->Fill(7., globalAndVtx * 1.);
    registry.get<TH1>(HIST("Stat"))->Fill(8., vtxAndGlobal * 1.);
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
    registry.get<TH1>(HIST("Stat"))->Fill(9., noAmbTracks * 1.);

    // check a given bc for possible ambiguous FwdTracks
    auto noAmbFwdTracks = isDGcandidate;
    for (auto const& bc : bcSlice) {
      if (afbcrs.isInRange(bc.globalIndex())) {
        noAmbFwdTracks = false;
        break;
      }
    }
    registry.get<TH1>(HIST("Stat"))->Fill(10., noAmbFwdTracks * 1.);

    // fraction of PV tracks with TOF hit
    isDGcandidate &= (rgtrwTOF >= diffCuts.minRgtrwTOF());
    registry.get<TH1>(HIST("Stat"))->Fill(11., isDGcandidate * 1.);

    // number of vertex tracks <= n
    isDGcandidate &= (collision.numContrib() >= diffCuts.minNTracks());
    registry.get<TH1>(HIST("Stat"))->Fill(12., isDGcandidate * 1.);
    isDGcandidate &= (collision.numContrib() <= diffCuts.maxNTracks());
    registry.get<TH1>(HIST("Stat"))->Fill(13., isDGcandidate * 1.);

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
        // update histogram rejectedTracks
        if (!track.isPVContributor()) {
          registry.get<TH1>(HIST("rejectedTracks"))->Fill(0., 1.);
          registry.get<TH1>(HIST("rejectedTracks"))->Fill(1., track.isGlobalTrackSDD() * 1.);
          registry.get<TH1>(HIST("rejectedTracks"))->Fill(2., track.passedTrackType() * 1.);
          registry.get<TH1>(HIST("rejectedTracks"))->Fill(3., track.passedPtRange() * 1.);
          registry.get<TH1>(HIST("rejectedTracks"))->Fill(4., track.passedEtaRange() * 1.);
          registry.get<TH1>(HIST("rejectedTracks"))->Fill(5., track.passedTPCNCls() * 1.);
          registry.get<TH1>(HIST("rejectedTracks"))->Fill(6., track.passedTPCCrossedRows() * 1.);
          registry.get<TH1>(HIST("rejectedTracks"))->Fill(7., track.passedTPCCrossedRowsOverNCls() * 1.);
          registry.get<TH1>(HIST("rejectedTracks"))->Fill(8., track.passedTPCChi2NDF() * 1.);
          registry.get<TH1>(HIST("rejectedTracks"))->Fill(9., track.passedTPCRefit() * 1.);
          registry.get<TH1>(HIST("rejectedTracks"))->Fill(10., track.passedITSNCls() * 1.);
          registry.get<TH1>(HIST("rejectedTracks"))->Fill(11., track.passedITSChi2NDF() * 1.);
          registry.get<TH1>(HIST("rejectedTracks"))->Fill(12., track.passedITSRefit() * 1.);
          registry.get<TH1>(HIST("rejectedTracks"))->Fill(13., track.passedITSHits() * 1.);
          registry.get<TH1>(HIST("rejectedTracks"))->Fill(14., track.passedGoldenChi2() * 1.);
          registry.get<TH1>(HIST("rejectedTracks"))->Fill(15., track.passedDCAxy() * 1.);
          registry.get<TH1>(HIST("rejectedTracks"))->Fill(16., track.passedDCAz() * 1.);
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
    registry.get<TH1>(HIST("Stat"))->Fill(14., isDGcandidate * 1.);
    isDGcandidate &= goodetas;
    registry.get<TH1>(HIST("Stat"))->Fill(15., isDGcandidate * 1.);
    auto netChargeValues = diffCuts.netCharges();
    if (std::find(netChargeValues.begin(), netChargeValues.end(), netCharge) == netChargeValues.end()) {
      goodnchs = false;
    }
    isDGcandidate &= goodnchs;
    registry.get<TH1>(HIST("Stat"))->Fill(16., isDGcandidate * 1.);
    registry.get<TH1>(HIST("Stat"))->Fill(17., isDGcandidate * 1.);
    isDGcandidate &= (ivm.M() >= diffCuts.minIVM());
    registry.get<TH1>(HIST("Stat"))->Fill(18., isDGcandidate * 1.);
    isDGcandidate &= (ivm.M() <= diffCuts.maxIVM());
    registry.get<TH1>(HIST("Stat"))->Fill(19., isDGcandidate * 1.);

    // update some DG histograms
    if (isDGcandidate) {
      // vertex position of DG events
      registry.get<TH2>(HIST("vtxPosxyDG"))->Fill(collision.posX(), collision.posY());
      registry.get<TH1>(HIST("vtxPoszDG"))->Fill(collision.posZ());
      registry.get<TH1>(HIST("netChargeDG"))->Fill(netCharge);
      registry.get<TH2>(HIST("IVMptSysDG"))->Fill(ivm.M(), ivm.Perp());

      // fill dEdx of DG event tracks
      for (auto const& track : tracks) {
        if (track.isPVContributor()) {
          LOGF(debug, "dEdx TPC %f TOF %i %f", track.tpcSignal(), track.hasTOF(), track.hasTOF() ? track.tofSignal() : 0.);
          registry.get<TH2>(HIST("dEdxTPCDG"))->Fill(track.p() * track.sign(), track.tpcSignal());
          registry.get<TH2>(HIST("etaptDG"))->Fill(track.eta(), track.pt());
          registry.get<TH2>(HIST("IVMptTrkDG"))->Fill(ivm.M(), track.pt());
          if (track.hasTOF()) {
            registry.get<TH2>(HIST("dEdxTOFDG"))->Fill(track.pt(), track.tofSignal());
          }
        }
      }
    }
  }
  PROCESS_SWITCH(DiffQA, processMain, "Process Main", true);

  // ...............................................................................................................
  // Distribution of number of PV contributors for all collisions and those with empty FT0
  void processFewProng(CC const& collision, BCs const& bct0s,
                       aod::FT0s const& ft0s, aod::FV0As const& fv0as, aod::FDDs const& fdds)
  {
    // count collisions
    registry.get<TH1>(HIST("fpStat"))->Fill(1., 1.);
    registry.get<TH1>(HIST("allPVC"))->Fill(collision.numContrib(), 1.);

    // check FT0 to be empty
    auto bc = collision.foundBC_as<BCs>();
    if (udhelpers::cleanFT0(bc, 0., 0.)) {
      // only collisions with empty FT0 arrive here
      registry.get<TH1>(HIST("fpStat"))->Fill(2., 1.);

      // update #PV contributors in collisions with empty FT0
      registry.get<TH1>(HIST("fpPVC"))->Fill(collision.numContrib(), 1.);
    }
  }
  PROCESS_SWITCH(DiffQA, processFewProng, "Process FewProng", true);

  // ...............................................................................................................
  // Fraction of collisions with empty FIT as function of NDtcoll
  void processCleanFIT1(CC const& collision, BCs const& bct0s,
                        aod::FT0s const& ft0s, aod::FV0As const& fv0as, aod::FDDs const& fdds)
  {
    LOGF(debug, "<CleanFit. Collision %d", collision.globalIndex());

    // test influence of BCrange width using a series of NDtcoll
    float ampFV0A, ampFT0A, ampFT0C, ampFDDA, ampFDDC;
    bool isDGcandidate = true;
    for (int NDtcoll = 0; NDtcoll < 20; NDtcoll++) {
      auto bcSlice = udhelpers::compatibleBCs(collision, NDtcoll, bct0s, 0);
      ampFV0A = ampFT0A = ampFT0C = ampFDDA = ampFDDC = 0.;
      isDGcandidate = true;
      for (auto const& bc : bcSlice) {
        isDGcandidate &= udhelpers::cleanFIT(bc, diffCuts.FITAmpLimits());

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
      registry.get<TH2>(HIST("cleanFIT1"))->Fill(NDtcoll, isDGcandidate * 1.);
      if (isDGcandidate) {
        registry.get<TH2>(HIST("cF1FV0Aamp"))->Fill(NDtcoll, ampFV0A);
        registry.get<TH2>(HIST("cF1FT0Aamp"))->Fill(NDtcoll, ampFT0A);
        registry.get<TH2>(HIST("cF1FT0Camp"))->Fill(NDtcoll, ampFT0C);
        registry.get<TH2>(HIST("cF1FDDAamp"))->Fill(NDtcoll, ampFDDA);
        registry.get<TH2>(HIST("cF1FDDCamp"))->Fill(NDtcoll, ampFDDC);
      }
    }
  }
  PROCESS_SWITCH(DiffQA, processCleanFIT1, "Process CleanFitTest1", true);

  // ...............................................................................................................
  void processCleanFIT2(CC const& collision, BCs const& bct0s,
                        aod::FT0s const& ft0s, aod::FV0As const& fv0as, aod::FDDs const& fdds)
  {
    LOGF(debug, "<CleanFit. Collision %d", collision.globalIndex());

    // test influence of BCrange width using a series of nMinBC
    float ampFV0A, ampFT0A, ampFT0C, ampFDDA, ampFDDC;
    bool isDGcandidate = true;
    for (int nMinBC = 0; nMinBC < 20; nMinBC++) {
      auto bcSlice = udhelpers::compatibleBCs(collision, 0, bct0s, nMinBC);
      ampFV0A = ampFT0A = ampFT0C = ampFDDA = ampFDDC = 0.;
      isDGcandidate = true;
      for (auto const& bc : bcSlice) {
        isDGcandidate &= udhelpers::cleanFIT(bc, diffCuts.FITAmpLimits());

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
      registry.get<TH2>(HIST("cleanFIT2"))->Fill(nMinBC, isDGcandidate * 1.);
      if (isDGcandidate) {
        registry.get<TH2>(HIST("cF2FV0Aamp"))->Fill(nMinBC, ampFV0A);
        registry.get<TH2>(HIST("cF2FT0Aamp"))->Fill(nMinBC, ampFT0A);
        registry.get<TH2>(HIST("cF2FT0Camp"))->Fill(nMinBC, ampFT0C);
        registry.get<TH2>(HIST("cF2FDDAamp"))->Fill(nMinBC, ampFDDA);
        registry.get<TH2>(HIST("cF2FDDCamp"))->Fill(nMinBC, ampFDDC);
      }
    }
  }
  PROCESS_SWITCH(DiffQA, processCleanFIT2, "Process CleanFitTest2", true);

  // ...............................................................................................................
  void processFV0(aod::FV0As const& fv0s, BCs const&)
  {
    LOGF(info, "<FV0Signals> %d", fv0s.size());
    if (fv0s.size() <= 0) {
      return;
    }

    int64_t lastBCwFV0 = fv0s.begin().bc_as<BCs>().globalBC();
    auto lastOrbit = lastBCwFV0 / nBCpOrbit;
    for (auto fv0 : fv0s) {

      // side A
      for (size_t ind = 0; ind < fv0.channel().size(); ind++) {
        registry.get<TH2>(HIST("FV0A"))->Fill((fv0.channel())[ind], (fv0.amplitude())[ind]);
      }

      // sequence of BCs
      auto bc = fv0.bc_as<BCs>();
      int64_t aBC = bc.globalBC();
      auto aOrbit = aBC / nBCpOrbit;
      auto ampA = udhelpers::FV0AmplitudeA(bc.foundFV0());
      if (ampA > 0.) {
        // require both BCs to be in same orbit
        if (aOrbit == lastOrbit)
          registry.get<TH1>(HIST("FV0BCNUM"))->Fill((bc.globalBC() - lastBCwFV0));
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
    float minAmpA = 15., fAmpA = 0.;
    float minAmpC = 15., fAmpC = 0.;

    int64_t lastBCwFT0 = ft0s.begin().bc_as<BCs>().globalBC();
    auto lastOrbit = lastBCwFT0 / nBCpOrbit;
    for (auto ft0 : ft0s) {

      // side A
      for (size_t ind = 0; ind < ft0.channelA().size(); ind++) {
        registry.get<TH2>(HIST("FT0A"))->Fill((ft0.channelA())[ind], (ft0.amplitudeA())[ind]);
      }

      // side C
      for (size_t ind = 0; ind < ft0.channelC().size(); ind++) {
        registry.get<TH2>(HIST("FT0C"))->Fill((ft0.channelC())[ind], (ft0.amplitudeC())[ind]);
      }

      // sequence of BCs
      auto bc = ft0.bc_as<BCs>();
      aBC = bc.globalBC();
      auto gBC = aBC % nBCpOrbit;
      auto ampA = udhelpers::FT0AmplitudeA(bc.foundFT0());
      auto ampC = udhelpers::FT0AmplitudeC(bc.foundFT0());

      // update FT0ABCNUM
      if (ampA > 0.) {
        registry.get<TH1>(HIST("FT0ABCNUM"))->Fill(gBC, 1.);
      }

      // update FT0CBCNUM
      if (ampC > 0.) {
        registry.get<TH1>(HIST("FT0CBCNUM"))->Fill(gBC, 1.);
      }

      // update dFT0BCNUM
      // require both BCs to be in same orbit
      auto aOrbit = aBC / nBCpOrbit;
      LOGF(debug, "lastOrbit %d aOrbit %d", lastOrbit, aOrbit);
      if (ampA > 0. || ampC > 0.) {
        if (aOrbit == lastOrbit)
          registry.get<TH1>(HIST("dFT0BCNUM"))->Fill((aBC - lastBCwFT0));
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
      if (dBC <= ns) {
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
    auto lastOrbit = lastBCwFDD / nBCpOrbit;
    for (auto fdd : fdds) {

      // side A
      for (auto ind = 0; ind < 8; ind++) {
        registry.get<TH2>(HIST("FDDA"))->Fill(ind, (fdd.chargeA())[ind]);
      }

      // side C
      for (auto ind = 0; ind < 8; ind++) {
        registry.get<TH2>(HIST("FDDC"))->Fill(ind, (fdd.chargeC())[ind]);
      }

      // sequence of BCs
      auto bc = fdd.bc_as<BCs>();
      auto aBC = bc.globalBC();
      int64_t aOrbit = aBC / nBCpOrbit;
      auto ampA = udhelpers::FDDAmplitudeA(bc.foundFDD());
      auto ampC = udhelpers::FDDAmplitudeC(bc.foundFDD());
      if (ampA > 0. || ampC > 0.) {
        // require both BCs to be in same orbit
        if (aOrbit == lastOrbit)
          registry.get<TH1>(HIST("FDDBCNUM"))->Fill((bc.globalBC() - lastBCwFDD));
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
  };
  PROCESS_SWITCH(DiffQA, processZDC, "Process ZDC", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<DiffQA>(cfgc, TaskName{"diffqa"}),
  };
}
