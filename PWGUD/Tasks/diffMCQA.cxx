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
///
///     options:
///           DiffCuts.mNDtcoll(4)
///           DiffCuts.mMinNBCs(7)
///           DiffCuts.mMinNTracks(0)
///           DiffCuts.mMaxNTracks(10000)
///           DiffCuts.mMinNetCharge(0)
///           DiffCuts.mMaxNetCharge(0)
///           DiffCuts.mPidHypo(211)
///           DiffCuts.mMinPosz(-1000.)
///           DiffCuts.mMaxPosz(1000.)
///           DiffCuts.mMinPt(0.)
///           DiffCuts.mMaxPt(1000.)
///           DiffCuts.mMinEta(-1.)
///           DiffCuts.mMaxEta(1.)
///           DiffCuts.mMinIVM(0.)
///           DiffCuts.mMaxIVM(1000.)
///           DiffCuts.mMaxnSigmaTPC(1000.)
///           DiffCuts.mMaxnSigmaTOF(1000.)
///           DiffCutsX.mFITAmpLimits({0., 0., 0., 0., 0.})
///
///     usage: copts="--configuration json://DiffQAConfig.json -b"
///
///           o2-analysis-timestamp $copts |
///           o2-analysis-track-propagation $copts |
///           o2-analysis-event-selection $copts |
///           o2-analysis-ft0-corrected-table $copts |
///           o2-analysis-trackextension $copts |
///           o2-analysis-trackselection $copts |
///           o2-analysis-ud-diff-mcqa $copts > diffQA.log
///
/// \author Paul Buehler, paul.buehler@oeaw.ac.at
/// \since  01.10.2021

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "TLorentzVector.h"
#include "ReconstructionDataFormats/BCRange.h"
#include "CommonConstants/PhysicsConstants.h"
#include "EventFiltering/PWGUD/DGHelpers.h"
#include "PWGUD/Core/DGMCHelpers.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct DiffMCQA {

  float maxdEdxTPC;
  float maxdEdxTOF;

  // get a DGCutparHolder
  DGCutparHolder diffCuts = DGCutparHolder();
  MutableConfigurable<DGCutparHolder> DGCuts{"DGCuts", {}, "DG event cuts"};
  Configurable<bool> withAmbTrackAnalysis{"ambiguousTracks", false, "with ambiguous tracks analysis"};
  Configurable<bool> withAmbFwdTrackAnalysis{"ambiguousFwdTracks", false, "with ambiguous forward tracks analysis"};

  // structures to hold information about the possible BCs the ambiguous tracks/FwdTracks belong to
  o2::dataformats::bcRanges abcrs = o2::dataformats::bcRanges("ambiguous_tracks");
  o2::dataformats::bcRanges afbcrs = o2::dataformats::bcRanges("ambiguous_fwdtracks");

  // define histograms
  // Stat:
  //  bin  0: all collisions
  //  bin  1: clean FIT
  //  bin  2: clean ZDC
  //  bin  3: clean Calo
  //  bin  4: no V0s
  //  bin  5: no Cascades
  //  bin  6: no FWD tracks
  //  bin  7: no global tracks which are no vtx tracks
  //  bin  8: no vtx tracks which are no global tracks
  //  bin  9: possible ambiguous tracks
  //  bin 10: possible ambiguous FwdTracks
  //  bin 11: fraction of PV tracks with TOF hit > minRgtrwTOF
  //  bin 12: number of tracks >= minimum number
  //  bin 13: number of tracks <= maximum number
  //  bin 14: minimum pt <= pt of vtx tracks <= maximum pt
  //  bin 15: minimum eta <= eta of vtx tracks <= maximum eta
  //  bin 16: net charge >= minimum net charge
  //  bin 17: net charge <= maximum net charge
  //  bin 18: IVM >= minimum IVM
  //  bin 19: IVM <= maximum IVM
  //
  // 3 diverent versions of histograms:
  //  Diff1: Pythia MBR
  //  Diff2: GRANIITTI
  //       : Rest (Pythia MB)
  //
  HistogramRegistry registry{
    "registry",
    {
      // non diffractive events
      {"Stat", "#Stat", {HistType::kTH1F, {{20, -0.5, 19.5}}}},
      {"cleanFIT", "#cleanFIT", {HistType::kTH2F, {{10, -0.5, 9.5}, {2, -0.5, 1.5}}}},
      {"Tracks", "#Tracks", {HistType::kTH1F, {{50, 0.5, 50.5}}}},
      {"vtxTracks", "#vtxTracks", {HistType::kTH1F, {{50, 0.5, 50.5}}}},
      {"globalTracks", "#globalTracks", {HistType::kTH1F, {{50, 0.5, 50.5}}}},
      {"rejectedTracks", "#rejectedTracks", {HistType::kTH1F, {{17, -0.5, 16.5}}}},
      {"tResvsrTOFTracks", "#tResvsrTOFTracks", {HistType::kTH2F, {{1000, 0., 1.E3}, {101, -0.01, 1.01}}}},
      {"vtxPosxy", "#vtxPosxy", {HistType::kTH2F, {{200, -2., 2.}, {200, -2., 2.}}}},
      {"vtxPosz", "#vtxPosz", {HistType::kTH1F, {{1000, -100., 100.}}}},
      {"etapt", "#etapt", {HistType::kTH2F, {{80, -2., 2.}, {100, 0., 5.}}}},
      {"dEdxTPC", "#dEdxTPC", {HistType::kTH2F, {{100, 0., 5.0}, {3000, 0., 30000.}}}},
      {"dEdxTOF", "#dEdxTOF", {HistType::kTH2F, {{100, 0., 5.0}, {1000, 0., 500000.}}}},
      {"vtxPosxyDG", "#vtxPosxyDG", {HistType::kTH2F, {{200, -2., 2.}, {200, -2., 2.}}}},
      {"vtxPoszDG", "#vtxPoszDG", {HistType::kTH1F, {{1000, -100., 100.}}}},
      {"etaptDG", "#etaptDG", {HistType::kTH2F, {{80, -2., 2.}, {100, 0., 5.}}}},
      {"dEdxTPCDG", "#dEdxTPCDG", {HistType::kTH2F, {{100, 0., 5.0}, {3000, 0., 30000.}}}},
      {"dEdxTOFDG", "#dEdxTOFDG", {HistType::kTH2F, {{100, 0., 5.0}, {1000, 0., 500000.}}}},
      {"netChargeDG", "#netChargeDG", {HistType::kTH1F, {{21, -10.5, 10.5}}}},
      {"IVMptSysDG", "#IVMptSysDG", {HistType::kTH2F, {{100, 0., 5.}, {350, 0., 3.5}}}},
      {"IVMptTrkDG", "#IVMptTrkDG", {HistType::kTH2F, {{100, 0., 5.}, {350, 0., 3.5}}}},
      // PYTHIA8 diffractive events
      {"StatDiff1", "#StatDiff1", {HistType::kTH1F, {{20, -0.5, 19.5}}}},
      {"cleanFITDiff1", "#cleanFITDiff1", {HistType::kTH2F, {{10, -0.5, 9.5}, {2, -0.5, 1.5}}}},
      {"TracksDiff1", "#TracksDiff1", {HistType::kTH1F, {{50, 0.5, 50.5}}}},
      {"vtxTracksDiff1", "#vtxTracksDiff1", {HistType::kTH1F, {{50, 0.5, 50.5}}}},
      {"globalTracksDiff1", "#globalTracksDiff1", {HistType::kTH1F, {{50, 0.5, 50.5}}}},
      {"rejectedTracksDiff1", "#rejectedTracksDiff1", {HistType::kTH1F, {{17, -0.5, 16.5}}}},
      {"tResvsrTOFTracksDiff1", "#tResvsrTOFTracksDiff1", {HistType::kTH2F, {{1000, 0., 1.E3}, {101, -0.01, 1.01}}}},
      {"vtxPosxyDiff1", "#vtxPosxyDiff1", {HistType::kTH2F, {{200, -2., 2.}, {200, -2., 2.}}}},
      {"vtxPoszDiff1", "#vtxPoszDiff1", {HistType::kTH1F, {{1000, -100., 100.}}}},
      {"etaptDiff1", "#etaptDiff1", {HistType::kTH2F, {{80, -2., 2.}, {100, 0., 5.}}}},
      {"dEdxTPCDiff1", "#dEdxTPCDiff1", {HistType::kTH2F, {{100, 0., 5.0}, {3000, 0., 30000.}}}},
      {"dEdxTOFDiff1", "#dEdxTOFDiff1", {HistType::kTH2F, {{100, 0., 5.0}, {1000, 0., 500000.}}}},
      {"vtxPosxyDGDiff1", "#vtxPosxyDGDiff1", {HistType::kTH2F, {{200, -2., 2.}, {200, -2., 2.}}}},
      {"vtxPoszDGDiff1", "#vtxPoszDGDiff1", {HistType::kTH1F, {{1000, -100., 100.}}}},
      {"etaptDGDiff1", "#etaptDGDiff1", {HistType::kTH2F, {{80, -2., 2.}, {100, 0., 5.}}}},
      {"dEdxTPCDGDiff1", "#dEdxTPCDGDiff1", {HistType::kTH2F, {{100, 0., 5.0}, {3000, 0., 30000.}}}},
      {"dEdxTOFDGDiff1", "#dEdxTOFDGDiff1", {HistType::kTH2F, {{100, 0., 5.0}, {1000, 0., 500000.}}}},
      {"netChargeDGDiff1", "#netChargeDGDiff1", {HistType::kTH1F, {{21, -10.5, 10.5}}}},
      {"IVMptSysDGDiff1", "#IVMptSysDGDiff1", {HistType::kTH2F, {{100, 0., 5.}, {350, 0., 3.5}}}},
      {"IVMptTrkDGDiff1", "#IVMptTrkDGDiff1", {HistType::kTH2F, {{100, 0., 5.}, {350, 0., 3.5}}}},
      // GRANIITTI diffractive events
      {"StatDiff2", "#StatDiff2", {HistType::kTH1F, {{20, -0.5, 19.5}}}},
      {"cleanFITDiff2", "#cleanFITDiff2", {HistType::kTH2F, {{10, -0.5, 9.5}, {2, -0.5, 1.5}}}},
      {"TracksDiff2", "#TracksDiff2", {HistType::kTH1F, {{50, 0.5, 50.5}}}},
      {"vtxTracksDiff2", "#vtxTracksDiff2", {HistType::kTH1F, {{50, 0.5, 50.5}}}},
      {"globalTracksDiff2", "#globalTracksDiff2", {HistType::kTH1F, {{50, 0.5, 50.5}}}},
      {"rejectedTracksDiff2", "#rejectedTracksDiff2", {HistType::kTH1F, {{17, -0.5, 16.5}}}},
      {"tResvsrTOFTracksDiff2", "#tResvsrTOFTracksDiff2", {HistType::kTH2F, {{1000, 0., 1.E3}, {101, -0.01, 1.01}}}},
      {"vtxPosxyDiff2", "#vtxPosxyDiff2", {HistType::kTH2F, {{200, -2., 2.}, {200, -2., 2.}}}},
      {"vtxPoszDiff2", "#vtxPoszDiff2", {HistType::kTH1F, {{1000, -100., 100.}}}},
      {"etaptDiff2", "#etaptDiff2", {HistType::kTH2F, {{80, -2., 2.}, {100, 0., 5.}}}},
      {"dEdxTPCDiff2", "#dEdxTPCDiff2", {HistType::kTH2F, {{100, 0., 5.0}, {3000, 0., 30000.}}}},
      {"dEdxTOFDiff2", "#dEdxTOFDiff2", {HistType::kTH2F, {{100, 0., 5.0}, {1000, 0., 500000.}}}},
      {"vtxPosxyDGDiff2", "#vtxPosxyDGDiff2", {HistType::kTH2F, {{200, -2., 2.}, {200, -2., 2.}}}},
      {"vtxPoszDGDiff2", "#vtxPoszDGDiff2", {HistType::kTH1F, {{1000, -100., 100.}}}},
      {"etaptDGDiff2", "#etaptDGDiff2", {HistType::kTH2F, {{80, -2., 2.}, {100, 0., 5.}}}},
      {"dEdxTPCDGDiff2", "#dEdxTPCDGDiff2", {HistType::kTH2F, {{100, 0., 5.0}, {3000, 0., 30000.}}}},
      {"dEdxTOFDGDiff2", "#dEdxTOFDGDiff2", {HistType::kTH2F, {{100, 0., 5.0}, {1000, 0., 500000.}}}},
      {"netChargeDGDiff2", "#netChargeDGDiff2", {HistType::kTH1F, {{21, -10.5, 10.5}}}},
      {"IVMptSysDGDiff2", "#IVMptSysDGDiff2", {HistType::kTH2F, {{100, 0., 5.}, {350, 0., 3.5}}}},
      {"IVMptTrkDGDiff2", "#IVMptTrkDGDiff2", {HistType::kTH2F, {{100, 0., 5.}, {350, 0., 3.5}}}},
    }};

  using CCs = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;
  using CC = CCs::iterator;
  using BCs = soa::Join<aod::BCs, aod::BcSels, aod::Run3MatchedToBCSparse>;
  using TCs = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TOFSignal>;
  using FWs = soa::Join<aod::FwdTracks, aod::McFwdTrackLabels>;
  using ATs = aod::AmbiguousTracks;
  using AFTs = aod::AmbiguousFwdTracks;

  void init(InitContext&)
  {
    maxdEdxTPC = 0.;
    maxdEdxTOF = 0.;
    diffCuts = (DGCutparHolder)DGCuts;
  }

  void run(ProcessingContext& pc)
  {
    // get ambiguous tracks table
    auto t1 = pc.inputs().get<TableConsumer>("BCs")->asArrowTable();
    auto t2 = pc.inputs().get<TableConsumer>("BcSels")->asArrowTable();
    auto t3 = pc.inputs().get<TableConsumer>("Run3MatchedToBCSparse")->asArrowTable();
    auto bcs = BCs({t1, t2, t3});

    if (withAmbFwdTrackAnalysis) {
      auto t4 = pc.inputs().get<TableConsumer>("AmbiguousTracks")->asArrowTable();
      auto ambtracks = ATs({t4});
      ambtracks.bindExternalIndices(&bcs);

      // make sorted list of BC ranges which are associated with an ambiguous track.
      // This is used to efficiently check whether a given BC is contained in one of these ranges
      abcrs.reset();
      LOGF(info, "<DiffMCQA> size of ambiguous tracks table %i", ambtracks.size());
      for (auto ambtrack : ambtracks) {
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
      for (auto ambfwdtrack : ambfwdtracks) {
        auto bcfirst = ambfwdtrack.bc().rawIteratorAt(0);
        auto bclast = ambfwdtrack.bc().rawIteratorAt(ambfwdtrack.bc().size() - 1);
        afbcrs.add(bcfirst.globalIndex(), bclast.globalIndex());
      }
      afbcrs.merge();
    }
    LOGF(info, "<DiffMCQA> Size of abcrs %i and afbcrs %i", abcrs.size(), afbcrs.size());
  }

  Preslice<aod::McParticles> perMcCollision = aod::mcparticle::mcCollisionId;

  void process(CC const& collision, BCs const& bct0s,
               TCs& tracks, FWs& fwdtracks, ATs& ambtracks, AFTs& ambfwdtracks,
               aod::FT0s& ft0s, aod::FV0As& fv0as, aod::FDDs& fdds,
               aod::Zdcs& zdcs, aod::Calos& calos,
               aod::V0s& v0s, aod::Cascades& cascades,
               aod::McCollisions& McCols, aod::McParticles const& McParts)
  {
    bool isDGcandidate = true;

    // is this a central diffractive event?
    // by default it is assumed to be a MB event
    bool isPythiaDiff = false;
    bool isGraniittiDiff = false;
    if (collision.has_mcCollision()) {
      auto MCCol = collision.mcCollision();
      auto MCPartSlice = McParts.sliceBy(perMcCollision, MCCol.globalIndex());
      isPythiaDiff = isPythiaCDE(MCPartSlice);
      isGraniittiDiff = isGraniittiCDE(MCPartSlice);
    }

    // global tracks
    Partition<TCs> goodTracks = requireGlobalTrackInFilter();
    goodTracks.bindTable(tracks);

    // update collision histograms
    if (isPythiaDiff) {
      registry.get<TH1>(HIST("StatDiff1"))->Fill(0., 1.);
      registry.get<TH2>(HIST("vtxPosxyDiff1"))->Fill(collision.posX(), collision.posY());
      registry.get<TH1>(HIST("vtxPoszDiff1"))->Fill(collision.posZ());
      registry.get<TH1>(HIST("TracksDiff1"))->Fill(tracks.size());
      registry.get<TH1>(HIST("vtxTracksDiff1"))->Fill(collision.numContrib());
      registry.get<TH1>(HIST("globalTracksDiff1"))->Fill(goodTracks.size());
    } else if (isGraniittiDiff) {
      registry.get<TH1>(HIST("StatDiff2"))->Fill(0., 1.);
      registry.get<TH2>(HIST("vtxPosxyDiff2"))->Fill(collision.posX(), collision.posY());
      registry.get<TH1>(HIST("vtxPoszDiff2"))->Fill(collision.posZ());
      registry.get<TH1>(HIST("TracksDiff2"))->Fill(tracks.size());
      registry.get<TH1>(HIST("vtxTracksDiff2"))->Fill(collision.numContrib());
      registry.get<TH1>(HIST("globalTracksDiff2"))->Fill(goodTracks.size());
    } else {
      registry.get<TH1>(HIST("Stat"))->Fill(0., 1.);
      registry.get<TH2>(HIST("vtxPosxy"))->Fill(collision.posX(), collision.posY());
      registry.get<TH1>(HIST("vtxPosz"))->Fill(collision.posZ());
      registry.get<TH1>(HIST("Tracks"))->Fill(tracks.size());
      registry.get<TH1>(HIST("vtxTracks"))->Fill(collision.numContrib());
      registry.get<TH1>(HIST("globalTracks"))->Fill(goodTracks.size());
    }

    // test influence of BCrange width
    for (int NDtcoll = 0; NDtcoll < 10; NDtcoll++) {
      auto bcSlice = MCcompatibleBCs(collision, NDtcoll, bct0s, 0);
      isDGcandidate = true;
      for (auto& bc : bcSlice) {
        isDGcandidate &= cleanFIT(bc, diffCuts.FITAmpLimits());
      }
      if (isPythiaDiff) {
        registry.get<TH2>(HIST("cleanFITDiff1"))->Fill(NDtcoll, isDGcandidate * 1.);
      } else if (isGraniittiDiff) {
        registry.get<TH2>(HIST("cleanFITDiff2"))->Fill(NDtcoll, isDGcandidate * 1.);
      } else {
        registry.get<TH2>(HIST("cleanFIT"))->Fill(NDtcoll, isDGcandidate * 1.);
      }
    }

    // number of vertex tracks with TOF hit
    float rgtrwTOF = 0.;
    for (auto& track : tracks) {
      // update eta vs pt and dEdx histograms histogram
      if (isPythiaDiff) {
        registry.get<TH2>(HIST("etaptDiff1"))->Fill(track.eta(), track.pt());
        registry.get<TH2>(HIST("dEdxTPCDiff1"))->Fill(track.p(), track.tpcSignal());
      } else if (isGraniittiDiff) {
        registry.get<TH2>(HIST("etaptDiff2"))->Fill(track.eta(), track.pt());
        registry.get<TH2>(HIST("dEdxTPCDiff2"))->Fill(track.p(), track.tpcSignal());
      } else {
        registry.get<TH2>(HIST("etapt"))->Fill(track.eta(), track.pt());
        registry.get<TH2>(HIST("dEdxTPC"))->Fill(track.p(), track.tpcSignal());
      }
      if (track.tpcSignal() > maxdEdxTPC) {
        maxdEdxTPC = track.tpcSignal();
        LOGF(info, "<DiffMCQA> New maxdEdx TPC %f", maxdEdxTPC);
      }

      // TOF hit?
      if (track.hasTOF()) {
        if (isPythiaDiff) {
          registry.get<TH2>(HIST("dEdxTOFDiff1"))->Fill(track.p(), track.tofSignal());
        } else if (isGraniittiDiff) {
          registry.get<TH2>(HIST("dEdxTOFDiff2"))->Fill(track.p(), track.tofSignal());
        } else {
          registry.get<TH2>(HIST("dEdxTOF"))->Fill(track.p(), track.tofSignal());
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
      registry.get<TH2>(HIST("tResvsrTOFTracksDiff1"))->Fill(collision.collisionTimeRes(), rgtrwTOF);
    } else if (isGraniittiDiff) {
      registry.get<TH2>(HIST("tResvsrTOFTracksDiff2"))->Fill(collision.collisionTimeRes(), rgtrwTOF);
    } else {
      registry.get<TH2>(HIST("tResvsrTOFTracks"))->Fill(collision.collisionTimeRes(), rgtrwTOF);
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
    auto bcSlice = MCcompatibleBCs(collision, diffCuts.NDtcoll(), bct0s, diffCuts.minNBCs());

    // no FIT signal in bcSlice
    for (auto& bc : bcSlice) {
      if (!cleanFIT(bc, diffCuts.FITAmpLimits())) {
        isDGcandidate = false;
        break;
      }
    }
    if (isPythiaDiff) {
      registry.get<TH1>(HIST("StatDiff1"))->Fill(1., isDGcandidate * 1.);
    } else if (isGraniittiDiff) {
      registry.get<TH1>(HIST("StatDiff2"))->Fill(1., isDGcandidate * 1.);
    } else {
      registry.get<TH1>(HIST("Stat"))->Fill(1., isDGcandidate * 1.);
    }

    // no Zdc signal in bcSlice
    std::vector<float> lims(10, 0.);
    for (auto& bc : bcSlice) {
      if (!cleanZDC(bc, zdcs, lims)) {
        isDGcandidate = false;
        break;
      }
    }
    if (isPythiaDiff) {
      registry.get<TH1>(HIST("StatDiff1"))->Fill(2., isDGcandidate * 1.);
    } else if (isGraniittiDiff) {
      registry.get<TH1>(HIST("StatDiff2"))->Fill(2., isDGcandidate * 1.);
    } else {
      registry.get<TH1>(HIST("Stat"))->Fill(2., isDGcandidate * 1.);
    }

    // no Calo signal in bcSlice
    for (auto& bc : bcSlice) {
      if (!cleanCalo(bc, calos, lims)) {
        isDGcandidate = false;
        break;
      }
    }
    if (isPythiaDiff) {
      registry.get<TH1>(HIST("StatDiff1"))->Fill(3., isDGcandidate * 1.);
    } else if (isGraniittiDiff) {
      registry.get<TH1>(HIST("StatDiff2"))->Fill(3., isDGcandidate * 1.);
    } else {
      registry.get<TH1>(HIST("Stat"))->Fill(3., isDGcandidate * 1.);
    }

    // no V0s
    isDGcandidate &= (v0s.size() == 0);
    if (isPythiaDiff) {
      registry.get<TH1>(HIST("StatDiff1"))->Fill(4., isDGcandidate * 1.);
    } else if (isGraniittiDiff) {
      registry.get<TH1>(HIST("StatDiff2"))->Fill(4., isDGcandidate * 1.);
    } else {
      registry.get<TH1>(HIST("Stat"))->Fill(4., isDGcandidate * 1.);
    }

    // no Cascades
    isDGcandidate &= (cascades.size() == 0);
    if (isPythiaDiff) {
      registry.get<TH1>(HIST("StatDiff1"))->Fill(5., isDGcandidate * 1.);
    } else if (isGraniittiDiff) {
      registry.get<TH1>(HIST("StatDiff2"))->Fill(5., isDGcandidate * 1.);
    } else {
      registry.get<TH1>(HIST("Stat"))->Fill(5., isDGcandidate * 1.);
    }

    // number of forward tracks = 0
    isDGcandidate &= (fwdtracks.size() == 0);
    if (isPythiaDiff) {
      registry.get<TH1>(HIST("StatDiff1"))->Fill(6., isDGcandidate * 1.);
    } else if (isGraniittiDiff) {
      registry.get<TH1>(HIST("StatDiff2"))->Fill(6., isDGcandidate * 1.);
    } else {
      registry.get<TH1>(HIST("Stat"))->Fill(6., isDGcandidate * 1.);
    }

    // no global tracks which are no vtx tracks
    bool globalAndVtx = isDGcandidate;
    bool vtxAndGlobal = isDGcandidate;
    for (auto& track : tracks) {
      if (track.isGlobalTrack() && !track.isPVContributor()) {
        globalAndVtx = false;
      }
      if (track.isPVContributor() && !track.isGlobalTrack()) {
        vtxAndGlobal = false;
      }
    }
    if (isPythiaDiff) {
      registry.get<TH1>(HIST("StatDiff1"))->Fill(7., globalAndVtx * 1.);
      registry.get<TH1>(HIST("StatDiff1"))->Fill(8., vtxAndGlobal * 1.);
    } else if (isGraniittiDiff) {
      registry.get<TH1>(HIST("StatDiff2"))->Fill(7., globalAndVtx * 1.);
      registry.get<TH1>(HIST("StatDiff2"))->Fill(8., vtxAndGlobal * 1.);
    } else {
      registry.get<TH1>(HIST("Stat"))->Fill(7., globalAndVtx * 1.);
      registry.get<TH1>(HIST("Stat"))->Fill(8., vtxAndGlobal * 1.);
    }
    isDGcandidate &= globalAndVtx;
    if (diffCuts.globalTracksOnly()) {
      isDGcandidate &= vtxAndGlobal;
    }

    // check a given bc for possible ambiguous Tracks
    auto noAmbTracks = isDGcandidate;
    for (auto& bc : bcSlice) {
      if (abcrs.isInRange(bc.globalIndex())) {
        noAmbTracks = false;
        break;
      }
    }
    if (isPythiaDiff) {
      registry.get<TH1>(HIST("StatDiff1"))->Fill(9., noAmbTracks * 1.);
    } else if (isGraniittiDiff) {
      registry.get<TH1>(HIST("StatDiff2"))->Fill(9., noAmbTracks * 1.);
    } else {
      registry.get<TH1>(HIST("Stat"))->Fill(9., noAmbTracks * 1.);
    }

    // check a given bc for possible ambiguous FwdTracks
    auto noAmbFwdTracks = isDGcandidate;
    for (auto& bc : bcSlice) {
      if (afbcrs.isInRange(bc.globalIndex())) {
        noAmbFwdTracks = false;
        break;
      }
    }
    if (isPythiaDiff) {
      registry.get<TH1>(HIST("StatDiff1"))->Fill(10., noAmbFwdTracks * 1.);
    } else if (isGraniittiDiff) {
      registry.get<TH1>(HIST("StatDiff2"))->Fill(10., noAmbFwdTracks * 1.);
    } else {
      registry.get<TH1>(HIST("Stat"))->Fill(10., noAmbFwdTracks * 1.);
    }

    // at least one vtx track with TOF hit
    isDGcandidate &= (rgtrwTOF >= diffCuts.minRgtrwTOF());
    if (isPythiaDiff) {
      registry.get<TH1>(HIST("StatDiff1"))->Fill(11., isDGcandidate * 1.);
    } else if (isGraniittiDiff) {
      registry.get<TH1>(HIST("StatDiff2"))->Fill(11., isDGcandidate * 1.);
    } else {
      registry.get<TH1>(HIST("Stat"))->Fill(11., isDGcandidate * 1.);
    }

    // number of vertex tracks <= n
    isDGcandidate &= (collision.numContrib() >= diffCuts.minNTracks());
    if (isPythiaDiff) {
      registry.get<TH1>(HIST("StatDiff1"))->Fill(12., isDGcandidate * 1.);
    } else if (isGraniittiDiff) {
      registry.get<TH1>(HIST("StatDiff2"))->Fill(12., isDGcandidate * 1.);
    } else {
      registry.get<TH1>(HIST("Stat"))->Fill(12., isDGcandidate * 1.);
    }
    isDGcandidate &= (collision.numContrib() <= diffCuts.maxNTracks());
    if (isPythiaDiff) {
      registry.get<TH1>(HIST("StatDiff1"))->Fill(13., isDGcandidate * 1.);
    } else if (isGraniittiDiff) {
      registry.get<TH1>(HIST("StatDiff2"))->Fill(13., isDGcandidate * 1.);
    } else {
      registry.get<TH1>(HIST("Stat"))->Fill(13., isDGcandidate * 1.);
    }

    // net charge and invariant mass
    bool goodetas = true;
    bool goodpts = true;
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
        // update histogram rejectedTracks
        if (!track.isPVContributor()) {
          if (isPythiaDiff) {
            registry.get<TH1>(HIST("rejectedTracksDiff1"))->Fill(0., 1.);
            registry.get<TH1>(HIST("rejectedTracksDiff1"))->Fill(1., track.isGlobalTrackSDD() * 1.);
            registry.get<TH1>(HIST("rejectedTracksDiff1"))->Fill(2., track.passedTrackType() * 1.);
            registry.get<TH1>(HIST("rejectedTracksDiff1"))->Fill(3., track.passedPtRange() * 1.);
            registry.get<TH1>(HIST("rejectedTracksDiff1"))->Fill(4., track.passedEtaRange() * 1.);
            registry.get<TH1>(HIST("rejectedTracksDiff1"))->Fill(5., track.passedTPCNCls() * 1.);
            registry.get<TH1>(HIST("rejectedTracksDiff1"))->Fill(6., track.passedTPCCrossedRows() * 1.);
            registry.get<TH1>(HIST("rejectedTracksDiff1"))->Fill(7., track.passedTPCCrossedRowsOverNCls() * 1.);
            registry.get<TH1>(HIST("rejectedTracksDiff1"))->Fill(8., track.passedTPCChi2NDF() * 1.);
            registry.get<TH1>(HIST("rejectedTracksDiff1"))->Fill(9., track.passedTPCRefit() * 1.);
            registry.get<TH1>(HIST("rejectedTracksDiff1"))->Fill(10., track.passedITSNCls() * 1.);
            registry.get<TH1>(HIST("rejectedTracksDiff1"))->Fill(11., track.passedITSChi2NDF() * 1.);
            registry.get<TH1>(HIST("rejectedTracksDiff1"))->Fill(12., track.passedITSRefit() * 1.);
            registry.get<TH1>(HIST("rejectedTracksDiff1"))->Fill(13., track.passedITSHits() * 1.);
            registry.get<TH1>(HIST("rejectedTracksDiff1"))->Fill(14., track.passedGoldenChi2() * 1.);
            registry.get<TH1>(HIST("rejectedTracksDiff1"))->Fill(15., track.passedDCAxy() * 1.);
            registry.get<TH1>(HIST("rejectedTracksDiff1"))->Fill(16., track.passedDCAz() * 1.);
          } else if (isGraniittiDiff) {
            registry.get<TH1>(HIST("rejectedTracksDiff2"))->Fill(0., 1.);
            registry.get<TH1>(HIST("rejectedTracksDiff2"))->Fill(1., track.isGlobalTrackSDD() * 1.);
            registry.get<TH1>(HIST("rejectedTracksDiff2"))->Fill(2., track.passedTrackType() * 1.);
            registry.get<TH1>(HIST("rejectedTracksDiff2"))->Fill(3., track.passedPtRange() * 1.);
            registry.get<TH1>(HIST("rejectedTracksDiff2"))->Fill(4., track.passedEtaRange() * 1.);
            registry.get<TH1>(HIST("rejectedTracksDiff2"))->Fill(5., track.passedTPCNCls() * 1.);
            registry.get<TH1>(HIST("rejectedTracksDiff2"))->Fill(6., track.passedTPCCrossedRows() * 1.);
            registry.get<TH1>(HIST("rejectedTracksDiff2"))->Fill(7., track.passedTPCCrossedRowsOverNCls() * 1.);
            registry.get<TH1>(HIST("rejectedTracksDiff2"))->Fill(8., track.passedTPCChi2NDF() * 1.);
            registry.get<TH1>(HIST("rejectedTracksDiff2"))->Fill(9., track.passedTPCRefit() * 1.);
            registry.get<TH1>(HIST("rejectedTracksDiff2"))->Fill(10., track.passedITSNCls() * 1.);
            registry.get<TH1>(HIST("rejectedTracksDiff2"))->Fill(11., track.passedITSChi2NDF() * 1.);
            registry.get<TH1>(HIST("rejectedTracksDiff2"))->Fill(12., track.passedITSRefit() * 1.);
            registry.get<TH1>(HIST("rejectedTracksDiff2"))->Fill(13., track.passedITSHits() * 1.);
            registry.get<TH1>(HIST("rejectedTracksDiff2"))->Fill(14., track.passedGoldenChi2() * 1.);
            registry.get<TH1>(HIST("rejectedTracksDiff2"))->Fill(15., track.passedDCAxy() * 1.);
            registry.get<TH1>(HIST("rejectedTracksDiff2"))->Fill(16., track.passedDCAz() * 1.);
          } else {
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
      registry.get<TH1>(HIST("StatDiff1"))->Fill(14., isDGcandidate * 1.);
    } else if (isGraniittiDiff) {
      registry.get<TH1>(HIST("StatDiff2"))->Fill(14., isDGcandidate * 1.);
    } else {
      registry.get<TH1>(HIST("Stat"))->Fill(14., isDGcandidate * 1.);
    }
    isDGcandidate &= goodetas;
    if (isPythiaDiff) {
      registry.get<TH1>(HIST("StatDiff1"))->Fill(15., isDGcandidate * 1.);
    } else if (isGraniittiDiff) {
      registry.get<TH1>(HIST("StatDiff2"))->Fill(15., isDGcandidate * 1.);
    } else {
      registry.get<TH1>(HIST("Stat"))->Fill(15., isDGcandidate * 1.);
    }
    isDGcandidate &= (netCharge >= diffCuts.minNetCharge());
    if (isPythiaDiff) {
      registry.get<TH1>(HIST("StatDiff1"))->Fill(16., isDGcandidate * 1.);
    } else if (isGraniittiDiff) {
      registry.get<TH1>(HIST("StatDiff2"))->Fill(16., isDGcandidate * 1.);
    } else {
      registry.get<TH1>(HIST("Stat"))->Fill(16., isDGcandidate * 1.);
    }
    isDGcandidate &= (netCharge <= diffCuts.maxNetCharge());
    if (isPythiaDiff) {
      registry.get<TH1>(HIST("StatDiff1"))->Fill(17., isDGcandidate * 1.);
    } else if (isGraniittiDiff) {
      registry.get<TH1>(HIST("StatDiff2"))->Fill(17., isDGcandidate * 1.);
    } else {
      registry.get<TH1>(HIST("Stat"))->Fill(17., isDGcandidate * 1.);
    }
    isDGcandidate &= (ivm.M() >= diffCuts.minIVM());
    if (isPythiaDiff) {
      registry.get<TH1>(HIST("StatDiff1"))->Fill(18., isDGcandidate * 1.);
    } else if (isGraniittiDiff) {
      registry.get<TH1>(HIST("StatDiff2"))->Fill(18., isDGcandidate * 1.);
    } else {
      registry.get<TH1>(HIST("Stat"))->Fill(18., isDGcandidate * 1.);
    }
    isDGcandidate &= (ivm.M() <= diffCuts.maxIVM());
    if (isPythiaDiff) {
      registry.get<TH1>(HIST("StatDiff1"))->Fill(19., isDGcandidate * 1.);
    } else if (isGraniittiDiff) {
      registry.get<TH1>(HIST("StatDiff2"))->Fill(19., isDGcandidate * 1.);
    } else {
      registry.get<TH1>(HIST("Stat"))->Fill(19., isDGcandidate * 1.);
    }

    // update some DG histograms
    if (isDGcandidate) {
      if (isPythiaDiff) {
        registry.get<TH2>(HIST("vtxPosxyDGDiff1"))->Fill(collision.posX(), collision.posY());
        registry.get<TH1>(HIST("vtxPoszDGDiff1"))->Fill(collision.posZ());
        registry.get<TH1>(HIST("netChargeDGDiff1"))->Fill(netCharge);
        registry.get<TH2>(HIST("IVMptSysDGDiff1"))->Fill(ivm.M(), ivm.Perp());
      } else if (isGraniittiDiff) {
        registry.get<TH2>(HIST("vtxPosxyDGDiff2"))->Fill(collision.posX(), collision.posY());
        registry.get<TH1>(HIST("vtxPoszDGDiff2"))->Fill(collision.posZ());
        registry.get<TH1>(HIST("netChargeDGDiff2"))->Fill(netCharge);
        registry.get<TH2>(HIST("IVMptSysDGDiff2"))->Fill(ivm.M(), ivm.Perp());
      } else {
        registry.get<TH2>(HIST("vtxPosxyDG"))->Fill(collision.posX(), collision.posY());
        registry.get<TH1>(HIST("vtxPoszDG"))->Fill(collision.posZ());
        registry.get<TH1>(HIST("netChargeDG"))->Fill(netCharge);
        registry.get<TH2>(HIST("IVMptSysDG"))->Fill(ivm.M(), ivm.Perp());
      }

      // fill dEdx of DG event tracks
      for (auto& track : tracks) {
        if (track.isPVContributor()) {
          if (isPythiaDiff) {
            registry.get<TH2>(HIST("etaptDGDiff1"))->Fill(track.eta(), track.pt());
            registry.get<TH2>(HIST("dEdxTPCDGDiff1"))->Fill(track.p(), track.tpcSignal());
            registry.get<TH2>(HIST("IVMptTrkDGDiff1"))->Fill(ivm.M(), track.pt());
          } else if (isGraniittiDiff) {
            registry.get<TH2>(HIST("etaptDGDiff2"))->Fill(track.eta(), track.pt());
            registry.get<TH2>(HIST("dEdxTPCDGDiff2"))->Fill(track.p(), track.tpcSignal());
            registry.get<TH2>(HIST("IVMptTrkDGDiff2"))->Fill(ivm.M(), track.pt());
          } else {
            registry.get<TH2>(HIST("etaptDG"))->Fill(track.eta(), track.pt());
            registry.get<TH2>(HIST("dEdxTPCDG"))->Fill(track.p(), track.tpcSignal());
            registry.get<TH2>(HIST("IVMptTrkDG"))->Fill(ivm.M(), track.pt());
          }
          if (track.hasTOF()) {
            if (isPythiaDiff) {
              registry.get<TH2>(HIST("dEdxTOFDGDiff1"))->Fill(track.p(), track.tofSignal());
            } else if (isGraniittiDiff) {
              registry.get<TH2>(HIST("dEdxTOFDGDiff2"))->Fill(track.p(), track.tofSignal());
            } else {
              registry.get<TH2>(HIST("dEdxTOFDG"))->Fill(track.p(), track.tofSignal());
            }
          }
        }
      }
    }
  };
};

struct FV0Signals {

  HistogramRegistry registry{
    "registry",
    {
      {"FV0A", "#FV0A", {HistType::kTH2F, {{48, -0.5, 47.5}, {1000, 0., 1000.}}}},
    }};

  void process(aod::FV0A const& fv0)
  {
    // side A
    for (size_t ind = 0; ind < fv0.channel().size(); ind++) {
      registry.get<TH2>(HIST("FV0A"))->Fill((fv0.channel())[ind], (fv0.amplitude())[ind]);
    }
  };
};

struct FT0Signals {

  HistogramRegistry registry{
    "registry",
    {
      {"FT0A", "#FT0A", {HistType::kTH2F, {{96, -0.5, 95.5}, {100, 0., 200.}}}},
      {"FT0C", "#FT0C", {HistType::kTH2F, {{112, -0.5, 111.5}, {100, 0., 200.}}}},
    }};

  void process(aod::FT0 const& ft0)
  {
    // side A
    for (size_t ind = 0; ind < ft0.channelA().size(); ind++) {
      registry.get<TH2>(HIST("FT0A"))->Fill((ft0.channelA())[ind], (ft0.amplitudeA())[ind]);
    }

    // side C
    for (size_t ind = 0; ind < ft0.channelC().size(); ind++) {
      registry.get<TH2>(HIST("FT0C"))->Fill((ft0.channelC())[ind], (ft0.amplitudeC())[ind]);
    }
  };
};

struct FDDSignals {

  HistogramRegistry registry{
    "registry",
    {
      {"FDDA", "#FDDA", {HistType::kTH2F, {{8, -0.5, 7.5}, {100, 0., 100.}}}},
      {"FDDC", "#FDDC", {HistType::kTH2F, {{8, -0.5, 7.5}, {100, 0., 100.}}}},
    }};

  void process(aod::FDD const& fdd)
  {
    // side A
    for (auto ind = 0; ind < 8; ind++) {
      registry.get<TH2>(HIST("FDDA"))->Fill(ind, (fdd.chargeA())[ind]);
    }

    // side C
    for (auto ind = 0; ind < 8; ind++) {
      registry.get<TH2>(HIST("FDDC"))->Fill(ind, (fdd.chargeC())[ind]);
    }
  };
};

struct ZDCSignals {

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
  HistogramRegistry registry{
    "registry",
    {
      {"ZdcEnergies", "#ZdcEnergies", {HistType::kTH2F, {{22, -0.5, 21.5}, {100, 0., 1000.}}}},
    }};

  void process(aod::Zdc const& zdc)
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
  };
};

struct CaloSignals {

  HistogramRegistry registry{
    "registry",
    {
      {"CaloCell", "#CaloCell", {HistType::kTH1I, {{18000, -0.5, 17999.5}}}},
      {"CaloAmplitude", "#CaloAmplitude", {HistType::kTH1F, {{100, 0, 10.}}}},
    }};

  void process(aod::Calo const& calo)
  {
    // cell number
    registry.get<TH1>(HIST("CaloCell"))->Fill(calo.cellNumber());

    // amplitude
    registry.get<TH1>(HIST("CaloAmplitude"))->Fill(calo.amplitude());
  };
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<DiffMCQA>(cfgc, TaskName{"diffmcqa"}),
    adaptAnalysisTask<FV0Signals>(cfgc, TaskName{"fv0signals"}),
    adaptAnalysisTask<FT0Signals>(cfgc, TaskName{"ft0signals"}),
    adaptAnalysisTask<FDDSignals>(cfgc, TaskName{"fddsignals"}),
    adaptAnalysisTask<ZDCSignals>(cfgc, TaskName{"zdcsignals"}),
    adaptAnalysisTask<CaloSignals>(cfgc, TaskName{"calosignals"}),
  };
}
