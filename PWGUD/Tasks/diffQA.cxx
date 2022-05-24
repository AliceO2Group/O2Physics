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
/// \brief A QA task for DG events events
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
/// \since  20.05.2022

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "TLorentzVector.h"
#include "ReconstructionDataFormats/BCRange.h"
#include "CommonConstants/PhysicsConstants.h"
#include "EventFiltering/PWGUD/diffHelpers.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct DiffQA {

  // get a cutHolder
  cutHolder diffCuts = cutHolder();
  MutableConfigurable<cutHolder> DGCuts{"DGCuts", {}, "DG event cuts"};

  // a structure to hold information about the possible BCs the ambiguous tracks belong to
  o2::dataformats::bcRanges abcrs = o2::dataformats::bcRanges("ambiguous_tracks");

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
  //  bin  8: possible ambiguous tracks
  //  bin  9: number of tracks >= minimum number
  //  bin 10: number of tracks <= maximum number
  //  bin 11: minimum pt <= pt of vtx tracks <= maximum pt
  //  bin 12: minimum eta <= eta of vtx tracks <= maximum eta
  //  bin 13: net charge >= minimum net charge
  //  bin 14: net charge <= maximum net charge
  //  bin 15: IVM >= minimum IVM
  //  bin 16: IVM <= maximum IVM
  HistogramRegistry registry{
    "registry",
    {
      {"Stat", "#Stat", {HistType::kTH1F, {{17, -0.5, 16.5}}}},
      {"cleanFIT", "#cleanFIT", {HistType::kTH2F, {{10, -0.5, 9.5}, {2, -0.5, 1.5}}}},
      {"Tracks", "#Tracks", {HistType::kTH1F, {{51, 0.5, 50.5}}}},
      {"vtxTracks", "#vtxTracks", {HistType::kTH1F, {{51, 0.5, 50.5}}}},
      {"tResvsrTOFTracks", "#tResvsrTOFTracks", {HistType::kTH2F, {{1000, 0., 1.E3}, {101, -0.01, 1.01}}}},
      {"vtxPosxy", "#vtxPosxy", {HistType::kTH2F, {{200, -2., 2.}, {200, -2., 2.}}}},
      {"vtxPosz", "#vtxPosz", {HistType::kTH1F, {{1000, -100., 100.}}}},
      {"dEdxTPC", "#dEdxTPC", {HistType::kTH2F, {{100, 0., 5.0}, {300, 0., 300.}}}},
      {"dEdxTOF", "#dEdxTOF", {HistType::kTH2F, {{100, 0., 5.0}, {100, 0., 50000.}}}},
      {"dEdxTPCDG", "#dEdxTPCDG", {HistType::kTH2F, {{100, 0., 5.0}, {300, 0., 300.}}}},
      {"dEdxTOFDG", "#dEdxTOFDG", {HistType::kTH2F, {{100, 0., 5.0}, {100, 0., 50000.}}}},
      {"vtxPosxyDG", "#vtxPosxyDG", {HistType::kTH2F, {{100, -10., 10.}, {100, -10., 10.}}}},
      {"vtxPoszDG", "#vtxPoszDG", {HistType::kTH1F, {{1000, -100., 100.}}}},
    }};

  using CCs = soa::Join<aod::Collisions, aod::EvSels>;
  using CC = CCs::iterator;
  using BCs = soa::Join<aod::BCs, aod::BcSels, aod::Run3MatchedToBCSparse>;
  using TCs = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TOFSignal>;
  using FWs = aod::FwdTracks;
  using ATs = aod::AmbiguousTracks;

  void init(InitContext&)
  {
    diffCuts = (cutHolder)DGCuts;
  }

  void run(ProcessingContext& pc)
  {
    // get ambiguous tracks table
    auto t1 = pc.inputs().get<TableConsumer>("BCs")->asArrowTable();
    auto t2 = pc.inputs().get<TableConsumer>("BcSels")->asArrowTable();
    auto t3 = pc.inputs().get<TableConsumer>("Run3MatchedToBCSparse")->asArrowTable();
    auto t4 = pc.inputs().get<TableConsumer>("AmbiguousTracks")->asArrowTable();
    auto bcs = BCs({t1, t2, t3});
    auto ambtracks = ATs({t4});
    ambtracks.bindExternalIndices(&bcs);

    // make sorted list of BC ranges, this is used to efficiently check whether a given BC
    // is contained in one of these ranges
    abcrs.reset();
    for (auto ambtrack : ambtracks) {
      auto bcfirst = ambtrack.bc().rawIteratorAt(0);
      auto bclast = ambtrack.bc().rawIteratorAt(ambtrack.bc().size() - 1);
      abcrs.add(bcfirst.globalIndex(), bclast.globalIndex());
    }
    abcrs.merge();
  }

  void process(CC const& collision, BCs const& bct0s,
               TCs& tracks, FWs& fwdtracks, ATs& ambtracks,
               aod::FT0s& ft0s, aod::FV0As& fv0as, aod::FDDs& fdds,
               aod::Zdcs& zdcs, aod::Calos& calos,
               aod::V0s& v0s, aod::Cascades& cascades)
  {
    LOGF(debug, "<DiffQA> Start %i", abcrs.size());
    bool isDGcandidate = true;
    registry.get<TH1>(HIST("Stat"))->Fill(0., isDGcandidate * 1.);

    // test influence of BCrange width
    for (int NDtcoll = 0; NDtcoll < 10; NDtcoll++) {
      auto bcSlice = compatibleBCs(collision, NDtcoll, bct0s, 0);
      isDGcandidate = true;
      for (auto& bc : bcSlice) {
        isDGcandidate &= cleanFIT(bc, diffCuts.FITAmpLimits());
      }
      registry.get<TH2>(HIST("cleanFIT"))->Fill(NDtcoll, isDGcandidate * 1.);
    }

    // vertex position
    registry.get<TH2>(HIST("vtxPosxy"))->Fill(collision.posX(), collision.posY());
    registry.get<TH1>(HIST("vtxPosz"))->Fill(collision.posZ());

    // check tracks
    registry.get<TH1>(HIST("Tracks"))->Fill(tracks.size());

    // global tracks
    Partition<TCs> goodTracks = requireGlobalTrackInFilter();
    goodTracks.bindTable(tracks);
    LOGF(debug, "<DiffQA> Number of good tracks: %i", goodTracks.size());

    // number of gobal tracks with TOF hit
    float rgtrwTOF = 0.;
    if (goodTracks.size() > 0) {
      for (auto& track : goodTracks) {
        // update dEdx histograms
        registry.get<TH2>(HIST("dEdxTPC"))->Fill(track.pt(), track.tpcSignal());
        if (track.hasTOF()) {
          registry.get<TH2>(HIST("dEdxTOF"))->Fill(track.pt(), track.tofSignal());

          // number of global tracks with TOF hit
          rgtrwTOF += 1.;
        }
      }
      rgtrwTOF /= goodTracks.size();
    }
    LOGF(debug, "<DiffQA> Good tracks with TOF: %f [1]", rgtrwTOF);
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
    auto bcSlice = compatibleBCs(collision, diffCuts.NDtcoll(), bct0s, diffCuts.minNBCs());

    // no FIT signal in bcSlice
    for (auto& bc : bcSlice) {
      if (!cleanFIT(bc, diffCuts.FITAmpLimits())) {
        isDGcandidate = false;
        break;
      }
    }
    registry.get<TH1>(HIST("Stat"))->Fill(1., isDGcandidate * 1.);

    // no Zdc signal in bcSlice
    std::vector<float> lims(10, 0.);
    for (auto& bc : bcSlice) {
      if (!cleanZDC(bc, zdcs, lims)) {
        isDGcandidate = false;
        break;
      }
    }
    registry.get<TH1>(HIST("Stat"))->Fill(2., isDGcandidate * 1.);

    // no Calo signal in bcSlice
    for (auto& bc : bcSlice) {
      if (!cleanCalo(bc, calos, lims)) {
        isDGcandidate = false;
        break;
      }
    }
    registry.get<TH1>(HIST("Stat"))->Fill(3., isDGcandidate * 1.);

    // no V0s
    auto colId = collision.globalIndex();
    const auto& V0Collision = v0s.sliceBy(aod::v0::collisionId, colId);
    isDGcandidate &= (V0Collision.size() == 0);
    registry.get<TH1>(HIST("Stat"))->Fill(4., isDGcandidate * 1.);

    // no Cascades
    const auto& CascadeCollision = cascades.sliceBy(aod::cascade::collisionId, colId);
    isDGcandidate &= (CascadeCollision.size() == 0);
    registry.get<TH1>(HIST("Stat"))->Fill(5., isDGcandidate * 1.);

    // number of forward tracks = 0
    isDGcandidate &= (fwdtracks.size() == 0);
    registry.get<TH1>(HIST("Stat"))->Fill(6., isDGcandidate * 1.);

    // no global tracks which are no vtx tracks
    for (auto& track : tracks) {
      if (track.isGlobalTrack() && !track.isPVContributor()) {
        isDGcandidate = false;
        break;
      }
    }
    registry.get<TH1>(HIST("Stat"))->Fill(7., isDGcandidate * 1.);

    // check a given bc for possible ambiguous tracks
    auto withAmbTracks = isDGcandidate;
    for (auto& bc : bcSlice) {
      if (abcrs.isInRange(bc.globalIndex())) {
        withAmbTracks = false;
        break;
      }
    }
    registry.get<TH1>(HIST("Stat"))->Fill(8., withAmbTracks * 1.);

    // number of vertex tracks <= n
    registry.get<TH1>(HIST("vtxTracks"))->Fill(collision.numContrib());
    isDGcandidate &= (collision.numContrib() >= diffCuts.minNTracks());
    registry.get<TH1>(HIST("Stat"))->Fill(9., isDGcandidate * 1.);
    isDGcandidate &= (collision.numContrib() <= diffCuts.maxNTracks());
    registry.get<TH1>(HIST("Stat"))->Fill(10., isDGcandidate * 1.);

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
        lvtmp.SetXYZM(track.px(), track.py(), track.pz(), mass2Use);
        if (lvtmp.Perp() <= diffCuts.minPt() || lvtmp.Perp() >= diffCuts.maxPt()) {
          goodpts = false;
        }
        if (lvtmp.Eta() <= diffCuts.minEta() || lvtmp.Eta() >= diffCuts.maxEta()) {
          goodetas = false;
        }
        netCharge += track.sign();
        ivm += lvtmp;
      }
    }
    isDGcandidate &= goodpts;
    registry.get<TH1>(HIST("Stat"))->Fill(11., isDGcandidate * 1.);
    isDGcandidate &= goodetas;
    registry.get<TH1>(HIST("Stat"))->Fill(12., isDGcandidate * 1.);
    isDGcandidate &= (netCharge >= diffCuts.minNetCharge());
    registry.get<TH1>(HIST("Stat"))->Fill(13., isDGcandidate * 1.);
    isDGcandidate &= (netCharge <= diffCuts.maxNetCharge());
    registry.get<TH1>(HIST("Stat"))->Fill(14., isDGcandidate * 1.);
    isDGcandidate &= (ivm.M() >= diffCuts.minIVM());
    registry.get<TH1>(HIST("Stat"))->Fill(15., isDGcandidate * 1.);
    isDGcandidate &= (ivm.M() <= diffCuts.maxIVM());
    registry.get<TH1>(HIST("Stat"))->Fill(16., isDGcandidate * 1.);

    // update some DG histograms
    if (isDGcandidate) {
      // vertex position of DG events
      registry.get<TH2>(HIST("vtxPosxy"))->Fill(collision.posX(), collision.posY());
      registry.get<TH1>(HIST("vtxPosz"))->Fill(collision.posZ());

      // fill dEdx of DG event tracks
      for (auto& track : tracks) {
        if (track.isPVContributor()) {
          LOGF(debug, "dEdx TPC %f TOF %i %f", track.tpcSignal(), track.hasTOF(), track.hasTOF() ? track.tofSignal() : 0.);
          registry.get<TH2>(HIST("dEdxTPCDG"))->Fill(track.pt(), track.tpcSignal());
          if (track.hasTOF()) {
            registry.get<TH2>(HIST("dEdxTOFDG"))->Fill(track.pt(), track.tofSignal());
          }
        }
      }
    }

    LOGF(debug, "<DiffQA> End");
  };
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<DiffQA>(cfgc, TaskName{"diffqa"}),
  };
}
