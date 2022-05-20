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
#include "CommonConstants/PhysicsConstants.h"
#include "EventFiltering/PWGUD/diffHelpers.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct DiffQA {

  // get a cutHolder
  cutHolder diffCuts = cutHolder();
  MutableConfigurable<cutHolder> DGCuts{"DGCuts", {}, "DG event cuts"};

  // define histograms
  HistogramRegistry registry{
    "registry",
    {
      {"Stat", "#Stat", {HistType::kTH1F, {{13, -0.5, 12.5}}}},
      {"NTracksStat", "#NTracksStat", {HistType::kTH1F, {{20, 0.5, 20.5}}}},
      {"cleanFIT", "#cleanFIT", {HistType::kTH2F, {{10, -0.5, 9.5}, {2, -0.5, 1.5}}}},
      {"tResvsrTOFTracks", "#tResvsrTOFTracks", {HistType::kTH2F, {{1000, 0., 1.E3}, {101, -0.01, 1.01}}}},
      {"dEdxTPC", "#dEdxTPC", {HistType::kTH2F, {{100, 0., 5.0}, {300, 0., 300.}}}},
      {"dEdxTOF", "#dEdxTOF", {HistType::kTH2F, {{100, 0., 5.0}, {100, 0., 2000.}}}},
    }};

  void init(InitContext&)
  {
    diffCuts = (cutHolder)DGCuts;
  }

  using CCs = soa::Join<aod::Collisions, aod::EvSels>;
  using CC = CCs::iterator;
  using BCs = soa::Join<aod::BCs, aod::BcSels, aod::Run3MatchedToBCSparse>;
  using TCs = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TOFSignal>;
  using FWs = aod::FwdTracks;

  void process(CC const& collision, BCs const& bct0s,
               TCs& tracks, FWs& fwdtracks, aod::FT0s& ft0s, aod::FV0As& fv0as, aod::FDDs& fdds)
  {
    LOGF(debug, "<DiffQA> Start");

    // test influence of BCrange width
    bool isDGcandidate;
    for (int NDtcoll = 0; NDtcoll < 10; NDtcoll++) {
      auto bcSlice = compatibleBCs(collision, NDtcoll, bct0s, 0);
      isDGcandidate = true;
      for (auto& bc : bcSlice) {
        isDGcandidate &= cleanFIT(bc, diffCuts.FITAmpLimits());
      }
      registry.get<TH2>(HIST("cleanFIT"))->Fill(NDtcoll, isDGcandidate * 1.);
    }

    // get BCrange to test for FIT signals
    auto bcSlice = compatibleBCs(collision, diffCuts.NDtcoll(), bct0s, diffCuts.minNBCs());

    // check tracks
    // global tracks
    Partition<TCs> goodTracks = requireGlobalTrackInFilter();
    goodTracks.bindTable(tracks);
    LOGF(debug, "<DiffQA> Number of good tracks: %i", goodTracks.size());

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
    registry.get<TH2>(HIST("tResvsrTOFTracks"))->Fill(collision.collisionTimeRes(), rgtrwTOF);

    // is it a DG candidate?
    // DG = no FIT signal in compatible BCs
    //    & number of forward tracks = 0
    //    & ntrMin <= number of vertex tracks <= ntrMax
    //    & no global track which is not a vertex track
    isDGcandidate = true;
    registry.get<TH1>(HIST("Stat"))->Fill(0., isDGcandidate * 1.);

    // no FIT signal in compatible BCs
    for (auto& bc : bcSlice) {
      if (!cleanFIT(bc, diffCuts.FITAmpLimits())) {
        isDGcandidate = false;
        continue;
      }
    }
    registry.get<TH1>(HIST("Stat"))->Fill(1., isDGcandidate * 1.);

    // number of forward tracks = 0
    isDGcandidate &= (fwdtracks.size() == 0);
    registry.get<TH1>(HIST("Stat"))->Fill(2., isDGcandidate * 1.);

    // no global tracks which are no vtx tracks
    for (auto& track : tracks) {
      if (track.isGlobalTrack() && !track.isPVContributor()) {
        isDGcandidate = false;
        continue;
      }
    }
    registry.get<TH1>(HIST("Stat"))->Fill(3., isDGcandidate * 1.);

    // number of vertex tracks <= n
    for (int ii = 1; ii <= 20; ii++) {
      registry.get<TH1>(HIST("NTracksStat"))->Fill(ii, isDGcandidate * (collision.numContrib() <= ii) * 1.);
    }
    isDGcandidate &= (collision.numContrib() >= diffCuts.minNTracks());
    registry.get<TH1>(HIST("Stat"))->Fill(4., isDGcandidate * 1.);
    isDGcandidate &= (collision.numContrib() <= diffCuts.maxNTracks());
    registry.get<TH1>(HIST("Stat"))->Fill(5., isDGcandidate * 1.);

    // invariant mass
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

      for (auto& track : tracks) {
        if (track.isPVContributor()) {
          LOGF(info, "dEdx TPC %f TOF %i %f", track.tpcSignal(), track.hasTOF(), track.hasTOF() ? track.tofSignal() : 0.);
          registry.get<TH2>(HIST("dEdxTPC"))->Fill(track.pt(), track.tpcSignal());
          if (track.hasTOF()) {
            registry.get<TH2>(HIST("dEdxTOF"))->Fill(track.pt(), track.tofSignal());
          }

          lvtmp.SetXYZM(track.px(), track.py(), track.pz(), mass2Use);
          if (lvtmp.Perp() < diffCuts.minPt() || lvtmp.Perp() > diffCuts.maxPt()) {
            goodpts = false;
          }
          if (lvtmp.Eta() < diffCuts.minEta() || lvtmp.Eta() > diffCuts.maxEta()) {
            goodetas = false;
          }
          netCharge += track.sign();
          ivm += lvtmp;
        }
      }
    }
    isDGcandidate &= goodpts;
    registry.get<TH1>(HIST("Stat"))->Fill(6., isDGcandidate * 1.);
    isDGcandidate &= goodetas;
    registry.get<TH1>(HIST("Stat"))->Fill(7., isDGcandidate * 1.);
    isDGcandidate &= (netCharge >= diffCuts.minNetCharge());
    registry.get<TH1>(HIST("Stat"))->Fill(8., isDGcandidate * 1.);
    isDGcandidate &= (netCharge <= diffCuts.maxNetCharge());
    registry.get<TH1>(HIST("Stat"))->Fill(9., isDGcandidate * 1.);
    isDGcandidate &= (ivm.M() >= diffCuts.minIVM());
    registry.get<TH1>(HIST("Stat"))->Fill(10., isDGcandidate * 1.);
    isDGcandidate &= (ivm.M() <= diffCuts.maxIVM());
    registry.get<TH1>(HIST("Stat"))->Fill(11., isDGcandidate * 1.);

    LOGF(debug, "<DiffQA> End");
  };
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<DiffQA>(cfgc, TaskName{"diffqa"}),
  };
}
