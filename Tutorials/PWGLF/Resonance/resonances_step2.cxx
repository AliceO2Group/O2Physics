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
/// \brief this is a starting point for the Resonances tutorial
/// \author
/// \since 02/11/2023

#include <TLorentzVector.h>

#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/runDataProcessing.h"
#include "PWGLF/DataModel/LFResonanceTables.h"
#include "PWGHF/Core/PDG.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// STEP 2
// Producing Mixed event invariant mass distribution
struct resonances_tutorial {
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Configurable for number of bins
  Configurable<int> nBins{"nBins", 100, "N bins in all histos"};

  // Configurable for min pT cut
  Configurable<double> cMinPtcut{"cMinPtcut", 0.15, "Track minium pt cut"};

  // Track selection
  // primary track condition
  Configurable<bool> cfgPrimaryTrack{"cfgPrimaryTrack", true, "Primary track selection"};                    // kGoldenChi2 | kDCAxy | kDCAz
  Configurable<bool> cfgGlobalWoDCATrack{"cfgGlobalWoDCATrack", true, "Global track selection without DCA"}; // kQualityTracks (kTrackType | kTPCNCls | kTPCCrossedRows | kTPCCrossedRowsOverNCls | kTPCChi2NDF | kTPCRefit | kITSNCls | kITSChi2NDF | kITSRefit | kITSHits) | kInAcceptanceTracks (kPtRange | kEtaRange)
  Configurable<bool> cfgPVContributor{"cfgPVContributor", true, "PV contributor track selection"};           // PV Contriuibutor
  // DCA Selections
  // DCAr to PV
  Configurable<double> cMaxDCArToPVcut{"cMaxDCArToPVcut", 0.5, "Track DCAr cut to PV Maximum"};
  // DCAz to PV
  Configurable<double> cMaxDCAzToPVcut{"cMaxDCAzToPVcut", 2.0, "Track DCAz cut to PV Maximum"};
  Configurable<double> cMinDCAzToPVcut{"cMinDCAzToPVcut", 0.0, "Track DCAz cut to PV Minimum"};

  // PID selection
  Configurable<float> nsigmaCutTPC{"nsigmacutTPC", 3.0, "Value of the TPC Nsigma cut"};
  Configurable<float> nsigmaCutCombined{"nsigmaCutCombined", 3.0, "Value of the TOF Nsigma cut"};

  // Event Mixing
  Configurable<int> nEvtMixing{"nEvtMixing", 5, "Number of events to mix"};
  ConfigurableAxis CfgVtxBins{"CfgVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis CfgMultBins{"CfgMultBins", {VARIABLE_WIDTH, 0., 1., 5., 10., 30., 50., 70., 100., 110.}, "Mixing bins - multiplicity"};

  // Initialize the ananlysis task
  void init(o2::framework::InitContext&)
  {
    // register histograms
    histos.add("hVertexZ", "hVertexZ", HistType::kTH1F, {{nBins, -15., 15.}});
    histos.add("InputTracks", "InputTracks", HistType::kTH1F, {{nBins, -5., 5.}});
    histos.add("hMultiplicityPercent", "Multiplicity Percentile", kTH1F, {{2000, 0.0f, 2000.0f}});
    histos.add("hEta", "Eta distribution", kTH1F, {{2000, -1.0f, 1.0f}});
    histos.add("hDcaxy", "Dcaxy distribution", kTH1F, {{2000, -1.0f, 1.0f}});
    histos.add("hDcaz", "Dcaz distribution", kTH1F, {{2000, -1.0f, 1.0f}});
    histos.add("hNsigmaKaonTPC", "NsigmaKaon TPC distribution", kTH1F, {{200, -10.0f, 10.0f}});
    histos.add("hNsigmaKaonTOF", "NsigmaKaon TOF distribution", kTH1F, {{200, -10.0f, 10.0f}});
    histos.add("h3PhiInvMassUnlikeSign", "Invariant mass of Phi meson Unlike Sign", kTH3F, {{2001, -0.5, 2000.5}, {100, 0.0f, 10.0f}, {300, 0.9, 1.2}});
    histos.add("h3PhiInvMassLikeSignPP", "Invariant mass of Phi meson Like Sign positive", kTH3F, {{2001, -0.5, 2000.5}, {100, 0.0f, 10.0f}, {300, 0.9, 1.2}});
    histos.add("h3PhiInvMassLikeSignMM", "Invariant mass of Phi meson Like Sign negative", kTH3F, {{2001, -0.5, 2000.5}, {100, 0.0f, 10.0f}, {300, 0.9, 1.2}});
    histos.add("h3PhiInvMassMixed", "Invariant mass of Phi meson Mixed", kTH3F, {{2001, -0.5, 2000.5}, {100, 0.0f, 10.0f}, {300, 0.9, 1.2}});
  }

  // variables
  double massKa = o2::analysis::pdg::MassKPlus;
  double rapidity, mass, pT, paircharge;
  TLorentzVector daughter1, daughter2, mother;

  // Track selection
  template <typename TrackType>
  bool trackCut(const TrackType track)
  {
    // basic track cuts
    if (std::abs(track.pt()) < cMinPtcut)
      return false;
    if (std::abs(track.dcaXY()) > cMaxDCArToPVcut)
      return false;
    if (std::abs(track.dcaZ()) > cMaxDCAzToPVcut)
      return false;
    if (cfgPrimaryTrack && !track.isPrimaryTrack())
      return false;
    if (cfgGlobalWoDCATrack && !track.isGlobalTrackWoDCA())
      return false;
    if (cfgPVContributor && !track.isPVContributor())
      return false;
    return true;
  }

  // Fill histograms (main function)
  template <bool IsMC, bool IsMix, typename CollisionType, typename TracksType>
  void fillHistograms(const CollisionType& collision, const TracksType& dTracks1, const TracksType& dTracks2)
  {
    for (const auto& trk : dTracks1) {
      if (!trackCut(trk))
        continue;
      histos.fill(HIST("InputTracks"), trk.pt() * trk.sign());
    }
  }

  // PID selection TPC +TOF Veto
  template <typename T>
  bool selectionPID(const T& candidate)
  {
    if (candidate.hasTOF() && (candidate.tofNSigmaKa() * candidate.tofNSigmaKa() + candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa()) < (2.0 * nsigmaCutCombined * nsigmaCutCombined)) {
      return true;
    } else if (std::abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC) {
      return true;
    }
    return false;
  }

  // Invariant mass distribution for different type of background
  template <typename T1, typename T2>
  void FillinvMass(const T1& candidate1, const T2& candidate2, float multiplicity, bool unlike, bool mix, bool likesign, float massd1, float massd2)
  {
    daughter1.SetXYZM(candidate1.px(), candidate1.py(), candidate1.pz(), massKa);
    daughter2.SetXYZM(candidate2.px(), candidate2.py(), candidate2.pz(), massKa);
    mother = daughter1 + daughter2;
    mass = mother.M();
    pT = mother.Pt();
    rapidity = mother.Rapidity();
    paircharge = candidate1.sign() * candidate2.sign();
    if (std::abs(rapidity) < 0.5 && paircharge < 0 && unlike) {
      histos.fill(HIST("h3PhiInvMassUnlikeSign"), multiplicity, pT, mass);
    } else if (std::abs(rapidity) < 0.5 && paircharge < 0 && mix) {
      histos.fill(HIST("h3PhiInvMassMixed"), multiplicity, pT, mass);
    } else if (std::abs(rapidity) < 0.5 && paircharge > 0 && likesign) {
      if (candidate1.sign() > 0 && candidate2.sign() > 0) {
        histos.fill(HIST("h3PhiInvMassLikeSignPP"), multiplicity, pT, mass);
      } else {
        histos.fill(HIST("h3PhiInvMassLikeSignMM"), multiplicity, pT, mass);
      }
    }
  }

  // Process the data in same event
  void process(aod::ResoCollision& collision, aod::ResoTracks const& resotracks)
  {
    // Fill the event counter
    histos.fill(HIST("hVertexZ"), collision.posZ());
    fillHistograms<false, false>(collision, resotracks, resotracks); // Fill histograms, no MC, no mixing
    auto multiplicity = collision.multV0M();
    histos.fill(HIST("hMultiplicityPercent"), multiplicity);
    for (auto track1 : resotracks) {
      if (!trackCut(track1) || !selectionPID(track1)) {
        continue;
      }
      histos.fill(HIST("hEta"), track1.eta());
      histos.fill(HIST("hDcaxy"), track1.dcaXY());
      histos.fill(HIST("hDcaz"), track1.dcaZ());
      histos.fill(HIST("hNsigmaKaonTPC"), track1.tpcNSigmaKa());
      if (track1.hasTOF()) {
        histos.fill(HIST("hNsigmaKaonTOF"), track1.tofNSigmaKa());
      }
      auto track1ID = track1.index();
      for (auto track2 : resotracks) {
        if (!trackCut(track2) || !selectionPID(track2)) {
          continue;
        }
        auto track2ID = track2.index();
        if (track2ID <= track1ID) {
          continue; // condition to avoid double counting of pair
        }
        bool unlike = true;
        bool mix = false;
        bool likesign = true;
        FillinvMass(track1, track2, multiplicity, unlike, mix, likesign, massKa, massKa);
      }
    }
  }

  // Processing Event Mixing
  using BinningTypeVtxZT0M = ColumnBinningPolicy<aod::collision::PosZ, aod::resocollision::MultV0M>;
  SliceCache cache;
  void processME(o2::aod::ResoCollisions& collisions, aod::ResoTracks const& resotracks)
  {
    auto tracksTuple = std::make_tuple(resotracks);
    BinningTypeVtxZT0M colBinning{{CfgVtxBins, CfgMultBins}, true};
    SameKindPair<aod::ResoCollisions, aod::ResoTracks, BinningTypeVtxZT0M> pairs{colBinning, nEvtMixing, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip

    for (auto& [collision1, tracks1, collision2, tracks2] : pairs) {
      auto multiplicity = collision1.multV0M();
      for (auto& [t1, t2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
        bool unlike = false;
        bool mix = true;
        bool likesign = false;
        if (!trackCut(t1) || !trackCut(t2)) {
          continue;
        }
        if (selectionPID(t1) && selectionPID(t2)) {
          FillinvMass(t1, t2, multiplicity, unlike, mix, likesign, massKa, massKa);
        }
      }
    }
  }

  PROCESS_SWITCH(resonances_tutorial, processME, "Process EventMixing for combinatorial background", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<resonances_tutorial>(cfgc)}; }
