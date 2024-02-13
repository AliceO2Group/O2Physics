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
/// \since 08/11/2023

#include <TLorentzVector.h>

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/runDataProcessing.h"
#include "PWGLF/DataModel/LFResonanceTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// STEP 2
// Producing Mixed event invariant mass distribution
struct resonances_tutorial {
  SliceCache cache;
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
  Configurable<float> nsigmacutTOF{"nsigmacutTOF", 3.0, "Value of the TOF Nsigma cut"};

  // Event Mixing
  Configurable<int> nEvtMixing{"nEvtMixing", 5, "Number of events to mix"};
  ConfigurableAxis CfgVtxBins{"CfgVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis CfgMultBins{"CfgMultBins", {VARIABLE_WIDTH, 0., 1., 5., 10., 30., 50., 70., 100., 110.}, "Mixing bins - multiplicity"};

  // variables
  double massKa = o2::constants::physics::MassKPlus;

  // Initialize the ananlysis task
  void init(o2::framework::InitContext&)
  {
    // register histograms
    histos.add("hVertexZ", "hVertexZ", HistType::kTH1F, {{nBins, -15., 15.}});
    histos.add("hMultiplicityPercent", "Multiplicity Percentile", kTH1F, {{120, 0.0f, 120.0f}});
    histos.add("hEta", "Eta distribution", kTH1F, {{200, -1.0f, 1.0f}});
    histos.add("hDcaxy", "Dcaxy distribution", kTH1F, {{200, -1.0f, 1.0f}});
    histos.add("hDcaz", "Dcaz distribution", kTH1F, {{200, -1.0f, 1.0f}});
    histos.add("hNsigmaKaonTPC", "NsigmaKaon TPC distribution", kTH1F, {{100, -10.0f, 10.0f}});
    histos.add("hNsigmaKaonTOF", "NsigmaKaon TOF distribution", kTH1F, {{100, -10.0f, 10.0f}});
    histos.add("h1PhiInvMassUnlikeSign", "Invariant mass of Phi meson Unlike Sign", kTH1F, {{300, 0.9, 1.2}});
    histos.add("h1PhiInvMassMixed", "Invariant mass of Phi meson Mixed", kTH1F, {{300, 0.9, 1.2}});
    histos.add("h1PhiInvMassLikeSignPP", "Invariant mass of Phi meson Like Sign positive", kTH1F, {{300, 0.9, 1.2}});
    histos.add("h1PhiInvMassLikeSignMM", "Invariant mass of Phi meson Like Sign negative", kTH1F, {{300, 0.9, 1.2}});
    histos.add("h3PhiInvMassUnlikeSign", "Invariant mass of Phi meson Unlike Sign", kTH3F, {{120, 0.0f, 120.0f}, {100, 0.0f, 10.0f}, {300, 0.9, 1.2}});
    histos.add("h3PhiInvMassMixed", "Invariant mass of Phi meson Mixed", kTH3F, {{120, 0.0f, 120.0f}, {100, 0.0f, 10.0f}, {300, 0.9, 1.2}});
    histos.add("h3PhiInvMassLikeSignPP", "Invariant mass of Phi meson Like Sign positive", kTH3F, {{120, 0.0f, 120.0f}, {100, 0.0f, 10.0f}, {300, 0.9, 1.2}});
    histos.add("h3PhiInvMassLikeSignMM", "Invariant mass of Phi meson Like Sign negative", kTH3F, {{120, 0.0f, 120.0f}, {100, 0.0f, 10.0f}, {300, 0.9, 1.2}});

    // Print output histograms statistics
    LOG(info) << "Size of the histograms in resonance tutorial step2:";
    histos.print();
  }

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

  // PID selection TPC +TOF Veto
  template <typename T>
  bool selectionPID(const T& candidate)
  {
    bool tpcPass = std::abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC;
    bool tofPass = (candidate.hasTOF()) ? std::abs(candidate.tofNSigmaKa()) < nsigmacutTOF : true;
    if (tpcPass && tofPass) {
      return true;
    }
    return false;
  }

  // Fill histograms (main function)
  double rapidity, mass, pT, paircharge;
  TLorentzVector daughter1, daughter2, mother;
  template <bool IsMC, bool IsMix, typename CollisionType, typename TracksType>
  void fillHistograms(const CollisionType& collision, const TracksType& dTracks1, const TracksType& dTracks2)
  {
    auto multiplicity = collision.cent();
    for (auto track1 : dTracks1) { // loop over all dTracks1
      if (!trackCut(track1) || !selectionPID(track1)) {
        continue; // track selection and PID selection
      }
      // QA plots
      histos.fill(HIST("hEta"), track1.eta());
      histos.fill(HIST("hDcaxy"), track1.dcaXY());
      histos.fill(HIST("hDcaz"), track1.dcaZ());
      histos.fill(HIST("hNsigmaKaonTPC"), track1.tpcNSigmaKa());
      if (track1.hasTOF()) {
        histos.fill(HIST("hNsigmaKaonTOF"), track1.tofNSigmaKa());
      }
      for (auto track2 : dTracks2) { // loop over all dTracks2
        if (!trackCut(track2) || !selectionPID(track2)) {
          continue; // track selection and PID selection
        }
        if (track2.index() <= track1.index()) {
          continue; // condition to avoid double counting of pair
        }
        daughter1.SetXYZM(track1.px(), track1.py(), track1.pz(), massKa); // set the daughter1 4-momentum
        daughter2.SetXYZM(track2.px(), track2.py(), track2.pz(), massKa); // set the daughter2 4-momentum
        mother = daughter1 + daughter2;                                   // calculate the mother 4-momentum
        mass = mother.M();
        pT = mother.Pt();
        rapidity = mother.Rapidity();
        paircharge = track1.sign() * track2.sign();

        if (std::abs(rapidity) > 0.5)
          continue; // rapidity cut

        if constexpr (!IsMix) { // Same event
          if (paircharge < 0) { // unlike sign
            histos.fill(HIST("h3PhiInvMassUnlikeSign"), multiplicity, pT, mass);
            histos.fill(HIST("h1PhiInvMassUnlikeSign"), mass);
          } else {                                        // like sign
            if (track1.sign() > 0 && track2.sign() > 0) { // positive
              histos.fill(HIST("h3PhiInvMassLikeSignPP"), multiplicity, pT, mass);
              histos.fill(HIST("h1PhiInvMassLikeSignPP"), mass);
            } else { // negative
              histos.fill(HIST("h3PhiInvMassLikeSignMM"), multiplicity, pT, mass);
              histos.fill(HIST("h1PhiInvMassLikeSignMM"), mass);
            }
          }
        } else {                // Mixed event
          if (paircharge < 0) { // unlike sign
            histos.fill(HIST("h3PhiInvMassMixed"), multiplicity, pT, mass);
            histos.fill(HIST("h1PhiInvMassMixed"), mass);
          }
        }
      }
    }
  }

  // Process the data
  void process(aod::ResoCollision& collision, aod::ResoTracks const& resotracks)
  {
    // Fill the event counter
    histos.fill(HIST("hVertexZ"), collision.posZ());
    histos.fill(HIST("hMultiplicityPercent"), collision.cent());

    fillHistograms<false, false>(collision, resotracks, resotracks); // Fill histograms, no MC, no mixing
  }
  PROCESS_SWITCH(resonances_tutorial, process, "Process event for data", true); // Basic processing

  // Processing Event Mixing
  using BinningTypeVtxZT0M = ColumnBinningPolicy<aod::collision::PosZ, aod::resocollision::Cent>;
  void processME(o2::aod::ResoCollisions& collisions, aod::ResoTracks const& resotracks)
  {
    auto tracksTuple = std::make_tuple(resotracks);
    BinningTypeVtxZT0M colBinning{{CfgVtxBins, CfgMultBins}, true};
    SameKindPair<aod::ResoCollisions, aod::ResoTracks, BinningTypeVtxZT0M> pairs{colBinning, nEvtMixing, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip

    for (auto& [collision1, tracks1, collision2, tracks2] : pairs) { // loop over all pairs
      fillHistograms<false, true>(collision1, tracks1, tracks2);     // Fill histograms, no MC, mixing
    }
  };
  PROCESS_SWITCH(resonances_tutorial, processME, "Process EventMixing for combinatorial background", false); // Event Mixing
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<resonances_tutorial>(cfgc)}; }
