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
/// \file resonancesCombine.cxx
/// \brief Resonance combination tutorial
/// \author Bong-Hwi Lim <bong-hwi.lim@cern.ch>
/// \since 13/12/2024

#include <TLorentzVector.h>
#include <TPDGCode.h>

#include "CommonConstants/PhysicsConstants.h"
#include "Common/DataModel/Qvectors.h"
#include "Common/Core/EventPlaneHelper.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/runDataProcessing.h"
#include "PWGLF/DataModel/LFResonanceTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// Extract STEP
// Combine Resonance tables into other O2 tables
struct ResonanceCombine {
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Configurable for number of bins
  Configurable<int> nBins{"nBins", 100, "N bins in all histos"};
  // Configurable for min pT cut
  Configurable<double> cMinPtcut{"cMinPtcut", 0.15, "Track minium pt cut"};
  // Configurable for event plane
  Configurable<int> cfgEvtPl{"cfgEvtPl", 40500, "Configuration of three subsystems for the event plane and its resolution, 10000*RefA + 100*RefB + S, where FT0C:0, FT0A:1, FT0M:2, FV0A:3, BPos:5, BNeg:6"};

  // Track selection
  // primary track condition
  Configurable<bool> cfgPrimaryTrack{"cfgPrimaryTrack", true, "Primary track selection"};
  Configurable<bool> cfgGlobalWoDCATrack{"cfgGlobalWoDCATrack", true, "Global track selection without DCA"};
  Configurable<bool> cfgPVContributor{"cfgPVContributor", true, "PV contributor track selection"};
  // DCA Selections
  // DCAr to PV
  Configurable<double> cMaxDCArToPVcut{"cMaxDCArToPVcut", 0.5, "Track DCAr cut to PV Maximum"};
  // DCAz to PV
  Configurable<double> cMaxDCAzToPVcut{"cMaxDCAzToPVcut", 2.0, "Track DCAz cut to PV Maximum"};
  Configurable<double> cMinDCAzToPVcut{"cMinDCAzToPVcut", 0.0, "Track DCAz cut to PV Minimum"};

  // PID selection
  Configurable<float> nSigmaCutTPC{"nSigmaCutTPC", 3.0, "Value of the TPC Nsigma cut"};
  Configurable<float> nSigmaCutTOF{"nSigmaCutTOF", 3.0, "Value of the TOF Nsigma cut"};

  double massKa = o2::constants::physics::MassKPlus;

  EventPlaneHelper helperEP;
  int evtPlRefAId = static_cast<int>(cfgEvtPl / 10000);
  int evtPlRefBId = static_cast<int>((cfgEvtPl - evtPlRefAId * 10000) / 100);
  int evtPlDetId = cfgEvtPl - evtPlRefAId * 10000 - evtPlRefBId * 100;

  void init(o2::framework::InitContext&)
  {
    histos.add("hVertexZ", "hVertexZ", HistType::kTH1F, {{nBins, -15., 15.}});
    histos.add("hMultiplicityPercent", "Multiplicity Percentile", kTH1F, {{120, 0.0f, 120.0f}});
    histos.add("hEvtPl", "Event Plane", kTH1F, {{100, -1.0f, 1.0f}});
    histos.add("hEta", "Eta distribution", kTH1F, {{200, -1.0f, 1.0f}});
    histos.add("hDcaxy", "Dcaxy distribution", kTH1F, {{200, -1.0f, 1.0f}});
    histos.add("hDcaz", "Dcaz distribution", kTH1F, {{200, -1.0f, 1.0f}});
    histos.add("hNsigmaKaonTPC", "NsigmaKaon TPC distribution", kTH1F, {{100, -10.0f, 10.0f}});
    histos.add("hNsigmaKaonTOF", "NsigmaKaon TOF distribution", kTH1F, {{100, -10.0f, 10.0f}});
    histos.add("hTiemResolutionTOF", "TOF time resolution", kTH1F, {{100, -10.0f, 10.0f}});
    histos.add("h1PhiInvMassUnlikeSign", "Invariant mass of Phi meson Unlike Sign", kTH1F, {{300, 0.9, 1.2}});
    histos.add("h1PhiInvMassLikeSignPP", "Invariant mass of Phi meson Like Sign positive", kTH1F, {{300, 0.9, 1.2}});
    histos.add("h1PhiInvMassLikeSignMM", "Invariant mass of Phi meson Like Sign negative", kTH1F, {{300, 0.9, 1.2}});
    histos.add("h3PhiInvMassUnlikeSign", "Invariant mass of Phi meson Unlike Sign", kTH3F, {{120, 0.0f, 120.0f}, {100, 0.0f, 10.0f}, {300, 0.9, 1.2}});
    histos.add("h3PhiInvMassLikeSignPP", "Invariant mass of Phi meson Like Sign positive", kTH3F, {{120, 0.0f, 120.0f}, {100, 0.0f, 10.0f}, {300, 0.9, 1.2}});
    histos.add("h3PhiInvMassLikeSignMM", "Invariant mass of Phi meson Like Sign negative", kTH3F, {{120, 0.0f, 120.0f}, {100, 0.0f, 10.0f}, {300, 0.9, 1.2}});

    LOG(info) << "Size of the histograms in resonance tutorial with table combination:";
    histos.print();
  }

  template <typename TrackType>
  bool trackCut(const TrackType track)
  {
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

  template <typename T>
  bool selectionPID(const T& candidate)
  {
    bool tpcPass = std::abs(candidate.tpcNSigmaKa()) < nSigmaCutTPC;
    bool tofPass = (candidate.hasTOF()) ? std::abs(candidate.tofNSigmaKa()) < nSigmaCutTOF : true;
    if (tpcPass && tofPass) {
      return true;
    }
    return false;
  }

  double rapidity, mass, pT, paircharge;
  TLorentzVector daughter1, daughter2, mother;
  template <bool IsMC, bool IsMix, typename CollisionType, typename TracksType>
  void fillHistograms(const CollisionType& collision, const TracksType& dTracks1, const TracksType& dTracks2)
  {
    auto multiplicity = collision.template collision_as<aod::ResoCollisions>().cent();
    for (auto const& track1 : dTracks1) {
      auto track1Reso = track1.template track_as<aod::ResoTracks>();
      auto track1FullPidExt = track1.template track_as<aod::Reso2TracksPIDExt>();
      if (!trackCut(track1Reso) || !selectionPID(track1Reso)) {
        continue;
      }
      histos.fill(HIST("hEta"), track1Reso.eta());
      histos.fill(HIST("hDcaxy"), track1Reso.dcaXY());
      histos.fill(HIST("hDcaz"), track1Reso.dcaZ());
      histos.fill(HIST("hNsigmaKaonTPC"), track1Reso.tpcNSigmaKa());
      if (track1Reso.hasTOF()) {
        histos.fill(HIST("hNsigmaKaonTOF"), track1Reso.tofNSigmaKa());
        histos.fill(HIST("hTiemResolutionTOF"), track1FullPidExt.trackTimeRes()); // TOF time resolution is not in the ResoTracks table (Important)
      }
      for (auto const& track2 : dTracks2) {
        auto track2Reso = track2.template track_as<aod::ResoTracks>();
        // auto track2FullPidExt = track2.template track_as<aod::Reso2TracksPIDExt>();

        if (!trackCut(track2Reso) || !selectionPID(track2Reso)) {
          continue;
        }
        if (track2Reso.index() <= track1Reso.index()) {
          continue;
        }
        daughter1.SetXYZM(track1Reso.px(), track1Reso.py(), track1Reso.pz(), massKa);
        daughter2.SetXYZM(track2Reso.px(), track2Reso.py(), track2Reso.pz(), massKa);
        mother = daughter1 + daughter2;
        mass = mother.M();
        pT = mother.Pt();
        rapidity = mother.Rapidity();
        paircharge = track1Reso.sign() * track2Reso.sign();

        if (std::abs(rapidity) > 0.5)
          continue;

        if (paircharge < 0) {
          histos.fill(HIST("h3PhiInvMassUnlikeSign"), multiplicity, pT, mass);
          histos.fill(HIST("h1PhiInvMassUnlikeSign"), mass);
        } else {
          if (track1Reso.sign() > 0 && track2Reso.sign() > 0) {
            histos.fill(HIST("h3PhiInvMassLikeSignPP"), multiplicity, pT, mass);
            histos.fill(HIST("h1PhiInvMassLikeSignPP"), mass);
          } else {
            histos.fill(HIST("h3PhiInvMassLikeSignMM"), multiplicity, pT, mass);
            histos.fill(HIST("h1PhiInvMassLikeSignMM"), mass);
          }
        }
      }
    }
  }

  void process(soa::Join<aod::ResoCollisions, aod::Qvectors>::iterator const& collision, soa::Join<aod::ResoTracks, aod::Reso2TracksPIDExt> const& resotracks)
  {
    histos.fill(HIST("hVertexZ"), collision.posZ());
    // Both resoCollisions and Qvectors have the same cent column, so we have to use "_as" to access it
    // Similarly, we can use "_as" to access the Qvectors table or other tables
    histos.fill(HIST("hMultiplicityPercent"), collision.collision_as<aod::ResoCollisions>().cent());
    // Event plane
    auto collisionQvec = collision.template collision_as<aod::Qvectors>(); // Qvectors table is not in the ResoCollisions table (Important)
    auto evtPl = -999.0;
    if (collisionQvec.qvecAmp()[evtPlDetId] > 1e-8)
      evtPl = helperEP.GetEventPlane(collisionQvec.qvecRe()[evtPlDetId * 4 + 3], collisionQvec.qvecIm()[evtPlDetId * 4 + 3], 2);
    if (evtPl > -999.0)
      histos.fill(HIST("hEvtPl"), evtPl);

    fillHistograms<false, false>(collision, resotracks, resotracks);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<ResonanceCombine>(cfgc)}; }
