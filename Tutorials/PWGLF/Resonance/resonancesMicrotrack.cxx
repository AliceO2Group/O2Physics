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
/// \file resonancesMicrotrack.cxx
/// \brief Resonance microtrack tutorial
/// \author Bong-Hwi Lim <bong-hwi.lim@cern.ch>
/// \since 07/03/2025

#include <TLorentzVector.h>
#include <TPDGCode.h>

#include "CommonConstants/MathConstants.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/runDataProcessing.h"
#include "PWGLF/DataModel/LFResonanceTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::math;
// Extract STEP
// Handle resomicrotracks
struct ResonancesMicrotrack {
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Configurable for number of bins
  Configurable<int> nBins{"nBins", 100, "N bins in all histos"};
  // Configurable for min pT cut
  Configurable<double> cMinPtcut{"cMinPtcut", 0.15, "Track minium pt cut"};
  // Configurable for event plane
  Configurable<int> cfgEvtPl{"cfgEvtPl", 40500, "Configuration of three subsystems for the event plane and its resolution, 10000*RefA + 100*RefB + S, where FT0C:0, FT0A:1, FT0M:2, FV0A:3, BPos:5, BNeg:6"};

  // Track selection
  // primary track condition
  Configurable<bool> cfgPrimaryTrack{"cfgPrimaryTrack", false, "Primary track selection"};
  Configurable<bool> cfgGlobalWoDCATrack{"cfgGlobalWoDCATrack", true, "Global track selection without DCA"};
  Configurable<bool> cfgPVContributor{"cfgPVContributor", false, "PV contributor track selection"};
  // DCA Selections
  // DCAr to PV
  Configurable<double> cMaxDCArToPVcut{"cMaxDCArToPVcut", 0.5, "Track DCAr cut to PV Maximum"};
  // DCAz to PV
  Configurable<double> cMaxDCAzToPVcut{"cMaxDCAzToPVcut", 1.0, "Track DCAz cut to PV Maximum"};
  Configurable<double> cMinDCAzToPVcut{"cMinDCAzToPVcut", 0.0, "Track DCAz cut to PV Minimum"};

  // PID selection
  Configurable<float> nSigmaCutTPC{"nSigmaCutTPC", 3.0, "Value of the TPC Nsigma cut"};
  Configurable<float> nSigmaCutTOF{"nSigmaCutTOF", 3.0, "Value of the TOF Nsigma cut"};

  void init(o2::framework::InitContext&)
  {
    histos.add("hVertexZ", "hVertexZ", HistType::kTH1F, {{nBins, -15., 15.}});
    histos.add("hMultiplicityPercent", "Multiplicity Percentile", kTH1F, {{120, 0.0f, 120.0f}});
    histos.add("hEta_ResoTracks", "Eta distribution", kTH1F, {{200, -1.0f, 1.0f}});
    histos.add("hEta_ResoMicroTracks", "Eta distribution", kTH1F, {{200, -1.0f, 1.0f}});
    histos.add("hPhi_ResoTracks", "Phi distribution", kTH1F, {{200, 0, TwoPI}});
    histos.add("hPhi_ResoMicroTracks", "Phi distribution", kTH1F, {{200, 0, TwoPI}});
    histos.add("hPt_ResoTracks", "Pt distribution", kTH1F, {{150, 0.0f, 15.0f}});
    histos.add("hPt_ResoMicroTracks", "Pt distribution", kTH1F, {{150, 0.0f, 15.0f}});
    histos.add("hPx_ResoTracks", "Px distribution", kTH1F, {{200, -10.0f, 10.0f}});
    histos.add("hPx_ResoMicroTracks", "Px distribution", kTH1F, {{200, -10.0f, 10.0f}});
    histos.add("hPy_ResoTracks", "Py distribution", kTH1F, {{200, -10.0f, 10.0f}});
    histos.add("hPy_ResoMicroTracks", "Py distribution", kTH1F, {{200, -10.0f, 10.0f}});
    histos.add("hPz_ResoTracks", "Pz distribution", kTH1F, {{200, -10.0f, 10.0f}});
    histos.add("hPz_ResoMicroTracks", "Pz distribution", kTH1F, {{200, -10.0f, 10.0f}});

    histos.add("hDcaxy_ResoTracks", "Dcaxy distribution", kTH1F, {{200, 0, 2.0f}});
    histos.add("hDcaz_ResoTracks", "Dcaz distribution", kTH1F, {{200, 0, 2.0f}});
    histos.add("hDcaxy_ResoMicroTracks", "Dcaxy distribution", kTH1F, {{200, 0, 2.0f}});
    histos.add("hDcaz_ResoMicroTracks", "Dcaz distribution", kTH1F, {{200, 0, 2.0f}});

    // PID
    histos.add("hNsigmaKaonTPC_ResoMicroTracks", "NsigmaKaonTPC distribution", kTH1F, {{240, -6.0f, 6.0f}});
    histos.add("hNsigmaKaonTOF_ResoMicroTracks", "NsigmaKaonTOF distribution", kTH1F, {{240, -6.0f, 6.0f}});
    histos.add("hNsigmaKaonTPC_ResoTracks", "NsigmaKaonTPC distribution", kTH1F, {{240, -6.0f, 6.0f}});
    histos.add("hNsigmaKaonTOF_ResoTracks", "NsigmaKaonTOF distribution", kTH1F, {{240, -6.0f, 6.0f}});

    LOG(info) << "Size of the histograms in resonance tutorial with table combination:";
    histos.print();
  }

  template <bool IsResoMicrotrack, typename TrackType>
  bool trackCut(const TrackType track)
  {
    if constexpr (!IsResoMicrotrack) {
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
    } else {
      if (std::abs(track.pt()) < cMinPtcut)
        return false;
      if (o2::aod::resomicrodaughter::ResoMicroTrackSelFlag::decodeDCAxy(track.trackSelectionFlags()) > cMaxDCArToPVcut - Epsilon)
        return false;
      if (o2::aod::resomicrodaughter::ResoMicroTrackSelFlag::decodeDCAz(track.trackSelectionFlags()) > cMaxDCAzToPVcut - Epsilon)
        return false;
      if (cfgPrimaryTrack && !track.isPrimaryTrack())
        return false;
      if (cfgGlobalWoDCATrack && !track.isGlobalTrackWoDCA())
        return false;
      if (cfgPVContributor && !track.isPVContributor())
        return false;
    }
    return true;
  }

  template <bool IsResoMicrotrack, typename T>
  bool selectionPID(const T& candidate)
  {
    if constexpr (!IsResoMicrotrack) {
      bool tpcPass = std::abs(candidate.tpcNSigmaPr()) < nSigmaCutTPC;
      bool tofPass = (candidate.hasTOF()) ? std::abs(candidate.tofNSigmaPr()) < nSigmaCutTOF : true;
      if (tpcPass && tofPass) {
        return true;
      }
      // return true;
    } else {
      bool tpcPass = std::abs(o2::aod::resomicrodaughter::PidNSigma::getTPCnSigma(candidate.pidNSigmaPrFlag())) < nSigmaCutTPC + Epsilon;
      bool tofPass = candidate.hasTOF() ? std::abs(o2::aod::resomicrodaughter::PidNSigma::getTOFnSigma(candidate.pidNSigmaPrFlag())) < nSigmaCutTOF + Epsilon : true;
      if (tpcPass && tofPass) {
        return true;
      }
      // return true;
    }
    return false;
  }

  template <bool IsMC, bool IsMix, bool IsResoMicrotrack, typename CollisionType, typename TracksType>
  void fillHistograms(const CollisionType& collision, const TracksType& dTracks)
  {
    auto multiplicity = collision.cent();
    histos.fill(HIST("hVertexZ"), collision.posZ());
    histos.fill(HIST("hMultiplicityPercent"), multiplicity);
    for (auto const& track : dTracks) {
      if (!trackCut<IsResoMicrotrack>(track) || !selectionPID<IsResoMicrotrack>(track)) {
        continue;
      }

      if constexpr (!IsResoMicrotrack) { // ResoTracks
        histos.fill(HIST("hEta_ResoTracks"), track.eta());
        histos.fill(HIST("hPhi_ResoTracks"), track.phi());
        histos.fill(HIST("hPt_ResoTracks"), track.pt());
        histos.fill(HIST("hPx_ResoTracks"), track.px());
        histos.fill(HIST("hPy_ResoTracks"), track.py());
        histos.fill(HIST("hPz_ResoTracks"), track.pz());
        histos.fill(HIST("hDcaxy_ResoTracks"), track.dcaXY());
        histos.fill(HIST("hDcaz_ResoTracks"), track.dcaZ());
        histos.fill(HIST("hNsigmaKaonTPC_ResoTracks"), track.tpcNSigmaPr());
        if (track.hasTOF()) {
          histos.fill(HIST("hNsigmaKaonTOF_ResoTracks"), track.tofNSigmaPr());
        }
      } else { // ResoMicroTracks
        histos.fill(HIST("hEta_ResoMicroTracks"), track.eta());
        histos.fill(HIST("hPhi_ResoMicroTracks"), track.phi());
        histos.fill(HIST("hPt_ResoMicroTracks"), track.pt());
        histos.fill(HIST("hPx_ResoMicroTracks"), track.px());
        histos.fill(HIST("hPy_ResoMicroTracks"), track.py());
        histos.fill(HIST("hPz_ResoMicroTracks"), track.pz());
        histos.fill(HIST("hDcaxy_ResoMicroTracks"), o2::aod::resomicrodaughter::ResoMicroTrackSelFlag::decodeDCAxy(track.trackSelectionFlags()));
        histos.fill(HIST("hDcaz_ResoMicroTracks"), o2::aod::resomicrodaughter::ResoMicroTrackSelFlag::decodeDCAz(track.trackSelectionFlags()));
        histos.fill(HIST("hNsigmaKaonTPC_ResoMicroTracks"), o2::aod::resomicrodaughter::PidNSigma::getTPCnSigma(track.pidNSigmaPrFlag()));
        if (track.hasTOF()) {
          histos.fill(HIST("hNsigmaKaonTOF_ResoMicroTracks"), o2::aod::resomicrodaughter::PidNSigma::getTOFnSigma(track.pidNSigmaPrFlag()));
        }
      }
    }
  }
  void processDummy(aod::ResoCollision const& /*collisions*/)
  {
  }
  PROCESS_SWITCH(ResonancesMicrotrack, processDummy, "Process Dummy", true);

  void processResoTracks(aod::ResoCollision const& collision, aod::ResoTracks const& resotracks)
  {
    fillHistograms<false, false, false>(collision, resotracks);
  }
  PROCESS_SWITCH(ResonancesMicrotrack, processResoTracks, "Process ResoTracks", false);

  void processResoMicroTracks(aod::ResoCollision const& collision, aod::ResoMicroTracks const& resomicrotracks)
  {
    fillHistograms<false, false, true>(collision, resomicrotracks);
  }
  PROCESS_SWITCH(ResonancesMicrotrack, processResoMicroTracks, "Process ResoMicroTracks", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<ResonancesMicrotrack>(cfgc)}; }
