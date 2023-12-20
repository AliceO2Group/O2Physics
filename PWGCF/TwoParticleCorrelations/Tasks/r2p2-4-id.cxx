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

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/PIDResponse.h"
#include "Framework/ASoAHelpers.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct r2p24ch {
  //-----Track&Event Selection----------------------------------------------
  Filter col = aod::evsel::sel8 == true;
  Filter collisionFilter = (nabs(aod::collision::posZ) < 10.f);
  Filter ptfilter = aod::track::pt > 0.2f && aod::track::pt < 2.0f;
  Filter globalfilter = requireGlobalTrackInFilter();
  //------------------------------------------------------------------------
  HistogramRegistry histos{"R2P2", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {
    //-----Defining Histograms------------------------------------------------
    const AxisSpec phi{36, 0, 2.0 * constants::math::PI, "phi"}, eta{24, -0.8, 0.8, "eta"}, etaphi1{864, 0, 864, "etaphi1"}, etaphi2{864, 0, 864, "etaphi2"};
    histos.add("h1d_n1_phi", "#phi distribution", kTH1D, {phi});
    histos.add("h1d_n1_eta", "#eta distribution", kTH1D, {eta});
    histos.add("h1d_n1_ptP", "p_T for +ve particles", kTH1D, {{100, 0, 6, "p_T"}});
    histos.add("h1d_n1_ptM", "p_T for -ve particles", kTH1D, {{100, 0, 6, "p_T"}});
    histos.add("h1i_n1_multPM", "Multiplicity", kTH1I, {{200, 0, 200, "Multiplicity"}});
    histos.add("h2d_n1_etaPhiP", "#rho_1 for +ve particles", kTH2D, {eta, phi});
    histos.add("h2d_n1_etaPhiM", "#rho_1 for -ve particles", kTH2D, {eta, phi});
    histos.add("h2d_n2_eta1Phi1Eta2Phi2PP", "#rho_2 for +ve+ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_n2_eta1Phi1Eta2Phi2PM", "#rho_2 for +ve-ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_n2_eta1Phi1Eta2Phi2MM", "#rho_2 for -ve-ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_pt_etaPhiP", "p_T for +ve particles", kTH2D, {eta, phi});
    histos.add("h2d_pt_etaPhiM", "p_T for -ve particles", kTH2D, {eta, phi});
    histos.add("h2d_ptpt_eta1Phi1Eta2Phi2PM", "p_Tp_T for +ve-ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_ptpt_eta1Phi1Eta2Phi2PP", "p_Tp_T for +ve+ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_ptpt_eta1Phi1Eta2Phi2MM", "p_Tp_T for -ve-ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_ptn_eta1Phi1Eta2Phi2PM", "p_Tn for +ve-ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_ptn_eta1Phi1Eta2Phi2PP", "p_Tn for +ve+ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_ptn_eta1Phi1Eta2Phi2MM", "p_Tn for -ve-ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_npt_eta1Phi1Eta2Phi2PM", "np_T for +ve-ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_npt_eta1Phi1Eta2Phi2PP", "np_T for +ve+ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_npt_eta1Phi1Eta2Phi2MM", "np_T for -ve-ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h1d_n1_pt", "p_T", kTH1D, {{100, 0, 5, "p_T"}});
    histos.add("h2d_n1_tpc1", "TPC", kTH2D, {{100, 0, 5, "p"}, {200, 0, 200, "#frac{dE}{dx}"}});
    histos.add("h2d_n1_tpc2", "TPC", kTH2D, {{100, 0, 5, "p"}, {200, 0, 200, "#frac{dE}{dx}"}});
    histos.add("h2d_n1_tof1", "TOF", kTH2D, {{100, 0, 5, "p"}, {50, 0.2, 1.2, "#beta"}});
    histos.add("h2d_n1_tof2", "TOF", kTH2D, {{100, 0, 5, "p"}, {50, 0.2, 1.2, "#beta"}});
    //------------------------------------------------------------------------
  }
  void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& filteredCollisions, soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::pidTPCPi, aod::pidTOFPi, aod::pidTPCEl, aod::pidTOFbeta, aod::TracksExtra>> const& tracks)
  {
    int mult = 0;
    int etabin1, etabin2, phibin1, phibin2;
    for (auto track1 : tracks) {
      //-----Single Particle Distribution---------------------------------------
      histos.fill(HIST("h1d_n1_pt"), track1.pt());
      histos.fill(HIST("h2d_n1_tof1"), track1.p(), track1.beta());
      histos.fill(HIST("h2d_n1_tpc1"), track1.p(), track1.tpcSignal());
      histos.fill(HIST("h2d_n1_tof2"), track1.p(), track1.beta());
      histos.fill(HIST("h2d_n1_tpc2"), track1.p(), track1.tpcSignal());

      histos.fill(HIST("h1d_n1_phi"), track1.phi());
      histos.fill(HIST("h1d_n1_eta"), track1.eta());
      if (track1.sign() == 1) //+ve particles
      {
        histos.fill(HIST("h2d_n1_etaPhiP"), track1.eta(), track1.phi());
        histos.fill(HIST("h1d_n1_ptP"), track1.pt(), 1.0 / (2.0 * constants::math::PI * track1.pt()));
        histos.fill(HIST("h2d_pt_etaPhiP"), track1.eta(), track1.phi(), track1.pt());
      } else { //-ve particles
        histos.fill(HIST("h2d_n1_etaPhiM"), track1.eta(), track1.phi());
        histos.fill(HIST("h2d_pt_etaPhiM"), track1.eta(), track1.phi(), track1.pt());
        histos.fill(HIST("h1d_n1_ptM"), track1.pt(), 1.0 / (2.0 * constants::math::PI * track1.pt()));
      }
      etabin1 = (track1.eta() + 0.8) * 15;
      phibin1 = 36 * track1.phi() / (2 * constants::math::PI);
      mult++;
      //------------------------------------------------------------------------
      for (auto track2 : tracks) {
        //-----Two Particle Distribution------------------------------------------
        if (track1.index() == track2.index())
          continue; // No auto-correlations

        etabin2 = (track2.eta() + 0.8) * 15;
        phibin2 = 36 * track2.phi() / (2 * constants::math::PI);
        if ((etabin1 >= 0) && (etabin2 >= 0) && (phibin1 >= 0) && (phibin2 >= 0) && (etabin1 < 24) && (etabin2 < 24) && (phibin1 < 36) && (phibin2 < 36)) {
          if (track1.sign() * track2.sign() == -1) // Unlike-Sign particles
          {
            histos.fill(HIST("h2d_n2_eta1Phi1Eta2Phi2PM"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5);
            histos.fill(HIST("h2d_npt_eta1Phi1Eta2Phi2PM"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track2.pt());
            histos.fill(HIST("h2d_ptn_eta1Phi1Eta2Phi2PM"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track1.pt());
            histos.fill(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PM"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track1.pt() * track2.pt());
          } else if (track1.sign() == 1) { //+ve particle pair
            histos.fill(HIST("h2d_n2_eta1Phi1Eta2Phi2PP"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5);
            histos.fill(HIST("h2d_npt_eta1Phi1Eta2Phi2PP"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track2.pt());
            histos.fill(HIST("h2d_ptn_eta1Phi1Eta2Phi2PP"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track1.pt());
            histos.fill(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PP"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track1.pt() * track2.pt());
          } else { //-ve particle pair
            histos.fill(HIST("h2d_n2_eta1Phi1Eta2Phi2MM"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5);
            histos.fill(HIST("h2d_npt_eta1Phi1Eta2Phi2MM"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track2.pt());
            histos.fill(HIST("h2d_ptn_eta1Phi1Eta2Phi2MM"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track1.pt());
            histos.fill(HIST("h2d_ptpt_eta1Phi1Eta2Phi2MM"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track1.pt() * track2.pt());
          }
        }
        //------------------------------------------------------------------------
      }
    }
    histos.fill(HIST("h1i_n1_multPM"), mult);
  }
};

struct r2p24pi {
  //-----Track&Event Selection----------------------------------------------
  Filter col = aod::evsel::sel8 == true;
  Filter collisionFilter = nabs(aod::collision::posZ) < 10.f;
  Filter ptfilter = aod::track::pt > 0.2f && aod::track::pt < 2.0f;
  Filter globalfilter = requireGlobalTrackInFilter();
  //------------------------------------------------------------------------
  HistogramRegistry histos{"R2P2", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {
    //-----Defining Histograms------------------------------------------------
    const AxisSpec phi{36, 0, 2.0 * constants::math::PI, "phi"}, eta{24, -0.8, 0.8, "eta"}, etaphi1{864, 0, 864, "etaphi1"}, etaphi2{864, 0, 864, "etaphi2"};
    histos.add("h1d_n1_phi", "#phi distribution", kTH1D, {phi});
    histos.add("h1d_n1_eta", "#eta distribution", kTH1D, {eta});
    histos.add("h1d_n1_ptP", "p_T for +ve particles", kTH1D, {{100, 0, 6, "p_T"}});
    histos.add("h1d_n1_ptM", "p_T for -ve particles", kTH1D, {{100, 0, 6, "p_T"}});
    histos.add("h1i_n1_multPM", "Multiplicity", kTH1I, {{200, 0, 200, "Multiplicity"}});
    histos.add("h2d_n1_etaPhiP", "#rho_1 for +ve particles", kTH2D, {eta, phi});
    histos.add("h2d_n1_etaPhiM", "#rho_1 for -ve particles", kTH2D, {eta, phi});
    histos.add("h2d_n2_eta1Phi1Eta2Phi2PP", "#rho_2 for +ve+ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_n2_eta1Phi1Eta2Phi2PM", "#rho_2 for +ve-ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_n2_eta1Phi1Eta2Phi2MM", "#rho_2 for -ve-ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_pt_etaPhiP", "p_T for +ve particles", kTH2D, {eta, phi});
    histos.add("h2d_pt_etaPhiM", "p_T for -ve particles", kTH2D, {eta, phi});
    histos.add("h2d_ptpt_eta1Phi1Eta2Phi2PM", "p_Tp_T for +ve-ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_ptpt_eta1Phi1Eta2Phi2PP", "p_Tp_T for +ve+ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_ptpt_eta1Phi1Eta2Phi2MM", "p_Tp_T for -ve-ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_ptn_eta1Phi1Eta2Phi2PM", "p_Tn for +ve-ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_ptn_eta1Phi1Eta2Phi2PP", "p_Tn for +ve+ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_ptn_eta1Phi1Eta2Phi2MM", "p_Tn for -ve-ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_npt_eta1Phi1Eta2Phi2PM", "np_T for +ve-ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_npt_eta1Phi1Eta2Phi2PP", "np_T for +ve+ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_npt_eta1Phi1Eta2Phi2MM", "np_T for -ve-ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h1d_n1_pt", "p_T", kTH1D, {{100, 0, 5, "p_T"}});
    histos.add("h2d_n1_tpc1", "TPC", kTH2D, {{100, 0, 5, "p"}, {200, 0, 200, "#frac{dE}{dx}"}});
    histos.add("h2d_n1_tpc2", "TPC", kTH2D, {{100, 0, 5, "p"}, {200, 0, 200, "#frac{dE}{dx}"}});
    histos.add("h2d_n1_tof1", "TOF", kTH2D, {{100, 0, 5, "p"}, {50, 0.2, 1.2, "#beta"}});
    histos.add("h2d_n1_tof2", "TOF", kTH2D, {{100, 0, 5, "p"}, {50, 0.2, 1.2, "#beta"}});
    histos.add("h2d_n1_nsigmatpc1", "N#sigma_{TPC} Before PID", kTH2D, {{100, 0, 5, "p_T"}, {25, -5, 5, "N#sigma"}});
    histos.add("h2d_n1_nsigmatpc2", "N#sigma_{TPC} After PID", kTH2D, {{100, 0, 5, "p_T"}, {25, -5, 5, "N#sigma"}});
    histos.add("h2d_n1_nsigmatof1", "N#sigma_{TOF} Before PID", kTH2D, {{100, 0, 5, "p_T"}, {25, -5, 5, "N#sigma"}});
    histos.add("h2d_n1_nsigmatof2", "N#sigma_{TOF} After PID", kTH2D, {{100, 0, 5, "p_T"}, {25, -5, 5, "N#sigma"}});
    histos.add("h1d_n1_tpcproj", "TPC Projection", kTH1D, {{25, -5, 5, "N#sigma"}});
    histos.add("h1d_n1_tofproj", "TOF Projection", kTH1D, {{25, -5, 5, "N#sigma"}});
    //------------------------------------------------------------------------
  }
  void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& filteredCollisions, soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::pidTPCPi, aod::pidTOFPi, aod::pidTOFbeta, aod::TracksExtra>> const& tracks)
  {
    int mult = 0;
    int etabin1, etabin2, phibin1, phibin2;
    float tpccut, tofcut;
    for (auto track1 : tracks) {
      //-----Single Particle Distribution---------------------------------------
      //-----PION PID----------------------------------------------------------------
      histos.fill(HIST("h1d_n1_pt"), track1.pt());
      histos.fill(HIST("h2d_n1_nsigmatpc1"), track1.p(), track1.tpcNSigmaPi());
      histos.fill(HIST("h2d_n1_nsigmatof1"), track1.p(), track1.tofNSigmaPi());
      histos.fill(HIST("h2d_n1_tof1"), track1.p(), track1.beta());
      histos.fill(HIST("h2d_n1_tpc1"), track1.p(), track1.tpcSignal());
      tpccut = 2.5;
      tofcut = 2.5;
      if (fabs(track1.tpcNSigmaPi()) >= tpccut)
        continue;
      if (track1.pt() >= 0.6) {
        if (track1.hasTOF()) {
          if (fabs(track1.tofNSigmaPi()) >= tofcut)
            continue;
        } else {
          continue;
        }
      }
      histos.fill(HIST("h2d_n1_tof2"), track1.p(), track1.beta());
      histos.fill(HIST("h2d_n1_nsigmatpc2"), track1.p(), track1.tpcNSigmaPi());
      histos.fill(HIST("h1d_n1_tpcproj"), track1.tpcNSigmaPi());
      histos.fill(HIST("h2d_n1_nsigmatof2"), track1.p(), track1.tofNSigmaPi());
      histos.fill(HIST("h1d_n1_tofproj"), track1.tofNSigmaPi());
      histos.fill(HIST("h2d_n1_tpc2"), track1.p(), track1.tpcSignal());
      //------------------------------------------------------------------------
      histos.fill(HIST("h1d_n1_phi"), track1.phi());
      histos.fill(HIST("h1d_n1_eta"), track1.eta());
      if (track1.sign() == 1) //+ve particles
      {
        histos.fill(HIST("h2d_n1_etaPhiP"), track1.eta(), track1.phi());
        histos.fill(HIST("h1d_n1_ptP"), track1.pt(), 1.0 / (2.0 * constants::math::PI * track1.pt()));
        histos.fill(HIST("h2d_pt_etaPhiP"), track1.eta(), track1.phi(), track1.pt());
      } else { //-ve particles
        histos.fill(HIST("h2d_n1_etaPhiM"), track1.eta(), track1.phi());
        histos.fill(HIST("h2d_pt_etaPhiM"), track1.eta(), track1.phi(), track1.pt());
        histos.fill(HIST("h1d_n1_ptM"), track1.pt(), 1.0 / (2.0 * constants::math::PI * track1.pt()));
      }
      etabin1 = (track1.eta() + 0.8) * 15;
      phibin1 = 36 * track1.phi() / (2 * constants::math::PI);
      mult++;
      //------------------------------------------------------------------------
      for (auto track2 : tracks) {
        //-----Two Particle Distribution------------------------------------------
        if (track1.index() == track2.index())
          continue; // No auto-correlations
        //-----PION PID-----------------------------------------------------------
        tpccut = 2.5;
        tofcut = 2.5;
        if (fabs(track2.tpcNSigmaPi()) >= tpccut)
          continue;
        if (track2.pt() >= 0.6) {
          if (track2.hasTOF()) {
            if (fabs(track2.tofNSigmaPi()) >= tofcut)
              continue;
          } else {
            continue;
          }
        }
        //------------------------------------------------------------------------
        etabin2 = (track2.eta() + 0.8) * 15;
        phibin2 = 36 * track2.phi() / (2 * constants::math::PI);
        if ((etabin1 >= 0) && (etabin2 >= 0) && (phibin1 >= 0) && (phibin2 >= 0) && (etabin1 < 24) && (etabin2 < 24) && (phibin1 < 36) && (phibin2 < 36)) {
          if (track1.sign() * track2.sign() == -1) // Unlike-sign particles
          {
            histos.fill(HIST("h2d_n2_eta1Phi1Eta2Phi2PM"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5);
            histos.fill(HIST("h2d_npt_eta1Phi1Eta2Phi2PM"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track2.pt());
            histos.fill(HIST("h2d_ptn_eta1Phi1Eta2Phi2PM"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track1.pt());
            histos.fill(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PM"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track1.pt() * track2.pt());
          } else if (track1.sign() == 1) { //+ve particle pair
            histos.fill(HIST("h2d_n2_eta1Phi1Eta2Phi2PP"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5);
            histos.fill(HIST("h2d_npt_eta1Phi1Eta2Phi2PP"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track2.pt());
            histos.fill(HIST("h2d_ptn_eta1Phi1Eta2Phi2PP"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track1.pt());
            histos.fill(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PP"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track1.pt() * track2.pt());
          } else { //-ve particle pair
            histos.fill(HIST("h2d_n2_eta1Phi1Eta2Phi2MM"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5);
            histos.fill(HIST("h2d_npt_eta1Phi1Eta2Phi2MM"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track2.pt());
            histos.fill(HIST("h2d_ptn_eta1Phi1Eta2Phi2MM"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track1.pt());
            histos.fill(HIST("h2d_ptpt_eta1Phi1Eta2Phi2MM"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track1.pt() * track2.pt());
          }
        }
        //------------------------------------------------------------------------
      }
    }
    histos.fill(HIST("h1i_n1_multPM"), mult);
  }
};

struct r2p24ka {
  //-----Track&Event Selection----------------------------------------------
  Filter col = aod::evsel::sel8 == true;
  Filter collisionFilter = nabs(aod::collision::posZ) < 10.f;
  Filter ptfilter = aod::track::pt > 0.2f && aod::track::pt < 2.0f;
  Filter globalfilter = requireGlobalTrackInFilter();
  //------------------------------------------------------------------------
  HistogramRegistry histos{"R2P2", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {
    //-----Defining Histograms------------------------------------------------
    const AxisSpec phi{36, 0, 2.0 * constants::math::PI, "phi"}, eta{24, -0.8, 0.8, "eta"}, etaphi1{864, 0, 864, "etaphi1"}, etaphi2{864, 0, 864, "etaphi2"};
    histos.add("h1d_n1_phi", "#phi distribution", kTH1D, {phi});
    histos.add("h1d_n1_eta", "#eta distribution", kTH1D, {eta});
    histos.add("h1d_n1_ptP", "p_T for +ve particles", kTH1D, {{100, 0, 6, "p_T"}});
    histos.add("h1d_n1_ptM", "p_T for -ve particles", kTH1D, {{100, 0, 6, "p_T"}});
    histos.add("h1i_n1_multPM", "Multiplicity", kTH1I, {{200, 0, 200, "Multiplicity"}});
    histos.add("h2d_n1_etaPhiP", "#rho_1 for +ve particles", kTH2D, {eta, phi});
    histos.add("h2d_n1_etaPhiM", "#rho_1 for -ve particles", kTH2D, {eta, phi});
    histos.add("h2d_n2_eta1Phi1Eta2Phi2PP", "#rho_2 for +ve+ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_n2_eta1Phi1Eta2Phi2PM", "#rho_2 for +ve-ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_n2_eta1Phi1Eta2Phi2MM", "#rho_2 for -ve-ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_pt_etaPhiP", "p_T for +ve particles", kTH2D, {eta, phi});
    histos.add("h2d_pt_etaPhiM", "p_T for -ve particles", kTH2D, {eta, phi});
    histos.add("h2d_ptpt_eta1Phi1Eta2Phi2PM", "p_Tp_T for +ve-ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_ptpt_eta1Phi1Eta2Phi2PP", "p_Tp_T for +ve+ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_ptpt_eta1Phi1Eta2Phi2MM", "p_Tp_T for -ve-ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_ptn_eta1Phi1Eta2Phi2PM", "p_Tn for +ve-ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_ptn_eta1Phi1Eta2Phi2PP", "p_Tn for +ve+ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_ptn_eta1Phi1Eta2Phi2MM", "p_Tn for -ve-ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_npt_eta1Phi1Eta2Phi2PM", "np_T for +ve-ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_npt_eta1Phi1Eta2Phi2PP", "np_T for +ve+ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_npt_eta1Phi1Eta2Phi2MM", "np_T for -ve-ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h1d_n1_pt", "p_T", kTH1D, {{100, 0, 5, "p_T"}});
    histos.add("h2d_n1_tpc1", "TPC", kTH2D, {{100, 0, 5, "p"}, {200, 0, 200, "#frac{dE}{dx}"}});
    histos.add("h2d_n1_tpc2", "TPC", kTH2D, {{100, 0, 5, "p"}, {200, 0, 200, "#frac{dE}{dx}"}});
    histos.add("h2d_n1_tof1", "TOF", kTH2D, {{100, 0, 5, "p"}, {50, 0.2, 1.2, "#beta"}});
    histos.add("h2d_n1_tof2", "TOF", kTH2D, {{100, 0, 5, "p"}, {50, 0.2, 1.2, "#beta"}});
    histos.add("h2d_n1_nsigmatpc1", "N#sigma_{TPC} Before PID", kTH2D, {{100, 0, 5, "p_T"}, {25, -5, 5, "N#sigma"}});
    histos.add("h2d_n1_nsigmatpc2", "N#sigma_{TPC} After PID", kTH2D, {{100, 0, 5, "p_T"}, {25, -5, 5, "N#sigma"}});
    histos.add("h2d_n1_nsigmatof1", "N#sigma_{TOF} Before PID", kTH2D, {{100, 0, 5, "p_T"}, {25, -5, 5, "N#sigma"}});
    histos.add("h2d_n1_nsigmatof2", "N#sigma_{TOF} After PID", kTH2D, {{100, 0, 5, "p_T"}, {25, -5, 5, "N#sigma"}});
    histos.add("h1d_n1_tpcproj", "TPC Projection", kTH1D, {{25, -5, 5, "N#sigma"}});
    histos.add("h1d_n1_tofproj", "TOF Projection", kTH1D, {{25, -5, 5, "N#sigma"}});
    //------------------------------------------------------------------------
  }
  void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& filteredCollisions, soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::pidTPCKa, aod::pidTOFKa, aod::pidTPCEl, aod::pidTOFbeta, aod::TracksExtra>> const& tracks)
  {
    int mult = 0;
    int etabin1, etabin2, phibin1, phibin2;
    float tpccut, tofcut;
    for (auto track1 : tracks) {
      //-----Single Particle Distribution--------------------------------------------
      histos.fill(HIST("h1d_n1_pt"), track1.pt());
      histos.fill(HIST("h2d_n1_nsigmatpc1"), track1.pt(), track1.tpcNSigmaKa());
      histos.fill(HIST("h2d_n1_nsigmatof1"), track1.pt(), track1.tofNSigmaKa());
      histos.fill(HIST("h2d_n1_tof1"), track1.p(), track1.beta());
      histos.fill(HIST("h2d_n1_tpc1"), track1.p(), track1.tpcSignal());
      //-----KAON PID-----------------------------------------------------------
      tpccut = 2;
      tofcut = 2;
      if (track1.pt() < 0.6) {
        if (track1.pt() < 0.45)
          tpccut = 2;
        else if (track1.pt() < 0.55)
          tpccut = 1;
        else
          tpccut = 0.6;
        if (fabs(track1.tpcNSigmaKa()) > tpccut)
          continue;
      } else if (track1.hasTOF()) {
        if ((fabs(track1.tpcNSigmaKa()) > tpccut) || (fabs(track1.tofNSigmaKa()) > tofcut))
          continue;
      } else {
        continue;
      }

      histos.fill(HIST("h2d_n1_tof2"), track1.p(), track1.beta());
      histos.fill(HIST("h2d_n1_nsigmatpc2"), track1.pt(), track1.tpcNSigmaKa());
      histos.fill(HIST("h1d_n1_tpcproj"), track1.tpcNSigmaKa());
      histos.fill(HIST("h2d_n1_nsigmatof2"), track1.pt(), track1.tofNSigmaKa());
      histos.fill(HIST("h1d_n1_tofproj"), track1.tofNSigmaKa());
      histos.fill(HIST("h2d_n1_tpc2"), track1.p(), track1.tpcSignal());

      //------------------------------------------------------------------------
      histos.fill(HIST("h1d_n1_phi"), track1.phi());
      histos.fill(HIST("h1d_n1_eta"), track1.eta());
      if (track1.sign() == 1) //+ve particle
      {
        histos.fill(HIST("h2d_n1_etaPhiP"), track1.eta(), track1.phi());
        histos.fill(HIST("h1d_n1_ptP"), track1.pt(), 1.0 / (2.0 * constants::math::PI * track1.pt()));
        histos.fill(HIST("h2d_pt_etaPhiP"), track1.eta(), track1.phi(), track1.pt());
      } else { //-ve particle
        histos.fill(HIST("h2d_n1_etaPhiM"), track1.eta(), track1.phi());
        histos.fill(HIST("h2d_pt_etaPhiM"), track1.eta(), track1.phi(), track1.pt());
        histos.fill(HIST("h1d_n1_ptM"), track1.pt(), 1.0 / (2.0 * constants::math::PI * track1.pt()));
      }

      etabin1 = (track1.eta() + 0.8) * 15;
      phibin1 = 36 * track1.phi() / (2 * constants::math::PI);
      mult++;
      //------------------------------------------------------------------------
      for (auto track2 : tracks) {
        //-----Two Particle Distribution------------------------------------------
        if (track1.index() == track2.index())
          continue; // No auto-correlations
        //-----KAON PID-----------------------------------------------------------
        tpccut = 2;
        tofcut = 2;
        if (track2.pt() < 0.6) {
          if (track2.pt() < 0.45)
            tpccut = 2;
          else if (track2.pt() < 0.55)
            tpccut = 1;
          else
            tpccut = 0.6;
          if (fabs(track2.tpcNSigmaKa()) > tpccut)
            continue;
        } else if (track2.hasTOF()) {
          if ((fabs(track2.tpcNSigmaKa()) > tpccut) || (fabs(track2.tofNSigmaKa()) > tofcut))
            continue;
        } else {
          continue;
        }
        //------------------------------------------------------------------------

        etabin2 = (track2.eta() + 0.8) * 15;
        phibin2 = 36 * track2.phi() / (2 * constants::math::PI);
        if ((etabin1 >= 0) && (etabin2 >= 0) && (phibin1 >= 0) && (phibin2 >= 0) && (etabin1 < 24) && (etabin2 < 24) && (phibin1 < 36) && (phibin2 < 36)) {
          if (track1.sign() * track2.sign() == -1) // Unlike-sign particles
          {
            histos.fill(HIST("h2d_n2_eta1Phi1Eta2Phi2PM"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5);
            histos.fill(HIST("h2d_npt_eta1Phi1Eta2Phi2PM"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track2.pt());
            histos.fill(HIST("h2d_ptn_eta1Phi1Eta2Phi2PM"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track1.pt());
            histos.fill(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PM"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track1.pt() * track2.pt());
          } else if (track1.sign() == 1) { //+ve particle pair
            histos.fill(HIST("h2d_n2_eta1Phi1Eta2Phi2PP"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5);
            histos.fill(HIST("h2d_npt_eta1Phi1Eta2Phi2PP"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track2.pt());
            histos.fill(HIST("h2d_ptn_eta1Phi1Eta2Phi2PP"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track1.pt());
            histos.fill(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PP"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track1.pt() * track2.pt());
          } else { //-ve particle pair
            histos.fill(HIST("h2d_n2_eta1Phi1Eta2Phi2MM"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5);
            histos.fill(HIST("h2d_npt_eta1Phi1Eta2Phi2MM"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track2.pt());
            histos.fill(HIST("h2d_ptn_eta1Phi1Eta2Phi2MM"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track1.pt());
            histos.fill(HIST("h2d_ptpt_eta1Phi1Eta2Phi2MM"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track1.pt() * track2.pt());
          }
        }
        //------------------------------------------------------------------------
      }
    }
    histos.fill(HIST("h1i_n1_multPM"), mult);
  }
};

struct r2p24pr {
  //-----Track&Event Selection----------------------------------------------
  Filter col = aod::evsel::sel8 == true;
  Filter collisionFilter = nabs(aod::collision::posZ) < 10.f;
  Filter ptfilter = aod::track::pt > 0.2f && aod::track::pt < 2.0f;
  Filter globalfilter = requireGlobalTrackInFilter();
  //------------------------------------------------------------------------
  HistogramRegistry histos{"R2P2", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {
    //-----Defining Histograms------------------------------------------------
    const AxisSpec phi{36, 0, 2.0 * constants::math::PI, "phi"}, eta{24, -0.8, 0.8, "eta"}, etaphi1{864, 0, 864, "etaphi1"}, etaphi2{864, 0, 864, "etaphi2"};
    histos.add("h1d_n1_phi", "#phi distribution", kTH1D, {phi});
    histos.add("h1d_n1_eta", "#eta distribution", kTH1D, {eta});
    histos.add("h1d_n1_ptP", "p_T for +ve particles", kTH1D, {{100, 0, 6, "p_T"}});
    histos.add("h1d_n1_ptM", "p_T for -ve particles", kTH1D, {{100, 0, 6, "p_T"}});
    histos.add("h1i_n1_multPM", "Multiplicity", kTH1I, {{200, 0, 200, "Multiplicity"}});
    histos.add("h2d_n1_etaPhiP", "#rho_1 for +ve particles", kTH2D, {eta, phi});
    histos.add("h2d_n1_etaPhiM", "#rho_1 for -ve particles", kTH2D, {eta, phi});
    histos.add("h2d_n2_eta1Phi1Eta2Phi2PP", "#rho_2 for +ve+ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_n2_eta1Phi1Eta2Phi2PM", "#rho_2 for +ve-ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_n2_eta1Phi1Eta2Phi2MM", "#rho_2 for -ve-ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_pt_etaPhiP", "p_T for +ve particles", kTH2D, {eta, phi});
    histos.add("h2d_pt_etaPhiM", "p_T for -ve particles", kTH2D, {eta, phi});
    histos.add("h2d_ptpt_eta1Phi1Eta2Phi2PM", "p_Tp_T for +ve-ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_ptpt_eta1Phi1Eta2Phi2PP", "p_Tp_T for +ve+ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_ptpt_eta1Phi1Eta2Phi2MM", "p_Tp_T for -ve-ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_ptn_eta1Phi1Eta2Phi2PM", "p_Tn for +ve-ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_ptn_eta1Phi1Eta2Phi2PP", "p_Tn for +ve+ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_ptn_eta1Phi1Eta2Phi2MM", "p_Tn for -ve-ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_npt_eta1Phi1Eta2Phi2PM", "np_T for +ve-ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_npt_eta1Phi1Eta2Phi2PP", "np_T for +ve+ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_npt_eta1Phi1Eta2Phi2MM", "np_T for -ve-ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h1d_n1_pt", "p_T", kTH1D, {{100, 0, 5, "p_T"}});
    histos.add("h2d_n1_tpc1", "TPC", kTH2D, {{100, 0, 5, "p"}, {200, 0, 200, "#frac{dE}{dx}"}});
    histos.add("h2d_n1_tpc2", "TPC", kTH2D, {{100, 0, 5, "p"}, {200, 0, 200, "#frac{dE}{dx}"}});
    histos.add("h2d_n1_tof1", "TOF", kTH2D, {{100, 0, 5, "p"}, {50, 0.2, 1.2, "#beta"}});
    histos.add("h2d_n1_tof2", "TOF", kTH2D, {{100, 0, 5, "p"}, {50, 0.2, 1.2, "#beta"}});
    histos.add("h2d_n1_nsigmatpc1", "N#sigma_{TPC} Before PID", kTH2D, {{100, 0, 5, "p_T"}, {25, -5, 5, "N#sigma"}});
    histos.add("h2d_n1_nsigmatpc2", "N#sigma_{TPC} After PID", kTH2D, {{100, 0, 5, "p_T"}, {25, -5, 5, "N#sigma"}});
    histos.add("h2d_n1_nsigmatof1", "N#sigma_{TOF} Before PID", kTH2D, {{100, 0, 5, "p_T"}, {25, -5, 5, "N#sigma"}});
    histos.add("h2d_n1_nsigmatof2", "N#sigma_{TOF} After PID", kTH2D, {{100, 0, 5, "p_T"}, {25, -5, 5, "N#sigma"}});
    histos.add("h1d_n1_tpcproj", "TPC Projection", kTH1D, {{25, -5, 5, "N#sigma"}});
    histos.add("h1d_n1_tofproj", "TOF Projection", kTH1D, {{25, -5, 5, "N#sigma"}});
    //------------------------------------------------------------------------
  }
  void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& filteredCollisions, soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::pidTPCPr, aod::pidTOFPr, aod::pidTPCEl, aod::pidTOFbeta, aod::TracksExtra>> const& tracks)
  {
    int mult = 0;
    int etabin1, etabin2, phibin1, phibin2;
    float tpccut, tofcut;
    for (auto track1 : tracks) {
      //-----Single Particle Distribution---------------------------------------
      //-----PROTON PID---------------------------------------------------------
      histos.fill(HIST("h1d_n1_pt"), track1.pt());
      histos.fill(HIST("h2d_n1_nsigmatpc1"), track1.pt(), track1.tpcNSigmaPr());
      histos.fill(HIST("h2d_n1_nsigmatof1"), track1.pt(), track1.tofNSigmaPr());
      histos.fill(HIST("h2d_n1_tof1"), track1.p(), track1.beta());
      histos.fill(HIST("h2d_n1_tpc1"), track1.p(), track1.tpcSignal());

      tpccut = 2.2;
      tofcut = 2;
      if (track1.pt() < 1.1) {
        if (track1.pt() < 0.85)
          tpccut = 2.2;
        else
          tpccut = 1;
        if (fabs(track1.tpcNSigmaPr()) > tpccut)
          continue;
      } else if (track1.hasTOF()) {
        if ((fabs(track1.tpcNSigmaPr()) > tpccut) || (fabs(track1.tofNSigmaPr()) > tofcut))
          continue;
      } else {
        continue;
      }

      histos.fill(HIST("h2d_n1_tof2"), track1.p(), track1.beta());
      histos.fill(HIST("h2d_n1_nsigmatpc2"), track1.pt(), track1.tpcNSigmaPr());
      histos.fill(HIST("h1d_n1_tpcproj"), track1.tpcNSigmaPr());
      histos.fill(HIST("h2d_n1_nsigmatof2"), track1.pt(), track1.tofNSigmaPr());
      histos.fill(HIST("h1d_n1_tofproj"), track1.tofNSigmaPr());
      histos.fill(HIST("h2d_n1_tpc2"), track1.p(), track1.tpcSignal());

      //------------------------------------------------------------------------
      histos.fill(HIST("h1d_n1_phi"), track1.phi());
      histos.fill(HIST("h1d_n1_eta"), track1.eta());
      if (track1.sign() == 1) //+ve particle
      {
        histos.fill(HIST("h2d_n1_etaPhiP"), track1.eta(), track1.phi());
        histos.fill(HIST("h1d_n1_ptP"), track1.pt(), 1.0 / (2.0 * constants::math::PI * track1.pt()));
        histos.fill(HIST("h2d_pt_etaPhiP"), track1.eta(), track1.phi(), track1.pt());
      } else { //-ve particle
        histos.fill(HIST("h2d_n1_etaPhiM"), track1.eta(), track1.phi());
        histos.fill(HIST("h2d_pt_etaPhiM"), track1.eta(), track1.phi(), track1.pt());
        histos.fill(HIST("h1d_n1_ptM"), track1.pt(), 1.0 / (2.0 * constants::math::PI * track1.pt()));
      }

      etabin1 = (track1.eta() + 0.8) * 15;
      phibin1 = 36 * track1.phi() / (2 * constants::math::PI);
      mult++;
      //------------------------------------------------------------------------
      for (auto track2 : tracks) {
        //-----Two Particle Distribution------------------------------------------
        if (track1.index() == track2.index())
          continue; // No auto-correlations
        //-----PROTON PID---------------------------------------------------------
        tpccut = 2.2;
        tofcut = 2;
        if (track2.pt() < 1.1) {
          if (track2.pt() < 0.85)
            tpccut = 2.2;
          else
            tpccut = 1;
          if (fabs(track2.tpcNSigmaPr()) > tpccut)
            continue;
        } else if (track2.hasTOF()) {
          if ((fabs(track2.tpcNSigmaPr()) > tpccut) || (fabs(track2.tofNSigmaPr()) > tofcut))
            continue;
        } else {
          continue;
        }

        //------------------------------------------------------------------------

        etabin2 = (track2.eta() + 0.8) * 15;
        phibin2 = 36 * track2.phi() / (2 * constants::math::PI);
        if ((etabin1 >= 0) && (etabin2 >= 0) && (phibin1 >= 0) && (phibin2 >= 0) && (etabin1 < 24) && (etabin2 < 24) && (phibin1 < 36) && (phibin2 < 36)) {
          if (track1.sign() * track2.sign() == -1) // Unlike-sign particle pair
          {
            histos.fill(HIST("h2d_n2_eta1Phi1Eta2Phi2PM"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5);
            histos.fill(HIST("h2d_npt_eta1Phi1Eta2Phi2PM"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track2.pt());
            histos.fill(HIST("h2d_ptn_eta1Phi1Eta2Phi2PM"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track1.pt());
            histos.fill(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PM"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track1.pt() * track2.pt());
          } else if (track1.sign() == 1) { //+ve particle pair
            histos.fill(HIST("h2d_n2_eta1Phi1Eta2Phi2PP"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5);
            histos.fill(HIST("h2d_npt_eta1Phi1Eta2Phi2PP"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track2.pt());
            histos.fill(HIST("h2d_ptn_eta1Phi1Eta2Phi2PP"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track1.pt());
            histos.fill(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PP"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track1.pt() * track2.pt());
          } else { //-ve particle pair
            histos.fill(HIST("h2d_n2_eta1Phi1Eta2Phi2MM"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5);
            histos.fill(HIST("h2d_npt_eta1Phi1Eta2Phi2MM"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track2.pt());
            histos.fill(HIST("h2d_ptn_eta1Phi1Eta2Phi2MM"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track1.pt());
            histos.fill(HIST("h2d_ptpt_eta1Phi1Eta2Phi2MM"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track1.pt() * track2.pt());
          }
        }
        //------------------------------------------------------------------------
      }
    }
    histos.fill(HIST("h1i_n1_multPM"), mult);
  }
};

struct crossr2p24pik {
  //-----Track&Event Selection----------------------------------------------
  Filter col = aod::evsel::sel8 == true;
  Filter collisionFilter = nabs(aod::collision::posZ) < 10.f;
  Filter ptfilter = aod::track::pt > 0.2f && aod::track::pt < 2.0f;
  Filter globalfilter = requireGlobalTrackInFilter();
  //------------------------------------------------------------------------
  HistogramRegistry histos{"R2P2", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {
    //-----Defining Histograms------------------------------------------------
    const AxisSpec phi{36, 0, 2.0 * constants::math::PI, "phi"}, eta{24, -0.8, 0.8, "eta"}, etaphi1{864, 0, 864, "etaphi1"}, etaphi2{864, 0, 864, "etaphi2"};
    histos.add("h1d_n1_phi1", "#phi distribution Particle 1", kTH1D, {phi});
    histos.add("h1d_n1_eta1", "#eta distribution Particle 1", kTH1D, {eta});
    histos.add("h1d_n1_phi2", "#phi distribution Particle 2", kTH1D, {phi});
    histos.add("h1d_n1_eta2", "#eta distribution Particle 2", kTH1D, {eta});
    histos.add("h1d_n1_ptP1", "p_T for +ve_1", kTH1D, {{30, 0, 6, "p_T"}});
    histos.add("h1d_n1_ptM1", "p_T for -ve_1", kTH1D, {{30, 0, 6, "p_T"}});
    histos.add("h1d_n1_ptP2", "p_T for +ve_2", kTH1D, {{30, 0, 6, "p_T"}});
    histos.add("h1d_n1_ptM2", "p_T for -ve_2", kTH1D, {{30, 0, 6, "p_T"}});
    histos.add("h1i_n1_multPM", "Multiplicity", kTH1I, {{200, 0, 200, "Multiplicity"}});
    histos.add("h2d_n1_etaPhiP1", "#rho_1 for +ve particle1", kTH2D, {eta, phi});
    histos.add("h2d_n1_etaPhiM1", "#rho_1 for -ve particle1", kTH2D, {eta, phi});
    histos.add("h2d_n1_etaPhiP2", "#rho_1 for +ve particle2", kTH2D, {eta, phi});
    histos.add("h2d_n1_etaPhiM2", "#rho_1 for -ve particle2", kTH2D, {eta, phi});
    histos.add("h2d_n2_eta1Phi1Eta2Phi2PP", "#rho_2 for +ve+ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_n2_eta1Phi1Eta2Phi2PM12", "#rho_2 for +ve_1 -ve_2 particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_n2_eta1Phi1Eta2Phi2PM21", "#rho_2 for +ve_2 -ve_1 particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_n2_eta1Phi1Eta2Phi2MM", "#rho_2 for -ve-ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_pt_etaPhiP1", "p_T for +ve_1", kTH2D, {eta, phi});
    histos.add("h2d_pt_etaPhiM1", "p_T for -ve_1", kTH2D, {eta, phi});
    histos.add("h2d_pt_etaPhiP2", "p_T for +ve_2", kTH2D, {eta, phi});
    histos.add("h2d_pt_etaPhiM2", "p_T for -ve_2", kTH2D, {eta, phi});
    histos.add("h2d_ptpt_eta1Phi1Eta2Phi2PM12", "p_Tp_T for +ve_1 -ve_2", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_ptpt_eta1Phi1Eta2Phi2PM21", "p_Tp_T for +ve_2 -ve_1", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_ptpt_eta1Phi1Eta2Phi2PP", "p_Tp_T for +ve+ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_ptpt_eta1Phi1Eta2Phi2MM", "p_Tp_T for -ve-ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_ptn_eta1Phi1Eta2Phi2PM12", "p_Tn for +ve_1 -ve_2", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_ptn_eta1Phi1Eta2Phi2PM21", "p_Tn for +ve_2 -ve_1", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_ptn_eta1Phi1Eta2Phi2PP", "p_Tn for +ve+ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_ptn_eta1Phi1Eta2Phi2MM", "p_Tn for -ve-ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_npt_eta1Phi1Eta2Phi2PM12", "np_T for +ve_1 -ve_2", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_npt_eta1Phi1Eta2Phi2PM21", "np_T for +ve_1 -ve_2", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_npt_eta1Phi1Eta2Phi2PP", "np_T for +ve+ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_npt_eta1Phi1Eta2Phi2MM", "np_T for -ve-ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h1d_n1_pt", "p_T", kTH1D, {{100, 0, 5, "p_T"}});
    histos.add("h2d_n1_tpc1", "TPC", kTH2D, {{100, 0, 5, "p_T"}, {200, 0, 200, "#frac{dE}{dx}"}});
    histos.add("h2d_n1_tpc2", "TPC", kTH2D, {{100, 0, 5, "p_T"}, {200, 0, 200, "#frac{dE}{dx}"}});
    histos.add("h2d_n1_tof1", "TOF", kTH2D, {{100, 0, 5, "p_T"}, {50, 0.2, 1.2, "#beta"}});
    histos.add("h2d_n1_tof2", "TOF", kTH2D, {{100, 0, 5, "p_T"}, {50, 0.2, 1.2, "#beta"}});
    histos.add("h2d_n1_nsigmatpc1", "N#sigma_{TPC} Before PID", kTH2D, {{100, 0, 5, "p_T"}, {25, -5, 5, "N#sigma"}});
    histos.add("h2d_n1_nsigmatpc2", "N#sigma_{TPC} After PID", kTH2D, {{100, 0, 5, "p_T"}, {25, -5, 5, "N#sigma"}});
    histos.add("h2d_n1_nsigmatof1", "N#sigma_{TOF} Before PID", kTH2D, {{100, 0, 5, "p_T"}, {25, -5, 5, "N#sigma"}});
    histos.add("h2d_n1_nsigmatof2", "N#sigma_{TOF} After PID", kTH2D, {{100, 0, 5, "p_T"}, {25, -5, 5, "N#sigma"}});
    histos.add("h1d_n1_tpcproj", "TPC Projection", kTH1D, {{25, -5, 5, "N#sigma"}});
    histos.add("h1d_n1_tofproj", "TOF Projection", kTH1D, {{25, -5, 5, "N#sigma"}});
    //------------------------------------------------------------------------
  }
  void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& filteredCollisions, soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::pidTPCPi, aod::pidTOFPi, aod::pidTPCPr, aod::pidTOFPr, aod::pidTPCKa, aod::pidTOFKa, aod::pidTPCEl, aod::pidTOFbeta, aod::TracksExtra>> const& tracks)
  {
    bool flag;
    int mult = 0;
    int etabin1, etabin2, phibin1, phibin2;
    float tpccut, tofcut;
    for (auto track1 : tracks) {
      //-----Single Particle Distribution---------------------------------------
      //-----KAON PID Particle2-------------------------------------------------
      flag = true;
      tpccut = 2;
      tofcut = 2;
      if (track1.pt() < 0.6) {
        if (track1.pt() < 0.45)
          tpccut = 2;
        else if (track1.pt() < 0.55)
          tpccut = 1;
        else
          tpccut = 0.6;
        if (fabs(track1.tpcNSigmaKa()) > tpccut)
          flag = false;
      } else if (track1.hasTOF()) {
        if ((fabs(track1.tpcNSigmaKa()) > tpccut) || (fabs(track1.tofNSigmaKa()) > tofcut))
          flag = false;
      } else {
        flag = false;
      }
      //------------------------------------------------------------------------
      if (flag == true) {
        histos.fill(HIST("h2d_n1_tof1"), track1.pt(), track1.beta());
        histos.fill(HIST("h2d_n1_tpc1"), track1.pt(), track1.tpcSignal());
        histos.fill(HIST("h1d_n1_phi1"), track1.phi());
        histos.fill(HIST("h1d_n1_eta1"), track1.eta());
        if (track1.sign() == 1) //+ve particle 2
        {
          histos.fill(HIST("h2d_n1_etaPhiP2"), track1.eta(), track1.phi());
          histos.fill(HIST("h1d_n1_ptP2"), track1.pt(), 1.0 / (2.0 * constants::math::PI * track1.pt()));
          histos.fill(HIST("h2d_pt_etaPhiP2"), track1.eta(), track1.phi(), track1.pt());
        } else { //-ve particle 2
          histos.fill(HIST("h2d_n1_etaPhiM2"), track1.eta(), track1.phi());
          histos.fill(HIST("h2d_pt_etaPhiM2"), track1.eta(), track1.phi(), track1.pt());
          histos.fill(HIST("h1d_n1_ptM2"), track1.pt(), 1.0 / (2.0 * constants::math::PI * track1.pt()));
        }
      }

      //-----PION PID Particle1-------------------------------------------------
      tpccut = 2.5;
      tofcut = 2.5;
      if (fabs(track1.tpcNSigmaPi()) >= tpccut)
        continue;
      if (track1.pt() >= 0.6) {
        if (track1.hasTOF()) {
          if (fabs(track1.tofNSigmaPi()) >= tofcut)
            continue;
        } else {
          continue;
        }
      }
      histos.fill(HIST("h2d_n1_tof2"), track1.pt(), track1.beta());
      histos.fill(HIST("h2d_n1_tpc2"), track1.pt(), track1.tpcSignal());
      histos.fill(HIST("h1d_n1_phi2"), track1.phi());
      histos.fill(HIST("h1d_n1_eta2"), track1.eta());
      if (track1.sign() == 1) //+ve particle1
      {
        histos.fill(HIST("h2d_n1_etaPhiP1"), track1.eta(), track1.phi());
        histos.fill(HIST("h1d_n1_ptP1"), track1.pt(), 1.0 / (2.0 * constants::math::PI * track1.pt()));
        histos.fill(HIST("h2d_pt_etaPhiP1"), track1.eta(), track1.phi(), track1.pt());
      } else { //-ve particle1
        histos.fill(HIST("h2d_n1_etaPhiM1"), track1.eta(), track1.phi());
        histos.fill(HIST("h2d_pt_etaPhiM1"), track1.eta(), track1.phi(), track1.pt());
        histos.fill(HIST("h1d_n1_ptM1"), track1.pt(), 1.0 / (2.0 * constants::math::PI * track1.pt()));
      }
      //------------------------------------------------------------------------
      etabin1 = (track1.eta() + 0.8) * 15;
      phibin1 = 36 * track1.phi() / (2 * constants::math::PI);
      mult++;
      //------------------------------------------------------------------------
      for (auto track2 : tracks) {
        //-----Two Particle Distribution------------------------------------------
        if (track1.index() == track2.index())
          continue; // No auto-correlations
        //-----KAON PID-----------------------------------------------------------
        tpccut = 2;
        tofcut = 2;
        if (track2.pt() < 0.6) {
          if (track2.pt() < 0.45)
            tpccut = 2;
          else if (track2.pt() < 0.55)
            tpccut = 1;
          else
            tpccut = 0.6;
          if (fabs(track2.tpcNSigmaKa()) > tpccut)
            continue;
        } else if (track2.hasTOF()) {
          if ((fabs(track2.tpcNSigmaKa()) > tpccut) || (fabs(track2.tofNSigmaKa()) > tofcut))
            continue;
        } else {
          continue;
        }
        //------------------------------------------------------------------------

        etabin2 = (track2.eta() + 0.8) * 15;
        phibin2 = 36 * track2.phi() / (2 * constants::math::PI);
        if ((etabin1 >= 0) && (etabin2 >= 0) && (phibin1 >= 0) && (phibin2 >= 0) && (etabin1 < 24) && (etabin2 < 24) && (phibin1 < 36) && (phibin2 < 36)) {
          if (track1.sign() * track2.sign() == -1) // Unlike-Sign particle pair
          {
            if (track1.sign() == 1) // Particle1 +ve
            {
              histos.fill(HIST("h2d_n2_eta1Phi1Eta2Phi2PM12"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5);
              histos.fill(HIST("h2d_npt_eta1Phi1Eta2Phi2PM12"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track2.pt());
              histos.fill(HIST("h2d_ptn_eta1Phi1Eta2Phi2PM12"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track1.pt());
              histos.fill(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PM12"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track1.pt() * track2.pt());
            } else { // Particle2 +ve
              histos.fill(HIST("h2d_n2_eta1Phi1Eta2Phi2PM21"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5);
              histos.fill(HIST("h2d_npt_eta1Phi1Eta2Phi2PM21"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track2.pt());
              histos.fill(HIST("h2d_ptn_eta1Phi1Eta2Phi2PM21"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track1.pt());
              histos.fill(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PM21"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track1.pt() * track2.pt());
            }
          } else if (track1.sign() == 1) { //+ve particle pair
            histos.fill(HIST("h2d_n2_eta1Phi1Eta2Phi2PP"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5);
            histos.fill(HIST("h2d_npt_eta1Phi1Eta2Phi2PP"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track2.pt());
            histos.fill(HIST("h2d_ptn_eta1Phi1Eta2Phi2PP"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track1.pt());
            histos.fill(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PP"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track1.pt() * track2.pt());
          } else { //-ve particle pair
            histos.fill(HIST("h2d_n2_eta1Phi1Eta2Phi2MM"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5);
            histos.fill(HIST("h2d_npt_eta1Phi1Eta2Phi2MM"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track2.pt());
            histos.fill(HIST("h2d_ptn_eta1Phi1Eta2Phi2MM"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track1.pt());
            histos.fill(HIST("h2d_ptpt_eta1Phi1Eta2Phi2MM"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track1.pt() * track2.pt());
          }
        }
        //------------------------------------------------------------------------
      }
    }
    histos.fill(HIST("h1i_n1_multPM"), mult);
  }
};

struct crossr2p24pip {
  //-----Track&Event Selection----------------------------------------------
  Filter col = aod::evsel::sel8 == true;
  Filter collisionFilter = nabs(aod::collision::posZ) < 10.f;
  Filter ptfilter = aod::track::pt > 0.2f && aod::track::pt < 2.0f;
  Filter globalfilter = requireGlobalTrackInFilter();
  //------------------------------------------------------------------------
  HistogramRegistry histos{"R2P2", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {
    //-----Defining Histograms------------------------------------------------
    const AxisSpec phi{36, 0, 2.0 * constants::math::PI, "phi"}, eta{24, -0.8, 0.8, "eta"}, etaphi1{864, 0, 864, "etaphi1"}, etaphi2{864, 0, 864, "etaphi2"};
    histos.add("h1d_n1_phi1", "#phi distribution Particle 1", kTH1D, {phi});
    histos.add("h1d_n1_eta1", "#eta distribution Particle 1", kTH1D, {eta});
    histos.add("h1d_n1_phi2", "#phi distribution Particle 2", kTH1D, {phi});
    histos.add("h1d_n1_eta2", "#eta distribution Particle 2", kTH1D, {eta});
    histos.add("h1d_n1_ptP1", "p_T for +ve_1", kTH1D, {{30, 0, 6, "p_T"}});
    histos.add("h1d_n1_ptM1", "p_T for -ve_1", kTH1D, {{30, 0, 6, "p_T"}});
    histos.add("h1d_n1_ptP2", "p_T for +ve_2", kTH1D, {{30, 0, 6, "p_T"}});
    histos.add("h1d_n1_ptM2", "p_T for -ve_2", kTH1D, {{30, 0, 6, "p_T"}});
    histos.add("h1i_n1_multPM", "Multiplicity", kTH1I, {{200, 0, 200, "Multiplicity"}});
    histos.add("h2d_n1_etaPhiP1", "#rho_1 for +ve particle1", kTH2D, {eta, phi});
    histos.add("h2d_n1_etaPhiM1", "#rho_1 for -ve particle1", kTH2D, {eta, phi});
    histos.add("h2d_n1_etaPhiP2", "#rho_1 for +ve particle2", kTH2D, {eta, phi});
    histos.add("h2d_n1_etaPhiM2", "#rho_1 for -ve particle2", kTH2D, {eta, phi});
    histos.add("h2d_n2_eta1Phi1Eta2Phi2PP", "#rho_2 for +ve+ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_n2_eta1Phi1Eta2Phi2PM12", "#rho_2 for +ve_1 -ve_2 particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_n2_eta1Phi1Eta2Phi2PM21", "#rho_2 for +ve_2 -ve_1 particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_n2_eta1Phi1Eta2Phi2MM", "#rho_2 for -ve-ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_pt_etaPhiP1", "p_T for +ve_1", kTH2D, {eta, phi});
    histos.add("h2d_pt_etaPhiM1", "p_T for -ve_1", kTH2D, {eta, phi});
    histos.add("h2d_pt_etaPhiP2", "p_T for +ve_2", kTH2D, {eta, phi});
    histos.add("h2d_pt_etaPhiM2", "p_T for -ve_2", kTH2D, {eta, phi});
    histos.add("h2d_ptpt_eta1Phi1Eta2Phi2PM12", "p_Tp_T for +ve_1 -ve_2", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_ptpt_eta1Phi1Eta2Phi2PM21", "p_Tp_T for +ve_2 -ve_1", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_ptpt_eta1Phi1Eta2Phi2PP", "p_Tp_T for +ve+ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_ptpt_eta1Phi1Eta2Phi2MM", "p_Tp_T for -ve-ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_ptn_eta1Phi1Eta2Phi2PM12", "p_Tn for +ve_1 -ve_2", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_ptn_eta1Phi1Eta2Phi2PM21", "p_Tn for +ve_2 -ve_1", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_ptn_eta1Phi1Eta2Phi2PP", "p_Tn for +ve+ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_ptn_eta1Phi1Eta2Phi2MM", "p_Tn for -ve-ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_npt_eta1Phi1Eta2Phi2PM12", "np_T for +ve_1 -ve_2", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_npt_eta1Phi1Eta2Phi2PM21", "np_T for +ve_1 -ve_2", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_npt_eta1Phi1Eta2Phi2PP", "np_T for +ve+ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_npt_eta1Phi1Eta2Phi2MM", "np_T for -ve-ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h1d_n1_pt", "p_T", kTH1D, {{100, 0, 5, "p_T"}});
    histos.add("h2d_n1_tpc1", "TPC", kTH2D, {{100, 0, 5, "p_T"}, {200, 0, 200, "#frac{dE}{dx}"}});
    histos.add("h2d_n1_tpc2", "TPC", kTH2D, {{100, 0, 5, "p_T"}, {200, 0, 200, "#frac{dE}{dx}"}});
    histos.add("h2d_n1_tof1", "TOF", kTH2D, {{100, 0, 5, "p_T"}, {50, 0.2, 1.2, "#beta"}});
    histos.add("h2d_n1_tof2", "TOF", kTH2D, {{100, 0, 5, "p_T"}, {50, 0.2, 1.2, "#beta"}});
    histos.add("h2d_n1_nsigmatpc1", "N#sigma_{TPC} Before PID", kTH2D, {{100, 0, 5, "p_T"}, {25, -5, 5, "N#sigma"}});
    histos.add("h2d_n1_nsigmatpc2", "N#sigma_{TPC} After PID", kTH2D, {{100, 0, 5, "p_T"}, {25, -5, 5, "N#sigma"}});
    histos.add("h2d_n1_nsigmatof1", "N#sigma_{TOF} Before PID", kTH2D, {{100, 0, 5, "p_T"}, {25, -5, 5, "N#sigma"}});
    histos.add("h2d_n1_nsigmatof2", "N#sigma_{TOF} After PID", kTH2D, {{100, 0, 5, "p_T"}, {25, -5, 5, "N#sigma"}});
    histos.add("h1d_n1_tpcproj", "TPC Projection", kTH1D, {{25, -5, 5, "N#sigma"}});
    histos.add("h1d_n1_tofproj", "TOF Projection", kTH1D, {{25, -5, 5, "N#sigma"}});
    //------------------------------------------------------------------------
  }
  void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& filteredCollisions, soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::pidTPCPi, aod::pidTOFPi, aod::pidTPCPr, aod::pidTOFPr, aod::pidTPCKa, aod::pidTOFKa, aod::pidTPCEl, aod::pidTOFbeta, aod::TracksExtra>> const& tracks)
  {
    bool flag;
    int mult = 0;
    int etabin1, etabin2, phibin1, phibin2;
    float tpccut, tofcut;
    for (auto track1 : tracks) {
      //-----Single Particle Distribution---------------------------------------
      //-----PROTON PID Particle2-----------------------------------------------
      flag = true;
      tpccut = 2.2;
      tofcut = 2;
      if (track1.pt() < 1.1) {
        if (track1.pt() < 0.85)
          tpccut = 2.2;
        else
          tpccut = 1;
        if (fabs(track1.tpcNSigmaPr()) > tpccut)
          flag = false;
      } else if (track1.hasTOF()) {
        if ((fabs(track1.tpcNSigmaPr()) > tpccut) || (fabs(track1.tofNSigmaPr()) > tofcut))
          flag = false;
      } else {
        flag = false;
      }

      if (flag == true) {
        histos.fill(HIST("h2d_n1_tof1"), track1.pt(), track1.beta());
        histos.fill(HIST("h2d_n1_tpc1"), track1.pt(), track1.tpcSignal());
        histos.fill(HIST("h1d_n1_phi1"), track1.phi());
        histos.fill(HIST("h1d_n1_eta1"), track1.eta());
        if (track1.sign() == 1) //+ve particle2
        {
          histos.fill(HIST("h2d_n1_etaPhiP2"), track1.eta(), track1.phi());
          histos.fill(HIST("h1d_n1_ptP2"), track1.pt(), 1.0 / (2.0 * constants::math::PI * track1.pt()));
          histos.fill(HIST("h2d_pt_etaPhiP2"), track1.eta(), track1.phi(), track1.pt());
        } else { //-ve particle2
          histos.fill(HIST("h2d_n1_etaPhiM2"), track1.eta(), track1.phi());
          histos.fill(HIST("h2d_pt_etaPhiM2"), track1.eta(), track1.phi(), track1.pt());
          histos.fill(HIST("h1d_n1_ptM2"), track1.pt(), 1.0 / (2.0 * constants::math::PI * track1.pt()));
        }
      }
      //------------------------------------------------------------------------
      //-----PION PID Partcile1-------------------------------------------------
      tpccut = 2.5;
      tofcut = 2.5;
      if (fabs(track1.tpcNSigmaPi()) >= tpccut)
        continue;
      if (track1.pt() >= 0.6) {
        if (track1.hasTOF()) {
          if (fabs(track1.tofNSigmaPi()) >= tofcut)
            continue;
        } else {
          continue;
        }
      }

      histos.fill(HIST("h2d_n1_tof2"), track1.pt(), track1.beta());
      histos.fill(HIST("h2d_n1_tpc2"), track1.pt(), track1.tpcSignal());
      histos.fill(HIST("h1d_n1_phi2"), track1.phi());
      histos.fill(HIST("h1d_n1_eta2"), track1.eta());
      if (track1.sign() == 1) //+ve particle 1
      {
        histos.fill(HIST("h2d_n1_etaPhiP1"), track1.eta(), track1.phi());
        histos.fill(HIST("h1d_n1_ptP1"), track1.pt(), 1.0 / (2.0 * constants::math::PI * track1.pt()));
        histos.fill(HIST("h2d_pt_etaPhiP1"), track1.eta(), track1.phi(), track1.pt());
      } else { //-ve partcile1
        histos.fill(HIST("h2d_n1_etaPhiM1"), track1.eta(), track1.phi());
        histos.fill(HIST("h2d_pt_etaPhiM1"), track1.eta(), track1.phi(), track1.pt());
        histos.fill(HIST("h1d_n1_ptM1"), track1.pt(), 1.0 / (2.0 * constants::math::PI * track1.pt()));
      }
      //------------------------------------------------------------------------
      etabin1 = (track1.eta() + 0.8) * 15;
      phibin1 = 36 * track1.phi() / (2 * constants::math::PI);
      mult++;
      //------------------------------------------------------------------------
      for (auto track2 : tracks) {
        //-----Two Particle Distribution------------------------------------------
        if (track1.index() == track2.index())
          continue; // No auto-correlations
        //-----PROTON PID Particle2-----------------------------------------------
        tpccut = 2.2;
        tofcut = 2;
        if (track2.pt() < 1.1) {
          if (track2.pt() < 0.85)
            tpccut = 2.2;
          else
            tpccut = 1;
          if (fabs(track2.tpcNSigmaPr()) > tpccut)
            continue;
        } else if (track2.hasTOF()) {
          if ((fabs(track2.tpcNSigmaPr()) > tpccut) || (fabs(track2.tofNSigmaPr()) > tofcut))
            continue;
        } else {
          continue;
        }

        //------------------------------------------------------------------------

        etabin2 = (track2.eta() + 0.8) * 15;
        phibin2 = 36 * track2.phi() / (2 * constants::math::PI);
        if ((etabin1 >= 0) && (etabin2 >= 0) && (phibin1 >= 0) && (phibin2 >= 0) && (etabin1 < 24) && (etabin2 < 24) && (phibin1 < 36) && (phibin2 < 36)) {
          if (track1.sign() * track2.sign() == -1) // Unlike-Sign particle pair
          {
            if (track1.sign() == 1) // Particle1 +ve
            {
              histos.fill(HIST("h2d_n2_eta1Phi1Eta2Phi2PM12"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5);
              histos.fill(HIST("h2d_npt_eta1Phi1Eta2Phi2PM12"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track2.pt());
              histos.fill(HIST("h2d_ptn_eta1Phi1Eta2Phi2PM12"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track1.pt());
              histos.fill(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PM12"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track1.pt() * track2.pt());
            } else { // Particle2 +ve
              histos.fill(HIST("h2d_n2_eta1Phi1Eta2Phi2PM21"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5);
              histos.fill(HIST("h2d_npt_eta1Phi1Eta2Phi2PM21"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track2.pt());
              histos.fill(HIST("h2d_ptn_eta1Phi1Eta2Phi2PM21"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track1.pt());
              histos.fill(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PM21"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track1.pt() * track2.pt());
            }
          } else if (track1.sign() == 1) { //+ve particle pair
            histos.fill(HIST("h2d_n2_eta1Phi1Eta2Phi2PP"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5);
            histos.fill(HIST("h2d_npt_eta1Phi1Eta2Phi2PP"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track2.pt());
            histos.fill(HIST("h2d_ptn_eta1Phi1Eta2Phi2PP"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track1.pt());
            histos.fill(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PP"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track1.pt() * track2.pt());
          } else { //-ve particle pair
            histos.fill(HIST("h2d_n2_eta1Phi1Eta2Phi2MM"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5);
            histos.fill(HIST("h2d_npt_eta1Phi1Eta2Phi2MM"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track2.pt());
            histos.fill(HIST("h2d_ptn_eta1Phi1Eta2Phi2MM"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track1.pt());
            histos.fill(HIST("h2d_ptpt_eta1Phi1Eta2Phi2MM"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track1.pt() * track2.pt());
          }
        }
        //------------------------------------------------------------------------
      }
    }
    histos.fill(HIST("h1i_n1_multPM"), mult);
  }
};

struct crossr2p24pk {
  //-----Track&Event Selection----------------------------------------------
  Filter col = aod::evsel::sel8 == true;
  Filter collisionFilter = nabs(aod::collision::posZ) < 10.f;
  Filter ptfilter = aod::track::pt > 0.2f && aod::track::pt < 2.0f;
  Filter globalfilter = requireGlobalTrackInFilter();
  //------------------------------------------------------------------------
  HistogramRegistry histos{"R2P2", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {
    //-----Defining Histograms------------------------------------------------
    const AxisSpec phi{36, 0, 2.0 * constants::math::PI, "phi"}, eta{24, -0.8, 0.8, "eta"}, etaphi1{864, 0, 864, "etaphi1"}, etaphi2{864, 0, 864, "etaphi2"};
    histos.add("h1d_n1_phi1", "#phi distribution Particle 1", kTH1D, {phi});
    histos.add("h1d_n1_eta1", "#eta distribution Particle 1", kTH1D, {eta});
    histos.add("h1d_n1_phi2", "#phi distribution Particle 2", kTH1D, {phi});
    histos.add("h1d_n1_eta2", "#eta distribution Particle 2", kTH1D, {eta});
    histos.add("h1d_n1_ptP1", "p_T for +ve_1", kTH1D, {{30, 0, 6, "p_T"}});
    histos.add("h1d_n1_ptM1", "p_T for -ve_1", kTH1D, {{30, 0, 6, "p_T"}});
    histos.add("h1d_n1_ptP2", "p_T for +ve_2", kTH1D, {{30, 0, 6, "p_T"}});
    histos.add("h1d_n1_ptM2", "p_T for -ve_2", kTH1D, {{30, 0, 6, "p_T"}});
    histos.add("h1i_n1_multPM", "Multiplicity", kTH1I, {{200, 0, 200, "Multiplicity"}});
    histos.add("h2d_n1_etaPhiP1", "#rho_1 for +ve particle1", kTH2D, {eta, phi});
    histos.add("h2d_n1_etaPhiM1", "#rho_1 for -ve particle1", kTH2D, {eta, phi});
    histos.add("h2d_n1_etaPhiP2", "#rho_1 for +ve particle2", kTH2D, {eta, phi});
    histos.add("h2d_n1_etaPhiM2", "#rho_1 for -ve particle2", kTH2D, {eta, phi});
    histos.add("h2d_n2_eta1Phi1Eta2Phi2PP", "#rho_2 for +ve+ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_n2_eta1Phi1Eta2Phi2PM12", "#rho_2 for +ve_1 -ve_2 particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_n2_eta1Phi1Eta2Phi2PM21", "#rho_2 for +ve_2 -ve_1 particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_n2_eta1Phi1Eta2Phi2MM", "#rho_2 for -ve-ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_pt_etaPhiP1", "p_T for +ve_1", kTH2D, {eta, phi});
    histos.add("h2d_pt_etaPhiM1", "p_T for -ve_1", kTH2D, {eta, phi});
    histos.add("h2d_pt_etaPhiP2", "p_T for +ve_2", kTH2D, {eta, phi});
    histos.add("h2d_pt_etaPhiM2", "p_T for -ve_2", kTH2D, {eta, phi});
    histos.add("h2d_ptpt_eta1Phi1Eta2Phi2PM12", "p_Tp_T for +ve_1 -ve_2", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_ptpt_eta1Phi1Eta2Phi2PM21", "p_Tp_T for +ve_2 -ve_1", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_ptpt_eta1Phi1Eta2Phi2PP", "p_Tp_T for +ve+ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_ptpt_eta1Phi1Eta2Phi2MM", "p_Tp_T for -ve-ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_ptn_eta1Phi1Eta2Phi2PM12", "p_Tn for +ve_1 -ve_2", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_ptn_eta1Phi1Eta2Phi2PM21", "p_Tn for +ve_2 -ve_1", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_ptn_eta1Phi1Eta2Phi2PP", "p_Tn for +ve+ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_ptn_eta1Phi1Eta2Phi2MM", "p_Tn for -ve-ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_npt_eta1Phi1Eta2Phi2PM12", "np_T for +ve_1 -ve_2", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_npt_eta1Phi1Eta2Phi2PM21", "np_T for +ve_1 -ve_2", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_npt_eta1Phi1Eta2Phi2PP", "np_T for +ve+ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h2d_npt_eta1Phi1Eta2Phi2MM", "np_T for -ve-ve particles", kTH2D, {etaphi1, etaphi2});
    histos.add("h1d_n1_pt", "p_T", kTH1D, {{100, 0, 5, "p_T"}});
    histos.add("h2d_n1_tpc1", "TPC", kTH2D, {{100, 0, 5, "p_T"}, {200, 0, 200, "#frac{dE}{dx}"}});
    histos.add("h2d_n1_tpc2", "TPC", kTH2D, {{100, 0, 5, "p_T"}, {200, 0, 200, "#frac{dE}{dx}"}});
    histos.add("h2d_n1_tof1", "TOF", kTH2D, {{100, 0, 5, "p_T"}, {50, 0.2, 1.2, "#beta"}});
    histos.add("h2d_n1_tof2", "TOF", kTH2D, {{100, 0, 5, "p_T"}, {50, 0.2, 1.2, "#beta"}});
    histos.add("h2d_n1_nsigmatpc1", "N#sigma_{TPC} Before PID", kTH2D, {{100, 0, 5, "p_T"}, {25, -5, 5, "N#sigma"}});
    histos.add("h2d_n1_nsigmatpc2", "N#sigma_{TPC} After PID", kTH2D, {{100, 0, 5, "p_T"}, {25, -5, 5, "N#sigma"}});
    histos.add("h2d_n1_nsigmatof1", "N#sigma_{TOF} Before PID", kTH2D, {{100, 0, 5, "p_T"}, {25, -5, 5, "N#sigma"}});
    histos.add("h2d_n1_nsigmatof2", "N#sigma_{TOF} After PID", kTH2D, {{100, 0, 5, "p_T"}, {25, -5, 5, "N#sigma"}});
    histos.add("h1d_n1_tpcproj", "TPC Projection", kTH1D, {{25, -5, 5, "N#sigma"}});
    histos.add("h1d_n1_tofproj", "TOF Projection", kTH1D, {{25, -5, 5, "N#sigma"}});
    //------------------------------------------------------------------------
  }
  void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& filteredCollisions, soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::pidTPCPi, aod::pidTOFPi, aod::pidTPCPr, aod::pidTOFPr, aod::pidTPCKa, aod::pidTOFKa, aod::pidTPCEl, aod::pidTOFbeta, aod::TracksExtra>> const& tracks)
  {
    bool flag;
    int mult = 0;
    int etabin1, etabin2, phibin1, phibin2;
    float tpccut, tofcut;
    for (auto track1 : tracks) {
      //-----Single Particle Distribution---------------------------------------
      //-----KAON PID Particle2-------------------------------------------------
      flag = true;
      tpccut = 2;
      tofcut = 2;
      if (track1.pt() < 0.6) {
        if (track1.pt() < 0.45)
          tpccut = 2;
        else if (track1.pt() < 0.55)
          tpccut = 1;
        else
          tpccut = 0.6;
        if (fabs(track1.tpcNSigmaKa()) > tpccut)
          flag = false;
      } else if (track1.hasTOF()) {
        if ((fabs(track1.tpcNSigmaKa()) > tpccut) || (fabs(track1.tofNSigmaKa()) > tofcut))
          flag = false;
      } else {
        flag = false;
      }

      if (flag == true) {
        histos.fill(HIST("h2d_n1_tof1"), track1.pt(), track1.beta());
        histos.fill(HIST("h2d_n1_tpc1"), track1.pt(), track1.tpcSignal());
        histos.fill(HIST("h1d_n1_phi1"), track1.phi());
        histos.fill(HIST("h1d_n1_eta1"), track1.eta());
        if (track1.sign() == 1) //+ve particle2
        {
          histos.fill(HIST("h2d_n1_etaPhiP2"), track1.eta(), track1.phi());
          histos.fill(HIST("h1d_n1_ptP2"), track1.pt(), 1.0 / (2.0 * constants::math::PI * track1.pt()));
          histos.fill(HIST("h2d_pt_etaPhiP2"), track1.eta(), track1.phi(), track1.pt());
        } else { //-ve particle1
          histos.fill(HIST("h2d_n1_etaPhiM2"), track1.eta(), track1.phi());
          histos.fill(HIST("h2d_pt_etaPhiM2"), track1.eta(), track1.phi(), track1.pt());
          histos.fill(HIST("h1d_n1_ptM2"), track1.pt(), 1.0 / (2.0 * constants::math::PI * track1.pt()));
        }
      }

      //------------------------------------------------------------------------
      //-----PROTON PID Particle1-----------------------------------------------
      tpccut = 2.2;
      tofcut = 2;
      if (track1.pt() < 1.1) {
        if (track1.pt() < 0.85)
          tpccut = 2.2;
        else
          tpccut = 1;
        if (fabs(track1.tpcNSigmaPr()) > tpccut)
          continue;
      } else if (track1.hasTOF()) {
        if ((fabs(track1.tpcNSigmaPr()) > tpccut) || (fabs(track1.tofNSigmaPr()) > tofcut))
          continue;
      } else {
        continue;
      }

      histos.fill(HIST("h2d_n1_tof2"), track1.pt(), track1.beta());
      histos.fill(HIST("h2d_n1_tpc2"), track1.pt(), track1.tpcSignal());
      histos.fill(HIST("h1d_n1_phi2"), track1.phi());
      histos.fill(HIST("h1d_n1_eta2"), track1.eta());
      if (track1.sign() == 1) //+ve particle 1
      {
        histos.fill(HIST("h2d_n1_etaPhiP1"), track1.eta(), track1.phi());
        histos.fill(HIST("h1d_n1_ptP1"), track1.pt(), 1.0 / (2.0 * constants::math::PI * track1.pt()));
        histos.fill(HIST("h2d_pt_etaPhiP1"), track1.eta(), track1.phi(), track1.pt());
      } else { //-ve particle 2
        histos.fill(HIST("h2d_n1_etaPhiM1"), track1.eta(), track1.phi());
        histos.fill(HIST("h2d_pt_etaPhiM1"), track1.eta(), track1.phi(), track1.pt());
        histos.fill(HIST("h1d_n1_ptM1"), track1.pt(), 1.0 / (2.0 * constants::math::PI * track1.pt()));
      }
      //------------------------------------------------------------------------
      etabin1 = (track1.eta() + 0.8) * 15;
      phibin1 = 36 * track1.phi() / (2 * constants::math::PI);
      mult++;
      //------------------------------------------------------------------------
      for (auto track2 : tracks) {
        //-----Two Particle Distribution------------------------------------------
        if (track1.index() == track2.index())
          continue; // No auto-correlations
        //-----KAON PID Particle 2------------------------------------------------
        tpccut = 2;
        tofcut = 2;
        if (track2.pt() < 0.6) {
          if (track2.pt() < 0.45)
            tpccut = 2;
          else if (track2.pt() < 0.55)
            tpccut = 1;
          else
            tpccut = 0.6;
          if (fabs(track2.tpcNSigmaKa()) > tpccut)
            continue;
        } else if (track2.hasTOF()) {
          if ((fabs(track2.tpcNSigmaKa()) > tpccut) || (fabs(track2.tofNSigmaKa()) > tofcut))
            continue;
        } else {
          continue;
        }
        //------------------------------------------------------------------------
        etabin2 = (track2.eta() + 0.8) * 15;
        phibin2 = 36 * track2.phi() / (2 * constants::math::PI);
        if ((etabin1 >= 0) && (etabin2 >= 0) && (phibin1 >= 0) && (phibin2 >= 0) && (etabin1 < 24) && (etabin2 < 24) && (phibin1 < 36) && (phibin2 < 36)) {
          if (track1.sign() * track2.sign() == -1) // Unlike-Sign particle pair
          {
            if (track1.sign() == 1) // Particle1 +ve
            {
              histos.fill(HIST("h2d_n2_eta1Phi1Eta2Phi2PM12"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5);
              histos.fill(HIST("h2d_npt_eta1Phi1Eta2Phi2PM12"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track2.pt());
              histos.fill(HIST("h2d_ptn_eta1Phi1Eta2Phi2PM12"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track1.pt());
              histos.fill(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PM12"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track1.pt() * track2.pt());
            } else { // Particle2 +ve
              histos.fill(HIST("h2d_n2_eta1Phi1Eta2Phi2PM21"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5);
              histos.fill(HIST("h2d_npt_eta1Phi1Eta2Phi2PM21"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track2.pt());
              histos.fill(HIST("h2d_ptn_eta1Phi1Eta2Phi2PM21"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track1.pt());
              histos.fill(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PM21"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track1.pt() * track2.pt());
            }
          } else if (track1.sign() == 1) { //+ve particle pair
            histos.fill(HIST("h2d_n2_eta1Phi1Eta2Phi2PP"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5);
            histos.fill(HIST("h2d_npt_eta1Phi1Eta2Phi2PP"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track2.pt());
            histos.fill(HIST("h2d_ptn_eta1Phi1Eta2Phi2PP"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track1.pt());
            histos.fill(HIST("h2d_ptpt_eta1Phi1Eta2Phi2PP"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track1.pt() * track2.pt());
          } else { //-ve particle pair
            histos.fill(HIST("h2d_n2_eta1Phi1Eta2Phi2MM"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5);
            histos.fill(HIST("h2d_npt_eta1Phi1Eta2Phi2MM"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track2.pt());
            histos.fill(HIST("h2d_ptn_eta1Phi1Eta2Phi2MM"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track1.pt());
            histos.fill(HIST("h2d_ptpt_eta1Phi1Eta2Phi2MM"), 36 * etabin1 + phibin1 + 0.5, 36 * etabin2 + phibin2 + 0.5, track1.pt() * track2.pt());
          }
        }
        //------------------------------------------------------------------------
      }
    }
    histos.fill(HIST("h1i_n1_multPM"), mult);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<r2p24ch>(cfgc),
    adaptAnalysisTask<r2p24pi>(cfgc),
    adaptAnalysisTask<r2p24ka>(cfgc),
    adaptAnalysisTask<r2p24pr>(cfgc),
    adaptAnalysisTask<crossr2p24pik>(cfgc),
    adaptAnalysisTask<crossr2p24pip>(cfgc),
    adaptAnalysisTask<crossr2p24pk>(cfgc),
  };
}
