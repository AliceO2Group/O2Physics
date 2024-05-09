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
// \Single Gap Event Analyzer
// \author Sandeep Dudi, sandeep.dudi3@gmail.com
// \since  May 2024

#include <cstdlib>
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "iostream"
#include "PWGUD/DataModel/UDTables.h"
#include "PWGUD/Core/SGSelector.h"
#include "PWGUD/Core/SGTrackSelector.h"
#include "Common/DataModel/PIDResponse.h"
#include <TString.h>
#include "TLorentzVector.h"
using namespace std;
using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
#define mpion 0.1396
#define mkaon 0.4937
#define mproton 0.9383
#define mmuon 0.1057

struct SGResonanceAnalyzer {
  SGSelector sgSelector;
  Configurable<float> FV0_cut{"FV0", 50., "FV0A threshold"};
  Configurable<float> FT0A_cut{"FT0A", 50., "FT0A threshold"};
  Configurable<float> FT0C_cut{"FT0C", 50., "FT0C threshold"};
  Configurable<float> FDDA_cut{"FDDA", 10000., "FDDA threshold"};
  Configurable<float> FDDC_cut{"FDDC", 10000., "FDDC threshold"};
  Configurable<float> ZDC_cut{"ZDC", 0., "ZDC threshold"};

  // Track Selections
  Configurable<float> PV_cut{"PV_cut", 1.0, "Use Only PV tracks"};
  Configurable<float> dcaZ_cut{"dcaZ_cut", 2.0, "dcaZ cut"};
  Configurable<float> dcaXY_cut{"dcaXY_cut", 0.0, "dcaXY cut (0 for Pt-function)"};
  Configurable<float> tpcChi2_cut{"tpcChi2_cut", 4, "Max tpcChi2NCl"};
  Configurable<float> tpcNClsFindable_cut{"tpcNClsFindable_cut", 70, "Min tpcNClsFindable"};
  Configurable<float> itsChi2_cut{"itsChi2_cut", 36, "Max itsChi2NCl"};
  Configurable<float> eta_cut{"eta_cut", 0.9, "Track Pseudorapidity"};

  Configurable<float> nsigmatpc_cut{"nsigmatpc", 3.0, "nsigma tpc cut"};
  Configurable<float> nsigmatof_cut{"nsigmatof", 9.0, "nsigma tof cut"};

  Configurable<int> mintrack{"min_track", 1, "min track"};
  Configurable<int> maxtrack{"max_track", 50, "max track"};

  Configurable<bool> use_tof{"Use_TOF", true, "TOF PID"};
  Configurable<bool> QA{"QA", true, ""};
  Configurable<bool> phi{"phi", true, ""};
  Configurable<bool> rho{"rho", true, ""};
  Configurable<bool> kstar{"kstar", true, ""};

  HistogramRegistry registry{
    "registry",
    {
      {"GapSide", "Gap Side; Entries", {HistType::kTH1F, {{4, -1.5, 2.5}}}},
      {"TrueGapSide", "Gap Side; Entries", {HistType::kTH1F, {{4, -1.5, 2.5}}}},
      {"os_KK_pT_0", "pt kaon pair", {HistType::kTH3F, {{100, 0.0, 10.0}, {200, -10.0, 10.0}, {220, 0.9, 1.12}}}},
      {"os_KK_pT_1", "pt kaon pair", {HistType::kTH3F, {{100, 0.0, 10.0}, {200, -10.0, 10.0}, {220, 0.9, 1.12}}}},
      {"os_KK_pT_2", "pt kaon pair", {HistType::kTH3F, {{100, 0.0, 10.0}, {200, -10.0, 10.0}, {220, 0.9, 1.12}}}},
      {"os_KK_ls_pT_0", "kaon pair like sign", {HistType::kTH3F, {{100, 0.0, 10.0}, {200, -10.0, 10.0}, {220, 0.9, 1.12}}}},
      {"os_KK_ls_pT_1", "kaon pair like sign", {HistType::kTH3F, {{100, 0.0, 10.0}, {200, -10.0, 10.0}, {220, 0.9, 1.12}}}},
      {"os_KK_ls_pT_2", "kaon pair like sign", {HistType::kTH3F, {{100, 0.0, 10.0}, {200, -10.0, 10.0}, {220, 0.9, 1.12}}}},
      {"os_pp_pT_0", "pt pion pair", {HistType::kTH3F, {{100, 0.0, 10.0}, {200, -10.0, 10.0}, {350, 0.0, 3.5}}}},
      {"os_pp_pT_1", "pt pion pair", {HistType::kTH3F, {{100, 0.0, 10.0}, {200, -10.0, 10.0}, {350, 0.0, 3.5}}}},
      {"os_pp_pT_2", "pt pion pair", {HistType::kTH3F, {{100, 0.0, 10.0}, {200, -10.0, 10.0}, {350, 0.0, 3.5}}}},
      {"os_pp_ls_pT_0", "pion pair like sign", {HistType::kTH3F, {{100, 0.0, 10.0}, {200, -10.0, 10.0}, {350, 0.0, 3.5}}}},
      {"os_pp_ls_pT_1", "pion pair like sign", {HistType::kTH3F, {{100, 0.0, 10.0}, {200, -10.0, 10.0}, {350, 0.0, 3.5}}}},
      {"os_pp_ls_pT_2", "pion pair like sign", {HistType::kTH3F, {{100, 0.0, 10.0}, {200, -10.0, 10.0}, {350, 0.0, 3.5}}}},
      {"os_pk_pT_0", "pion-kaon pair", {HistType::kTH3F, {{100, 0.0, 10.0}, {200, -10.0, 10.0}, {400, 0.0, 2.0}}}},
      {"os_pk_pT_1", "pion-kaon pair", {HistType::kTH3F, {{100, 0.0, 10.0}, {200, -10.0, 10.0}, {400, 0.0, 2.0}}}},
      {"os_pk_pT_2", "pion-kaon pair", {HistType::kTH3F, {{100, 0.0, 10.0}, {200, -10.0, 10.0}, {400, 0.0, 2.0}}}},
      {"os_pk_ls_pT_0", "pion-kaon pair like sign", {HistType::kTH3F, {{100, 0.0, 10.0}, {200, -10.0, 10.0}, {400, 0.0, 2.0}}}},
      {"os_pk_ls_pT_1", "pion-kaon like sign", {HistType::kTH3F, {{100, 0.0, 10.0}, {200, -10.0, 10.0}, {400, 0.0, 2.0}}}},
      {"os_pk_ls_pT_2", "pion-kaon like sign", {HistType::kTH3F, {{100, 0.0, 10.0}, {200, -10.0, 10.0}, {400, 0.0, 2.0}}}},
      // QA plots

      {"tpc_dedx", "p vs dE/dx", {HistType::kTH2F, {{100, 0.0, 10.0}, {5000, 0.0, 5000.0}}}},
      {"tpc_dedx_kaon", "p#k dE/dx", {HistType::kTH2F, {{100, 0.0, 10.0}, {5000, 0.0, 5000.0}}}},
      {"tpc_dedx_pion", "p#pi dE/dx", {HistType::kTH2F, {{100, 0.0, 10.0}, {5000, 0.0, 5000.0}}}},
      {"tpc_dedx_kaon_1", "tpc+tof pid cut p#k dE/dx", {HistType::kTH2F, {{100, 0.0, 10.0}, {5000, 0.0, 5000.0}}}},
      {"tpc_dedx_kaon_2", "tpc+tof pid cut1 p#k dE/dx", {HistType::kTH2F, {{100, 0.0, 10.0}, {5000, 0.0, 5000.0}}}},
      {"tpc_dedx_pion_1", "tpc+tof pid cut p#pi dE/dx", {HistType::kTH2F, {{100, 0.0, 10.0}, {5000, 0.0, 5000.0}}}},
      {"tpc_nsigma_kaon", "p#k n#sigma", {HistType::kTH2F, {{100, 0.0, 10.0}, {20, -10.0, 10.0}}}},
      {"tpc_nsigma_pion", "p#pi n#sigma", {HistType::kTH2F, {{100, 0.0, 10.0}, {20, -10.0, 10.0}}}},
      {"FT0A", "T0A amplitude", {HistType::kTH1F, {{500, 0.0, 500.0}}}},
      {"FT0A_0", "T0A amplitude", {HistType::kTH1F, {{500, 0.0, 500.0}}}},
      {"FT0A_1", "T0A amplitude", {HistType::kTH1F, {{20000, 0.0, 20000.0}}}},
      {"FT0C", "T0C amplitude", {HistType::kTH1F, {{500, 0.0, 500.0}}}},
      {"FT0C_0", "T0C amplitude", {HistType::kTH1F, {{20000, 0.0, 20000.0}}}},
      {"FT0C_1", "T0C amplitude", {HistType::kTH1F, {{500, 0.0, 500.0}}}},
      {"ZDC_A", "ZDC amplitude", {HistType::kTH1F, {{2000, 0.0, 1000.0}}}},
      {"ZDC_A_0", "ZDC amplitude", {HistType::kTH1F, {{2000, 0.0, 1000.0}}}},
      {"ZDC_A_1", "ZDC amplitude", {HistType::kTH1F, {{2000, 0.0, 1000.0}}}},
      {"ZDC_C", "ZDC amplitude", {HistType::kTH1F, {{2000, 0.0, 1000.0}}}},
      {"ZDC_C_0", "ZDC amplitude", {HistType::kTH1F, {{2000, 0.0, 1000.0}}}},
      {"ZDC_C_1", "ZDC amplitude", {HistType::kTH1F, {{2000, 0.0, 1000.0}}}},
      {"V0A", "V0A amplitude", {HistType::kTH1F, {{1000, 0.0, 1000.0}}}},
      {"V0A_0", "V0A amplitude", {HistType::kTH1F, {{1000, 0.0, 1000.0}}}},
      {"V0A_1", "V0A amplitude", {HistType::kTH1F, {{1000, 0.0, 1000.0}}}},
      {"gap_mult0", "Mult", {HistType::kTH1F, {{100, 0.0, 100.0}}}},
      {"gap_mult1", "Mult", {HistType::kTH1F, {{100, 0.0, 100.0}}}},
      {"gap_mult2", "Mult", {HistType::kTH1F, {{100, 0.0, 100.0}}}},

    }};
  using udtracks = soa::Join<aod::UDTracks, aod::UDTracksExtra, aod::UDTracksPID>;
  using udtracksfull = soa::Join<aod::UDTracks, aod::UDTracksPID, aod::UDTracksExtra, aod::UDTracksFlags, aod::UDTracksDCA>;
  using UDCollisionsFull = soa::Join<aod::UDCollisions, aod::SGCollisions, aod::UDCollisionsSels, aod::UDZdcsReduced>; //
  using UDCollisionFull = UDCollisionsFull::iterator;

  void process(UDCollisionFull const& collision, udtracksfull const& tracks)
  {
    TLorentzVector v0;
    TLorentzVector v1;
    TLorentzVector v01;
    int gapSide = collision.gapSide();
    float FIT_cut[5] = {FV0_cut, FT0A_cut, FT0C_cut, FDDA_cut, FDDC_cut};
    std::vector<float> parameters = {PV_cut, dcaZ_cut, dcaXY_cut, tpcChi2_cut, tpcNClsFindable_cut, itsChi2_cut, eta_cut};
    int truegapSide = sgSelector.trueGap(collision, FIT_cut[0], FIT_cut[1], FIT_cut[2], ZDC_cut);
    registry.fill(HIST("GapSide"), gapSide);
    registry.fill(HIST("TrueGapSide"), truegapSide);
    gapSide = truegapSide;
    if (gapSide < 0 || gapSide > 2)
      return;
    Int_t mult = collision.numContrib();
    if (gapSide == 0) {
      registry.get<TH1>(HIST("gap_mult0"))->Fill(mult);
    }
    if (gapSide == 1) {
      registry.get<TH1>(HIST("gap_mult1"))->Fill(mult);
    }
    if (gapSide == 2) {
      registry.get<TH1>(HIST("gap_mult2"))->Fill(mult);
    }
    if (mult < mintrack || mult > maxtrack)
      return;
    for (auto track1 : tracks) {
      if (!trackselector(track1, parameters))
        continue;
      v0.SetXYZM(track1.px(), track1.py(), track1.pz(), o2::constants::physics::MassPionCharged);
      if (QA) {
        registry.get<TH2>(HIST("tpc_dedx"))->Fill(v0.P(), track1.tpcSignal());
        registry.get<TH2>(HIST("tpc_nsigma_kaon"))->Fill(v0.Pt(), track1.tpcNSigmaKa());
        registry.get<TH2>(HIST("tpc_nsigma_pion"))->Fill(v0.Pt(), track1.tpcNSigmaPi());
        if (std::abs(track1.tpcNSigmaKa()) < 3.0) {
          registry.get<TH2>(HIST("tpc_dedx_kaon"))->Fill(v0.P(), track1.tpcSignal());
        } else if (std::abs(track1.tpcNSigmaPi()) < 3.0) {
          registry.get<TH2>(HIST("tpc_dedx_pion"))->Fill(v0.P(), track1.tpcSignal());
        }
        if (selectionPIDKaon(track1, use_tof, nsigmatpc_cut, nsigmatof_cut)) {
          registry.get<TH2>(HIST("tpc_dedx_kaon_1"))->Fill(v0.P(), track1.tpcSignal());
        }
        if (selectionPIDKaon(track1, use_tof, nsigmatpc_cut, nsigmatof_cut) && std::abs(track1.tpcNSigmaPi()) > 3.0) {
          registry.get<TH2>(HIST("tpc_dedx_kaon_2"))->Fill(v0.P(), track1.tpcSignal());
        }
        if (selectionPIDPion(track1, use_tof, nsigmatpc_cut, nsigmatof_cut)) {
          registry.get<TH2>(HIST("tpc_dedx_pion_1"))->Fill(v0.P(), track1.tpcSignal());
        }
      }
    }
    if (QA) {
      if (gapSide == 0) {
        registry.get<TH1>(HIST("V0A_0"))->Fill(collision.totalFV0AmplitudeA());
        registry.get<TH1>(HIST("FT0A_0"))->Fill(collision.totalFT0AmplitudeA());
        registry.get<TH1>(HIST("FT0C_0"))->Fill(collision.totalFT0AmplitudeC());
        registry.get<TH1>(HIST("ZDC_A_0"))->Fill(collision.energyCommonZNA());
        registry.get<TH1>(HIST("ZDC_C_0"))->Fill(collision.energyCommonZNC());
      }
      if (gapSide == 1) {
        registry.get<TH1>(HIST("V0A_1"))->Fill(collision.totalFV0AmplitudeA());
        registry.get<TH1>(HIST("FT0A_1"))->Fill(collision.totalFT0AmplitudeA());
        registry.get<TH1>(HIST("FT0C_1"))->Fill(collision.totalFT0AmplitudeC());
        registry.get<TH1>(HIST("ZDC_A_1"))->Fill(collision.energyCommonZNA());
        registry.get<TH1>(HIST("ZDC_C_1"))->Fill(collision.energyCommonZNC());
      }
      if (gapSide == 2) {
        registry.get<TH1>(HIST("V0A"))->Fill(collision.totalFV0AmplitudeA());
        registry.get<TH1>(HIST("FT0A"))->Fill(collision.totalFT0AmplitudeA());
        registry.get<TH1>(HIST("FT0C"))->Fill(collision.totalFT0AmplitudeC());
        registry.get<TH1>(HIST("ZDC_A"))->Fill(collision.energyCommonZNA());
        registry.get<TH1>(HIST("ZDC_C"))->Fill(collision.energyCommonZNC());
      }
    }
    for (auto& [t0, t1] : combinations(tracks, tracks)) {
      if (!trackselector(t0, parameters) && !trackselector(t1, parameters))
        continue;
      if (phi && selectionPIDKaon(t0, use_tof, nsigmatpc_cut, nsigmatof_cut) && selectionPIDKaon(t1, use_tof, nsigmatpc_cut, nsigmatof_cut)) {
        // Apply kaon hypothesis and create pairs
        v0.SetXYZM(t0.px(), t0.py(), t0.pz(), o2::constants::physics::MassKaonCharged);
        v1.SetXYZM(t1.px(), t1.py(), t1.pz(), o2::constants::physics::MassKaonCharged);
        v01 = v0 + v1;
        // Opposite sign pairs
        if (t0.sign() != t1.sign()) {
          if (gapSide == 0) {
            registry.fill(HIST("os_KK_pT_0"), v01.Pt(), v01.Rapidity(), v01.M());
          }
          if (gapSide == 1) {
            registry.fill(HIST("os_KK_pT_1"), v01.Pt(), v01.Rapidity(), v01.M());
          }
          if (gapSide == 2) {
            registry.fill(HIST("os_KK_pT_2"), v01.Pt(), v01.Rapidity(), v01.M());
          }
        }
        // samesignpair
        if (t0.sign() == t1.sign()) {
          if (gapSide == 0) {
            registry.fill(HIST("os_KK_ls_pT_0"), v01.Pt(), v01.Rapidity(), v01.M());
          }
          if (gapSide == 1) {
            registry.fill(HIST("os_KK_ls_pT_1"), v01.Pt(), v01.Rapidity(), v01.M());
          }
          if (gapSide == 2) {
            registry.fill(HIST("os_KK_ls_pT_2"), v01.Pt(), v01.Rapidity(), v01.M());
          }
        }
      }
      if (rho && selectionPIDPion(t0, use_tof, nsigmatpc_cut, nsigmatof_cut) && selectionPIDPion(t1, use_tof, nsigmatpc_cut, nsigmatof_cut)) {
        v0.SetXYZM(t0.px(), t0.py(), t0.pz(), o2::constants::physics::MassPionCharged);
        v1.SetXYZM(t1.px(), t1.py(), t1.pz(), o2::constants::physics::MassPionCharged);
        v01 = v0 + v1;
        // Opposite sign pairs
        if (t0.sign() != t1.sign()) {
          if (gapSide == 0) {
            registry.fill(HIST("os_pp_pT_0"), v01.Pt(), v01.Rapidity(), v01.M());
          }
          if (gapSide == 1) {
            registry.fill(HIST("os_pp_pT_1"), v01.Pt(), v01.Rapidity(), v01.M());
          }
          if (gapSide == 2) {
            registry.fill(HIST("os_pp_pT_2"), v01.Pt(), v01.Rapidity(), v01.M());
          }
        } // same sign pair
        if (t0.sign() == t1.sign()) {
          if (gapSide == 0) {
            registry.fill(HIST("os_pp_ls_pT_0"), v01.Pt(), v01.Rapidity(), v01.M());
          }
          if (gapSide == 1) {
            registry.fill(HIST("os_pp_ls_pT_1"), v01.Pt(), v01.Rapidity(), v01.M());
          }
          if (gapSide == 2) {
            registry.fill(HIST("os_pp_ls_pT_2"), v01.Pt(), v01.Rapidity(), v01.M());
          }
        }
      }
      if (kstar && selectionPIDKaon(t0, use_tof, nsigmatpc_cut, nsigmatof_cut) && std::abs(t0.tpcNSigmaPi()) > 3.0 && selectionPIDPion(t1, use_tof, nsigmatpc_cut, nsigmatof_cut)) {
        v0.SetXYZM(t0.px(), t0.py(), t0.pz(), o2::constants::physics::MassKaonCharged);
        v1.SetXYZM(t1.px(), t1.py(), t1.pz(), o2::constants::physics::MassPionCharged);
        v01 = v0 + v1;
        // Opposite sign pairs
        if (t0.sign() != t1.sign()) {
          if (gapSide == 0) {
            registry.fill(HIST("os_pk_pT_0"), v01.Pt(), v01.Rapidity(), v01.M());
          }
          if (gapSide == 1) {
            registry.fill(HIST("os_pk_pT_1"), v01.Pt(), v01.Rapidity(), v01.M());
          }
          if (gapSide == 2) {
            registry.fill(HIST("os_pk_pT_2"), v01.Pt(), v01.Rapidity(), v01.M());
          }
        } // same sign pair
        if (t0.sign() == t1.sign()) {
          if (gapSide == 0) {
            registry.fill(HIST("os_pk_ls_pT_0"), v01.Pt(), v01.Rapidity(), v01.M());
          }
          if (gapSide == 1) {
            registry.fill(HIST("os_pk_ls_pT_1"), v01.Pt(), v01.Rapidity(), v01.M());
          }
          if (gapSide == 2) {
            registry.fill(HIST("os_pk_ls_pT_2"), v01.Pt(), v01.Rapidity(), v01.M());
          }
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<SGResonanceAnalyzer>(cfgc)};
}
