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

#include <cstdlib>
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "iostream"
#include "PWGUD/DataModel/UDTables.h"
#include "PWGUD/Core/SGSelector.h"
// #include "Common/DataModel/PIDResponse.h"
// #include "PWGUD/Core/RLhelper.h"
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

struct SGSpectraAnalyzer {
  SGSelector sgSelector;
  Configurable<float> FV0_cut{"FV0", 100., "FV0A threshold"};
  Configurable<float> FT0A_cut{"FT0A", 100., "FT0A threshold"};
  Configurable<float> FT0C_cut{"FT0C", 100., "FT0C threshold"};
  Configurable<float> ZDC_cut{"ZDC", 0., "ZDC threshold"};
  Configurable<float> eta_cut{"Eta", 0.9, "Eta cut"};
  Configurable<float> pt_cut{"Pt", 0.15, "Pt cut"};
  Configurable<float> nsigmatpc_cut{"nsigmatpc", 3.0, "nsigma tpc cut"};
  Configurable<float> nsigmatof_cut{"nsigmatof", 9.0, "nsigma tof cut"};

  Configurable<int> mintrack{"min_track", 1, "min track"};
  Configurable<int> maxtrack{"max_track", 50, "max track"};

  Configurable<float> dcaz_cut{"dcaz", 0.1, "dcaz cut"};
  Configurable<float> dcaxy_cut{"dcaxy", 0.1, "dcaxy cut"};
  Configurable<float> nclstpc_cut{"ncls_tpc", 50.0, "ncls tpc cut"};
  Configurable<float> itschi2_cut{"itschi2", 36.0, "its chi2 cut"};
  Configurable<float> tpcchi2_cut{"tpcchi2", 4.0, "tpc chi2 cut"};
  Configurable<bool> use_tof{"Use_TOF", true, "TOF PID"};
  Configurable<bool> QA{"QA", true, ""};
  Configurable<bool> phi{"phi", true, ""};
  Configurable<bool> rho{"rho", true, ""};
  Configurable<bool> kstar{"kstar", false, ""};

  HistogramRegistry registry{
    "registry",
    {
      // Pion histograms for each eta bin and gapSide
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
      {"fit_3d_0", "FIT amplitude", {HistType::kTH3F, {{1000, 0.0, 1000.0}, {500, 0.0, 500.0}, {80, 0.0, 80.0}}}},
      {"fit_3d_1", "FIT amplitude", {HistType::kTH3F, {{1000, 0.0, 1000.0}, {500, 0.0, 500.0}, {80, 0.0, 80.0}}}},
      {"V0A", "V0A amplitude", {HistType::kTH1F, {{1000, 0.0, 1000.0}}}},
      {"V0A_0", "V0A amplitude", {HistType::kTH1F, {{1000, 0.0, 1000.0}}}},
      {"V0A_1", "V0A amplitude", {HistType::kTH1F, {{1000, 0.0, 1000.0}}}},
      {"qt_pt_gap0", "qt vs pT", {HistType::kTH2F, {{200, 0.0, 20.0}, {200, 0.0, 20.0}}}},
      {"qt_pt_gap1", "qt vs  pT", {HistType::kTH2F, {{200, 0.0, 20.0}, {200, 0.0, 20.0}}}},
      {"qt_pt_gap2", "qt vs pT", {HistType::kTH2F, {{200, 0.0, 20.0}, {200, 0.0, 20.0}}}},
      {"track_pt_gap0", "Track pT", {HistType::kTH1F, {{200, 0.0, 20.0}}}},
      {"track_pt_gap1", "Track pT", {HistType::kTH1F, {{200, 0.0, 20.0}}}},
      {"track_pt_gap2", "Track pT", {HistType::kTH1F, {{200, 0.0, 20.0}}}},
      {"track_dpt_gap0", "Track pT", {HistType::kTH1F, {{200, 0.0, 20.0}}}},
      {"track_dpt_gap1", "Track pT", {HistType::kTH1F, {{200, 0.0, 20.0}}}},
      {"track_dpt_gap2", "Track pT", {HistType::kTH1F, {{200, 0.0, 20.0}}}},

      {"gap_mult0", "Mult", {HistType::kTH1F, {{100, 0.0, 100.0}}}},
      {"gap_mult1", "Mult", {HistType::kTH1F, {{100, 0.0, 100.0}}}},
      {"gap_mult2", "Mult", {HistType::kTH1F, {{100, 0.0, 100.0}}}},
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

      // QA plots
      {"tpc_dedx", "p vs dE/dx", {HistType::kTH2F, {{100, 0.0, 10.0}, {5000, 0.0, 5000.0}}}},
      {"tpc_dedx_kaon", "p#k dE/dx", {HistType::kTH2F, {{100, 0.0, 10.0}, {5000, 0.0, 5000.0}}}},
      {"tpc_dedx_pion", "p#pi dE/dx", {HistType::kTH2F, {{100, 0.0, 10.0}, {5000, 0.0, 5000.0}}}},
      {"tpc_nsigma_kaon", "p#k n#sigma", {HistType::kTH2F, {{100, 0.0, 10.0}, {20, -10.0, 10.0}}}},
      {"tpc_nsigma_pion", "p#pi n#sigma", {HistType::kTH2F, {{100, 0.0, 10.0}, {20, -10.0, 10.0}}}},
    }};
  using udtracks = soa::Join<aod::UDTracks, aod::UDTracksExtra, aod::UDTracksPID>;
  using udtracksfull = soa::Join<aod::UDTracks, aod::UDTracksPID, aod::UDTracksExtra, aod::UDTracksFlags, aod::UDTracksDCA>;
  using UDCollisionsFull = soa::Join<aod::UDCollisions, aod::SGCollisions, aod::UDCollisionsSels, aod::UDZdcsReduced>; //
  using UDCollisionFull = UDCollisionsFull::iterator;

  template <typename T>
  int trackselector(const T& track, bool use_tof)
  {
    TLorentzVector a;
    a.SetXYZM(track.px(), track.py(), track.pz(), mpion);
    if (std::abs(track.dcaZ()) > dcaz_cut)
      return 0;
    //    if (std::abs(track.dcaXY()) > .0105 + .035 / pow(a.Pt(), 1.1))
    if (std::abs(track.dcaXY()) > dcaxy_cut)
      return 0;
    if (track.tpcChi2NCl() > tpcchi2_cut)
      return 0;
    if (track.tpcNClsFindable() < nclstpc_cut)
      return 0;
    if (track.itsChi2NCl() > itschi2_cut)
      return 0;
    if (TMath::Sqrt(track.px() * track.px() + track.py() * track.py()) < pt_cut)
      return 0;

    return 1;
  }
  template <typename T>
  bool selectionPIDKaon(const T& candidate)
  {

    if (use_tof && candidate.hasTOF() && (candidate.tofNSigmaKa() * candidate.tofNSigmaKa() + candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa()) < nsigmatof_cut) {
      return true;
    }
    if (use_tof && !candidate.hasTOF() && std::abs(candidate.tpcNSigmaKa()) < nsigmatpc_cut) {
      return true;
    }
    if (!use_tof && std::abs(candidate.tpcNSigmaKa()) < nsigmatpc_cut) {
      return true;
    }
    return false;
  }

  template <typename T>
  bool selectionPIDPion(const T& candidate)
  {
    if (use_tof && candidate.hasTOF() && (candidate.tofNSigmaPi() * candidate.tofNSigmaPi() + candidate.tpcNSigmaPi() * candidate.tpcNSigmaPi()) < nsigmatof_cut) {
      return true;
    }
    if (use_tof && !candidate.hasTOF() && std::abs(candidate.tpcNSigmaPi()) < nsigmatpc_cut) {
      return true;
    }
    if (!use_tof && std::abs(candidate.tpcNSigmaPi()) < nsigmatpc_cut) {
      return true;
    }
    return false;
  }
  template <typename T>
  bool selectionPIDMuon(const T& candidate)
  {
    if (use_tof && candidate.hasTOF() && (candidate.tofNSigmaMu() * candidate.tofNSigmaMu() + candidate.tpcNSigmaMu() * candidate.tpcNSigmaMu()) < nsigmatof_cut) {
      return true;
    }
    if (use_tof && !candidate.hasTOF() && std::abs(candidate.tpcNSigmaMu()) < nsigmatpc_cut) {
      return true;
    }
    if (!use_tof && std::abs(candidate.tpcNSigmaMu()) < nsigmatpc_cut) {
      return true;
    }
    return false;
  }

  void process(UDCollisionFull const& collision, udtracksfull const& tracks)
  {
    TLorentzVector a;
    TLorentzVector v0;
    TLorentzVector v1;
    TLorentzVector v01;
    TLorentzVector v02, vpt, vqt;
    int gapSide = collision.gapSide();
    int truegapSide = sgSelector.trueGap(collision, FV0_cut, ZDC_cut);
    gapSide = truegapSide;
    if (gapSide < 0 || gapSide > 2)
      return;
    Int_t mult = collision.numContrib();
    Double_t t0a_amp = collision.totalFT0AmplitudeA();
    Double_t t0c_amp = collision.totalFT0AmplitudeC();
    if (gapSide == 0 && t0a_amp > FT0A_cut)
      return;
    if (gapSide == 1 && t0c_amp > FT0C_cut)
      return;
    if (gapSide == 2 && (t0c_amp > FT0A_cut || t0a_amp > FT0C_cut))
      return;
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
    TLorentzVector vpt_sum, vpt_diff;
    Double_t x, y;
    Double_t z = 0.0;
    for (auto track1 : tracks) {
      if (!trackselector(track1, use_tof))
        continue;
      v0.SetXYZM(track1.px(), track1.py(), track1.pz(), o2::constants::physics::MassPionCharged);
      auto track1ID = track1.index();
      for (auto track2 : tracks) {
        if (!trackselector(track2, use_tof))
          continue;
        auto track2ID = track2.index();
        if (track2ID <= track1ID)
          continue;
        v1.SetXYZM(track2.px(), track2.py(), track2.pz(), o2::constants::physics::MassPionCharged);
        x = v0.Pt();
        y = v1.Pt();
        if ((x + y) > z) {
          z = x + y;
          vpt_sum = v0 + v1;
          vpt_diff = 0.5 * (v0 - v1);
        }
      }
      if (QA) {
        registry.get<TH2>(HIST("tpc_dedx"))->Fill(v0.P(), track1.tpcSignal());
        registry.get<TH2>(HIST("tpc_nsigma_kaon"))->Fill(v0.Pt(), track1.tpcNSigmaKa());
        registry.get<TH2>(HIST("tpc_nsigma_pion"))->Fill(v0.Pt(), track1.tpcNSigmaPi());
        if (std::abs(track1.tpcNSigmaKa()) < 3.0) {
          registry.get<TH2>(HIST("tpc_dedx_kaon"))->Fill(v0.P(), track1.tpcSignal());
        } else if (std::abs(track1.tpcNSigmaPi()) < 3.0) {
          registry.get<TH2>(HIST("tpc_dedx_pion"))->Fill(v0.P(), track1.tpcSignal());
        }
      }
    }
    if (gapSide == 0) {
      registry.get<TH2>(HIST("qt_pt_gap0"))->Fill(vpt_sum.Pt(), vpt_diff.Pt());
      registry.get<TH1>(HIST("track_pt_gap0"))->Fill(vpt_sum.Pt());
      registry.get<TH1>(HIST("track_dpt_gap0"))->Fill(vpt_diff.Pt());
      registry.get<TH1>(HIST("V0A_0"))->Fill(collision.totalFV0AmplitudeA());
      registry.get<TH1>(HIST("FT0A_0"))->Fill(collision.totalFT0AmplitudeA());
      registry.get<TH1>(HIST("FT0C_0"))->Fill(collision.totalFT0AmplitudeC());
      registry.get<TH1>(HIST("ZDC_A_0"))->Fill(collision.energyCommonZNA());
      registry.get<TH1>(HIST("ZDC_C_0"))->Fill(collision.energyCommonZNC());
      registry.get<TH3>(HIST("fit_3d_0"))->Fill(collision.energyCommonZNA(), collision.totalFT0AmplitudeA(), mult);
    }
    if (gapSide == 1) {
      registry.get<TH2>(HIST("qt_pt_gap1"))->Fill(vpt_sum.Pt(), vpt_diff.Pt());
      registry.get<TH1>(HIST("track_pt_gap1"))->Fill(vpt_sum.Pt());
      registry.get<TH1>(HIST("track_dpt_gap1"))->Fill(vpt_diff.Pt());
      registry.get<TH1>(HIST("V0A_1"))->Fill(collision.totalFV0AmplitudeA());
      registry.get<TH1>(HIST("FT0A_1"))->Fill(collision.totalFT0AmplitudeA());
      registry.get<TH1>(HIST("FT0C_1"))->Fill(collision.totalFT0AmplitudeC());
      registry.get<TH1>(HIST("ZDC_A_1"))->Fill(collision.energyCommonZNA());
      registry.get<TH1>(HIST("ZDC_C_1"))->Fill(collision.energyCommonZNC());
      registry.get<TH3>(HIST("fit_3d_1"))->Fill(collision.energyCommonZNC(), collision.totalFT0AmplitudeC(), mult);
    }
    if (gapSide == 2) {
      registry.get<TH2>(HIST("qt_pt_gap2"))->Fill(vpt_sum.Pt(), vpt_diff.Pt());
      registry.get<TH1>(HIST("track_pt_gap2"))->Fill(vpt_sum.Pt());
      registry.get<TH1>(HIST("track_dpt_gap2"))->Fill(vpt_diff.Pt());
      registry.get<TH1>(HIST("V0A"))->Fill(collision.totalFV0AmplitudeA());
      registry.get<TH1>(HIST("FT0A"))->Fill(collision.totalFT0AmplitudeA());
      registry.get<TH1>(HIST("FT0C"))->Fill(collision.totalFT0AmplitudeC());
      registry.get<TH1>(HIST("ZDC_A"))->Fill(collision.energyCommonZNA());
      registry.get<TH1>(HIST("ZDC_C"))->Fill(collision.energyCommonZNC());
    }
    for (auto& [t0, t1] : combinations(tracks, tracks)) {
      if (!(trackselector(t0, use_tof) && trackselector(t1, use_tof)))
        continue;

      if (phi && selectionPIDKaon(t0) && selectionPIDKaon(t1)) {
        // Apply kaon hypothesis and create pairs
        v0.SetXYZM(t0.px(), t0.py(), t0.pz(), o2::constants::physics::MassKaonCharged);
        v1.SetXYZM(t1.px(), t1.py(), t1.pz(), o2::constants::physics::MassKaonCharged);
        if (v0.Eta() > std::abs(1.0) || v1.Eta() > std::abs(1.0))
          continue;
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
      if (rho && selectionPIDPion(t0) && selectionPIDPion(t1)) {
        v0.SetXYZM(t0.px(), t0.py(), t0.pz(), o2::constants::physics::MassPionCharged);
        v1.SetXYZM(t1.px(), t1.py(), t1.pz(), o2::constants::physics::MassPionCharged);
        if (v0.Eta() > std::abs(1.0) || v1.Eta() > std::abs(1.0))
          continue;
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
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<SGSpectraAnalyzer>(cfgc)};
}
