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
  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

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
  Configurable<float> pt_cut{"pt_cut", 0.15, "Track pt cut"};

  Configurable<float> EtaGapMin{"EtaGapMin", 0.0, "Track eta min"};
  Configurable<float> EtaGapMax{"EtaGapMax", 0.9, "Track eta min"};
  Configurable<float> EtaDG{"EtaDG", 0.5, "Track eta DG"};

  Configurable<float> nsigmatpc_cut{"nsigmatpc", 3.0, "nsigma tpc cut"};
  Configurable<float> nsigmatof_cut{"nsigmatof", 9.0, "nsigma tof cut"};

  Configurable<int> mintrack{"min_track", 1, "min track"};
  Configurable<int> maxtrack{"max_track", 50, "max track"};

  Configurable<bool> use_tof{"Use_TOF", true, "TOF PID"};
  Configurable<bool> QA{"QA", true, ""};
  Configurable<bool> rapidity_gap{"rapidity_gap", true, ""};

  Configurable<bool> phi{"phi", true, ""};
  Configurable<bool> rho{"rho", true, ""};
  Configurable<bool> kstar{"kstar", true, ""};
  void init(InitContext const&)
  {
    registry.add("GapSide", "Gap Side; Entries", kTH1F, {{4, -1.5, 2.5}});
    registry.add("TrueGapSide", "Gap Side; Entries", kTH1F, {{4, -1.5, 2.5}});
    if (phi) {
      registry.add("os_KK_pT_0", "pt kaon pair", kTH3F, {{100, 0.0, 10.0}, {80, -2.0, 2.0}, {220, 0.9, 1.12}});
      registry.add("os_KK_pT_1", "pt kaon pair", kTH3F, {{100, 0.0, 10.0}, {80, -2.0, 2.0}, {220, 0.9, 1.12}});
      registry.add("os_KK_pT_2", "pt kaon pair", kTH3F, {{100, 0.0, 10.0}, {80, -2.0, 2.0}, {220, 0.9, 1.12}});
      registry.add("os_KK_ls_pT_0", "kaon pair like sign", kTH3F, {{100, 0.0, 10.0}, {80, -2.0, 2.0}, {220, 0.9, 1.12}});
      registry.add("os_KK_ls_pT_1", "kaon pair like sign", kTH3F, {{100, 0.0, 10.0}, {80, -2.0, 2.0}, {220, 0.9, 1.12}});
      registry.add("os_KK_ls_pT_2", "kaon pair like sign", kTH3F, {{100, 0.0, 10.0}, {80, -2.0, 2.0}, {220, 0.9, 1.12}});
    }
    if (rho) {
      registry.add("os_pp_pT_0", "pt pion pair", kTH3F, {{100, 0.0, 10.0}, {80, -2.0, 2.0}, {350, 0.0, 3.5}});
      registry.add("os_pp_pT_1", "pt pion pair", kTH3F, {{100, 0.0, 10.0}, {80, -2.0, 2.0}, {350, 0.0, 3.5}});
      registry.add("os_pp_pT_2", "pt pion pair", kTH3F, {{100, 0.0, 10.0}, {80, -2.0, 2.0}, {350, 0.0, 3.5}});
      registry.add("os_pp_ls_pT_0", "pion pair like sign", kTH3F, {{100, 0.0, 10.0}, {80, -2.0, 2.0}, {350, 0.0, 3.5}});
      registry.add("os_pp_ls_pT_1", "pion pair like sign", kTH3F, {{100, 0.0, 10.0}, {80, -2.0, 2.0}, {350, 0.0, 3.5}});
      registry.add("os_pp_ls_pT_2", "pion pair like sign", kTH3F, {{100, 0.0, 10.0}, {80, -2.0, 2.0}, {350, 0.0, 3.5}});
    }
    if (kstar) {
      registry.add("os_pk_pT_0", "pion-kaon pair", kTH3F, {{100, 0.0, 10.0}, {80, -2.0, 2.0}, {400, 0.0, 2.0}});
      registry.add("os_pk_pT_1", "pion-kaon pair", kTH3F, {{100, 0.0, 10.0}, {80, -2.0, 2.0}, {400, 0.0, 2.0}});
      registry.add("os_pk_pT_2", "pion-kaon pair", kTH3F, {{100, 0.0, 10.0}, {80, -2.0, 2.0}, {400, 0.0, 2.0}});
      registry.add("os_pk_ls_pT_0", "pion-kaon pair like sign", kTH3F, {{100, 0.0, 10.0}, {80, -2.0, 2.0}, {400, 0.0, 2.0}});
      registry.add("os_pk_ls_pT_1", "pion-kaon like sign", kTH3F, {{100, 0.0, 10.0}, {80, -2.0, 2.0}, {400, 0.0, 2.0}});
      registry.add("os_pk_ls_pT_2", "pion-kaon like sign", kTH3F, {{100, 0.0, 10.0}, {80, -2.0, 2.0}, {400, 0.0, 2.0}});
    }
    // QA plots
    if (QA) {
      registry.add("tpc_dedx", "p vs dE/dx", kTH2F, {{100, 0.0, 10.0}, {5000, 0.0, 5000.0}});
      registry.add("tpc_dedx_kaon", "p#k dE/dx", kTH2F, {{100, 0.0, 10.0}, {5000, 0.0, 5000.0}});
      registry.add("tpc_dedx_pion", "p#pi dE/dx", kTH2F, {{100, 0.0, 10.0}, {5000, 0.0, 5000.0}});
      registry.add("tpc_dedx_kaon_1", "tpc+tof pid cut p#k dE/dx", kTH2F, {{100, 0.0, 10.0}, {5000, 0.0, 5000.0}});
      registry.add("tpc_dedx_kaon_2", "tpc+tof pid cut1 p#k dE/dx", kTH2F, {{100, 0.0, 10.0}, {5000, 0.0, 5000.0}});
      registry.add("tpc_dedx_pion_1", "tpc+tof pid cut p#pi dE/dx", kTH2F, {{100, 0.0, 10.0}, {5000, 0.0, 5000.0}});
      registry.add("tpc_nsigma_kaon", "p#k n#sigma", kTH2F, {{100, 0.0, 10.0}, {100, -10.0, 10.0}});
      registry.add("tpc_nsigma_pion", "p#pi n#sigma", kTH2F, {{100, 0.0, 10.0}, {100, -10.0, 10.0}});
      registry.add("tpc_tof_nsigma_kaon", "p#k n#sigma TPC vs TOF", kTH2F, {{100, -10.0, 10.0}, {100, -10.0, 10.0}});
      registry.add("tpc_tof_nsigma_pion", "p#pi n#sigma TPC vs TOF", kTH2F, {{100, -10.0, 10.0}, {100, -10.0, 10.0}});

      registry.add("FT0A", "T0A amplitude", kTH1F, {{500, 0.0, 500.0}});
      registry.add("FT0A_0", "T0A amplitude", kTH1F, {{500, 0.0, 500.0}});
      registry.add("FT0A_1", "T0A amplitude", kTH1F, {{20000, 0.0, 20000.0}});
      registry.add("FT0C", "T0C amplitude", kTH1F, {{500, 0.0, 500.0}});
      registry.add("FT0C_0", "T0C amplitude", kTH1F, {{20000, 0.0, 20000.0}});
      registry.add("FT0C_1", "T0C amplitude", kTH1F, {{500, 0.0, 500.0}});
      registry.add("ZDC_A", "ZDC amplitude", kTH1F, {{2000, 0.0, 1000.0}});
      registry.add("ZDC_A_0", "ZDC amplitude", kTH1F, {{2000, 0.0, 1000.0}});
      registry.add("ZDC_A_1", "ZDC amplitude", kTH1F, {{2000, 0.0, 1000.0}});
      registry.add("ZDC_C", "ZDC amplitude", kTH1F, {{2000, 0.0, 1000.0}});
      registry.add("ZDC_C_0", "ZDC amplitude", kTH1F, {{2000, 0.0, 1000.0}});
      registry.add("ZDC_C_1", "ZDC amplitude", kTH1F, {{2000, 0.0, 1000.0}});
      registry.add("V0A", "V0A amplitude", kTH1F, {{1000, 0.0, 1000.0}});
      registry.add("V0A_0", "V0A amplitude", kTH1F, {{1000, 0.0, 1000.0}});
      registry.add("V0A_1", "V0A amplitude", kTH1F, {{1000, 0.0, 1000.0}});
      if (rapidity_gap) {
        registry.add("mult_0", "mult0", kTH1F, {{150, 0, 150}});
        registry.add("mult_1", "mult1", kTH1F, {{150, 0, 150}});
        registry.add("mult_2", "mult2", kTH1F, {{150, 0, 150}});
        registry.add("mult_0_pt", "mult0_pt", kTH1F, {{150, 0, 150}});
        registry.add("mult_1_pt", "mult1_pt", kTH1F, {{150, 0, 150}});
        registry.add("mult_2_pt", "mult2_pt", kTH1F, {{150, 0, 150}});
        registry.add("mult_0_pt1", "mult0_pt1", kTH1F, {{150, 0, 150}});
        registry.add("mult_1_pt1", "mult1_pt1", kTH1F, {{150, 0, 150}});
        registry.add("mult_2_pt1", "mult2_pt1", kTH1F, {{150, 0, 150}});
        registry.add("mult_0_pt2", "mult0_pt2", kTH1F, {{150, 0, 150}});
        registry.add("mult_1_pt2", "mult1_pt2", kTH1F, {{150, 0, 150}});
        registry.add("mult_2_pt2", "mult2_pt2", kTH1F, {{150, 0, 150}});
        registry.add("event_rap_gap", "rap_gap", kTH1F, {{15, 0, 15.0}});
        registry.add("rap_mult1", "rap_mult1", kTH1F, {{150, 0, 150}});
        registry.add("rap_mult2", "rap_mult2", kTH1F, {{150, 0, 150}});
        registry.add("rap_mult3", "rap_mult3", kTH1F, {{150, 0, 150}});
        registry.add("rap1_mult1", "rap1_mult1", kTH1F, {{150, 0, 150}});
        registry.add("rap1_mult2", "rap1_mult2", kTH1F, {{150, 0, 150}});
        registry.add("rap1_mult3", "rap1_mult3", kTH1F, {{150, 0, 150}});
        registry.add("rap2_mult1", "rap2_mult1", kTH1F, {{150, 0, 150}});
        registry.add("rap2_mult2", "rap2_mult2", kTH1F, {{150, 0, 150}});
        registry.add("rap2_mult3", "rap2_mult3", kTH1F, {{150, 0, 150}});
      }
    }
    registry.add("gap_mult0", "Mult 0", kTH1F, {{100, 0.0, 100.0}});
    registry.add("gap_mult1", "Mult 1", kTH1F, {{100, 0.0, 100.0}});
    registry.add("gap_mult2", "Mult 2", kTH1F, {{100, 0.0, 100.0}});
    // Multiplicity plot
    if (rapidity_gap && phi) {
      registry.add("os_kk_mass_rap", "phi mass", kTH3F, {{220, 0.98, 1.12}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_kk_mass_rap1", "phi mass", kTH3F, {{220, 0.98, 1.12}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_kk_mass_rap2", "phi mass", kTH3F, {{220, 0.98, 1.12}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_kk_mass1_rap", "phi mass gap1", kTH3F, {{220, 0.98, 1.12}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_kk_mass1_rap1", "phi mass gap1", kTH3F, {{220, 0.98, 1.12}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_kk_mass1_rap2", "phi mass gap1", kTH3F, {{220, 0.98, 1.12}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_kk_mass2_rap", "phi mass DG", kTH3F, {{220, 0.98, 1.12}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_kk_mass2_rap1", "phi mass DG", kTH3F, {{220, 0.98, 1.12}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_kk_mass2_rap2", "phi mass DG", kTH3F, {{220, 0.98, 1.12}, {80, -2.0, 2.0}, {100, 0, 10}});
    }

    if (rapidity_gap && kstar) {
      registry.add("os_kp_mass_rap", "kstar mass", kTH3F, {{400, 0.0, 2.0}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_kp_mass_rap1", "kstar mass", kTH3F, {{400, 0.0, 2.0}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_kp_mass_rap2", "kstar mass", kTH3F, {{400, 0.0, 2.0}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_kp_mass1_rap", "kstar mass gap1", kTH3F, {{400, 0.0, 2.0}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_kp_mass1_rap1", "kstar mass gap1", kTH3F, {{400, 0.0, 2.0}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_kp_mass1_rap2", "kstar mass gap1", kTH3F, {{400, 0.0, 2.0}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_kp_mass2_rap", "kstar mass DG", kTH3F, {{400, 0.0, 2.0}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_kp_mass2_rap1", "kstar mass DG", kTH3F, {{400, 0.0, 2.0}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_kp_mass2_rap2", "kstar mass DG", kTH3F, {{400, 0.0, 2.0}, {80, -2.0, 2.0}, {100, 0, 10}});
    }
  }
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
    std::vector<float> parameters = {PV_cut, dcaZ_cut, dcaXY_cut, tpcChi2_cut, tpcNClsFindable_cut, itsChi2_cut, eta_cut, pt_cut};
    int truegapSide = sgSelector.trueGap(collision, FIT_cut[0], FIT_cut[1], FIT_cut[2], ZDC_cut);
    registry.fill(HIST("GapSide"), gapSide);
    registry.fill(HIST("TrueGapSide"), truegapSide);
    gapSide = truegapSide;
    if (gapSide < 0 || gapSide > 2)
      return;
    Int_t mult = collision.numContrib();
    if (gapSide == 0) {
      registry.fill(HIST("gap_mult0"), mult);
    }
    if (gapSide == 1) {
      registry.fill(HIST("gap_mult1"), mult);
    }
    if (gapSide == 2) {
      registry.fill(HIST("gap_mult2"), mult);
    }
    if (mult < mintrack || mult > maxtrack)
      return;
    Int_t mult0 = 0;
    Int_t mult1 = 0;
    Int_t mult2 = 0;
    Int_t mult_pt[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    Int_t trackgapA = 0;
    Int_t trackgapC = 0;
    Int_t trackDG = 0;
    Int_t trackextra = 0;
    Int_t trackextraDG = 0;
    for (auto track1 : tracks) {
      if (!trackselector(track1, parameters))
        continue;
      v0.SetXYZM(track1.px(), track1.py(), track1.pz(), o2::constants::physics::MassPionCharged);
      if (gapSide == 0) {
        mult0++;
        if (v0.Pt() < 1.0) {
          mult_pt[0]++;
        }
        if (v0.Pt() >= 1.0 && v0.Pt() < 3.0) {
          mult_pt[1]++;
        }
        if (v0.Pt() >= 3.0) {
          mult_pt[2]++;
        }
      }
      if (gapSide == 1) {
        mult1++;
        if (v0.Pt() < 1.0) {
          mult_pt[3]++;
        }
        if (v0.Pt() >= 1.0 && v0.Pt() < 3.0) {
          mult_pt[4]++;
        }
        if (v0.Pt() >= 3.0) {
          mult_pt[5]++;
        }
      }
      if (gapSide == 2) {
        mult2++;
        if (v0.Pt() < 1.0) {
          mult_pt[6]++;
        }
        if (v0.Pt() >= 1.0 && v0.Pt() < 3.0) {
          mult_pt[7]++;
        }
        if (v0.Pt() >= 3.0) {
          mult_pt[8]++;
        }
      }
      if (TMath::Abs(v0.Eta()) < EtaDG) {
        trackDG++;
      }
      if (v0.Eta() > EtaGapMin && v0.Eta() < EtaGapMax) {
        trackgapA++;
      }
      if (v0.Eta() < EtaGapMin && v0.Eta() > -EtaGapMax) {
        trackgapC++;
      }
      if (TMath::Abs(v0.Eta()) > EtaGapMax || TMath::Abs(v0.Eta()) < EtaGapMin) {
        trackextra++;
      }
      if (TMath::Abs(v0.Eta()) > EtaDG) {
        trackextraDG++;
      }

      if (QA) {
        registry.fill(HIST("tpc_dedx"), v0.P(), track1.tpcSignal());
        if (std::abs(track1.tpcNSigmaKa()) < 3.0) {
          registry.fill(HIST("tpc_dedx_kaon"), v0.P(), track1.tpcSignal());
        } else if (std::abs(track1.tpcNSigmaPi()) < 3.0) {
          registry.fill(HIST("tpc_dedx_pion"), v0.P(), track1.tpcSignal());
        }
        if (selectionPIDKaon(track1, use_tof, nsigmatpc_cut, nsigmatof_cut)) {
          registry.fill(HIST("tpc_dedx_kaon_1"), v0.P(), track1.tpcSignal());
          registry.fill(HIST("tpc_nsigma_kaon"), v0.Pt(), track1.tpcNSigmaKa());
          registry.fill(HIST("tpc_tof_nsigma_kaon"), track1.tpcNSigmaKa(), track1.tofNSigmaKa());
        }
        if (selectionPIDKaon(track1, use_tof, nsigmatpc_cut, nsigmatof_cut) && std::abs(track1.tpcNSigmaPi()) > 3.0) {
          registry.fill(HIST("tpc_dedx_kaon_2"), v0.P(), track1.tpcSignal());
        }
        if (selectionPIDPion(track1, use_tof, nsigmatpc_cut, nsigmatof_cut)) {
          registry.fill(HIST("tpc_dedx_pion_1"), v0.P(), track1.tpcSignal());
          registry.fill(HIST("tpc_nsigma_pion"), v0.Pt(), track1.tpcNSigmaPi());
          registry.fill(HIST("tpc_tof_nsigma_pion"), track1.tpcNSigmaPi(), track1.tofNSigmaPi());
        }
      }
    }
    if (QA) {
      if (gapSide == 0) {
        registry.fill(HIST("V0A_0"), collision.totalFV0AmplitudeA());
        registry.fill(HIST("FT0A_0"), collision.totalFT0AmplitudeA());
        registry.fill(HIST("FT0C_0"), collision.totalFT0AmplitudeC());
        registry.fill(HIST("ZDC_A_0"), collision.energyCommonZNA());
        registry.fill(HIST("ZDC_C_0"), collision.energyCommonZNC());
        registry.fill(HIST("mult_0"), mult0);
        registry.fill(HIST("mult_0_pt"), mult_pt[0]);
        registry.fill(HIST("mult_1_pt"), mult_pt[1]);
        registry.fill(HIST("mult_2_pt"), mult_pt[2]);
      }
      if (gapSide == 1) {
        registry.fill(HIST("V0A_1"), collision.totalFV0AmplitudeA());
        registry.fill(HIST("FT0A_1"), collision.totalFT0AmplitudeA());
        registry.fill(HIST("FT0C_1"), collision.totalFT0AmplitudeC());
        registry.fill(HIST("ZDC_A_1"), collision.energyCommonZNA());
        registry.fill(HIST("ZDC_C_1"), collision.energyCommonZNC());
        registry.fill(HIST("mult_1"), mult1);
        registry.fill(HIST("mult_0_pt1"), mult_pt[3]);
        registry.fill(HIST("mult_1_pt1"), mult_pt[4]);
        registry.fill(HIST("mult_2_pt1"), mult_pt[5]);
      }
      if (gapSide == 2) {
        registry.fill(HIST("V0A"), collision.totalFV0AmplitudeA());
        registry.fill(HIST("FT0A"), collision.totalFT0AmplitudeA());
        registry.fill(HIST("FT0C"), collision.totalFT0AmplitudeC());
        registry.fill(HIST("ZDC_A"), collision.energyCommonZNA());
        registry.fill(HIST("ZDC_C"), collision.energyCommonZNC());
        registry.fill(HIST("mult_2"), mult2);
        registry.fill(HIST("mult_0_pt2"), mult_pt[6]);
        registry.fill(HIST("mult_1_pt2"), mult_pt[7]);
        registry.fill(HIST("mult_2_pt2"), mult_pt[8]);
      }
      if (rapidity_gap) {
        if (trackgapC > 0 && trackgapA == 0 && trackextra == 0) {
          if (gapSide == 0) {
            registry.fill(HIST("event_rap_gap"), 1);
            registry.fill(HIST("rap_mult1"), trackgapC);
          }
          if (gapSide == 1) {
            registry.fill(HIST("event_rap_gap"), 4);
            registry.fill(HIST("rap1_mult1"), trackgapC);
          }
          if (gapSide == 2) {
            registry.fill(HIST("event_rap_gap"), 7);
            registry.fill(HIST("rap2_mult1"), trackgapC);
          }
        }
        if (trackgapC == 0 && trackgapA > 0 && trackextra == 0) {
          if (gapSide == 0) {
            registry.fill(HIST("event_rap_gap"), 2);
            registry.fill(HIST("rap_mult2"), trackgapA);
          }
          if (gapSide == 1) {
            registry.fill(HIST("event_rap_gap"), 5);
            registry.fill(HIST("rap1_mult2"), trackgapA);
          }
          if (gapSide == 2) {
            registry.fill(HIST("event_rap_gap"), 8);
            registry.fill(HIST("rap2_mult2"), trackgapA);
          }
        }
        if (trackDG > 0 && trackextraDG == 0) {
          if (gapSide == 0) {
            registry.fill(HIST("event_rap_gap"), 3);
            registry.fill(HIST("rap_mult3"), trackDG);
          }
          if (gapSide == 1) {
            registry.fill(HIST("event_rap_gap"), 6);
            registry.fill(HIST("rap1_mult3"), trackDG);
          }
          if (gapSide == 2) {
            registry.fill(HIST("event_rap_gap"), 9);
            registry.fill(HIST("rap2_mult3"), trackDG);
          }
        }
      }
    }

    if (rapidity_gap) {
      if (trackgapC > 0 && trackgapA == 0 && trackextra == 0) {
        for (auto& [t0, t1] : combinations(tracks, tracks)) {
          if (!trackselector(t0, parameters) || !trackselector(t1, parameters))
            continue;
          if (phi && selectionPIDKaon(t0, use_tof, nsigmatpc_cut, nsigmatof_cut) && selectionPIDKaon(t1, use_tof, nsigmatpc_cut, nsigmatof_cut)) {
            // Apply kaon hypothesis and create pairs
            v0.SetXYZM(t0.px(), t0.py(), t0.pz(), o2::constants::physics::MassKaonCharged);
            v1.SetXYZM(t1.px(), t1.py(), t1.pz(), o2::constants::physics::MassKaonCharged);
            v01 = v0 + v1;
            // Opposite sign pairs
            if (t0.sign() != t1.sign()) {
              if (gapSide == 0) {
                registry.fill(HIST("os_kk_mass_rap"), v01.M(), v01.Rapidity(), v01.Pt());
              }
              if (gapSide == 1) {
                registry.fill(HIST("os_kk_mass1_rap"), v01.M(), v01.Rapidity(), v01.Pt());
              }
              if (gapSide == 2) {
                registry.fill(HIST("os_kk_mass2_rap"), v01.M(), v01.Rapidity(), v01.Pt());
              }
            }
          }
          if (kstar && selectionPIDKaon(t0, use_tof, nsigmatpc_cut, nsigmatof_cut) && std::abs(t0.tpcNSigmaPi()) > 3.0 && selectionPIDPion(t1, use_tof, nsigmatpc_cut, nsigmatof_cut) && std::abs(t1.tpcNSigmaKa()) > 3.0) {
            // Apply kaon hypothesis and create pairs
            v0.SetXYZM(t0.px(), t0.py(), t0.pz(), o2::constants::physics::MassKaonCharged);
            v1.SetXYZM(t1.px(), t1.py(), t1.pz(), o2::constants::physics::MassPionCharged);
            v01 = v0 + v1;
            // Opposite sign pairs
            if (t0.sign() != t1.sign()) {
              if (gapSide == 0) {
                registry.fill(HIST("os_kp_mass_rap"), v01.M(), v01.Rapidity(), v01.Pt());
              }
              if (gapSide == 1) {
                registry.fill(HIST("os_kp_mass1_rap"), v01.M(), v01.Rapidity(), v01.Pt());
              }
              if (gapSide == 2) {
                registry.fill(HIST("os_kp_mass2_rap"), v01.M(), v01.Rapidity(), v01.Pt());
              }
            }
          }
        }
      }

      if (trackgapC == 0 && trackgapA > 0 && trackextra == 0) {
        for (auto& [t0, t1] : combinations(tracks, tracks)) {
          if (!trackselector(t0, parameters) || !trackselector(t1, parameters))
            continue;
          if (phi && selectionPIDKaon(t0, use_tof, nsigmatpc_cut, nsigmatof_cut) && selectionPIDKaon(t1, use_tof, nsigmatpc_cut, nsigmatof_cut)) {
            // Apply kaon hypothesis and create pairs
            v0.SetXYZM(t0.px(), t0.py(), t0.pz(), o2::constants::physics::MassKaonCharged);
            v1.SetXYZM(t1.px(), t1.py(), t1.pz(), o2::constants::physics::MassKaonCharged);
            v01 = v0 + v1;
            // Opposite sign pairs
            if (t0.sign() != t1.sign()) {
              if (gapSide == 0) {
                registry.fill(HIST("os_kk_mass_rap1"), v01.M(), v01.Rapidity(), v01.Pt());
              }
              if (gapSide == 1) {
                registry.fill(HIST("os_kk_mass1_rap1"), v01.M(), v01.Rapidity(), v01.Pt());
              }
              if (gapSide == 2) {
                registry.fill(HIST("os_kk_mass2_rap1"), v01.M(), v01.Rapidity(), v01.Pt());
              }
            }
          }
          if (kstar && selectionPIDKaon(t0, use_tof, nsigmatpc_cut, nsigmatof_cut) && std::abs(t0.tpcNSigmaPi()) > 3.0 && selectionPIDPion(t1, use_tof, nsigmatpc_cut, nsigmatof_cut) && std::abs(t1.tpcNSigmaKa()) > 3.0) {
            // Apply kaon hypothesis and create pairs
            v0.SetXYZM(t0.px(), t0.py(), t0.pz(), o2::constants::physics::MassKaonCharged);
            v1.SetXYZM(t1.px(), t1.py(), t1.pz(), o2::constants::physics::MassPionCharged);
            v01 = v0 + v1;
            // Opposite sign pairs
            if (t0.sign() != t1.sign()) {
              if (gapSide == 0) {
                registry.fill(HIST("os_kp_mass_rap1"), v01.M(), v01.Rapidity(), v01.Pt());
              }
              if (gapSide == 1) {
                registry.fill(HIST("os_kp_mass1_rap1"), v01.M(), v01.Rapidity(), v01.Pt());
              }
              if (gapSide == 2) {
                registry.fill(HIST("os_kp_mass2_rap1"), v01.M(), v01.Rapidity(), v01.Pt());
              }
            }
          }
        }
      }
      if (trackDG > 0 && trackextraDG == 0) {
        for (auto& [t0, t1] : combinations(tracks, tracks)) {
          if (!trackselector(t0, parameters) || !trackselector(t1, parameters))
            continue;
          if (phi && selectionPIDKaon(t0, use_tof, nsigmatpc_cut, nsigmatof_cut) && selectionPIDKaon(t1, use_tof, nsigmatpc_cut, nsigmatof_cut)) {
            // Apply kaon hypothesis and create pairs
            v0.SetXYZM(t0.px(), t0.py(), t0.pz(), o2::constants::physics::MassKaonCharged);
            v1.SetXYZM(t1.px(), t1.py(), t1.pz(), o2::constants::physics::MassKaonCharged);
            v01 = v0 + v1;
            // Opposite sign pairs
            if (t0.sign() != t1.sign()) {
              if (gapSide == 0) {
                registry.fill(HIST("os_kk_mass_rap2"), v01.M(), v01.Rapidity(), v01.Pt());
              }
              if (gapSide == 1) {
                registry.fill(HIST("os_kk_mass1_rap2"), v01.M(), v01.Rapidity(), v01.Pt());
              }
              if (gapSide == 2) {
                registry.fill(HIST("os_kk_mass2_rap2"), v01.M(), v01.Rapidity(), v01.Pt());
              }
            }
          }
          if (kstar && selectionPIDKaon(t0, use_tof, nsigmatpc_cut, nsigmatof_cut) && std::abs(t0.tpcNSigmaPi()) > 3.0 && selectionPIDPion(t1, use_tof, nsigmatpc_cut, nsigmatof_cut) && std::abs(t1.tpcNSigmaKa()) > 3.0) {
            // Apply kaon hypothesis and create pairs
            v0.SetXYZM(t0.px(), t0.py(), t0.pz(), o2::constants::physics::MassKaonCharged);
            v1.SetXYZM(t1.px(), t1.py(), t1.pz(), o2::constants::physics::MassPionCharged);
            v01 = v0 + v1;
            // Opposite sign pairs
            if (t0.sign() != t1.sign()) {
              if (gapSide == 0) {
                registry.fill(HIST("os_kp_mass_rap2"), v01.M(), v01.Rapidity(), v01.Pt());
              }
              if (gapSide == 1) {
                registry.fill(HIST("os_kp_mass1_rap2"), v01.M(), v01.Rapidity(), v01.Pt());
              }
              if (gapSide == 2) {
                registry.fill(HIST("os_kp_mass2_rap2"), v01.M(), v01.Rapidity(), v01.Pt());
              }
            }
          }
        }
      }
    }

    for (auto& [t0, t1] : combinations(tracks, tracks)) {
      if (!trackselector(t0, parameters) || !trackselector(t1, parameters))
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
      if (kstar && selectionPIDKaon(t0, use_tof, nsigmatpc_cut, nsigmatof_cut) && std::abs(t0.tpcNSigmaPi()) > 3.0 && selectionPIDPion(t1, use_tof, nsigmatpc_cut, nsigmatof_cut) && std::abs(t1.tpcNSigmaKa()) > 3.0) {
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
