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
#include <vector>
#include <TString.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include "Math/Vector4D.h"
#include "Math/Vector3D.h"
#include "Math/GenVector/Boost.h"

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Vertex.h"

#include "PWGUD/DataModel/UDTables.h"
#include "PWGUD/Core/SGSelector.h"
#include "PWGUD/Core/SGTrackSelector.h"
#include "PWGUD/Core/UPCHelpers.h"

#include "Common/DataModel/PIDResponse.h"

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
  Configurable<float> FT0A_cut{"FT0A", 100., "FT0A threshold"};
  Configurable<float> FT0C_cut{"FT0C", 50., "FT0C threshold"};
  Configurable<float> FDDA_cut{"FDDA", 10000., "FDDA threshold"};
  Configurable<float> FDDC_cut{"FDDC", 10000., "FDDC threshold"};
  Configurable<float> ZDC_cut{"ZDC", 0., "ZDC threshold"};
  Configurable<float> Vz_cut{"Vz_cut", 10., "Vz position"};
  Configurable<float> OccT_cut{"OccT", 1000., "Occupancy cut"};

  // Track Selections
  Configurable<float> PV_cut{"PV_cut", 1.0, "Use Only PV tracks"};
  Configurable<float> dcaZ_cut{"dcaZ_cut", 2.0, "dcaZ cut"};
  Configurable<float> dcaXY_cut{"dcaXY_cut", 0.0, "dcaXY cut (0 for Pt-function)"};
  Configurable<float> tpcChi2_cut{"tpcChi2_cut", 4, "Max tpcChi2NCl"};
  Configurable<float> tpcNClsFindable_cut{"tpcNClsFindable_cut", 70, "Min tpcNClsFindable"};
  Configurable<float> itsChi2_cut{"itsChi2_cut", 36, "Max itsChi2NCl"};
  Configurable<float> eta_cut{"eta_cut", 0.9, "Track Pseudorapidity"};
  Configurable<float> pt_cut{"pt_cut", 0.15, "Track pt cut"};
  Configurable<float> pt1{"pt1", 0.3, "pid selection pt1"};
  Configurable<float> pt2{"pt2", 0.4, "pid selection pt2"};
  Configurable<float> pt3{"pt3", 0.5, "pid selection pt3"};

  Configurable<float> EtaGapMin{"EtaGapMin", 0.0, "Track eta min"};
  Configurable<float> EtaGapMax{"EtaGapMax", 0.9, "Track eta min"};
  Configurable<float> EtaDG{"EtaDG", 0.5, "Track eta DG"};

  Configurable<float> nsigmatpc_cut1{"nsigmatpc1", 3.0, "nsigma tpc cut1"};
  Configurable<float> nsigmatpc_cut2{"nsigmatpc2", 3.0, "nsigma tpc cut2"};
  Configurable<float> nsigmatpc_cut3{"nsigmatpc3", 3.0, "nsigma tpc cut3"};
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
  Configurable<bool> fourpion{"fourpion", true, ""};

  Configurable<int> cfgNoMixedEvents{"cfgNoMixedEvents", 1, "Number of mixed events per event"};
  Configurable<int> nBkgRotations{"nBkgRotations", 9, "Number of rotated copies (background) per each original candidate"};
  Configurable<bool> fillRotation{"fillRotation", true, "fill rotation"};
  Configurable<float> confMinRot{"confMinRot", 5.0 * TMath::Pi() / 6.0, "Minimum of rotation"};
  Configurable<float> confMaxRot{"confMaxRot", 7.0 * TMath::Pi() / 6.0, "Maximum of rotation"};

  void init(InitContext const&)
  {
    registry.add("GapSide", "Gap Side; Entries", kTH1F, {{4, -1.5, 2.5}});
    registry.add("TrueGapSide", "Gap Side; Entries", kTH1F, {{4, -1.5, 2.5}});
    if (phi) {
      registry.add("os_KK_pT_0", "pt kaon pair", kTH3F, {{220, 0.98, 1.2}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_KK_pT_1", "pt kaon pair", kTH3F, {{220, 0.98, 1.2}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_KK_pT_2", "pt kaon pair", kTH3F, {{220, 0.98, 1.2}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_KK_ls_pT_0", "kaon pair like sign", kTH3F, {{220, 0.98, 1.2}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_KK_ls_pT_1", "kaon pair like sign", kTH3F, {{220, 0.98, 1.2}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_KK_ls_pT_2", "kaon pair like sign", kTH3F, {{220, 0.98, 1.2}, {80, -2.0, 2.0}, {100, 0, 10}});

      registry.add("os_KK_mix_pT_0", "kaon pair mix event", kTH3F, {{220, 0.98, 1.2}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_KK_mix_pT_1", "kaon pair mix event", kTH3F, {{220, 0.98, 1.2}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_KK_mix_pT_2", "kaon pair mix event", kTH3F, {{220, 0.98, 1.2}, {80, -2.0, 2.0}, {100, 0, 10}});

      registry.add("os_KK_rot_pT_0", "kaon pair mix event", kTH3F, {{220, 0.98, 1.2}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_KK_rot_pT_1", "kaon pair mix event", kTH3F, {{220, 0.98, 1.2}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_KK_rot_pT_2", "kaon pair mix event", kTH3F, {{220, 0.98, 1.2}, {80, -2.0, 2.0}, {100, 0, 10}});
    }
    if (rho) {
      registry.add("os_pp_pT_0", "pt pion pair", kTH3F, {{120, 1.44, 2.04}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_pp_pT_1", "pt pion pair", kTH3F, {{120, 1.44, 2.04}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_pp_pT_2", "pt pion pair", kTH3F, {{120, 1.44, 2.04}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_pp_ls_pT_0", "pion pair like sign", kTH3F, {{120, 1.44, 2.04}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_pp_ls_pT_1", "pion pair like sign", kTH3F, {{120, 1.44, 2.04}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_pp_ls_pT_2", "pion pair like sign", kTH3F, {{120, 1.44, 2.04}, {80, -2.0, 2.0}, {100, 0, 10}});
    }
    if (kstar) {
      registry.add("os_pk_pT_0", "pion-kaon pair", kTH3F, {{400, 0.0, 2.0}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_pk_pT_1", "pion-kaon pair", kTH3F, {{400, 0.0, 2.0}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_pk_pT_2", "pion-kaon pair", kTH3F, {{400, 0.0, 2.0}, {80, -2.0, 2.0}, {100, 0, 10}});

      registry.add("os_pk_mix_pT_0", "pion-kaon mix pair", kTH3F, {{400, 0.0, 2.0}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_pk_mix_pT_1", "pion-kaon mix pair", kTH3F, {{400, 0.0, 2.0}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_pk_mix_pT_2", "pion-kaon mix pair", kTH3F, {{400, 0.0, 2.0}, {80, -2.0, 2.0}, {100, 0, 10}});

      registry.add("os_pk_rot_pT_0", "pion-kaon rotional pair", kTH3F, {{400, 0.0, 2.0}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_pk_rot_pT_1", "pion-kaon rotional pair", kTH3F, {{400, 0.0, 2.0}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_pk_rot_pT_2", "pion-kaon rotional pair", kTH3F, {{400, 0.0, 2.0}, {80, -2.0, 2.0}, {100, 0, 10}});

      registry.add("os_pk_ls_pT_0", "pion-kaon pair like sign", kTH3F, {{400, 0.0, 2.0}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_pk_ls_pT_1", "pion-kaon like sign", kTH3F, {{400, 0.0, 2.0}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_pk_ls_pT_2", "pion-kaon like sign", kTH3F, {{400, 0.0, 2.0}, {80, -2.0, 2.0}, {100, 0, 10}});

      registry.add("hRotation", "hRotation", kTH1F, {{360, 0.0, 2.0 * TMath::Pi()}});
    }
    // QA plots
    if (QA) {
      registry.add("tpc_dedx", "p vs dE/dx", kTH2F, {{100, 0.0, 10.0}, {5000, 0.0, 5000.0}});
      registry.add("tof_beta", "p vs beta", kTH2F, {{100, 0.0, 10.0}, {5000, 0.0, 5000.0}});

      registry.add("tpc_dedx_kaon", "p#k dE/dx", kTH2F, {{100, 0.0, 10.0}, {5000, 0.0, 5000.0}});
      registry.add("tpc_dedx_pion", "p#pi dE/dx", kTH2F, {{100, 0.0, 10.0}, {5000, 0.0, 5000.0}});
      registry.add("tpc_dedx_kaon_1", "tpc+tof pid cut p#k dE/dx", kTH2F, {{100, 0.0, 10.0}, {5000, 0.0, 5000.0}});
      registry.add("tpc_dedx_kaon_2", "tpc+tof pid cut1 p#k dE/dx", kTH2F, {{100, 0.0, 10.0}, {5000, 0.0, 5000.0}});
      registry.add("tpc_dedx_pion_1", "tpc+tof pid cut p#pi dE/dx", kTH2F, {{100, 0.0, 10.0}, {5000, 0.0, 5000.0}});
      registry.add("tpc_nsigma_kaon", "p#k n#sigma", kTH2F, {{100, 0.0, 10.0}, {100, -10.0, 10.0}});
      registry.add("tpc_nsigma_pion", "p#pi n#sigma", kTH2F, {{100, 0.0, 10.0}, {100, -10.0, 10.0}});
      registry.add("tpc_tof_nsigma_kaon", "p#k n#sigma TPC vs TOF", kTH2F, {{100, -10.0, 10.0}, {100, -10.0, 10.0}});
      registry.add("tpc_tof_nsigma_pion", "p#pi n#sigma TPC vs TOF", kTH2F, {{100, -10.0, 10.0}, {100, -10.0, 10.0}});

      registry.add("tof_nsigma_kaon", "p#k n#sigma", kTH2F, {{100, 0.0, 10.0}, {100, -10.0, 10.0}});
      registry.add("tof_nsigma_pion", "p#pi n#sigma", kTH2F, {{100, 0.0, 10.0}, {100, -10.0, 10.0}});

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
      registry.add("os_kk_mass_rap", "phi mass1", kTH3F, {{220, 0.98, 1.2}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_kk_mass_rap1", "phi mass2", kTH3F, {{220, 0.98, 1.2}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_kk_mass_rap2", "phi mass3", kTH3F, {{220, 0.98, 1.2}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_kk_mass1_rap", "phi mass1 gap1", kTH3F, {{220, 0.98, 1.2}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_kk_mass1_rap1", "phi mass2 gap1", kTH3F, {{220, 0.98, 1.2}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_kk_mass1_rap2", "phi mass3 gap1", kTH3F, {{220, 0.98, 1.2}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_kk_mass2_rap", "phi mass1 DG", kTH3F, {{220, 0.98, 1.2}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_kk_mass2_rap1", "phi mass2 DG", kTH3F, {{220, 0.98, 1.2}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_kk_mass2_rap2", "phi mass3 DG", kTH3F, {{220, 0.98, 1.2}, {80, -2.0, 2.0}, {100, 0, 10}});

      // like sign bkg
      registry.add("os_kk_ls_mass_rap", "phi ls mass1", kTH3F, {{220, 0.98, 1.2}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_kk_ls_mass_rap1", "phi ls mass2", kTH3F, {{220, 0.98, 1.2}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_kk_ls_mass_rap2", "phi ls mass3", kTH3F, {{220, 0.98, 1.2}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_kk_ls_mass1_rap", "phi ls mass1 gap1", kTH3F, {{220, 0.98, 1.2}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_kk_ls_mass1_rap1", "phi ls mass2 gap1", kTH3F, {{220, 0.98, 1.2}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_kk_ls_mass1_rap2", "phi ls mass3 gap1", kTH3F, {{220, 0.98, 1.2}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_kk_ls_mass2_rap", "phi ls mass1 DG", kTH3F, {{220, 0.98, 1.2}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_kk_ls_mass2_rap1", "phi ls mass2 DG", kTH3F, {{220, 0.98, 1.2}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_kk_ls_mass2_rap2", "phi ls mass3 DG", kTH3F, {{220, 0.98, 1.2}, {80, -2.0, 2.0}, {100, 0, 10}});
    }

    if (rapidity_gap && kstar) {
      registry.add("os_kp_mass_rap", "kstar mass1", kTH3F, {{400, 0.0, 2.0}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_kp_mass_rap1", "kstar mass2", kTH3F, {{400, 0.0, 2.0}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_kp_mass_rap2", "kstar mass3", kTH3F, {{400, 0.0, 2.0}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_kp_mass1_rap", "kstar mass1 gap1", kTH3F, {{400, 0.0, 2.0}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_kp_mass1_rap1", "kstar mass2 gap1", kTH3F, {{400, 0.0, 2.0}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_kp_mass1_rap2", "kstar mass3 gap1", kTH3F, {{400, 0.0, 2.0}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_kp_mass2_rap", "kstar mass1 DG", kTH3F, {{400, 0.0, 2.0}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_kp_mass2_rap1", "kstar mass2 DG", kTH3F, {{400, 0.0, 2.0}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_kp_mass2_rap2", "kstar mass3 DG", kTH3F, {{400, 0.0, 2.0}, {80, -2.0, 2.0}, {100, 0, 10}});

      // like sign bkg

      registry.add("os_kp_ls_mass_rap", "kstar ls mass1", kTH3F, {{400, 0.0, 2.0}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_kp_ls_mass_rap1", "kstar ls mass2", kTH3F, {{400, 0.0, 2.0}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_kp_ls_mass_rap2", "kstar ls mass3", kTH3F, {{400, 0.0, 2.0}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_kp_ls_mass1_rap", "kstar ls mass1 gap1", kTH3F, {{400, 0.0, 2.0}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_kp_ls_mass1_rap1", "kstar ls mass2 gap1", kTH3F, {{400, 0.0, 2.0}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_kp_ls_mass1_rap2", "kstar ls mass3 gap1", kTH3F, {{400, 0.0, 2.0}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_kp_ls_mass2_rap", "kstar ls mass1 DG", kTH3F, {{400, 0.0, 2.0}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_kp_ls_mass2_rap1", "kstar ls mass2 DG", kTH3F, {{400, 0.0, 2.0}, {80, -2.0, 2.0}, {100, 0, 10}});
      registry.add("os_kp_ls_mass2_rap2", "kstar ls mass3 DG", kTH3F, {{400, 0.0, 2.0}, {80, -2.0, 2.0}, {100, 0, 10}});
    }
    if (fourpion) {
      registry.add("os_pppp_pT_2", "4 pion pair", kTH3F, {{800, 0.5, 4.5}, {250, 0.0, 5.0}, {30, -1.5, 1.5}});
      registry.add("os_pppp_pT_2_ls", "4 pion pair", kTH3F, {{800, 0.5, 4.5}, {250, 0.0, 5.0}, {30, -1.5, 1.5}});
      registry.add("os_pp_vs_pp_mass", "pair1 vd pair2 ", kTH2F, {{800, 0.5, 4.5}, {800, 0.5, 4.5}});
      registry.add("os_pp_vs_pp_pt", "pair1 pt vs pair2 pt", kTH2F, {{250, 0.0, 5.0}, {250, 0.0, 5.0}});
      registry.add("os_pp_vs_pp_mass1", "pair3 vd pair4 ", kTH2F, {{800, 0.5, 4.5}, {800, 0.5, 4.5}});
      registry.add("os_pp_vs_pp_pt1", "pair3 pt vs pair4 pt", kTH2F, {{250, 0.0, 5.0}, {250, 0.0, 5.0}});
      registry.add("phi_dis", "phi_dis", kTH1F, {{360, 0, 6.28}});
      registry.add("costheta_dis", "costheta_dis", kTH1F, {{40, -1.0, 1.0}});
      registry.add("costheta_vs_phi", "costheta_vs_phi", kTH2F, {{40, -1.0, 1.0}, {360, 0.0, 6.28}});
      registry.add("phi_dis1", "phi_dis1", kTH1F, {{360, 0, 6.28}});
      registry.add("costheta_dis1", "costheta_dis1", kTH1F, {{40, -1.0, 1.0}});
      registry.add("costheta_vs_phi1", "costheta_vs_phi1", kTH2F, {{40, -1.0, 1.0}, {360, 0.0, 6.28}});
    }
  }

  //_____________________________________________________________________________
  Double_t CosThetaCollinsSoperFrame(ROOT::Math::PtEtaPhiMVector pair1,
                                     ROOT::Math::PtEtaPhiMVector pair2,
                                     ROOT::Math::PtEtaPhiMVector fourpion)
  {
    Double_t HalfSqrtSnn = 2680.;
    Double_t MassOfLead208 = 193.6823;
    Double_t MomentumBeam = TMath::Sqrt(HalfSqrtSnn * HalfSqrtSnn * 208 * 208 - MassOfLead208 * MassOfLead208);

    TLorentzVector pProjCM(0., 0., -MomentumBeam, HalfSqrtSnn * 208); // projectile
    TLorentzVector pTargCM(0., 0., MomentumBeam, HalfSqrtSnn * 208);  // target

    //  TVector3 beta = (-1. / fourpion.E()) * fourpion.Vect();
    ROOT::Math::PtEtaPhiMVector v1 = pair1;
    ROOT::Math::PtEtaPhiMVector v2 = pair2;
    ROOT::Math::PtEtaPhiMVector v12 = fourpion;

    // Boost to center of mass frame
    ROOT::Math::Boost boostv12{v12.BoostToCM()};
    ROOT::Math::XYZVectorF v1_CM{(boostv12(v1).Vect()).Unit()};
    ROOT::Math::XYZVectorF v2_CM{(boostv12(v2).Vect()).Unit()};
    ROOT::Math::XYZVectorF Beam1_CM{(boostv12(pProjCM).Vect()).Unit()};
    ROOT::Math::XYZVectorF Beam2_CM{(boostv12(pTargCM).Vect()).Unit()};

    // Axes
    ROOT::Math::XYZVectorF zaxis_CS{((Beam1_CM.Unit() - Beam2_CM.Unit()).Unit())};

    Double_t CosThetaCS = zaxis_CS.Dot((v1_CM));
    return CosThetaCS;
  }

  template <typename T>
  bool selectionPIDKaon1(const T& candidate)
  {
    auto pt = TMath::Sqrt(candidate.px() * candidate.px() + candidate.py() * candidate.py());
    //  float pt1, pt2, pt3 , nsigmatpc_cut1, nsigmatpc_cut2, nsigmatpc_cut3;
    if (use_tof && pt < pt1 && std::abs(candidate.tpcNSigmaKa()) < nsigmatpc_cut1) {
      return true;
    }
    if (use_tof && pt >= pt1 && pt < pt2 && std::abs(candidate.tpcNSigmaKa()) < nsigmatpc_cut2) {
      return true;
    }
    if (use_tof && pt >= pt2 && pt < pt3 && std::abs(candidate.tpcNSigmaKa()) < nsigmatpc_cut3) {
      return true;
    }
    if (use_tof && pt >= pt3 && !candidate.hasTOF() && std::abs(candidate.tpcNSigmaKa()) < nsigmatpc_cut) {
      return true;
    }
    if (use_tof && candidate.hasTOF() && (candidate.tofNSigmaKa() * candidate.tofNSigmaKa() + candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa()) < nsigmatof_cut) {
      return true;
    }
    if (!use_tof && std::abs(candidate.tpcNSigmaKa()) < nsigmatpc_cut) {
      return true;
    }
    return false;
  }

  template <typename T>
  bool selectionPIDPion1(const T& candidate)
  {
    auto pt = TMath::Sqrt(candidate.px() * candidate.px() + candidate.py() * candidate.py());

    if (use_tof && pt < pt1 && std::abs(candidate.tpcNSigmaPi()) < nsigmatpc_cut1) {
      return true;
    }
    if (use_tof && pt >= pt1 && pt < pt2 && std::abs(candidate.tpcNSigmaPi()) < nsigmatpc_cut2) {
      return true;
    }
    if (use_tof && pt >= pt2 && pt < pt3 && std::abs(candidate.tpcNSigmaPi()) < nsigmatpc_cut3) {
      return true;
    }
    if (use_tof && pt >= pt3 && !candidate.hasTOF() && std::abs(candidate.tpcNSigmaPi()) < nsigmatpc_cut) {
      return true;
    }
    if (use_tof && candidate.hasTOF() && (candidate.tofNSigmaPi() * candidate.tofNSigmaPi() + candidate.tpcNSigmaPi() * candidate.tpcNSigmaPi()) < nsigmatof_cut) {
      return true;
    }
    if (!use_tof && std::abs(candidate.tpcNSigmaPi()) < nsigmatpc_cut) {
      return true;
    }
    return false;
  }

  //------------------------------------------------------------------------------------------------------
  Double_t PhiCollinsSoperFrame(ROOT::Math::PtEtaPhiMVector pair1, ROOT::Math::PtEtaPhiMVector pair2, ROOT::Math::PtEtaPhiMVector fourpion)
  {
    // Half of the energy per pair of the colliding nucleons.
    Double_t HalfSqrtSnn = 2680.;
    Double_t MassOfLead208 = 193.6823;
    Double_t MomentumBeam = TMath::Sqrt(HalfSqrtSnn * HalfSqrtSnn * 208 * 208 - MassOfLead208 * MassOfLead208);

    TLorentzVector pProjCM(0., 0., -MomentumBeam, HalfSqrtSnn * 208); // projectile
    TLorentzVector pTargCM(0., 0., MomentumBeam, HalfSqrtSnn * 208);  // target
    ROOT::Math::PtEtaPhiMVector v1 = pair1;
    ROOT::Math::PtEtaPhiMVector v2 = pair2;
    ROOT::Math::PtEtaPhiMVector v12 = fourpion;
    // Boost to center of mass frame
    ROOT::Math::Boost boostv12{v12.BoostToCM()};
    ROOT::Math::XYZVectorF v1_CM{(boostv12(v1).Vect()).Unit()};
    ROOT::Math::XYZVectorF v2_CM{(boostv12(v2).Vect()).Unit()};
    ROOT::Math::XYZVectorF Beam1_CM{(boostv12(pProjCM).Vect()).Unit()};
    ROOT::Math::XYZVectorF Beam2_CM{(boostv12(pTargCM).Vect()).Unit()};
    // Axes
    ROOT::Math::XYZVectorF zaxis_CS{((Beam1_CM.Unit() - Beam2_CM.Unit()).Unit())};
    ROOT::Math::XYZVectorF yaxis_CS{(Beam1_CM.Cross(Beam2_CM)).Unit()};
    ROOT::Math::XYZVectorF xaxis_CS{(yaxis_CS.Cross(zaxis_CS)).Unit()};

    Double_t phi = TMath::ATan2(yaxis_CS.Dot(v1_CM), xaxis_CS.Dot(v1_CM));
    return phi;
  }

  using udtracks = soa::Join<aod::UDTracks, aod::UDTracksExtra, aod::UDTracksPID>;
  using udtracksfull = soa::Join<aod::UDTracks, aod::UDTracksPID, aod::UDTracksExtra, aod::UDTracksFlags, aod::UDTracksDCA>;
  using UDCollisionsFull = soa::Join<aod::UDCollisions, aod::SGCollisions, aod::UDCollisionSelExtras, aod::UDCollisionsSels, aod::UDZdcsReduced>; //
  using UDCollisionFull = UDCollisionsFull::iterator;

  void process(UDCollisionFull const& collision, udtracksfull const& tracks)
  {
    TLorentzVector v0;
    TLorentzVector v1;
    TLorentzVector v01;
    TLorentzVector v0_1;
    int gapSide = collision.gapSide();
    float FIT_cut[5] = {FV0_cut, FT0A_cut, FT0C_cut, FDDA_cut, FDDC_cut};
    std::vector<float> parameters = {PV_cut, dcaZ_cut, dcaXY_cut, tpcChi2_cut, tpcNClsFindable_cut, itsChi2_cut, eta_cut, pt_cut};
    int truegapSide = sgSelector.trueGap(collision, FIT_cut[0], FIT_cut[1], FIT_cut[2], ZDC_cut);

    ROOT::Math::PtEtaPhiMVector phiv;
    ROOT::Math::PtEtaPhiMVector phiv1;

    std::vector<ROOT::Math::PtEtaPhiMVector> onlyPionTracks_p;
    std::vector<decltype(tracks.begin())> rawPionTracks_p;

    std::vector<ROOT::Math::PtEtaPhiMVector> onlyPionTracks_pm;
    std::vector<decltype(tracks.begin())> rawPionTracks_pm;

    std::vector<ROOT::Math::PtEtaPhiMVector> onlyPionTracks_n;
    std::vector<decltype(tracks.begin())> rawPionTracks_n;

    registry.fill(HIST("GapSide"), gapSide);
    registry.fill(HIST("TrueGapSide"), truegapSide);
    gapSide = truegapSide;
    if (gapSide < 0 || gapSide > 2)
      return;
    if (std::abs(collision.posZ()) > Vz_cut)
      return;
    if (std::abs(collision.occupancyInTime()) > OccT_cut)
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
    Int_t trackgapA = 0;
    Int_t trackgapC = 0;
    Int_t trackDG = 0;
    Int_t trackextra = 0;
    Int_t trackextraDG = 0;
    for (auto track1 : tracks) {
      if (!trackselector(track1, parameters))
        continue;
      v0.SetXYZM(track1.px(), track1.py(), track1.pz(), o2::constants::physics::MassPionCharged);
      ROOT::Math::PtEtaPhiMVector vv1(v0.Pt(), v0.Eta(), v0.Phi(), o2::constants::physics::MassPionCharged);
      if (selectionPIDPion1(track1)) {
        onlyPionTracks_pm.push_back(vv1);
        rawPionTracks_pm.push_back(track1);
        if (track1.sign() == 1) {
          onlyPionTracks_p.push_back(vv1);
          rawPionTracks_p.push_back(track1);
        }
        if (track1.sign() == -1) {
          onlyPionTracks_n.push_back(vv1);
          rawPionTracks_n.push_back(track1);
        }
      }

      if (gapSide == 0) {
        mult0++;
      }
      if (gapSide == 1) {
        mult1++;
      }
      if (gapSide == 2) {
        mult2++;
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
        registry.fill(HIST("tof_beta"), v0.P(), track1.beta());
        registry.fill(HIST("tof_nsigma_kaon"), v0.Pt(), track1.tofNSigmaKa());
        registry.fill(HIST("tof_nsigma_pion"), v0.Pt(), track1.tofNSigmaPi());

        if (std::abs(track1.tpcNSigmaKa()) < 3.0) {
          registry.fill(HIST("tpc_dedx_kaon"), v0.P(), track1.tpcSignal());
        } else if (std::abs(track1.tpcNSigmaPi()) < 3.0) {
          registry.fill(HIST("tpc_dedx_pion"), v0.P(), track1.tpcSignal());
        }
        if (selectionPIDKaon1(track1)) {
          registry.fill(HIST("tpc_dedx_kaon_1"), v0.P(), track1.tpcSignal());
          registry.fill(HIST("tpc_nsigma_kaon"), v0.Pt(), track1.tpcNSigmaKa());
          registry.fill(HIST("tpc_tof_nsigma_kaon"), track1.tpcNSigmaKa(), track1.tofNSigmaKa());
        }
        if (selectionPIDKaon1(track1) && std::abs(track1.tpcNSigmaPi()) > 3.0) {
          registry.fill(HIST("tpc_dedx_kaon_2"), v0.P(), track1.tpcSignal());
        }
        if (selectionPIDPion1(track1)) {
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
      }
      if (gapSide == 1) {
        registry.fill(HIST("V0A_1"), collision.totalFV0AmplitudeA());
        registry.fill(HIST("FT0A_1"), collision.totalFT0AmplitudeA());
        registry.fill(HIST("FT0C_1"), collision.totalFT0AmplitudeC());
        registry.fill(HIST("ZDC_A_1"), collision.energyCommonZNA());
        registry.fill(HIST("ZDC_C_1"), collision.energyCommonZNC());
        registry.fill(HIST("mult_1"), mult1);
      }
      if (gapSide == 2) {
        registry.fill(HIST("V0A"), collision.totalFV0AmplitudeA());
        registry.fill(HIST("FT0A"), collision.totalFT0AmplitudeA());
        registry.fill(HIST("FT0C"), collision.totalFT0AmplitudeC());
        registry.fill(HIST("ZDC_A"), collision.energyCommonZNA());
        registry.fill(HIST("ZDC_C"), collision.energyCommonZNC());
        registry.fill(HIST("mult_2"), mult2);
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
          if (phi && selectionPIDKaon1(t0) && selectionPIDKaon1(t1)) {
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
            if (t0.sign() == t1.sign()) {
              if (gapSide == 0) {
                registry.fill(HIST("os_kk_ls_mass_rap"), v01.M(), v01.Rapidity(), v01.Pt());
              }
              if (gapSide == 1) {
                registry.fill(HIST("os_kk_ls_mass1_rap"), v01.M(), v01.Rapidity(), v01.Pt());
              }
              if (gapSide == 2) {
                registry.fill(HIST("os_kk_ls_mass2_rap"), v01.M(), v01.Rapidity(), v01.Pt());
              }
            }
          }
        }
        for (auto& [t0, t1] : combinations(o2::soa::CombinationsFullIndexPolicy(tracks, tracks))) {
          if (!trackselector(t0, parameters) || !trackselector(t1, parameters))
            continue;
          if (t0.globalIndex() == t1.globalIndex())
            continue;
          if (kstar && selectionPIDKaon1(t0) && selectionPIDPion1(t1)) {
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
            if (t0.sign() == t1.sign()) {
              if (gapSide == 0) {
                registry.fill(HIST("os_kp_ls_mass_rap"), v01.M(), v01.Rapidity(), v01.Pt());
              }
              if (gapSide == 1) {
                registry.fill(HIST("os_kp_ls_mass1_rap"), v01.M(), v01.Rapidity(), v01.Pt());
              }
              if (gapSide == 2) {
                registry.fill(HIST("os_kp_ls_mass2_rap"), v01.M(), v01.Rapidity(), v01.Pt());
              }
            }
          }
        }
      }

      if (trackgapC == 0 && trackgapA > 0 && trackextra == 0) {
        for (auto& [t0, t1] : combinations(tracks, tracks)) {
          if (!trackselector(t0, parameters) || !trackselector(t1, parameters))
            continue;
          if (phi && selectionPIDKaon1(t0) && selectionPIDKaon1(t1)) {
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
            if (t0.sign() == t1.sign()) {
              if (gapSide == 0) {
                registry.fill(HIST("os_kk_ls_mass_rap1"), v01.M(), v01.Rapidity(), v01.Pt());
              }
              if (gapSide == 1) {
                registry.fill(HIST("os_kk_ls_mass1_rap1"), v01.M(), v01.Rapidity(), v01.Pt());
              }
              if (gapSide == 2) {
                registry.fill(HIST("os_kk_ls_mass2_rap1"), v01.M(), v01.Rapidity(), v01.Pt());
              }
            }
          }
        }
        for (auto& [t0, t1] : combinations(o2::soa::CombinationsFullIndexPolicy(tracks, tracks))) {
          if (!trackselector(t0, parameters) || !trackselector(t1, parameters))
            continue;
          if (t0.globalIndex() == t1.globalIndex())
            continue;
          if (kstar && selectionPIDKaon1(t0) && selectionPIDPion1(t1)) {
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
            if (t0.sign() == t1.sign()) {
              if (gapSide == 0) {
                registry.fill(HIST("os_kp_ls_mass_rap1"), v01.M(), v01.Rapidity(), v01.Pt());
              }
              if (gapSide == 1) {
                registry.fill(HIST("os_kp_ls_mass1_rap1"), v01.M(), v01.Rapidity(), v01.Pt());
              }
              if (gapSide == 2) {
                registry.fill(HIST("os_kp_ls_mass2_rap1"), v01.M(), v01.Rapidity(), v01.Pt());
              }
            }
          }
        }
      }
      if (trackDG > 0 && trackextraDG == 0) {
        for (auto& [t0, t1] : combinations(tracks, tracks)) {
          if (!trackselector(t0, parameters) || !trackselector(t1, parameters))
            continue;
          if (phi && selectionPIDKaon1(t0) && selectionPIDKaon1(t1)) {
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
            if (t0.sign() == t1.sign()) {
              if (gapSide == 0) {
                registry.fill(HIST("os_kk_ls_mass_rap2"), v01.M(), v01.Rapidity(), v01.Pt());
              }
              if (gapSide == 1) {
                registry.fill(HIST("os_kk_ls_mass1_rap2"), v01.M(), v01.Rapidity(), v01.Pt());
              }
              if (gapSide == 2) {
                registry.fill(HIST("os_kk_ls_mass2_rap2"), v01.M(), v01.Rapidity(), v01.Pt());
              }
            }
          }
        }
        for (auto& [t0, t1] : combinations(o2::soa::CombinationsFullIndexPolicy(tracks, tracks))) {
          if (!trackselector(t0, parameters) || !trackselector(t1, parameters))
            continue;
          if (t0.globalIndex() == t1.globalIndex())
            continue;
          if (kstar && selectionPIDKaon1(t0) && selectionPIDPion1(t1)) {
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
            if (t0.sign() == t1.sign()) {
              if (gapSide == 0) {
                registry.fill(HIST("os_kp_ls_mass_rap2"), v01.M(), v01.Rapidity(), v01.Pt());
              }
              if (gapSide == 1) {
                registry.fill(HIST("os_kp_ls_mass1_rap2"), v01.M(), v01.Rapidity(), v01.Pt());
              }
              if (gapSide == 2) {
                registry.fill(HIST("os_kp_ls_mass2_rap2"), v01.M(), v01.Rapidity(), v01.Pt());
              }
            }
          }
        }
      }
    }

    for (auto& [t0, t1] : combinations(tracks, tracks)) {
      if (!trackselector(t0, parameters) || !trackselector(t1, parameters))
        continue;

      if (phi && selectionPIDKaon1(t0) && selectionPIDKaon1(t1)) {
        // Apply kaon hypothesis and create pairs
        v0.SetXYZM(t0.px(), t0.py(), t0.pz(), o2::constants::physics::MassKaonCharged);
        v1.SetXYZM(t1.px(), t1.py(), t1.pz(), o2::constants::physics::MassKaonCharged);
        v01 = v0 + v1;
        // Opposite sign pairs
        if (t0.sign() != t1.sign()) {
          if (gapSide == 0) {
            registry.fill(HIST("os_KK_pT_0"), v01.M(), v01.Rapidity(), v01.Pt());
          }
          if (gapSide == 1) {
            registry.fill(HIST("os_KK_pT_1"), v01.M(), v01.Rapidity(), v01.Pt());
          }
          if (gapSide == 2) {
            registry.fill(HIST("os_KK_pT_2"), v01.M(), v01.Rapidity(), v01.Pt());
          }
        }
        // samesignpair
        if (t0.sign() == t1.sign()) {
          if (gapSide == 0) {
            registry.fill(HIST("os_KK_ls_pT_0"), v01.M(), v01.Rapidity(), v01.Pt());
          }
          if (gapSide == 1) {
            registry.fill(HIST("os_KK_ls_pT_1"), v01.M(), v01.Rapidity(), v01.Pt());
          }
          if (gapSide == 2) {
            registry.fill(HIST("os_KK_ls_pT_2"), v01.M(), v01.Rapidity(), v01.Pt());
          }
        }

        if (fillRotation) {
          for (int nrotbkg = 0; nrotbkg < nBkgRotations; nrotbkg++) {
            auto anglestart = confMinRot;
            auto angleend = confMaxRot;
            auto anglestep = (angleend - anglestart) / (1.0 * (nBkgRotations - 1));
            auto rotangle = anglestart + nrotbkg * anglestep;
            registry.fill(HIST("hRotation"), rotangle);

            auto rotkaonPx = t0.px() * std::cos(rotangle) - t0.py() * std::sin(rotangle);
            auto rotkaonPy = t0.px() * std::sin(rotangle) + t0.py() * std::cos(rotangle);

            v0.SetXYZM(rotkaonPx, rotkaonPy, t0.pz(), o2::constants::physics::MassKaonCharged);
            v1.SetXYZM(t1.px(), t1.py(), t1.pz(), o2::constants::physics::MassKaonCharged);
            v01 = v0 + v1;
            if (t0.sign() != t1.sign()) {
              if (gapSide == 0) {
                registry.fill(HIST("os_KK_rot_pT_0"), v01.M(), v01.Rapidity(), v01.Pt());
              }
              if (gapSide == 1) {
                registry.fill(HIST("os_KK_rot_pT_1"), v01.M(), v01.Rapidity(), v01.Pt());
              }
              if (gapSide == 2) {
                registry.fill(HIST("os_KK_rot_pT_2"), v01.M(), v01.Rapidity(), v01.Pt());
              }
            }
          }
        }
      }
    }
    for (auto& [t0, t1] : combinations(o2::soa::CombinationsFullIndexPolicy(tracks, tracks))) {
      if (!trackselector(t0, parameters) || !trackselector(t1, parameters))
        continue;
      if (t0.globalIndex() == t1.globalIndex())
        continue;
      if (rho && selectionPIDProton(t0, use_tof, nsigmatpc_cut, nsigmatof_cut) && selectionPIDPion1(t1)) {
        v0.SetXYZM(t0.px(), t0.py(), t0.pz(), mproton);
        v1.SetXYZM(t1.px(), t1.py(), t1.pz(), o2::constants::physics::MassPionCharged);
        v01 = v0 + v1;
        // Opposite sign pairs
        if (t0.sign() != t1.sign()) {
          if (gapSide == 0) {
            registry.fill(HIST("os_pp_pT_0"), v01.M(), v01.Rapidity(), v01.Pt());
          }
          if (gapSide == 1) {
            registry.fill(HIST("os_pp_pT_1"), v01.M(), v01.Rapidity(), v01.Pt());
          }
          if (gapSide == 2) {
            registry.fill(HIST("os_pp_pT_2"), v01.M(), v01.Rapidity(), v01.Pt());
          }
        } // same sign pair
        if (t0.sign() == t1.sign()) {
          if (gapSide == 0) {
            registry.fill(HIST("os_pp_ls_pT_0"), v01.M(), v01.Rapidity(), v01.Pt());
          }
          if (gapSide == 1) {
            registry.fill(HIST("os_pp_ls_pT_1"), v01.M(), v01.Rapidity(), v01.Pt());
          }
          if (gapSide == 2) {
            registry.fill(HIST("os_pp_ls_pT_2"), v01.M(), v01.Rapidity(), v01.Pt());
          }
        }
      }
      if (kstar && selectionPIDKaon1(t0) && selectionPIDPion1(t1)) {
        v0.SetXYZM(t0.px(), t0.py(), t0.pz(), o2::constants::physics::MassKaonCharged);
        v1.SetXYZM(t1.px(), t1.py(), t1.pz(), o2::constants::physics::MassPionCharged);
        v01 = v0 + v1;
        // Opposite sign pairs
        if (t0.sign() != t1.sign()) {
          if (gapSide == 0) {
            registry.fill(HIST("os_pk_pT_0"), v01.M(), v01.Rapidity(), v01.Pt());
          }
          if (gapSide == 1) {
            registry.fill(HIST("os_pk_pT_1"), v01.M(), v01.Rapidity(), v01.Pt());
          }
          if (gapSide == 2) {
            registry.fill(HIST("os_pk_pT_2"), v01.M(), v01.Rapidity(), v01.Pt());
          }
        } // same sign pair
        if (t0.sign() == t1.sign()) {
          if (gapSide == 0) {
            registry.fill(HIST("os_pk_ls_pT_0"), v01.M(), v01.Rapidity(), v01.Pt());
          }
          if (gapSide == 1) {
            registry.fill(HIST("os_pk_ls_pT_1"), v01.M(), v01.Rapidity(), v01.Pt());
          }
          if (gapSide == 2) {
            registry.fill(HIST("os_pk_ls_pT_2"), v01.M(), v01.Rapidity(), v01.Pt());
          }
        }
        if (fillRotation) {
          for (int nrotbkg = 0; nrotbkg < nBkgRotations; nrotbkg++) {
            auto anglestart = confMinRot;
            auto angleend = confMaxRot;
            auto anglestep = (angleend - anglestart) / (1.0 * (nBkgRotations - 1));
            auto rotangle = anglestart + nrotbkg * anglestep;
            registry.fill(HIST("hRotation"), rotangle);

            auto rotkaonPx = t0.px() * std::cos(rotangle) - t0.py() * std::sin(rotangle);
            auto rotkaonPy = t0.px() * std::sin(rotangle) + t0.py() * std::cos(rotangle);

            v0.SetXYZM(rotkaonPx, rotkaonPy, t0.pz(), o2::constants::physics::MassKaonCharged);
            v1.SetXYZM(t1.px(), t1.py(), t1.pz(), o2::constants::physics::MassPionCharged);
            v01 = v0 + v1;
            if (t0.sign() != t1.sign()) {
              if (gapSide == 0) {
                registry.fill(HIST("os_pk_rot_pT_0"), v01.M(), v01.Rapidity(), v01.Pt());
              }
              if (gapSide == 1) {
                registry.fill(HIST("os_pk_rot_pT_1"), v01.M(), v01.Rapidity(), v01.Pt());
              }
              if (gapSide == 2) {
                registry.fill(HIST("os_pk_rot_pT_2"), v01.M(), v01.Rapidity(), v01.Pt());
              }
            }
          }
        }
      }
    }
    if (fourpion) {
      if (gapSide == 2 && mult2 == 4) {

        ROOT::Math::PtEtaPhiMVector pair1, pair2, pair3, pair4;
        if (onlyPionTracks_p.size() == 2 && onlyPionTracks_n.size() == 2) {
          ROOT::Math::PtEtaPhiMVector k1 = onlyPionTracks_p.at(0);
          ROOT::Math::PtEtaPhiMVector k2 = onlyPionTracks_p.at(1);
          ROOT::Math::PtEtaPhiMVector k3 = onlyPionTracks_n.at(0);
          ROOT::Math::PtEtaPhiMVector k4 = onlyPionTracks_n.at(1);
          phiv = k1 + k2 + k3 + k4;
          pair1 = k1 + k3;
          pair2 = k2 + k4;
          pair3 = k1 + k4;
          pair4 = k2 + k3;
          registry.fill(HIST("os_pppp_pT_2"), phiv.M(), phiv.Pt(), phiv.Rapidity());
          registry.fill(HIST("os_pp_vs_pp_mass"), pair1.M(), pair2.M());
          registry.fill(HIST("os_pp_vs_pp_pt"), pair1.Pt(), pair2.Pt());
          auto costhetaPair = CosThetaCollinsSoperFrame(pair1, pair2, phiv);
          auto phiPair = 1. * TMath::Pi() + PhiCollinsSoperFrame(pair1, pair2, phiv);
          registry.fill(HIST("phi_dis"), phiPair);
          registry.fill(HIST("costheta_dis"), costhetaPair);
          registry.fill(HIST("costheta_vs_phi"), costhetaPair, phiPair);
          registry.fill(HIST("os_pp_vs_pp_mass1"), pair3.M(), pair4.M());
          registry.fill(HIST("os_pp_vs_pp_pt1"), pair3.Pt(), pair4.Pt());
          auto costhetaPair1 = CosThetaCollinsSoperFrame(pair3, pair4, phiv);
          auto phiPair1 = 1. * TMath::Pi() + PhiCollinsSoperFrame(pair3, pair4, phiv);
          registry.fill(HIST("phi_dis1"), phiPair1);
          registry.fill(HIST("costheta_dis1"), costhetaPair1);
          registry.fill(HIST("costheta_vs_phi1"), costhetaPair1, phiPair1);
        }
        if (onlyPionTracks_p.size() != 2 && onlyPionTracks_n.size() != 2) {
          if (onlyPionTracks_p.size() + onlyPionTracks_n.size() != 4)
            return;
          ROOT::Math::PtEtaPhiMVector l1 = onlyPionTracks_pm.at(0);
          ROOT::Math::PtEtaPhiMVector l2 = onlyPionTracks_pm.at(1);
          ROOT::Math::PtEtaPhiMVector l3 = onlyPionTracks_pm.at(2);
          ROOT::Math::PtEtaPhiMVector l4 = onlyPionTracks_pm.at(3);
          phiv1 = l1 + l2 + l3 + l4;
          registry.fill(HIST("os_pppp_pT_2_ls"), phiv1.M(), phiv1.Pt(), phiv1.Rapidity());
        }
      }
    }
  }
  PROCESS_SWITCH(SGResonanceAnalyzer, process, "Process unlike event", false);

  using UDCollisionsFull1 = soa::Join<aod::UDCollisions, aod::SGCollisions, aod::UDCollisionsSels, aod::UDZdcsReduced>; //
  SliceCache cache;
  Partition<udtracksfull> posTracks = aod::udtrack::sign > 0;
  Partition<udtracksfull> negTracks = aod::udtrack::sign < 0;

  ConfigurableAxis axisVertex{"axisVertex", {10, -10, 10}, "vertex axis for bin"};
  ConfigurableAxis axisMultiplicityClass{"axisMultiplicityClass", {10, 0, 100}, "multiplicity percentile for bin"};
  using BinningTypeVertexContributor = ColumnBinningPolicy<aod::collision::PosZ, aod::collision::NumContrib>;
  void mixprocess(UDCollisionsFull1 const& collisions, udtracksfull const& /*track*/)
  {
    TLorentzVector v0;
    TLorentzVector v1;
    TLorentzVector v01;
    float FIT_cut[5] = {FV0_cut, FT0A_cut, FT0C_cut, FDDA_cut, FDDC_cut};
    std::vector<float> parameters = {PV_cut, dcaZ_cut, dcaXY_cut, tpcChi2_cut, tpcNClsFindable_cut, itsChi2_cut, eta_cut, pt_cut};
    BinningTypeVertexContributor binningOnPositions{{axisVertex, axisMultiplicityClass}, true};
    for (auto const& [collision1, collision2] : o2::soa::selfCombinations(binningOnPositions, cfgNoMixedEvents, -1, collisions, collisions)) {
      int truegapSide1 = sgSelector.trueGap(collision1, FIT_cut[0], FIT_cut[1], FIT_cut[2], ZDC_cut);
      int truegapSide2 = sgSelector.trueGap(collision2, FIT_cut[0], FIT_cut[1], FIT_cut[2], ZDC_cut);
      if (truegapSide1 != truegapSide2)
        continue;
      if (truegapSide1 == -1)
        continue;
      auto posThisColl = posTracks->sliceByCached(aod::udtrack::udCollisionId, collision1.globalIndex(), cache);
      auto negThisColl = negTracks->sliceByCached(aod::udtrack::udCollisionId, collision2.globalIndex(), cache);
      //      for (auto& [track1, track2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(posThisColl, negThisColl))) {
      for (auto& [track1, track2] : o2::soa::combinations(posThisColl, negThisColl)) {
        if (!trackselector(track1, parameters) || !trackselector(track2, parameters))
          continue;
        if (selectionPIDKaon1(track1) && selectionPIDKaon1(track2)) {
          v0.SetXYZM(track1.px(), track1.py(), track1.pz(), o2::constants::physics::MassKaonCharged);
          v1.SetXYZM(track2.px(), track2.py(), track2.pz(), o2::constants::physics::MassKaonCharged);
          v01 = v0 + v1;
          // Opposite sign pairs
          if (track1.sign() != track2.sign()) {
            if (truegapSide1 == 0) {
              registry.fill(HIST("os_KK_mix_pT_0"), v01.M(), v01.Rapidity(), v01.Pt());
            }
            if (truegapSide1 == 1) {
              registry.fill(HIST("os_KK_mix_pT_1"), v01.M(), v01.Rapidity(), v01.Pt());
            }
            if (truegapSide1 == 2) {
              registry.fill(HIST("os_KK_mix_pT_2"), v01.M(), v01.Rapidity(), v01.Pt());
            }
          }
        }
      }
      for (auto& [track1, track2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(posThisColl, negThisColl))) {
        if (!trackselector(track1, parameters) || !trackselector(track2, parameters))
          continue;
        if (track1.globalIndex() == track2.globalIndex())
          continue;
        if (selectionPIDKaon1(track1) && selectionPIDPion1(track2)) {
          v0.SetXYZM(track1.px(), track1.py(), track1.pz(), o2::constants::physics::MassKaonCharged);
          v1.SetXYZM(track2.px(), track2.py(), track2.pz(), o2::constants::physics::MassPionCharged);
          v01 = v0 + v1;
          // Opposite sign pairs
          if (track1.sign() != track2.sign()) {
            if (truegapSide1 == 0) {
              registry.fill(HIST("os_pk_mix_pT_0"), v01.M(), v01.Rapidity(), v01.Pt());
            }
            if (truegapSide1 == 1) {
              registry.fill(HIST("os_pk_mix_pT_1"), v01.M(), v01.Rapidity(), v01.Pt());
            }
            if (truegapSide1 == 2) {
              registry.fill(HIST("os_pk_mix_pT_2"), v01.M(), v01.Rapidity(), v01.Pt());
            }
          }
        }
      }
    }
  }
  PROCESS_SWITCH(SGResonanceAnalyzer, mixprocess, "Process Mixed event", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<SGResonanceAnalyzer>(cfgc)};
}
