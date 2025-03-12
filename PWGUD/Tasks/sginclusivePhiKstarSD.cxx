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
#include "TDatabasePDG.h"

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
  // a pdg object
  TDatabasePDG* pdg = nullptr;
  // get a DGCutparHolder
  Configurable<int> verbosity{"Verbosity", 0, "Determines level of verbosity"};

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
  Configurable<float> nsigmatpc_cut4{"nsigmatpc4", 3.0, "nsigma tpc cut4"};
  Configurable<float> nsigmatpc_cut{"nsigmatpc", 3.0, "nsigma tpc cut"};
  Configurable<float> nsigmatof_cut{"nsigmatof", 9.0, "nsigma tpc+tof cut"};
  Configurable<float> nsigmatof_cut1{"nsigmatof1", 3.0, "nsigma tof cut"};
  Configurable<float> pionnsigmacut{"pionnsigmacut", 3.0, "nsigma tpc cut for kaon"};

  Configurable<int> mintrack{"min_track", 1, "min track"};
  Configurable<int> maxtrack{"max_track", 50, "max track"};
  Configurable<bool> use_tof{"Use_TOF", true, "TOF PID"};
  Configurable<bool> ccut{"ccut", true, "TPC + TOF PID"};
  Configurable<bool> kaoncut{"kaoncut", true, " kaon slection cut for kstar "};

  Configurable<bool> QA{"QA", true, ""};
  Configurable<bool> rapidity_gap{"rapidity_gap", true, ""};

  Configurable<bool> phi{"phi", true, ""};
  Configurable<bool> rho{"rho", true, ""};
  Configurable<bool> kstar{"kstar", true, ""};
  Configurable<bool> fourpion{"fourpion", true, ""};
  Configurable<bool> MC{"MC", true, ""};
  Configurable<int> gapsideMC{"gapsideMC", 1, "gapside MC"};

  Configurable<int> cfgNoMixedEvents{"cfgNoMixedEvents", 1, "Number of mixed events per event"};
  Configurable<int> nBkgRotations{"nBkgRotations", 9, "Number of rotated copies (background) per each original candidate"};
  Configurable<bool> fillRotation{"fillRotation", true, "fill rotation"};
  Configurable<float> confMinRot{"confMinRot", 5.0 * TMath::Pi() / 6.0, "Minimum of rotation"};
  Configurable<float> confMaxRot{"confMaxRot", 7.0 * TMath::Pi() / 6.0, "Maximum of rotation"};

  void init(InitContext const& context)
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
      registry.add("tof_nsigma_kaon_f", "p#k n#sigma", kTH2F, {{100, 0.0, 10.0}, {100, -10.0, 10.0}});
      registry.add("tof_nsigma_pion_f", "p#pi n#sigma", kTH2F, {{100, 0.0, 10.0}, {100, -10.0, 10.0}});
      registry.add("tpc_nsigma_kaon_f", "p#k n#sigma", kTH2F, {{100, 0.0, 10.0}, {100, -10.0, 10.0}});
      registry.add("tpc_nsigma_pion_f", "p#pi n#sigma", kTH2F, {{100, 0.0, 10.0}, {100, -10.0, 10.0}});

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
    if (MC) {
      pdg = TDatabasePDG::Instance();
      // add histograms for the different process functions
      if (context.mOptions.get<bool>("processMCTruth")) {
        registry.add("MC/Stat", "Count generated events; ; Entries", {HistType::kTH1F, {{1, -0.5, 0.5}}});
        registry.add("MC/recCols", "Number of reconstructed collisions; Number of reconstructed collisions; Entries", {HistType::kTH1F, {{31, -0.5, 30.5}}});
        registry.add("MC/nParts", "Number of McParticles per collision; Number of McParticles; Entries", {HistType::kTH1F, {{1001, -0.5, 1000.5}}});
        registry.add("MC/nRecTracks", "Number of reconstructed tracks per McParticle; Number of reconstructed tracks per McParticle; Entries", {HistType::kTH1F, {{11, -0.5, 10.5}}});
        registry.add("MC/genEtaPt", "Generated events; eta (1); Pt (GeV/c)", {HistType::kTH2F, {{60, -1.5, 1.5}, {250, 0.0, 5.0}}});
        registry.add("MC/genRap", "Generated events; Rapidity (1)", {HistType::kTH1F, {{60, -1.5, 1.5}}});
        registry.add("MC/genMPt", "Generated events; Mass (GeV/c^2); Pt (GeV/c)", {HistType::kTH2F, {{220, 0.98, 1.2}, {200, 0.0, 10.0}}});
        registry.add("MC/genM", "Generated events; Mass (GeV/c^2)", {HistType::kTH1F, {
                                                                                        {220, 0.98, 1.2},
                                                                                      }});
        registry.add("MC/genM_1", "Generated events; Mass (GeV/c^2)", {HistType::kTH1F, {
                                                                                          {220, 0.98, 1.2},
                                                                                        }});

        registry.add("MC/accMPtRap_phi_G", "Generated Phi; Mass (GeV/c^2); Pt (GeV/c)", {HistType::kTH3F, {{220, 0.98, 1.20}, {200, 0.0, 10.0}, {60, -1.5, 1.5}}});
        registry.add("MC/accEtaPt", "Generated events in acceptance; eta (1); Pt (GeV/c)", {HistType::kTH2F, {{60, -1.5, 1.5}, {250, 0.0, 5.0}}});
        registry.add("MC/accRap", "Generated events in acceptance; Rapidity (1)", {HistType::kTH1F, {{60, -1.5, 1.5}}});
        registry.add("MC/accMPt", "Generated events in acceptance; Mass (GeV/c^2); Pt (GeV/c)", {HistType::kTH2F, {{220, 0.98, 1.20}, {200, 0.0, 10.0}}});
        registry.add("MC/accMPtRap", "Generated events in acceptance; Mass (GeV/c^2); Pt (GeV/c)", {HistType::kTH3F, {{220, 0.98, 1.20}, {200, 0.0, 10.0}, {60, -1.5, 1.5}}});
        registry.add("MC/accM", "Generated events in acceptance; Mass (GeV/c^2)", {HistType::kTH1F, {{220, 0.98, 1.20}}});
        registry.add("MC/selEtaPt", "Selected events in acceptance; eta (1); Pt (GeV/c)", {HistType::kTH2F, {{60, -1.5, 1.5}, {250, 0.0, 5.0}}});
        registry.add("MC/selRap", "Selected events in acceptance; Rapidity (1)", {HistType::kTH1F, {{60, -1.5, 1.5}}});
        registry.add("MC/selMPt", "Selected events in acceptance; Mass (GeV/c^2); Pt (GeV/c)", {HistType::kTH2F, {{250, 2.5, 5.0}, {100, 0.0, 1.0}}});
        registry.add("MC/pDiff", "McTruth - reconstructed track momentum; McTruth - reconstructed track momentum; Entries", {HistType::kTH2F, {{240, -6., 6.}, {3, -1.5, 1.5}}});
        // K*0
        registry.add("MC/accMPtRap_kstar_G", "Generated K*0; Mass (GeV/c^2); Pt (GeV/c)", {HistType::kTH3F, {{400, 0., 2.0}, {200, 0.0, 10.0}, {60, -1.5, 1.5}}});
        registry.add("MC/genEtaPt_k", "Generated events; eta (1); Pt (GeV/c)", {HistType::kTH2F, {{60, -1.5, 1.5}, {250, 0.0, 5.0}}});
        registry.add("MC/genRap_k", "Generated events; Rapidity (1)", {HistType::kTH1F, {{60, -1.5, 1.5}}});
        registry.add("MC/genMPt_k", "Generated events; Mass (GeV/c^2); Pt (GeV/c)", {HistType::kTH2F, {{400, 0., 2.0}, {200, 0.0, 10.0}}});
        registry.add("MC/genM_k", "Generated events; Mass (GeV/c^2)", {HistType::kTH1F, {
                                                                                          {400, 0., 2.0},
                                                                                        }});
        registry.add("MC/genM_1_k", "Generated events; Mass (GeV/c^2)", {HistType::kTH1F, {
                                                                                            {400, 0., 2.0},
                                                                                          }});

        registry.add("MC/accEtaPt_k", "Generated events in acceptance; eta (1); Pt (GeV/c)", {HistType::kTH2F, {{60, -1.5, 1.5}, {250, 0.0, 5.0}}});
        registry.add("MC/accRap_k", "Generated events in acceptance; Rapidity (1)", {HistType::kTH1F, {{60, -1.5, 1.5}}});
        registry.add("MC/accMPt_k", "Generated events in acceptance; Mass (GeV/c^2); Pt (GeV/c)", {HistType::kTH2F, {{400, 0., 2.0}, {200, 0.0, 10.0}}});
        registry.add("MC/accMPtRap_k", "Generated events in acceptance; Mass (GeV/c^2); Pt (GeV/c)", {HistType::kTH3F, {{400, 0., 2.0}, {200, 0.0, 10.0}, {60, -1.5, 1.5}}});
        registry.add("MC/accM_k", "Generated events in acceptance; Mass (GeV/c^2)", {HistType::kTH1F, {{400, 0., 2.0}}});
      }
      if (context.mOptions.get<bool>("processReco")) {
        registry.add("Reco/Stat", "Count reconstruted events; ; Entries", {HistType::kTH1F, {{5, -0.5, 4.5}}});
        // registry.add("Reco/nTracks", "Number of reconstructed tracks per collision; Number of reconstructed tracks; Entries", {HistType::kTH1F, {{51, -0.5, 50.5}}});
        registry.add("Reco/nPVContributors", "Number of PV contributors per collision; Number of PV contributors; Entries", {HistType::kTH1F, {{51, -0.5, 50.5}}});
        registry.add("Reco/selEtaPt", "Selected events in acceptance; eta (1); Pt (GeV/c)", {HistType::kTH2F, {{300, -1.5, 1.5}, {250, 0.0, 5.0}}});
        registry.add("Reco/selRap", "Selected events in acceptance; Rapidity (1)", {HistType::kTH1F, {{60, -1.5, 1.5}}});
        registry.add("Reco/selMPt", "Reconstructed events in acceptance; Mass (GeV/c^2); Pt (GeV/c)", {HistType::kTH2F, {{220, 0.98, 1.20}, {200, 0.0, 10.0}}});
        registry.add("Reco/selMPtRap", "Reconstructed events in acceptance; Mass (GeV/c^2); Pt (GeV/c)", {HistType::kTH3F, {{220, 0.98, 1.20}, {200, 0.0, 10.0}, {60, -1.5, 1.5}}});
        registry.add("Reco/selMPtRap_gen", "Reconstructed(gen) events in acceptance; Mass (GeV/c^2); Pt (GeV/c)", {HistType::kTH3F, {{220, 0.98, 1.20}, {200, 0.0, 10.0}, {60, -1.5, 1.5}}});
        registry.add("Reco/selPt", "Reconstructed events in acceptance;Pt (GeV/c)", {HistType::kTH1F, {{200, 0.0, 10.0}}});
        registry.add("Reco/selM", "Reconstructed events in acceptance; Mass (GeV/c^2); ", {HistType::kTH1F, {{220, 0.98, 1.20}}});
        registry.add("Reco/mcEtaPt", "Generated events in acceptance; eta (1); Pt (GeV/c)", {HistType::kTH2F, {{60, -1.5, 1.5}, {250, 0.0, 5.0}}});
        registry.add("Reco/mcRap", "Generated events in acceptance; Rapidity (1)", {HistType::kTH1F, {{60, -1.5, 1.5}}});
        registry.add("Reco/mcMPt", "Generated events in acceptance; Mass (GeV/c^2); Pt (GeV/c)", {HistType::kTH2F, {{250, 2.5, 5.0}, {100, 0.0, 1.0}}});
        registry.add("Reco/pDiff", "McTruth - reconstructed track momentum; McTruth - reconstructed track momentum; Entries", {HistType::kTH2F, {{240, -6., 6.}, {3, -1.5, 1.5}}});

        registry.add("Reco/selEtaPt_k", "Selected events in acceptance; eta (1); Pt (GeV/c)", {HistType::kTH2F, {{300, -1.5, 1.5}, {250, 0.0, 5.0}}});
        registry.add("Reco/selRap_k", "Selected events in acceptance; Rapidity (1)", {HistType::kTH1F, {{60, -1.5, 1.5}}});
        registry.add("Reco/selMPt_k", "Reconstructed events in acceptance; Mass (GeV/c^2); Pt (GeV/c)", {HistType::kTH2F, {{400, 0., 2.0}, {200, 0.0, 10.0}}});
        registry.add("Reco/selMPtRap_k", "Reconstructed events in acceptance; Mass (GeV/c^2); Pt (GeV/c)", {HistType::kTH3F, {{400, 0., 2.0}, {200, 0.0, 10.0}, {60, -1.5, 1.5}}});
        registry.add("Reco/selMPtRap_k_gen", "Reconstructed(gen) events in acceptance; Mass (GeV/c^2); Pt (GeV/c)", {HistType::kTH3F, {{400, 0., 2.0}, {200, 0.0, 10.0}, {60, -1.5, 1.5}}});
        registry.add("Reco/selPt_k", "Reconstructed events in acceptance;Pt (GeV/c)", {HistType::kTH1F, {{200, 0.0, 10.0}}});
        registry.add("Reco/selM_k", "Reconstructed events in acceptance; Mass (GeV/c^2); ", {HistType::kTH1F, {{400, 0., 2.0}}});
        registry.add("Reco/selM_k_K", "Reconstructed events in acceptance; Mass (GeV/c^2); ", {HistType::kTH1F, {{400, 0., 2.0}}});

        registry.add("Reco/nTracks", "Number of reconstructed tracks per collision; Number of reconstructed tracks; Entries", {HistType::kTH1F, {{101, -0.5, 100.5}}});
        registry.add("Reco/treta_k", "track kaon eta", {HistType::kTH1F, {{200, -5.0, 5.0}}});
        registry.add("Reco/trpt_k", "rec kaon track pt", {HistType::kTH1F, {{200, 0.0, 10.0}}});
        registry.add("Reco/trpt", "rec track pt", {HistType::kTH1F, {{200, 0.0, 10.0}}});

        registry.add("Reco/tr_dcaz_1", "dcaz-", {HistType::kTH1F, {{1000, -5.0, 5.0}}});
        registry.add("Reco/tr_dcaxy_1", "dcaxy-", {HistType::kTH1F, {{1000, -5.0, 5.0}}});
        registry.add("Reco/tr_chi2ncl_1", "chi2ncl-", {HistType::kTH1F, {{100, 0.0, 100.0}}});
        registry.add("Reco/tr_tpcnclfind_1", "tpcnclfind-", {HistType::kTH1F, {{300, 0.0, 300.0}}});
        registry.add("Reco/tr_itsChi2NCl_1", "itsChi2NCl-", {HistType::kTH1F, {{200, 0.0, 200.0}}});
        registry.add("Reco/tr_Eta_1", " eta -", {HistType::kTH1F, {{300, 1.5, 1.5}}});

        registry.add("Reco/tr_dcaz_2", "dcaz", {HistType::kTH1F, {{1000, -5.0, 5.0}}});
        registry.add("Reco/tr_dcaxy_2", "dcaxy", {HistType::kTH1F, {{1000, -5.0, 5.0}}});
        registry.add("Reco/tr_chi2ncl_2", "chi2ncl", {HistType::kTH1F, {{100, 0.0, 100.0}}});
        registry.add("Reco/tr_tpcnclfind_2", "tpcnclfind", {HistType::kTH1F, {{300, 0.0, 300.0}}});
        registry.add("Reco/tr_itsChi2NCl_2", "itsChi2NCl", {HistType::kTH1F, {{200, 0.0, 200.0}}});
      }
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
    if (use_tof && candidate.hasTOF() && ccut && (candidate.tofNSigmaKa() * candidate.tofNSigmaKa() + candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa()) < nsigmatof_cut) {
      return true;
    }
    if (use_tof && candidate.hasTOF() && !ccut && std::abs(candidate.tpcNSigmaKa()) < nsigmatpc_cut4 && std::abs(candidate.tofNSigmaKa()) < nsigmatof_cut1) {
      return true;
    }

    if (!use_tof && std::abs(candidate.tpcNSigmaKa()) < nsigmatpc_cut) {
      return true;
    }
    return false;
  }

  float particleMass(TDatabasePDG* pdg, int pid)
  {
    auto mass = 0.;
    TParticlePDG* pdgparticle = pdg->GetParticle(pid);
    if (pdgparticle != nullptr) {
      mass = pdgparticle->Mass();
    }
    return mass;
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
    if (use_tof && candidate.hasTOF() && ccut && (candidate.tofNSigmaPi() * candidate.tofNSigmaPi() + candidate.tpcNSigmaPi() * candidate.tpcNSigmaPi()) < nsigmatof_cut) {
      return true;
    }
    if (use_tof && candidate.hasTOF() && !ccut && std::abs(candidate.tpcNSigmaKa()) < nsigmatpc_cut4 && std::abs(candidate.tofNSigmaKa()) < nsigmatof_cut1) {
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
        registry.fill(HIST("tof_nsigma_kaon_f"), v0.Pt(), track1.tofNSigmaKa());
        registry.fill(HIST("tof_nsigma_pion_f"), v0.Pt(), track1.tofNSigmaPi());
        registry.fill(HIST("tpc_nsigma_kaon_f"), v0.Pt(), track1.tpcNSigmaKa());
        registry.fill(HIST("tpc_nsigma_pion_f"), v0.Pt(), track1.tpcNSigmaPi());

        if (selectionPIDKaon1(track1)) {
          registry.fill(HIST("tpc_dedx_kaon_1"), v0.P(), track1.tpcSignal());
          registry.fill(HIST("tpc_nsigma_kaon"), v0.Pt(), track1.tpcNSigmaKa());
          registry.fill(HIST("tof_nsigma_kaon"), v0.Pt(), track1.tofNSigmaKa());
          registry.fill(HIST("tpc_tof_nsigma_kaon"), track1.tpcNSigmaKa(), track1.tofNSigmaKa());
        }

        if (selectionPIDPion1(track1)) {
          registry.fill(HIST("tpc_dedx_pion_1"), v0.P(), track1.tpcSignal());
          registry.fill(HIST("tpc_nsigma_pion"), v0.Pt(), track1.tpcNSigmaPi());
          registry.fill(HIST("tof_nsigma_pion"), v0.Pt(), track1.tofNSigmaPi());
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
        if (kaoncut && t0.tpcNSigmaPi() < pionnsigmacut)
          continue;
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

  // define abbreviations  , aod::UDCollisions_001,
  using CCs = soa::Join<aod::UDCollisions, aod::SGCollisions, aod::UDCollisionSelExtras, aod::UDCollisionsSels, aod::UDMcCollsLabels, aod::UDZdcsReduced>;
  using CC = CCs::iterator;
  // using TCs = soa::Join<aod::UDTracks>;
  using TCs = soa::Join<aod::UDTracks, aod::UDTracksExtra, aod::UDTracksDCA, aod::UDTracksFlags, aod::UDTracksPID, aod::UDMcTrackLabels>;
  using TC = TCs::iterator;

  PresliceUnsorted<aod::UDMcParticles> partPerMcCollision = aod::udmcparticle::udMcCollisionId;
  PresliceUnsorted<CCs> colPerMcCollision = aod::udcollision::udMcCollisionId;
  PresliceUnsorted<TCs> trackPerMcParticle = aod::udmctracklabel::udMcParticleId;

  void processMCTruth(aod::UDMcCollisions const& mccollisions, CCs const& collisions, aod::UDMcParticles const& McParts, TCs const& tracks)
  {
    // number of McCollisions in DF
    if (verbosity > 0) {
      LOGF(info, "Number of MC collisions %d", mccollisions.size());
    }
    bool Kstar = true;
    TLorentzVector v0;
    TLorentzVector v1;
    TLorentzVector v01;
    TLorentzVector v_kstar;
    TLorentzVector v_phi;

    // some variables
    TLorentzVector* lv1_rec = new TLorentzVector();
    TLorentzVector* lv2_rec = new TLorentzVector();
    TLorentzVector* lv_rec = new TLorentzVector();
    int y = 0;
    // loop over all generated collisions
    int kk1 = 0;
    for (auto mccollision : mccollisions) {
      registry.get<TH1>(HIST("MC/Stat"))->Fill(0., 1.);
      y++;
      // get reconstructed collision which belongs to mccollision
      auto colSlice = collisions.sliceBy(colPerMcCollision, mccollision.globalIndex());
      registry.get<TH1>(HIST("MC/recCols"))->Fill(colSlice.size(), 1.);

      // get McParticles which belong to mccollision
      auto partSlice = McParts.sliceBy(partPerMcCollision, mccollision.globalIndex());
      registry.get<TH1>(HIST("MC/nParts"))->Fill(partSlice.size(), 1.);
      if (verbosity > 0) {
        LOGF(info, "Number of McParts %d", partSlice.size());
      }
      for (auto& [tr1, tr2] : combinations(partSlice, partSlice)) {
        if ((tr1.pdgCode() == 321 && tr2.pdgCode() == -211) || (tr1.pdgCode() == -321 && tr2.pdgCode() == 211) || (tr1.pdgCode() == 211 && tr2.pdgCode() == -321) || (tr1.pdgCode() == -211 && tr2.pdgCode() == 321)) {
          auto m1 = particleMass(pdg, 321);
          auto m2 = particleMass(pdg, 211);
          auto m3 = particleMass(pdg, 313);
          if (std::abs(tr1.pdgCode()) == 321) {
            v0.SetXYZM(tr1.px(), tr1.py(), tr1.pz(), m1);
            v1.SetXYZM(tr2.px(), tr2.py(), tr2.pz(), m2);
          } else {
            v0.SetXYZM(tr1.px(), tr1.py(), tr1.pz(), m2);
            v1.SetXYZM(tr2.px(), tr2.py(), tr2.pz(), m1);
          }
          if (!tr1.isPhysicalPrimary() || !tr2.isPhysicalPrimary())
            continue;
          v01 = v0 + v1;
          if (tr1.globalIndex() + 1 != tr2.globalIndex()) {
            registry.get<TH1>(HIST("MC/genM_1_k"))->Fill(v01.M(), 1.);
          }
          if (std::abs(tr1.globalIndex() - tr2.globalIndex()) != 1)
            continue;
          bool flag = false;
          bool flag1 = false;
          if (tr1.has_mothers() && tr2.has_mothers()) {
            for (const auto& mother : tr1.mothers_as<aod::UDMcParticles>()) {
              if (std::abs(mother.pdgCode()) == 313) {
                v_kstar.SetXYZM(mother.px(), mother.py(), mother.pz(), m3);
                registry.get<TH3>(HIST("MC/accMPtRap_kstar_G"))->Fill(v_kstar.M(), v_kstar.Pt(), v_kstar.Rapidity(), 1.);
                flag = true;
              }
            }
            for (const auto& mother1 : tr2.mothers_as<aod::UDMcParticles>()) {
              if (std::abs(mother1.pdgCode()) == 313) {
                flag1 = true;
              }
            }
          }
          if (flag && flag1) {
            registry.get<TH1>(HIST("MC/genRap_k"))->Fill(v01.Rapidity(), 1.);
            registry.get<TH2>(HIST("MC/genMPt_k"))->Fill(v01.M(), v01.Pt(), 1.);
            registry.get<TH1>(HIST("MC/genM_k"))->Fill(v01.M(), 1.);
            if (std::abs(v0.Eta()) < 0.8 && std::abs(v1.Eta()) < 0.8 && v0.Pt() > 0.15 && v1.Pt() > 0.15) {
              registry.get<TH1>(HIST("MC/accRap_k"))->Fill(v01.Rapidity(), 1.);
              registry.get<TH2>(HIST("MC/accMPt_k"))->Fill(v01.M(), v01.Pt(), 1.);
              registry.get<TH1>(HIST("MC/accM_k"))->Fill(v01.M(), 1.);
              registry.get<TH3>(HIST("MC/accMPtRap_k"))->Fill(v01.M(), v01.Pt(), v01.Rapidity(), 1.);
            }
          }
        }
        if (std::abs(tr1.pdgCode()) != 321 || std::abs(tr2.pdgCode()) != 321)
          continue;
        auto m1 = particleMass(pdg, 321);
        v0.SetXYZM(tr1.px(), tr1.py(), tr1.pz(), m1);
        v1.SetXYZM(tr2.px(), tr2.py(), tr2.pz(), m1);
        if (tr1.pdgCode() == tr2.pdgCode())
          continue;
        v01 = v0 + v1;
        if (!tr1.isPhysicalPrimary() || !tr2.isPhysicalPrimary())
          continue;
        if (tr1.globalIndex() + 1 != tr2.globalIndex()) {
          registry.get<TH1>(HIST("MC/genM_1"))->Fill(v01.M(), 1.);
        }
        if (std::abs(tr1.globalIndex() - tr2.globalIndex()) != 1)
          continue;
        bool flag = false;
        bool flag1 = false;
        if (tr1.has_mothers() && tr2.has_mothers()) {
          for (const auto& mother : tr1.mothers_as<aod::UDMcParticles>()) {
            auto m4 = particleMass(pdg, 333);
            if (std::abs(mother.pdgCode()) == 333) {
              v_phi.SetXYZM(mother.px(), mother.py(), mother.pz(), m4);
              registry.get<TH3>(HIST("MC/accMPtRap_phi_G"))->Fill(v_phi.M(), v_phi.Pt(), v_phi.Rapidity(), 1.);
              flag = true;
            }
          }
          for (const auto& mother1 : tr2.mothers_as<aod::UDMcParticles>()) {
            if (std::abs(mother1.pdgCode()) == 333) {
              flag1 = true;
            }
          }
        }
        if (flag && flag1) {
          registry.get<TH1>(HIST("MC/genRap"))->Fill(v01.Rapidity(), 1.);
          registry.get<TH2>(HIST("MC/genMPt"))->Fill(v01.M(), v01.Pt(), 1.);
          registry.get<TH1>(HIST("MC/genM"))->Fill(v01.M(), 1.);
          if (std::abs(v0.Eta()) < 0.8 && std::abs(v1.Eta()) < 0.8 && v0.Pt() > 0.15 && v1.Pt() > 0.15) {
            registry.get<TH1>(HIST("MC/accRap"))->Fill(v01.Rapidity(), 1.);
            registry.get<TH2>(HIST("MC/accMPt"))->Fill(v01.M(), v01.Pt(), 1.);
            registry.get<TH1>(HIST("MC/accM"))->Fill(v01.M(), 1.);
            registry.get<TH3>(HIST("MC/accMPtRap"))->Fill(v01.M(), v01.Pt(), v01.Rapidity(), 1.);
          }
        }
      }
      // compute the difference between generated and reconstructed particle momentum
      for (auto McPart : partSlice) {
        // get track which corresponds to McPart
        auto trackSlice = tracks.sliceBy(trackPerMcParticle, McPart.globalIndex());
        registry.get<TH1>(HIST("MC/nRecTracks"))->Fill(trackSlice.size(), 1.);
        // compute momentum difference between MCTruth and Reconstruction
        if (trackSlice.size() > 0) {
          for (auto track : trackSlice) {
            // std::cout<<trackSlice.size()<<std::endl;
            auto pTrack = sqrt(track.px() * track.px() + track.py() * track.py() + track.pz() * track.pz());
            auto pPart = sqrt(McPart.px() * McPart.px() + McPart.py() * McPart.py() + McPart.pz() * McPart.pz());
            auto pDiff = pTrack - pPart;
            registry.get<TH2>(HIST("MC/pDiff"))->Fill(pDiff, track.isPVContributor(), 1.);
            if (verbosity > 0) {
              LOGF(info, "  PID: %d Generated: %d Process: %d PV contributor: %d dP: %f", McPart.pdgCode(), McPart.producedByGenerator(), McPart.getProcess(), track.isPVContributor(), pDiff);
            }
          }
        } else {
          registry.get<TH2>(HIST("MC/pDiff"))->Fill(-5.9, -1, 1.);
          if (verbosity > 0) {
            LOGF(info, "  PID: %d Generated: %d Process: %d PV contributor: No dP: nan", McPart.pdgCode(), McPart.producedByGenerator(), McPart.getProcess());
          }
        }
      }
      if (verbosity > 0) {
        LOGF(info, "");
      }
    }
  }
  PROCESS_SWITCH(SGResonanceAnalyzer, processMCTruth, "Process MC truth", true);
  // ...............................................................................................................
  void processReco(CC const& collision, TCs const& tracks, aod::UDMcCollisions const& /*mccollisions*/, aod::UDMcParticles const& McParts)
  {
    // number of McCollisions in DF
    if (verbosity > 0) {
      LOGF(info, "Number of MC collisions %d", collision.size());
    }
    float FIT_cut[5] = {FV0_cut, FT0A_cut, FT0C_cut, FDDA_cut, FDDC_cut};
    std::vector<float> parameters = {PV_cut, dcaZ_cut, dcaXY_cut, tpcChi2_cut, tpcNClsFindable_cut, itsChi2_cut, eta_cut, pt_cut};
    int truegapSide = sgSelector.trueGap(collision, FIT_cut[0], FIT_cut[1], FIT_cut[2], ZDC_cut);
    registry.get<TH1>(HIST("Reco/Stat"))->Fill(4.0, 1.);
    Partition<TCs> PVContributors = aod::udtrack::isPVContributor == true;
    PVContributors.bindTable(tracks);
    if (std::abs(collision.posZ()) > 10.0)
      return;
    if (std::abs(collision.occupancyInTime()) > 1000.0)
      return;
    registry.get<TH1>(HIST("Reco/Stat"))->Fill(truegapSide, 1.);
    if (truegapSide != gapsideMC)
      return;
    registry.get<TH1>(HIST("Reco/nPVContributors"))->Fill(PVContributors.size(), 1.);
    TLorentzVector v_phi;
    TLorentzVector v_kstar;
    TLorentzVector v0;
    TLorentzVector vr0;
    TLorentzVector vr1;
    TLorentzVector vr01;
    TLorentzVector vr0_g;
    TLorentzVector vr1_g;
    TLorentzVector vr01_g;
    int t1 = 0;
    if (QA) {
      if (truegapSide == 0) {
        registry.fill(HIST("V0A_0"), collision.totalFV0AmplitudeA());
        registry.fill(HIST("FT0A_0"), collision.totalFT0AmplitudeA());
        registry.fill(HIST("FT0C_0"), collision.totalFT0AmplitudeC());
        registry.fill(HIST("ZDC_A_0"), collision.energyCommonZNA());
        registry.fill(HIST("ZDC_C_0"), collision.energyCommonZNC());
      }
      if (truegapSide == 1) {
        registry.fill(HIST("V0A_1"), collision.totalFV0AmplitudeA());
        registry.fill(HIST("FT0A_1"), collision.totalFT0AmplitudeA());
        registry.fill(HIST("FT0C_1"), collision.totalFT0AmplitudeC());
        registry.fill(HIST("ZDC_A_1"), collision.energyCommonZNA());
        registry.fill(HIST("ZDC_C_1"), collision.energyCommonZNC());
      }
    }
    for (auto tr1 : tracks) {
      if (!tr1.has_udMcParticle())
        continue;
      auto McPart_1 = tr1.udMcParticle();
      registry.get<TH1>(HIST("Reco/tr_dcaz_1"))->Fill(tr1.dcaZ(), 1.);
      registry.get<TH1>(HIST("Reco/tr_dcaxy_1"))->Fill(tr1.dcaXY(), 1.);
      registry.get<TH1>(HIST("Reco/tr_chi2ncl_1"))->Fill(tr1.tpcChi2NCl(), 1.);
      registry.get<TH1>(HIST("Reco/tr_tpcnclfind_1"))->Fill(tr1.tpcNClsFindable(), 1.);
      registry.get<TH1>(HIST("Reco/tr_itsChi2NCl_1"))->Fill(tr1.itsChi2NCl(), 1.);

      if (!trackselector(tr1, parameters))
        continue;

      registry.get<TH1>(HIST("Reco/tr_dcaz_2"))->Fill(tr1.dcaZ(), 1.);
      registry.get<TH1>(HIST("Reco/tr_dcaxy_2"))->Fill(tr1.dcaXY(), 1.);
      registry.get<TH1>(HIST("Reco/tr_chi2ncl_2"))->Fill(tr1.tpcChi2NCl(), 1.);
      registry.get<TH1>(HIST("Reco/tr_tpcnclfind_2"))->Fill(tr1.tpcNClsFindable(), 1.);
      registry.get<TH1>(HIST("Reco/tr_itsChi2NCl_2"))->Fill(tr1.itsChi2NCl(), 1.);
      v0.SetXYZM(tr1.px(), tr1.py(), tr1.pz(), o2::constants::physics::MassPionCharged);
      if (QA) {
        registry.fill(HIST("tpc_dedx"), v0.P(), tr1.tpcSignal());
        registry.fill(HIST("tof_beta"), v0.P(), tr1.beta());
        registry.fill(HIST("tof_nsigma_kaon_f"), v0.Pt(), tr1.tofNSigmaKa());
        registry.fill(HIST("tof_nsigma_pion_f"), v0.Pt(), tr1.tofNSigmaPi());
        registry.fill(HIST("tpc_nsigma_kaon_f"), v0.Pt(), tr1.tpcNSigmaKa());
        registry.fill(HIST("tpc_nsigma_pion_f"), v0.Pt(), tr1.tpcNSigmaPi());
        if (selectionPIDKaon1(tr1)) {
          registry.fill(HIST("tpc_dedx_kaon_1"), v0.P(), tr1.tpcSignal());
          registry.fill(HIST("tpc_nsigma_kaon"), v0.Pt(), tr1.tpcNSigmaKa());
          registry.fill(HIST("tof_nsigma_kaon"), v0.Pt(), tr1.tofNSigmaKa());
          registry.fill(HIST("tpc_tof_nsigma_kaon"), tr1.tpcNSigmaKa(), tr1.tofNSigmaKa());
        }
        if (selectionPIDPion1(tr1)) {
          registry.fill(HIST("tpc_dedx_pion_1"), v0.P(), tr1.tpcSignal());
          registry.fill(HIST("tpc_nsigma_pion"), v0.Pt(), tr1.tpcNSigmaPi());
          registry.fill(HIST("tof_nsigma_pion"), v0.Pt(), tr1.tofNSigmaPi());
          registry.fill(HIST("tpc_tof_nsigma_pion"), tr1.tpcNSigmaPi(), tr1.tofNSigmaPi());
        }
      }
      t1++;
      auto m1 = particleMass(pdg, 321);
      auto m4 = particleMass(pdg, 333);
      vr0.SetXYZM(tr1.px(), tr1.py(), tr1.pz(), m1);
      registry.get<TH1>(HIST("Reco/trpt"))->Fill(vr0.Pt(), 1.);
      registry.get<TH1>(HIST("Reco/treta_k"))->Fill(vr0.Eta(), 1.);
      if (!selectionPIDKaon1(tr1))
        continue;
      registry.get<TH1>(HIST("Reco/trpt_k"))->Fill(vr0.Pt(), 1.);
      int t2 = 0;
      for (auto tr2 : tracks) {
        if (!tr2.has_udMcParticle())
          continue;
        if (!trackselector(tr2, parameters))
          continue;
        t2++;
        if (t2 > t1) {
          if (!selectionPIDKaon1(tr2))
            continue;
          vr1.SetXYZM(tr2.px(), tr2.py(), tr2.pz(), m1);
          auto McPart_2 = tr2.udMcParticle();
          if (std::abs(McPart_2.globalIndex() - McPart_1.globalIndex()) != 1)
            continue;
          if (std::abs(McPart_1.pdgCode()) != 321 || std::abs(McPart_2.pdgCode()) != 321)
            continue;
          if (McPart_1.pdgCode() == McPart_2.pdgCode())
            continue;
          if (!McPart_1.isPhysicalPrimary() || !McPart_2.isPhysicalPrimary())
            continue;
          bool flag = false;
          bool flag1 = false;
          if (McPart_1.has_mothers() && McPart_2.has_mothers()) {
            for (const auto& mother : McPart_1.mothers_as<aod::UDMcParticles>()) {
              if (std::abs(mother.pdgCode()) == 333) {
                v_phi.SetXYZM(mother.px(), mother.py(), mother.pz(), m4);
                registry.get<TH3>(HIST("MC/accMPtRap_phi_T"))->Fill(v_phi.M(), v_phi.Pt(), v_phi.Rapidity(), 1.);
                flag = true;
              }
            }
            for (const auto& mother1 : McPart_1.mothers_as<aod::UDMcParticles>()) {
              if (std::abs(mother1.pdgCode()) == 333) {
                flag1 = true;
              }
            }
          }
          if (flag && flag1) {
            vr0_g.SetXYZM(McPart_1.px(), McPart_1.py(), McPart_1.pz(), m1);
            vr1_g.SetXYZM(McPart_2.px(), McPart_2.py(), McPart_2.pz(), m1);
            vr01_g = vr0_g + vr1_g;
            vr01 = vr0 + vr1;
            registry.get<TH1>(HIST("Reco/selRap"))->Fill(vr01.Rapidity(), 1.);
            registry.get<TH2>(HIST("Reco/selMPt"))->Fill(vr01.M(), vr01.Pt(), 1.);
            registry.get<TH3>(HIST("Reco/selMPtRap"))->Fill(vr01.M(), vr01.Pt(), vr01.Rapidity(), 1.);
            registry.get<TH1>(HIST("Reco/selPt"))->Fill(vr01.Pt(), 1.);
            registry.get<TH1>(HIST("Reco/selM"))->Fill(vr01.M(), 1.);
            registry.get<TH3>(HIST("Reco/selMPtRap_gen"))->Fill(vr01_g.M(), vr01_g.Pt(), vr01_g.Rapidity(), 1.);
          }
        }
      }
    }
    // KStar
    for (auto& [tr1, tr2] : combinations(o2::soa::CombinationsFullIndexPolicy(tracks, tracks))) {
      if (!tr1.has_udMcParticle() || !tr2.has_udMcParticle())
        continue;
      auto m1 = particleMass(pdg, 211);
      auto m2 = particleMass(pdg, 321);
      auto m3 = particleMass(pdg, 313);
      if (!selectionPIDPion1(tr1) || !selectionPIDKaon1(tr2))
        continue;
      //  if (tr1.index() == tr2.index())
      // continue; // We need to run (0,1), (1,0) pairs as well. but same id pairs are not needed.
      auto McPart_1 = tr1.udMcParticle();
      auto McPart_2 = tr2.udMcParticle();
      if (std::abs(McPart_1.pdgCode()) != 211 || std::abs(McPart_2.pdgCode()) != 321)
        continue;
      if (!McPart_1.isPhysicalPrimary() || !McPart_2.isPhysicalPrimary())
        continue;
      if (std::abs(McPart_2.globalIndex() - McPart_1.globalIndex()) != 1)
        continue;

      if (tr1.sign() * tr2.sign() > 0)
        continue;

      vr0.SetXYZM(tr1.px(), tr1.py(), tr1.pz(), m1);
      vr1.SetXYZM(tr2.px(), tr2.py(), tr2.pz(), m2);
      vr0_g.SetXYZM(McPart_1.px(), McPart_1.py(), McPart_1.pz(), m1);
      vr1_g.SetXYZM(McPart_2.px(), McPart_2.py(), McPart_2.pz(), m2);
      vr01_g = vr0_g + vr1_g;
      vr01 = vr0 + vr1;
      if (!trackselector(tr1, parameters) || !trackselector(tr2, parameters)) {
        registry.get<TH1>(HIST("Reco/selM_k"))->Fill(vr01.M(), 1.);
      }
      if (trackselector(tr1, parameters) && trackselector(tr2, parameters)) {
        bool flag = false;
        bool flag1 = false;
        if (McPart_1.has_mothers() && McPart_2.has_mothers()) {
          for (const auto& mother : McPart_1.mothers_as<aod::UDMcParticles>()) {
            if (std::abs(mother.pdgCode()) == 313) {
              v_kstar.SetXYZM(mother.px(), mother.py(), mother.pz(), m3);
              registry.get<TH3>(HIST("MC/accMPtRap_kstar_T"))->Fill(v_kstar.M(), v_kstar.Pt(), v_kstar.Rapidity(), 1.);
              flag = true;
            }
          }
          for (const auto& mother1 : McPart_1.mothers_as<aod::UDMcParticles>()) {
            if (std::abs(mother1.pdgCode()) == 313) {
              flag1 = true;
            }
          }
        }
        if (flag && flag1) {
          registry.get<TH1>(HIST("Reco/selM_k_K"))->Fill(vr01.M(), 1.);
          registry.get<TH1>(HIST("Reco/selRap_k"))->Fill(vr01.Rapidity(), 1.);
          registry.get<TH2>(HIST("Reco/selMPt_k"))->Fill(vr01.M(), vr01.Pt(), 1.);
          registry.get<TH3>(HIST("Reco/selMPtRap_k"))->Fill(vr01.M(), vr01.Pt(), vr01.Rapidity(), 1.);
          registry.get<TH1>(HIST("Reco/selPt_k"))->Fill(vr01.Pt(), 1.);
          // registry.get<TH1>(HIST("Reco/selM_k"))->Fill(vr01_g.M(), 1.);
          registry.get<TH3>(HIST("Reco/selMPtRap_k_gen"))->Fill(vr01_g.M(), vr01_g.Pt(), vr01_g.Rapidity(), 1.);
        }
      }
    }
    registry.get<TH1>(HIST("Reco/nTracks"))->Fill(t1, 1.);
    // now access the McTruth information
    // get McCollision belonging to collision
    if (collision.has_udMcCollision()) {
      // auto mccollision = collision.udMcCollision();
      registry.get<TH1>(HIST("Reco/Stat"))->Fill(3., 1.);
    } else {
      if (verbosity > 0) {
        LOGF(info, "This collision has no associated McCollision");
      }
    }
    // compute the difference between generated and reconstructed momentum
    for (auto track : tracks) {
      // is there an associated McParticle?
      if (track.has_udMcParticle()) {
        auto pTrack = sqrt(track.px() * track.px() + track.py() * track.py() + track.pz() * track.pz());
        auto McPart = track.udMcParticle();
        auto pPart = sqrt(McPart.px() * McPart.px() + McPart.py() * McPart.py() + McPart.pz() * McPart.pz());
        auto pDiff = pTrack - pPart;
        if (std::abs(McPart.pdgCode()) == 321) {
          //	  std::cout<< McPart.pdgCode()<<" ============----------------- "<< McPart.globalIndex()<<std::endl;
        }
        registry.get<TH2>(HIST("Reco/pDiff"))->Fill(pDiff, track.isPVContributor(), 1.);
        if (verbosity > 0) {
          LOGF(info, "  PID: %d Generated: %d Process: %d PV contributor: %d dP: %f", McPart.pdgCode(), McPart.producedByGenerator(), McPart.getProcess(), track.isPVContributor(), pDiff);
        }
      } else {
        registry.get<TH2>(HIST("Reco/pDiff"))->Fill(-5.9, -1, 1.);
        if (verbosity > 0) {
          LOGF(info, "  This track has no associated McParticle");
        }
      }
    }
    if (verbosity > 0) {
      LOGF(info, "");
    }
  }
  PROCESS_SWITCH(SGResonanceAnalyzer, processReco, "Process reconstructed data", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<SGResonanceAnalyzer>(cfgc)};
}
