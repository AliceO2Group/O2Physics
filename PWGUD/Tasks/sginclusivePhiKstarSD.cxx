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
//
/// \file sginclusivePhiKstarSD.cxx
/// \brief Single Gap Event Analyzer for phi and Kstar
/// \author Sandeep Dudi, sandeep.dudi3@gmail.com
/// \since  May 2024

#include <cstdlib>
#include <vector>
#include <TString.h>
#include <TMath.h>
#include "Math/Vector4D.h"
#include "Math/Vector3D.h"
#include "Math/GenVector/Boost.h"
#include "TPDGCode.h"

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/O2DatabasePDGPlugin.h"
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
using namespace o2::constants::physics;

struct SginclusivePhiKstarSD {
  SGSelector sgSelector;
  Service<o2::framework::O2DatabasePDG> pdg;

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  Configurable<float> fv0Cut{"fv0Cut", 50., "FV0A threshold"};
  Configurable<float> ft0aCut{"ft0aCut", 100., "FT0A threshold"};
  Configurable<float> ft0cCut{"ft0cCut", 50., "FT0C threshold"};
  Configurable<float> fddaCut{"fddaCut", 10000., "FDDA threshold"};
  Configurable<float> fddcCut{"fddcCut", 10000., "FDDC threshold"};
  Configurable<float> zdcCut{"zdcCut", 0., "ZDC threshold"};
  Configurable<float> vzCut{"vzCut", 10., "Vz position"};
  Configurable<float> occCut{"occCut", 1000., "Occupancy cut"};
  Configurable<float> hadronicRate{"hadronicRate", 1000., "hadronicRate cut"};
  Configurable<int> useTrs{"useTrs", -1, "kNoCollInTimeRangeStandard cut"};
  Configurable<int> useTrofs{"useTrofs", -1, "kNoCollInRofStandard cut"};
  Configurable<int> useHmpr{"useHmpr", -1, "kNoHighMultCollInPrevRof cut"};
  Configurable<int> useTfb{"useTfb", -1, "kNoTimeFrameBorder cut"};
  Configurable<int> useItsrofb{"useItsrofb", -1, "kNoITSROFrameBorder cut"};
  Configurable<int> useSbp{"useSbp", -1, "kNoSameBunchPileup cut"};
  Configurable<int> useZvtxftovpv{"useZvtxftovpv", -1, "kIsGoodZvtxFT0vsPV cut"};
  Configurable<int> useVtxItsTpc{"useVtxItsTpc", -1, "kIsVertexITSTPC cut"};

  // Track Selections
  Configurable<float> pvCut{"pvCut", 1.0, "Use Only PV tracks"};
  Configurable<float> dcazCut{"dcazCut", 2.0, "dcaZ cut"};
  Configurable<float> dcaxyCut{"dcaxyCut", 0.0, "dcaXY cut (0 for Pt-function)"};
  Configurable<float> tpcChi2Cut{"tpcChi2Cut", 4, "Max tpcChi2NCl"};
  Configurable<float> tpcNClsFindableCut{"tpcNClsFindableCut", 70, "Min tpcNClsFindable"};
  Configurable<float> itsChi2Cut{"itsChi2Cut", 36, "Max itsChi2NCl"};
  Configurable<float> etaCut{"etaCut", 0.9, "Track Pseudorapidity"};
  Configurable<float> ptCut{"ptCut", 0.15, "Track pt cut"};
  Configurable<float> pt1{"pt1", 0.3, "pid selection pt1"};
  Configurable<float> pt2{"pt2", 0.4, "pid selection pt2"};
  Configurable<float> pt3{"pt3", 0.5, "pid selection pt3"};

  Configurable<float> etaGapMin{"etaGapMin", 0.0, "Track eta min"};
  Configurable<float> etaGapMax{"etaGapMax", 0.9, "Track eta max"};
  Configurable<float> etaDG{"etaDG", 0.5, "Track eta DG"};

  Configurable<float> nsigmaTpcCut1{"nsigmaTpcCut1", 3.0, "nsigma tpc cut1"};
  Configurable<float> nsigmaTpcCut2{"nsigmaTpcCut2", 3.0, "nsigma tpc cut2"};
  Configurable<float> nsigmaTpcCut3{"nsigmaTpcCut3", 3.0, "nsigma tpc cut3"};
  Configurable<float> nsigmaTpcCut4{"nsigmaTpcCut4", 3.0, "nsigma tpc cut4"};
  Configurable<float> nsigmaTpcCut{"nsigmaTpcCut", 3.0, "nsigma tpc cut"};
  Configurable<float> nsigmaTofCut{"nsigmaTofCut", 9.0, "nsigma tpc+tof cut"};
  Configurable<float> nsigmaTofCut1{"nsigmaTofCut1", 3.0, "nsigma tof cut"};
  Configurable<float> pionNsigmaCut{"pionNsigmaCut", 3.0, "nsigma tpc cut for kaon"};

  Configurable<int> mintrack{"mintrack", 1, "min track"};
  Configurable<int> maxtrack{"maxtrack", 50, "max track"};
  Configurable<bool> useTof{"useTof", true, "TOF PID"};
  Configurable<bool> ccut{"ccut", true, "TPC + TOF PID"};
  Configurable<bool> kaoncut{"kaoncut", true, " kaon slection cut for kstar "};

  Configurable<bool> qa{"qa", true, ""};
  Configurable<bool> rapidityGap{"rapidityGap", true, ""};

  Configurable<bool> phi{"phi", true, ""};
  Configurable<bool> rho{"rho", true, ""};
  Configurable<bool> kstar{"kstar", true, ""};
  Configurable<bool> fourpion{"fourpion", true, ""};
  Configurable<bool> mc{"mc", true, ""};
  Configurable<int> gapsideMC{"gapsideMC", 1, "gapside MC"};

  Configurable<int> cfgNoMixedEvents{"cfgNoMixedEvents", 1, "Number of mixed events per event"};
  Configurable<int> nBkgRotations{"nBkgRotations", 9, "Number of rotated copies (background) per each original candidate"};
  Configurable<bool> fillRotation{"fillRotation", true, "fill rotation"};
  Configurable<float> confMinRot{"confMinRot", 5.0 * o2::constants::math::PI / 6.0, "Minimum of rotation"};
  Configurable<float> confMaxRot{"confMaxRot", 7.0 * o2::constants::math::PI / 6.0, "Maximum of rotation"};
  //
  Configurable<bool> reconstruction{"reconstruction", true, ""};
  Configurable<int> generatedId{"generatedId", 31, ""};

  void init(InitContext const& context)
  {
    registry.add("GapSide", "Gap Side; Entries", kTH1F, {{4, -1.5, 2.5}});
    registry.add("TrueGapSide", "Gap Side; Entries", kTH1F, {{4, -1.5, 2.5}});
    registry.add("nPVContributors_data", "Multiplicity_dist_before track cut gap A", kTH1F, {{110, 0, 110}});
    registry.add("nPVContributors_data_1", "Multiplicity_dist_before track cut gap C", kTH1F, {{110, 0, 110}});

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

      registry.add("hRotation", "hRotation", kTH1F, {{360, 0.0, o2::constants::math::TwoPI}});
    }
    // qa plots
    if (qa) {
      registry.add("tpc_dedx", "p vs dE/dx", kTH2F, {{100, 0.0, 10.0}, {10000, 0.0, 2500.0}});
      registry.add("tof_beta", "p vs beta", kTH2F, {{100, 0.0, 10.0}, {100, 0.0, 5.0}});

      registry.add("tpc_dedx_kaon", "p#k dE/dx", kTH2F, {{100, 0.0, 10.0}, {10000, 0.0, 2500.0}});
      registry.add("tpc_dedx_pion", "p#pi dE/dx", kTH2F, {{100, 0.0, 10.0}, {10000, 0.0, 2500.0}});
      registry.add("tpc_dedx_kaon_1", "tpc+tof pid cut p#k dE/dx", kTH2F, {{100, 0.0, 10.0}, {10000, 0.0, 2500.0}});
      registry.add("tpc_dedx_kaon_2", "tpc+tof pid cut1 p#k dE/dx", kTH2F, {{100, 0.0, 10.0}, {10000, 0.0, 2500.0}});
      registry.add("tpc_dedx_pion_1", "tpc+tof pid cut p#pi dE/dx", kTH2F, {{100, 0.0, 10.0}, {10000, 0.0, 25000.0}});
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

      if (rapidityGap) {
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
    if (rapidityGap && phi) {
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

    if (rapidityGap && kstar) {
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
    if (mc) {
      // add histograms for the different process functions
      if (context.mOptions.get<bool>("processMCTruth")) {
        registry.add("MC/Stat", "Count generated events; ; Entries", {HistType::kTH1F, {{1, -0.5, 0.5}}});
        registry.add("MC/recCols", "Number of reconstructed collisions; Number of reconstructed collisions; Entries", {HistType::kTH1F, {{31, -0.5, 30.5}}});
        registry.add("MC/nParts", "Number of McParticles per collision; Number of McParticles; Entries", {HistType::kTH1F, {{1001, -0.5, 1000.5}}});
        registry.add("MC/nRecTracks", "Number of reconstructed tracks per McParticle; Number of reconstructed tracks per McParticle; Entries", {HistType::kTH1F, {{11, -0.5, 10.5}}});
        registry.add("MC/genEtaPt", "Generated events; eta (1); Pt (GeV/c)", {HistType::kTH2F, {{60, -1.5, 1.5}, {250, 0.0, 5.0}}});
        registry.add("MC/genRap", "Generated events; Rapidity (1)", {HistType::kTH1F, {{60, -1.5, 1.5}}});
        registry.add("MC/genMPt", "Generated events; Mass (GeV/c^2); Pt (GeV/c)", {HistType::kTH2F, {{220, 0.98, 1.2}, {200, 0.0, 10.0}}});
        registry.add("MC/genM", "Generated events; Mass (GeV/c^2)", {HistType::kTH1F, {{220, 0.98, 1.2}}});
        registry.add("MC/genM_1", "Generated events; Mass (GeV/c^2)", {HistType::kTH1F, {{220, 0.98, 1.2}}});

        registry.add("MC/accMPtRap_phi_G", "Generated Phi; Mass (GeV/c^2); Pt (GeV/c)", {HistType::kTH3F, {{220, 0.98, 1.20}, {200, 0.0, 10.0}, {60, -1.5, 1.5}}});

        registry.add("MC/accEtaPt", "Generated events in acceptance; eta (1); Pt (GeV/c)", {HistType::kTH2F, {{60, -1.5, 1.5}, {250, 0.0, 5.0}}});
        registry.add("MC/accRap", "Generated events in acceptance; Rapidity (1)", {HistType::kTH1F, {{60, -1.5, 1.5}}});
        registry.add("MC/accMPt", "Generated events in acceptance; Mass (GeV/c^2); Pt (GeV/c)", {HistType::kTH2F, {{220, 0.98, 1.20}, {200, 0.0, 10.0}}});
        registry.add("MC/accMPtRap", "Generated events in acceptance; Mass (GeV/c^2); Pt (GeV/c)", {HistType::kTH3F, {{220, 0.98, 1.20}, {200, 0.0, 10.0}, {60, -1.5, 1.5}}});
        registry.add("MC/accM", "Generated events in acceptance; Mass (GeV/c^2)", {HistType::kTH1F, {{220, 0.98, 1.20}}});
        registry.add("MC/selRap", "Selected events in acceptance; Rapidity (1)", {HistType::kTH1F, {{60, -1.5, 1.5}}});
        registry.add("MC/selMPt", "Selected events in acceptance; Mass (GeV/c^2); Pt (GeV/c)", {HistType::kTH2F, {{250, 2.5, 5.0}, {100, 0.0, 1.0}}});
        registry.add("MC/pDiff", "McTruth - reconstructed track momentum; McTruth - reconstructed track momentum; Entries", {HistType::kTH2F, {{240, -6., 6.}, {3, -1.5, 1.5}}});
        // K*0
        registry.add("MC/accMPtRap_kstar_G", "Generated K*0; Mass (GeV/c^2); Pt (GeV/c)", {HistType::kTH3F, {{400, 0., 2.0}, {200, 0.0, 10.0}, {60, -1.5, 1.5}}});
        registry.add("MC/genEtaPt_k", "Generated events; eta (1); Pt (GeV/c)", {HistType::kTH2F, {{60, -1.5, 1.5}, {250, 0.0, 5.0}}});
        registry.add("MC/genRap_k", "Generated events; Rapidity (1)", {HistType::kTH1F, {{60, -1.5, 1.5}}});
        registry.add("MC/genMPt_k", "Generated events; Mass (GeV/c^2); Pt (GeV/c)", {HistType::kTH2F, {{400, 0., 2.0}, {200, 0.0, 10.0}}});
        registry.add("MC/genM_k", "Generated events; Mass (GeV/c^2)", {HistType::kTH1F, {{400, 0., 2.0}}});
        registry.add("MC/genM_1_k", "Generated events; Mass (GeV/c^2)", {HistType::kTH1F, {{400, 0., 2.0}}});

        registry.add("MC/accEtaPt_k", "Generated events in acceptance; eta (1); Pt (GeV/c)", {HistType::kTH2F, {{60, -1.5, 1.5}, {250, 0.0, 5.0}}});
        registry.add("MC/accRap_k", "Generated events in acceptance; Rapidity (1)", {HistType::kTH1F, {{60, -1.5, 1.5}}});
        registry.add("MC/accMPt_k", "Generated events in acceptance; Mass (GeV/c^2); Pt (GeV/c)", {HistType::kTH2F, {{400, 0., 2.0}, {200, 0.0, 10.0}}});
        registry.add("MC/accMPtRap_k", "Generated events in acceptance; Mass (GeV/c^2); Pt (GeV/c)", {HistType::kTH3F, {{400, 0., 2.0}, {200, 0.0, 10.0}, {60, -1.5, 1.5}}});
        registry.add("MC/accM_k", "Generated events in acceptance; Mass (GeV/c^2)", {HistType::kTH1F, {{400, 0., 2.0}}});
      }
      if (context.mOptions.get<bool>("processReco")) {
        registry.add("Reco/Stat", "Count reconstruted events; ; Entries", {HistType::kTH1F, {{5, -0.5, 4.5}}});
        registry.add("Reco/nPVContributors", "Number of PV contributors per collision; Number of PV contributors; Entries", {HistType::kTH1F, {{51, -0.5, 50.5}}});
        registry.add("Reco/selRap", "Selected events in acceptance; Rapidity (1)", {HistType::kTH1F, {{60, -1.5, 1.5}}});
        registry.add("Reco/selMPt", "Reconstructed events in acceptance; Mass (GeV/c^2); Pt (GeV/c)", {HistType::kTH2F, {{220, 0.98, 1.20}, {200, 0.0, 10.0}}});
        registry.add("Reco/selMPtRap", "Reconstructed events in acceptance; Mass (GeV/c^2); Pt (GeV/c)", {HistType::kTH3F, {{220, 0.98, 1.20}, {200, 0.0, 10.0}, {60, -1.5, 1.5}}});
        registry.add("Reco/selMPtRap_gen", "Reconstructed(gen) events in acceptance; Mass (GeV/c^2); Pt (GeV/c)", {HistType::kTH3F, {{220, 0.98, 1.20}, {200, 0.0, 10.0}, {60, -1.5, 1.5}}});
        registry.add("MC/accMPtRap_phi_T", "Reconstrcted Phi; Mass (GeV/c^2); Pt (GeV/c)", {HistType::kTH3F, {{220, 0.98, 1.20}, {200, 0.0, 10.0}, {60, -1.5, 1.5}}});
        registry.add("MC/accMPtRap_kstar_T", "Reconstructed K*0; Mass (GeV/c^2); Pt (GeV/c)", {HistType::kTH3F, {{400, 0., 2.0}, {200, 0.0, 10.0}, {60, -1.5, 1.5}}});

        registry.add("Reco/selPt", "Reconstructed events in acceptance;Pt (GeV/c)", {HistType::kTH1F, {{200, 0.0, 10.0}}});
        registry.add("Reco/selM", "Reconstructed events in acceptance; Mass (GeV/c^2); ", {HistType::kTH1F, {{220, 0.98, 1.20}}});
        registry.add("Reco/mcEtaPt", "Generated events in acceptance; eta (1); Pt (GeV/c)", {HistType::kTH2F, {{60, -1.5, 1.5}, {250, 0.0, 5.0}}});
        registry.add("Reco/mcRap", "Generated events in acceptance; Rapidity (1)", {HistType::kTH1F, {{60, -1.5, 1.5}}});
        registry.add("Reco/mcMPt", "Generated events in acceptance; Mass (GeV/c^2); Pt (GeV/c)", {HistType::kTH2F, {{250, 2.5, 5.0}, {100, 0.0, 1.0}}});
        registry.add("Reco/pDiff", "McTruth - reconstructed track momentum; McTruth - reconstructed track momentum; Entries", {HistType::kTH2F, {{240, -6., 6.}, {3, -1.5, 1.5}}});

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

        // qa
        registry.add("tpc_dedx_mc", "p vs dE/dx", kTH2F, {{100, 0.0, 10.0}, {5000, 0.0, 5000.0}});
        registry.add("tof_beta_mc", "p vs beta", kTH2F, {{100, 0.0, 10.0}, {5000, 0.0, 5000.0}});

        registry.add("tpc_dedx_kaon_mc", "p#k dE/dx", kTH2F, {{100, 0.0, 10.0}, {5000, 0.0, 5000.0}});
        registry.add("tpc_dedx_pion_mc", "p#pi dE/dx", kTH2F, {{100, 0.0, 10.0}, {5000, 0.0, 5000.0}});
        registry.add("tpc_nsigma_kaon_mc", "p#k n#sigma", kTH2F, {{100, 0.0, 10.0}, {100, -10.0, 10.0}});
        registry.add("tpc_nsigma_pion_mc", "p#pi n#sigma", kTH2F, {{100, 0.0, 10.0}, {100, -10.0, 10.0}});

        registry.add("tpc_tof_nsigma_kaon_mc", "p#k n#sigma TPC vs TOF", kTH2F, {{100, -10.0, 10.0}, {100, -10.0, 10.0}});
        registry.add("tpc_tof_nsigma_pion_mc", "p#pi n#sigma TPC vs TOF", kTH2F, {{100, -10.0, 10.0}, {100, -10.0, 10.0}});

        registry.add("tof_nsigma_kaon_mc", "p#k n#sigma", kTH2F, {{100, 0.0, 10.0}, {100, -10.0, 10.0}});
        registry.add("tof_nsigma_pion_mc", "p#pi n#sigma", kTH2F, {{100, 0.0, 10.0}, {100, -10.0, 10.0}});
        registry.add("tof_nsigma_kaon_f_mc", "p#k n#sigma", kTH2F, {{100, 0.0, 10.0}, {100, -10.0, 10.0}});
        registry.add("tof_nsigma_pion_f_mc", "p#pi n#sigma", kTH2F, {{100, 0.0, 10.0}, {100, -10.0, 10.0}});

        registry.add("tpc_nsigma_kaon_f_mc", "p#k n#sigma", kTH2F, {{100, 0.0, 10.0}, {100, -10.0, 10.0}});
        registry.add("tpc_nsigma_pion_f_mc", "p#pi n#sigma", kTH2F, {{100, 0.0, 10.0}, {100, -10.0, 10.0}});

        registry.add("FT0A_0_mc", "T0A amplitude", kTH1F, {{500, 0.0, 500.0}});
        registry.add("FT0A_1_mc", "T0A amplitude", kTH1F, {{20000, 0.0, 20000.0}});
        registry.add("FT0C_0_mc", "T0C amplitude", kTH1F, {{20000, 0.0, 20000.0}});
        registry.add("FT0C_1_mc", "T0C amplitude", kTH1F, {{500, 0.0, 500.0}});
        registry.add("ZDC_A_0_mc", "ZDC amplitude", kTH1F, {{2000, 0.0, 1000.0}});
        registry.add("ZDC_A_1_mc", "ZDC amplitude", kTH1F, {{2000, 0.0, 1000.0}});
        registry.add("ZDC_C_0_mc", "ZDC amplitude", kTH1F, {{2000, 0.0, 1000.0}});
        registry.add("ZDC_C_1_mc", "ZDC amplitude", kTH1F, {{2000, 0.0, 1000.0}});
        registry.add("V0A_0_mc", "V0A amplitude", kTH1F, {{1000, 0.0, 1000.0}});
        registry.add("V0A_1_mc", "V0A amplitude", kTH1F, {{1000, 0.0, 1000.0}});
      }
    }
  }

  //_____________________________________________________________________________
  double cosThetaCollinsSoperFrame(ROOT::Math::PxPyPzMVector pair1,
                                   ROOT::Math::PxPyPzMVector pair2,
                                   ROOT::Math::PxPyPzMVector fourpion)
  {
    double halfSqrtSnn = 2680.;
    double massOfLead208 = 193.6823;
    double momentumBeam = std::sqrt(halfSqrtSnn * halfSqrtSnn * 208 * 208 - massOfLead208 * massOfLead208);

    ROOT::Math::PxPyPzEVector pProjCM(0., 0., -momentumBeam, halfSqrtSnn * 208); // projectile
    ROOT::Math::PxPyPzEVector pTargCM(0., 0., momentumBeam, halfSqrtSnn * 208);  // target

    ROOT::Math::PxPyPzMVector v1 = ROOT::Math::PxPyPzMVector(pair1.Px(), pair1.Py(), pair1.Pz(), pair1.M());
    ROOT::Math::PxPyPzMVector v2 = ROOT::Math::PxPyPzMVector(pair2.Px(), pair2.Py(), pair2.Pz(), pair2.M());
    ROOT::Math::PxPyPzMVector v12 = ROOT::Math::PxPyPzMVector(fourpion.Px(), fourpion.Py(), fourpion.Pz(), fourpion.M());

    // Boost to center of mass frame
    ROOT::Math::Boost boostv12{v12.BoostToCM()};
    ROOT::Math::XYZVectorF v1Cm{(boostv12(v1).Vect()).Unit()};
    ROOT::Math::XYZVectorF v2Cm{(boostv12(v2).Vect()).Unit()};
    ROOT::Math::XYZVectorF beam1Cm{(boostv12(pProjCM).Vect()).Unit()};
    ROOT::Math::XYZVectorF beam2Cm{(boostv12(pTargCM).Vect()).Unit()};
    // Axes
    ROOT::Math::XYZVectorF zaxisCs{((beam1Cm.Unit() - beam2Cm.Unit()).Unit())};

    double cosThetaCs = zaxisCs.Dot((v1Cm));
    return cosThetaCs;
  }

  template <typename T>
  bool selectionPIDKaon1(const T& candidate)
  {
    auto pt = std::sqrt(candidate.px() * candidate.px() + candidate.py() * candidate.py());
    //  float pt1, pt2, pt3 , nsigmatpc_cut1, nsigmatpc_cut2, nsigmatpc_cut3;
    if (useTof && pt < pt1 && std::abs(candidate.tpcNSigmaKa()) < nsigmaTpcCut1) {
      return true;
    }
    if (useTof && pt >= pt1 && pt < pt2 && std::abs(candidate.tpcNSigmaKa()) < nsigmaTpcCut2) {
      return true;
    }
    if (useTof && pt >= pt2 && pt < pt3 && std::abs(candidate.tpcNSigmaKa()) < nsigmaTpcCut3) {
      return true;
    }
    if (useTof && pt >= pt3 && !candidate.hasTOF() && std::abs(candidate.tpcNSigmaKa()) < nsigmaTpcCut) {
      return true;
    }
    if (useTof && candidate.hasTOF() && ccut && (candidate.tofNSigmaKa() * candidate.tofNSigmaKa() + candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa()) < nsigmaTofCut) {
      return true;
    }
    if (useTof && candidate.hasTOF() && !ccut && std::abs(candidate.tpcNSigmaKa()) < nsigmaTpcCut && std::abs(candidate.tofNSigmaKa()) < nsigmaTofCut1) {
      return true;
    }

    if (!useTof && std::abs(candidate.tpcNSigmaKa()) < nsigmaTpcCut) {
      return true;
    }
    return false;
  }

  template <typename T>
  bool selectionPIDPion1(const T& candidate)
  {
    auto pt = std::sqrt(candidate.px() * candidate.px() + candidate.py() * candidate.py());

    if (useTof && pt < pt1 && std::abs(candidate.tpcNSigmaPi()) < nsigmaTpcCut1) {
      return true;
    }
    if (useTof && pt >= pt1 && pt < pt2 && std::abs(candidate.tpcNSigmaPi()) < nsigmaTpcCut2) {
      return true;
    }
    if (useTof && pt >= pt2 && pt < pt3 && std::abs(candidate.tpcNSigmaPi()) < nsigmaTpcCut3) {
      return true;
    }
    if (useTof && pt >= pt3 && !candidate.hasTOF() && std::abs(candidate.tpcNSigmaPi()) < nsigmaTpcCut) {
      return true;
    }
    if (useTof && candidate.hasTOF() && ccut && (candidate.tofNSigmaPi() * candidate.tofNSigmaPi() + candidate.tpcNSigmaPi() * candidate.tpcNSigmaPi()) < nsigmaTofCut) {
      return true;
    }
    if (useTof && candidate.hasTOF() && !ccut && std::abs(candidate.tpcNSigmaPi()) < nsigmaTpcCut && std::abs(candidate.tofNSigmaPi()) < nsigmaTofCut1) {
      return true;
    }
    if (!useTof && std::abs(candidate.tpcNSigmaPi()) < nsigmaTpcCut) {
      return true;
    }
    return false;
  }

  //------------------------------------------------------------------------------------------------------
  double phiCollinsSoperFrame(ROOT::Math::PxPyPzMVector pair1, ROOT::Math::PxPyPzMVector pair2, ROOT::Math::PxPyPzMVector fourpion)
  {
    // Half of the energy per pair of the colliding nucleons.
    double halfSqrtSnn = 2680.;
    double massOfLead208 = 193.6823;
    double momentumBeam = std::sqrt(halfSqrtSnn * halfSqrtSnn * 208 * 208 - massOfLead208 * massOfLead208);

    ROOT::Math::PxPyPzEVector pProjCM(0., 0., -momentumBeam, halfSqrtSnn * 208); // projectile
    ROOT::Math::PxPyPzEVector pTargCM(0., 0., momentumBeam, halfSqrtSnn * 208);  // target

    ROOT::Math::PxPyPzMVector v1 = ROOT::Math::PxPyPzMVector(pair1.Px(), pair1.Py(), pair1.Pz(), pair1.M());
    ROOT::Math::PxPyPzMVector v2 = ROOT::Math::PxPyPzMVector(pair2.Px(), pair2.Py(), pair2.Pz(), pair2.M());
    ROOT::Math::PxPyPzMVector v12 = ROOT::Math::PxPyPzMVector(fourpion.Px(), fourpion.Py(), fourpion.Pz(), fourpion.M());

    // Boost to center of mass frame
    ROOT::Math::Boost boostv12{v12.BoostToCM()};
    ROOT::Math::XYZVectorF v1Cm{(boostv12(v1).Vect()).Unit()};
    ROOT::Math::XYZVectorF v2Cm{(boostv12(v2).Vect()).Unit()};
    ROOT::Math::XYZVectorF beam1Cm{(boostv12(pProjCM).Vect()).Unit()};
    ROOT::Math::XYZVectorF beam2Cm{(boostv12(pTargCM).Vect()).Unit()};
    // Axes
    ROOT::Math::XYZVectorF zaxisCs{((beam1Cm.Unit() - beam2Cm.Unit()).Unit())};
    ROOT::Math::XYZVectorF yaxisCs{(beam1Cm.Cross(beam2Cm)).Unit()};
    ROOT::Math::XYZVectorF xaxisCs{(yaxisCs.Cross(zaxisCs)).Unit()};

    double phi = std::atan2(yaxisCs.Dot(v1Cm), xaxisCs.Dot(v1Cm));
    return phi;
  }

  using UDtracksfull = soa::Join<aod::UDTracks, aod::UDTracksPID, aod::UDTracksExtra, aod::UDTracksFlags, aod::UDTracksDCA>;
  using UDCollisionsFull = soa::Join<aod::UDCollisions, aod::SGCollisions, aod::UDCollisionSelExtras, aod::UDCollisionsSels, aod::UDZdcsReduced>; //
  using UDCollisionFull = UDCollisionsFull::iterator;

  void process(UDCollisionFull const& collision, UDtracksfull const& tracks)
  {
    ROOT::Math::PxPyPzMVector v0;
    ROOT::Math::PxPyPzMVector v1;
    ROOT::Math::PxPyPzMVector v01;
    int gapSide = collision.gapSide();
    float fitCut[5] = {fv0Cut, ft0aCut, ft0cCut, fddaCut, fddcCut};
    std::vector<float> parameters = {pvCut, dcazCut, dcaxyCut, tpcChi2Cut, tpcNClsFindableCut, itsChi2Cut, etaCut, ptCut};
    int truegapSide = sgSelector.trueGap(collision, fitCut[0], fitCut[1], fitCut[2], zdcCut);

    ROOT::Math::PxPyPzMVector phiv;
    ROOT::Math::PxPyPzMVector phiv1;

    std::vector<ROOT::Math::PxPyPzMVector> onlyPionTracksp;
    std::vector<decltype(tracks.begin())> rawPionTracksp;

    std::vector<ROOT::Math::PxPyPzMVector> onlyPionTrackspm;
    std::vector<decltype(tracks.begin())> rawPionTrackspm;

    std::vector<ROOT::Math::PxPyPzMVector> onlyPionTracksn;
    std::vector<decltype(tracks.begin())> rawPionTracksn;

    registry.fill(HIST("GapSide"), gapSide);
    registry.fill(HIST("TrueGapSide"), truegapSide);
    gapSide = truegapSide;
    if (gapSide < 0 || gapSide > 2)
      return;
    if (std::abs(collision.posZ()) > vzCut)
      return;
    if (std::abs(collision.occupancyInTime()) > occCut)
      return;
    if (std::abs(collision.hadronicRate()) > hadronicRate)
      return;

    if (useTrs != -1 && collision.trs() != useTrs)
      return;
    if (useTrofs != -1 && collision.trofs() != useTrofs)
      return;
    if (useHmpr != -1 && collision.hmpr() != useHmpr)
      return;
    if (useTfb != -1 && collision.tfb() != useTfb)
      return;
    if (useItsrofb != -1 && collision.itsROFb() != useItsrofb)
      return;
    if (useSbp != -1 && collision.sbp() != useSbp)
      return;
    if (useZvtxftovpv != -1 && collision.zVtxFT0vPV() != useZvtxftovpv)
      return;
    if (useVtxItsTpc != -1 && collision.vtxITSTPC() != useVtxItsTpc)
      return;

    int mult = collision.numContrib();
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
    int mult0 = 0;
    int mult1 = 0;
    int mult2 = 0;
    int trackgapA = 0;
    int trackgapC = 0;
    int trackDG = 0;
    int trackextra = 0;
    int trackextraDG = 0;

    /*   Partition<UDtracksfull> pvContributors1 = aod::udtrack::isPVContributor == true;
   pvContributors1.bindTable(tracks);
   if (gapSide == 0) {
   registry.get<TH1>(HIST("nPVContributors_data"))->Fill(pvContributors1.size(), 1.);
   }
   if (gapSide == 1) {
   registry.get<TH1>(HIST("nPVContributors_data_1"))->Fill(pvContributors1.size(), 1.);
   }
    */
    for (const auto& track1 : tracks) {
      if (!trackselector(track1, parameters))
        continue;
      v0.SetCoordinates(track1.px(), track1.py(), track1.pz(), o2::constants::physics::MassPionCharged);
      if (selectionPIDPion1(track1)) {
        onlyPionTrackspm.push_back(v0);
        rawPionTrackspm.push_back(track1);
        if (track1.sign() == 1) {
          onlyPionTracksp.push_back(v0);
          rawPionTracksp.push_back(track1);
        }
        if (track1.sign() == -1) {
          onlyPionTracksn.push_back(v0);
          rawPionTracksn.push_back(track1);
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
      if (std::abs(v0.Eta()) < etaDG) {
        trackDG++;
      }
      if (v0.Eta() > etaGapMin && v0.Eta() < etaGapMax) {
        trackgapA++;
      }
      if (v0.Eta() < etaGapMin && v0.Eta() > -etaGapMax) {
        trackgapC++;
      }
      if (std::abs(v0.Eta()) > etaGapMax || std::abs(v0.Eta()) < etaGapMin) {
        trackextra++;
      }
      if (std::abs(v0.Eta()) > etaDG) {
        trackextraDG++;
      }

      if (qa) {
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
    if (qa) {
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
      if (rapidityGap) {
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

    if (rapidityGap) {
      if (trackgapC > 0 && trackgapA == 0 && trackextra == 0) {
        for (const auto& [t0, t1] : combinations(tracks, tracks)) {
          if (!trackselector(t0, parameters) || !trackselector(t1, parameters))
            continue;
          if (phi && selectionPIDKaon1(t0) && selectionPIDKaon1(t1)) {
            // Apply kaon hypothesis and create pairs
            v0.SetCoordinates(t0.px(), t0.py(), t0.pz(), o2::constants::physics::MassKaonCharged);
            v1.SetCoordinates(t1.px(), t1.py(), t1.pz(), o2::constants::physics::MassKaonCharged);
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
        for (const auto& [t0, t1] : combinations(o2::soa::CombinationsFullIndexPolicy(tracks, tracks))) {
          if (!trackselector(t0, parameters) || !trackselector(t1, parameters))
            continue;
          if (t0.globalIndex() == t1.globalIndex())
            continue;
          if (kstar && selectionPIDKaon1(t0) && selectionPIDPion1(t1)) {
            // Apply kaon hypothesis and create pairs
            v0.SetCoordinates(t0.px(), t0.py(), t0.pz(), o2::constants::physics::MassKaonCharged);
            v1.SetCoordinates(t1.px(), t1.py(), t1.pz(), o2::constants::physics::MassPionCharged);
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
        for (const auto& [t0, t1] : combinations(tracks, tracks)) {
          if (!trackselector(t0, parameters) || !trackselector(t1, parameters))
            continue;
          if (phi && selectionPIDKaon1(t0) && selectionPIDKaon1(t1)) {
            // Apply kaon hypothesis and create pairs
            v0.SetCoordinates(t0.px(), t0.py(), t0.pz(), o2::constants::physics::MassKaonCharged);
            v1.SetCoordinates(t1.px(), t1.py(), t1.pz(), o2::constants::physics::MassKaonCharged);
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
        for (const auto& [t0, t1] : combinations(o2::soa::CombinationsFullIndexPolicy(tracks, tracks))) {
          if (!trackselector(t0, parameters) || !trackselector(t1, parameters))
            continue;
          if (t0.globalIndex() == t1.globalIndex())
            continue;
          if (kstar && selectionPIDKaon1(t0) && selectionPIDPion1(t1)) {
            // Apply kaon hypothesis and create pairs
            v0.SetCoordinates(t0.px(), t0.py(), t0.pz(), o2::constants::physics::MassKaonCharged);
            v1.SetCoordinates(t1.px(), t1.py(), t1.pz(), o2::constants::physics::MassPionCharged);
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
        for (const auto& [t0, t1] : combinations(tracks, tracks)) {
          if (!trackselector(t0, parameters) || !trackselector(t1, parameters))
            continue;
          if (phi && selectionPIDKaon1(t0) && selectionPIDKaon1(t1)) {
            // Apply kaon hypothesis and create pairs
            v0.SetCoordinates(t0.px(), t0.py(), t0.pz(), o2::constants::physics::MassKaonCharged);
            v1.SetCoordinates(t1.px(), t1.py(), t1.pz(), o2::constants::physics::MassKaonCharged);
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
        for (const auto& [t0, t1] : combinations(o2::soa::CombinationsFullIndexPolicy(tracks, tracks))) {
          if (!trackselector(t0, parameters) || !trackselector(t1, parameters))
            continue;
          if (t0.globalIndex() == t1.globalIndex())
            continue;
          if (kstar && selectionPIDKaon1(t0) && selectionPIDPion1(t1)) {
            // Apply kaon hypothesis and create pairs
            v0.SetCoordinates(t0.px(), t0.py(), t0.pz(), o2::constants::physics::MassKaonCharged);
            v1.SetCoordinates(t1.px(), t1.py(), t1.pz(), o2::constants::physics::MassPionCharged);
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

    for (const auto& [t0, t1] : combinations(tracks, tracks)) {
      if (!trackselector(t0, parameters) || !trackselector(t1, parameters))
        continue;

      if (phi && selectionPIDKaon1(t0) && selectionPIDKaon1(t1)) {
        // Apply kaon hypothesis and create pairs
        v0.SetCoordinates(t0.px(), t0.py(), t0.pz(), o2::constants::physics::MassKaonCharged);
        v1.SetCoordinates(t1.px(), t1.py(), t1.pz(), o2::constants::physics::MassKaonCharged);
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

            v0.SetCoordinates(rotkaonPx, rotkaonPy, t0.pz(), o2::constants::physics::MassKaonCharged);
            v1.SetCoordinates(t1.px(), t1.py(), t1.pz(), o2::constants::physics::MassKaonCharged);
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
    for (const auto& [t0, t1] : combinations(o2::soa::CombinationsFullIndexPolicy(tracks, tracks))) {
      if (!trackselector(t0, parameters) || !trackselector(t1, parameters))
        continue;
      if (t0.globalIndex() == t1.globalIndex())
        continue;
      if (rho && selectionPIDProton(t0, useTof, nsigmaTpcCut, nsigmaTofCut) && selectionPIDPion1(t1)) {
        v0.SetCoordinates(t0.px(), t0.py(), t0.pz(), o2::constants::physics::MassProton);
        v1.SetCoordinates(t1.px(), t1.py(), t1.pz(), o2::constants::physics::MassPionCharged);
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
        if (kaoncut && t0.tpcNSigmaPi() < pionNsigmaCut)
          continue;
        v0.SetCoordinates(t0.px(), t0.py(), t0.pz(), o2::constants::physics::MassKaonCharged);
        v1.SetCoordinates(t1.px(), t1.py(), t1.pz(), o2::constants::physics::MassPionCharged);
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

            v0.SetCoordinates(rotkaonPx, rotkaonPy, t0.pz(), o2::constants::physics::MassKaonCharged);
            v1.SetCoordinates(t1.px(), t1.py(), t1.pz(), o2::constants::physics::MassPionCharged);
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
        ROOT::Math::PxPyPzMVector pair1, pair2, pair3, pair4;
        if (onlyPionTracksp.size() == 2 && onlyPionTracksn.size() == 2) {
          ROOT::Math::PxPyPzMVector k1 = onlyPionTracksp.at(0);
          ROOT::Math::PxPyPzMVector k2 = onlyPionTracksp.at(1);
          ROOT::Math::PxPyPzMVector k3 = onlyPionTracksn.at(0);
          ROOT::Math::PxPyPzMVector k4 = onlyPionTracksn.at(1);
          phiv = k1 + k2 + k3 + k4;
          pair1 = k1 + k3;
          pair2 = k2 + k4;
          pair3 = k1 + k4;
          pair4 = k2 + k3;
          registry.fill(HIST("os_pppp_pT_2"), phiv.M(), phiv.Pt(), phiv.Rapidity());
          registry.fill(HIST("os_pp_vs_pp_mass"), pair1.M(), pair2.M());
          registry.fill(HIST("os_pp_vs_pp_pt"), pair1.Pt(), pair2.Pt());
          auto costhetaPair = cosThetaCollinsSoperFrame(pair1, pair2, phiv);
          auto phiPair = 1. * o2::constants::math::PI + phiCollinsSoperFrame(pair1, pair2, phiv);
          registry.fill(HIST("phi_dis"), phiPair);
          registry.fill(HIST("costheta_dis"), costhetaPair);
          registry.fill(HIST("costheta_vs_phi"), costhetaPair, phiPair);
          registry.fill(HIST("os_pp_vs_pp_mass1"), pair3.M(), pair4.M());
          registry.fill(HIST("os_pp_vs_pp_pt1"), pair3.Pt(), pair4.Pt());
          auto costhetaPair1 = cosThetaCollinsSoperFrame(pair3, pair4, phiv);
          auto phiPair1 = 1. * o2::constants::math::PI + phiCollinsSoperFrame(pair3, pair4, phiv);
          registry.fill(HIST("phi_dis1"), phiPair1);
          registry.fill(HIST("costheta_dis1"), costhetaPair1);
          registry.fill(HIST("costheta_vs_phi1"), costhetaPair1, phiPair1);
        }
        if (onlyPionTracksp.size() != 2 && onlyPionTracksn.size() != 2) {
          if (onlyPionTracksp.size() + onlyPionTracksn.size() != 4)
            return;
          ROOT::Math::PxPyPzMVector l1 = onlyPionTrackspm.at(0);
          ROOT::Math::PxPyPzMVector l2 = onlyPionTrackspm.at(1);
          ROOT::Math::PxPyPzMVector l3 = onlyPionTrackspm.at(2);
          ROOT::Math::PxPyPzMVector l4 = onlyPionTrackspm.at(3);
          phiv1 = l1 + l2 + l3 + l4;
          registry.fill(HIST("os_pppp_pT_2_ls"), phiv1.M(), phiv1.Pt(), phiv1.Rapidity());
        }
      }
    }
  }
  PROCESS_SWITCH(SginclusivePhiKstarSD, process, "Process unlike event", false);

  using UDCollisionsFull1 = soa::Join<aod::UDCollisions, aod::SGCollisions, aod::UDCollisionSelExtras, aod::UDCollisionsSels, aod::UDZdcsReduced>; //
  SliceCache cache;
  Partition<UDtracksfull> posTracks = aod::udtrack::sign > 0;
  Partition<UDtracksfull> negTracks = aod::udtrack::sign < 0;

  ConfigurableAxis axisVertex{"axisVertex", {10, -10, 10}, "vertex axis for bin"};
  ConfigurableAxis axisMultiplicityClass{"axisMultiplicityClass", {10, 0, 100}, "multiplicity percentile for bin"};
  using BinningTypeVertexContributor = ColumnBinningPolicy<aod::collision::PosZ, aod::collision::NumContrib>;
  void mixprocess(UDCollisionsFull1 const& collisions, UDtracksfull const& /*track*/)
  {
    ROOT::Math::PxPyPzMVector v0;
    ROOT::Math::PxPyPzMVector v1;
    ROOT::Math::PxPyPzMVector v01;
    float fitCut[5] = {fv0Cut, ft0aCut, ft0cCut, fddaCut, fddcCut};
    std::vector<float> parameters = {pvCut, dcazCut, dcaxyCut, tpcChi2Cut, tpcNClsFindableCut, itsChi2Cut, etaCut, ptCut};
    BinningTypeVertexContributor binningOnPositions{{axisVertex, axisMultiplicityClass}, true};
    for (auto const& [collision1, collision2] : o2::soa::selfCombinations(binningOnPositions, cfgNoMixedEvents, -1, collisions, collisions)) {
      int truegapSide1 = sgSelector.trueGap(collision1, fitCut[0], fitCut[1], fitCut[2], zdcCut);
      int truegapSide2 = sgSelector.trueGap(collision2, fitCut[0], fitCut[1], fitCut[2], zdcCut);
      if (truegapSide1 != truegapSide2)
        continue;
      if (truegapSide1 == -1)
        continue;
      if (std::abs(collision1.posZ()) > vzCut || std::abs(collision2.posZ()) > vzCut)
        continue;
      if (std::abs(collision1.occupancyInTime()) > occCut || std::abs(collision2.occupancyInTime()) > occCut)
        continue;
      auto posThisColl = posTracks->sliceByCached(aod::udtrack::udCollisionId, collision1.globalIndex(), cache);
      auto negThisColl = negTracks->sliceByCached(aod::udtrack::udCollisionId, collision2.globalIndex(), cache);
      //      for (auto& [track1, track2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(posThisColl, negThisColl))) {
      for (const auto& [track1, track2] : o2::soa::combinations(posThisColl, negThisColl)) {
        if (!trackselector(track1, parameters) || !trackselector(track2, parameters))
          continue;
        if (phi && selectionPIDKaon1(track1) && selectionPIDKaon1(track2)) {
          v0.SetCoordinates(track1.px(), track1.py(), track1.pz(), o2::constants::physics::MassKaonCharged);
          v1.SetCoordinates(track2.px(), track2.py(), track2.pz(), o2::constants::physics::MassKaonCharged);
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
      for (const auto& [track1, track2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(posThisColl, negThisColl))) {
        if (!trackselector(track1, parameters) || !trackselector(track2, parameters))
          continue;
        if (track1.globalIndex() == track2.globalIndex())
          continue;
        if (kstar && selectionPIDKaon1(track1) && selectionPIDPion1(track2)) {
          v0.SetCoordinates(track1.px(), track1.py(), track1.pz(), o2::constants::physics::MassKaonCharged);
          v1.SetCoordinates(track2.px(), track2.py(), track2.pz(), o2::constants::physics::MassPionCharged);
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
  PROCESS_SWITCH(SginclusivePhiKstarSD, mixprocess, "Process Mixed event", false);

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
    ROOT::Math::PxPyPzMVector v0;
    ROOT::Math::PxPyPzMVector v1;
    ROOT::Math::PxPyPzMVector v01;
    ROOT::Math::PxPyPzMVector vkstar;
    ROOT::Math::PxPyPzMVector vphi;
    for (const auto& mccollision : mccollisions) {
      if (mccollision.generatorsID() != generatedId)
        continue;
      registry.get<TH1>(HIST("MC/Stat"))->Fill(0., 1.);
      // get reconstructed collision which belongs to mccollision
      auto colSlice = collisions.sliceBy(colPerMcCollision, mccollision.globalIndex());
      registry.get<TH1>(HIST("MC/recCols"))->Fill(colSlice.size(), 1.);
      if (reconstruction && colSlice.size() < 1)
        continue;
      // get McParticles which belong to mccollision
      auto partSlice = McParts.sliceBy(partPerMcCollision, mccollision.globalIndex());
      registry.get<TH1>(HIST("MC/nParts"))->Fill(partSlice.size(), 1.);
      for (const auto& [tr1, tr2] : combinations(partSlice, partSlice)) {
        if ((tr1.pdgCode() == kKPlus && tr2.pdgCode() == kPiMinus) || (tr1.pdgCode() == kKMinus && tr2.pdgCode() == kPiPlus) || (tr1.pdgCode() == kPiPlus && tr2.pdgCode() == kKMinus) || (tr1.pdgCode() == kPiMinus && tr2.pdgCode() == kKPlus)) {
          if (std::abs(tr1.pdgCode()) == kKPlus) {
            v0.SetCoordinates(tr1.px(), tr1.py(), tr1.pz(), o2::constants::physics::MassKaonCharged);
            v1.SetCoordinates(tr2.px(), tr2.py(), tr2.pz(), o2::constants::physics::MassPionCharged);
          } else {
            v0.SetCoordinates(tr1.px(), tr1.py(), tr1.pz(), o2::constants::physics::MassPionCharged);
            v1.SetCoordinates(tr2.px(), tr2.py(), tr2.pz(), o2::constants::physics::MassKaonCharged);
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
              if (std::abs(mother.pdgCode()) == o2::constants::physics::Pdg::kK0Star892) {
                vkstar.SetCoordinates(mother.px(), mother.py(), mother.pz(), o2::constants::physics::MassK0Star892);
                registry.get<TH3>(HIST("MC/accMPtRap_kstar_G"))->Fill(vkstar.M(), vkstar.Pt(), vkstar.Rapidity(), 1.);
                flag = true;
              }
            }
            for (const auto& mother1 : tr2.mothers_as<aod::UDMcParticles>()) {
              if (std::abs(mother1.pdgCode()) == o2::constants::physics::Pdg::kK0Star892) {
                flag1 = true;
              }
            }
          }
          if (flag && flag1) {
            registry.get<TH2>(HIST("MC/genMPt_k"))->Fill(v01.M(), v01.Pt(), 1.);
          }
          registry.get<TH1>(HIST("MC/genRap_k"))->Fill(v01.Rapidity(), 1.);
          registry.get<TH1>(HIST("MC/genM_k"))->Fill(v01.M(), 1.);
          registry.get<TH1>(HIST("MC/accRap_k"))->Fill(v01.Rapidity(), 1.);
          registry.get<TH2>(HIST("MC/accMPt_k"))->Fill(v01.M(), v01.Pt(), 1.);
          registry.get<TH1>(HIST("MC/accM_k"))->Fill(v01.M(), 1.);
          registry.get<TH3>(HIST("MC/accMPtRap_k"))->Fill(v01.M(), v01.Pt(), v01.Rapidity(), 1.);
        }

        if (std::abs(tr1.pdgCode()) != kKPlus || std::abs(tr2.pdgCode()) != kKPlus)
          continue;
        v0.SetCoordinates(tr1.px(), tr1.py(), tr1.pz(), o2::constants::physics::MassKaonCharged);
        v1.SetCoordinates(tr2.px(), tr2.py(), tr2.pz(), o2::constants::physics::MassKaonCharged);
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
            if (std::abs(mother.pdgCode()) == o2::constants::physics::Pdg::kPhi) {
              vphi.SetCoordinates(mother.px(), mother.py(), mother.pz(), o2::constants::physics::MassPhi);
              registry.get<TH3>(HIST("MC/accMPtRap_phi_G"))->Fill(vphi.M(), vphi.Pt(), vphi.Rapidity(), 1.);
              flag = true;
            }
          }
          for (const auto& mother1 : tr2.mothers_as<aod::UDMcParticles>()) {
            if (std::abs(mother1.pdgCode()) == o2::constants::physics::Pdg::kPhi) {
              flag1 = true;
            }
          }
        }
        if (flag && flag1) {
          registry.get<TH2>(HIST("MC/genMPt"))->Fill(v01.M(), v01.Pt(), 1.);
        }
        // registry.get<TH1>(HIST("MC/genRap"))->Fill(v01.Rapidity(), 1.);
        // registry.get<TH1>(HIST("MC/genM"))->Fill(v01.M(), 1.);
        registry.get<TH1>(HIST("MC/accRap"))->Fill(v01.Rapidity(), 1.);
        registry.get<TH2>(HIST("MC/accMPt"))->Fill(v01.M(), v01.Pt(), 1.);
        registry.get<TH1>(HIST("MC/accM"))->Fill(v01.M(), 1.);
        registry.get<TH3>(HIST("MC/accMPtRap"))->Fill(v01.M(), v01.Pt(), v01.Rapidity(), 1.);
      }
      // compute the difference between generated and reconstructed particle momentum
      for (const auto& McPart : partSlice) {
        auto trackSlice = tracks.sliceBy(trackPerMcParticle, McPart.globalIndex());
        registry.get<TH1>(HIST("MC/nRecTracks"))->Fill(trackSlice.size(), 1.);
        // compute momentum difference between MCTruth and Reconstruction
        if (trackSlice.size() > 0) {
          for (const auto& track : trackSlice) {
            // std::cout<<trackSlice.size()<<std::endl;
            auto pTrack = std::sqrt(track.px() * track.px() + track.py() * track.py() + track.pz() * track.pz());
            auto pPart = std::sqrt(McPart.px() * McPart.px() + McPart.py() * McPart.py() + McPart.pz() * McPart.pz());
            auto pDiff = pTrack - pPart;
            registry.get<TH2>(HIST("MC/pDiff"))->Fill(pDiff, track.isPVContributor(), 1.);
          }
        } else {
          registry.get<TH2>(HIST("MC/pDiff"))->Fill(-5.9, -1, 1.);
        }
      }
    }
  }
  PROCESS_SWITCH(SginclusivePhiKstarSD, processMCTruth, "Process MC truth", true);
  // ...............................................................................................................
  void processReco(CC const& collision, TCs const& tracks, aod::UDMcCollisions const& /*mccollisions*/, aod::UDMcParticles const& /*McParts*/)
  {
    if (!collision.has_udMcCollision())
      return;
    auto mccoll = collision.udMcCollision();
    if (mccoll.generatorsID() != generatedId)
      return;
    float fitCut[5] = {fv0Cut, ft0aCut, ft0cCut, fddaCut, fddcCut};
    std::vector<float> parameters = {pvCut, dcazCut, dcaxyCut, tpcChi2Cut, tpcNClsFindableCut, itsChi2Cut, etaCut, ptCut};
    int truegapSide = sgSelector.trueGap(collision, fitCut[0], fitCut[1], fitCut[2], zdcCut);
    registry.get<TH1>(HIST("Reco/Stat"))->Fill(4.0, 1.);
    //    Partition<TCs> pvContributors = aod::udtrack::isPVContributor == true;
    // pvContributors.bindTable(tracks);
    if (std::abs(collision.posZ()) > vzCut)
      return;
    if (std::abs(collision.occupancyInTime()) > occCut)
      return;
    registry.get<TH1>(HIST("Reco/Stat"))->Fill(truegapSide, 1.);
    if (truegapSide != gapsideMC)
      return;
    // registry.get<TH1>(HIST("Reco/nPVContributors"))->Fill(pvContributors.size(), 1.);
    ROOT::Math::PxPyPzMVector vphi;
    ROOT::Math::PxPyPzMVector vkstar;
    ROOT::Math::PxPyPzMVector v0;
    ROOT::Math::PxPyPzMVector vr0;
    ROOT::Math::PxPyPzMVector vr1;
    ROOT::Math::PxPyPzMVector vr01;
    ROOT::Math::PxPyPzMVector vr0g;
    ROOT::Math::PxPyPzMVector vr1g;
    ROOT::Math::PxPyPzMVector vr01g;
    int t1 = 0;
    if (truegapSide == 0) {
      registry.fill(HIST("V0A_0_mc"), collision.totalFV0AmplitudeA());
      registry.fill(HIST("FT0A_0_mc"), collision.totalFT0AmplitudeA());
      registry.fill(HIST("FT0C_0_mc"), collision.totalFT0AmplitudeC());
      registry.fill(HIST("ZDC_A_0_mc"), collision.energyCommonZNA());
      registry.fill(HIST("ZDC_C_0_mc"), collision.energyCommonZNC());
    }
    if (truegapSide == 1) {
      registry.fill(HIST("V0A_1_mc"), collision.totalFV0AmplitudeA());
      registry.fill(HIST("FT0A_1_mc"), collision.totalFT0AmplitudeA());
      registry.fill(HIST("FT0C_1_mc"), collision.totalFT0AmplitudeC());
      registry.fill(HIST("ZDC_A_1_mc"), collision.energyCommonZNA());
      registry.fill(HIST("ZDC_C_1_mc"), collision.energyCommonZNC());
    }
    for (const auto& tr1 : tracks) {
      if (!tr1.has_udMcParticle())
        continue;
      auto mcPart1 = tr1.udMcParticle();
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
      v0.SetCoordinates(tr1.px(), tr1.py(), tr1.pz(), o2::constants::physics::MassPionCharged);

      registry.fill(HIST("tpc_dedx_mc"), v0.P(), tr1.tpcSignal());
      registry.fill(HIST("tof_beta_mc"), v0.P(), tr1.beta());
      registry.fill(HIST("tof_nsigma_kaon_f_mc"), v0.Pt(), tr1.tofNSigmaKa());
      registry.fill(HIST("tof_nsigma_pion_f_mc"), v0.Pt(), tr1.tofNSigmaPi());
      registry.fill(HIST("tpc_nsigma_kaon_f_mc"), v0.Pt(), tr1.tpcNSigmaKa());
      registry.fill(HIST("tpc_nsigma_pion_f_mc"), v0.Pt(), tr1.tpcNSigmaPi());
      if (selectionPIDKaon1(tr1)) {
        registry.fill(HIST("tpc_nsigma_kaon_mc"), v0.Pt(), tr1.tpcNSigmaKa());
        registry.fill(HIST("tof_nsigma_kaon_mc"), v0.Pt(), tr1.tofNSigmaKa());
        registry.fill(HIST("tpc_tof_nsigma_kaon_mc"), tr1.tpcNSigmaKa(), tr1.tofNSigmaKa());
      }
      if (selectionPIDPion1(tr1)) {
        registry.fill(HIST("tpc_nsigma_pion_mc"), v0.Pt(), tr1.tpcNSigmaPi());
        registry.fill(HIST("tof_nsigma_pion_mc"), v0.Pt(), tr1.tofNSigmaPi());
        registry.fill(HIST("tpc_tof_nsigma_pion_mc"), tr1.tpcNSigmaPi(), tr1.tofNSigmaPi());
      }

      t1++;
      vr0.SetCoordinates(tr1.px(), tr1.py(), tr1.pz(), o2::constants::physics::MassKaonCharged);
      registry.get<TH1>(HIST("Reco/trpt"))->Fill(vr0.Pt(), 1.);
      registry.get<TH1>(HIST("Reco/treta_k"))->Fill(vr0.Eta(), 1.);
      if (!selectionPIDKaon1(tr1))
        continue;
      registry.get<TH1>(HIST("Reco/trpt_k"))->Fill(vr0.Pt(), 1.);
      int t2 = 0;
      for (const auto& tr2 : tracks) {
        if (!tr2.has_udMcParticle())
          continue;
        if (!trackselector(tr2, parameters))
          continue;
        t2++;
        if (t2 > t1) {
          if (!selectionPIDKaon1(tr2))
            continue;
          vr1.SetCoordinates(tr2.px(), tr2.py(), tr2.pz(), o2::constants::physics::MassKaonCharged);
          auto mcPart2 = tr2.udMcParticle();
          if (std::abs(mcPart2.globalIndex() - mcPart1.globalIndex()) != 1)
            continue;
          if (std::abs(mcPart1.pdgCode()) != kKPlus || std::abs(mcPart2.pdgCode()) != kKPlus)
            continue;
          if (mcPart1.pdgCode() == mcPart2.pdgCode())
            continue;
          if (!mcPart1.isPhysicalPrimary() || !mcPart2.isPhysicalPrimary())
            continue;
          bool flag = false;
          bool flag1 = false;
          if (mcPart1.has_mothers() && mcPart2.has_mothers()) {
            for (const auto& mother : mcPart1.mothers_as<aod::UDMcParticles>()) {
              if (std::abs(mother.pdgCode()) == o2::constants::physics::Pdg::kPhi) {
                vphi.SetCoordinates(mother.px(), mother.py(), mother.pz(), o2::constants::physics::MassPhi);
                registry.get<TH3>(HIST("MC/accMPtRap_phi_T"))->Fill(vphi.M(), vphi.Pt(), vphi.Rapidity(), 1.);
                flag = true;
              }
            }
            for (const auto& mother1 : mcPart1.mothers_as<aod::UDMcParticles>()) {
              if (std::abs(mother1.pdgCode()) == o2::constants::physics::Pdg::kPhi) {
                flag1 = true;
              }
            }
          }
          vr0g.SetCoordinates(mcPart1.px(), mcPart1.py(), mcPart1.pz(), o2::constants::physics::MassKaonCharged);
          vr1g.SetCoordinates(mcPart2.px(), mcPart2.py(), mcPart2.pz(), o2::constants::physics::MassKaonCharged);
          vr01g = vr0g + vr1g;
          vr01 = vr0 + vr1;

          if (flag && flag1) {
            registry.get<TH2>(HIST("Reco/selMPt"))->Fill(vr01.M(), vr01.Pt(), 1.);
          }
          registry.get<TH1>(HIST("Reco/selRap"))->Fill(vr01.Rapidity(), 1.);
          registry.get<TH3>(HIST("Reco/selMPtRap"))->Fill(vr01.M(), vr01.Pt(), vr01.Rapidity(), 1.);
          registry.get<TH1>(HIST("Reco/selPt"))->Fill(vr01.Pt(), 1.);
          registry.get<TH1>(HIST("Reco/selM"))->Fill(vr01.M(), 1.);
          registry.get<TH3>(HIST("Reco/selMPtRap_gen"))->Fill(vr01g.M(), vr01g.Pt(), vr01g.Rapidity(), 1.);
        }
      }
    }
    // KStar
    for (const auto& [tr1, tr2] : combinations(o2::soa::CombinationsFullIndexPolicy(tracks, tracks))) {
      if (!tr1.has_udMcParticle() || !tr2.has_udMcParticle())
        continue;
      if (!selectionPIDPion1(tr1) || !selectionPIDKaon1(tr2))
        continue;
      //  if (tr1.index() == tr2.index())
      // continue; // We need to run (0,1), (1,0) pairs as well. but same id pairs are not needed.
      auto mcPart1 = tr1.udMcParticle();
      auto mcPart2 = tr2.udMcParticle();
      if (std::abs(mcPart1.pdgCode()) != kPiPlus || std::abs(mcPart2.pdgCode()) != kKPlus)
        continue;
      if (!mcPart1.isPhysicalPrimary() || !mcPart2.isPhysicalPrimary())
        continue;
      if (std::abs(mcPart2.globalIndex() - mcPart1.globalIndex()) != 1)
        continue;

      if (tr1.sign() * tr2.sign() > 0)
        continue;

      vr0.SetCoordinates(tr1.px(), tr1.py(), tr1.pz(), o2::constants::physics::MassPionCharged);
      vr1.SetCoordinates(tr2.px(), tr2.py(), tr2.pz(), o2::constants::physics::MassKaonCharged);
      vr0g.SetCoordinates(mcPart1.px(), mcPart1.py(), mcPart1.pz(), o2::constants::physics::MassPionCharged);
      vr1g.SetCoordinates(mcPart2.px(), mcPart2.py(), mcPart2.pz(), o2::constants::physics::MassKaonCharged);
      vr01g = vr0g + vr1g;
      vr01 = vr0 + vr1;
      if (!trackselector(tr1, parameters) || !trackselector(tr2, parameters)) {
        registry.get<TH1>(HIST("Reco/selM_k"))->Fill(vr01.M(), 1.);
      }
      if (trackselector(tr1, parameters) && trackselector(tr2, parameters)) {
        bool flag = false;
        bool flag1 = false;
        if (mcPart1.has_mothers() && mcPart2.has_mothers()) {
          for (const auto& mother : mcPart1.mothers_as<aod::UDMcParticles>()) {
            if (std::abs(mother.pdgCode()) == o2::constants::physics::Pdg::kK0Star892) {
              vkstar.SetCoordinates(mother.px(), mother.py(), mother.pz(), o2::constants::physics::MassK0Star892);
              registry.get<TH3>(HIST("MC/accMPtRap_kstar_T"))->Fill(vkstar.M(), vkstar.Pt(), vkstar.Rapidity(), 1.);
              flag = true;
            }
          }
          for (const auto& mother1 : mcPart2.mothers_as<aod::UDMcParticles>()) {
            if (std::abs(mother1.pdgCode()) == o2::constants::physics::Pdg::kK0Star892) {
              flag1 = true;
            }
          }
        }
        if (flag && flag1) {
          registry.get<TH2>(HIST("Reco/selMPt_k"))->Fill(vr01.M(), vr01.Pt(), 1.);
        }
        registry.get<TH1>(HIST("Reco/selM_k_K"))->Fill(vr01.M(), 1.);
        registry.get<TH1>(HIST("Reco/selRap_k"))->Fill(vr01.Rapidity(), 1.);
        registry.get<TH3>(HIST("Reco/selMPtRap_k"))->Fill(vr01.M(), vr01.Pt(), vr01.Rapidity(), 1.);
        registry.get<TH1>(HIST("Reco/selPt_k"))->Fill(vr01.Pt(), 1.);
        registry.get<TH3>(HIST("Reco/selMPtRap_k_gen"))->Fill(vr01g.M(), vr01g.Pt(), vr01g.Rapidity(), 1.);
      }
    }
    registry.get<TH1>(HIST("Reco/nTracks"))->Fill(t1, 1.);
    // now access the McTruth information
    // get McCollision belonging to collision
    if (collision.has_udMcCollision()) {
      // auto mccollision = collision.udMcCollision();
      registry.get<TH1>(HIST("Reco/Stat"))->Fill(3., 1.);
    }
    // compute the difference between generated and reconstructed momentum
    for (const auto& track : tracks) {
      // is there an associated McParticle?
      if (track.has_udMcParticle()) {
        auto pTrack = std::sqrt(track.px() * track.px() + track.py() * track.py() + track.pz() * track.pz());
        auto mcPart = track.udMcParticle();
        auto pPart = std::sqrt(mcPart.px() * mcPart.px() + mcPart.py() * mcPart.py() + mcPart.pz() * mcPart.pz());
        auto pDiff = pTrack - pPart;
        registry.get<TH2>(HIST("Reco/pDiff"))->Fill(pDiff, track.isPVContributor(), 1.);
      } else {
        registry.get<TH2>(HIST("Reco/pDiff"))->Fill(-5.9, -1, 1.);
      }
    }
  }
  PROCESS_SWITCH(SginclusivePhiKstarSD, processReco, "Process reconstructed data", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<SginclusivePhiKstarSD>(cfgc)};
}
