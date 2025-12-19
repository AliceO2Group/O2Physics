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
// \Single Gap Event Analyzer
// \author Sasha Bylinkin, alexander.bylinkin@gmail.com
// \since  April 2023

#include "PWGUD/Core/SGSelector.h"
#include "PWGUD/Core/SGTrackSelector.h"
#include "PWGUD/DataModel/UDTables.h"

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"

#include <iostream>
// #include "PWGUD/Core/RLhelper.h"
#include "TLorentzVector.h"
#include <TString.h>
using namespace std;
using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
namespace excl_fs
{
DECLARE_SOA_COLUMN(GS, gs, int);
DECLARE_SOA_COLUMN(PV, pv, int);
DECLARE_SOA_COLUMN(ZA, za, int);
DECLARE_SOA_COLUMN(ZC, zc, int);
DECLARE_SOA_COLUMN(SIGN, sign, std::vector<int>);
DECLARE_SOA_COLUMN(PX, px, std::vector<float>);
DECLARE_SOA_COLUMN(PY, py, std::vector<float>);
DECLARE_SOA_COLUMN(PZ, pz, std::vector<float>);
DECLARE_SOA_COLUMN(ISELEC, iselec, std::vector<int>);
DECLARE_SOA_COLUMN(ISMUON, ismuon, std::vector<int>);
DECLARE_SOA_COLUMN(ISPION, ispion, std::vector<int>);
DECLARE_SOA_COLUMN(ISKAON, iskaon, std::vector<int>);
DECLARE_SOA_COLUMN(ISPROTON, isproton, std::vector<int>);
} // namespace excl_fs
namespace o2::aod
{
DECLARE_SOA_TABLE(Excl_fs, "AOD", "EXCL_FS",
                  excl_fs::GS, excl_fs::PV, excl_fs::ZA, excl_fs::ZC, excl_fs::SIGN, excl_fs::PX, excl_fs::PY, excl_fs::PZ, excl_fs::ISELEC, excl_fs::ISMUON, excl_fs::ISPION, excl_fs::ISKAON, excl_fs::ISPROTON);
}
struct SGSpectraAnalyzer {
  SGSelector sgSelector;
  Service<o2::framework::O2DatabasePDG> pdg;
  Configurable<float> FV0_cut{"FV0", 100., "FV0A threshold"};
  Configurable<float> FT0A_cut{"FT0A", 100., "FT0A threshold"};
  Configurable<float> FT0C_cut{"FT0C", 50., "FT0C threshold"};
  Configurable<float> FDDA_cut{"FDDA", 10000., "FDDA threshold"};
  Configurable<float> FDDC_cut{"FDDC", 10000., "FDDC threshold"};
  Configurable<float> ZDC_cut{"ZDC", 10., "ZDC threshold"};
  Configurable<bool> use_tof{"Use_TOF", true, "TOF PID"};
  // Track Selections
  Configurable<float> PV_cut{"PV_cut", 1.0, "Use Only PV tracks"};
  Configurable<float> dcaZ_cut{"dcaZ_cut", 2.0, "dcaZ cut"};
  Configurable<float> dcaXY_cut{"dcaXY_cut", 0.0, "dcaXY cut (0 for Pt-function)"};
  Configurable<float> tpcChi2_cut{"tpcChi2_cut", 4, "Max tpcChi2NCl"};
  Configurable<float> tpcNClsFindable_cut{"tpcNClsFindable_cut", 70, "Min tpcNClsFindable"};
  Configurable<float> itsChi2_cut{"itsChi2_cut", 36, "Max itsChi2NCl"};
  Configurable<float> eta_cut{"Eta", 0.9, "Eta cut"};
  Configurable<float> pt_cut{"pt_cut", 0.1, "Track Pt"};
  HistogramRegistry registry{
    "registry",
    {// Pion histograms for each eta bin and gapSide
     {"E_mult_SG", "Photon-Ion c.m.s. Energy vs Multiplicity", {HistType::kTH2F, {{50, -.5, 49.5}, {250, 10, 300}}}},
     {"E_mult_DG", "Photon-Ion c.m.s. Energy vs Multiplicity", {HistType::kTH2F, {{50, -.5, 49.5}, {250, 10, 300}}}},
     {"E_M_SG", "Photon-Ion c.m.s. Energy vs Mass", {HistType::kTH2F, {{1000, .05, 25.05}, {250, 10, 300}}}},
     {"E_M_DG", "Photon-Ion c.m.s. Energy vs Mass", {HistType::kTH2F, {{1000, .05, 25.05}, {250, 10, 300}}}},
     {"E_Y_SG", "Photon-Ion c.m.s. Energy vs Rapidity", {HistType::kTH2F, {{200, -1., 1.}, {250, 10, 300}}}},
     {"E_Y_DG", "Photon-Ion c.m.s. Energy vs Rapidity", {HistType::kTH2F, {{200, -1., 1.}, {250, 10, 300}}}},

     {"ITS_Cluster_nonPV", "ITS Cluster Size", {HistType::kTH1F, {{140, -.5, 139.5}}}},
     {"all_tracks_SG", "All Tracks SG", {HistType::kTH1F, {{50, -.5, 49.5}}}},
     {"all_tracks_DG", "All Tracks DG", {HistType::kTH1F, {{50, -.5, 49.5}}}},
     {"good_tracks_SG", "Good Tracks SG", {HistType::kTH1F, {{50, -.5, 49.5}}}},
     {"good_tracks_DG", "Good Tracks DG", {HistType::kTH1F, {{50, -.5, 49.5}}}},
     {"ITS_Cluster_PV", "ITS Cluster Size", {HistType::kTH1F, {{140, -.5, 139.5}}}},
     {"ITS_Chi2_PV", "ITS Chi2", {HistType::kTH1F, {{10000, -999.5, 999.5}}}},
     {"ITS_Chi2_nonPV", "ITS Chi2", {HistType::kTH1F, {{10000, -999.5, 999.5}}}},
     {"TPC_Chi2_PV", "TPC Chi2", {HistType::kTH1F, {{10000, -999.5, 999.5}}}},
     {"TPC_Chi2_nonPV", "TPC Chi2", {HistType::kTH1F, {{10000, -999.5, 999.5}}}},
     {"Length_nonPV", "Length", {HistType::kTH1F, {{10000, -1199.5, 999.5}}}},
     {"Length_PV", "Length", {HistType::kTH1F, {{10000, -1199.5, 999.5}}}},
     {"DcaZ_PV", "Dca Z", {HistType::kTH1F, {{10000, -19.5, 19.5}}}},
     {"DcaXY_PV", "Dca Z", {HistType::kTH1F, {{10000, -9.5, 9.5}}}},
     {"DcaZ_nonPV", "Dca Z", {HistType::kTH1F, {{10000, -19.5, 19.5}}}},
     {"DcaXY_nonPV", "Dca Z", {HistType::kTH1F, {{10000, -9.5, 9.5}}}},
     {"TPC_Cluster_PV", "TPC Cluster Size", {HistType::kTH1F, {{240, -.5, 239.5}}}},
     {"TPC_Cluster_nonPV", "TPC Cluster Size", {HistType::kTH1F, {{240, -.5, 239.5}}}},
     {"TPC_ClusterS_PV", "TPC Cluster Shared", {HistType::kTH1F, {{240, -.5, 239.5}}}},
     {"TPC_ClusterS_nonPV", "TPC Cluster Shared", {HistType::kTH1F, {{240, -.5, 239.5}}}},
     {"TPC_ClustermF_PV", "TPC Cluster Size minus Found", {HistType::kTH1F, {{240, -.5, 239.5}}}},
     {"TPC_ClustermF_nonPV", "TPC Cluster Size minus Found", {HistType::kTH1F, {{240, -.5, 239.5}}}},
     {"TPC_ClustermCR_PV", "TPC Cluster Size minus Crossed Rows", {HistType::kTH1F, {{240, -60.5, 179.5}}}},
     {"TPC_ClustermCR_nonPV", "TPC Cluster Size minus Crossed Rows", {HistType::kTH1F, {{240, -60.5, 179.5}}}},
     {"TPC_IP_PV", "TPC Inner Param", {HistType::kTH1F, {{200, -.5, 4999.5}}}},
     {"TPC_IP_nonPV", "TPC Inner Param", {HistType::kTH1F, {{200, -.5, 4999.5}}}},
     {"hPtPion_gap0", "Pion pT in eta [-1, 1] Gap 0; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtPion_gap1", "Pion pT in eta [-1, 1] Gap 1; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtPion_gap2", "Pion pT in eta [-1, 1] Gap 2; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hEtaPion_gap0", "Pion eta Gap 0; eta; Entries", {HistType::kTH1F, {{100, -1, 1}}}},
     {"hEtaPion_gap1", "Pion eta Gap 1; eta; Entries", {HistType::kTH1F, {{100, -1, 1}}}},
     {"hEtaPion_gap2", "Pion eta Gap 2; eta; Entries", {HistType::kTH1F, {{100, -1, 1}}}},
     {"hPtPion_etaBin1_gap0", "Pion pT in eta [-1, -0.5] Gap 0; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtPion_etaBin2_gap0", "Pion pT in eta [-0.5, 0] Gap 0; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtPion_etaBin3_gap0", "Pion pT in eta [0, 0.5] Gap 0; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtPion_etaBin4_gap0", "Pion pT in eta [0.5, 1] Gap 0; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtPion_etaBin1_gap1", "Pion pT in eta [-1, -0.5] Gap 1; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtPion_etaBin2_gap1", "Pion pT in eta [-0.5, 0] Gap 1; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtPion_etaBin3_gap1", "Pion pT in eta [0, 0.5] Gap 1; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtPion_etaBin4_gap1", "Pion pT in eta [0.5, 1] Gap 1; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtPion_etaBin1_gap2", "Pion pT in eta [-1, -0.5] Gap 2; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtPion_etaBin2_gap2", "Pion pT in eta [-0.5, 0] Gap 2; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtPion_etaBin3_gap2", "Pion pT in eta [0, 0.5] Gap 2; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtPion_etaBin4_gap2", "Pion pT in eta [0.5, 1] Gap 2; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},

     // Kaon histograms for each eta bin and gapSide
     {"hPtKaon_gap0", "Kaon pT in eta [-1, 1] Gap 0; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtKaon_gap1", "Kaon pT in eta [-1, 1] Gap 1; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtKaon_gap2", "Kaon pT in eta [-1, 1] Gap 2; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hEtaKaon_gap0", "Kaon eta Gap 0; eta; Entries", {HistType::kTH1F, {{100, -1, 1}}}},
     {"hEtaKaon_gap1", "Kaon eta Gap 1; eta; Entries", {HistType::kTH1F, {{100, -1, 1}}}},
     {"hEtaKaon_gap2", "Kaon eta Gap 2; eta; Entries", {HistType::kTH1F, {{100, -1, 1}}}},
     {"hPtKaon_etaBin1_gap0", "Kaon pT in eta [-1, -0.5] Gap 0; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtKaon_etaBin2_gap0", "Kaon pT in eta [-0.5, 0] Gap 0; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtKaon_etaBin3_gap0", "Kaon pT in eta [0, 0.5] Gap 0; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtKaon_etaBin4_gap0", "Kaon pT in eta [0.5, 1] Gap 0; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtKaon_etaBin1_gap1", "Kaon pT in eta [-1, -0.5] Gap 1; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtKaon_etaBin2_gap1", "Kaon pT in eta [-0.5, 0] Gap 1; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtKaon_etaBin3_gap1", "Kaon pT in eta [0, 0.5] Gap 1; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtKaon_etaBin4_gap1", "Kaon pT in eta [0.5, 1] Gap 1; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtKaon_etaBin1_gap2", "Kaon pT in eta [-1, -0.5] Gap 2; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtKaon_etaBin2_gap2", "Kaon pT in eta [-0.5, 0] Gap 2; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtKaon_etaBin3_gap2", "Kaon pT in eta [0, 0.5] Gap 2; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtKaon_etaBin4_gap2", "Kaon pT in eta [0.5, 1] Gap 2; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},

     // Proton histograms for each eta bin and gapSide
     {"hPtProton_gap0", "Proton pT in eta [-1, 1] Gap 0; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtProton_gap1", "Proton pT in eta [-1, 1] Gap 1; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtProton_gap2", "Proton pT in eta [-1, 1] Gap 2; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hEtaProton_gap0", "Proton eta Gap 0; eta; Entries", {HistType::kTH1F, {{100, -1, 1}}}},
     {"hEtaProton_gap1", "Proton eta Gap 1; eta; Entries", {HistType::kTH1F, {{100, -1, 1}}}},
     {"hEtaProton_gap2", "Proton eta Gap 2; eta; Entries", {HistType::kTH1F, {{100, -1, 1}}}},
     {"hPtProton_etaBin1_gap0", "Proton pT in eta [-1, -0.5] Gap 0; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtProton_etaBin2_gap0", "Proton pT in eta [-0.5, 0] Gap 0; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtProton_etaBin3_gap0", "Proton pT in eta [0, 0.5] Gap 0; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtProton_etaBin4_gap0", "Proton pT in eta [0.5, 1] Gap 0; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtProton_etaBin1_gap1", "Proton pT in eta [-1, -0.5] Gap 1; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtProton_etaBin2_gap1", "Proton pT in eta [-0.5, 0] Gap 1; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtProton_etaBin3_gap1", "Proton pT in eta [0, 0.5] Gap 1; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtProton_etaBin4_gap1", "Proton pT in eta [0.5, 1] Gap 1; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtProton_etaBin1_gap2", "Proton pT in eta [-1, -0.5] Gap 2; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtProton_etaBin2_gap2", "Proton pT in eta [-0.5, 0] Gap 2; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtProton_etaBin3_gap2", "Proton pT in eta [0, 0.5] Gap 2; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}},
     {"hPtProton_etaBin4_gap2", "Proton pT in eta [0.5, 1] Gap 2; pT (GeV/c); Entries", {HistType::kTH1F, {{100, 0, 10}}}}}};
  using udtracks = soa::Join<aod::UDTracks, aod::UDTracksExtra, aod::UDTracksPID>;
  // using udtracks = soa::Join<aod::UDTracks,aod::UDTracksExtra,aod::UDTracksPID,aod::UDTracksDCA>;
  using udtracksfull = soa::Join<aod::UDTracks, aod::UDTracksPID, aod::UDTracksExtra, aod::UDTracksFlags, aod::UDTracksDCA>;
  // using UDCollisionsFull = soa::Join<aod::UDCollisions, aod::SGCollisions, aod::UDCollisionsSels>;
  using UDCollisionsFull = soa::Join<aod::UDCollisions, aod::SGCollisions, aod::UDCollisionsSels, aod::UDZdcsReduced>; //
  // using UDCollisionsFull = soa::Join<aod::UDCollisions, aod::SGCollisions, aod::UDCollisionsSels, aod::UDZdcs>;
  using UDCollisionFull = UDCollisionsFull::iterator;

  /*
  template <typename T>
  bool ispion(const T& track, bool use_tof){
    int pid = 0;
    if (std::abs(track.tpcNSigmaPi() < 3)){
      if (use_tof && track.tofChi2()>-1){
        if (std::abs(track.tofNSigmaPi() < 3)) pid = 1;
        else pid = 0;
      } else pid = 1;
    }
          if (pid ==1) return true;
    else return false;
  }
  template <typename T>
  bool iskaon(const T& track, bool use_tof){
    int pid = 0;
    if (std::abs(track.tpcNSigmaKa() < 3)){
      if (use_tof && track.tofChi2()>-1){
        if (std::abs(track.tofNSigmaKa() < 3)) pid = 2;
        else pid =  0;
      } else pid = 2;
    }
          if (pid ==2) return true;
    else return false;
  }
  template <typename T>
  bool isproton(const T& track, bool use_tof){
    int pid = 0;
          if (std::abs(track.tpcNSigmaPr() < 3)){
      if (use_tof && track.tofChi2()>-1){
        if (std::abs(track.tofNSigmaPr() < 3)) pid = 3;
        else pid =  0;
      } else pid = 3;
    }
          if (pid ==3) return true;
    else return false;
  }
  */
  void process(UDCollisionFull const& collision, udtracksfull const& tracks)
  {
    const float mpion = pdg->Mass(211);
    const float mkaon = pdg->Mass(321);
    const float mproton = pdg->Mass(2212);
    TLorentzVector a;
    TLorentzVector am;
    TLorentzVector sum;
    int goodtracks = 0;
    int alltracks = 0;
    sum.SetXYZM(0, 0, 0, 0);
    int gapSide = collision.gapSide();
    float FIT_cut[5] = {FV0_cut, FT0A_cut, FT0C_cut, FDDA_cut, FDDC_cut};
    std::vector<float> parameters = {PV_cut, dcaZ_cut, dcaXY_cut, tpcChi2_cut, tpcNClsFindable_cut, itsChi2_cut, eta_cut, pt_cut};
    int truegapSide = sgSelector.trueGap(collision, FIT_cut[0], FIT_cut[1], FIT_cut[2], ZDC_cut);
    gapSide = truegapSide;
    if (gapSide < 0 || gapSide > 2)
      return;
    for (auto& track : tracks) {
      if (!track.isPVContributor()) {
        registry.get<TH1>(HIST("ITS_Cluster_nonPV"))->Fill(track.itsClusterSizes());
        registry.get<TH1>(HIST("ITS_Chi2_nonPV"))->Fill(track.itsChi2NCl());
        registry.get<TH1>(HIST("TPC_Chi2_nonPV"))->Fill(track.tpcChi2NCl());
        registry.get<TH1>(HIST("Length_nonPV"))->Fill(track.length());
        registry.get<TH1>(HIST("TPC_Cluster_nonPV"))->Fill(track.tpcNClsFindable());
        registry.get<TH1>(HIST("TPC_ClusterS_nonPV"))->Fill(track.tpcNClsShared());
        registry.get<TH1>(HIST("TPC_ClustermF_nonPV"))->Fill(track.tpcNClsFindableMinusFound());
        registry.get<TH1>(HIST("TPC_ClustermCR_nonPV"))->Fill(track.tpcNClsFindableMinusCrossedRows());
        registry.get<TH1>(HIST("TPC_IP_nonPV"))->Fill(track.tpcInnerParam());
        registry.get<TH1>(HIST("DcaZ_nonPV"))->Fill(track.dcaZ());
        registry.get<TH1>(HIST("DcaXY_nonPV"))->Fill(track.dcaXY());
      } else {
        registry.get<TH1>(HIST("ITS_Cluster_PV"))->Fill(track.itsClusterSizes());
        registry.get<TH1>(HIST("ITS_Chi2_PV"))->Fill(track.itsChi2NCl());
        registry.get<TH1>(HIST("TPC_Chi2_PV"))->Fill(track.tpcChi2NCl());
        registry.get<TH1>(HIST("Length_PV"))->Fill(track.length());
        registry.get<TH1>(HIST("TPC_Cluster_PV"))->Fill(track.tpcNClsFindable());
        registry.get<TH1>(HIST("TPC_ClusterS_PV"))->Fill(track.tpcNClsShared());
        registry.get<TH1>(HIST("TPC_ClustermF_PV"))->Fill(track.tpcNClsFindableMinusFound());
        registry.get<TH1>(HIST("TPC_ClustermCR_PV"))->Fill(track.tpcNClsFindableMinusCrossedRows());
        registry.get<TH1>(HIST("TPC_IP_PV"))->Fill(track.tpcInnerParam());
        registry.get<TH1>(HIST("DcaZ_PV"))->Fill(track.dcaZ());
        registry.get<TH1>(HIST("DcaXY_PV"))->Fill(track.dcaXY());
        alltracks++;
        if (trackselector(track, parameters)) {
          int track_pid = trackpid(track, use_tof);
          // if (ispion(track, use_tof)) {
          if (track_pid <= 1) {
            a.SetXYZM(track.px(), track.py(), track.pz(), mpion);
            am.SetXYZM(track.px(), track.py(), -track.pz(), mpion);
            if (std::abs(a.Eta()) < eta_cut) {
              goodtracks++;
              fillHistograms("Pion", a.Pt(), a.Eta(), gapSide);
              if (gapSide == 0)
                sum = sum + a;
              else
                sum = sum + am;
            }
          }
          //      if (iskaon(track, use_tof)) {
          if (track_pid == 2) {
            a.SetXYZM(track.px(), track.py(), track.pz(), mkaon);
            am.SetXYZM(track.px(), track.py(), -track.pz(), mkaon);
            if (std::abs(a.Eta()) < eta_cut) {
              goodtracks++;
              fillHistograms("Kaon", a.Pt(), a.Eta(), gapSide);
              if (gapSide == 0)
                sum = sum + a;
              else
                sum = sum + am;
            }
          }
          //    if (isproton(track, use_tof)) {
          if (track_pid == 3) {
            a.SetXYZM(track.px(), track.py(), track.pz(), mproton);
            am.SetXYZM(track.px(), track.py(), -track.pz(), mproton);
            if (std::abs(a.Eta()) < eta_cut) {
              goodtracks++;
              fillHistograms("Proton", a.Pt(), a.Eta(), gapSide);
              if (gapSide == 0)
                sum = sum + a;
              else
                sum = sum + am;
            }
          }
        }
      }
    }
    if (goodtracks > 1) {
      float W_gPb = TMath::Sqrt(2 * 2680 * sum.M() * TMath::Exp(sum.Rapidity()));
      if (sum.M() < .2)
        std::cout << goodtracks << "\t" << sum.M() << "\t" << sum.Pt() << std::endl;
      if (gapSide < 2) {
        registry.get<TH2>(HIST("E_mult_SG"))->Fill(goodtracks, W_gPb);
        registry.get<TH2>(HIST("E_M_SG"))->Fill(sum.M(), W_gPb);
        registry.get<TH2>(HIST("E_Y_SG"))->Fill(sum.Rapidity(), W_gPb);
        registry.get<TH1>(HIST("all_tracks_SG"))->Fill(alltracks);
        registry.get<TH1>(HIST("good_tracks_SG"))->Fill(goodtracks);
      } else {
        registry.get<TH2>(HIST("E_mult_DG"))->Fill(goodtracks, W_gPb);
        registry.get<TH2>(HIST("E_M_DG"))->Fill(sum.M(), W_gPb);
        registry.get<TH2>(HIST("E_Y_DG"))->Fill(sum.Rapidity(), W_gPb);
        registry.get<TH1>(HIST("all_tracks_DG"))->Fill(alltracks);
        registry.get<TH1>(HIST("good_tracks_DG"))->Fill(goodtracks);
      }
    }
  }

  void fillHistograms(const std::string& particleType, float pt, float eta, int gapSide)
  {
    if (particleType == "Pion") {
      if (gapSide == 0) {
        registry.get<TH1>(HIST("hPtPion_gap0"))->Fill(pt);
        registry.get<TH1>(HIST("hEtaPion_gap0"))->Fill(eta);
      } else if (gapSide == 1) {
        registry.get<TH1>(HIST("hPtPion_gap1"))->Fill(pt);
        registry.get<TH1>(HIST("hEtaPion_gap1"))->Fill(eta);
      } else if (gapSide == 2) {
        registry.get<TH1>(HIST("hPtPion_gap2"))->Fill(pt);
        registry.get<TH1>(HIST("hEtaPion_gap2"))->Fill(eta);
      }
    } else if (particleType == "Kaon") {
      if (gapSide == 0) {
        registry.get<TH1>(HIST("hPtKaon_gap0"))->Fill(pt);
        registry.get<TH1>(HIST("hEtaKaon_gap0"))->Fill(eta);
      } else if (gapSide == 1) {
        registry.get<TH1>(HIST("hPtKaon_gap1"))->Fill(pt);
        registry.get<TH1>(HIST("hEtaKaon_gap1"))->Fill(eta);
      } else if (gapSide == 2) {
        registry.get<TH1>(HIST("hPtKaon_gap2"))->Fill(pt);
        registry.get<TH1>(HIST("hEtaKaon_gap2"))->Fill(eta);
      }
    } else if (particleType == "Proton") {
      if (gapSide == 0) {
        registry.get<TH1>(HIST("hPtProton_gap0"))->Fill(pt);
        registry.get<TH1>(HIST("hEtaProton_gap0"))->Fill(eta);
      } else if (gapSide == 1) {
        registry.get<TH1>(HIST("hPtProton_gap1"))->Fill(pt);
        registry.get<TH1>(HIST("hEtaProton_gap1"))->Fill(eta);
      } else if (gapSide == 2) {
        registry.get<TH1>(HIST("hPtProton_gap2"))->Fill(pt);
        registry.get<TH1>(HIST("hEtaProton_gap2"))->Fill(eta);
      }
    }
    std::vector<std::pair<float, float>> etaBins = {{-1, -0.5}, {-0.5, 0}, {0, 0.5}, {0.5, 1}};
    for (int i = 0; i < 4; ++i) {
      if (eta > etaBins[i].first && eta <= etaBins[i].second) {
        if (i == 0) {
          if (particleType == "Pion") {
            if (gapSide == 0) {
              registry.get<TH1>(HIST("hPtPion_etaBin1_gap0"))->Fill(pt);
            } else if (gapSide == 1) {
              registry.get<TH1>(HIST("hPtPion_etaBin1_gap1"))->Fill(pt);
            } else if (gapSide == 2) {
              registry.get<TH1>(HIST("hPtPion_etaBin1_gap2"))->Fill(pt);
            }
          } else if (particleType == "Kaon") {
            if (gapSide == 0) {
              registry.get<TH1>(HIST("hPtKaon_etaBin1_gap0"))->Fill(pt);
            } else if (gapSide == 1) {
              registry.get<TH1>(HIST("hPtKaon_etaBin1_gap1"))->Fill(pt);
            } else if (gapSide == 2) {
              registry.get<TH1>(HIST("hPtKaon_etaBin1_gap2"))->Fill(pt);
            }
          } else if (particleType == "Proton") {
            if (gapSide == 0) {
              registry.get<TH1>(HIST("hPtProton_etaBin1_gap0"))->Fill(pt);
            } else if (gapSide == 1) {
              registry.get<TH1>(HIST("hPtProton_etaBin1_gap1"))->Fill(pt);
            } else if (gapSide == 2) {
              registry.get<TH1>(HIST("hPtProton_etaBin1_gap2"))->Fill(pt);
            }
          }
        }
        if (i == 1) {
          if (particleType == "Pion") {
            if (gapSide == 0) {
              registry.get<TH1>(HIST("hPtPion_etaBin2_gap0"))->Fill(pt);
            } else if (gapSide == 1) {
              registry.get<TH1>(HIST("hPtPion_etaBin2_gap1"))->Fill(pt);
            } else if (gapSide == 2) {
              registry.get<TH1>(HIST("hPtPion_etaBin2_gap2"))->Fill(pt);
            }
          } else if (particleType == "Kaon") {
            if (gapSide == 0) {
              registry.get<TH1>(HIST("hPtKaon_etaBin2_gap0"))->Fill(pt);
            } else if (gapSide == 1) {
              registry.get<TH1>(HIST("hPtKaon_etaBin2_gap1"))->Fill(pt);
            } else if (gapSide == 2) {
              registry.get<TH1>(HIST("hPtKaon_etaBin2_gap2"))->Fill(pt);
            }
          } else if (particleType == "Proton") {
            if (gapSide == 0) {
              registry.get<TH1>(HIST("hPtProton_etaBin2_gap0"))->Fill(pt);
            } else if (gapSide == 1) {
              registry.get<TH1>(HIST("hPtProton_etaBin2_gap1"))->Fill(pt);
            } else if (gapSide == 2) {
              registry.get<TH1>(HIST("hPtProton_etaBin2_gap2"))->Fill(pt);
            }
          }
        }
        if (i == 2) {
          if (particleType == "Pion") {
            if (gapSide == 0) {
              registry.get<TH1>(HIST("hPtPion_etaBin3_gap0"))->Fill(pt);
            } else if (gapSide == 1) {
              registry.get<TH1>(HIST("hPtPion_etaBin3_gap1"))->Fill(pt);
            } else if (gapSide == 2) {
              registry.get<TH1>(HIST("hPtPion_etaBin3_gap2"))->Fill(pt);
            }
          } else if (particleType == "Kaon") {
            if (gapSide == 0) {
              registry.get<TH1>(HIST("hPtKaon_etaBin3_gap0"))->Fill(pt);
            } else if (gapSide == 1) {
              registry.get<TH1>(HIST("hPtKaon_etaBin3_gap1"))->Fill(pt);
            } else if (gapSide == 2) {
              registry.get<TH1>(HIST("hPtKaon_etaBin3_gap2"))->Fill(pt);
            }
          } else if (particleType == "Proton") {
            if (gapSide == 0) {
              registry.get<TH1>(HIST("hPtProton_etaBin3_gap0"))->Fill(pt);
            } else if (gapSide == 1) {
              registry.get<TH1>(HIST("hPtProton_etaBin3_gap1"))->Fill(pt);
            } else if (gapSide == 2) {
              registry.get<TH1>(HIST("hPtProton_etaBin3_gap2"))->Fill(pt);
            }
          }
        }
        if (i == 3) {
          if (particleType == "Pion") {
            if (gapSide == 0) {
              registry.get<TH1>(HIST("hPtPion_etaBin4_gap0"))->Fill(pt);
            } else if (gapSide == 1) {
              registry.get<TH1>(HIST("hPtPion_etaBin4_gap1"))->Fill(pt);
            } else if (gapSide == 2) {
              registry.get<TH1>(HIST("hPtPion_etaBin4_gap2"))->Fill(pt);
            }
          } else if (particleType == "Kaon") {
            if (gapSide == 0) {
              registry.get<TH1>(HIST("hPtKaon_etaBin4_gap0"))->Fill(pt);
            } else if (gapSide == 1) {
              registry.get<TH1>(HIST("hPtKaon_etaBin4_gap1"))->Fill(pt);
            } else if (gapSide == 2) {
              registry.get<TH1>(HIST("hPtKaon_etaBin4_gap2"))->Fill(pt);
            }
          } else if (particleType == "Proton") {
            if (gapSide == 0) {
              registry.get<TH1>(HIST("hPtProton_etaBin4_gap0"))->Fill(pt);
            } else if (gapSide == 1) {
              registry.get<TH1>(HIST("hPtProton_etaBin4_gap1"))->Fill(pt);
            } else if (gapSide == 2) {
              registry.get<TH1>(HIST("hPtProton_etaBin4_gap2"))->Fill(pt);
            }
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
