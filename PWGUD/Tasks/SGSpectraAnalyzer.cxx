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

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "iostream"
#include "PWGUD/DataModel/UDTables.h"
#include "PWGUD/Core/SGSelector.h"
//#include "Common/DataModel/PIDResponse.h"
//#include "PWGUD/Core/RLhelper.h"
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
struct SGSpectraAnalyzer {
  SGSelector sgSelector;
  Configurable<float> FV0_cut{"FV0", 100., "FV0A threshold"};
  Configurable<float> ZDC_cut{"ZDC", 10., "ZDC threshold"};
  Configurable<float> eta_cut{"Eta", 0.9, "Eta cut"};
  Configurable<bool> use_tof{"Use_TOF", true, "TOF PID"};
  HistogramRegistry registry{
    "registry",
    {// Pion histograms for each eta bin and gapSide
     {"ITS_Cluster_nonPV", "ITS Cluster Size", {HistType::kTH1F, {{140, -.5, 139.5}}}},
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

  template <typename T>
  int trackselector(const T& track, bool use_tof)
  {
    TLorentzVector a;
    a.SetXYZM(track.px(), track.py(), track.pz(), mpion);
    if (std::abs(track.dcaZ()) > 2.)
      return 0;
    if (std::abs(track.dcaXY()) > .0105 + .035 / pow(a.Pt(), 1.1))
      return 0;
    if (track.tpcChi2NCl() > 4)
      return 0;
    if (track.tpcNClsFindable() < 70)
      return 0;
    if (track.itsChi2NCl() > 36)
      return 0;
    return 1;
  }
  template <typename T>
  int trackpid(const T& track, bool use_tof)
  {
    int pid = 0;
    float pi, ka, pr;
    float tpi, tka, tpr;
    pi = std::abs(track.tpcNSigmaPi());
    ka = std::abs(track.tpcNSigmaKa());
    pr = std::abs(track.tpcNSigmaPr());
    if (pi < 1. && pi < ka && pi < pr)
      pid = 1;
    else if (ka < 1. && ka < pi && ka < pr)
      pid = 2;
    else if (pr < 1. && pr < pi && pr < ka)
      pid = 3;
    if (use_tof && track.tofChi2() > -1) {
      tpi = std::abs(track.tofNSigmaPi());
      tka = std::abs(track.tofNSigmaKa());
      tpr = std::abs(track.tofNSigmaPr());
      if (std::sqrt(pi * pi + tpi * tpi) < 2 && std::sqrt(pi * pi + tpi * tpi) < std::sqrt(ka * ka + tka * tka) && std::sqrt(pi * pi + tpi * tpi) < std::sqrt(pr * pr + tpr * tpr))
        pid = 1;
      else if (std::sqrt(ka * ka + tka * tka) < 2 && std::sqrt(pi * pi + tpi * tpi) > std::sqrt(ka * ka + tka * tka) && std::sqrt(ka * ka + tka * tka) < std::sqrt(pr * pr + tpr * tpr))
        pid = 2;
      else if (std::sqrt(pr * pr + tpr * tpr) < 2 && std::sqrt(pr * pr + tpr * tpr) < std::sqrt(ka * ka + tka * tka) && std::sqrt(pi * pi + tpi * tpi) > std::sqrt(pr * pr + tpr * tpr))
        pid = 3;
    }
    return pid;
  }
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
    TLorentzVector a;
    int gapSide = collision.gapSide();
    int truegapSide = sgSelector.trueGap(collision, FV0_cut, ZDC_cut);
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
        if (trackselector(track, use_tof)) {
          // if (ispion(track, use_tof)) {
          if (trackpid(track, use_tof) == 1) {
            a.SetXYZM(track.px(), track.py(), track.pz(), mpion);
            if (std::abs(a.Eta()) < eta_cut)
              fillHistograms("Pion", a.Pt(), a.Eta(), gapSide);
          }
          //      if (iskaon(track, use_tof)) {
          if (trackpid(track, use_tof) == 2) {
            a.SetXYZM(track.px(), track.py(), track.pz(), mkaon);
            if (std::abs(a.Eta()) < eta_cut)
              fillHistograms("Kaon", a.Pt(), a.Eta(), gapSide);
          }
          //    if (isproton(track, use_tof)) {
          if (trackpid(track, use_tof) == 3) {
            a.SetXYZM(track.px(), track.py(), track.pz(), mproton);
            if (std::abs(a.Eta()) < eta_cut)
              fillHistograms("Proton", a.Pt(), a.Eta(), gapSide);
          }
        }
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
