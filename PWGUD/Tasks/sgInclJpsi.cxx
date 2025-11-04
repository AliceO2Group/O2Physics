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

#include "TLorentzVector.h"
#include <TString.h>

#include <iostream>
using namespace std;
using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
struct SGInclJpsi {
  SGSelector sgSelector;
  Service<o2::framework::O2DatabasePDG> pdg;
  Configurable<float> FV0_cut{"FV0", 100., "FV0A threshold"};
  Configurable<float> FT0A_cut{"FT0A", 100., "FT0A threshold"};
  Configurable<float> FT0C_cut{"FT0C", 50., "FT0C threshold"};
  Configurable<float> FDDA_cut{"FDDA", 10000., "FDDA threshold"};
  Configurable<float> FDDC_cut{"FDDC", 10000., "FDDC threshold"};
  Configurable<float> ZDC_cut{"ZDC", 1., "ZDC threshold"};
  // Track Selections
  Configurable<float> PV_cut{"PV_cut", 1.0, "Use Only PV tracks"};
  Configurable<float> dcaZ_cut{"dcaZ_cut", 2.0, "dcaZ cut"};
  Configurable<float> dcaXY_cut{"dcaXY_cut", 0, "dcaXY cut (0 for Pt-function)"};
  Configurable<float> tpcChi2_cut{"tpcChi2_cut", 4, "Max tpcChi2NCl"};
  Configurable<float> tpcNClsFindable_cut{"tpcNClsFindable_cut", 70, "Min tpcNClsFindable"};
  Configurable<float> itsChi2_cut{"itsChi2_cut", 36, "Max itsChi2NCl"};
  Configurable<float> eta_cut{"eta_cut", 0.9, "Track Pseudorapidity"};
  Configurable<float> pt_cut{"pt_cut", 0.1, "Track Pt"};
  // JPsi Specific Cuts
  Configurable<float> nsigmatpc_cut{"nsigmatpc", 3.0, "nsigma tpc cut"};
  Configurable<float> nsigmatof_cut{"nsigmatof", 9.0, "nsigma tof cut"};
  Configurable<bool> use_tof{"Use_TOF", true, "TOF PID"};
  HistogramRegistry registry{
    "registry",
    {

      {"GapSide", "Gap Side; Entries", {HistType::kTH1F, {{4, -1.5, 2.5}}}},
      {"TrueGapSide", "Gap Side; Entries", {HistType::kTH1F, {{4, -1.5, 2.5}}}},
      {"Ntr", "Ntracks", {HistType::kTH1F, {{50, -.5, 49.5}}}},
      {"Ntr_0", "Ntracks", {HistType::kTH1F, {{50, -.5, 49.5}}}},
      {"Ntr_1", "Ntracks", {HistType::kTH1F, {{50, -.5, 49.5}}}},
      {"Ntr_2", "Ntracks", {HistType::kTH1F, {{50, -.5, 49.5}}}},
      {"Ntr_3", "Ntracks", {HistType::kTH1F, {{50, -.5, 49.5}}}},
      {"os_mm_pT", "pT (GeV/c); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"os_mm_eTa", "eTa (GeV/c); Entries", {HistType::kTH1F, {{100, -1., 1.}}}},
      {"os_mm_invm", "Mass (GeV/c^2); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"ss_mm_pT", "pT (GeV/c); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"ss_mm_eTa", "eTa (GeV/c); Entries", {HistType::kTH1F, {{100, -1., 1.}}}},
      {"ss_mm_invm", "Mass (GeV/c^2); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"os_mm_pT_0", "pT (GeV/c); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"os_mm_eTa_0", "eTa (GeV/c); Entries", {HistType::kTH1F, {{100, -1., 1.}}}},
      {"os_mm_invm_0", "Mass (GeV/c^2); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"ss_mm_pT_0", "pT (GeV/c); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"ss_mm_eTa_0", "eTa (GeV/c); Entries", {HistType::kTH1F, {{100, -1., 1.}}}},
      {"ss_mm_invm_0", "Mass (GeV/c^2); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"os_mm_pT_1", "pT (GeV/c); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"os_mm_eTa_1", "eTa (GeV/c); Entries", {HistType::kTH1F, {{100, -1., 1.}}}},
      {"os_mm_invm_1", "Mass (GeV/c^2); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"ss_mm_pT_1", "pT (GeV/c); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"ss_mm_eTa_1", "eTa (GeV/c); Entries", {HistType::kTH1F, {{100, -1., 1.}}}},
      {"ss_mm_invm_1", "Mass (GeV/c^2); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"os_mm_pT_2", "pT (GeV/c); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"os_mm_eTa_2", "eTa (GeV/c); Entries", {HistType::kTH1F, {{100, -1., 1.}}}},
      {"os_mm_invm_2", "Mass (GeV/c^2); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"ss_mm_pT_2", "pT (GeV/c); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"ss_mm_eTa_2", "eTa (GeV/c); Entries", {HistType::kTH1F, {{100, -1., 1.}}}},
      {"ss_mm_invm_2", "Mass (GeV/c^2); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"os_mm_pT_3", "pT (GeV/c); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"os_mm_eTa_3", "eTa (GeV/c); Entries", {HistType::kTH1F, {{100, -1., 1.}}}},
      {"os_mm_invm_3", "Mass (GeV/c^2); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"ss_mm_pT_3", "pT (GeV/c); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"ss_mm_eTa_3", "eTa (GeV/c); Entries", {HistType::kTH1F, {{100, -1., 1.}}}},
      {"ss_mm_invm_3", "Mass (GeV/c^2); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"os_Ntr_mm_invm_0", "N tracks vs #mu#mu Mass", {HistType::kTH2F, {{48, 1.5, 49.5}, {600, 2., 5.}}}},
      {"os_Ntr_mm_invm_1", "N tracks vs #mu#mu Mass", {HistType::kTH2F, {{48, 1.5, 49.5}, {600, 2., 5.}}}},
      {"os_Ntr_mm_invm_2", "N tracks vs #mu#mu Mass", {HistType::kTH2F, {{48, 1.5, 49.5}, {600, 2., 5.}}}},
      {"os_Ntr_mm_invm_3", "N tracks vs #mu#mu Mass", {HistType::kTH2F, {{48, 1.5, 49.5}, {600, 2., 5.}}}},
      {"ss_Ntr_mm_invm_0", "N tracks vs #mu#mu Mass", {HistType::kTH2F, {{48, 1.5, 49.5}, {600, 2., 5.}}}},
      {"ss_Ntr_mm_invm_1", "N tracks vs #mu#mu Mass", {HistType::kTH2F, {{48, 1.5, 49.5}, {600, 2., 5.}}}},
      {"ss_Ntr_mm_invm_2", "N tracks vs #mu#mu Mass", {HistType::kTH2F, {{48, 1.5, 49.5}, {600, 2., 5.}}}},
      {"ss_Ntr_mm_invm_3", "N tracks vs #mu#mu Mass", {HistType::kTH2F, {{48, 1.5, 49.5}, {600, 2., 5.}}}},
      {"sss_Ntr_mm_invm_0", "N tracks vs #mu#mu Mass", {HistType::kTH2F, {{48, 1.5, 49.5}, {600, 2., 5.}}}},
      {"sss_Ntr_mm_invm_1", "N tracks vs #mu#mu Mass", {HistType::kTH2F, {{48, 1.5, 49.5}, {600, 2., 5.}}}},
      {"sss_Ntr_mm_invm_2", "N tracks vs #mu#mu Mass", {HistType::kTH2F, {{48, 1.5, 49.5}, {600, 2., 5.}}}},
      {"sss_Ntr_mm_invm_3", "N tracks vs #mu#mu Mass", {HistType::kTH2F, {{48, 1.5, 49.5}, {600, 2., 5.}}}},
      {"os_ee_pT", "pT (GeV/c); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"os_ee_eTa", "eTa (GeV/c); Entries", {HistType::kTH1F, {{100, -1., 1.}}}},
      {"os_ee_invm", "Mass (GeV/c^2); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"ss_ee_pT", "pT (GeV/c); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"ss_ee_eTa", "eTa (GeV/c); Entries", {HistType::kTH1F, {{100, -1., 1.}}}},
      {"ss_ee_invm", "Mass (GeV/c^2); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"os_ee_pT_0", "pT (GeV/c); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"os_ee_eTa_0", "eTa (GeV/c); Entries", {HistType::kTH1F, {{100, -1., 1.}}}},
      {"os_ee_invm_0", "Mass (GeV/c^2); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"ss_ee_pT_0", "pT (GeV/c); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"ss_ee_eTa_0", "eTa (GeV/c); Entries", {HistType::kTH1F, {{100, -1., 1.}}}},
      {"ss_ee_invm_0", "Mass (GeV/c^2); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"os_ee_pT_1", "pT (GeV/c); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"os_ee_eTa_1", "eTa (GeV/c); Entries", {HistType::kTH1F, {{100, -1., 1.}}}},
      {"os_ee_invm_1", "Mass (GeV/c^2); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"ss_ee_pT_1", "pT (GeV/c); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"ss_ee_eTa_1", "eTa (GeV/c); Entries", {HistType::kTH1F, {{100, -1., 1.}}}},
      {"ss_ee_invm_1", "Mass (GeV/c^2); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"os_ee_pT_2", "pT (GeV/c); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"os_ee_eTa_2", "eTa (GeV/c); Entries", {HistType::kTH1F, {{100, -1., 1.}}}},
      {"os_ee_invm_2", "Mass (GeV/c^2); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"ss_ee_pT_2", "pT (GeV/c); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"ss_ee_eTa_2", "eTa (GeV/c); Entries", {HistType::kTH1F, {{100, -1., 1.}}}},
      {"ss_ee_invm_2", "Mass (GeV/c^2); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"os_ee_pT_3", "pT (GeV/c); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"os_ee_eTa_3", "eTa (GeV/c); Entries", {HistType::kTH1F, {{100, -1., 1.}}}},
      {"os_ee_invm_3", "Mass (GeV/c^2); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"ss_ee_pT_3", "pT (GeV/c); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"ss_ee_eTa_3", "eTa (GeV/c); Entries", {HistType::kTH1F, {{100, -1., 1.}}}},
      {"ss_ee_invm_3", "Mass (GeV/c^2); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"os_ee_pt_invm_0", "N tracks vs ee Mass", {HistType::kTH2F, {{500, 0, 10}, {600, 2., 5.}}}},
      {"os_ee_pt_invm_1", "N tracks vs ee Mass", {HistType::kTH2F, {{500, 0, 10}, {600, 2., 5.}}}},
      {"os_ee_pt_invm_2", "N tracks vs ee Mass", {HistType::kTH2F, {{500, 0, 10}, {600, 2., 5.}}}},
      {"os_ee_pt_invm_3", "N tracks vs ee Mass", {HistType::kTH2F, {{500, 0, 10}, {600, 2., 5.}}}},
      {"os_mm_pt_invm_0", "N tracks vs ee Mass", {HistType::kTH2F, {{500, 0, 10}, {600, 2., 5.}}}},
      {"os_mm_pt_invm_1", "N tracks vs ee Mass", {HistType::kTH2F, {{500, 0, 10}, {600, 2., 5.}}}},
      {"os_mm_pt_invm_2", "N tracks vs ee Mass", {HistType::kTH2F, {{500, 0, 10}, {600, 2., 5.}}}},
      {"os_mm_pt_invm_3", "N tracks vs ee Mass", {HistType::kTH2F, {{500, 0, 10}, {600, 2., 5.}}}},
      {"os_Ntr_ee_invm_0", "N tracks vs ee Mass", {HistType::kTH2F, {{48, 1.5, 49.5}, {600, 2., 5.}}}},
      {"os_Ntr_ee_invm_1", "N tracks vs ee Mass", {HistType::kTH2F, {{48, 1.5, 49.5}, {600, 2., 5.}}}},
      {"os_Ntr_ee_invm_2", "N tracks vs ee Mass", {HistType::kTH2F, {{48, 1.5, 49.5}, {600, 2., 5.}}}},
      {"os_Ntr_ee_invm_3", "N tracks vs ee Mass", {HistType::kTH2F, {{48, 1.5, 49.5}, {600, 2., 5.}}}},
      {"ss_ee_pt_invm_0", "N tracks vs ee Mass", {HistType::kTH2F, {{500, 0, 10}, {600, 2., 5.}}}},
      {"ss_ee_pt_invm_1", "N tracks vs ee Mass", {HistType::kTH2F, {{500, 0, 10}, {600, 2., 5.}}}},
      {"ss_ee_pt_invm_2", "N tracks vs ee Mass", {HistType::kTH2F, {{500, 0, 10}, {600, 2., 5.}}}},
      {"ss_ee_pt_invm_3", "N tracks vs ee Mass", {HistType::kTH2F, {{500, 0, 10}, {600, 2., 5.}}}},
      {"ss_mm_pt_invm_0", "N tracks vs ee Mass", {HistType::kTH2F, {{500, 0, 10}, {600, 2., 5.}}}},
      {"ss_mm_pt_invm_1", "N tracks vs ee Mass", {HistType::kTH2F, {{500, 0, 10}, {600, 2., 5.}}}},
      {"ss_mm_pt_invm_2", "N tracks vs ee Mass", {HistType::kTH2F, {{500, 0, 10}, {600, 2., 5.}}}},
      {"ss_mm_pt_invm_3", "N tracks vs ee Mass", {HistType::kTH2F, {{500, 0, 10}, {600, 2., 5.}}}},
      {"ss_Ntr_ee_invm_0", "N tracks vs ee Mass", {HistType::kTH2F, {{48, 1.5, 49.5}, {600, 2., 5.}}}},
      {"ss_Ntr_ee_invm_1", "N tracks vs ee Mass", {HistType::kTH2F, {{48, 1.5, 49.5}, {600, 2., 5.}}}},
      {"ss_Ntr_ee_invm_2", "N tracks vs ee Mass", {HistType::kTH2F, {{48, 1.5, 49.5}, {600, 2., 5.}}}},
      {"ss_Ntr_ee_invm_3", "N tracks vs ee Mass", {HistType::kTH2F, {{48, 1.5, 49.5}, {600, 2., 5.}}}},
      {"sss_ee_pt_invm_0", "N tracks vs ee Mass", {HistType::kTH2F, {{500, 0, 10}, {600, 2., 5.}}}},
      {"sss_ee_pt_invm_1", "N tracks vs ee Mass", {HistType::kTH2F, {{500, 0, 10}, {600, 2., 5.}}}},
      {"sss_ee_pt_invm_2", "N tracks vs ee Mass", {HistType::kTH2F, {{500, 0, 10}, {600, 2., 5.}}}},
      {"sss_ee_pt_invm_3", "N tracks vs ee Mass", {HistType::kTH2F, {{500, 0, 10}, {600, 2., 5.}}}},
      {"sss_mm_pt_invm_0", "N tracks vs ee Mass", {HistType::kTH2F, {{500, 0, 10}, {600, 2., 5.}}}},
      {"sss_mm_pt_invm_1", "N tracks vs ee Mass", {HistType::kTH2F, {{500, 0, 10}, {600, 2., 5.}}}},
      {"sss_mm_pt_invm_2", "N tracks vs ee Mass", {HistType::kTH2F, {{500, 0, 10}, {600, 2., 5.}}}},
      {"sss_mm_pt_invm_3", "N tracks vs ee Mass", {HistType::kTH2F, {{500, 0, 10}, {600, 2., 5.}}}},
      {"sss_Ntr_ee_invm_0", "N tracks vs ee Mass", {HistType::kTH2F, {{48, 1.5, 49.5}, {600, 2., 5.}}}},
      {"sss_Ntr_ee_invm_1", "N tracks vs ee Mass", {HistType::kTH2F, {{48, 1.5, 49.5}, {600, 2., 5.}}}},
      {"sss_Ntr_ee_invm_2", "N tracks vs ee Mass", {HistType::kTH2F, {{48, 1.5, 49.5}, {600, 2., 5.}}}},
      {"sss_Ntr_ee_invm_3", "N tracks vs ee Mass", {HistType::kTH2F, {{48, 1.5, 49.5}, {600, 2., 5.}}}},
    }};
  using udtracks = soa::Join<aod::UDTracks, aod::UDTracksExtra, aod::UDTracksPID>;
  using udtracksfull = soa::Join<aod::UDTracks, aod::UDTracksPID, aod::UDTracksExtra, aod::UDTracksFlags, aod::UDTracksDCA>;
  using UDCollisionsFull = soa::Join<aod::UDCollisions, aod::SGCollisions, aod::UDCollisionsSels, aod::UDZdcsReduced>; //
  using UDCollisionFull = UDCollisionsFull::iterator;

  void process(UDCollisionFull const& collision, udtracksfull const& tracks)
  {
    TLorentzVector a;
    int gapSide = collision.gapSide();
    if (gapSide < 0 || gapSide > 2)
      return;
    // Single gap either side
    TLorentzVector v0;
    TLorentzVector v1;
    TLorentzVector v01;
    //  int truegapSide = sgSelector.trueGap(collision);
    // int truegapSide = sgSelector.trueGap(collision, FV0_cut, ZDC_cut);
    float FIT_cut[5] = {FV0_cut, FT0A_cut, FT0C_cut, FDDA_cut, FDDC_cut};
    std::vector<float> parameters = {PV_cut, dcaZ_cut, dcaXY_cut, tpcChi2_cut, tpcNClsFindable_cut, itsChi2_cut, eta_cut, pt_cut};
    // int truegapSide = sgSelector.trueGap(collision, *FIT_cut, ZDC_cut);
    int truegapSide = sgSelector.trueGap(collision, FIT_cut[0], FIT_cut[1], FIT_cut[2], ZDC_cut);
    registry.fill(HIST("GapSide"), gapSide);
    registry.fill(HIST("TrueGapSide"), truegapSide);
    gapSide = truegapSide;
    int pvtracks = 0;
    for (auto& t0 : tracks) {
      if (trackselector(t0, parameters) && t0.isPVContributor())
        pvtracks++;
    }
    registry.fill(HIST("Ntr"), pvtracks);
    if (gapSide == 0)
      registry.fill(HIST("Ntr_0"), pvtracks);
    if (gapSide == 1)
      registry.fill(HIST("Ntr_1"), pvtracks);
    if (gapSide == 2)
      registry.fill(HIST("Ntr_2"), pvtracks);
    if (gapSide == -1)
      registry.fill(HIST("Ntr_3"), pvtracks);
    // Look for D0 and D0bar
    for (auto& [t0, t1] : combinations(tracks, tracks)) {
      // PID cut - t0=K, t1=pi
      if (!trackselector(t0, parameters) || !trackselector(t1, parameters))
        continue;
      if (selectionPIDMuon(t0, use_tof, nsigmatpc_cut, nsigmatof_cut) && selectionPIDMuon(t1, use_tof, nsigmatpc_cut, nsigmatof_cut)) {
        // Apply pion hypothesis and create pairs
        v0.SetXYZM(t0.px(), t0.py(), t0.pz(), o2::constants::physics::MassMuon);
        v1.SetXYZM(t1.px(), t1.py(), t1.pz(), o2::constants::physics::MassMuon);
        v01 = v0 + v1;
        // Opposite sign pairs
        if (t0.sign() != t1.sign()) {
          registry.fill(HIST("os_mm_pT"), v01.Pt());
          registry.fill(HIST("os_mm_eTa"), v01.Eta());
          registry.fill(HIST("os_mm_invm"), v01.M());
          if (gapSide == 0) {
            registry.fill(HIST("os_mm_pT_0"), v01.Pt());
            registry.fill(HIST("os_mm_eTa_0"), v01.Eta());
            registry.fill(HIST("os_mm_invm_0"), v01.M());
            registry.fill(HIST("os_Ntr_mm_invm_0"), tracks.size(), v01.M());
            registry.fill(HIST("os_mm_pt_invm_0"), v01.Pt(), v01.M());
            registry.fill(HIST("sss_Ntr_mm_invm_0"), tracks.size(), v01.M());
            registry.fill(HIST("sss_mm_pt_invm_0"), v01.Pt(), v01.M());
          } else if (gapSide == 1) {
            registry.fill(HIST("os_mm_pT_1"), v01.Pt());
            registry.fill(HIST("os_mm_eTa_1"), v01.Eta());
            registry.fill(HIST("os_mm_invm_1"), v01.M());
            registry.fill(HIST("os_Ntr_mm_invm_1"), tracks.size(), v01.M());
            registry.fill(HIST("os_mm_pt_invm_1"), v01.Pt(), v01.M());
            registry.fill(HIST("sss_Ntr_mm_invm_1"), tracks.size(), v01.M());
            registry.fill(HIST("sss_mm_pt_invm_1"), v01.Pt(), v01.M());
          } else if (gapSide == 2) {
            registry.fill(HIST("os_mm_pT_2"), v01.Pt());
            registry.fill(HIST("os_mm_eTa_2"), v01.Eta());
            registry.fill(HIST("os_mm_invm_2"), v01.M());
            registry.fill(HIST("os_Ntr_mm_invm_2"), tracks.size(), v01.M());
            registry.fill(HIST("os_mm_pt_invm_2"), v01.Pt(), v01.M());
            registry.fill(HIST("sss_Ntr_mm_invm_2"), tracks.size(), v01.M());
            registry.fill(HIST("sss_mm_pt_invm_2"), v01.Pt(), v01.M());
          } else {
            registry.fill(HIST("os_mm_pT_3"), v01.Pt());
            registry.fill(HIST("os_mm_eTa_3"), v01.Eta());
            registry.fill(HIST("os_mm_invm_3"), v01.M());
            registry.fill(HIST("os_Ntr_mm_invm_3"), tracks.size(), v01.M());
            registry.fill(HIST("os_mm_pt_invm_3"), v01.Pt(), v01.M());
            registry.fill(HIST("sss_Ntr_mm_invm_3"), tracks.size(), v01.M());
            registry.fill(HIST("sss_mm_pt_invm_3"), v01.Pt(), v01.M());
          }
        } else if (t0.sign() == t1.sign()) {
          registry.fill(HIST("ss_mm_pT"), v01.Pt());
          registry.fill(HIST("ss_mm_eTa"), v01.Eta());
          registry.fill(HIST("ss_mm_invm"), v01.M());
          if (gapSide == 0) {
            registry.fill(HIST("ss_mm_pT_0"), v01.Pt());
            registry.fill(HIST("ss_mm_eTa_0"), v01.Eta());
            registry.fill(HIST("ss_mm_invm_0"), v01.M());
            registry.fill(HIST("ss_Ntr_mm_invm_0"), tracks.size(), v01.M());
            registry.fill(HIST("ss_mm_pt_invm_0"), v01.Pt(), v01.M());
            registry.fill(HIST("sss_Ntr_mm_invm_0"), tracks.size(), v01.M(), -1);
            registry.fill(HIST("sss_mm_pt_invm_0"), v01.Pt(), v01.M(), -1);
          }
          if (gapSide == 1) {
            registry.fill(HIST("ss_mm_pT_1"), v01.Pt());
            registry.fill(HIST("ss_mm_eTa_1"), v01.Eta());
            registry.fill(HIST("ss_mm_invm_1"), v01.M());
            registry.fill(HIST("ss_Ntr_mm_invm_1"), tracks.size(), v01.M());
            registry.fill(HIST("ss_mm_pt_invm_1"), v01.Pt(), v01.M());
            registry.fill(HIST("sss_Ntr_mm_invm_1"), tracks.size(), v01.M(), -1);
            registry.fill(HIST("sss_mm_pt_invm_1"), v01.Pt(), v01.M(), -1);
          }
          if (gapSide == 2) {
            registry.fill(HIST("ss_mm_pT_2"), v01.Pt());
            registry.fill(HIST("ss_mm_eTa_2"), v01.Eta());
            registry.fill(HIST("ss_mm_invm_2"), v01.M());
            registry.fill(HIST("ss_Ntr_mm_invm_2"), tracks.size(), v01.M());
            registry.fill(HIST("ss_mm_pt_invm_2"), v01.Pt(), v01.M());
            registry.fill(HIST("sss_Ntr_mm_invm_2"), tracks.size(), v01.M(), -1);
            registry.fill(HIST("sss_mm_pt_invm_2"), v01.Pt(), v01.M(), -1);
          }
          if (gapSide == -1) {
            registry.fill(HIST("ss_mm_pT_3"), v01.Pt());
            registry.fill(HIST("ss_mm_eTa_3"), v01.Eta());
            registry.fill(HIST("ss_mm_invm_3"), v01.M());
            registry.fill(HIST("ss_Ntr_mm_invm_3"), tracks.size(), v01.M());
            registry.fill(HIST("ss_mm_pt_invm_3"), v01.Pt(), v01.M());
            registry.fill(HIST("sss_Ntr_mm_invm_3"), tracks.size(), v01.M(), -1);
            registry.fill(HIST("sss_mm_pt_invm_3"), v01.Pt(), v01.M(), -1);
          }
        }
      }
      if (selectionPIDElec(t0, use_tof, nsigmatpc_cut, nsigmatof_cut) && selectionPIDElec(t1, use_tof, nsigmatpc_cut, nsigmatof_cut)) {
        // Apply pion hypothesis and create pairs
        v0.SetXYZM(t0.px(), t0.py(), t0.pz(), o2::constants::physics::MassElectron);
        v1.SetXYZM(t1.px(), t1.py(), t1.pz(), o2::constants::physics::MassElectron);
        v01 = v0 + v1;
        // Opposite sign pairs
        if (t0.sign() != t1.sign()) {
          registry.fill(HIST("os_ee_pT"), v01.Pt());
          registry.fill(HIST("os_ee_eTa"), v01.Eta());
          registry.fill(HIST("os_ee_invm"), v01.M());
          if (gapSide == 0) {
            registry.fill(HIST("os_ee_pT_0"), v01.Pt());
            registry.fill(HIST("os_ee_eTa_0"), v01.Eta());
            registry.fill(HIST("os_ee_invm_0"), v01.M());
            registry.fill(HIST("os_Ntr_ee_invm_0"), tracks.size(), v01.M());
            registry.fill(HIST("os_ee_pt_invm_0"), v01.Pt(), v01.M());
            registry.fill(HIST("sss_ee_pt_invm_0"), v01.Pt(), v01.M());
            registry.fill(HIST("sss_Ntr_ee_invm_0"), tracks.size(), v01.M());
          } else if (gapSide == 1) {
            registry.fill(HIST("os_ee_pT_1"), v01.Pt());
            registry.fill(HIST("os_ee_eTa_1"), v01.Eta());
            registry.fill(HIST("os_ee_invm_1"), v01.M());
            registry.fill(HIST("os_Ntr_ee_invm_1"), tracks.size(), v01.M());
            registry.fill(HIST("os_ee_pt_invm_1"), v01.Pt(), v01.M());
            registry.fill(HIST("sss_Ntr_ee_invm_1"), tracks.size(), v01.M());
            registry.fill(HIST("sss_ee_pt_invm_1"), v01.Pt(), v01.M());
          } else if (gapSide == 2) {
            registry.fill(HIST("os_ee_pT_2"), v01.Pt());
            registry.fill(HIST("os_ee_eTa_2"), v01.Eta());
            registry.fill(HIST("os_ee_invm_2"), v01.M());
            registry.fill(HIST("os_Ntr_ee_invm_2"), tracks.size(), v01.M());
            registry.fill(HIST("os_ee_pt_invm_2"), v01.Pt(), v01.M());
            registry.fill(HIST("sss_Ntr_ee_invm_2"), tracks.size(), v01.M());
            registry.fill(HIST("sss_ee_pt_invm_2"), v01.Pt(), v01.M());
          } else {
            registry.fill(HIST("os_ee_pT_3"), v01.Pt());
            registry.fill(HIST("os_ee_eTa_3"), v01.Eta());
            registry.fill(HIST("os_ee_invm_3"), v01.M());
            registry.fill(HIST("os_Ntr_ee_invm_3"), tracks.size(), v01.M());
            registry.fill(HIST("os_ee_pt_invm_3"), v01.Pt(), v01.M());
            registry.fill(HIST("sss_Ntr_ee_invm_3"), tracks.size(), v01.M());
            registry.fill(HIST("sss_ee_pt_invm_3"), v01.Pt(), v01.M());
          }
        } else if (t0.sign() == t1.sign()) {
          registry.fill(HIST("ss_ee_pT"), v01.Pt());
          registry.fill(HIST("ss_ee_eTa"), v01.Eta());
          registry.fill(HIST("ss_ee_invm"), v01.M());
          if (gapSide == 0) {
            registry.fill(HIST("ss_ee_pT_0"), v01.Pt());
            registry.fill(HIST("ss_ee_eTa_0"), v01.Eta());
            registry.fill(HIST("ss_ee_invm_0"), v01.M());
            registry.fill(HIST("ss_Ntr_ee_invm_0"), tracks.size(), v01.M());
            registry.fill(HIST("ss_ee_pt_invm_0"), v01.Pt(), v01.M());
            registry.fill(HIST("sss_Ntr_ee_invm_0"), tracks.size(), v01.M(), -1);
            registry.fill(HIST("sss_ee_pt_invm_0"), v01.Pt(), v01.M(), -1);
          }
          if (gapSide == 1) {
            registry.fill(HIST("ss_ee_pT_1"), v01.Pt());
            registry.fill(HIST("ss_ee_eTa_1"), v01.Eta());
            registry.fill(HIST("ss_ee_invm_1"), v01.M());
            registry.fill(HIST("ss_Ntr_ee_invm_1"), tracks.size(), v01.M());
            registry.fill(HIST("ss_ee_pt_invm_1"), v01.Pt(), v01.M());
            registry.fill(HIST("sss_Ntr_ee_invm_1"), tracks.size(), v01.M(), -1);
            registry.fill(HIST("sss_ee_pt_invm_1"), v01.Pt(), v01.M(), -1);
          }
          if (gapSide == 2) {
            registry.fill(HIST("ss_ee_pT_2"), v01.Pt());
            registry.fill(HIST("ss_ee_eTa_2"), v01.Eta());
            registry.fill(HIST("ss_ee_invm_2"), v01.M());
            registry.fill(HIST("ss_Ntr_ee_invm_2"), tracks.size(), v01.M());
            registry.fill(HIST("ss_ee_pt_invm_2"), v01.Pt(), v01.M());
            registry.fill(HIST("sss_Ntr_ee_invm_2"), tracks.size(), v01.M(), -1);
            registry.fill(HIST("sss_ee_pt_invm_2"), v01.Pt(), v01.M(), -1);
          }
          if (gapSide == -1) {
            registry.fill(HIST("ss_ee_pT_3"), v01.Pt());
            registry.fill(HIST("ss_ee_eTa_3"), v01.Eta());
            registry.fill(HIST("ss_ee_invm_3"), v01.M());
            registry.fill(HIST("ss_Ntr_ee_invm_3"), tracks.size(), v01.M());
            registry.fill(HIST("ss_ee_pt_invm_3"), v01.Pt(), v01.M());
            registry.fill(HIST("sss_Ntr_ee_invm_3"), tracks.size(), v01.M(), -1);
            registry.fill(HIST("sss_ee_pt_invm_3"), v01.Pt(), v01.M(), -1);
          }
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<SGInclJpsi>(cfgc)};
}
