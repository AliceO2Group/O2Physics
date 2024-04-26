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
#include "Common/DataModel/PIDResponse.h"
#include "PWGUD/Tasks/SGTrackSelector.h"
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
struct SGD0Analyzer {
  SGSelector sgSelector;
  Configurable<float> FV0_cut{"FV0", 100., "FV0A threshold"};
  Configurable<float> FT0A_cut{"FT0A", 100., "FT0A threshold"};
  Configurable<float> FT0C_cut{"FT0C", 50., "FT0C threshold"};
  Configurable<float> FDDA_cut{"FDDA", 10000., "FDDA threshold"};
  Configurable<float> FDDC_cut{"FDDC", 10000., "FDDC threshold"};
  Configurable<float> ZDC_cut{"ZDC", 10., "ZDC threshold"};
  // Track Selections
  Configurable<float> PV_cut{"PV_cut", 0.0, "Use Only PV tracks"};
  Configurable<float> dcaZ_cut{"dcaZ_cut", 2.0, "dcaZ cut"};
  Configurable<float> dcaXY_cut{"dcaXY_cut", 2.0, "dcaXY cut (0 for Pt-function)"};
  Configurable<float> tpcChi2_cut{"tpcChi2_cut", 4, "Max tpcChi2NCl"};
  Configurable<float> tpcNClsFindable_cut{"tpcNClsFindable_cut", 70, "Min tpcNClsFindable"};
  Configurable<float> itsChi2_cut{"itsChi2_cut", 36, "Max itsChi2NCl"};
  Configurable<float> eta_cut{"eta_cut", 0.9, "Track Pseudorapidity"};
  // D0 Specific Cuts
  Configurable<float> Ntr_min{"Ntr_min", 2., "Minimum Number of Tracks"};
  Configurable<float> Ntr_max{"Ntr_max", 50., "Maximum Number of Tracks"};
  HistogramRegistry registry{
    "registry",
    {

      {"GapSide", "Gap Side; Entries", {HistType::kTH1F, {{4, -1.5, 2.5}}}},
      {"TrueGapSide", "Gap Side; Entries", {HistType::kTH1F, {{4, -1.5, 2.5}}}},
      {"os_KPi_pT", "K#pi pT (GeV/c); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"os_KPi_eTa", "K#pi eTa (GeV/c); Entries", {HistType::kTH1F, {{100, -1., 1.}}}},
      {"os_KPi_invm", "K#pi Mass (GeV/c^2); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"ss_KPi_pT", "K#pi pT (GeV/c); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"ss_KPi_eTa", "K#pi eTa (GeV/c); Entries", {HistType::kTH1F, {{100, -1., 1.}}}},
      {"ss_KPi_invm", "K#pi Mass (GeV/c^2); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"os_KPi_pT_1", "K#pi pT (GeV/c); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"os_KPi_eTa_1", "K#pi eTa (GeV/c); Entries", {HistType::kTH1F, {{100, -1., 1.}}}},
      {"os_KPi_invm_1", "K#pi Mass (GeV/c^2); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"ss_KPi_pT_1", "K#pi pT (GeV/c); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"ss_KPi_eTa_1", "K#pi eTa (GeV/c); Entries", {HistType::kTH1F, {{100, -1., 1.}}}},
      {"ss_KPi_invm_1", "K#pi Mass (GeV/c^2); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"os_KPi_pT_0", "K#pi pT (GeV/c); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"os_KPi_eTa_0", "K#pi eTa (GeV/c); Entries", {HistType::kTH1F, {{100, -1., 1.}}}},
      {"os_KPi_invm_0", "K#pi Mass (GeV/c^2); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"ss_KPi_pT_0", "K#pi pT (GeV/c); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"ss_KPi_eTa_0", "K#pi eTa (GeV/c); Entries", {HistType::kTH1F, {{100, -1., 1.}}}},
      {"ss_KPi_invm_0", "K#pi Mass (GeV/c^2); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"os_KPi_pT_2", "K#pi pT (GeV/c); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"os_KPi_eTa_2", "K#pi eTa (GeV/c); Entries", {HistType::kTH1F, {{100, -1., 1.}}}},
      {"os_KPi_invm_2", "K#pi Mass (GeV/c^2); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"ss_KPi_pT_2", "K#pi pT (GeV/c); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"ss_KPi_eTa_2", "K#pi eTa (GeV/c); Entries", {HistType::kTH1F, {{100, -1., 1.}}}},
      {"ss_KPi_invm_2", "K#pi Mass (GeV/c^2); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"os_Ntr_KPi_invm_0", "N tracks vs K#pi Mass", {HistType::kTH2F, {{48, 1.5, 49.5}, {750, .7, 2.2}}}},
      {"os_Ntr_KPi_invm_1", "N tracks vs K#pi Mass", {HistType::kTH2F, {{48, 1.5, 49.5}, {750, .7, 2.2}}}},
      {"os_Ntr_KPi_invm_2", "N tracks vs K#pi Mass", {HistType::kTH2F, {{48, 1.5, 49.5}, {750, .7, 2.2}}}}}};
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
    std::vector<float> parameters = {PV_cut, dcaZ_cut, dcaXY_cut, tpcChi2_cut, tpcNClsFindable_cut, itsChi2_cut, eta_cut};
    // int truegapSide = sgSelector.trueGap(collision, *FIT_cut, ZDC_cut);
    int truegapSide = sgSelector.trueGap(collision, FIT_cut[0], FIT_cut[1], FIT_cut[3], ZDC_cut);
    registry.fill(HIST("GapSide"), gapSide);
    registry.fill(HIST("TrueGapSide"), truegapSide);
    gapSide = truegapSide;
    // Look for D0 and D0bar
    if (tracks.size() < Ntr_min || tracks.size() > Ntr_max)
      return;
    for (auto& [t0, t1] : combinations(tracks, tracks)) {
      // PID cut - t0=K, t1=pi
      if (!trackselector(t0, parameters) || !trackselector(t1, parameters))
        continue;
      if (std::abs(t0.tpcNSigmaKa()) < 3 && std::abs(t1.tpcNSigmaPi()) < 3 && std::abs(t0.tofNSigmaKa()) < 3 && std::abs(t1.tofNSigmaPi()) < 3) {
        // Apply pion hypothesis and create pairs
        v0.SetXYZM(t0.px(), t0.py(), t0.pz(), o2::constants::physics::MassKaonCharged);
        v1.SetXYZM(t1.px(), t1.py(), t1.pz(), o2::constants::physics::MassPionCharged);
        v01 = v0 + v1;
        // Opposite sign pairs
        if (t0.sign() != t1.sign()) {
          registry.fill(HIST("os_KPi_pT"), v01.Pt());
          registry.fill(HIST("os_KPi_eTa"), v01.Eta());
          registry.fill(HIST("os_KPi_invm"), v01.M());
          if (gapSide == 0) {
            registry.fill(HIST("os_KPi_pT_0"), v01.Pt());
            registry.fill(HIST("os_KPi_eTa_0"), v01.Eta());
            registry.fill(HIST("os_KPi_invm_0"), v01.M());
            registry.fill(HIST("os_Ntr_KPi_invm_0"), tracks.size(), v01.M());
          }
          if (gapSide == 1) {
            registry.fill(HIST("os_KPi_pT_1"), v01.Pt());
            registry.fill(HIST("os_KPi_eTa_1"), v01.Eta());
            registry.fill(HIST("os_KPi_invm_1"), v01.M());
            registry.fill(HIST("os_Ntr_KPi_invm_1"), tracks.size(), v01.M());
          }
          if (gapSide == 2) {
            registry.fill(HIST("os_KPi_pT_2"), v01.Pt());
            registry.fill(HIST("os_KPi_eTa_2"), v01.Eta());
            registry.fill(HIST("os_KPi_invm_2"), v01.M());
            registry.fill(HIST("os_Ntr_KPi_invm_2"), tracks.size(), v01.M());
          }
        } else if (t0.sign() == t1.sign()) {
          registry.fill(HIST("ss_KPi_pT"), v01.Pt());
          registry.fill(HIST("ss_KPi_eTa"), v01.Eta());
          registry.fill(HIST("ss_KPi_invm"), v01.M());
          if (gapSide == 0) {
            registry.fill(HIST("ss_KPi_pT_0"), v01.Pt());
            registry.fill(HIST("ss_KPi_eTa_0"), v01.Eta());
            registry.fill(HIST("ss_KPi_invm_0"), v01.M());
          }
          if (gapSide == 1) {
            registry.fill(HIST("ss_KPi_pT_1"), v01.Pt());
            registry.fill(HIST("ss_KPi_eTa_1"), v01.Eta());
            registry.fill(HIST("ss_KPi_invm_1"), v01.M());
          }
          if (gapSide == 2) {
            registry.fill(HIST("ss_KPi_pT_2"), v01.Pt());
            registry.fill(HIST("ss_KPi_eTa_2"), v01.Eta());
            registry.fill(HIST("ss_KPi_invm_2"), v01.M());
          }
        }
      } else if (std::abs(t1.tpcNSigmaKa()) < 3 && std::abs(t0.tpcNSigmaPi()) < 3 && std::abs(t1.tofNSigmaKa()) < 3 && std::abs(t0.tofNSigmaPi()) < 3) {
        // Apply pion hypothesis and create pairs
        v0.SetXYZM(t0.px(), t0.py(), t0.pz(), o2::constants::physics::MassPionCharged);
        v1.SetXYZM(t1.px(), t1.py(), t1.pz(), o2::constants::physics::MassKaonCharged);
        v01 = v0 + v1;
        // Opposite sign pairs
        if (t0.sign() != t1.sign()) {
          registry.fill(HIST("os_KPi_pT"), v01.Pt());
          registry.fill(HIST("os_KPi_eTa"), v01.Eta());
          registry.fill(HIST("os_KPi_invm"), v01.M());
          if (gapSide == 0) {
            registry.fill(HIST("os_KPi_pT_0"), v01.Pt());
            registry.fill(HIST("os_KPi_eTa_0"), v01.Eta());
            registry.fill(HIST("os_KPi_invm_0"), v01.M());
            registry.fill(HIST("os_Ntr_KPi_invm_0"), tracks.size(), v01.M());
          }
          if (gapSide == 1) {
            registry.fill(HIST("os_KPi_pT_1"), v01.Pt());
            registry.fill(HIST("os_KPi_eTa_1"), v01.Eta());
            registry.fill(HIST("os_KPi_invm_1"), v01.M());
            registry.fill(HIST("os_Ntr_KPi_invm_1"), tracks.size(), v01.M());
          }
          if (gapSide == 2) {
            registry.fill(HIST("os_KPi_pT_2"), v01.Pt());
            registry.fill(HIST("os_KPi_eTa_2"), v01.Eta());
            registry.fill(HIST("os_KPi_invm_2"), v01.M());
            registry.fill(HIST("os_Ntr_KPi_invm_2"), tracks.size(), v01.M());
          }
        } else if (t0.sign() == t1.sign()) {
          registry.fill(HIST("ss_KPi_pT"), v01.Pt());
          registry.fill(HIST("ss_KPi_eTa"), v01.Eta());
          registry.fill(HIST("ss_KPi_invm"), v01.M());
          if (gapSide == 0) {
            registry.fill(HIST("ss_KPi_pT_0"), v01.Pt());
            registry.fill(HIST("ss_KPi_eTa_0"), v01.Eta());
            registry.fill(HIST("ss_KPi_invm_0"), v01.M());
          }
          if (gapSide == 1) {
            registry.fill(HIST("ss_KPi_pT_1"), v01.Pt());
            registry.fill(HIST("ss_KPi_eTa_1"), v01.Eta());
            registry.fill(HIST("ss_KPi_invm_1"), v01.M());
          }
          if (gapSide == 2) {
            registry.fill(HIST("ss_KPi_pT_2"), v01.Pt());
            registry.fill(HIST("ss_KPi_eTa_2"), v01.Eta());
            registry.fill(HIST("ss_KPi_invm_2"), v01.M());
          }
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<SGD0Analyzer>(cfgc)};
}
