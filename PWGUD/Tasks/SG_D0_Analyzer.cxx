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
struct SG_D0_Analyzer {
  SGSelector sgSelector;
  Configurable<float> FV0_cut{"FV0", 100., "FV0A threshold"};
  Configurable<float> FT0A_cut{"FT0A", 100., "FT0A threshold"};
  Configurable<float> FT0C_cut{"FT0C", 50., "FT0C threshold"};
  Configurable<float> FDDA_cut{"FDDA", 10000., "FDDA threshold"};
  Configurable<float> FDDC_cut{"FDDC", 10000., "FDDC threshold"};
  Configurable<float> ZDC_cut{"ZDC", 10., "ZDC threshold"};
  HistogramRegistry registry{
    "registry",
    {

      {"GapSide", "Gap Side; Entries", {HistType::kTH1F, {{4, -1.5, 2.5}}}},
      {"TrueGapSide", "Gap Side; Entries", {HistType::kTH1F, {{4, -1.5, 2.5}}}},
      {"os_KPi_pT", "#K#Pi pT (GeV/c); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"os_KPi_eTa", "#K#Pi eTa (GeV/c); Entries", {HistType::kTH1F, {{100, -1., 1.}}}},
      {"os_KPi_invm", "#K#Pi Mass (GeV/c^2); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"ss_KPi_pT", "#K#Pi pT (GeV/c); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"ss_KPi_eTa", "#K#Pi eTa (GeV/c); Entries", {HistType::kTH1F, {{100, -1., 1.}}}},
      {"ss_KPi_invm", "#K#Pi Mass (GeV/c^2); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"os_KPi_pT_1", "#K#Pi pT (GeV/c); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"os_KPi_eTa_1", "#K#Pi eTa (GeV/c); Entries", {HistType::kTH1F, {{100, -1., 1.}}}},
      {"os_KPi_invm_1", "#K#Pi Mass (GeV/c^2); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"ss_KPi_pT_1", "#K#Pi pT (GeV/c); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"ss_KPi_eTa_1", "#K#Pi eTa (GeV/c); Entries", {HistType::kTH1F, {{100, -1., 1.}}}},
      {"ss_KPi_invm_1", "#K#Pi Mass (GeV/c^2); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"os_KPi_pT_0", "#K#Pi pT (GeV/c); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"os_KPi_eTa_0", "#K#Pi eTa (GeV/c); Entries", {HistType::kTH1F, {{100, -1., 1.}}}},
      {"os_KPi_invm_0", "#K#Pi Mass (GeV/c^2); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"ss_KPi_pT_0", "#K#Pi pT (GeV/c); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"ss_KPi_eTa_0", "#K#Pi eTa (GeV/c); Entries", {HistType::kTH1F, {{100, -1., 1.}}}},
      {"ss_KPi_invm_0", "#K#Pi Mass (GeV/c^2); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"os_KPi_pT_2", "#K#Pi pT (GeV/c); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"os_KPi_eTa_2", "#K#Pi eTa (GeV/c); Entries", {HistType::kTH1F, {{100, -1., 1.}}}},
      {"os_KPi_invm_2", "#K#Pi Mass (GeV/c^2); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"ss_KPi_pT_2", "#K#Pi pT (GeV/c); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"ss_KPi_eTa_2", "#K#Pi eTa (GeV/c); Entries", {HistType::kTH1F, {{100, -1., 1.}}}},
      {"ss_KPi_invm_2", "#K#Pi Mass (GeV/c^2); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
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
    // int truegapSide = sgSelector.trueGap(collision, *FIT_cut, ZDC_cut);
    int truegapSide = sgSelector.trueGap(collision, FIT_cut[0], FIT_cut[1], FIT_cut[3], ZDC_cut);
    registry.fill(HIST("GapSide"), gapSide);
    registry.fill(HIST("TrueGapSide"), truegapSide);
    gapSide = truegapSide;
    // Look for D0 and D0bar
    for (auto& [t0, t1] : combinations(tracks, tracks)) {
      // PID cut - t0=K, t1=pi
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
          }
          if (gapSide == 1) {
            registry.fill(HIST("os_KPi_pT_1"), v01.Pt());
            registry.fill(HIST("os_KPi_eTa_1"), v01.Eta());
            registry.fill(HIST("os_KPi_invm_1"), v01.M());
          }
          if (gapSide == 2) {
            registry.fill(HIST("os_KPi_pT_2"), v01.Pt());
            registry.fill(HIST("os_KPi_eTa_2"), v01.Eta());
            registry.fill(HIST("os_KPi_invm_2"), v01.M());
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
          }
          if (gapSide == 1) {
            registry.fill(HIST("os_KPi_pT_1"), v01.Pt());
            registry.fill(HIST("os_KPi_eTa_1"), v01.Eta());
            registry.fill(HIST("os_KPi_invm_1"), v01.M());
          }
          if (gapSide == 2) {
            registry.fill(HIST("os_KPi_pT_2"), v01.Pt());
            registry.fill(HIST("os_KPi_eTa_2"), v01.Eta());
            registry.fill(HIST("os_KPi_invm_2"), v01.M());
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
    adaptAnalysisTask<SG_D0_Analyzer>(cfgc)};
}
