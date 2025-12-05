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
struct SGD0Analyzer {
  SGSelector sgSelector;
  Service<o2::framework::O2DatabasePDG> pdg;
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
  Configurable<float> pt_cut{"pt_cut", 0.1, "Track Pt"};
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
      {"os_Ntr_KPi_invm_2", "N tracks vs K#pi Mass", {HistType::kTH2F, {{48, 1.5, 49.5}, {750, .7, 2.2}}}},
      {"D0_FT0A", "FIT Amplitude", {HistType::kTH1F, {{500, 0, 1000}}}},
      {"D0_FT0C", "FIT Amplitude", {HistType::kTH1F, {{500, 0, 1000}}}},
      {"D0_FV0A", "FIT Amplitude", {HistType::kTH1F, {{500, 0, 1000}}}},
      {"D0_FDDA", "FIT Amplitude", {HistType::kTH1F, {{500, 0, 1000}}}},
      {"D0_FDDC", "FIT Amplitude", {HistType::kTH1F, {{500, 0, 1000}}}},
      {"D0_ZNA", "ZDC Energy", {HistType::kTH1F, {{500, 0, 200}}}},
      {"D0_ZNC", "ZDC Energy", {HistType::kTH1F, {{500, 0, 200}}}},
      {"D0_Ntr_0", "Tracks", {HistType::kTH1F, {{100, 0.5, 100.5}}}},
      {"D0_Ntr_1", "Tracks", {HistType::kTH1F, {{100, 0.5, 100.5}}}},
      {"D0_Ntr_2", "Tracks", {HistType::kTH1F, {{100, 0.5, 100.5}}}},
      {"D0_NtrPV_0", "Tracks", {HistType::kTH1F, {{50, 0.5, 50.5}}}},
      {"D0_NtrPV_1", "Tracks", {HistType::kTH1F, {{50, 0.5, 50.5}}}},
      {"D0_NtrPV_2", "Tracks", {HistType::kTH1F, {{50, 0.5, 50.5}}}}}};
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
            if (v01.M() > 1.8 && v01.M() < 1.9) {
              registry.fill(HIST("D0_Ntr_0"), tracks.size());
              registry.fill(HIST("D0_NtrPV_0"), pvtracks);
            }
          }
          if (gapSide == 1) {
            registry.fill(HIST("os_KPi_pT_1"), v01.Pt());
            registry.fill(HIST("os_KPi_eTa_1"), v01.Eta());
            registry.fill(HIST("os_KPi_invm_1"), v01.M());
            registry.fill(HIST("os_Ntr_KPi_invm_1"), tracks.size(), v01.M());
            if (v01.M() > 1.8 && v01.M() < 1.9) {
              registry.fill(HIST("D0_Ntr_1"), tracks.size());
              registry.fill(HIST("D0_NtrPV_1"), pvtracks);
            }
          }
          if (gapSide == 2) {
            registry.fill(HIST("os_KPi_pT_2"), v01.Pt());
            registry.fill(HIST("os_KPi_eTa_2"), v01.Eta());
            registry.fill(HIST("os_KPi_invm_2"), v01.M());
            registry.fill(HIST("os_Ntr_KPi_invm_2"), tracks.size(), v01.M());
            if (v01.M() > 1.8 && v01.M() < 1.9) {
              registry.fill(HIST("D0_Ntr_2"), tracks.size());
              registry.fill(HIST("D0_NtrPV_2"), pvtracks);
              registry.fill(HIST("D0_FT0A"), collision.totalFT0AmplitudeA());
              registry.fill(HIST("D0_FT0C"), collision.totalFT0AmplitudeC());
              registry.fill(HIST("D0_FV0A"), collision.totalFV0AmplitudeA());
              registry.fill(HIST("D0_FDDA"), collision.totalFDDAmplitudeA());
              registry.fill(HIST("D0_FDDC"), collision.totalFDDAmplitudeC());
              registry.fill(HIST("D0_ZNA"), collision.energyCommonZNA());
              registry.fill(HIST("D0_ZNC"), collision.energyCommonZNC());
            }
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
            if (v01.M() > 1.8 && v01.M() < 1.9) {
              registry.fill(HIST("D0_Ntr_0"), tracks.size());
              registry.fill(HIST("D0_NtrPV_0"), pvtracks);
            }
          }
          if (gapSide == 1) {
            registry.fill(HIST("os_KPi_pT_1"), v01.Pt());
            registry.fill(HIST("os_KPi_eTa_1"), v01.Eta());
            registry.fill(HIST("os_KPi_invm_1"), v01.M());
            registry.fill(HIST("os_Ntr_KPi_invm_1"), tracks.size(), v01.M());
            if (v01.M() > 1.8 && v01.M() < 1.9) {
              registry.fill(HIST("D0_Ntr_1"), tracks.size());
              registry.fill(HIST("D0_NtrPV_1"), pvtracks);
            }
          }
          if (gapSide == 2) {
            registry.fill(HIST("os_KPi_pT_2"), v01.Pt());
            registry.fill(HIST("os_KPi_eTa_2"), v01.Eta());
            registry.fill(HIST("os_KPi_invm_2"), v01.M());
            registry.fill(HIST("os_Ntr_KPi_invm_2"), tracks.size(), v01.M());
            if (v01.M() > 1.8 && v01.M() < 1.9) {
              registry.fill(HIST("D0_Ntr_2"), tracks.size());
              registry.fill(HIST("D0_NtrPV_2"), pvtracks);
              registry.fill(HIST("D0_FT0A"), collision.totalFT0AmplitudeA());
              registry.fill(HIST("D0_FT0C"), collision.totalFT0AmplitudeC());
              registry.fill(HIST("D0_FV0A"), collision.totalFV0AmplitudeA());
              registry.fill(HIST("D0_FDDA"), collision.totalFDDAmplitudeA());
              registry.fill(HIST("D0_FDDC"), collision.totalFDDAmplitudeC());
              registry.fill(HIST("D0_ZNA"), collision.energyCommonZNA());
              registry.fill(HIST("D0_ZNC"), collision.energyCommonZNC());
            }
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
