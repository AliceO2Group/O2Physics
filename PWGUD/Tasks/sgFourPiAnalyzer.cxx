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

#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/DataTypes.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
#include "MathUtils/Utils.h"

#include "TLorentzVector.h"
#include <TString.h>

#include <iostream>
using namespace std;
using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
struct SGFourPiAnalyzer {
  SGSelector sgSelector;
  Service<o2::framework::O2DatabasePDG> pdg;
  // Adjusted Gap thresholds
  Configurable<float> FV0_cut{"FV0", 50., "FV0A threshold"};
  Configurable<float> FT0A_cut{"FT0A", 150., "FT0A threshold"};
  Configurable<float> FT0C_cut{"FT0C", 50., "FT0C threshold"};
  Configurable<float> FDDA_cut{"FDDA", 10000., "FDDA threshold"};
  Configurable<float> FDDC_cut{"FDDC", 10000., "FDDC threshold"};
  Configurable<float> ZDC_cut{"ZDC", 10., "ZDC threshold"};
  // Track Selections
  Configurable<float> PV_cut{"PV_cut", 1.0, "Use Only PV tracks"};
  Configurable<float> dcaZ_cut{"dcaZ_cut", 2.0, "dcaZ cut"};
  Configurable<float> dcaXY_cut{"dcaXY_cut", 0.0, "dcaXY cut (0 for Pt-function)"};
  Configurable<float> tpcChi2_cut{"tpcChi2_cut", 4, "Max tpcChi2NCl"};
  Configurable<float> tpcNClsFindable_cut{"tpcNClsFindable_cut", 70, "Min tpcNClsFindable"};
  Configurable<float> itsChi2_cut{"itsChi2_cut", 36, "Max itsChi2NCl"};
  Configurable<float> eta_cut{"eta_cut", 0.9, "Track Pseudorapidity"};
  Configurable<float> pt_cut{"pt_cut", 0.1, "Track Pt"};
  HistogramRegistry registry{
    "registry",
    {

      {"GapSide", "Gap Side; Entries", {HistType::kTH1F, {{4, -1.5, 2.5}}}},
      {"TrueGapSide", "Gap Side; Entries", {HistType::kTH1F, {{4, -1.5, 2.5}}}},
      {"ITSNCls", "ITS Clusters", {HistType::kTH1F, {{10, -.5, 9.5}}}},
      {"os_4Pi_pT", "#K#Pi pT (GeV/c); Entries", {HistType::kTH1F, {{1000, 0, 10}}}},
      {"os_4Pi_eTa", "#K#Pi eTa (GeV/c); Entries", {HistType::kTH1F, {{100, -2., 2.}}}},
      {"os_4Pi_invm", "#K#Pi Mass (GeV/c^2); Entries", {HistType::kTH1F, {{1000, 0, 10}}}},
      {"ss_4Pi_pT", "#K#Pi pT (GeV/c); Entries", {HistType::kTH1F, {{1000, 0, 10}}}},
      {"ss_4Pi_eTa", "#K#Pi eTa (GeV/c); Entries", {HistType::kTH1F, {{100, -2., 2.}}}},
      {"ss_4Pi_invm", "#K#Pi Mass (GeV/c^2); Entries", {HistType::kTH1F, {{1000, 0, 10}}}},
      {"os_4Pi_pT_1", "#K#Pi pT (GeV/c); Entries", {HistType::kTH1F, {{1000, 0, 10}}}},
      {"os_4Pi_eTa_1", "#K#Pi eTa (GeV/c); Entries", {HistType::kTH1F, {{100, -2., 2.}}}},
      {"os_4Pi_invm_1", "#K#Pi Mass (GeV/c^2); Entries", {HistType::kTH1F, {{1000, 0, 10}}}},
      {"ss_4Pi_pT_1", "#K#Pi pT (GeV/c); Entries", {HistType::kTH1F, {{1000, 0, 10}}}},
      {"ss_4Pi_eTa_1", "#K#Pi eTa (GeV/c); Entries", {HistType::kTH1F, {{100, -2., 2.}}}},
      {"ss_4Pi_invm_1", "#K#Pi Mass (GeV/c^2); Entries", {HistType::kTH1F, {{1000, 0, 10}}}},
      {"os_4Pi_pT_0", "#K#Pi pT (GeV/c); Entries", {HistType::kTH1F, {{1000, 0, 10}}}},
      {"os_4Pi_eTa_0", "#K#Pi eTa (GeV/c); Entries", {HistType::kTH1F, {{100, -2., 2.}}}},
      {"os_4Pi_invm_0", "#K#Pi Mass (GeV/c^2); Entries", {HistType::kTH1F, {{1000, 0, 10}}}},
      {"ss_4Pi_pT_0", "#K#Pi pT (GeV/c); Entries", {HistType::kTH1F, {{1000, 0, 10}}}},
      {"ss_4Pi_eTa_0", "#K#Pi eTa (GeV/c); Entries", {HistType::kTH1F, {{100, -2., 2.}}}},
      {"ss_4Pi_invm_0", "#K#Pi Mass (GeV/c^2); Entries", {HistType::kTH1F, {{1000, 0, 10}}}},
      {"os_4Pi_pT_2", "#K#Pi pT (GeV/c); Entries", {HistType::kTH1F, {{1000, 0, 10}}}},
      {"os_4Pi_eTa_2", "#K#Pi eTa (GeV/c); Entries", {HistType::kTH1F, {{100, -2., 2.}}}},
      {"os_4Pi_invm_2", "#K#Pi Mass (GeV/c^2); Entries", {HistType::kTH1F, {{1000, 0, 10}}}},
      {"ss_4Pi_pT_2", "#K#Pi pT (GeV/c); Entries", {HistType::kTH1F, {{1000, 0, 10}}}},
      {"ss_4Pi_eTa_2", "#K#Pi eTa (GeV/c); Entries", {HistType::kTH1F, {{100, -2., 2.}}}},
      {"ss_4Pi_invm_2", "#K#Pi Mass (GeV/c^2); Entries", {HistType::kTH1F, {{1000, 0, 10}}}},
    }};
  using udtracks = soa::Join<aod::UDTracks, aod::UDTracksExtra, aod::UDTracksPID>;
  using udtracksfull = soa::Join<aod::UDTracks, aod::UDTracksPID, aod::UDTracksExtra, aod::UDTracksFlags, aod::UDTracksDCA>;
  using UDCollisionsFull = soa::Join<aod::UDCollisions, aod::SGCollisions, aod::UDCollisionsSels, aod::UDZdcsReduced>; //
  using UDCollisionFull = UDCollisionsFull::iterator;

  void process(UDCollisionFull const& collision, udtracksfull const& tracks)
  {
    const float mpion = pdg->Mass(211);
    TLorentzVector a;
    int gapSide = collision.gapSide();
    if (gapSide < 0 || gapSide > 2)
      return;
    // Single gap either side
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
    std::vector<TLorentzVector> goodTracks;
    float sign = 0;
    for (auto t : tracks) {
      int itsNCls = t.itsNCls();
      // if (itsNCls) {
      registry.fill(HIST("ITSNCls"), itsNCls);
      //}
      TLorentzVector a;
      a.SetXYZM(t.px(), t.py(), t.pz(), mpion);
      if (trackselector(t, parameters)) {
        sign += t.sign();
        goodTracks.push_back(a);
      }
    }
    //    std::cout << goodTracks.size()<<std::endl;
    if (goodTracks.size() == 4) {
      for (auto pion : goodTracks) {
        v01 += pion;
      }
      // Apply pion hypothesis and create pairs
      // Opposite sign pairs

      if (sign == 0) {
        registry.fill(HIST("os_4Pi_pT"), v01.Pt());
        registry.fill(HIST("os_4Pi_eTa"), v01.Eta());
        if (TMath::Abs(v01.Eta() < 0.9) && v01.Pt() < .15)
          registry.fill(HIST("os_4Pi_invm"), v01.M());
        if (gapSide == 0) {
          registry.fill(HIST("os_4Pi_pT_0"), v01.Pt());
          registry.fill(HIST("os_4Pi_eTa_0"), v01.Eta());
          if (TMath::Abs(v01.Eta() < 0.9) && v01.Pt() < .15)
            registry.fill(HIST("os_4Pi_invm_0"), v01.M());
        }
        if (gapSide == 1) {
          registry.fill(HIST("os_4Pi_pT_1"), v01.Pt());
          registry.fill(HIST("os_4Pi_eTa_1"), v01.Eta());
          if (TMath::Abs(v01.Eta() < 0.9) && v01.Pt() < .15)
            registry.fill(HIST("os_4Pi_invm_1"), v01.M());
        }
        if (gapSide == 2) {
          registry.fill(HIST("os_4Pi_pT_2"), v01.Pt());
          registry.fill(HIST("os_4Pi_eTa_2"), v01.Eta());
          if (TMath::Abs(v01.Eta() < 0.9) && v01.Pt() < .15)
            registry.fill(HIST("os_4Pi_invm_2"), v01.M());
        }
      } else {
        registry.fill(HIST("ss_4Pi_pT"), v01.Pt());
        registry.fill(HIST("ss_4Pi_eTa"), v01.Eta());
        if (TMath::Abs(v01.Eta() < 0.9) && v01.Pt() < .15)
          registry.fill(HIST("ss_4Pi_invm"), v01.M());
        if (gapSide == 0) {
          registry.fill(HIST("ss_4Pi_pT_0"), v01.Pt());
          registry.fill(HIST("ss_4Pi_eTa_0"), v01.Eta());
          if (TMath::Abs(v01.Eta() < 0.9) && v01.Pt() < .15)
            registry.fill(HIST("ss_4Pi_invm_0"), v01.M());
        }
        if (gapSide == 1) {
          registry.fill(HIST("ss_4Pi_pT_1"), v01.Pt());
          registry.fill(HIST("ss_4Pi_eTa_1"), v01.Eta());
          if (TMath::Abs(v01.Eta() < 0.9) && v01.Pt() < .15)
            registry.fill(HIST("ss_4Pi_invm_1"), v01.M());
        }
        if (gapSide == 2) {
          registry.fill(HIST("ss_4Pi_pT_2"), v01.Pt());
          registry.fill(HIST("ss_4Pi_eTa_2"), v01.Eta());
          if (TMath::Abs(v01.Eta() < 0.9) && v01.Pt() < .15)
            registry.fill(HIST("ss_4Pi_invm_2"), v01.M());
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<SGFourPiAnalyzer>(cfgc)};
}
