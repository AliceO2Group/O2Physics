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
struct SGExclOmega {
  SGSelector sgSelector;
  Service<o2::framework::O2DatabasePDG> pdg;
  Configurable<float> FV0_cut{"FV0", 100., "FV0A threshold"};
  Configurable<float> FT0A_cut{"FT0A", 100., "FT0A threshold"};
  Configurable<float> FT0C_cut{"FT0C", 50., "FT0C threshold"};
  Configurable<float> FDDA_cut{"FDDA", 10000., "FDDA threshold"};
  Configurable<float> FDDC_cut{"FDDC", 10000., "FDDC threshold"};
  Configurable<float> ZDC_cut{"ZDC", 10., "ZDC threshold"};
  // Track Selections
  Configurable<float> PV_cut{"PV_cut", 1.0, "Use Only PV tracks"};
  Configurable<float> dcaZ_cut{"dcaZ_cut", 2.0, "dcaZ cut"};
  Configurable<float> dcaXY_cut{"dcaXY_cut", 2.0, "dcaXY cut (0 for Pt-function)"};
  Configurable<float> tpcChi2_cut{"tpcChi2_cut", 4, "Max tpcChi2NCl"};
  Configurable<float> tpcNClsFindable_cut{"tpcNClsFindable_cut", 70, "Min tpcNClsFindable"};
  Configurable<float> itsChi2_cut{"itsChi2_cut", 36, "Max itsChi2NCl"};
  Configurable<float> eta_cut{"eta_cut", 0.9, "Track Pseudorapidity"};
  Configurable<float> pt_cut{"pt_cut", 0.1, "Track Pt"};
  // D0 Specific Cuts
  Configurable<float> nsigmatpc_cut{"nsigmatpc", 3.0, "nsigma tpc cut"};
  Configurable<float> nsigmatof_cut{"nsigmatof", 9.0, "nsigma tof cut"};
  Configurable<bool> use_tof{"Use_TOF", true, "TOF PID"};
  HistogramRegistry registry{
    "registry",
    {

      {"GapSide", "Gap Side; Entries", {HistType::kTH1F, {{4, -1.5, 2.5}}}},
      {"TrueGapSide", "Gap Side; Entries", {HistType::kTH1F, {{4, -1.5, 2.5}}}},
      {"os_O_pT", "pT (GeV/c); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"os_O_eTa", "eTa (GeV/c); Entries", {HistType::kTH1F, {{100, -1., 1.}}}},
      {"os_O_invm", "Mass (GeV/c^2); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"ss_O_pT", "pT (GeV/c); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"ss_O_eTa", "eTa (GeV/c); Entries", {HistType::kTH1F, {{100, -1., 1.}}}},
      {"ss_O_invm", "Mass (GeV/c^2); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"pi0_invm", "Mass (GeV/c^2); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"os_O_pt_invm", " Nt vs Mass", {HistType::kTH2F, {{500, 0, 10.}, {500, 0, 2.5}}}},
      {"os_O_pT_pid", "pT (GeV/c); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"os_O_eTa_pid", "eTa (GeV/c); Entries", {HistType::kTH1F, {{100, -1., 1.}}}},
      {"os_O_invm_pid", "Mass (GeV/c^2); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"ss_O_pT_pid", "pT (GeV/c); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"ss_O_eTa_pid", "eTa (GeV/c); Entries", {HistType::kTH1F, {{100, -1., 1.}}}},
      {"ss_O_invm_pid", "Mass (GeV/c^2); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"os_O_pt_invm_pid", "pt vs Mass", {HistType::kTH2F, {{500, 0, 10.}, {500, 0, 2.5}}}},
      {"pi0_invm_pid", "Mass (GeV/c^2); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"os_O_pT_pid_pi0", "pT (GeV/c); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"os_O_eTa_pid_pi0", "eTa (GeV/c); Entries", {HistType::kTH1F, {{100, -1., 1.}}}},
      {"os_O_invm_pid_pi0", "Mass (GeV/c^2); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"ss_O_pT_pid_pi0", "pT (GeV/c); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"ss_O_eTa_pid_pi0", "eTa (GeV/c); Entries", {HistType::kTH1F, {{100, -1., 1.}}}},
      {"ss_O_invm_pid_pi0", "Mass (GeV/c^2); Entries", {HistType::kTH1F, {{5000, 0, 10}}}},
      {"os_O_pt_invm_pid_pi0", "pt vs Mass", {HistType::kTH2F, {{500, 0, 10.}, {500, 0, 2.5}}}},
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
    TLorentzVector v2;
    TLorentzVector v3;
    TLorentzVector v4;
    TLorentzVector v5;
    TLorentzVector v01;
    std::vector<TLorentzVector> els, pis;
    TLorentzVector v00;
    //  int truegapSide = sgSelector.trueGap(collision);
    // int truegapSide = sgSelector.trueGap(collision, FV0_cut, ZDC_cut);
    float FIT_cut[5] = {FV0_cut, FT0A_cut, FT0C_cut, FDDA_cut, FDDC_cut};
    std::vector<float> parameters = {PV_cut, dcaZ_cut, dcaXY_cut, tpcChi2_cut, tpcNClsFindable_cut, itsChi2_cut, eta_cut, pt_cut};
    // int truegapSide = sgSelector.trueGap(collision, *FIT_cut, ZDC_cut);
    int truegapSide = sgSelector.trueGap(collision, FIT_cut[0], FIT_cut[1], FIT_cut[2], ZDC_cut);
    registry.fill(HIST("GapSide"), gapSide);
    registry.fill(HIST("TrueGapSide"), truegapSide);
    gapSide = truegapSide;
    // if (gapSide!=2) return;
    int pvtracks = 0;
    int esign = 0;
    int nElec = 0;
    for (auto& t0 : tracks) {
      if (trackselector(t0, parameters) && t0.isPVContributor()) {
        pvtracks++;
        if (selectionPIDElec(t0, use_tof, nsigmatpc_cut, nsigmatof_cut)) {
          nElec++;
          esign += t0.sign();
          a.SetXYZM(t0.px(), t0.py(), t0.pz(), o2::constants::physics::MassElectron);
          els.push_back(a);
        } else {
          a.SetXYZM(t0.px(), t0.py(), t0.pz(), o2::constants::physics::MassPionCharged);
          pis.push_back(a);
        }
      }
    }
    // Look for D0 and D0bar
    // if (pvtracks != 6) return;
    if (nElec != 4 || esign != 0 || pvtracks < 6)
      return;
    // Apply pion hypothesis and create pairs
    v00 = els[0] + els[1] + els[2] + els[3];
    for (auto& [t0, t1] : combinations(tracks, tracks)) {
      if (selectionPIDElec(t0, use_tof, nsigmatpc_cut, nsigmatof_cut) || selectionPIDElec(t1, use_tof, nsigmatpc_cut, nsigmatof_cut))
        continue;
      v0.SetXYZM(t0.px(), t0.py(), t0.pz(), o2::constants::physics::MassPionCharged);
      v1.SetXYZM(t0.px(), t0.py(), t0.pz(), o2::constants::physics::MassPionCharged);
      v01 = v0 + v1 + v00;
      // Opposite sign pairs
      if (t0.sign() != t1.sign()) {
        registry.fill(HIST("os_O_pT"), v01.Pt());
        registry.fill(HIST("os_O_eTa"), v01.Eta());
        registry.fill(HIST("os_O_invm"), v01.M());
        registry.fill(HIST("os_O_pt_invm"), v01.Pt(), v01.M());
        registry.fill(HIST("pi0_invm"), v00.M());
      } else {
        registry.fill(HIST("ss_O_pT"), v01.Pt());
        registry.fill(HIST("ss_O_eTa"), v01.Eta());
        registry.fill(HIST("ss_O_invm"), v01.M());
      }
      if (selectionPIDPion(t0, use_tof, nsigmatpc_cut, nsigmatof_cut) && selectionPIDPion(t1, use_tof, nsigmatpc_cut, nsigmatof_cut)) {
        if (t0.sign() != t1.sign()) {
          registry.fill(HIST("os_O_pT_pid"), v01.Pt());
          registry.fill(HIST("os_O_eTa_pid"), v01.Eta());
          registry.fill(HIST("os_O_invm_pid"), v01.M());
          registry.fill(HIST("os_O_pt_invm_pid"), v01.Pt(), v01.M());
        } else {
          registry.fill(HIST("ss_O_pT_pid_pi0"), v01.Pt());
          registry.fill(HIST("ss_O_eTa_pid_pi0"), v01.Eta());
          registry.fill(HIST("ss_O_invm_pid_pi0"), v01.M());
        }

        if (abs(v00.M() - o2::constants::physics::MassPionNeutral) < 0.1) {
          if (t0.sign() != t1.sign()) {
            registry.fill(HIST("os_O_pT_pid_pi0"), v01.Pt());
            registry.fill(HIST("os_O_eTa_pid_pi0"), v01.Eta());
            registry.fill(HIST("os_O_invm_pid_pi0"), v01.M());
            registry.fill(HIST("os_O_pt_invm_pid_pi0"), v01.Pt(), v01.M());
          } else {
            registry.fill(HIST("ss_O_pT_pid_pi0"), v01.Pt());
            registry.fill(HIST("ss_O_eTa_pid_pi0"), v01.Eta());
            registry.fill(HIST("ss_O_invm_pid_pi0"), v01.M());
          }
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<SGExclOmega>(cfgc)};
}
