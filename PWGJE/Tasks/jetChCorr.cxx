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

// jet analysis tasks (subscribing to jet finder task)
//
/// Mriganka Mouli Mondal <mriganka.mouli.mondal@cern.ch>    originally modified from  Nima Zardoshti <nima.zardoshti@cern.ch>
//

#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"

#include <TMath.h>
#include <TMathBase.h>
#include <TVector3.h>

#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/PseudoJet.hh"
#include <fastjet/JetDefinition.hh>

#include <cmath>
#include <utility>
#include <vector>

// #include "PWGLF/DataModel/LFResonanceTables.h"

#include <Framework/AnalysisHelpers.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

using namespace std;
using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

int trackSelection = -1;
int trackL = -1;
int n_trackL = -1;

TVector3 v1, v2, vR; // leading and next to leading vector
int ch_l = 0;        // charge leading
int ch_nl = 0;       // charge next to leading

int ch_mult = 0; // charge next to leading

struct JetChCorr {
  // Produces<aod::CJetSSs> jetSubstructureDataTable;
  // Produces<aod::CMCDJetSSs> jetSubstructureMCDTable;
  // Produces<aod::CMCPJetSSs> jetSubstructureMCPTable;
  // Produces<aod::CEWSJetSSs> jetSubstructureDataSubTable;

  Configurable<float> zCut{"zCut", 0.05, "soft drop z cut"};
  Configurable<float> beta{"beta", 0.1, "soft drop beta"};
  Configurable<float> jetPtTh{"jetPtTh", 15.0, "jet transverse momentum cut"};

  Service<o2::framework::O2DatabasePDG> pdg;
  std::vector<fastjet::PseudoJet> jetConstituents;
  std::vector<fastjet::PseudoJet> jetReclustered;
  JetFinder jetReclusterer;

  HistogramRegistry registry;

  void init(InitContext const&)
  {
    registry.add("h2_jet_pt_jet_zg", ";#it{p}_{T,jet} (GeV/#it{c});#it{z}_{g}", {HistType::kTH2F, {{200, 0., 200.}, {22, 0.0, 1.1}}});
    registry.add("h2_jet_pt_jet_rg", ";#it{p}_{T,jet} (GeV/#it{c});#it{R}_{g}", {HistType::kTH2F, {{200, 0., 200.}, {22, 0.0, 1.1}}});
    registry.add("h2_jet_pt_jet_nsd", ";#it{p}_{T,jet} (GeV/#it{c});#it{n}_{SD}", {HistType::kTH2F, {{200, 0., 200.}, {15, -0.5, 14.5}}});

    registry.add("h2_jet_pt_part_jet_zg_part", ";#it{p}_{T,jet}^{part} (GeV/#it{c});#it{z}_{g}^{part}", {HistType::kTH2F, {{200, 0., 200.}, {22, 0.0, 1.1}}});
    registry.add("h2_jet_pt_part_jet_rg_part", ";#it{p}_{T,jet}^{part} (GeV/#it{c});#it{R}_{g}^{part}", {HistType::kTH2F, {{200, 0., 200.}, {22, 0.0, 1.1}}});
    registry.add("h2_jet_pt_part_jet_nsd_part", ";#it{p}_{T,jet}^{part} (GeV/#it{c});#it{n}_{SD}^{part}", {HistType::kTH2F, {{200, 0., 200.}, {15, -0.5, 14.5}}});

    registry.add("h2_jet_pt_jet_zg_eventwiseconstituentsubtracted", ";#it{p}_{T,jet} (GeV/#it{c});#it{z}_{g}", {HistType::kTH2F, {{200, 0., 200.}, {22, 0.0, 1.1}}});
    registry.add("h2_jet_pt_jet_rg_eventwiseconstituentsubtracted", ";#it{p}_{T,jet} (GeV/#it{c});#it{R}_{g}", {HistType::kTH2F, {{200, 0., 200.}, {22, 0.0, 1.1}}});
    registry.add("h2_jet_pt_jet_nsd_eventwiseconstituentsubtracted", ";#it{p}_{T,jet} (GeV/#it{c});#it{n}_{SD}", {HistType::kTH2F, {{200, 0., 200.}, {15, -0.5, 14.5}}});

    ////  For L and NL particles

    ////  For L and NL particles
    registry.add("h_ch_s_pt", ";#it{p}_{T,pair} (GeV/#it{c});N_{ch} (hh)", {HistType::kTH1F, {{50, 0., 50.}}});
    registry.add("h_ch_d_pt", ";#it{p}_{T,pair} (GeV/#it{c});N_{ch} (h#bar{h})", {HistType::kTH1F, {{50, 0., 50.}}});
    registry.add("h_ch_s_kt", ";#it{k}_{T,pair} (GeV/#it{c});N_{ch} (hh)", {HistType::kTH1F, {{50, 0., 5.}}});
    registry.add("h_ch_d_kt", ";#it{k}_{T,pair} (GeV/#it{c});N_{ch} (h#bar{h})", {HistType::kTH1F, {{50, 0., 5.}}});
    registry.add("h_ch_s_z", ";#it{z}_{pair} ;N_{ch} (hh)", {HistType::kTH1F, {{20, 0., .5}}});
    registry.add("h_ch_d_z", ";#it{z}_{pair} ;ch_{h#bar{h}}", {HistType::kTH1F, {{20, 0., .5}}});
    registry.add("h_ch_s_ft", ";#it{fm}_{time} (GeV/#it{c});N_{ch} (hh)", {HistType::kTH1F, {{20, -1., 4.}}});
    registry.add("h_ch_d_ft", ";#it{fm}_{time} (GeV/#it{c});N_{ch} (h#bar{h})", {HistType::kTH1F, {{20, -1., 4.}}});
    registry.add("h_ch_s_mult", ";#it{N}_{ch,jet} (GeV/#it{c});N_{ch} (hh)", {HistType::kTH1F, {{15, 0., 30.}}});
    registry.add("h_ch_d_mult", ";#it{N}_{ch,jet} (GeV/#it{c});N_{ch} (h#bar{h})", {HistType::kTH1F, {{15, 0., 30.}}});

    registry.add("hr1_ch_s_pt", ";#it{p}_{T,pair} (GeV/#it{c}) {nsd1};N_{ch} (hh)", {HistType::kTH1F, {{50, 0., 50.}}});
    registry.add("hr1_ch_d_pt", ";#it{p}_{T,pair} (GeV/#it{c}) {nsd1};N_{ch} (h#bar{h})", {HistType::kTH1F, {{50, 0., 50.}}});
    registry.add("hr1_ch_s_kt", ";#it{k}_{T,pair} (GeV/#it{c}) {nsd1};N_{ch} (hh)", {HistType::kTH1F, {{50, 0., 5.}}});
    registry.add("hr1_ch_d_kt", ";#it{k}_{T,pair} (GeV/#it{c}) {nsd1};N_{ch} (h#bar{h})", {HistType::kTH1F, {{50, 0., 5.}}});
    registry.add("hr1_ch_s_z", ";#it{z}_{pair} ;N_{ch} (hh) {nsd1}", {HistType::kTH1F, {{20, 0., .5}}});
    registry.add("hr1_ch_d_z", ";#it{z}_{pair} ;ch_{h#bar{h}} {nsd1}", {HistType::kTH1F, {{20, 0., .5}}});
    registry.add("hr1_ch_s_ft", ";#it{fm}_{time} (GeV/#it{c}) {nsd1};N_{ch} (hh)", {HistType::kTH1F, {{20, -1., 4.}}});
    registry.add("hr1_ch_d_ft", ";#it{fm}_{time} (GeV/#it{c}) {nsd1};N_{ch} (h#bar{h})", {HistType::kTH1F, {{20, -1., 4.}}});
    registry.add("hr1_ch_s_mult", ";#it{N}_{ch,jet} (GeV/#it{c}) {nsd1};N_{ch} (hh)", {HistType::kTH1F, {{15, 0., 30.}}});
    registry.add("hr1_ch_d_mult", ";#it{N}_{ch,jet} (GeV/#it{c}) {nsd1};N_{ch} (h#bar{h})", {HistType::kTH1F, {{15, 0., 30.}}});
    registry.add("hr1_ch_s_zg", ";#it{z}_{g} (GeV/#it{c}) {nsd1};N_{ch} (hh)", {HistType::kTH1F, {{20, 0., .5}}});
    registry.add("hr1_ch_d_zg", ";#it{z}_{g} (GeV/#it{c}) {nsd1};N_{ch} (h#bar{h})", {HistType::kTH1F, {{20, 0., .5}}});
    registry.add("hr1_ch_s_rg", ";#it{r}_{g} (GeV/#it{c}) {nsd1};N_{ch} (hh)", {HistType::kTH1F, {{20, -2., 0.}}});
    registry.add("hr1_ch_d_rg", ";#it{r}_{g} (GeV/#it{c}) {nsd1};N_{ch} (h#bar{h})", {HistType::kTH1F, {{20, -2., 0.}}});
    registry.add("hr1_ch_s_lu", ";#it{z}_{g} {nsd1};#it{r}_{g}", {HistType::kTH2F, {{20, 0., 5.}, {28, -2., 5.}}});
    registry.add("hr1_ch_d_lu", ";#it{z}_{g} {nsd1};#it{r}_{g}", {HistType::kTH2F, {{20, 0., 5.}, {28, -2., 5.}}});

    registry.add("hr2_ch_s_pt", ";#it{p}_{T,pair} (GeV/#it{c}) {nsd2};N_{ch} (hh)", {HistType::kTH1F, {{50, 0., 50.}}});
    registry.add("hr2_ch_d_pt", ";#it{p}_{T,pair} (GeV/#it{c}) {nsd2};N_{ch} (h#bar{h})", {HistType::kTH1F, {{50, 0., 50.}}});
    registry.add("hr2_ch_s_kt", ";#it{k}_{T,pair} (GeV/#it{c}) {nsd2};N_{ch} (hh)", {HistType::kTH1F, {{50, 0., 5.}}});
    registry.add("hr2_ch_d_kt", ";#it{k}_{T,pair} (GeV/#it{c}) {nsd2};N_{ch} (h#bar{h})", {HistType::kTH1F, {{50, 0., 5.}}});
    registry.add("hr2_ch_s_z", ";#it{z}_{pair} ;N_{ch} (hh) {nsd2}", {HistType::kTH1F, {{20, 0., .5}}});
    registry.add("hr2_ch_d_z", ";#it{z}_{pair} ;ch_{h#bar{h}} {nsd2}", {HistType::kTH1F, {{20, 0., .5}}});
    registry.add("hr2_ch_s_ft", ";#it{fm}_{time} (GeV/#it{c}) {nsd2};N_{ch} (hh)", {HistType::kTH1F, {{20, -1., 4.}}});
    registry.add("hr2_ch_d_ft", ";#it{fm}_{time} (GeV/#it{c}) {nsd2};N_{ch} (h#bar{h})", {HistType::kTH1F, {{20, -1., 4.}}});
    registry.add("hr2_ch_s_mult", ";#it{N}_{ch,jet} (GeV/#it{c}) {nsd2};N_{ch} (hh)", {HistType::kTH1F, {{15, 0., 30.}}});
    registry.add("hr2_ch_d_mult", ";#it{N}_{ch,jet} (GeV/#it{c}) {nsd2};N_{ch} (h#bar{h})", {HistType::kTH1F, {{15, 0., 30.}}});
    registry.add("hr2_ch_s_zg", ";#it{z}_{g} (GeV/#it{c}) {nsd2};N_{ch} (hh)", {HistType::kTH1F, {{20, 0., .5}}});
    registry.add("hr2_ch_d_zg", ";#it{z}_{g} (GeV/#it{c}) {nsd2};N_{ch} (h#bar{h})", {HistType::kTH1F, {{20, 0., .5}}});
    registry.add("hr2_ch_s_rg", ";#it{r}_{g} (GeV/#it{c}) {nsd2};N_{ch} (hh)", {HistType::kTH1F, {{20, -2., 0.}}});
    registry.add("hr2_ch_d_rg", ";#it{r}_{g} (GeV/#it{c}) {nsd2};N_{ch} (h#bar{h})", {HistType::kTH1F, {{20, -2., 0.}}});
    registry.add("hr2_ch_s_lu", ";#it{z}_{g} {nsd2};#it{r}_{g}", {HistType::kTH2F, {{20, 0., 5.}, {28, -2., 5.}}});
    registry.add("hr2_ch_d_lu", ";#it{z}_{g} {nsd2};#it{r}_{g}", {HistType::kTH2F, {{20, 0., 5.}, {28, -2., 5.}}});

    registry.add("hr3_ch_s_pt", ";#it{p}_{T,pair} (GeV/#it{c}) {nsd3};N_{ch} (hh)", {HistType::kTH1F, {{50, 0., 50.}}});
    registry.add("hr3_ch_d_pt", ";#it{p}_{T,pair} (GeV/#it{c}) {nsd3};N_{ch} (h#bar{h})", {HistType::kTH1F, {{50, 0., 50.}}});
    registry.add("hr3_ch_s_kt", ";#it{k}_{T,pair} (GeV/#it{c}) {nsd3};N_{ch} (hh)", {HistType::kTH1F, {{50, 0., 5.}}});
    registry.add("hr3_ch_d_kt", ";#it{k}_{T,pair} (GeV/#it{c}) {nsd3};N_{ch} (h#bar{h})", {HistType::kTH1F, {{50, 0., 5.}}});
    registry.add("hr3_ch_s_z", ";#it{z}_{pair} ;N_{ch} (hh) {nsd3}", {HistType::kTH1F, {{20, 0., .5}}});
    registry.add("hr3_ch_d_z", ";#it{z}_{pair} ;ch_{h#bar{h}} {nsd3}", {HistType::kTH1F, {{20, 0., .5}}});
    registry.add("hr3_ch_s_ft", ";#it{fm}_{time} (GeV/#it{c}) {nsd3};N_{ch} (hh)", {HistType::kTH1F, {{20, -1., 4.}}});
    registry.add("hr3_ch_d_ft", ";#it{fm}_{time} (GeV/#it{c}) {nsd3};N_{ch} (h#bar{h})", {HistType::kTH1F, {{20, -1., 4.}}});
    registry.add("hr3_ch_s_mult", ";#it{N}_{ch,jet} (GeV/#it{c}) {nsd3};N_{ch} (hh)", {HistType::kTH1F, {{15, 0., 30.}}});
    registry.add("hr3_ch_d_mult", ";#it{N}_{ch,jet} (GeV/#it{c}) {nsd3};N_{ch} (h#bar{h})", {HistType::kTH1F, {{15, 0., 30.}}});
    registry.add("hr3_ch_s_zg", ";#it{z}_{g} (GeV/#it{c}) {nsd3};N_{ch} (hh)", {HistType::kTH1F, {{20, 0., .5}}});
    registry.add("hr3_ch_d_zg", ";#it{z}_{g} (GeV/#it{c}) {nsd3};N_{ch} (h#bar{h})", {HistType::kTH1F, {{20, 0., .5}}});
    registry.add("hr3_ch_s_rg", ";#it{r}_{g} (GeV/#it{c}) {nsd3};N_{ch} (hh)", {HistType::kTH1F, {{20, -2., 0.}}});
    registry.add("hr3_ch_d_rg", ";#it{r}_{g} (GeV/#it{c}) {nsd3};N_{ch} (h#bar{h})", {HistType::kTH1F, {{20, -2., 0.}}});
    registry.add("hr3_ch_s_lu", ";#it{z}_{g} {nsd3};#it{r}_{g}", {HistType::kTH2F, {{20, 0., 5.}, {28, -2., 5.}}});
    registry.add("hr3_ch_d_lu", ";#it{z}_{g} {nsd3};#it{r}_{g}", {HistType::kTH2F, {{20, 0., 5.}, {28, -2., 5.}}});

    jetReclusterer.isReclustering = true;
    jetReclusterer.algorithm = fastjet::JetAlgorithm::cambridge_algorithm;
  }

  // template <typename JetTable, typename JetTableMCP, typename SubstructureTable>
  //  template <bool isMCP, bool isSubtracted, typename T, typename U>
  //  void jetReclustering(T const& jet, U& outputTable)
  template <bool isMCP, bool isSubtracted, typename T>
  void jetReclustering(T const& jet)
  {
    jetReclustered.clear();
    fastjet::ClusterSequenceArea clusterSeq(jetReclusterer.findJets(jetConstituents, jetReclustered));
    jetReclustered = sorted_by_pt(jetReclustered);
    fastjet::PseudoJet daughterSubJet = jetReclustered[0];
    fastjet::PseudoJet parentSubJet1;
    fastjet::PseudoJet parentSubJet2;
    bool softDropped = false;
    auto nsd = 0.0;

    // std::vector<float> energyMotherVec;
    // std::vector<float> ptLeadingVec;
    // std::vector<float> ptSubLeadingVec;
    // std::vector<float> thetaVec;

    float ptJet = daughterSubJet.pt();

    while (daughterSubJet.has_parents(parentSubJet1, parentSubJet2)) {
      if (parentSubJet1.perp() < parentSubJet2.perp()) {
        std::swap(parentSubJet1, parentSubJet2);
      }

      vector<fastjet::PseudoJet> constituents1 = parentSubJet1.constituents();
      vector<fastjet::PseudoJet> constituents2 = parentSubJet2.constituents();
      // cout<<"in subJET1 ********************************************* "<<endl;
      bool found1 = false;
      bool found2 = false;

      for (unsigned int j = 0; j < constituents1.size(); j++) {
        // cout<<constituents1[j].user_index()-1<<", ";
        if ((n_trackL == (constituents1[j].user_index() - 1)) || (trackL == (constituents1[j].user_index() - 1)))
          found1 = true;
      }
      // cout<<endl;
      // cout<<"in subJET2 ********************************************* "<<endl;
      for (unsigned int j = 0; j < constituents2.size(); j++) {
        // cout<<constituents2[j].user_index()-1<<", ";
        if ((n_trackL == constituents2[j].user_index() - 1) || (trackL == constituents2[j].user_index() - 1))
          found2 = true;
      }
      // cout<<endl;

      auto z = parentSubJet2.perp() / (parentSubJet1.perp() + parentSubJet2.perp());
      auto theta = parentSubJet1.delta_R(parentSubJet2);
      // energyMotherVec.push_back(daughterSubJet.e());
      // ptLeadingVec.push_back(parentSubJet1.pt());
      // ptSubLeadingVec.push_back(parentSubJet2.pt());
      // thetaVec.push_back(theta);

      if (z >= zCut * TMath::Power(theta / (jet.r() / 100.f), beta)) {
        if (found1 == true && found2 == true) { // found leading and next-to-leading in seperate prongs
          if (!softDropped) {
            auto zg = z;
            auto rg = theta;
            if constexpr (!isSubtracted && !isMCP) {
              registry.fill(HIST("h2_jet_pt_jet_zg"), jet.pt(), zg);
              registry.fill(HIST("h2_jet_pt_jet_rg"), jet.pt(), rg);
            }
            if constexpr (!isSubtracted && isMCP) {
              registry.fill(HIST("h2_jet_pt_part_jet_zg_part"), jet.pt(), zg);
              registry.fill(HIST("h2_jet_pt_part_jet_rg_part"), jet.pt(), rg);
            }
            if constexpr (isSubtracted && !isMCP) {
              registry.fill(HIST("h2_jet_pt_jet_zg_eventwiseconstituentsubtracted"), jet.pt(), zg);
              registry.fill(HIST("h2_jet_pt_jet_rg_eventwiseconstituentsubtracted"), jet.pt(), rg);
            }
            softDropped = true;

            v1.SetXYZ(parentSubJet1.px(), parentSubJet1.py(), parentSubJet1.pz());
            v2.SetXYZ(parentSubJet2.px(), parentSubJet2.py(), parentSubJet2.pz());

            vR = v1 + v2;
            float z = v2.Perp(vR.Orthogonal()) / (v1.Perp(vR.Orthogonal()) + v2.Perp(vR.Orthogonal())); // you redefine z here, please check!
            float fT = ((2. * z * (1 - z) * vR.Mag()) / v1.Perp2(vR)) / 6.;
            float kt_p = v1.Perp(vR);

            ////////  for first matched split
            if (nsd == 0) {
              if (ch_l == ch_nl) {
                registry.fill(HIST("hr1_ch_s_pt"), ptJet);
                registry.fill(HIST("hr1_ch_s_kt"), kt_p);
                registry.fill(HIST("hr1_ch_s_z"), z);
                registry.fill(HIST("hr1_ch_s_ft"), log10(fT));
                registry.fill(HIST("hr1_ch_s_mult"), ch_mult);
                registry.fill(HIST("hr1_ch_s_zg"), zg);
                registry.fill(HIST("hr1_ch_s_rg"), log(rg));
                registry.fill(HIST("hr1_ch_s_lu"), log(1 / rg), log(kt_p));
              }
              if (ch_l != ch_nl) {
                registry.fill(HIST("hr1_ch_d_pt"), ptJet);
                registry.fill(HIST("hr1_ch_d_kt"), kt_p);
                registry.fill(HIST("hr1_ch_d_z"), z);
                registry.fill(HIST("hr1_ch_d_ft"), log10(fT));
                registry.fill(HIST("hr1_ch_d_mult"), ch_mult);
                registry.fill(HIST("hr1_ch_d_zg"), zg);
                registry.fill(HIST("hr1_ch_d_rg"), log(rg));
                registry.fill(HIST("hr1_ch_d_lu"), log(1 / rg), log(kt_p));
              }
            }
            ////////  for 2nd matched split
            if (nsd == 1) {
              if (ch_l == ch_nl) {
                registry.fill(HIST("hr2_ch_s_pt"), ptJet);
                registry.fill(HIST("hr2_ch_s_kt"), kt_p);
                registry.fill(HIST("hr2_ch_s_z"), z);
                registry.fill(HIST("hr2_ch_s_ft"), log10(fT));
                registry.fill(HIST("hr2_ch_s_mult"), ch_mult);
                registry.fill(HIST("hr2_ch_s_zg"), zg);
                registry.fill(HIST("hr2_ch_s_rg"), log(rg));
                registry.fill(HIST("hr2_ch_s_lu"), log(1 / rg), log(kt_p));
              }
              if (ch_l != ch_nl) {
                registry.fill(HIST("hr2_ch_d_pt"), ptJet);
                registry.fill(HIST("hr2_ch_d_kt"), kt_p);
                registry.fill(HIST("hr2_ch_d_z"), z);
                registry.fill(HIST("hr2_ch_d_ft"), log10(fT));
                registry.fill(HIST("hr2_ch_d_mult"), ch_mult);
                registry.fill(HIST("hr2_ch_d_zg"), zg);
                registry.fill(HIST("hr2_ch_d_rg"), log(rg));
                registry.fill(HIST("hr2_ch_d_lu"), log(1 / rg), log(kt_p));
              }
            }
            ////////  for 3rd marched split
            if (nsd >= 2) {
              if (ch_l == ch_nl) {
                registry.fill(HIST("hr3_ch_s_pt"), ptJet);
                registry.fill(HIST("hr3_ch_s_kt"), kt_p);
                registry.fill(HIST("hr3_ch_s_z"), z);
                registry.fill(HIST("hr3_ch_s_ft"), log10(fT));
                registry.fill(HIST("hr3_ch_s_mult"), ch_mult);
                registry.fill(HIST("hr3_ch_s_zg"), zg);
                registry.fill(HIST("hr3_ch_s_rg"), log(rg));
                registry.fill(HIST("hr3_ch_s_lu"), log(1 / rg), log(kt_p));
              }
              if (ch_l != ch_nl) {
                registry.fill(HIST("hr3_ch_d_pt"), ptJet);
                registry.fill(HIST("hr3_ch_d_kt"), kt_p);
                registry.fill(HIST("hr3_ch_d_z"), z);
                registry.fill(HIST("hr3_ch_d_ft"), log10(fT));
                registry.fill(HIST("hr3_ch_d_mult"), ch_mult);
                registry.fill(HIST("hr3_ch_d_zg"), zg);
                registry.fill(HIST("hr3_ch_d_rg"), log(rg));
                registry.fill(HIST("hr3_ch_d_lu"), log(1 / rg), log(kt_p));
              }
            }
          }
          break;
        }
        nsd++;
      }
      daughterSubJet = parentSubJet1;
    }
    if constexpr (!isSubtracted && !isMCP) {
      registry.fill(HIST("h2_jet_pt_jet_nsd"), jet.pt(), nsd);
    }
    if constexpr (!isSubtracted && isMCP) {
      registry.fill(HIST("h2_jet_pt_part_jet_nsd_part"), jet.pt(), nsd);
    }
    if constexpr (isSubtracted && !isMCP) {
      registry.fill(HIST("h2_jet_pt_jet_nsd_eventwiseconstituentsubtracted"), jet.pt(), nsd);
    }
  }

  template <bool isSubtracted, typename T, typename U>
  void analyseCharged(T const& jet, U const& /*tracks*/)
  {
    if (jet.pt() > jetPtTh) {
      jetConstituents.clear();

      int nn = 0;
      int iord[50];
      float ptc[50];
      for (auto& jetConstituent : jet.template tracks_as<aod::JTracks>()) {
        fastjetutilities::fillTracks(jetConstituent, jetConstituents, jetConstituent.globalIndex());
        ptc[nn] = jetConstituent.pt();
        nn++;
      }
      TMath::Sort(nn, ptc, iord);

      nn = 0;
      ch_mult = jet.tracksIds().size();
      for (auto& jetConstituent : jet.template tracks_as<aod::JTracks>()) {
        if (iord[nn] > 1)
          continue;

        if (iord[nn] == 0) {
          trackL = jetConstituent.globalIndex();
          ch_l = jetConstituent.sign();
          v1.SetXYZ(jetConstituent.pt(), jetConstituent.py(), jetConstituent.pz());
        }

        if (iord[nn] == 1) {
          n_trackL = jetConstituent.globalIndex();
          ch_nl = jetConstituent.sign();
          v2.SetXYZ(jetConstituent.px(), jetConstituent.py(), jetConstituent.pz());

          vR = v1 + v2;
          float z = v2.Perp(vR.Orthogonal()) / (v1.Perp(vR.Orthogonal()) + v2.Perp(vR.Orthogonal()));
          float fT = ((2. * z * (1 - z) * vR.Mag()) / v1.Perp2(vR)) / 6.;
          float kt_p = v1.Perp(vR);

          if (ch_l == ch_nl) {
            registry.fill(HIST("h_ch_s_pt"), jet.pt());
            registry.fill(HIST("h_ch_s_kt"), kt_p);
            registry.fill(HIST("h_ch_s_z"), z);
            registry.fill(HIST("h_ch_s_ft"), log10(fT));
            registry.fill(HIST("h_ch_s_mult"), ch_mult);
          }
          if (ch_l != ch_nl) {
            registry.fill(HIST("h_ch_d_pt"), jet.pt());
            registry.fill(HIST("h_ch_d_kt"), kt_p);
            registry.fill(HIST("h_ch_d_z"), z);
            registry.fill(HIST("h_ch_d_ft"), log10(fT));
            registry.fill(HIST("h_ch_d_mult"), ch_mult);
          }
        }

        nn++;
      }
      if (nn > 1)
        jetReclustering<false, isSubtracted>(jet);
    }
  }

  void processDummy(aod::JetTracks const&)
  {
  }
  PROCESS_SWITCH(JetChCorr, processDummy, "Dummy process function turned on by default", true);

  void processChargedJetsData(soa::Join<aod::ChargedJets, aod::ChargedJetConstituents>::iterator const& jet, aod::JetTracks const& tracks)
  {
    analyseCharged<false>(jet, tracks);
  }
  PROCESS_SWITCH(JetChCorr, processChargedJetsData, "charged jet substructure", false);

  void processChargedJetsEventWiseSubData(soa::Join<aod::ChargedEventWiseSubtractedJets, aod::ChargedEventWiseSubtractedJetConstituents>::iterator const& jet,
                                          aod::JetTracksSub const& tracks)
  {
    analyseCharged<true>(jet, tracks);
  }
  PROCESS_SWITCH(JetChCorr, processChargedJetsEventWiseSubData, "eventwise-constituent subtracted charged jet substructure", false);

  void processChargedJetsMCD(typename soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents>::iterator const& jet,
                             aod::JetTracks const& tracks)
  {
    analyseCharged<false>(jet, tracks);
  }
  PROCESS_SWITCH(JetChCorr, processChargedJetsMCD, "charged jet substructure", false);

  void processChargedJetsMCP(typename soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents>::iterator const& jet,
                             aod::JetParticles const&)
  {
    jetConstituents.clear();
    for (auto& jetConstituent : jet.template tracks_as<aod::JetParticles>()) {
      fastjetutilities::fillTracks(jetConstituent, jetConstituents, jetConstituent.globalIndex(), static_cast<int>(JetConstituentStatus::track), pdg->Mass(jetConstituent.pdgCode()));
    }
    jetReclustering<true, false>(jet);
  }
  PROCESS_SWITCH(JetChCorr, processChargedJetsMCP, "charged jet substructure on MC particle level", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{

  return WorkflowSpec{adaptAnalysisTask<JetChCorr>(
    cfgc, TaskName{"jet-ch-corr"})};
}
