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

/// \file rsn_analysis.cxx
/// \brief  Analysis task for the measurement of resonances
/// \author Nicola Rubini <nrubini@cern.ch>

// O2 includes
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "TLorentzVector.h"

using namespace std;
using namespace o2;
using namespace o2::framework;

enum EventSelection { kNoSelection = 0,
                      kTVXselection = 1,
                      kVertexCut = 2 };

using kCollisionsTable = soa::Join<aod::Collisions, aod::EvSels>;
using kTracksTable = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection,
                               aod::TracksDCA, aod::pidTOFFullKa, aod::pidTPCFullKa,
                               aod::pidTOFFullPi, aod::pidTPCFullPi>;
using kCollisionsTableMC = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;
using kTracksTableMC = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection,
                                 aod::TracksDCA, aod::pidTOFFullKa, aod::pidTPCFullKa,
                                 aod::pidTOFFullPi, aod::pidTPCFullPi, o2::aod::McTrackLabels>;

struct rsn_analysis {
  //

  //  Configurables
  //  --- General
  Configurable<bool> kProcessData{"kProcessData", true, "ProcessData"};
  //  --- Event
  Configurable<float> kVertexZ{"kVertexCut", 10., "Vertex Z Cut"};
  //  --- Track
  Configurable<float> kTrackEta{"kTrackEta", .80, "Track eta cut"};
  Configurable<float> kTrackPT{"kTrackPT", .15, "Track Pt cut"};
  //  --- PID
  //  --- --- Kaons
  Configurable<float> kKaonsTOFveto{"kKaonsTOFveto", 3.,
                                    "TOF Veto NSigma for Kaons"};
  Configurable<float> kKaonsTPCveto{"kKaonsTPCveto", 5.,
                                    "TPC NSigma for Kaons w/ TOF Veto"};
  Configurable<float> kKaonsTPCstda{"kKaonsTPCstda", 3.,
                                    "TPC NSigma for Kaons Standalone"};
  //  --- --- Pions
  Configurable<float> kPionsTOFveto{"kPionsTOFveto", 3.,
                                    "TOF Veto NSigma for Pions"};
  Configurable<float> kPionsTPCveto{"kPionsTPCveto", 5.,
                                    "TPC NSigma for Pions w/ TOF Veto"};
  Configurable<float> kPionsTPCstda{"kPionsTPCstda", 3.,
                                    "TPC NSigma for Pions Standalone"};
  //
  HistogramRegistry uHistograms{
    "Histograms",
    {},
    OutputObjHandlingPolicy::AnalysisObject};
  //
  // Create Output Objects
  void init(o2::framework::InitContext&)
  {
    // Histogram is added to the output registry
    //
    //  Event QA
    uHistograms.add("QA/Event/EnumEvents",
                    "Event selection;"
                    ";"
                    "Selected Events",
                    kTH1F, {{10, 0, 10}});
    uHistograms.add("QA/Event/VertexZ",
                    "Event selection;"
                    "Vertex Z (cm);"
                    "Selected Events",
                    kTH1F, {{4000, -20, 20}});
    uHistograms.add("QA/Event/Selected/VertexZ",
                    "Event selection;"
                    "Vertex Z (cm);"
                    "Selected Events",
                    kTH1F, {{4000, -20, 20}});
    //
    //  Track QA
    uHistograms.add("QA/Track/Momentum",
                    "Momentum distribution of selected Tracks;"
                    "#it{p} (GeV/#it{c});"
                    "Selected Tracks",
                    kTH1F, {{1000, 0, 20}});
    uHistograms.add("QA/Track/TMomentum",
                    "Transverse momentum distribution of selected Tracks;"
                    "#it{p}_{T} (GeV/#it{c});"
                    "Selected Tracks",
                    kTH1F, {{1000, 0, 20}});
    uHistograms.add("QA/Track/dcaZ",
                    "DCA_{Z} distribution of selected Tracks;"
                    "DCA_{Z} (cm);"
                    "#it{p}_{T} (GeV/#it{c});",
                    kTH2F, {{1000, -5, 5}, {1000, 0, 20}});
    uHistograms.add("QA/Track/dcaXY",
                    "DCA_{XY} momentum distribution of selected Tracks;"
                    "DCA_{XY} (cm);"
                    "#it{p}_{T} (GeV/#it{c});",
                    kTH2F, {{1000, -5, 5}, {1000, 0, 20}});
    uHistograms.add("QA/Track/TPC_CR",
                    "Transverse momentum distribution of selected Tracks;"
                    "#it{p}_{T} (GeV/#it{c});"
                    "Selected Tracks",
                    kTH1F, {{1000, 0, 300}});
    uHistograms.add("QA/Track/TPC_CRoverFnd",
                    "Transverse momentum distribution of selected Tracks;"
                    "#it{p}_{T} (GeV/#it{c});"
                    "Selected Tracks",
                    kTH1F, {{1000, 0, 2}});
    uHistograms.add("QA/Track/TPC_Chi2overCls",
                    "Transverse momentum distribution of selected Tracks;"
                    "#it{p}_{T} (GeV/#it{c});"
                    "Selected Tracks",
                    kTH1F, {{1000, 0, 10}});
    uHistograms.add("QA/Track/ITS_Chi2overCls",
                    "Transverse momentum distribution of selected Tracks;"
                    "#it{p}_{T} (GeV/#it{c});"
                    "Selected Tracks",
                    kTH1F, {{1000, 0, 10}});
    //
    //  PID QA
    //  --- Kaon
    uHistograms.add("QA/Kaon/TOF_TPC_Map",
                    "TOF + TPC Combined PID for Kaons;"
                    "#sigma_{TOF}^{Kaon};"
                    "#sigma_{TPC}^{Kaon}",
                    kTH2F, {{1000, -10, 10}, {1000, -10, 10}});
    uHistograms.add("QA/Kaon/TOF_Nsigma",
                    "TOF NSigma for Kaons;"
                    "#it{p}_{T} (GeV/#it{c});"
                    "#sigma_{TOF}^{Kaon};",
                    kTH2F, {{1000, 0, 20}, {1000, -10, 10}});
    uHistograms.add("QA/Kaon/TPC_Nsigma",
                    "TPC NSigma for Kaons;"
                    "#it{p}_{T} (GeV/#it{c});"
                    "#sigma_{TPC}^{Kaon};",
                    kTH2F, {{1000, 0, 20}, {1000, -10, 10}});
    uHistograms.add("QA/Kaon/Selected/TOF_TPC_Map",
                    "TOF + TPC Combined PID for Kaons;"
                    "#sigma_{TPC}^{Kaon};"
                    "#sigma_{TPC}^{Kaon}",
                    kTH2F, {{1000, -10, 10}, {1000, -10, 10}});
    uHistograms.add("QA/Kaon/Selected/TOF_Nsigma",
                    "TOF NSigma for Kaons;"
                    "#it{p}_{T} (GeV/#it{c});"
                    "#sigma_{TOF}^{Kaon};",
                    kTH2F, {{1000, 0, 20}, {1000, -10, 10}});
    uHistograms.add("QA/Kaon/Selected/TPC_Nsigma",
                    "TPC NSigma for Kaons;"
                    "#it{p}_{T} (GeV/#it{c});"
                    "#sigma_{TPC}^{Kaon};",
                    kTH2F, {{1000, 0, 20}, {1000, -10, 10}});
    //  --- Pion
    uHistograms.add("QA/Pion/TOF_TPC_Map",
                    "TOF + TPC Combined PID for Pions;"
                    "#sigma_{TOF}^{Pion};"
                    "#sigma_{TPC}^{Pion}",
                    kTH2F, {{1000, -10, 10}, {1000, -10, 10}});
    uHistograms.add("QA/Pion/TOF_Nsigma",
                    "TOF NSigma for Pions;"
                    "#it{p}_{T} (GeV/#it{c});"
                    "#sigma_{TOF}^{Pion};",
                    kTH2F, {{1000, 0, 20}, {1000, -10, 10}});
    uHistograms.add("QA/Pion/TPC_Nsigma",
                    "TPC NSigma for Pions;"
                    "#it{p}_{T} (GeV/#it{c});"
                    "#sigma_{TPC}^{Pion};",
                    kTH2F, {{1000, 0, 20}, {1000, -10, 10}});
    uHistograms.add("QA/Pion/Selected/TOF_TPC_Map",
                    "TOF + TPC Combined PID for Pions;"
                    "#sigma_{TPC}^{Pion};"
                    "#sigma_{TPC}^{Pion}",
                    kTH2F, {{1000, -10, 10}, {1000, -10, 10}});
    uHistograms.add("QA/Pion/Selected/TOF_Nsigma",
                    "TOF NSigma for Pions;"
                    "#it{p}_{T} (GeV/#it{c});"
                    "#sigma_{TOF}^{Pion};",
                    kTH2F, {{1000, 0, 20}, {1000, -10, 10}});
    uHistograms.add("QA/Pion/Selected/TPC_Nsigma",
                    "TPC NSigma for Pions;"
                    "#it{p}_{T} (GeV/#it{c});"
                    "#sigma_{TPC}^{Pion};",
                    kTH2F, {{1000, 0, 20}, {1000, -10, 10}});
    //
    //  Analysis Plots
    //  --- Phi
    uHistograms.add("Analysis/Phi/FullInvariantMass",
                    "K^{+}K^{-} Invariant Mass;"
                    "M_{K^{+}K^{-}};"
                    "Counts;",
                    kTH1F, {{1000, 0.99, 1.10}});
    uHistograms.add("Analysis/Phi/PTInvariantMass",
                    "#it{p}_{T} (GeV/#it{c});"
                    "K^{+}K^{-} Invariant Mass;"
                    "#it{p}_{T} dependent K^{+}K^{-} Invariant Mass",
                    kTH2F, {{1000, 0., 20.}, {1000, 0.99, 1.10}});
    //  --- K*
    uHistograms.add("Analysis/Kstar/FullInvariantMass",
                    "K^{#mp}#pi^{#pm} Invariant Mass;"
                    "M_{K^{#mp}#pi^{#pm}};"
                    "Counts;",
                    kTH1F, {{1000, 0.70, 1.10}});
    uHistograms.add("Analysis/Kstar/PTInvariantMass",
                    "#it{p}_{T} (GeV/#it{c});"
                    "K^{#mp}#pi^{#pm} Invariant Mass;"
                    "#it{p}_{T} dependent K^{#mp}#pi^{#pm} Invariant Mass",
                    kTH2F, {{1000, 0., 20.}, {1000, 0.70, 1.10}});
    if (!kProcessData) {
      //  Resolution Plots
      //  --- Phi
      uHistograms.add("Analysis/Phi/MassResolution",
                      "K^{+}K^{-} Invariant Mass Resolution;"
                      "(M_{K^{+}K^{-}}_{Rec}-M_{K^{+}K^{-}}_{True})/M_{K^{+}K^{-}}_{True};"
                      "#it{p}_{T} (GeV/#it{c});",
                      kTH2F, {{1000, -1.00, 1.00}, {1000, 0., 20.}});
      //  --- K*
      uHistograms.add("Analysis/Kstar/MassResolution",
                      "K^{#mp}#pi^{#pm} Invariant Mass Resolution;"
                      "(M_{K^{#mp}#pi^{#pm}}_{Rec}-M_{K^{#mp}#pi^{#pm}}_{True})/M_{K^{#mp}#pi^{#pm}}_{True};"
                      "#it{p}_{T} (GeV/#it{c});",
                      kTH2F, {{200, -1.00, 1.00}, {1000, 0., 20.}});
      //  Efficiency Plots
      //  --- Phi
      uHistograms.add("Analysis/Phi/Reconstructed",
                      "K^{+}K^{-} Invariant Mass Resolution;"
                      "#it{p}_{T} (GeV/#it{c});"
                      "Counts",
                      kTH1F, {{1000, 0., 20.}});
      uHistograms.add("Analysis/Phi/Generated",
                      "K^{+}K^{-} Invariant Mass Resolution;"
                      "#it{p}_{T} (GeV/#it{c});"
                      "Counts",
                      kTH1F, {{1000, 0., 20.}});
      uHistograms.add("Analysis/Phi/True",
                      "K^{+}K^{-} Invariant Mass Resolution;"
                      "#it{p}_{T} (GeV/#it{c});"
                      "Counts",
                      kTH1F, {{1000, 0., 20.}});
      //  --- K*
      uHistograms.add("Analysis/Kstar/Reconstructed",
                      "K^{#mp}#pi^{#pm} Invariant Mass Resolution;"
                      "#it{p}_{T} (GeV/#it{c});"
                      "Counts",
                      kTH1F, {{1000, 0., 20.}});
      uHistograms.add("Analysis/Kstar/Generated",
                      "K^{#mp}#pi^{#pm} Invariant Mass Resolution;"
                      "#it{p}_{T} (GeV/#it{c});"
                      "Counts",
                      kTH1F, {{1000, 0., 20.}});
      uHistograms.add("Analysis/Kstar/True",
                      "K^{#mp}#pi^{#pm} Invariant Mass Resolution;"
                      "#it{p}_{T} (GeV/#it{c});"
                      "Counts",
                      kTH1F, {{1000, 0., 20.}});
    } else {
      //  Background Plots
      //  --- Phi
      uHistograms.add("Analysis/Phi/BKG_FullInvariantMass",
                      "K^{#pm}K^{#pm} Invariant Mass;"
                      "M_{K^{#pm}K^{#pm}};"
                      "Counts;",
                      kTH1F, {{1000, 0.99, 1.10}});
      uHistograms.add("Analysis/Phi/BKG_PTInvariantMass",
                      "#it{p}_{T} (GeV/#it{c});"
                      "M_{K^{#pm}K^{#pm}};"
                      "#it{p}_{T} dependent K^{#pm}K^{#pm} Invariant Mass",
                      kTH2F, {{1000, 0., 20.}, {1000, 0.99, 1.10}});
      //  --- K*
      uHistograms.add("Analysis/Kstar/BKG_FullInvariantMass",
                      "K^{#pm}#pi^{#pm} Invariant Mass;"
                      "M_{K^{#pm}#pi^{#pm}};"
                      "Counts;",
                      kTH1F, {{1000, 0.70, 1.10}});
      uHistograms.add("Analysis/Kstar/BKG_PTInvariantMass",
                      "#it{p}_{T} (GeV/#it{c});"
                      "K^{#pm}#pi^{#pm} Invariant Mass;"
                      "#it{p}_{T} dependent K^{#pm}#pi^{#pm} Invariant Mass",
                      kTH2F, {{1000, 0., 20.}, {1000, 0.70, 1.10}});
    }
  }
  //
  //  Processing A02D
  void processData(kCollisionsTable::iterator const& kCurrentCollision,
                   kTracksTable const& kTracks)
  {
    //
    //  Collision QA
    uHistograms.fill(HIST("QA/Event/EnumEvents"), EventSelection::kNoSelection);
    //
    //  --- Event selection: Trigger based on TVX
    if (!kCurrentCollision.sel8())
      return;
    uHistograms.fill(HIST("QA/Event/EnumEvents"),
                     EventSelection::kTVXselection);
    //
    //  --- Event selection: Vertex position
    uHistograms.fill(HIST("QA/Event/VertexZ"), kCurrentCollision.posZ());
    if (!(fabs(kCurrentCollision.posZ()) < kVertexZ))
      return;
    uHistograms.fill(HIST("QA/Event/EnumEvents"), EventSelection::kVertexCut);
    uHistograms.fill(HIST("QA/Event/Selected/VertexZ"), kCurrentCollision.posZ());
    //
    //  Storage for Kaons
    //  --- PX  PY  PZ  HasPositiveCharge
    std::vector<std::tuple<Float_t, Float_t, Float_t>> kPosSelectedKaons;
    std::vector<std::tuple<Float_t, Float_t, Float_t>> kNegSelectedKaons;
    std::vector<std::tuple<Float_t, Float_t, Float_t>> kPosSelectedPions;
    std::vector<std::tuple<Float_t, Float_t, Float_t>> kNegSelectedPions;
    //
    //  Loop on Tracks
    for (auto kCurrentTrack : kTracks) {
      //
      //  Track Selection
      if (!uIsTrackSelected(kCurrentTrack))
        continue;
      //    PID QA
      uFillFullPIDQA(kCurrentTrack);
      //
      if (uIsKaonSelected(kCurrentTrack)) {
        if (kCurrentTrack.sign() > 0)
          kPosSelectedKaons.push_back(std::tuple<Float_t, Float_t, Float_t>(
            kCurrentTrack.px(), kCurrentTrack.py(), kCurrentTrack.pz()));
        else
          kNegSelectedKaons.push_back(std::tuple<Float_t, Float_t, Float_t>(
            kCurrentTrack.px(), kCurrentTrack.py(), kCurrentTrack.pz()));
      }
      //
      if (uIsPionSelected(kCurrentTrack)) {
        if (kCurrentTrack.sign() > 0)
          kPosSelectedPions.push_back(std::tuple<Float_t, Float_t, Float_t>(
            kCurrentTrack.px(), kCurrentTrack.py(), kCurrentTrack.pz()));
        else
          kNegSelectedPions.push_back(std::tuple<Float_t, Float_t, Float_t>(
            kCurrentTrack.px(), kCurrentTrack.py(), kCurrentTrack.pz()));
      }
    }
    //
    //  Invariant Mass for Phi
    TLorentzVector lPositiveDecayDaughter;
    TLorentzVector lNegativeDecayDaughter;
    TLorentzVector lResonanceCandidate;
    Int_t iUtility = 0;
    Int_t jUtility = 0;
    for (auto kPosKaon : kPosSelectedKaons) {
      iUtility++;
      for (auto kNegKaon : kNegSelectedKaons) {
        jUtility++;
        if (iUtility == jUtility)
          continue;
        //
        lPositiveDecayDaughter.SetXYZM(get<0>(kPosKaon), get<1>(kPosKaon),
                                       get<2>(kPosKaon), .493677);
        lNegativeDecayDaughter.SetXYZM(get<0>(kNegKaon), get<1>(kNegKaon),
                                       get<2>(kNegKaon), .493677);
        lResonanceCandidate = lPositiveDecayDaughter + lNegativeDecayDaughter;
        //
        if (lResonanceCandidate.Mag() < 0.90 ||
            lResonanceCandidate.Mag() > 1.10)
          continue;
        if (fabs(lResonanceCandidate.Rapidity()) > 0.5)
          continue;
        //
        uHistograms.fill(HIST("Analysis/Phi/FullInvariantMass"),
                         lResonanceCandidate.Mag());
        uHistograms.fill(HIST("Analysis/Phi/PTInvariantMass"),
                         lResonanceCandidate.Pt(),
                         lResonanceCandidate.Mag());
      }
    }
    iUtility = 0;
    jUtility = 0;
    for (auto kPosKaon : kPosSelectedKaons) {
      iUtility++;
      for (auto kNegKaon : kPosSelectedKaons) {
        jUtility++;
        if (iUtility == jUtility)
          continue;
        //
        lPositiveDecayDaughter.SetXYZM(get<0>(kPosKaon), get<1>(kPosKaon),
                                       get<2>(kPosKaon), .493677);
        lNegativeDecayDaughter.SetXYZM(get<0>(kNegKaon), get<1>(kNegKaon),
                                       get<2>(kNegKaon), .493677);
        lResonanceCandidate = lPositiveDecayDaughter + lNegativeDecayDaughter;
        //
        if (lResonanceCandidate.Mag() < 0.90 ||
            lResonanceCandidate.Mag() > 1.10)
          continue;
        if (fabs(lResonanceCandidate.Rapidity()) > 0.5)
          continue;
        //
        uHistograms.fill(HIST("Analysis/Phi/BKG_FullInvariantMass"),
                         lResonanceCandidate.Mag());
        uHistograms.fill(HIST("Analysis/Phi/BKG_PTInvariantMass"),
                         lResonanceCandidate.Pt(),
                         lResonanceCandidate.Mag());
      }
    }
    iUtility = 0;
    jUtility = 0;
    for (auto kPosKaon : kNegSelectedKaons) {
      iUtility++;
      for (auto kNegKaon : kNegSelectedKaons) {
        jUtility++;
        if (iUtility == jUtility)
          continue;
        //
        lPositiveDecayDaughter.SetXYZM(get<0>(kPosKaon), get<1>(kPosKaon),
                                       get<2>(kPosKaon), .493677);
        lNegativeDecayDaughter.SetXYZM(get<0>(kNegKaon), get<1>(kNegKaon),
                                       get<2>(kNegKaon), .493677);
        lResonanceCandidate = lPositiveDecayDaughter + lNegativeDecayDaughter;
        //
        if (lResonanceCandidate.Mag() < 0.90 ||
            lResonanceCandidate.Mag() > 1.10)
          continue;
        if (fabs(lResonanceCandidate.Rapidity()) > 0.5)
          continue;
        //
        uHistograms.fill(HIST("Analysis/Phi/BKG_FullInvariantMass"),
                         lResonanceCandidate.Mag());
        uHistograms.fill(HIST("Analysis/Phi/BKG_PTInvariantMass"),
                         lResonanceCandidate.Pt(),
                         lResonanceCandidate.Mag());
      }
    }
    //  Invariant Mass for K*
    for (auto kPosKaon : kPosSelectedKaons) {
      for (auto kNegPion : kNegSelectedPions) {
        //
        lPositiveDecayDaughter.SetXYZM(get<0>(kPosKaon), get<1>(kPosKaon),
                                       get<2>(kPosKaon), .493677);
        lNegativeDecayDaughter.SetXYZM(get<0>(kNegPion), get<1>(kNegPion),
                                       get<2>(kNegPion), .139570);
        lResonanceCandidate = lPositiveDecayDaughter + lNegativeDecayDaughter;
        //
        if (lResonanceCandidate.Mag() < 0.70 ||
            lResonanceCandidate.Mag() > 1.10)
          continue;
        if (fabs(lResonanceCandidate.Rapidity()) > 0.5)
          continue;
        //
        uHistograms.fill(HIST("Analysis/Kstar/FullInvariantMass"),
                         lResonanceCandidate.Mag());
        uHistograms.fill(HIST("Analysis/Kstar/PTInvariantMass"),
                         lResonanceCandidate.Pt(),
                         lResonanceCandidate.Mag());
      }
    }
    for (auto kPosPion : kPosSelectedPions) {
      for (auto kNegKaon : kNegSelectedKaons) {
        //
        lPositiveDecayDaughter.SetXYZM(get<0>(kPosPion), get<1>(kPosPion),
                                       get<2>(kPosPion), .139570);
        lNegativeDecayDaughter.SetXYZM(get<0>(kNegKaon), get<1>(kNegKaon),
                                       get<2>(kNegKaon), .493677);
        lResonanceCandidate = lPositiveDecayDaughter + lNegativeDecayDaughter;
        //
        if (lResonanceCandidate.Mag() < 0.70 ||
            lResonanceCandidate.Mag() > 1.10)
          continue;
        if (fabs(lResonanceCandidate.Rapidity()) > 0.5)
          continue;
        //
        uHistograms.fill(HIST("Analysis/Kstar/FullInvariantMass"),
                         lResonanceCandidate.Mag());
        uHistograms.fill(HIST("Analysis/Kstar/PTInvariantMass"),
                         lResonanceCandidate.Pt(),
                         lResonanceCandidate.Mag());
      }
    }
    for (auto kPosKaon : kPosSelectedKaons) {
      for (auto kNegPion : kPosSelectedPions) {
        //
        lPositiveDecayDaughter.SetXYZM(get<0>(kPosKaon), get<1>(kPosKaon),
                                       get<2>(kPosKaon), .493677);
        lNegativeDecayDaughter.SetXYZM(get<0>(kNegPion), get<1>(kNegPion),
                                       get<2>(kNegPion), .139570);
        lResonanceCandidate = lPositiveDecayDaughter + lNegativeDecayDaughter;
        //
        if (lResonanceCandidate.Mag() < 0.70 ||
            lResonanceCandidate.Mag() > 1.10)
          continue;
        if (fabs(lResonanceCandidate.Rapidity()) > 0.5)
          continue;
        //
        uHistograms.fill(HIST("Analysis/Kstar/BKG_FullInvariantMass"),
                         lResonanceCandidate.Mag());
        uHistograms.fill(HIST("Analysis/Kstar/BKG_PTInvariantMass"),
                         lResonanceCandidate.Pt(),
                         lResonanceCandidate.Mag());
      }
    }
    for (auto kPosPion : kNegSelectedPions) {
      for (auto kNegKaon : kNegSelectedKaons) {
        //
        lPositiveDecayDaughter.SetXYZM(get<0>(kPosPion), get<1>(kPosPion),
                                       get<2>(kPosPion), .139570);
        lNegativeDecayDaughter.SetXYZM(get<0>(kNegKaon), get<1>(kNegKaon),
                                       get<2>(kNegKaon), .493677);
        lResonanceCandidate = lPositiveDecayDaughter + lNegativeDecayDaughter;
        //
        if (lResonanceCandidate.Mag() < 0.70 ||
            lResonanceCandidate.Mag() > 1.10)
          continue;
        if (fabs(lResonanceCandidate.Rapidity()) > 0.5)
          continue;
        //
        uHistograms.fill(HIST("Analysis/Kstar/BKG_FullInvariantMass"),
                         lResonanceCandidate.Mag());
        uHistograms.fill(HIST("Analysis/Kstar/BKG_PTInvariantMass"),
                         lResonanceCandidate.Pt(),
                         lResonanceCandidate.Mag());
      }
    }
  }
  PROCESS_SWITCH(rsn_analysis, processData, "Process Data", true);
  //
  //  Processing Monte Carlo A02D
  void processMC(kCollisionsTableMC::iterator const& kCurrentCollision,
                 kTracksTableMC const& kTracks,
                 aod::McParticles const& /*mcParticles*/)
  {
    //
    //  Collision QA
    uHistograms.fill(HIST("QA/Event/EnumEvents"), EventSelection::kNoSelection);
    //
    //  --- Event selection: Trigger based on TVX
    if (!kCurrentCollision.sel8())
      return;
    uHistograms.fill(HIST("QA/Event/EnumEvents"),
                     EventSelection::kTVXselection);
    //
    //  --- Event selection: Vertex position
    uHistograms.fill(HIST("QA/Event/VertexZ"), kCurrentCollision.posZ());
    if (!(fabs(kCurrentCollision.posZ()) < kVertexZ))
      return;
    uHistograms.fill(HIST("QA/Event/EnumEvents"), EventSelection::kVertexCut);
    uHistograms.fill(HIST("QA/Event/Selected/VertexZ"), kCurrentCollision.posZ());
    //
    //  Storage for Kaons
    //  --- PX  PY  PZ  HasPositiveCharge
    std::vector<std::tuple<Float_t, Float_t, Float_t, o2::aod::McParticles::iterator>> kPosSelectedKaons;
    std::vector<std::tuple<Float_t, Float_t, Float_t, o2::aod::McParticles::iterator>> kNegSelectedKaons;
    std::vector<std::tuple<Float_t, Float_t, Float_t, o2::aod::McParticles::iterator>> kPosSelectedPions;
    std::vector<std::tuple<Float_t, Float_t, Float_t, o2::aod::McParticles::iterator>> kNegSelectedPions;
    //
    //  Loop on Tracks
    for (auto kCurrentTrack : kTracks) {
      //
      //  Track Selection
      if (!uIsTrackSelected(kCurrentTrack))
        continue;
      //    PID QA
      uFillFullPIDQA(kCurrentTrack);
      //
      if (!kCurrentTrack.has_mcParticle())
        continue;
      //
      if (uIsKaonSelected(kCurrentTrack)) {
        if (kCurrentTrack.sign() > 0)
          kPosSelectedKaons.push_back(std::tuple<Float_t, Float_t, Float_t, o2::aod::McParticles::iterator>(
            kCurrentTrack.px(), kCurrentTrack.py(), kCurrentTrack.pz(), kCurrentTrack.mcParticle()));
        else
          kNegSelectedKaons.push_back(std::tuple<Float_t, Float_t, Float_t, o2::aod::McParticles::iterator>(
            kCurrentTrack.px(), kCurrentTrack.py(), kCurrentTrack.pz(), kCurrentTrack.mcParticle()));
      }
      //
      if (uIsPionSelected(kCurrentTrack)) {
        if (kCurrentTrack.sign() > 0)
          kPosSelectedPions.push_back(std::tuple<Float_t, Float_t, Float_t, o2::aod::McParticles::iterator>(
            kCurrentTrack.px(), kCurrentTrack.py(), kCurrentTrack.pz(), kCurrentTrack.mcParticle()));
        else
          kNegSelectedPions.push_back(std::tuple<Float_t, Float_t, Float_t, o2::aod::McParticles::iterator>(
            kCurrentTrack.px(), kCurrentTrack.py(), kCurrentTrack.pz(), kCurrentTrack.mcParticle()));
      }
    }
    //
    //  Invariant Mass for Phi
    TLorentzVector lPositiveDecayDaughter;
    TLorentzVector lNegativeDecayDaughter;
    TLorentzVector lPositiveDecayDaughterMC;
    TLorentzVector lNegativeDecayDaughterMC;
    TLorentzVector lResonanceCandidate;
    TLorentzVector lResonanceCandidateMC;
    Int_t iUtility = 0;
    Int_t jUtility = 0;
    //
    for (auto kPosKaon : kPosSelectedKaons) {
      iUtility++;
      for (auto kNegKaon : kNegSelectedKaons) {
        jUtility++;
        if (iUtility == jUtility)
          continue;
        //
        lPositiveDecayDaughter.SetXYZM(get<0>(kPosKaon), get<1>(kPosKaon),
                                       get<2>(kPosKaon), .493677);
        lNegativeDecayDaughter.SetXYZM(get<0>(kNegKaon), get<1>(kNegKaon),
                                       get<2>(kNegKaon), .493677);
        lResonanceCandidate = lPositiveDecayDaughter + lNegativeDecayDaughter;
        //
        if (lResonanceCandidate.Mag() < 0.90 ||
            lResonanceCandidate.Mag() > 1.10)
          continue;
        if (fabs(lResonanceCandidate.Rapidity()) > 0.5)
          continue;
        //
        uHistograms.fill(HIST("Analysis/Phi/FullInvariantMass"),
                         lResonanceCandidate.Mag());
        uHistograms.fill(HIST("Analysis/Phi/PTInvariantMass"),
                         lResonanceCandidate.Pt(),
                         lResonanceCandidate.Mag());
        auto kPositiveDecayDaughterMC = get<3>(kPosKaon);
        auto kNegativeDecayDaughterMC = get<3>(kNegKaon);
        //
        //  Check they both have mothers
        if (!kPositiveDecayDaughterMC.has_mothers() || !kNegativeDecayDaughterMC.has_mothers())
          continue;
        //  Check they both have 1 mother
        if ((kPositiveDecayDaughterMC.mothersIds().size() != 1) || (kNegativeDecayDaughterMC.mothersIds().size() != 1))
          continue;
        //  Check they both have the same mother
        if (*kPositiveDecayDaughterMC.mothersIds().begin() == *kNegativeDecayDaughterMC.mothersIds().begin())
          continue;
        //  Check this mother is a phi meson within rapidity
        auto kResonanceMCTruthParticle = *kPositiveDecayDaughterMC.mothers_as<aod::McParticles>().begin();
        if (kResonanceMCTruthParticle.pdgCode() != 333)
          continue;
        if (fabs(kResonanceMCTruthParticle.y()) > 0.5)
          continue;
        //
        lPositiveDecayDaughterMC.SetXYZM(kPositiveDecayDaughterMC.px(), kPositiveDecayDaughterMC.py(), kPositiveDecayDaughterMC.pz(), .493677);
        lNegativeDecayDaughterMC.SetXYZM(kNegativeDecayDaughterMC.px(), kNegativeDecayDaughterMC.py(), kNegativeDecayDaughterMC.pz(), .493677);
        lResonanceCandidateMC = lPositiveDecayDaughterMC + lNegativeDecayDaughterMC;
        //
        uHistograms.fill(HIST("Analysis/Phi/MassResolution"),
                         lResonanceCandidate.Mag() - lResonanceCandidateMC.Mag(),
                         lResonanceCandidate.Pt());
        //
        uHistograms.fill(HIST("Analysis/Phi/Reconstructed"),
                         kResonanceMCTruthParticle.pt());
      }
    }
    //  Invariant Mass for K*
    for (auto kPosKaon : kPosSelectedKaons) {
      for (auto kNegPion : kNegSelectedPions) {
        //
        lPositiveDecayDaughter.SetXYZM(get<0>(kPosKaon), get<1>(kPosKaon),
                                       get<2>(kPosKaon), .493677);
        lNegativeDecayDaughter.SetXYZM(get<0>(kNegPion), get<1>(kNegPion),
                                       get<2>(kNegPion), .139570);
        lResonanceCandidate = lPositiveDecayDaughter + lNegativeDecayDaughter;
        //
        if (lResonanceCandidate.Mag() < 0.70 ||
            lResonanceCandidate.Mag() > 1.10)
          continue;
        if (fabs(lResonanceCandidate.Rapidity()) > 0.5)
          continue;
        //
        uHistograms.fill(HIST("Analysis/Kstar/FullInvariantMass"),
                         lResonanceCandidate.Mag());
        uHistograms.fill(HIST("Analysis/Kstar/PTInvariantMass"),
                         lResonanceCandidate.Pt(),
                         lResonanceCandidate.Mag());
        auto kPositiveDecayDaughterMC = get<3>(kPosKaon);
        auto kNegativeDecayDaughterMC = get<3>(kNegPion);
        //
        //  Check they both have mothers
        if (!kPositiveDecayDaughterMC.has_mothers() || !kNegativeDecayDaughterMC.has_mothers())
          continue;
        //  Check they both have 1 mother
        if ((kPositiveDecayDaughterMC.mothersIds().size() != 1) || (kNegativeDecayDaughterMC.mothersIds().size() != 1))
          continue;
        //  Check they both have the same mother
        if (*kPositiveDecayDaughterMC.mothersIds().begin() == *kNegativeDecayDaughterMC.mothersIds().begin())
          continue;
        //  Check this mother is a phi meson within rapidity
        auto kResonanceMCTruthParticle = *kPositiveDecayDaughterMC.mothers_as<aod::McParticles>().begin();
        if (fabs(kResonanceMCTruthParticle.pdgCode()) != 313)
          continue;
        if (fabs(kResonanceMCTruthParticle.y()) > 0.5)
          continue;
        //
        lPositiveDecayDaughterMC.SetXYZM(kPositiveDecayDaughterMC.px(), kPositiveDecayDaughterMC.py(), kPositiveDecayDaughterMC.pz(), .493677);
        lNegativeDecayDaughterMC.SetXYZM(kNegativeDecayDaughterMC.px(), kNegativeDecayDaughterMC.py(), kNegativeDecayDaughterMC.pz(), .139570);
        lResonanceCandidateMC = lPositiveDecayDaughterMC + lNegativeDecayDaughterMC;
        //
        uHistograms.fill(HIST("Analysis/Kstar/MassResolution"),
                         lResonanceCandidate.Mag() - lResonanceCandidateMC.Mag(),
                         lResonanceCandidate.Pt());
        uHistograms.fill(HIST("Analysis/Kstar/Reconstructed"),
                         kResonanceMCTruthParticle.pt());
      }
    }
  }
  PROCESS_SWITCH(rsn_analysis, processMC, "Process Monte Carlo", false);
  //
  void processMCTruth(aod::McCollision const& /*mcCollision*/,
                      aod::McParticles const& mcParticles)
  {

    // Loop on all mc particles
    for (auto kCurrentParticle : mcParticles) {
      //
      if (!kCurrentParticle.producedByGenerator())
        continue;
      if (fabs(kCurrentParticle.y()) > 0.5)
        continue;
      //
      if (abs(kCurrentParticle.pdgCode()) == 333) {
        uHistograms.fill(HIST("Analysis/Phi/True"),
                         kCurrentParticle.pt());
        auto kDaughters = kCurrentParticle.daughters_as<aod::McParticles>();
        auto kHasKaonp = false;
        auto kHasKaonm = false;
        if (kDaughters.size() == 2) {
          for (auto kCurrentDaughter : kDaughters) {
            if (kCurrentDaughter.pdgCode() == +321)
              kHasKaonp = true;
            if (kCurrentDaughter.pdgCode() == -321)
              kHasKaonm = true;
          }
        }
        if (kHasKaonp && kHasKaonm)
          uHistograms.fill(HIST("Analysis/Phi/Generated"),
                           kCurrentParticle.pt());
      }
      //
      if (abs(kCurrentParticle.pdgCode()) == 313) {
        uHistograms.fill(HIST("Analysis/Kstar/True"),
                         kCurrentParticle.pt());
        auto kDaughters = kCurrentParticle.daughters_as<aod::McParticles>();
        auto kHasKaonp = false;
        auto kHasKaonm = false;
        if (kDaughters.size() == 2) {
          for (auto kCurrentDaughter : kDaughters) {
            if (fabs(kCurrentDaughter.pdgCode()) == 321)
              kHasKaonp = true;
            if (fabs(kCurrentDaughter.pdgCode()) == 211)
              kHasKaonm = true;
          }
        }
        if (kHasKaonp && kHasKaonm)
          uHistograms.fill(HIST("Analysis/Kstar/Generated"),
                           kCurrentParticle.pt());
      }
    }
  }
  PROCESS_SWITCH(rsn_analysis, processMCTruth, "Process Monte Carlo Truth", false);
  //
  template <typename kCurrentTrackType>
  bool uIsKaonSelected(kCurrentTrackType const& kCurrentTrack)
  {
    //
    //  Utility variables
    Bool_t kHasTPC = kCurrentTrack.hasTPC();
    Double_t kTPCSigma = kCurrentTrack.tpcNSigmaKa();
    Bool_t kHasTOF = kCurrentTrack.hasTOF();
    Double_t kTOFSigma = kCurrentTrack.tofNSigmaKa();
    Double_t kTMomentum = kCurrentTrack.pt();
    //
    //  Selection
    if (!kHasTPC || fabs(kTPCSigma) > kKaonsTPCveto)
      return false;
    if (!kHasTOF && fabs(kTPCSigma) > kKaonsTPCstda)
      return false;
    if (kHasTOF && fabs(kTOFSigma) > kKaonsTOFveto)
      return false;
    //
    uHistograms.fill(HIST("QA/Kaon/Selected/TOF_Nsigma"), kTMomentum,
                     kTOFSigma);
    uHistograms.fill(HIST("QA/Kaon/Selected/TPC_Nsigma"), kTMomentum,
                     kTPCSigma);
    uHistograms.fill(HIST("QA/Kaon/Selected/TOF_TPC_Map"), kTOFSigma,
                     kTPCSigma);
    return true;
  }
  //
  template <typename kCurrentTrackType>
  bool uIsPionSelected(kCurrentTrackType const& kCurrentTrack)
  {
    //
    //  Utility variables
    Bool_t kHasTPC = kCurrentTrack.hasTPC();
    Double_t kTPCSigma = kCurrentTrack.tpcNSigmaPi();
    Bool_t kHasTOF = kCurrentTrack.hasTOF();
    Double_t kTOFSigma = kCurrentTrack.tofNSigmaPi();
    Double_t kTMomentum = kCurrentTrack.pt();
    //
    //  Selection
    if (!kHasTPC || fabs(kTPCSigma) > kPionsTPCveto)
      return false;
    if (!kHasTOF && fabs(kTPCSigma) > kPionsTPCstda)
      return false;
    if (kHasTOF && fabs(kTOFSigma) > kPionsTOFveto)
      return false;
    //
    uHistograms.fill(HIST("QA/Pion/Selected/TOF_Nsigma"), kTMomentum,
                     kTOFSigma);
    uHistograms.fill(HIST("QA/Pion/Selected/TPC_Nsigma"), kTMomentum,
                     kTPCSigma);
    uHistograms.fill(HIST("QA/Pion/Selected/TOF_TPC_Map"), kTOFSigma,
                     kTPCSigma);
    return true;
  }
  //
  template <typename kCurrentTrackType>
  bool uIsTrackSelected(kCurrentTrackType const& kCurrentTrack)
  {
    //
    //    Selection
    if (!kCurrentTrack.isGlobalTrack())
      return false;
    //
    //    Custom Selection
    // if ( fabs(kCurrentTrack.eta()) > kTrackEta ) return false;
    // if ( fabs(kCurrentTrack.pt()) > kTrackPT ) return false;
    //
    //    QA Plots
    //    --- General
    uHistograms.fill(HIST("QA/Track/Momentum"), kCurrentTrack.p());
    uHistograms.fill(HIST("QA/Track/TMomentum"), kCurrentTrack.pt());
    uHistograms.fill(HIST("QA/Track/dcaZ"), kCurrentTrack.dcaZ(), kCurrentTrack.pt());
    uHistograms.fill(HIST("QA/Track/dcaXY"), kCurrentTrack.dcaXY(), kCurrentTrack.pt());
    //    --- TPC
    uHistograms.fill(HIST("QA/Track/TPC_CR"), kCurrentTrack.tpcNClsCrossedRows());
    uHistograms.fill(HIST("QA/Track/TPC_CRoverFnd"), kCurrentTrack.tpcCrossedRowsOverFindableCls());
    uHistograms.fill(HIST("QA/Track/TPC_Chi2overCls"), kCurrentTrack.tpcChi2NCl());
    //    --- ITS
    uHistograms.fill(HIST("QA/Track/ITS_Chi2overCls"), kCurrentTrack.itsChi2NCl());
    return true;
  }
  //
  template <typename kCurrentTrackType>
  void uFillFullPIDQA(kCurrentTrackType const& kCurrentTrack)
  {
    //
    //  PID Selection
    //  --- PID QA Kaons
    uHistograms.fill(HIST("QA/Kaon/TOF_Nsigma"), kCurrentTrack.pt(),
                     kCurrentTrack.tofNSigmaKa());
    uHistograms.fill(HIST("QA/Kaon/TPC_Nsigma"), kCurrentTrack.pt(),
                     kCurrentTrack.tpcNSigmaKa());
    uHistograms.fill(HIST("QA/Kaon/TOF_TPC_Map"), kCurrentTrack.tofNSigmaKa(),
                     kCurrentTrack.tpcNSigmaKa());
    //  --- PID QA Pions
    uHistograms.fill(HIST("QA/Pion/TOF_Nsigma"), kCurrentTrack.pt(),
                     kCurrentTrack.tofNSigmaPi());
    uHistograms.fill(HIST("QA/Pion/TPC_Nsigma"), kCurrentTrack.pt(),
                     kCurrentTrack.tpcNSigmaPi());
    uHistograms.fill(HIST("QA/Pion/TOF_TPC_Map"), kCurrentTrack.tofNSigmaPi(),
                     kCurrentTrack.tpcNSigmaPi());
    //
  }
};

WorkflowSpec defineDataProcessing(
  ConfigContext const& cfgc) // This puts your task in the DPL workflow
{
  // Equivalent to the AddTask in AliPhysics
  WorkflowSpec workflow{adaptAnalysisTask<rsn_analysis>(cfgc)};
  return workflow;
}
