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
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
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
using kTracksTable =
  soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection,
            aod::TracksDCA, aod::pidTOFFullKa, aod::pidTPCFullKa,
            aod::pidTOFFullPi, aod::pidTPCFullPi>;

struct rsn_analysis {
  HistogramRegistry uHistograms{
    "Histograms",
    {},
    OutputObjHandlingPolicy::AnalysisObject};
  //
  // Create Output Objects
  void init(o2::framework::InitContext&)
  {
    // Histogram is added to the ouput registry
    //
    //  Event QA
    uHistograms.add("QA/Event/EnumEvents",
                    "Event selection;                                       ;  "
                    "                         Selected Events",
                    kTH1F, {{10, 0, 10}});
    uHistograms.add("QA/Event/VertexZ",
                    "Event selection;                                       "
                    "Vertex Z (cm);              Selected Events",
                    kTH1F, {{4000, -20, 20}});
    uHistograms.add("QA/Event/Selected/VertexZ",
                    "Event selection;                                       "
                    "Vertex Z (cm);              Selected Events",
                    kTH1F, {{4000, -20, 20}});
    //
    //  Track QA
    uHistograms.add("QA/Track/Momentum",
                    "Momentum distribution of selected Tracks;              "
                    "#it{p} (GeV/#it{c});        Selected Tracks",
                    kTH1F, {{1000, 0, 20}});
    uHistograms.add("QA/Track/TMomentum",
                    "Transverse momentum distribution of selected Tracks;   "
                    "#it{p}_{T} (GeV/#it{c});    Selected Tracks",
                    kTH1F, {{1000, 0, 20}});
    //
    //  PID QA
    //  --- Kaon
    uHistograms.add("QA/Kaon/TOF_TPC_Map",
                    "TOF + TPC Combined PID for Kaons;                      "
                    "#sigma_{TOF}^{Kaon};        #sigma_{TPC}^{Kaon}",
                    kTH2F, {{1000, -10, 10}, {1000, -10, 10}});
    uHistograms.add("QA/Kaon/TOF_Nsigma",
                    "TOF NSigma for Kaons;                                  "
                    "#it{p}_{T} (GeV/#it{c});    #sigma_{TOF}^{Kaon};",
                    kTH2F, {{1000, 0, 20}, {1000, -10, 10}});
    uHistograms.add("QA/Kaon/TPC_Nsigma",
                    "TPC NSigma for Kaons;                                  "
                    "#it{p}_{T} (GeV/#it{c});    #sigma_{TPC}^{Kaon};",
                    kTH2F, {{1000, 0, 20}, {1000, -10, 10}});
    uHistograms.add("QA/Kaon/Selected/TOF_TPC_Map",
                    "TOF + TPC Combined PID for Kaons;                      "
                    "#sigma_{TPC}^{Kaon};        #sigma_{TPC}^{Kaon}",
                    kTH2F, {{1000, -10, 10}, {1000, -10, 10}});
    uHistograms.add("QA/Kaon/Selected/TOF_Nsigma",
                    "TOF NSigma for Kaons;                                  "
                    "#it{p}_{T} (GeV/#it{c});    #sigma_{TOF}^{Kaon};",
                    kTH2F, {{1000, 0, 20}, {1000, -10, 10}});
    uHistograms.add("QA/Kaon/Selected/TPC_Nsigma",
                    "TPC NSigma for Kaons;                                  "
                    "#it{p}_{T} (GeV/#it{c});    #sigma_{TPC}^{Kaon};",
                    kTH2F, {{1000, 0, 20}, {1000, -10, 10}});
    //  --- Pion
    uHistograms.add("QA/Pion/TOF_TPC_Map",
                    "TOF + TPC Combined PID for Pions;                      "
                    "#sigma_{TOF}^{Pion};        #sigma_{TPC}^{Pion}",
                    kTH2F, {{1000, -10, 10}, {1000, -10, 10}});
    uHistograms.add("QA/Pion/TOF_Nsigma",
                    "TOF NSigma for Pions;                                  "
                    "#it{p}_{T} (GeV/#it{c});    #sigma_{TOF}^{Pion};",
                    kTH2F, {{1000, 0, 20}, {1000, -10, 10}});
    uHistograms.add("QA/Pion/TPC_Nsigma",
                    "TPC NSigma for Pions;                                  "
                    "#it{p}_{T} (GeV/#it{c});    #sigma_{TPC}^{Pion};",
                    kTH2F, {{1000, 0, 20}, {1000, -10, 10}});
    uHistograms.add("QA/Pion/Selected/TOF_TPC_Map",
                    "TOF + TPC Combined PID for Pions;                      "
                    "#sigma_{TPC}^{Pion};        #sigma_{TPC}^{Pion}",
                    kTH2F, {{1000, -10, 10}, {1000, -10, 10}});
    uHistograms.add("QA/Pion/Selected/TOF_Nsigma",
                    "TOF NSigma for Pions;                                  "
                    "#it{p}_{T} (GeV/#it{c});    #sigma_{TOF}^{Pion};",
                    kTH2F, {{1000, 0, 20}, {1000, -10, 10}});
    uHistograms.add("QA/Pion/Selected/TPC_Nsigma",
                    "TPC NSigma for Pions;                                  "
                    "#it{p}_{T} (GeV/#it{c});    #sigma_{TPC}^{Pion};",
                    kTH2F, {{1000, 0, 20}, {1000, -10, 10}});
    //
    //  Analysis Plots
    //  --- Phi
    uHistograms.add("Analysis/Phi/FullInvariantMass",
                    "K^{+}K^{-} Invariant Mass;                             "
                    "M_{K^{+}K^{-}};             Counts;",
                    kTH1F, {{1000, 0.99, 1.10}});
    uHistograms.add("Analysis/Phi/PTInvariantMass",
                    "#it{p}_{T} (GeV/#it{c}); K^{+}K^{-} Invariant Mass;       "
                    "                      M_{K^{+}K^{-}};",
                    kTH2F, {{1000, 0., 20.}, {1000, 0.99, 1.10}});
    //  --- K*
    uHistograms.add("Analysis/Kstar/FullInvariantMass",
                    "K^{#pm}#pi^{#pm} Invariant Mass;                       "
                    "M_{K^{#pm}#pi^{#pm}};       Counts;",
                    kTH1F, {{1000, 0.84, 0.95}});
    uHistograms.add("Analysis/Kstar/PTInvariantMass",
                    "#it{p}_{T} (GeV/#it{c}); K^{#pm}#pi^{#pm} Invariant Mass; "
                    "                      M_{K^{#pm}#pi^{#pm}};",
                    kTH2F, {{1000, 0., 20.}, {1000, 0.84, 0.95}});
  }
  //
  //  Configurables
  //  --- Event
  Configurable<double> kVertexZ{"kVertexCut", 10., "Vertex Z Cut"};
  //  --- PID
  //  --- --- Kaons
  Configurable<double> kKaonsTOFveto{"kKaonsTOFveto", 3.,
                                     "TOF Veto NSigma for Kaons"};
  Configurable<double> kKaonsTPCveto{"kKaonsTPCveto", 5.,
                                     "TPC NSigma for Kaons w/ TOF Veto"};
  Configurable<double> kKaonsTPCstda{"kKaonsTPCstda", 3.,
                                     "TPC NSigma for Kaons Standalone"};
  //  --- --- Pions
  Configurable<double> kPionsTOFveto{"kPionsTOFveto", 3.,
                                     "TOF Veto NSigma for Pions"};
  Configurable<double> kPionsTPCveto{"kPionsTPCveto", 5.,
                                     "TPC NSigma for Pions w/ TOF Veto"};
  Configurable<double> kPionsTPCstda{"kPionsTPCstda", 3.,
                                     "TPC NSigma for Pions Standalone"};
  //
  //  Processing A02D
  void process(kCollisionsTable::iterator const& kCollision,
               kTracksTable const& kTracks)
  {
    //
    //  Collision QA
    uHistograms.fill(HIST("QA/Event/EnumEvents"), EventSelection::kNoSelection);
    //
    //  --- Event selection: Trigger based on TVX
    if (!kCollision.sel8())
      return;
    uHistograms.fill(HIST("QA/Event/EnumEvents"),
                     EventSelection::kTVXselection);
    //
    //  --- Event selection: Vertex position
    uHistograms.fill(HIST("QA/Event/VertexZ"), kCollision.posZ());
    if (!(fabs(kCollision.posZ()) < kVertexZ))
      return;
    uHistograms.fill(HIST("QA/Event/EnumEvents"), EventSelection::kVertexCut);
    uHistograms.fill(HIST("QA/Event/Selected/VertexZ"), kCollision.posZ());
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
      if (!kCurrentTrack.isGlobalTrack())
        continue;
      uHistograms.fill(HIST("QA/Track/Momentum"), kCurrentTrack.dcaZ());
      uHistograms.fill(HIST("QA/Track/TMomentum"), kCurrentTrack.pt());
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
    TLorentzVector kPositiveDecayDaughter;
    TLorentzVector kNegativeDecayDaughter;
    TLorentzVector kResonanceCandidate;
    for (auto kPosKaon : kPosSelectedKaons) {
      for (auto kNegKaon : kNegSelectedKaons) {
        //
        kPositiveDecayDaughter.SetXYZM(get<0>(kPosKaon), get<1>(kPosKaon),
                                       get<2>(kPosKaon), .493677);
        kNegativeDecayDaughter.SetXYZM(get<0>(kNegKaon), get<1>(kNegKaon),
                                       get<2>(kNegKaon), .493677);
        kResonanceCandidate = kPositiveDecayDaughter + kNegativeDecayDaughter;
        //
        if (kResonanceCandidate.Mag() < 0.90 ||
            kResonanceCandidate.Mag() > 1.10)
          continue;
        if (fabs(kResonanceCandidate.Rapidity()) > 0.5)
          continue;
        //
        uHistograms.fill(HIST("Analysis/Phi/FullInvariantMass"),
                         kResonanceCandidate.Mag());
      }
    }
    //  Invariant Mass for K*
    for (auto kPosKaon : kPosSelectedKaons) {
      for (auto kNegKaon : kNegSelectedPions) {
        //
        kPositiveDecayDaughter.SetXYZM(get<0>(kPosKaon), get<1>(kPosKaon),
                                       get<2>(kPosKaon), .493677);
        kNegativeDecayDaughter.SetXYZM(get<0>(kNegKaon), get<1>(kNegKaon),
                                       get<2>(kNegKaon), .139570);
        kResonanceCandidate = kPositiveDecayDaughter + kNegativeDecayDaughter;
        //
        if (kResonanceCandidate.Mag() < 0.84 ||
            kResonanceCandidate.Mag() > 0.95)
          continue;
        if (fabs(kResonanceCandidate.Rapidity()) > 0.5)
          continue;
        //
        uHistograms.fill(HIST("Analysis/Kstar/FullInvariantMass"),
                         kResonanceCandidate.Mag());
      }
    }
    for (auto kPosKaon : kPosSelectedPions) {
      for (auto kNegKaon : kNegSelectedKaons) {
        //
        kPositiveDecayDaughter.SetXYZM(get<0>(kPosKaon), get<1>(kPosKaon),
                                       get<2>(kPosKaon), .139570);
        kNegativeDecayDaughter.SetXYZM(get<0>(kNegKaon), get<1>(kNegKaon),
                                       get<2>(kNegKaon), .493677);
        kResonanceCandidate = kPositiveDecayDaughter + kNegativeDecayDaughter;
        //
        if (kResonanceCandidate.Mag() < 0.84 ||
            kResonanceCandidate.Mag() > 0.95)
          continue;
        if (fabs(kResonanceCandidate.Rapidity()) > 0.5)
          continue;
        //
        uHistograms.fill(HIST("Analysis/Kstar/FullInvariantMass"),
                         kResonanceCandidate.Mag());
      }
    }
  }
  //
  bool uIsKaonSelected(kTracksTable::iterator const& kCurrentTrack)
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
  bool uIsPionSelected(kTracksTable::iterator const& kCurrentTrack)
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
  bool uIsTrackSelected() { return true; }
};

WorkflowSpec defineDataProcessing(
  ConfigContext const& cfgc) // This puts your task in the DPL workflow
{
  // Equivalent to the AddTask in AliPhysics
  WorkflowSpec workflow{adaptAnalysisTask<rsn_analysis>(cfgc)};
  return workflow;
}
