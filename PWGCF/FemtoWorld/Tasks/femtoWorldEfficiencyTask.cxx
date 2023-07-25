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

/// \file femtoWorldEfficiencyTask.cxx
/// \author Lukasz Graczykowski, WUT Warsaw, lgraczyk@cern.ch
/// \author Malgorzata Janik, WUT Warsaw, majanik@cern.ch

// O2 includes
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"
#include "PWGLF/DataModel/LFResonanceTables.h"
#include "PWGCF/FemtoWorld/Core/FemtoWorldCollisionSelection.h"

#include "TPDGCode.h"

using namespace o2;
using namespace o2::analysis::femtoWorld;
using namespace o2::track;
using namespace o2::framework;
using namespace o2::framework::expressions;

// using TracksPID = aod::FullTracks; // This is okay.
using TracksPID = soa::Join<aod::FullTracks, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>; // for helper task with "full"
// using TracksPID = soa::Join<aod::FullTracks, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr, aod::pidTOFPi, aod::pidTOFKa, aod::pidTOFPr>; // This is okay for "no full"

using CollisionsEvSel = soa::Join<aod::Collisions, aod::EvSels, aod::Mults>;

// Femto World Efficiency task
struct femtoWorldEficiencyTask {
  // histogram registries produced
  HistogramRegistry registryQAevent{"QAHistosEvent", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry registryQAtrack{"QAHistosTrack", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry registryPID{"PIDHistos", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry registryPDG{"PDGHistos", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry registryPri{"PriHistos", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};

  HistogramRegistry registryMCtruth{"MCtruthHistos", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};

  // configurables
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "Eta range for tracks"};
  Configurable<float> pidnSigmaCut{"pidnSigmaCut", 3.0f, "TPC and TOF PID cut "};
  Configurable<float> tofPtCut{"tofPtCut", 0.5f, "From what pT TOF is used"};
  Configurable<bool> ConfIsRun3{"ConfIsRun3", false, "Running on Run 3 data"}; // Choose if running on converted data or  run 3 data
  /// Event cuts
  o2::analysis::femtoWorld::FemtoWorldCollisionSelection colCuts;
  Configurable<float> ConfEvtZvtx{"ConfEvtZvtx", 10.f, "Evt sel: Max. z-Vertex (cm)"};
  Configurable<bool> ConfEvtTriggerCheck{"ConfEvtTriggerCheck", true, "Evt sel: check for trigger"};
  Configurable<int> ConfEvtTriggerSel{"ConfEvtTriggerSel", kINT7, "Evt sel: trigger"};
  Configurable<bool> ConfEvtOfflineCheck{"ConfEvtOfflineCheck", false, "Evt sel: check for offline selection"};

  void init(InitContext&)
  {
    colCuts.setCuts(ConfEvtZvtx, ConfEvtTriggerCheck, ConfEvtTriggerSel, ConfEvtOfflineCheck, ConfIsRun3);
    colCuts.init(&registryQAevent);

    // event cuts - already done in FemtoWorldCollisionSelection.h
    // registryQAevent.add("before/reco/zvtx", "vtx_{#it{z}}", kTH1F, {{300, -15.0, 15.0, "vtx_{#it{z}} (cm)"}});
    // registryQAevent.add("before/reco/multiplicity", "V0M multiplicity class", kTH1F, {{100, 0.0, 100.0, "V0M multiplicity (%)"}});

    // track cuts
    registryQAtrack.add("after/all/plus/pt", "Charged particles #it{p}_{T}", kTH1F, {{150, 0.0, 15.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryQAtrack.add("after/all/plus/eta", "Charged particles #eta", kTH1F, {{400, -1.0, 1.0, "#eta"}});
    registryQAtrack.add("after/all/plus/phi", "Charged particles #varphi", kTH1F, {{360, 0.0, constants::math::TwoPI, "#varphi"}});
    registryQAtrack.add("after/all/plus/etaphi", "#eta - #varphi;#eta;#varphi", {HistType::kTH2F, {{200, -1, 1}, {200, 0, 2 * TMath::Pi()}}});
    registryQAtrack.addClone("after/all/plus/", "after/all/minus/");
    registryQAtrack.addClone("after/all/", "after/pion/");
    registryQAtrack.addClone("after/all/", "after/kaon/");
    registryQAtrack.addClone("after/all/", "after/proton/");

    // pid cuts

    // nsigmas
    registryPID.add("pid/plus/TOF_TPC_Map", "TOF + TPC Combined PID for all;#sigma_{TOF}^{all};#sigma_{TPC}^{all}", {HistType::kTH2F, {{100, -5, 5}, {100, -5, 5}}});
    registryPID.addClone("pid/plus/", "pid/minus/");

    registryPID.add("pid/kaon/plus/TOF_TPC_Map", "TOF + TPC Combined PID;#sigma_{TOF};#sigma_{TPC}", {HistType::kTH2F, {{100, -5, 5}, {100, -5, 5}}});
    registryPID.add("pid/kaon/plus/TOF_Nsigma", "TOF NSigma;#it{p}_{T} (GeV/#it{c});#sigma_{TOF};", {HistType::kTH2F, {{100, 0, 10}, {100, -5, 5}}});
    registryPID.add("pid/kaon/plus/TPC_Nsigma", "TPC NSigma;#it{p}_{T} (GeV/#it{c});#sigma_{TPC};", {HistType::kTH2F, {{100, 0, 10}, {100, -5, 5}}});
    registryPID.addClone("pid/kaon/plus/", "pid/kaon/minus/");

    registryPID.addClone("pid/kaon/", "pid/pion/");
    registryPID.addClone("pid/kaon/", "pid/proton/");

    // PDG
    registryPDG.add("plus/PDGPi", "PDGPi;#it{p}_{T} (GeV/c); PDG", {HistType::kTH2F, {{500, 0, 10}, {8001, -4000.5, 4000.5}}});
    registryPDG.add("plus/PDGKa", "PDGKa;#it{p}_{T} (GeV/c); PDG", {HistType::kTH2F, {{500, 0, 10}, {8001, -4000.5, 4000.5}}});
    registryPDG.add("plus/PDGPr", "PDGPr;#it{p}_{T} (GeV/c); PDG", {HistType::kTH2F, {{500, 0, 10}, {8001, -4000.5, 4000.5}}});

    registryPDG.add("minus/PDGPi", "PDGPi;#it{p}_{T} (GeV/c); PDG", {HistType::kTH2F, {{500, 0, 10}, {8001, -4000.5, 4000.5}}});
    registryPDG.add("minus/PDGKa", "PDGKa;#it{p}_{T} (GeV/c); PDG", {HistType::kTH2F, {{500, 0, 10}, {8001, -4000.5, 4000.5}}});
    registryPDG.add("minus/PDGPr", "PDGPr;#it{p}_{T} (GeV/c); PDG", {HistType::kTH2F, {{500, 0, 10}, {8001, -4000.5, 4000.5}}});

    registryPri.add("plus/PiPri", "PiPri;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 10}, {400, -1.0, 1.0}}});
    registryPri.add("plus/KaPri", "KaPri;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 10}, {400, -1.0, 1.0}}});
    registryPri.add("plus/PrPri", "PrPri;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 10}, {400, -1.0, 1.0}}});
    registryPri.add("plus/AllPri", "AllPri;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 10}, {400, -1.0, 1.0}}});

    registryPri.add("minus/PiPri", "PiPri;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 10}, {400, -1.0, 1.0}}});
    registryPri.add("minus/KaPri", "KaPri;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 10}, {400, -1.0, 1.0}}});
    registryPri.add("minus/PrPri", "PrPri;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 10}, {400, -1.0, 1.0}}});
    registryPri.add("minus/AllPri", "AllPri;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 10}, {400, -1.0, 1.0}}});

    registryPri.add("plus/PiPriPt", "PiPri;#it{p}_{T} (GeV/c)", {HistType::kTH1F, {{500, 0, 10}}});
    registryPri.add("plus/KaPriPt", "KaPri;#it{p}_{T} (GeV/c)", {HistType::kTH1F, {{500, 0, 10}}});
    registryPri.add("plus/PrPriPt", "PrPri;#it{p}_{T} (GeV/c)", {HistType::kTH1F, {{500, 0, 5}}});
    registryPri.add("plus/AllPriPt", "AllPri;#it{p}_{T} (GeV/c)", {HistType::kTH1F, {{500, 0, 10}}});

    registryPri.add("minus/PiPriPt", "PiPri;#it{p}_{T} (GeV/c)", {HistType::kTH1F, {{500, 0, 10}}});
    registryPri.add("minus/KaPriPt", "KaPri;#it{p}_{T} (GeV/c)", {HistType::kTH1F, {{500, 0, 10}}});
    registryPri.add("minus/PrPriPt", "PrPri;#it{p}_{T} (GeV/c)", {HistType::kTH1F, {{500, 0, 5}}});
    registryPri.add("minus/AllPriPt", "AllPri;#it{p}_{T} (GeV/c)", {HistType::kTH1F, {{500, 0, 10}}});

    // MC truth
    registryMCtruth.add("plus/MCtruthPi", "MC truth pions;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 10}, {400, -1.0, 1.0}}});
    registryMCtruth.add("plus/MCtruthKa", "MC truth kaons;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 10}, {400, -1.0, 1.0}}});
    registryMCtruth.add("plus/MCtruthPr", "MC truth protons;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 10}, {400, -1.0, 1.0}}});

    registryMCtruth.add("minus/MCtruthPi", "MC truth pions;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 10}, {400, -1.0, 1.0}}});
    registryMCtruth.add("minus/MCtruthKa", "MC truth kaons;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 10}, {400, -1.0, 1.0}}});
    registryMCtruth.add("minus/MCtruthPr", "MC truth protons;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 10}, {400, -1.0, 1.0}}});

    registryMCtruth.add("plus/MCtruthPiPt", "MC truth pions;#it{p}_{T} (GeV/c)", {HistType::kTH1F, {{500, 0, 10}}});
    registryMCtruth.add("plus/MCtruthKaPt", "MC truth kaons;#it{p}_{T} (GeV/c)", {HistType::kTH1F, {{500, 0, 10}}});
    registryMCtruth.add("plus/MCtruthPrPt", "MC truth protons;#it{p}_{T} (GeV/c)", {HistType::kTH1F, {{500, 0, 5}}});
    registryMCtruth.add("plus/MCtruthAllPt", "MC truth all;#it{p}_{T} (GeV/c)", {HistType::kTH1F, {{500, 0, 10}}});

    registryMCtruth.add("minus/MCtruthPiPt", "MC truth pions;#it{p}_{T} (GeV/c)", {HistType::kTH1F, {{500, 0, 10}}});
    registryMCtruth.add("minus/MCtruthKaPt", "MC truth kaons;#it{p}_{T} (GeV/c)", {HistType::kTH1F, {{500, 0, 10}}});
    registryMCtruth.add("minus/MCtruthPrPt", "MC truth protons;#it{p}_{T} (GeV/c)", {HistType::kTH1F, {{500, 0, 5}}});
    registryMCtruth.add("minus/MCtruthAllPt", "MC truth all;#it{p}_{T} (GeV/c)", {HistType::kTH1F, {{500, 0, 10}}});
  }

  bool IsNSigmaAccept(float nsigmaTPC, float nsigmaTOF, float mom)
  {

    if (mom > tofPtCut) {
      if (TMath::Hypot(nsigmaTOF, nsigmaTPC) < pidnSigmaCut)
        return true;
    } else {
      if (TMath::Abs(nsigmaTPC) < pidnSigmaCut)
        return true;
    }

    return false;
  }

  Filter trackFilter = nabs(aod::track::eta) < cfgCutEta; // Eta cut
  Filter trackCutFilter = requireGlobalTrackInFilter();   // Global track cuts
  void processReco(const CollisionsEvSel::iterator& collision,
                   soa::Filtered<TracksPID> const& tracks /*, aod::BCsWithTimestamps const&*/)
  {
    // auto bc = collision.bc_as<aod::BCsWithTimestamps>(); /// adding timestamp to access magnetic field later
    //  Default event selection
    colCuts.printCuts();
    if (!colCuts.isSelected(collision))
      return;
    colCuts.fillQA(collision);

    // Loop over tracks
    for (auto& track : tracks) {
      // Tracks are already filtered by the pre-filters

      // track.hasTOF() - to check if TOF info available

      // cuts on tracks:
      if (track.pt() > tofPtCut && !track.hasTOF())
        continue; // if no TOF information above tofPtCut reject such track

      // no PID histograms
      if (track.sign() > 0) {
        registryQAtrack.fill(HIST("after/all/plus/etaphi"), track.eta(), track.phi());
        registryQAtrack.fill(HIST("after/all/plus/pt"), track.pt());
        registryQAtrack.fill(HIST("after/all/plus/eta"), track.eta());
        registryQAtrack.fill(HIST("after/all/plus/phi"), track.phi());
      }
      if (track.sign() < 0) {
        registryQAtrack.fill(HIST("after/all/minus/etaphi"), track.eta(), track.phi());
        registryQAtrack.fill(HIST("after/all/minus/pt"), track.pt());
        registryQAtrack.fill(HIST("after/all/minus/eta"), track.eta());
        registryQAtrack.fill(HIST("after/all/minus/phi"), track.phi());
      }

      // Add PID selection criteria here
      if (IsNSigmaAccept(std::abs(track.tpcNSigmaPi()), std::abs(track.tofNSigmaPi()), track.pt())) {
        if (track.sign() > 0) {
          registryPID.fill(HIST("pid/pion/plus/TOF_TPC_Map"), track.tofNSigmaPi(), track.tpcNSigmaPi());
          registryPID.fill(HIST("pid/pion/plus/TOF_Nsigma"), track.pt(), track.tofNSigmaPi());
          registryPID.fill(HIST("pid/pion/plus/TPC_Nsigma"), track.pt(), track.tpcNSigmaPi());

          registryQAtrack.fill(HIST("after/pion/plus/etaphi"), track.eta(), track.phi());
          registryQAtrack.fill(HIST("after/pion/plus/pt"), track.pt());
          registryQAtrack.fill(HIST("after/pion/plus/eta"), track.eta());
          registryQAtrack.fill(HIST("after/pion/plus/phi"), track.phi());
        }
        if (track.sign() < 0) {
          registryPID.fill(HIST("pid/pion/minus/TOF_TPC_Map"), track.tofNSigmaPi(), track.tpcNSigmaPi());
          registryPID.fill(HIST("pid/pion/minus/TOF_Nsigma"), track.pt(), track.tofNSigmaPi());
          registryPID.fill(HIST("pid/pion/minus/TPC_Nsigma"), track.pt(), track.tpcNSigmaPi());

          registryQAtrack.fill(HIST("after/pion/minus/etaphi"), track.eta(), track.phi());
          registryQAtrack.fill(HIST("after/pion/minus/pt"), track.pt());
          registryQAtrack.fill(HIST("after/pion/minus/eta"), track.eta());
          registryQAtrack.fill(HIST("after/pion/minus/phi"), track.phi());
        }
      }
      if (IsNSigmaAccept(std::abs(track.tpcNSigmaKa()), std::abs(track.tofNSigmaKa()), track.pt())) {
        if (track.sign() > 0) {
          registryPID.fill(HIST("pid/kaon/plus/TOF_TPC_Map"), track.tofNSigmaKa(), track.tpcNSigmaKa());
          registryPID.fill(HIST("pid/kaon/plus/TOF_Nsigma"), track.pt(), track.tofNSigmaKa());
          registryPID.fill(HIST("pid/kaon/plus/TPC_Nsigma"), track.pt(), track.tpcNSigmaKa());
          registryQAtrack.fill(HIST("after/kaon/plus/etaphi"), track.eta(), track.phi());
          registryQAtrack.fill(HIST("after/kaon/plus/pt"), track.pt());
          registryQAtrack.fill(HIST("after/kaon/plus/eta"), track.eta());
          registryQAtrack.fill(HIST("after/kaon/plus/phi"), track.phi());
        }
        if (track.sign() < 0) {
          registryPID.fill(HIST("pid/kaon/minus/TOF_TPC_Map"), track.tofNSigmaKa(), track.tpcNSigmaKa());
          registryPID.fill(HIST("pid/kaon/minus/TOF_Nsigma"), track.pt(), track.tofNSigmaKa());
          registryPID.fill(HIST("pid/kaon/minus/TPC_Nsigma"), track.pt(), track.tpcNSigmaKa());
          registryQAtrack.fill(HIST("after/kaon/minus/etaphi"), track.eta(), track.phi());
          registryQAtrack.fill(HIST("after/kaon/minus/pt"), track.pt());
          registryQAtrack.fill(HIST("after/kaon/minus/eta"), track.eta());
          registryQAtrack.fill(HIST("after/kaon/minus/phi"), track.phi());
        }
      }
      if (IsNSigmaAccept(std::abs(track.tpcNSigmaPr()), std::abs(track.tofNSigmaPr()), track.pt())) {
        if (track.sign() > 0) {
          registryPID.fill(HIST("pid/proton/plus/TOF_TPC_Map"), track.tofNSigmaPr(), track.tpcNSigmaPr());
          registryPID.fill(HIST("pid/proton/plus/TOF_Nsigma"), track.pt(), track.tofNSigmaPr());
          registryPID.fill(HIST("pid/proton/plus/TPC_Nsigma"), track.pt(), track.tpcNSigmaPr());
          registryQAtrack.fill(HIST("after/proton/plus/etaphi"), track.eta(), track.phi());
          registryQAtrack.fill(HIST("after/proton/plus/pt"), track.pt());
          registryQAtrack.fill(HIST("after/proton/plus/eta"), track.eta());
          registryQAtrack.fill(HIST("after/proton/plus/phi"), track.phi());
        }
        if (track.sign() < 0) {
          registryPID.fill(HIST("pid/proton/minus/TOF_TPC_Map"), track.tofNSigmaPr(), track.tpcNSigmaPr());
          registryPID.fill(HIST("pid/proton/minus/TOF_Nsigma"), track.pt(), track.tofNSigmaPr());
          registryPID.fill(HIST("pid/proton/minus/TPC_Nsigma"), track.pt(), track.tpcNSigmaPr());
          registryQAtrack.fill(HIST("after/proton/minus/etaphi"), track.eta(), track.phi());
          registryQAtrack.fill(HIST("after/proton/minus/pt"), track.pt());
          registryQAtrack.fill(HIST("after/proton/minus/eta"), track.eta());
          registryQAtrack.fill(HIST("after/proton/minus/phi"), track.phi());
        }
      }
    }
  }
  PROCESS_SWITCH(femtoWorldEficiencyTask, processReco, "Process reconstructed data", true);

  using BigTracksMC = soa::Join<TracksPID, aod::McTrackLabels>;
  Preslice<BigTracksMC> perCollisionID = aod::track::collisionId;

  void processMCTruth(aod::McCollision const& collision, soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions>> const& collisions, aod::McParticles const& mcparticles, soa::Filtered<BigTracksMC> const& tracks)
  {
    // Loop over reconstructed collisions corresponding to MC collision
    for (auto& collision : collisions) {

      // Group tracks belonging to collision
      auto groupedTracks = tracks.sliceBy(perCollisionID, collision.globalIndex());

      // Loop over tracks
      for (auto& track : groupedTracks) {
        if (!track.has_mcParticle()) {
          continue;
        }
        const auto mcParticle = track.mcParticle();

        if (track.sign() > 0) {
          if (IsNSigmaAccept(std::abs(track.tpcNSigmaPi()), std::abs(track.tofNSigmaPi()), track.pt())) {
            registryPDG.fill(HIST("plus/PDGPi"), track.pt(), mcParticle.pdgCode());
          }
          if (IsNSigmaAccept(std::abs(track.tpcNSigmaKa()), std::abs(track.tofNSigmaKa()), track.pt())) {
            registryPDG.fill(HIST("plus/PDGKa"), track.pt(), mcParticle.pdgCode());
          }
          if (IsNSigmaAccept(std::abs(track.tpcNSigmaPr()), std::abs(track.tofNSigmaPr()), track.pt())) {
            registryPDG.fill(HIST("plus/PDGPr"), track.pt(), mcParticle.pdgCode());
          }
        }
        if (track.sign() < 0) {
          if (IsNSigmaAccept(std::abs(track.tpcNSigmaPi()), std::abs(track.tofNSigmaPi()), track.pt())) {
            registryPDG.fill(HIST("minus/PDGPi"), track.pt(), mcParticle.pdgCode());
          }
          if (IsNSigmaAccept(std::abs(track.tpcNSigmaKa()), std::abs(track.tofNSigmaKa()), track.pt())) {
            registryPDG.fill(HIST("minus/PDGKa"), track.pt(), mcParticle.pdgCode());
          }
          if (IsNSigmaAccept(std::abs(track.tpcNSigmaPr()), std::abs(track.tofNSigmaPr()), track.pt())) {
            registryPDG.fill(HIST("minus/PDGPr"), track.pt(), mcParticle.pdgCode());
          }
        }

        if (mcParticle.isPhysicalPrimary() && TMath::Abs(mcParticle.eta()) < cfgCutEta) {
          if (track.sign() > 0) {
            registryPri.fill(HIST("plus/AllPri"), track.pt(), track.eta());
            registryPri.fill(HIST("plus/AllPriPt"), track.pt());
            if (IsNSigmaAccept(std::abs(track.tpcNSigmaPi()), std::abs(track.tofNSigmaPi()), track.pt()) && mcParticle.pdgCode() == 211) {
              registryPri.fill(HIST("plus/PiPri"), track.pt(), track.eta());
              registryPri.fill(HIST("plus/PiPriPt"), mcParticle.pt());
            }
            if (IsNSigmaAccept(std::abs(track.tpcNSigmaKa()), std::abs(track.tofNSigmaKa()), track.pt()) && mcParticle.pdgCode() == 321) {
              registryPri.fill(HIST("plus/KaPri"), track.pt(), track.eta());
              registryPri.fill(HIST("plus/KaPriPt"), mcParticle.pt());
            }
            if (IsNSigmaAccept(std::abs(track.tpcNSigmaPr()), std::abs(track.tofNSigmaPr()), track.pt()) && mcParticle.pdgCode() == 2212) {
              registryPri.fill(HIST("plus/PrPri"), track.pt(), track.eta());
              registryPri.fill(HIST("plus/PrPriPt"), mcParticle.pt());
            }
          }
          if (track.sign() < 0) {
            registryPri.fill(HIST("minus/AllPri"), track.pt(), track.eta());
            registryPri.fill(HIST("minus/AllPriPt"), mcParticle.pt());
            if (IsNSigmaAccept(std::abs(track.tpcNSigmaPi()), std::abs(track.tofNSigmaPi()), track.pt()) && mcParticle.pdgCode() == -211) {
              registryPri.fill(HIST("minus/PiPri"), track.pt(), track.eta());
              registryPri.fill(HIST("minus/PiPriPt"), mcParticle.pt());
            }
            if (IsNSigmaAccept(std::abs(track.tpcNSigmaKa()), std::abs(track.tofNSigmaKa()), track.pt()) && mcParticle.pdgCode() == -321) {
              registryPri.fill(HIST("minus/KaPri"), track.pt(), track.eta());
              registryPri.fill(HIST("minus/KaPriPt"), mcParticle.pt());
            }
            if (IsNSigmaAccept(std::abs(track.tpcNSigmaPr()), std::abs(track.tofNSigmaPr()), track.pt()) && mcParticle.pdgCode() == -2212) {
              registryPri.fill(HIST("minus/PrPri"), track.pt(), track.eta());
              registryPri.fill(HIST("minus/PrPriPt"), mcParticle.pt());
            }
          }
        }
      }
    }
    // loop over MC particles
    for (auto& mcparticle : mcparticles) {
      if (!mcparticle.isPhysicalPrimary() || TMath::Abs(mcparticle.eta()) > cfgCutEta)
        continue;
      if (mcparticle.pdgCode() > 0) {
        registryMCtruth.fill(HIST("plus/MCtruthAllPt"), mcparticle.pt());
      }
      if (mcparticle.pdgCode() == 211) {
        registryMCtruth.fill(HIST("plus/MCtruthPi"), mcparticle.pt(), mcparticle.eta());
        registryMCtruth.fill(HIST("plus/MCtruthPiPt"), mcparticle.pt());
      }
      if (mcparticle.pdgCode() == 321) {
        registryMCtruth.fill(HIST("plus/MCtruthKa"), mcparticle.pt(), mcparticle.eta());
        registryMCtruth.fill(HIST("plus/MCtruthKaPt"), mcparticle.pt());
      }
      if (mcparticle.pdgCode() == 2212) {
        registryMCtruth.fill(HIST("plus/MCtruthPr"), mcparticle.pt(), mcparticle.eta());
        registryMCtruth.fill(HIST("plus/MCtruthPrPt"), mcparticle.pt());
      }

      if (mcparticle.pdgCode() < 0) {
        registryMCtruth.fill(HIST("minus/MCtruthAllPt"), mcparticle.pt());
      }
      if (mcparticle.pdgCode() == -211) {
        registryMCtruth.fill(HIST("minus/MCtruthPi"), mcparticle.pt(), mcparticle.eta());
        registryMCtruth.fill(HIST("minus/MCtruthPiPt"), mcparticle.pt());
      }
      if (mcparticle.pdgCode() == -321) {
        registryMCtruth.fill(HIST("minus/MCtruthKa"), mcparticle.pt(), mcparticle.eta());
        registryMCtruth.fill(HIST("minus/MCtruthKaPt"), mcparticle.pt());
      }
      if (mcparticle.pdgCode() == -2212) {
        registryMCtruth.fill(HIST("minus/MCtruthPr"), mcparticle.pt(), mcparticle.eta());
        registryMCtruth.fill(HIST("minus/MCtruthPrPt"), mcparticle.pt());
      }
    }
  }
  PROCESS_SWITCH(femtoWorldEficiencyTask, processMCTruth, "Process MC truth data", true);

}; // end of spectra task
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<femtoWorldEficiencyTask>(cfgc, TaskName{"femtoWorldEficiencyTask"}),
  };
}
