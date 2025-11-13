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
/// \author Alicja Plachta, WUT Warsaw, alicja.plachta@cern.ch
/// \author Barbara Chytla, WUT Warsaw, barbara.chytla@cern.ch

// O2 includes
#include "PWGCF/FemtoWorld/Core/FemtoWorldCollisionSelection.h"
#include "PWGLF/DataModel/LFResonanceTables.h"

#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include "TPDGCode.h"

using namespace o2;
using namespace o2::analysis::femtoWorld;
using namespace o2::track;
using namespace o2::framework;
using namespace o2::framework::expressions;

using TracksPID = soa::Join<aod::FullTracks, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullDe, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFFullDe>; // for helper task with "full"

using CollisionsEvSel = soa::Join<aod::Collisions, aod::EvSels, aod::Mults>;

// Femto World Efficiency task
struct femtoWorldEficiencyTaskDe {

  Service<o2::framework::O2DatabasePDG> pdgDB;

  // histogram registries produced
  HistogramRegistry registryQAevent{"QAHistosEvent", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry registryQAtrack{"QAHistosTrack", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry registryPID{"PIDHistos", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry registryPDG{"PDGHistos", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry registryPri{"PriHistos", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry registryPriCuts{"PriHistosCuts", {}, OutputObjHandlingPolicy::AnalysisObject, false, true}; // for tracking efficiency = cuts only
  HistogramRegistry registryGlobal{"GlobalTrackHistos", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};

  HistogramRegistry registryMCtruth{"MCtruthHistos", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};

  // configurables
  Configurable<float> pidnSigmaCut{"pidnSigmaCut", 5.0f, "TPC and TOF PID cut"};
  Configurable<float> tofPtCut{"tofPtCut", 0.5f, "From what pT TOF is used"};
  Configurable<bool> ConfIsRun3{"ConfIsRun3", false, "Running on Run 3 data"}; // Choose if running on converted data or  run 3 data
  // track cuts
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "Eta range for tracks"};
  Configurable<float> cfgPtLow{"cfgPtLow", 0.2, "Lower limit for Pt"};
  Configurable<float> cfgPtHigh{"cfgPtHigh", 4., "Higher limit for Pt"};
  Configurable<float> cfgDcaXY{"cfgDcaXY", 2.4, "Value of max. DCA_XY"};
  Configurable<float> cfgDcaZ{"cfgDcaZ", 3.2, "Value of max. DCA_Z"};
  // Event cuts
  o2::analysis::femtoWorld::FemtoWorldCollisionSelection colCuts;
  Configurable<float> ConfEvtZvtx{"ConfEvtZvtx", 10.f, "Evt sel: Max. z-Vertex (cm)"};
  Configurable<bool> ConfEvtTriggerCheck{"ConfEvtTriggerCheck", true, "Evt sel: check for trigger"};
  Configurable<int> ConfEvtTriggerSel{"ConfEvtTriggerSel", kINT7, "Evt sel: trigger"};
  Configurable<bool> ConfEvtOfflineCheck{"ConfEvtOfflineCheck", false, "Evt sel: check for offline selection"};

  Configurable<bool> ConfPIDnoTOF{"ConfPIDnoTOF", false, "Use only TPC in PID, no TOF"};

  void init(InitContext&)
  {
    colCuts.setCuts(ConfEvtZvtx, ConfEvtTriggerCheck, ConfEvtTriggerSel, ConfEvtOfflineCheck, ConfIsRun3);
    colCuts.init(&registryQAevent);

    // event cuts - already done in FemtoWorldCollisionSelection.h

    // track cuts
    registryQAtrack.add("after/all/plus/pt", "Charged particles #it{p}_{T}", kTH1F, {{150, 0.0, 15.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryQAtrack.add("after/all/plus/eta", "Charged particles #eta", kTH1F, {{400, -1.0, 1.0, "#eta"}});
    registryQAtrack.add("after/all/plus/phi", "Charged particles #varphi", kTH1F, {{360, 0.0, constants::math::TwoPI, "#varphi"}});
    registryQAtrack.add("after/all/plus/etaphi", "#eta - #varphi;#eta;#varphi", {HistType::kTH2F, {{200, -1, 1}, {200, 0, 2 * TMath::Pi()}}});
    registryQAtrack.add("after/all/plus/DCAxy", "DCAyx; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", kTH2F, {{1000, 0, 10}, {1000, -5, 5}});
    registryQAtrack.add("after/all/plus/DCAz", "DCAz; #it{p}_{T} (GeV/#it{c}); DCA_{z} (cm)", kTH2F, {{1000, 0, 10}, {1000, -5, 5}});
    registryQAtrack.addClone("after/all/plus/", "after/all/minus/");
    registryQAtrack.addClone("after/all/", "after/pion/");
    registryQAtrack.addClone("after/all/", "after/kaon/");
    registryQAtrack.addClone("after/all/", "after/deuteron/");

    // global track histos
    registryGlobal.add("crossedRows", "Number of crossed rows TPC", kTH1F, {{159, 0, 158, "N crossed rows"}});
    registryGlobal.add("RowsOverClustersTPC", "Ratio of crossed rows over findable TPC clusters", kTH1F, {{100, 0.5, 2., "N crossed rows"}});
    registryGlobal.add("chi2TPC", "Chi2 per TPC cluster", kTH1F, {{500, 0, 5, "Chi2 per TPC cluster"}});
    registryGlobal.add("chi2ITS", "Chi2 per ITS cluster", kTH1F, {{500, 0, 40, "Chi2 per ITS cluster"}});
    registryGlobal.add("dcaZ", "DCA to z vertex", kTH1F, {{500, 0, 4., "DCA to z vertex"}});
    registryGlobal.add("dcaXY", "DCA to xy vertex", kTH1F, {{500, 0, 4., "DCA to xy vertex"}});
    registryGlobal.add("clustersITS", "Number of ITS clusters", kTH1F, {{8, 0, 7, "N crossed rows"}}); // perhaps change to itsClusterMap
    registryGlobal.add("pt", "#it{p}_{T}", kTH1F, {{150, 0.0, 5.0, "#it{p}_{T} (GeV/#it{c})"}});
    registryGlobal.add("eta", "#eta", kTH1F, {{400, -1.0, 1.0, "#eta"}});

    // nsigmas
    registryPID.add("pid/plus/TOF_TPC_Map", "TOF + TPC Combined PID for all;#sigma_{TOF}^{all};#sigma_{TPC}^{all}", {HistType::kTH2F, {{100, -5, 5}, {100, -5, 5}}});
    registryPID.addClone("pid/plus/", "pid/minus/");

    registryPID.add("pid/kaon/plus/TOF_TPC_Map", "TOF + TPC Combined PID;#sigma_{TOF};#sigma_{TPC}", {HistType::kTH2F, {{100, -5, 5}, {100, -5, 5}}});
    registryPID.add("pid/kaon/plus/TOF_Nsigma", "TOF NSigma;#it{p}_{T} (GeV/#it{c});#sigma_{TOF};", {HistType::kTH2F, {{100, 0, 10}, {100, -5, 5}}});
    registryPID.add("pid/kaon/plus/TPC_Nsigma", "TPC NSigma;#it{p}_{T} (GeV/#it{c});#sigma_{TPC};", {HistType::kTH2F, {{100, 0, 10}, {100, -5, 5}}});
    registryPID.add("pid/kaon/plus/TPC_dEdx", "TPC dE/dx;#it{p}_{T} (GeV/#it{c});#it{dE/dx};", {HistType::kTH2F, {{100, 0, 10}, {100, -50, 200}}});
    registryPID.addClone("pid/kaon/plus/", "pid/kaon/minus/");

    registryPID.addClone("pid/kaon/", "pid/pion/");
    registryPID.addClone("pid/kaon/", "pid/deuteron/");

    // PDG
    registryPDG.add("plus/PDGPi", "PDGPi;#it{p}_{T} (GeV/c); PDG", {HistType::kTH2F, {{500, 0, 5}, {8001, -4000.5, 4000.5}}});
    registryPDG.add("plus/PDGKa", "PDGKa;#it{p}_{T} (GeV/c); PDG", {HistType::kTH2F, {{500, 0, 5}, {8001, -4000.5, 4000.5}}});
    registryPDG.add("plus/PDGDe", "PDGDe;#it{p}_{T} (GeV/c); PDG", {HistType::kTH2F, {{500, 0, 5}, {8001, -4000.5, 4000.5}}});

    registryPDG.add("minus/PDGPi", "PDGPi;#it{p}_{T} (GeV/c); PDG", {HistType::kTH2F, {{500, 0, 5}, {8001, -4000.5, 4000.5}}});
    registryPDG.add("minus/PDGKa", "PDGKa;#it{p}_{T} (GeV/c); PDG", {HistType::kTH2F, {{500, 0, 5}, {8001, -4000.5, 4000.5}}});
    registryPDG.add("minus/PDGDe", "PDGDe;#it{p}_{T} (GeV/c); PDG", {HistType::kTH2F, {{500, 0, 5}, {8001, -4000.5, 4000.5}}});

    // Pri
    registryPri.add("plus/PiPri", "PiPri;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 5}, {400, -1.0, 1.0}}});
    registryPri.add("plus/KaPri", "KaPri;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 5}, {400, -1.0, 1.0}}});
    registryPri.add("plus/DePri", "DePri;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 5}, {400, -1.0, 1.0}}});
    registryPri.add("plus/AllPri", "AllPri;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 5}, {400, -1.0, 1.0}}});

    registryPri.add("minus/PiPri", "PiPri;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 5}, {400, -1.0, 1.0}}});
    registryPri.add("minus/KaPri", "KaPri;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 5}, {400, -1.0, 1.0}}});
    registryPri.add("minus/DePri", "DePri;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 5}, {400, -1.0, 1.0}}});
    registryPri.add("minus/AllPri", "AllPri;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 5}, {400, -1.0, 1.0}}});

    registryPri.add("plus/PiPriPt", "PiPri;#it{p}_{T} (GeV/c)", {HistType::kTH1F, {{500, 0, 5}}});
    registryPri.add("plus/KaPriPt", "KaPri;#it{p}_{T} (GeV/c)", {HistType::kTH1F, {{500, 0, 5}}});
    registryPri.add("plus/DePriPt", "DePri;#it{p}_{T} (GeV/c)", {HistType::kTH1F, {{500, 0, 5}}});
    registryPri.add("plus/AllPriPt", "AllPri;#it{p}_{T} (GeV/c)", {HistType::kTH1F, {{500, 0, 5}}});

    registryPri.add("minus/PiPriPt", "PiPri;#it{p}_{T} (GeV/c)", {HistType::kTH1F, {{500, 0, 5}}});
    registryPri.add("minus/KaPriPt", "KaPri;#it{p}_{T} (GeV/c)", {HistType::kTH1F, {{500, 0, 5}}});
    registryPri.add("minus/DePriPt", "DePri;#it{p}_{T} (GeV/c)", {HistType::kTH1F, {{500, 0, 5}}});
    registryPri.add("minus/AllPriPt", "AllPri;#it{p}_{T} (GeV/c)", {HistType::kTH1F, {{500, 0, 5}}});

    registryPri.add("plus/TOFmatchingAll", ";#it{p}_{T} (GeV/c)", {HistType::kTH1F, {{500, 0, 5}}});
    registryPri.add("minus/TOFmatchingAll", ";#it{p}_{T} (GeV/c)", {HistType::kTH1F, {{500, 0, 5}}});

    // Pri our tracking cuts only
    registryPriCuts.add("plus/PiPriPt", "PiPri;#it{p}_{T} (GeV/c)", {HistType::kTH1F, {{500, 0, 5}}});
    registryPriCuts.add("plus/KaPriPt", "KaPri;#it{p}_{T} (GeV/c)", {HistType::kTH1F, {{500, 0, 5}}});
    registryPriCuts.add("plus/DePriPt", "DePri;#it{p}_{T} (GeV/c)", {HistType::kTH1F, {{500, 0, 5}}});

    registryPriCuts.add("minus/PiPriPt", "PiPri;#it{p}_{T} (GeV/c)", {HistType::kTH1F, {{500, 0, 5}}});
    registryPriCuts.add("minus/KaPriPt", "KaPri;#it{p}_{T} (GeV/c)", {HistType::kTH1F, {{500, 0, 5}}});
    registryPriCuts.add("minus/DePriPt", "DePri;#it{p}_{T} (GeV/c)", {HistType::kTH1F, {{500, 0, 5}}});

    // MC truth
    registryMCtruth.add("plus/MCtruthPi", "MC truth pions;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 5}, {400, -1.0, 1.0}}});
    registryMCtruth.add("plus/MCtruthKa", "MC truth kaons;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 5}, {400, -1.0, 1.0}}});
    registryMCtruth.add("plus/MCtruthDe", "MC truth deuterons;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 5}, {400, -1.0, 1.0}}});

    registryMCtruth.add("minus/MCtruthPi", "MC truth pions;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 5}, {400, -1.0, 1.0}}});
    registryMCtruth.add("minus/MCtruthKa", "MC truth kaons;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 5}, {400, -1.0, 1.0}}});
    registryMCtruth.add("minus/MCtruthDe", "MC truth deuterons;#it{p}_{T} (GeV/c); #eta", {HistType::kTH2F, {{500, 0, 5}, {400, -1.0, 1.0}}});

    registryMCtruth.add("plus/MCtruthPiPt", "MC truth pions;#it{p}_{T} (GeV/c)", {HistType::kTH1F, {{500, 0, 5}}});
    registryMCtruth.add("plus/MCtruthKaPt", "MC truth kaons;#it{p}_{T} (GeV/c)", {HistType::kTH1F, {{500, 0, 5}}});
    registryMCtruth.add("plus/MCtruthDePt", "MC truth deuterons;#it{p}_{T} (GeV/c)", {HistType::kTH1F, {{500, 0, 5}}});
    registryMCtruth.add("plus/MCtruthAllPt", "MC truth all;#it{p}_{T} (GeV/c)", {HistType::kTH1F, {{500, 0, 5}}});

    registryMCtruth.add("minus/MCtruthPiPt", "MC truth pions;#it{p}_{T} (GeV/c)", {HistType::kTH1F, {{500, 0, 5}}});
    registryMCtruth.add("minus/MCtruthKaPt", "MC truth kaons;#it{p}_{T} (GeV/c)", {HistType::kTH1F, {{500, 0, 5}}});
    registryMCtruth.add("minus/MCtruthDePt", "MC truth deuterons;#it{p}_{T} (GeV/c)", {HistType::kTH1F, {{500, 0, 5}}});
    registryMCtruth.add("minus/MCtruthAllPt", "MC truth all;#it{p}_{T} (GeV/c)", {HistType::kTH1F, {{500, 0, 5}}});
  }

  bool IsNSigmaAccept(float nsigmaTPC, float nsigmaTOF, float mom)
  {
    if (ConfPIDnoTOF) {
      if (TMath::Abs(nsigmaTPC) < pidnSigmaCut) {
        return true;
      } else {
        return false;
      }
    } else {
      if (mom > tofPtCut) {
        if (TMath::Hypot(nsigmaTOF, nsigmaTPC) < pidnSigmaCut)
          return true;
      } else {
        if (TMath::Abs(nsigmaTPC) < pidnSigmaCut)
          return true;
      }
      return false;
    }
  }

  Filter trackFilter = nabs(aod::track::eta) < cfgCutEta && requireGlobalTrackInFilter() && aod::track::pt > cfgPtLow&& aod::track::pt < cfgPtHigh;

  void processReco(const CollisionsEvSel::iterator& collision,
                   soa::Filtered<soa::Join<TracksPID, aod::TrackSelection>> const& tracks /*, aod::BCsWithTimestamps const&*/)
  {
    if (ConfIsRun3) {
      if (!collision.sel8())
        return;
    }

    // Loop over tracks
    for (auto& track : tracks) {
      // Tracks are already filtered by the pre-filters

      registryGlobal.fill(HIST("crossedRows"), track.tpcNClsCrossedRows());
      registryGlobal.fill(HIST("RowsOverClustersTPC"), track.tpcCrossedRowsOverFindableCls());
      registryGlobal.fill(HIST("chi2TPC"), track.tpcChi2NCl());
      registryGlobal.fill(HIST("chi2ITS"), track.itsChi2NCl());
      registryGlobal.fill(HIST("dcaZ"), track.dcaZ());
      registryGlobal.fill(HIST("dcaXY"), track.dcaXY());
      registryGlobal.fill(HIST("clustersITS"), track.itsNCls());
      registryGlobal.fill(HIST("pt"), track.pt());
      registryGlobal.fill(HIST("eta"), track.eta());

      // no PID histograms
      if (track.sign() > 0) {
        registryQAtrack.fill(HIST("after/all/plus/etaphi"), track.eta(), track.phi());
        registryQAtrack.fill(HIST("after/all/plus/pt"), track.pt());
        registryQAtrack.fill(HIST("after/all/plus/eta"), track.eta());
        registryQAtrack.fill(HIST("after/all/plus/phi"), track.phi());
        registryQAtrack.fill(HIST("after/all/plus/DCAxy"), track.pt(), track.dcaXY());
        registryQAtrack.fill(HIST("after/all/plus/DCAz"), track.pt(), track.dcaZ());
      }
      if (track.sign() < 0) {
        registryQAtrack.fill(HIST("after/all/minus/etaphi"), track.eta(), track.phi());
        registryQAtrack.fill(HIST("after/all/minus/pt"), track.pt());
        registryQAtrack.fill(HIST("after/all/minus/eta"), track.eta());
        registryQAtrack.fill(HIST("after/all/minus/phi"), track.phi());
        registryQAtrack.fill(HIST("after/all/minus/DCAxy"), track.pt(), track.dcaXY());
        registryQAtrack.fill(HIST("after/all/minus/DCAz"), track.pt(), track.dcaZ());
      }

      // Add PID selection criteria here
      if (IsNSigmaAccept(std::abs(track.tpcNSigmaPi()), std::abs(track.tofNSigmaPi()), track.pt())) {
        if (track.sign() > 0) {
          registryPID.fill(HIST("pid/pion/plus/TOF_TPC_Map"), track.tofNSigmaPi(), track.tpcNSigmaPi());
          registryPID.fill(HIST("pid/pion/plus/TOF_Nsigma"), track.pt(), track.tofNSigmaPi());
          registryPID.fill(HIST("pid/pion/plus/TPC_Nsigma"), track.pt(), track.tpcNSigmaPi());
          registryPID.fill(HIST("pid/pion/plus/TPC_dEdx"), track.pt(), track.tpcSignal());

          registryQAtrack.fill(HIST("after/pion/plus/etaphi"), track.eta(), track.phi());
          registryQAtrack.fill(HIST("after/pion/plus/pt"), track.pt());
          registryQAtrack.fill(HIST("after/pion/plus/eta"), track.eta());
          registryQAtrack.fill(HIST("after/pion/plus/phi"), track.phi());
          registryQAtrack.fill(HIST("after/pion/plus/DCAxy"), track.pt(), track.dcaXY());
          registryQAtrack.fill(HIST("after/pion/plus/DCAz"), track.pt(), track.dcaZ());
        }
        if (track.sign() < 0) {
          registryPID.fill(HIST("pid/pion/minus/TOF_TPC_Map"), track.tofNSigmaPi(), track.tpcNSigmaPi());
          registryPID.fill(HIST("pid/pion/minus/TOF_Nsigma"), track.pt(), track.tofNSigmaPi());
          registryPID.fill(HIST("pid/pion/minus/TPC_Nsigma"), track.pt(), track.tpcNSigmaPi());
          registryPID.fill(HIST("pid/pion/minus/TPC_dEdx"), track.pt(), track.tpcSignal());

          registryQAtrack.fill(HIST("after/pion/minus/etaphi"), track.eta(), track.phi());
          registryQAtrack.fill(HIST("after/pion/minus/pt"), track.pt());
          registryQAtrack.fill(HIST("after/pion/minus/eta"), track.eta());
          registryQAtrack.fill(HIST("after/pion/minus/phi"), track.phi());
          registryQAtrack.fill(HIST("after/pion/minus/DCAxy"), track.pt(), track.dcaXY());
          registryQAtrack.fill(HIST("after/pion/minus/DCAz"), track.pt(), track.dcaZ());
        }
      }
      if (IsNSigmaAccept(std::abs(track.tpcNSigmaKa()), std::abs(track.tofNSigmaKa()), track.pt())) {
        if (track.sign() > 0) {
          registryPID.fill(HIST("pid/kaon/plus/TOF_TPC_Map"), track.tofNSigmaKa(), track.tpcNSigmaKa());
          registryPID.fill(HIST("pid/kaon/plus/TOF_Nsigma"), track.pt(), track.tofNSigmaKa());
          registryPID.fill(HIST("pid/kaon/plus/TPC_Nsigma"), track.pt(), track.tpcNSigmaKa());
          registryPID.fill(HIST("pid/kaon/plus/TPC_dEdx"), track.pt(), track.tpcSignal());

          registryQAtrack.fill(HIST("after/kaon/plus/etaphi"), track.eta(), track.phi());
          registryQAtrack.fill(HIST("after/kaon/plus/pt"), track.pt());
          registryQAtrack.fill(HIST("after/kaon/plus/eta"), track.eta());
          registryQAtrack.fill(HIST("after/kaon/plus/phi"), track.phi());
          registryQAtrack.fill(HIST("after/kaon/plus/DCAxy"), track.pt(), track.dcaXY());
          registryQAtrack.fill(HIST("after/kaon/plus/DCAz"), track.pt(), track.dcaZ());
        }
        if (track.sign() < 0) {
          registryPID.fill(HIST("pid/kaon/minus/TOF_TPC_Map"), track.tofNSigmaKa(), track.tpcNSigmaKa());
          registryPID.fill(HIST("pid/kaon/minus/TOF_Nsigma"), track.pt(), track.tofNSigmaKa());
          registryPID.fill(HIST("pid/kaon/minus/TPC_Nsigma"), track.pt(), track.tpcNSigmaKa());
          registryPID.fill(HIST("pid/kaon/minus/TPC_dEdx"), track.pt(), track.tpcSignal());

          registryQAtrack.fill(HIST("after/kaon/minus/etaphi"), track.eta(), track.phi());
          registryQAtrack.fill(HIST("after/kaon/minus/pt"), track.pt());
          registryQAtrack.fill(HIST("after/kaon/minus/eta"), track.eta());
          registryQAtrack.fill(HIST("after/kaon/minus/phi"), track.phi());
          registryQAtrack.fill(HIST("after/kaon/minus/DCAxy"), track.pt(), track.dcaXY());
          registryQAtrack.fill(HIST("after/kaon/minus/DCAz"), track.pt(), track.dcaZ());
        }
      }
      if (IsNSigmaAccept(std::abs(track.tpcNSigmaDe()), std::abs(track.tofNSigmaDe()), track.pt())) {
        if (track.sign() > 0) {
          registryPID.fill(HIST("pid/deuteron/plus/TOF_TPC_Map"), track.tofNSigmaDe(), track.tpcNSigmaDe());
          registryPID.fill(HIST("pid/deuteron/plus/TOF_Nsigma"), track.pt(), track.tofNSigmaDe());
          registryPID.fill(HIST("pid/deuteron/plus/TPC_Nsigma"), track.pt(), track.tpcNSigmaDe());
          registryPID.fill(HIST("pid/deuteron/plus/TPC_dEdx"), track.pt(), track.tpcSignal());

          registryQAtrack.fill(HIST("after/deuteron/plus/etaphi"), track.eta(), track.phi());
          registryQAtrack.fill(HIST("after/deuteron/plus/pt"), track.pt());
          registryQAtrack.fill(HIST("after/deuteron/plus/eta"), track.eta());
          registryQAtrack.fill(HIST("after/deuteron/plus/phi"), track.phi());
          registryQAtrack.fill(HIST("after/deuteron/plus/DCAxy"), track.pt(), track.dcaXY());
          registryQAtrack.fill(HIST("after/deuteron/plus/DCAz"), track.pt(), track.dcaZ());
        }
        if (track.sign() < 0) {
          registryPID.fill(HIST("pid/deuteron/minus/TOF_TPC_Map"), track.tofNSigmaDe(), track.tpcNSigmaDe());
          registryPID.fill(HIST("pid/deuteron/minus/TOF_Nsigma"), track.pt(), track.tofNSigmaDe());
          registryPID.fill(HIST("pid/deuteron/minus/TPC_Nsigma"), track.pt(), track.tpcNSigmaDe());
          registryPID.fill(HIST("pid/deuteron/minus/TPC_dEdx"), track.pt(), track.tpcSignal());

          registryQAtrack.fill(HIST("after/deuteron/minus/etaphi"), track.eta(), track.phi());
          registryQAtrack.fill(HIST("after/deuteron/minus/pt"), track.pt());
          registryQAtrack.fill(HIST("after/deuteron/minus/eta"), track.eta());
          registryQAtrack.fill(HIST("after/deuteron/minus/phi"), track.phi());
          registryQAtrack.fill(HIST("after/deuteron/minus/DCAxy"), track.pt(), track.dcaXY());
          registryQAtrack.fill(HIST("after/deuteron/minus/DCAz"), track.pt(), track.dcaZ());
        }
      }
    }
  }
  PROCESS_SWITCH(femtoWorldEficiencyTaskDe, processReco, "Process reconstructed data", true);

  using BigTracksMC = soa::Join<TracksPID, aod::McTrackLabels, aod::TrackSelection>;
  Preslice<BigTracksMC> perCollisionID = aod::track::collisionId;
  void processMCTruth(aod::McCollision const&, soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions>> const& collisions, aod::McParticles const& mcparticles, soa::Filtered<BigTracksMC> const& tracks)
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

          int pdg = mcParticle.pdgCode();
          if (pdg == 1000010020) {
            pdg = 999;
          }

          registryPDG.fill(HIST("plus/PDGDe"), track.pt(), pdg);
        }
        if (track.sign() < 0) {
          if (IsNSigmaAccept(std::abs(track.tpcNSigmaPi()), std::abs(track.tofNSigmaPi()), track.pt())) {
            registryPDG.fill(HIST("minus/PDGPi"), track.pt(), mcParticle.pdgCode());
          }
          if (IsNSigmaAccept(std::abs(track.tpcNSigmaKa()), std::abs(track.tofNSigmaKa()), track.pt())) {
            registryPDG.fill(HIST("minus/PDGKa"), track.pt(), mcParticle.pdgCode());
          }
          int pdg = mcParticle.pdgCode();
          if (pdg == -1000010020) {
            pdg = -999;
          }
          registryPDG.fill(HIST("minus/PDGDe"), track.pt(), pdg);
        }

        if (mcParticle.isPhysicalPrimary()) {
          if (track.sign() > 0) {
            // PID only
            registryPri.fill(HIST("plus/AllPri"), track.pt(), track.eta());
            registryPri.fill(HIST("plus/AllPriPt"), mcParticle.pt());
            // histogram pt TOF matching
            if (track.hasTOF()) {
              registryPri.fill(HIST("plus/TOFmatchingAll"), mcParticle.pt());
            }
            if (IsNSigmaAccept(std::abs(track.tpcNSigmaPi()), std::abs(track.tofNSigmaPi()), track.pt()) && mcParticle.pdgCode() == 211) {
              registryPri.fill(HIST("plus/PiPri"), track.pt(), track.eta());
              registryPri.fill(HIST("plus/PiPriPt"), mcParticle.pt());
            }
            if (IsNSigmaAccept(std::abs(track.tpcNSigmaKa()), std::abs(track.tofNSigmaKa()), track.pt()) && mcParticle.pdgCode() == 321) {
              registryPri.fill(HIST("plus/KaPri"), track.pt(), track.eta());
              registryPri.fill(HIST("plus/KaPriPt"), mcParticle.pt());
            }
            if (IsNSigmaAccept(std::abs(track.tpcNSigmaDe()), std::abs(track.tofNSigmaDe()), track.pt()) && mcParticle.pdgCode() == 1000010020) {
              registryPri.fill(HIST("plus/DePri"), track.pt(), track.eta());
              registryPri.fill(HIST("plus/DePriPt"), mcParticle.pt());
            }
            // tracking efficiency only
            if (mcParticle.pdgCode() == 211) {
              registryPriCuts.fill(HIST("plus/PiPriPt"), mcParticle.pt());
            }
            if (mcParticle.pdgCode() == 321) {
              registryPriCuts.fill(HIST("plus/KaPriPt"), mcParticle.pt());
            }
            if (mcParticle.pdgCode() == 1000010020) {
              registryPriCuts.fill(HIST("plus/DePriPt"), mcParticle.pt());
            }
          }
          if (track.sign() < 0) {
            // PID only
            registryPri.fill(HIST("minus/AllPri"), track.pt(), track.eta());
            registryPri.fill(HIST("minus/AllPriPt"), mcParticle.pt());
            // histogram pt TOF matching
            if (track.hasTOF()) {
              registryPri.fill(HIST("minus/TOFmatchingAll"), mcParticle.pt());
            }
            if (IsNSigmaAccept(std::abs(track.tpcNSigmaPi()), std::abs(track.tofNSigmaPi()), track.pt()) && mcParticle.pdgCode() == -211) {
              registryPri.fill(HIST("minus/PiPri"), track.pt(), track.eta());
              registryPri.fill(HIST("minus/PiPriPt"), mcParticle.pt());
            }
            if (IsNSigmaAccept(std::abs(track.tpcNSigmaKa()), std::abs(track.tofNSigmaKa()), track.pt()) && mcParticle.pdgCode() == -321) {
              registryPri.fill(HIST("minus/KaPri"), track.pt(), track.eta());
              registryPri.fill(HIST("minus/KaPriPt"), mcParticle.pt());
            }
            if (IsNSigmaAccept(std::abs(track.tpcNSigmaDe()), std::abs(track.tofNSigmaDe()), track.pt()) && mcParticle.pdgCode() == -1000010020) {
              registryPri.fill(HIST("minus/DePri"), track.pt(), track.eta());
              registryPri.fill(HIST("minus/DePriPt"), mcParticle.pt());
            }
            // tracking efficiency only
            if (mcParticle.pdgCode() == -211) {
              registryPriCuts.fill(HIST("minus/PiPriPt"), mcParticle.pt());
            }
            if (mcParticle.pdgCode() == -321) {
              registryPriCuts.fill(HIST("minus/KaPriPt"), mcParticle.pt());
            }
            if (mcParticle.pdgCode() == -1000010020) {
              registryPriCuts.fill(HIST("minus/DePriPt"), mcParticle.pt());
            }
          }
        }
      }
    }
    // loop over MC particles
    for (auto& mcparticle : mcparticles) {
      if (!mcparticle.isPhysicalPrimary() || TMath::Abs(mcparticle.eta()) > cfgCutEta || mcparticle.pt() < cfgPtLow || mcparticle.pt() > cfgPtHigh)
        continue;

      const auto& pdgParticle = pdgDB->GetParticle(mcparticle.pdgCode());
      if (!pdgParticle) {
        continue;
      }

      if (pdgParticle->Charge() > 0) {
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
      if (mcparticle.pdgCode() == 1000010020) {
        registryMCtruth.fill(HIST("plus/MCtruthDe"), mcparticle.pt(), mcparticle.eta());
        registryMCtruth.fill(HIST("plus/MCtruthDePt"), mcparticle.pt());
      }

      if (pdgParticle->Charge() < 0) {
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
      if (mcparticle.pdgCode() == -1000010020) {
        registryMCtruth.fill(HIST("minus/MCtruthDe"), mcparticle.pt(), mcparticle.eta());
        registryMCtruth.fill(HIST("minus/MCtruthDePt"), mcparticle.pt());
      }
    }
  }
  PROCESS_SWITCH(femtoWorldEficiencyTaskDe, processMCTruth, "Process MC truth data", true);

}; // end of spectra task
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<femtoWorldEficiencyTaskDe>(cfgc, TaskName{"femtoWorldEficiencyTaskDe"}),
  };
}
