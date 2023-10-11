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

// \author
// Alice Caluisi   -   alice.caluisi@cern.ch
// \since July 2023

//
// Task producing jet tracking qa histograms
//
#include <iostream>

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "Framework/ASoA.h"

#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/TrackJetQa.h"
#include "PWGJE/TableProducer/jetfinder.h"

using namespace o2;
using namespace o2::track;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct TrackJetQa {
  HistogramRegistry histos{"JetQAHistograms"};
  // Configurable<int> selectedTracks{"select", 1, "Choice of track selection. 0 = no selection, 1 = globalTracks"}; --not in use
  Configurable<bool> enable{"selectTrack", true, "false = disable track selection, true = enable track selection"};
  Configurable<int> nBins{"nBins", 200, "N bins in histos"};

  Configurable<double> ValVtx{"ValVtx", 10, "Value of the vertex position"};
  Configurable<float> ValCutEta{"ValCutEta", 0.8f, "Eta range for tracks"};
  Configurable<float> ValCutY{"ValCutY", 0.5f, "Y range for tracks"};

  // Custom track cuts for the cut variation study
  TrackSelection customTrackCuts;
  Configurable<int> itsPattern{"itsPattern", 1, "0 = Run3ITSibAny, 1 = Run3ITSallAny, 2 = Run3ITSall7Layers, 3 = Run3ITSibTwo"};
  Configurable<bool> requireITS{"requireITS", true, "Additional cut on the ITS requirement"};
  Configurable<bool> requireTPC{"requireTPC", true, "Additional cut on the TPC requirement"};
  Configurable<bool> requireGoldenChi2{"requireGoldenChi2", true, "Additional cut on the GoldenChi2"};
  Configurable<float> minNCrossedRowsTPC{"minNCrossedRowsTPC", 60.f, "Additional cut on the minimum number of crossed rows in the TPC"};
  Configurable<float> minNCrossedRowsOverFindableClustersTPC{"minNCrossedRowsOverFindableClustersTPC", 0.7f, "Additional cut on the minimum value of the ratio between crossed rows and findable clusters in the TPC"};
  Configurable<float> maxChi2PerClusterTPC{"maxChi2PerClusterTPC", 7.f, "Additional cut on the maximum value of the chi2 per cluster in the TPC"};
  Configurable<float> maxChi2PerClusterITS{"maxChi2PerClusterITS", 36.f, "Additional cut on the maximum value of the chi2 per cluster in the ITS"};
  Configurable<float> maxDcaXYFactor{"maxDcaXYFactor", 1.f, "Additional cut on the maximum value of the DCA xy (multiplicative factor)"};
  Configurable<float> maxDcaZ{"maxDcaZ", 3.f, "Additional cut on the maximum value of the DCA z"};
  Configurable<float> minTPCNClsFound{"minTPCNClsFound", 0.f, "Additional cut on the minimum value of the number of found clusters in the TPC"};

  void init(o2::framework::InitContext&)
  {
    // Custom track cuts - to adjust cuts for the tree analysis via configurables
    LOG(info) << "Using custom track cuts from values:";
    LOG(info) << "\trequireITS=" << requireITS.value;
    LOG(info) << "\trequireTPC=" << requireTPC.value;
    LOG(info) << "\trequireGoldenChi2=" << requireGoldenChi2.value;
    LOG(info) << "\tmaxChi2PerClusterTPC=" << maxChi2PerClusterTPC.value;
    LOG(info) << "\tminNCrossedRowsTPC=" << minNCrossedRowsTPC.value;
    LOG(info) << "\tminTPCNClsFound=" << minTPCNClsFound.value;
    LOG(info) << "\tmaxChi2PerClusterITS=" << maxChi2PerClusterITS.value;
    LOG(info) << "\tmaxDcaZ=" << maxDcaZ.value;

    customTrackCuts = getGlobalTrackSelectionRun3ITSMatch(itsPattern.value);
    LOG(info) << "Customizing track cuts:";
    customTrackCuts.SetRequireITSRefit(requireITS.value);
    customTrackCuts.SetRequireTPCRefit(requireTPC.value);
    customTrackCuts.SetRequireGoldenChi2(requireGoldenChi2.value);
    customTrackCuts.SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC.value);
    customTrackCuts.SetMaxChi2PerClusterITS(maxChi2PerClusterITS.value);
    customTrackCuts.SetMinNCrossedRowsTPC(minNCrossedRowsTPC.value);
    customTrackCuts.SetMinNClustersTPC(minTPCNClsFound.value);
    customTrackCuts.SetMinNCrossedRowsOverFindableClustersTPC(minNCrossedRowsOverFindableClustersTPC.value);
    customTrackCuts.SetMaxDcaXYPtDep([](float pt) { return 10.f; }); // No DCAxy cut will be used, this is done via the member function of the task
    customTrackCuts.SetMaxDcaZ(maxDcaZ.value);
    customTrackCuts.print();
    // kinetic histograms
    histos.add("Kine/pt", "#it{p}_{T};#it{p}_{T} [GeV/c];number of entries", HistType::kTH1F, {{nBins, 0, 200}});
    histos.add("Kine/pt_TRD", "#it{p}_{T} if track has a TRD match;#it{p}_{T} [GeV/c];number of entries", HistType::kTH1F, {{nBins, 0, 200}});
    histos.add("Kine/eta", "#eta;#it{p}_{T} [GeV/c];#eta", {HistType::kTH2F, {{nBins, 0, 200}, {180, -0.9, 0.9}}});
    histos.add("Kine/phi", "#phi;#it{p}_{T} [GeV/c];#phi [rad]", {HistType::kTH2F, {{nBins, 0, 200}, {180, 0., 2 * M_PI}}});
    histos.add("Kine/etaVSphi", "#eta VS phi;#eta;#phi [rad]", {HistType::kTH2F, {{180, -0.9, 0.9}, {180, 0., 2 * M_PI}}});
    histos.add("Kine/etaVSphi0", "#eta VS phi in pT range {0, 1};#eta;#phi [rad]", {HistType::kTH2F, {{180, -0.9, 0.9}, {180, 0., 2 * M_PI}}});
    histos.add("Kine/etaVSphi1", "#eta VS phi in pT range {1, 2};#eta;#phi [rad]", {HistType::kTH2F, {{180, -0.9, 0.9}, {180, 0., 2 * M_PI}}});
    histos.add("Kine/etaVSphi2", "#eta VS phi in pT range {2, 5};#eta;#phi [rad]", {HistType::kTH2F, {{180, -0.9, 0.9}, {180, 0., 2 * M_PI}}});
    histos.add("Kine/etaVSphi3", "#eta VS phi in pT range {5, 10};#eta;#phi [rad]", {HistType::kTH2F, {{180, -0.9, 0.9}, {180, 0., 2 * M_PI}}});
    histos.add("Kine/etaVSphi4", "#eta VS phi in pT range {10, 20};#eta;#phi [rad]", {HistType::kTH2F, {{180, -0.9, 0.9}, {180, 0., 2 * M_PI}}});
    histos.add("Kine/etaVSphi5", "#eta VS phi in pT range {20, 50};#eta;#phi [rad]", {HistType::kTH2F, {{180, -0.9, 0.9}, {180, 0., 2 * M_PI}}});
    histos.add("Kine/etaVSphi6", "#eta VS phi in pT range {50, 200};#eta;#phi [rad]", {HistType::kTH2F, {{180, -0.9, 0.9}, {180, 0., 2 * M_PI}}});

    // track parameter histograms
    histos.add("TrackPar/x", "track #it{x} position at dca in local coordinate system;#it{p}_{T} [GeV/c];#it{x} [cm]", {HistType::kTH2F, {{nBins, 0, 200}, {200, -0.36, 0.36}}});
    histos.add("TrackPar/y", "track #it{y} position at dca in local coordinate system;#it{p}_{T} [GeV/c];#it{y} [cm]", {HistType::kTH2F, {{nBins, 0, 200}, {200, -0.5, 0.5}}});
    histos.add("TrackPar/z", "track #it{z} position at dca in local coordinate system;#it{p}_{T} [GeV/c];#it{z} [cm]", {HistType::kTH2F, {{nBins, 0, 200}, {200, -11., 11.}}});
    histos.add("TrackPar/alpha", "rotation angle of local wrt. global coordinate system;#it{p}_{T} [GeV/c];#alpha [rad]", {HistType::kTH2F, {{nBins, 0, 200}, {36, -M_PI, M_PI}}});
    histos.add("TrackPar/signed1Pt", "track signed 1/#it{p}_{T};#it{p}_{T} [GeV/c];#it{q}/#it{p}_{T}", {HistType::kTH2F, {{nBins, 0, 200}, {200, -8, 8}}});
    histos.add("TrackPar/snp", "sinus of track momentum azimuthal angle;#it{p}_{T} [GeV/c];snp", {HistType::kTH2F, {{nBins, 0, 200}, {11, -0.1, 0.1}}});
    histos.add("TrackPar/tgl", "tangent of the track momentum dip angle;#it{p}_{T} [GeV/c];tgl;", {HistType::kTH2F, {{nBins, 0, 200}, {200, -1., 1.}}});
    histos.add("TrackPar/flags", "track flag;#it{p}_{T} [GeV/c];flag bit", {HistType::kTH2F, {{nBins, 0, 200}, {64, -0.5, 63.5}}});
    histos.add("TrackPar/dcaXY", "distance of closest approach in #it{xy} plane;#it{p}_{T} [GeV/c];#it{dcaXY} [cm];", {HistType::kTH2F, {{nBins, 0, 200}, {200, -0.15, 0.15}}});
    histos.add("TrackPar/dcaZ", "distance of closest approach in #it{z};#it{p}_{T} [GeV/c];#it{dcaZ} [cm];", {HistType::kTH2F, {{nBins, 0, 200}, {200, -0.15, 0.15}}});
    histos.add("TrackPar/length", "track length in cm;#it{p}_{T} [GeV/c];#it{Length} [cm];", {HistType::kTH2F, {{nBins, 0, 200}, {200, 0, 1000}}});
    histos.add("TrackPar/Sigma1Pt", "uncertainty over #it{p}_{T};#it{p}_{T} [GeV/c];#it{p}_{T}*#it{sigma1}{p}_{T};", {HistType::kTH2F, {{nBins, 0, 200}, {100, 0, 1}}});
    histos.add("TrackPar/Sigma1Pt_Layer1", "uncertainty over #it{p}_{T} with only 1st ITS layer active;#it{p}_{T} [GeV/c];#it{p}_{T}*#it{sigma1}{p}_{T};", {HistType::kTH2F, {{nBins, 0, 200}, {100, 0, 1}}});
    histos.add("TrackPar/Sigma1Pt_Layer2", "uncertainty over #it{p}_{T} with only 2nd ITS layer active;#it{p}_{T} [GeV/c];#it{p}_{T}*#it{sigma1}{p}_{T};", {HistType::kTH2F, {{nBins, 0, 200}, {100, 0, 1}}});
    histos.add("TrackPar/Sigma1Pt_Layers12", "uncertainty over #it{p}_{T} with only 1st and 2nd ITS layers active;#it{p}_{T} [GeV/c];#it{p}_{T}*#it{sigma1}{p}_{T};", {HistType::kTH2F, {{nBins, 0, 200}, {100, 0, 1}}});
    histos.add("TrackPar/Sigma1Pt_Layer4", "uncertainty over #it{p}_{T} with only 4th ITS layer active;#it{p}_{T} [GeV/c];#it{p}_{T}*#it{sigma1}{p}_{T};", {HistType::kTH2F, {{nBins, 0, 200}, {100, 0, 1}}});
    histos.add("TrackPar/Sigma1Pt_Layer5", "uncertainty over #it{p}_{T} with only 5th ITS layer active;#it{p}_{T} [GeV/c];#it{p}_{T}*#it{sigma1}{p}_{T};", {HistType::kTH2F, {{nBins, 0, 200}, {100, 0, 1}}});
    histos.add("TrackPar/Sigma1Pt_Layer6", "uncertainty over #it{p}_{T} with only 6th ITS layer active;#it{p}_{T} [GeV/c];#it{p}_{T}*#it{sigma1}{p}_{T};", {HistType::kTH2F, {{nBins, 0, 200}, {100, 0, 1}}});
    histos.add("TrackPar/Sigma1Pt_Layers45", "uncertainty over #it{p}_{T} with only 4th and 5th ITS layers active;#it{p}_{T} [GeV/c];#it{p}_{T}*#it{sigma1}{p}_{T};", {HistType::kTH2F, {{nBins, 0, 200}, {100, 0, 1}}});
    histos.add("TrackPar/Sigma1Pt_Layers56", "uncertainty over #it{p}_{T} with only 5th and 6th ITS layers active;#it{p}_{T} [GeV/c];#it{p}_{T}*#it{sigma1}{p}_{T};", {HistType::kTH2F, {{nBins, 0, 200}, {100, 0, 1}}});
    histos.add("TrackPar/Sigma1Pt_Layers46", "uncertainty over #it{p}_{T} with only 4th and 6th ITS layers active;#it{p}_{T} [GeV/c];#it{p}_{T}*#it{sigma1}{p}_{T};", {HistType::kTH2F, {{nBins, 0, 200}, {100, 0, 1}}});
    histos.add("TrackPar/Sigma1Pt_Layers456", "uncertainty over #it{p}_{T} with only 4th, 5th and 6th ITS layers active;#it{p}_{T} [GeV/c];#it{p}_{T}*#it{sigma1}{p}_{T};", {HistType::kTH2F, {{nBins, 0, 200}, {100, 0, 1}}});

    // event property histograms
    histos.add("EventProp/collisionVtxZ", "Collsion Vertex Z;#it{Vtx}_{z} [cm];number of entries", HistType::kTH1F, {{nBins, -20, 20}});
    histos.add("EventProp/collisionVtxZnoSel", "Collsion Vertex Z without event selection;#it{Vtx}_{z} [cm];number of entries", HistType::kTH1F, {{nBins, -20, 20}});
    histos.add("EventProp/collisionVtxZSel8", "Collsion Vertex Z with event selection;#it{Vtx}_{z} [cm];number of entries", HistType::kTH1F, {{nBins, -20, 20}});

    // ITS histograms
    histos.add("ITS/itsNCls", "number of found ITS clusters;#it{p}_{T} [GeV/c];# clusters ITS", {HistType::kTH2F, {{nBins, 0, 200}, {8, -0.5, 7.5}}});
    histos.add("ITS/itsChi2NCl", "chi2 per ITS cluster;#it{p}_{T} [GeV/c];chi2 / cluster ITS", {HistType::kTH2F, {{nBins, 0, 200}, {100, 0, 40}}});
    histos.add("ITS/itsHits", "hitmap ITS;#it{p}_{T} [GeV/c];layer ITS", {HistType::kTH2F, {{nBins, 0, 200}, {7, -0.5, 6.5}}});

    // TPC histograms
    histos.add("TPC/tpcNClsFindable", "number of findable TPC clusters;#it{p}_{T} [GeV/c];# findable clusters TPC", {HistType::kTH2F, {{nBins, 0, 200}, {165, -0.5, 164.5}}});
    histos.add("TPC/tpcNClsFound", "number of found TPC clusters;#it{p}_{T} [GeV/c];# clusters TPC", {HistType::kTH2F, {{nBins, 0, 200}, {165, -0.5, 164.5}}});
    histos.add("TPC/tpcNClsShared", "number of shared TPC clusters;#it{p}_{T} [GeV/c];# shared clusters TPC", {HistType::kTH2F, {{nBins, 0, 200}, {165, -0.5, 164.5}}});
    histos.add("TPC/tpcNClsCrossedRows", "number of crossed TPC rows;#it{p}_{T} [GeV/c];# crossed rows TPC", {HistType::kTH2F, {{nBins, 0, 200}, {165, -0.5, 164.5}}});
    histos.add("TPC/tpcFractionSharedCls", "fraction of shared TPC clusters;#it{p}_{T} [GeV/c];fraction shared clusters TPC", {HistType::kTH2F, {{nBins, 0, 200}, {100, 0., 1.}}});
    histos.add("TPC/tpcCrossedRowsOverFindableCls", "crossed TPC rows over findable clusters;#it{p}_{T} [GeV/c];crossed rows / findable clusters TPC", {HistType::kTH2F, {{nBins, 0, 200}, {120, 0.0, 1.2}}});
    histos.add("TPC/tpcChi2NCl", "chi2 per cluster in TPC;#it{p}_{T} [GeV/c];chi2 / cluster TPC", {HistType::kTH2F, {{nBins, 0, 200}, {100, 0, 10}}});

    histos.print();
  }

  // @Alice, please make these blocks where you fill histograms templates that we can call in both process functions. Thank you !
  void processFull(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                   soa::Join<aod::FullTracks, aod::TracksDCA, aod::TrackSelection, aod::TracksCov> const& tracks)
  {

    for (auto& collision : collisions) {

      // fill event property variables
      histos.fill(HIST("EventProp/collisionVtxZnoSel"), collision.posZ());

      if (!collision.sel8()) {
        return;
      }
      histos.fill(HIST("EventProp/collisionVtxZSel8"), collision.posZ());

      if (fabs(collision.posZ()) > ValVtx) {
        return;
      }
      histos.fill(HIST("EventProp/collisionVtxZ"), collision.posZ());

      Partition<soa::Join<aod::FullTracks, aod::TracksDCA, aod::TrackSelection, aod::TracksCov>> groupedTracks = aod::track::collisionId == collision.globalIndex();
      groupedTracks.bindTable(tracks);
      for (auto& track : groupedTracks) {
        // Track selection
        if (enable && !track.isGlobalTrackWoPtEta()) {
          continue;
        }

        // fill kinematic variables
        histos.fill(HIST("Kine/pt"), track.pt());
        histos.fill(HIST("Kine/eta"), track.pt(), track.eta());
        histos.fill(HIST("Kine/phi"), track.pt(), track.phi());
        histos.fill(HIST("Kine/etaVSphi"), track.eta(), track.phi());
        if (track.hasTRD()) {
          histos.fill(HIST("Kine/pt_TRD"), track.pt());
        }

        //// eta VS phi for different pT ranges
        double pt = track.pt();
        if (pt >= 0.0 && pt < 1.0) {
          histos.fill(HIST("Kine/etaVSphi0"), track.eta(), track.phi());
        } else if (pt >= 1.0 && pt < 2.0) {
          histos.fill(HIST("Kine/etaVSphi1"), track.eta(), track.phi());
        } else if (pt >= 2.0 && pt <= 5.0) {
          histos.fill(HIST("Kine/etaVSphi2"), track.eta(), track.phi());
        } else if (pt >= 5.0 && pt < 10.0) {
          histos.fill(HIST("Kine/etaVSphi3"), track.eta(), track.phi());
        } else if (pt >= 10.0 && pt <= 20.0) {
          histos.fill(HIST("Kine/etaVSphi4"), track.eta(), track.phi());
        } else if (pt >= 20.0 && pt <= 50.0) {
          histos.fill(HIST("Kine/etaVSphi5"), track.eta(), track.phi());
        } else if (pt >= 50.0 && pt <= 200.0) {
          histos.fill(HIST("Kine/etaVSphi6"), track.eta(), track.phi());
        }

        // fill track parameter variables
        histos.fill(HIST("TrackPar/alpha"), track.pt(), track.alpha());
        histos.fill(HIST("TrackPar/x"), track.pt(), track.x());
        histos.fill(HIST("TrackPar/y"), track.pt(), track.y());
        histos.fill(HIST("TrackPar/z"), track.pt(), track.z());
        histos.fill(HIST("TrackPar/signed1Pt"), track.pt(), track.signed1Pt());
        histos.fill(HIST("TrackPar/snp"), track.pt(), track.snp());
        histos.fill(HIST("TrackPar/tgl"), track.pt(), track.tgl());
        for (unsigned int i = 0; i < 64; i++) {
          if (track.flags() & (1 << i)) {
            histos.fill(HIST("TrackPar/flags"), track.pt(), i);
          }
        }
        histos.fill(HIST("TrackPar/dcaXY"), track.pt(), track.dcaXY());
        histos.fill(HIST("TrackPar/dcaZ"), track.pt(), track.dcaZ());
        histos.fill(HIST("TrackPar/length"), track.pt(), track.length());
        histos.fill(HIST("TrackPar/Sigma1Pt"), track.pt(), track.sigma1Pt() * track.pt());

        //// check the uncertainty over pT activating several ITS layers
        bool firstLayerActive = track.itsClusterMap() & (1 << 0);
        bool secondLayerActive = track.itsClusterMap() & (1 << 1);
        bool fourthLayerActive = track.itsClusterMap() & (1 << 3);
        bool fifthLayerActive = track.itsClusterMap() & (1 << 4);
        bool sixthLayerActive = track.itsClusterMap() & (1 << 5);
        if (firstLayerActive) {
          histos.fill(HIST("TrackPar/Sigma1Pt_Layer1"), track.pt(), track.sigma1Pt() * track.pt());
        }
        if (secondLayerActive) {
          histos.fill(HIST("TrackPar/Sigma1Pt_Layer2"), track.pt(), track.sigma1Pt() * track.pt());
        }
        if (firstLayerActive && secondLayerActive) {
          histos.fill(HIST("TrackPar/Sigma1Pt_Layers12"), track.pt(), track.sigma1Pt() * track.pt());
        }
        if (fourthLayerActive) {
          histos.fill(HIST("TrackPar/Sigma1Pt_Layer4"), track.pt(), track.sigma1Pt() * track.pt());
        }
        if (fifthLayerActive) {
          histos.fill(HIST("TrackPar/Sigma1Pt_Layer5"), track.pt(), track.sigma1Pt() * track.pt());
        }
        if (sixthLayerActive) {
          histos.fill(HIST("TrackPar/Sigma1Pt_Layer6"), track.pt(), track.sigma1Pt() * track.pt());
        }
        if (fourthLayerActive && fifthLayerActive) {
          histos.fill(HIST("TrackPar/Sigma1Pt_Layers45"), track.pt(), track.sigma1Pt() * track.pt());
        }
        if (fifthLayerActive && sixthLayerActive) {
          histos.fill(HIST("TrackPar/Sigma1Pt_Layers56"), track.pt(), track.sigma1Pt() * track.pt());
        }
        if (fourthLayerActive && sixthLayerActive) {
          histos.fill(HIST("TrackPar/Sigma1Pt_Layers46"), track.pt(), track.sigma1Pt() * track.pt());
        }
        if (fourthLayerActive && fifthLayerActive && sixthLayerActive) {
          histos.fill(HIST("TrackPar/Sigma1Pt_Layers456"), track.pt(), track.sigma1Pt() * track.pt());
        }

        // fill ITS variables
        histos.fill(HIST("ITS/itsNCls"), track.pt(), track.itsNCls());
        histos.fill(HIST("ITS/itsChi2NCl"), track.pt(), track.itsChi2NCl());
        for (unsigned int i = 0; i < 7; i++) {
          if (track.itsClusterMap() & (1 << i)) {
            histos.fill(HIST("ITS/itsHits"), track.pt(), i);
          }
        }

        // fill TPC variables
        histos.fill(HIST("TPC/tpcNClsFindable"), track.pt(), track.tpcNClsFindable());
        histos.fill(HIST("TPC/tpcNClsFound"), track.pt(), track.tpcNClsFound());
        histos.fill(HIST("TPC/tpcNClsShared"), track.pt(), track.tpcNClsShared());
        histos.fill(HIST("TPC/tpcNClsCrossedRows"), track.pt(), track.tpcNClsCrossedRows());
        histos.fill(HIST("TPC/tpcCrossedRowsOverFindableCls"), track.pt(), track.tpcCrossedRowsOverFindableCls());
        histos.fill(HIST("TPC/tpcFractionSharedCls"), track.pt(), track.tpcFractionSharedCls());
        histos.fill(HIST("TPC/tpcChi2NCl"), track.pt(), track.tpcChi2NCl());
      }
    }
  }
  PROCESS_SWITCH(TrackJetQa, processFull, "Standard data processor", true);

  Preslice<aod::SpTracks> spPerCol = aod::spectra::collisionId;
  SliceCache cacheTrk;
  void processDerived(aod::SpColls const& collisions,
                      aod::SpTracks const& tracks)
  {
    for (const auto& collision : collisions) {
      // fill event property variables - i put here a duplicate of your code above to test. Please make templates to shorten the code !
      histos.fill(HIST("EventProp/collisionVtxZnoSel"), collision.posZ());

      if (!collision.sel8()) {
        return;
      }
      histos.fill(HIST("EventProp/collisionVtxZSel8"), collision.posZ());

      if (fabs(collision.posZ()) > ValVtx) {
        return;
      }
      histos.fill(HIST("EventProp/collisionVtxZ"), collision.posZ());

      // Partition<soa::Join<aod::FullTracks, aod::TracksDCA, aod::TrackSelection, aod::TracksCov>> groupedTracks = aod::track::collisionId == collision.globalIndex();
      // groupedTracks.bindTable(tracks);

      const auto& tracksInCollision = tracks.sliceByCached(aod::spectra::collisionId, collision.globalIndex(), cacheTrk);
      for (const auto& track : tracksInCollision) {
        if (enable && !track.isGlobalTrackWoPtEta()) {
          continue;
        }
        // filling all your histos.... please put them into templates and call the templates in the process functions.
        //  else you really just duplicate code and we want to save line here :)
      }
    }
  } // end of the process function
  PROCESS_SWITCH(TrackJetQa, processDerived, "Derived data processor", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow;
  workflow.push_back(adaptAnalysisTask<TrackJetQa>(cfgc));
  return workflow;
}
