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

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "Framework/ASoA.h"

#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/TableProducer/jetfinder.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

//****************************************************************************************
/**
 * QA histograms for track quantities.
 */
//****************************************************************************************
struct TrackJetQa 
{
  HistogramRegistry histos{"JetQAHistograms"};

  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};

  std::string trackSelection;

  // function that performs track selections on each track
  using JetTracks = soa::Join<aod::Tracks, aod::TrackSelection>;
  template <typename T>
  bool selectTrack(T const& track, std::string trackSelection)
  {
    if (trackSelection == "globalTracks") {
      return track.isGlobalTrackWoPtEta();
    }
  }

  Filter trackFilter = requireGlobalTrackInFilter();
  
  Configurable<int> nBins{"nBins", 200, "N bins in histos"}; 
  
  Configurable<double> ValVtx{"ValVtx", 10, "Value of the vertex position"};

  std::vector<double> ptBinning = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
                                     1.1, 1.2, 1.3, 1.4, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 
                                     6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0};

  void init(o2::framework::InitContext&)
  {
    // set track selections
    trackSelection = static_cast<std::string>(trackSelections);

    // kinetic histograms
    histos.add("Kine/pt", "#it{p}_{T};#it{p}_{T} [GeV/c];number of entries", HistType::kTH1F, {{ptBinning}});
    histos.add("Kine/eta", "#eta;#it{p}_{T} [GeV/c];#eta", {HistType::kTH2F,  {{nBins, 0, 20}, {180, -0.9, 0.9}}});
    histos.add("Kine/phi", "#phi;#it{p}_{T} [GeV/c];#phi [rad]", {HistType::kTH2F,  {{nBins, 0, 20}, {180, 0., 2 * M_PI}}});
    histos.add("Kine/etaVSphi", "#eta VS phi;#eta;#phi [rad]", {HistType::kTH2F,  {{180, -0.9, 0.9}, {180, 0., 2 * M_PI}}});

    // track 2D histograms
    histos.add("TrackPar/x", "track #it{x} position at dca in local coordinate system;#it{p}_{T} [GeV/c];#it{x} [cm]", {HistType::kTH2F, {{nBins, 0, 20}, {200, -0.36, 0.36}}});
    histos.add("TrackPar/y", "track #it{y} position at dca in local coordinate system;#it{p}_{T} [GeV/c];#it{y} [cm]", {HistType::kTH2F, {{nBins, 0, 20}, {200, -0.5, 0.5}}});
    histos.add("TrackPar/z", "track #it{z} position at dca in local coordinate system;#it{p}_{T} [GeV/c];#it{z} [cm]", {HistType::kTH2F, {{nBins, 0, 20}, {200, -11., 11.}}});
    histos.add("TrackPar/alpha", "rotation angle of local wrt. global coordinate system;#it{p}_{T} [GeV/c];#alpha [rad]", {HistType::kTH2F, {{nBins, 0, 20}, {36, -M_PI, M_PI}}});
    histos.add("TrackPar/signed1Pt", "track signed 1/#it{p}_{T};#it{p}_{T} [GeV/c];#it{q}/#it{p}_{T}", {HistType::kTH2F, {{nBins, 0, 20}, {200, -8, 8}}});
    histos.add("TrackPar/snp", "sinus of track momentum azimuthal angle;#it{p}_{T} [GeV/c];snp", {HistType::kTH2F, {{nBins, 0, 20}, {11, -0.1, 0.1}}});
    histos.add("TrackPar/tgl", "tangent of the track momentum dip angle;#it{p}_{T} [GeV/c];tgl;", {HistType::kTH2F, {{nBins, 0, 20}, {200, -1., 1.}}});
    histos.add("TrackPar/flags", "track flag;#it{p}_{T} [GeV/c];flag bit", {HistType::kTH2F, {{nBins, 0, 20}, {64, -0.5, 63.5}}});
    histos.add("TrackPar/dcaXY", "distance of closest approach in #it{xy} plane;#it{p}_{T} [GeV/c];#it{dcaXY} [cm];", {HistType::kTH2F, {{nBins, 0, 20}, {200, -0.15, 0.15}}});
    histos.add("TrackPar/dcaZ", "distance of closest approach in #it{z};#it{p}_{T} [GeV/c];#it{dcaZ} [cm];", {HistType::kTH2F, {{nBins, 0, 20}, {200, -0.15, 0.15}}});
    histos.add("TrackPar/length", "track length in cm;#it{p}_{T} [GeV/c];#it{Length} [cm];", {HistType::kTH2F, {{nBins, 0, 20}, {200, 0, 1000}}});
    histos.add("TrackPar/Sigma1Pt", "uncertainty over #it{p}_{T};#it{p}_{T} [GeV/c];#it{p}_{T}*#it{p}_{T}*#it{sigma1}{p}_{T};", {HistType::kTH2F, {{400, 0, 100}, {100, 0, 1}}});
    
    // event property histograms
    histos.add("EventProp/collisionVtxZ", "Collsion Vertex Z;#it{Vtx}_{z} [cm];number of entries", HistType::kTH1F, {{nBins, -20, 20}});
    histos.add("EventProp/collisionVtxZnoSel", "Collsion Vertex Z without event selection;#it{Vtx}_{z} [cm];number of entries", HistType::kTH1F, {{nBins, -20, 20}});
    histos.add("EventProp/collisionVtxZSel8", "Collsion Vertex Z with event selection;#it{Vtx}_{z} [cm];number of entries", HistType::kTH1F, {{nBins, -20, 20}});

    // ITS histograms
    histos.add("ITS/itsNCls", "number of found ITS clusters;#it{p}_{T} [GeV/c];# clusters ITS", {HistType::kTH2F, {{nBins, 0, 20}, {8, -0.5, 7.5}}});
    histos.add("ITS/itsChi2NCl", "chi2 per ITS cluster;#it{p}_{T} [GeV/c];chi2 / cluster ITS", {HistType::kTH2F, {{nBins, 0, 20}, {100, 0, 40}}});
    histos.add("ITS/itsHits", "hitmap ITS;#it{p}_{T} [GeV/c];layer ITS", {HistType::kTH2F, {{nBins, 0, 20}, {7, -0.5, 6.5}}});

    // TPC histograms
    histos.add("TPC/tpcNClsFindable", "number of findable TPC clusters;#it{p}_{T} [GeV/c];# findable clusters TPC", {HistType::kTH2F, {{nBins, 0, 20}, {165, -0.5, 164.5}}});
    histos.add("TPC/tpcNClsFound", "number of found TPC clusters;#it{p}_{T} [GeV/c];# clusters TPC", {HistType::kTH2F, {{nBins, 0, 20}, {165, -0.5, 164.5}}});
    histos.add("TPC/tpcNClsShared", "number of shared TPC clusters;#it{p}_{T} [GeV/c];# shared clusters TPC", {HistType::kTH2F, {{nBins, 0, 20}, {165, -0.5, 164.5}}});
    histos.add("TPC/tpcNClsCrossedRows", "number of crossed TPC rows;#it{p}_{T} [GeV/c];# crossed rows TPC", {HistType::kTH2F, {{nBins, 0, 20}, {165, -0.5, 164.5}}});
    histos.add("TPC/tpcFractionSharedCls", "fraction of shared TPC clusters;#it{p}_{T} [GeV/c];fraction shared clusters TPC", {HistType::kTH2F, {{nBins, 0, 20}, {100, 0., 1.}}});
    histos.add("TPC/tpcCrossedRowsOverFindableCls", "crossed TPC rows over findable clusters;#it{p}_{T} [GeV/c];crossed rows / findable clusters TPC", {HistType::kTH2F, {{nBins, 0, 20}, {120, 0.0, 1.2}}});
    histos.add("TPC/tpcChi2NCl", "chi2 per cluster in TPC;#it{p}_{T} [GeV/c];chi2 / cluster TPC", {HistType::kTH2F, {{nBins, 0, 20}, {100, 0, 10}}});

    histos.print();
  }

  void process(soa::Join<aod::Collisions, aod::EvSels> const& collisions, 
               soa::Filtered<soa::Join<aod::FullTracks, aod::TracksDCA, aod::TrackSelection, aod::TracksCov>> const& tracks)
  {
    // Loop over collisions
    //int collisionNumber = 0;
    //int maxCollisions = 10;

    for (const auto& collision : collisions) {
      /*if (collisionNumber >= maxCollisions) {
            std::cout << "Reached the maximum number of collisions to print." << std::endl;
            break;
      }*/

      // fill event property variables
      histos.fill(HIST("EventProp/collisionVtxZnoSel"), collision.posZ());
        
      if (!collision.sel8()){
        return;
      }
      histos.fill(HIST("EventProp/collisionVtxZSel8"), collision.posZ());

      if (fabs(collision.posZ()) > ValVtx) {
        return;
      }
      histos.fill(HIST("EventProp/collisionVtxZ"), collision.posZ());
      
      // Loop over tracks
      //int trackNumber = 0;
      //int maxTracks = 10;

      for (const auto& track : tracks) {
        /*if (trackNumber >= maxTracks) {
              std::cout << "Reached the maximum number of tracks to print." << std::endl;
              break;
        }
        
        // Print the track number and its properties
        std::cout << "Track " << trackNumber << ": "
                  << "pt = " << track.pt() << " GeV, "
                  << "phi = " << track.phi() << " rad, "
                  << "eta = " << track.eta() << std::endl;
        trackNumber++;*/
        if (track.collisionId() == collision.globalIndex()){
      
      // Track selection
      //if (selectTrack(track, trackSelection) == true)) {
        /*if (!selectTrack(track, trackSelection)){
          return;
        }*/
     
      // fill kinematic variables
      histos.fill(HIST("Kine/pt"), track.pt());
      histos.fill(HIST("Kine/eta"), track.pt(), track.eta());
      histos.fill(HIST("Kine/phi"), track.pt(), track.phi());
      histos.fill(HIST("Kine/etaVSphi"), track.eta(), track.phi());

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
      histos.fill(HIST("TrackPar/Sigma1Pt"), track.pt(), track.sigma1Pt()*track.pt()*track.pt());
      
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
      //collisionNumber++;
      }
    }
  }
};

//****************************************************************************************
/**
 * Workflow definition.
 */
//****************************************************************************************
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow;
  workflow.push_back(adaptAnalysisTask<TrackJetQa>(cfgc));
  return workflow;
}
