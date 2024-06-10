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

///
/// \file   qaTrackSplitting.cxx
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \brief  Task to analyse the numbers of particles reconstructed more than once
///

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/TrackSelectionTables.h"

using namespace o2::framework;

struct qaTrackSplitting {
  Configurable<int> pdg{"pdg", 2212, "PDG code of the particle to be analysed"};
  struct : ConfigurableGroup {
    Configurable<bool> enableTrackCuts{"enableTrackCuts", false, "Enable the custom track cuts"};
    Configurable<int> itsPattern{"itsPattern", 0, "0 = Run3ITSibAny, 1 = Run3ITSallAny, 2 = Run3ITSall7Layers, 3 = Run3ITSibTwo"};
    Configurable<bool> requireITS{"requireITS", true, "Additional cut on the ITS requirement"};
    Configurable<bool> requireTPC{"requireTPC", true, "Additional cut on the TPC requirement"};
    Configurable<bool> requireGoldenChi2{"requireGoldenChi2", true, "Additional cut on the GoldenChi2"};
    Configurable<int> minITScl{"minITScl", 4, "Additional cut on the ITS cluster"};
    Configurable<float> minNCrossedRowsTPC{"minNCrossedRowsTPC", 70.f, "Additional cut on the minimum number of crossed rows in the TPC"};
    Configurable<float> minNCrossedRowsOverFindableClustersTPC{"minNCrossedRowsOverFindableClustersTPC", 0.8f, "Additional cut on the minimum value of the ratio between crossed rows and findable clusters in the TPC"};
    Configurable<float> maxChi2PerClusterTPC{"maxChi2PerClusterTPC", 4.f, "Additional cut on the maximum value of the chi2 per cluster in the TPC"};
    Configurable<float> maxChi2PerClusterITS{"maxChi2PerClusterITS", 36.f, "Additional cut on the maximum value of the chi2 per cluster in the ITS"};
    Configurable<float> maxDcaXY{"maxDcaXY", 10000.f, "Additional cut on the maximum abs value of the DCA xy"};
    Configurable<float> maxDcaZ{"maxDcaZ", 2.f, "Additional cut on the maximum abs value of the DCA z"};
    Configurable<float> minTPCNClsFound{"minTPCNClsFound", 0.f, "Additional cut on the minimum value of the number of found clusters in the TPC"};
  } cfgCustomTrackCuts;

  // Histograms
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  TrackSelection customTrackCuts;

  void init(InitContext&)
  {
    histos.add("tracks", "tracsk", kTH1D, {{10, -0.5, 9.5, "Track selection"}});
    histos.add("numberOfRecoed", "recoed", kTH1D, {{10, -0.5, 9.5, "Number of tracks associated to a particle"}});
    histos.add("map", "map", kTH3D, {{100, -1, 1, "#Delta #eta"}, {100, -1, 1, "#Delta #varphi"}, {100, -1, 1, "#Delta #it{p}_{T}"}});
    histos.add("mapMC", "mapMC", kTH3D, {{100, -1, 1, "#Delta #eta"}, {100, -1, 1, "#Delta #varphi"}, {100, -1, 1, "#Delta #it{p}_{T}"}});

    customTrackCuts = getGlobalTrackSelectionRun3ITSMatch(cfgCustomTrackCuts.itsPattern);
    LOG(info) << "Customizing track cuts:";
    customTrackCuts.SetRequireITSRefit(cfgCustomTrackCuts.requireITS);
    customTrackCuts.SetRequireTPCRefit(cfgCustomTrackCuts.requireTPC);
    customTrackCuts.SetRequireGoldenChi2(cfgCustomTrackCuts.requireGoldenChi2);
    customTrackCuts.SetRequireHitsInITSLayers(cfgCustomTrackCuts.minITScl.value, {0, 1, 2, 3, 4, 5, 6});
    customTrackCuts.SetMaxChi2PerClusterTPC(cfgCustomTrackCuts.maxChi2PerClusterTPC);
    customTrackCuts.SetMaxChi2PerClusterITS(cfgCustomTrackCuts.maxChi2PerClusterITS);
    customTrackCuts.SetMinNCrossedRowsTPC(cfgCustomTrackCuts.minNCrossedRowsTPC);
    customTrackCuts.SetMinNClustersTPC(cfgCustomTrackCuts.minTPCNClsFound);
    customTrackCuts.SetMinNCrossedRowsOverFindableClustersTPC(cfgCustomTrackCuts.minNCrossedRowsOverFindableClustersTPC);
    customTrackCuts.SetMaxDcaXYPtDep([&](float /*pt*/) { return cfgCustomTrackCuts.maxDcaXY; }); // No DCAxy cut will be used, this is done via the member function of the task
    customTrackCuts.SetMaxDcaZ(cfgCustomTrackCuts.maxDcaZ);
    customTrackCuts.print();
  }

  // Global process
  using TrackCandidates = o2::soa::Join<o2::aod::Tracks, o2::aod::TracksExtra, o2::aod::TracksDCA, o2::aod::McTrackLabels>;
  void process(o2::soa::Join<o2::aod::Collisions, o2::aod::McCollisionLabels, o2::aod::EvSels>::iterator const& collision,
               TrackCandidates const& tracks,
               o2::aod::McParticles const&)
  {
    if (!collision.sel8()) {
      return;
    }
    typedef std::shared_ptr<TrackCandidates::iterator> trkType;

    std::map<int64_t, std::vector<trkType>> particleUsageCounter;
    for (auto track : tracks) {
      histos.fill(HIST("tracks"), 0);
      if (!track.has_mcParticle()) {
        continue;
      }
      histos.fill(HIST("tracks"), 1);
      if (track.mcParticle().pdgCode() != pdg) {
        continue;
      }
      histos.fill(HIST("tracks"), 2);
      if (!track.mcParticle().isPhysicalPrimary()) {
        continue;
      }
      histos.fill(HIST("tracks"), 3);
      if (cfgCustomTrackCuts.enableTrackCuts.value && !customTrackCuts.IsSelected(track)) {
        continue;
      }
      histos.fill(HIST("tracks"), 4);
      particleUsageCounter[track.mcParticleId()].push_back(std::make_shared<decltype(track)>(track));
    }
    for (const auto& [mcId, tracks] : particleUsageCounter) {
      histos.fill(HIST("numberOfRecoed"), tracks.size());
      if (tracks.size() > 1) {
        bool isFirst = true;
        for (const auto& track : tracks) {
          if (isFirst) {
            isFirst = false;
            histos.fill(HIST("mapMC"), track->eta() - track->mcParticle().eta(), track->phi() - track->mcParticle().phi(), track->pt() - track->mcParticle().pt());
            continue;
          }
          histos.fill(HIST("map"), track->eta() - tracks[0]->eta(), track->phi() - tracks[0]->phi(), track->pt() - tracks[0]->pt());
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<qaTrackSplitting>(cfgc)}; }
