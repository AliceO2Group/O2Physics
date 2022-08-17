// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Centrality.h"

#include "Common/DataModel/PIDResponse.h"

#include <TH1F.h>
#include <TParameter.h>

#include <cmath>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct SimpleHistogramTask {
  HistogramRegistry registry{
    "registry",
    {
      // name, title, hist type, vector of axes
      // here each axis is defined as: {nBins, min, max, optional title, optional name}
      // alternatively, for variable binning: {{binEdge1, binEdge2, ...}, optional title, optional name}
      {"eta", "#eta", {HistType::kTH1F, {{102, -2.01, 2.01}}}},                                               //
      {"phi", "#varphi", {HistType::kTH1F, {{100, 0., 2. * M_PI}}}},                                          //
      {"pt", "pt", {HistType::kTH1F, {{100, -0.01, 10.01, "#it{p}_{T} (GeV/c)"}}}},                           //
      {"ptVariableBinning", "pt", {HistType::kTH1F, {{                                                        //
                                                      {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, //
                                                       1.1, 1.2, 1.3, 1.4, 1.5, 2.0, 5.0, 10.0, 20.0, 50.0},
                                                      "#it{p}_{T} (GeV/c)"}}}} //
    }                                                                          //
  };

  // aod::Track is an iterator over tracks table
  void process(aod::Track const& track)
  {

    registry.get<TH1>(HIST("eta"))->Fill(track.eta());
    registry.get<TH1>(HIST("phi"))->Fill(track.phi());
    registry.get<TH1>(HIST("pt"))->Fill(track.pt());
    registry.get<TH1>(HIST("ptVariableBinning"))->Fill(track.pt());
  }
};

//start TPCvsPT ZCH
struct HistogramsWithTpc {
  HistogramRegistry registry{
    "registry",
    {
      {"TPCSignal", "TPC signal", {HistType::kTH2F, {{500, 0, 10, "pt"}, {1000, 0, 600, "tpc"}}}} //
    }                                                                                             //
  };

  // If we fill histograms with filters, we need to provide full tables
  // aod::Tracks instead of aod::Track

  Filter trackFilter = requireGlobalTrackInFilter();
  using FilteredTrack = soa::Filtered<soa::Join<aod::FullTracks, aod::TrackSelection>>::iterator;
  //using pidTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::McTrackLabels, aod::TracksDCA, aod::TrackSelection, aod::pidTOFbeta, aod::TOFSignal>>::iterator;
  void process(FilteredTrack const& track)
  {
    registry.fill(HIST("TPCSignal"), track.pt(), track.tpcSignal());
  }
};
//end TPCvsPT ZCH

//start TOFvsPT ZCH

//using BigTrack = soa::Filtered<soa::Join<aod::FullTracks, aod::TracksDCA, aod::pidTOFbeta, aod::pidTPCFullEl, aod::pidTOFFullEl, aod::pidTPCFullMu, aod::pidTOFFullMu, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::TrackSelection, aod::TOFSignal>>::iterator;

struct HistogramsWithTof {
  HistogramRegistry registry{
    "registry",
    {
      {"TOFSignal", "TOF signal", {HistType::kTH2F, {{500, 0, 10, "pt"}, {5000, 0, 50000, "tof"}}}}, //
      {"TOFbeta", "TOF beta", {HistType::kTH2F, {{500, 0, 10, "pt"}, {500, 0, 2, "beta"}}}}          //
    }                                                                                                //
  };

  // If we fill histograms with filters, we need to provide full tables
  // aod::Tracks instead of aod::Track
  //void process(BigTrack const& track)

  Filter trackFilter = requireGlobalTrackInFilter();
  using FilteredTrackWithTOF = soa::Filtered<soa::Join<aod::FullTracks, aod::pidTOFbeta, aod::TOFSignal, aod::TrackSelection>>::iterator;

  void process(FilteredTrackWithTOF const& track)
  {
    registry.fill(HIST("TOFSignal"), track.pt(), track.tofSignal());
    registry.fill(HIST("TOFbeta"), track.pt(), track.beta());
  }
};
//end TOFvsPT ZCH

//start TRDvsPT ZCH
struct HistogramsWithTrd {
  HistogramRegistry registry{
    "registry",
    {
      {"TRDSignal", "TRD signal", {HistType::kTH2F, {{500, 0, 10, "pt"}, {2500, 0, 100, "trd"}}}} //
    }                                                                                             //
  };

  // If we fill histograms with filters, we need to provide full tables
  // aod::Tracks instead of aod::Track

  Filter trackFilter = requireGlobalTrackInFilter();
  using FilteredTracks = soa::Filtered<soa::Join<aod::FullTracks, aod::TrackSelection>>;

  void process(aod::Collision const& col, FilteredTracks const& tracks)
  {
    for (const auto& track : tracks) {
      registry.fill(HIST("TRDSignal"), track.pt(), track.trdSignal());
    }
  }
};
//end TOFvsPT ZCH

//start ETAvsPHI
struct MultAndCentrality {
  HistogramRegistry registry{
    "registry",
    {{"mult", "multiplicity", {HistType::kTH1F, {{102, 0, 1000}}}}} //
  };

  // If we fill histograms with filters, we need to provide full tables
  // aod::Tracks instead of aod::Track
  //using CollisionInstance = aod::Collisions::iterator;
  using CollisionInstance = soa::Join<aod::Collisions, aod::CentRun2V0Ms>::iterator;
  void process(CollisionInstance const& collision, aod::Tracks const& tracks)
  //void process(aod::FullTrack const& track)
  {
    registry.fill(HIST("mult"), 1);
    registry.get<TH1>(HIST("mult"))->Fill(2);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<SimpleHistogramTask>(cfgc),
    //adaptAnalysisTask<MultAndCentrality>(cfgc),
    adaptAnalysisTask<HistogramsWithTof>(cfgc),
    adaptAnalysisTask<HistogramsWithTrd>(cfgc),
    adaptAnalysisTask<HistogramsWithTpc>(cfgc)

  };
}
