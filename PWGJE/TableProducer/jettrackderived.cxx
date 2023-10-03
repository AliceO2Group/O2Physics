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
/// \file   for derived process JetTrackQa.h and .cxx
/// \author Johanna Lömker
/// \since  2023-10-02
/// \brief  Header for the trackJetQa task for the analysis of the tracks for jets..
///

// O2 includes
#include "ReconstructionDataFormats/Track.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/AnalysisDataModel.h"
#include "PWGJE/DataModel/Jet.h"
#include "Framework/StaticFor.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "PWGJE/DataModel/TrackJetQa.h"

using namespace o2;
using namespace o2::track;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct spectraDerivedMaker {
  Configurable<int> nBins{"nBins", 200, "N bins in histos"};

  Configurable<double> ValVtx{"ValVtx", 10, "Value of the vertex position"};
  Configurable<float> ValCutEta{"ValCutEta", 0.8f, "Eta range for tracks"};
  Configurable<float> ValCutY{"ValCutY", 0.5f, "Y range for tracks"};
  Configurable<float> fractionOfEvents{"fractionOfEvents", 2.f, "Downsampling factor for the events for derived data"};
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

  // Histograms
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext&)
  {
    // Custom track cuts
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

    // event property histograms
    histos.add("EventProp/collisionVtxZ", "Collsion Vertex Z;#it{Vtx}_{z} [cm];number of entries", HistType::kTH1F, {{nBins, -20, 20}});
    histos.add("EventProp/collisionVtxZnoSel", "Collsion Vertex Z without event selection;#it{Vtx}_{z} [cm];number of entries", HistType::kTH1F, {{nBins, -20, 20}});
    histos.add("EventProp/collisionVtxZSel8", "Collsion Vertex Z with event selection;#it{Vtx}_{z} [cm];number of entries", HistType::kTH1F, {{nBins, -20, 20}});
    histos.add("EventProp/sampledvertexz", "Sampled collsion Vertex Z with event selection;#it{Vtx}_{z} [cm];number of entries", HistType::kTH1F, {{nBins, -20, 20}});
  }

  template <typename CollisionType, typename TrackType>
  bool isEventSelected(CollisionType const& collision, TrackType const& tracks)
  {
    // here we could already fill some histos fo cross checks
    if (!collision.sel8()) {
      return false;
    }

    if (abs(collision.posZ()) > ValVtx) {
      return false;
    }
    histos.fill(HIST("EventProp/collisionVtxZ"), collision.posZ());                                                              // test fill
                                                                                                                                 //  Last thing, check the sampling
    if (fractionOfEvents < 1.f && (static_cast<float>(rand_r(&randomSeed)) / static_cast<float>(RAND_MAX)) > fractionOfEvents) { // Skip events that are not sampled
      return false;
    }
    histos.fill(HIST("EventProp/sampledvertexz"), collision.posZ());
    return true;
  }

  template <typename TrackType>
  bool isTrackSelected(TrackType const& track) // add trackselections and corresponding histos for cross checks to derived table
  {
    if (!track.isGlobalTrackWoPtEta()) { // in principle we would liek to check all these cuts
      return false;
    }
    return true;
  }

  using CollisionCandidate = soa::Join<aod::Collisions, aod::EvSels>;
  using TrackCandidates = soa::Join<aod::FullTracks, aod::TracksDCA, aod::TrackSelection, aod::TracksCov>;

  Produces<o2::aod::SpColls> tableColl;
  Produces<o2::aod::SpTracks> tableTrack;
  unsigned int randomSeed = 0;
  void processData(CollisionCandidate::iterator const& collision,
                   TrackCandidates const& tracks,
                   aod::BCs const&)
  {
    if (!isEventSelected(collision, tracks)) {
      return;
    }

    tableColl(collision.numContrib(),
              collision.posX(),
              collision.posY(),
              collision.posZ(),
              collision.sel8(),
              collision.bc().runNumber());

    tableTrack.reserve(tracks.size());
    for (const auto& trk : tracks) {
      if (!isTrackSelected(trk)) {
        return;
      }

      tableTrack(tableColl.lastIndex(),
                 trk.pt() * trk.sign(), trk.eta(), trk.phi(),
                 trk.sigma1Pt(),
                 o2::aod::spectra::packInTable<o2::aod::spectra::binningDCA>(trk.dcaXY()),
                 o2::aod::spectra::packInTable<o2::aod::spectra::binningDCA>(trk.dcaZ()),
                 trk.length(),
                 trk.tpcSignal(),
                 trk.tpcChi2NCl(), trk.itsChi2NCl(), trk.tofChi2(),
                 trk.tpcNClsShared(),
                 trk.tpcNClsFindable(),
                 trk.tpcNClsFindableMinusFound(),
                 trk.tpcNClsFindableMinusCrossedRows(),
                 trk.isPVContributor(),
                 trk.itsClusterMap(),
                 trk.hasTRD(),
                 trk.isGlobalTrack(),
                 trk.isGlobalTrackWoDCA(),
                 trk.isGlobalTrackWoPtEta());
    }
  }
  PROCESS_SWITCH(spectraDerivedMaker, processData, "Process data for derived dataset production", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<spectraDerivedMaker>(cfgc)}; }
