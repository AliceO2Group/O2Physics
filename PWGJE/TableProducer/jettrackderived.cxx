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
/// \author Johanna Lömker <johanna.lomker@cern.ch>
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
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"

using namespace o2;
using namespace o2::track;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct jetspectraDerivedMaker {
  Configurable<float> fractionOfEvents{"fractionOfEvents", 1.f, "Downsampling factor for the events for derived data"};

  Configurable<double> ValVtx{"ValVtx", 10, "Value of the vertex position"};
  Configurable<float> ValCutEta{"ValCutEta", 0.8f, "Eta range for tracks"};
  Configurable<float> minPt{"minPt", 0.15f, "minimum pT for tracks"};
  Configurable<float> maxPt{"maxPt", 10e10, "maximum pT for tracks"};
  Configurable<bool> fillMultiplicity{"fillMultiplicity", true, "To fill multiplicity and centrality histograms"};

  // Custom track cuts for the cut variation study
  TrackSelection customTrackCuts;
  Configurable<int> itsPattern{"itsPattern", 2, "0 = Run3ITSibAny, 1 = Run3ITSibTwo, 2 = Run3ITSallAny, 3 = Run3ITSall7Layers"};
  Configurable<bool> requireITS{"requireITS", true, "Additional cut on the ITS requirement"};
  Configurable<bool> requireTPC{"requireTPC", true, "Additional cut on the TPC requirement"};
  Configurable<bool> requireGoldenChi2{"requireGoldenChi2", true, "Additional cut on the GoldenChi2"};
  Configurable<float> minNCrossedRowsTPC{"minNCrossedRowsTPC", 60.f, "Additional cut on the minimum number of crossed rows in the TPC"};
  Configurable<float> minNCrossedRowsOverFindableClustersTPC{"minNCrossedRowsOverFindableClustersTPC", 0.7f, "Additional cut on the minimum value of the ratio between crossed rows and findable clusters in the TPC"};
  Configurable<float> maxChi2PerClusterTPC{"maxChi2PerClusterTPC", 7.f, "Additional cut on the maximum value of the chi2 per cluster in the TPC"};
  Configurable<float> maxChi2PerClusterITS{"maxChi2PerClusterITS", 36.f, "Additional cut on the maximum value of the chi2 per cluster in the ITS"};
  Configurable<float> maxDcaXY{"maxDcaXY", 0.25f, "Cut on the maximum value of the DCA xy "};
  Configurable<float> maxDcaZ{"maxDcaZ", 3.f, "Additional cut on the maximum value of the DCA z"};
  Configurable<float> minTPCNClsFound{"minTPCNClsFound", 0.f, "Additional cut on the minimum value of the number of found clusters in the TPC"};

  // Histograms
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Configurable<int> nBins{"nBins", 200, "N bins in histos"};
  ConfigurableAxis binsMultiplicity{"binsMultiplicity", {100, 0, 100}, "Binning for multiplicity"};
  ConfigurableAxis binsPercentile{"binsPercentile", {100, 0, 100}, "Binning for percentiles"};

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
    LOG(info) << "\tRequireHitsInITSLayers=" << maxChi2PerClusterITS.value;
    LOG(info) << "\tmaxDcaXY=" << maxDcaXY.value;
    LOG(info) << "\tmaxDcaZ=" << maxDcaZ.value;
    LOG(info) << "\tminPt=" << minPt.value;
    LOG(info) << "\tmaxPt=" << maxPt.value;
    LOG(info) << "\tmaxEta=" << ValCutEta.value;

    customTrackCuts = getGlobalTrackSelectionRun3ITSMatch(itsPattern.value);
    LOG(info) << "Customizing track cuts:";
    customTrackCuts.SetEtaRange(-ValCutEta.value, ValCutEta.value);
    customTrackCuts.SetPtRange(minPt.value, maxPt.value);
    customTrackCuts.SetRequireITSRefit(requireITS.value);
    customTrackCuts.SetRequireTPCRefit(requireTPC.value);
    customTrackCuts.SetRequireGoldenChi2(requireGoldenChi2.value);
    customTrackCuts.SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC.value);
    customTrackCuts.SetMaxChi2PerClusterITS(maxChi2PerClusterITS.value);
    customTrackCuts.SetMinNCrossedRowsTPC(minNCrossedRowsTPC.value);
    customTrackCuts.SetMinNClustersTPC(minTPCNClsFound.value);
    // customTrackCuts.SetRequireHitsInITSLayers(nHits.value, {0, 1}); // one hit in any SPD layer (#hits, {layer0, layer1,...})
    customTrackCuts.SetMinNCrossedRowsOverFindableClustersTPC(minNCrossedRowsOverFindableClustersTPC.value);
    customTrackCuts.SetMaxDcaXY(maxDcaXY.value);
    customTrackCuts.SetMaxDcaZ(maxDcaZ.value);
    customTrackCuts.print();

    // event property histograms
    histos.add("EventProp/collisionVtxZ", "Collsion Vertex Z;#it{Vtx}_{z} [cm];number of entries", HistType::kTH1F, {{nBins, -20, 20}});
    histos.add("EventProp/sampledvertexz", "Sampled collsion Vertex Z with event (sel8) selection;#it{Vtx}_{z} [cm];number of entries", HistType::kTH1F, {{nBins, -20, 20}});
    histos.add("EventProp/NumContrib", "Number of contributors to vertex of collision; number of contributors to vtx; number of entries", HistType::kTH1F, {{nBins, 0, 600}});
    histos.add("EventProp/rejectedCollId", "CollisionId of collisions that did not pass the event selection; collisionId; number of entries", HistType::kTH1F, {{10, 0, 5}});

    const AxisSpec axisPercentileFT0A{binsPercentile, "Centrality FT0A"};
    const AxisSpec axisPercentileFT0C{binsPercentile, "Centrality FT0C"};
    const AxisSpec axisMultiplicityPV{binsMultiplicity, "MultNTracksPV"};
    const AxisSpec axisMultiplicityFT0A{binsMultiplicity, "Multiplicity FT0A"};
    const AxisSpec axisMultiplicityFT0C{binsMultiplicity, "Multiplicity FT0C"};

    histos.add("Centrality/FT0A", "CentFT0A", HistType::kTH1D, {axisPercentileFT0A});
    histos.add("Centrality/FT0C", "CentFT0C", HistType::kTH1D, {axisPercentileFT0C});
    histos.add("Mult/NTracksPV", "MultNTracksPV", HistType::kTH1D, {axisMultiplicityPV});
    histos.add("Mult/FT0A", "MultFT0A", HistType::kTH1D, {axisMultiplicityFT0A});
    histos.add("Mult/FT0C", "MultFT0C", HistType::kTH1D, {axisMultiplicityFT0C});
  }

  template <typename CollisionType>
  bool isEventSelected(CollisionType const& collision)
  {
    // here we already fill some event histos for cross checks
    if (!collision.sel8()) {
      histos.fill(HIST("EventProp/rejectedCollId"), 2);
      return false;
    }

    if (abs(collision.posZ()) > ValVtx) {
      histos.fill(HIST("EventProp/rejectedCollId"), 3);
      return false;
    }
    histos.fill(HIST("EventProp/collisionVtxZ"), collision.posZ());                                                              // test fill
                                                                                                                                 //  Last thing, check the sampling
    if (fractionOfEvents < 1.f && (static_cast<float>(rand_r(&randomSeed)) / static_cast<float>(RAND_MAX)) > fractionOfEvents) { // Skip events that are not sampled
      histos.fill(HIST("EventProp/rejectedCollId"), 4);
      return false;
    }
    histos.fill(HIST("EventProp/sampledvertexz"), collision.posZ());
    histos.fill(HIST("EventProp/NumContrib"), collision.numContrib());
    if (fillMultiplicity == true) {
      histos.fill(HIST("Centrality/FT0A"), collision.centFT0A());
      histos.fill(HIST("Centrality/FT0C"), collision.centFT0C());
      histos.fill(HIST("Mult/FT0C"), collision.multFT0C());
      histos.fill(HIST("Mult/FT0A"), collision.multFT0A());
      histos.fill(HIST("Mult/NTracksPV"), collision.multNTracksPV());
    }
    return true;
  }

  Preslice<aod::Track> trackPerColl = aod::track::collisionId;
  Produces<o2::aod::JeTracks> tableTrack;
  Produces<o2::aod::JeColls> tableColl;
  using CollisionCandidate = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>;
  using TrackCandidates = soa::Join<aod::FullTracks, aod::TracksDCA, aod::TrackSelection, aod::TracksCov>;
  unsigned int randomSeed = 0;
  void processData(CollisionCandidate const& collisions,
                   TrackCandidates const& tracks, aod::BCs const&)
  {
    for (const auto& collision : collisions) {
      if (!isEventSelected(collision)) {
        histos.fill(HIST("EventProp/rejectedCollId"), 1);
        continue;
      } else {
        auto tracksInCollision = tracks.sliceBy(trackPerColl, collision.globalIndex());
        tableColl(collision.globalIndex(),
                  collision.collisionTime(),
                  collision.numContrib(),
                  collision.posX(),
                  collision.posY(),
                  collision.posZ(),
                  collision.sel8(),
                  tracksInCollision.size(),
                  collision.multNTracksPV(),
                  collision.multFT0A(),
                  collision.multFT0C(),
                  collision.centFT0A(),
                  collision.centFT0C(),
                  collision.bc().runNumber());
        tableTrack.reserve(tracks.size());
        for (const auto& trk : tracksInCollision) {
          if (!customTrackCuts.IsSelected(trk)) {
            continue;
          } else {
            tableTrack(trk.collisionId(),
                       trk.trackTime(),
                       trk.signed1Pt(), trk.eta(), trk.phi(), trk.pt(),
                       trk.sigma1Pt(),
                       trk.alpha(),
                       trk.x(), trk.y(), trk.z(),
                       trk.snp(),
                       trk.tgl(),
                       trk.isPVContributor(),
                       trk.hasTRD(),
                       trk.hasITS(),
                       trk.hasTPC(),
                       trk.isGlobalTrack(),
                       trk.isGlobalTrackWoDCA(),
                       trk.isGlobalTrackWoPtEta(),
                       trk.flags(),
                       trk.trackType(),
                       trk.length(),
                       trk.tpcChi2NCl(), trk.itsChi2NCl(), trk.tofChi2(),
                       trk.tpcNClsShared(),
                       trk.tpcNClsFindable(),
                       trk.tpcNClsFindableMinusFound(),
                       trk.tpcNClsFindableMinusCrossedRows(),
                       trk.itsClusterMap(),
                       trk.itsNCls(),
                       trk.tpcFractionSharedCls(),
                       trk.tpcNClsFound(),
                       trk.tpcNClsCrossedRows(),
                       trk.tpcCrossedRowsOverFindableCls(),
                       trk.tpcFoundOverFindableCls(),
                       trk.dcaXY(),
                       trk.dcaZ());
          }
        }
      }
    }
  }

  PROCESS_SWITCH(jetspectraDerivedMaker, processData, "Process collision data for derived dataset production", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<jetspectraDerivedMaker>(cfgc)};
}
