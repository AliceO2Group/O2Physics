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

// \brief   Task for the calculation of SPC with filtered data.
// \author  Maxim Virta (maxim.virta@cern.ch), Cindy Mordasini (cindy.mordasini@cern.ch)

// Standard headers.
#include <chrono>
#include <string>
#include <vector>
#include <TRandom3.h>

// O2 headers. //
// The first two are mandatory.
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"

#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"

// O2 Physics headers. //
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGCF/DataModel/CorrelationsDerived.h"
#include "PWGCF/JCorran/DataModel/JCatalyst.h"
#include "PWGCF/JCorran/Core/FlowJSPCAnalysis.h"
#include "PWGCF/JCorran/Core/FlowJSPCObservables.h"
#include "PWGCF/JCorran/Core/FlowJHistManager.h"

// Namespaces and definitions.
using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults,
                               aod::FT0sCorrected, aod::CentFT0Ms,
                               aod::CentFT0As, aod::CentFT0Cs, aod::CentFV0As,
                               aod::CentFDDMs, aod::CentNTPVs>;

using MyTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::JWeights>;

struct flowJSPCAnalysis {
  HistogramRegistry spcHistograms{"SPCResults", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  FlowJSPCAnalysis spcAnalysis;
  FlowJSPCAnalysis::JQVectorsT jqvecs;
  template <class T>
  using HasWeightNUA = decltype(std::declval<T&>().weightNUA());

  HistogramRegistry qaHistRegistry{"qaHistRegistry", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  FlowJHistManager histManager;

  FlowJSPCObservables spcObservables;

  // Set Configurables here
  Configurable<bool> cfgFillQA{"cfgFillQA", true, "Fill QA plots"};

  Configurable<int> cfgWhichSPC{"cfgWhichSPC", 0, "Which SPC observables to compute."};

  struct : ConfigurableGroup {
    Configurable<float> cfgPtMin{"cfgPtMin", 0.2f, "Minimum pT used for track selection."};
    Configurable<float> cfgPtMax{"cfgPtMax", 5.0f, "Maximum pT used for track selection."};
    Configurable<float> cfgEtaMax{"cfgEtaMax", 0.8f, "Maximum eta used for track selection."};
  } cfgTrackCuts;

  struct : ConfigurableGroup {
    Configurable<int> cfgCentEst{"cfgCentEst", 2, "Centrality estimator."};
    Configurable<float> cfgZvtxMax{"cfgZvtxMax", 10.0f, "Maximum primary vertex cut applied for the events."};
    Configurable<int> cfgMultMin{"cfgMultMin", 10, "Minimum number of particles required for the event to have."};
  } cfgEventCuts;

  // // Filters to be applied to the received data.
  // // The analysis assumes the data has been subjected to a QA of its selection,
  // // and thus only the final distributions of the data for analysis are saved.
  Filter collFilter = (nabs(aod::collision::posZ) < cfgEventCuts.cfgZvtxMax);
  Filter trackFilter = (aod::track::pt > cfgTrackCuts.cfgPtMin) && (aod::track::pt < cfgTrackCuts.cfgPtMax) && (nabs(aod::track::eta) < cfgTrackCuts.cfgEtaMax);
  Filter cftrackFilter = (aod::cftrack::pt > cfgTrackCuts.cfgPtMin) && (aod::cftrack::pt < cfgTrackCuts.cfgPtMax); // eta cuts done by jfluc

  void init(InitContext const&)
  {
    // Add histomanager here
    spcAnalysis.setHistRegistry(&spcHistograms);
    spcAnalysis.createHistos();

    spcObservables.setSPCObservables(cfgWhichSPC);
    spcAnalysis.setFullCorrSet(spcObservables.harmonicArray);

    histManager.setHistRegistryQA(&qaHistRegistry);
    histManager.setDebugLog(false);
    histManager.createHistQA();
  }

  template <class CollisionT, class TrackT>
  void analyze(CollisionT const& collision, TrackT const& tracks)
  // void process(soa::Filtered<MyCollisions>::iterator const& coll, soa::Filtered<soa::Join<aod::MyTracks, aod::JWeights>> const& tracks)
  {
    if (tracks.size() < cfgEventCuts.cfgMultMin)
      return;

    float cent = collision.multiplicity();
    if (cent < 0. || cent > 100.) {
      return;
    }
    int cBin = histManager.getCentBin(cent);
    spcHistograms.fill(HIST("FullCentrality"), cent);
    int nTracks = tracks.size();
    for (const auto& track : tracks) {
      if (cfgFillQA) {
        // histManager.FillTrackQA<0>(track, cBin, collision.posZ());

        using JInputClassIter = typename TrackT::iterator;
        if constexpr (std::experimental::is_detected<HasWeightNUA, const JInputClassIter>::value) {
          spcAnalysis.fillQAHistograms(cBin, track.phi(), 1. / track.weightNUA());
        }
      }
    }

    if (cfgFillQA)
      histManager.fillEventQA<1>(collision, cBin, cent, nTracks);

    jqvecs.Calculate(tracks, 0.0, cfgTrackCuts.cfgEtaMax);
    spcAnalysis.setQvectors(&jqvecs);
    spcAnalysis.calculateCorrelators(cBin);
  }

  void processJDerived(aod::JCollision const& collision, soa::Filtered<aod::JTracks> const& tracks)
  {
    analyze(collision, tracks);
  }
  PROCESS_SWITCH(flowJSPCAnalysis, processJDerived, "Process derived data", false);

  void processJDerivedCorrected(aod::JCollision const& collision, soa::Filtered<soa::Join<aod::JTracks, aod::JWeights>> const& tracks)
  {
    analyze(collision, tracks);
  }
  PROCESS_SWITCH(flowJSPCAnalysis, processJDerivedCorrected, "Process derived data with corrections", false);

  void processCFDerived(aod::CFCollision const& collision, soa::Filtered<aod::CFTracks> const& tracks)
  {
    analyze(collision, tracks);
  }
  PROCESS_SWITCH(flowJSPCAnalysis, processCFDerived, "Process CF derived data", false);

  void processCFDerivedCorrected(aod::CFCollision const& collision, soa::Filtered<soa::Join<aod::CFTracks, aod::JWeights>> const& tracks)
  {
    analyze(collision, tracks);
  }
  PROCESS_SWITCH(flowJSPCAnalysis, processCFDerivedCorrected, "Process CF derived data with corrections", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<flowJSPCAnalysis>(cfgc)};
}
