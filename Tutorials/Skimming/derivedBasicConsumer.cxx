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
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>, CERN

// O2 includes
#include "ReconstructionDataFormats/Track.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "DataModel/DerivedExampleTable.h"
#include "PWGJE/DataModel/Jet.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#include "Framework/runDataProcessing.h"

struct DerivedBasicConsumer {
  SliceCache cache;

  /// Function to aid in calculating delta-phi
  /// \param phi1 first phi value
  /// \param phi2 second phi value
  Double_t ComputeDeltaPhi(Double_t phi1, Double_t phi2)
  {
    Double_t deltaPhi = phi1 - phi2;
    if (deltaPhi < -TMath::Pi() / 2.) {
      deltaPhi += 2. * TMath::Pi();
    }
    if (deltaPhi > 3 * TMath::Pi() / 2.) {
      deltaPhi -= 2. * TMath::Pi();
    }
    return deltaPhi;
  }

  Configurable<float> minPtAssoc{"minPtAssoc", 2.0, "min pT associated"};
  Configurable<float> maxPtAssoc{"maxPtAssoc", 4.0, "max pT associated"};
  Configurable<float> minPtTrig{"minPtTrig", 4.0, "min pT trigger"};

  ConfigurableAxis axisPhi{"axisPhi", {72, 0, 2 * TMath::Pi()}, "#phi"};
  ConfigurableAxis axisEta{"axisEta", {80, -0.8, +0.8}, "#eta"};
  ConfigurableAxis axisDeltaPhi{"axisDeltaPhi", {36, -o2::constants::math::PIHalf, o2::constants::math::PIHalf * 3}, "delta #varphi axis for histograms"};

  // Histogram registry: an object to hold your histograms
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Partition<aod::DerivedTracks> triggerTracks = aod::track::pt > minPtTrig;
  Partition<aod::DerivedTracks> assocTracks = aod::track::pt > minPtAssoc&& aod::track::pt < maxPtAssoc;

  Preslice<aod::Tracks> perCollision = aod::track::collisionId;

  void init(InitContext const&)
  {
    // define axes you want to use
    const AxisSpec axisCounter{1, 0, +1, ""};
    histos.add("eventCounter", "eventCounter", kTH1F, {axisCounter});

    // for correlation study
    histos.add("etaHistogramTrigger", "etaHistogramTrigger", kTH1F, {axisEta});
    histos.add("etaHistogramAssoc", "etaHistogramAssoc", kTH1F, {axisEta});
    histos.add("phiHistogramTrigger", "phiHistogramTrigger", kTH1F, {axisPhi});
    histos.add("phiHistogramAssoc", "phiHistogramAssoc", kTH1F, {axisPhi});
    histos.add("correlationFunction", "correlationFunction", kTH1F, {axisDeltaPhi});
  }

  void process(aod::DerivedCollision const& collision, aod::DerivedTracks const& tracks)
  {
    histos.fill(HIST("eventCounter"), 0.5);

    // partitions are not grouped by default
    auto triggerTracksGrouped = triggerTracks->sliceByCached(aod::exampleTrackSpace::derivedCollisionId, collision.globalIndex(), cache);
    auto assocTracksGrouped = assocTracks->sliceByCached(aod::exampleTrackSpace::derivedCollisionId, collision.globalIndex(), cache);

    // Inspect the trigger and associated populations
    for (auto& track : triggerTracksGrouped) {                         //<- only for a subset
      histos.get<TH1>(HIST("etaHistogramTrigger"))->Fill(track.eta()); //<- this should show the selection
      histos.get<TH1>(HIST("ptHistogramTrigger"))->Fill(track.pt());
    }
    for (auto& track : assocTracksGrouped) {                         //<- only for a subset
      histos.get<TH1>(HIST("etaHistogramAssoc"))->Fill(track.eta()); //<- this should show the selection
      histos.get<TH1>(HIST("ptHistogramAssoc"))->Fill(track.pt());
    }

    // Now we do two-particle correlations, using "combinations"
    for (auto& [trackTrigger, trackAssoc] : combinations(o2::soa::CombinationsFullIndexPolicy(triggerTracksGrouped, assocTracksGrouped))) {
      histos.get<TH1>(HIST("correlationFunction"))->Fill(ComputeDeltaPhi(trackTrigger.phi(), trackAssoc.phi()));
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<DerivedBasicConsumer>(cfgc, TaskName{"derived-basic-consumer"})};
  return workflow;
}
