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

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#include "Framework/runDataProcessing.h"

struct DerivedBasicConsumer {
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

  Configurable<float> associatedMinPt{"associatedMinPt", 4.0f, "NSassociatedMinPt"};
  Configurable<float> associatedMaxPt{"associatedMaxPt", 6.0f, "associatedMaxPt"};
  Configurable<float> triggerMinPt{"triggerMinPt", 6.0f, "triggerMinPt"};
  ConfigurableAxis axisPt{"axisPt", {200,0.0f,20.0f}, "pt axis"};

  SliceCache cache;

  // define partitions
  Partition<aod::DrTracks> associatedTracks = aod::exampleTrackSpace::pt < associatedMaxPt && aod::exampleTrackSpace::pt > associatedMinPt;
  Partition<aod::DrTracks> triggerTracks = aod::exampleTrackSpace::pt > triggerMinPt;

  Filter collZfilter = nabs(aod::collision::posZ) < 10.0f;

  // Histogram registry: an object to hold your histograms
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {
    // define axes you want to use
    const AxisSpec axisCounter{1, 0, +1, ""};
    const AxisSpec axisPVz{300, -15.0f, +15.0f, ""};
    const AxisSpec axisDeltaPhi{100, -0.5*o2::constants::math::PI, +1.5*o2::constants::math::PI, "#Delta#phi"};
    const AxisSpec axisDeltaEta{100, -1.0, +1.0, "#Delta#eta"};

    histos.add("eventCounter", "eventCounter", kTH1D, {axisCounter});
    histos.add("hEventPVz", "hEventPVz", kTH1D, {axisPVz});

    histos.add("ptAssoHistogram", "ptAssoHistogram", kTH1D, {axisPt});
    histos.add("ptTrigHistogram", "ptTrigHistogram", kTH1D, {axisPt});

    histos.add("correlationFunction", "correlationFunction", kTH1D, {axisDeltaPhi});
    histos.add("correlationFunctionO2", "correlationFunctionO2", kTH1D, {axisDeltaPhi});

    histos.add("correlationFunction2d", "correlationFunction2d", kTH2F, {axisDeltaPhi, axisDeltaEta});
  }

  void process(soa::Filtered<aod::DrCollisions>::iterator const& collision, aod::DrTracks const&)
  {
    histos.fill(HIST("eventCounter"), 0.5);
    histos.fill(HIST("hEventPVz"), collision.posZ());

    auto assoTracksThisCollision = associatedTracks->sliceByCached(aod::exampleTrackSpace::drCollisionId, collision.globalIndex(), cache);
    auto trigTracksThisCollision = triggerTracks->sliceByCached(aod::exampleTrackSpace::drCollisionId, collision.globalIndex(), cache);

    for (auto& track : assoTracksThisCollision)
      histos.fill(HIST("ptAssoHistogram"), track.pt());
    for (auto& track : trigTracksThisCollision)
      histos.fill(HIST("ptTrigHistogram"), track.pt());

    for (auto& trigger : trigTracksThisCollision){
      for (auto& associated : assoTracksThisCollision){
        histos.fill(HIST("correlationFunction"), ComputeDeltaPhi(trigger.phi(),associated.phi()));
      }
    }

    for (auto& [trigger, associated] : combinations(o2::soa::CombinationsFullIndexPolicy(trigTracksThisCollision, assoTracksThisCollision))) {
      histos.fill(HIST("correlationFunctionO2"), ComputeDeltaPhi(trigger.phi(),associated.phi()));
      histos.fill(HIST("correlationFunction2d"), ComputeDeltaPhi(trigger.phi(),associated.phi()), trigger.eta() - associated.eta());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<DerivedBasicConsumer>(cfgc, TaskName{"derived-basic-consumer"})};
  return workflow;
}