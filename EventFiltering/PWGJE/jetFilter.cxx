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

// Author: Filip Krizek

#include <TMath.h>
#include <cmath>
#include <string>

#include "Framework/ASoA.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/DataModel/EMCALClusters.h"
#include "PWGJE/DataModel/Jet.h"

#include "../filterTables.h"

#include "Framework/HistogramRegistry.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

static const std::vector<std::string> highPtObjectsNames{"JetChHighPt"};

struct jetFilter {
  enum { kJetChLowPt = 0,
         kJetChHighPt = 1,
         kHighPtObjects };

  Produces<aod::JetFilters> tags;

  Configurable<float> cfgJetR{"cfgJetR", 0.6,
                              "jet resolution parameter"}; // jet cone radius

  Configurable<float> jetPtLowThreshold{"jetPtLowThreshold", 30.0, "threshold for low pT jet trigger"};
  Configurable<float> jetPtHighThreshold{"jetPtHighThreshold", 50.0, "threshold for high pT jet trigger"};

  HistogramRegistry spectra{
    "spectra",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    true,
    true};

  void init(o2::framework::InitContext&)
  {

    spectra.add("fCollZpos", "collision z position", HistType::kTH1F,
                {{200, -20., +20., "#it{z}_{vtx} position (cm)"}});

    spectra.add("ptphiJetChSelected_lowptjettrigger", "pT of selected low pT charged jet trigger vs phi", HistType::kTH2F,
                {{150, 0., +150., "charged jet #it{p}_{T} (GeV/#it{c})"},
                 {60, 0, TMath::TwoPi()}});

    spectra.add("ptphiJetChSelected_highptjettrigger", "pT of selected high pT charged jet trigger vs phi", HistType::kTH2F,
                {{150, 0., +150., "charged jet #it{p}_{T} (GeV/#it{c})"},
                 {60, 0, TMath::TwoPi()}});

    spectra.add("ptetaJetChSelected_lowptjettrigger", "pT of selected low pT charged jet trigger vs eta", HistType::kTH2F,
                {{150, 0., +150., "charged jet #it{p}_{T} (GeV/#it{c})"},
                 {40, -1.0, 1.0}});

    spectra.add("ptetaJetChSelected_highptjettrigger", "pT of selected high pT charged jet trigger vs eta", HistType::kTH2F,
                {{150, 0., +150., "charged jet #it{p}_{T} (GeV/#it{c})"},
                 {40, -1.0, 1.0}});

    auto scalers{std::get<std::shared_ptr<TH1>>(spectra.add(
      "fProcessedEvents", ";;Number of filtered events", HistType::kTH1F,
      {{kHighPtObjects, -0.5, kHighPtObjects - 0.5}}))};
    for (uint32_t iS{1}; iS <= highPtObjectsNames.size(); ++iS) {
      scalers->GetXaxis()->SetBinLabel(iS, highPtObjectsNames[iS - 1].data());
    }
  }

  // declare filters on tracks
  // Filter collisionFilter = nabs(aod::jcollision::posZ) < cfgVertexCut;

  Filter jetRadiusSelection = o2::aod::jet::r == nround(cfgJetR.node() * 100.0f);
  using filteredJets = o2::soa::Filtered<o2::aod::ChargedJets>;

  void process(aod::JCollision const& collision, filteredJets const& jets)
  {
    // collision process loop
    bool keepEvent[kHighPtObjects]{false};
    spectra.fill(HIST("fCollZpos"), collision.posZ());

    for (const auto& jet : jets) { // jets are ordered by pT

      if (jet.pt() >= jetPtLowThreshold) {
        spectra.fill(HIST("ptphiJetChSelected_lowptjettrigger"), jet.pt(), jet.phi()); // charged jet pT vs phi
        spectra.fill(HIST("ptetaJetChSelected_lowptjettrigger"), jet.pt(), jet.eta()); // charged jet pT vs eta
        keepEvent[kJetChLowPt] = true;
      }
      if (jet.pt() >= jetPtHighThreshold) {
        spectra.fill(HIST("ptphiJetChSelected_highptjettrigger"), jet.pt(), jet.phi()); // charged jet pT vs phi
        spectra.fill(HIST("ptetaJetChSelected_highptjettrigger"), jet.pt(), jet.eta()); // charged jet pT vs eta
        keepEvent[kJetChHighPt] = true;
      }
      break; // only looks at the highest pT jet in the event
    }

    for (int iDecision{0}; iDecision < kHighPtObjects; ++iDecision) {
      if (keepEvent[iDecision]) {
        spectra.fill(HIST("fProcessedEvents"), iDecision);
      }
    }
    tags(keepEvent[kJetChLowPt], keepEvent[kJetChHighPt]);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{

  return WorkflowSpec{adaptAnalysisTask<jetFilter>(cfg)};
}
