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
#include "DataModel/JEDerived.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/Core/JetFinder.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#include "Framework/runDataProcessing.h"

struct JetProviderTask {

  Produces<aod::JEJets> outputJets;
  Produces<aod::JEConstituents> outputConstituents;

  Configurable<float> jetPtMin{"jetPtMin", 0.0, "minimum jet pT cut"};
  Configurable<bool> keepConstituents{"keepConstituents", true, "Constituent table is filled"};
  Configurable<bool> DoConstSub{"DoConstSub", false, "do constituent subtraction"};
  Filter jetCuts = aod::jet::pt > jetPtMin;

  void process(soa::Filtered<soa::Join<aod::ChargedJets, aod::ChargedJetConstituents>>::iterator const& jet,
               aod::Tracks const& tracks)
  {
    outputJets(jet.pt(), jet.eta(), jet.phi(), jet.energy(), jet.mass(), jet.area());
    if (keepConstituents) {
      outputConstituents.reserve(jet.tracks().size());
      for (const auto& constituent : jet.tracks_as<aod::Tracks>()) {
        outputConstituents(outputJets.lastIndex(), constituent.pt(), constituent.eta(), constituent.phi());
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<JetProviderTask>(cfgc, TaskName{"jet-task-skim-provider"})};
  return workflow;
}
