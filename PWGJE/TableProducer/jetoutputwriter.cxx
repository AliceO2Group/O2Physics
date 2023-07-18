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

/// \file jetoutputwriter.cxx
/// \brief Task to skim jet framework tables (JetCollisions, JetTracks, JetClusters, ...)
/// while adjusting indices accordingly
///
/// \author Jochen Klein <jochen.klein@cern.ch>
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/runDataProcessing.h"

#include "PWGJE/DataModel/Jet.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct JetOutputWriter {

  Produces<JetCollisionsMarked> jetCollisionsOutput;
  Produces<JetCollisionsMarked> jetTracksOutput;

  bool acceptCollision(JetCollision const &collision) {
    return true;
  }

  template <typename Jets>
  void process(JetCollision const &collision, Jets const& jets) {
    bool keepCollision = false;
    for (const auto &jet : jets) {
      // TODO: replace with meangingful condition
      keepCollision = jet.pt() > 10.;
    }
  }

  // add specializations for all jets

  // how can we ensure running order?

  // how to template for all possible jets?
  void process(JetCollisions const& collision, JetTracks const& tracks, Jets const &jets) {
    if (acceptCollision(collision)) {
      jetCollisionsOutput(collision);

      for (const auto &track : tracks) {
        jetTracksOutput(jetCollisionsOutput.lastIndex(), ...);
      }

      // do we want to modify the jets table at this stage?
      // or do we filter the jets at the level of the jet finder already?
      // the latter would avoid having to rewrite the jet tables
      // but we anyway have to rewrite the consituent tables
      for (const auto &jet : jets) {
        jetsOutput();
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;

  tasks.emplace_back(adaptAnalysisTask<JetOutputWriter>(cfgc, TaskName{"jet-output-writer"}));

  return WorkflowSpec{tasks};
}