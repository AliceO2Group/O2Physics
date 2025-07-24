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
/// \brief A task to create a flat tree for ML model training as an input
///        for ML tutorial task
/// \author
/// \since

/*
 * command to run:
 * o2-analysis-track-propagation --aod-writer-keep="AOD/TRTR/0" |\
 * o2-analysis-timestamp --aod-writer-keep="AOD/TRTR/0" |\
 * o2-analysis-multiplicity-table --aod-writer-keep="AOD/TRTR/0" |\
 * o2-analysis-event-selection --aod-writer-keep="AOD/TRTR/0" |\
 * o2-analysistutorial-tree-creator --aod-writer-keep="AOD/TRTR/0"
 */

/// Note that output in AnalysisResults_trees.root will be split
/// by dataframe, use o2-aod-merger to merge them into one`
/// use train_model.pynb to train the model

// This workflow is used to create a flat tree for model training
// Use o2-aod-merger to combine dataframes in output AnalysisResults_trees.root

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/Multiplicity.h"
#include "TrainingTree.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct CreateTree {
  Configurable<float> centralEtaCut{"centralEtaCut", 0.8, "central eta limit"};

  Produces<aod::TrainingTree> tt;
  Filter centralTracks = nabs(aod::track::eta) < centralEtaCut;

  void processRun2(soa::Join<aod::Collisions, aod::Mults>::iterator const& collision, soa::Filtered<aod::Tracks> const& tracks)
  {
    tt(collision.posZ(), collision.posX(), collision.posY(), analysis::meanPt(tracks), tracks.size(), collision.multFT0M(), collision.multFV0M());
  }

  PROCESS_SWITCH(CreateTree, processRun2, "Use Run 2 parameters", false);

  void processRun3(soa::Join<aod::Collisions, aod::Mults>::iterator const& collision, soa::Filtered<aod::Tracks> const& tracks)
  {
    tt(collision.posZ(), collision.posX(), collision.posY(), analysis::meanPt(tracks), tracks.size(), collision.multFT0M(), collision.multFV0A());
  }

  PROCESS_SWITCH(CreateTree, processRun3, "Use Run 3 parameters", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<CreateTree>(cfgc)};
}
