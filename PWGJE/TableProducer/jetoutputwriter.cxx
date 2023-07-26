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
#include "PWGJE/DataModel/JetReducedData.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct JetOutputWriter {

  Configurable<float> jetPtMin{"jetPtMin", 1.f, "Minimum jet pt to accept event"};

  Produces<o2::aod::StoredJCollisions> storedJetCollisions;
  Produces<o2::aod::StoredJTracks> storedJetTracks;

  std::vector<bool> collFlag;

  bool acceptCollision(o2::aod::JCollision const& collision)
  {
    return true;
  }

  // explicit process function used to reset acceptance flags
  void process(o2::aod::JCollisions const& collisions)
  {
    collFlag.reserve(collisions.size());
    std::fill(collFlag.begin(), collFlag.end(), false);
  }

  template <typename Jets>
  void processJets(Jets& jets)
  {
    for (const auto& jet : jets) {
      if (jet.pt() > jetPtMin) {
        collFlag[jet.collisionId()] = true;
      }
    }
  }
// TODO: replace with PROCESS_SWITCH_FULL when available
#define PROCESS_SWITCH_JKL(_Class_, _Method_, _Name_, _Help_, _Default_) \
  decltype(ProcessConfigurable{&_Class_ ::_Method_, #_Name_, _Default_, _Help_}) do##_Name_ = ProcessConfigurable{&_Class_ ::_Method_, #_Name_, _Default_, _Help_};
  PROCESS_SWITCH_JKL(JetOutputWriter, processJets<o2::aod::ChargedJets>, processChargedJets, "process charged jets", true);
  PROCESS_SWITCH_JKL(JetOutputWriter, processJets<o2::aod::NeutralJets>, processNeutralJets, "process neutral jets", true);
  PROCESS_SWITCH_JKL(JetOutputWriter, processJets<o2::aod::D0ChargedJets>, processD0ChargedJets, "process D0 charged jets", true);
  PROCESS_SWITCH_JKL(JetOutputWriter, processJets<o2::aod::LcChargedJets>, processLcChargedJets, "process Lc charged jets", true);

  void processCollisions(o2::aod::JCollision const& collision, o2::aod::JTracks const& tracks)
  {
    if (collFlag[collision.globalIndex()]) {
      storedJetCollisions(collision.posZ(), collision.eventSel());

      for (const auto &track : tracks) {
        storedJetTracks(storedJetCollisions.lastIndex(), track.pt(), track.eta(), track.phi(), track.energy(), track.trackSel());
      }
    }
  }
  // process switch for output writing must be last
  // to run after all jet selections
  PROCESS_SWITCH(JetOutputWriter, processCollisions, "write collisions and tracks to output tables", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;

  tasks.emplace_back(adaptAnalysisTask<JetOutputWriter>(cfgc, TaskName{"jet-output-writer"}));

  return WorkflowSpec{tasks};
}