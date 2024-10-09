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

// task to produce a self contained data format for trigger tables in jet analyses from the full AO2D
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#include <vector>
#include <string>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/O2DatabasePDGPlugin.h"

#include "EventFiltering/filterTables.h"

#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct JetDerivedDataTriggerProducerTask {
  Produces<aod::JChTrigSels> jChargedTriggerSelsTable;
  Produces<aod::JFullTrigSels> jFullTriggerSelsTable;
  Produces<aod::JChHFTrigSels> jChargedHFTriggerSelsTable;

  void init(InitContext const&)
  {
  }

  void processChargedJetTriggers(soa::Join<aod::Collisions, aod::JetFilters>::iterator const& collision)
  {
    jChargedTriggerSelsTable(jetderiveddatautilities::setChargedTriggerSelectionBit(collision));
  }
  PROCESS_SWITCH(JetDerivedDataTriggerProducerTask, processChargedJetTriggers, "produces derived charged trigger table", false);

  void processFullJetTriggers(soa::Join<aod::Collisions, aod::FullJetFilters>::iterator const& collision)
  {
    jFullTriggerSelsTable(jetderiveddatautilities::setFullTriggerSelectionBit(collision));
  }
  PROCESS_SWITCH(JetDerivedDataTriggerProducerTask, processFullJetTriggers, "produces derived full trigger table", false);

  void processChargedHFJetTriggers(soa::Join<aod::Collisions, aod::JetHFFilters>::iterator const& collision)
  {
    jChargedHFTriggerSelsTable(jetderiveddatautilities::setChargedHFTriggerSelectionBit(collision));
  }
  PROCESS_SWITCH(JetDerivedDataTriggerProducerTask, processChargedHFJetTriggers, "produces derived charged hf trigger table", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<JetDerivedDataTriggerProducerTask>(cfgc, TaskName{"jet-deriveddata-trigger-producer"})};
}
