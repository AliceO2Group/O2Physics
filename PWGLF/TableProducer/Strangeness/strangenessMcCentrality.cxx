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
/// \file   mcCentrality.cxx
/// \author Romain Schotter romain.schotter@cern.ch
/// \brief  Task to produce the MC centrality table for strangeness derived data
///

#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/Utils/mcCentralityModule.h"

#include <CCDB/BasicCCDBManager.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <string>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

/// Task to produce the response table
struct StrangenessMcCentrality {
  // Input parameters
  o2::framework::Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  o2::pwglf::mccentrality::products products;
  o2::pwglf::mccentrality::coreConfigurables baseOpts;
  o2::pwglf::mccentrality::BuilderModule mcCentralityBuilderModule;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext& initContext)
  {
    // Set up the CCDB
    ccdb->setURL(ccdbUrl.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    mcCentralityBuilderModule.init(baseOpts, histos, initContext);
  }

  // Full tables (independent on central calibrations)
  void process(aod::StraMCCollMults const& mcCollisions,
               aod::StraStamps const& bcs)
  {
    mcCentralityBuilderModule.dataProcess(ccdb, histos, bcs, mcCollisions, products);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<StrangenessMcCentrality>(cfgc)};
}
