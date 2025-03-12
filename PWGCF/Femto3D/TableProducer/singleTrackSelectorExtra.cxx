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
//
/// \brief create a table applying some basic cuts on the ITS and DCA.
/// \author Sofia Tomassini, Gleb Romanenko, Nicolò Jacazio
/// \since 31 May 2023

// this task produces a dummy "SingleCollExtras" table that is now required in the analysis tasks. Needed to have a compatibility with old der. data

#include <fairlogger/Logger.h>
#include <Framework/AnalysisDataModel.h>

#include "PWGCF/Femto3D/DataModel/singletrackselector.h"

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::track;
using namespace o2::aod;
//::singletrackselector; // the namespace defined in .h

struct singleTrackSelectorDummy {

  Produces<o2::aod::SingleCollExtras> tableRowCollExtra;

  void process(aod::SingleCollSels::iterator const&)
  {
    uint64_t selection = 0;
    tableRowCollExtra(selection,
                      0.0,
                      0);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<singleTrackSelectorDummy>(cfgc)};
}
