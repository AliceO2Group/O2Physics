// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
#include "dndeta.h"
#include <cmath>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::pwgmm::multiplicity;

void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  ConfigParamSpec optionDoMC{"doMC", VariantType::Bool, false, {"Use MC information"}};
  workflowOptions.push_back(optionDoMC);
}
// always should be after customize() function
#include "Framework/runDataProcessing.h"

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  const bool doMC = cfgc.options().get<bool>("doMC");
  if (!doMC) {
    return WorkflowSpec{adaptAnalysisTask<PseudorapidityDensity<aod::track::TrackTypeEnum::Run2Tracklet>>(cfgc, TaskName{"pseudorapidity-density"})};
  }
  return WorkflowSpec{adaptAnalysisTask<PseudorapidityDensity<aod::track::TrackTypeEnum::Run2Tracklet>>(cfgc, TaskName{"pseudorapidity-density"}, SetDefaultProcesses{{{"processGen", true}, {"processMatching", true}}})};
}
