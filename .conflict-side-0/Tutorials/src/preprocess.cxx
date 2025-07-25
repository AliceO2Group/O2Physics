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
/// \brief Preprocess task creates intermediate tables for Consume task
/// and automatically adjust its configuration by inspecting the requirements
/// of the Consume task
/// \author
/// \since

#include "IntermediateTables.h"

#include "Framework/AnalysisTask.h"
#include "Framework/ConfigParamSpec.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;

struct PreProcess {
  Produces<aod::Decisions> decisions;

  Configurable<float> ptCutOff{"ptcutoff", 2.f, "Pt cut off"};
  Configurable<float> etaRange{"etarange", 0.8f, "eta range"};

  void init(InitContext& context)
  {
    auto& rw = context.services().get<RunningWorkflowInfo const>();
    auto consumer = std::find_if(rw.devices.begin(), rw.devices.end(), [](auto const& device) { return device.name.find("consume") != std::string::npos; });
    if (consumer != rw.devices.end()) {
      auto locate = std::find_if(consumer->inputs.begin(), consumer->inputs.end(), [](auto const& i) { return i.matcher.binding == "Decisions"; });
      if (locate != consumer->inputs.end()) {
        auto ptoption = std::find_if(consumer->options.begin(), consumer->options.end(), [](auto const& o) { return o.name == "ptcutoff"; });
        if (ptoption != consumer->options.end()) {
          auto value = ptoption->defaultValue.get<float>();
          LOG(info) << "Resetting option pt cut off from " << static_cast<float>(ptCutOff) << " to " << value;
          ptCutOff.value = value;
        }
        auto etaoption = std::find_if(consumer->options.begin(), consumer->options.end(), [](auto const& o) { return o.name == "etarange"; });
        if (etaoption != consumer->options.end()) {
          auto value = etaoption->defaultValue.get<float>();
          LOG(info) << "Resetting option eta range from " << static_cast<float>(etaRange) << " to " << value;
          etaRange.value = value;
        }
      }
    }
  }

  void process(aod::Tracks const& tracks)
  {
    for (auto& track : tracks) {
      decisions(track.pt() < ptCutOff, std::abs(track.eta()) < etaRange);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& context)
{
  return {
    adaptAnalysisTask<PreProcess>(context),
  };
}
