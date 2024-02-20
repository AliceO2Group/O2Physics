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
#include "Framework/ConfigParamSpec.h"

using namespace o2;
using namespace o2::framework;

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"
#include "Common/DataModel/EventSelection.h"
#include <CommonConstants/LHCConstants.h>

struct ZdcRateTask {
  SliceCache cache;
  using BCsRun3 = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;
  HistogramRegistry registry{"Histos"};
  void init(InitContext& context)
  {
    registry.add("rate", "rate;ZDC rate (Hz)", kTH1F, {{10000, 1.f, 1000001.f}});
    registry.add("energyCommon", "energyCommon;energyCommonZNA;energyCommonZNC", kTH2F, {{400, -0.5f, 399.5f}, {400, -0.5f, 399.5f}});
  }
  double zdcRate = 0;
  double zdcWindow = 0;
  static constexpr double windowDuration = (o2::constants::lhc::LHCOrbitNS * 1e-9) * 32;

  void process(BCsRun3 const& bcs,
               aod::Zdcs const& zdcs)
  {
    int lastBc = 0;
    for (const auto& bc : bcs) {
      const int bcNumberInWindow = bc.globalBC() % (3564 * 32);
      if (lastBc < bcNumberInWindow) {
        lastBc = bcNumberInWindow;
        LOG(info) << "BC " << bc.globalBC() << " in window " << bcNumberInWindow << " rate is " << zdcRate << " / " << windowDuration << " = " << zdcRate / windowDuration << " Hz";
        zdcRate /= windowDuration;
        registry.fill(HIST("rate"), zdcRate / 20);
        zdcRate = 0;
      }
      if (!bc.has_zdc()) {
        continue;
      }
      registry.fill(HIST("energyCommon"), bc.zdc().energyCommonZNA(), bc.zdc().energyCommonZNC());
      if (bc.zdc().energyCommonZNA() > 0.f || bc.zdc().energyCommonZNC() > 0.f) {
        zdcRate += 1.0;
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<ZdcRateTask>(cfgc)}; }
