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
/// \file lumiStabilityLightIons.cxx
/// \brief Analysis over BCs to study the luminosity stability along time
///
/// \author Nicolas Strangmann (nicolas.strangmann@cern.ch) - Goethe University Frankfurt
/// \author Stefanie Mrozinski (stefanie.mrozinski@cern.ch) - Goethe University Frankfurt
/// \author Lorenzo Mattei (lorenzo.mattei@cern.ch) - Turin University

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/Core/MetadataHelper.h"
#include "Common/DataModel/EventSelection.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/LHCConstants.h>
#include <DataFormatsParameters/AggregatedRunInfo.h>
#include <DataFormatsParameters/GRPLHCIFData.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TH2.h>

#include <array>
#include <bitset>
#include <cstdint>
#include <cstdlib>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

o2::common::core::MetadataHelper metadataInfo;

namespace o2::aod
{
namespace myBc_aod
{
DECLARE_SOA_COLUMN(Timestamp, timestamp, uint64_t);
DECLARE_SOA_COLUMN(BCid, bcId, int);
DECLARE_SOA_COLUMN(TimeZNA, timeZNA, float);
DECLARE_SOA_COLUMN(TimeZNC, timeZNC, float);
DECLARE_SOA_COLUMN(AmplitudeZNA, amplitudeZNA, float);
DECLARE_SOA_COLUMN(AmplitudeZNC, amplitudeZNC, float);
} // namespace myBc_aod
DECLARE_SOA_TABLE(MyBCaod, "AOD", "MYBCAOD",
                  myBc_aod::Timestamp,
                  myBc_aod::BCid,
                  myBc_aod::TimeZNA,
                  myBc_aod::TimeZNC,
                  myBc_aod::AmplitudeZNA,
                  myBc_aod::AmplitudeZNC);
} // namespace o2::aod

using MyBCs = soa::Join<aod::BCs, aod::BcSels, aod::Timestamps, aod::Run3MatchedToBCSparse>;

struct LumiZdcQa {
  Produces<aod::MyBCaod> BCaod;

  Configurable<bool> cfgRequireZDCTriggerForZDCQA{"cfgRequireZDCTriggerForZDCQA", true, "Require ZDC trigger (1ZNC) for filling QA histograms"};
  Configurable<bool> cfgRequireTVXTriggerForZDCQA{"cfgRequireTVXTriggerForZDCQA", true, "Require FT0 vertex trigger (MTVX) for filling ZDC QA histograms"};
  Configurable<bool> cfgRequireZEDTriggerForZDCQA{"cfgRequireZEDTriggerForZDCQA", true, "Require ZED trigger (1ZNC||1ZNA) for filling QA histograms"};

  Configurable<bool> cfgFillBCao2d{"cfgFillBCao2d", false, "Fill BC ao2d with timestamps and ZDC times"};
  Configurable<uint64_t> cfgTstampStartFillingBCao2d{"cfgTstampStartFillingBCao2d", 0, "Minimum value of timestamp for output bc ao2d to be filled"};
  Configurable<uint64_t> cfgTstampEndFillingBCao2d{"cfgTstampEndFillingBCao2d", 0, "Maximum value of timestamp for output bc ao2d to be filled"};

  const int nBCsPerOrbit = 3564;

  HistogramRegistry mHistManager{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  void init(InitContext&)
  {
    AxisSpec zdcTimeAxis{200, -50., 50.};
    mHistManager.add("ZDCQA/BCHasZDC", "Does the BC have ZDC?;BC has ZDC;Has ZNC according to CTP;#bf{#it{N}_{BC}}", HistType::kTH2D, {{2, -0.5, 1.5}, {2, -0.5, 1.5}});
    mHistManager.get<TH2>(HIST("ZDCQA/BCHasZDC")).get()->GetYaxis()->SetBinLabel(1, "No CTP trigger");
    mHistManager.get<TH2>(HIST("ZDCQA/BCHasZDC")).get()->GetYaxis()->SetBinLabel(2, "CTP triggered");
    mHistManager.get<TH2>(HIST("ZDCQA/BCHasZDC")).get()->GetXaxis()->SetBinLabel(1, "No found ZDC");
    mHistManager.get<TH2>(HIST("ZDCQA/BCHasZDC")).get()->GetXaxis()->SetBinLabel(2, "Good ZDC");
    mHistManager.add("ZDCQA/ZNCTimeVsEnergy", "ZDC properties in BCs with found ZDC;Energy;#bf{ZNC arrival time (ns)};#bf{#it{N}_{BC}}", HistType::kTH2D, {{1501, -10, 1.5E4}, zdcTimeAxis});
    mHistManager.add("ZDCQA/ZDCTimes", "Correlation between ZNA and ZNC timing;#bf{ZNC arrival time (ns)};#bf{ZNA arrival time (ns)}", HistType::kTH2D, {zdcTimeAxis, zdcTimeAxis});
    mHistManager.add("ZDCQA/ZNATime", "Time of the ZNA signal;#bf{ZNA arrival time (ns)};#bf{#it{N}_{BC}}", HistType::kTH1D, {zdcTimeAxis});
    mHistManager.add("ZDCQA/ZNCTime", "Time of the ZNC signal;#bf{ZNC arrival time (ns)};#bf{#it{N}_{BC}}", HistType::kTH1D, {zdcTimeAxis});
  }

  void process(MyBCs const& bcs, aod::Zdcs const&)
  {
    const int maxTimeZDC = 50;
    const float dummyZDCTime = 42.f;

    for (const auto& bc : bcs) {
      std::bitset<64> ctpInputMask(bc.inputMask());

      if (cfgRequireTVXTriggerForZDCQA && !(ctpInputMask.test(2))) {
        continue;
      }
      if (cfgRequireZDCTriggerForZDCQA && !(ctpInputMask.test(25))) {
        continue;
      }
      if (cfgRequireZEDTriggerForZDCQA && !(ctpInputMask.test(24))) {
        continue;
      }

      bool zdcHit = !bc.has_zdc() ? 0 : ((bc.zdc().energyCommonZNC() > -1 && std::abs(bc.zdc().timeZNC()) < 1E5) ? 1 : 0);
      mHistManager.fill(HIST("ZDCQA/BCHasZDC"), zdcHit, ctpInputMask.test(25) ? 1 : 0);

      if (!bc.has_zdc()) {
        continue;
      }

      mHistManager.fill(HIST("ZDCQA/ZNCTimeVsEnergy"),
                        bc.zdc().energyCommonZNC() > -1 ? bc.zdc().energyCommonZNC() : -1,
                        std::abs(bc.zdc().timeZNC()) < maxTimeZDC ? bc.zdc().timeZNC() : dummyZDCTime);

      float timeZNA = bc.zdc().timeZNA();
      float timeZNC = bc.zdc().timeZNC();

      if (std::abs(timeZNA) > maxTimeZDC) {
        timeZNA = dummyZDCTime;
        mHistManager.fill(HIST("ZDCQA/ZNCTime"), timeZNC);
      }
      if (std::abs(timeZNC) > maxTimeZDC) {
        timeZNC = dummyZDCTime;
        if (timeZNA != dummyZDCTime) {
          mHistManager.fill(HIST("ZDCQA/ZNATime"), timeZNA);
        }
      }

      mHistManager.fill(HIST("ZDCQA/ZDCTimes"), timeZNA, timeZNC);

      uint64_t timestamp = bc.timestamp();
      int64_t globalBC = bc.globalBC();
      int localBC = globalBC % nBCsPerOrbit;
      float amplitudeZNA = bc.zdc().amplitudeZNA();
      float amplitudeZNC = bc.zdc().amplitudeZNC();

      if (cfgFillBCao2d && timestamp >= cfgTstampStartFillingBCao2d && timestamp <= cfgTstampEndFillingBCao2d) {
        BCaod(timestamp, localBC, timeZNA, timeZNC, amplitudeZNA, amplitudeZNC);
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  metadataInfo.initMetadata(cfgc);
  return WorkflowSpec{adaptAnalysisTask<LumiZdcQa>(cfgc)};
}
