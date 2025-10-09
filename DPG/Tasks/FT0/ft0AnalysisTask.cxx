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

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/DataModel/EventSelection.h"

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TAxis.h>
#include <TH2.h>

#include <bitset>
#include <cstddef>
#include <map>
#include <memory>
#include <numeric>
#include <string>

using namespace o2;
using namespace o2::framework;
using namespace o2::aod::evsel;

struct ft0AnalysisTask {
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  constexpr static int sNtriggers = 5; // Number of produced triggers
  enum ETriggers {
    kNoZDC,
    kZNA,
    kZNC,
    kZNAorZNC,
    kZNAandZNC
  };
  using TriggerWord_t = std::bitset<sNtriggers>;
  void init(InitContext&)
  {
    const AxisSpec axisSumAmpA{2300, 0., 230000., "SumAmpA [ADC]"};
    const AxisSpec axisSumAmpC{2300, 0., 230000., "SumAmpC [ADC]"};
    const AxisSpec axisSumAmp{2500, 0., 250000., "SumAmp [ADC]"};

    const AxisSpec axisLowSumAmpA{2000, 0., 2000., "SumAmpA [ADC]"};
    const AxisSpec axisLowSumAmpC{2000, 0., 2000., "SumAmpC [ADC]"};
    const AxisSpec axisLowSumAmp{2000, 0., 2000., "SumAmp [ADC]"};

    std::map<unsigned int, std::string> mapTriggers = {
      {kNoZDC, "No ZDC"},
      {kZNA, "ZNA"},
      {kZNC, "ZNC"},
      {kZNAorZNC, "ZNA or ZNC"},
      {kZNAandZNC, "ZNA and ZNC"}};
    const AxisSpec axisTriggers{static_cast<int>(mapTriggers.size()), 0., static_cast<double>(mapTriggers.size()), "Triggers"};
    int axisPerTrg = 1;
    // temporary solution
    auto makeBinLabels = [&mapLabels = mapTriggers, &axisPos = axisPerTrg](auto histPtrVariant) {
      auto hist = std::get<std::shared_ptr<TH2>>(histPtrVariant);
      TAxis* axis = nullptr;
      if (axisPos == 0) {
        axis = hist->GetXaxis();
      } else if (axisPos == 1) {
        axis = hist->GetYaxis();
      } else if (axisPos == 2) {
        axis = hist->GetZaxis();
      } else {
        return;
      }
      for (const auto& [binIdx, label] : mapLabels) {
        axis->SetBinLabel(binIdx + 1, label.c_str());
      }
    };

    makeBinLabels(histos.add("hSumAmpA_perTriggers", "Sum amp FT0A per triggers", kTH2F, {axisSumAmpA, axisTriggers}));
    makeBinLabels(histos.add("hSumAmpC_perTriggers", "Sum amp FT0C per triggers", kTH2F, {axisSumAmpC, axisTriggers}));
    makeBinLabels(histos.add("hSumAmp_perTriggers", "Sum amp FT0 per triggers", kTH2F, {axisSumAmp, axisTriggers}));

    makeBinLabels(histos.add("hLowSumAmpA_perTriggers", "Sum amp FT0A per triggers", kTH2F, {axisLowSumAmpA, axisTriggers}));
    makeBinLabels(histos.add("hLowSumAmpC_perTriggers", "Sum amp FT0C per triggers", kTH2F, {axisLowSumAmpC, axisTriggers}));
    makeBinLabels(histos.add("hLowSumAmp_perTriggers", "Sum amp FT0 per triggers", kTH2F, {axisLowSumAmp, axisTriggers}));

    makeBinLabels(histos.add("hSumAmpA_perNotTriggers", "Sum amp FT0A per NOT triggers", kTH2F, {axisSumAmpA, axisTriggers}));
    makeBinLabels(histos.add("hSumAmpC_perNotTriggers", "Sum amp FT0C per NOT triggers", kTH2F, {axisSumAmpC, axisTriggers}));
    makeBinLabels(histos.add("hSumAmp_perNotTriggers", "Sum amp FT0 per NOT triggers", kTH2F, {axisSumAmp, axisTriggers}));

    makeBinLabels(histos.add("hLowSumAmpA_perNotTriggers", "Sum amp FT0A per NOT triggers", kTH2F, {axisLowSumAmpA, axisTriggers}));
    makeBinLabels(histos.add("hLowSumAmpC_perNotTriggers", "Sum amp FT0C per NOT triggers", kTH2F, {axisLowSumAmpC, axisTriggers}));
    makeBinLabels(histos.add("hLowSumAmp_perNotTriggers", "Sum amp FT0 per NOT triggers", kTH2F, {axisLowSumAmp, axisTriggers}));
  }

  void process(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, aod::FT0s const&)
  {
    if (!collision.has_foundFT0()) { // Analysis is fully based on FT0 presence
      return;
    }
    const auto& ft0 = collision.foundFT0();
    TriggerWord_t trgWord = (((!collision.has_foundZDC()) << kNoZDC) |
                             (collision.selection_bit(kIsBBZNA) << kZNA) |
                             (collision.selection_bit(kIsBBZNC) << kZNC) |
                             ((collision.selection_bit(kIsBBZNA) || collision.selection_bit(kIsBBZNC)) << kZNAorZNC) |
                             ((collision.selection_bit(kIsBBZNA) && collision.selection_bit(kIsBBZNC)) << kZNAandZNC));
    auto sumCalc = [](auto&& sum, auto&& curr) { return sum + (curr > 0 ? curr : 0); };
    constexpr float dummyValue = -5.e3;
    auto checkDummy = [&dummyValue](auto&& val) { return val > 0 ? val : dummyValue; }; // zero values are moves to dummy
    const float sumAmpA = checkDummy(std::accumulate(ft0.amplitudeA().begin(), ft0.amplitudeA().end(), 0.f, sumCalc));
    const float sumAmpC = checkDummy(std::accumulate(ft0.amplitudeC().begin(), ft0.amplitudeC().end(), 0.f, sumCalc));
    const float sumAmp = (sumAmpA > dummyValue && sumAmpC > dummyValue) ? (sumAmpA + sumAmpC) : dummyValue;

    for (std::size_t iTrgBit = 0; iTrgBit < trgWord.size(); iTrgBit++) {
      if (trgWord.test(iTrgBit)) {
        histos.fill(HIST("hSumAmpA_perTriggers"), sumAmpA, iTrgBit);
        histos.fill(HIST("hLowSumAmpA_perTriggers"), sumAmpA, iTrgBit);
        histos.fill(HIST("hSumAmpC_perTriggers"), sumAmpC, iTrgBit);
        histos.fill(HIST("hLowSumAmpC_perTriggers"), sumAmpC, iTrgBit);
        histos.fill(HIST("hSumAmp_perTriggers"), sumAmp, iTrgBit);
        histos.fill(HIST("hLowSumAmp_perTriggers"), sumAmp, iTrgBit);
      } else { // NOT case
        histos.fill(HIST("hSumAmpA_perNotTriggers"), sumAmpA, iTrgBit);
        histos.fill(HIST("hLowSumAmpA_perNotTriggers"), sumAmpA, iTrgBit);
        histos.fill(HIST("hSumAmpC_perNotTriggers"), sumAmpC, iTrgBit);
        histos.fill(HIST("hLowSumAmpC_perNotTriggers"), sumAmpC, iTrgBit);
        histos.fill(HIST("hSumAmp_perNotTriggers"), sumAmp, iTrgBit);
        histos.fill(HIST("hLowSumAmp_perNotTriggers"), sumAmp, iTrgBit);
      }
    }
  } // end of processCollsions()

  PROCESS_SWITCH(ft0AnalysisTask, process, "collision analysis", true);
}; // end of struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<ft0AnalysisTask>(cfgc, TaskName{"ft0-task"})};
}
