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

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/ctpRateFetcher.h"
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
#include <TString.h>

#include <array>
#include <bitset>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

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

struct LumiStabilityLightIons {
  Produces<aod::MyBCaod> BCaod;

  Configurable<bool> cfgDoFT0Vtx{"cfgDoFT0Vtx", true, "Create and fill histograms for the FT0 vertex trigger"};
  Configurable<bool> cfgDoFT0CE{"cfgDoFT0CE", true, "Create and fill histograms for the FT0 centrality trigger"};
  Configurable<bool> cfgDoFDD{"cfgDoFDD", true, "Create and fill histograms for the FDD trigger"};
  Configurable<bool> cfgDo1ZNC{"cfgDo1ZNC", true, "Create and fill histograms for the 1ZNC trigger"};

  Configurable<bool> cfgDoBCA{"cfgDoBCA", true, "Create and fill histograms for the BCs of type A"};
  Configurable<bool> cfgDoBCB{"cfgDoBCB", true, "Create and fill histograms for the BCs of type B"};
  Configurable<bool> cfgDoBCC{"cfgDoBCC", true, "Create and fill histograms for the BCs of type C"};
  Configurable<bool> cfgDoBCE{"cfgDoBCE", true, "Create and fill histograms for the BCs of type E"};
  Configurable<bool> cfgDoBCL{"cfgDoBCL", true, "Create and fill histograms for leading BCs of type B (non-B BCs before)"};
  Configurable<bool> cfgDoBCLE{"cfgDoBCLE", true, "Create and fill histograms for leading BCs of type B (strictly empty BCs before)"};
  Configurable<bool> cfgDoBCNL{"cfgDoBCNL", true, "Create and fill histograms for non-leading BCs of type B (complement of BCL)"};
  Configurable<bool> cfgDoBCNLE{"cfgDoBCNLE", true, "Create and fill histograms for non-leading BCs of type B (complement of BCLE)"};
  Configurable<bool> cfgDoBCSLFDD{"cfgDoBCSLFDD", true, "Create and fill histograms for super-leading BCs w.r.t. FDD activity"};
  Configurable<bool> cfgDoBCSLFT0{"cfgDoBCSLFT0", true, "Create and fill histograms for super-leading BCs w.r.t. FT0 activity"};
  Configurable<bool> cfgDoBCNSLFDD{"cfgDoBCNSLFDD", true, "Create and fill histograms for non-super-leading BCs w.r.t. FDD activity"};
  Configurable<bool> cfgDoBCNSLFT0{"cfgDoBCNSLFT0", true, "Create and fill histograms for non-super-leading BCs w.r.t. FT0 activity"};

  Configurable<bool> cfgRequireZDCTriggerForZDCQA{"cfgRequireZDCTriggerForZDCQA", true, "Require ZDC trigger (1ZNC) for filling QA histograms"};
  Configurable<bool> cfgRequireTVXTriggerForZDCQA{"cfgRequireTVXTriggerForZDCQA", true, "Require FT0 vertex trigger (MTVX) for filling ZDC QA histograms"};
  Configurable<bool> cfgRequireZEDTriggerForZDCQA{"cfgRequireZEDTriggerForZDCQA", true, "Require ZED trigger (1ZNC||1ZNA) for filling QA histograms"};

  Configurable<bool> cfgRequireNoT0ForSLBC{"cfgRequireNoT0ForSLBC", false, "Require no T0 signal for definition of super leading BC (otherwise only no FDD)"};

  Configurable<int> cfgEmptyBCsBeforeLeadingBC{"cfgEmptyBCsBeforeLeadingBC", 5, "Minimum number of non-B BCs before a BCL leading BC"};
  Configurable<int> cfgEmptyBCsBeforeLeadingBCLE{"cfgEmptyBCsBeforeLeadingBCLE", 5, "Minimum number of strictly empty (E-type) BCs before a BCLE leading BC"};
  Configurable<int> cfgBCsBeforeSuperLeading{"cfgBCsBeforeSuperLeading", 5, "Minimum number of BCs without FDD/FT0 activity before a super-leading BC"};

  Configurable<bool> cfgFillBCao2d{"cfgFillBCao2d", false, "Fill BC ao2d with timestamps and ZDC times"};
  Configurable<uint64_t> cfgTstampStartFillingBCao2d{"cfgTstampStartFillingBCao2d", 0, "Minimum value of timestamp for output bc ao2d to be filled"};
  Configurable<uint64_t> cfgTstampEndFillingBCao2d{"cfgTstampEndFillingBCao2d", 0, "Maximum value of timestamp for output bc ao2d to be filled"};

  Configurable<int> cfgBcShiftFDDForData2023{"cfgBcShiftFDDForData2023", 7, "Number of BCs to shift FDD, applied for 2023 data only"};

  std::bitset<o2::constants::lhc::LHCMaxBunches> beamPatternA, beamPatternC;
  std::bitset<o2::constants::lhc::LHCMaxBunches> bcPatternA, bcPatternC, bcPatternB, bcPatternE;
  std::bitset<o2::constants::lhc::LHCMaxBunches> bcPatternL;
  std::bitset<o2::constants::lhc::LHCMaxBunches> bcPatternLE;

  std::string strLPMProductionTag = "";
  const int nBCsPerOrbit = 3564;

  parameters::GRPLHCIFData* mLHCIFdata = nullptr;
  int mRunNumber = -1;
  bool isData23 = false;
  int mBcShiftFDD = 0;
  ctpRateFetcher mRateFetcher;

  HistogramRegistry mHistManager{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  static constexpr int nTriggers = 5;
  enum TriggerAliases {
    kAllBCs = 0,
    kFT0Vtx = 1,
    kFT0CE = 2,
    kFDD = 3,
    k1ZNC = 4
  };

  static constexpr int nBCCategories = 12;
  enum BCCategories {
    kBCA = 0,
    kBCB = 1,
    kBCC = 2,
    kBCE = 3,
    kBCL = 4,
    kBCLE = 5,
    kBCNL = 6,
    kBCNLE = 7,
    kBCSLFDD = 8,
    kBCSLFT0 = 9,
    kBCNSLFDD = 10,
    kBCNSLFT0 = 11
  };

  static constexpr std::string_view NBCsVsTimeHistNames[6][12] = {
    {"AllBCs/BC_A/nBCsVsTime", "AllBCs/BC_B/nBCsVsTime", "AllBCs/BC_C/nBCsVsTime", "AllBCs/BC_E/nBCsVsTime", "AllBCs/BC_L/nBCsVsTime", "AllBCs/BC_LE/nBCsVsTime", "AllBCs/BC_NL/nBCsVsTime", "AllBCs/BC_NLE/nBCsVsTime", "AllBCs/BC_SL_FDD/nBCsVsTime", "AllBCs/BC_SL_FT0/nBCsVsTime", "AllBCs/BC_NSL_FDD/nBCsVsTime", "AllBCs/BC_NSL_FT0/nBCsVsTime"},
    {"FT0VTx/BC_A/nBCsVsTime", "FT0VTx/BC_B/nBCsVsTime", "FT0VTx/BC_C/nBCsVsTime", "FT0VTx/BC_E/nBCsVsTime", "FT0VTx/BC_L/nBCsVsTime", "FT0VTx/BC_LE/nBCsVsTime", "FT0VTx/BC_NL/nBCsVsTime", "FT0VTx/BC_NLE/nBCsVsTime", "FT0VTx/BC_SL_FDD/nBCsVsTime", "FT0VTx/BC_SL_FT0/nBCsVsTime", "FT0VTx/BC_NSL_FDD/nBCsVsTime", "FT0VTx/BC_NSL_FT0/nBCsVsTime"},
    {"FT0CE/BC_A/nBCsVsTime", "FT0CE/BC_B/nBCsVsTime", "FT0CE/BC_C/nBCsVsTime", "FT0CE/BC_E/nBCsVsTime", "FT0CE/BC_L/nBCsVsTime", "FT0CE/BC_LE/nBCsVsTime", "FT0CE/BC_NL/nBCsVsTime", "FT0CE/BC_NLE/nBCsVsTime", "FT0CE/BC_SL_FDD/nBCsVsTime", "FT0CE/BC_SL_FT0/nBCsVsTime", "FT0CE/BC_NSL_FDD/nBCsVsTime", "FT0CE/BC_NSL_FT0/nBCsVsTime"},
    {"FDD/BC_A/nBCsVsTime", "FDD/BC_B/nBCsVsTime", "FDD/BC_C/nBCsVsTime", "FDD/BC_E/nBCsVsTime", "FDD/BC_L/nBCsVsTime", "FDD/BC_LE/nBCsVsTime", "FDD/BC_NL/nBCsVsTime", "FDD/BC_NLE/nBCsVsTime", "FDD/BC_SL_FDD/nBCsVsTime", "FDD/BC_SL_FT0/nBCsVsTime", "FDD/BC_NSL_FDD/nBCsVsTime", "FDD/BC_NSL_FT0/nBCsVsTime"},
    {"1ZNC/BC_A/nBCsVsTime", "1ZNC/BC_B/nBCsVsTime", "1ZNC/BC_C/nBCsVsTime", "1ZNC/BC_E/nBCsVsTime", "1ZNC/BC_L/nBCsVsTime", "1ZNC/BC_LE/nBCsVsTime", "1ZNC/BC_NL/nBCsVsTime", "1ZNC/BC_NLE/nBCsVsTime", "1ZNC/BC_SL_FDD/nBCsVsTime", "1ZNC/BC_SL_FT0/nBCsVsTime", "1ZNC/BC_NSL_FDD/nBCsVsTime", "1ZNC/BC_NSL_FT0/nBCsVsTime"}};

  static constexpr std::string_view NBCsVsBCIDHistNames[5][12] = {
    {"AllBCs/BC_A/nBCsVsBCID", "AllBCs/BC_B/nBCsVsBCID", "AllBCs/BC_C/nBCsVsBCID", "AllBCs/BC_E/nBCsVsBCID", "AllBCs/BC_L/nBCsVsBCID", "AllBCs/BC_LE/nBCsVsBCID", "AllBCs/BC_NL/nBCsVsBCID", "AllBCs/BC_NLE/nBCsVsBCID", "AllBCs/BC_SL_FDD/nBCsVsBCID", "AllBCs/BC_SL_FT0/nBCsVsBCID", "AllBCs/BC_NSL_FDD/nBCsVsBCID", "AllBCs/BC_NSL_FT0/nBCsVsBCID"},
    {"FT0VTx/BC_A/nBCsVsBCID", "FT0VTx/BC_B/nBCsVsBCID", "FT0VTx/BC_C/nBCsVsBCID", "FT0VTx/BC_E/nBCsVsBCID", "FT0VTx/BC_L/nBCsVsBCID", "FT0VTx/BC_LE/nBCsVsBCID", "FT0VTx/BC_NL/nBCsVsBCID", "FT0VTx/BC_NLE/nBCsVsBCID", "FT0VTx/BC_SL_FDD/nBCsVsBCID", "FT0VTx/BC_SL_FT0/nBCsVsBCID", "FT0VTx/BC_NSL_FDD/nBCsVsBCID", "FT0VTx/BC_NSL_FT0/nBCsVsBCID"},
    {"FT0CE/BC_A/nBCsVsBCID", "FT0CE/BC_B/nBCsVsBCID", "FT0CE/BC_C/nBCsVsBCID", "FT0CE/BC_E/nBCsVsBCID", "FT0CE/BC_L/nBCsVsBCID", "FT0CE/BC_LE/nBCsVsBCID", "FT0CE/BC_NL/nBCsVsBCID", "FT0CE/BC_NLE/nBCsVsBCID", "FT0CE/BC_SL_FDD/nBCsVsBCID", "FT0CE/BC_SL_FT0/nBCsVsBCID", "FT0CE/BC_NSL_FDD/nBCsVsBCID", "FT0CE/BC_NSL_FT0/nBCsVsBCID"},
    {"FDD/BC_A/nBCsVsBCID", "FDD/BC_B/nBCsVsBCID", "FDD/BC_C/nBCsVsBCID", "FDD/BC_E/nBCsVsBCID", "FDD/BC_L/nBCsVsBCID", "FDD/BC_LE/nBCsVsBCID", "FDD/BC_NL/nBCsVsBCID", "FDD/BC_NLE/nBCsVsBCID", "FDD/BC_SL_FDD/nBCsVsBCID", "FDD/BC_SL_FT0/nBCsVsBCID", "FDD/BC_NSL_FDD/nBCsVsBCID", "FDD/BC_NSL_FT0/nBCsVsBCID"},
    {"1ZNC/BC_A/nBCsVsBCID", "1ZNC/BC_B/nBCsVsBCID", "1ZNC/BC_C/nBCsVsBCID", "1ZNC/BC_E/nBCsVsBCID", "1ZNC/BC_L/nBCsVsBCID", "1ZNC/BC_LE/nBCsVsBCID", "1ZNC/BC_NL/nBCsVsBCID", "1ZNC/BC_NLE/nBCsVsBCID", "1ZNC/BC_SL_FDD/nBCsVsBCID", "1ZNC/BC_SL_FT0/nBCsVsBCID", "1ZNC/BC_NSL_FDD/nBCsVsBCID", "1ZNC/BC_NSL_FT0/nBCsVsBCID"}};

  static constexpr std::string_view NBCsInspectedVsBCIDHistNames[5][12] = {
    {"AllBCs/BC_A/nBCsInspectedVsBCID", "AllBCs/BC_B/nBCsInspectedVsBCID", "AllBCs/BC_C/nBCsInspectedVsBCID", "AllBCs/BC_E/nBCsInspectedVsBCID", "AllBCs/BC_L/nBCsInspectedVsBCID", "AllBCs/BC_LE/nBCsInspectedVsBCID", "AllBCs/BC_NL/nBCsInspectedVsBCID", "AllBCs/BC_NLE/nBCsInspectedVsBCID", "AllBCs/BC_SL_FDD/nBCsInspectedVsBCID", "AllBCs/BC_SL_FT0/nBCsInspectedVsBCID", "AllBCs/BC_NSL_FDD/nBCsInspectedVsBCID", "AllBCs/BC_NSL_FT0/nBCsInspectedVsBCID"},
    {"FT0VTx/BC_A/nBCsInspectedVsBCID", "FT0VTx/BC_B/nBCsInspectedVsBCID", "FT0VTx/BC_C/nBCsInspectedVsBCID", "FT0VTx/BC_E/nBCsInspectedVsBCID", "FT0VTx/BC_L/nBCsInspectedVsBCID", "FT0VTx/BC_LE/nBCsInspectedVsBCID", "FT0VTx/BC_NL/nBCsInspectedVsBCID", "FT0VTx/BC_NLE/nBCsInspectedVsBCID", "FT0VTx/BC_SL_FDD/nBCsInspectedVsBCID", "FT0VTx/BC_SL_FT0/nBCsInspectedVsBCID", "FT0VTx/BC_NSL_FDD/nBCsInspectedVsBCID", "FT0VTx/BC_NSL_FT0/nBCsInspectedVsBCID"},
    {"FT0CE/BC_A/nBCsInspectedVsBCID", "FT0CE/BC_B/nBCsInspectedVsBCID", "FT0CE/BC_C/nBCsInspectedVsBCID", "FT0CE/BC_E/nBCsInspectedVsBCID", "FT0CE/BC_L/nBCsInspectedVsBCID", "FT0CE/BC_LE/nBCsInspectedVsBCID", "FT0CE/BC_NL/nBCsInspectedVsBCID", "FT0CE/BC_NLE/nBCsInspectedVsBCID", "FT0CE/BC_SL_FDD/nBCsInspectedVsBCID", "FT0CE/BC_SL_FT0/nBCsInspectedVsBCID", "FT0CE/BC_NSL_FDD/nBCsInspectedVsBCID", "FT0CE/BC_NSL_FT0/nBCsInspectedVsBCID"},
    {"FDD/BC_A/nBCsInspectedVsBCID", "FDD/BC_B/nBCsInspectedVsBCID", "FDD/BC_C/nBCsInspectedVsBCID", "FDD/BC_E/nBCsInspectedVsBCID", "FDD/BC_L/nBCsInspectedVsBCID", "FDD/BC_LE/nBCsInspectedVsBCID", "FDD/BC_NL/nBCsInspectedVsBCID", "FDD/BC_NLE/nBCsInspectedVsBCID", "FDD/BC_SL_FDD/nBCsInspectedVsBCID", "FDD/BC_SL_FT0/nBCsInspectedVsBCID", "FDD/BC_NSL_FDD/nBCsInspectedVsBCID", "FDD/BC_NSL_FT0/nBCsInspectedVsBCID"},
    {"1ZNC/BC_A/nBCsInspectedVsBCID", "1ZNC/BC_B/nBCsInspectedVsBCID", "1ZNC/BC_C/nBCsInspectedVsBCID", "1ZNC/BC_E/nBCsInspectedVsBCID", "1ZNC/BC_L/nBCsInspectedVsBCID", "1ZNC/BC_LE/nBCsInspectedVsBCID", "1ZNC/BC_NL/nBCsInspectedVsBCID", "1ZNC/BC_NLE/nBCsInspectedVsBCID", "1ZNC/BC_SL_FDD/nBCsInspectedVsBCID", "1ZNC/BC_SL_FT0/nBCsInspectedVsBCID", "1ZNC/BC_NSL_FDD/nBCsInspectedVsBCID", "1ZNC/BC_NSL_FT0/nBCsInspectedVsBCID"}};

  std::array<std::array<std::shared_ptr<TH1>, nBCCategories>, nTriggers> mInspectedHistos{};

  int64_t bcSOR = 0;
  int nBCsPerTF = 0;
  int64_t currentTFid = -1;

  int64_t globalBCIdOfLastBCWithActivityFDD{std::numeric_limits<int64_t>::min() / 2};
  int64_t globalBCIdOfLastBCWithActivityFT0{std::numeric_limits<int64_t>::min() / 2};
  int64_t globalBCLastInspectedBC{-1};

  using DenomCounter = std::vector<std::array<std::array<int, nBCCategories>, nTriggers>>;

  bool hasAnyFDDTrigger(const std::bitset<64>& ctpInputMask) const
  {
    return ctpInputMask.test(12) || ctpInputMask.test(14) || ctpInputMask.test(15) ||
           ctpInputMask.test(16) || ctpInputMask.test(17);
  }

  bool hasAnyFT0Trigger(const std::bitset<64>& ctpInputMask) const
  {
    return ctpInputMask.test(0) || ctpInputMask.test(1) || ctpInputMask.test(2) ||
           ctpInputMask.test(3) || ctpInputMask.test(4);
  }

  void init(InitContext&)
  {
    strLPMProductionTag = metadataInfo.get("LPMProductionTag");

    LOG(info) << "strLPMProductionTag: " << strLPMProductionTag;

    AxisSpec timeAxis{1440, 0., 1440., "#bf{t-t_{SOF} (min)}"};
    AxisSpec bcIDAxis{3600, 0., 3600., "#bf{BC ID in orbit}"};

    for (int iTrigger = 0; iTrigger < nTriggers; iTrigger++) {
      if ((iTrigger == kAllBCs) || (iTrigger == kFT0Vtx && cfgDoFT0Vtx) || (iTrigger == kFT0CE && cfgDoFT0CE) || (iTrigger == kFDD && cfgDoFDD) || (iTrigger == k1ZNC && cfgDo1ZNC)) {
        for (int iBCCategory = 0; iBCCategory < nBCCategories; iBCCategory++) {
          if ((iBCCategory == kBCA && cfgDoBCA) || (iBCCategory == kBCB && cfgDoBCB) || (iBCCategory == kBCC && cfgDoBCC) || (iBCCategory == kBCE && cfgDoBCE) ||
              (iBCCategory == kBCL && cfgDoBCL) || (iBCCategory == kBCLE && cfgDoBCLE) ||
              (iBCCategory == kBCNL && cfgDoBCNL) || (iBCCategory == kBCNLE && cfgDoBCNLE) ||
              (iBCCategory == kBCSLFDD && cfgDoBCSLFDD) || (iBCCategory == kBCSLFT0 && cfgDoBCSLFT0) ||
              (iBCCategory == kBCNSLFDD && cfgDoBCNSLFDD) || (iBCCategory == kBCNSLFT0 && cfgDoBCNSLFT0)) {
            mHistManager.add(Form("%s", std::string(NBCsVsTimeHistNames[iTrigger][iBCCategory]).c_str()), "Time of triggered BCs since the start of fill;#bf{t-t_{SOF} (min)};#bf{#it{N}_{BC}}", HistType::kTH1D, {timeAxis});
            mHistManager.add(Form("%s", std::string(NBCsVsBCIDHistNames[iTrigger][iBCCategory]).c_str()), "BC ID of triggered BCs;#bf{BC ID in orbit};#bf{#it{N}_{BC}}", HistType::kTH1D, {bcIDAxis});
            mInspectedHistos[iTrigger][iBCCategory] = mHistManager.add<TH1>(
              Form("%s", std::string(NBCsInspectedVsBCIDHistNames[iTrigger][iBCCategory]).c_str()),
              "Inspected BC ID (denominator for mu);#bf{BC ID in orbit};#bf{#it{N}_{BC}}",
              HistType::kTH1D, {bcIDAxis});
          }
        }
      }
    }

    if (cfgDoBCSLFDD || cfgDoBCSLFT0 || cfgDoBCNSLFDD || cfgDoBCNSLFT0) {
      mHistManager.add("FITQA/BCHasFT0", "Does the BC have FT0?;BC has FT0;TVX triggered according to CTP;#bf{#it{N}_{BC}}", HistType::kTH2D, {{2, -0.5, 1.5}, {2, -0.5, 1.5}});
      mHistManager.get<TH2>(HIST("FITQA/BCHasFT0")).get()->GetYaxis()->SetBinLabel(1, "No CTP trigger");
      mHistManager.get<TH2>(HIST("FITQA/BCHasFT0")).get()->GetYaxis()->SetBinLabel(2, "CTP triggered");
      mHistManager.get<TH2>(HIST("FITQA/BCHasFT0")).get()->GetXaxis()->SetBinLabel(1, "No found FT0");
      mHistManager.get<TH2>(HIST("FITQA/BCHasFT0")).get()->GetXaxis()->SetBinLabel(2, "Found FT0");

      mHistManager.add("FITQA/BCHasFDD", "Does the BC have FDD?;BC has FDD;FDD triggered according to CTP;#bf{#it{N}_{BC}}", HistType::kTH2D, {{2, -0.5, 1.5}, {2, -0.5, 1.5}});
      mHistManager.get<TH2>(HIST("FITQA/BCHasFDD")).get()->GetYaxis()->SetBinLabel(1, "No CTP trigger");
      mHistManager.get<TH2>(HIST("FITQA/BCHasFDD")).get()->GetYaxis()->SetBinLabel(2, "CTP triggered");
      mHistManager.get<TH2>(HIST("FITQA/BCHasFDD")).get()->GetXaxis()->SetBinLabel(1, "No found FDD");
      mHistManager.get<TH2>(HIST("FITQA/BCHasFDD")).get()->GetXaxis()->SetBinLabel(2, "Found FDD");
    }

    mHistManager.add("FT0Vtx_EvSel/nBCsVsTime", "Time of TVX triggered BCs since the start of fill;;#bf{#it{N}_{BC}}", HistType::kTH1D, {timeAxis});
    mHistManager.add("nBCsVsBCID", "Time of TVX triggered BCs since the start of fill;#bf{t-t_{SOF} (min)};#bf{#it{N}_{BC}}", HistType::kTH1D, {bcIDAxis});
    mHistManager.add("TFsPerMinute", "TFs seen in this minute (to account for failed jobs);#bf{t-t_{SOF} (min)};#bf{#it{N}_{TFs}}", HistType::kTH1D, {timeAxis});

    if (cfgDo1ZNC) {
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
  }

  void setLHCIFData(const auto& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }

    auto& ccdbMgr = o2::ccdb::BasicCCDBManager::instance();
    uint64_t timeStamp = bc.timestamp();

    const int runStart2023{535069};
    const int runStop2023{543113};
    isData23 = (bc.runNumber() >= runStart2023 && bc.runNumber() <= runStop2023);
    mBcShiftFDD = isData23 ? static_cast<int>(cfgBcShiftFDDForData2023) : 0;

    std::map<std::string, std::string> metadata;
    mLHCIFdata = ccdbMgr.getSpecific<o2::parameters::GRPLHCIFData>("GLO/Config/GRPLHCIF", timeStamp, metadata);
    if (mLHCIFdata == nullptr) {
      LOG(fatal) << "GRPLHCIFData not in database, timestamp:" << timeStamp;
    }

    mRunNumber = bc.runNumber();
    LOG(info) << "LHCIF data fetched for run " << mRunNumber << " and timestamp " << timeStamp
              << " (isData23=" << isData23 << ", bcShiftFDD=" << mBcShiftFDD << ")";

    beamPatternA = mLHCIFdata->getBunchFilling().getBeamPattern(0);
    beamPatternC = mLHCIFdata->getBunchFilling().getBeamPattern(1);
    bcPatternA = beamPatternA & ~beamPatternC;
    bcPatternC = ~beamPatternA & beamPatternC;
    bcPatternB = beamPatternA & beamPatternC;
    bcPatternE = ~beamPatternA & ~beamPatternC;

    bcPatternL.reset();
    bcPatternLE.reset();

    int totalLeadingBCsL = 0;
    int totalLeadingBCsLE = 0;

    for (int iBC = 0; iBC < o2::constants::lhc::LHCMaxBunches; iBC++) {
      if (!bcPatternB[iBC]) {
        continue;
      }

      int nonBBefore = 0;
      int emptyBefore = 0;

      for (int j = 1; j <= cfgEmptyBCsBeforeLeadingBC; j++) {
        int prevBC = (iBC - j + o2::constants::lhc::LHCMaxBunches) % o2::constants::lhc::LHCMaxBunches;
        if (!bcPatternB[prevBC]) {
          nonBBefore++;
        } else {
          break;
        }
      }

      for (int j = 1; j <= cfgEmptyBCsBeforeLeadingBCLE; j++) {
        int prevBC = (iBC - j + o2::constants::lhc::LHCMaxBunches) % o2::constants::lhc::LHCMaxBunches;
        if (bcPatternE[prevBC]) {
          emptyBefore++;
        } else {
          break;
        }
      }

      if (nonBBefore >= cfgEmptyBCsBeforeLeadingBC) {
        bcPatternL[iBC] = true;
        totalLeadingBCsL++;
      }
      if (emptyBefore >= cfgEmptyBCsBeforeLeadingBCLE) {
        bcPatternLE[iBC] = true;
        totalLeadingBCsLE++;
      }
    }

    LOG(info) << "bcPatternL (non-B before) complete. Leading BCs found: " << totalLeadingBCsL;
    LOG(info) << "bcPatternLE (empty before) complete. Leading BCs found: " << totalLeadingBCsLE;

    auto runInfo = o2::parameters::AggregatedRunInfo::buildAggregatedRunInfo(o2::ccdb::BasicCCDBManager::instance(), mRunNumber, strLPMProductionTag);
    bcSOR = runInfo.orbitSOR * nBCsPerOrbit;
    LOG(info) << "BC SOR: " << bcSOR << " (orbit SOR: " << runInfo.orbitSOR << ") NBCs per orbit: " << nBCsPerOrbit;
    nBCsPerTF = runInfo.orbitsPerTF * nBCsPerOrbit;
  }

  float getTimeSinceSOF(const auto& bc)
  {
    return (bc.timestamp() - mLHCIFdata->getFillNumberTime()) / 1e3 / 60.f;
  }

  template <int iTrigger, int iBCCategory>
  void fillHistograms(float timeSinceSOF, int64_t localBC)
  {
    mHistManager.fill(HIST(NBCsVsTimeHistNames[iTrigger][iBCCategory]), timeSinceSOF);
    mHistManager.fill(HIST(NBCsVsBCIDHistNames[iTrigger][iBCCategory]), localBC);
  }

  template <int iTrigger>
  void countInspectedBC(DenomCounter& nBCsPerBcId,
                        int iLBC,
                        int iLBCFDD,
                        int64_t iGBC,
                        int64_t lastFT0ActivityBC,
                        int64_t lastFDDActivityBC)
  {
    if constexpr (iTrigger == kFDD) {
      if (bcPatternA[iLBCFDD]) {
        nBCsPerBcId[iLBCFDD][iTrigger][kBCA]++;
      }
      if (bcPatternB[iLBCFDD]) {
        nBCsPerBcId[iLBCFDD][iTrigger][kBCB]++;
      }
      if (bcPatternC[iLBCFDD]) {
        nBCsPerBcId[iLBCFDD][iTrigger][kBCC]++;
      }
      if (bcPatternE[iLBCFDD]) {
        nBCsPerBcId[iLBCFDD][iTrigger][kBCE]++;
      }
      if (bcPatternL[iLBCFDD]) {
        nBCsPerBcId[iLBCFDD][iTrigger][kBCL]++;
      }
      if (bcPatternLE[iLBCFDD]) {
        nBCsPerBcId[iLBCFDD][iTrigger][kBCLE]++;
      }
      if (bcPatternB[iLBCFDD] && !bcPatternL[iLBCFDD]) {
        nBCsPerBcId[iLBCFDD][iTrigger][kBCNL]++;
      }
      if (bcPatternB[iLBCFDD] && !bcPatternLE[iLBCFDD]) {
        nBCsPerBcId[iLBCFDD][iTrigger][kBCNLE]++;
      }
      if (bcPatternB[iLBCFDD]) {
        const bool slFDD = ((iGBC + mBcShiftFDD) - lastFDDActivityBC >= cfgBCsBeforeSuperLeading);
        nBCsPerBcId[iLBCFDD][iTrigger][slFDD ? kBCSLFDD : kBCNSLFDD]++;

        const bool slFT0 = (iGBC - lastFT0ActivityBC >= cfgBCsBeforeSuperLeading);
        nBCsPerBcId[iLBCFDD][iTrigger][slFT0 ? kBCSLFT0 : kBCNSLFT0]++;
      }
    } else {
      if (bcPatternA[iLBC]) {
        nBCsPerBcId[iLBC][iTrigger][kBCA]++;
      }
      if (bcPatternB[iLBC]) {
        nBCsPerBcId[iLBC][iTrigger][kBCB]++;
      }
      if (bcPatternC[iLBC]) {
        nBCsPerBcId[iLBC][iTrigger][kBCC]++;
      }
      if (bcPatternE[iLBC]) {
        nBCsPerBcId[iLBC][iTrigger][kBCE]++;
      }
      if (bcPatternL[iLBC]) {
        nBCsPerBcId[iLBC][iTrigger][kBCL]++;
      }
      if (bcPatternLE[iLBC]) {
        nBCsPerBcId[iLBC][iTrigger][kBCLE]++;
      }
      if (bcPatternB[iLBC] && !bcPatternL[iLBC]) {
        nBCsPerBcId[iLBC][iTrigger][kBCNL]++;
      }
      if (bcPatternB[iLBC] && !bcPatternLE[iLBC]) {
        nBCsPerBcId[iLBC][iTrigger][kBCNLE]++;
      }
      if (bcPatternB[iLBC]) {
        const bool slFT0 = (iGBC - lastFT0ActivityBC >= cfgBCsBeforeSuperLeading);
        nBCsPerBcId[iLBC][iTrigger][slFT0 ? kBCSLFT0 : kBCNSLFT0]++;
      }
      if (bcPatternB[iLBCFDD]) {
        const bool slFDD = ((iGBC + mBcShiftFDD) - lastFDDActivityBC >= cfgBCsBeforeSuperLeading);
        nBCsPerBcId[iLBCFDD][iTrigger][slFDD ? kBCSLFDD : kBCNSLFDD]++;
      }
    }
  }

  void processZDCQA(MyBCs const& bcs, aod::Zdcs const&)
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
  PROCESS_SWITCH(LumiStabilityLightIons, processZDCQA, "process QA for the ZDC triggers (light ions and PbPb)", false);

  void process(MyBCs const& bcs, aod::FT0s const&, aod::FDDs const&)
  {
    DenomCounter nBCsPerBcId(nBCsPerOrbit);
    for (auto& triggerArr : nBCsPerBcId) {
      for (auto& catArr : triggerArr) {
        catArr.fill(0);
      }
    }

    for (const auto& bc : bcs) {
      if (bc.timestamp() == 0) {
        continue;
      }

      setLHCIFData(bc);

      float timeSinceSOF = getTimeSinceSOF(bc);

      if (bc.selection_bit(aod::evsel::kIsTriggerTVX)) {
        mHistManager.fill(HIST("FT0Vtx_EvSel/nBCsVsTime"), timeSinceSOF);
      }

      int64_t globalBC = bc.globalBC();
      int localBC = static_cast<int>(globalBC % nBCsPerOrbit);

      int64_t globalBCFDD = globalBC + mBcShiftFDD;
      int localBCFDD = static_cast<int>((globalBCFDD % nBCsPerOrbit + nBCsPerOrbit) % nBCsPerOrbit);

      int64_t thisTFid = (globalBC - bcSOR) / nBCsPerTF;
      if (thisTFid != currentTFid) {
        currentTFid = thisTFid;
        mHistManager.fill(HIST("TFsPerMinute"), timeSinceSOF);
      }

      std::bitset<64> ctpInputMask(bc.inputMask());

      const bool anyFT0Trigger = hasAnyFT0Trigger(ctpInputMask);
      const bool anyFDDTrigger = hasAnyFDDTrigger(ctpInputMask);

      bool isSuperLeadingBcFDD = bcPatternB[localBCFDD] &&
                                 (globalBCFDD - globalBCIdOfLastBCWithActivityFDD >= cfgBCsBeforeSuperLeading);

      bool isSuperLeadingBcFT0 = bcPatternB[localBC] &&
                                 (globalBC - globalBCIdOfLastBCWithActivityFT0 >= cfgBCsBeforeSuperLeading);

      if (cfgDoBCSLFDD || cfgDoBCSLFT0 || cfgDoBCNSLFDD || cfgDoBCNSLFT0) {
        mHistManager.fill(HIST("FITQA/BCHasFT0"), bc.has_ft0(), ctpInputMask.test(2));
        mHistManager.fill(HIST("FITQA/BCHasFDD"), bc.has_fdd(), anyFDDTrigger);
      }

      int64_t globalBCStart = (globalBCLastInspectedBC >= 0 && globalBCLastInspectedBC < globalBC) ? globalBCLastInspectedBC + 1 : globalBC;
      const int64_t maxBcGap = 2LL * nBCsPerOrbit;
      if (globalBC - globalBCStart > maxBcGap) {
        globalBCStart = globalBC;
      }

      for (int64_t iGBC = globalBCStart; iGBC <= globalBC; ++iGBC) {
        const int iLBC = static_cast<int>((iGBC % nBCsPerOrbit + nBCsPerOrbit) % nBCsPerOrbit);
        const int iLBCFDD = static_cast<int>(((iGBC + mBcShiftFDD) % nBCsPerOrbit + nBCsPerOrbit) % nBCsPerOrbit);

        countInspectedBC<kAllBCs>(nBCsPerBcId, iLBC, iLBCFDD, iGBC,
                                  globalBCIdOfLastBCWithActivityFT0,
                                  globalBCIdOfLastBCWithActivityFDD);

        if (cfgDoFT0Vtx && ctpInputMask.test(2)) {
          countInspectedBC<kFT0Vtx>(nBCsPerBcId, iLBC, iLBCFDD, iGBC,
                                    globalBCIdOfLastBCWithActivityFT0,
                                    globalBCIdOfLastBCWithActivityFDD);
        }

        if (cfgDoFT0CE && ctpInputMask.test(4)) {
          countInspectedBC<kFT0CE>(nBCsPerBcId, iLBC, iLBCFDD, iGBC,
                                   globalBCIdOfLastBCWithActivityFT0,
                                   globalBCIdOfLastBCWithActivityFDD);
        }

        if (cfgDoFDD && anyFDDTrigger) {
          countInspectedBC<kFDD>(nBCsPerBcId, iLBC, iLBCFDD, iGBC,
                                 globalBCIdOfLastBCWithActivityFT0,
                                 globalBCIdOfLastBCWithActivityFDD);
        }

        if (cfgDo1ZNC && ctpInputMask.test(25)) {
          countInspectedBC<k1ZNC>(nBCsPerBcId, iLBC, iLBCFDD, iGBC,
                                  globalBCIdOfLastBCWithActivityFT0,
                                  globalBCIdOfLastBCWithActivityFDD);
        }
      }

      if (anyFDDTrigger) {
        globalBCIdOfLastBCWithActivityFDD = globalBCFDD;
      }
      if (anyFT0Trigger) {
        globalBCIdOfLastBCWithActivityFT0 = globalBC;
      }

      globalBCLastInspectedBC = globalBC;

      if (cfgDoBCA && bcPatternA[localBC]) {
        fillHistograms<kAllBCs, kBCA>(timeSinceSOF, localBC);
      }
      if (cfgDoBCB && bcPatternB[localBC]) {
        fillHistograms<kAllBCs, kBCB>(timeSinceSOF, localBC);
      }
      if (cfgDoBCC && bcPatternC[localBC]) {
        fillHistograms<kAllBCs, kBCC>(timeSinceSOF, localBC);
      }
      if (cfgDoBCE && bcPatternE[localBC]) {
        fillHistograms<kAllBCs, kBCE>(timeSinceSOF, localBC);
      }
      if (cfgDoBCL && bcPatternL[localBC]) {
        fillHistograms<kAllBCs, kBCL>(timeSinceSOF, localBC);
      }
      if (cfgDoBCLE && bcPatternLE[localBC]) {
        fillHistograms<kAllBCs, kBCLE>(timeSinceSOF, localBC);
      }
      if (cfgDoBCNL && bcPatternB[localBC] && !bcPatternL[localBC]) {
        fillHistograms<kAllBCs, kBCNL>(timeSinceSOF, localBC);
      }
      if (cfgDoBCNLE && bcPatternB[localBC] && !bcPatternLE[localBC]) {
        fillHistograms<kAllBCs, kBCNLE>(timeSinceSOF, localBC);
      }
      if (cfgDoBCSLFDD && isSuperLeadingBcFDD) {
        fillHistograms<kAllBCs, kBCSLFDD>(timeSinceSOF, localBCFDD);
      }
      if (cfgDoBCSLFT0 && isSuperLeadingBcFT0) {
        fillHistograms<kAllBCs, kBCSLFT0>(timeSinceSOF, localBC);
      }
      if (cfgDoBCNSLFDD && bcPatternB[localBCFDD] && !isSuperLeadingBcFDD) {
        fillHistograms<kAllBCs, kBCNSLFDD>(timeSinceSOF, localBCFDD);
      }
      if (cfgDoBCNSLFT0 && bcPatternB[localBC] && !isSuperLeadingBcFT0) {
        fillHistograms<kAllBCs, kBCNSLFT0>(timeSinceSOF, localBC);
      }

      if (cfgDoFT0Vtx && ctpInputMask.test(2)) {
        if (cfgDoBCA && bcPatternA[localBC])
          fillHistograms<kFT0Vtx, kBCA>(timeSinceSOF, localBC);
        if (cfgDoBCB && bcPatternB[localBC])
          fillHistograms<kFT0Vtx, kBCB>(timeSinceSOF, localBC);
        if (cfgDoBCC && bcPatternC[localBC])
          fillHistograms<kFT0Vtx, kBCC>(timeSinceSOF, localBC);
        if (cfgDoBCE && bcPatternE[localBC])
          fillHistograms<kFT0Vtx, kBCE>(timeSinceSOF, localBC);
        if (cfgDoBCL && bcPatternL[localBC])
          fillHistograms<kFT0Vtx, kBCL>(timeSinceSOF, localBC);
        if (cfgDoBCLE && bcPatternLE[localBC])
          fillHistograms<kFT0Vtx, kBCLE>(timeSinceSOF, localBC);
        if (cfgDoBCNL && bcPatternB[localBC] && !bcPatternL[localBC])
          fillHistograms<kFT0Vtx, kBCNL>(timeSinceSOF, localBC);
        if (cfgDoBCNLE && bcPatternB[localBC] && !bcPatternLE[localBC])
          fillHistograms<kFT0Vtx, kBCNLE>(timeSinceSOF, localBC);
        if (cfgDoBCSLFDD && isSuperLeadingBcFDD)
          fillHistograms<kFT0Vtx, kBCSLFDD>(timeSinceSOF, localBCFDD);
        if (cfgDoBCSLFT0 && isSuperLeadingBcFT0)
          fillHistograms<kFT0Vtx, kBCSLFT0>(timeSinceSOF, localBC);
        if (cfgDoBCNSLFDD && bcPatternB[localBCFDD] && !isSuperLeadingBcFDD)
          fillHistograms<kFT0Vtx, kBCNSLFDD>(timeSinceSOF, localBCFDD);
        if (cfgDoBCNSLFT0 && bcPatternB[localBC] && !isSuperLeadingBcFT0)
          fillHistograms<kFT0Vtx, kBCNSLFT0>(timeSinceSOF, localBC);
      }

      if (cfgDoFT0CE && ctpInputMask.test(4)) {
        if (cfgDoBCA && bcPatternA[localBC])
          fillHistograms<kFT0CE, kBCA>(timeSinceSOF, localBC);
        if (cfgDoBCB && bcPatternB[localBC])
          fillHistograms<kFT0CE, kBCB>(timeSinceSOF, localBC);
        if (cfgDoBCC && bcPatternC[localBC])
          fillHistograms<kFT0CE, kBCC>(timeSinceSOF, localBC);
        if (cfgDoBCE && bcPatternE[localBC])
          fillHistograms<kFT0CE, kBCE>(timeSinceSOF, localBC);
        if (cfgDoBCL && bcPatternL[localBC])
          fillHistograms<kFT0CE, kBCL>(timeSinceSOF, localBC);
        if (cfgDoBCLE && bcPatternLE[localBC])
          fillHistograms<kFT0CE, kBCLE>(timeSinceSOF, localBC);
        if (cfgDoBCNL && bcPatternB[localBC] && !bcPatternL[localBC])
          fillHistograms<kFT0CE, kBCNL>(timeSinceSOF, localBC);
        if (cfgDoBCNLE && bcPatternB[localBC] && !bcPatternLE[localBC])
          fillHistograms<kFT0CE, kBCNLE>(timeSinceSOF, localBC);
        if (cfgDoBCSLFDD && isSuperLeadingBcFDD)
          fillHistograms<kFT0CE, kBCSLFDD>(timeSinceSOF, localBCFDD);
        if (cfgDoBCSLFT0 && isSuperLeadingBcFT0)
          fillHistograms<kFT0CE, kBCSLFT0>(timeSinceSOF, localBC);
        if (cfgDoBCNSLFDD && bcPatternB[localBCFDD] && !isSuperLeadingBcFDD)
          fillHistograms<kFT0CE, kBCNSLFDD>(timeSinceSOF, localBCFDD);
        if (cfgDoBCNSLFT0 && bcPatternB[localBC] && !isSuperLeadingBcFT0)
          fillHistograms<kFT0CE, kBCNSLFT0>(timeSinceSOF, localBC);
      }

      if (cfgDoFDD && anyFDDTrigger) {
        if (cfgDoBCA && bcPatternA[localBCFDD])
          fillHistograms<kFDD, kBCA>(timeSinceSOF, localBCFDD);
        if (cfgDoBCB && bcPatternB[localBCFDD])
          fillHistograms<kFDD, kBCB>(timeSinceSOF, localBCFDD);
        if (cfgDoBCC && bcPatternC[localBCFDD])
          fillHistograms<kFDD, kBCC>(timeSinceSOF, localBCFDD);
        if (cfgDoBCE && bcPatternE[localBCFDD])
          fillHistograms<kFDD, kBCE>(timeSinceSOF, localBCFDD);
        if (cfgDoBCL && bcPatternL[localBCFDD])
          fillHistograms<kFDD, kBCL>(timeSinceSOF, localBCFDD);
        if (cfgDoBCLE && bcPatternLE[localBCFDD])
          fillHistograms<kFDD, kBCLE>(timeSinceSOF, localBCFDD);
        if (cfgDoBCNL && bcPatternB[localBCFDD] && !bcPatternL[localBCFDD])
          fillHistograms<kFDD, kBCNL>(timeSinceSOF, localBCFDD);
        if (cfgDoBCNLE && bcPatternB[localBCFDD] && !bcPatternLE[localBCFDD])
          fillHistograms<kFDD, kBCNLE>(timeSinceSOF, localBCFDD);
        if (cfgDoBCSLFDD && isSuperLeadingBcFDD)
          fillHistograms<kFDD, kBCSLFDD>(timeSinceSOF, localBCFDD);
        if (cfgDoBCSLFT0 && bcPatternB[localBCFDD] && isSuperLeadingBcFT0)
          fillHistograms<kFDD, kBCSLFT0>(timeSinceSOF, localBCFDD);
        if (cfgDoBCNSLFDD && bcPatternB[localBCFDD] && !isSuperLeadingBcFDD)
          fillHistograms<kFDD, kBCNSLFDD>(timeSinceSOF, localBCFDD);
        if (cfgDoBCNSLFT0 && bcPatternB[localBCFDD] && !isSuperLeadingBcFT0)
          fillHistograms<kFDD, kBCNSLFT0>(timeSinceSOF, localBCFDD);
      }

      if (cfgDo1ZNC && ctpInputMask.test(25)) {
        if (cfgDoBCA && bcPatternA[localBC])
          fillHistograms<k1ZNC, kBCA>(timeSinceSOF, localBC);
        if (cfgDoBCB && bcPatternB[localBC])
          fillHistograms<k1ZNC, kBCB>(timeSinceSOF, localBC);
        if (cfgDoBCC && bcPatternC[localBC])
          fillHistograms<k1ZNC, kBCC>(timeSinceSOF, localBC);
        if (cfgDoBCE && bcPatternE[localBC])
          fillHistograms<k1ZNC, kBCE>(timeSinceSOF, localBC);
        if (cfgDoBCL && bcPatternL[localBC])
          fillHistograms<k1ZNC, kBCL>(timeSinceSOF, localBC);
        if (cfgDoBCLE && bcPatternLE[localBC])
          fillHistograms<k1ZNC, kBCLE>(timeSinceSOF, localBC);
        if (cfgDoBCNL && bcPatternB[localBC] && !bcPatternL[localBC])
          fillHistograms<k1ZNC, kBCNL>(timeSinceSOF, localBC);
        if (cfgDoBCNLE && bcPatternB[localBC] && !bcPatternLE[localBC])
          fillHistograms<k1ZNC, kBCNLE>(timeSinceSOF, localBC);
        if (cfgDoBCSLFDD && isSuperLeadingBcFDD)
          fillHistograms<k1ZNC, kBCSLFDD>(timeSinceSOF, localBCFDD);
        if (cfgDoBCSLFT0 && isSuperLeadingBcFT0)
          fillHistograms<k1ZNC, kBCSLFT0>(timeSinceSOF, localBC);
        if (cfgDoBCNSLFDD && bcPatternB[localBCFDD] && !isSuperLeadingBcFDD)
          fillHistograms<k1ZNC, kBCNSLFDD>(timeSinceSOF, localBCFDD);
        if (cfgDoBCNSLFT0 && bcPatternB[localBC] && !isSuperLeadingBcFT0)
          fillHistograms<k1ZNC, kBCNSLFT0>(timeSinceSOF, localBC);
      }

      mHistManager.fill(HIST("nBCsVsBCID"), localBC);
    }

    for (int iT = 0; iT < nTriggers; ++iT) {
      for (int iC = 0; iC < nBCCategories; ++iC) {
        if (!mInspectedHistos[iT][iC]) {
          continue;
        }
        for (int iBcId = 0; iBcId < nBCsPerOrbit; ++iBcId) {
          const int value = nBCsPerBcId[iBcId][iT][iC];
          if (value > 0) {
            mInspectedHistos[iT][iC]->Fill(iBcId, value);
          }
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  metadataInfo.initMetadata(cfgc);
  return WorkflowSpec{adaptAnalysisTask<LumiStabilityLightIons>(cfgc)};
}
