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
/// \author Nicolas Strangmann (nicolas.strangmann@cern.ch) - Goethe University Frankfurt, Stefanie Mrozinski (stefanie.mrozinski@cern.ch) - Goethe University Frankfurt

#include "Common/CCDB/ctpRateFetcher.h"
#include "Common/Core/MetadataHelper.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/AggregatedRunInfo.h"
#include "DataFormatsParameters/GRPLHCIFData.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"

#include <limits>
#include <map>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

o2::common::core::MetadataHelper metadataInfo; // Metadata helper

namespace o2::aod {
namespace myBc_aod {
DECLARE_SOA_COLUMN(Timestamp, timestamp, uint64_t);
DECLARE_SOA_COLUMN(TimeZNA, timeZNA, float);
DECLARE_SOA_COLUMN(TimeZNC, timeZNC, float);
} //namespace myBc_aod
DECLARE_SOA_TABLE(MyBCaod, "AOD", "MYBCAOD", myBc_aod::Timestamp, myBc_aod::TimeZNA, myBc_aod::TimeZNC);
} //namespace o2::aod

using MyBCs = soa::Join<aod::BCs, aod::BcSels, aod::Timestamps, aod::Run3MatchedToBCSparse>;

struct LumiStabilityLightIons {
  Produces<aod::MyBCaod> BCaod;

  Configurable<bool> cfgDoFT0Vtx{"cfgDoFT0Vtx", true, "Create and fill histograms for the FT0 vertex trigger"};
  Configurable<bool> cfgDoFT0CE{"cfgDoFT0CE", true, "Create and fill histograms for the FT0 centrality trigger"};
  Configurable<bool> cfgDoFDD{"cfgDoFDD", true, "Create and fill histograms for the FDD trigger"};
  Configurable<bool> cfgDo1ZNC{"cfgDo1ZNC", true, "Create and fill histograms for the 1ZNC trigger"};

  Configurable<bool> cfgDoBCA{"cfgDoBCA", false, "Create and fill histograms for the BCs of type A"};
  Configurable<bool> cfgDoBCB{"cfgDoBCB", true, "Create and fill histograms for the BCs of type B"};
  Configurable<bool> cfgDoBCC{"cfgDoBCC", false, "Create and fill histograms for the BCs of type C"};
  Configurable<bool> cfgDoBCE{"cfgDoBCE", false, "Create and fill histograms for the BCs of type E"};
  Configurable<bool> cfgDoBCL{"cfgDoBCL", false, "Create and fill histograms for leading BCs of type B"};
  Configurable<bool> cfgDoBCSL{"cfgDoBCSL", false, "Create and fill histograms for super-leading BCs (no preceding FT0/FDD activity) of type B"};

  Configurable<bool> cfgRequireNoT0ForSLBC{"cfgRequireNoT0ForSLBC", false, "Require no T0 signal for definition of super leading BC (otherwise only no FDD)"};

  Configurable<int> cfgEmptyBCsBeforeLeadingBC{"cfgEmptyBCsBeforeLeadingBC", 5, "Minimum number of empty BCs before a leading BC to identify it as such"};

  //Configurables specific to VdM analysis: output ao2d with timestamps and ZDC times
  Configurable<bool> cfgFillBCao2d{"cfgFillBCao2d", false, "Fill BC ao2d with timestamps and ZDC times"};
  Configurable<uint64_t> cfgTstampStartFillingBCao2d{"cfgTstampStartFillingBCao2d", 0, "Minimum value of timestamp for output bc ao2d to be filled"};
  Configurable<uint64_t> cfgTstampEndFillingBCao2d{"cfgTstampEndFillingBCao2d", 0, "Maximum value of timestamp for output bc ao2d to be filled"};

  std::bitset<o2::constants::lhc::LHCMaxBunches> beamPatternA, beamPatternC;
  std::bitset<o2::constants::lhc::LHCMaxBunches> bcPatternA, bcPatternC, bcPatternB, bcPatternE, bcPatternL;

  std::string strLPMProductionTag = ""; // MC production tag to be retrieved from AO2D metadata

  const int nBCsPerOrbit = 3564;

  parameters::GRPLHCIFData* mLHCIFdata = nullptr;
  int mRunNumber = -1;
  ctpRateFetcher mRateFetcher;
  bool isLeadingBC = false;

  HistogramRegistry mHistManager{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  const int nTriggers = 5;
  enum TriggerAliases { kAllBCs = 0,
                        kFT0Vtx = 1,
                        kFT0CE = 2,
                        kFDD = 3,
                        k1ZNC = 4 };

  const int nBCCategories = 6;
  enum BCCategories { kBCA = 0,    // A side BCs (bunch-crossings that had beam only from A side)
                      kBCB = 1,    // B type BCs (bunch-crossings that had beam from both sides)
                      kBCC = 2,    // C side BCs (bunch-crossings that had beam only from C side)
                      kBCE = 3,    // empty BCs (bunch-crossings that did not have beam from either side)
                      kBCL = 4,    // leading BCs (bunch-crossings that did not have interacting bunches for a configurable number of preceding BCs)
                      kBCSL = 5 }; // super-leading BCs (bunch-crossings that did not have FDD/FT0 activity for a configurable number of preceding BCs)

  static constexpr std::string_view NBCsVsTimeHistNames[5][6] =
    {{"AllBCs/BC_A/nBCsVsTime", "AllBCs/BC_B/nBCsVsTime", "AllBCs/BC_C/nBCsVsTime", "AllBCs/BC_E/nBCsVsTime", "AllBCs/BC_L/nBCsVsTime", "AllBCs/BC_SL/nBCsVsTime"},
     {"FT0VTx/BC_A/nBCsVsTime", "FT0VTx/BC_B/nBCsVsTime", "FT0VTx/BC_C/nBCsVsTime", "FT0VTx/BC_E/nBCsVsTime", "FT0VTx/BC_L/nBCsVsTime", "FT0VTx/BC_SL/nBCsVsTime"},
     {"FT0CE/BC_A/nBCsVsTime", "FT0CE/BC_B/nBCsVsTime", "FT0CE/BC_C/nBCsVsTime", "FT0CE/BC_E/nBCsVsTime", "FT0CE/BC_L/nBCsVsTime", "FT0CE/BC_SL/nBCsVsTime"},
     {"FDD/BC_A/nBCsVsTime", "FDD/BC_B/nBCsVsTime", "FDD/BC_C/nBCsVsTime", "FDD/BC_E/nBCsVsTime", "FDD/BC_L/nBCsVsTime", "FDD/BC_SL/nBCsVsTime"},
     {"1ZNC/BC_A/nBCsVsTime", "1ZNC/BC_B/nBCsVsTime", "1ZNC/BC_C/nBCsVsTime", "1ZNC/BC_E/nBCsVsTime", "1ZNC/BC_L/nBCsVsTime", "1ZNC/BC_SL/nBCsVsTime"}};

  static constexpr std::string_view NBCsVsBCIDHistNames[5][6] =
    {{"AllBCs/BC_A/nBCsVsBCID", "AllBCs/BC_B/nBCsVsBCID", "AllBCs/BC_C/nBCsVsBCID", "AllBCs/BC_E/nBCsVsBCID", "AllBCs/BC_L/nBCsVsBCID", "AllBCs/BC_SL/nBCsVsBCID"},
     {"FT0VTx/BC_A/nBCsVsBCID", "FT0VTx/BC_B/nBCsVsBCID", "FT0VTx/BC_C/nBCsVsBCID", "FT0VTx/BC_E/nBCsVsBCID", "FT0VTx/BC_L/nBCsVsBCID", "FT0VTx/BC_SL/nBCsVsBCID"},
     {"FT0CE/BC_A/nBCsVsBCID", "FT0CE/BC_B/nBCsVsBCID", "FT0CE/BC_C/nBCsVsBCID", "FT0CE/BC_E/nBCsVsBCID", "FT0CE/BC_L/nBCsVsBCID", "FT0CE/BC_SL/nBCsVsBCID"},
     {"FDD/BC_A/nBCsVsBCID", "FDD/BC_B/nBCsVsBCID", "FDD/BC_C/nBCsVsBCID", "FDD/BC_E/nBCsVsBCID", "FDD/BC_L/nBCsVsBCID", "FDD/BC_SL/nBCsVsBCID"},
     {"1ZNC/BC_A/nBCsVsBCID", "1ZNC/BC_B/nBCsVsBCID", "1ZNC/BC_C/nBCsVsBCID", "1ZNC/BC_E/nBCsVsBCID", "1ZNC/BC_L/nBCsVsBCID", "1ZNC/BC_SL/nBCsVsBCID"}};

  int64_t bcSOR;
  int nBCsPerTF;
  int64_t currentTFid = -1;

  void init(InitContext&)
  {
    strLPMProductionTag = metadataInfo.get("LPMProductionTag"); // to extract info from ccdb by the tag

    LOG(info) << "strLPMProductionTag: " << strLPMProductionTag;

    AxisSpec timeAxis{1440, 0., 1440., "#bf{t-t_{SOF} (min)}"}, bcIDAxis{3600, 0., 3600., "#bf{BC ID in orbit}"};

    for (int iTrigger = 0; iTrigger < nTriggers; iTrigger++) {
      if ((iTrigger == kAllBCs) || (iTrigger == kFT0Vtx && cfgDoFT0Vtx) || (iTrigger == kFT0CE && cfgDoFT0CE) || (iTrigger == kFDD && cfgDoFDD) || (iTrigger == k1ZNC && cfgDo1ZNC)) {
        for (int iBCCategory = 0; iBCCategory < nBCCategories; iBCCategory++) {
          if ((iBCCategory == kBCA && cfgDoBCA) || (iBCCategory == kBCB && cfgDoBCB) || (iBCCategory == kBCC && cfgDoBCC) || (iBCCategory == kBCE && cfgDoBCE) || (iBCCategory == kBCL && cfgDoBCL)) {
            mHistManager.add(Form("%s", std::string(NBCsVsTimeHistNames[iTrigger][iBCCategory]).c_str()), "Time of triggered BCs since the start of fill;#bf{t-t_{SOF} (min)};#bf{#it{N}_{BC}}", HistType::kTH1D, {timeAxis});
            mHistManager.add(Form("%s", std::string(NBCsVsBCIDHistNames[iTrigger][iBCCategory]).c_str()), "BC ID of triggered BCs;#bf{BC ID in orbit};#bf{#it{N}_{BC}}", HistType::kTH1D, {bcIDAxis});
          }
        }
        if (cfgDoBCSL && (iTrigger == kFT0Vtx || iTrigger == kFDD)) { // only for FT0Vtx and FDD fill super-leading BC histograms
          mHistManager.add(Form("%s", std::string(NBCsVsTimeHistNames[iTrigger][5]).c_str()), "Time of triggered BCs since the start of fill;#bf{t-t_{SOF} (min)};#bf{#it{N}_{BC}}", HistType::kTH1D, {timeAxis});
          mHistManager.add(Form("%s", std::string(NBCsVsBCIDHistNames[iTrigger][5]).c_str()), "BC ID of triggered BCs;#bf{BC ID in orbit};#bf{#it{N}_{BC}}", HistType::kTH1D, {bcIDAxis});
        }
      }
    }

    if (cfgDoBCSL) {
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
    if (mRunNumber == bc.runNumber())
      return;

    auto& ccdbMgr = o2::ccdb::BasicCCDBManager::instance();
    uint64_t timeStamp = bc.timestamp();

    std::map<std::string, std::string> metadata;
    mLHCIFdata = ccdbMgr.getSpecific<o2::parameters::GRPLHCIFData>("GLO/Config/GRPLHCIF", timeStamp, metadata);
    if (mLHCIFdata == nullptr)
      LOG(fatal) << "GRPLHCIFData not in database, timestamp:" << timeStamp;
    mRunNumber = bc.runNumber();
    LOG(info) << "LHCIF data fetched for run " << mRunNumber << " and timestamp " << timeStamp;

    beamPatternA = mLHCIFdata->getBunchFilling().getBeamPattern(0);
    beamPatternC = mLHCIFdata->getBunchFilling().getBeamPattern(1);
    bcPatternA = beamPatternA & ~beamPatternC;
    bcPatternC = ~beamPatternA & beamPatternC;
    bcPatternB = beamPatternA & beamPatternC;
    bcPatternE = ~beamPatternA & ~beamPatternC;

    // Create bcPatternL: leading BCs of type B that follow at least "cfgEmptyBCsBeforeLeadingBC" empty BCs
    bcPatternL.reset(); // Initialize all bits to false
    LOG(info) << "Starting to create bcPatternL from bcPatternB";
    LOG(info) << "Total number of BCs to check: " << o2::constants::lhc::LHCMaxBunches;

    int totalLeadingBCs = 0;
    for (int iBC = 0; iBC < o2::constants::lhc::LHCMaxBunches; iBC++) {
      if (bcPatternB[iBC]) {    // Check if current BC is of type B
        int emptyBCsBefore = 0; // Count how many consecutive BCs before this one are NOT type B
        for (int j = 1; j <= cfgEmptyBCsBeforeLeadingBC; j++) {
          int prevBC = (iBC - j + o2::constants::lhc::LHCMaxBunches) % o2::constants::lhc::LHCMaxBunches; // Protection for BCs at small indices to check the end of the orbit
          if (!bcPatternB[prevBC]) {
            emptyBCsBefore++;
          } else {
            break; // Stop counting if we hit a type B BC
          }
        }
        if (emptyBCsBefore >= cfgEmptyBCsBeforeLeadingBC) { // If we found at least cfgEmptyBCsBeforeLeadingBC empty BCs before this one, mark it as leading
          bcPatternL[iBC] = true;
          totalLeadingBCs++;
        }
      }
    }
    LOG(info) << "bcPatternL creation complete. Total leading BCs found: " << totalLeadingBCs;

    auto runInfo = o2::parameters::AggregatedRunInfo::buildAggregatedRunInfo(o2::ccdb::BasicCCDBManager::instance(), mRunNumber, strLPMProductionTag);
    bcSOR = runInfo.orbitSOR * nBCsPerOrbit; // first bc of the first orbit
    LOG(info) << "BC SOR: " << bcSOR << " (orbit SOR: " << runInfo.orbitSOR << ") NBCs per orbit: " << nBCsPerOrbit;
    nBCsPerTF = runInfo.orbitsPerTF * nBCsPerOrbit; // duration of TF in bcs

    return;
  }

  float getTimeSinceSOF(const auto& bc)
  {
    return (bc.timestamp() - mLHCIFdata->getFillNumberTime()) / 1e3 / 60; // Convert to minutes
  }

  template <int iTrigger, int iBCCategory>
  void fillHistograms(float timeSinceSOF, int64_t localBC)
  {
    mHistManager.fill(HIST(NBCsVsTimeHistNames[iTrigger][iBCCategory]), timeSinceSOF);
    mHistManager.fill(HIST(NBCsVsBCIDHistNames[iTrigger][iBCCategory]), localBC);
  }

  void processZDCQA(MyBCs const& bcs, aod::Zdcs const&)
  {
    const int maxTimeZDC = 50;     // Maximum time the histogram allows before setting a dummy value
    const int dummyZDCTime = 42.f; // Time value to indicate missing ZDC time
    for (const auto& bc : bcs) {

      std::bitset<64> ctpInputMask(bc.inputMask());

      bool zdcHit = !bc.has_zdc() ? 0 : ((bc.zdc().energyCommonZNC() > -1 && std::abs(bc.zdc().timeZNC()) < 1E5) ? 1 : 0);
      mHistManager.fill(HIST("ZDCQA/BCHasZDC"), zdcHit, ctpInputMask.test(25) ? 1 : 0);
      if (!bc.has_zdc())
        continue;

      mHistManager.fill(HIST("ZDCQA/ZNCTimeVsEnergy"), bc.zdc().energyCommonZNC() > -1 ? bc.zdc().energyCommonZNC() : -1, std::abs(bc.zdc().timeZNC()) < maxTimeZDC ? bc.zdc().timeZNC() : dummyZDCTime);

      float timeZNA = bc.zdc().timeZNA();
      float timeZNC = bc.zdc().timeZNC();

      if (std::abs(timeZNA) > maxTimeZDC) {
        timeZNA = dummyZDCTime; // set dummy value for missing ZDC times to be able to plot them
        mHistManager.fill(HIST("ZDCQA/ZNCTime"), timeZNC);
      }
      if (std::abs(timeZNC) > maxTimeZDC) {
        timeZNC = dummyZDCTime;      // set dummy value for missing ZDC times to be able to plot them
        if (timeZNA != dummyZDCTime) // If ZNA and ZNC are both missing, do not fill the ZNA histogram with the dummy value
          mHistManager.fill(HIST("ZDCQA/ZNATime"), timeZNA);
      }

      mHistManager.fill(HIST("ZDCQA/ZDCTimes"), timeZNA, timeZNC);

      //For VdM analysis: fill timestamps and ZDC times in output tree, if enabled
      uint64_t timestamp = bc.timestamp();
      if(cfgFillBCao2d && timestamp>=cfgTstampStartFillingBCao2d && timestamp<=cfgTstampEndFillingBCao2d) {
        BCaod(timestamp, timeZNA, timeZNC);
      }
    }
  }
  PROCESS_SWITCH(LumiStabilityLightIons, processZDCQA, "process QA for the ZDC triggers (light ions and PbPb)", false);

  void processSLBunches(MyBCs const& bcs, aod::FT0s const&, aod::FDDs const&)
  {
    int64_t globalBCIdOfLastBCWithActivity = 0;
    for (const auto& bc : bcs) {
      if (bc.timestamp() == 0)
        continue;

      setLHCIFData(bc);

      std::bitset<64> ctpInputMask(bc.inputMask());

      mHistManager.fill(HIST("FITQA/BCHasFT0"), bc.has_ft0(), ctpInputMask.test(2));
      mHistManager.fill(HIST("FITQA/BCHasFDD"), bc.has_fdd(), ctpInputMask.test(15));

      int64_t globalBC = bc.globalBC();

      if (globalBC - globalBCIdOfLastBCWithActivity < cfgEmptyBCsBeforeLeadingBC)
        continue; // not a super-leading BC

      if (bc.has_fdd() || (cfgRequireNoT0ForSLBC && bc.has_ft0()))
        globalBCIdOfLastBCWithActivity = globalBC;

      float timeSinceSOF = getTimeSinceSOF(bc);

      int localBC = globalBC % nBCsPerOrbit;

      if (!bcPatternB[localBC])
        continue;

      if (ctpInputMask.test(2))
        fillHistograms<kFT0Vtx, kBCSL>(timeSinceSOF, localBC);
      if (ctpInputMask.test(15))
        fillHistograms<kFDD, kBCSL>(timeSinceSOF, localBC);
    }
  }
  PROCESS_SWITCH(LumiStabilityLightIons, processSLBunches, "process trigger counting of TVX and FDD for bunches without preceding single-arm activity", false);

  void process(MyBCs const& bcs, aod::FT0s const&)
  {
    for (const auto& bc : bcs) {

      if (bc.timestamp() == 0)
        continue;

      setLHCIFData(bc);

      float timeSinceSOF = getTimeSinceSOF(bc);

      if (bc.selection_bit(aod::evsel::kIsTriggerTVX))
        mHistManager.fill(HIST("FT0Vtx_EvSel/nBCsVsTime"), timeSinceSOF);

      int64_t globalBC = bc.globalBC();
      int localBC = globalBC % nBCsPerOrbit;

      int64_t thisTFid = (globalBC - bcSOR) / nBCsPerTF;

      if (thisTFid != currentTFid) {
        currentTFid = thisTFid;
        mHistManager.fill(HIST("TFsPerMinute"), timeSinceSOF);
      }

      std::bitset<64> ctpInputMask(bc.inputMask());

      for (int iTrigger = 0; iTrigger < nTriggers; iTrigger++) {
        if ((iTrigger == kAllBCs) || (iTrigger == kFT0Vtx && cfgDoFT0Vtx) || (iTrigger == kFT0CE && cfgDoFT0CE) || (iTrigger == kFDD && cfgDoFDD) || (iTrigger == k1ZNC && cfgDo1ZNC)) {
          for (int iBCCategory = 0; iBCCategory < nBCCategories - 1; iBCCategory++) { // Don't do SL BCs here
            if ((iBCCategory == kBCA && cfgDoBCA) || (iBCCategory == kBCB && cfgDoBCB) || (iBCCategory == kBCC && cfgDoBCC) || (iBCCategory == kBCE && cfgDoBCE) || (iBCCategory == kBCL && cfgDoBCL)) {
              if (iTrigger == kAllBCs) {
                if (iBCCategory == kBCA && bcPatternA[localBC])
                  fillHistograms<kAllBCs, kBCA>(timeSinceSOF, localBC);
                if (iBCCategory == kBCB && bcPatternB[localBC])
                  fillHistograms<kAllBCs, kBCB>(timeSinceSOF, localBC);
                if (iBCCategory == kBCC && bcPatternC[localBC])
                  fillHistograms<kAllBCs, kBCC>(timeSinceSOF, localBC);
                if (iBCCategory == kBCE && bcPatternE[localBC])
                  fillHistograms<kAllBCs, kBCE>(timeSinceSOF, localBC);
                if (iBCCategory == kBCL && bcPatternL[localBC])
                  fillHistograms<kAllBCs, kBCL>(timeSinceSOF, localBC);
              }
              if (iTrigger == kFT0Vtx && ctpInputMask.test(2)) {
                if (iBCCategory == kBCA && bcPatternA[localBC])
                  fillHistograms<kFT0Vtx, kBCA>(timeSinceSOF, localBC);
                if (iBCCategory == kBCB && bcPatternB[localBC])
                  fillHistograms<kFT0Vtx, kBCB>(timeSinceSOF, localBC);
                if (iBCCategory == kBCC && bcPatternC[localBC])
                  fillHistograms<kFT0Vtx, kBCC>(timeSinceSOF, localBC);
                if (iBCCategory == kBCE && bcPatternE[localBC])
                  fillHistograms<kFT0Vtx, kBCE>(timeSinceSOF, localBC);
                if (iBCCategory == kBCL && bcPatternL[localBC])
                  fillHistograms<kFT0Vtx, kBCL>(timeSinceSOF, localBC);
              }
              if (iTrigger == kFT0CE && ctpInputMask.test(4)) {
                if (iBCCategory == kBCA && bcPatternA[localBC])
                  fillHistograms<kFT0CE, kBCA>(timeSinceSOF, localBC);
                if (iBCCategory == kBCB && bcPatternB[localBC])
                  fillHistograms<kFT0CE, kBCB>(timeSinceSOF, localBC);
                if (iBCCategory == kBCC && bcPatternC[localBC])
                  fillHistograms<kFT0CE, kBCC>(timeSinceSOF, localBC);
                if (iBCCategory == kBCE && bcPatternE[localBC])
                  fillHistograms<kFT0CE, kBCE>(timeSinceSOF, localBC);
                if (iBCCategory == kBCL && bcPatternL[localBC])
                  fillHistograms<kFT0CE, kBCL>(timeSinceSOF, localBC);
              }
              if (iTrigger == kFDD && ctpInputMask.test(15)) {
                if (iBCCategory == kBCA && bcPatternA[localBC])
                  fillHistograms<kFDD, kBCA>(timeSinceSOF, localBC);
                if (iBCCategory == kBCB && bcPatternB[localBC])
                  fillHistograms<kFDD, kBCB>(timeSinceSOF, localBC);
                if (iBCCategory == kBCC && bcPatternC[localBC])
                  fillHistograms<kFDD, kBCC>(timeSinceSOF, localBC);
                if (iBCCategory == kBCE && bcPatternE[localBC])
                  fillHistograms<kFDD, kBCE>(timeSinceSOF, localBC);
                if (iBCCategory == kBCL && bcPatternL[localBC])
                  fillHistograms<kFDD, kBCL>(timeSinceSOF, localBC);
              }
              if (iTrigger == k1ZNC && ctpInputMask.test(25)) {
                if (iBCCategory == kBCA && bcPatternA[localBC])
                  fillHistograms<k1ZNC, kBCA>(timeSinceSOF, localBC);
                if (iBCCategory == kBCB && bcPatternB[localBC])
                  fillHistograms<k1ZNC, kBCB>(timeSinceSOF, localBC);
                if (iBCCategory == kBCC && bcPatternC[localBC])
                  fillHistograms<k1ZNC, kBCC>(timeSinceSOF, localBC);
                if (iBCCategory == kBCE && bcPatternE[localBC])
                  fillHistograms<k1ZNC, kBCE>(timeSinceSOF, localBC);
                if (iBCCategory == kBCL && bcPatternL[localBC])
                  fillHistograms<k1ZNC, kBCL>(timeSinceSOF, localBC);
              }
            }
          }
        }
      }
      mHistManager.fill(HIST("nBCsVsBCID"), localBC);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  metadataInfo.initMetadata(cfgc);
  return WorkflowSpec{adaptAnalysisTask<LumiStabilityLightIons>(cfgc)};
}
