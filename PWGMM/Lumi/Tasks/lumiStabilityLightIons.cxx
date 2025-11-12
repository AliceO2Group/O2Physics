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

#include <limits>
#include <map>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

o2::common::core::MetadataHelper metadataInfo; // Metadata helper

using MyBCs = soa::Join<aod::BCs, aod::BcSels, aod::Timestamps, aod::Run3MatchedToBCSparse>;

struct LumiStabilityLightIons {
  Configurable<bool> cfgDoFT0Vtx{"cfgDoFT0Vtx", true, "Create and fill histograms for the FT0 vertex trigger"};
  Configurable<bool> cfgDoFT0CE{"cfgDoFT0CE", true, "Create and fill histograms for the FT0 centrality trigger"};
  Configurable<bool> cfgDoFDD{"cfgDoFDD", true, "Create and fill histograms for the FDD trigger"};
  Configurable<bool> cfgDo1ZNC{"cfgDo1ZNC", true, "Create and fill histograms for the 1ZNC trigger"};

  Configurable<bool> cfgDoBCA{"cfgDoBCA", false, "Create and fill histograms for the BCs of type A"};
  Configurable<bool> cfgDoBCB{"cfgDoBCB", true, "Create and fill histograms for the BCs of type B"};
  Configurable<bool> cfgDoBCC{"cfgDoBCC", false, "Create and fill histograms for the BCs of type C"};
  Configurable<bool> cfgDoBCE{"cfgDoBCE", false, "Create and fill histograms for the BCs of type E"};
  Configurable<bool> cfgDoBCL{"cfgDoBCL", false, "Create and fill histograms for leading BCs of type B"};

  Configurable<int> cfgEmptyBCsBeforeLeadingBC{"cfgEmptyBCsBeforeLeadingBC", 5, "Minimum number of empty BCs before a leading BC to identify it as such"};

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
  const int nBCCategories = 5;
  enum BCCategories { kBCA = 0,
                      kBCB = 1,
                      kBCC = 2,
                      kBCE = 3,
                      kBCL = 4 };

  static constexpr std::string_view NBCsVsTimeHistNames[5][5] =
    {{"AllBCs/BC_A/nBCsVsTime", "AllBCs/BC_B/nBCsVsTime", "AllBCs/BC_C/nBCsVsTime", "AllBCs/BC_E/nBCsVsTime", "AllBCs/BC_L/nBCsVsTime"},
     {"FT0VTx/BC_A/nBCsVsTime", "FT0VTx/BC_B/nBCsVsTime", "FT0VTx/BC_C/nBCsVsTime", "FT0VTx/BC_E/nBCsVsTime", "FT0VTx/BC_L/nBCsVsTime"},
     {"FT0CE/BC_A/nBCsVsTime", "FT0CE/BC_B/nBCsVsTime", "FT0CE/BC_C/nBCsVsTime", "FT0CE/BC_E/nBCsVsTime", "FT0CE/BC_L/nBCsVsTime"},
     {"FDD/BC_A/nBCsVsTime", "FDD/BC_B/nBCsVsTime", "FDD/BC_C/nBCsVsTime", "FDD/BC_E/nBCsVsTime", "FDD/BC_L/nBCsVsTime"},
     {"1ZNC/BC_A/nBCsVsTime", "1ZNC/BC_B/nBCsVsTime", "1ZNC/BC_C/nBCsVsTime", "1ZNC/BC_E/nBCsVsTime", "1ZNC/BC_L/nBCsVsTime"}};

  static constexpr std::string_view NBCsVsBCIDHistNames[5][5] =
    {{"AllBCs/BC_A/nBCsVsBCID", "AllBCs/BC_B/nBCsVsBCID", "AllBCs/BC_C/nBCsVsBCID", "AllBCs/BC_E/nBCsVsBCID", "AllBCs/BC_L/nBCsVsBCID"},
     {"FT0VTx/BC_A/nBCsVsBCID", "FT0VTx/BC_B/nBCsVsBCID", "FT0VTx/BC_C/nBCsVsBCID", "FT0VTx/BC_E/nBCsVsBCID", "FT0VTx/BC_L/nBCsVsBCID"},
     {"FT0CE/BC_A/nBCsVsBCID", "FT0CE/BC_B/nBCsVsBCID", "FT0CE/BC_C/nBCsVsBCID", "FT0CE/BC_E/nBCsVsBCID", "FT0CE/BC_L/nBCsVsBCID"},
     {"FDD/BC_A/nBCsVsBCID", "FDD/BC_B/nBCsVsBCID", "FDD/BC_C/nBCsVsBCID", "FDD/BC_E/nBCsVsBCID", "FDD/BC_L/nBCsVsBCID"},
     {"1ZNC/BC_A/nBCsVsBCID", "1ZNC/BC_B/nBCsVsBCID", "1ZNC/BC_C/nBCsVsBCID", "1ZNC/BC_E/nBCsVsBCID", "1ZNC/BC_L/nBCsVsBCID"}};

  int64_t bcSOR;
  int nBCsPerTF;
  int64_t currentTFid = -1;

  void init(InitContext&)
  {
    strLPMProductionTag = metadataInfo.get("LPMProductionTag"); // to extract info from ccdb by the tag

    LOG(info) << "strLPMProductionTag: " << strLPMProductionTag;

    AxisSpec timeAxis{1200, 0., 1200., "#bf{t-t_{SOF} (min)}"}, bcIDAxis{3600, 0., 3600., "#bf{BC ID in orbit}"};

    for (int iTrigger = 0; iTrigger < nTriggers; iTrigger++) {
      if ((iTrigger == kAllBCs) || (iTrigger == kFT0Vtx && cfgDoFT0Vtx) || (iTrigger == kFT0CE && cfgDoFT0CE) || (iTrigger == kFDD && cfgDoFDD) || (iTrigger == k1ZNC && cfgDo1ZNC)) {
        for (int iBCCategory = 0; iBCCategory < nBCCategories; iBCCategory++) {
          if ((iBCCategory == kBCA && cfgDoBCA) || (iBCCategory == kBCB && cfgDoBCB) || (iBCCategory == kBCC && cfgDoBCC) || (iBCCategory == kBCE && cfgDoBCE) || (iBCCategory == kBCL && cfgDoBCL)) {
            mHistManager.add(Form("%s", std::string(NBCsVsTimeHistNames[iTrigger][iBCCategory]).c_str()), "Time of triggered BCs since the start of fill;#bf{t-t_{SOF} (min)};#bf{#it{N}_{BC}}", HistType::kTH1F, {timeAxis});
            mHistManager.add(Form("%s", std::string(NBCsVsBCIDHistNames[iTrigger][iBCCategory]).c_str()), "BC ID of triggered BCs;#bf{BC ID in orbit};#bf{#it{N}_{BC}}", HistType::kTH1F, {bcIDAxis});
          }
        }
      }
    }

    mHistManager.add("FT0Vtx_EvSel/nBCsVsTime", "Time of TVX triggered BCs since the start of fill;;#bf{#it{N}_{BC}}", HistType::kTH1F, {timeAxis});
    mHistManager.add("nBCsVsBCID", "Time of TVX triggered BCs since the start of fill;#bf{t-t_{SOF} (min)};#bf{#it{N}_{BC}}", HistType::kTH1F, {bcIDAxis});
    mHistManager.add("TFsPerMinute", "TFs seen in this minute (to account for failed jobs);#bf{t-t_{SOF} (min)};#bf{#it{N}_{TFs}}", HistType::kTH1F, {timeAxis});
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
          for (int iBCCategory = 0; iBCCategory < nBCCategories; iBCCategory++) {
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
