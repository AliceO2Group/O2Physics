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
/// \file lumiStabilityPP.cxx
/// \brief Analysis over BCs to study the luminosity stability along time for pp collisions
///
/// \author Fabrizio Grosa (fabrizio.grosa@cern.ch), CERN

#include "Common/CCDB/ctpRateFetcher.h"
#include "Common/Core/MetadataHelper.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FT0Corrected.h"

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsFT0/Digit.h"
#include "DataFormatsParameters/AggregatedRunInfo.h"
#include "DataFormatsParameters/GRPLHCIFData.h"
#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include <limits>
#include <map>
#include <memory>
#include <string>
#include <vector>

o2::common::core::MetadataHelper metadataInfo; // Metadata helper

namespace o2
{
namespace lumi
{
enum TriggerAliases { AllBCs = 0,
                      FT0Vtx = 1,
                      FT0CE = 2,
                      FDD = 3,
                      NTriggerAliases };
enum BCCategories { BCA = 0,  // A side BCs (bunch-crossings that had beam only from A side)
                    BCB = 1,  // B type BCs (bunch-crossings that had beam from both sides)
                    BCC = 2,  // C side BCs (bunch-crossings that had beam only from C side)
                    BCE = 3,  // empty BCs (bunch-crossings that did not have beam from either side)
                    BCL = 4,  // leading BCs (bunch-crossings that did not have interacting bunches for a configurable number of preceding BCs)
                    BCSL = 5, // super-leading BCs (bunch-crossings that did not have FDD/FT0 activity for a configurable number of preceding BCs)
                    NBCCategories };
} // namespace lumi
namespace aod
{
// Columns to store the information about the presence of FT0 and FDD signals associated to a given BC
DECLARE_SOA_TABLE(BcDetectorInfo, "AOD", "BCDETECTORINFO", //!
                  indices::FT0Id,
                  indices::FDDId);
} // namespace aod
} // namespace o2

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::lumi;

using BCsWithTimeStamps = soa::Join<aod::BCs, aod::Timestamps, aod::BcDetectorInfo>;

struct BuildBcFlagTable {

  Produces<aod::BcDetectorInfo> bcFlags;

  void init(InitContext&) {}

  void process(aod::BC const& bc,
               aod::FT0s const& ft0s,
               aod::FDDs const& fdds)
  {
    int64_t idxFT0{-1}, idxFDD{-1};
    for (const auto& ft0 : ft0s) {
      if (ft0.bcId() == bc.globalIndex()) {
        idxFT0 = ft0.globalIndex();
        break;
      }
    }
    for (const auto& fdd : fdds) {
      if (fdd.bcId() == bc.globalIndex()) {
        idxFDD = fdd.globalIndex();
        break;
      }
    }
    bcFlags(idxFT0, idxFDD);
  }
};

struct LumiStabilityPP {

  Configurable<bool> doBCA{"doBCA", false, "Create and fill histograms for the BCs of type A"};
  Configurable<bool> doBCB{"doBCB", true, "Create and fill histograms for the BCs of type B"};
  Configurable<bool> doBCC{"doBCC", false, "Create and fill histograms for the BCs of type C"};
  Configurable<bool> doBCE{"doBCE", false, "Create and fill histograms for the BCs of type E"};
  Configurable<bool> doBCL{"doBCL", false, "Create and fill histograms for leading BCs of type B"};
  Configurable<bool> doBCSL{"doBCSL", false, "Create and fill histograms for super-leading BCs (no preceding FT0/FDD activity) of type B"};
  Configurable<int> numEmptyBCsBeforeLeadingBC{"numEmptyBCsBeforeLeadingBC", 5, "Number of empty BCs before a leading BC"};
  Configurable<bool> requireNoT0ForSLBC{"requireNoT0ForSLBC", false, "Require no T0 signal for definition of super leading BC (otherwise only no FDD)"};

  std::bitset<o2::constants::lhc::LHCMaxBunches> beamPatternA, beamPatternC;
  std::bitset<o2::constants::lhc::LHCMaxBunches> bcPatternA, bcPatternC, bcPatternB, bcPatternE, bcPatternL;
  const int nBCsPerOrbit = o2::constants::lhc::LHCMaxBunches;

  parameters::GRPLHCIFData* mLHCIFdata = nullptr;
  int runNumber{-1};
  ctpRateFetcher mRateFetcher;
  std::string injectionScheme;

  HistogramRegistry registry{"registry"};

  std::array<std::array<std::map<int, std::shared_ptr<TH1>>, NBCCategories>, NTriggerAliases> histBcVsTime;
  std::array<std::array<std::map<int, std::shared_ptr<TH1>>, NBCCategories>, NTriggerAliases> histBcVsBcId;
  std::map<int, std::shared_ptr<TH1>> histNBcsVsTime;
  std::map<int, std::shared_ptr<TH1>> histNBcsVsBcId;
  std::map<int, std::shared_ptr<TH1>> histTfPerMin;
  std::map<int, std::shared_ptr<TH1>> histBcHasFT0;
  std::map<int, std::shared_ptr<TH1>> histBcHasFDD;
  std::map<int, std::shared_ptr<TH1>> histFillingScheme;

  static constexpr std::string_view NBCsVsTimeHistNames[NTriggerAliases][NBCCategories] =
    {{"AllBCs/BC_A/nBCsVsTime", "AllBCs/BC_B/nBCsVsTime", "AllBCs/BC_C/nBCsVsTime", "AllBCs/BC_E/nBCsVsTime", "AllBCs/BC_L/nBCsVsTime", "AllBCs/BC_SL/nBCsVsTime"},
     {"FT0VTx/BC_A/nBCsVsTime", "FT0VTx/BC_B/nBCsVsTime", "FT0VTx/BC_C/nBCsVsTime", "FT0VTx/BC_E/nBCsVsTime", "FT0VTx/BC_L/nBCsVsTime", "FT0VTx/BC_SL/nBCsVsTime"},
     {"FT0CE/BC_A/nBCsVsTime", "FT0CE/BC_B/nBCsVsTime", "FT0CE/BC_C/nBCsVsTime", "FT0CE/BC_E/nBCsVsTime", "FT0CE/BC_L/nBCsVsTime", "FT0CE/BC_SL/nBCsVsTime"},
     {"FDD/BC_A/nBCsVsTime", "FDD/BC_B/nBCsVsTime", "FDD/BC_C/nBCsVsTime", "FDD/BC_E/nBCsVsTime", "FDD/BC_L/nBCsVsTime", "FDD/BC_SL/nBCsVsTime"}};

  static constexpr std::string_view NBCsVsBCIDHistNames[NTriggerAliases][NBCCategories] =
    {{"AllBCs/BC_A/nBCsVsBCID", "AllBCs/BC_B/nBCsVsBCID", "AllBCs/BC_C/nBCsVsBCID", "AllBCs/BC_E/nBCsVsBCID", "AllBCs/BC_L/nBCsVsBCID", "AllBCs/BC_SL/nBCsVsBCID"},
     {"FT0VTx/BC_A/nBCsVsBCID", "FT0VTx/BC_B/nBCsVsBCID", "FT0VTx/BC_C/nBCsVsBCID", "FT0VTx/BC_E/nBCsVsBCID", "FT0VTx/BC_L/nBCsVsBCID", "FT0VTx/BC_SL/nBCsVsBCID"},
     {"FT0CE/BC_A/nBCsVsBCID", "FT0CE/BC_B/nBCsVsBCID", "FT0CE/BC_C/nBCsVsBCID", "FT0CE/BC_E/nBCsVsBCID", "FT0CE/BC_L/nBCsVsBCID", "FT0CE/BC_SL/nBCsVsBCID"},
     {"FDD/BC_A/nBCsVsBCID", "FDD/BC_B/nBCsVsBCID", "FDD/BC_C/nBCsVsBCID", "FDD/BC_E/nBCsVsBCID", "FDD/BC_L/nBCsVsBCID", "FDD/BC_SL/nBCsVsBCID"}};

  const AxisSpec timeAxis{2880, 0., 2880., "#bf{t-t_{SOF} (min)}"}, bcIDAxis{nBCsPerOrbit, -0.5, static_cast<float>(nBCsPerOrbit) - 0.5, "#bf{BC ID in orbit}"};

  int64_t bcSOR;
  int nBCsPerTF;
  int64_t currentTFid = -1;

  void init(InitContext&) {}

  void createHistograms()
  {
    if (histNBcsVsTime[runNumber]) { // histograms for this run already there
      return;
    }

    histNBcsVsTime[runNumber] = registry.add<TH1>(Form("%d/FT0Vtx_EvSel/nBCsVsTime", runNumber), "Time of TVX triggered BCs since the start of fill;;#bf{#it{N}_{BC}}", HistType::kTH1D, {timeAxis});
    histNBcsVsBcId[runNumber] = registry.add<TH1>(Form("%d/nBCsVsBCID", runNumber), "Time of TVX triggered BCs since the start of fill;#bf{t-t_{SOF} (min)};#bf{#it{N}_{BC}}", HistType::kTH1D, {bcIDAxis});
    histTfPerMin[runNumber] = registry.add<TH1>(Form("%d/TFsPerMinute", runNumber), "TFs seen in this minute (to account for failed jobs);#bf{t-t_{SOF} (min)};#bf{#it{N}_{TFs}}", HistType::kTH1D, {timeAxis});
    histFillingScheme[runNumber] = registry.add<TH1>(Form("%d/FillingScheme", runNumber), "Filling Scheme;Filling Scheme;", HistType::kTH1D, {{1, 0, 1}});

    histBcHasFT0[runNumber] = registry.add<TH2>(Form("%d/FITQA/BCHasFT0", runNumber), "Does the BC have FT0?;BC has FT0;TVX triggered according to CTP;#bf{#it{N}_{BC}}", HistType::kTH2D, {{2, -0.5, 1.5}, {2, -0.5, 1.5}});
    histBcHasFT0[runNumber]->GetYaxis()->SetBinLabel(1, "No CTP trigger");
    histBcHasFT0[runNumber]->GetYaxis()->SetBinLabel(2, "CTP triggered");
    histBcHasFT0[runNumber]->GetXaxis()->SetBinLabel(1, "No found FT0");
    histBcHasFT0[runNumber]->GetXaxis()->SetBinLabel(2, "Found FT0");
    histBcHasFDD[runNumber] = registry.add<TH2>(Form("%d/FITQA/BCHasFDD", runNumber), "Does the BC have FDD?;BC has FDD;FDD triggered according to CTP;#bf{#it{N}_{BC}}", HistType::kTH2D, {{2, -0.5, 1.5}, {2, -0.5, 1.5}});
    histBcHasFDD[runNumber]->GetYaxis()->SetBinLabel(1, "No CTP trigger");
    histBcHasFDD[runNumber]->GetYaxis()->SetBinLabel(2, "CTP triggered");
    histBcHasFDD[runNumber]->GetXaxis()->SetBinLabel(1, "No found FDD");
    histBcHasFDD[runNumber]->GetXaxis()->SetBinLabel(2, "Found FDD");

    for (int iTrigger{0}; iTrigger < NTriggerAliases; ++iTrigger) {
      for (int iBCCategory{0}; iBCCategory < NBCCategories; ++iBCCategory) { // Don't do SL BCs here
        if ((iBCCategory == BCA && doBCA) || (iBCCategory == BCB && doBCB) || (iBCCategory == BCC && doBCC) || (iBCCategory == BCE && doBCE) || (iBCCategory == BCL && doBCL) || (iBCCategory == BCSL && doBCSL)) {
          histBcVsTime[iTrigger][iBCCategory][runNumber] = registry.add<TH1>(Form("%d/%s", runNumber, std::string(NBCsVsTimeHistNames[iTrigger][iBCCategory]).c_str()), "Time of triggered BCs since the start of fill;#bf{t-t_{SOF} (min)};#bf{#it{N}_{BC}}", HistType::kTH1D, {timeAxis});
          histBcVsBcId[iTrigger][iBCCategory][runNumber] = registry.add<TH1>(Form("%d/%s", runNumber, std::string(NBCsVsBCIDHistNames[iTrigger][iBCCategory]).c_str()), "BC ID of triggered BCs;#bf{BC ID in orbit};#bf{#it{N}_{BC}}", HistType::kTH1D, {bcIDAxis});
        }
      }
    }
  }

  void setLHCIFData(const auto& bc)
  {
    if (runNumber == bc.runNumber()) {
      return;
    }

    auto& ccdbMgr = o2::ccdb::BasicCCDBManager::instance();
    uint64_t timeStamp = bc.timestamp();

    std::map<std::string, std::string> metadata;
    mLHCIFdata = ccdbMgr.getSpecific<o2::parameters::GRPLHCIFData>("GLO/Config/GRPLHCIF", timeStamp, metadata);
    if (mLHCIFdata == nullptr) {
      LOG(fatal) << "GRPLHCIFData not in database, timestamp:" << timeStamp;
    }

    runNumber = bc.runNumber();
    LOG(info) << "LHCIF data fetched for run " << runNumber << " and timestamp " << timeStamp;
    createHistograms();

    histFillingScheme[runNumber]->Fill(mLHCIFdata->getInjectionScheme().c_str(), 0);

    beamPatternA = mLHCIFdata->getBunchFilling().getBeamPattern(0);
    beamPatternC = mLHCIFdata->getBunchFilling().getBeamPattern(1);
    bcPatternA = beamPatternA & ~beamPatternC;
    bcPatternC = ~beamPatternA & beamPatternC;
    bcPatternB = beamPatternA & beamPatternC;
    bcPatternE = ~beamPatternA & ~beamPatternC;

    // Create bcPatternL: leading BCs of type B that follow at least "numEmptyBCsBeforeLeadingBC" empty BCs
    bcPatternL.reset(); // Initialize all bits to false
    LOG(info) << "Starting to create bcPatternL from bcPatternB";
    LOG(info) << "Total number of BCs to check: " << o2::constants::lhc::LHCMaxBunches;

    int totalLeadingBCs = 0;
    for (int iBC = 0; iBC < o2::constants::lhc::LHCMaxBunches; iBC++) {
      if (bcPatternB[iBC]) {    // Check if current BC is of type B
        int emptyBCsBefore = 0; // Count how many consecutive BCs before this one are NOT type B
        for (int j = 1; j <= numEmptyBCsBeforeLeadingBC; j++) {
          int prevBC = (iBC - j + o2::constants::lhc::LHCMaxBunches) % o2::constants::lhc::LHCMaxBunches; // Protection for BCs at small indices to check the end of the orbit
          if (!bcPatternB[prevBC]) {
            emptyBCsBefore++;
          } else {
            break; // Stop counting if we hit a type B BC
          }
        }
        if (emptyBCsBefore >= numEmptyBCsBeforeLeadingBC) { // If we found at least numEmptyBCsBeforeLeadingBC empty BCs before this one, mark it as leading
          bcPatternL[iBC] = true;
          totalLeadingBCs++;
        }
      }
    }
    LOG(info) << "bcPatternL creation complete. Total leading BCs found: " << totalLeadingBCs;

    auto runInfo = o2::parameters::AggregatedRunInfo::buildAggregatedRunInfo(o2::ccdb::BasicCCDBManager::instance(), runNumber, metadataInfo.get("LPMProductionTag"));
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
    histBcVsTime[iTrigger][iBCCategory][runNumber]->Fill(timeSinceSOF);
    histBcVsBcId[iTrigger][iBCCategory][runNumber]->Fill(localBC);
  }

  void process(BCsWithTimeStamps const& bcs,
               aod::FT0s const&,
               aod::FDDs const&)
  {
    int64_t globalBCIdOfLastBCWithActivity = 0;
    for (const auto& bc : bcs) {

      if (bc.timestamp() == 0) {
        continue;
      }

      setLHCIFData(bc);

      float timeSinceSOF = getTimeSinceSOF(bc);
      bool isTriggerTVX = (bc.has_ft0() ? TESTBIT(bc.ft0().triggerMask(), o2::ft0::Triggers::bitVertex) : false);

      if (isTriggerTVX) {
        histNBcsVsTime[runNumber]->Fill(timeSinceSOF);
      }

      int64_t globalBC = bc.globalBC();
      int localBC = globalBC % nBCsPerOrbit;

      bool isSuperLeadingBc{true};
      if (globalBC - globalBCIdOfLastBCWithActivity < numEmptyBCsBeforeLeadingBC) {
        isSuperLeadingBc = false; // not a super-leading BC
      }

      if (bc.has_fdd() || (requireNoT0ForSLBC && bc.has_ft0())) {
        globalBCIdOfLastBCWithActivity = globalBC;
      }

      if (!bcPatternB[localBC]) {
        isSuperLeadingBc = false; // not a super-leading BC
      }

      int64_t thisTFid = (globalBC - bcSOR) / nBCsPerTF;

      if (thisTFid != currentTFid) {
        currentTFid = thisTFid;
        histTfPerMin[runNumber]->Fill(timeSinceSOF);
      }

      std::bitset<64> ctpInputMask(bc.inputMask());

      histBcHasFT0[runNumber]->Fill(bc.has_ft0(), ctpInputMask.test(2));
      histBcHasFDD[runNumber]->Fill(bc.has_fdd(), ctpInputMask.test(15));

      for (int iTrigger{0}; iTrigger < NTriggerAliases; ++iTrigger) {
        for (int iBCCategory{0}; iBCCategory < NBCCategories; ++iBCCategory) { // Don't do SL BCs here
          if ((iBCCategory == BCA && doBCA) || (iBCCategory == BCB && doBCB) || (iBCCategory == BCC && doBCC) || (iBCCategory == BCE && doBCE) || (iBCCategory == BCL && doBCL)) {
            if (iTrigger == AllBCs) {
              if (iBCCategory == BCA && bcPatternA[localBC])
                fillHistograms<AllBCs, BCA>(timeSinceSOF, localBC);
              if (iBCCategory == BCB && bcPatternB[localBC])
                fillHistograms<AllBCs, BCB>(timeSinceSOF, localBC);
              if (iBCCategory == BCC && bcPatternC[localBC])
                fillHistograms<AllBCs, BCC>(timeSinceSOF, localBC);
              if (iBCCategory == BCE && bcPatternE[localBC])
                fillHistograms<AllBCs, BCE>(timeSinceSOF, localBC);
              if (iBCCategory == BCL && bcPatternL[localBC])
                fillHistograms<AllBCs, BCL>(timeSinceSOF, localBC);
              if (iBCCategory == BCSL && isSuperLeadingBc)
                fillHistograms<AllBCs, BCSL>(timeSinceSOF, localBC);
            }
            if (iTrigger == FT0Vtx && ctpInputMask.test(2)) {
              if (iBCCategory == BCA && bcPatternA[localBC])
                fillHistograms<FT0Vtx, BCA>(timeSinceSOF, localBC);
              if (iBCCategory == BCB && bcPatternB[localBC])
                fillHistograms<FT0Vtx, BCB>(timeSinceSOF, localBC);
              if (iBCCategory == BCC && bcPatternC[localBC])
                fillHistograms<FT0Vtx, BCC>(timeSinceSOF, localBC);
              if (iBCCategory == BCE && bcPatternE[localBC])
                fillHistograms<FT0Vtx, BCE>(timeSinceSOF, localBC);
              if (iBCCategory == BCL && bcPatternL[localBC])
                fillHistograms<FT0Vtx, BCL>(timeSinceSOF, localBC);
              if (iBCCategory == BCSL && isSuperLeadingBc)
                fillHistograms<FT0Vtx, BCSL>(timeSinceSOF, localBC);
            }
            if (iTrigger == FT0CE && ctpInputMask.test(4)) {
              if (iBCCategory == BCA && bcPatternA[localBC])
                fillHistograms<FT0CE, BCA>(timeSinceSOF, localBC);
              if (iBCCategory == BCB && bcPatternB[localBC])
                fillHistograms<FT0CE, BCB>(timeSinceSOF, localBC);
              if (iBCCategory == BCC && bcPatternC[localBC])
                fillHistograms<FT0CE, BCC>(timeSinceSOF, localBC);
              if (iBCCategory == BCE && bcPatternE[localBC])
                fillHistograms<FT0CE, BCE>(timeSinceSOF, localBC);
              if (iBCCategory == BCL && bcPatternL[localBC])
                fillHistograms<FT0CE, BCL>(timeSinceSOF, localBC);
            }
            if (iTrigger == FDD && ctpInputMask.test(15)) {
              if (iBCCategory == BCA && bcPatternA[localBC])
                fillHistograms<FDD, BCA>(timeSinceSOF, localBC);
              if (iBCCategory == BCB && bcPatternB[localBC])
                fillHistograms<FDD, BCB>(timeSinceSOF, localBC);
              if (iBCCategory == BCC && bcPatternC[localBC])
                fillHistograms<FDD, BCC>(timeSinceSOF, localBC);
              if (iBCCategory == BCE && bcPatternE[localBC])
                fillHistograms<FDD, BCE>(timeSinceSOF, localBC);
              if (iBCCategory == BCL && bcPatternL[localBC])
                fillHistograms<FDD, BCL>(timeSinceSOF, localBC);
              if (iBCCategory == BCSL && isSuperLeadingBc)
                fillHistograms<FDD, BCSL>(timeSinceSOF, localBC);
            }
          }
        }
      }
      histNBcsVsBcId[runNumber]->Fill(localBC);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  metadataInfo.initMetadata(cfgc);
  return WorkflowSpec{adaptAnalysisTask<BuildBcFlagTable>(cfgc), adaptAnalysisTask<LumiStabilityPP>(cfgc)};
}
