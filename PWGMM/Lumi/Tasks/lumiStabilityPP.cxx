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
/// \author Fabrizio Chinu (fabrizio.chinu@cern.ch), INFN and University of Turin

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
#include <Framework/Array2D.h>
#include <Framework/Configurable.h>

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

// SL types must be after all the others
enum BCCategories { BCA = 0,  // A side BCs (bunch-crossings that had beam only from A side)
                    BCB,      // B type BCs (bunch-crossings that had beam from both sides)
                    BCC,      // C side BCs (bunch-crossings that had beam only from C side)
                    BCE,      // empty BCs (bunch-crossings that did not have beam from either side)
                    BCL,      // leading BCs (bunch-crossings that have not-B BCs for a configurable number of preceding BCs)
                    BCLE,     // leading BCs (bunch-crossings that did not have interacting bunches for a configurable number of preceding BCs)
                    BCNL,     // non-leading BCs of type B (bunch-crossings that come after a BCL and are of type B)
                    BCNLE,    // non-leading BCs of type B (bunch-crossings that come after a BCLE and are of type B)
                    BCSLFDD,  // super-leading BCs for FDD (bunch-crossings that had beam from both sides but did not have FDD activity for a configurable number of preceding BCs)
                    BCSLFT0,  // super-leading BCs for FT0 (bunch-crossings that had beam from both sides but did not have FT0 activity for a configurable number of preceding BCs)
                    BCNSLFDD, // non-super-leading BCs for FDD of type B (bunch-crossings that had beam from both sides but are not SL for FDD activity for a configurable number of preceding BCs)
                    BCNSLFT0, // non-super-leading BCs for FT0 of type B (bunch-crossings that had beam from both sides but are not SL for FT0 activity for a configurable number of preceding BCs)
                    NBCCategories };
} // namespace lumi
} // namespace o2

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::lumi;

using BCsWithTimeStamps = soa::Join<aod::BCs, aod::Timestamps>;

struct LumiStabilityPP {

  static constexpr int defaulFlags[1][NBCCategories] = {{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}};
  Configurable<LabeledArray<int>> doTypeBC{"doTypeBC", {defaulFlags[0], NBCCategories, {"BCA", "BCB", "BCC", "BCE", "BCL", "BCLE", "BCSLFDD", "BCSLFT0", "BCNL", "BCNLE", "BCNSLFDD", "BCNSLFT0"}}, "Create and fill histograms for different BC types"};

  static constexpr int defaulNumBCsBeforeLeadingBC[1][3] = {{5, 5, 5}};
  Configurable<LabeledArray<int>> numEmptyBCsBeforeLeadingBC{"numEmptyBCsBeforeLeadingBC", {defaulNumBCsBeforeLeadingBC[0], 3, {"BCL", "BCLE", "BCSL"}}, "Create and fill histograms for different BC types"};
  Configurable<int> bcShiftFDDForData2023{"bcShiftFDDForData2023", 7, "Number of bc to shift for FDD to be applied for 2023 data only"};

  std::bitset<o2::constants::lhc::LHCMaxBunches> beamPatternA, beamPatternC;
  std::bitset<o2::constants::lhc::LHCMaxBunches> bcPatternA, bcPatternC, bcPatternB, bcPatternE, bcPatternL, bcPatternLE;
  const int nBCsPerOrbit = o2::constants::lhc::LHCMaxBunches;

  o2::framework::Service<o2::ccdb::BasicCCDBManager> ccdb;
  parameters::GRPLHCIFData* mLHCIFdata = nullptr;
  int runNumber{-1};
  bool isData23{false};
  ctpRateFetcher mRateFetcher;
  int nBunchesFillingScheme;

  HistogramRegistry registry{"registry"};

  std::array<std::array<std::map<int, std::shared_ptr<TH1>>, NBCCategories>, NTriggerAliases> histBcVsTime;
  std::array<std::array<std::map<int, std::shared_ptr<TH1>>, NBCCategories>, NTriggerAliases> histBcVsBcId;
  std::array<std::array<std::map<int, std::shared_ptr<TH1>>, NBCCategories>, NTriggerAliases> histBcInspectVsBcId;
  std::array<std::map<int, std::shared_ptr<TH1>>, BCSLFDD> histBcPattern; // undefined for BC(N)SL
  std::map<int, std::shared_ptr<TH1>> histNBcsVsTime;
  std::map<int, std::shared_ptr<TH1>> histNBcsVsBcId;
  std::map<int, std::shared_ptr<TH1>> histTfPerMin;
  std::map<int, std::shared_ptr<TH1>> histFillingScheme;
  std::map<int, std::shared_ptr<TH1>> histFillTime;
  std::map<int, std::shared_ptr<TH1>> histInteractionRate;

  static constexpr std::string_view NBCsVsTimeHistNames[NTriggerAliases][NBCCategories] =
    {{"AllBCs/BC_A/nBCsVsTime", "AllBCs/BC_B/nBCsVsTime", "AllBCs/BC_C/nBCsVsTime", "AllBCs/BC_E/nBCsVsTime", "AllBCs/BC_L/nBCsVsTime", "AllBCs/BC_LE/nBCsVsTime", "AllBCs/BC_NL/nBCsVsTime", "AllBCs/BC_NLE/nBCsVsTime", "AllBCs/BC_SL_FDD/nBCsVsTime", "AllBCs/BC_SL_FT0/nBCsVsTime", "AllBCs/BC_NSL_FT0/nBCsVsTime", "AllBCs/BC_NSL_FDD/nBCsVsTime"},
     {"FT0VTx/BC_A/nBCsVsTime", "FT0VTx/BC_B/nBCsVsTime", "FT0VTx/BC_C/nBCsVsTime", "FT0VTx/BC_E/nBCsVsTime", "FT0VTx/BC_L/nBCsVsTime", "FT0VTx/BC_LE/nBCsVsTime", "FT0VTx/BC_NL/nBCsVsTime", "FT0VTx/BC_NLE/nBCsVsTime", "FT0VTx/BC_SL_FDD/nBCsVsTime", "FT0VTx/BC_SL_FT0/nBCsVsTime", "FT0VTx/BC_NSL_FT0/nBCsVsTime", "FT0VTx/BC_NSL_FDD/nBCsVsTime"},
     {"FT0CE/BC_A/nBCsVsTime", "FT0CE/BC_B/nBCsVsTime", "FT0CE/BC_C/nBCsVsTime", "FT0CE/BC_E/nBCsVsTime", "FT0CE/BC_L/nBCsVsTime", "FT0CE/BC_LE/nBCsVsTime", "FT0CE/BC_NL/nBCsVsTime", "FT0CE/BC_NLE/nBCsVsTime", "FT0CE/BC_SL_FDD/nBCsVsTime", "FT0CE/BC_SL_FT0/nBCsVsTime", "FT0CE/BC_NSL_FT0/nBCsVsTime", "FT0CE/BC_NSL_FDD/nBCsVsTime"},
     {"FDD/BC_A/nBCsVsTime", "FDD/BC_B/nBCsVsTime", "FDD/BC_C/nBCsVsTime", "FDD/BC_E/nBCsVsTime", "FDD/BC_L/nBCsVsTime", "FDD/BC_LE/nBCsVsTime", "FDD/BC_NL/nBCsVsTime", "FDD/BC_NLE/nBCsVsTime", "FDD/BC_SL_FDD/nBCsVsTime", "FDD/BC_SL_FT0/nBCsVsTime", "FDD/BC_NSL_FT0/nBCsVsTime", "FDD/BC_NSL_FDD/nBCsVsTime"}};

  static constexpr std::string_view NBCsVsBCIDHistNames[NTriggerAliases][NBCCategories] =
    {{"AllBCs/BC_A/nBCsVsBCID", "AllBCs/BC_B/nBCsVsBCID", "AllBCs/BC_C/nBCsVsBCID", "AllBCs/BC_E/nBCsVsBCID", "AllBCs/BC_L/nBCsVsBCID", "AllBCs/BC_LE/nBCsVsBCID", "AllBCs/BC_NL/nBCsVsBCID", "AllBCs/BC_NLE/nBCsVsBCID", "AllBCs/BC_SL_FDD/nBCsVsBCID", "AllBCs/BC_SL_FT0/nBCsVsBCID", "AllBCs/BC_NSL_FT0/nBCsVsBCID", "AllBCs/BC_NSL_FDD/nBCsVsBCID"},
     {"FT0VTx/BC_A/nBCsVsBCID", "FT0VTx/BC_B/nBCsVsBCID", "FT0VTx/BC_C/nBCsVsBCID", "FT0VTx/BC_E/nBCsVsBCID", "FT0VTx/BC_L/nBCsVsBCID", "FT0VTx/BC_LE/nBCsVsBCID", "FT0VTx/BC_NL/nBCsVsBCID", "FT0VTx/BC_NLE/nBCsVsBCID", "FT0VTx/BC_SL_FDD/nBCsVsBCID", "FT0VTx/BC_SL_FT0/nBCsVsBCID", "FT0VTx/BC_NSL_FT0/nBCsVsBCID", "FT0VTx/BC_NSL_FDD/nBCsVsBCID"},
     {"FT0CE/BC_A/nBCsVsBCID", "FT0CE/BC_B/nBCsVsBCID", "FT0CE/BC_C/nBCsVsBCID", "FT0CE/BC_E/nBCsVsBCID", "FT0CE/BC_L/nBCsVsBCID", "FT0CE/BC_LE/nBCsVsBCID", "FT0CE/BC_NL/nBCsVsBCID", "FT0CE/BC_NLE/nBCsVsBCID", "FT0CE/BC_SL_FDD/nBCsVsBCID", "FT0CE/BC_SL_FT0/nBCsVsBCID", "FT0CE/BC_NSL_FT0/nBCsVsBCID", "FT0CE/BC_NSL_FDD/nBCsVsBCID"},
     {"FDD/BC_A/nBCsVsBCID", "FDD/BC_B/nBCsVsBCID", "FDD/BC_C/nBCsVsBCID", "FDD/BC_E/nBCsVsBCID", "FDD/BC_L/nBCsVsBCID", "FDD/BC_LE/nBCsVsBCID", "FDD/BC_NL/nBCsVsBCID", "FDD/BC_NLE/nBCsVsBCID", "FDD/BC_SL_FDD/nBCsVsBCID", "FDD/BC_SL_FT0/nBCsVsBCID", "FDD/BC_NSL_FT0/nBCsVsBCID", "FDD/BC_NSL_FDD/nBCsVsBCID"}};

  static constexpr std::string_view NBCsInspectVsBCIDHistNames[NTriggerAliases][NBCCategories] =
    {{"AllBCs/BC_A/nBCsInspectedVsBCID", "AllBCs/BC_B/nBCsInspectedVsBCID", "AllBCs/BC_C/nBCsInspectedVsBCID", "AllBCs/BC_E/nBCsInspectedVsBCID", "AllBCs/BC_L/nBCsInspectedVsBCID", "AllBCs/BC_LE/nBCsInspectedVsBCID", "AllBCs/BC_NL/nBCsInspectedVsBCID", "AllBCs/BC_NLE/nBCsInspectedVsBCID", "AllBCs/BC_SL_FDD/nBCsInspectedVsBCID", "AllBCs/BC_SL_FT0/nBCsInspectedVsBCID", "AllBCs/BC_NSL_FT0/nBCsInspectedVsBCID", "AllBCs/BC_NSL_FDD/nBCsInspectedVsBCID"},
     {"FT0VTx/BC_A/nBCsInspectedVsBCID", "FT0VTx/BC_B/nBCsInspectedVsBCID", "FT0VTx/BC_C/nBCsInspectedVsBCID", "FT0VTx/BC_E/nBCsInspectedVsBCID", "FT0VTx/BC_L/nBCsInspectedVsBCID", "FT0VTx/BC_LE/nBCsInspectedVsBCID", "FT0VTx/BC_NL/nBCsInspectedVsBCID", "FT0VTx/BC_NLE/nBCsInspectedVsBCID", "FT0VTx/BC_SL_FDD/nBCsInspectedVsBCID", "FT0VTx/BC_SL_FT0/nBCsInspectedVsBCID", "FT0VTx/BC_NSL_FT0/nBCsInspectedVsBCID", "FT0VTx/BC_NSL_FDD/nBCsInspectedVsBCID"},
     {"FT0CE/BC_A/nBCsInspectedVsBCID", "FT0CE/BC_B/nBCsInspectedVsBCID", "FT0CE/BC_C/nBCsInspectedVsBCID", "FT0CE/BC_E/nBCsInspectedVsBCID", "FT0CE/BC_L/nBCsInspectedVsBCID", "FT0CE/BC_LE/nBCsInspectedVsBCID", "FT0CE/BC_NL/nBCsInspectedVsBCID", "FT0CE/BC_NLE/nBCsInspectedVsBCID", "FT0CE/BC_SL_FDD/nBCsInspectedVsBCID", "FT0CE/BC_SL_FT0/nBCsInspectedVsBCID", "FT0CE/BC_NSL_FT0/nBCsInspectedVsBCID", "FT0CE/BC_NSL_FDD/nBCsInspectedVsBCID"},
     {"FDD/BC_A/nBCsInspectedVsBCID", "FDD/BC_B/nBCsInspectedVsBCID", "FDD/BC_C/nBCsInspectedVsBCID", "FDD/BC_E/nBCsInspectedVsBCID", "FDD/BC_L/nBCsInspectedVsBCID", "FDD/BC_LE/nBCsInspectedVsBCID", "FDD/BC_NL/nBCsInspectedVsBCID", "FDD/BC_NLE/nBCsInspectedVsBCID", "FDD/BC_SL_FDD/nBCsInspectedVsBCID", "FDD/BC_SL_FT0/nBCsInspectedVsBCID", "FDD/BC_NSL_FT0/nBCsInspectedVsBCID", "FDD/BC_NSL_FDD/nBCsInspectedVsBCID"}};

  static constexpr std::string_view PatternHistNames[BCSLFDD] = {"BC_A/BcPattern", "BC_B/BcPattern", "BC_C/BcPattern", "BC_E/BcPattern", "BC_L/BcPattern", "BC_LE/BcPattern", "BC_NL/BcPattern", "BC_NLE/BcPattern"};

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
    histFillTime[runNumber] = registry.add<TH1>(Form("%d/FillTime", runNumber), "Fill time;Fill time;", HistType::kTH1D, {{1, 0, 1}});
    histInteractionRate[runNumber] = registry.add<TH1>(Form("%d/InteractionRate", runNumber), "Interaction rate (kHz);Interaction rate (kHz);", HistType::kTH1D, {{3000, 0., 3000.}});

    for (int iTrigger{0}; iTrigger < NTriggerAliases; ++iTrigger) {
      for (int iBCCategory{0}; iBCCategory < NBCCategories; ++iBCCategory) {
        if (doTypeBC->get(0u, iBCCategory)) {
          histBcVsTime[iTrigger][iBCCategory][runNumber] = registry.add<TH1>(Form("%d/%s", runNumber, std::string(NBCsVsTimeHistNames[iTrigger][iBCCategory]).c_str()), "Time of triggered BCs since the start of fill;#bf{t-t_{SOF} (min)};#bf{#it{N}_{BC}}", HistType::kTH1D, {timeAxis});
          histBcVsBcId[iTrigger][iBCCategory][runNumber] = registry.add<TH1>(Form("%d/%s", runNumber, std::string(NBCsVsBCIDHistNames[iTrigger][iBCCategory]).c_str()), "BC ID of triggered BCs;#bf{BC ID in orbit};#bf{#it{N}_{BC}}", HistType::kTH1D, {bcIDAxis});
          histBcInspectVsBcId[iTrigger][iBCCategory][runNumber] = registry.add<TH1>(Form("%d/%s", runNumber, std::string(NBCsInspectVsBCIDHistNames[iTrigger][iBCCategory]).c_str()), "BC ID of inspecred BCs;#bf{BC ID in orbit};#bf{#it{N}_{BC}}", HistType::kTH1D, {bcIDAxis});
        }
      }
    }

    for (int iBCCategory{0}; iBCCategory < BCSLFDD; ++iBCCategory) {
      histBcPattern[iBCCategory][runNumber] = registry.add<TH1>(Form("%d/%s", runNumber, std::string(PatternHistNames[iBCCategory]).c_str()), "BC Pattern;#bf{BC ID in orbit};", HistType::kTH1D, {bcIDAxis});
    }
  }

  void setLHCIFData(const auto& bc)
  {

    if (runNumber == bc.runNumber()) {
      return;
    }

    const int runStart2023{535069};
    const int runStop2023{543113};
    if (bc.runNumber() >= runStart2023 && bc.runNumber() <= runStop2023) {
      isData23 = true;
    }

    uint64_t timeStamp = bc.timestamp();

    std::map<std::string, std::string> metadata;
    mLHCIFdata = ccdb.service->getSpecific<o2::parameters::GRPLHCIFData>("GLO/Config/GRPLHCIF", timeStamp, metadata);
    if (mLHCIFdata == nullptr) {
      LOG(fatal) << "GRPLHCIFData not in database, timestamp:" << timeStamp;
    }

    runNumber = bc.runNumber();
    LOG(info) << "LHCIF data fetched for run " << runNumber << " and timestamp " << timeStamp;
    createHistograms();

    std::string_view injectionScheme = mLHCIFdata->getInjectionScheme();
    size_t underScorePos = injectionScheme.find('_');
    size_t bPos = injectionScheme.find('b', underScorePos);
    if (underScorePos != std::string_view::npos && bPos != std::string_view::npos && bPos > underScorePos) {
      std::string_view nBunchesFillingSchemeStr = injectionScheme.substr(underScorePos + 1, bPos - (underScorePos + 1));
      nBunchesFillingScheme = std::stoi(std::string(nBunchesFillingSchemeStr));
    }
    histFillingScheme[runNumber]->Fill(std::string(injectionScheme).c_str(), 0);
    histFillTime[runNumber]->Fill(0.5, mLHCIFdata->getFillNumberTime());

    beamPatternA = mLHCIFdata->getBunchFilling().getBeamPattern(0);
    beamPatternC = mLHCIFdata->getBunchFilling().getBeamPattern(1);
    bcPatternA = beamPatternA & ~beamPatternC;
    bcPatternC = ~beamPatternA & beamPatternC;
    bcPatternB = beamPatternA & beamPatternC;
    bcPatternE = ~beamPatternA & ~beamPatternC;

    // Create bcPatternL: leading BCs of type B that follow at least "numEmptyBCsBeforeLeadingBC" non-B BCs
    bcPatternL.reset(); // Initialize all bits to false
    bcPatternLE.reset();
    LOG(info) << "Starting to create bcPatternL from bcPatternB";
    LOG(info) << "Total number of BCs to check: " << o2::constants::lhc::LHCMaxBunches;

    std::array<int, 2> totalLeadingBCs = {0, 0};
    for (int iBC = 0; iBC < o2::constants::lhc::LHCMaxBunches; iBC++) {
      if (bcPatternB[iBC]) {    // Check if current BC is of type B
        int nonBtypeBCsBefore{0}, emptyBCsBefore{0}; // Count how many consecutive BCs before this one are non-B
        for (int j = 1; j <= numEmptyBCsBeforeLeadingBC->get(0u, 0u); j++) {
          int prevBC = (iBC - j + o2::constants::lhc::LHCMaxBunches) % o2::constants::lhc::LHCMaxBunches; // Protection for BCs at small indices to check the end of the orbit
          if (!bcPatternB[prevBC]) {
            nonBtypeBCsBefore++;
          } else {
            break; // Stop counting if we hit a BCB
          }
        }
        for (int j = 1; j <= numEmptyBCsBeforeLeadingBC->get(0u, 1u); j++) {
          int prevBC = (iBC - j + o2::constants::lhc::LHCMaxBunches) % o2::constants::lhc::LHCMaxBunches; // Protection for BCs at small indices to check the end of the orbit
          if (bcPatternE[prevBC]) {
            emptyBCsBefore++;
          } else {
            break; // Stop counting if we hit a non BCE
          }
        }
        if (nonBtypeBCsBefore >= numEmptyBCsBeforeLeadingBC->get(0u, 0u)) { // If we found at least numEmptyBCsBeforeLeadingBC[0] non-B BCs before this one, mark it as leading
          bcPatternL[iBC] = true;
          totalLeadingBCs[0]++;
        }
        if (emptyBCsBefore >= numEmptyBCsBeforeLeadingBC->get(0u, 1u)) { // If we found at least numEmptyBCsBeforeLeadingBC[1] empty BCs before this one, mark it as leading
          bcPatternLE[iBC] = true;
          totalLeadingBCs[1]++;
        }
      }
      if (bcPatternA[iBC]) {
        histBcPattern[BCA][runNumber]->Fill(iBC);
      }
      if (bcPatternB[iBC]) {
        histBcPattern[BCB][runNumber]->Fill(iBC);
        if (!bcPatternL[iBC]) {
          histBcPattern[BCNL][runNumber]->Fill(iBC);
        }
        if (!bcPatternLE[iBC]) {
          histBcPattern[BCNLE][runNumber]->Fill(iBC);
        }
      }
      if (bcPatternC[iBC]) {
        histBcPattern[BCC][runNumber]->Fill(iBC);
      }
      if (bcPatternE[iBC]) {
        histBcPattern[BCE][runNumber]->Fill(iBC);
      }
      if (bcPatternL[iBC]) {
        histBcPattern[BCL][runNumber]->Fill(iBC);
      }
      if (bcPatternLE[iBC]) {
        histBcPattern[BCLE][runNumber]->Fill(iBC);
      }
    }
    LOG(info) << "bcPatternL creation complete. Total leading BCs found: " << totalLeadingBCs[0];
    LOG(info) << "bcPatternLE creation complete. Total leading BCs found: " << totalLeadingBCs[1];

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
    int64_t globalBCIdOfLastBCWithActivityFDD{0}, globalBCIdOfLastBCWithActivityFT0{0}, globalBCLastInspectedBC{-1};
    std::vector<std::array<int, NBCCategories>> nBCsPerBcId;
    nBCsPerBcId.resize(nBCsPerOrbit);
    std::fill(&nBCsPerBcId[0][0], &nBCsPerBcId[0][0] + (static_cast<int>(nBCsPerOrbit) * static_cast<int>(NBCCategories)), 0); // Initialize to 0

    double rate{-1.};
    for (const auto& bc : bcs) {

      if (bc.timestamp() == 0) {
        continue;
      }

      setLHCIFData(bc);
      int bcShiftFDD{0};
      if (isData23) {
        bcShiftFDD = bcShiftFDDForData2023;
      } else {
        bcShiftFDD = 0;
      }
      float timeSinceSOF = getTimeSinceSOF(bc);

      std::bitset<64> ctpInputMask(bc.inputMask());
      if (ctpInputMask.test(2)) {
        histNBcsVsTime[runNumber]->Fill(timeSinceSOF);
        int runVdM23Start{542757};
        int runVdM23Stop{542768};
        if (runNumber < runVdM23Start || runNumber > runVdM23Stop) {
          rate = mRateFetcher.fetch(ccdb.service, bc.timestamp(), bc.runNumber(), std::string("T0VTX"), true) * 1.e-3; // kHz
        }
        histInteractionRate[runNumber]->Fill(rate);
      }

      int64_t globalBC = bc.globalBC();
      int64_t globalBCFDD = bc.globalBC() + bcShiftFDD;
      int localBC = globalBC % nBCsPerOrbit;
      int localBCFDD = globalBCFDD % nBCsPerOrbit;

      bool isSuperLeadingBcFDD{true}, isSuperLeadingBcFT0{true};
      if (globalBCFDD - globalBCIdOfLastBCWithActivityFDD < numEmptyBCsBeforeLeadingBC->get(0u, 2u)) {
        isSuperLeadingBcFDD = false; // not a super-leading BC for FDD
      }
      if (globalBC - globalBCIdOfLastBCWithActivityFT0 < numEmptyBCsBeforeLeadingBC->get(0u, 2u)) {
        isSuperLeadingBcFT0 = false; // not a super-leading BC for FT0
      }

      if (ctpInputMask.test(13) || ctpInputMask.test(15) || ctpInputMask.test(16) || ctpInputMask.test(17) || ctpInputMask.test(18)) { // 5 FDD triggers
        globalBCIdOfLastBCWithActivityFDD = globalBC;
      }
      if (ctpInputMask.test(1) || ctpInputMask.test(2) || ctpInputMask.test(3) || ctpInputMask.test(4) || ctpInputMask.test(5)) { // 5 FT0 triggers
        globalBCIdOfLastBCWithActivityFT0 = globalBCFDD;
      }

      if (!bcPatternB[localBC]) {
        isSuperLeadingBcFT0 = false; // not a super-leading BC
      }
      if (!bcPatternB[localBCFDD]) {
        isSuperLeadingBcFDD = false; // not a super-leading BC
      }

      int64_t globalBCStart = (globalBCLastInspectedBC >= 0 && globalBCLastInspectedBC < globalBC) ? globalBCLastInspectedBC + 1 : globalBC;
      int64_t maxBcDiff = (rate > 0) ? 10 * static_cast<int>(nBunchesFillingScheme * constants::lhc::LHCRevFreq / rate / 1.e3) : 1500;
      if (globalBC - globalBCStart > maxBcDiff) { // we changed fill, we should not count all BCs between the current and the previous one
        globalBCStart = globalBC;
      }
      for (int64_t iGlobalBC{globalBCStart}; iGlobalBC <= globalBC; ++iGlobalBC) { // we count all BCs in between one and another stored in the AO2Ds
        int iLocalBC = iGlobalBC % nBCsPerOrbit;
        if (bcPatternA[iLocalBC]) {
          nBCsPerBcId[iLocalBC][BCA]++;
        }
        if (bcPatternB[iLocalBC]) {
          nBCsPerBcId[iLocalBC][BCB]++;
          if (iGlobalBC - globalBCIdOfLastBCWithActivityFDD > numEmptyBCsBeforeLeadingBC->get(0u, 2u)) {
            nBCsPerBcId[iLocalBC][BCSLFDD]++;
          } else {
            nBCsPerBcId[iLocalBC][BCNSLFDD]++;
          }
          if (iGlobalBC - globalBCIdOfLastBCWithActivityFT0 > numEmptyBCsBeforeLeadingBC->get(0u, 2u)) {
            nBCsPerBcId[iLocalBC][BCSLFT0]++;
          } else {
            nBCsPerBcId[iLocalBC][BCNSLFT0]++;
          }
          if (!bcPatternL[iLocalBC]) {
            nBCsPerBcId[iLocalBC][BCNL]++;
          }
          if (!bcPatternLE[iLocalBC]) {
            nBCsPerBcId[iLocalBC][BCNLE]++;
          }
        }
        if (bcPatternC[iLocalBC]) {
          nBCsPerBcId[iLocalBC][BCC]++;
        }
        if (bcPatternE[iLocalBC]) {
          nBCsPerBcId[iLocalBC][BCE]++;
        }
        if (bcPatternL[iLocalBC]) {
          nBCsPerBcId[iLocalBC][BCL]++;
        }
        if (bcPatternLE[iLocalBC]) {
          nBCsPerBcId[iLocalBC][BCLE]++;
        }
      }

      int64_t thisTFid = (globalBC - bcSOR) / nBCsPerTF;

      if (thisTFid != currentTFid) {
        currentTFid = thisTFid;
        histTfPerMin[runNumber]->Fill(timeSinceSOF);
      }

      for (int iTrigger{0}; iTrigger < NTriggerAliases; ++iTrigger) {
        for (int iBCCategory{0}; iBCCategory < NBCCategories; ++iBCCategory) {
          if (doTypeBC->get(0u, iBCCategory)) {
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
              if (iBCCategory == BCLE && bcPatternLE[localBC])
                fillHistograms<AllBCs, BCLE>(timeSinceSOF, localBC);
              if (iBCCategory == BCSLFDD && isSuperLeadingBcFDD)
                fillHistograms<AllBCs, BCSLFDD>(timeSinceSOF, localBC);
              if (iBCCategory == BCSLFT0 && isSuperLeadingBcFT0)
                fillHistograms<AllBCs, BCSLFT0>(timeSinceSOF, localBC);
              if (iBCCategory == BCNL && !bcPatternL[localBC] && bcPatternB[localBC])
                fillHistograms<AllBCs, BCNL>(timeSinceSOF, localBC);
              if (iBCCategory == BCNLE && !bcPatternLE[localBC] && bcPatternB[localBC])
                fillHistograms<AllBCs, BCNLE>(timeSinceSOF, localBC);
              if (iBCCategory == BCNSLFDD && !isSuperLeadingBcFDD && bcPatternB[localBC])
                fillHistograms<AllBCs, BCNSLFDD>(timeSinceSOF, localBC);
              if (iBCCategory == BCNSLFT0 && !isSuperLeadingBcFT0 && bcPatternB[localBC])
                fillHistograms<AllBCs, BCNSLFT0>(timeSinceSOF, localBC);
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
              if (iBCCategory == BCLE && bcPatternLE[localBC])
                fillHistograms<FT0Vtx, BCLE>(timeSinceSOF, localBC);
              if (iBCCategory == BCSLFDD && isSuperLeadingBcFDD)
                fillHistograms<FT0Vtx, BCSLFDD>(timeSinceSOF, localBC);
              if (iBCCategory == BCSLFT0 && isSuperLeadingBcFT0)
                fillHistograms<FT0Vtx, BCSLFT0>(timeSinceSOF, localBC);
              if (iBCCategory == BCNL && !bcPatternL[localBC] && bcPatternB[localBC])
                fillHistograms<FT0Vtx, BCNL>(timeSinceSOF, localBC);
              if (iBCCategory == BCNLE && !bcPatternLE[localBC] && bcPatternB[localBC])
                fillHistograms<FT0Vtx, BCNLE>(timeSinceSOF, localBC);
              if (iBCCategory == BCNSLFDD && !isSuperLeadingBcFDD && bcPatternB[localBC])
                fillHistograms<FT0Vtx, BCNSLFDD>(timeSinceSOF, localBC);
              if (iBCCategory == BCNSLFT0 && !isSuperLeadingBcFT0 && bcPatternB[localBC])
                fillHistograms<FT0Vtx, BCNSLFT0>(timeSinceSOF, localBC);
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
              if (iBCCategory == BCLE && bcPatternLE[localBC])
                fillHistograms<FT0CE, BCLE>(timeSinceSOF, localBC);
              if (iBCCategory == BCSLFDD && isSuperLeadingBcFDD)
                fillHistograms<FT0CE, BCSLFDD>(timeSinceSOF, localBC);
              if (iBCCategory == BCSLFT0 && isSuperLeadingBcFT0)
                fillHistograms<FT0CE, BCSLFT0>(timeSinceSOF, localBC);
              if (iBCCategory == BCNL && !bcPatternL[localBC] && bcPatternB[localBC])
                fillHistograms<FT0CE, BCNL>(timeSinceSOF, localBC);
              if (iBCCategory == BCNLE && !bcPatternLE[localBC] && bcPatternB[localBC])
                fillHistograms<FT0CE, BCNLE>(timeSinceSOF, localBC);
              if (iBCCategory == BCNSLFDD && !isSuperLeadingBcFDD && bcPatternB[localBC])
                fillHistograms<FT0CE, BCNSLFDD>(timeSinceSOF, localBC);
              if (iBCCategory == BCNSLFT0 && !isSuperLeadingBcFT0 && bcPatternB[localBC])
                fillHistograms<FT0CE, BCNSLFT0>(timeSinceSOF, localBC);
            }
            if (iTrigger == FDD && ctpInputMask.test(15)) {
              if (iBCCategory == BCA && bcPatternA[localBCFDD])
                fillHistograms<FDD, BCA>(timeSinceSOF, localBCFDD);
              if (iBCCategory == BCB && bcPatternB[localBCFDD])
                fillHistograms<FDD, BCB>(timeSinceSOF, localBCFDD);
              if (iBCCategory == BCC && bcPatternC[localBCFDD])
                fillHistograms<FDD, BCC>(timeSinceSOF, localBCFDD);
              if (iBCCategory == BCE && bcPatternE[localBCFDD])
                fillHistograms<FDD, BCE>(timeSinceSOF, localBCFDD);
              if (iBCCategory == BCL && bcPatternL[localBCFDD])
                fillHistograms<FDD, BCL>(timeSinceSOF, localBCFDD);
              if (iBCCategory == BCLE && bcPatternLE[localBCFDD])
                fillHistograms<FDD, BCLE>(timeSinceSOF, localBCFDD);
              if (iBCCategory == BCSLFDD && isSuperLeadingBcFDD)
                fillHistograms<FDD, BCSLFDD>(timeSinceSOF, localBCFDD);
              if (iBCCategory == BCSLFT0 && isSuperLeadingBcFT0)
                fillHistograms<FDD, BCSLFT0>(timeSinceSOF, localBCFDD);
              if (iBCCategory == BCNL && !bcPatternL[localBCFDD] && bcPatternB[localBCFDD])
                fillHistograms<FDD, BCNL>(timeSinceSOF, localBCFDD);
              if (iBCCategory == BCNLE && !bcPatternLE[localBCFDD] && bcPatternB[localBCFDD])
                fillHistograms<FDD, BCNLE>(timeSinceSOF, localBCFDD);
              if (iBCCategory == BCNSLFDD && !isSuperLeadingBcFDD && bcPatternB[localBCFDD])
                fillHistograms<FDD, BCNSLFDD>(timeSinceSOF, localBCFDD);
              if (iBCCategory == BCNSLFT0 && !isSuperLeadingBcFT0 && bcPatternB[localBCFDD])
                fillHistograms<FDD, BCNSLFT0>(timeSinceSOF, localBCFDD);
            }
          }
        }
      }
      histNBcsVsBcId[runNumber]->Fill(localBC);
      if (globalBCLastInspectedBC < globalBC) {
        globalBCLastInspectedBC = globalBC;
      } else {
        globalBCLastInspectedBC = -1;
      }
    }
    // fill histogram for mu
    for (int iTrigger{0}; iTrigger < NTriggerAliases; ++iTrigger) {
      for (int iBCCategory{0}; iBCCategory < NBCCategories; ++iBCCategory) {
        if (doTypeBC->get(0u, iBCCategory)) {
          for (int iBcId{0}; iBcId < nBCsPerOrbit; ++iBcId) {
            histBcInspectVsBcId[iTrigger][iBCCategory][runNumber]->Fill(iBcId, nBCsPerBcId[iBcId][iBCCategory]);
          }
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  metadataInfo.initMetadata(cfgc);
  return WorkflowSpec{adaptAnalysisTask<LumiStabilityPP>(cfgc)};
}
